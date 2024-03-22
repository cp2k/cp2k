!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      11.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_util
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_indxg2l,&
                                              cp_fm_indxg2p,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE rpa_communication,               ONLY: communicate_buffer
   USE util,                            ONLY: sort,&
                                              sort_unique
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_util'

   PUBLIC :: mult_B_with_W, fm_general_add_bse, truncate_fm, fm_write_thresh, print_BSE_start_flag, &
             deallocate_matrices_bse, comp_eigvec_coeff_BSE, sort_excitations, &
             estimate_BSE_resources, filter_eigvec_contrib, truncate_BSE_matrices

CONTAINS

! **************************************************************************************************
!> \brief Multiplies B-matrix (RI-3c-Integrals) with W (screening) to obtain \bar{B}
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE mult_B_with_W(fm_mat_S_ij_bse, fm_mat_S_ia_bse, fm_mat_S_bar_ia_bse, &
                            fm_mat_S_bar_ij_bse, fm_mat_Q_static_bse_gemm, &
                            dimen_RI, homo, virtual)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ij_bse, fm_mat_S_ia_bse
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q_static_bse_gemm
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_B_with_W'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_S_bar_ia_bse, fm_mat_S_ia_bse%matrix_struct)
      CALL cp_fm_set_all(fm_mat_S_bar_ia_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_S_bar_ij_bse, fm_mat_S_ij_bse%matrix_struct)
      CALL cp_fm_set_all(fm_mat_S_bar_ij_bse, 0.0_dp)

      CALL cp_fm_create(fm_work, fm_mat_Q_static_bse_gemm%matrix_struct)
      CALL cp_fm_set_all(fm_work, 0.0_dp)

      ! get info of fm_mat_Q_static_bse and compute ((1+Q(0))^-1-1)
      CALL cp_fm_get_info(matrix=fm_mat_Q_static_bse_gemm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) = fm_mat_Q_static_bse_gemm%local_data(iiB, jjB) + 1.0_dp
            END IF
         END DO
      END DO

      ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q_static_bse_gemm, n=dimen_RI, info_out=info_chol)

      CPASSERT(info_chol == 0)

      ! calculate [1+Q(i0)]^-1
      CALL cp_fm_cholesky_invert(fm_mat_Q_static_bse_gemm)

      ! symmetrize the result
      CALL cp_fm_upper_to_full(fm_mat_Q_static_bse_gemm, fm_work)

      CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=homo**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_Q_static_bse_gemm, matrix_b=fm_mat_S_ij_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_S_bar_ij_bse)

      ! fm_mat_S_bar_ia_bse has a different blacs_env as fm_mat_S_ij_bse since we take
      ! fm_mat_S_ia_bse from RPA. Therefore, we also need a different fm_mat_Q_static_bse_gemm
      CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=homo*virtual, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_Q_static_bse_gemm, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_S_bar_ia_bse)

      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Adds and reorders full matrices with a combined index structure, e.g. adding W_ij,ab
!> to A_ia, which needs MPI communication.
!> \param fm_out ...
!> \param fm_in ...
!> \param beta ...
!> \param nrow_secidx_in ...
!> \param ncol_secidx_in ...
!> \param nrow_secidx_out ...
!> \param ncol_secidx_out ...
!> \param unit_nr ...
!> \param reordering ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE fm_general_add_bse(fm_out, fm_in, beta, nrow_secidx_in, ncol_secidx_in, &
                                 nrow_secidx_out, ncol_secidx_out, unit_nr, reordering, mp2_env)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_out, fm_in
      REAL(kind=dp)                                      :: beta
      INTEGER, INTENT(IN)                                :: nrow_secidx_in, ncol_secidx_in, &
                                                            nrow_secidx_out, ncol_secidx_out
      INTEGER                                            :: unit_nr
      INTEGER, DIMENSION(4)                              :: reordering
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_general_add_bse'

      INTEGER :: col_idx_loc, dummy, handle, handle2, i_entry_rec, idx_col_out, idx_row_out, ii, &
         iproc, jj, ncol_block_in, ncol_block_out, ncol_local_in, ncol_local_out, npcol, nprocs, &
         nprow, nrow_block_in, nrow_block_out, nrow_local_in, nrow_local_out, proc_send, &
         row_idx_loc, send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(4)                              :: indices_in
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_in, col_indices_out, &
                                                            row_indices_in, row_indices_out
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env_out
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

      CALL timeset(routineN, handle)
      CALL timeset(routineN//"_1_setup", handle2)

      para_env_out => fm_out%matrix_struct%para_env
      ! A_iajb
      ! We start by moving data from local parts of W_ijab to the full matrix A_iajb using buffers
      CALL cp_fm_get_info(matrix=fm_out, &
                          nrow_local=nrow_local_out, &
                          ncol_local=ncol_local_out, &
                          row_indices=row_indices_out, &
                          col_indices=col_indices_out, &
                          nrow_block=nrow_block_out, &
                          ncol_block=ncol_block_out)

      nprow = fm_out%matrix_struct%context%num_pe(1)
      npcol = fm_out%matrix_struct%context%num_pe(2)

      ALLOCATE (num_entries_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env_out%num_pe - 1))

      num_entries_rec(:) = 0
      num_entries_send(:) = 0

      dummy = 0

      CALL cp_fm_get_info(matrix=fm_in, &
                          nrow_local=nrow_local_in, &
                          ncol_local=ncol_local_in, &
                          row_indices=row_indices_in, &
                          col_indices=col_indices_in, &
                          nrow_block=nrow_block_in, &
                          ncol_block=ncol_block_in)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A14,A10,T71,I10)') 'BSE|DEBUG|', 'Row number of ', fm_out%name, &
            fm_out%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A10,T13,A17,A10,T71,I10)') 'BSE|DEBUG|', 'Column number of ', fm_out%name, &
            fm_out%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A10,T13,A18,A10,T71,I10)') 'BSE|DEBUG|', 'Row block size of ', fm_out%name, nrow_block_out
         WRITE (unit_nr, '(T2,A10,T13,A21,A10,T71,I10)') 'BSE|DEBUG|', 'Column block size of ', fm_out%name, ncol_block_out

         WRITE (unit_nr, '(T2,A10,T13,A14,A10,T71,I10)') 'BSE|DEBUG|', 'Row number of ', fm_in%name, &
            fm_in%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A10,T13,A17,A10,T71,I10)') 'BSE|DEBUG|', 'Column number of ', fm_in%name, &
            fm_in%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A10,T13,A18,A10,T71,I10)') 'BSE|DEBUG|', 'Row block size of ', fm_in%name, nrow_block_in
         WRITE (unit_nr, '(T2,A10,T13,A21,A10,T71,I10)') 'BSE|DEBUG|', 'Column block size of ', fm_in%name, ncol_block_in
      END IF

      ! Use scalapack wrapper to find process index in fm_out
      ! To that end, we obtain the global index in fm_out from the level indices
      indices_in(:) = 0
      DO row_idx_loc = 1, nrow_local_in
         indices_in(1) = (row_indices_in(row_idx_loc) - 1)/nrow_secidx_in + 1
         indices_in(2) = MOD(row_indices_in(row_idx_loc) - 1, nrow_secidx_in) + 1
         DO col_idx_loc = 1, ncol_local_in
            indices_in(3) = (col_indices_in(col_idx_loc) - 1)/ncol_secidx_in + 1
            indices_in(4) = MOD(col_indices_in(col_idx_loc) - 1, ncol_secidx_in) + 1

            idx_row_out = indices_in(reordering(2)) + (indices_in(reordering(1)) - 1)*nrow_secidx_out
            idx_col_out = indices_in(reordering(4)) + (indices_in(reordering(3)) - 1)*ncol_secidx_out

            send_prow = cp_fm_indxg2p(idx_row_out, nrow_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A27)') 'BSE|DEBUG|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A18)') 'BSE|DEBUG|', 'Allocating buffers'
      END IF

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_out%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
         buffer_rec(iproc)%indx = 0

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%indx(num_entries_send(iproc), 2))
         buffer_send(iproc)%indx = 0

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_4_buf_from_fmin_"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A18,A10,A13)') 'BSE|DEBUG|', 'Writing data from ', fm_in%name, ' into buffers'
      END IF

      ALLOCATE (entry_counter(0:para_env_out%num_pe - 1))
      entry_counter(:) = 0

      ! Now we can write the actual data and indices to the send-buffer
      DO row_idx_loc = 1, nrow_local_in
         indices_in(1) = (row_indices_in(row_idx_loc) - 1)/nrow_secidx_in + 1
         indices_in(2) = MOD(row_indices_in(row_idx_loc) - 1, nrow_secidx_in) + 1
         DO col_idx_loc = 1, ncol_local_in
            indices_in(3) = (col_indices_in(col_idx_loc) - 1)/ncol_secidx_in + 1
            indices_in(4) = MOD(col_indices_in(col_idx_loc) - 1, ncol_secidx_in) + 1

            idx_row_out = indices_in(reordering(2)) + (indices_in(reordering(1)) - 1)*nrow_secidx_out
            idx_col_out = indices_in(reordering(4)) + (indices_in(reordering(3)) - 1)*ncol_secidx_out

            send_prow = cp_fm_indxg2p(idx_row_out, nrow_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               fm_in%local_data(row_idx_loc, col_idx_loc)

            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = idx_row_out
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = idx_col_out

         END DO
      END DO

      ALLOCATE (req_array(1:para_env_out%num_pe, 4))

      CALL timestop(handle2)

      CALL timeset(routineN//"_5_comm_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A21)') 'BSE|DEBUG|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A24,A10)') 'BSE|DEBUG|', 'Writing from buffers to ', fm_out%name
      END IF

      ! fill fm_out with the entries from buffer_rec, i.e. buffer_rec are parts of fm_in
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_out, nprocs, nrow_block_out, ncol_block_out, &
!$OMP          num_entries_rec, buffer_rec, beta, dummy, nprow, npcol) &
!$OMP PRIVATE(iproc, i_entry_rec, ii, jj)
      DO iproc = 0, nprocs - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            ii = cp_fm_indxg2l(buffer_rec(iproc)%indx(i_entry_rec, 1), nrow_block_out, &
                               dummy, dummy, nprow)
            jj = cp_fm_indxg2l(buffer_rec(iproc)%indx(i_entry_rec, 2), ncol_block_out, &
                               dummy, dummy, npcol)

            fm_out%local_data(ii, jj) = fm_out%local_data(ii, jj) + beta*buffer_rec(iproc)%msg(i_entry_rec)
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A41)') 'BSE|DEBUG|', 'Starting cleanup of communication buffers'
      END IF

      !Clean up all the arrays from the communication process
      DO iproc = 0, para_env_out%num_pe - 1
         DEALLOCATE (buffer_rec(iproc)%msg)
         DEALLOCATE (buffer_rec(iproc)%indx)
         DEALLOCATE (buffer_send(iproc)%msg)
         DEALLOCATE (buffer_send(iproc)%indx)
      END DO
      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (req_array)
      DEALLOCATE (entry_counter)
      DEALLOCATE (num_entries_rec, num_entries_send)

      CALL timestop(handle2)
      CALL timestop(handle)

   END SUBROUTINE fm_general_add_bse

! **************************************************************************************************
!> \brief Routine for truncating a full matrix as given by the energy cutoffs in the input file.
!>  Logic: Matrices have some dimension dimen_RI x nrow_in*ncol_in  for the incoming (untruncated) matrix
!>  and dimen_RI x nrow_out*ncol_out for the truncated matrix. The truncation is done by resorting the indices
!>  via parallel communication.
!> \param fm_out ...
!> \param fm_in ...
!> \param ncol_in ...
!> \param nrow_out ...
!> \param ncol_out ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param nrow_offset ...
!> \param ncol_offset ...
! **************************************************************************************************
   SUBROUTINE truncate_fm(fm_out, fm_in, ncol_in, &
                          nrow_out, ncol_out, unit_nr, mp2_env, &
                          nrow_offset, ncol_offset)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_out
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_in
      INTEGER                                            :: ncol_in, nrow_out, ncol_out, unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nrow_offset, ncol_offset

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'truncate_fm'

      INTEGER :: col_idx_loc, dummy, handle, handle2, i_entry_rec, idx_col_first, idx_col_in, &
         idx_col_out, idx_col_sec, idx_row_in, ii, iproc, jj, ncol_block_in, ncol_block_out, &
         ncol_local_in, ncol_local_out, npcol, nprocs, nprow, nrow_block_in, nrow_block_out, &
         nrow_local_in, nrow_local_out, proc_send, row_idx_loc, send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_in, col_indices_out, &
                                                            row_indices_in, row_indices_out
      LOGICAL                                            :: correct_ncol, correct_nrow
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env_out
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

      CALL timeset(routineN, handle)
      CALL timeset(routineN//"_1_setup", handle2)

      correct_nrow = .FALSE.
      correct_ncol = .FALSE.
      !In case of truncation in the occupied space, we need to correct the interval of indices
      IF (PRESENT(nrow_offset)) THEN
         correct_nrow = .TRUE.
      END IF
      IF (PRESENT(ncol_offset)) THEN
         correct_ncol = .TRUE.
      END IF

      para_env_out => fm_out%matrix_struct%para_env

      CALL cp_fm_get_info(matrix=fm_out, &
                          nrow_local=nrow_local_out, &
                          ncol_local=ncol_local_out, &
                          row_indices=row_indices_out, &
                          col_indices=col_indices_out, &
                          nrow_block=nrow_block_out, &
                          ncol_block=ncol_block_out)

      nprow = fm_out%matrix_struct%context%num_pe(1)
      npcol = fm_out%matrix_struct%context%num_pe(2)

      ALLOCATE (num_entries_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env_out%num_pe - 1))

      num_entries_rec(:) = 0
      num_entries_send(:) = 0

      dummy = 0

      CALL cp_fm_get_info(matrix=fm_in, &
                          nrow_local=nrow_local_in, &
                          ncol_local=ncol_local_in, &
                          row_indices=row_indices_in, &
                          col_indices=col_indices_in, &
                          nrow_block=nrow_block_in, &
                          ncol_block=ncol_block_in)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A14,A10,T71,I10)') 'BSE|DEBUG|', 'Row number of ', fm_out%name, &
            fm_out%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A10,T13,A17,A10,T71,I10)') 'BSE|DEBUG|', 'Column number of ', fm_out%name, &
            fm_out%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A10,T13,A18,A10,T71,I10)') 'BSE|DEBUG|', 'Row block size of ', fm_out%name, nrow_block_out
         WRITE (unit_nr, '(T2,A10,T13,A21,A10,T71,I10)') 'BSE|DEBUG|', 'Column block size of ', fm_out%name, ncol_block_out

         WRITE (unit_nr, '(T2,A10,T13,A14,A10,T71,I10)') 'BSE|DEBUG|', 'Row number of ', fm_in%name, &
            fm_in%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A10,T13,A17,A10,T71,I10)') 'BSE|DEBUG|', 'Column number of ', fm_in%name, &
            fm_in%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A10,T13,A18,A10,T71,I10)') 'BSE|DEBUG|', 'Row block size of ', fm_in%name, nrow_block_in
         WRITE (unit_nr, '(T2,A10,T13,A21,A10,T71,I10)') 'BSE|DEBUG|', 'Column block size of ', fm_in%name, ncol_block_in
      END IF

      ! We find global indices in S with nrow_in and ncol_in for truncation
      DO col_idx_loc = 1, ncol_local_in
         idx_col_in = col_indices_in(col_idx_loc)

         idx_col_first = (idx_col_in - 1)/ncol_in + 1
         idx_col_sec = MOD(idx_col_in - 1, ncol_in) + 1

         ! If occupied orbitals are included, these have to be handled differently
         ! due to their reversed indexing
         IF (correct_nrow) THEN
            idx_col_first = idx_col_first - nrow_offset + 1
            IF (idx_col_first .LE. 0) CYCLE
         ELSE
            IF (idx_col_first > nrow_out) EXIT
         END IF
         IF (correct_ncol) THEN
            idx_col_sec = idx_col_sec - ncol_offset + 1
            IF (idx_col_sec .LE. 0) CYCLE
         ELSE
            IF (idx_col_sec > ncol_out) CYCLE
         END IF

         idx_col_out = idx_col_sec + (idx_col_first - 1)*ncol_out

         DO row_idx_loc = 1, nrow_local_in
            idx_row_in = row_indices_in(row_idx_loc)

            send_prow = cp_fm_indxg2p(idx_row_in, nrow_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A27)') 'BSE|DEBUG|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A18)') 'BSE|DEBUG|', 'Allocating buffers'
      END IF

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_out%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
         buffer_rec(iproc)%indx = 0

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%indx(num_entries_send(iproc), 2))
         buffer_send(iproc)%indx = 0

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_4_buf_from_fmin_"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A18,A10,A13)') 'BSE|DEBUG|', 'Writing data from ', fm_in%name, ' into buffers'
      END IF

      ALLOCATE (entry_counter(0:para_env_out%num_pe - 1))
      entry_counter(:) = 0

      ! Now we can write the actual data and indices to the send-buffer
      DO col_idx_loc = 1, ncol_local_in
         idx_col_in = col_indices_in(col_idx_loc)

         idx_col_first = (idx_col_in - 1)/ncol_in + 1
         idx_col_sec = MOD(idx_col_in - 1, ncol_in) + 1

         ! If occupied orbitals are included, these have to be handled differently
         ! due to their reversed indexing
         IF (correct_nrow) THEN
            idx_col_first = idx_col_first - nrow_offset + 1
            IF (idx_col_first .LE. 0) CYCLE
         ELSE
            IF (idx_col_first > nrow_out) EXIT
         END IF
         IF (correct_ncol) THEN
            idx_col_sec = idx_col_sec - ncol_offset + 1
            IF (idx_col_sec .LE. 0) CYCLE
         ELSE
            IF (idx_col_sec > ncol_out) CYCLE
         END IF

         idx_col_out = idx_col_sec + (idx_col_first - 1)*ncol_out

         DO row_idx_loc = 1, nrow_local_in
            idx_row_in = row_indices_in(row_idx_loc)

            send_prow = cp_fm_indxg2p(idx_row_in, nrow_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy, &
                                      fm_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               fm_in%local_data(row_idx_loc, col_idx_loc)
            !No need to create row_out, since it is identical to incoming
            !We dont change the RI index for any fm_mat_XX_BSE
            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = idx_row_in
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = idx_col_out

         END DO
      END DO

      ALLOCATE (req_array(1:para_env_out%num_pe, 4))

      CALL timestop(handle2)

      CALL timeset(routineN//"_5_comm_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A21)') 'BSE|DEBUG|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A24,A10)') 'BSE|DEBUG|', 'Writing from buffers to ', fm_out%name
      END IF

      ! fill fm_out with the entries from buffer_rec, i.e. buffer_rec are parts of fm_in
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_out, nprocs, nrow_block_out, ncol_block_out, &
!$OMP          num_entries_rec, buffer_rec, nprow, npcol, dummy) &
!$OMP PRIVATE(iproc, i_entry_rec, ii, jj)
      DO iproc = 0, nprocs - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            ii = cp_fm_indxg2l(buffer_rec(iproc)%indx(i_entry_rec, 1), nrow_block_out, &
                               dummy, dummy, nprow)
            jj = cp_fm_indxg2l(buffer_rec(iproc)%indx(i_entry_rec, 2), ncol_block_out, &
                               dummy, dummy, npcol)

            fm_out%local_data(ii, jj) = fm_out%local_data(ii, jj) + buffer_rec(iproc)%msg(i_entry_rec)
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A41)') 'BSE|DEBUG|', 'Starting cleanup of communication buffers'
      END IF

      !Clean up all the arrays from the communication process
      DO iproc = 0, para_env_out%num_pe - 1
         DEALLOCATE (buffer_rec(iproc)%msg)
         DEALLOCATE (buffer_rec(iproc)%indx)
         DEALLOCATE (buffer_send(iproc)%msg)
         DEALLOCATE (buffer_send(iproc)%indx)
      END DO
      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (req_array)
      DEALLOCATE (entry_counter)
      DEALLOCATE (num_entries_rec, num_entries_send)

      CALL timestop(handle2)
      CALL timestop(handle)

   END SUBROUTINE truncate_fm

! **************************************************************************************************
!> \brief Debug function to write elements of a full matrix to file, if they are larger than a given threshold
!> \param fm ...
!> \param thresh ...
!> \param header ...
!> \param unit_nr ...
!> \param abs_vals ...
! **************************************************************************************************
   SUBROUTINE fm_write_thresh(fm, thresh, header, unit_nr, abs_vals)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm
      REAL(KIND=dp), INTENT(IN)                          :: thresh
      CHARACTER(LEN=*), INTENT(IN)                       :: header
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, OPTIONAL                                  :: abs_vals

      CHARACTER(LEN=*), PARAMETER :: my_footer = " | ENDING WRITING OF MATRIX", &
         routineN = 'fm_write_thresh'

      INTEGER                                            :: handle, i, j, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_abs_vals

      CALL timeset(routineN, handle)

      IF (PRESENT(abs_vals)) THEN
         my_abs_vals = abs_vals
      ELSE
         my_abs_vals = .FALSE.
      END IF

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) header
      END IF
      IF (my_abs_vals) THEN
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm%local_data(i, j)) > thresh) THEN
                  WRITE (unit_nr, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     ABS(fm%local_data(i, j))
               END IF
            END DO
         END DO
      ELSE
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm%local_data(i, j)) > thresh) THEN
                  WRITE (unit_nr, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     fm%local_data(i, j)
               END IF
            END DO
         END DO
      END IF
      CALL fm%matrix_struct%para_env%sync()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) my_footer
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param bse_tda ...
!> \param bse_abba ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_BSE_start_flag(bse_tda, bse_abba, unit_nr)

      LOGICAL, INTENT(IN)                                :: bse_tda, bse_abba
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_BSE_start_flag'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '**           Bethe Salpeter equation (BSE) for excitation energies           **'
         IF (bse_tda .AND. bse_abba) THEN
            WRITE (unit_nr, '(T2,A79)') '**          solved with and without Tamm-Dancoff approximation (TDA)         **'
         ELSE IF (bse_tda) THEN
            WRITE (unit_nr, '(T2,A79)') '**                solved with Tamm-Dancoff approximation (TDA)               **'
         ELSE
            WRITE (unit_nr, '(T2,A79)') '**               solved without Tamm-Dancoff approximation (TDA)             **'
         END IF

         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, *) ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_trunc ...
!> \param fm_mat_S_ij_trunc ...
!> \param fm_mat_S_ab_trunc ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_bse(fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                      fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                      fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)

      TYPE(cp_fm_type), INTENT(INOUT) :: fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_trunc, &
         fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_bse'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_release(fm_mat_S_bar_ia_bse)
      CALL cp_fm_release(fm_mat_S_bar_ij_bse)
      CALL cp_fm_release(fm_mat_S_trunc)
      CALL cp_fm_release(fm_mat_S_ij_trunc)
      CALL cp_fm_release(fm_mat_S_ab_trunc)
      CALL cp_fm_release(fm_mat_Q_static_bse)
      CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

      CALL timestop(handle)

   END SUBROUTINE deallocate_matrices_bse

! **************************************************************************************************
!> \brief Routine for computing the coefficients of the eigenvectors of the BSE matrix from a
!>  multiplication with the eigenvalues
!> \param fm_work ...
!> \param eig_vals ...
!> \param beta ...
!> \param gamma ...
!> \param do_transpose ...
! **************************************************************************************************
   SUBROUTINE comp_eigvec_coeff_BSE(fm_work, eig_vals, beta, gamma, do_transpose)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: eig_vals
      REAL(KIND=dp), INTENT(IN)                          :: beta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: gamma
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_transpose

      CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_eigvec_coeff_BSE'

      INTEGER                                            :: handle, i_row_global, ii, j_col_global, &
                                                            jj, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_transpose
      REAL(KIND=dp)                                      :: coeff, my_gamma

      CALL timeset(routineN, handle)

      IF (PRESENT(gamma)) THEN
         my_gamma = gamma
      ELSE
         my_gamma = 2.0_dp
      END IF

      IF (PRESENT(do_transpose)) THEN
         my_do_transpose = do_transpose
      ELSE
         my_do_transpose = .FALSE.
      END IF

      CALL cp_fm_get_info(matrix=fm_work, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF (my_do_transpose) THEN
         DO jj = 1, ncol_local
            j_col_global = col_indices(jj)
            DO ii = 1, nrow_local
               coeff = (eig_vals(j_col_global)**beta)/my_gamma
               fm_work%local_data(ii, jj) = fm_work%local_data(ii, jj)*coeff
            END DO
         END DO
      ELSE
         DO jj = 1, ncol_local
            DO ii = 1, nrow_local
               i_row_global = row_indices(ii)
               coeff = (eig_vals(i_row_global)**beta)/my_gamma
               fm_work%local_data(ii, jj) = fm_work%local_data(ii, jj)*coeff
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param idx_prim ...
!> \param idx_sec ...
!> \param eigvec_entries ...
! **************************************************************************************************
   SUBROUTINE sort_excitations(idx_prim, idx_sec, eigvec_entries)

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_prim, idx_sec
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sort_excitations'

      INTEGER                                            :: handle, ii, kk, num_entries, num_mults
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_prim_work, idx_sec_work, tmp_index
      LOGICAL                                            :: unique_entries
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries_work

      CALL timeset(routineN, handle)

      num_entries = SIZE(idx_prim)

      ALLOCATE (tmp_index(num_entries))

      CALL sort(idx_prim, num_entries, tmp_index)

      ALLOCATE (idx_sec_work(num_entries))
      ALLOCATE (eigvec_entries_work(num_entries))

      DO ii = 1, num_entries
         idx_sec_work(ii) = idx_sec(tmp_index(ii))
         eigvec_entries_work(ii) = eigvec_entries(tmp_index(ii))
      END DO

      DEALLOCATE (tmp_index)
      DEALLOCATE (idx_sec)
      DEALLOCATE (eigvec_entries)

      CALL MOVE_ALLOC(idx_sec_work, idx_sec)
      CALL MOVE_ALLOC(eigvec_entries_work, eigvec_entries)

      !Now check for multiple entries in first idx to check necessity of sorting in second idx
      CALL sort_unique(idx_prim, unique_entries)
      IF (.NOT. unique_entries) THEN
         ALLOCATE (idx_prim_work(num_entries))
         idx_prim_work(:) = idx_prim(:)
         ! Find duplicate entries in idx_prim
         DO ii = 1, num_entries
            IF (idx_prim_work(ii) == 0) CYCLE
            num_mults = COUNT(idx_prim_work == idx_prim_work(ii))
            IF (num_mults > 1) THEN
               !Set all duplicate entries to 0
               idx_prim_work(ii:ii + num_mults - 1) = 0
               !Start sorting in secondary index
               ALLOCATE (idx_sec_work(num_mults))
               ALLOCATE (eigvec_entries_work(num_mults))
               idx_sec_work(:) = idx_sec(ii:ii + num_mults - 1)
               eigvec_entries_work(:) = eigvec_entries(ii:ii + num_mults - 1)
               ALLOCATE (tmp_index(num_mults))
               CALL sort(idx_sec_work, num_mults, tmp_index)

               !Now write newly sorted indices to original arrays
               DO kk = ii, ii + num_mults - 1
                  idx_sec(kk) = idx_sec_work(kk - ii + 1)
                  eigvec_entries(kk) = eigvec_entries_work(tmp_index(kk - ii + 1))
               END DO
               !Deallocate work arrays
               DEALLOCATE (tmp_index)
               DEALLOCATE (idx_sec_work)
               DEALLOCATE (eigvec_entries_work)
            END IF
            idx_prim_work(ii) = idx_prim(ii)
         END DO
         DEALLOCATE (idx_prim_work)
      END IF

      CALL timestop(handle)

   END SUBROUTINE sort_excitations

! **************************************************************************************************
!> \brief Roughly estimates the needed runtime and memory during the BSE run
!> \param homo_red ...
!> \param virtual_red ...
!> \param unit_nr ...
!> \param bse_abba ...
!> \param para_env ...
!> \param diag_runtime_est ...
! **************************************************************************************************
   SUBROUTINE estimate_BSE_resources(homo_red, virtual_red, unit_nr, bse_abba, &
                                     para_env, diag_runtime_est)

      INTEGER                                            :: homo_red, virtual_red, unit_nr
      LOGICAL                                            :: bse_abba
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp)                                      :: diag_runtime_est

      CHARACTER(LEN=*), PARAMETER :: routineN = 'estimate_BSE_resources'

      INTEGER                                            :: handle, num_BSE_matrices
      REAL(KIND=dp)                                      :: mem_est

      CALL timeset(routineN, handle)

      ! Number of matrices with size of A in TDA is 2 (A itself and W_ijab)
      num_BSE_matrices = 2
      ! With the full diagonalization of ABBA, we need several auxiliary matrices in the process
      ! The maximum number is 2 + 2 + 6 (additional B and C matrix as well as 6 matrices to create C)
      IF (bse_abba) THEN
         num_BSE_matrices = 10
      END IF
      mem_est = REAL(REAL(8*REAL(homo_red**2*virtual_red**2)*num_BSE_matrices)/(1024**3))
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A40,T68,F13.3)') 'BSE|', 'Total peak memory estimate from BSE [GB]', &
            mem_est
         WRITE (unit_nr, '(T2,A4,T7,A47,T68,F13.3)') 'BSE|', 'Peak memory estimate per MPI rank from BSE [GB]', &
            REAL(mem_est/para_env%num_pe)
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF
      ! Rough estimation of diagonalization runtimes. Baseline was a full BSE Naphthalene
      ! run with 11000x11000 entries in A/B/C, which took 10s on 32 ranks
      diag_runtime_est = REAL(REAL(homo_red*virtual_red/11000)**3*10*32/para_env%num_pe)

      CALL timestop(handle)

   END SUBROUTINE estimate_BSE_resources

! **************************************************************************************************
!> \brief Filters eigenvector entries above a given threshold to describe excitations in the
!> singleparticle basis
!> \param fm_eigvec ...
!> \param idx_homo ...
!> \param idx_virt ...
!> \param eigvec_entries ...
!> \param i_exc ...
!> \param virtual ...
!> \param num_entries ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE filter_eigvec_contrib(fm_eigvec, idx_homo, idx_virt, eigvec_entries, &
                                    i_exc, virtual, num_entries, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_homo, idx_virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries
      INTEGER                                            :: i_exc, virtual, num_entries
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'filter_eigvec_contrib'

      INTEGER                                            :: eigvec_idx, handle, ii, iproc, jj, kk, &
                                                            ncol_local, nrow_local, &
                                                            num_entries_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_entries_to_comm
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigvec_entry
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_entries
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      para_env => fm_eigvec%matrix_struct%para_env

      CALL cp_fm_get_info(matrix=fm_eigvec, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ALLOCATE (num_entries_to_comm(0:para_env%num_pe - 1))
      num_entries_to_comm(:) = 0

      DO jj = 1, ncol_local
         !First check if i is localized on this proc
         IF (col_indices(jj) /= i_exc) THEN
            CYCLE
         END IF
         DO ii = 1, nrow_local
            eigvec_idx = row_indices(ii)
            eigvec_entry = fm_eigvec%local_data(ii, jj)/SQRT(2.0_dp)
            IF (ABS(eigvec_entry) > mp2_env%ri_g0w0%eps_x) THEN
               num_entries_to_comm(para_env%mepos) = num_entries_to_comm(para_env%mepos) + 1
            END IF
         END DO
      END DO

      !Gather number of entries of other processes
      CALL para_env%sum(num_entries_to_comm)

      num_entries_local = num_entries_to_comm(para_env%mepos)

      ALLOCATE (buffer_entries(0:para_env%num_pe - 1))

      DO iproc = 0, para_env%num_pe - 1
         ALLOCATE (buffer_entries(iproc)%msg(num_entries_to_comm(iproc)))
         ALLOCATE (buffer_entries(iproc)%indx(num_entries_to_comm(iproc), 2))
         buffer_entries(iproc)%msg = 0.0_dp
         buffer_entries(iproc)%indx = 0
      END DO

      kk = 1
      DO jj = 1, ncol_local
         !First check if i is localized on this proc
         IF (col_indices(jj) /= i_exc) THEN
            CYCLE
         END IF
         DO ii = 1, nrow_local
            eigvec_idx = row_indices(ii)
            eigvec_entry = fm_eigvec%local_data(ii, jj)/SQRT(2.0_dp)
            IF (ABS(eigvec_entry) > mp2_env%ri_g0w0%eps_x) THEN
               buffer_entries(para_env%mepos)%indx(kk, 1) = (eigvec_idx - 1)/virtual + 1
               buffer_entries(para_env%mepos)%indx(kk, 2) = MOD(eigvec_idx - 1, virtual) + 1
               buffer_entries(para_env%mepos)%msg(kk) = eigvec_entry
               kk = kk + 1
            END IF
         END DO
      END DO

      DO iproc = 0, para_env%num_pe - 1
         CALL para_env%sum(buffer_entries(iproc)%msg)
         CALL para_env%sum(buffer_entries(iproc)%indx)
      END DO

      !Now sum up gathered information
      num_entries = SUM(num_entries_to_comm)
      ALLOCATE (idx_homo(num_entries))
      ALLOCATE (idx_virt(num_entries))
      ALLOCATE (eigvec_entries(num_entries))

      kk = 1
      DO iproc = 0, para_env%num_pe - 1
         IF (num_entries_to_comm(iproc) /= 0) THEN
            DO ii = 1, num_entries_to_comm(iproc)
               idx_homo(kk) = buffer_entries(iproc)%indx(ii, 1)
               idx_virt(kk) = buffer_entries(iproc)%indx(ii, 2)
               eigvec_entries(kk) = buffer_entries(iproc)%msg(ii)
               kk = kk + 1
            END DO
         END IF
      END DO

      !Deallocate all the used arrays
      DO iproc = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_entries(iproc)%msg)
         DEALLOCATE (buffer_entries(iproc)%indx)
      END DO
      DEALLOCATE (buffer_entries)
      DEALLOCATE (num_entries_to_comm)
      NULLIFY (row_indices)
      NULLIFY (col_indices)

      !Now sort the results according to the involved singleparticle orbitals
      ! (homo first, then virtual)
      CALL sort_excitations(idx_homo, idx_virt, eigvec_entries)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Determines indices within the given energy cutoffs and truncates Eigenvalues and matrices
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S_trunc ...
!> \param fm_mat_S_ij_trunc ...
!> \param fm_mat_S_ab_trunc ...
!> \param Eigenval ...
!> \param Eigenval_reduced ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param homo_red ...
!> \param virt_red ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE truncate_BSE_matrices(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                    Eigenval, Eigenval_reduced, &
                                    homo, virtual, dimen_RI, unit_nr, &
                                    homo_red, virt_red, &
                                    mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ia_bse, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_S_trunc, fm_mat_S_ij_trunc, &
                                                            fm_mat_S_ab_trunc
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, unit_nr
      INTEGER, INTENT(OUT)                               :: homo_red, virt_red
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'truncate_BSE_matrices'

      INTEGER                                            :: handle, homo_incl, i_homo, j_virt, &
                                                            virt_incl
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ab, fm_struct_ia, fm_struct_ij
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Determine index in homo and virtual for truncation
      ! Uses indices of outermost orbitals within energy range (-mp2_env%ri_g0w0%bse_cutoff_occ,mp2_env%ri_g0w0%bse_cutoff_virt)
      IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0 .OR. mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
         IF (-mp2_env%ri_g0w0%bse_cutoff_occ .LT. Eigenval(1) - Eigenval(homo) &
             .OR. mp2_env%ri_g0w0%bse_cutoff_occ < 0) THEN
            homo_red = homo
            homo_incl = 1
         ELSE
            homo_incl = 1
            DO i_homo = 1, homo
               IF (Eigenval(i_homo) - Eigenval(homo) .GT. -mp2_env%ri_g0w0%bse_cutoff_occ) THEN
                  homo_incl = i_homo
                  EXIT
               END IF
            END DO
            homo_red = homo - homo_incl + 1
         END IF

         IF (mp2_env%ri_g0w0%bse_cutoff_virt .GT. Eigenval(homo + virtual) - Eigenval(homo + 1) &
             .OR. mp2_env%ri_g0w0%bse_cutoff_virt < 0) THEN
            virt_red = virtual
            virt_incl = virtual
         ELSE
            virt_incl = homo + 1
            DO j_virt = 1, virtual
               IF (Eigenval(homo + j_virt) - Eigenval(homo + 1) .GT. mp2_env%ri_g0w0%bse_cutoff_virt) THEN
                  virt_incl = j_virt - 1
                  EXIT
               END IF
            END DO
            virt_red = virt_incl
         END IF
      ELSE
         homo_red = homo
         virt_red = virtual
         homo_incl = 1
         virt_incl = virtual
      END IF
      IF (unit_nr > 0) THEN
         IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A29,T71,F10.3)') 'BSE|', 'Cutoff occupied orbitals [eV]', &
               mp2_env%ri_g0w0%bse_cutoff_occ*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No cutoff given for occupied orbitals'
         END IF
         IF (mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A28,T71,F10.3)') 'BSE|', 'Cutoff virtual orbitals [eV]', &
               mp2_env%ri_g0w0%bse_cutoff_virt*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A36)') 'BSE|', 'No cutoff given for virtual orbitals'
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A20,T71,I10)') 'BSE|', 'First occupied index', homo_incl
         WRITE (unit_nr, '(T2,A4,T7,A34,T71,I10)') 'BSE|', 'Last virtual index (not MO index!)', virt_incl
         WRITE (unit_nr, '(T2,A4,T7,A35,T71,F10.3)') 'BSE|', 'Energy of first occupied index [eV]', Eigenval(homo_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A33,T71,F10.3)') 'BSE|', 'Energy of last virtual index [eV]', Eigenval(homo + virt_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A54,T71,F10.3)') 'BSE|', 'Energy difference of first occupied index to HOMO [eV]', &
            -(Eigenval(homo_incl) - Eigenval(homo))*evolt
         WRITE (unit_nr, '(T2,A4,T7,A52,T71,F10.3)') 'BSE|', 'Energy difference of last virtual index to LUMO [eV]', &
            (Eigenval(homo + virt_incl) - Eigenval(homo + 1))*evolt
         WRITE (unit_nr, '(T2,A4,T7,A35,T71,I10)') 'BSE|', 'Number of GW-corrected occupied MOs', mp2_env%ri_g0w0%corr_mos_occ
         WRITE (unit_nr, '(T2,A4,T7,A34,T71,I10)') 'BSE|', 'Number of GW-corrected virtual MOs', mp2_env%ri_g0w0%corr_mos_virt
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF
      IF (unit_nr > 0) THEN
         IF (homo - homo_incl + 1 > mp2_env%ri_g0w0%corr_mos_occ) THEN
            CPABORT("Number of GW-corrected occupied MOs too small for chosen BSE cutoff")
         END IF
         IF (virt_incl > mp2_env%ri_g0w0%corr_mos_virt) THEN
            CPABORT("Number of GW-corrected virtual MOs too small for chosen BSE cutoff")
         END IF
      END IF
      !Truncate full fm_S matrices
      !Allocate new truncated matrices of proper size
      para_env => fm_mat_S_ia_bse%matrix_struct%para_env
      context => fm_mat_S_ia_bse%matrix_struct%context

      CALL cp_fm_struct_create(fm_struct_ia, para_env, context, dimen_RI, homo_red*virt_red)
      CALL cp_fm_struct_create(fm_struct_ij, para_env, context, dimen_RI, homo_red*homo_red)
      CALL cp_fm_struct_create(fm_struct_ab, para_env, context, dimen_RI, virt_red*virt_red)

      CALL cp_fm_create(fm_mat_S_trunc, fm_struct_ia, "fm_S_trunc")
      CALL cp_fm_create(fm_mat_S_ij_trunc, fm_struct_ij, "fm_S_ij_trunc")
      CALL cp_fm_create(fm_mat_S_ab_trunc, fm_struct_ab, "fm_S_ab_trunc")

      !Copy parts of original matrices to truncated ones
      IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0 .OR. mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
         !Truncate eigenvals
         ALLOCATE (Eigenval_reduced(homo_red + virt_red))
         Eigenval_reduced(:) = Eigenval(homo_incl:homo + virt_incl)

         CALL truncate_fm(fm_mat_S_trunc, fm_mat_S_ia_bse, virtual, &
                          homo_red, virt_red, unit_nr, mp2_env, &
                          nrow_offset=homo_incl)
         CALL truncate_fm(fm_mat_S_ij_trunc, fm_mat_S_ij_bse, homo, &
                          homo_red, homo_red, unit_nr, mp2_env, &
                          homo_incl, homo_incl)
         CALL truncate_fm(fm_mat_S_ab_trunc, fm_mat_S_ab_bse, virtual, &
                          virt_red, virt_red, unit_nr, mp2_env)

      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No truncation of BSE matrices applied'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF
         ALLOCATE (Eigenval_reduced(homo_red + virt_red))
         Eigenval_reduced(:) = Eigenval(:)
         CALL cp_fm_to_fm_submat_general(fm_mat_S_ia_bse, fm_mat_S_trunc, dimen_RI, homo_red*virt_red, &
                                         1, 1, 1, 1, context)
         CALL cp_fm_to_fm_submat_general(fm_mat_S_ij_bse, fm_mat_S_ij_trunc, dimen_RI, homo_red*homo_red, &
                                         1, 1, 1, 1, context)
         CALL cp_fm_to_fm_submat_general(fm_mat_S_ab_bse, fm_mat_S_ab_trunc, dimen_RI, virt_red*virt_red, &
                                         1, 1, 1, 1, context)
      END IF

      CALL cp_fm_struct_release(fm_struct_ia)
      CALL cp_fm_struct_release(fm_struct_ij)
      CALL cp_fm_struct_release(fm_struct_ab)

      NULLIFY (para_env)
      NULLIFY (context)

      CALL timestop(handle)

   END SUBROUTINE truncate_BSE_matrices

END MODULE
