!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      11.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_util
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_constants,                 ONLY: bse_screening_alpha,&
                                              bse_screening_rpa,&
                                              bse_screening_tdhf,&
                                              use_mom_ref_coac
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE moments_utils,                   ONLY: get_reference_point
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE rpa_communication,               ONLY: communicate_buffer
   USE util,                            ONLY: sort,&
                                              sort_unique
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_util'

   PUBLIC :: mult_B_with_W, fm_general_add_bse, truncate_fm, &
             deallocate_matrices_bse, comp_eigvec_coeff_BSE, sort_excitations, &
             estimate_BSE_resources, filter_eigvec_contrib, truncate_BSE_matrices, &
             determine_cutoff_indices, adapt_BSE_input_params, get_multipoles_mo, &
             reshuffle_eigvec, print_bse_nto_cubes, trace_exciton_descr

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

      IF (info_chol /= 0) THEN
         CALL cp_abort(__LOCATION__, 'Cholesky decomposition failed for static polarization in BSE')
      END IF

      ! calculate [1+Q(i0)]^-1
      CALL cp_fm_cholesky_invert(fm_mat_Q_static_bse_gemm)

      ! symmetrize the result
      CALL cp_fm_uplo_to_full(fm_mat_Q_static_bse_gemm, fm_work)

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
!> to A_ia, jb which needs MPI communication.
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

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_out
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_in
      REAL(kind=dp)                                      :: beta
      INTEGER, INTENT(IN)                                :: nrow_secidx_in, ncol_secidx_in, &
                                                            nrow_secidx_out, ncol_secidx_out
      INTEGER                                            :: unit_nr
      INTEGER, DIMENSION(4)                              :: reordering
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_general_add_bse'

      INTEGER :: col_idx_loc, dummy, handle, handle2, i_entry_rec, idx_col_out, idx_row_out, ii, &
         iproc, jj, ncol_block_in, ncol_block_out, ncol_local_in, ncol_local_out, nprocs, &
         nrow_block_in, nrow_block_out, nrow_local_in, nrow_local_out, proc_send, row_idx_loc, &
         send_pcol, send_prow
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

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

            send_prow = fm_out%matrix_struct%g2p_row(idx_row_out)
            send_pcol = fm_out%matrix_struct%g2p_col(idx_col_out)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A27)') 'BSE|DEBUG|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

            send_prow = fm_out%matrix_struct%g2p_row(idx_row_out)
            send_pcol = fm_out%matrix_struct%g2p_col(idx_col_out)

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
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A21)') 'BSE|DEBUG|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A24,A10)') 'BSE|DEBUG|', 'Writing from buffers to ', fm_out%name
      END IF

      ! fill fm_out with the entries from buffer_rec, i.e. buffer_rec are parts of fm_in
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_out, nprocs, num_entries_rec, buffer_rec, beta) &
!$OMP PRIVATE(iproc, i_entry_rec, ii, jj)
      DO iproc = 0, nprocs - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            ii = fm_out%matrix_struct%g2l_row(buffer_rec(iproc)%indx(i_entry_rec, 1))
            jj = fm_out%matrix_struct%g2l_col(buffer_rec(iproc)%indx(i_entry_rec, 2))

            fm_out%local_data(ii, jj) = fm_out%local_data(ii, jj) + beta*buffer_rec(iproc)%msg(i_entry_rec)
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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
         ncol_local_in, ncol_local_out, nprocs, nrow_block_in, nrow_block_out, nrow_local_in, &
         nrow_local_out, proc_send, row_idx_loc, send_pcol, send_prow
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

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

            send_prow = fm_out%matrix_struct%g2p_row(idx_row_in)
            send_pcol = fm_out%matrix_struct%g2p_col(idx_col_out)

            proc_send = fm_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A27)') 'BSE|DEBUG|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

            send_prow = fm_out%matrix_struct%g2p_row(idx_row_in)

            send_pcol = fm_out%matrix_struct%g2p_col(idx_col_out)

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
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A21)') 'BSE|DEBUG|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A24,A10)') 'BSE|DEBUG|', 'Writing from buffers to ', fm_out%name
      END IF

      ! fill fm_out with the entries from buffer_rec, i.e. buffer_rec are parts of fm_in
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_out, nprocs, num_entries_rec, buffer_rec) &
!$OMP PRIVATE(iproc, i_entry_rec, ii, jj)
      DO iproc = 0, nprocs - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            ii = fm_out%matrix_struct%g2l_row(buffer_rec(iproc)%indx(i_entry_rec, 1))
            jj = fm_out%matrix_struct%g2l_col(buffer_rec(iproc)%indx(i_entry_rec, 2))

            fm_out%local_data(ii, jj) = fm_out%local_data(ii, jj) + buffer_rec(iproc)%msg(i_entry_rec)
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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
!> \brief ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_trunc ...
!> \param fm_mat_S_ij_trunc ...
!> \param fm_mat_S_ab_trunc ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_bse(fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                      fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                      fm_mat_Q_static_bse_gemm, mp2_env)

      TYPE(cp_fm_type), INTENT(INOUT) :: fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_trunc, &
         fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, fm_mat_Q_static_bse_gemm
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_bse'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_release(fm_mat_S_bar_ia_bse)
      CALL cp_fm_release(fm_mat_S_bar_ij_bse)
      CALL cp_fm_release(fm_mat_S_trunc)
      CALL cp_fm_release(fm_mat_S_ij_trunc)
      CALL cp_fm_release(fm_mat_S_ab_trunc)
      CALL cp_fm_release(fm_mat_Q_static_bse_gemm)
      IF (mp2_env%bse%do_nto_analysis) THEN
         DEALLOCATE (mp2_env%bse%bse_nto_state_list_final)
      END IF

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
      INTEGER(KIND=int_8)                                :: full_dim
      REAL(KIND=dp)                                      :: mem_est, mem_est_per_rank

      CALL timeset(routineN, handle)

      ! Number of matrices with size of A in TDA is 2 (A itself and W_ijab)
      num_BSE_matrices = 2
      ! With the full diagonalization of ABBA, we need several auxiliary matrices in the process
      ! The maximum number is 2 + 2 + 6 (additional B and C matrix as well as 6 matrices to create C)
      IF (bse_abba) THEN
         num_BSE_matrices = 10
      END IF

      full_dim = (INT(homo_red, KIND=int_8)**2*INT(virtual_red, KIND=int_8)**2)*INT(num_BSE_matrices, KIND=int_8)
      mem_est = REAL(8*full_dim, KIND=dp)/REAL(1024**3, KIND=dp)
      mem_est_per_rank = REAL(mem_est/para_env%num_pe, KIND=dp)

      IF (unit_nr > 0) THEN
         ! WRITE (unit_nr, '(T2,A4,T7,A40,T68,F13.3)') 'BSE|', 'Total peak memory estimate from BSE [GB]', &
         !    mem_est
         WRITE (unit_nr, '(T2,A4,T7,A40,T68,ES13.3)') 'BSE|', 'Total peak memory estimate from BSE [GB]', &
            mem_est
         WRITE (unit_nr, '(T2,A4,T7,A47,T68,F13.3)') 'BSE|', 'Peak memory estimate per MPI rank from BSE [GB]', &
            mem_est_per_rank
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF
      ! Rough estimation of diagonalization runtimes. Baseline was a full BSE Naphthalene
      ! run with 11000x11000 entries in A/B/C, which took 10s on 32 ranks
      diag_runtime_est = REAL(INT(homo_red, KIND=int_8)*INT(virtual_red, KIND=int_8)/11000_int_8, KIND=dp)**3* &
                         10*32/REAL(para_env%num_pe, KIND=dp)

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
            IF (ABS(eigvec_entry) > mp2_env%bse%eps_x) THEN
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
            IF (ABS(eigvec_entry) > mp2_env%bse%eps_x) THEN
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
!> \brief Reads cutoffs for BSE from mp2_env and compares to energies in Eigenval to extract
!>        reduced homo/virtual and
!> \param Eigenval array (1d) with energies, can be e.g. from GW or DFT
!> \param homo Total number of occupied orbitals
!> \param virtual Total number of unoccupied orbitals
!> \param homo_red Total number of occupied orbitals to include after cutoff
!> \param virt_red Total number of unoccupied orbitals to include after ctuoff
!> \param homo_incl First occupied index to include after cutoff
!> \param virt_incl Last unoccupied index to include after cutoff
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE determine_cutoff_indices(Eigenval, &
                                       homo, virtual, &
                                       homo_red, virt_red, &
                                       homo_incl, virt_incl, &
                                       mp2_env)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: homo, virtual
      INTEGER, INTENT(OUT)                               :: homo_red, virt_red, homo_incl, virt_incl
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'determine_cutoff_indices'

      INTEGER                                            :: handle, i_homo, j_virt

      CALL timeset(routineN, handle)
      ! Determine index in homo and virtual for truncation
      ! Uses indices of outermost orbitals within energy range (-mp2_env%bse%bse_cutoff_occ,mp2_env%bse%bse_cutoff_empty)
      IF (mp2_env%bse%bse_cutoff_occ > 0 .OR. mp2_env%bse%bse_cutoff_empty > 0) THEN
         IF (-mp2_env%bse%bse_cutoff_occ .LT. Eigenval(1) - Eigenval(homo) &
             .OR. mp2_env%bse%bse_cutoff_occ < 0) THEN
            homo_red = homo
            homo_incl = 1
         ELSE
            homo_incl = 1
            DO i_homo = 1, homo
               IF (Eigenval(i_homo) - Eigenval(homo) .GT. -mp2_env%bse%bse_cutoff_occ) THEN
                  homo_incl = i_homo
                  EXIT
               END IF
            END DO
            homo_red = homo - homo_incl + 1
         END IF

         IF (mp2_env%bse%bse_cutoff_empty .GT. Eigenval(homo + virtual) - Eigenval(homo + 1) &
             .OR. mp2_env%bse%bse_cutoff_empty < 0) THEN
            virt_red = virtual
            virt_incl = virtual
         ELSE
            virt_incl = homo + 1
            DO j_virt = 1, virtual
               IF (Eigenval(homo + j_virt) - Eigenval(homo + 1) .GT. mp2_env%bse%bse_cutoff_empty) THEN
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
!> \param Eigenval_scf ...
!> \param Eigenval ...
!> \param Eigenval_reduced ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param bse_lev_virt ...
!> \param homo_red ...
!> \param virt_red ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE truncate_BSE_matrices(fm_mat_S_ia_bse, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                    Eigenval_scf, Eigenval, Eigenval_reduced, &
                                    homo, virtual, dimen_RI, unit_nr, &
                                    bse_lev_virt, &
                                    homo_red, virt_red, &
                                    mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ia_bse, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_S_trunc, fm_mat_S_ij_trunc, &
                                                            fm_mat_S_ab_trunc
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval_scf, Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, unit_nr, &
                                                            bse_lev_virt
      INTEGER, INTENT(OUT)                               :: homo_red, virt_red
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'truncate_BSE_matrices'

      INTEGER                                            :: handle, homo_incl, virt_incl
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ab, fm_struct_ia, fm_struct_ij
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Determine index in homo and virtual for truncation
      ! Uses indices of outermost orbitals within energy range (-mp2_env%bse%bse_cutoff_occ,mp2_env%bse%bse_cutoff_empty)

      CALL determine_cutoff_indices(Eigenval_scf, &
                                    homo, virtual, &
                                    homo_red, virt_red, &
                                    homo_incl, virt_incl, &
                                    mp2_env)

      IF (unit_nr > 0) THEN
         IF (mp2_env%bse%bse_cutoff_occ > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A29,T71,F10.3)') 'BSE|', 'Cutoff occupied orbitals [eV]', &
               mp2_env%bse%bse_cutoff_occ*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No cutoff given for occupied orbitals'
         END IF
         IF (mp2_env%bse%bse_cutoff_empty > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A26,T71,F10.3)') 'BSE|', 'Cutoff empty orbitals [eV]', &
               mp2_env%bse%bse_cutoff_empty*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A34)') 'BSE|', 'No cutoff given for empty orbitals'
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A20,T71,I10)') 'BSE|', 'First occupied index', homo_incl
         WRITE (unit_nr, '(T2,A4,T7,A32,T71,I10)') 'BSE|', 'Last empty index (not MO index!)', virt_incl
         WRITE (unit_nr, '(T2,A4,T7,A35,T71,F10.3)') 'BSE|', 'Energy of first occupied index [eV]', Eigenval(homo_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A31,T71,F10.3)') 'BSE|', 'Energy of last empty index [eV]', Eigenval(homo + virt_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A54,T71,F10.3)') 'BSE|', 'Energy difference of first occupied index to HOMO [eV]', &
            -(Eigenval(homo_incl) - Eigenval(homo))*evolt
         WRITE (unit_nr, '(T2,A4,T7,A50,T71,F10.3)') 'BSE|', 'Energy difference of last empty index to LUMO [eV]', &
            (Eigenval(homo + virt_incl) - Eigenval(homo + 1))*evolt
         WRITE (unit_nr, '(T2,A4,T7,A35,T71,I10)') 'BSE|', 'Number of GW-corrected occupied MOs', mp2_env%ri_g0w0%corr_mos_occ
         WRITE (unit_nr, '(T2,A4,T7,A32,T71,I10)') 'BSE|', 'Number of GW-corrected empty MOs', mp2_env%ri_g0w0%corr_mos_virt
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
      IF (mp2_env%bse%bse_cutoff_occ > 0 .OR. mp2_env%bse%bse_cutoff_empty > 0) THEN
         !Truncate eigenvals
         ALLOCATE (Eigenval_reduced(homo_red + virt_red))
         ! Include USE_KS_ENERGIES input
         IF (mp2_env%bse%use_ks_energies) THEN
            Eigenval_reduced(:) = Eigenval_scf(homo_incl:homo + virt_incl)
         ELSE
            Eigenval_reduced(:) = Eigenval(homo_incl:homo + virt_incl)
         END IF

         CALL truncate_fm(fm_mat_S_trunc, fm_mat_S_ia_bse, virtual, &
                          homo_red, virt_red, unit_nr, mp2_env, &
                          nrow_offset=homo_incl)
         CALL truncate_fm(fm_mat_S_ij_trunc, fm_mat_S_ij_bse, homo, &
                          homo_red, homo_red, unit_nr, mp2_env, &
                          homo_incl, homo_incl)
         CALL truncate_fm(fm_mat_S_ab_trunc, fm_mat_S_ab_bse, bse_lev_virt, &
                          virt_red, virt_red, unit_nr, mp2_env)

      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No truncation of BSE matrices applied'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF
         ALLOCATE (Eigenval_reduced(homo_red + virt_red))
         ! Include USE_KS_ENERGIES input
         IF (mp2_env%bse%use_ks_energies) THEN
            Eigenval_reduced(:) = Eigenval_scf(:)
         ELSE
            Eigenval_reduced(:) = Eigenval(:)
         END IF
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

! **************************************************************************************************
!> \brief ...
!> \param fm_eigvec ...
!> \param fm_eigvec_reshuffled ...
!> \param homo ...
!> \param virtual ...
!> \param n_exc ...
!> \param do_transpose ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE reshuffle_eigvec(fm_eigvec, fm_eigvec_reshuffled, homo, virtual, n_exc, do_transpose, &
                               unit_nr, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_eigvec_reshuffled
      INTEGER, INTENT(IN)                                :: homo, virtual, n_exc
      LOGICAL, INTENT(IN)                                :: do_transpose
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reshuffle_eigvec'

      INTEGER                                            :: handle, my_m_col, my_n_row
      INTEGER, DIMENSION(4)                              :: reordering
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_eigvec_col, &
                                                            fm_struct_eigvec_reshuffled
      TYPE(cp_fm_type)                                   :: fm_eigvec_col

      CALL timeset(routineN, handle)

      ! Define reordering:
      ! (ia,11) to (a1,i1) for transposition
      ! (ia,11) to (i1,a1) for default
      IF (do_transpose) THEN
         reordering = (/2, 3, 1, 4/)
         my_n_row = virtual
         my_m_col = homo
      ELSE
         reordering = (/1, 3, 2, 4/)
         my_n_row = homo
         my_m_col = virtual
      END IF

      CALL cp_fm_struct_create(fm_struct_eigvec_col, &
                               fm_eigvec%matrix_struct%para_env, fm_eigvec%matrix_struct%context, &
                               homo*virtual, 1)
      CALL cp_fm_struct_create(fm_struct_eigvec_reshuffled, &
                               fm_eigvec%matrix_struct%para_env, fm_eigvec%matrix_struct%context, &
                               my_n_row, my_m_col)

      ! Resort indices
      CALL cp_fm_create(fm_eigvec_col, fm_struct_eigvec_col, name="BSE_column_vector")
      CALL cp_fm_set_all(fm_eigvec_col, 0.0_dp)
      CALL cp_fm_create(fm_eigvec_reshuffled, fm_struct_eigvec_reshuffled, name="BSE_reshuffled_eigenvector")
      CALL cp_fm_set_all(fm_eigvec_reshuffled, 0.0_dp)
      ! Fill matrix
      CALL cp_fm_to_fm_submat(fm_eigvec, fm_eigvec_col, &
                              homo*virtual, 1, &
                              1, n_exc, &
                              1, 1)
      ! Reshuffle
      CALL fm_general_add_bse(fm_eigvec_reshuffled, fm_eigvec_col, 1.0_dp, &
                              virtual, 1, &
                              1, 1, &
                              unit_nr, reordering, mp2_env)

      CALL cp_fm_release(fm_eigvec_col)
      CALL cp_fm_struct_release(fm_struct_eigvec_col)
      CALL cp_fm_struct_release(fm_struct_eigvec_reshuffled)

      CALL timestop(handle)

   END SUBROUTINE reshuffle_eigvec

! **************************************************************************************************
!> \brief Borrowed from the tddfpt module with slight adaptions
!> \param qs_env ...
!> \param mos ...
!> \param istate ...
!> \param info_approximation ...
!> \param stride ...
!> \param append_cube ...
!> \param print_section ...
! **************************************************************************************************
   SUBROUTINE print_bse_nto_cubes(qs_env, mos, istate, info_approximation, &
                                  stride, append_cube, print_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      INTEGER, INTENT(IN)                                :: istate
      CHARACTER(LEN=10)                                  :: info_approximation
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL, INTENT(IN)                                :: append_cube
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_bse_nto_cubes'

      CHARACTER(LEN=default_path_length)                 :: filename, info_approx_trunc, &
                                                            my_pos_cube, title
      INTEGER                                            :: handle, i, iset, nmo, unit_nr_cube
      LOGICAL                                            :: mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: wf_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL auxbas_pw_pool%create_pw(wf_r)
      CALL auxbas_pw_pool%create_pw(wf_g)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      my_pos_cube = "REWIND"
      IF (append_cube) THEN
         my_pos_cube = "APPEND"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set)

      DO iset = 1, 2
         CALL get_mo_set(mo_set=mos(iset), mo_coeff=mo_coeff, nmo=nmo)
         DO i = 1, nmo
            CALL calculate_wavefunction(mo_coeff, i, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                        cell, dft_control, particle_set, pw_env)
            IF (iset == 1) THEN
               WRITE (filename, '(A6,I3.3,A5,I2.2,a11)') "_NEXC_", istate, "_NTO_", i, "_Hole_State"
            ELSEIF (iset == 2) THEN
               WRITE (filename, '(A6,I3.3,A5,I2.2,a15)') "_NEXC_", istate, "_NTO_", i, "_Particle_State"
            END IF
            info_approx_trunc = TRIM(ADJUSTL(info_approximation))
            info_approx_trunc = info_approx_trunc(2:LEN_TRIM(info_approx_trunc) - 1)
            filename = TRIM(info_approx_trunc)//TRIM(filename)
            mpi_io = .TRUE.
            unit_nr_cube = cp_print_key_unit_nr(logger, print_section, '', extension=".cube", &
                                                middle_name=TRIM(filename), file_position=my_pos_cube, &
                                                log_filename=.FALSE., ignore_should_output=.TRUE., mpi_io=mpi_io)
            IF (iset == 1) THEN
               WRITE (title, *) "Natural Transition Orbital Hole State", i
            ELSEIF (iset == 2) THEN
               WRITE (title, *) "Natural Transition Orbital Particle State", i
            END IF
            CALL cp_pw_to_cube(wf_r, unit_nr_cube, title, particles=particles, stride=stride, mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr_cube, logger, print_section, '', &
                                              ignore_should_output=.TRUE., mpi_io=mpi_io)
         END DO
      END DO

      CALL auxbas_pw_pool%give_back_pw(wf_g)
      CALL auxbas_pw_pool%give_back_pw(wf_r)

      CALL timestop(handle)
   END SUBROUTINE print_bse_nto_cubes

! **************************************************************************************************
!> \brief Checks BSE input section and adapts them if necessary
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE adapt_BSE_input_params(homo, virtual, unit_nr, mp2_env, qs_env)

      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'adapt_BSE_input_params'

      INTEGER                                            :: handle, i, j, n, ndim_periodic_cell, &
                                                            ndim_periodic_poisson, &
                                                            num_state_list_exceptions
      TYPE(cell_type), POINTER                           :: cell_ref
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env

      CALL timeset(routineN, handle)
      ! Get environment infos for later usage
      NULLIFY (pw_env, cell_ref, poisson_env)
      CALL get_qs_env(qs_env, pw_env=pw_env, cell_ref=cell_ref)
      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      ndim_periodic_poisson = COUNT(poisson_env%parameters%periodic == 1)
      ndim_periodic_cell = SUM(cell_ref%perd(1:3)) ! Borrowed from cell_methods.F/write_cell_low

      ! Handle negative NUM_PRINT_EXC
      IF (mp2_env%bse%num_print_exc < 0 .OR. &
          mp2_env%bse%num_print_exc > homo*virtual) THEN
         mp2_env%bse%num_print_exc = homo*virtual
         IF (unit_nr > 0) THEN
            CALL cp_hint(__LOCATION__, &
                         "Keyword NUM_PRINT_EXC is either negative or too large. "// &
                         "Printing all computed excitations.")
         END IF
      END IF

      ! Default to NUM_PRINT_EXC if too large or negative,
      ! but only if NTOs are called - would be confusing for the user otherwise
      ! Prepare and adapt user inputs for NTO analysis
      ! Logic: Explicit state list overrides NUM_PRINT_EXC_NTOS
      !        If only NUM_PRINT_EXC_NTOS is given, we write the array 1,...,NUM_PRINT_EXC_NTOS to
      !        bse_nto_state_list
      IF (mp2_env%bse%do_nto_analysis) THEN
         IF (mp2_env%bse%explicit_nto_list) THEN
            IF (mp2_env%bse%num_print_exc_ntos > 0) THEN
               IF (unit_nr > 0) THEN
                  CALL cp_hint(__LOCATION__, &
                               "Keywords NUM_PRINT_EXC_NTOS and STATE_LIST are both given in input. "// &
                               "Overriding NUM_PRINT_EXC_NTOS.")
               END IF
            END IF
            ! Check if all states are within the range
            ! Count them and initialize new array afterwards
            num_state_list_exceptions = 0
            DO i = 1, SIZE(mp2_env%bse%bse_nto_state_list)
               IF (mp2_env%bse%bse_nto_state_list(i) < 1 .OR. &
                   mp2_env%bse%bse_nto_state_list(i) > mp2_env%bse%num_print_exc) THEN
                  num_state_list_exceptions = num_state_list_exceptions + 1
               END IF
            END DO
            IF (num_state_list_exceptions > 0) THEN
               IF (unit_nr > 0) THEN
                  CALL cp_hint(__LOCATION__, &
                               "STATE_LIST contains indices outside the range of included excitation levels. "// &
                               "Ignoring these states.")
               END IF
            END IF
            n = SIZE(mp2_env%bse%bse_nto_state_list) - num_state_list_exceptions
            ALLOCATE (mp2_env%bse%bse_nto_state_list_final(n))
            mp2_env%bse%bse_nto_state_list_final(:) = 0
            i = 1
            DO j = 1, SIZE(mp2_env%bse%bse_nto_state_list)
               IF (mp2_env%bse%bse_nto_state_list(j) >= 1 .AND. &
                   mp2_env%bse%bse_nto_state_list(j) <= mp2_env%bse%num_print_exc) THEN
                  mp2_env%bse%bse_nto_state_list_final(i) = mp2_env%bse%bse_nto_state_list(j)
                  i = i + 1
               END IF
            END DO

            mp2_env%bse%num_print_exc_ntos = SIZE(mp2_env%bse%bse_nto_state_list_final)
         ELSE
            IF (mp2_env%bse%num_print_exc_ntos > mp2_env%bse%num_print_exc .OR. &
                mp2_env%bse%num_print_exc_ntos < 0) THEN
               mp2_env%bse%num_print_exc_ntos = mp2_env%bse%num_print_exc
            END IF
            ALLOCATE (mp2_env%bse%bse_nto_state_list_final(mp2_env%bse%num_print_exc_ntos))
            DO i = 1, mp2_env%bse%num_print_exc_ntos
               mp2_env%bse%bse_nto_state_list_final(i) = i
            END DO
         END IF
      END IF

      ! Takes care of triplet states, when oscillator strengths are 0
      IF (mp2_env%bse%bse_spin_config /= 0 .AND. &
          mp2_env%bse%eps_nto_osc_str > 0) THEN
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Cannot apply EPS_OSC_STR for Triplet excitations. "// &
                         "Resetting EPS_OSC_STR to default.")
         END IF
         mp2_env%bse%eps_nto_osc_str = -1.0_dp
      END IF

      ! Take care of number for computed exciton descriptors
      IF (mp2_env%bse%num_print_exc_descr < 0 .OR. &
          mp2_env%bse%num_print_exc_descr > mp2_env%bse%num_print_exc) THEN
         IF (unit_nr > 0) THEN
            CALL cp_hint(__LOCATION__, &
                         "Keyword NUM_PRINT_EXC_DESCR is either negative or too large. "// &
                         "Printing exciton descriptors up to NUM_PRINT_EXC.")
         END IF
         mp2_env%bse%num_print_exc_descr = mp2_env%bse%num_print_exc
      END IF

      ! Handle screening factor options
      IF (mp2_env%BSE%screening_factor > 0.0_dp) THEN
         IF (mp2_env%BSE%screening_method /= bse_screening_alpha) THEN
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "Screening factor is only supported for &SCREENING_IN_W ALPHA. "// &
                            "Resetting SCREENING_IN_W to ALPHA.")
            END IF
            mp2_env%BSE%screening_method = bse_screening_alpha
         END IF
         IF (mp2_env%BSE%screening_factor > 1.0_dp) THEN
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "Screening factor is larger than 1.0. ")
            END IF
         END IF
      END IF

      IF (mp2_env%BSE%screening_factor < 0.0_dp .AND. &
          mp2_env%BSE%screening_method == bse_screening_alpha) THEN
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Screening factor is negative. Defaulting to 0.25")
         END IF
         mp2_env%BSE%screening_factor = 0.25_dp
      END IF

      IF (mp2_env%BSE%screening_factor == 0.0_dp) THEN
         ! Use RPA internally in this case
         mp2_env%BSE%screening_method = bse_screening_rpa
      END IF
      IF (mp2_env%BSE%screening_factor == 1.0_dp) THEN
         ! Use TDHF internally in this case
         mp2_env%BSE%screening_method = bse_screening_tdhf
      END IF

      ! Add warning for usage of KS energies
      IF (mp2_env%bse%use_ks_energies) THEN
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Using KS energies for BSE calculations. Therefore, no quantities "// &
                         "of the preceeding GW calculation enter the BSE.")
         END IF
      END IF

      ! Add warning if periodic calculation is invoked
      IF (ndim_periodic_poisson /= 0) THEN
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Poisson solver should be invoked by PERIODIC NONE. "// &
                         "The applied length gauge might give misleading results for "// &
                         "oscillator strengths.")
         END IF
      END IF
      IF (ndim_periodic_cell /= 0) THEN
         IF (unit_nr > 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "CELL in SUBSYS should be invoked with PERIODIC NONE. "// &
                         "The applied length gauge might give misleading results for "// &
                         "oscillator strengths.")
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE adapt_BSE_input_params

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param fm_multipole_ai_trunc ...
!> \param fm_multipole_ij_trunc ...
!> \param fm_multipole_ab_trunc ...
!> \param qs_env ...
!> \param mo_coeff ...
!> \param rpoint ...
!> \param n_moments ...
!> \param homo_red ...
!> \param virtual_red ...
!> \param context_BSE ...
! **************************************************************************************************
   SUBROUTINE get_multipoles_mo(fm_multipole_ai_trunc, fm_multipole_ij_trunc, fm_multipole_ab_trunc, &
                                qs_env, mo_coeff, rpoint, n_moments, &
                                homo_red, virtual_red, context_BSE)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: fm_multipole_ai_trunc, &
                                                            fm_multipole_ij_trunc, &
                                                            fm_multipole_ab_trunc
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: rpoint
      INTEGER, INTENT(IN)                                :: n_moments, homo_red, virtual_red
      TYPE(cp_blacs_env_type), POINTER                   :: context_BSE

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_multipoles_mo'

      INTEGER                                            :: handle, idir, n_multipole, n_occ, &
                                                            n_virt, nao
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_struct_type), POINTER :: fm_struct_mp_ab_trunc, fm_struct_mp_ai_trunc, &
         fm_struct_mp_ij_trunc, fm_struct_multipoles_ao, fm_struct_nao_nao
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_multipole_per_dir
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_multipole, matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env_BSE
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      !First, we calculate the AO dipoles
      NULLIFY (sab_orb, matrix_s)
      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      matrix_s=matrix_s, &
                      sab_orb=sab_orb)

      ! Use the same blacs environment as for the MO coefficients to ensure correct multiplication dbcsr x fm later on
      fm_struct_multipoles_ao => mos(1)%mo_coeff%matrix_struct
      ! BSE has different contexts and blacsenvs
      para_env_BSE => context_BSE%para_env
      ! Get size of multipole tensor
      n_multipole = (6 + 11*n_moments + 6*n_moments**2 + n_moments**3)/6 - 1
      NULLIFY (matrix_multipole)
      CALL dbcsr_allocate_matrix_set(matrix_multipole, n_multipole)
      ALLOCATE (fm_multipole_per_dir(n_multipole))
      DO idir = 1, n_multipole
         CALL dbcsr_init_p(matrix_multipole(idir)%matrix)
         CALL dbcsr_create(matrix_multipole(idir)%matrix, name="ao_multipole", &
                           template=matrix_s(1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_multipole(idir)%matrix, sab_orb)
         CALL dbcsr_set(matrix_multipole(idir)%matrix, 0._dp)
      END DO

      CALL get_reference_point(rpoint, qs_env=qs_env, reference=use_mom_ref_coac, ref_point=ref_point)

      CALL build_local_moment_matrix(qs_env, matrix_multipole, n_moments, ref_point=rpoint)

      NULLIFY (sab_orb)

      ! Now we transform them to MO
      ! n_occ is the number of occupied MOs, nao the number of all AOs
      ! Writing homo to n_occ instead if nmo,
      ! takes care of ADDED_MOS, which would overwrite nmo, if invoked
      CALL get_mo_set(mo_set=mos(1), homo=n_occ, nao=nao)
      n_virt = nao - n_occ

      ! At the end, we need four different layouts of matrices in this multiplication, e.g. for a dipole:
      ! D_pq = full multipole matrix for occupied and unoccupied
      ! Final result:D_pq= C_{mu p}        <\mu|\vec{r}|\nu>        C_{\nu q}              EQ.I
      !                   \_______/         \___________/          \______/
      !                    fm_coeff            matrix_multipole              fm_coeff
      !                    (EQ.Ia)             (EQ.Ib)              (EQ.Ia)
      ! Intermediate work matrices:
      ! fm_work =                 <\mu|\vec{r}|\nu>        C_{\nu q}              EQ.II

      ! Struct for the full multipole matrix
      CALL cp_fm_struct_create(fm_struct_nao_nao, &
                               fm_struct_multipoles_ao%para_env, fm_struct_multipoles_ao%context, &
                               nao, nao)

      ! At the very end, we copy the multipoles corresponding to truncated BSE indices in i and a
      CALL cp_fm_struct_create(fm_struct_mp_ai_trunc, para_env_BSE, &
                               context_BSE, virtual_red, homo_red)
      CALL cp_fm_struct_create(fm_struct_mp_ij_trunc, para_env_BSE, &
                               context_BSE, homo_red, homo_red)
      CALL cp_fm_struct_create(fm_struct_mp_ab_trunc, para_env_BSE, &
                               context_BSE, virtual_red, virtual_red)
      DO idir = 1, n_multipole
         CALL cp_fm_create(fm_multipole_ai_trunc(idir), matrix_struct=fm_struct_mp_ai_trunc, &
                           name="dipoles_mo_ai_trunc")
         CALL cp_fm_set_all(fm_multipole_ai_trunc(idir), 0.0_dp)
         CALL cp_fm_create(fm_multipole_ij_trunc(idir), matrix_struct=fm_struct_mp_ij_trunc, &
                           name="dipoles_mo_ij_trunc")
         CALL cp_fm_set_all(fm_multipole_ij_trunc(idir), 0.0_dp)
         CALL cp_fm_create(fm_multipole_ab_trunc(idir), matrix_struct=fm_struct_mp_ab_trunc, &
                           name="dipoles_mo_ab_trunc")
         CALL cp_fm_set_all(fm_multipole_ab_trunc(idir), 0.0_dp)
      END DO

      ! Need another temporary matrix to store intermediate result from right multiplication
      ! D = C_{mu a}        <\mu|\vec{r}|\nu>        C_{\nu i}
      CALL cp_fm_create(fm_work, matrix_struct=fm_struct_nao_nao, name="multipole_work")
      CALL cp_fm_set_all(fm_work, 0.0_dp)

      DO idir = 1, n_multipole
         ! Create the full multipole matrix per direction
         CALL cp_fm_create(fm_multipole_per_dir(idir), matrix_struct=fm_struct_nao_nao, name="multipoles_mo")
         CALL cp_fm_set_all(fm_multipole_per_dir(idir), 0.0_dp)
         ! Fill final (MO) multipole matrix
         CALL cp_dbcsr_sm_fm_multiply(matrix_multipole(idir)%matrix, mo_coeff(1), &
                                      fm_work, ncol=nao)
         ! Now obtain the multipoles by the final multiplication;
         ! We do that inside the loop to obtain multipoles per axis for print
         CALL parallel_gemm('T', 'N', nao, nao, nao, 1.0_dp, mo_coeff(1), fm_work, 0.0_dp, fm_multipole_per_dir(idir))
         ! Truncate full matrix to the BSE indices
         ! D_ai
         CALL cp_fm_to_fm_submat_general(fm_multipole_per_dir(idir), &
                                         fm_multipole_ai_trunc(idir), &
                                         virtual_red, &
                                         homo_red, &
                                         n_occ + 1, &
                                         n_occ - homo_red + 1, &
                                         1, &
                                         1, &
                                         fm_multipole_per_dir(idir)%matrix_struct%context)
         ! D_ij
         CALL cp_fm_to_fm_submat_general(fm_multipole_per_dir(idir), &
                                         fm_multipole_ij_trunc(idir), &
                                         homo_red, &
                                         homo_red, &
                                         n_occ - homo_red + 1, &
                                         n_occ - homo_red + 1, &
                                         1, &
                                         1, &
                                         fm_multipole_per_dir(idir)%matrix_struct%context)
         ! D_ab
         CALL cp_fm_to_fm_submat_general(fm_multipole_per_dir(idir), &
                                         fm_multipole_ab_trunc(idir), &
                                         virtual_red, &
                                         virtual_red, &
                                         n_occ + 1, &
                                         n_occ + 1, &
                                         1, &
                                         1, &
                                         fm_multipole_per_dir(idir)%matrix_struct%context)
      END DO

      !Release matrices and structs
      NULLIFY (fm_struct_multipoles_ao)
      CALL cp_fm_struct_release(fm_struct_mp_ai_trunc)
      CALL cp_fm_struct_release(fm_struct_mp_ij_trunc)
      CALL cp_fm_struct_release(fm_struct_mp_ab_trunc)
      CALL cp_fm_struct_release(fm_struct_nao_nao)
      DO idir = 1, n_multipole
         CALL cp_fm_release(fm_multipole_per_dir(idir))
      END DO
      DEALLOCATE (fm_multipole_per_dir)
      CALL cp_fm_release(fm_work)
      CALL dbcsr_deallocate_matrix_set(matrix_multipole)

      CALL timestop(handle)

   END SUBROUTINE get_multipoles_mo

! **************************************************************************************************
!> \brief Computes trace of form Tr{A^T B C} for exciton descriptors
!> \param fm_A Full Matrix, typically X or Y, in format homo x virtual
!> \param fm_B ...
!> \param fm_C ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE trace_exciton_descr(fm_A, fm_B, fm_C, alpha)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A, fm_B, fm_C
      REAL(KIND=dp), INTENT(OUT)                         :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trace_exciton_descr'

      INTEGER                                            :: handle, ncol_A, ncol_B, ncol_C, nrow_A, &
                                                            nrow_B, nrow_C
      TYPE(cp_fm_type)                                   :: fm_work_ia

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_work_ia, fm_A%matrix_struct)
      CALL cp_fm_get_info(fm_A, nrow_global=nrow_A, ncol_global=ncol_A)
      CALL cp_fm_get_info(fm_B, nrow_global=nrow_B, ncol_global=ncol_B)
      CALL cp_fm_get_info(fm_C, nrow_global=nrow_C, ncol_global=ncol_C)

      ! Check matrix sizes
      CPASSERT(nrow_A == nrow_B .AND. ncol_A == ncol_C .AND. ncol_B == nrow_C)

      CALL cp_fm_set_all(fm_work_ia, 0.0_dp)

      CALL parallel_gemm("N", "N", nrow_A, ncol_A, nrow_C, 1.0_dp, &
                         fm_B, fm_C, 0.0_dp, fm_work_ia)

      CALL cp_fm_trace(fm_A, fm_work_ia, alpha)

      CALL cp_fm_release(fm_work_ia)

      CALL timestop(handle)

   END SUBROUTINE trace_exciton_descr

END MODULE
