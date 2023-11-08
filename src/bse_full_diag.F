!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2017 created [Jan Wilhelm]
! **************************************************************************************************
MODULE bse_full_diag
   !MG check if all are necessary; copied from bse.F

   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_type
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_indxg2p,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
!MG
   USE input_constants,                 ONLY: bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_full_diag'

   PUBLIC :: create_and_diagonalize_A

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param blacs_env_BSE ...
! **************************************************************************************************
   SUBROUTINE create_and_diagonalize_A(fm_mat_S, fm_mat_S_bar_ij_bse, fm_mat_S_ab_bse, &
                                       Eigenval, &
                                       homo, virtual, dimen_RI, bse_spin_config, unit_nr, &
                                       para_env, para_env_RPA, blacs_env_BSE)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_bar_ij_bse, &
                                                            fm_mat_S_ab_bse
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, bse_spin_config
      INTEGER                                            :: unit_nr
      TYPE(mp_para_env_type), INTENT(INOUT)              :: para_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_RPA
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_BSE

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_and_diagonalize_A'

      INTEGER :: a_virt_row, ab_loc, b_virt_col, dummy_proc, handle, i_entry_rec, i_occ_row, &
         i_row_global, ii, ij_loc, iproc, j_col_global, j_occ_col, jj, ncol_block_A, ncol_block_W, &
         ncol_local_A, ncol_local_W, npcol, nprow, nrow_block_A, nrow_block_W, nrow_local_A, &
         nrow_local_W, proc_send, send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_A, col_indices_W, &
                                                            row_indices_A, row_indices_W
      REAL(KIND=dp)                                      :: alpha, eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_A, fm_struct_eigvec, &
                                                            fm_struct_W
      TYPE(cp_fm_type)                                   :: fm_mat_A, fm_mat_eigvec, fm_mat_W
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

!MG for debugging

      CALL timeset(routineN, handle)

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      ! create the blacs env for BSE matrices (NOT fm_mat_S related parallel_gemms!)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      !We have to use the same blacs_env for A as for the matrices fm_mat_S from RPA
      !Logic: A_iajb = (eps_a-eps_i) delta_ij delta _ab + alpha * v_iajb - W_ijab
      ! We create v_iajb and W_ijab, then we communicate entries from local W_ijab
      ! to the full matrix v_iajb. By adding these and the energy diffenences: v_iajb -> A_iajb
      ! We use the A matrix already from the start instead of v
      CALL cp_fm_struct_create(fm_struct_a, context=blacs_env_BSE, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=para_env_RPA)
      CALL cp_fm_create(fm_mat_A, fm_struct_a, name="fm_mat_A_iajb_BSE")
      CALL cp_fm_set_all(fm_mat_A, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct_W, context=blacs_env, nrow_global=homo**2, &
                               ncol_global=virtual**2, para_env=para_env)
      CALL cp_fm_create(fm_mat_W, fm_struct_W, name="fm_mat_W_ijab_BSE")
      CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

      !Create A matrix from GW Energies, v_iajb and W_ijab (different blacs_env!)
      !v_iajb, which is directly initialized as A
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_A)

      !CALL cp_fm_write_formatted(fm_mat_A,unit_nr,"MG: print full A in fulldiag")

      !CALL cp_fm_write_formatted(fm_mat_S_bar_ij_bse,unit_nr,"MG: print full S_ij_bar in fulldiag")
      !CALL cp_fm_write_formatted(fm_mat_S_ab_bse,unit_nr,"MG: print full fm_mat_S_ab_bse in fulldiag")
      !W_ijab
      CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_W)

      !CALL cp_fm_write_formatted(fm_mat_W,unit_nr,"MG: print full W in fulldiag")

      ! A_iajb
      ! We start by moving data from local parts of W_ijab to the full matrix A_iajb using buffers
      CALL cp_fm_get_info(matrix=fm_mat_A, &
                          nrow_local=nrow_local_A, &
                          ncol_local=ncol_local_A, &
                          row_indices=row_indices_A, &
                          col_indices=col_indices_A, &
                          nrow_block=nrow_block_A, &
                          ncol_block=ncol_block_A)

      nprow = fm_mat_A%matrix_struct%context%num_pe(1)
      npcol = fm_mat_A%matrix_struct%context%num_pe(2)

      ALLOCATE (num_entries_rec(0:para_env_RPA%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env_RPA%num_pe - 1))

      num_entries_rec(:) = 0
      num_entries_send(:) = 0

      dummy_proc = 0

      CALL cp_fm_get_info(matrix=fm_mat_W, &
                          nrow_local=nrow_local_W, &
                          ncol_local=ncol_local_W, &
                          row_indices=row_indices_W, &
                          col_indices=col_indices_W, &
                          nrow_block=nrow_block_W, &
                          ncol_block=ncol_block_W)

      !MG to del: cf. mp2_r_2c.F and rpa_communication.F
      ! Use scalapack wrapper to find process index in A, where ijab -> iajb is needed
      ! To that end, we obtain the global index in A from the level indices
      DO ij_loc = 1, nrow_local_W
         i_occ_row = (row_indices_W(ij_loc) - 1)/virtual + 1
         j_occ_col = MOD(row_indices_W(ij_loc) - 1, virtual) + 1

         DO ab_loc = 1, ncol_local_W
            a_virt_row = (col_indices_W(ab_loc) - 1)/virtual + 1
            b_virt_col = MOD(col_indices_W(ab_loc) - 1, virtual) + 1

            send_prow = cp_fm_indxg2p(i_occ_row*a_virt_row, nrow_block_A, dummy_proc, &
                                      fm_mat_A%matrix_struct%first_p_pos(1), nprow)
            send_pcol = cp_fm_indxg2p(j_occ_col*b_virt_col, ncol_block_A, dummy_proc, &
                                      fm_mat_A%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_A%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL para_env_RPA%alltoall(num_entries_send, num_entries_rec, 1)

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_RPA%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_RPA%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_RPA%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_RPA%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_RPA%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
         buffer_rec(iproc)%indx = 0

      END DO

      DO iproc = 0, para_env_RPA%num_pe - 1

         ALLOCATE (buffer_send(iproc)%indx(num_entries_send(iproc), 2))
         buffer_send(iproc)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env_RPA%num_pe - 1))
      entry_counter(:) = 0

      ! Now we can write the actual data and indices to the send-buffer
      DO ij_loc = 1, nrow_local_W
         i_occ_row = (row_indices_W(ij_loc) - 1)/virtual + 1
         j_occ_col = MOD(row_indices_W(ij_loc) - 1, virtual) + 1

         DO ab_loc = 1, ncol_local_W
            a_virt_row = (col_indices_W(ab_loc) - 1)/virtual + 1
            b_virt_col = MOD(col_indices_W(ab_loc) - 1, virtual) + 1

            send_prow = cp_fm_indxg2p(i_occ_row*a_virt_row, nrow_block_A, dummy_proc, &
                                      fm_mat_A%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(j_occ_col*b_virt_col, ncol_block_A, dummy_proc, &
                                      fm_mat_A%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_A%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               fm_mat_W%local_data(ij_loc, ab_loc)

            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = i_occ_row*a_virt_row
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = j_occ_col*b_virt_col
         END DO
      END DO

      ALLOCATE (req_array(1:para_env_RPA%num_pe, 4))
      ! communicate the buffer
      CALL communicate_buffer(para_env_RPA, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      ! fill fm_mat_A with the entries from buffer_rec, i.e. buffer_rec are parts of W_ijab
      DO iproc = 0, para_env_RPA%num_pe - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            DO ii = 1, nrow_local_A

               i_row_global = row_indices_A(ii)

               DO jj = 1, ncol_local_A

                  j_col_global = col_indices_A(jj)

                  IF (i_row_global == buffer_rec(iproc)%indx(i_entry_rec, 1) .AND. &
                      j_col_global == buffer_rec(iproc)%indx(i_entry_rec, 2)) THEN
                     !Subtract buffer, since A = Delta E + alpha V - W
                     fm_mat_A%local_data(ii, jj) = fm_mat_A%local_data(ii, jj) - buffer_rec(iproc)%msg(i_entry_rec)

                  END IF
               END DO
            END DO
         END DO
      END DO

      !Clean up all the arrays from the communication process W_ijab -> A_iajb
      DO iproc = 0, para_env_RPA%num_pe - 1
         DEALLOCATE (buffer_rec(iproc)%msg)
         DEALLOCATE (buffer_rec(iproc)%indx)
         DEALLOCATE (buffer_send(iproc)%msg)
         DEALLOCATE (buffer_send(iproc)%indx)
      END DO
      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (req_array)
      DEALLOCATE (entry_counter)
      DEALLOCATE (num_entries_rec, num_entries_send)

      !full matrix W is not needed anymore, release it to save memory
      CALL cp_fm_struct_release(fm_struct_W)
      CALL cp_fm_release(fm_mat_W)

      !Now add the energy differences to A
      DO ii = 1, nrow_local_A

         i_row_global = row_indices_A(ii)

         DO jj = 1, ncol_local_A

            j_col_global = col_indices_A(jj)

            IF (i_row_global == j_col_global) THEN
               i_occ_row = (i_row_global - 1)/virtual + 1
               a_virt_row = MOD(i_row_global - 1, virtual) + 1
               eigen_diff = Eigenval(a_virt_row + homo) - Eigenval(i_occ_row)
               fm_mat_A%local_data(ii, jj) = fm_mat_A%local_data(ii, jj) + eigen_diff

            END IF
         END DO
      END DO

      !We have now the full matrix A_iajb, distributed over all ranks
      !Now: Diagonalize it
      CALL cp_fm_struct_create(fmstruct=fm_struct_eigvec, &
                               context=blacs_env_BSE, &
                               nrow_global=homo*virtual, &
                               ncol_global=homo*virtual)

      CALL cp_fm_create(fm_mat_eigvec, fm_struct_eigvec)

      ALLOCATE (Exc_ens(homo*virtual))
      CALL cp_fm_syevd(fm_mat_A, fm_mat_eigvec, Exc_ens)

      CALL success_message(Exc_ens, bse_spin_config, unit_nr, homo*virtual)

      !Release all local matrices and
      CALL cp_fm_struct_release(fm_struct_A)
      CALL cp_fm_release(fm_mat_A)
      CALL cp_fm_struct_release(fm_struct_eigvec)
      CALL cp_fm_release(fm_mat_eigvec)
      DEALLOCATE (Exc_ens)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Exc_ens ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param num_exc_en ...
! **************************************************************************************************
   SUBROUTINE success_message(Exc_ens, bse_spin_config, unit_nr, num_exc_en)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      INTEGER                                            :: bse_spin_config, unit_nr, num_exc_en

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'success_message'

      CHARACTER(LEN=10)                                  :: multiplet
      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: alpha

      CALL timeset(routineN, handle)

      !Prepare variables for printing
      IF (bse_spin_config == 0) THEN
         multiplet = "Singlet"
         alpha = 2.0_dp
      ELSE
         multiplet = "Triplet"
         alpha = 0.0_dp
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A)') '**                        BSE-TDA EXCITONIC ENERGIES                        **'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' The excitation energies are calculated by fully diagonalizing: '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') '    A_iajb   =  (E_a-E_i) delta_ij delta_ab   +  alpha * v_iajb   -  W_ijab   '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A48,A7,A12,F3.1)') &
            ' The spin-dependent factor is for the requested ', multiplet, " is alpha = ", alpha
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A16,T50,A22)') &
            ' Excitonic level', 'Excitation energy (eV)'
         !prints actual energies values
         DO i = 1, num_exc_en
            WRITE (unit_nr, '(T3,I16,T50,F22.3)') i, Exc_ens(i)*evolt
         END DO
         WRITE (unit_nr, '(T3,A)') ' '
      END IF

      CALL timestop(handle)
   END SUBROUTINE

END MODULE bse_full_diag
