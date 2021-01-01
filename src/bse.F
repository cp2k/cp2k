!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2017 created [Jan Wilhelm]
! **************************************************************************************************
MODULE bse
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE group_dist_types,                ONLY: get_group_dist,&
                                              group_dist_d1_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_sum
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse'

   PUBLIC :: mult_B_with_W_and_fill_local_3c_arrays, do_subspace_iterations

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param max_iter ...
!> \param threshold_min_trans ...
!> \param Eigenval ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo, virtual, num_Z_vectors, &
                                     max_iter, threshold_min_trans, Eigenval, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER                                            :: homo, virtual, num_Z_vectors, max_iter
      REAL(KIND=dp)                                      :: threshold_min_trans
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'do_subspace_iterations'

      INTEGER                                            :: handle, i_iter, local_RI_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_ia_tmp, M_ji_tmp, RI_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, BZ, Z_vectors

      CALL timeset(routineN, handle)

      ! JW hack 2del
      threshold_min_trans = 0.01_dp

      ALLOCATE (Z_vectors(homo, virtual, num_Z_vectors))
      Z_vectors = 0.0_dp

      ALLOCATE (AZ(homo, virtual, num_Z_vectors))
      AZ = 0.0_dp

      ALLOCATE (BZ(homo, virtual, num_Z_vectors))
      BZ = 0.0_dp

      local_RI_size = SIZE(B_iaQ_bse_local, 3)

      ALLOCATE (M_ia_tmp(homo, virtual))
      M_ia_tmp = 0.0_dp

      ALLOCATE (M_ji_tmp(homo, homo))
      M_ji_tmp = 0.0_dp

      ALLOCATE (RI_vector(local_RI_size, num_Z_vectors))
      RI_vector = 0.0_dp

      CALL initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)

      DO i_iter = 1, max_iter

         CALL compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, &
                         M_ia_tmp, RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, &
                         para_env)

         CALL compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
                         M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, &
                         para_env)

      END DO

      DEALLOCATE (AZ, BZ, Z_vectors, M_ia_tmp, M_ji_tmp, RI_vector)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param BZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param M_ji_tmp ...
!> \param RI_vector ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_BZ(BZ, Z_vectors, B_iaQ_bse_local, B_bar_iaQ_bse_local, &
                         M_ji_tmp, RI_vector, homo, virtual, num_Z_vectors, local_RI_size, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ji_tmp, RI_vector
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i_Z_vector, LLL

      BZ(:, :, :) = 0.0_dp

      CALL compute_v_ia_jb_part(BZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                num_Z_vectors, homo, virtual)

      DO i_Z_vector = 1, num_Z_vectors

         DO LLL = 1, local_RI_size

            ! M_ji^P = sum_b Z_jb*B_bi^P
            CALL DGEMM("N", "T", homo, homo, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_iaQ_bse_local(:, :, LLL), homo, 0.0_dp, M_ji_tmp, homo)
            ! (BZ)_ia = sum_jP M_ij^P*B^bar_ja^P
            CALL DGEMM("T", "N", homo, virtual, homo, 1.0_dp, M_ji_tmp, homo, &
                       B_bar_iaQ_bse_local, homo, 1.0_dp, BZ(:, :, i_Z_vector), homo)

         END DO

      END DO

      ! we make the mp_sum to sum over all RI basis functions
      CALL mp_sum(BZ, para_env%group)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param M_ia_tmp ...
!> \param RI_vector ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param num_Z_vectors ...
!> \param local_RI_size ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_AZ(AZ, Z_vectors, B_iaQ_bse_local, B_bar_ijQ_bse_local, B_abQ_bse_local, M_ia_tmp, &
                         RI_vector, Eigenval, homo, virtual, num_Z_vectors, local_RI_size, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, Z_vectors, B_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_abQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: M_ia_tmp, RI_vector
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo, virtual, num_Z_vectors, &
                                                            local_RI_size
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_AZ'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, LLL
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      AZ(:, :, :) = 0.0_dp

      CALL compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                num_Z_vectors, homo, virtual)

      DO i_Z_vector = 1, num_Z_vectors

         ! JW TO DO: OMP PARALLELIZATION
         DO LLL = 1, local_RI_size

            ! M_ja^P = sum_j Z_jb*B_ba^P
            CALL DGEMM("N", "N", homo, virtual, virtual, 1.0_dp, Z_vectors(:, :, i_Z_vector), homo, &
                       B_abQ_bse_local(:, :, LLL), virtual, 0.0_dp, M_ia_tmp, homo)

            ! (AZ)_ia = sum_jP B_bar_ij^P*M_ja^P
            CALL DGEMM("N", "N", homo, virtual, homo, 1.0_dp, B_bar_ijQ_bse_local(:, :, LLL), homo, &
                       M_ia_tmp, homo, 1.0_dp, AZ(:, :, i_Z_vector), homo)

         END DO

      END DO

      ! we make the mp_sum to sum over all RI basis functions
      CALL mp_sum(AZ, para_env%group)

      ! add (e_a-e_i)*Z_ia
      DO i_occ = 1, homo
         DO a_virt = 1, virtual

            eigen_diff = Eigenval(a_virt + homo) - Eigenval(i_occ)

            AZ(i_occ, a_virt, :) = AZ(i_occ, a_virt, :) + Z_vectors(i_occ, a_virt, :)*eigen_diff

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param AZ ...
!> \param Z_vectors ...
!> \param B_iaQ_bse_local ...
!> \param RI_vector ...
!> \param local_RI_size ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE compute_v_ia_jb_part(AZ, Z_vectors, B_iaQ_bse_local, RI_vector, local_RI_size, &
                                   num_Z_vectors, homo, virtual)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: AZ, Z_vectors, B_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:, :)                     :: RI_vector
      INTEGER                                            :: local_RI_size, num_Z_vectors, homo, &
                                                            virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_v_ia_jb_part'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, LLL

      CALL timeset(routineN, handle)

      RI_vector = 0.0_dp

      ! v_P = sum_jb B_jb^P Z_jb
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual

                  RI_vector(LLL, i_Z_vector) = RI_vector(LLL, i_Z_vector) + &
                                               Z_vectors(i_occ, a_virt, i_Z_vector)* &
                                               B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO

      ! AZ = sum_P B_ia^P*v_P + ...
      DO LLL = 1, local_RI_size
         DO i_Z_vector = 1, num_Z_vectors
            DO i_occ = 1, homo
               DO a_virt = 1, virtual

                  AZ(i_occ, a_virt, i_Z_vector) = AZ(i_occ, a_virt, i_Z_vector) + &
                                                  RI_vector(LLL, i_Z_vector)* &
                                                  B_iaQ_bse_local(i_occ, a_virt, LLL)

               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Z_vectors ...
!> \param Eigenval ...
!> \param num_Z_vectors ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE initial_guess_Z_vectors(Z_vectors, Eigenval, num_Z_vectors, homo, virtual)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Z_vectors
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: num_Z_vectors, homo, virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initial_guess_Z_vectors'

      INTEGER                                            :: a_virt, handle, i_occ, i_Z_vector, &
                                                            min_loc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigen_diff_ia

      CALL timeset(routineN, handle)

      ALLOCATE (eigen_diff_ia(homo, virtual))

      DO i_occ = 1, homo
         DO a_virt = 1, virtual
            eigen_diff_ia(i_occ, a_virt) = Eigenval(a_virt + homo) - Eigenval(i_occ)
         END DO
      END DO

      DO i_Z_vector = 1, num_Z_vectors

         min_loc = MINLOC(eigen_diff_ia)

         Z_vectors(min_loc(1), min_loc(2), i_Z_vector) = 1.0_dp

         eigen_diff_ia(min_loc(1), min_loc(2)) = 1.0E20_dp

      END DO

      DEALLOCATE (eigen_diff_ia)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param B_bar_ijQ_bse_local ...
!> \param B_abQ_bse_local ...
!> \param B_bar_iaQ_bse_local ...
!> \param B_iaQ_bse_local ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
!> \param gd_array ...
!> \param color_sub ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE mult_B_with_W_and_fill_local_3c_arrays(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S, fm_mat_Q_static_bse, &
                                                     fm_mat_Q_static_bse_gemm, &
                                                     B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                     B_iaQ_bse_local, dimen_RI, homo, virtual, dimen_ia, &
                                                     gd_array, color_sub, para_env)

      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                                            fm_mat_S, fm_mat_Q_static_bse, &
                                                            fm_mat_Q_static_bse_gemm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                            B_bar_iaQ_bse_local, B_iaQ_bse_local
      INTEGER                                            :: dimen_RI, homo, virtual, dimen_ia
      TYPE(group_dist_d1_type)                           :: gd_array
      INTEGER                                            :: color_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mult_B_with_W_and_fill_local_3c_arrays'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_bar_ia_bse, &
                                                            fm_mat_S_bar_ij_bse, fm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_S_bar_ia_bse, fm_mat_S%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S, fm_mat_S_bar_ia_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ia_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_S_bar_ij_bse, fm_mat_S_ij_bse%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_S_ij_bse, fm_mat_S_bar_ij_bse)
      CALL cp_fm_set_all(fm_mat_S_bar_ij_bse, 0.0_dp)

      CALL cp_fm_create(fm_mat_work, fm_mat_Q_static_bse_gemm%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_Q_static_bse_gemm, fm_mat_work)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

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
      CALL cp_fm_upper_to_full(fm_mat_Q_static_bse_gemm, fm_mat_work)

      CALL cp_gemm(transa="N", transb="N", m=homo**2, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S_ij_bse, matrix_b=fm_mat_Q_static_bse, beta=0.0_dp, &
                   matrix_c=fm_mat_S_bar_ij_bse)

      ! fm_mat_S_bar_ia_bse has a different blacs_env as fm_mat_S_ij_bse since we take
      ! fm_mat_S from RPA. Therefore, we also need a different fm_mat_Q_static_bse_gemm
      CALL cp_gemm(transa="N", transb="N", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S, matrix_b=fm_mat_Q_static_bse_gemm, beta=0.0_dp, &
                   matrix_c=fm_mat_S_bar_ia_bse)

      CALL allocate_and_fill_local_array(B_iaQ_bse_local, fm_mat_S, gd_array, color_sub, homo, virtual, dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_iaQ_bse_local, fm_mat_S_bar_ia_bse, gd_array, color_sub, homo, virtual, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_bar_ijQ_bse_local, fm_mat_S_bar_ij_bse, gd_array, color_sub, homo, homo, &
                                         dimen_RI, para_env)

      CALL allocate_and_fill_local_array(B_abQ_bse_local, fm_mat_S_ab_bse, gd_array, color_sub, virtual, virtual, &
                                         dimen_RI, para_env)

      CALL cp_fm_release(fm_mat_S_bar_ia_bse)
      CALL cp_fm_release(fm_mat_S_bar_ij_bse)
      CALL cp_fm_release(fm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param B_local ...
!> \param fm_mat_S ...
!> \param gd_array ...
!> \param color_sub ...
!> \param small_size ...
!> \param big_size ...
!> \param dimen_RI ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_local_array(B_local, fm_mat_S, gd_array, &
                                            color_sub, small_size, big_size, dimen_RI, para_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_local
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      TYPE(group_dist_d1_type)                           :: gd_array
      INTEGER                                            :: color_sub, small_size, big_size, dimen_RI
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_fill_local_array'

      INTEGER :: combi_index, end_RI, handle, handle1, i_comm, i_entry, iiB, imepos, jjB, &
         level_big_size, level_small_size, ncol_local, nrow_local, num_comm_cycles, RI_index, &
         size_RI, start_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, mepos_from_RI_index, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp)                                      :: matrix_el
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      ALLOCATE (mepos_from_RI_index(dimen_RI))
      mepos_from_RI_index = 0

      DO imepos = 0, para_env%num_pe - 1

         CALL get_group_dist(gd_array, pos=imepos, starts=start_RI, ends=end_RI)

         mepos_from_RI_index(start_RI:end_RI) = imepos

      END DO

      ! color_sub is automatically the number of the process since every subgroup has only one MPI rank
      CALL get_group_dist(gd_array, color_sub, start_RI, end_RI, size_RI)

      ALLOCATE (B_local(small_size, big_size, 1:size_RI))

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))

      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      num_comm_cycles = 10

      ! communicate not all due to huge memory overhead, since for every number in fm_mat_S, we store
      ! three additional ones (RI index, first MO index, second MO index!!)
      DO i_comm = 0, num_comm_cycles - 1

         num_entries_send = 0
         num_entries_rec = 0

         ! loop over RI index to get the number of sent entries
         DO jjB = 1, ncol_local

            RI_index = col_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            num_entries_send(imepos) = num_entries_send(imepos) + nrow_local

         END DO

         CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

         ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
         ALLOCATE (buffer_send(0:para_env%num_pe - 1))

         ! allocate data message and corresponding indices
         DO imepos = 0, para_env%num_pe - 1

            ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
            buffer_rec(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
            buffer_send(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
            buffer_rec(imepos)%indx = 0

            ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
            buffer_send(imepos)%indx = 0

         END DO

         entry_counter(:) = 0

         ! loop over RI index for filling the send-buffer
         DO jjB = 1, ncol_local

            RI_index = col_indices(jjB)

            IF (MODULO(RI_index, num_comm_cycles) /= i_comm) CYCLE

            imepos = mepos_from_RI_index(RI_index)

            DO iiB = 1, nrow_local

               combi_index = row_indices(iiB)
               level_small_size = MAX(1, combi_index - 1)/big_size + 1
               level_big_size = combi_index - (level_small_size - 1)*big_size

               entry_counter(imepos) = entry_counter(imepos) + 1

               buffer_send(imepos)%msg(entry_counter(imepos)) = fm_mat_S%local_data(iiB, jjB)

               buffer_send(imepos)%indx(entry_counter(imepos), 1) = RI_index
               buffer_send(imepos)%indx(entry_counter(imepos), 2) = level_small_size
               buffer_send(imepos)%indx(entry_counter(imepos), 3) = level_big_size

            END DO

         END DO

         CALL timeset("BSE_comm_data", handle1)

         CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

         CALL timestop(handle1)

         ! fill B_local
         DO imepos = 0, para_env%num_pe - 1

            DO i_entry = 1, num_entries_rec(imepos)

               RI_index = buffer_rec(imepos)%indx(i_entry, 1) - start_RI + 1
               level_small_size = buffer_rec(imepos)%indx(i_entry, 2)
               level_big_size = buffer_rec(imepos)%indx(i_entry, 3)

               matrix_el = buffer_rec(imepos)%msg(i_entry)

               B_local(level_small_size, level_big_size, RI_index) = matrix_el

            END DO

         END DO

         DO imepos = 0, para_env%num_pe - 1
            DEALLOCATE (buffer_send(imepos)%msg)
            DEALLOCATE (buffer_send(imepos)%indx)
            DEALLOCATE (buffer_rec(imepos)%msg)
            DEALLOCATE (buffer_rec(imepos)%indx)
         END DO

         DEALLOCATE (buffer_rec, buffer_send)

      END DO

      DEALLOCATE (num_entries_send, num_entries_rec)

      DEALLOCATE (mepos_from_RI_index)

      DEALLOCATE (entry_counter, req_array)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE bse
