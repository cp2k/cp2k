!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 08.2023
! **************************************************************************************************
MODULE gw_utils_communication
   USE cp_dbcsr_api,                    ONLY: dbcsr_filter,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_dbcsr_redistribution,         ONLY: cp_dbcsr_to_subgroups,&
                                              cp_dbcsr_from_subgroups
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_copy,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_communication'

   PUBLIC :: local_dbt_to_global_mat, fm_to_local_tensor, local_dbt_to_global_fm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_global ...
!> \param mat_global ...
!> \param mat_local ...
!> \param tensor ...
!> \param bs_env Band-structure environment containing GW parameters.
!> \param atom_ranges ...
! **************************************************************************************************
   SUBROUTINE fm_to_local_tensor(fm_global, mat_global, mat_local, tensor, bs_env, atom_ranges)

      TYPE(cp_fm_type)                                   :: fm_global
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(dbt_type)                                     :: tensor
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:, :), OPTIONAL                 :: atom_ranges

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_to_local_tensor'

      INTEGER                                            :: handle
      TYPE(dbt_type)                                     :: tensor_tmp

      CALL timeset(routineN, handle)

      CALL dbt_clear(tensor)
      CALL copy_fm_to_dbcsr(fm_global, mat_global, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_global, bs_env%eps_filter)
      IF (PRESENT(atom_ranges)) THEN
         CALL cp_dbcsr_to_subgroups(mat_global, mat_local, bs_env%para_env, &
                                            bs_env%para_env_tensor%num_pe, atom_ranges)
      ELSE
         CALL cp_dbcsr_to_subgroups(mat_global, mat_local, bs_env%para_env, &
                                            bs_env%para_env_tensor%num_pe)
      END IF
      CALL dbt_create(mat_local, tensor_tmp)
      CALL dbt_copy_matrix_to_tensor(mat_local, tensor_tmp)
      CALL dbt_copy(tensor_tmp, tensor, move_data=.TRUE.)
      CALL dbt_destroy(tensor_tmp)
      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)

      CALL timestop(handle)

   END SUBROUTINE fm_to_local_tensor

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param mat_tensor ...
!> \param mat_global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE local_dbt_to_global_mat(tensor, mat_tensor, mat_global, para_env)

      TYPE(dbt_type)                                     :: tensor
      TYPE(dbcsr_type)                                   :: mat_tensor, mat_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_dbt_to_global_mat'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_copy_tensor_to_matrix(tensor, mat_tensor)
      CALL dbt_clear(tensor)
      ! the next para_env%sync is not mandatory, but it makes the timing output
      ! of cp_dbcsr_from_subgroups correct
      CALL para_env%sync()
      CALL cp_dbcsr_from_subgroups(mat_tensor, mat_global, para_env)

      CALL timestop(handle)

   END SUBROUTINE local_dbt_to_global_mat





! **************************************************************************************************
!> \brief ...
!> \param t_R ...
!> \param fm_R ...
!> \param mat_global ...
!> \param mat_local ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE local_dbt_to_global_fm(t_R, fm_R, mat_global, mat_local, bs_env)
      TYPE(dbt_type), DIMENSION(:)                       :: t_R
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_R
      TYPE(dbcsr_p_type)                                 :: mat_global, mat_local
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_dbt_to_global_fm'

      INTEGER                                            :: handle, i_cell, n_images

      CALL timeset(routineN, handle)

      n_images = SIZE(t_R)

      CPASSERT(n_images == SIZE(fm_R))

      DO i_cell = 1, n_images
         CALL dbcsr_set(mat_global%matrix, 0.0_dp)
         CALL dbcsr_set(mat_local%matrix, 0.0_dp)
         CALL local_dbt_to_global_mat(t_R(i_cell), mat_local%matrix, mat_global%matrix, &
                                      bs_env%para_env)
         CALL copy_dbcsr_to_fm(mat_global%matrix, fm_R(i_cell))
      END DO

      CALL timestop(handle)

   END SUBROUTINE local_dbt_to_global_fm

END MODULE gw_utils_communication
