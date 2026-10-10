!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Redistribute matrices between full-matrix and tensor process layouts.
! **************************************************************************************************
MODULE dbt_matrix_redistribution
   USE cp_dbcsr_api,                    ONLY: dbcsr_filter,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_dbcsr_redistribution,         ONLY: cp_dbcsr_from_subgroups,&
                                              cp_dbcsr_to_subgroups
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_copy,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_matrix_redistribution'
   PUBLIC :: cp_fm_to_dbt, cp_dbt_to_dbcsr, cp_dbt_to_fm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_global ...
!> \param mat_global ...
!> \param mat_local ...
!> \param tensor ...
!> \param para_env communicator covering the full matrix and all subgroups
!> \param num_pe_sub number of ranks per subgroup
!> \param eps_filter DBCSR filtering threshold before subgroup distribution
!> \param block_row_ranges ...
! **************************************************************************************************
   SUBROUTINE cp_fm_to_dbt(fm_global, mat_global, mat_local, tensor, para_env, num_pe_sub, eps_filter, block_row_ranges)

      TYPE(cp_fm_type)                                   :: fm_global
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(dbt_type)                                     :: tensor
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: num_pe_sub
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, DIMENSION(:, :), OPTIONAL                 :: block_row_ranges

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_fm_to_dbt'

      INTEGER                                            :: handle
      TYPE(dbt_type)                                     :: tensor_tmp

      CALL timeset(routineN, handle)

      CALL dbt_clear(tensor)
      CALL copy_fm_to_dbcsr(fm_global, mat_global, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_global, eps_filter)
      CALL cp_dbcsr_to_subgroups(mat_global, mat_local, para_env, num_pe_sub, block_row_ranges)
      CALL dbt_create(mat_local, tensor_tmp)
      CALL dbt_copy_matrix_to_tensor(mat_local, tensor_tmp)
      CALL dbt_copy(tensor_tmp, tensor, move_data=.TRUE.)
      CALL dbt_destroy(tensor_tmp)
      CALL dbcsr_set(mat_local, 0.0_dp)
      CALL dbcsr_filter(mat_local, 1.0_dp)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_to_dbt

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param mat_tensor ...
!> \param mat_global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cp_dbt_to_dbcsr(tensor, mat_tensor, mat_global, para_env)

      TYPE(dbt_type)                                     :: tensor
      TYPE(dbcsr_type)                                   :: mat_tensor, mat_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_dbt_to_dbcsr'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_copy_tensor_to_matrix(tensor, mat_tensor)
      CALL dbt_clear(tensor)
      ! the next para_env%sync is not mandatory, but it makes the timing output
      ! of cp_dbcsr_from_subgroups correct
      CALL para_env%sync()
      CALL cp_dbcsr_from_subgroups(mat_tensor, mat_global, para_env)

      CALL timestop(handle)

   END SUBROUTINE cp_dbt_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param t_R ...
!> \param fm_R ...
!> \param mat_global ...
!> \param mat_local ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cp_dbt_to_fm(t_R, fm_R, mat_global, mat_local, para_env)
      TYPE(dbt_type), DIMENSION(:)                       :: t_R
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_R
      TYPE(dbcsr_p_type)                                 :: mat_global, mat_local
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_dbt_to_fm'

      INTEGER                                            :: handle, i_cell, n_images

      CALL timeset(routineN, handle)

      n_images = SIZE(t_R)

      CPASSERT(n_images == SIZE(fm_R))

      DO i_cell = 1, n_images
         CALL cp_dbt_to_dbcsr(t_R(i_cell), mat_local%matrix, mat_global%matrix, &
                              para_env)
         CALL copy_dbcsr_to_fm(mat_global%matrix, fm_R(i_cell))
      END DO

      CALL timestop(handle)

   END SUBROUTINE cp_dbt_to_fm

END MODULE dbt_matrix_redistribution
