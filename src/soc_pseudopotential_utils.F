!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE soc_pseudopotential_utils
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'soc_pseudopotential_utils'

   PUBLIC :: add_dbcsr_submat, cfm_add_on_diag, add_fm_submat, get_cfm_submat, create_cfm_double, &
             add_cfm_submat

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_target ...
!> \param mat_source ...
!> \param fm_struct_source ...
!> \param nstart_row ...
!> \param nstart_col ...
!> \param factor ...
!> \param add_also_herm_conj ...
! **************************************************************************************************
   SUBROUTINE add_dbcsr_submat(cfm_mat_target, mat_source, fm_struct_source, &
                               nstart_row, nstart_col, factor, add_also_herm_conj)
      TYPE(cp_cfm_type)                                  :: cfm_mat_target
      TYPE(dbcsr_type)                                   :: mat_source
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_source
      INTEGER                                            :: nstart_row, nstart_col
      COMPLEX(KIND=dp)                                   :: factor
      LOGICAL                                            :: add_also_herm_conj

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_dbcsr_submat'

      INTEGER                                            :: handle, nao
      TYPE(cp_cfm_type)                                  :: cfm_mat_work_double, &
                                                            cfm_mat_work_double_2
      TYPE(cp_fm_type)                                   :: fm_mat_work_double_im, fm_mat_work_im

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_work_double_im, cfm_mat_target%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_double_im, 0.0_dp)

      CALL cp_cfm_create(cfm_mat_work_double, cfm_mat_target%matrix_struct)
      CALL cp_cfm_create(cfm_mat_work_double_2, cfm_mat_target%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work_double, z_zero)
      CALL cp_cfm_set_all(cfm_mat_work_double_2, z_zero)

      CALL cp_fm_create(fm_mat_work_im, fm_struct_source)

      CALL copy_dbcsr_to_fm(mat_source, fm_mat_work_im)

      CALL cp_fm_get_info(fm_mat_work_im, nrow_global=nao)

      CALL cp_fm_to_fm_submat(msource=fm_mat_work_im, mtarget=fm_mat_work_double_im, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=1, s_firstcol=1, &
                              t_firstrow=nstart_row, t_firstcol=nstart_col)

      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_work_double, gaussi, fm_mat_work_double_im)

      CALL cp_cfm_scale(factor, cfm_mat_work_double)

      CALL cp_cfm_scale_and_add(z_one, cfm_mat_target, z_one, cfm_mat_work_double)

      IF (add_also_herm_conj) THEN
         CALL cp_cfm_transpose(cfm_mat_work_double, 'C', cfm_mat_work_double_2)
         CALL cp_cfm_scale_and_add(z_one, cfm_mat_target, z_one, cfm_mat_work_double_2)
      END IF

      CALL cp_fm_release(fm_mat_work_double_im)
      CALL cp_cfm_release(cfm_mat_work_double)
      CALL cp_cfm_release(cfm_mat_work_double_2)
      CALL cp_fm_release(fm_mat_work_im)

      CALL timestop(handle)

   END SUBROUTINE add_dbcsr_submat

! **************************************************************************************************
!> \brief ...
!> \param cfm ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE cfm_add_on_diag(cfm, alpha)

      TYPE(cp_cfm_type)                                  :: cfm
      REAL(KIND=dp), DIMENSION(:)                        :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cfm_add_on_diag'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, nao, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      nao = SIZE(alpha)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            i_global = row_indices(i_row)
            IF (j_global == i_global) THEN
               IF (i_global .LE. nao) THEN
                  cfm%local_data(i_row, j_col) = cfm%local_data(i_row, j_col) + alpha(i_global)*z_one
               ELSE
                  cfm%local_data(i_row, j_col) = cfm%local_data(i_row, j_col) + alpha(i_global - nao)*z_one
               END IF
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE cfm_add_on_diag

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_target ...
!> \param fm_mat_source ...
!> \param nstart_row ...
!> \param nstart_col ...
! **************************************************************************************************
   SUBROUTINE add_fm_submat(cfm_mat_target, fm_mat_source, nstart_row, nstart_col)

      TYPE(cp_cfm_type)                                  :: cfm_mat_target
      TYPE(cp_fm_type)                                   :: fm_mat_source
      INTEGER                                            :: nstart_row, nstart_col

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_fm_submat'

      INTEGER                                            :: handle, nao
      TYPE(cp_fm_type)                                   :: fm_mat_work_double_re

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_work_double_re, cfm_mat_target%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_double_re, 0.0_dp)

      CALL cp_fm_get_info(fm_mat_source, nrow_global=nao)

      CALL cp_fm_to_fm_submat(msource=fm_mat_source, mtarget=fm_mat_work_double_re, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=1, s_firstcol=1, &
                              t_firstrow=nstart_row, t_firstcol=nstart_col)

      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_target, z_one, fm_mat_work_double_re)

      CALL cp_fm_release(fm_mat_work_double_re)

      CALL timestop(handle)

   END SUBROUTINE add_fm_submat

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_target ...
!> \param cfm_mat_source ...
!> \param nstart_row ...
!> \param nstart_col ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE add_cfm_submat(cfm_mat_target, cfm_mat_source, nstart_row, nstart_col, factor)

      TYPE(cp_cfm_type)                                  :: cfm_mat_target, cfm_mat_source
      INTEGER                                            :: nstart_row, nstart_col
      COMPLEX(KIND=dp), OPTIONAL                         :: factor

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_cfm_submat'

      COMPLEX(KIND=dp)                                   :: factor_im, factor_re
      INTEGER                                            :: handle, nao
      TYPE(cp_fm_type)                                   :: fm_mat_source_im, fm_mat_source_re, &
                                                            fm_mat_work_double_im, &
                                                            fm_mat_work_double_re

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_work_double_re, cfm_mat_target%matrix_struct)
      CALL cp_fm_create(fm_mat_work_double_im, cfm_mat_target%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_double_re, 0.0_dp)
      CALL cp_fm_set_all(fm_mat_work_double_im, 0.0_dp)

      CALL cp_fm_create(fm_mat_source_re, cfm_mat_source%matrix_struct)
      CALL cp_fm_create(fm_mat_source_im, cfm_mat_source%matrix_struct)
      CALL cp_cfm_to_fm(cfm_mat_source, fm_mat_source_re, fm_mat_source_im)

      CALL cp_cfm_get_info(cfm_mat_source, nrow_global=nao)

      CALL cp_fm_to_fm_submat(msource=fm_mat_source_re, mtarget=fm_mat_work_double_re, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=1, s_firstcol=1, &
                              t_firstrow=nstart_row, t_firstcol=nstart_col)

      CALL cp_fm_to_fm_submat(msource=fm_mat_source_im, mtarget=fm_mat_work_double_im, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=1, s_firstcol=1, &
                              t_firstrow=nstart_row, t_firstcol=nstart_col)

      IF (PRESENT(factor)) THEN
         factor_re = factor
         factor_im = gaussi*factor
      ELSE
         factor_re = z_one
         factor_im = gaussi
      END IF

      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_target, factor_re, fm_mat_work_double_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_target, factor_im, fm_mat_work_double_im)

      CALL cp_fm_release(fm_mat_work_double_re)
      CALL cp_fm_release(fm_mat_work_double_im)
      CALL cp_fm_release(fm_mat_source_re)
      CALL cp_fm_release(fm_mat_source_im)

      CALL timestop(handle)

   END SUBROUTINE add_cfm_submat

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_target ...
!> \param cfm_mat_source ...
!> \param nstart_row ...
!> \param nstart_col ...
! **************************************************************************************************
   SUBROUTINE get_cfm_submat(cfm_mat_target, cfm_mat_source, nstart_row, nstart_col)

      TYPE(cp_cfm_type)                                  :: cfm_mat_target, cfm_mat_source
      INTEGER                                            :: nstart_row, nstart_col

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_cfm_submat'

      INTEGER                                            :: handle, nao
      TYPE(cp_fm_type)                                   :: fm_mat_source_double_im, &
                                                            fm_mat_source_double_re, &
                                                            fm_mat_work_im, fm_mat_work_re

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mat_source_double_re, cfm_mat_source%matrix_struct)
      CALL cp_fm_create(fm_mat_source_double_im, cfm_mat_source%matrix_struct)
      CALL cp_cfm_to_fm(cfm_mat_source, fm_mat_source_double_re, fm_mat_source_double_im)

      CALL cp_fm_create(fm_mat_work_re, cfm_mat_target%matrix_struct)
      CALL cp_fm_create(fm_mat_work_im, cfm_mat_target%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_re, 0.0_dp)
      CALL cp_fm_set_all(fm_mat_work_im, 0.0_dp)

      CALL cp_cfm_get_info(cfm_mat_target, nrow_global=nao)

      CALL cp_fm_to_fm_submat(msource=fm_mat_source_double_re, mtarget=fm_mat_work_re, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=nstart_row, s_firstcol=nstart_col, &
                              t_firstrow=1, t_firstcol=1)

      CALL cp_fm_to_fm_submat(msource=fm_mat_source_double_im, mtarget=fm_mat_work_im, &
                              nrow=nao, ncol=nao, &
                              s_firstrow=nstart_row, s_firstcol=nstart_col, &
                              t_firstrow=1, t_firstcol=1)

      CALL cp_fm_to_cfm(fm_mat_work_re, fm_mat_work_im, cfm_mat_target)

      CALL cp_fm_release(fm_mat_work_re)
      CALL cp_fm_release(fm_mat_work_im)
      CALL cp_fm_release(fm_mat_source_double_re)
      CALL cp_fm_release(fm_mat_source_double_im)

      CALL timestop(handle)

   END SUBROUTINE get_cfm_submat

! **************************************************************************************************
!> \brief ...
!> \param fm_orig ...
!> \param cfm_double ...
! **************************************************************************************************
   SUBROUTINE create_cfm_double(fm_orig, cfm_double)
      TYPE(cp_fm_type)                                   :: fm_orig
      TYPE(cp_cfm_type)                                  :: cfm_double

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_cfm_double'

      INTEGER                                            :: handle, ncol_global_orig, &
                                                            nrow_global_orig
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_double

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_orig, nrow_global=nrow_global_orig, ncol_global=ncol_global_orig)

      CALL cp_fm_struct_create(fm_struct_double, &
                               nrow_global=2*nrow_global_orig, &
                               ncol_global=2*ncol_global_orig, &
                               template_fmstruct=fm_orig%matrix_struct)

      CALL cp_cfm_create(cfm_double, fm_struct_double)

      CALL cp_fm_struct_release(fm_struct_double)

      CALL timestop(handle)

   END SUBROUTINE create_cfm_double

END MODULE soc_pseudopotential_utils
