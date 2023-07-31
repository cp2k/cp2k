!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Interface to (sca)lapack for the Cholesky based procedures
!> \author  VW
!> \date    2009-09-08
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-09-08
! **************************************************************************************************
MODULE cp_dbcsr_cholesky
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_cholesky_restore,&
                                              cp_fm_potrf,&
                                              cp_fm_potri,&
                                              cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_get_info,&
                                              dbcsr_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_cholesky'

   PUBLIC :: cp_dbcsr_cholesky_decompose, cp_dbcsr_cholesky_invert, &
             cp_dbcsr_cholesky_restore

   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief used to replace a symmetric positive def. matrix M with its cholesky
!>      decomposition U: M = U^T * U, with U upper triangular
!> \param matrix the matrix to replace with its cholesky decomposition
!> \param n the number of row (and columns) of the matrix &
!>        (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \par History
!>      05.2002 created [JVdV]
!>      12.2002 updated, added n optional parm [fawzi]
!> \author Joost
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_cholesky_decompose(matrix, n, para_env, blacs_env)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_cholesky_decompose'

      INTEGER                                            :: handle, my_n, nfullcols_total, &
                                                            nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_matrix

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)

      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix, fm_struct, name="fm_matrix")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fm_matrix)

      my_n = MIN(fm_matrix%matrix_struct%nrow_global, &
                 fm_matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      CALL cp_fm_potrf(fm_matrix, my_n)

      CALL copy_fm_to_dbcsr(fm_matrix, matrix)

      CALL cp_fm_release(fm_matrix)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_cholesky_decompose

! **************************************************************************************************
!> \brief used to replace the cholesky decomposition by the inverse
!> \param matrix the matrix to invert (must be an upper triangular matrix)
!> \param n size of the matrix to invert (defaults to the min(size(matrix)))
!> \param para_env ...
!> \param blacs_env ...
!> \param upper_to_full ...
!> \par History
!>      05.2002 created [JVdV]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_cholesky_invert(matrix, n, para_env, blacs_env, upper_to_full)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      LOGICAL, INTENT(IN)                                :: upper_to_full

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_cholesky_invert'

      INTEGER                                            :: handle, my_n, nfullcols_total, &
                                                            nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_matrix, fm_matrix_tmp

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)

      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullrows_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix, fm_struct, name="fm_matrix")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fm_matrix)

      my_n = MIN(fm_matrix%matrix_struct%nrow_global, &
                 fm_matrix%matrix_struct%ncol_global)
      IF (PRESENT(n)) THEN
         CPASSERT(n <= my_n)
         my_n = n
      END IF

      CALL cp_fm_potri(fm_matrix, my_n)

      IF (upper_to_full) THEN
         CALL cp_fm_create(fm_matrix_tmp, fm_matrix%matrix_struct, name="fm_matrix_tmp")
         CALL cp_fm_upper_to_full(fm_matrix, fm_matrix_tmp)
         CALL cp_fm_release(fm_matrix_tmp)
      END IF

      CALL copy_fm_to_dbcsr(fm_matrix, matrix)

      CALL cp_fm_release(fm_matrix)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_cholesky_invert

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param neig ...
!> \param matrixb ...
!> \param matrixout ...
!> \param op ...
!> \param pos ...
!> \param transa ...
!> \param para_env ...
!> \param blacs_env ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_cholesky_restore(matrix, neig, matrixb, matrixout, op, pos, transa, &
                                        para_env, blacs_env)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER, INTENT(IN)                                :: neig
      TYPE(dbcsr_type)                                   :: matrixb, matrixout
      CHARACTER(LEN=*), INTENT(IN)                       :: op
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: pos, transa
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_cholesky_restore'

      CHARACTER                                          :: chol_pos, chol_transa
      INTEGER                                            :: handle, nfullcols_total, nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_matrix, fm_matrixb, fm_matrixout

      CALL timeset(routineN, handle)

      NULLIFY (fm_struct)

      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix, fm_struct, name="fm_matrix")
      CALL cp_fm_struct_release(fm_struct)

      CALL dbcsr_get_info(matrixb, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrixb, fm_struct, name="fm_matrixb")
      CALL cp_fm_struct_release(fm_struct)

      CALL dbcsr_get_info(matrixout, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrixout, fm_struct, name="fm_matrixout")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fm_matrix)
      CALL copy_dbcsr_to_fm(matrixb, fm_matrixb)
      !CALL copy_dbcsr_to_fm(matrixout, fm_matrixout)

      IF (op /= "SOLVE" .AND. op /= "MULTIPLY") &
         CPABORT("wrong argument op")

      IF (PRESENT(pos)) THEN
         SELECT CASE (pos)
         CASE ("LEFT")
            chol_pos = 'L'
         CASE ("RIGHT")
            chol_pos = 'R'
         CASE DEFAULT
            CPABORT("wrong argument pos")
         END SELECT
      ELSE
         chol_pos = 'L'
      END IF

      chol_transa = 'N'
      IF (PRESENT(transa)) chol_transa = transa

      IF ((fm_matrix%use_sp .NEQV. fm_matrixb%use_sp) .OR. (fm_matrix%use_sp .NEQV. fm_matrixout%use_sp)) &
         CPABORT("not the same precision")

      CALL cp_fm_cholesky_restore(fm_matrix, neig, fm_matrixb, fm_matrixout, op, chol_pos, chol_transa)

      CALL copy_fm_to_dbcsr(fm_matrixout, matrixout)

      CALL cp_fm_release(fm_matrix)
      CALL cp_fm_release(fm_matrixb)
      CALL cp_fm_release(fm_matrixout)

      CALL timestop(handle)

   END SUBROUTINE cp_dbcsr_cholesky_restore

END MODULE cp_dbcsr_cholesky

