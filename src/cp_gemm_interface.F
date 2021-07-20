!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief basic linear algebra operations for full matrixes
!> \par History
!>      08.2002 splitted out of qs_blacs [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_gemm_interface
   USE ISO_C_BINDING, ONLY: C_CHAR, &
                            C_DOUBLE, &
                            C_INT, &
                            C_LOC, &
                            C_PTR
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm_bc, &
                                  copy_fm_to_dbcsr_bc
   USE cp_fm_basic_linalg, ONLY: cp_fm_gemm
   USE cp_fm_types, ONLY: cp_fm_get_info, &
                          cp_fm_get_mm_type, &
                          cp_fm_type
   USE dbcsr_api, ONLY: dbcsr_multiply, &
                        dbcsr_release, &
                        dbcsr_type
   USE input_constants, ONLY: do_cosma, &
                              do_dbcsr, &
                              do_scalapack
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_min
   USE offload_api, ONLY: offload_set_device
   USE string_utilities, ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_gemm_interface'

   PUBLIC :: cp_gemm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param a_first_col ...
!> \param a_first_row ...
!> \param b_first_col ...
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
! **************************************************************************************************
   SUBROUTINE cp_gemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, &
                      matrix_c, a_first_col, a_first_row, b_first_col, b_first_row, &
                      c_first_col, c_first_row)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      INTEGER, INTENT(IN)                                :: m, n, k
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                          :: beta
      TYPE(cp_fm_type), POINTER                          :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                      :: a_first_col, a_first_row, b_first_col, &
                                                            b_first_row, c_first_col, c_first_row

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_gemm'

      CHARACTER(LEN=1)                                   :: my_trans
      INTEGER                                            :: handle, handle1, my_multi
      INTEGER, PARAMETER                                 :: &
#if defined(__COSMA)
         my_multi_fallback = do_cosma
#else
      my_multi_fallback = do_scalapack
#endif
      INTEGER, DIMENSION(:), POINTER                     :: a_col_loc, a_row_loc, b_col_loc, &
                                                            b_row_loc, c_col_loc, c_row_loc
      TYPE(dbcsr_type)                                   :: a_db, b_db, c_db

      CALL timeset(routineN, handle)

      my_multi = cp_fm_get_mm_type()

      ! catch the special case that matrices have different blocking
      ! SCALAPACK/COSMA can deal with it but dbcsr doesn't like it
      CALL cp_fm_get_info(matrix_a, nrow_locals=a_row_loc, ncol_locals=a_col_loc)
      CALL cp_fm_get_info(matrix_b, nrow_locals=b_row_loc, ncol_locals=b_col_loc)
      CALL cp_fm_get_info(matrix_c, nrow_locals=c_row_loc, ncol_locals=c_col_loc)
      IF (PRESENT(a_first_row)) my_multi = my_multi_fallback
      IF (PRESENT(a_first_col)) my_multi = my_multi_fallback
      IF (PRESENT(b_first_row)) my_multi = my_multi_fallback
      IF (PRESENT(b_first_col)) my_multi = my_multi_fallback
      IF (PRESENT(c_first_row)) my_multi = my_multi_fallback
      IF (PRESENT(c_first_col)) my_multi = my_multi_fallback

      my_trans = transa; CALL uppercase(my_trans)
      IF (my_trans == 'T') THEN
         CALL cp_fm_get_info(matrix_a, nrow_locals=a_col_loc, ncol_locals=a_row_loc)
      END IF

      my_trans = transb; CALL uppercase(my_trans)
      IF (my_trans == 'T') THEN
         CALL cp_fm_get_info(matrix_b, nrow_locals=b_col_loc, ncol_locals=b_row_loc)
      END IF

      IF (my_multi .NE. do_scalapack .AND. my_multi .NE. do_cosma) THEN
         IF (SIZE(a_row_loc) == SIZE(c_row_loc)) THEN
            IF (ANY(a_row_loc - c_row_loc .NE. 0)) my_multi = my_multi_fallback
         ELSE
            my_multi = my_multi_fallback
         END IF
      END IF
      IF (my_multi .NE. do_scalapack .AND. my_multi .NE. do_cosma) THEN
         IF (SIZE(b_col_loc) == SIZE(c_col_loc)) THEN
            IF (ANY(b_col_loc - c_col_loc .NE. 0)) my_multi = my_multi_fallback
         ELSE
            my_multi = my_multi_fallback
         END IF
      END IF
      IF (my_multi .NE. do_scalapack .AND. my_multi .NE. do_cosma) THEN
         IF (SIZE(a_col_loc) == SIZE(b_row_loc)) THEN
            IF (ANY(a_col_loc - b_row_loc .NE. 0)) my_multi = my_multi_fallback
         ELSE
            my_multi = my_multi_fallback
         END IF
      END IF

      ! IMPORTANT do_scalapack is lowest value. If one processor has it set make all use it.
      IF (cp_fm_get_mm_type() .NE. do_scalapack .AND. &
          cp_fm_get_mm_type() .NE. do_cosma) CALL mp_min(my_multi, matrix_a%matrix_struct%para_env%group)

      SELECT CASE (my_multi)
      CASE (do_scalapack)
         CALL timeset("cp_gemm_fm_gemm", handle1)
         CALL cp_fm_gemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, matrix_c, &
                         a_first_col=a_first_col, &
                         a_first_row=a_first_row, &
                         b_first_col=b_first_col, &
                         b_first_row=b_first_row, &
                         c_first_col=c_first_col, &
                         c_first_row=c_first_row)
         CALL timestop(handle1)
      CASE (do_cosma)
#if defined(__COSMA)
         CALL timeset("cp_gemm_cosma", handle1)
         CALL offload_set_device()
         CALL cosma_pdgemm(transa=transa, transb=transb, m=m, n=n, k=k, alpha=alpha, &
                           matrix_a=matrix_a, matrix_b=matrix_b, beta=beta, matrix_c=matrix_c)
         CALL timestop(handle1)
#else
         CPABORT("CP2K compiled without the COSMA library.")
#endif
      CASE (do_dbcsr)
         CALL timeset("cp_gemm_dbcsr_mm", handle1)
         CALL copy_fm_to_dbcsr_bc(matrix_a, a_db)
         CALL copy_fm_to_dbcsr_bc(matrix_b, b_db)
         CALL copy_fm_to_dbcsr_bc(matrix_c, c_db)

         CALL dbcsr_multiply(transa, transb, alpha, a_db, b_db, beta, c_db, last_k=k)

         CALL copy_dbcsr_to_fm_bc(c_db, matrix_c)
         CALL dbcsr_release(a_db)
         CALL dbcsr_release(b_db)
         CALL dbcsr_release(c_db)
         CALL timestop(handle1)
      END SELECT
      CALL timestop(handle)

   END SUBROUTINE cp_gemm

#if defined(__COSMA)
! **************************************************************************************************
!> \brief Fortran wrapper for cosma_pdgemm.
!> \param transa ...
!> \param transb ...
!> \param m ...
!> \param n ...
!> \param k ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cosma_pdgemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, matrix_c)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      INTEGER, INTENT(IN)                                :: m, n, k
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                          :: beta
      TYPE(cp_fm_type), POINTER                          :: matrix_c

      INTERFACE
         SUBROUTINE cosma_pdgemm_c(transa, transb, m, n, k, alpha, a, ia, ja, desca, &
                                   b, ib, jb, descb, beta, c, ic, jc, descc) &
            BIND(C, name="cosma_pdgemm")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_CHAR
            CHARACTER(KIND=C_CHAR)                    :: transa
            CHARACTER(KIND=C_CHAR)                    :: transb
            INTEGER(KIND=C_INT)                       :: m
            INTEGER(KIND=C_INT)                       :: n
            INTEGER(KIND=C_INT)                       :: k
            REAL(KIND=C_DOUBLE)                       :: alpha
            TYPE(C_PTR), VALUE                        :: a
            INTEGER(KIND=C_INT)                       :: ia
            INTEGER(KIND=C_INT)                       :: ja
            TYPE(C_PTR), VALUE                        :: desca
            TYPE(C_PTR), VALUE                        :: b
            INTEGER(KIND=C_INT)                       :: ib
            INTEGER(KIND=C_INT)                       :: jb
            TYPE(C_PTR), VALUE                        :: descb
            REAL(KIND=C_DOUBLE)                       :: beta
            TYPE(C_PTR), VALUE                        :: c
            INTEGER(KIND=C_INT)                       :: ic
            INTEGER(KIND=C_INT)                       :: jc
            TYPE(C_PTR), VALUE                        :: descc
         END SUBROUTINE cosma_pdgemm_c
      END INTERFACE

      CALL cosma_pdgemm_c(transa=transa, transb=transb, m=m, n=n, k=k, &
                          alpha=alpha, &
                          a=C_LOC(matrix_a%local_data(1, 1)), ia=1, ja=1, &
                          desca=C_LOC(matrix_a%matrix_struct%descriptor(1)), &
                          b=C_LOC(matrix_b%local_data(1, 1)), ib=1, jb=1, &
                          descb=C_LOC(matrix_b%matrix_struct%descriptor(1)), &
                          beta=beta, &
                          c=C_LOC(matrix_c%local_data(1, 1)), ic=1, jc=1, &
                          descc=C_LOC(matrix_c%matrix_struct%descriptor(1)))

   END SUBROUTINE cosma_pdgemm
#endif

END MODULE cp_gemm_interface
