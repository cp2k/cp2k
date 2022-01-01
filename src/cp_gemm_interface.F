!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
   USE cp_fm_types,                     ONLY: cp_fm_get_mm_type,&
                                              cp_fm_type
   USE input_constants,                 ONLY: do_cosma,&
                                              do_scalapack
   USE kinds,                           ONLY: dp
   USE offload_api,                     ONLY: offload_set_device
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

      INTEGER                                            :: handle, handle1, my_multi

      CALL timeset(routineN, handle)

      my_multi = cp_fm_get_mm_type()

      SELECT CASE (my_multi)
      CASE (do_scalapack)
         CALL timeset(routineN//"_fm_gemm", handle1)
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
         CALL timeset(routineN//"_cosma", handle1)
         CALL offload_set_device()
         CALL cosma_pdgemm(transa=transa, transb=transb, m=m, n=n, k=k, alpha=alpha, &
                           matrix_a=matrix_a, matrix_b=matrix_b, beta=beta, matrix_c=matrix_c, &
                           a_first_col=a_first_col, &
                           a_first_row=a_first_row, &
                           b_first_col=b_first_col, &
                           b_first_row=b_first_row, &
                           c_first_col=c_first_col, &
                           c_first_row=c_first_row)
         CALL timestop(handle1)
#else
         CPABORT("CP2K compiled without the COSMA library.")
#endif
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
!> \param a_first_col ...
!> \param a_first_row ...
!> \param b_first_col ...
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cosma_pdgemm(transa, transb, m, n, k, alpha, matrix_a, matrix_b, beta, matrix_c, &
                           a_first_col, a_first_row, b_first_col, b_first_row, &
                           c_first_col, c_first_row)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      INTEGER, INTENT(IN)                                :: m, n, k
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: matrix_a, matrix_b
      REAL(KIND=dp), INTENT(IN)                          :: beta
      TYPE(cp_fm_type), POINTER                          :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                      :: a_first_col, a_first_row, b_first_col, &
                                                            b_first_row, c_first_col, c_first_row

      INTEGER                                            :: i_a, i_b, i_c, j_a, j_b, j_c
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

      IF (PRESENT(a_first_row)) THEN
         i_a = a_first_row
      ELSE
         i_a = 1
      END IF
      IF (PRESENT(a_first_col)) THEN
         j_a = a_first_col
      ELSE
         j_a = 1
      END IF
      IF (PRESENT(b_first_row)) THEN
         i_b = b_first_row
      ELSE
         i_b = 1
      END IF
      IF (PRESENT(b_first_col)) THEN
         j_b = b_first_col
      ELSE
         j_b = 1
      END IF
      IF (PRESENT(c_first_row)) THEN
         i_c = c_first_row
      ELSE
         i_c = 1
      END IF
      IF (PRESENT(c_first_col)) THEN
         j_c = c_first_col
      ELSE
         j_c = 1
      END IF

      CALL cosma_pdgemm_c(transa=transa, transb=transb, m=m, n=n, k=k, &
                          alpha=alpha, &
                          a=C_LOC(matrix_a%local_data(1, 1)), ia=i_a, ja=j_a, &
                          desca=C_LOC(matrix_a%matrix_struct%descriptor(1)), &
                          b=C_LOC(matrix_b%local_data(1, 1)), ib=i_b, jb=j_b, &
                          descb=C_LOC(matrix_b%matrix_struct%descriptor(1)), &
                          beta=beta, &
                          c=C_LOC(matrix_c%local_data(1, 1)), ic=i_c, jc=j_c, &
                          descc=C_LOC(matrix_c%matrix_struct%descriptor(1)))

   END SUBROUTINE cosma_pdgemm
#endif

END MODULE cp_gemm_interface
