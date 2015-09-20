!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic linear algebra operations for full matrixes
!> \par History
!>      08.2002 splitted out of qs_blacs [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_gemm_interface
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_multiply,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm_bc,&
                                             copy_fm_to_dbcsr_bc
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_mm_type,&
                                             cp_fm_type
  USE input_constants,                 ONLY: do_dbcsr,&
                                             do_pdgemm
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_min
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_gemm_interface'

  PUBLIC :: cp_gemm

CONTAINS

! *****************************************************************************
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
!> \param b_first_col ...
!> \param a_first_row ...
!> \param b_first_row ...
!> \param c_first_col ...
!> \param c_first_row ...
! *****************************************************************************
  SUBROUTINE cp_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
       matrix_c,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    INTEGER, INTENT(IN)                      :: m, n, k
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_fm_type), POINTER                :: matrix_a, matrix_b
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_fm_type), POINTER                :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: b_first_col, a_first_row, &
                                                b_first_row, c_first_col, &
                                                c_first_row

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_gemm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans
    INTEGER                                  :: handle, handle1, my_multi
    INTEGER, DIMENSION(:), POINTER           :: a_col_loc, a_row_loc, &
                                                b_col_loc, b_row_loc, &
                                                c_col_loc, c_row_loc
    TYPE(cp_dbcsr_type)                      :: a_db, b_db, c_db

    CALL timeset(routineN,handle)

    my_multi=cp_fm_get_mm_type()
    ! catch the special case that matrices have different blocking
    ! SCALAPACK can deal with it but dbcsr doesn't like it
    CALL cp_fm_get_info(matrix_a,nrow_locals=a_row_loc,ncol_locals=a_col_loc)
    CALL cp_fm_get_info(matrix_b,nrow_locals=b_row_loc,ncol_locals=b_col_loc)
    CALL cp_fm_get_info(matrix_c,nrow_locals=c_row_loc,ncol_locals=c_col_loc)
    IF(PRESENT(b_first_col))my_multi=do_pdgemm
    IF(PRESENT(a_first_row))my_multi=do_pdgemm
    IF(PRESENT(b_first_row))my_multi=do_pdgemm
    IF(PRESENT(c_first_col))my_multi=do_pdgemm
    IF(PRESENT(c_first_row))my_multi=do_pdgemm

    my_trans=transa; CALL uppercase(my_trans)
    IF(my_trans=='T')THEN
       CALL cp_fm_get_info(matrix_a,nrow_locals=a_col_loc,ncol_locals=a_row_loc)
    END IF 

    my_trans=transb; CALL uppercase(my_trans)
    IF(my_trans=='T')THEN
       CALL cp_fm_get_info(matrix_b,nrow_locals=b_col_loc,ncol_locals=b_row_loc)
    END IF

    IF(my_multi.NE.do_pdgemm)THEN
       IF(SIZE(a_row_loc)==SIZE(c_row_loc))THEN
          IF(ANY(a_row_loc-c_row_loc.NE.0)) my_multi=do_pdgemm
       ELSE
          my_multi=do_pdgemm
       END IF
    END IF
    IF(my_multi.NE.do_pdgemm)THEN
       IF(SIZE(b_col_loc)==SIZE(c_col_loc))THEN
          IF(ANY(b_col_loc-c_col_loc.NE.0)) my_multi=do_pdgemm
       ELSE
          my_multi=do_pdgemm
       END IF
    END IF
    IF(my_multi.NE.do_pdgemm)THEN
       IF(SIZE(a_col_loc)==SIZE(b_row_loc))THEN
          IF(ANY(a_col_loc-b_row_loc.NE.0)) my_multi=do_pdgemm
       ELSE
          my_multi=do_pdgemm
       END IF
    END IF

    ! IMPORTANT do_pdgemm is lowest value. If one processor has it set make all do pdgemm
    IF(cp_fm_get_mm_type().NE.do_pdgemm) CALL mp_min(my_multi,matrix_a%matrix_struct%para_env%group)  
    
    SELECT CASE(my_multi)
    CASE (do_pdgemm)
        CALL timeset("cp_gemm_fm_gemm",handle1)
       CALL cp_fm_gemm(transa,transb,m,n,k,alpha,matrix_a,matrix_b,beta,&
        matrix_c,b_first_col,a_first_row,b_first_row,c_first_col,c_first_row)
       CALL timestop(handle1)
    CASE (do_dbcsr)
        CALL timeset("cp_gemm_dbcsr_mm",handle1)
       CALL copy_fm_to_dbcsr_bc(matrix_a,a_db)
       CALL copy_fm_to_dbcsr_bc(matrix_b,b_db)
       CALL copy_fm_to_dbcsr_bc(matrix_c,c_db)
   
       CALL cp_dbcsr_multiply(transa,transb,alpha, a_db, b_db, beta, c_db, last_k=k)
   
       CALL copy_dbcsr_to_fm_bc(c_db, matrix_c)
       CALL cp_dbcsr_release (a_db)
       CALL cp_dbcsr_release (b_db)
       CALL cp_dbcsr_release (c_db)
       CALL timestop(handle1)
    END SELECT
    CALL timestop(handle)

  END SUBROUTINE cp_gemm

END MODULE cp_gemm_interface
