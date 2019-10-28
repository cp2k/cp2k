!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate general matrix functions using (sca)lapack diagonalization
! **************************************************************************************************

MODULE cp_dbcsr_func
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add_on_diag,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_multiply,&
                                              dbcsr_release,&
                                              dbcsr_set_diag,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_func'

   PUBLIC :: cp_dbcsr_matrix_function

CONTAINS

! **************************************************************************************************
!> \brief Calculate general matrix functions using (sca)lapack diagonalization
!> \param A matrix
!> \param f function f
!> \param f_A f(A)
!> \param para_env ...
!> \param blacs_env ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_matrix_function(A, f, f_A, para_env, blacs_env)

      TYPE(dbcsr_type), INTENT(IN)                       :: A
      INTERFACE
         FUNCTION f(values)
            IMPORT :: dp
            REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: values
            REAL(KIND=dp), DIMENSION(SIZE(values)) :: f
         END FUNCTION
      END INTERFACE
      TYPE(dbcsr_type), INTENT(INOUT)                    :: f_A
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_matrix_function', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigvals
      TYPE(dbcsr_type)                                   :: diag, eigvec, tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_create(eigvec, template=A, matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_get_info(A, nfullrows_total=n)
      ALLOCATE (eigvals(n))
      CALL cp_dbcsr_syevd(A, eigvec, eigvals, para_env, blacs_env)

      CALL dbcsr_create(diag, template=A)
      CALL dbcsr_add_on_diag(diag, 1.0_dp)
      eigvals(:) = f(eigvals)
      CALL dbcsr_set_diag(diag, eigvals)
      DEALLOCATE (eigvals)

      CALL dbcsr_create(tmp, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_multiply("N", "N", 1.0_dp, eigvec, diag, 0.0_dp, tmp)
      CALL dbcsr_release(diag)
      CALL dbcsr_multiply("N", "T", 1.0_dp, tmp, eigvec, 0.0_dp, f_A)
      CALL dbcsr_release(eigvec)
      CALL dbcsr_release(tmp)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE
