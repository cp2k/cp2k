!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE spfft_lib

   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_CHAR, &
                            C_DOUBLE, &
                            C_DOUBLE_COMPLEX, &
                            C_INT, &
                            C_PTR
#if defined(__SPFFT) && defined(__parallel)
   USE ISO_C_BINDING, ONLY: &
      C_FLOAT, &
      C_FLOAT_COMPLEX, &
      C_FUNPTR, &
      C_INT32_T, &
      C_INTPTR_T, &
      C_LOC, &
      C_NULL_CHAR, &
      C_SIZE_T, C_F_POINTER
   USE spfft, ONLY: &
      spfft_grid_create, spfft_grid_destroy, spfft_transform_create, SPFFT_PU_HOST, SPFFT_PU_GPU, SPFFT_SUCCESS, &
      SPFFT_TRANS_C2C, SPFFT_INDEX_TRIPLETS, spfft_transform_destroy
   USE omp_lib, ONLY: omp_get_thread_num
#endif
   USE cp_files, ONLY: get_unit_number
   USE fft_kinds, ONLY: dp, int_8
   USE fft_plan, ONLY: fft_plan_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#if defined(__SPFFT) && defined(__parallel)
   INTEGER(C_INT), PARAMETER :: spfft_processing_unit = SPFFT_PU_HOST + SPFFT_PU_GPU
#endif

   PUBLIC :: spfft_create_plan_3d, spfft_create_plan_1dm, spfft_destroy_plan

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_create_plan_3d(plan)

      TYPE(fft_plan_type), INTENT(INOUT)              :: plan
#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER :: number_of_threads, i, j, k
      INTEGER(KIND=int_8) :: counter
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: indices
      TYPE(C_PTR) :: grid

      number_of_threads = 1
!$    number_of_threads = omp_get_thread_num()
      IF (number_of_threads > 1) &
         CPABORT("SPFFT can currently not be called with more than one thread!")

      errorcode = spfft_grid_create(grid, plan%n_3d(1), plan%n_3d(2), plan%n_3d(3), plan%n_3d(2)*plan%n_3d(3), &
                                    spfft_processing_unit, number_of_threads)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")

      ! Setup index map for SPFFT
      counter = 1
      ALLOCATE (indices(3, plan%n_3d(1), plan%n_3d(2), plan%n_3d(3)))
      DO k = 1, plan%n_3d(3)
      DO j = 1, plan%n_3d(2)
      DO i = 1, plan%n_3d(1)
         indices(1, i, j, k) = i - 1
         indices(2, i, j, k) = j - 1
         indices(3, i, j, k) = k - 1
      END DO
      END DO
      END DO

      ! TODO:
      ! Use ngpts instead of n (before SPFFT_INDEX_TRIPLETS)
      errorcode = spfft_transform_create(plan%spfft_transform, grid, spfft_processing_unit, SPFFT_TRANS_C2C, &
                                         plan%n_3d(1), plan%n_3d(2), plan%n_3d(3), plan%n_3d(3), &
                                         plan%n, SPFFT_INDEX_TRIPLETS, indices)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")

      errorcode = spfft_grid_destroy(grid)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_create_plan_3d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_create_plan_1dm(plan)

      TYPE(fft_plan_type), INTENT(INOUT)              :: plan
#if defined (__SPFFT) && defined(__parallel)
      MARK_USED(plan)
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_create_plan_1dm

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_destroy_plan(plan)

      TYPE(fft_plan_type), INTENT(INOUT)   :: plan

#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER(KIND=C_INT) :: errorCode

      errorCode = spfft_transform_destroy(plan%spfft_transform)
      IF (errorCode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_destroy_plan

END MODULE
