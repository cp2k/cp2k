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
      spfft_transform_create_independent, SPFFT_PU_HOST, SPFFT_PU_GPU, SPFFT_SUCCESS, &
      SPFFT_TRANS_C2C, SPFFT_INDEX_TRIPLETS, spfft_transform_destroy, SPFFT_NO_SCALING, spfft_transform_get_space_domain, &
      spfft_transform_forward
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
   INTEGER(C_INT), PARAMETER :: spfft_processing_unit = SPFFT_PU_HOST !+ SPFFT_PU_GPU
#endif

   PUBLIC :: spfft_create_plan_3d, spfft_destroy_plan, spfft_3d

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_create_plan_3d(plan)

      TYPE(fft_plan_type), INTENT(INOUT)              :: plan
#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER :: number_of_threads, i, j, k
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: indices

      number_of_threads = 1
      !$OMP PARALLEL DEFAULT(NONE) SHARED(number_of_threads)
      !$OMP MASTER
!$    number_of_threads = omp_get_thread_num()
!$OMP END MASTER
!$OMP END PARALLEL

      ! Setup index map for SPFFT
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

      errorcode = spfft_transform_create_independent(plan%spfft_transform, number_of_threads, spfft_processing_unit, &
                                                     SPFFT_TRANS_C2C, plan%n_3d(1), plan%n_3d(2), plan%n_3d(3), &
                                                     PRODUCT(plan%n_3d), SPFFT_INDEX_TRIPLETS, indices)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_create_plan_3d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param scale ...
!> \param zin ...
!> \param zout ...
!> \param stat ...
! **************************************************************************************************
   SUBROUTINE spfft_3d(plan, scale, zin, zout, stat)

      TYPE(fft_plan_type), INTENT(IN)                      :: plan
      REAL(KIND=dp), INTENT(IN)                            :: scale
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), CONTIGUOUS, TARGET :: zin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), CONTIGUOUS, TARGET :: zout
      INTEGER, INTENT(OUT)                                 :: stat
#if defined (__SPFFT) && defined(__parallel)
      COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:, :, :), POINTER, CONTIGUOUS :: complex_ptr
      INTEGER(KIND=C_INT) :: errorcode
      TYPE(C_PTR) :: spfft_ptr

      ! set space domain array to use memory allocted by the library
      errorCode = spfft_transform_get_space_domain(plan%spfft_transform, SPFFT_PU_HOST, spfft_ptr)
      IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

      CALL C_F_POINTER(spfft_ptr, complex_ptr, plan%n_3d)

      errorCode = spfft_transform_forward(plan%spfft_transform, SPFFT_PU_HOST, zin, SPFFT_NO_SCALING)
      IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zout, scale, complex_ptr)
      zout = scale*CMPLX(complex_ptr, KIND=dp)
!$OMP END PARALLEL WORKSHARE
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif
      stat = 0

   END SUBROUTINE spfft_3d

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
