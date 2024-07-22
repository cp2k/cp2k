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
      spfft_transform_forward, spfft_multi_transform_forward_ptr, spfft_transform_backward, spfft_multi_transform_backward_ptr
   USE omp_lib, ONLY: omp_get_thread_num
#endif
   USE cp_files, ONLY: get_unit_number
   USE fft_kinds, ONLY: dp, int_8
   USE fft_plan, ONLY: fft_plan_type, FWFFT, BWFFT

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#if defined(__SPFFT) && defined(__parallel)
   INTEGER(C_INT), PARAMETER :: spfft_processing_unit = SPFFT_PU_HOST !+ SPFFT_PU_GPU
#endif

   PUBLIC :: spfft_create_plan_3d, spfft_3d
   PUBLIC :: spfft_create_plan_2d, spfft_2d
   PUBLIC :: spfft_create_plan_1d, spfft_1d
   PUBLIC :: spfft_destroy_plan

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
! **************************************************************************************************
   SUBROUTINE spfft_create_plan_2d(plan)

      TYPE(fft_plan_type), INTENT(INOUT)              :: plan
#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER :: number_of_threads, i, j
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: indices

      number_of_threads = 1
      !$OMP PARALLEL DEFAULT(NONE) SHARED(number_of_threads)
      !$OMP MASTER
!$    number_of_threads = omp_get_thread_num()
!$OMP END MASTER
!$OMP END PARALLEL

      ! Setup index map for SPFFT
      ALLOCATE (indices(3, plan%n_3d(1), plan%n_3d(2)))
      DO j = 1, plan%n_3d(2)
      DO i = 1, plan%n_3d(1)
         indices(1, i, j) = i - 1
         indices(2, i, j) = j - 1
         indices(3, i, j) = 0
      END DO
      END DO

      errorcode = spfft_transform_create_independent(plan%spfft_transform, number_of_threads, spfft_processing_unit, &
                                                     SPFFT_TRANS_C2C, plan%n_3d(1), plan%n_3d(2), 1, &
                                                     PRODUCT(plan%n_3d(1:2)), SPFFT_INDEX_TRIPLETS, indices)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_create_plan_2d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_create_plan_1d(plan)

      TYPE(fft_plan_type), INTENT(INOUT)              :: plan
#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER :: number_of_threads, i
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: indices

      number_of_threads = 1
      !$OMP PARALLEL DEFAULT(NONE) SHARED(number_of_threads)
      !$OMP MASTER
!$    number_of_threads = omp_get_thread_num()
!$OMP END MASTER
!$OMP END PARALLEL

      ! Setup index map for SPFFT
      ! We map 1D FFTs to 3D FFTs
      ALLOCATE (indices(3, plan%n))
      DO i = 1, plan%n
         indices(1, i) = i - 1
         indices(2, i) = 0
         indices(3, i) = 0
      END DO

      errorcode = spfft_transform_create_independent(plan%spfft_transform, number_of_threads, spfft_processing_unit, &
                                                     SPFFT_TRANS_C2C, plan%n, 1, 1, plan%n, SPFFT_INDEX_TRIPLETS, indices)
      IF (errorcode /= SPFFT_SUCCESS) CPABORT("Error in SPFFT!")
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_create_plan_1d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param scale ...
!> \param zin ...
!> \param zout ...
! **************************************************************************************************
   SUBROUTINE spfft_3d(plan, fsign, scale, zin, zout)

      TYPE(fft_plan_type), INTENT(IN)                      :: plan
      REAL(KIND=dp), INTENT(IN)                            :: scale
      INTEGER, INTENT(IN)                                  :: fsign
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(plan%n_3d)), INTENT(INOUT) :: zin, zout
#if defined (__SPFFT) && defined(__parallel)
      COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:), POINTER, CONTIGUOUS :: complex_ptr
      INTEGER(KIND=C_INT) :: errorcode
      TYPE(C_PTR) :: spfft_ptr

      ! set space domain array to use memory allocted by the library
      errorCode = spfft_transform_get_space_domain(plan%spfft_transform, SPFFT_PU_HOST, spfft_ptr)
      IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

      CALL C_F_POINTER(spfft_ptr, complex_ptr, [PRODUCT(plan%n_3d)])

      IF (plan%fft_in_place) THEN
         IF (fsign == FWFFT) THEN
            complex_ptr = zin
            errorCode = spfft_transform_forward(plan%spfft_transform, SPFFT_PU_HOST, zin, SPFFT_NO_SCALING)
         ELSE
            errorCode = spfft_transform_backward(plan%spfft_transform, zin, SPFFT_PU_HOST)
         END IF
      ELSE
         IF (fsign == FWFFT) THEN
            complex_ptr = zin
            errorCode = spfft_transform_forward(plan%spfft_transform, SPFFT_PU_HOST, zout, SPFFT_NO_SCALING)
         ELSE
            errorCode = spfft_transform_backward(plan%spfft_transform, zin, SPFFT_PU_HOST)
         END IF
      END IF
      IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

      IF (plan%fft_in_place) THEN
         IF (fsign == FWFFT) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zin, scale)
            zin = scale*CMPLX(zin, KIND=dp)
!$OMP END PARALLEL WORKSHARE
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zin, scale, complex_ptr)
            zin = scale*CMPLX(complex_ptr, KIND=dp)
!$OMP END PARALLEL WORKSHARE
         END IF
      ELSE
         IF (fsign == FWFFT) THEN
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zout, scale)
            zout = scale*CMPLX(zout, KIND=dp)
!$OMP END PARALLEL WORKSHARE
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zout, scale, complex_ptr)
            zout = scale*CMPLX(complex_ptr, KIND=dp)
!$OMP END PARALLEL WORKSHARE
         END IF
      END IF
#else
      MARK_USED(plan)
      MARK_USED(zin)
      MARK_USED(zout)
      MARK_USED(scale)
      MARK_USED(fsign)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_3d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param zin ...
!> \param zintermed ...
!> \param zout ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE spfft_2d(plan, fsign, zin, zintermed, zout, scale)

      TYPE(fft_plan_type), INTENT(IN)                      :: plan
      INTEGER, INTENT(IN)                                  :: fsign
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(plan%n_3d)), INTENT(INOUT), TARGET :: zin, zintermed, zout
      REAL(KIND=dp), INTENT(IN)                            :: scale
#if defined (__SPFFT) && defined(__parallel)
      INTEGER :: i, j
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER(KIND=C_INT), DIMENSION(:), ALLOCATABLE, TARGET :: scalingTypes
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: inputLocations, outputLocations

      ALLOCATE (scalingTypes(plan%n_3d(3)), inputLocations(plan%n_3d(3)), outputLocations(plan%n_3d(3)))
      scalingTypes = SPFFT_NO_SCALING

      IF (fsign == FWFFT) THEN
         IF (plan%trans) THEN
            ! Copy input to intermediate buffer
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(zin,scale,plan,zintermed) PRIVATE(i, j) COLLAPSE(2)
            DO i = 1, plan%n_3d(1)*plan%n_3d(2)
            DO j = 1, plan%n_3d(3)
               zintermed(i + (j - 1)*plan%n_3d(1)) = scale*zin(j + (i - 1)*plan%n_3d(3))
            END DO
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zin,scale,zintermed)
            zintermed = scale*zin
!$OMP END PARALLEL WORKSHARE
         END IF

         ! Compile an array with the starting points of each FFT
         DO i = 1, plan%n_3d(3)
            inputLocations(i) = C_LOC(zintermed(i*plan%n_3d(1)*plan%n_3d(2)))
            outputLocations(i) = C_LOC(zout(i*plan%n_3d(1)*plan%n_3d(2)))
         END DO

         ! Perform the actual FFTs
         errorCode = spfft_multi_transform_forward_ptr(plan%n_3d(3), plan%spfft_transform, &
                                                       inputLocations, outputLocations, C_LOC(scalingTypes))
         IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")
      ELSE
         DO i = 1, plan%n_3d(3)
            inputLocations(i) = C_LOC(zin(i*plan%n_3d(1)*plan%n_3d(2)))
            outputLocations(i) = C_LOC(zintermed(i*plan%n_3d(1)*plan%n_3d(2)))
         END DO

         errorCode = spfft_multi_transform_backward_ptr(plan%n_3d(3), plan%spfft_transform, &
                                                        inputLocations, outputLocations)
         IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

         IF (plan%trans) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(zout,scale,plan,zintermed) PRIVATE(i, j) COLLAPSE(2)
            DO i = 1, plan%n_3d(1)*plan%n_3d(2)
            DO j = 1, plan%n_3d(3)
               zout(j + (i - 1)*plan%n_3d(3)) = scale*zintermed(i + (j - 1)*plan%n_3d(3))
            END DO
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zout,scale,zintermed)
            zout = scale*zintermed
!$OMP END PARALLEL WORKSHARE
         END IF
      END IF

      DEALLOCATE (scalingTypes, inputLocations, outputLocations)
#else
      MARK_USED(plan)
      MARK_USED(zin)
      MARK_USED(zout)
      MARK_USED(scale)
      MARK_USED(fsign)
      MARK_USED(zintermed)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_2d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
!> \param fsign ...
!> \param zin ...
!> \param zout ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE spfft_1d(plan, fsign, zin, zout, scale)

      TYPE(fft_plan_type), INTENT(IN)                      :: plan
      INTEGER, INTENT(IN)                                  :: fsign
      COMPLEX(KIND=dp), DIMENSION(plan%n*plan%m), INTENT(INOUT), TARGET :: zin, zout
      REAL(KIND=dp), INTENT(IN)                            :: scale
#if defined (__SPFFT) && defined(__parallel)
      COMPLEX(KIND=C_DOUBLE_COMPLEX), DIMENSION(:, :), ALLOCATABLE, TARGET :: complex_array
      INTEGER :: i, j
      INTEGER(KIND=C_INT) :: errorcode
      INTEGER(KIND=C_INT), DIMENSION(:), ALLOCATABLE, TARGET :: scalingTypes
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: inputLocations, outputLocations, transforms

      ! Allocate a temporary array to ensure enough memory being available
      ALLOCATE (complex_array(plan%n, plan%m))

      ALLOCATE (scalingTypes(plan%m), inputLocations(plan%m), outputLocations(plan%m), transforms(plan%m))
      scalingTypes = SPFFT_NO_SCALING
      transforms = plan%spfft_transform

      IF (fsign == FWFFT) THEN
         ! Scale here to prevent additional copies between the buffer and the output array
         IF (plan%trans) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(zin,complex_array,scale,plan) PRIVATE(i,j) COLLAPSE(2)
            DO i = 1, plan%n
            DO j = 1, plan%m
               complex_array(i, j) = scale*zin(j + (i - 1)*plan%m)
            END DO
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zin,scale,complex_array,plan)
            complex_array(:, :) = scale*RESHAPE(zin, [plan%n, plan%m])
!$OMP END PARALLEL WORKSHARE
         END IF

         DO i = 1, plan%m
            inputLocations(i) = C_LOC(complex_array(:, i))
            outputLocations(i) = C_LOC(zout(i*plan%n))
         END DO

         errorCode = spfft_multi_transform_forward_ptr(plan%m, transforms, &
                                                       inputLocations, outputLocations, C_LOC(scalingTypes))
         IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

      ELSE
         ! Just setup pointers
         DO i = 1, plan%m
            inputLocations(i) = C_LOC(zin(i*plan%n))
            outputLocations(i) = C_LOC(complex_array(:, i))
         END DO
         errorCode = spfft_multi_transform_backward_ptr(plan%m, transforms, &
                                                        inputLocations, outputLocations)
         IF (errorCode /= SPFFT_SUCCESS) CPABORT("FFT with SpFFT failed!")

         IF (plan%trans) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(zout,complex_array,scale,plan) PRIVATE(i,j) COLLAPSE(2)
            DO i = 1, plan%n
            DO j = 1, plan%m
               zout(j - 1 + (i - 1)*plan%m) = scale*complex_array(i, j)
            END DO
            END DO
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(zout,scale,plan,complex_array)
            zout = scale*RESHAPE(complex_array, [plan%m*plan%n])
!$OMP END PARALLEL WORKSHARE
         END IF
      END IF

      DEALLOCATE (complex_array, inputLocations, outputLocations, transforms)
#else
      MARK_USED(plan)
      MARK_USED(zin)
      MARK_USED(zout)
      MARK_USED(scale)
      MARK_USED(fsign)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_1d

! **************************************************************************************************
!> \brief ...
!> \param plan ...
! **************************************************************************************************
   SUBROUTINE spfft_destroy_plan(plan)

      TYPE(fft_plan_type), INTENT(INOUT)   :: plan

#if defined ( __SPFFT ) && defined(__parallel)
      INTEGER(KIND=C_INT) :: errorCode

      errorCode = spfft_transform_destroy(plan%spfft_transform)
#else
      MARK_USED(plan)
      CPABORT("In order to use SPFFT, you have to install it!")
#endif

   END SUBROUTINE spfft_destroy_plan

END MODULE
