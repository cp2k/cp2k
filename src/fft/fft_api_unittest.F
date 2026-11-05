!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!
PROGRAM fft_api_unittest
   USE fft_api,                         ONLY: &
        FFT_LIBRARY_BACKEND_DEFAULT, FFT_LIBRARY_FFTW_MODE_ESTIMATE, fft_1d_local, &
        fft_1d_local_inplace, fft_2d_local, fft_2d_local_inplace, fft_3d_local, &
        fft_3d_local_inplace, fft_library_finalize, fft_library_init, fft_register_1d_local, &
        fft_register_1d_local_inplace, fft_register_2d_local, fft_register_2d_local_inplace, &
        fft_register_3d_local, fft_register_3d_local_inplace
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE offload_api,                     ONLY: offload_get_device_count,&
                                              offload_set_chosen_device

   IMPLICIT NONE

   TYPE(mp_comm_type)                                  :: mp_comm
   INTEGER                                             :: errors, my_rank, status

   CALL mp_world_init(mp_comm)
   my_rank = mp_comm%mepos
   IF (offload_get_device_count() > 0) THEN
      CALL offload_set_chosen_device(MOD(my_rank, offload_get_device_count()))
   END IF

   CALL fft_library_init(FFT_LIBRARY_BACKEND_DEFAULT, FFT_LIBRARY_FFTW_MODE_ESTIMATE, &
                         .TRUE., "", status)
   IF (status /= 0) ERROR STOP "FFT library initialization failed"

   errors = 0
   errors = errors + test_1d_c2c(15, 13, .TRUE., .TRUE.)
   errors = errors + test_1d_c2c(9, 11, .TRUE., .FALSE.)
   errors = errors + test_1d_c2c(10, 14, .FALSE., .TRUE.)
   errors = errors + test_1d_c2c(14, 13, .FALSE., .FALSE.)

   errors = errors + test_1d_c2c_inplace(15, 13, .TRUE.)
   errors = errors + test_1d_c2c_inplace(14, 13, .FALSE.)

   errors = errors + test_2d_c2c([10, 10], 10, .TRUE., .TRUE.)
   errors = errors + test_2d_c2c([8, 9], 10, .TRUE., .FALSE.)
   errors = errors + test_2d_c2c([7, 10], 11, .FALSE., .TRUE.)
   errors = errors + test_2d_c2c([6, 7], 10, .FALSE., .FALSE.)

   errors = errors + test_2d_c2c_inplace([8, 9], 10, .TRUE.)
   errors = errors + test_2d_c2c_inplace([7, 10], 11, .FALSE.)

   errors = errors + test_3d_c2c([8, 8, 8], .FALSE., .FALSE.)
   errors = errors + test_3d_c2c([7, 5, 3], .TRUE., .TRUE.)

   errors = errors + test_3d_c2c_inplace([8, 8, 8], .FALSE.)
   errors = errors + test_3d_c2c_inplace([7, 5, 3], .TRUE.)

   CALL fft_library_finalize("")
   CALL mp_world_finalize()

   IF (my_rank == 0) THEN
      IF (errors == 0) THEN
         WRITE (*, '(A)') "All Fortran FFT-API tests have passed :-)"
      ELSE
         WRITE (*, '(A,I0)') "Fortran FFT-API tests failed: ", errors
      END IF
   END IF
   IF (errors /= 0) ERROR STOP

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_1d_c2c(fft_size, number_of_ffts, transpose_in, transpose_out)
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out

      COMPLEX(KIND=dp), ALLOCATABLE                      :: input(:), output(:), reference(:)
      INTEGER                                            :: i, ld_in, ld_out
      REAL(KIND=dp)                                      :: max_error

      ld_in = MERGE(number_of_ffts, fft_size, transpose_in)
      ld_out = MERGE(number_of_ffts, fft_size, transpose_out)
      ALLOCATE (input(fft_size*number_of_ffts), output(fft_size*number_of_ffts), &
                reference(fft_size*number_of_ffts))
      DO i = 1, SIZE(input)
         input(i) = CMPLX(SIN(0.1_dp*i), COS(0.07_dp*i), KIND=dp)
      END DO
      reference = input

      CALL fft_register_1d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, &
                                 ld_in, ld_out, input, output)
      CALL fft_1d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, &
                        ld_in, ld_out, input, output)
      CALL fft_register_1d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, &
                                 ld_out, ld_in, output, input)
      CALL fft_1d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, &
                        ld_out, ld_in, output, input)

      max_error = MAXVAL(ABS(input - REAL(fft_size, KIND=dp)*reference))
      test_1d_c2c = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,2I4,2L2)') &
            "1D C2C roundtrip error ", max_error, " size/batch ", fft_size, number_of_ffts, &
            transpose_in, transpose_out
         test_1d_c2c = 1
      END IF
      DEALLOCATE (input, output, reference)
   END FUNCTION test_1d_c2c

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_1d_c2c_inplace(fft_size, number_of_ffts, transpose)
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose

      COMPLEX(KIND=dp), ALLOCATABLE                      :: grid(:), reference(:)
      INTEGER                                            :: i, ld
      REAL(KIND=dp)                                      :: max_error

      ld = MERGE(number_of_ffts, fft_size, transpose)
      ALLOCATE (grid(fft_size*number_of_ffts), reference(fft_size*number_of_ffts))
      DO i = 1, SIZE(grid)
         grid(i) = CMPLX(SIN(0.11_dp*i), COS(0.05_dp*i), KIND=dp)
      END DO
      reference = grid
      CALL fft_register_1d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, ld, grid)
      CALL fft_1d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, ld, grid)
      CALL fft_register_1d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, ld, grid)
      CALL fft_1d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, ld, grid)

      max_error = MAXVAL(ABS(grid - REAL(fft_size, KIND=dp)*reference))
      test_1d_c2c_inplace = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,2I4,L2)') &
            "1D C2C inplace roundtrip error ", max_error, " size/batch ", fft_size, &
            number_of_ffts, transpose
         test_1d_c2c_inplace = 1
      END IF
      DEALLOCATE (grid, reference)
   END FUNCTION test_1d_c2c_inplace

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_2d_c2c(fft_size, number_of_ffts, transpose_in, transpose_out)
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out

      COMPLEX(KIND=dp), ALLOCATABLE                      :: input(:), output(:), reference(:)
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: max_error

      ALLOCATE (input(PRODUCT(fft_size)*number_of_ffts), &
                output(PRODUCT(fft_size)*number_of_ffts), &
                reference(PRODUCT(fft_size)*number_of_ffts))
      DO i = 1, SIZE(input)
         input(i) = CMPLX(SIN(0.09_dp*i), COS(0.04_dp*i), KIND=dp)
      END DO
      reference = input

      CALL fft_register_2d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, &
                                 input, output)
      CALL fft_2d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, input, output)
      CALL fft_register_2d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, &
                                 output, input)
      CALL fft_2d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, output, input)

      max_error = MAXVAL(ABS(input - REAL(PRODUCT(fft_size), KIND=dp)*reference))
      test_2d_c2c = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,3I4,2L2)') &
            "2D C2C roundtrip error ", max_error, " size/batch ", fft_size, number_of_ffts, &
            transpose_in, transpose_out
         test_2d_c2c = 1
      END IF
      DEALLOCATE (input, output, reference)
   END FUNCTION test_2d_c2c

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_2d_c2c_inplace(fft_size, number_of_ffts, transpose)
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose

      COMPLEX(KIND=dp), ALLOCATABLE                      :: grid(:), reference(:)
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: max_error

      ALLOCATE (grid(PRODUCT(fft_size)*number_of_ffts), &
                reference(PRODUCT(fft_size)*number_of_ffts))
      DO i = 1, SIZE(grid)
         grid(i) = CMPLX(SIN(0.08_dp*i), COS(0.03_dp*i), KIND=dp)
      END DO
      reference = grid
      CALL fft_register_2d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_2d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_register_2d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_2d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, grid)

      max_error = MAXVAL(ABS(grid - REAL(PRODUCT(fft_size), KIND=dp)*reference))
      test_2d_c2c_inplace = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,3I4,L2)') &
            "2D C2C inplace roundtrip error ", max_error, " size/batch ", fft_size, &
            number_of_ffts, transpose
         test_2d_c2c_inplace = 1
      END IF
      DEALLOCATE (grid, reference)
   END FUNCTION test_2d_c2c_inplace

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_3d_c2c(fft_size, transpose_in, transpose_out)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out

      COMPLEX(KIND=dp), ALLOCATABLE                      :: input(:), output(:), reference(:)
      INTEGER                                            :: i, number_of_ffts
      REAL(KIND=dp)                                      :: max_error

      number_of_ffts = PRODUCT(fft_size)
      ALLOCATE (input(PRODUCT(fft_size)*number_of_ffts), &
                output(PRODUCT(fft_size)*number_of_ffts), &
                reference(PRODUCT(fft_size)*number_of_ffts))
      DO i = 1, SIZE(input)
         input(i) = CMPLX(SIN(0.07_dp*i), COS(0.02_dp*i), KIND=dp)
      END DO
      reference = input

      CALL fft_register_3d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, &
                                 input, output)
      CALL fft_3d_local(.TRUE., fft_size, number_of_ffts, transpose_in, transpose_out, input, output)
      CALL fft_register_3d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, &
                                 output, input)
      CALL fft_3d_local(.FALSE., fft_size, number_of_ffts, transpose_out, transpose_in, output, input)

      max_error = MAXVAL(ABS(input - REAL(PRODUCT(fft_size), KIND=dp)*reference))
      test_3d_c2c = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,4I4,2L2)') &
            "3D C2C roundtrip error ", max_error, " size/batch ", fft_size, number_of_ffts, &
            transpose_in, transpose_out
         test_3d_c2c = 1
      END IF
      DEALLOCATE (input, output, reference)
   END FUNCTION test_3d_c2c

! **************************************************************************************************
!> \brief ...
!> \param fft_size ...
!> \param transpose ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION test_3d_c2c_inplace(fft_size, transpose)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      LOGICAL, INTENT(IN)                                :: transpose

      COMPLEX(KIND=dp), ALLOCATABLE                      :: grid(:), reference(:)
      INTEGER                                            :: i, number_of_ffts
      REAL(KIND=dp)                                      :: max_error

      number_of_ffts = PRODUCT(fft_size)
      ALLOCATE (grid(PRODUCT(fft_size)*number_of_ffts), &
                reference(PRODUCT(fft_size)*number_of_ffts))
      DO i = 1, SIZE(grid)
         grid(i) = CMPLX(SIN(0.06_dp*i), COS(0.015_dp*i), KIND=dp)
      END DO
      reference = grid
      CALL fft_register_3d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_3d_local_inplace(.TRUE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_register_3d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, grid)
      CALL fft_3d_local_inplace(.FALSE., fft_size, number_of_ffts, transpose, grid)

      max_error = MAXVAL(ABS(grid - REAL(PRODUCT(fft_size), KIND=dp)*reference))
      test_3d_c2c_inplace = 0
      IF (max_error > 1.0E-8_dp) THEN
         IF (my_rank == 0) WRITE (*, '(A,F12.3,A,4I4,L2)') &
            "3D C2C inplace roundtrip error ", max_error, " size/batch ", fft_size, &
            number_of_ffts, transpose
         test_3d_c2c_inplace = 1
      END IF
      DEALLOCATE (grid, reference)
   END FUNCTION test_3d_c2c_inplace

END PROGRAM fft_api_unittest
