!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fortran API for the FFT package, which is written in C.
!> \author Frederick Stein
! **************************************************************************************************
MODULE fft_api
   USE ISO_C_BINDING, ONLY: &
      C_ASSOCIATED, C_BOOL, C_CHAR, C_DOUBLE, C_DOUBLE_COMPLEX, C_F_POINTER, C_INT, C_LOC, C_NULL_CHAR, &
      C_NULL_PTR, C_PTR
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_api'

   PUBLIC :: fft_library_init, fft_library_finalize
   PUBLIC :: fft_print_timing_report
   PUBLIC :: fft_backend_in_use, fft_has_guru_interface
   PUBLIC :: fft_allocate, fft_deallocate
   PUBLIC :: fft_register_1d_local, fft_register_2d_local, fft_register_3d_local
   PUBLIC :: fft_register_1d_local_inplace, fft_register_2d_local_inplace, fft_register_3d_local_inplace
   PUBLIC :: fft_1d_local, fft_2d_local, fft_3d_local
   PUBLIC :: fft_1d_local_inplace, fft_2d_local_inplace, fft_3d_local_inplace

   ! Keep it synchronized with fft_api_F.c
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_BACKEND_DEFAULT = 1
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_BACKEND_FFTW = 2
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_BACKEND_GPU = 3
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_FFTW_MODE_ESTIMATE = 11
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_FFTW_MODE_MEASURE = 12
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_FFTW_MODE_PATIENT = 13
   INTEGER, PARAMETER, PUBLIC :: FFT_LIBRARY_FFTW_MODE_EXHAUSTIVE = 14

   INTERFACE fft_allocate
      MODULE PROCEDURE :: fft_allocate_real_1
      MODULE PROCEDURE :: fft_allocate_real_2
      MODULE PROCEDURE :: fft_allocate_real_3
      MODULE PROCEDURE :: fft_allocate_complex_1
      MODULE PROCEDURE :: fft_allocate_complex_2
      MODULE PROCEDURE :: fft_allocate_complex_3
   END INTERFACE fft_allocate

   INTERFACE fft_deallocate
      MODULE PROCEDURE :: fft_free_real_1
      MODULE PROCEDURE :: fft_free_real_2
      MODULE PROCEDURE :: fft_free_real_3
      MODULE PROCEDURE :: fft_free_complex_1
      MODULE PROCEDURE :: fft_free_complex_2
      MODULE PROCEDURE :: fft_free_complex_3
   END INTERFACE fft_deallocate

   INTERFACE
      SUBROUTINE fft_allocate_double_c(length, buffer) BIND(C, name="fft_allocate_double_F")
         IMPORT :: C_PTR, C_INT
         INTEGER(C_INT), VALUE :: length
         TYPE(C_PTR) :: buffer
      END SUBROUTINE fft_allocate_double_c
   END INTERFACE

   INTERFACE
      SUBROUTINE fft_allocate_complex_c(length, buffer) BIND(C, name="fft_allocate_complex_F")
         IMPORT :: C_PTR, C_INT
         INTEGER(C_INT), VALUE :: length
         TYPE(C_PTR) :: buffer
      END SUBROUTINE fft_allocate_complex_c
   END INTERFACE

   INTERFACE
      SUBROUTINE fft_free_double_c(buffer) BIND(C, name="fft_free_double_F")
         IMPORT :: C_PTR
         TYPE(C_PTR), VALUE :: buffer
      END SUBROUTINE fft_free_double_c
   END INTERFACE

   INTERFACE
      SUBROUTINE fft_free_complex_c(buffer) BIND(C, name="fft_free_complex_F")
         IMPORT :: C_PTR
         TYPE(C_PTR), VALUE :: buffer
      END SUBROUTINE fft_free_complex_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Initialize FFT library
!> \param backend ...
!> \param planning_mode ...
!> \param use_guru_interface ...
!> \param wisdom_file ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_library_init(backend, planning_mode, use_guru_interface, wisdom_file, status)
      INTEGER, INTENT(IN)                                :: backend, planning_mode
      LOGICAL, INTENT(IN)                                :: use_guru_interface
      CHARACTER(LEN=*), INTENT(IN)                       :: wisdom_file
      INTEGER, INTENT(OUT)                               :: status

      CHARACTER(KIND=C_CHAR, LEN=1), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: wisdom_file_c
      INTEGER                                            :: i
      INTERFACE
         INTEGER FUNCTION fft_library_init_c(backend, planning_mode, use_guru_interface, &
                                             wisdom_file) BIND(C, name="fft_library_init_F")
            IMPORT :: C_INT, C_BOOL, C_PTR
            INTEGER(C_INT), VALUE :: backend, planning_mode
            LOGICAL(C_BOOL), VALUE :: use_guru_interface
            TYPE(C_PTR), VALUE :: wisdom_file

         END FUNCTION fft_library_init_c
      END INTERFACE

      ALLOCATE (wisdom_file_c(LEN_TRIM(wisdom_file) + 1))
      DO i = 1, LEN_TRIM(wisdom_file)
         wisdom_file_c(i) = wisdom_file(i:i)
      END DO
      wisdom_file_c(LEN_TRIM(wisdom_file) + 1) = C_NULL_CHAR

      status = INT(fft_library_init_c(backend, planning_mode, &
                                      LOGICAL(use_guru_interface, C_BOOL), C_LOC(wisdom_file_c)))

      DEALLOCATE (wisdom_file_c)

   END SUBROUTINE fft_library_init

! **************************************************************************************************
!> \brief Finalize FFT library
!> \param wisdom_file ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_library_finalize(wisdom_file)
      CHARACTER(LEN=*), INTENT(IN)                       :: wisdom_file

      CHARACTER(KIND=C_CHAR, LEN=1), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: wisdom_file_c
      INTEGER                                            :: i
      INTERFACE
         SUBROUTINE fft_library_finalize_c(wisdom_file) BIND(C, name="fft_library_finalize_F")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: wisdom_file
         END SUBROUTINE fft_library_finalize_c
      END INTERFACE

      ALLOCATE (wisdom_file_c(LEN_TRIM(wisdom_file) + 1))
      DO i = 1, LEN_TRIM(wisdom_file)
         wisdom_file_c(i) = wisdom_file(i:i)
      END DO
      wisdom_file_c(LEN_TRIM(wisdom_file) + 1) = C_NULL_CHAR

      CALL fft_library_finalize_c(C_LOC(wisdom_file_c))

      DEALLOCATE (wisdom_file_c)

   END SUBROUTINE fft_library_finalize

! **************************************************************************************************
!> \brief Print FFT library statistics
!> \param threshold ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_print_timing_report(threshold)
      REAL(KIND=dp), INTENT(IN)                          :: threshold

      INTERFACE
         SUBROUTINE fft_print_timing_report_c(threshold) &
            BIND(C, name="fft_print_timing_report_F")
            IMPORT :: C_DOUBLE
            REAL(C_DOUBLE), VALUE :: threshold
         END SUBROUTINE fft_print_timing_report_c
      END INTERFACE

      CALL fft_print_timing_report_c(threshold)

   END SUBROUTINE fft_print_timing_report

! **************************************************************************************************
!> \brief Inquire the library backend in use
!> \return ...
!> \author Frederick Stein
! **************************************************************************************************
   INTEGER FUNCTION fft_backend_in_use()

      INTERFACE
         INTEGER(C_INT) FUNCTION fft_backend_in_use_c() &
            BIND(C, name="fft_backend_in_use_F")
            IMPORT :: C_INT
         END FUNCTION fft_backend_in_use_c
      END INTERFACE

      fft_backend_in_use = fft_backend_in_use_c()

   END FUNCTION fft_backend_in_use

! **************************************************************************************************
!> \brief Inquire whether a guru interface is available
!> \return ...
!> \author Frederick Stein
! **************************************************************************************************
   LOGICAL FUNCTION fft_has_guru_interface()

      INTERFACE
         LOGICAL(C_BOOL) FUNCTION fft_has_guru_interface_c() &
            BIND(C, name="fft_has_guru_interface_F")
            IMPORT :: C_BOOL
         END FUNCTION fft_has_guru_interface_c
      END INTERFACE

      fft_has_guru_interface = fft_has_guru_interface_c()

   END FUNCTION fft_has_guru_interface

   #:set maxdim = 3
   #:for dim in range(1, maxdim+1)
      #:for type in ["COMPLEX", "REAL"]
! **************************************************************************************************
!> \brief Allocate an array for use with the FFT library
!> \param buffer ...
!> \param length ...
!> \author Frederick Stein
! **************************************************************************************************
         SUBROUTINE fft_allocate_${type}$_${dim}$ (buffer, length)
            ${type}$ (KIND=dp), CONTIGUOUS, DIMENSION(:${", :"*(dim-1)}$), &
               POINTER                                         :: buffer
            INTEGER, DIMENSION(${dim}$), INTENT(IN)                                :: length

            TYPE(C_PTR)                                        :: buffer_c

            buffer_c = C_NULL_PTR

            CALL fft_allocate_complex_c(INT(PRODUCT(length), C_INT), buffer_c)

            CALL C_F_POINTER(buffer_c, buffer, length)

         END SUBROUTINE fft_allocate_${type}$_${dim}$
      #:endfor
   #:endfor

! **************************************************************************************************
!> \brief Free a complex array allocated with fft_allocate_complex
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_complex_1(buffer)
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:), &
         POINTER                                         :: buffer

      CALL fft_free_complex_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_complex_1

! **************************************************************************************************
!> \brief Free a complex array allocated with fft_allocate_complex
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_complex_2(buffer)
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: buffer

      CALL fft_free_complex_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_complex_2

! **************************************************************************************************
!> \brief Free a complex array allocated with fft_allocate_complex
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_complex_3(buffer)
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: buffer

      CALL fft_free_complex_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_complex_3

! **************************************************************************************************
!> \brief Free a double array allocated with fft_allocate_double
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_real_1(buffer)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:), POINTER   :: buffer

      CALL fft_free_double_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_real_1

! **************************************************************************************************
!> \brief Free a double array allocated with fft_allocate_double
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_real_2(buffer)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: buffer

      CALL fft_free_double_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_real_2

! **************************************************************************************************
!> \brief Free a double array allocated with fft_allocate_double
!> \param buffer ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_free_real_3(buffer)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: buffer

      CALL fft_free_double_c(C_LOC(buffer))

      NULLIFY (buffer)

   END SUBROUTINE fft_free_real_3

! **************************************************************************************************
!> \brief Registers a complex-to-complex 1D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \param leading_dimension_in ...
!> \param leading_dimension_out ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_1d_local(direction, fft_size, number_of_ffts, transpose_in, transpose_out, &
                                    leading_dimension_in, leading_dimension_out, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out
      INTEGER, INTENT(IN)                                :: leading_dimension_in, leading_dimension_out
      COMPLEX(KIND=dp), DIMENSION(leading_dimension_in*MERGE(fft_size, number_of_ffts, transpose_in)), &
         INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(leading_dimension_out*MERGE(fft_size, number_of_ffts, transpose_out)), &
         INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_register_1d_local_c(direction, fft_size, number_of_ffts, transpose_in, transpose_out, &
                                            leading_dimension_in, leading_dimension_out, grid_in, grid_out) &
            BIND(C, name="fft_register_1d_local_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), VALUE :: fft_size, number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose_in, transpose_out
            INTEGER(C_INT), VALUE :: leading_dimension_in, leading_dimension_out
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_register_1d_local_c
      END INTERFACE

      CALL fft_register_1d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                   INT(number_of_ffts, C_INT), LOGICAL(transpose_in, C_BOOL), &
                                   LOGICAL(transpose_out, C_BOOL), leading_dimension_in, leading_dimension_out, &
                                   C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_register_1d_local

! **************************************************************************************************
!> \brief Registers a complex-to-complex 2D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_2d_local(direction, fft_size, number_of_ffts, &
                                    transpose_in, transpose_out, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_register_2d_local_c(direction, fft_size, number_of_ffts, &
                                            transpose_in, transpose_out, grid_in, grid_out) &
            BIND(C, name="fft_register_2d_local_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(2) :: fft_size
            INTEGER(C_INT), VALUE :: number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose_in, transpose_out
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_register_2d_local_c
      END INTERFACE

      CALL fft_register_2d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                   INT(number_of_ffts, C_INT), LOGICAL(transpose_in, C_BOOL), &
                                   LOGICAL(transpose_out, C_BOOL), C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_register_2d_local

! **************************************************************************************************
!> \brief Registers a complex-to-complex 3D FFT
!> \param direction ...
!> \param fft_size ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_3d_local(direction, fft_size, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_register_3d_local_c(direction, fft_size, grid_in, grid_out) &
            BIND(C, name="fft_register_3d_local_F")
            IMPORT :: C_INT, C_DOUBLE_COMPLEX, C_BOOL, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(3), INTENT(IN) :: fft_size
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_register_3d_local_c
      END INTERFACE

      CALL fft_register_3d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                   C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_register_3d_local

! **************************************************************************************************
!> \brief Registers a complex-to-complex 1D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \param leading_dimension ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_1d_local_inplace(direction, fft_size, number_of_ffts, &
                                            transpose, leading_dimension, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose
      INTEGER, INTENT(IN)                                :: leading_dimension
      COMPLEX(KIND=dp), DIMENSION(leading_dimension*MERGE(fft_size, number_of_ffts, transpose)), INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_register_1d_local_inplace_c(direction, fft_size, number_of_ffts, &
                                                    transpose, leading_dimension, grid) &
            BIND(C, name="fft_register_1d_local_inplace_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), VALUE :: fft_size, number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose
            INTEGER(C_INT), VALUE :: leading_dimension
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_register_1d_local_inplace_c
      END INTERFACE

      CALL fft_register_1d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                           INT(number_of_ffts, C_INT), &
                                           LOGICAL(transpose, C_BOOL), leading_dimension, C_LOC(grid))

   END SUBROUTINE fft_register_1d_local_inplace

! **************************************************************************************************
!> \brief Registers a complex-to-complex 2D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_2d_local_inplace(direction, fft_size, number_of_ffts, transpose, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_register_2d_local_inplace_c(direction, fft_size, number_of_ffts, transpose, grid) &
            BIND(C, name="fft_register_2d_local_inplace_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(2) :: fft_size
            INTEGER(C_INT), VALUE :: number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_register_2d_local_inplace_c
      END INTERFACE

      CALL fft_register_2d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                           INT(number_of_ffts, C_INT), &
                                           LOGICAL(transpose, C_BOOL), C_LOC(grid))

   END SUBROUTINE fft_register_2d_local_inplace

! **************************************************************************************************
!> \brief Registers a complex-to-complex 3D FFT
!> \param direction ...
!> \param fft_size ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_register_3d_local_inplace(direction, fft_size, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_register_3d_local_inplace_c(direction, fft_size, grid) &
            BIND(C, name="fft_register_3d_local_inplace_F")
            IMPORT :: C_INT, C_DOUBLE_COMPLEX, C_BOOL, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(3), INTENT(IN) :: fft_size
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_register_3d_local_inplace_c
      END INTERFACE

      CALL fft_register_3d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), C_LOC(grid))

   END SUBROUTINE fft_register_3d_local_inplace

! **************************************************************************************************
!> \brief Performs a complex-to-complex 1D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \param leading_dimension_in ...
!> \param leading_dimension_out ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_1d_local(direction, fft_size, number_of_ffts, transpose_in, transpose_out, &
                           leading_dimension_in, leading_dimension_out, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out
      INTEGER, INTENT(IN)                                :: leading_dimension_in, leading_dimension_out
      COMPLEX(KIND=dp), DIMENSION(leading_dimension_in*MERGE(fft_size, number_of_ffts, transpose_in)), &
         INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(leading_dimension_out*MERGE(fft_size, number_of_ffts, transpose_out)), &
         INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_1d_local_c(direction, fft_size, number_of_ffts, transpose_in, transpose_out, &
                                   leading_dimension_in, leading_dimension_out, grid_in, grid_out) &
            BIND(C, name="fft_1d_local_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), VALUE :: fft_size, number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose_in, transpose_out
            INTEGER(C_INT), VALUE :: leading_dimension_in, leading_dimension_out
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_1d_local_c
      END INTERFACE

      CALL fft_1d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                          INT(number_of_ffts, C_INT), LOGICAL(transpose_in, C_BOOL), &
                          LOGICAL(transpose_out, C_BOOL), leading_dimension_in, leading_dimension_out, &
                          C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_1d_local

! **************************************************************************************************
!> \brief Performs a complex-to-complex 2D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose_in ...
!> \param transpose_out ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_2d_local(direction, fft_size, number_of_ffts, transpose_in, transpose_out, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose_in, transpose_out
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_2d_local_c(direction, fft_size, number_of_ffts, transpose_in, transpose_out, grid_in, grid_out) &
            BIND(C, name="fft_2d_local_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(2) :: fft_size
            INTEGER(C_INT), VALUE :: number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose_in, transpose_out
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_2d_local_c
      END INTERFACE

      CALL fft_2d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                          INT(number_of_ffts, C_INT), LOGICAL(transpose_in, C_BOOL), &
                          LOGICAL(transpose_out, C_BOOL), C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_2d_local

! **************************************************************************************************
!> \brief Performs a complex-to-complex 3D FFT
!> \param direction ...
!> \param fft_size ...
!> \param grid_in ...
!> \param grid_out ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_3d_local(direction, fft_size, grid_in, grid_out)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET      :: grid_in
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(OUT), TARGET        :: grid_out

      INTERFACE
         SUBROUTINE fft_3d_local_c(direction, fft_size, grid_in, grid_out) &
            BIND(C, name="fft_3d_local_F")
            IMPORT :: C_INT, C_DOUBLE_COMPLEX, C_BOOL, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(3), INTENT(IN) :: fft_size
            TYPE(C_PTR), VALUE :: grid_in, grid_out
         END SUBROUTINE fft_3d_local_c
      END INTERFACE

      CALL fft_3d_local_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), C_LOC(grid_in), C_LOC(grid_out))

   END SUBROUTINE fft_3d_local

! **************************************************************************************************
!> \brief Performs a complex-to-complex 1D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \param leading_dimension ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_1d_local_inplace(direction, fft_size, number_of_ffts, transpose, leading_dimension, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, INTENT(IN)                                :: fft_size, number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose
      INTEGER, INTENT(IN)                                :: leading_dimension
      COMPLEX(KIND=dp), DIMENSION(leading_dimension*MERGE(fft_size, number_of_ffts, transpose)), &
         INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_1d_local_inplace_c(direction, fft_size, number_of_ffts, transpose, leading_dimension, grid) &
            BIND(C, name="fft_1d_local_inplace_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), VALUE :: fft_size, number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose
            INTEGER(C_INT), VALUE :: leading_dimension
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_1d_local_inplace_c
      END INTERFACE

      CALL fft_1d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                  INT(number_of_ffts, C_INT), &
                                  LOGICAL(transpose, C_BOOL), leading_dimension, C_LOC(grid))

   END SUBROUTINE fft_1d_local_inplace

! **************************************************************************************************
!> \brief Performs a complex-to-complex 2D FFT
!> \param direction ...
!> \param fft_size ...
!> \param number_of_ffts ...
!> \param transpose ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_2d_local_inplace(direction, fft_size, number_of_ffts, transpose, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(2), INTENT(IN)                  :: fft_size
      INTEGER, INTENT(IN)                                :: number_of_ffts
      LOGICAL, INTENT(IN)                                :: transpose
      COMPLEX(KIND=dp), DIMENSION(PRODUCT(fft_size)*number_of_ffts), INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_2d_local_inplace_c(direction, fft_size, number_of_ffts, transpose, grid) &
            BIND(C, name="fft_2d_local_inplace_F")
            IMPORT :: C_INT, C_BOOL, C_DOUBLE_COMPLEX, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(2) :: fft_size
            INTEGER(C_INT), VALUE :: number_of_ffts
            LOGICAL(C_BOOL), VALUE :: transpose
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_2d_local_inplace_c
      END INTERFACE

      CALL fft_2d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), &
                                  INT(number_of_ffts, C_INT), LOGICAL(transpose, C_BOOL), C_LOC(grid))

   END SUBROUTINE fft_2d_local_inplace

! **************************************************************************************************
!> \brief Performs a complex-to-complex 3D FFT
!> \param direction ...
!> \param fft_size ...
!> \param grid ...
!> \author Frederick Stein
! **************************************************************************************************
   SUBROUTINE fft_3d_local_inplace(direction, fft_size, grid)
      LOGICAL, INTENT(IN)                                :: direction
      INTEGER, DIMENSION(3), INTENT(IN)                  :: fft_size
      COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET      :: grid

      INTERFACE
         SUBROUTINE fft_3d_local_inplace_c(direction, fft_size, grid) &
            BIND(C, name="fft_3d_local_inplace_F")
            IMPORT :: C_INT, C_DOUBLE_COMPLEX, C_BOOL, C_PTR
            LOGICAL(C_BOOL), VALUE :: direction
            INTEGER(C_INT), DIMENSION(3), INTENT(IN) :: fft_size
            TYPE(C_PTR), VALUE :: grid
         END SUBROUTINE fft_3d_local_inplace_c
      END INTERFACE

      CALL fft_3d_local_inplace_c(LOGICAL(direction, C_BOOL), INT(fft_size, C_INT), C_LOC(grid))

   END SUBROUTINE fft_3d_local_inplace

END MODULE fft_api
