!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on INTEL
!> FPGAs using OpenCL. It depends at execution time on the board support
!> packages of the specific FPGA
!> \author Arjun Ramaswami
!> \author Robert Schade
! **************************************************************************************************

MODULE pw_fpga
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_CHAR,&
                                              C_DOUBLE_COMPLEX,&
                                              C_FLOAT_COMPLEX,&
                                              C_INT
   USE cp_files,                        ONLY: get_data_dir
   USE kinds,                           ONLY: dp,&
                                              sp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_fpga_init, pw_fpga_finalize
   PUBLIC :: pw_fpga_check_bitstream
   PUBLIC :: pw_fpga_r3dc1d_3d_sp, pw_fpga_c1dr3d_3d_sp
   PUBLIC :: pw_fpga_r3dc1d_3d_dp, pw_fpga_c1dr3d_3d_dp

   INTERFACE
! **************************************************************************************************
!> \brief Initialize FPGA
!> \retval status if the routine failed or not
! **************************************************************************************************
      FUNCTION pw_fpga_initialize() RESULT(stat) &
         BIND(C, name="pw_fpga_initialize_")
         IMPORT
         INTEGER(KIND=C_INT)                    :: stat
      END FUNCTION pw_fpga_initialize

! **************************************************************************************************
!> \brief Destroy FPGA
! **************************************************************************************************
      SUBROUTINE pw_fpga_final() &
         BIND(C, name="pw_fpga_final_")
      END SUBROUTINE pw_fpga_final

   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief Check whether an fpga bitstream for the given FFT3d size is present
!> \param npts ...
!> \return res ...
! **************************************************************************************************
      FUNCTION pw_fpga_check_bitstream(n) RESULT(res) &
         BIND(C, name="pw_fpga_check_bitstream_")
         IMPORT
         INTEGER(KIND=C_INT), DIMENSION(*), &
            INTENT(IN)                            :: n
         LOGICAL(KIND=C_BOOL)                     :: res

      END FUNCTION pw_fpga_check_bitstream

   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief single precision FFT3d using FPGA
!> \param data_path_len - length of the path to the data directory
!> \param data_path - path to the data directory
!> \param dir - direction of FFT3d
!> \param npts - dimensions of FFT3d
!> \param single precision c_in...
! **************************************************************************************************
      SUBROUTINE pw_fpga_fft3d_sp(data_path_len, data_path, dir, n, c_in_sp) &
         BIND(C, name="pw_fpga_fft3d_sp_")
         IMPORT
         INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: data_path_len
         CHARACTER(len=1, kind=C_CHAR), TARGET    :: data_path(data_path_len)
         INTEGER(KIND=C_INT), VALUE, INTENT(IN)   :: dir
         INTEGER(KIND=C_INT), DIMENSION(*), &
            INTENT(IN)                            :: n
         COMPLEX(KIND=C_FLOAT_COMPLEX)           :: c_in_sp(n(1), n(2), n(3))
      END SUBROUTINE pw_fpga_fft3d_sp
   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief double precision FFT3d using FPGA
!> \param data_path_len - length of the path to the data directory
!> \param data_path - path to the data directory
!> \param dir - direction of FFT3d
!> \param npts - dimensions of FFT3d
!> \param double precision c_in...
! **************************************************************************************************
      SUBROUTINE pw_fpga_fft3d_dp(data_path_len, data_path, dir, n, c_in_dp) &
         BIND(C, name="pw_fpga_fft3d_dp_")
         IMPORT
         INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: data_path_len
         CHARACTER(len=1, kind=C_CHAR), TARGET    :: data_path(data_path_len)
         INTEGER(KIND=C_INT), VALUE, INTENT(IN)   :: dir
         INTEGER(KIND=C_INT), DIMENSION(*), &
            INTENT(IN)                            :: n
         COMPLEX(KIND=C_DOUBLE_COMPLEX)           :: c_in_dp(n(1), n(2), n(3))
      END SUBROUTINE pw_fpga_fft3d_dp
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocates resources on the fpga device
! **************************************************************************************************
   SUBROUTINE pw_fpga_init()
#if defined (__PW_FPGA)
      INTEGER :: stat

#if defined (__PW_CUDA)
#error "CUDA and FPGA cannot be configured concurrently!"
#endif

      stat = pw_fpga_initialize()
      IF (stat /= 0) &
         CPABORT("pw_fpga_init: failed")
#endif
   END SUBROUTINE pw_fpga_init

! **************************************************************************************************
!> \brief Releases resources on the fpga device
! **************************************************************************************************
   SUBROUTINE pw_fpga_finalize()
#if defined (__PW_FPGA)
      CALL pw_fpga_final()
#endif
   END SUBROUTINE pw_fpga_finalize

! **************************************************************************************************
!> \brief perform an in-place double precision fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_r3dc1d_3d_dp(n, c_out)
      INTEGER, DIMENSION(:), INTENT(IN)          :: n
      COMPLEX(KIND=dp), INTENT(INOUT)            :: c_out(n(1), n(2), n(3))

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)
#else
      INTEGER                                     :: handle3
      CHARACTER(len=100)                          :: data_path
      INTEGER                                     :: data_path_len

      CHARACTER(len=*), PARAMETER :: routineX = 'fw_fft_fpga_r3dc1d_dp'

      data_path = TRIM(get_data_dir())
      data_path_len = LEN_TRIM(data_path)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_dp(data_path_len, data_path, +1, n, c_out)
      CALL timestop(handle3)

#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief perform an in-place double precision inverse fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_c1dr3d_3d_dp(n, c_out)
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n
      COMPLEX(KIND=dp), INTENT(INOUT)    :: c_out(n(1), n(2), n(3))

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)
#else
      INTEGER                                          :: handle3
      CHARACTER(len=100)                               :: data_path
      INTEGER                                          :: data_path_len

      CHARACTER(len=*), PARAMETER :: routineX = 'bw_fft_fpga_c1dr3d_dp'

      data_path = TRIM(get_data_dir())
      data_path_len = LEN_TRIM(data_path)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_dp(data_path_len, data_path, -1, n, c_out)
      CALL timestop(handle3)

#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief perform an in-place single precision fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_r3dc1d_3d_sp(n, c_out)
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n
      COMPLEX(KIND=dp), INTENT(INOUT)             :: c_out(n(1), n(2), n(3))

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)
#else
      COMPLEX, DIMENSION(:, :, :), POINTER             :: c_in_sp
      INTEGER                                          :: handle3
      CHARACTER(len=100)                               :: data_path
      INTEGER                                          :: data_path_len

      CHARACTER(len=*), PARAMETER :: routineX = 'fw_fft_fpga_r3dc1d_sp'

      data_path = TRIM(get_data_dir())
      data_path_len = LEN_TRIM(data_path)

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp = CMPLX(c_out, KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(data_path_len, data_path, +1, n, c_in_sp)
      CALL timestop(handle3)

      ! typecast sp back to dp
      !c_out = CMPLX(real(c_in_sp), 0.0_dp, KIND=dp)
      c_out = CMPLX(c_in_sp, KIND=dp)

      DEALLOCATE (c_in_sp)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief perform an in-place single precision inverse fft3d on the FPGA
!> \param n ...
!> \param c_out  ...
! **************************************************************************************************
   SUBROUTINE pw_fpga_c1dr3d_3d_sp(n, c_out)
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n
      COMPLEX(KIND=dp), INTENT(INOUT)             :: c_out(n(1), n(2), n(3))

#if ! defined (__PW_FPGA)
      MARK_USED(c_out)
      MARK_USED(n)

#else
      COMPLEX, DIMENSION(:, :, :), POINTER             :: c_in_sp
      INTEGER                                          :: handle3
      CHARACTER(len=100)                               :: data_path
      INTEGER                                          :: data_path_len

      CHARACTER(len=*), PARAMETER :: routineX = 'bw_fft_fpga_c1dr3d_sp'

      data_path = TRIM(get_data_dir())
      data_path_len = LEN_TRIM(data_path)

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp = CMPLX(c_out, KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(data_path_len, data_path, -1, n, c_in_sp)
      CALL timestop(handle3)

      ! typecast sp back to dp
      c_out = CMPLX(c_in_sp, KIND=dp)

      DEALLOCATE (c_in_sp)
#endif
   END SUBROUTINE

END MODULE pw_fpga

