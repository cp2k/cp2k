!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE_COMPLEX,&
                                              C_FLOAT_COMPLEX,&
                                              C_INT,&
                                              C_NULL_CHAR
   USE cp_files,                        ONLY: get_data_dir
   USE kinds,                           ONLY: dp,&
                                              sp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_fpga_init, pw_fpga_finalize
   PUBLIC :: pw_fpga_init_bitstream
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
!> \brief Check whether an fpga bitstream for the given FFT3d size is present & load binary if needed
!> \param data_path - path to the data directory
!> \param npts - fft3d size
!> \return res - true if fft3d size supported
! **************************************************************************************************
      FUNCTION pw_fpga_check_bitstream(data_path, n) RESULT(res) &
         BIND(C, name="pw_fpga_check_bitstream_")
         IMPORT
         CHARACTER(KIND=C_CHAR)        :: data_path(*)
         INTEGER(KIND=C_INT)           :: n(3)
         INTEGER(KIND=C_INT)           :: res
      END FUNCTION pw_fpga_check_bitstream

   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief single precision FFT3d using FPGA
!> \param dir - direction of FFT3d
!> \param npts - dimensions of FFT3d
!> \param single precision c_in...
! **************************************************************************************************
      SUBROUTINE pw_fpga_fft3d_sp(dir, n, c_in_sp) &
         BIND(C, name="pw_fpga_fft3d_sp_")
         IMPORT
         INTEGER(KIND=C_INT), VALUE              :: dir
         INTEGER(KIND=C_INT)                     :: n(3)
         COMPLEX(KIND=C_FLOAT_COMPLEX)           :: c_in_sp(n(1), n(2), n(3))
      END SUBROUTINE pw_fpga_fft3d_sp
   END INTERFACE

   INTERFACE
! **************************************************************************************************
!> \brief double precision FFT3d using FPGA
!> \param dir - direction of FFT3d
!> \param npts - dimensions of FFT3d
!> \param double precision c_in...
! **************************************************************************************************
      SUBROUTINE pw_fpga_fft3d_dp(dir, n, c_in_dp) &
         BIND(C, name="pw_fpga_fft3d_dp_")
         IMPORT
         INTEGER(KIND=C_INT), VALUE              :: dir
         INTEGER(KIND=C_INT)                     :: n(3)
         COMPLEX(KIND=C_DOUBLE_COMPLEX)          :: c_in_dp(n(1), n(2), n(3))
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
      CPABORT("CUDA and FPGA cannot be configured concurrently!")
#endif
      stat = pw_fpga_initialize()
      IF (stat /= 0) &
         CPABORT("pw_fpga_init: failed")
#endif

#if (__PW_FPGA_SP && !(__PW_FPGA))
#error "Define both __PW_FPGA_SP and __PW_FPGA"
      CPABORT("Define both __PW_FPGA_SP and __PW_FPGA")
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

      CHARACTER(len=*), PARAMETER :: routineX = 'fw_fft_fpga_r3dc1d_dp'

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_dp(+1, n, c_out)
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

      CHARACTER(len=*), PARAMETER :: routineX = 'bw_fft_fpga_c1dr3d_dp'

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_dp(-1, n, c_out)
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

      CHARACTER(len=*), PARAMETER :: routineX = 'fw_fft_fpga_r3dc1d_sp'

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp = CMPLX(c_out, KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(+1, n, c_in_sp)
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

      CHARACTER(len=*), PARAMETER :: routineX = 'bw_fft_fpga_c1dr3d_sp'

      ALLOCATE (c_in_sp(n(1), n(2), n(3)))
      ! pointer to single precision complex array
      c_in_sp = CMPLX(c_out, KIND=sp)

      CALL timeset(routineX, handle3)
      CALL pw_fpga_fft3d_sp(-1, n, c_in_sp)
      CALL timestop(handle3)

      ! typecast sp back to dp
      c_out = CMPLX(c_in_sp, KIND=dp)

      DEALLOCATE (c_in_sp)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief  Invoke the pw_fpga_check_bitstream C function passing the path to the data dir
!> \param  n   - fft3d size
!> \return ...
!> \retval res - true if fft size found and initialized else false
! **************************************************************************************************
   FUNCTION pw_fpga_init_bitstream(n) RESULT(res)
      INTEGER, DIMENSION(:), INTENT(IN)                 :: n
      INTEGER                                           :: res

#if ! defined (__PW_FPGA)
      res = 0
      MARK_USED(n)
      MARK_USED(res)
#else
      CHARACTER(len=100)                               :: data_path
      INTEGER                                          :: data_path_len

      data_path = TRIM(get_data_dir())//C_NULL_CHAR
      data_path_len = LEN_TRIM(data_path)

      res = pw_fpga_check_bitstream(data_path, n)
#endif
   END FUNCTION

END MODULE pw_fpga

