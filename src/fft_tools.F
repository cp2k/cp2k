!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!
! How to add a new FFT library:
!  - create a new interface library : fftXX_lib with the entries
!    fft3d and fft_get_lengths, see fftw_lib for a template
!  - add in this file the entries to the new library; in each
!    subroutine there will be an additional CASE
!

MODULE fft_tools
  
  USE kinds, ONLY: dbl, sgl
  USE stop_program, ONLY : stop_prg
  USE fftw_lib, ONLY : w_fft3d => fft3d, &
                       w_fft_get_lengths => fft_get_lengths
  USE fftsg_lib, ONLY : sg_fft3d => fft3d, &
                        sg_fft_get_lengths => fft_get_lengths
  
  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: init_fft, get_fft_library, fft3d
  PUBLIC :: fft_radix_operations
  PUBLIC :: FWFFT, BWFFT
  PUBLIC :: FFT_RADIX_CLOSEST, FFT_RADIX_NEXT, FFT_RADIX_ALLOWED
  PUBLIC :: FFT_RADIX_DISALLOWED
  
  INTEGER, PARAMETER :: FWFFT = +1, BWFFT = -1
  INTEGER, PARAMETER :: FFT_RADIX_CLOSEST = 493, FFT_RADIX_NEXT = 494
  INTEGER, PARAMETER :: FFT_RADIX_ALLOWED = 495, FFT_RADIX_DISALLOWED = 496
  
  INTEGER :: fft_type = 0

!******************************************************************************

CONTAINS
  
!******************************************************************************

SUBROUTINE init_fft ( fftlib, fftnum )

  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: fftlib
  INTEGER, INTENT ( IN ), OPTIONAL :: fftnum

  INTEGER :: n(3) = 4
  COMPLEX ( dbl ) :: zz ( 4, 4, 4 )
  INTEGER :: stat, i

  IF ( PRESENT ( fftlib ) ) THEN
    SELECT CASE ( fftlib )
    CASE DEFAULT
      CALL stop_prg ("init_fft","Unknown FFT library")
    CASE ( "FFTSG" )
      fft_type = 1
    CASE ( "FFTW" )
      fft_type = 2
    END SELECT
  ELSE IF ( PRESENT ( fftnum ) ) THEN
    fft_type = fftnum
    zz = 0.1_dbl
    CALL fft3d ( 1, n, zz, status=stat )
    IF ( stat /= 0 ) call stop_prg ("init_fft","FFT library not available")
  ELSE
    zz = 0.1_dbl
    DO i = 1, 100
      fft_type = i
      CALL fft3d ( 1, n, zz, status=stat )
      IF ( stat == 0 ) EXIT
      fft_type = 0
    END DO
    IF (fft_type == 0 ) call stop_prg ("init_fft","No FFT library available")
  END IF

  IF ( PRESENT ( fftnum ) ) THEN
    IF ( fft_type /= fftnum ) call stop_prg ("init_fft",&
       " Inconsistent Arguments ")
  END IF

END SUBROUTINE init_fft

!******************************************************************************

SUBROUTINE get_fft_library ( fft_handle, library )

  INTEGER, INTENT ( OUT ) :: fft_handle
  CHARACTER ( LEN = * ), INTENT ( OUT ) :: library

  SELECT CASE ( fft_type )
  CASE DEFAULT
    library = " No FFT library initialized "
    fft_handle = 0
  CASE ( 1 )
    library = " FFTsg library initialized "
    fft_handle = 1
  CASE ( 2 )
    library = " FFTw library initialized "
    fft_handle = 2
  END SELECT

END SUBROUTINE get_fft_library

!******************************************************************************

SUBROUTINE fft_radix_operations ( radix_in, radix_out, operation )
!! Determine the allowed lengths of FFT's   '''
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: radix_in, operation
  INTEGER, INTENT ( OUT ) :: radix_out
  
! Locals
  INTEGER :: i, iloc
  INTEGER, DIMENSION ( : ), POINTER :: data
  
!------------------------------------------------------------------------------
  
  SELECT CASE ( fft_type )
  CASE ( 1 )
    CALL sg_fft_get_lengths ( data, max_length = 1024 )
  CASE ( 2 )
    CALL w_fft_get_lengths ( data, max_length = 1024 )
  CASE DEFAULT
    CALL stop_prg ("fft3d","Unknown FFT library")
  END SELECT
  
  iloc = 0
  DO i = 1, SIZE ( data )
     IF ( data ( i ) == radix_in ) THEN
        iloc = i
        EXIT
     ELSE
        IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
           CYCLE
        ELSE IF ( data ( i ) > radix_in ) THEN
           iloc = i
           EXIT
        END IF
     END IF
  END DO
  
  IF ( iloc == 0 ) THEN
     CALL stop_prg ( "fft_radix_operations", "index to radix array not found" )
  END IF
  
  IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
     IF ( data ( iloc ) == radix_in ) THEN
        radix_out = FFT_RADIX_ALLOWED
     ELSE
        radix_out = FFT_RADIX_DISALLOWED
     END IF
     
  ELSE IF ( OPERATION == FFT_RADIX_CLOSEST ) THEN
     IF ( data ( iloc ) == radix_in ) THEN
        radix_out = data ( iloc )
     ELSE
        IF ( ABS ( data ( iloc - 1 ) - radix_in ) <= &
             ABS ( data ( iloc ) - radix_in ) ) THEN
           radix_out = data ( iloc - 1 )
        ELSE
           radix_out = data ( iloc )
        END IF
     END IF
     
  ELSE IF ( OPERATION == FFT_RADIX_NEXT ) THEN
     radix_out = data ( iloc )
     
  ELSE
     CALL stop_prg ( "fft_radix_operations", "disallowed radix operation" )
  END IF
  
END SUBROUTINE fft_radix_operations

!******************************************************************************

SUBROUTINE fft3d ( fsign, n, zg, zg_out, scale, status )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: fsign
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: zg
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ), OPTIONAL :: zg_out
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  INTEGER, INTENT ( OUT ), OPTIONAL :: status

  REAL ( dbl ) :: norm
  INTEGER :: sign

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1._dbl
  END IF

  sign = fsign

  SELECT CASE ( fft_type )
  CASE ( 1 )
    IF ( PRESENT ( zg_out ) ) THEN
      CALL sg_fft3d ( sign, norm, n, zg, zg_out )
    ELSE
      CALL sg_fft3d ( sign, norm, n, zg )
    END IF
  CASE ( 2 )
    IF ( PRESENT ( zg_out ) ) THEN
      CALL w_fft3d ( sign, norm, n, zg, zg_out )
    ELSE
      CALL w_fft3d ( sign, norm, n, zg )
    END IF
  CASE DEFAULT
    CALL stop_prg ("fft3d","Unknown FFT library")
  END SELECT

  IF ( PRESENT ( status ) ) THEN
    IF ( sign == 0 ) THEN
      status = 1
    ELSE
      status = 0
    END IF
  END IF

END SUBROUTINE fft3d 

!******************************************************************************

END MODULE fft_tools

!******************************************************************************
