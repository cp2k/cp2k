!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fft_tools [1.0] *
!!
!!   NAME
!!     fft_tools
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (5-Jan-2001): Added SGI library FFT
!!
!!   NOTES
!!     How to add a new FFT library:
!!      - create a new interface library : fftXX_lib with the entries
!!        fft3d and fft_get_lengths, see fftw_lib for a template
!!      - add in this file the entries to the new library; in each
!!        subroutine there will be an additional CASE
!!
!!   SOURCE
!******************************************************************************

MODULE fft_tools
  
  USE kinds, ONLY: dbl, sgl
  USE stop_program, ONLY : stop_prg, stop_memory
  USE fftessl_lib, ONLY : i_fft3d => fft3d, &
                       i_fft_get_lengths => fft_get_lengths
  USE fftw_lib, ONLY : w_fft3d => fft3d, &
                       w_fft_get_lengths => fft_get_lengths
  USE fftsg_lib, ONLY : sg_fft3d => fft3d, sg_mltfft => mltfft, &
                        sg_fft_get_lengths => fft_get_lengths
  USE fftsgi_lib, ONLY : sgi_fft3d => fft3d, &
                         sgi_fft_get_lengths => fft_get_lengths
  USE message_passing, ONLY : mp_environ, mp_comm_compare
  
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

  INTERFACE fft3d
     MODULE PROCEDURE fft3d_s, fft3d_ps, fft3d_pb
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS
  
!******************************************************************************
!!****** fft_tools/init_fft [1.0] *
!!
!!   NAME
!!     init_fft
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     30-Nov-2000 (JGH) ESSL FFT Library added
!!
!!   NOTES
!!     Label 1 : FFTsg
!!     Label 2 : FFTw 
!!     Label 3 : FFT ESSL
!!
!!   SOURCE
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
      CALL stop_prg ("init_fft","Unknown FFT library",fftlib)
    CASE ( "FFTSG" )
      fft_type = 1
    CASE ( "FFTW" )
      fft_type = 2
    CASE ( "FFTESSL" )
      fft_type = 3
    CASE ( "FFTSGI" )
      fft_type = 4
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

!!*****
!******************************************************************************
!!****** fft_tools/get_fft_library [1.0] *
!!
!!   NAME
!!     get_fft_library
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     30-Nov-2000 (JGH) ESSL FFT Library added
!!
!!   SOURCE
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
  CASE ( 3 )
    library = " ESSL FFT library initialized "
    fft_handle = 3
  CASE ( 4 )
    library = " SGI FFT library initialized "
    fft_handle = 4
  END SELECT

END SUBROUTINE get_fft_library

!!*****
!******************************************************************************
!!****** fft_tools/fft_radix_operations [1.0] *
!!
!!   NAME
!!     fft_radix_operations
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Ari Seitsonen
!!
!!   MODIFICATION HISTORY
!!     new library structure (JGH)
!!
!!   SOURCE
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
  CASE ( 3 )
    CALL i_fft_get_lengths ( data, max_length = 1024 )
  CASE ( 4 )
    CALL sgi_fft_get_lengths ( data, max_length = 1024 )
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

!!*****
!******************************************************************************
!!****** fft_tools/fft3d [1.0] *
!!
!!   NAME
!!     fft3d
!!
!!   FUNCTION
!!     Calls the 3D-FFT function from the initialized library
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     30-Nov-2000 (JGH) ESSL FFT Library added
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft3d_s ( fsign, n, zin, zout, scale, status )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: fsign
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: zin
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ), OPTIONAL :: zout
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  INTEGER, INTENT ( OUT ), OPTIONAL :: status

  REAL ( dbl ) :: norm
  INTEGER :: sign

!------------------------------------------------------------------------------

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1._dbl
  END IF

  sign = fsign

  SELECT CASE ( fft_type )
  CASE ( 1 )
    IF ( PRESENT ( zout ) ) THEN
      CALL sg_fft3d ( sign, norm, n, zin, zout )
    ELSE
      CALL sg_fft3d ( sign, norm, n, zin )
    END IF
  CASE ( 2 )
    IF ( PRESENT ( zout ) ) THEN
      CALL w_fft3d ( sign, norm, n, zin, zout )
    ELSE
      CALL w_fft3d ( sign, norm, n, zin )
    END IF
  CASE ( 3 )
    IF ( PRESENT ( zout ) ) THEN
      CALL i_fft3d ( sign, norm, n, zin, zout )
    ELSE
      CALL i_fft3d ( sign, norm, n, zin )
    END IF
  CASE ( 4 )
    IF ( PRESENT ( zout ) ) THEN
      CALL sgi_fft3d ( sign, norm, n, zin, zout )
    ELSE
      CALL sgi_fft3d ( sign, norm, n, zin )
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

END SUBROUTINE fft3d_s

!******************************************************************************

SUBROUTINE fft3d_ps ( fsign, n, cin, gin, rs_group, gs_group, yzp, bounds, &
                      scale, status )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: fsign
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: cin
  COMPLEX ( dbl ), DIMENSION(:,:), INTENT ( INOUT ) :: gin
  INTEGER, INTENT ( IN ) :: rs_group, gs_group
  INTEGER, DIMENSION ( :, :, 0: ), INTENT ( IN ) :: yzp
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: bounds
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  INTEGER, INTENT ( OUT ), OPTIONAL :: status

! Locals
  INTEGER :: numtask, numtask_g, numtask_r, g_pos, r_pos(2), r_dim(2), iout
  INTEGER :: sign, ierr, lg, mg, ng
  COMPLEX ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: sbuf
  REAL ( dbl ) :: norm

!------------------------------------------------------------------------------

  CALL mp_environ ( numtask_g, g_pos, gs_group )
  CALL mp_environ ( numtask_r, r_dim, r_pos, rs_group )
  IF ( numtask_g /= numtask_r ) THEN
    CALL stop_prg ( "fft3d", "Real space and G space groups are different" )
  END IF
  numtask = numtask_r
  CALL mp_comm_compare ( rs_group, gs_group, iout )
  IF ( iout /= 0 ) THEN
    CALL stop_prg ( "fft3d", "Real space and G space groups are different" )
  END IF

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1._dbl
  END IF

  sign = fsign

  IF ( r_dim ( 2 ) > 1 ) THEN

    !
    ! real space is distributed over x and y coordinate
    ! we have to stages of communication
    !

  ELSE

    !
    ! real space is only distributed over x coordinate
    ! we have one stage of communication, after the transform of
    ! direction x
    !

    IF ( sign == FWFFT ) THEN
      ! gin -> cin

      lg = SIZE ( gin ( :, 1 ) )
      mg = SIZE ( gin ( 1, : ) )

      ALLOCATE ( sbuf ( mg, ng ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "sbuf", mg*ng )

      ! FFT along x
      SELECT CASE ( fft_type )
      CASE ( 1 )
        CALL sg_mltfft ( sign, "T", norm, ng, mg, gin, sbuf )
      CASE ( 2 )
!        CALL w_mltfft ( )
      CASE ( 3 )
!        CALL i_mltfft ( )
      CASE ( 4 )
!        CALL sgi_mltfft ( )
      CASE DEFAULT
        CALL stop_prg ("fft3d","Unknown FFT library")
      END SELECT

      ! Exchange data ( transpose of matrix )
      ! Sort data
      ! FFT along y and z

      DEALLOCATE ( sbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "sbuf" )

    ELSE IF ( sign == BWFFT ) THEN
      ! cin -> gin

      SELECT CASE ( fft_type )
      CASE ( 1 )
        CALL sg_mltfft ( sign, "T", norm, ng, mg, gin, sbuf )
      CASE ( 2 )
!        CALL w_mltfft ( sign, norm, n, gin, sbuf )
      CASE ( 3 )
!        CALL i_mltfft ( sign, norm, n, gin, sbuf )
      CASE ( 4 )
!        CALL sgi_mltfft ( sign, norm, n, gin, sbuf )
      CASE DEFAULT
        CALL stop_prg ("fft3d","Unknown FFT library")
      END SELECT

    ELSE
      CALL stop_prg ( "fft3d", "Illegal fsign parameter" )
    ENDIF
  
  ENDIF
  
  IF ( PRESENT ( status ) ) THEN
    IF ( sign == 0 ) THEN
      status = 1
    ELSE
      status = 0
    END IF
  END IF

END SUBROUTINE fft3d_ps

!******************************************************************************

SUBROUTINE fft3d_pb ( fsign )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: fsign

!------------------------------------------------------------------------------


END SUBROUTINE fft3d_pb

!!*****
!******************************************************************************

END MODULE fft_tools

!******************************************************************************
