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
!!     JGH (30-Nov-2000: ESSL FFT Library added
!!     JGH (5-Jan-2001): Added SGI library FFT
!!     JGH (14-Jan-2001): Added parallel 3d FFT
!!
!!   NOTES
!!     How to add a new FFT library:
!!      - create a new interface library : fftXX_lib with the entries
!!        fft3d, mfft2d, mltfft, and fft_get_lengths, 
!!        see fftsg_lib for a template
!!      - add in this file the entries to the new library; in each
!!        subroutine there will be an additional CASE
!!
!!   SOURCE
!******************************************************************************

MODULE fft_tools

  USE fftessl_lib,                     ONLY: i_fft3d => fft3d,&
                                             i_fft_get_lengths => fft_get_lengths,&
                                             i_mfft2d => mfft2d,&
                                             i_mltfft => mltfft
  USE fftsg_lib,                       ONLY: sg_fft3d => fft3d,&
                                             sg_fft_get_lengths => fft_get_lengths,&
                                             sg_mfft2d => mfft2d,&
                                             sg_mltfft => mltfft
  USE fftsgi_lib,                      ONLY: sgi_fft3d => fft3d,&
                                             sgi_fft_get_lengths => fft_get_lengths,&
                                             sgi_mfft2d => mfft2d,&
                                             sgi_mltfft => mltfft
  USE fftw_lib,                        ONLY: w_fft3d => fft3d,&
                                             w_fft_get_lengths => fft_get_lengths,&
                                             w_mfft2d => mfft2d,&
                                             w_mltfft => mltfft
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_alltoall,&
                                             mp_cart_coords,&
                                             mp_cart_rank,&
                                             mp_cart_sub,&
                                             mp_comm_compare,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_rank_compare,&
                                             mp_sum
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: init_fft, get_fft_library, fft1dm, fft2dm, fft3d
  PUBLIC :: fft_radix_operations
  PUBLIC :: FWFFT, BWFFT
  PUBLIC :: FFT_RADIX_CLOSEST, FFT_RADIX_NEXT, FFT_RADIX_ALLOWED
  PUBLIC :: FFT_RADIX_DISALLOWED, FFT_RADIX_NEXT_ODD
  
  INTEGER, PARAMETER :: FWFFT = +1, BWFFT = -1
  INTEGER, PARAMETER :: FFT_RADIX_CLOSEST = 493, FFT_RADIX_NEXT = 494
  INTEGER, PARAMETER :: FFT_RADIX_ALLOWED = 495, FFT_RADIX_DISALLOWED = 496
  INTEGER, PARAMETER :: FFT_RADIX_NEXT_ODD = 497
  
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
!!     Label 4 : FFT SGI
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_fft ( fftlib, fftnum )

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: fftlib
    INTEGER, INTENT(IN), OPTIONAL            :: fftnum

    COMPLEX(KIND=dp)                             :: zz( 4, 4, 4 )
    INTEGER                                  :: i, n(3) = 4, stat

  IF ( PRESENT ( fftlib ) ) THEN
    SELECT CASE ( fftlib )
    CASE DEFAULT
      CALL stop_program ("init_fft","Unknown FFT library: "//TRIM(fftlib))
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
    zz = 0.1_dp
    CALL fft3d ( 1, n, zz, status=stat )
    IF ( stat /= 0 ) call stop_program ("init_fft","FFT library not available")
  ELSE
    zz = 0.1_dp
    DO i = 1, 100
      fft_type = i
      CALL fft3d ( 1, n, zz, status=stat )
      IF ( stat == 0 ) EXIT
      fft_type = 0
    END DO
    IF (fft_type == 0 ) call stop_program ("init_fft",&
                                           "No FFT library available")
  END IF

  IF ( PRESENT ( fftnum ) ) THEN
    IF ( fft_type /= fftnum ) call stop_program ("init_fft",&
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

    INTEGER, INTENT(OUT)                     :: fft_handle
    CHARACTER(LEN=*), INTENT(OUT)            :: library

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


    INTEGER, INTENT(IN)                      :: radix_in
    INTEGER, INTENT(OUT)                     :: radix_out
    INTEGER, INTENT(IN)                      :: operation

    INTEGER                                  :: i, iloc
    INTEGER, DIMENSION(:), POINTER           :: data

!------------------------------------------------------------------------------

  SELECT CASE ( fft_type )
  CASE ( 1 )
    CALL sg_fft_get_lengths ( data, max_length = 1024 )
  CASE ( 2 )
    CALL w_fft_get_lengths ( data, max_length = 5096 )
  CASE ( 3 )
    CALL i_fft_get_lengths ( data, max_length = 5096 )
  CASE ( 4 )
    CALL sgi_fft_get_lengths ( data, max_length = 1024 )
  CASE DEFAULT
    CALL stop_program ("fft3d","Unknown FFT library")
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
     CALL stop_program ( "fft_radix_operations",&
                         "index to radix array not found" )
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

  ELSE IF ( OPERATION == FFT_RADIX_NEXT_ODD ) THEN
     DO i = iloc, SIZE ( data )
       IF ( MOD ( data ( i ), 2 ) == 1 ) THEN
         radix_out = data ( i )
         EXIT
       END IF
     END DO
     IF ( MOD ( radix_out, 2 ) == 0 ) THEN
        CALL stop_program ( "fft_radix_operations", "no odd radix found" )
     END IF
     
  ELSE
     CALL stop_program ( "fft_radix_operations",&
                         "disallowed radix operation" )
  END IF

  ! seems to be allocated in fft_get_lengths
  DEALLOCATE(data) 
  
END SUBROUTINE fft_radix_operations

!!*****
!******************************************************************************
!!****** fft_tools/fft1dm [1.0] *
!!
!!   NAME
!!     fft1dm
!!
!!   FUNCTION
!!     Calls the multiple 1D-FFT function from the initialized library
!!
!!   AUTHOR
!!     JGH (8-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft1dm ( fsign, tin, tout, n, m, zin, zout, scale )


    INTEGER, INTENT(IN)                      :: fsign
    CHARACTER(LEN=*), INTENT(IN)             :: tin, tout
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: zout
    REAL(KIND=dp), INTENT(IN)                    :: scale

    INTEGER                                  :: sign

!------------------------------------------------------------------------------

  IF ( m == 0 ) RETURN

! CALL timeset("fft1dm","I","Mflops",handle)

  sign = fsign

  SELECT CASE ( fft_type )
  CASE ( 1 )
    CALL sg_mltfft ( tin, tout, sign, scale, n, m, zin, zout )
  CASE ( 2 )
    CALL w_mltfft ( tin, tout, sign, scale, n, m, zin, zout )
  CASE ( 3 )
    CALL i_mltfft ( tin, tout, sign, scale, n, m, zin, zout )
  CASE ( 4 )
    CALL sgi_mltfft ( tin, tout, sign, scale, n, m, zin, zout )
  CASE DEFAULT
    CALL stop_program ("fft1dm","Unknown FFT library")
  END SELECT

! flops = REAL ( m*n,KIND=dp) * 5.e-6_dp * LOG ( REAL ( n,KIND=dp) )
! CALL timestop(flops,handle)

END SUBROUTINE fft1dm

!!*****
!******************************************************************************
!!****** fft_tools/fft2dm [1.0] *
!!
!!   NAME
!!     fft2dm
!!
!!   FUNCTION
!!     Calls the multiple 2D-FFT function from the initialized library
!!
!!   AUTHOR
!!     JGH (9-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft2dm ( fsign, tin, tout, na, nb, m, zin, zout, scale )


    INTEGER, INTENT(IN)                      :: fsign
    CHARACTER(LEN=*), INTENT(IN)             :: tin, tout
    INTEGER, INTENT(IN)                      :: na, nb, m
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zout
    REAL(KIND=dp), INTENT(IN)                    :: scale

    INTEGER                                  :: sign

!------------------------------------------------------------------------------

  IF ( m == 0 ) RETURN

! CALL timeset("fft2dm","I","Mflops",handle)

  sign = fsign

  SELECT CASE ( fft_type )
  CASE ( 1 )
    CALL sg_mfft2d ( tin, tout, sign, scale, na, nb, m, zin, zout )
  CASE ( 2 )
    CALL w_mfft2d ( tin, tout, sign, scale, na, nb, m, zin, zout )
  CASE ( 3 )
    CALL i_mfft2d ( tin, tout, sign, scale, na, nb, m, zin, zout )
  CASE ( 4 )
    CALL sgi_mfft2d ( tin, tout, sign, scale, na, nb, m, zin, zout )
  CASE DEFAULT
    CALL stop_program ("fft2dm","Unknown FFT library")
  END SELECT

! flops = 5.e-6_dp * REAL ( m*na*nb,KIND=dp) * ( LOG ( REAL ( na*nb,KIND=dp) ) )
! CALL timestop(flops,handle)

END SUBROUTINE fft2dm

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

SUBROUTINE fft3d_s ( fsign, n, zin, zout, scale, status, debug )


    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zout
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    INTEGER                                  :: handle, ld(3), lo(3), sign
    LOGICAL                                  :: test
    REAL(KIND=dp)                                :: flops, in_sum, norm, nx, &
                                                out_sum

!------------------------------------------------------------------------------

  CALL timeset("fft3d_s","I","Mflops",handle)

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1.0_dp
  END IF

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

  IF ( test ) THEN
    ld ( 1 ) = SIZE ( zin ,1 )
    ld ( 2 ) = SIZE ( zin ,2 )
    ld ( 3 ) = SIZE ( zin ,3 )
    in_sum = SUM ( ABS ( zin ) )
  ENDIF

  sign = fsign

  SELECT CASE ( fft_type )
  CASE ( 1 )
    CALL sg_fft3d ( sign, norm, n, zin, zout )
  CASE ( 2 )
    CALL w_fft3d ( sign, norm, n, zin, zout )
  CASE ( 3 )
    CALL i_fft3d ( sign, norm, n, zin, zout )
  CASE ( 4 )
    CALL sgi_fft3d ( sign, norm, n, zin, zout )
  CASE DEFAULT
    CALL stop_program ("fft3d", "Unknown FFT library")
  END SELECT

  IF ( PRESENT ( status ) ) THEN
    IF ( sign == 0 ) THEN
      status = 1
    ELSE
      status = 0
    END IF
  END IF

  IF ( test ) THEN
    IF ( PRESENT ( zout ) ) THEN
      out_sum = SUM ( ABS ( zout ) )
      lo ( 1 ) = SIZE ( zout ,1 )
      lo ( 2 ) = SIZE ( zout ,2 )
      lo ( 3 ) = SIZE ( zout ,3 )
      WRITE ( *, '(A)') "  Out of place 3D FFT (local)  : fft3d_s"
      WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
      WRITE ( *, '(A,T60,3I7)') "     Input array dimensions ",ld
      WRITE ( *, '(A,T60,3I7)') "     Output array dimensions ",lo
      WRITE ( *, '(A,T61,E20.14)') "     Sum of input data ",in_sum
      WRITE ( *, '(A,T61,E20.14)') "     Sum of output data ",out_sum
    ELSE
      out_sum = SUM ( ABS ( zin ) )
      WRITE ( *, '(A)') "  In place 3D FFT (local)  : fft3d_s"
      WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
      WRITE ( *, '(A,T60,3I7)') "     Input/output array dimensions ",ld
      WRITE ( *, '(A,T61,E20.14)') "     Sum of input data ",in_sum
      WRITE ( *, '(A,T61,E20.14)') "     Sum of output data ",out_sum
    END IF
  END IF

  nx = 5.e-6_dp * REAL ( n(1)*n(2)*n(3),KIND=dp)
  flops = nx * LOG ( REAL ( n(1)*n(2)*n(3),KIND=dp) ) 
  CALL timestop(flops,handle)

END SUBROUTINE fft3d_s

!******************************************************************************

SUBROUTINE fft3d_ps ( fsign, n, cin, gin, gs_group, rs_group, yzp, nyzray, &
                      bo, scale, status, debug )


    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: cin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: gin
    INTEGER, INTENT(IN)                      :: gs_group, rs_group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nyzray
    INTEGER, DIMENSION(:, :, 0:, :), &
      INTENT(IN)                             :: bo
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuf, qbuf, rbuf, sbuf
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: tbuf
    INTEGER :: g_pos, handle, ierr, iout, lg, lmax, mg, mx1, mx2, my1, mz2, &
      numtask, numtask_g, numtask_r, nx, ny, nz, r_dim(2), r_pos(2), rp, sign
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: p2p
    LOGICAL                                  :: test
    REAL(KIND=dp)                                :: flops, nf, norm, sum_data

!------------------------------------------------------------------------------

  CALL timeset("fft3d_ps","I","Gflops",handle)

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

  CALL mp_environ ( numtask_g, g_pos, gs_group )
  CALL mp_environ ( numtask_r, r_dim, r_pos, rs_group )
  IF ( numtask_g /= numtask_r ) THEN
    CALL stop_program ( "fft3d", "Real space and G space groups are different")
  END IF
  numtask = numtask_r
  CALL mp_comm_compare ( rs_group, gs_group, iout )
  IF ( iout /= 0 ) THEN
    CALL stop_program("fft3d","Real space and G space groups are different")
  END IF

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1.0_dp
  END IF

  sign = fsign

  lg = SIZE ( gin ,1 )
  mg = SIZE ( gin ,2 )

  nx = SIZE ( cin ,1 )
  ny = SIZE ( cin ,2 )
  nz = SIZE ( cin ,3 )

  ALLOCATE ( p2p ( 0 : numtask - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "p2p", numtask )
    
  CALL mp_rank_compare ( gs_group, rs_group, p2p )

  rp = p2p ( g_pos )
  mx1 = bo ( 2, 1, rp, 1 ) - bo ( 1, 1, rp, 1 ) + 1
  my1 = bo ( 2, 2, rp, 1 ) - bo ( 1, 2, rp, 1 ) + 1
  mx2 = bo ( 2, 1, rp, 2 ) - bo ( 1, 1, rp, 2 ) + 1
  mz2 = bo ( 2, 3, rp, 2 ) - bo ( 1, 3, rp, 2 ) + 1

  IF ( test ) THEN
    IF ( g_pos == 0 ) THEN
      WRITE ( *, '(A)') "  Parallel 3D FFT : fft3d_ps"
      WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
      WRITE ( *, '(A,T67,2I7)') "     Array dimensions (gin) ",lg,mg
      WRITE ( *, '(A,T60,3I7)') "     Array dimensions (cin) ",nx,ny,nz
    END IF
  END IF

  IF ( r_dim ( 2 ) > 1 ) THEN

    !
    ! real space is distributed over x and y coordinate
    ! we have two stages of communication
    !

    IF ( r_dim ( 1 ) == 1 ) THEN
      CALL stop_program ( "fft3d", "This processor distribution is not supported" )
    END IF

    IF ( sign == FWFFT ) THEN
      ! cin -> gin

      IF ( test ) THEN
        sum_data = ABS ( SUM ( cin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Two step communication algorithm "
          WRITE ( *, '(A,T60,3I7)') "     Transform Z ",n(3),mx1*my1
          WRITE ( *, '(A,T60,3I7)') "     Transform Y ",n(2),mx2*mz2
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      ALLOCATE ( pbuf ( mx1*my1, n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf", n(3)*mx1*my1 )
      ALLOCATE ( qbuf ( n(3), mx1*my1 ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf", n(3)*mx1*my1 )

      ! FFT along z
      CALL zcopy ( mx1*my1*n(3), cin, 1 , pbuf, 1 )
      CALL fft1dm ( sign, "T", "N", n ( 3 ), mx1*my1, pbuf, qbuf, norm )

      DEALLOCATE ( pbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf" )
      ALLOCATE ( rbuf ( mx2*mz2, n(2) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "rbuf", n(2)*mx2*mz2 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( qbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
        END IF
      END IF

     ! Exchange data ( transpose of matrix )
      CALL cube_transpose_2 ( qbuf, rs_group, bo ( :, :, :, 1 ), &
                            bo ( :, :, :, 2 ), rbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( rbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) T",sum_data
        END IF
      END IF

      DEALLOCATE ( qbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf" )
      ALLOCATE ( pbuf ( n(2), mx2*mz2 ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf", mx2*mz2*n(2) )

      ! FFT along y
      CALL fft1dm ( sign, "T", "N", n ( 2 ), mx2*mz2, rbuf, pbuf, 1.0_dp )

      DEALLOCATE ( rbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "rbuf" )
      ALLOCATE ( qbuf ( nyzray ( g_pos ), n(1) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf", nyzray ( g_pos )*n(1) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( pbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) TS",sum_data
        END IF
      END IF

      ! Exchange data ( transpose of matrix ) and sort
      CALL xz_to_yz ( pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
                      bo ( :, : , 0: , 2 ), qbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( qbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) TS",sum_data
        END IF
      END IF

      DEALLOCATE ( pbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf" )

      ! FFT along x
      CALL fft1dm ( sign, "T", "N", n ( 1 ), nyzray ( g_pos ), qbuf, gin, 1.0_dp )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
        END IF
      END IF

      DEALLOCATE ( qbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf" )

    ELSE IF ( sign == BWFFT ) THEN
      ! gin -> cin

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Two step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
          WRITE ( *, '(A,T60,3I7)') "     Transform Y ",n(2),mx2*mz2
          WRITE ( *, '(A,T60,3I7)') "     Transform Z ",n(3),mx1*my1
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      ALLOCATE ( pbuf ( mg, lg ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf", lg*mg )
      ALLOCATE ( qbuf ( n ( 2 ), mx2*mz2 ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf", n(2)*mx2*mz2 )

      ! FFT along x
      CALL fft1dm ( sign, "N", "T", n ( 1 ), nyzray ( g_pos ), gin, pbuf, norm )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( pbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
        END IF
      END IF

      ! Exchange data ( transpose of matrix ) and sort
      CALL yz_to_xz ( pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
                      bo ( :, : , : , 2 ), qbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( qbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
        END IF
      END IF

      DEALLOCATE ( pbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf" )
      ALLOCATE ( rbuf ( mx2*mz2, n ( 2 ) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "rbuf", n ( 2 ) * mx2*mz2 )

      ! FFT along y
      CALL fft1dm ( sign, "N", "T", n ( 2 ), mx2*mz2, qbuf, rbuf, 1.0_dp )

      DEALLOCATE ( qbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf" )
      ALLOCATE ( pbuf ( n(3), mx1*my1 ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf", n(3)*mx1*my1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( rbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
        END IF
      END IF

      ! Exchange data ( transpose of matrix )
      CALL cube_transpose_1 ( rbuf, rs_group, bo ( :, :, :, 2 ), &
                            bo ( :, :, :, 1 ), pbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( pbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) T",sum_data
        END IF
      END IF

      DEALLOCATE ( rbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "rbuf" )
      ALLOCATE ( qbuf ( mx1*my1, n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf", mx1*my1*n(3) )

      ! FFT along z
      CALL fft1dm ( sign, "N", "T", n ( 3 ), mx1*my1, pbuf, qbuf, 1.0_dp )
      CALL zcopy ( mx1*my1*n(3), qbuf, 1 , cin, 1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( cin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
        END IF
      END IF

      DEALLOCATE ( qbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "qbuf" )
      DEALLOCATE ( pbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "pbuf" )

    ELSE

      CALL stop_program ( "fft3d", "Illegal fsign parameter" )

    ENDIF

  ELSE

    !
    ! real space is only distributed over x coordinate
    ! we have one stage of communication, after the transform of
    ! direction x
    !

    IF ( mg == 0 ) THEN
      lmax = lg
    ELSE
      lmax = MAX ( lg, (nx*ny*nz)/mg + 1 )
    END IF
    ALLOCATE ( sbuf ( mg, lmax ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "sbuf", mg*lmax )

    ALLOCATE ( tbuf ( ny, nz, nx ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "tbuf", nx*ny*nz )

    IF ( sign == FWFFT ) THEN
      ! cin -> gin

      IF ( test ) THEN
        sum_data = ABS ( SUM ( cin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A)') "     One step communication algorithm "
          WRITE ( *, '(A,T60,3I7)') "     Transform YZ ",n(2),n(3),nx
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      ! FFT along y and z
      CALL fft2dm ( sign, "T", "N", n ( 2 ), n ( 3 ), nx, cin, tbuf, 1.0_dp )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( tbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
        END IF
      END IF

      ! Exchange data ( transpose of matrix ) and sort
      CALL yz_to_x ( tbuf, gs_group, g_pos, p2p, yzp, nyzray, &
                     bo ( :, :, :, 2 ), sbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( sbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
        END IF
      END IF

      ! FFT along x
      CALL fft1dm ( sign, "T", "N", n ( 1 ), nyzray ( g_pos ), sbuf, gin, norm )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
        END IF
      END IF

    ELSE IF ( sign == BWFFT ) THEN
      ! gin -> cin

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A)') "  One step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),nyzray(g_pos)
          WRITE ( *, '(A,T60,3I7)') "     Transform YZ ",n(2),n(3),nx
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      ! FFT along x
      CALL fft1dm ( sign, "N", "T", n ( 1 ), nyzray ( g_pos ), gin, sbuf, norm )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( sbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) TS",sum_data
        END IF
      END IF

      ! Exchange data ( transpose of matrix ) and sort
      CALL x_to_yz ( sbuf, gs_group, g_pos, p2p, yzp, nyzray, &
                     bo ( :, :, :, 2 ), tbuf )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( tbuf ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) TS",sum_data
        END IF
      END IF

      ! FFT along y and z
      CALL fft2dm ( sign, "N", "T", n ( 2 ), n ( 3 ), nx, tbuf, cin, 1.0_dp )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( cin ) )
        CALL mp_sum ( sum_data, gs_group )
        IF ( g_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
        END IF
      END IF

    ELSE
      CALL stop_program ( "fft3d", "Illegal fsign parameter" )
    ENDIF

    DEALLOCATE ( tbuf, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "tbuf" ) 
    DEALLOCATE ( sbuf, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "sbuf" )
  
  ENDIF

  DEALLOCATE ( p2p, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "p2p" )

  IF ( PRESENT ( status ) ) THEN
    IF ( sign == 0 ) THEN
      status = 1
    ELSE
      status = 0
    END IF
  END IF


  nf = 5.e-6_dp * REAL ( n(1)*n(2)*n(3),KIND=dp)
  flops = nf * LOG ( REAL ( n(1)*n(2)*n(3),KIND=dp) )
  flops = 1.0e-3_dp*flops/numtask_r  ! per cpu
  CALL timestop(flops,handle)

END SUBROUTINE fft3d_ps

!******************************************************************************

SUBROUTINE fft3d_pb ( fsign, n, zin, gin, group, bo, scale, status, debug )


    INTEGER, INTENT(IN)                      :: fsign
    INTEGER, DIMENSION(3), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: gin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:, :), &
      INTENT(IN)                             :: bo
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: scale
    INTEGER, INTENT(OUT), OPTIONAL           :: status
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: abuf, bbuf
    INTEGER                                  :: handle, ierr, lg(2), lz(3), &
                                                mx, my, my_pos, mz, np, sign
    INTEGER, DIMENSION(2)                    :: dim, pos
    LOGICAL                                  :: test
    REAL(KIND=dp)                                :: flops, norm, nx, sum_data

!------------------------------------------------------------------------------
! "Real Space"  1) xyZ      or      1) xYZ
!               2) xYz      or         not used
! "G Space"     3) Xyz      or      3) XYz
!
! There is one communicator (2-dimensional) for all distributions
! np = n1 * n2, where np is the total number of processors
! If n2 = 1, we have the second case and only one transpose step is needed
! 
! Assignment of dimensions to axis for different steps
! First case: 1) n1=x; n2=y
!             2) n1=x; n2=z
!             3) n1=y; n2=z
! Second case 1) n1=x
!             3) n1=z
!
! The more general case with two communicators for the initial and final
! distribution is not covered.
!------------------------------------------------------------------------------

  CALL timeset("fft3d_pb","I","Mflops",handle)

  CALL mp_environ ( np, dim, pos, group )
  CALL mp_cart_rank ( group, pos, my_pos )

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

  IF ( PRESENT ( scale ) ) THEN
    norm = scale
  ELSE
    norm = 1.0_dp
  END IF

  sign = fsign

  IF ( test ) THEN
    lg ( 1 ) = SIZE ( gin ,1 )
    lg ( 2 ) = SIZE ( gin ,2 )
    lz ( 1 ) = SIZE ( zin ,1 )
    lz ( 2 ) = SIZE ( zin ,2 )
    lz ( 3 ) = SIZE ( zin ,3 )
    IF ( my_pos == 0 ) THEN
      WRITE ( *, '(A)') "  Parallel 3D FFT : fft3d_pb"
      WRITE ( *, '(A,T60,3I7)') "     Transform lengths ",n
      WRITE ( *, '(A,T67,2I7)') "     Array dimensions (gin) ",lg
      WRITE ( *, '(A,T60,3I7)') "     Array dimensions (cin) ",lz
    END IF
  END IF

  IF ( dim ( 1 ) > 1 .AND. dim ( 2 ) > 1 ) THEN

    !
    ! First case; two stages of communication
    !

    IF ( sign == FWFFT ) THEN
      ! Stage 1 -> 3

      mx = bo ( 2, 1, my_pos, 1 ) - bo ( 1, 1, my_pos, 1 ) + 1
      my = bo ( 2, 2, my_pos, 1 ) - bo ( 1, 2, my_pos, 1 ) + 1
      ALLOCATE ( abuf ( mx*my, n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*my*n(3) )
      ALLOCATE ( bbuf ( n(3), mx*my ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(3)*mx*my )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( zin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Two step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx*my
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      CALL zcopy ( n(3)*mx*my, zin, 1, abuf, 1 )

      ! FFT along z
      CALL fft1dm ( sign, "T", "N", n ( 3 ), mx*my, abuf, bbuf, norm )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      mx = bo ( 2, 1, my_pos, 2 ) - bo ( 1, 1, my_pos, 2 ) + 1
      mz = bo ( 2, 3, my_pos, 2 ) - bo ( 1, 3, my_pos, 2 ) + 1
      ALLOCATE ( abuf ( mx*mz, n(2) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*mz*n(2) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
        END IF
      END IF

      CALL cube_transpose_2 ( bbuf, group, bo(:,:,:,1), bo(:,:,:,2), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( n(2), mx*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*mz*n(2) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
        END IF
      END IF

      ! FFT along y
      CALL fft1dm ( sign, "T", "N", n ( 2 ), mx*mz, abuf, bbuf, 1.0_dp )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      my = bo ( 2, 2, my_pos, 3 ) - bo ( 1, 2, my_pos, 3 ) + 1
      mz = bo ( 2, 3, my_pos, 3 ) - bo ( 1, 3, my_pos, 3 ) + 1
      ALLOCATE ( abuf ( my*mz, n(1) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", my*mz*n(1) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
        END IF
      END IF

      CALL cube_transpose_4 ( bbuf, group, bo(:,:,:,2), bo(:,:,:,3), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( n(1), my*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*my*n(1) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) ",sum_data
        END IF
      END IF

      ! FFT along x
      CALL fft1dm ( sign, "T", "N", n ( 1 ), my*mz, abuf, bbuf, 1.0_dp )

      CALL zcopy ( n ( 1 ) * my * mz, bbuf, 1, gin, 1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
        END IF
      END IF

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )

    ELSEIF ( sign == BWFFT ) THEN
      ! Stage 3 -> 1

      my = bo ( 2, 2, my_pos, 3 ) - bo ( 1, 2, my_pos, 3 ) + 1
      mz = bo ( 2, 3, my_pos, 3 ) - bo ( 1, 3, my_pos, 3 ) + 1
      ALLOCATE ( abuf ( n(1), my*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(1)*my*mz )
      ALLOCATE ( bbuf ( my*mz, n(1) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(1)*my*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A)') "  Two step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      CALL zcopy ( n(1)*my*mz, gin, 1, abuf, 1 )

      ! FFT along x
      CALL fft1dm ( sign, "N", "T", n ( 1 ), my*mz, abuf, bbuf, 1.0_dp )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      mx = bo ( 2, 1, my_pos, 2 ) - bo ( 1, 1, my_pos, 2 ) + 1
      mz = bo ( 2, 3, my_pos, 2 ) - bo ( 1, 3, my_pos, 2 ) + 1
      ALLOCATE ( abuf ( n(2), mx*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(2)*mx*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
        END IF
      END IF

      CALL cube_transpose_3 ( bbuf, group, bo(:,:,:,3), bo(:,:,:,2), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( mx*mz, n(2) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(2)*mx*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
        END IF
      END IF

      ! FFT along y
      CALL fft1dm ( sign, "N", "T", n ( 2 ), mx*mz, abuf, bbuf, 1.0_dp )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      mx = bo ( 2, 1, my_pos, 1 ) - bo ( 1, 1, my_pos, 1 ) + 1
      my = bo ( 2, 2, my_pos, 1 ) - bo ( 1, 2, my_pos, 1 ) + 1
      ALLOCATE ( abuf ( n(3), mx*my ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(3)*mx*my )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) T",sum_data
        END IF
      END IF

      CALL cube_transpose_1 ( bbuf, group, bo(:,:,:,2), bo(:,:,:,1), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( mx*my, n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(3)*mx*my )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx*my
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(5) ",sum_data
        END IF
      END IF

      ! FFT along z
      CALL fft1dm ( sign, "N", "T", n ( 3 ), mx*my, abuf, bbuf, norm )

      CALL zcopy ( n(3)*mx*my, bbuf, 1, zin, 1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( zin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(6) ",sum_data
        END IF
      END IF

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )

    ELSE
      CALL stop_program ( "fft3d", "Illegal fsign parameter" )
    ENDIF

  ELSEIF ( dim ( 2 ) == 1 ) THEN

    !
    ! Second case; one stage of communication
    !

    IF ( sign == FWFFT ) THEN
      ! Stage 1 -> 3

      mx = bo ( 2, 1, my_pos, 1 ) - bo ( 1, 1, my_pos, 1 ) + 1
      my = bo ( 2, 2, my_pos, 1 ) - bo ( 1, 2, my_pos, 1 ) + 1
      mz = bo ( 2, 3, my_pos, 1 ) - bo ( 1, 3, my_pos, 1 ) + 1

      ALLOCATE ( abuf ( mx*my, n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*my*n(3) )
      ALLOCATE ( bbuf ( n(3), mx*my ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(3)*mx*my )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( zin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A)') "  one step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx*my
          WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      CALL zcopy ( n(3)*mx*my, zin, 1, abuf, 1 )

      ! FFT along z
      CALL fft1dm ( sign, "T", "N", n ( 3 ), mx*my, abuf, bbuf, norm )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      ALLOCATE ( abuf ( mx*mz, n(2) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", mx*mz*n(2) )

      call zcopy ( n ( 3 ) * mx*my, bbuf,1,abuf,1)

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( n(2), mx*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", mx*mz*n(2) )

      ! FFT along y
      CALL fft1dm ( sign, "T", "N", n ( 2 ), mx*mz, abuf, bbuf, 1.0_dp )

      mx = bo ( 2, 1, my_pos, 3 ) - bo ( 1, 1, my_pos, 3 ) + 1
      my = bo ( 2, 2, my_pos, 3 ) - bo ( 1, 2, my_pos, 3 ) + 1
      mz = bo ( 2, 3, my_pos, 3 ) - bo ( 1, 3, my_pos, 3 ) + 1

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( abuf ( my*mz, n ( 1 ) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", my*mz*n(1) )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
        END IF
      END IF

      CALL cube_transpose_6 ( bbuf, group, bo(:,:,:,1), bo(:,:,:,3), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( n ( 1 ), my*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(1)*my*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
        END IF
      END IF

      ! FFT along x
      CALL fft1dm ( sign, "T", "N", n ( 1 ), my*mz, abuf, bbuf, 1.0_dp )

      CALL zcopy ( n ( 1 ) * my * mz, bbuf, 1, gin, 1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
        END IF
      END IF

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )

    ELSEIF ( sign == BWFFT ) THEN
      ! Stage 3 -> 1

      mx = bo ( 2, 1, my_pos, 3 ) - bo ( 1, 1, my_pos, 3 ) + 1
      my = bo ( 2, 2, my_pos, 3 ) - bo ( 1, 2, my_pos, 3 ) + 1
      mz = bo ( 2, 3, my_pos, 3 ) - bo ( 1, 3, my_pos, 3 ) + 1
      ALLOCATE ( abuf ( n(1), my*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(1)*my*mz )
      ALLOCATE ( bbuf ( my*mz, n(1) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(1)*my*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( gin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A)') "  one step communication algorithm "
          WRITE ( *, '(A,T67,2I7)') "     Transform X ",n(1),my*mz
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(1) ",sum_data
        END IF
      END IF

      CALL zcopy ( n(1)*my*mz, gin, 1, abuf, 1 )

      ! FFT along x
      CALL fft1dm ( sign, "N", "T", n ( 1 ), my*mz, abuf, bbuf, 1.0_dp )

      mx = bo ( 2, 1, my_pos, 1 ) - bo ( 1, 1, my_pos, 1 ) + 1
      my = bo ( 2, 2, my_pos, 1 ) - bo ( 1, 2, my_pos, 1 ) + 1
      mz = bo ( 2, 3, my_pos, 1 ) - bo ( 1, 3, my_pos, 1 ) + 1

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      ALLOCATE ( abuf ( n(2), mx*mz ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(2)*mx*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( bbuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(2) T",sum_data
        END IF
      END IF

      CALL cube_transpose_5 ( bbuf, group, bo(:,:,:,3), bo(:,:,:,1), abuf )

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( mz*mx, n(2) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(2)*mx*mz )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( abuf ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T67,2I7)') "     Transform Y ",n(2),mx*mz
          WRITE ( *, '(A,T67,2I7)') "     Transform Z ",n(3),mx*my
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(3) ",sum_data
        END IF
      END IF

      ! FFT along y
      CALL fft1dm ( sign, "N", "T", n ( 2 ), mx*mz, abuf, bbuf, 1.0_dp )

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      ALLOCATE ( abuf ( n(3), mx*my ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf", n(3)*mx*my )

      CALL zcopy(n(3)*mx*my,bbuf,1,abuf,1)

      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )
      ALLOCATE ( bbuf ( mx*my,n(3) ), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf", n(3)*mx*my )

      ! FFT along z
      CALL fft1dm ( sign, "N", "T", n ( 3 ), mx*my, abuf, bbuf, norm )

      CALL zcopy ( n(3)*mx*my, bbuf, 1, zin, 1 )

      IF ( test ) THEN
        sum_data = ABS ( SUM ( zin ) )
        CALL mp_sum ( sum_data, group )
        IF ( my_pos == 0 ) THEN
          WRITE ( *, '(A,T61,E20.14)') "     Sum of data(4) ",sum_data
        END IF
      END IF

      DEALLOCATE ( abuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "abuf" )
      DEALLOCATE ( bbuf, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "fft3d", "bbuf" )

    ELSE
      CALL stop_program ( "fft3d", "Illegal fsign parameter" )
    ENDIF

  ELSE
     
    CALL stop_program ( "fft3d", "This partition not implemented" )

  END IF

  IF ( PRESENT ( status ) ) THEN
    IF ( sign == 0 ) THEN
      status = 1
    ELSE
      status = 0
    END IF
  END IF

  nx = 5.e-6_dp * REAL ( n(1)*n(2)*n(3),KIND=dp)
  flops = nx * LOG ( REAL ( n(1)*n(2)*n(3),KIND=dp) ) 
  CALL timestop(flops,handle)

END SUBROUTINE fft3d_pb

!!*****
!******************************************************************************
!!****** fft_tools/x_to_yz [1.0] *
!!
!!   NAME
!!     x_to_yz
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (14-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE x_to_yz ( sb, group, my_pos, p2p, yzp, nray, bo, tb )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group, my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: tb

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rr
    INTEGER                                  :: handle, ierr, ip, ir, ix, &
                                                ixx, iy, iz, mpr, nm, np, nr, &
                                                nx
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
    REAL(KIND=dp)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset("x_to_yz","I","Mcopy",handle)

  np = SIZE ( p2p )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "x_to_yz", "scount", 4*np )

  mpr = p2p ( my_pos )
  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  nm = MAXVAL ( nray ( 0 : np - 1 ) )
  ALLOCATE ( rr ( nm * nx, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "x_to_yz", "rr", np*nm*nx )
  nr = SIZE ( rr )
  CALL dcopy ( 2*nr, 0.0_dp, 0, rr, 1 )
  
  nr = nray ( my_pos )
  DO ip = 0, np - 1
    ix = p2p ( ip )
    nx = bo ( 2, 1, ix ) - bo ( 1, 1, ix ) + 1
    scount ( ip ) = nr * nx
    sdispl ( ip ) = nr * ( bo ( 1, 1, ix ) - 1 ) 
  END DO
  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  DO ip = 0, np - 1
    nr = nray ( ip )
    rcount ( ip ) = nr * nx
    rdispl ( ip ) = nm * nx * ip 
  END DO
  flop = REAL ( nm * nx * np,KIND=dp)

  CALL mp_alltoall ( sb, scount, sdispl, rr, rcount, rdispl, group )

  nx = SIZE ( tb )
  CALL dcopy ( 2*nx, 0.0_dp, 0, tb, 1 )
  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  DO ip = 0, np - 1
    nr = nray ( ip )
    DO ix = 1, nx
      ixx = nr * ( ix - 1 )
      DO ir = 1, nray ( ip )
        iy = yzp ( 1, ir, ip )
        iz = yzp ( 2, ir, ip )
        tb ( iy, iz, ix ) = rr ( ir + ixx, ip )
      END DO
    END DO
  END DO

  DEALLOCATE ( rr, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "x_to_yz", "rr" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "x_to_yz", "scount" )

  flop = flop * 1.e-6_dp
  CALL timestop(flop,handle)

END SUBROUTINE x_to_yz

!!*****
!******************************************************************************
!!****** fft_tools/yz_to_x [1.0] *
!!
!!   NAME
!!     yz_to_x
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (14-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE yz_to_x ( tb, group, my_pos, p2p, yzp, nray, bo, sb )


    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: tb
    INTEGER, INTENT(IN)                      :: group, my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sb

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rr
    INTEGER                                  :: handle, ierr, ip, ir, ix, &
                                                ixx, iy, iz, mpr, nm, np, nr, &
                                                nx
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
    REAL(KIND=dp)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset("yz_to_x","I","Mcopy",handle)

  np = SIZE ( p2p )
  mpr = p2p ( my_pos )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_x", "scount", 4*np )

  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  nm = MAXVAL ( nray ( 0 : np - 1 ) )
  ALLOCATE ( rr ( nm * nx, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_x", "rr", np*nm*nx )
  nx = SIZE ( rr )
  CALL dcopy ( 2*nx, 0.0_dp, 0, rr, 1 )
  
  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  DO ip = 0, np - 1
    nr = nray ( ip )
    DO ix = 1, nx
      ixx = nr * ( ix - 1 )
      DO ir = 1, nray ( ip )
        iy = yzp ( 1, ir, ip )
        iz = yzp ( 2, ir, ip )
        rr ( ir + ixx, ip ) = tb ( iy, iz, ix )
      END DO
    END DO
  END DO

  nr = nray ( my_pos )
  DO ip = 0, np - 1
    ix = p2p ( ip )
    nx = bo ( 2, 1, ix ) - bo ( 1, 1, ix ) + 1
    rcount ( ip ) = nr * nx
    rdispl ( ip ) = nr * ( bo ( 1, 1, ix ) - 1 )
  END DO
  nx = bo ( 2, 1, mpr ) - bo ( 1, 1, mpr ) + 1
  DO ip = 0, np - 1
    nr = nray ( ip )
    scount ( ip ) = nr * nx
    sdispl ( ip ) = nm * nx * ip
  END DO
  flop = REAL ( nm * nx * np,KIND=dp)

  CALL mp_alltoall ( rr, scount, sdispl, sb, rcount, rdispl, group )

  DEALLOCATE ( rr, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_x", "rr" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_x", "scount" )

  flop = flop * 1.e-6
  CALL timestop(flop,handle)

END SUBROUTINE yz_to_x

!!*****
!******************************************************************************
!!****** fft_tools/yz_to_xz [1.0] *
!!
!!   NAME
!!     yz_to_xz
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (18-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE yz_to_xz ( sb, group, dims, my_pos, p2p, yzp, nray, bo, tb )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(2), INTENT(IN)        :: dims
    INTEGER, INTENT(IN)                      :: my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: tb

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rq, rs
    INTEGER                                  :: cpos( 2 ), handle, ierr, ip, &
                                                ipr, ir, ix, iz, jj, jx, jy, &
                                                jz, myx, myz, nmax, np, npx, &
                                                npz, nx, nz, rs_pos
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, &
                                                sdispl, xcor, zcor
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: pgrid
    REAL(KIND=dp)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset("yz_to_xz","I","Mcopy",handle)

  np = SIZE ( p2p )
  nx = MAXVAL ( bo ( 2, 1, : ) - bo ( 1, 1, : ) + 1 )
  nmax = nx * MAXVAL ( nray ( : ) )
  rs_pos = p2p ( my_pos )

  ALLOCATE ( rs ( nmax, 0:np-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rs", nmax*np )
  rs = 0.0_dp
  ALLOCATE ( rq ( nmax, 0:np-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rq", nmax*np )
  npx = dims ( 1 )
  npz = dims ( 2 )
  ALLOCATE ( pgrid ( 0:npx-1, 0:npz-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "pgrid", npx*npz )
  nx = MAXVAL ( bo ( 2, 1, : ) )
  ALLOCATE ( xcor ( 1:nx ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "xcor", nx )
  nz = MAXVAL ( bo ( 2, 3, : ) )
  ALLOCATE ( zcor ( 1:nz ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "zcor", nz )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "scount", 4*np )


  DO ix = 0, npx - 1
    DO iz = 0, npz - 1
      cpos ( 1 ) = ix
      cpos ( 2 ) = iz
      CALL mp_cart_rank ( group, cpos, pgrid ( ix, iz ) )
    END DO
  END DO
  DO ix = 0, npx - 1
    ip = pgrid ( ix, 0 )
    xcor ( bo ( 1, 1, ip ) : bo ( 2, 1, ip ) ) = ix
  END DO
  DO iz = 0, npz - 1
    ip = pgrid ( 0, iz )
    zcor ( bo ( 1, 3, ip ) : bo ( 2, 3, ip ) ) = iz
  END DO

  scount = 0
  nx = SIZE ( sb ,2 )
  DO jx = 1, nx
    DO ir = 1, nray ( my_pos )
      jy = yzp ( 1, ir, my_pos )
      jz = yzp ( 2, ir, my_pos )
      ip = pgrid ( xcor ( jx ), zcor ( jz ) )
      scount ( ip ) = scount ( ip ) + 1
      rs ( scount ( ip ), ip ) = sb ( ir, jx )
    END DO
  END DO

  DO ip = 0, np - 1
    sdispl ( ip ) = nmax * ip
    rdispl ( ip ) = nmax * ip
  END DO
  flop = REAL ( nmax * np,KIND=dp)

  CALL mp_alltoall ( scount, rcount, 1, group )
  CALL mp_alltoall ( rs, scount, sdispl, rq, rcount, rdispl, group )

  CALL mp_cart_coords ( group, rs_pos, cpos )
  myx = cpos ( 1 )
  myz = cpos ( 2 )
  nz = bo ( 2, 3, rs_pos ) - bo ( 1, 3, rs_pos ) + 1
  tb = 0.0_dp
  DO ip = 0, np - 1
    ipr = p2p ( ip )
    jj = 0
    DO jx = 0, bo ( 2, 1, rs_pos ) - bo ( 1, 1, rs_pos )
      DO ir = 1, nray ( ip )
        jz = yzp ( 2, ir, ip )
        if ( zcor ( jz ) == myz ) THEN
          jj = jj + 1
          jy = yzp ( 1, ir, ip )
          jz = jz - bo ( 1, 3, rs_pos ) + 1
          tb ( jy, jz + jx * nz ) = rq ( jj, ipr )
        END IF
      END DO
    END DO
  END DO

  DEALLOCATE ( rs, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rs" )
  DEALLOCATE ( rq, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rq" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "pgrid" )
  DEALLOCATE ( xcor, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "xcor" )
  DEALLOCATE ( zcor, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "zcor" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "scount" )

  flop = flop * 1.e-6
  CALL timestop(flop,handle)

END SUBROUTINE yz_to_xz

!!*****
!******************************************************************************
!!****** fft_tools/xz_to_yz [1.0] *
!!
!!   NAME
!!     xz_to_yz
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (19-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE xz_to_yz ( sb, group, dims, my_pos, p2p, yzp, nray, bo, tb )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(2), INTENT(IN)        :: dims
    INTEGER, INTENT(IN)                      :: my_pos
    INTEGER, DIMENSION(0:), INTENT(IN)       :: p2p
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: yzp
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nray
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: bo
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: tb

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rq, rs
    INTEGER                                  :: cpos( 2 ), handle, ierr, ip, &
                                                ipl, ir, ix, ixx, iz, jj, jx, &
                                                jy, jz, mp, myx, myz, nmax, &
                                                np, npx, npz, nx, ny, nz
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, &
                                                sdispl, xcor, zcor
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: pgrid
    REAL(KIND=dp)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset("xz_to_yz","I","Mcopy",handle)

  np = SIZE ( p2p )
  ny = MAXVAL ( bo ( 2, 2, : ) - bo ( 1, 2, : ) + 1 )
  nmax = (2*ny)/np * SIZE ( sb ,2 )

  ALLOCATE ( rs ( nmax, 0:np-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "rs", nmax*np )
  rs = 0.0_dp
  ALLOCATE ( rq ( nmax, 0:np-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "rq", nmax*np )
  rq = 0.0_dp
  npx = dims ( 1 )
  npz = dims ( 2 )
  ALLOCATE ( pgrid ( 0:npx-1, 0:npz-1 ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "pgrid", npx*npz )
  nx = MAXVAL ( bo ( 2, 1, : ) )
  ALLOCATE ( xcor ( 1:nx ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "xcor", nx )
  nz = MAXVAL ( bo ( 2, 3, : ) )
  ALLOCATE ( zcor ( 1:nz ) , STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "zcor", nz )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "xz_to_yz", "scount", 4*np )

  DO ix = 0, npx - 1
    DO iz = 0, npz - 1
      cpos ( 1 ) = ix
      cpos ( 2 ) = iz
      CALL mp_cart_rank ( group, cpos, pgrid ( ix, iz ) )
    END DO
  END DO
  DO ix = 0, npx - 1
    ip = pgrid ( ix, 0 )
    xcor ( bo ( 1, 1, ip ) : bo ( 2, 1, ip ) ) = ix
  END DO
  DO iz = 0, npz - 1
    ip = pgrid ( 0, iz )
    zcor ( bo ( 1, 3, ip ) : bo ( 2, 3, ip ) ) = iz
  END DO

  rcount = 0
  nx = MAXVAL ( bo ( 2, 1, : ) )
  DO jx = 1, nx
    DO ir = 1, nray ( my_pos )
      jy = yzp ( 1, ir, my_pos )
      jz = yzp ( 2, ir, my_pos )
      ip = pgrid ( xcor ( jx ), zcor ( jz ) )
      rcount ( ip ) = rcount ( ip ) + 1
    END DO
  END DO

  CALL mp_alltoall ( rcount, scount, 1, group )

  CALL mp_cart_coords ( group, p2p ( my_pos ), cpos )
  myx = cpos ( 1 )
  myz = cpos ( 2 )
  mp = p2p ( my_pos )
  nz = bo ( 2, 3, mp ) - bo ( 1, 3, mp ) + 1
  nx = bo ( 2, 1, mp ) - bo ( 1, 1, mp ) + 1
  DO ip = 0, np - 1
    jj = 0
    ipl = p2p ( ip )
    DO ir = 1, nray ( ip )
      jz = yzp ( 2, ir, ip )
      if ( zcor ( jz ) == myz ) THEN
        jj = jj + 1
        jy = yzp ( 1, ir, ip )
        jz = yzp ( 2, ir, ip ) - bo ( 1, 3, mp ) + 1
        DO jx = 0, nx - 1
          ixx = jj + jx * scount ( ipl )/nx
          rq ( ixx, ipl ) = sb ( jy, jz + jx * nz ) 
        END DO
      END IF
    END DO
  END DO

  DO ip = 0, np - 1
    sdispl ( ip ) = nmax * ip
    rdispl ( ip ) = nmax * ip
  END DO
  flop = REAL ( nmax * np,KIND=dp)

  CALL mp_alltoall ( rq, scount, sdispl, rs, rcount, rdispl, group )

  rcount = 0
  nx = MAXVAL ( bo ( 2, 1, : ) )
  DO jx = 1, nx
    DO ir = 1, nray ( my_pos )
      jy = yzp ( 1, ir, my_pos )
      jz = yzp ( 2, ir, my_pos )
      ip = pgrid ( xcor ( jx ), zcor ( jz ) )
      rcount ( ip ) = rcount ( ip ) + 1
      tb ( ir, jx ) = rs ( rcount ( ip ), ip )
    END DO
  END DO

  DEALLOCATE ( rs, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rs" )
  DEALLOCATE ( rq, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "rq" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "pgrid" )
  DEALLOCATE ( xcor, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "xcor" )
  DEALLOCATE ( zcor, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "zcor" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "yz_to_xz", "scount" )

  flop = flop * 1.e-6
  CALL timestop(flop,handle)

END SUBROUTINE xz_to_yz

!!*****
!******************************************************************************
!!****** fft_tools/cube_transpose [1.0] *
!!
!!   NAME
!!     cube_transpose
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (20-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE cube_transpose_1 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ipl, ir, &
                                                is, ixy, iz, mip, mz, n, np, &
                                                nx, ny, nz, sub_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pgrid, rcount, rdispl, &
                                                scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos, pp
    LOGICAL, DIMENSION(2)                    :: rdim

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_1","I","",handle)

  rdim ( 1 ) = .FALSE.
  rdim ( 2 ) = .TRUE.
  CALL mp_environ ( n, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )
  CALL mp_cart_sub ( group, rdim, sub_group )


  np = dim ( 2 )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )
  ALLOCATE ( pgrid ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid", np )

  pp = pos
  DO ip = 0, np - 1
    pp ( 2 ) = ip
    CALL mp_cart_rank ( group, pp, pgrid ( ip ) )
  END DO

  nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
  nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
    scount ( ip ) = nx * nz * ny
    sdispl ( ip ) = nx * nz * ( boout ( 1, 2, ipl ) - 1 )
  END DO
  ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
  mz = MAXVAL ( boin ( 2, 3, : ) - boin ( 1, 3, : ) + 1 )
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    nz = boin ( 2, 3, ipl ) - boin ( 1, 3, ipl ) + 1
    rcount ( ip ) = nx * nz * ny
    rdispl ( ip ) = nx * ny * mz * ip
  END DO

  ALLOCATE ( rbuf ( mz*nx*ny, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", mz*nx*ny*np )
  rbuf = 0.0_dp

  CALL mp_alltoall ( sin, scount, sdispl, rbuf, rcount, rdispl, sub_group )

  DO ixy = 1, nx * ny
    DO ip = 0, np - 1
      ipl = pgrid ( ip )
      nz = boin ( 2, 3, ipl ) - boin ( 1, 3, ipl ) + 1
      DO iz = 1, nz
        is = boin ( 1, 3, ipl ) + iz - 1
        ir = iz + nz * ( ixy - 1 )
        sout ( is, ixy ) = rbuf ( ir, ip )
      END DO
    END DO
  END DO

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid" )

  CALL mp_comm_free ( sub_group )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_1

!******************************************************************************

SUBROUTINE cube_transpose_2 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ipl, ir, &
                                                ixy, iz, mip, mz, n, np, nx, &
                                                ny, nz, sub_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pgrid, rcount, rdispl, &
                                                scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos, pp
    LOGICAL, DIMENSION(2)                    :: rdim

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_2","I","",handle)

  rdim ( 1 ) = .FALSE.
  rdim ( 2 ) = .TRUE.
  CALL mp_environ ( n, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )
  CALL mp_cart_sub ( group, rdim, sub_group )

  np = dim ( 2 )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )
  ALLOCATE ( pgrid ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid", np )

  pp = pos
  DO ip = 0, np - 1
    pp ( 2 ) = ip
    CALL mp_cart_rank ( group, pp, pgrid ( ip ) )
  END DO

  nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
  ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
  mz = MAXVAL ( boin ( 2, 3, : ) - boin ( 1, 3, : ) + 1 )

  ALLOCATE ( rbuf ( mz*nx*ny, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", mz*nx*ny*np )
  rbuf = 0.0_dp

  DO ixy = 1, nx * ny
    DO ip = 0, np - 1
      ipl = pgrid ( ip )
      nz = boout ( 2, 3, ipl ) - boout ( 1, 3, ipl ) + 1
      DO iz = boout ( 1, 3, ipl ), boout ( 2, 3, ipl )
        ir = iz - boout ( 1, 3, ipl ) + 1 + ( ixy - 1 ) * nz
        rbuf ( ir, ip ) = sin ( iz, ixy )
      END DO
    END DO
  END DO

  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    nz = boout ( 2, 3, ipl ) - boout ( 1, 3, ipl ) + 1
    scount ( ip ) = nx * ny * nz
    sdispl ( ip ) = nx * ny * mz * ip
  END DO
  nz = boout ( 2, 3, mip ) - boout ( 1, 3, mip ) + 1
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
    rcount ( ip ) = nx * ny * nz
    rdispl ( ip ) = nx * nz * ( boin ( 1, 2, ipl ) - 1 )
  END DO

  sout = 0.0_dp
  CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, sub_group )

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid" )

  CALL mp_comm_free ( sub_group )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_2

!******************************************************************************

SUBROUTINE cube_transpose_3 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ipl, ir, &
                                                is, ixz, iy, mip, my, n, np, &
                                                nx, ny, nz, sub_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pgrid, rcount, rdispl, &
                                                scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos, pp
    LOGICAL, DIMENSION(2)                    :: rdim

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_3","I","",handle)

  rdim ( 1 ) = .TRUE. 
  rdim ( 2 ) = .FALSE.
  CALL mp_environ ( n, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )
  CALL mp_cart_sub ( group, rdim, sub_group )

  np = dim ( 1 )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )
  ALLOCATE ( pgrid ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid", np )

  pp = pos
  DO ip = 0, np - 1
    pp ( 1 ) = ip
    CALL mp_cart_rank ( group, pp, pgrid ( ip ) )
  END DO

  ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
  nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    nx = boout ( 2, 1, ipl ) - boout ( 1, 1, ipl ) + 1
    scount ( ip ) = nx * nz * ny
    sdispl ( ip ) = ny * nz * ( boout ( 1, 1, ipl ) - 1 )
  END DO
  nx = boout ( 2, 1, mip ) - boout ( 1, 1, mip ) + 1
  my = MAXVAL ( boin ( 2, 2, : ) - boin ( 1, 2, : ) + 1 )
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
    rcount ( ip ) = nx * nz * ny
    rdispl ( ip ) = nx * my * nz * ip
  END DO

  ALLOCATE ( rbuf ( nz*nx*my, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", nz*nx*my*np )
  rbuf = 0.0_dp

  CALL mp_alltoall ( sin, scount, sdispl, rbuf, rcount, rdispl, sub_group )

  DO ixz = 1, nx * nz
    DO ip = 0, np - 1
      ipl = pgrid ( ip )
      ny = boin ( 2, 2, ipl ) - boin ( 1, 2, ipl ) + 1
      DO iy = 1, ny
        is = boin ( 1, 2, ipl ) + iy - 1
        ir = iy + ny * ( ixz - 1 )
        sout ( is, ixz ) = rbuf ( ir, ip )
      END DO
    END DO
  END DO

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid" )

  CALL mp_comm_free ( sub_group )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_3

!******************************************************************************

SUBROUTINE cube_transpose_4 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ipl, ir, &
                                                iy, izx, mip, my, n, np, nx, &
                                                ny, nz, sub_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: pgrid, rcount, rdispl, &
                                                scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos, pp
    LOGICAL, DIMENSION(2)                    :: rdim

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_4","I","",handle)

  rdim ( 1 ) = .TRUE.
  rdim ( 2 ) = .FALSE.
  CALL mp_environ ( n, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )
  CALL mp_cart_sub ( group, rdim, sub_group )

  np = dim ( 1 )
  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )
  ALLOCATE ( pgrid ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid", np )

  pp = pos
  DO ip = 0, np - 1
    pp ( 1 ) = ip
    CALL mp_cart_rank ( group, pp, pgrid ( ip ) )
  END DO

  nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
  nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
  my = MAXVAL ( boout ( 2, 2, : ) - boout ( 1, 2, : ) + 1 )

  ALLOCATE ( rbuf ( nz*nx*my, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", nz*nx*my*np )
  rbuf = 0.0_dp

  DO izx = 1, nz * nx
    DO ip = 0, np - 1
      ipl = pgrid ( ip )
      ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
      DO iy = boout ( 1, 2, ipl ), boout ( 2, 2, ipl )
        ir = iy - boout ( 1, 2, ipl ) + 1 + ( izx - 1 ) * ny
        rbuf ( ir, ip ) = sin ( iy, izx )
      END DO
    END DO
  END DO

  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    ny = boout ( 2, 2, ipl ) - boout ( 1, 2, ipl ) + 1
    scount ( ip ) = nx * ny * nz
    sdispl ( ip ) = nx * nz * my * ip
  END DO
  ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
  DO ip = 0, np - 1
    ipl = pgrid ( ip )
    nx = boin ( 2, 1, ipl ) - boin ( 1, 1, ipl ) + 1
    rcount ( ip ) = nx * ny * nz
    rdispl ( ip ) = ny * nz * ( boin ( 1, 1, ipl ) - 1 )
  END DO

  CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, sub_group )

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )
  DEALLOCATE ( pgrid, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "pgrid" )

  CALL mp_comm_free ( sub_group )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_4

!******************************************************************************

SUBROUTINE cube_transpose_5 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ir, is, &
                                                ixz, iy, mip, my, np, nx, ny, &
                                                nz
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_5","I","",handle)

  CALL mp_environ ( np, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )

  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )

  ny = boin ( 2, 2, mip ) - boin ( 1, 2, mip ) + 1
  nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
  DO ip = 0, np - 1
    nx = boout ( 2, 1, ip ) - boout ( 1, 1, ip ) + 1
    scount ( ip ) = nx * nz * ny
    sdispl ( ip ) = ny * nz * ( boout ( 1, 1, ip ) - 1 )
  END DO
  nx = boout ( 2, 1, mip ) - boout ( 1, 1, mip ) + 1
  my = MAXVAL ( boin ( 2, 2, : ) - boin ( 1, 2, : ) + 1 )
  DO ip = 0, np - 1
    ny = boin ( 2, 2, ip ) - boin ( 1, 2, ip ) + 1
    rcount ( ip ) = nx * nz * ny
    rdispl ( ip ) = nx * my * nz * ip
  END DO

  ALLOCATE ( rbuf ( nz*nx*my, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", nz*nx*my*np )
  rbuf = 0.0_dp

  CALL mp_alltoall ( sin, scount, sdispl, rbuf, rcount, rdispl, group )

  DO ixz = 1, nx * nz
    DO ip = 0, np - 1
      ny = boin ( 2, 2, ip ) - boin ( 1, 2, ip ) + 1
      DO iy = 1, ny
        is = boin ( 1, 2, ip ) + iy - 1
        ir = iy + ny * ( ixz - 1 )
        sout ( is, ixz ) = rbuf ( ir, ip )
      END DO
    END DO
  END DO

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_5

!******************************************************************************

SUBROUTINE cube_transpose_6 ( sin, group, boin, boout, sout )


    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sin
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: boin, boout
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    INTEGER                                  :: handle, ierr, ip, ir, iy, &
                                                izx, mip, my, np, nx, ny, nz
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
    INTEGER, DIMENSION(2)                    :: dim, pos

!------------------------------------------------------------------------------

  CALL timeset("cube_transpose_6","I","",handle)

  CALL mp_environ ( np, dim, pos, group )
  CALL mp_cart_rank ( group, pos, mip )

  ALLOCATE ( scount ( 0 : np - 1 ), sdispl ( 0 : np - 1 ), &
             rcount ( 0 : np - 1 ), rdispl ( 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount", 4*np )

  nx = boin ( 2, 1, mip ) - boin ( 1, 1, mip ) + 1
  nz = boin ( 2, 3, mip ) - boin ( 1, 3, mip ) + 1
  my = MAXVAL ( boout ( 2, 2, : ) - boout ( 1, 2, : ) + 1 )

  ALLOCATE ( rbuf ( nz*nx*my, 0 : np - 1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf", nz*nx*my*np )
  rbuf = 0.0_dp

  DO izx = 1, nz * nx
    DO ip = 0, np - 1
      ny = boout ( 2, 2, ip ) - boout ( 1, 2, ip ) + 1
      DO iy = boout ( 1, 2, ip ), boout ( 2, 2, ip )
        ir = iy - boout ( 1, 2, ip ) + 1 + ( izx - 1 ) * ny
        rbuf ( ir, ip ) = sin ( iy, izx )
      END DO
    END DO
  END DO

  DO ip = 0, np - 1
    ny = boout ( 2, 2, ip ) - boout ( 1, 2, ip ) + 1
    scount ( ip ) = nx * ny * nz
    sdispl ( ip ) = nx * nz * my * ip
  END DO
  ny = boout ( 2, 2, mip ) - boout ( 1, 2, mip ) + 1
  DO ip = 0, np - 1
    nx = boin ( 2, 1, ip ) - boin ( 1, 1, ip ) + 1
    rcount ( ip ) = nx * ny * nz
    rdispl ( ip ) = ny * nz * ( boin ( 1, 1, ip ) - 1 )
  END DO

  CALL mp_alltoall ( rbuf, scount, sdispl, sout, rcount, rdispl, group )

  DEALLOCATE ( rbuf, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "rbuf" )
  DEALLOCATE ( scount, sdispl, rcount, rdispl, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "cube_transpose", "scount" )

  CALL timestop(0.0_dp,handle)

END SUBROUTINE cube_transpose_6

!!*****
!******************************************************************************

END MODULE fft_tools

!******************************************************************************
