!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/pw_types [1.0] *
!!
!!   NAME
!!     pw_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (29-12-2000) : Changes for parallel use
!!
!!   NOTES
!!     If parallel mode is distributed certain combination of
!!     "in_use" and "in_space" can not be used.
!!     For performance reasons it would be better to have the loops
!!     over g-vectros in the gather/scatter routines in new subprograms
!!     with the actual arrays (also the adressing) in the parameter list
!!
!!
!!   SOURCE
!******************************************************************************

MODULE pw_types
  
  USE kinds, ONLY : dbl
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE, FULLSPACE
  USE stop_program, ONLY : stop_prg, stop_memory
  USE fft_tools, ONLY : fft3d, FWFFT, BWFFT
  USE timings, ONLY : timeset, timestop
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_type
  PUBLIC :: pw_allocate, pw_deallocate
  PUBLIC :: pw_zero, pw_gather, pw_scatter
  PUBLIC :: pw_fft_wrap
  PUBLIC :: REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D
  PUBLIC :: NOSPACE, REALSPACE, RECIPROCALSPACE
  
  TYPE pw_type
     REAL ( dbl ), DIMENSION ( : ), POINTER :: cr
     REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: cr3d
     COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: cc
     COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: cc3d

     INTEGER :: in_use              ! Which data is used [r1d/c1d/r3d/c3d]
     INTEGER :: in_space            ! Real/Reciprocal space
     
     TYPE ( pw_grid_type ), POINTER :: pw_grid
  END TYPE pw_type
  
  ! Flags for the structure member 'in_use'
  INTEGER, PARAMETER :: REALDATA1D = 301, COMPLEXDATA1D = 302
  INTEGER, PARAMETER :: REALDATA3D = 303, COMPLEXDATA3D = 304
  
  ! Flags for the structure member 'in_space'
  INTEGER, PARAMETER :: NOSPACE = 371, REALSPACE = 372, RECIPROCALSPACE = 373
  
  INTERFACE pw_gather
     MODULE PROCEDURE pw_gather_s, pw_gather_p
  END INTERFACE

  INTERFACE pw_scatter
     MODULE PROCEDURE pw_scatter_s, pw_scatter_p
  END INTERFACE

  INTERFACE pw_fft_wrap
     MODULE PROCEDURE fft_wrap_pw1, fft_wrap_pw1pw2
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE pw_allocate ( pw, pw_grid, use_data )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: use_data
  TYPE ( pw_grid_type ), INTENT ( IN ), TARGET :: pw_grid
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw
  
! Locals
  INTEGER :: allocstat
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  pw % in_use = use_data
  pw % pw_grid => pw_grid
  pw % in_space = NOSPACE
  bounds => pw % pw_grid % bounds_local
  
  NULLIFY ( pw % cr, pw % cc, pw % cr3d, pw % cc3d )
  
  IF ( use_data == REALDATA1D ) THEN
     ALLOCATE ( pw % cr ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cr", pw % pw_grid % ngpts_cut_local )
     END IF
     
  ELSE IF ( use_data == COMPLEXDATA1D ) THEN
     ALLOCATE ( pw % cc ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cc", pw % pw_grid % ngpts_cut_local )
     END IF
     
  ELSE IF ( use_data == REALDATA3D ) THEN
     ALLOCATE ( pw % cr3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cr3d", 0 )
     END IF
     
  ELSE IF ( use_data == COMPLEXDATA3D ) THEN
     ALLOCATE ( pw % cc3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cc3d", 0 )
     END IF
     
  ELSE
     CALL stop_prg ( "pw_allocate", "no possible data field" )
  END IF
  
END SUBROUTINE pw_allocate

!******************************************************************************

SUBROUTINE pw_deallocate ( pw )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  
! Locals
  INTEGER :: allocstat
  
!------------------------------------------------------------------------------
  
  IF ( pw % in_use == REALDATA1D ) THEN
     DEALLOCATE ( pw % cr, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cr" )
     END IF
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     DEALLOCATE ( pw % cc, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cc" )
     END IF
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     DEALLOCATE ( pw % cr3d, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cr3d" )
     END IF
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     DEALLOCATE ( pw % cc3d, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cc3d" )
     END IF
  ELSE
     CALL stop_prg ( "pw_deallocate", "no possible data field" )
  END IF
  
END SUBROUTINE pw_deallocate

!******************************************************************************

SUBROUTINE pw_zero ( pw )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw

! Local
  INTEGER :: ns, handle
  REAL ( dbl ) :: zr
  
!------------------------------------------------------------------------------

  CALL timeset ( 'PW_ZERO', 'I', 'Mzero', handle )
  
  IF ( pw % in_use == REALDATA1D ) THEN
     ns = SIZE ( pw % cr ( : ) )
     pw % cr = 0.0_dbl
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     ns = 2 * SIZE ( pw % cc ( : ) )
     pw % cc = 0.0_dbl
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     ns = SIZE ( pw % cr3d ( :, :, : ) )
     pw % cr3d = 0.0_dbl
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     ns = 2 * SIZE ( pw % cc3d ( :, :, : ) )
     pw % cc3d = 0.0_dbl
  ELSE
     CALL stop_prg ( "pw_zero", "no possible data field" )
  END IF
  
  zr = REAL ( ns, dbl ) * 1.e-6_dbl
  CALL timestop ( zr, handle )

END SUBROUTINE pw_zero

!******************************************************************************

SUBROUTINE pw_gather_s ( pw, c, scale )
 
  IMPLICIT NONE
 
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  COMPLEX (dbl ), DIMENSION ( :, :, : ), INTENT ( IN ) :: c
  COMPLEX ( dbl ), OPTIONAL, INTENT ( IN ) :: scale
 
! Local
  REAL ( dbl ) :: cpy
  INTEGER :: l, m, n, ngpts, gpt,handle
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn
  INTEGER, DIMENSION ( :, : ), POINTER :: ghat

!------------------------------------------------------------------------------

  CALL timeset ( 'PW_GATHER', 'I', 'Mcopy', handle )

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_prg ( "pw_gather_s", "Data field has to be COMPLEXDATA1D" )
  ENDIF

! after the gather we are in g-space
  pw % in_space = RECIPROCALSPACE

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos
 
  ngpts = SIZE ( pw % pw_grid % gsq ( : ) )

  ghat => pw % pw_grid % g_hat

  IF ( PRESENT ( scale ) ) THEN

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      pw % cc ( gpt ) = pw % cc ( gpt ) + scale * c ( l, m, n )

    END DO
 
  ELSE

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      pw % cc ( gpt ) = c ( l, m, n )

    END DO
 
  END IF

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_gather_s

!******************************************************************************

SUBROUTINE pw_gather_p ( pw, c, scale )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw
  COMPLEX (dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: c
  COMPLEX ( dbl ), OPTIONAL, INTENT ( IN ) :: scale

! Local
  REAL ( dbl ) :: cpy
  INTEGER :: l, m, n, mn, ngpts, gpt, handle
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn
  INTEGER, DIMENSION ( :, : ), POINTER :: yzq
  INTEGER, DIMENSION ( :, : ), POINTER :: ghat
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds

!------------------------------------------------------------------------------

  CALL timeset ( 'PW_GATHER', 'I', 'Mcopy', handle )

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_prg ( "pw_gather_p", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % pw_grid % para % mode /= 1 ) THEN
     CALL stop_prg ( "pw_gather_p", "This grid type is not distributed")
  ENDIF

! after the gather we are in g-space
  pw % in_space = RECIPROCALSPACE

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ngpts = SIZE ( pw % pw_grid % gsq ( : ) )

  ghat => pw % pw_grid % g_hat
  yzq => pw % pw_grid % para % yzq
  bounds => pw % pw_grid % bounds

  IF ( PRESENT ( scale ) ) THEN

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
      n = ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
      mn = yzq ( m, n )
      pw % cc ( gpt ) = pw % cc ( gpt ) + scale * c ( l, mn )

    END DO

  ELSE

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
      n = ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
      mn = yzq ( m, n )
      pw % cc ( gpt ) = c ( l, mn )

    END DO

  END IF

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_gather_p

!******************************************************************************

SUBROUTINE pw_scatter_s ( pw, c, scale )
 
  IMPLICIT NONE
 
! Arguments
  TYPE ( pw_type ), INTENT ( IN ) :: pw
  COMPLEX (dbl ), DIMENSION ( :, :, : ), INTENT ( INOUT ) :: c
  COMPLEX ( dbl ), OPTIONAL, INTENT ( IN ) :: scale
 
! Local
  REAL ( dbl ) :: cpy
  INTEGER :: l, m, n, ngpts, gpt, handle
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn
  INTEGER, DIMENSION ( :, : ), POINTER :: ghat

!------------------------------------------------------------------------------

  CALL timeset ( 'PW_SCATTER', 'I', 'Mcopy', handle )

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_prg ( "pw_scatter_s", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % in_space /= RECIPROCALSPACE ) THEN
     CALL stop_prg ( "pw_scatter_s", "Data has to be in RECIPROCALSPACE" )
  ENDIF

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ghat => pw % pw_grid % g_hat

  ngpts = SIZE ( pw % pw_grid % gsq ( : ) )

  IF ( .NOT. PRESENT ( scale ) ) c = 0._dbl
 
  IF ( PRESENT ( scale ) ) THEN

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      c ( l, m, n ) = c ( l, m, n ) + scale * pw % cc ( gpt )

    END DO

  ELSE

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      c ( l, m, n ) = pw % cc ( gpt )

    END DO

  END IF

  IF ( pw % pw_grid % grid_span == HALFSPACE ) THEN

    mapl => pw % pw_grid % mapl % neg
    mapm => pw % pw_grid % mapm % neg
    mapn => pw % pw_grid % mapn % neg

    IF ( PRESENT ( scale ) ) THEN

      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        c ( l, m, n ) = c ( l, m, n ) + scale * CONJG ( pw % cc ( gpt ) )

      END DO

    ELSE

      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        c ( l, m, n ) = CONJG ( pw % cc ( gpt ) )

      END DO

    END IF

  END IF

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_scatter_s

!******************************************************************************

SUBROUTINE pw_scatter_p ( pw, c, scale )
 
  IMPLICIT NONE
 
! Arguments
  TYPE ( pw_type ), INTENT ( IN ), TARGET :: pw
  COMPLEX (dbl ), DIMENSION ( :, : ), INTENT ( INOUT ) :: c
  COMPLEX ( dbl ), OPTIONAL, INTENT ( IN ) :: scale

! Local
  INTEGER :: l, m, n, mn, ngpts, gpt, handle
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn
  INTEGER, DIMENSION ( :, : ), POINTER :: ghat
  INTEGER, DIMENSION ( :, : ), POINTER :: yzq
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  REAL ( dbl ) :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'PW_SCATTER', 'I', 'Mcopy', handle )
 
  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_prg ( "pw_scatter_p", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % in_space /= RECIPROCALSPACE ) THEN
     CALL stop_prg ( "pw_scatter_p", "Data has to be in RECIPROCALSPACE" )
  ENDIF

  IF ( pw % pw_grid % para % mode /= 1 ) THEN
     CALL stop_prg ( "pw_scatter_p", "This grid type is not distributed")
  ENDIF

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ghat => pw % pw_grid % g_hat
  yzq => pw % pw_grid % para % yzq
  bounds => pw % pw_grid % bounds

  ngpts = SIZE ( pw % pw_grid % gsq ( : ) )

  IF ( .NOT. PRESENT ( scale ) ) c = 0._dbl

  IF ( PRESENT ( scale ) ) THEN

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
      n = ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
      mn = yzq ( m, n )
      c ( l, mn ) = c ( l, mn ) + scale * pw % cc ( gpt )

    END DO

  ELSE

    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
      n = ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
      mn = yzq ( m, n )
      c ( l, mn ) = pw % cc ( gpt )

    END DO

  END IF

  IF ( pw % pw_grid % grid_span == HALFSPACE ) THEN

    mapm => pw % pw_grid % mapm % neg
    mapn => pw % pw_grid % mapn % neg
    mapl => pw % pw_grid % mapl % neg

    IF ( PRESENT ( scale ) ) THEN

      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = -ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
        n = -ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
        mn = yzq ( m, n )
        c ( l, mn ) = c ( l, mn ) + scale * CONJG ( pw % cc ( gpt ) )

      END DO

    ELSE

      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = -ghat ( 2, gpt ) - bounds ( 1, 2 ) + 1
        n = -ghat ( 3, gpt ) - bounds ( 1, 3 ) + 1
        mn = yzq ( m, n )
        c ( l, mn ) = CONJG ( pw % cc ( gpt ) )

      END DO

    END IF

  END IF

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_scatter_p

!******************************************************************************
!!****** pw_types/pw_fft_wrap [1.1] *
!!
!!   NAME
!!     pw_fft_wrap
!!
!!   FUNCTION
!!     Generic function for 3d FFT of a coefficient_type or pw_type
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (30-12-2000): New setup of functions and adaptation to parallelism
!!     JGH (04-01-2001): Moved routine from pws to this module, only covers
!!                       pw_types, no more coefficient types
!!
!!   NOTES
!!     The following functions are covered
!!      fft_wrap_pw1, fft_wrap_pw1pw2
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft_wrap_pw1 ( pw1 )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw1

! Local
  INTEGER :: dir, out_space
  REAL ( dbl ) :: norm
!------------------------------------------------------------------------------

  IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
       pw1 % pw_grid % para % mode == 0 ) THEN
!..what dirction will the transform be?
    IF ( pw1 % in_space == REALSPACE ) THEN
       dir = FWFFT
       norm = pw1 % pw_grid % dvol
       out_space = RECIPROCALSPACE
    ELSEIF ( pw1 % in_space == RECIPROCALSPACE ) THEN
       dir = BWFFT
       norm = 1._dbl
       out_space = REALSPACE
    ELSE
       CALL stop_prg ( "fft_wrap_pw1", "PW structure is missing a", &
                       "proper tag to identidy its space" )
    END IF

    CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, scale = norm )

!..tag new data with correct space
    pw1 % in_space = out_space
  ELSE
    CALL stop_prg ( "fft_wrap_pw1", "In place FFT only possible for", &
                    "replicated data with COMPLEXDATA3D structure" )
  END IF

END SUBROUTINE fft_wrap_pw1


!******************************************************************************

SUBROUTINE fft_wrap_pw1pw2 ( pw1, pw2 )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( IN ), TARGET :: pw1
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw2
  
! Locals
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: c_in
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: c_out
  COMPLEX ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: grays
  REAL ( dbl ) :: norm
  INTEGER, DIMENSION ( : ), POINTER :: n
  INTEGER :: dir, out_space, ierr, nrays, my_pos
  CHARACTER ( LEN = 9 ) :: mode
  
!------------------------------------------------------------------------------
  
  NULLIFY ( c_in )
  NULLIFY ( c_out )

!..check if grids are compatible
  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) THEN
    IF ( pw1 % pw_grid % dvol /= pw2 % pw_grid % dvol ) THEN
      CALL stop_prg ( "fft_wrap_pw1pw2", "PW grids not compatible" )
    END IF
    IF ( pw1 % pw_grid % para %group /= pw2 % pw_grid % para % group) THEN
      CALL stop_prg ( "fft_wrap_pw1pw2", &
                      "PW grids have not compatible MPI groups" )
    END IF
  END IF

!..prepare input
  IF ( pw1 % in_space == REALSPACE ) THEN
    dir = FWFFT
    norm = 1._dbl / pw1 % pw_grid % ngpts
    out_space = RECIPROCALSPACE
  ELSE IF ( pw1 % in_space == RECIPROCALSPACE ) THEN
    dir = BWFFT
    norm = 1._dbl
    out_space = REALSPACE
  ELSE
    CALL stop_prg ( "fft_wrap_pw1pw2", "Error in space tag" )
  END IF

  n => pw1 % pw_grid % npts

  mode = fftselect ( pw1 % in_use, pw2 % in_use, pw1 % in_space )

  IF ( pw1 % pw_grid % para % mode == 0 ) THEN

!
!..replicated data, use local FFT
!

    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_prg ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space" )
      CASE ( "FW_C3DC3D" )
        c_in => pw1 % cc3d
        c_out => pw2 % cc3d
      CASE ( "FW_R3DC3D" )
        pw2 % cc3d = CMPLX ( pw1 % cr3d )
        c_out => pw2 % cc3d
      CASE ( "FW_C3DC1D" )
        c_in => pw1 % cc3d
        ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
      CASE ( "FW_R3DC1D" )
        ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
        c_out = CMPLX ( pw1 % cr3d )
      CASE ( "BW_C3DC3D" )
        c_in => pw1 % cc3d
        c_out => pw2 % cc3d
      CASE ( "BW_C3DR3D" )
        c_in => pw1 % cc3d
        ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
      CASE ( "BW_C1DC3D" )
        c_out => pw2 % cc3d
        CALL pw_scatter ( pw1, c_out )
      CASE ( "BW_C1DR3D" )
        ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
        CALL pw_scatter ( pw1, c_out )
    END SELECT

!..transform
    IF ( ASSOCIATED ( c_in ) ) THEN
      CALL fft3d ( dir, n, c_in, c_out, scale = norm )
    ELSE
      CALL fft3d ( dir, n, c_out, scale = norm )
    ENDIF
     
!..prepare output
    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_prg ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space" )
      CASE ( "FW_C3DC3D" )
        !     nothing to do
      CASE ( "FW_R3DC3D" )
        !     nothing to do
      CASE ( "FW_C3DC1D" )
        !     gather results
        CALL pw_gather ( pw2, c_out )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "FW_R3DC1D" )
        !     gather results
        CALL pw_gather ( pw2, c_out )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "BW_C3DC3D" )
        !     nothing to do
      CASE ( "BW_C3DR3D" )
        !     use real part only
        pw2 % cr3d = REAL ( c_out, dbl )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "BW_C1DC3D" )
        !     nothing to do
      CASE ( "BW_C1DR3D" )
        !     use real part only
        pw2 % cr3d = REAL ( c_out, dbl )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
    END SELECT

  ELSE
!
!..parallel FFT
!

    my_pos = pw1 % pw_grid % para % my_pos
    nrays = pw1 % pw_grid % para % nyzray ( my_pos )
    ALLOCATE ( grays ( n(1), nrays ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", &
      "grays", n(1) * nrays )
    grays = 0._dbl

    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_prg ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space", &
                        "in parallel 3d FFT" )
      CASE ( "FW_C3DC1D" )
        c_in => pw1 % cc3d
      CASE ( "FW_R3DC1D" )
        n => pw1 % pw_grid % npts_local
        ALLOCATE ( c_in( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in", &
            pw1 % pw_grid % ngpts_local )
        c_in = CMPLX ( pw1 % cr3d )
      CASE ( "BW_C1DC3D" )
        CALL pw_scatter ( pw1, grays )
        c_in => pw2 % cc3d
      CASE ( "BW_C1DR3D" )
        CALL pw_scatter ( pw1, grays )
        n => pw2 % pw_grid % npts_local
        ALLOCATE ( c_in( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in", &
            pw2 % pw_grid % ngpts_local )
    END SELECT

!..transform
!debCALL pfft3d ( dir, n, c_in, grays, gs_group, rs_group, &
!deb              pw1 % pw_grid % bounds_local, &
!deb              pw1 % pw_grid % para % yzp, scale = norm )

!..prepare output
    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_prg ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space", &
                        "in parallel 3d FFT" )
      CASE ( "FW_C3DC1D" )
        CALL pw_gather ( pw2, grays )
      CASE ( "FW_R3DC1D" )
        CALL pw_gather ( pw2, grays )
        DEALLOCATE ( c_in, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in" )
      CASE ( "BW_C1DC3D" )
        ! nothing to do
      CASE ( "BW_C1DR3D" )
        pw2 % cr3d = REAL ( c_in, dbl )
        DEALLOCATE ( c_in, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in" )
    END SELECT

    DEALLOCATE ( grays, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "grays" )

  END IF

!..update the space tag for pw2
  pw2 % in_space = out_space

END SUBROUTINE fft_wrap_pw1pw2

!******************************************************************************

FUNCTION fftselect ( use1, use2, space1 ) RESULT ( mode ) 
 
  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: use1, use2, space1 

  CHARACTER ( LEN = 9 ) :: mode
  
!------------------------------------------------------------------------------
  
  IF ( space1 == REALSPACE ) THEN
    mode ( 1 : 3 ) = "FW_"
  ELSE IF ( space1 == RECIPROCALSPACE ) THEN
    mode ( 1 : 3 ) = "BW_"
  ELSE
    CALL stop_prg ( "fftselect", "Error in space tag" )
  END IF

  SELECT CASE ( use1 )
    CASE ( COMPLEXDATA3D )
      mode ( 4 : 6 ) = "C3D"
    CASE ( REALDATA3D )
      mode ( 4 : 6 ) = "R3D"
    CASE ( COMPLEXDATA1D )
      mode ( 4 : 6 ) = "C1D"
    CASE ( REALDATA1D )
      mode ( 4 : 6 ) = "R1D"
    CASE DEFAULT
      CALL stop_prg ( "fftselect", "Error in use1 tag" )
  END SELECT

  SELECT CASE ( use2 )
    CASE ( COMPLEXDATA3D )
      mode ( 7 : 9 ) = "C3D"
    CASE ( REALDATA3D )
      mode ( 7 : 9 ) = "R3D"
    CASE ( COMPLEXDATA1D )
      mode ( 7 : 9 ) = "C1D"
    CASE ( REALDATA1D )
      mode ( 7 : 9 ) = "R1D"
    CASE DEFAULT
      CALL stop_prg ( "fftselect", "Error in use1 tag" )
  END SELECT

END FUNCTION fftselect

!!*****
!******************************************************************************
  
END MODULE pw_types

!******************************************************************************
