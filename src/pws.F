!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE pws
  
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, coeff_deallocate, &
       SQUARE, SQUAREROOT
  USE fft_tools, ONLY : fft3d, &
       FWFFTu => FWFFT, BWFFTu => BWFFT, FWFFTw => FWFFT, BWFFTw => BWFFT
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : fourpi
  USE message_passing, ONLY : mp_sum
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE
  USE pw_types, ONLY : pw_type, pw_allocate, pw_deallocate, &
       REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D
  USE stop_program, ONLY : stop_prg
  USE util, ONLY : dotprod_3d
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_poisson_solver, pw_kinetic_evaluate
  PUBLIC :: pw_integral_a, pw_integral_ab, pw_h_psi, fft_wrap
  
  INTERFACE fft_wrap
     MODULE PROCEDURE fft_wrap_c1c2, fft_wrap_c1pw2, &
          fft_wrap_pw1c2, fft_wrap_pw1pw2
  END INTERFACE
  
CONTAINS

!******************************************************************************

SUBROUTINE pw_poisson_solver ( n, vhartree, dvhartree )
  
!! obtains electrostatic potential and its derivatives with respect to r
!!    from the density
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ), TARGET :: n
  TYPE ( coeff_type ), INTENT ( OUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( OUT ), OPTIONAL :: dvhartree
  
! Locals
  INTEGER :: gpt, dir
  REAL ( dbl ) :: gf
  COMPLEX ( dbl ) :: i_cmplx
  INTEGER :: lp, mp, np, ln, mn, nn
  INTEGER, DIMENSION ( : ), POINTER :: ghat
  
  INTEGER, DIMENSION ( :, : ), POINTER :: g_hat
  TYPE ( pw_grid_type ), POINTER :: pw_grid
  TYPE ( pw_type ) :: rho_g, cfft
  
!------------------------------------------------------------------------------
  
  i_cmplx = ( 0._dbl, 1._dbl )
  
  pw_grid => n % pw % pw_grid
  g_hat => pw_grid % g_hat
  
  CALL pw_allocate ( rho_g, pw_grid, use_data = COMPLEXDATA3D )
  
  CALL pw_allocate ( cfft, pw_grid, use_data = COMPLEXDATA3D )
  
  CALL fft_wrap ( n, rho_g, dir = FWFFTu )
  
  IF ( PRESENT ( vhartree ) ) THEN
     cfft % cc3d = 0._dbl
     
     DO gpt = 1, pw_grid % ngpts_cut
        
        ! compute Green's function  '''
        CALL pw_greens_function ( gf, pw_grid % g ( :, gpt ) )
        
        ghat => g_hat ( :, gpt )
        
        lp = pw_grid % mapl % pos ( ghat ( 1 ) )
        ln = pw_grid % mapl % neg ( ghat ( 1 ) )
        mp = pw_grid % mapm % pos ( ghat ( 2 ) )
        mn = pw_grid % mapm % neg ( ghat ( 2 ) )
        np = pw_grid % mapn % pos ( ghat ( 3 ) )
        nn = pw_grid % mapn % neg ( ghat ( 3 ) )
        
        lp = lp + LBOUND ( cfft % cc3d, 1 )
        mp = mp + LBOUND ( cfft % cc3d, 2 )
        np = np + LBOUND ( cfft % cc3d, 3 )
        ln = ln + LBOUND ( cfft % cc3d, 1 )
        mn = mn + LBOUND ( cfft % cc3d, 2 )
        nn = nn + LBOUND ( cfft % cc3d, 3 )
        
        cfft % cc3d ( lp, mp, np ) = gf * rho_g % cc3d ( lp, mp, np )
        cfft % cc3d ( ln, mn, nn ) = gf * rho_g % cc3d ( ln, mn, nn )
     END DO
     
! get phi_r
     CALL fft_wrap ( cfft, vhartree, dir = BWFFTu )
  END IF
  
! multiply by the Green's function   '''
  
  IF ( PRESENT ( dvhartree ) ) THEN
     DO dir = 1, 3
        cfft % cc3d = 0.0_dbl
        
        DO gpt = 1, pw_grid % ngpts_cut
           
           ! compute Green's function   '''
           CALL pw_greens_function ( gf, pw_grid % g ( :, gpt ) )
           
           ghat => g_hat ( :, gpt )
           
           lp = pw_grid % mapl % pos ( ghat ( 1 ) )
           ln = pw_grid % mapl % neg ( ghat ( 1 ) )
           mp = pw_grid % mapm % pos ( ghat ( 2 ) )
           mn = pw_grid % mapm % neg ( ghat ( 2 ) )
           np = pw_grid % mapn % pos ( ghat ( 3 ) )
           nn = pw_grid % mapn % neg ( ghat ( 3 ) )
            
           lp = lp + LBOUND ( cfft % cc3d, 1 )
           mp = mp + LBOUND ( cfft % cc3d, 2 )
           np = np + LBOUND ( cfft % cc3d, 3 )
           ln = ln + LBOUND ( cfft % cc3d, 1 )
           mn = mn + LBOUND ( cfft % cc3d, 2 )
           nn = nn + LBOUND ( cfft % cc3d, 3 )
           
           cfft % cc3d ( lp, mp, np ) = &
                  i_cmplx * pw_grid % g ( dir, gpt ) * gf &
                * rho_g % cc3d ( lp, mp, np )
           cfft % cc3d ( ln, mn, nn ) = &
                - i_cmplx * pw_grid % g ( dir, gpt ) * gf &
                * rho_g % cc3d ( ln, mn, nn )
        END DO
        
! get derivative of phi_r
        CALL fft_wrap ( cfft, dvhartree ( dir ), dir = BWFFTu )
     END DO
  END IF
  
  CALL pw_deallocate ( cfft )
  
  CALL pw_deallocate ( rho_g )
  
END SUBROUTINE pw_poisson_solver

!******************************************************************************

SUBROUTINE pw_greens_function ( gf, g )
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: g
  REAL ( dbl ), INTENT ( OUT ) :: gf
  
! Locals
  INTEGER :: dimensionality  !TMPTMPTMP
  REAL ( dbl ) :: g2
  
!------------------------------------------------------------------------------
  
  dimensionality = SIZE ( g )
  
  SELECT CASE ( dimensionality )
  CASE ( 3 )                   ! 3-dim
     g2 = dotprod_3d ( g, g )
     IF ( g2 > 1.0E-10_dbl ) THEN
        gf = fourpi / g2
     ELSE
        gf = 0.0_dbl
     END IF
  CASE ( 2 )                   ! 2-dim
     stop "greens_function, no 2-dim function implemented yet"
  CASE ( 1 )                   ! 1-dim
     stop "greens_function, no 1-dim function implemented yet"
  CASE ( 0 )                   ! 0-dim = cluster
     stop "greens_function, no 0-dim function implemented yet"
  END SELECT
  
END SUBROUTINE pw_greens_function

!******************************************************************************

SUBROUTINE pw_allocate_coeff ( pw_out, pw_model, use_data )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( IN ), TARGET :: pw_model
  TYPE ( pw_type ), INTENT ( OUT ) :: pw_out
  INTEGER, INTENT ( IN ), OPTIONAL :: use_data
  
! Locals
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( use_data ) ) THEN
     pw_out % in_use = use_data
  ELSE
     pw_out % in_use = pw_model % in_use
  END IF
  pw_out % pw_grid => pw_model % pw_grid
  bounds => pw_model % pw_grid % bounds
  
  CALL pw_allocate ( pw_out, pw_out % pw_grid, use_data = use_data )
  
END SUBROUTINE pw_allocate_coeff

!******************************************************************************

SUBROUTINE pw_deallocate_coeff ( pw )
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  
!------------------------------------------------------------------------------
  
  IF ( pw % in_use == REALDATA1D ) THEN
     DEALLOCATE ( pw % cr )
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     DEALLOCATE ( pw % cc )
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     DEALLOCATE ( pw % cr3d )
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     DEALLOCATE ( pw % cc3d )
  ELSE 
     CALL stop_prg ( "pw_deallocate_coeff", "error de-allocating data" )
  END IF
  
END SUBROUTINE pw_deallocate_coeff

!******************************************************************************

SUBROUTINE fft_wrap_c1c2 ( c1, c2, dir, scale )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: c1
  TYPE ( coeff_type ), INTENT ( INOUT ), OPTIONAL, TARGET :: c2
  INTEGER, INTENT ( IN ), OPTIONAL :: dir
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  
! Locals
  INTEGER :: gpt, ghat ( 3 )
  INTEGER, POINTER :: bds ( :, : )
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: cfft
  TYPE ( pw_type ), POINTER :: pw1, pw2
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( c2 ) ) THEN
     pw1 => c1 % pw
     pw2 => c2 % pw
     
     IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        IF ( .NOT. PRESENT ( dir ) ) THEN
           CALL stop_prg ( "fft_wrap_c1c2", &
                 "direction must be given in 3D => 3D transform" )
        ELSE
           CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, pw2 % cc3d, &
                scale = scale )
        END IF
        
     ELSE IF ( pw1 % in_use == REALDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        
        pw2 % cc3d = pw1 % cr3d / REAL ( pw1 % pw_grid % ngpts, dbl )
        CALL fft3d ( FWFFTu, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == REALDATA3D ) THEN
        
        bds => pw1 % pw_grid % bounds
        ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
        cfft = pw1 % cc3d
        
        CALL fft3d ( BWFFTu, pw1 % pw_grid % npts, cfft, scale = scale )
        pw2 % cr3d = REAL ( cfft )
        DEALLOCATE ( cfft )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        
        pw2 % cc3d = 0._dbl
        
        DO gpt = 1, pw1 % pw_grid % ngpts_cut
           ghat = pw2 % pw_grid % g_hat ( :, gpt )
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
!TMPTMPTMP                pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
!TMPTMPTMP                pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           END IF
        END DO
        
        CALL fft3d ( BWFFTw, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
        
        bds => pw1 % pw_grid % bounds
        ALLOCATE (cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)) )
        cfft = pw1 % cc3d / REAL ( pw1 % pw_grid % ngpts, dbl )
        
        CALL fft3d ( FWFFTw, pw1 % pw_grid % npts, cfft, scale = scale )
        
        DO gpt = 1, pw2 % pw_grid % ngpts_cut
           ghat = pw1 % pw_grid % g_hat ( :, gpt )
           pw2 % cc ( gpt ) = cfft ( &
                LBOUND ( cfft, 1 ) + pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
                LBOUND ( cfft, 2 ) + pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
                LBOUND ( cfft, 3 ) + pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
!TMPTMPTMP                pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
!TMPTMPTMP                pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
!TMPTMPTMP                pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
        END DO
        
        DEALLOCATE ( cfft )
        
     ELSE
        CALL stop_prg ( "fft_wrap_c1c2", "no suitable FFT found" )
     END IF
  ELSE
     pw1 => c1 % pw
     IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        IF ( PRESENT ( dir ) ) THEN
           CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, scale = scale )
        ELSE
           CALL stop_prg ( "fft_wrap_c1c2", &
                "no direction specified for in-place FFT" )
        END IF
     ELSE
        CALL stop_prg ( "fft_wrap_c1c2", "wrong data type to be FFTed" )
     END IF
  END IF
  
END SUBROUTINE fft_wrap_c1c2

!******************************************************************************

SUBROUTINE fft_wrap_c1pw2 ( c1, pw2, dir, scale )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ), TARGET :: c1
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw2
  INTEGER, INTENT ( IN ), OPTIONAL :: dir
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  
! Locals
  INTEGER :: gpt, ghat ( 3 )
  INTEGER, POINTER :: bds ( :, : )
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: cfft
  TYPE ( pw_type ), POINTER :: pw1
  
!------------------------------------------------------------------------------
  
  pw1 => c1 % pw
  
  IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     IF ( .NOT. PRESENT ( dir ) ) THEN
        CALL stop_prg ( "fft_wrap_c1pw2", &
           "direction must be given in 3D => 3D transform" )
     ELSE
        CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, pw2 % cc3d, &
             scale = scale )
     END IF
     
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     
     pw2 % cc3d = pw1 % cr3d / REAL ( pw1 % pw_grid % ngpts, dbl )
     CALL fft3d ( FWFFTu, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     bds => pw1 % pw_grid % bounds
     ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
     cfft = pw1 % cc3d
     CALL fft3d ( BWFFTu, pw1 % pw_grid % npts, cfft, scale = scale )
     pw2 % cr3d = REAL ( cfft )
     DEALLOCATE ( cfft )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     
     pw2 % cc3d = 0._dbl
     
     DO gpt = 1, pw1 % pw_grid % ngpts_cut
        ghat = pw2 % pw_grid % g_hat ( :, gpt )
        pw2 % cc3d ( &
             LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
             LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
             LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        !TMPTMPTMP                pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
        !TMPTMPTMP                pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
        !TMPTMPTMP                pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        
        IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           !TMPTMPTMP                pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
           !TMPTMPTMP                pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
           !TMPTMPTMP                pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        END IF
     END DO
     
     CALL fft3d ( BWFFTw, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     
     bds => pw1 % pw_grid % bounds
     ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
     cfft = pw1 % cc3d / REAL ( pw1 % pw_grid % ngpts, dbl )
     
     CALL fft3d ( FWFFTw, pw1 % pw_grid % npts, cfft, scale = scale )
     
     DO gpt = 1, pw2 % pw_grid % ngpts_cut
        ghat = pw1 % pw_grid % g_hat ( :, gpt )
        pw2 % cc ( gpt ) = cfft ( &
             LBOUND ( cfft, 1 ) + pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
             LBOUND ( cfft, 2 ) + pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
             LBOUND ( cfft, 3 ) + pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
        !TMPTMPTMP                pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
        !TMPTMPTMP                pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
        !TMPTMPTMP                pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
     END DO
     
     DEALLOCATE ( cfft )
     
  ELSE
     CALL stop_prg ( "fft_wrap_c1pw2", "no suitable FFT found" )
  END IF
  
END SUBROUTINE fft_wrap_c1pw2

!******************************************************************************

SUBROUTINE fft_wrap_pw1c2 ( pw1, c2, dir, scale )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw1
  TYPE ( coeff_type ), INTENT ( INOUT ), TARGET :: c2
  INTEGER, INTENT ( IN ), OPTIONAL :: dir
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  
! Locals
  INTEGER :: gpt, ghat ( 3 )
  INTEGER, POINTER :: bds ( :, : )
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: cfft
  TYPE ( pw_type ), POINTER :: pw2
  
!------------------------------------------------------------------------------
  
  pw2 => c2 % pw
  
  IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     IF ( .NOT. PRESENT ( dir ) ) THEN
        CALL stop_prg ( "fft_wrap_pw1c2", &
           "direction must be given in 3D => 3D transform" )
     ELSE
        CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, pw2 % cc3d, &
             scale = scale )
     END IF
     
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     
     pw2 % cc3d = pw1 % cr3d / REAL ( pw1 % pw_grid % ngpts, dbl )
     CALL fft3d ( FWFFTu, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     bds => pw1 % pw_grid % bounds
     ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
     cfft = pw1 % cc3d
     CALL fft3d ( BWFFTu, pw1 % pw_grid % npts, cfft, scale = scale )
     pw2 % cr3d = REAL ( cfft )
     DEALLOCATE ( cfft )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     
     pw2 % cc3d = 0._dbl
     
     DO gpt = 1, pw1 % pw_grid % ngpts_cut
        ghat = pw2 % pw_grid % g_hat ( :, gpt )
        pw2 % cc3d ( &
             LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
             LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
             LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        !TMPTMPTMP                pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
        !TMPTMPTMP                pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
        !TMPTMPTMP                pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           !TMPTMPTMP                pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
           !TMPTMPTMP                pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
           !TMPTMPTMP                pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
        END IF
     END DO
     
     CALL fft3d ( BWFFTw, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
     
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     
     bds => pw1 % pw_grid % bounds
     ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
     cfft = pw1 % cc3d / REAL ( pw1 % pw_grid % ngpts, dbl )
     
     CALL fft3d ( FWFFTw, pw1 % pw_grid % npts, cfft, scale = scale )
     
     DO gpt = 1, pw2 % pw_grid % ngpts_cut
        ghat = pw1 % pw_grid % g_hat ( :, gpt )
        pw2 % cc ( gpt ) = cfft ( &
             LBOUND ( cfft, 1 ) + pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
             LBOUND ( cfft, 2 ) + pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
             LBOUND ( cfft, 3 ) + pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
        !TMPTMPTMP                pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
        !TMPTMPTMP                pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
        !TMPTMPTMP                pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
     END DO
     
     DEALLOCATE ( cfft )
     
  ELSE
     CALL stop_prg ( "fft_wrap_pw1c2", "no suitable FFT found" )
  END IF
  
END SUBROUTINE fft_wrap_pw1c2

!******************************************************************************

SUBROUTINE fft_wrap_pw1pw2 ( pw1, pw2, dir, scale )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( OUT ), TARGET :: pw1
  TYPE ( pw_type ), INTENT ( INOUT ), OPTIONAL, TARGET :: pw2
  INTEGER, INTENT ( IN ), OPTIONAL :: dir
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  
! Locals
  INTEGER :: gpt, ghat ( 3 )
  INTEGER, POINTER :: bds ( :, : )
  COMPLEX ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: cfft
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( pw2 ) ) THEN
     
     IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        IF ( .NOT. PRESENT ( dir ) ) THEN
           CALL stop_prg ( "fft_wrap_pw1pw2", &
                 "direction must be given in 3D => 3D transform" )
        ELSE
           CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, pw2 % cc3d, &
                scale = scale )
        END IF
        
     ELSE IF ( pw1 % in_use == REALDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        
        pw2 % cc3d = pw1 % cr3d / REAL ( pw1 % pw_grid % ngpts, dbl )
        CALL fft3d ( FWFFTu, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == REALDATA3D ) THEN
        bds => pw1 % pw_grid % bounds
        ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
        cfft = pw1 % cc3d
        CALL fft3d ( BWFFTu, pw1 % pw_grid % npts, cfft, scale = scale )
        pw2 % cr3d = REAL ( cfft )
        DEALLOCATE ( cfft )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
          .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
        
        pw2 % cc3d = 0._dbl
        
        DO gpt = 1, pw1 % pw_grid % ngpts_cut
           ghat = pw2 % pw_grid % g_hat ( :, gpt )
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
!TMPTMPTMP                pw2 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapn % pos ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
           pw2 % cc3d ( &
                LBOUND ( pw2 % cc3d, 1 ) + pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
                LBOUND ( pw2 % cc3d, 2 ) + pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
                LBOUND ( pw2 % cc3d, 3 ) + pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
!TMPTMPTMP                pw2 % pw_grid % mapl % neg ( ghat ( 1 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapm % neg ( ghat ( 2 ) ), &
!TMPTMPTMP                pw2 % pw_grid % mapn % neg ( ghat ( 3 ) ) ) = pw1 % cc ( gpt )
           END IF
        END DO
        
        CALL fft3d ( BWFFTw, pw2 % pw_grid % npts, pw2 % cc3d, scale = scale )
        
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
          .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
        
        bds => pw1 % pw_grid % bounds
        ALLOCATE(cfft(bds(1,1):bds(2,1),bds(1,2):bds(2,2),bds(1,3):bds(2,3)))
        cfft = pw1 % cc3d / REAL ( pw1 % pw_grid % ngpts, dbl )
        
        CALL fft3d ( FWFFTw, pw1 % pw_grid % npts, cfft, scale = scale )
        
        DO gpt = 1, pw2 % pw_grid % ngpts_cut
           ghat = pw1 % pw_grid % g_hat ( :, gpt )
           pw2 % cc ( gpt ) = cfft ( &
                LBOUND ( cfft, 1 ) + pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
                LBOUND ( cfft, 2 ) + pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
                LBOUND ( cfft, 3 ) + pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
!TMPTMPTMP                pw1 % pw_grid % mapl % pos ( ghat ( 1 ) ), &
!TMPTMPTMP                pw1 % pw_grid % mapm % pos ( ghat ( 2 ) ), &
!TMPTMPTMP                pw1 % pw_grid % mapn % pos ( ghat ( 3 ) ) )
        END DO
        
        DEALLOCATE ( cfft )
        
     ELSE
        CALL stop_prg ( "fft_wrap_pw1pw2", "no suitable FFT found" )
     END IF
     
  ELSE
     IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        IF ( PRESENT ( dir ) ) THEN
           CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, scale = scale )
        ELSE
           CALL stop_prg ( "fft_wrap_pw1pw2", &
                "no direction specified for in-place FFT" )
        END IF
     ELSE
        CALL stop_prg ( "fft_wrap_pw1pw2", "wrong data type to be FFTed" )
     END IF
  END IF
  
END SUBROUTINE fft_wrap_pw1pw2

!******************************************************************************

FUNCTION pw_integral_ab ( pw1, pw2, flag ) RESULT ( integral_value )
  
  TYPE ( pw_type ), INTENT ( IN ) :: pw1, pw2
  INTEGER, INTENT ( IN ), OPTIONAL :: flag
  
  REAL ( dbl ) :: integral_value
  
!------------------------------------------------------------------------------
  
  IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) ) 
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cc3d ( :, :, : ) ) 
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cc3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) ) 
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
          * pw2 % cc3d ( :, :, : ) ) 
     
  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cr ( : ) ) 
  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cc ( : ) ) 
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cc ( : ) * pw2 % cr ( : ) ) 
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) * pw2 % cc ( : ) ) 
  ELSE
     CALL stop_prg ( "integral_ab", "no possible data" )
  END IF
  
  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
     integral_value = integral_value
  ELSE
     CALL stop_prg ( "pw_integral_ab", "volume element not implemented" )
  END IF
  
  IF ( pw1 % group_id >= 0 ) CALL mp_sum ( integral_value, pw1 % group_id )
  
END FUNCTION pw_integral_ab

!******************************************************************************

FUNCTION pw_integral_a ( pw1, flag ) RESULT ( integral_value )
  
  TYPE ( pw_type ), INTENT ( IN ) :: pw1
  INTEGER, INTENT ( IN ), OPTIONAL :: flag
  
  REAL ( dbl ) :: integral_value
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( flag ) ) THEN
     IF ( flag == SQUARE ) THEN
        IF ( pw1 % in_use == REALDATA3D ) THEN
           integral_value = SUM ( pw1 % cr3d ( :, :, : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
                * pw1 % cc3d ( :, :, : ) )
        ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
           integral_value = SUM ( pw1 % cr ( : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) &
                * pw1 % cc ( : ) )
        ELSE
           CALL stop_prg ( "integral_aa", "no possible SQUARE data" )
        END IF
        
     ELSE IF ( flag == SQUAREROOT ) THEN
        CALL stop_prg ( "integral_aa", "no SQUAREROOT defined" )
     END IF
     
  ELSE
     IF ( pw1 % in_use == REALDATA3D ) THEN
        integral_value = SUM ( pw1 % cr3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        integral_value = SUM ( pw1 % cc3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
        integral_value = SUM ( pw1 % cr ( : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
        integral_value = SUM ( pw1 % cc ( : ) )
     ELSE
        CALL stop_prg ( "integral_aa", "no possible data" )
     END IF
  END IF
  
  IF ( pw1 % in_use == REALDATA3D ) &
       integral_value = integral_value * pw1 % pw_grid % dvol
  
  IF ( pw1 % group_id >= 0 ) CALL mp_sum ( integral_value, pw1 % group_id )
  
END FUNCTION pw_integral_a

!******************************************************************************

SUBROUTINE pw_kinetic_evaluate ( psi, xkpt, ekinetic, kinetic_psi )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: psi
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  
  REAL ( dbl ), INTENT ( OUT ) :: ekinetic
  TYPE ( coeff_type ), INTENT ( INOUT ) :: kinetic_psi
  
! Local
  INTEGER :: gpt
  
!------------------------------------------------------------------------------
  
  IF ( psi % pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_prg ( "pw_calculate_kinetic", "wrong data in use" )
  END IF
  
  !OpenMP: parallel do
  DO gpt = 1, psi % pw % pw_grid % ngpts_cut
     kinetic_psi % pw % cc ( gpt ) = kinetic_psi % pw % cc ( gpt ) &
          + 0.5_dbl &
          * SUM ( ( psi % pw % pw_grid % g ( :, gpt ) + xkpt ( : ) ) ** 2 ) &
          * psi % pw % cc ( gpt )
  END DO
  
  ekinetic = pw_integral_ab ( psi % pw, kinetic_psi % pw )
  
END SUBROUTINE pw_kinetic_evaluate

!******************************************************************************

SUBROUTINE pw_h_psi ( psi_iks, xkpt, vloc, h_psi, psi_iks_r )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: psi_iks
  TYPE ( coeff_type ), INTENT ( IN ), OPTIONAL :: psi_iks_r
  TYPE ( coeff_type ), INTENT ( IN ) :: vloc
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  
  TYPE ( coeff_type ), INTENT ( INOUT ) :: h_psi
  
! Locals
  INTEGER :: gpt, gpt2, i1, i2, i3, nhits !TMPTMPTMP
  REAL ( dbl ) :: ekinetic !TMPTMPTMP
!  COMPLEX ( dbl ) :: foo !TMPTMPTMP
  TYPE ( pw_type ) :: cfft
  
!------------------------------------------------------------------------------
  
  CALL pw_allocate ( cfft, psi_iks % pw % pw_grid, &
       use_data = COMPLEXDATA3D )
  
#if defined ( DEBUG )
  nhits = 0
  OPEN ( UNIT = 99, FILE = "psik", FORM = "FORMATTED" )
  DO gpt = 1, psi_iks % grid % ngpts
     READ ( UNIT = 99, FMT = * ) i1, i2, i3, foo
     DO gpt2 = 1, psi_iks % grid % ngpts
        IF ( psi_iks % grid % g_hat ( 1, gpt2 ) == i1 .AND. &
             psi_iks % grid % g_hat ( 2, gpt2 ) == i2 .AND. &
             psi_iks % grid % g_hat ( 3, gpt2 ) == i3 ) THEN
           psi_iks % cc ( gpt2 ) = foo
           nhits = nhits + 1
           EXIT
        END IF
     END DO
  END DO
  CLOSE ( UNIT = 99 )
  WRITE ( 6, * ) "DDD: #hits = ", nhits
#endif
  
  IF ( PRESENT ( psi_iks_r ) ) THEN
     cfft % cr3d = vloc % pw % cr3d * psi_iks_r % pw % cc3d
  ELSE
     CALL fft_wrap ( psi_iks, cfft )
     
     cfft % cc3d = cfft % cc3d * vloc % pw % cr3d
  END IF
  
  CALL fft_wrap ( cfft, h_psi )
  
  !! This must be done after the local potential
  CALL pw_kinetic_evaluate ( psi_iks, xkpt, ekinetic, h_psi )
  
  CALL pw_deallocate ( cfft )
  
END SUBROUTINE pw_h_psi

!******************************************************************************

END MODULE pws
