!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pw_grids [1.0] *
!!
!!   NAME
!!     pw_grids
!!
!!   FUNCTION
!!     This module defines the grid data type and some basic operations on it
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pw_grids
  
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : twopi
  USE pw_grid_types, ONLY : pw_grid_type, map_pn, HALFSPACE, FULLSPACE
  USE simulation_cell, ONLY : cell_type
  USE stop_program, ONLY : stop_prg
  USE util, ONLY : dotprod_3d, sort
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_grid_setup, pw_grid_change, pw_find_cutoff
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE pw_grid_setup ( cell, pw_grid, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: cutoff
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( cutoff ) &
       .AND. SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) THEN
     CALL pw_grid_find_bounds ( pw_grid, cell, cutoff )
  ELSE
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) THEN
        CALL stop_prg ( "grid_setup", &
             "provide initial values for bounds yourself" )
     END IF
  END IF
  
  pw_grid % npts ( : ) = &
       pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1
  
  ! Find the number of grid points
  CALL pw_grid_count ( cell, pw_grid, cutoff )
  
  pw_grid % ngpts = PRODUCT ( pw_grid % npts )
  
! IF ( PRESENT ( cutoff ) ) THEN
!    write ( 6, '( A, I8 )' ) " PW_grid_setup| #gpts = ", pw_grid % ngpts_cut
! ELSE
!    write ( 6, '( A, I8 )' ) " PW_grid_setup| #gpts = ", pw_grid % ngpts
! END IF
  
  ! Allocate the grid fields
  CALL pw_grid_allocate ( pw_grid )
  
  ! Fill in the grid structure
  CALL pw_grid_assign ( cell, pw_grid, cutoff )
  
  CALL pw_grid_sort ( pw_grid )
  
  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts, dbl )
  pw_grid % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ) ** 2, 1 ) ) &
       / REAL ( pw_grid % npts ( : ), dbl )
  
END SUBROUTINE pw_grid_setup

!******************************************************************************

SUBROUTINE pw_grid_count ( cell, pw_grid, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: cutoff
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER :: gpt, l, m, n, n_upperlimit
  REAL ( dbl ) :: length_x, length_y, length_z, length
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  bounds => pw_grid % bounds
  h_inv => cell % h_inv
  
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     n_upperlimit = 0
  ELSE IF ( pw_grid % grid_span == FULLSPACE ) THEN
     n_upperlimit = bounds ( 2, 3 )
  ELSE
     CALL stop_prg ( "pw_grid_count", "no type set for the grid" )
  END IF
  
! finds valid g-points within grid
  gpt = 0
  DO n = bounds ( 1, 3 ), n_upperlimit
     DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
        DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
           IF ( pw_grid % grid_span == HALFSPACE .AND. n == 0 ) THEN
              IF ( ( m == 0 .AND. l > 0 ) .OR. ( m > 0 ) ) CYCLE
           END IF
           
           IF ( PRESENT ( cutoff ) ) THEN
              length_x &
                   = REAL(l,dbl) * h_inv(1,1) &
                   + REAL(m,dbl) * h_inv(2,1) &
                   + REAL(n,dbl) * h_inv(3,1)
              length_y &
                   = REAL(l,dbl) * h_inv(1,2) &
                   + REAL(m,dbl) * h_inv(2,2) &
                   + REAL(n,dbl) * h_inv(3,2)
              length_z &
                   = REAL(l,dbl) * h_inv(1,3) &
                   + REAL(m,dbl) * h_inv(2,3) &
                   + REAL(n,dbl) * h_inv(3,3)
              
              length = SQRT ( length_x ** 2 + length_y ** 2 + length_z ** 2 )
              length = twopi * length
              
              IF ( length <= cutoff ) THEN
                 gpt = gpt + 1
              END IF
           ELSE
              gpt = gpt + 1
           END IF
           
        END DO
     END DO
  END DO
  
! number of g-vectors for grid
  pw_grid % ngpts_cut = gpt

END SUBROUTINE pw_grid_count

!******************************************************************************

SUBROUTINE pw_grid_assign ( cell, pw_grid, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: cutoff
  TYPE ( cell_type ), INTENT ( IN ), OPTIONAL, TARGET :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER :: gpt,l,m,n, n_upperlimit
  REAL ( dbl ) :: length_x, length_y, length_z, length
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  bounds => pw_grid % bounds
  h_inv => cell % h_inv
  
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     n_upperlimit = 0
  ELSE IF ( pw_grid % grid_span == FULLSPACE ) THEN
     n_upperlimit = bounds ( 2, 3 )
  ELSE
     CALL stop_prg ( "pw_grid_assign", "no type set for the grid" )
  END IF
  
! finds valid g-points within grid
  gpt = 0
  DO n = bounds ( 1, 3 ), n_upperlimit
     DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
        DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
           IF ( pw_grid % grid_span == HALFSPACE .AND. n == 0 ) THEN
              IF ( ( m == 0 .AND. l > 0 ) .OR. ( m > 0 ) ) CYCLE
           END IF
           
           length_x &
                = REAL(l,dbl) * h_inv(1,1) &
                + REAL(m,dbl) * h_inv(2,1) &
                + REAL(n,dbl) * h_inv(3,1)
           length_y &
                = REAL(l,dbl) * h_inv(1,2) &
                + REAL(m,dbl) * h_inv(2,2) &
                + REAL(n,dbl) * h_inv(3,2)
           length_z &
                = REAL(l,dbl) * h_inv(1,3) &
                + REAL(m,dbl) * h_inv(2,3) &
                + REAL(n,dbl) * h_inv(3,3)
           
           length = SQRT ( length_x ** 2 + length_y ** 2 + length_z ** 2 )
           length = twopi * length
           
           IF ( PRESENT ( cutoff ) ) THEN
              IF ( length <= cutoff ) THEN
                 gpt = gpt + 1
                 pw_grid % g ( 1, gpt ) = twopi * length_x
                 pw_grid % g ( 2, gpt ) = twopi * length_y
                 pw_grid % g ( 3, gpt ) = twopi * length_z
                 pw_grid % gsq ( gpt ) = length ** 2
                 pw_grid % g_hat ( 1, gpt ) = l
                 pw_grid % g_hat ( 2, gpt ) = m
                 pw_grid % g_hat ( 3, gpt ) = n
              END IF
           ELSE
              gpt = gpt + 1
              pw_grid % g ( 1, gpt ) = twopi * length_x
              pw_grid % g ( 2, gpt ) = twopi * length_y
              pw_grid % g ( 3, gpt ) = twopi * length_z
              pw_grid % gsq ( gpt ) = length ** 2
              pw_grid % g_hat ( 1, gpt ) = l
              pw_grid % g_hat ( 2, gpt ) = m
              pw_grid % g_hat ( 3, gpt ) = n
           END IF
           
        END DO
     END DO
  END DO
  
! Check the number of g-vectors for grid
  IF ( pw_grid % ngpts_cut /= gpt ) THEN
     CALL stop_prg ( "pw_grid_assign", "error re-counting the vectors" )
  END IF
  
  pw_grid % have_g0 = .FALSE.
  pw_grid % first_gne0 = 2
  DO gpt = 1, pw_grid % ngpts_cut
     IF ( ALL ( pw_grid % g_hat ( :, gpt ) == 0 ) ) THEN
        pw_grid % have_g0 = .TRUE.
        pw_grid % first_gne0 = 1
        EXIT
     END IF
  END DO
  
  CALL pw_grid_set_maps ( pw_grid % g_hat, pw_grid % ngpts_cut, &
       pw_grid % mapl, pw_grid % mapm, pw_grid % mapn, pw_grid % npts )
  
END SUBROUTINE pw_grid_assign

!******************************************************************************

SUBROUTINE pw_grid_set_maps ( g_hat, ngpts_cut, mapl, mapm, mapn, npts )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: ngpts_cut
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: g_hat
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapl
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapm
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapn
  
! Locals
  INTEGER :: l, m, n, gpt
  
!------------------------------------------------------------------------------
  
  DO gpt = 1, ngpts_cut
     l = g_hat ( 1, gpt )
     m = g_hat ( 2, gpt )
     n = g_hat ( 3, gpt )
     IF ( l < 0 ) THEN
        mapl % pos ( g_hat ( 1, gpt ) ) = g_hat ( 1, gpt ) + npts ( 1 )
     ELSE
        mapl % pos ( g_hat ( 1, gpt ) ) = g_hat ( 1, gpt )
     END IF
     IF ( m < 0 ) THEN
        mapm % pos ( g_hat ( 2, gpt ) ) = g_hat ( 2, gpt ) + npts ( 2 )
     ELSE
        mapm % pos ( g_hat ( 2, gpt ) ) = g_hat ( 2, gpt )
     END IF
     IF ( n < 0 ) THEN
        mapn % pos ( g_hat ( 3, gpt ) ) = g_hat ( 3, gpt ) + npts ( 3 )
     ELSE
        mapn % pos ( g_hat ( 3, gpt ) ) = g_hat ( 3, gpt )
     END IF
     
! Generating the maps to the full 3-d space
     
     IF ( l <= 0 ) THEN
        mapl % neg ( g_hat ( 1, gpt ) ) = - g_hat ( 1, gpt )
     ELSE
        mapl % neg ( g_hat ( 1, gpt ) ) = npts ( 1 ) - g_hat ( 1, gpt )
     END IF
     IF ( m <= 0 ) THEN
        mapm % neg ( g_hat ( 2, gpt ) ) = - g_hat ( 2, gpt )
     ELSE
        mapm % neg ( g_hat ( 2, gpt ) ) = npts ( 2 ) - g_hat ( 2, gpt )
     END IF
     IF ( n <= 0 ) THEN
        mapn % neg ( g_hat ( 3, gpt ) ) = - g_hat ( 3, gpt )
     ELSE
        mapn % neg ( g_hat ( 3, gpt ) ) = npts ( 3 ) - g_hat ( 3, gpt )
     END IF
  END DO
  
END SUBROUTINE pw_grid_set_maps

!******************************************************************************

SUBROUTINE pw_grid_find_bounds ( pw_grid, cell, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  REAL ( dbl ), INTENT ( IN ) :: cutoff
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: cell
  
! Locals
  INTEGER :: l, m, n, test_max = 50
  REAL ( dbl ) :: length_x, length_y, length_z, length
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  bounds => pw_grid % bounds
  h_inv => cell % h_inv
  
  bounds ( :, : ) = 0
  
! finds valid g-points within grid
  DO n = 0, +test_max
     DO m = -test_max, +test_max
        DO l = -test_max, +test_max
           length_x &
                = REAL(l,dbl) * h_inv(1,1) &
                + REAL(m,dbl) * h_inv(2,1) &
                + REAL(n,dbl) * h_inv(3,1)
           length_y &
                = REAL(l,dbl) * h_inv(1,2) &
                + REAL(m,dbl) * h_inv(2,2) &
                + REAL(n,dbl) * h_inv(3,2)
           length_z &
                = REAL(l,dbl) * h_inv(1,3) &
                + REAL(m,dbl) * h_inv(2,3) &
                + REAL(n,dbl) * h_inv(3,3)
           
           length = SQRT ( length_x ** 2 + length_y ** 2 + length_z ** 2 )
           length = twopi * length
           
           IF ( length <= cutoff ) THEN
              bounds ( 2, 1 ) = MAX ( bounds ( 2, 1 ), l )
              bounds ( 2, 2 ) = MAX ( bounds ( 2, 2 ), m )
              bounds ( 2, 3 ) = MAX ( bounds ( 2, 3 ), n )
           END IF
           
        END DO
     END DO
  END DO
  
  IF ( bounds ( 2, 1 ) == test_max .OR. &
       bounds ( 2, 2 ) == test_max .OR. &
       bounds ( 2, 3 ) == test_max ) THEN
     WRITE ( 6, * ) "grid_find_bounds: WARNING"
     WRITE ( 6, * ) "  initial bounds maybe too small:"
     WRITE ( 6, * ) bounds ( 2, : ), test_max
  END IF
  
  write(6,*) "apsi: automatic grid bounds ", bounds(2,:)
  
  bounds ( 1, : ) = - bounds ( 2, : )
  
END SUBROUTINE pw_grid_find_bounds

!****************************************************************************

SUBROUTINE pw_grid_allocate ( pw_grid )
  IMPLICIT NONE
  
! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER, POINTER :: ngpts_cut
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  ngpts_cut => pw_grid % ngpts_cut
  bounds => pw_grid % bounds
  
  ALLOCATE ( pw_grid % g ( 3, ngpts_cut ) )
  ALLOCATE ( pw_grid % gsq ( ngpts_cut ) )
  ALLOCATE ( pw_grid % g_hat ( 3, ngpts_cut ) )
  ALLOCATE ( pw_grid % mapl % pos ( bounds ( 1, 1 ):bounds ( 2, 1 ) ) )
  ALLOCATE ( pw_grid % mapl % neg ( bounds ( 1, 1 ):bounds ( 2, 1 ) ) )
  ALLOCATE ( pw_grid % mapm % pos ( bounds ( 1, 2 ):bounds ( 2, 2 ) ) )
  ALLOCATE ( pw_grid % mapm % neg ( bounds ( 1, 2 ):bounds ( 2, 2 ) ) )
  ALLOCATE ( pw_grid % mapn % pos ( bounds ( 1, 3 ):bounds ( 2, 3 ) ) )
  ALLOCATE ( pw_grid % mapn % neg ( bounds ( 1, 3 ):bounds ( 2, 3 ) ) )
  
END SUBROUTINE pw_grid_allocate

!****************************************************************************

SUBROUTINE pw_grid_sort ( pw_grid )
  
  IMPLICIT NONE
  
! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER :: gpt, gpt_curr, gpt_tmp, int_tmp ( 3 )
  INTEGER :: idx ( SIZE ( pw_grid % gsq ) )
  REAL ( dbl ) :: real_tmp ( 3 )
  
!------------------------------------------------------------------------------
  
  IF ( .NOT. ASSOCIATED ( pw_grid % gsq ) ) THEN
     CALL stop_prg ( "pw_grid_sort", &
          "works currently only if lenghts are filled" )
  END IF
  
  CALL sort ( pw_grid % gsq, pw_grid % ngpts_cut, idx )
  
  DO gpt = 1, pw_grid % ngpts_cut - 1
     gpt_curr = gpt
     
     DO
        IF ( idx ( gpt_curr ) == gpt ) EXIT
        
        real_tmp ( : ) = pw_grid % g ( :, gpt_curr )
        pw_grid % g ( :, gpt_curr ) = pw_grid % g ( :, idx ( gpt_curr ) )
        pw_grid % g ( :, idx ( gpt_curr ) ) = real_tmp ( : )
        
        int_tmp ( : ) = pw_grid % g_hat ( :, gpt_curr )
        pw_grid % g_hat ( :, gpt_curr ) &
             = pw_grid % g_hat ( :, idx ( gpt_curr ) )
        pw_grid % g_hat ( :, idx ( gpt_curr ) ) = int_tmp ( : )
        
        gpt_tmp = gpt_curr
        gpt_curr = idx ( gpt_curr )
        idx ( gpt_tmp ) = gpt_tmp
        
        IF ( idx ( gpt_curr ) == gpt ) then
           idx ( gpt_curr ) = gpt_curr
           EXIT
        END IF
     END DO
  END DO
  
END SUBROUTINE pw_grid_sort

!****************************************************************************

SUBROUTINE pw_grid_change ( cell, pw_grid )
  IMPLICIT NONE
  
! Argument
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER :: gpt, l, m, n
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: g
  
!------------------------------------------------------------------------------
  
  h_inv => cell % h_inv
  g => pw_grid % g
  
  DO gpt = 1, pw_grid % ngpts_cut
     l = pw_grid % g_hat ( 1, gpt )
     m = pw_grid % g_hat ( 2, gpt )
     n = pw_grid % g_hat ( 3, gpt )
     
     g ( 1, gpt ) &
          = REAL(l,dbl) * h_inv(1,1) &
          + REAL(m,dbl) * h_inv(2,1) &
          + REAL(n,dbl) * h_inv(3,1)
     g ( 2, gpt ) &
          = REAL(l,dbl) * h_inv(1,2) &
          + REAL(m,dbl) * h_inv(2,2) &
          + REAL(n,dbl) * h_inv(3,2)
     g ( 3, gpt ) &
          = REAL(l,dbl) * h_inv(1,3) &
          + REAL(m,dbl) * h_inv(2,3) &
          + REAL(n,dbl) * h_inv(3,3)
     
     IF ( ASSOCIATED ( pw_grid % gsq ) ) &
          pw_grid % gsq ( gpt ) = SUM ( g ( :, gpt ) ** 2 )
  END DO
  
END SUBROUTINE pw_grid_change

!****************************************************************************

SUBROUTINE pw_find_cutoff ( npts, box, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( OUT ) :: cutoff
  
! Locals
  INTEGER :: i
  REAL ( dbl ) :: gdum ( 3 ), length, gcut
  
!------------------------------------------------------------------------------
  
! compute 2*pi*h_inv^t*g  where g = (npts[1],0,0)
  gdum ( : ) = twopi * box % h_inv ( 1, : ) &
       * REAL ( npts ( 1 ), dbl ) / 2.0_dbl
!*apsi:old       * REAL ( npts ( 1 ) - 1, dbl ) / 2.0_dbl
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = length
  
! compute 2*pi*h_inv^t*g  where g = (0,npts[2],0)
  gdum ( : ) = twopi * box % h_inv ( 2, : ) &
       * REAL ( npts ( 2 ), dbl ) / 2.0_dbl
!*apsi:old       * REAL ( npts ( 2 ) - 1, dbl ) / 2.0_dbl
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )
  
! compute 2*pi*h_inv^t*g  where g = (0,0,npts[3])
  gdum ( : ) = twopi * box % h_inv ( 3, : ) &
       * REAL ( npts ( 3 ), dbl ) / 2.0_dbl
!*apsi:old       * REAL ( npts ( 3 ) - 1, dbl ) / 2.0_dbl
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )
  
! Because there is a HUGE OZONE hole over FINLAND that
! effects SOME  co-workers...We subtract epsilon.
  
  cutoff = gcut - 1.0E-5_dbl
!*apsi:old  cutoff = gcut + 1e-5_dbl
  
END SUBROUTINE pw_find_cutoff

!******************************************************************************

END MODULE pw_grids
