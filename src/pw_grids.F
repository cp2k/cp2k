!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pw_grids [1.1] *
!!
!!   NAME
!!     pw_grids
!!
!!   FUNCTION
!!     This module defines the grid data type and some basic operations on it
!!
!!   AUTHOR
!!     apsi
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Adapted for parallel use
!!
!!   NOTES
!!     pw_grid_setup  : main routine to set up a grid
!!          input: cell (the box for the grid)
!!                 pw_grid (the grid; pw_grid%grid_span has to be set)
!!                 cutoff (optional, if not given pw_grid%bounds has to be set)
!!                 pe_group (optional, if not given we have a local grid)
!!
!!                 if no cutoff or a negative cutoff is given, all g-vectors
!!                 in the box are included (no spherical cutoff)
!!          output: pw_grid
!!
!!     pw_grid_change : updates g-vectors after a change of the box
!!
!!     pw_find_cutoff : Calculates the cutoff for given box and points
!!
!!   SOURCE
!******************************************************************************

MODULE pw_grids
  
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : twopi
  USE message_passing, ONLY : mp_environ, mp_sum, mp_cart_create, mp_group, &
                              mp_max, mp_min
  USE pw_grid_types, ONLY : pw_grid_type, map_pn, HALFSPACE, FULLSPACE
  USE output_utilities, ONLY : print_warning
  USE simulation_cell, ONLY : cell_type
  USE stop_program, ONLY : stop_prg, stop_memory
  USE util, ONLY : matvec_3x3, dotprod_3d, sort, get_limit
  USE fft_tools, ONLY : fft_radix_operations, FFT_RADIX_NEXT
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_grid_setup, pw_grid_change, pw_find_cutoff
  
  INTEGER :: grid_tag = 0
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pw_grids/pw_grid_setup [1.1] *
!!
!!   NAME
!!     pw_grid_setup
!!
!!   FUNCTION
!!     Defines the PW grid
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Adapted for parallel use
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_setup ( cell, pw_grid, cutoff, pe_group, info )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: cutoff
  INTEGER, OPTIONAL :: pe_group          ! parent group for this grid
  INTEGER, OPTIONAL :: info              ! output unit for information on grid

! Locals
  REAL ( dbl ) :: ecut
  INTEGER, DIMENSION ( :, : ), ALLOCATABLE :: yz_mask
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: group_list
  INTEGER :: allocstat, group_size, my_pos, i, n(3)
  REAL ( dbl ) :: rv(3,3)
  INTEGER :: aver, amax, amin 
  
!------------------------------------------------------------------------------
  ! assign a unique tag to this grid
  grid_tag = grid_tag + 1
  pw_grid % identifier = grid_tag

  ! parallel info
  IF (PRESENT ( pe_group ) ) THEN
    CALL mp_environ ( group_size, my_pos, pe_group )
    IF ( group_size > 1 ) THEN
      ALLOCATE ( group_list ( group_size ), STAT = allocstat )
      IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_setup", "group_list", &
          group_size )
      group_list = (/( I, I = 0, group_size - 1 )/)
      CALL mp_group ( group_list, group_size, pe_group, &
                      pw_grid % para % group )
      CALL mp_environ ( pw_grid % para % group_size, &
                        pw_grid % para % my_pos, &
                        pw_grid % para % group )
      pw_grid % para % group_head_id = 0
      pw_grid % para % group_head = &
                ( pw_grid % para % group_head_id == pw_grid % para % my_pos )
      pw_grid % para % mode = 1
      DEALLOCATE ( group_list, STAT = allocstat )
      IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_setup", "group_list" )
    ELSE
      pw_grid % para % mode = 0
    END IF
  ELSE
    pw_grid % para % mode = 0
  END IF

  IF ( PRESENT ( cutoff ) ) THEN
     ecut = cutoff
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) &
       CALL pw_grid_find_bounds ( pw_grid % bounds, cell % h_inv, ecut )
     pw_grid % spherical = .TRUE.
     pw_grid % cutoff = 0.5_dbl * ABS ( ecut )
  ELSE
     ecut = 1.e10_dbl              ! all g-vectors in the box will be included
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) THEN
        CALL stop_prg ( "grid_setup", &
             "provide initial values for bounds yourself" )
     END IF
     pw_grid % spherical = .FALSE.
     pw_grid % cutoff = 0._dbl
!    ! without a cutoff and HALFSPACE we have to be sure that there is
!    ! a negative counterpart to every g vector (-> odd number of grid points)
!    IF ( pw_grid % grid_span == HALFSPACE ) THEN
!      n = pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1
!      IF ( MOD ( n ( 1 ), 2 ) == 0 .OR. &
!           MOD ( n ( 2 ), 2 ) == 0 .OR. &
!           MOD ( n ( 3 ), 2 ) == 0 ) &
!         CALL stop_prg ( "grid_setup", &
!            "No spherical cutoff and HALFSPACE option", &
!            "Number of grid points in each direction", "has to be odd" )
!    END IF
  END IF
  
  pw_grid % npts ( : ) = &
       pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1

  n ( : )  = pw_grid % npts ( : )

  ! Find the number of grid points
  ! yz_mask counts the number of g-vectors orthogonal to the yz plane
  ! the indices in yz_mask are from -n/2 .. n/2 shifted by n/2 + 1
  ! these are not mapped indices !
  ALLOCATE ( yz_mask ( n(2), n(3) ), STAT = allocstat )
  IF ( allocstat /= 0 ) call stop_memory ( "pw_grid_setup", "yz_mask", 0 )
  CALL pw_grid_count ( cell % h_inv, pw_grid, ecut, yz_mask )
  
  ! Distribute grid 
  CALL pw_grid_distribute ( pw_grid, yz_mask )

  ! Allocate the grid fields
  CALL pw_grid_allocate ( pw_grid, pw_grid % ngpts_cut_local, &
                          pw_grid % bounds )
  
  ! Fill in the grid structure
  CALL pw_grid_assign ( cell % h_inv, pw_grid, ecut )
  
  DEALLOCATE ( yz_mask )
  IF ( allocstat /= 0 ) call stop_memory ( "pw_grid_setup", "yz_mask", 0 )

  CALL pw_grid_sort ( pw_grid )
  
  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts, dbl )
  pw_grid % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ) ** 2, 1 ) ) &
       / REAL ( pw_grid % npts ( : ), dbl )

!
! Output: All the information of this grid type
!

  IF ( PRESENT ( info ) ) THEN
    IF ( pw_grid % para % mode == 0 ) THEN
      WRITE ( info, '(/,A,T71,I10)' ) &
        " PW_GRID: Information for grid number ", pw_grid % identifier
      IF ( pw_grid % spherical ) THEN
        WRITE ( info, '(A,T35,A,T66,F10.1,T77,A)' ) " PW_GRID: ", &
          "spherical cutoff", pw_grid % cutoff, "a.u."
        WRITE ( info, '(A,T71,I10)' ) " PW_GRID: Grid points within cutoff", &
          pw_grid % ngpts_cut
      END IF
      DO i = 1, 3
        WRITE ( info, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " PW_GRID:   Bounds ", &
         i, pw_grid % bounds ( 1, I ), pw_grid % bounds ( 2, I ), &
         "Points:",pw_grid % npts ( I )
      END DO
      WRITE ( info, '(A,G12.4,T50,A,T67,F14.4)' ) &
         " PW_GRID: Volume element (a.u.^3)", &
         pw_grid % dvol," Volume (a.u.^3) :",pw_grid % vol
      IF ( pw_grid % grid_span == HALFSPACE ) THEN
        WRITE ( info, '(A,T72,A)' ) " PW_GRID: Grid span","HALFSPACE"
      ELSE
        WRITE ( info, '(A,T72,A)' ) " PW_GRID: Grid span","FULLSPACE"
      END IF
      WRITE ( info, '(/)' ) 
    ELSE
      IF ( pw_grid % para % group_head ) THEN
        WRITE ( info, '(/,A,T71,I10)' ) &
          " PW_GRID: Information for grid number ", pw_grid % identifier
        WRITE ( info, '(A,T60,I10,A)' ) &
          " PW_GRID: Grid distributed over ", pw_grid % para % group_size, &
          " processors"
        WRITE ( info, '(A,T71,2I5)' ) &
          " PW_GRID: Real space group dimensions ", pw_grid % para % rs_dims
        IF ( pw_grid % spherical ) THEN
          WRITE ( info, '(A,T35,A,T66,F10.1,T77,A)' ) " PW_GRID: ", &
            "spherical cutoff", pw_grid % cutoff, "a.u."
          WRITE ( info, '(A,T71,I10)' ) " PW_GRID: Grid points within cutoff", &
            pw_grid % ngpts_cut
        END IF
        DO i = 1, 3
          WRITE ( info, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " PW_GRID:   Bounds ", &
           i, pw_grid % bounds ( 1, I ), pw_grid % bounds ( 2, I ), &
           "Points:",pw_grid % npts ( I )
        END DO
        WRITE ( info, '(A,G12.4,T50,A,T67,F14.4)' ) &
           " PW_GRID: Volume element (a.u.^3)", &
           pw_grid % dvol," Volume (a.u.^3) :",pw_grid % vol
        IF ( pw_grid % grid_span == HALFSPACE ) THEN
          WRITE ( info, '(A,T72,A)' ) " PW_GRID: Grid span","HALFSPACE"
        ELSE
          WRITE ( info, '(A,T72,A)' ) " PW_GRID: Grid span","FULLSPACE"
        END IF
      END IF
      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_sum ( n(1:2), pw_grid % para % group )
      n ( 3 ) = SUM ( pw_grid % para % nyzray )
      rv ( :, 1 ) = REAL ( n, dbl ) / REAL ( group_size, dbl )
      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_max ( n(1:2), pw_grid % para % group )
      n ( 3 ) = MAXVAL ( pw_grid % para % nyzray )
      rv ( :, 2 ) = REAL ( n, dbl )
      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_min ( n(1:2), pw_grid % para % group )
      n ( 3 ) = MInVAL ( pw_grid % para % nyzray )
      rv ( :, 3 ) = REAL ( n, dbl )
      IF ( pw_grid % para % group_head ) THEN
        WRITE ( info, '(A,T48,A)' ) " PW_GRID:   Distribution", &
             "  Average         Max         Min"
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   G-Vectors", &
             rv ( 1, 1 ), NINT ( rv ( 1, 2 ) ), NINT ( rv ( 1, 3 ) )
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   G-Rays   ", &
             rv ( 3, 1 ), NINT ( rv ( 3, 2 ) ), NINT ( rv ( 3, 3 ) )
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   Real Space Points", &
             rv ( 2, 1 ), NINT ( rv ( 2, 2 ) ), NINT ( rv ( 2, 3 ) )
      END IF
    END IF
  END IF

END SUBROUTINE pw_grid_setup

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_distribute [1.0] *
!!
!!   NAME
!!     pw_grid_distribute
!!
!!   FUNCTION
!!     Distribute grids in real and Fourier Space to the processors in group
!!
!!   AUTHOR
!!     JGH (22-12-2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_distribute ( pw_grid, yz_mask )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  INTEGER, DIMENSION ( :, : ), INTENT ( INOUT ) :: yz_mask

! Locals
  INTEGER :: gmin, gmax, nx, ny, nz, np, ns, i, j, k, m, n, npx_real, npz_g
  INTEGER :: lo(2), ip, ierr, gneg, lby, lbz
  REAL ( dbl ) :: tfun, tt
  
!------------------------------------------------------------------------------

  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )

  pw_grid % ngpts = PRODUCT ( pw_grid % npts )

  IF ( pw_grid % para % mode == 0) THEN

    pw_grid % bounds_local = pw_grid % bounds
    pw_grid % npts_local = pw_grid % npts
    pw_grid % ngpts_cut_local = pw_grid % ngpts_cut

  ELSE

!..find the real space distribution
    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )

    np = pw_grid % para % group_size

    ns = INT ( SQRT ( REAL ( np, dbl ) ) )

    tfun = 1.e20_dbl
    DO i = ns, 1, -1
      IF ( MOD ( np, i ) == 0 ) THEN
        j = np/i
        k = nx*ny - i*(nx/i) * j*(ny/j)
        IF ( k == 0 ) THEN
          tfun = 0._dbl
          m = i
        ELSE
          tt = 1._dbl - REAL ( k, dbl ) / REAL ( i*j, dbl )
          IF ( tt < tfun ) THEN
            tfun = tt
            m = i
          END IF
        END IF
      END IF
    END DO
    k = ( nx - np*(nx/np)) * ny
    IF ( k == 0 ) THEN
      m = np
    ELSE
      tt = 1._dbl - REAL ( k, dbl ) / REAL ( np, dbl )
      IF ( tt < tfun ) m = np
    END IF
!..create group for real space distribution
    pw_grid % para % rs_dims ( 1 ) = m
    pw_grid % para % rs_dims ( 2 ) = np/m
    CALL mp_cart_create ( pw_grid % para % group, 2, &
                          pw_grid % para % rs_dims, &
                          pw_grid % para % rs_pos, &
                          pw_grid % para % rs_group )
    lo = get_limit ( nx, pw_grid % para % rs_dims ( 1 ), &
                     pw_grid % para % rs_pos ( 1 ) )
    pw_grid % bounds_local ( :, 1 ) = lo + pw_grid % bounds ( 1, 1 ) - 1
    lo = get_limit ( nx, pw_grid % para % rs_dims ( 2 ), &
                     pw_grid % para % rs_pos ( 2 ) )
    pw_grid % bounds_local ( :, 2 ) = lo + pw_grid % bounds ( 1, 2 ) - 1
    pw_grid % bounds_local ( :, 3 ) = pw_grid % bounds ( :, 3 )
    pw_grid % npts_local ( : ) = pw_grid % bounds_local ( 2, : ) &
                                 - pw_grid % bounds_local ( 1, : ) + 1

!..find the g space distribution
    pw_grid % ngpts_cut_local = 0

    ALLOCATE ( pw_grid % para % yzq ( ny, nz ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
       "pw_grid % para % yzq", ny*nz )
    ALLOCATE ( pw_grid % para % nyzray ( 0: np-1 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
          "pw_grid % para % nyzray", np )

    ! Round Robin distribution 
    ! Processors 0 .. NP-1, NP-1 .. 0  get the largest remaining batch
    ! of g vectors in turn

    pw_grid % para % yzq = -1
    DO
      DO i = 1, 2*np
        IF ( i > np ) THEN
          ip = 2*np - i
        ELSE
          ip = i - 1
        ENDIF
        lo = MAXLOC ( yz_mask )
        gmax = yz_mask ( lo(1), lo(2) )
        IF ( gmax == 0 ) EXIT
        yz_mask ( lo(1), lo(2) ) = 0
        IF ( ip == pw_grid % para % my_pos ) THEN
          pw_grid % ngpts_cut_local = pw_grid % ngpts_cut_local + gmax
        END IF
        pw_grid % para % yzq ( lo(1), lo(2) ) = ip
        IF ( pw_grid % grid_span == HALFSPACE ) THEN
          m = -lo(1) - 2*lby + 2
          n = -lo(2) - 2*lbz + 2
          pw_grid % para % yzq ( m, n ) = ip
          yz_mask ( m, n ) = 0
        END IF
      END DO
      IF ( gmax == 0 ) EXIT
    END DO

    ! Count the total number of rays on each processor
    pw_grid % para % nyzray = 0
    DO i = 1, nz
      DO j = 1, ny
        ip = pw_grid % para % yzq ( j, i )
        IF ( ip >= 0 ) pw_grid % para % nyzray ( ip ) = &
                       pw_grid % para % nyzray ( ip ) + 1
      END DO
    END DO

    ! Allocate mapping array (y:z, nray, nproc)
    ns = MAXVAL ( pw_grid % para % nyzray ( 0: np-1 ) )
    ALLOCATE ( pw_grid % para % yzp ( 2, ns, 0: np-1 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
          "pw_grid % para % yzp", 2*ns*np )

    ! Fill mapping array, recalculate nyzray for convenience
    pw_grid % para % nyzray = 0
    DO i = 1, nz
      DO j = 1, ny
        ip = pw_grid % para % yzq ( j, i )
        IF ( ip >= 0 ) THEN
          pw_grid % para % nyzray ( ip ) = &
                       pw_grid % para % nyzray ( ip ) + 1
          ns = pw_grid % para % nyzray ( ip )
          pw_grid % para % yzp ( 1, ns, ip ) = j + lby - 1
          pw_grid % para % yzp ( 2, ns, ip ) = i + lbz - 1
          IF ( ip == pw_grid % para % my_pos ) THEN
            pw_grid % para % yzq ( j, i ) = ns
          ELSE
            pw_grid % para % yzq ( j, i ) = -1
          END IF
        ELSE
          pw_grid % para % yzq ( j, i ) = -2
        END IF
      END DO
    END DO

  END IF

  pw_grid % ngpts_local = PRODUCT ( pw_grid % npts_local )

END SUBROUTINE pw_grid_distribute

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_count [1.1] *
!!
!!   NAME
!!     pw_grid_count
!!
!!   FUNCTION
!!     Count total number of g vectors
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (22-12-2000) : Adapted for parallel use
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_count ( h_inv, pw_grid, cutoff, yz_mask )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: h_inv
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  REAL ( dbl ), INTENT ( IN ) :: cutoff
  INTEGER, DIMENSION ( :, : ), INTENT ( OUT ) :: yz_mask
  
! Locals
  INTEGER :: gpt, l, m, n, n_upperlimit, nlim ( 2 ), mm, nn
  REAL ( dbl ) :: gmat ( 3, 3 ), gi ( 3 )
  REAL ( dbl ) :: length
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  bounds => pw_grid % bounds
  
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     n_upperlimit = 0
  ELSE IF ( pw_grid % grid_span == FULLSPACE ) THEN
     n_upperlimit = bounds ( 2, 3 )
  ELSE
     CALL stop_prg ( "pw_grid_count", "no type set for the grid" )
  END IF
  
! finds valid g-points within grid
  gmat = MATMUL ( h_inv, TRANSPOSE ( h_inv ) )
  gpt = 0
  IF ( pw_grid % para % mode == 0 ) THEN
    nlim ( 1 ) = bounds ( 1, 3 )
    nlim ( 2 ) = n_upperlimit
  ELSE IF ( pw_grid % para % mode == 1 ) THEN
    n = n_upperlimit - bounds ( 1, 3 ) + 1
    nlim = get_limit ( n, pw_grid % para % group_size, pw_grid % para % my_pos )
    nlim = nlim + bounds ( 1, 3 ) - 1
  ELSE
     CALL stop_prg ( "pw_grid_count", "para % mode not specified" )
  END IF

  yz_mask = 0
  DO n = nlim ( 1 ), nlim ( 2 )
     gi ( 3 ) = REAL(n,dbl)
     nn = n - bounds ( 1, 3) + 1
     DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
        gi ( 2 ) = REAL(m,dbl)
        mm = m - bounds ( 1, 2) + 1
        DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
           IF ( pw_grid % grid_span == HALFSPACE .AND. n == 0 ) THEN
              IF ( ( m == 0 .AND. l > 0 ) .OR. ( m > 0 ) ) CYCLE
           END IF

           gi ( 1 ) = REAL(l,dbl)
           length = twopi * twopi * dotprod_3d ( matvec_3x3 ( gmat, gi ), gi )
           IF ( length <= cutoff * cutoff ) THEN
             gpt = gpt + 1
             yz_mask ( mm, nn ) = yz_mask ( mm, nn ) + 1
           END IF
           
        END DO
     END DO
  END DO
  
! number of g-vectors for grid
  IF ( pw_grid % para % mode == 1 ) THEN
    CALL mp_sum ( gpt, pw_grid % para % group )
    CALL mp_sum ( yz_mask, pw_grid % para % group )
  ENDIF
  pw_grid % ngpts_cut = gpt

END SUBROUTINE pw_grid_count

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_assign [1.1] *
!!
!!   NAME
!!     pw_grid_assign
!!
!!   FUNCTION
!!     Setup maps from 1d to 3d space
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (29-12-2000) : Adapted for parallel use
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_assign ( h_inv, pw_grid, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: h_inv
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  REAL ( dbl ), INTENT ( IN ) :: cutoff
  
! Locals
  INTEGER :: gpt, i, l, m, n, n_upperlimit, ip
  REAL ( dbl ) :: length_x, length_y, length_z, length
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  bounds => pw_grid % bounds
  
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     n_upperlimit = 0
  ELSE IF ( pw_grid % grid_span == FULLSPACE ) THEN
     n_upperlimit = bounds ( 2, 3 )
  ELSE
     CALL stop_prg ( "pw_grid_assign", "no type set for the grid" )
  END IF
  
! finds valid g-points within grid
  IF ( pw_grid % para % mode == 0 ) THEN
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
               
             length = length_x ** 2 + length_y ** 2 + length_z ** 2
             length = twopi * twopi * length
           
             IF ( length <= cutoff * cutoff ) THEN
                gpt = gpt + 1
                pw_grid % g ( 1, gpt ) = twopi * length_x
                pw_grid % g ( 2, gpt ) = twopi * length_y
                pw_grid % g ( 3, gpt ) = twopi * length_z
                pw_grid % gsq ( gpt ) = length
                pw_grid % g_hat ( 1, gpt ) = l
                pw_grid % g_hat ( 2, gpt ) = m
                pw_grid % g_hat ( 3, gpt ) = n
             END IF
             
          END DO
       END DO
    END DO

  ELSE

    gpt = 0
    ip = pw_grid % para % my_pos
    DO i = 1, pw_grid % para % nyzray ( ip )
       n = pw_grid % para % yzp ( 2, i, ip )
       m = pw_grid % para % yzp ( 1, i, ip )
       IF ( n > n_upperlimit ) CYCLE
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
            
          length = length_x ** 2 + length_y ** 2 + length_z ** 2
          length = twopi * twopi * length
          
          IF ( length <= cutoff * cutoff ) THEN
             gpt = gpt + 1
             pw_grid % g ( 1, gpt ) = twopi * length_x
             pw_grid % g ( 2, gpt ) = twopi * length_y
             pw_grid % g ( 3, gpt ) = twopi * length_z
             pw_grid % gsq ( gpt ) = length
             pw_grid % g_hat ( 1, gpt ) = l
             pw_grid % g_hat ( 2, gpt ) = m
             pw_grid % g_hat ( 3, gpt ) = n
          END IF
             
       END DO
    END DO

  END IF
  
! Check the number of g-vectors for grid
  IF ( pw_grid % ngpts_cut_local /= gpt ) THEN
     CALL stop_prg ( "pw_grid_assign", "error re-counting the vectors" )
  END IF
  IF ( pw_grid % para % mode == 1 ) THEN
     CALL mp_sum ( gpt, pw_grid % para % group )
     IF ( pw_grid % ngpts_cut /= gpt ) &
     CALL stop_prg ( "pw_grid_assign", " sum on all processors", &
                     "error re-counting the vectors" )
  ENDIF
  
  pw_grid % have_g0 = .FALSE.
  pw_grid % first_gne0 = 2
  DO gpt = 1, pw_grid % ngpts_cut
     IF ( ALL ( pw_grid % g_hat ( :, gpt ) == 0 ) ) THEN
        pw_grid % have_g0 = .TRUE.
        pw_grid % first_gne0 = 1
        EXIT
     END IF
  END DO
  
  CALL pw_grid_set_maps ( pw_grid % grid_span, pw_grid % g_hat, &
       pw_grid % mapl, pw_grid % mapm, pw_grid % mapn, pw_grid % npts )
  
END SUBROUTINE pw_grid_assign

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_set_maps [1.1] *
!!
!!   NAME
!!     pw_grid_set_maps
!!
!!   FUNCTION
!!     Setup maps from 1d to 3d space
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (21-12-2000) : Size of g_hat locally determined
!!
!!   NOTES
!!     Maps are to full 3D space (not distributed)
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_set_maps ( grid_span, g_hat, mapl, mapm, mapn, npts )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: grid_span
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: g_hat
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: npts
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapl
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapm
  TYPE ( map_pn ), INTENT ( INOUT ) :: mapn
  
! Locals
  INTEGER :: l, m, n, gpt, ng
  
!------------------------------------------------------------------------------
  
  ng = SIZE ( g_hat ( 1, : ) )

  DO gpt = 1, ng
     l = g_hat ( 1, gpt )
     m = g_hat ( 2, gpt )
     n = g_hat ( 3, gpt )
     IF ( l < 0 ) THEN
        mapl % pos ( l ) = l + npts ( 1 )
     ELSE
        mapl % pos ( l ) = l
     END IF
     IF ( m < 0 ) THEN
        mapm % pos ( m ) = m + npts ( 2 )
     ELSE
        mapm % pos ( m ) = m
     END IF
     IF ( n < 0 ) THEN
        mapn % pos ( n ) = n + npts ( 3 )
     ELSE
        mapn % pos ( n ) = n
     END IF
     
! Generating the maps to the full 3-d space

     IF ( grid_span == HALFSPACE ) THEN
     
       IF ( l <= 0 ) THEN
          mapl % neg ( l ) = - l
       ELSE
          mapl % neg ( l ) = npts ( 1 ) - l
       END IF
       IF ( m <= 0 ) THEN
          mapm % neg ( m ) = - m
       ELSE
          mapm % neg ( m ) = npts ( 2 ) - m
       END IF
       IF ( n <= 0 ) THEN
          mapn % neg ( n ) = - n
       ELSE
          mapn % neg ( n ) = npts ( 3 ) - n
       END IF
     
     END IF

  END DO
  
END SUBROUTINE pw_grid_set_maps

!!*****
!****************************************************************************
!!****** pw_grids/pw_grid_find_bounds [1.1] *
!!
!!   NAME
!!     pw_grid_find_bounds
!!
!!   FUNCTION
!!     Find bounds of grid array
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (21-12-2000) : Simplify parameter list, bounds will be global
!!     JGH ( 8-01-2001) : Add check to FFT allowd grids (this now depends
!!                        on the FFT library.
!!                        Should the pw_grid_type have a reference to the FFT
!!                        library ?
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_find_bounds ( bounds, h_inv, cutoff )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, DIMENSION ( 2, 3 ), INTENT ( OUT ) :: bounds
  REAL ( dbl ), INTENT ( IN ) :: cutoff
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( IN ) :: h_inv
  
! Locals
  INTEGER :: l, m, n, test_max, rin(3), rout(3)
  REAL ( dbl ) :: length
  REAL ( dbl ) :: gmat ( 3, 3 ), gi ( 3 )
  
!------------------------------------------------------------------------------
  
  bounds ( :, : ) = 0

  gmat = MATMUL ( h_inv, TRANSPOSE ( h_inv ) )
  gi = 1._dbl
  length = twopi * SQRT ( dotprod_3d ( matvec_3x3 ( gmat, gi ), gi ) )
  test_max = 2 * NINT ( cutoff / length )
  
! finds valid g-points within grid
  DO n = 0, test_max
     gi ( 3 ) = REAL ( n, dbl )
     DO m = 0, test_max
        gi ( 2 ) = REAL ( m, dbl )
        DO l = 0, test_max
           gi ( 1 ) = REAL ( l, dbl )
           length = twopi * twopi * dotprod_3d ( matvec_3x3 ( gmat, gi ), gi )
           
           IF ( length <= cutoff * cutoff ) THEN
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
       CALL print_warning ( "pw_grid_find_bounds", &
                            "initial bounds maybe too small:" )
  END IF
  
  rin(1) = 2 * bounds ( 2, 1 ) + 1
  CALL fft_radix_operations ( rin(1), rout(1), FFT_RADIX_NEXT )
  bounds ( 1, 1 ) = -rout(1)/2
  rin(2) = 2 * bounds ( 2, 2 ) + 1
  CALL fft_radix_operations ( rin(2), rout(2), FFT_RADIX_NEXT )
  bounds ( 1, 2 ) = -rout(2)/2
  rin(3) = 2 * bounds ( 2, 3 ) + 1
  CALL fft_radix_operations ( rin(3), rout(3), FFT_RADIX_NEXT )
  bounds ( 1, 3 ) = -rout(3)/2

  bounds ( 2, : ) = bounds ( 1, : ) + rout - 1
  
END SUBROUTINE pw_grid_find_bounds

!!*****
!****************************************************************************
!!****** pw_grids/pw_grid_allocate [1.1] *
!!
!!   NAME
!!     pw_grid_allocate
!!
!!   FUNCTION
!!     Allocate all (Pointer) Arrays in pw_grid
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Added status variable 
!!                        Bounds of arrays now from calling routine, this
!!                        makes it independent from parallel setup
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_allocate ( pw_grid, ng, bounds )
  IMPLICIT NONE
  
! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  INTEGER, INTENT ( IN ) :: ng
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: bounds

! LOCALS
  INTEGER :: allocstat
  
!------------------------------------------------------------------------------
  
  ALLOCATE ( pw_grid % g ( 3, ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_allocate", "g", 3*ng )
  ALLOCATE ( pw_grid % gsq ( ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_allocate", "gsq", ng )
  ALLOCATE ( pw_grid % g_hat ( 3, ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_allocate", "g_hat", 3*ng )

  ALLOCATE ( pw_grid % mapl % pos ( bounds ( 1, 1 ):bounds ( 2, 1 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapl % pos", 0 )
  ALLOCATE ( pw_grid % mapl % neg ( bounds ( 1, 1 ):bounds ( 2, 1 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapl % neg", 0 )

  ALLOCATE ( pw_grid % mapm % pos ( bounds ( 1, 2 ):bounds ( 2, 2 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapm % pos", 0 )
  ALLOCATE ( pw_grid % mapm % neg ( bounds ( 1, 2 ):bounds ( 2, 2 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapm % neg", 0 )

  ALLOCATE ( pw_grid % mapn % pos ( bounds ( 1, 3 ):bounds ( 2, 3 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapn % pos", 0 )
  ALLOCATE ( pw_grid % mapn % neg ( bounds ( 1, 3 ):bounds ( 2, 3 ) ), &
             STAT = allocstat )
  IF ( allocstat /= 0 ) &
       CALL stop_memory ( "pw_grid_allocate", "mapn % neg", 0 )
  
END SUBROUTINE pw_grid_allocate

!!*****
!****************************************************************************
!!****** pw_grids/pw_grid_sort [1.1] *
!!
!!   NAME
!!     pw_grid_sort
!!
!!   FUNCTION
!!     Sort g-vectors according to length
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : allocate idx, ng = SIZE ( pw_grid % gsq ) the
!!                        sorting is local and independent from parallelisation
!!                        WARNING: Global ordering depends now on the number
!!                                 of cpus.
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_sort ( pw_grid )
  
  IMPLICIT NONE
  
! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  
! Locals
  INTEGER :: gpt, gpt_curr, gpt_tmp, int_tmp ( 3 ), allocstat, ng
  INTEGER, ALLOCATABLE, DIMENSION ( : ) :: idx 
  REAL ( dbl ) :: real_tmp ( 3 )
  
!------------------------------------------------------------------------------
  
  ng = SIZE ( pw_grid % gsq ( : ) )
  ALLOCATE ( idx ( ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "idx", ng )

  CALL sort ( pw_grid % gsq, ng, idx )
  
  DO gpt = 1, ng - 1
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

  DEALLOCATE ( idx, STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "idx" )
  
END SUBROUTINE pw_grid_sort

!!*****
!****************************************************************************
!!****** pw_grids/pw_grid_change [1.1] *
!!
!!   NAME
!!     pw_grid_change
!!
!!   FUNCTION
!!     Recalculate the g-vectors after a change of the box
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : get local grid size from definition of g.
!!                        Assume that gsq is allocated.
!!                        Local routine, no information on distribution of
!!                        PW required.
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_change ( cell, pw_grid )
  IMPLICIT NONE
  
! Argument
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  
! Locals
  INTEGER :: gpt
  REAL ( dbl ) :: lr, mr, nr
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: g
  
!------------------------------------------------------------------------------
  
  h_inv => cell % h_inv
  g => pw_grid % g
  
  DO gpt = 1, SIZE ( g ( 1, : ) )

     lr = REAL ( pw_grid % g_hat ( 1, gpt ), dbl ) * twopi
     mr = REAL ( pw_grid % g_hat ( 2, gpt ), dbl ) * twopi
     nr = REAL ( pw_grid % g_hat ( 3, gpt ), dbl ) * twopi
     
     g ( 1, gpt ) &
          = lr * h_inv(1,1) &
          + mr * h_inv(2,1) &
          + nr * h_inv(3,1)
     g ( 2, gpt ) &
          = lr * h_inv(1,2) &
          + mr * h_inv(2,2) &
          + nr * h_inv(3,2)
     g ( 3, gpt ) &
          = lr * h_inv(1,3) &
          + mr * h_inv(2,3) &
          + nr * h_inv(3,3)
     
     pw_grid % gsq ( gpt ) = SUM ( g ( :, gpt ) ** 2 )

  END DO
  
END SUBROUTINE pw_grid_change

!!*****
!****************************************************************************
!!****** pw_grids/pw_find_cutoff [1.1] *
!!
!!   NAME
!!     pw_find_cutoff
!!
!!   FUNCTION
!!     Given a grid and a box, calculate the corresponding cutoff
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Deleted some strange comments
!!
!!   NOTES
!!     This routine is local. It works independent from the distribution
!!     of PW on processors.
!!     npts is the grid size for the full box.
!!
!!   SOURCE
!******************************************************************************

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
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = length
  
! compute 2*pi*h_inv^t*g  where g = (0,npts[2],0)
  gdum ( : ) = twopi * box % h_inv ( 2, : ) &
       * REAL ( npts ( 2 ), dbl ) / 2.0_dbl
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )
  
! compute 2*pi*h_inv^t*g  where g = (0,0,npts[3])
  gdum ( : ) = twopi * box % h_inv ( 3, : ) &
       * REAL ( npts ( 3 ), dbl ) / 2.0_dbl
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )
  
  cutoff = gcut - 1.0E-5_dbl
  
END SUBROUTINE pw_find_cutoff

!!*****
!******************************************************************************

END MODULE pw_grids

!******************************************************************************
