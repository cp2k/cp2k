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
!!     JGH (07-02-2001) : Added constructor and destructor routines
!!
!!   NOTES
!!     pw_grid_construct : set the defaults
!!     pw_grid_destruct : release all memory connected to type
!!     pw_grid_setup  : main routine to set up a grid
!!          input: cell (the box for the grid)
!!                 pw_grid (the grid; pw_grid%grid_span has to be set)
!!                 cutoff (optional, if not given pw_grid%bounds has to be set)
!!                 pe_group (optional, if not given we have a local grid)
!!
!!                 if no cutoff or a negative cutoff is given, all g-vectors
!!                 in the box are included (no spherical cutoff)
!!
!!                 for a distributed setup the array in para rs_dims has to
!!                 be initialized
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
                              mp_max, mp_min, mp_cart_coords, mp_dims_create, &
                              mp_cart_rank, mp_sync, mp_comm_free, mp_comm_dup
  USE pw_grid_types, ONLY : pw_grid_type, map_pn, HALFSPACE, FULLSPACE
  USE output_utilities, ONLY : print_warning
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program
  USE util, ONLY : matvec_3x3, dotprod_3d, sort, get_limit
  USE fft_tools, ONLY : fft_radix_operations, FFT_RADIX_NEXT, &
                        FFT_RADIX_NEXT_ODD
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_grid_construct, pw_grid_destruct, pw_grid_compare
  PUBLIC :: pw_grid_setup, pw_grid_change, pw_find_cutoff

  INTEGER :: grid_tag = 0
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pw_grids/pw_grid_construct [1.0] *
!!
!!   NAME
!!     pw_grid_construct
!!
!!   FUNCTION
!!     Initialize a PW grid with all defaults
!!
!!   AUTHOR
!!     JGH (7-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_construct ( pw_grid )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid

  pw_grid % bounds = 0
  pw_grid % cutoff = 0._dbl
  pw_grid % grid_span = FULLSPACE
  pw_grid % para % rs_dims = 0

END SUBROUTINE pw_grid_construct

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_destruct [1.0] *
!!
!!   NAME
!!     pw_grid_destruct
!!
!!   FUNCTION
!!     Deallocate all pointers within a pw_grid_type
!!
!!   AUTHOR
!!     JGH (7-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Mar-2001) : release group communicators
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_destruct ( pw_grid )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid

! Local
  INTEGER :: ierr

  IF ( ASSOCIATED ( pw_grid % g ) ) THEN
    DEALLOCATE ( pw_grid % g, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", "pw_grid % g" )
  END IF
  IF ( ASSOCIATED ( pw_grid % gsq ) ) THEN
    DEALLOCATE ( pw_grid % gsq, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", "pw_grid % gsq" )
  END IF
  IF ( ASSOCIATED ( pw_grid % g_hat ) ) THEN
    DEALLOCATE ( pw_grid % g_hat, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", "pw_grid % g_hat" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapl % pos ) ) THEN
    DEALLOCATE ( pw_grid % mapl % pos, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapl % pos" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapm % pos ) ) THEN
    DEALLOCATE ( pw_grid % mapm % pos, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapm % pos" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapn % pos ) ) THEN
    DEALLOCATE ( pw_grid % mapn % pos, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapn % pos" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapl % neg ) ) THEN
    DEALLOCATE ( pw_grid % mapl % neg, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapl % neg" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapm % neg ) ) THEN
    DEALLOCATE ( pw_grid % mapm % neg, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapm % neg" )
  END IF
  IF ( ASSOCIATED ( pw_grid % mapn % neg ) ) THEN
    DEALLOCATE ( pw_grid % mapn % neg, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
       "pw_grid % mapn % neg" )
  END IF
  IF ( pw_grid % para % mode == 1 ) THEN
    IF ( ASSOCIATED ( pw_grid % para % yzp ) ) THEN
      DEALLOCATE ( pw_grid % para % yzp, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
        "pw_grid % para % yzp" )
    END IF
    IF ( ASSOCIATED ( pw_grid % para % yzq ) ) THEN
      DEALLOCATE ( pw_grid % para % yzq, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
        "pw_grid % para % yzq" )
    END IF
    IF ( ASSOCIATED ( pw_grid % para % nyzray ) ) THEN
      DEALLOCATE ( pw_grid % para % nyzray, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
        "pw_grid % para % nyzray" )
    END IF
    IF ( ASSOCIATED ( pw_grid % para % bo ) ) THEN
      DEALLOCATE ( pw_grid % para % bo, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_destruct", &
        "pw_grid % para % bo" )
    END IF

    ! also release groups
    CALL mp_comm_free ( pw_grid % para % group )
    CALL mp_comm_free ( pw_grid % para % rs_group )

  END IF
  

END SUBROUTINE pw_grid_destruct

!!*****
!******************************************************************************
!!****** pw_grids/pw_grid_compare [1.0] *
!!
!!   NAME
!!     pw_grid_compare
!!
!!   FUNCTION
!!     Check if two pw_grids are equal
!!
!!   AUTHOR
!!     JGH (14-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

FUNCTION pw_grid_compare ( grida, gridb ) RESULT ( equal )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( IN ) :: grida
  TYPE ( pw_grid_type ), INTENT ( IN ) :: gridb
  LOGICAL :: equal

!------------------------------------------------------------------------------

  IF ( grida % identifier == gridb % identifier ) THEN
    equal = .TRUE.
  ELSE
! for the moment all grids with different identifiers are considered as not equal
! later we can get this more relaxed
    equal = .FALSE.
  END IF

END function pw_grid_compare

!!*****
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
!!     JGH (20-Dec-2000) : Adapted for parallel use
!!     JGH (28-Feb-2001) : New optional argument fft_usage
!!     JGH (21-Mar-2001) : Reference grid code
!!     JGH (21-Mar-2001) : New optional argument symm_usage
!!     JGH (22-Mar-2001) : Simplify group assignment (mp_comm_dup)
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_setup ( cell, pw_grid, cutoff, pe_group, info, fft_usage, &
                           symm_usage, ref_grid, orthorhombic )

  IMPLICIT NONE

! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: cell
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: cutoff ! in atomic units
  INTEGER, OPTIONAL, INTENT ( IN ) :: pe_group ! parent group for this grid
  INTEGER, OPTIONAL, INTENT ( IN ) :: info  ! output unit for information on grid
  LOGICAL, OPTIONAL, INTENT ( IN ) :: fft_usage  ! has the grid size to be 
                                                 ! compatible with the FFT
  LOGICAL, OPTIONAL, INTENT ( IN ) :: symm_usage ! has the grid size to be 
                                                 ! symmetric (g <-> -g)    
  LOGICAL, OPTIONAL, INTENT ( IN ) :: orthorhombic
  TYPE ( pw_grid_type ), OPTIONAL, INTENT ( IN ) :: ref_grid

! Locals
  REAL ( dbl ) :: ecut
  INTEGER, DIMENSION ( :, : ), ALLOCATABLE :: yz_mask
  INTEGER :: allocstat, group_size, my_pos, i, n(3), ntest(3)
  REAL ( dbl ) :: rv(3,3)
  INTEGER :: aver, amax, amin
  LOGICAL :: fft, symmetry

!------------------------------------------------------------------------------
  ! default is to use only fft compatible grids
  IF ( PRESENT ( fft_usage ) ) THEN
    fft = fft_usage
  ELSE
    fft = .TRUE.
  END IF
  ! default is to use only symmetric grids
  IF ( PRESENT ( symm_usage ) ) THEN
    symmetry = symm_usage
  ELSE
    symmetry = .TRUE.
  END IF

  ! assign a unique tag to this grid
  grid_tag = grid_tag + 1
  pw_grid % identifier = grid_tag

  ! parallel info
  IF (PRESENT ( pe_group ) ) THEN
    CALL mp_environ ( group_size, my_pos, pe_group )
    IF ( group_size > 1 ) THEN
      CALL mp_comm_dup ( pe_group, pw_grid % para % group )
      CALL mp_environ ( pw_grid % para % group_size, &
                        pw_grid % para % my_pos, &
                        pw_grid % para % group )
      pw_grid % para % group_head_id = 0
      pw_grid % para % group_head = &
                ( pw_grid % para % group_head_id == pw_grid % para % my_pos )
      pw_grid % para % mode = 1
    ELSE
      pw_grid % para % mode = 0
    END IF
  ELSE
    pw_grid % para % mode = 0
  END IF

  IF ( PRESENT ( cutoff ) ) THEN
     ecut = cutoff
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) &
               CALL pw_grid_find_bounds &
               ( pw_grid % bounds, cell % h_inv, ecut, fft, orthorhombic )
     pw_grid % spherical = .TRUE.
     pw_grid % cutoff = ABS ( ecut )
  ELSE
     ecut = 1.e10_dbl              ! all g-vectors in the box will be included
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) THEN
        CALL stop_program ( "grid_setup", &
                            "provide initial values for bounds yourself" )
     END IF
     pw_grid % spherical = .FALSE.
     pw_grid % cutoff = 0._dbl

     IF ( fft ) THEN

       ! select FFT compatible bounds
       ntest = pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1
       ! without a cutoff and HALFSPACE we have to be sure that there is
       ! a negative counterpart to every g vector (-> odd number of grid points)
       IF ( pw_grid % grid_span == HALFSPACE .AND. symmetry ) THEN

         CALL fft_radix_operations ( ntest(1), n(1), FFT_RADIX_NEXT_ODD )
         CALL fft_radix_operations ( ntest(2), n(2), FFT_RADIX_NEXT_ODD )
         CALL fft_radix_operations ( ntest(3), n(3), FFT_RADIX_NEXT_ODD )

       ELSE

         CALL fft_radix_operations ( ntest(1), n(1), FFT_RADIX_NEXT )
         CALL fft_radix_operations ( ntest(2), n(2), FFT_RADIX_NEXT )
         CALL fft_radix_operations ( ntest(3), n(3), FFT_RADIX_NEXT )

       END IF

     ELSE

       n = pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1
       ! without a cutoff and HALFSPACE we have to be sure that there is
       ! a negative counterpart to every g vector (-> odd number of grid points)
       IF ( pw_grid % grid_span == HALFSPACE .AND. symmetry ) &
          n = n + MOD ( n + 1, 2 )

     END IF

     pw_grid % bounds ( 1, : ) = - n / 2
     pw_grid % bounds ( 2, : ) = pw_grid % bounds ( 1, : ) + n - 1
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

  ! Check if reference grid is compatible
  IF ( PRESENT ( ref_grid ) ) THEN
    IF ( pw_grid % para % mode /= ref_grid % para % mode ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible parallelisation scheme" )
    END IF
    IF ( pw_grid % para % mode == 1 ) THEN
      IF ( pw_grid % para % group /= ref_grid % para % group ) THEN
        CALL stop_program ( "pw_grid_setup", "Incompatible MPI groups" )
      END IF
    END IF
    IF ( pw_grid % grid_span /= ref_grid % grid_span ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible grid types" )
    END IF
    IF ( pw_grid % spherical .EQV. ref_grid % spherical ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible cutoff schemes" )
    END IF
  END IF

  ! Distribute grid
  CALL pw_grid_distribute ( pw_grid, yz_mask, ref_grid )

  ! Allocate the grid fields
  CALL pw_grid_allocate ( pw_grid, pw_grid % ngpts_cut_local, &
                          pw_grid % bounds )

  ! Fill in the grid structure
  CALL pw_grid_assign ( cell % h_inv, pw_grid, ecut )
  
  ! Sort g vector wrt length (only local for each processor)
  CALL pw_grid_sort ( pw_grid, ref_grid )
  
  CALL pw_grid_remap ( pw_grid, yz_mask )

  DEALLOCATE ( yz_mask )
  IF ( allocstat /= 0 ) call stop_memory ( "pw_grid_setup", "yz_mask" )

  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts, dbl )
  pw_grid % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ) ** 2, 1 ) ) &
       / REAL ( pw_grid % npts ( : ), dbl )

!
! Output: All the information of this grid type
!

  IF ( PRESENT ( info ) ) THEN
    IF ( pw_grid % para % mode == 0  .AND. info >= 0 ) THEN
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
    ELSEIF ( info >= 0 ) THEN
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
!!     JGH (01-Mar-2001) optional reference grid 
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_distribute ( pw_grid, yz_mask, ref_grid )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  INTEGER, DIMENSION ( :, : ), INTENT ( INOUT ) :: yz_mask
  TYPE ( pw_grid_type ), OPTIONAL, INTENT ( IN ) :: ref_grid

! Locals
  INTEGER :: gmin, gmax, nx, ny, nz, np, ns, i, j, k, l, m, n
  INTEGER :: lo ( 2 ), ip, ierr, gneg, lby, lbz, coor ( 2 ), rsd ( 2 ), ipl
  REAL ( dbl ) :: tfun, tt

!------------------------------------------------------------------------------

  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )

  pw_grid % ngpts = PRODUCT ( pw_grid % npts )

  IF ( pw_grid % para % mode == 0) THEN

    pw_grid % bounds_local = pw_grid % bounds
    pw_grid % npts_local = pw_grid % npts
    pw_grid % ngpts_cut_local = pw_grid % ngpts_cut
    pw_grid % ngpts_local = PRODUCT ( pw_grid % npts_local )

  ELSE

!..find the real space distribution
    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )

    np = pw_grid % para % group_size

    IF ( PRODUCT ( pw_grid % para % rs_dims ) == 0 ) THEN

      ns = INT ( SQRT ( REAL ( np, dbl ) ) )

      tfun = 1.e20_dbl
      DO i = ns, 2, -1
        IF ( MOD ( np, i ) == 0 ) THEN
          j = np/i
          k = nx*ny - i*(nx/i) * j*(ny/j)
          IF ( k == 0 ) THEN
            tfun = 0._dbl
            m = i
          ELSE
            tt = 1._dbl - REAL ( k, dbl ) / REAL ( np, dbl )
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
      pw_grid % para % rs_dims ( 1 ) = m
      pw_grid % para % rs_dims ( 2 ) = np/m

    ELSEIF ( PRODUCT ( pw_grid % para % rs_dims ) /= np ) THEN

      pw_grid % para % rs_dims = 0
      CALL mp_dims_create ( np, pw_grid % para % rs_dims )

    END IF
!..create group for real space distribution
    CALL mp_cart_create ( pw_grid % para % group, 2, &
                          pw_grid % para % rs_dims, &
                          pw_grid % para % rs_pos, &
                          pw_grid % para % rs_group )
    CALL mp_cart_rank ( pw_grid % para % rs_group, &
                        pw_grid % para % rs_pos, &
                        pw_grid % para % rs_mpo )
    lo = get_limit ( nx, pw_grid % para % rs_dims ( 1 ), &
                     pw_grid % para % rs_pos ( 1 ) )
    pw_grid % bounds_local ( :, 1 ) = lo + pw_grid % bounds ( 1, 1 ) - 1
    lo = get_limit ( ny, pw_grid % para % rs_dims ( 2 ), &
                     pw_grid % para % rs_pos ( 2 ) )
    pw_grid % bounds_local ( :, 2 ) = lo + pw_grid % bounds ( 1, 2 ) - 1
    pw_grid % bounds_local ( :, 3 ) = pw_grid % bounds ( :, 3 )
    pw_grid % npts_local ( : ) = pw_grid % bounds_local ( 2, : ) &
                                 - pw_grid % bounds_local ( 1, : ) + 1

!..the third distribution is needed for the second step in the FFT
    ALLOCATE ( pw_grid % para % bo ( 2, 3, 0:np-1, 3 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
       "pw_grid % para % bo", 2*np )
    rsd = pw_grid % para % rs_dims
    DO ip = 0, np - 1
      CALL mp_cart_coords ( pw_grid % para % rs_group, ip, coor )
      ! distribution xyZ
      pw_grid % para % bo ( 1:2, 1, ip, 1 ) = get_limit ( nx, rsd ( 1 ), coor ( 1 ) )
      pw_grid % para % bo ( 1:2, 2, ip, 1 ) = get_limit ( ny, rsd ( 2 ), coor ( 2 ) )
      pw_grid % para % bo ( 1, 3, ip, 1 ) = 1
      pw_grid % para % bo ( 2, 3, ip, 1 ) = nz
      ! distribution xYz
      pw_grid % para % bo ( 1:2, 1, ip, 2 ) = get_limit ( nx, rsd ( 1 ), coor ( 1 ) )
      pw_grid % para % bo ( 1, 2, ip, 2 ) = 1
      pw_grid % para % bo ( 2, 2, ip, 2 ) = ny
      pw_grid % para % bo ( 1:2, 3, ip, 2 ) = get_limit ( nz, rsd ( 2 ), coor ( 2 ) )
      ! distribution Xyz
      pw_grid % para % bo ( 1, 1, ip, 3 ) = 1
      pw_grid % para % bo ( 2, 1, ip, 3 ) = nx
      pw_grid % para % bo ( 1:2, 2, ip, 3 ) = get_limit ( ny, rsd ( 1 ), coor ( 1 ) )
      pw_grid % para % bo ( 1:2, 3, ip, 3 ) = get_limit ( nz, rsd ( 2 ), coor ( 2 ) )
    END DO

!..find the g space distribution
    pw_grid % ngpts_cut_local = 0

    ALLOCATE ( pw_grid % para % nyzray ( 0: np-1 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
          "pw_grid % para % nyzray", np )

    ALLOCATE ( pw_grid % para % yzq ( ny, nz ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
       "pw_grid % para % yzq", ny*nz )

    IF ( pw_grid % spherical .OR. pw_grid % grid_span == HALFSPACE ) THEN

      pw_grid % para % ray_distribution = .TRUE.

      pw_grid % para % yzq = -1
      IF ( PRESENT ( ref_grid ) ) THEN
        ! tag all vectors from the reference grid
        CALL pre_tag ( pw_grid, yz_mask, ref_grid )
      END IF

      ! Round Robin distribution 
      ! Processors 0 .. NP-1, NP-1 .. 0  get the largest remaining batch
      ! of g vectors in turn
  
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
            pw_grid % para % yzp ( 1, ns, ip ) = j
            pw_grid % para % yzp ( 2, ns, ip ) = i
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

      pw_grid % ngpts_local = PRODUCT ( pw_grid % npts_local )

    ELSE
      !
      !  block distribution of g vectors, we do not have a spherical cutoff
      !

      pw_grid % para % ray_distribution = .FALSE.

      DO ip = 0, np - 1
        m = pw_grid % para % bo ( 2, 2, ip, 3 ) - &
            pw_grid % para % bo ( 1, 2, ip, 3 ) + 1
        n = pw_grid % para % bo ( 2, 3, ip, 3 ) - &
            pw_grid % para % bo ( 1, 3, ip, 3 ) + 1
        pw_grid % para % nyzray ( ip ) = n*m
      END DO

      ipl = pw_grid % para % rs_mpo
      l = pw_grid % para % bo ( 2, 1, ipl, 3 ) - &
          pw_grid % para % bo ( 1, 1, ipl, 3 ) + 1
      m = pw_grid % para % bo ( 2, 2, ipl, 3 ) - &
          pw_grid % para % bo ( 1, 2, ipl, 3 ) + 1
      n = pw_grid % para % bo ( 2, 3, ipl, 3 ) - &
          pw_grid % para % bo ( 1, 3, ipl, 3 ) + 1
      pw_grid % ngpts_cut_local = l * m * n
      pw_grid % ngpts_local = pw_grid % ngpts_cut_local

      pw_grid % para % yzq = 0
      ny = pw_grid % para % bo ( 2, 2, ipl, 3 ) - &
           pw_grid % para % bo ( 1, 2, ipl, 3 ) + 1
      DO n = pw_grid % para % bo ( 1, 3, ipl, 3 ), &
             pw_grid % para % bo ( 2, 3, ipl, 3 )
        i = n - pw_grid % para % bo ( 1, 3, ipl, 3 )
        DO m = pw_grid % para % bo ( 1, 2, ipl, 3 ), &
               pw_grid % para % bo ( 2, 2, ipl, 3 )
          j = m - pw_grid % para % bo ( 1, 2, ipl, 3 ) + 1
          pw_grid % para % yzq ( m, n ) = j + i * ny
        END DO
      END DO

    END IF

  END IF

END SUBROUTINE pw_grid_distribute

!******************************************************************************

SUBROUTINE pre_tag ( pw_grid, yz_mask, ref_grid )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  INTEGER, DIMENSION ( :, : ), INTENT ( INOUT ) :: yz_mask
  TYPE ( pw_grid_type ), INTENT ( IN ) :: ref_grid

! Locals
  INTEGER :: ip, ig, y, z, ny, nz, lby, lbz, uby, ubz, gmax

!------------------------------------------------------------------------------

  ny = ref_grid % npts ( 2 )
  nz = ref_grid % npts ( 3 )
  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )
  uby = pw_grid % bounds ( 2, 2 )
  ubz = pw_grid % bounds ( 2, 3 )
 
  ! loop over all processors and all g vectors yz lines on this processor
  DO ip = 0, ref_grid % para % group_size
    DO ig = 1, ref_grid % para % nyzray ( ip )
      ! go from mapped coordinates to original coordinates
      ! 0 .. N-1 -> -n/2 .. (n+1)/2
      y = ref_grid % para % yzp ( 1, ig, ip ) - 1
      IF ( y > ny/2 ) y = y - ny 
      z = ref_grid % para % yzp ( 2, ig, ip ) - 1
      IF ( z > nz/2 ) z = z - nz 
      ! check if this is inside the realm of the new grid
      IF ( y < lby .OR. y > uby .OR. z < lbz .OR. z > ubz ) CYCLE
      ! go to shifted coordinates
      y = y - lby + 1
      z = z - lbz + 1
      gmax = yz_mask ( y, z )
      yz_mask ( y, z ) = 0
      IF ( ip == pw_grid % para % my_pos ) THEN
         pw_grid % ngpts_cut_local = pw_grid % ngpts_cut_local + gmax
      END IF
      pw_grid % para % yzq ( y, z ) = ip
    END DO
  END DO

END SUBROUTINE pre_tag

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
     CALL stop_program ( "pw_grid_count", "no type set for the grid" )
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
     CALL stop_program ( "pw_grid_count", "para % mode not specified" )
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
           IF ( 0.5_dbl * length <= cutoff ) THEN
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
  INTEGER :: gpt, i, l, m, n, n_upperlimit, ip, lby, lbz, nn, mm, ll
  REAL ( dbl ) :: length_x, length_y, length_z, length
  INTEGER, DIMENSION ( 2, 3 ) :: bounds, bol
  
!------------------------------------------------------------------------------
  
  bounds = pw_grid % bounds
  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )
  
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     n_upperlimit = 0
  ELSE IF ( pw_grid % grid_span == FULLSPACE ) THEN
     n_upperlimit = bounds ( 2, 3 )
  ELSE
     CALL stop_program ( "pw_grid_assign", "no type set for the grid" )
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

             IF ( 0.5_dbl * length <= cutoff ) THEN
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

    IF ( pw_grid % para % ray_distribution ) THEN

      gpt = 0
      ip = pw_grid % para % my_pos
      DO i = 1, pw_grid % para % nyzray ( ip )
         n = pw_grid % para % yzp ( 2, i, ip ) + lbz - 1
         m = pw_grid % para % yzp ( 1, i, ip ) + lby - 1
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
          
            IF ( 0.5_dbl * length <= cutoff ) THEN
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

    ELSE

      bol = pw_grid % para % bo ( :, :, pw_grid % para % rs_mpo, 3 )
      gpt = 0
      DO n = bounds ( 1, 3 ), bounds ( 2, 3 )
         IF ( n < 0 ) THEN
            nn = n + pw_grid % npts ( 3 ) + 1
         ELSE
            nn = n + 1
         END IF
         IF ( nn < bol ( 1, 3 ) .OR. nn > bol ( 2, 3 ) ) CYCLE
         DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
            IF ( m < 0 ) THEN
               mm = m + pw_grid % npts ( 2 ) + 1
            ELSE
               mm = m + 1
            END IF
            IF ( mm < bol ( 1, 2 ) .OR. mm > bol ( 2, 2 ) ) CYCLE
            DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
               IF ( l < 0 ) THEN
                  ll = l + pw_grid % npts ( 1 ) + 1
               ELSE
                  ll = l + 1
               END IF
               IF ( ll < bol ( 1, 1 ) .OR. ll > bol ( 2, 1 ) ) CYCLE

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
 
               gpt = gpt + 1
               pw_grid % g ( 1, gpt ) = twopi * length_x
               pw_grid % g ( 2, gpt ) = twopi * length_y
               pw_grid % g ( 3, gpt ) = twopi * length_z
               pw_grid % gsq ( gpt ) = length
               pw_grid % g_hat ( 1, gpt ) = l
               pw_grid % g_hat ( 2, gpt ) = m
               pw_grid % g_hat ( 3, gpt ) = n

            END DO
         END DO
      END DO

    END IF

  END IF

! Check the number of g-vectors for grid
  IF ( pw_grid % ngpts_cut_local /= gpt ) THEN
     CALL stop_program ( "pw_grid_assign", "error re-counting the vectors" )
  END IF
  IF ( pw_grid % para % mode == 1 ) THEN
     CALL mp_sum ( gpt, pw_grid % para % group )
     IF ( pw_grid % ngpts_cut /= gpt ) &
     CALL stop_program ( "pw_grid_assign", " sum on all processors"//&
                         "error re-counting the vectors" )
  ENDIF

  pw_grid % have_g0 = .FALSE.
  pw_grid % first_gne0 = 1
  DO gpt = 1, pw_grid % ngpts_cut_local
     IF ( ALL ( pw_grid % g_hat ( :, gpt ) == 0 ) ) THEN
        pw_grid % have_g0 = .TRUE.
        pw_grid % first_gne0 = 2
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
!!     JGH (28-02-2001) : Only do conditional check for FFT
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_find_bounds ( bounds, h_inv, cutoff, fft, orthorhombic )

  IMPLICIT NONE

! Arguments
  INTEGER, DIMENSION ( 2, 3 ), INTENT ( OUT ) :: bounds
  REAL ( dbl ), INTENT ( IN ) :: cutoff
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( IN ) :: h_inv
  LOGICAL, INTENT ( IN ) :: fft
  LOGICAL, OPTIONAL, INTENT ( IN ) :: orthorhombic

! Locals
  INTEGER :: l, m, n, test_max, rin(3), rout(3)
  REAL ( dbl ) :: length
  REAL ( dbl ) :: gmat ( 3, 3 ), gi ( 3 )

!------------------------------------------------------------------------------

  bounds ( :, : ) = 0

  IF (PRESENT(orthorhombic)) THEN

    bounds(2,1) = FLOOR(SQRT(2.0_dbl*cutoff)/(twopi*h_inv(1,1)))
    bounds(2,2) = FLOOR(SQRT(2.0_dbl*cutoff)/(twopi*h_inv(2,2)))
    bounds(2,3) = FLOOR(SQRT(2.0_dbl*cutoff)/(twopi*h_inv(3,3)))

  ELSE

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
             length = 0.5_dbl * twopi * twopi * &
                      dotprod_3d ( matvec_3x3 ( gmat, gi ), gi )
           
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
         CALL print_warning ( "pw_grid_find_bounds", &
                              "initial bounds maybe too small:" )
    END IF

  END IF

  IF ( fft ) THEN
    rin(1) = 2 * bounds ( 2, 1 ) + 1
    CALL fft_radix_operations ( rin(1), rout(1), FFT_RADIX_NEXT )
    bounds ( 1, 1 ) = -rout(1)/2
    rin(2) = 2 * bounds ( 2, 2 ) + 1
    CALL fft_radix_operations ( rin(2), rout(2), FFT_RADIX_NEXT )
    bounds ( 1, 2 ) = -rout(2)/2
    rin(3) = 2 * bounds ( 2, 3 ) + 1
    CALL fft_radix_operations ( rin(3), rout(3), FFT_RADIX_NEXT )
    bounds ( 1, 3 ) = -rout(3)/2
  ELSE
    rout = rin
  END IF

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
!!     JGH (28-02-2001) : check for ordering against reference grid
!!     JGH (01-05-2001) : sort spherical cutoff grids also within shells
!!                        reference grids for non-spherical cutoffs
!!     JGH (20-06-2001) : do not sort non-spherical grids
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_sort ( pw_grid, ref_grid )

  IMPLICIT NONE

! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ) :: pw_grid
  TYPE ( pw_grid_type ), OPTIONAL, INTENT ( IN ) :: ref_grid

! Locals
  INTEGER :: ig, allocstat, ng, ngr
  INTEGER, ALLOCATABLE, DIMENSION ( : ) :: idx, int_tmp
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( : ) :: real_tmp

!------------------------------------------------------------------------------

  IF ( pw_grid % spherical ) THEN

    ! spherical grids are (locally) ordered by length of G-vectors
    ng = SIZE ( pw_grid % gsq ( : ) )
    ALLOCATE ( idx ( ng ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "idx", ng )
    ALLOCATE ( int_tmp ( ng ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "int_tmp", ng )
    ALLOCATE ( real_tmp ( ng ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "real_tmp", ng )

    CALL sort ( pw_grid % gsq, ng, idx )
    CALL sort_shells ( pw_grid % gsq, pw_grid % g_hat, idx )

    real_tmp ( 1:ng ) = pw_grid % g ( 1, 1:ng )
    pw_grid % g ( 1, 1:ng ) = real_tmp ( idx ( 1:ng ) )
    real_tmp ( 1:ng ) = pw_grid % g ( 2, 1:ng )
    pw_grid % g ( 2, 1:ng ) = real_tmp ( idx ( 1:ng ) )
    real_tmp ( 1:ng ) = pw_grid % g ( 3, 1:ng )
    pw_grid % g ( 3, 1:ng ) = real_tmp ( idx ( 1:ng ) )
    int_tmp ( 1:ng ) = pw_grid % g_hat ( 1, 1:ng )
    pw_grid % g_hat ( 1, 1:ng ) = int_tmp ( idx ( 1:ng ) )
    int_tmp ( 1:ng ) = pw_grid % g_hat ( 2, 1:ng )
    pw_grid % g_hat ( 2, 1:ng ) = int_tmp ( idx ( 1:ng ) )
    int_tmp ( 1:ng ) = pw_grid % g_hat ( 3, 1:ng )
    pw_grid % g_hat ( 3, 1:ng ) = int_tmp ( idx ( 1:ng ) )

    DEALLOCATE ( idx, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "idx" )
    DEALLOCATE ( int_tmp, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "int_tmp" )
    DEALLOCATE ( real_tmp, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "real_tmp" )

  ELSE

    ! We do not sort non-spherical grids, but we have to make sure that
    ! possible reference grids are compatible
    ! We have to test if the new grid is larger then the reference grid
    IF ( PRESENT ( ref_grid ) ) THEN
      IF ( ALL ( pw_grid % npts >= ref_grid % npts ) ) THEN
        idx = 0
        ngr = SIZE ( ref_grid % gsq ( : ) )
        ! stop here as this code is not yet working
        CALL stop_program ( "pw_grid_sort", "Reference grid code for"//&
            " non-spherical cutoff not yet working")
      ELSE
        CALL stop_program ( "pw_grid_sort", "Reference grid has to"//&
            " be smaller then new grid for non-spherical cutoff")
      END IF
    END IF

  END IF

  ! check if ordering is compatible to reference grid
  IF ( PRESENT ( ref_grid ) ) THEN
    ngr = SIZE ( ref_grid % gsq ( : ) )
    ngr = MIN ( ng, ngr )
    IF ( .NOT. ALL ( pw_grid % g_hat ( 1:3, 1:ngr ) &
                  == ref_grid % g_hat ( 1:3, 1:ngr ) ) ) THEN
      CALL stop_program ( "pw_grid_sort", "G space sorting not compatible" )
    END IF
  END IF

END SUBROUTINE pw_grid_sort

SUBROUTINE sort_shells ( gsq, g_hat, idx )

  IMPLICIT NONE

! Argument
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: gsq
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: g_hat
  INTEGER, DIMENSION ( : ), INTENT ( INOUT ) :: idx

! Locals
  REAL ( dbl ), PARAMETER :: small = 1.e-12_dbl
  REAL ( dbl ) :: s_begin
  INTEGER :: ng, ig, s1, s2

  ng = SIZE ( gsq )
  s_begin = -1._dbl
  s1 = 0
  s2 = 0
  ig = 1
  DO ig = 1, ng
    IF ( ABS ( gsq ( ig ) - s_begin ) < small ) THEN
      s2 = ig
    ELSE
      CALL redist ( g_hat, idx, s1, s2)
      s_begin = gsq ( ig )
      s1 = ig
      s2 = ig
    END IF
  END DO
  CALL redist ( g_hat, idx, s1, s2 )

END SUBROUTINE sort_shells

SUBROUTINE redist ( g_hat, idx, s1, s2 )

  IMPLICIT NONE

! Argument
  INTEGER, DIMENSION ( :, : ), INTENT ( IN ) :: g_hat
  INTEGER, DIMENSION ( : ), INTENT ( INOUT ) :: idx
  INTEGER, INTENT ( IN ) :: s1, s2

! Locals
  INTEGER :: ns, info, i, ii
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: indl
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: slen

  IF ( s2 <= s1 ) RETURN
  ns = s2 - s1 + 1
  ALLOCATE ( indl ( ns ), STAT = info )
  IF ( info /= 0 ) call stop_memory ( "redist", "indl", ns )
  ALLOCATE ( slen ( ns ), STAT = info )
  IF ( info /= 0 ) call stop_memory ( "redist", "slen", ns )

  DO i = s1, s2
    ii = idx ( i )
    slen ( i - s1 + 1 ) = 1000._dbl * REAL(g_hat(1,ii),dbl) + &
       REAL(g_hat(2,ii),dbl) + 0.001_dbl * REAL(g_hat(3,ii),dbl)
  END DO
  CALL sort ( slen, ns, indl )
  DO i = 1, ns
    ii = indl ( i ) + s1 - 1
    indl ( i ) = idx ( ii )
  END DO
  idx ( s1:s2 ) = indl ( 1:ns )

  DEALLOCATE ( indl, STAT = info )
  IF ( info /= 0 ) call stop_memory ( "redist", "indl" )
  DEALLOCATE ( slen, STAT = info )
  IF ( info /= 0 ) call stop_memory ( "redist", "slen" )

END SUBROUTINE redist

!!*****
!****************************************************************************
!!****** pw_grids/pw_grid_remap [1.0] *
!!
!!   NAME
!!     pw_grid_remap
!!
!!   FUNCTION
!!     Reorder yzq and yzp arrays for parallel FFT according to FFT mapping
!!
!!   AUTHOR
!!     JGH (17-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_grid_remap ( pw_grid, yz )
  IMPLICIT NONE

! Argument
  TYPE ( pw_grid_type ), INTENT ( INOUT ), TARGET :: pw_grid
  INTEGER, DIMENSION ( :, : ), INTENT ( OUT ) :: yz

! Locals
  INTEGER :: i, j, is, ip, m, n, gpt, ny, nz
  INTEGER, DIMENSION ( : ), POINTER :: mapl, mapm, mapn

!------------------------------------------------------------------------------

  IF ( pw_grid % para % mode == 0 ) RETURN
  IF ( .NOT. pw_grid % para % ray_distribution ) RETURN

  ny = pw_grid % npts ( 2 )
  nz = pw_grid % npts ( 3 )

  yz = 0
  pw_grid % para % yzp = 0
  pw_grid % para % yzq = 0

  mapl => pw_grid % mapl % pos
  mapm => pw_grid % mapm % pos
  mapn => pw_grid % mapn % pos

  DO gpt = 1, SIZE ( pw_grid % gsq ( : ) )
    m = mapm ( pw_grid % g_hat ( 2, gpt ) ) + 1
    n = mapn ( pw_grid % g_hat ( 3, gpt ) ) + 1
    yz ( m, n ) = yz ( m, n ) + 1
  END DO
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     mapl => pw_grid % mapl % neg
     mapm => pw_grid % mapm % neg
     mapn => pw_grid % mapn % neg
     DO gpt = 1, SIZE ( pw_grid % gsq ( : ) )
       m = mapm ( pw_grid % g_hat ( 2, gpt ) ) + 1
       n = mapn ( pw_grid % g_hat ( 3, gpt ) ) + 1
       yz ( m, n ) = yz ( m, n ) + 1
     END DO
  END IF

  ip =  pw_grid % para % my_pos
  is = 0
  DO i = 1, nz
    DO j = 1, ny
      IF ( yz ( j, i ) > 0 ) THEN
        is = is + 1
        pw_grid % para % yzp ( 1, is, ip ) = j
        pw_grid % para % yzp ( 2, is, ip ) = i
        pw_grid % para % yzq ( j, i ) = is
      END IF
    END DO
  END DO
  IF ( is /= pw_grid % para % nyzray ( ip ) ) THEN
    CALL stop_program ( "pw_grid_remap", "recount of yz ray failed" )
  END IF
  CALL mp_sum ( pw_grid % para % yzp, pw_grid % para % group )

END SUBROUTINE pw_grid_remap

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
!!     JGH (8-Mar-2001) : also update information on volume and grid spaceing
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
  REAL ( dbl ) :: l, m, n
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: h_inv
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: g

!------------------------------------------------------------------------------

  h_inv => cell % h_inv
  g => pw_grid % g

  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts, dbl )
  pw_grid % dr ( : ) = SQRT ( SUM ( cell % hmat ( :, : ) ** 2, 1 ) ) &
       / REAL ( pw_grid % npts ( : ), dbl )

  DO gpt = 1, SIZE ( g ( 1, : ) )

     l = twopi * REAL ( pw_grid % g_hat ( 1, gpt ), dbl )
     m = twopi * REAL ( pw_grid % g_hat ( 2, gpt ), dbl )
     n = twopi * REAL ( pw_grid % g_hat ( 3, gpt ), dbl )
     
     g ( 1, gpt ) = l * h_inv(1,1) + m * h_inv(2,1) + n * h_inv(3,1)
     g ( 2, gpt ) = l * h_inv(1,2) + m * h_inv(2,2) + n * h_inv(3,2)
     g ( 3, gpt ) = l * h_inv(1,3) + m * h_inv(2,3) + n * h_inv(3,3)
     
     pw_grid % gsq ( gpt ) = g ( 1, gpt ) * g ( 1, gpt ) &
                           + g ( 2, gpt ) * g ( 2, gpt ) &
                           + g ( 3, gpt ) * g ( 3, gpt )

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
!!     This routine calculates the cutoff in momentum units!
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
       * REAL ( ( npts ( 1 ) - 1 ) / 2, dbl )
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = length

! compute 2*pi*h_inv^t*g  where g = (0,npts[2],0)
  gdum ( : ) = twopi * box % h_inv ( 2, : ) &
       * REAL ( ( npts ( 2 ) - 1 ) / 2, dbl )
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )

! compute 2*pi*h_inv^t*g  where g = (0,0,npts[3])
  gdum ( : ) = twopi * box % h_inv ( 3, : ) &
       * REAL ( ( npts ( 3 ) - 1 ) / 2, dbl )
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )

  cutoff = gcut - 1.0E-5_dbl

END SUBROUTINE pw_find_cutoff

!!*****
!******************************************************************************

END MODULE pw_grids

!******************************************************************************
