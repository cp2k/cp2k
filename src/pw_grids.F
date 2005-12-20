!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
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
!!     JGH (21-02-2003) : Generalized reference grid concept
!!
!!   NOTES
!!     pw_grid_create : set the defaults
!!     pw_grid_release : release all memory connected to type
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
  USE fft_tools,                       ONLY: FFT_RADIX_NEXT,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft_radix_operations
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_cart_coords,&
                                             mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_compare,&
                                             mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_dims_create,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sum,&
                                             MPI_COMM_SELF
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             map_pn,&
                                             pw_grid_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: dotprod_3d,&
                                             get_limit,&
                                             matvec_3x3,&
                                             sort
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_grid_create, pw_grid_compare
  PUBLIC :: pw_grid_setup, pw_grid_change, pw_find_cutoff, pw_grid_find_bounds
  PUBLIC :: pw_grid_retain, pw_grid_release, pw_grid_create_copy_no_pbc
  PUBLIC :: create_gvectors

  INTEGER :: grid_tag = 0
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_grids'
!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pw_grids/pw_grid_create [1.0] *
!!
!!   NAME
!!     pw_grid_create
!!
!!   FUNCTION
!!     Initialize a PW grid with all defaults
!!
!!   AUTHOR
!!     JGH (7-Feb-2001) & fawzi
!!
!!   MODIFICATION HISTORY
!!     JGH (21-Feb-2003) : initialize pw_grid%reference
!!
!!*** *************************************************************************

  SUBROUTINE pw_grid_create ( pw_grid, pe_group, error )

    TYPE(pw_grid_type), pointer        :: pw_grid
    INTEGER, intent(in) :: pe_group
    TYPE(cp_error_type), INTENT(inout), optional :: error

    LOGICAL :: failure
    integer :: stat
    CHARACTER(len=*), PARAMETER :: routineN="pw_grid_create",&
         routineP=moduleN//':'//routineN

    failure=.false.
    CPPrecondition(.NOT.ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(pw_grid,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       pw_grid % bounds = 0
       pw_grid % cutoff = 0.0_dp
       pw_grid % grid_span = FULLSPACE
       pw_grid % para % mode = PW_MODE_LOCAL
       pw_grid % para % rs_dims = 0
       pw_grid % reference = 0
       pw_grid % ref_count = 1
       NULLIFY ( pw_grid % g )
       NULLIFY ( pw_grid % gsq )
       NULLIFY ( pw_grid % g_hat )
       NULLIFY ( pw_grid % gidx )
       NULLIFY ( pw_grid % mapl % pos )
       NULLIFY ( pw_grid % mapl % neg )
       NULLIFY ( pw_grid % mapm % pos )
       NULLIFY ( pw_grid % mapm % neg )
       NULLIFY ( pw_grid % mapn % pos )
       NULLIFY ( pw_grid % mapn % neg )
       NULLIFY ( pw_grid % para % yzp )
       NULLIFY ( pw_grid % para % yzq )
       NULLIFY ( pw_grid % para % nyzray )
       NULLIFY ( pw_grid % para % bo )
       NULLIFY ( pw_grid % para % pos_of_x )

       ! assign a unique tag to this grid
       grid_tag = grid_tag + 1
       pw_grid %id_nr = grid_tag

       ! parallel info
       CALL mp_comm_dup ( pe_group, pw_grid % para % group )
       CALL mp_environ ( pw_grid % para % group_size, &
            pw_grid % para % my_pos, &
            pw_grid % para % group )
       pw_grid % para % group_head_id = 0
       pw_grid % para % group_head = &
            ( pw_grid % para % group_head_id == pw_grid % para % my_pos )
       IF (pw_grid % para % group_size > 1) THEN
          pw_grid % para % mode = PW_MODE_DISTRIBUTED
       ELSE
          pw_grid % para % mode = PW_MODE_LOCAL
       END IF
    END IF

  END SUBROUTINE pw_grid_create

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
!!*** *************************************************************************

FUNCTION pw_grid_compare ( grida, gridb ) RESULT ( equal )


    TYPE(pw_grid_type), INTENT(IN)           :: grida, gridb
    LOGICAL                                  :: equal

!------------------------------------------------------------------------------

  IF ( grida %id_nr == gridb %id_nr ) THEN
    equal = .TRUE.
  ELSE
! for the moment all grids with different identifiers are considered as not equal
! later we can get this more relaxed
    equal = .FALSE.
  END IF

END FUNCTION pw_grid_compare

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
!!     JGH (21-May-2002) : Remove orthorhombic keyword (code is fast enough)
!!     JGH (19-Feb-2003) : Negative cutoff can be used for non-spheric grids
!!     Joost VandeVondele (Feb-2004) : optionally generate pw grids that are commensurate in rs
!!
!!*** *************************************************************************

SUBROUTINE pw_grid_setup ( cell, pw_grid, cutoff, info, fft_usage, &
                           symm_usage, blocked, ref_grid, ncommensurate,&
                           icommensurate, rs_dims )


    TYPE(cell_type), INTENT(IN)              :: cell
    TYPE(pw_grid_type), POINTER              :: pw_grid
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cutoff
    INTEGER, INTENT(IN), OPTIONAL            :: info
    LOGICAL, INTENT(IN), OPTIONAL            :: fft_usage, symm_usage, blocked
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid
    INTEGER, INTENT(IN), OPTIONAL            :: ncommensurate, icommensurate
    INTEGER, DIMENSION(2), INTENT(in), optional :: rs_dims

    INTEGER :: allocstat, handle, i, ires, my_icommensurate, &
      my_ncommensurate, n(3), nlowest, nlowest_new, ntest(3)
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: yz_mask
    LOGICAL                                  :: blocking, fft, symmetry
    REAL(KIND=dp)                            :: ecut, rv(3,3)

! cutoff could be inout, and set to the effective cutoff of the given grid
! in atomic units
! parent group for this grid
! output unit for information on grid
! has the grid size to be 
! compatible with the FFT
! has the grid size to be 
! symmetric (g <-> -g)    
! block or ray-distribution
!------------------------------------------------------------------------------

  CALL timeset("pw_grid_setup","I","",handle)

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

  ! default is to allow block usage
  IF ( PRESENT ( blocked ) ) THEN
    blocking = blocked
  ELSE
    blocking = .TRUE.
  END IF
 
  ! ncommensurate is the number of commensurate grids 
  ! in order to have non-commensurate grids ncommensurate must be 0
  ! icommensurte  is the level number of communensurate grids
  ! this implies that the number of grid points in each direction
  ! is k*2**(ncommensurate-icommensurate) 
  IF ( PRESENT ( ncommensurate ) .AND. PRESENT( icommensurate ) ) THEN
     my_ncommensurate=ncommensurate
     IF (my_ncommensurate .GT. 0 ) THEN
         my_icommensurate=icommensurate
     ELSE
         my_icommensurate=0
     ENDIF
     IF (my_icommensurate > my_ncommensurate ) THEN
        CALL stop_program ( "grid_setup", &
                            "my_icommensurate > my_ncommensurate" )
     ENDIF
     IF (my_icommensurate<=0 .AND. my_ncommensurate > 1) THEN
        CALL stop_program ( "grid_setup", &
                            " my_incommensurate<=0 .AND. my_ncommensurate > 1 " )
     ENDIF
     IF (my_ncommensurate < 0 ) THEN
        CALL stop_program ( "grid_setup", &
                            "my_ncommensurate < 0 " )
     ENDIF
     IF (my_icommensurate > 1 .AND. .NOT. PRESENT(ref_grid)  ) THEN
        CALL stop_program ( "grid_setup", &
                            "my_icommensurate > 1 .AND. .NOT. PRESENT(ref_grid)" )
     ENDIF
  ELSE
     IF (PRESENT( ncommensurate ) .OR. PRESENT( icommensurate ) ) THEN
        CALL stop_program ( "grid_setup", &
                            "both ncommensurate and icommensurate needed" )
     ELSE
        my_ncommensurate=0
        my_icommensurate=0
     ENDIF
  ENDIF

  ! set pointer to possible reference grid
  IF ( PRESENT ( ref_grid ) ) THEN
    pw_grid % reference = ref_grid %id_nr
  END IF

  IF ( PRESENT ( cutoff ) ) THEN
     pw_grid % cutoff = ABS ( cutoff )
     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) &
          CALL pw_grid_find_bounds ( pw_grid % bounds, cell % h_inv, &
                                     pw_grid % cutoff, fft )
     IF ( cutoff < 0.0_dp ) THEN
        pw_grid % spherical = .FALSE.
        ecut = 1.e10_dp
     ELSE
        pw_grid % spherical = .TRUE.
        ecut = cutoff
     END IF
  ELSE
     ecut = 1.e10_dp              ! all g-vectors in the box will be included
     pw_grid % spherical = .FALSE.
     pw_grid % cutoff = 0.0_dp
  END IF

  IF ( my_ncommensurate > 0 ) THEN
     IF (my_icommensurate > 1) THEN ! we should use the reference grid to find out the number of gridpoints
         pw_grid % bounds ( 1, : ) = - ref_grid % npts  / 2**my_icommensurate
         pw_grid % bounds ( 2, : ) = pw_grid % bounds ( 1, : ) +  ref_grid % npts / (2**(my_icommensurate-1)) - 1
     ENDIF
  ENDIF

  IF ( .NOT. pw_grid % spherical ) THEN

     IF ( SUM ( ABS ( pw_grid % bounds ( :, : ) ) ) == 0 ) THEN
        CALL stop_program ( "grid_setup", &
                            "provide initial values for bounds" )
     END IF

     ntest = pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1

     IF ( fft ) THEN
 
         ! select FFT compatible bounds
         ! without a cutoff and HALFSPACE we have to be sure that there is
         ! a negative counterpart to every g vector (-> odd number of grid points)
         IF ( pw_grid % grid_span == HALFSPACE .AND. symmetry ) THEN
  
           CALL fft_radix_operations ( ntest(1), n(1), FFT_RADIX_NEXT_ODD )
           CALL fft_radix_operations ( ntest(2), n(2), FFT_RADIX_NEXT_ODD )
           CALL fft_radix_operations ( ntest(3), n(3), FFT_RADIX_NEXT_ODD )
  
         ELSE
           ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( ntest(1), n(1), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=n(1)/2**(my_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
             IF (nlowest==nlowest_new .AND. MODULO(n(1),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                ntest(1)=n(1)+1
             ENDIF
           ENDDO
           ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( ntest(2), n(2), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=n(2)/2**(my_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
             IF (nlowest==nlowest_new .AND. MODULO(n(2),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                ntest(2)=n(2)+1
             ENDIF
           ENDDO
           ! keep looping to find the right one
           DO
             CALL fft_radix_operations ( ntest(3), n(3), FFT_RADIX_NEXT )
             ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
             nlowest=n(3)/2**(my_ncommensurate-my_icommensurate)
             CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
             IF (nlowest==nlowest_new .AND. MODULO(n(3),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                EXIT
             ELSE
                ntest(3)=n(3)+1
             ENDIF
           ENDDO

         END IF

     ELSE

       ! without a cutoff and HALFSPACE we have to be sure that there is
       ! a negative counterpart to every g vector (-> odd number of grid points)
       IF ( pw_grid % grid_span == HALFSPACE .AND. symmetry ) &
            n = ntest + MOD ( ntest + 1, 2 )

     END IF

     pw_grid % bounds ( 1, : ) = - n / 2
     pw_grid % bounds ( 2, : ) = pw_grid % bounds ( 1, : ) + n - 1

  END IF

  pw_grid % npts ( : ) = &
       pw_grid % bounds ( 2, : ) - pw_grid % bounds ( 1, : ) + 1

  ! final check if all went fine ...
  IF (my_ncommensurate>0) THEN
     IF ( ANY( MODULO(pw_grid % npts,2**(my_ncommensurate-my_icommensurate)).NE.0 ) ) THEN ! nope, sorry
        CALL stop_program ( "grid_setup", &
            "commensurate option failed (I) ... maybe not yet programmed for this combination of options ?" )
     END IF
     IF ( my_icommensurate > 1 ) THEN 
        IF ( ANY(pw_grid % npts * 2 ** (my_icommensurate-1) .NE. ref_grid % npts ) ) THEN
          CALL stop_program ( "grid_setup", &
              "commensurate option failed (II) ... maybe not yet programmed for this combination of options ?" )
        ENDIF
     END IF
  ENDIF

  n ( : )  = pw_grid % npts ( : )

  ! Find the number of grid points
  ! yz_mask counts the number of g-vectors orthogonal to the yz plane
  ! the indices in yz_mask are from -n/2 .. n/2 shifted by n/2 + 1
  ! these are not mapped indices !
  ALLOCATE ( yz_mask ( n(2), n(3) ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_setup", "yz_mask", 0 )
  CALL pw_grid_count ( cell % h_inv, pw_grid, ecut, yz_mask )

  ! Check if reference grid is compatible
  IF ( PRESENT ( ref_grid ) ) THEN
    IF ( pw_grid % para % mode /= ref_grid % para % mode ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible parallelisation scheme" )
    END IF
    IF ( pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
      CALL mp_comm_compare(pw_grid % para % group, ref_grid % para % group, ires )
      IF ( ires >2 ) & !FM make it >3 ?
        CALL stop_program ( "pw_grid_setup", "Incompatible MPI groups" )
    END IF
    IF ( pw_grid % grid_span /= ref_grid % grid_span ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible grid types" )
    END IF
    IF ( pw_grid % spherical .NEQV. ref_grid % spherical ) THEN
      CALL stop_program ( "pw_grid_setup", "Incompatible cutoff schemes" )
    END IF
  END IF

  ! Distribute grid
  CALL pw_grid_distribute ( pw_grid, yz_mask, ref_grid, blocking, &
       rs_dims=rs_dims )

  ! Allocate the grid fields
  CALL pw_grid_allocate ( pw_grid, pw_grid % ngpts_cut_local, &
                          pw_grid % bounds )

  ! Fill in the grid structure
  CALL pw_grid_assign ( cell % h_inv, pw_grid, ecut )
  
  ! Sort g vector wrt length (only local for each processor)
  CALL pw_grid_sort ( pw_grid, ref_grid )
  
  CALL pw_grid_remap ( pw_grid, yz_mask )

  DEALLOCATE ( yz_mask )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_setup", "yz_mask" )

  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts,KIND=dp)
  pw_grid % dr ( 1 ) = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 1 ),KIND=dp)
  pw_grid % dr ( 2 ) = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 2 ),KIND=dp)
  pw_grid % dr ( 3 ) = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 3 ),KIND=dp)

!
! Output: All the information of this grid type
!

  IF ( PRESENT ( info ) ) THEN
    IF ( pw_grid % para % mode == PW_MODE_LOCAL) THEN
       IF (info >= 0 ) THEN
          WRITE ( info, '(/,A,T71,I10)' ) &
               " PW_GRID: Information for grid number ", pw_grid %id_nr
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
       END IF
    ELSE

      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_sum ( n(1:2), pw_grid % para % group )
      n ( 3 ) = SUM ( pw_grid % para % nyzray )
      rv ( :, 1 ) = REAL ( n,KIND=dp) / REAL ( pw_grid % para % group_size,KIND=dp)
      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_max ( n(1:2), pw_grid % para % group )
      n ( 3 ) = MAXVAL ( pw_grid % para % nyzray )
      rv ( :, 2 ) = REAL ( n,KIND=dp)
      n ( 1 ) = pw_grid % ngpts_cut_local
      n ( 2 ) = pw_grid % ngpts_local
      CALL mp_min ( n(1:2), pw_grid % para % group )
      n ( 3 ) = MINVAL ( pw_grid % para % nyzray )
      rv ( :, 3 ) = REAL ( n,KIND=dp)

      IF ( pw_grid % para % group_head .AND. info>=0) THEN
        WRITE ( info, '(/,A,T71,I10)' ) &
          " PW_GRID: Information for grid number ", pw_grid %id_nr
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
        WRITE ( info, '(A,T48,A)' ) " PW_GRID:   Distribution", &
             "  Average         Max         Min"
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   G-Vectors", &
             rv ( 1, 1 ), NINT ( rv ( 1, 2 ) ), NINT ( rv ( 1, 3 ) )
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   G-Rays   ", &
             rv ( 3, 1 ), NINT ( rv ( 3, 2 ) ), NINT ( rv ( 3, 3 ) )
        WRITE ( info, '(A,T45,F12.1,2I12)' ) " PW_GRID:   Real Space Points", &
             rv ( 2, 1 ), NINT ( rv ( 2, 2 ) ), NINT ( rv ( 2, 3 ) )
      END IF ! group head
    END IF ! local
  END IF ! present info

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_setup

!******************************************************************************

 SUBROUTINE create_gvectors(pw_grid,cell,ecut,blocking,ref_grid)
  ! Find the number of grid points
  ! yz_mask counts the number of g-vectors orthogonal to the yz plane
  ! the indices in yz_mask are from -n/2 .. n/2 shifted by n/2 + 1
  ! these are not mapped indices !

    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp)                            :: ecut
    LOGICAL :: blocking
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid

    INTEGER :: istat, n(3)
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: yz_mask

    n ( : )  = pw_grid % npts ( : )

    ALLOCATE ( yz_mask ( n(2), n(3) ), STAT = istat )
    IF ( istat /= 0 ) CALL stop_memory ( "create_gvectors", "yz_mask", 0 )
    CALL pw_grid_count ( cell % h_inv, pw_grid, ecut, yz_mask )

    ! Distribute grid
    CALL pw_grid_distribute ( pw_grid, yz_mask, ref_grid, blocking )

    ! Allocate the grid fields
    CALL pw_grid_allocate ( pw_grid, pw_grid % ngpts_cut_local, &
                            pw_grid % bounds )

    ! Fill in the grid structure
    CALL pw_grid_assign ( cell % h_inv, pw_grid, ecut )
  
    ! Sort g vector wrt length (only local for each processor)
    CALL pw_grid_sort ( pw_grid, ref_grid )
  
    CALL pw_grid_remap ( pw_grid, yz_mask )
 
    DEALLOCATE ( yz_mask, STAT=istat )
    IF ( istat /= 0 ) CALL stop_memory ( "create_gvectors", "yz_mask" )

  END SUBROUTINE create_gvectors


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
!!     JGH (22-May-2002) bug fix for pre_tag and HALFSPACE grids
!!     JGH (09-Sep-2003) reduce scaling for distribution
!!
!!*** *************************************************************************

SUBROUTINE pw_grid_distribute ( pw_grid, yz_mask, ref_grid, blocking, rs_dims,&
     error)


    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, DIMENSION(:, :), INTENT(INOUT)  :: yz_mask
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid
    LOGICAL, INTENT(IN)                      :: blocking
    INTEGER, DIMENSION(2), INTENT(in), optional :: rs_dims
    TYPE(cp_error_type), INTENT(inout), optional :: error

    INTEGER                                  :: coor( 2 ), gmax, handle, i, &
                                                i1, i2, ierr, ip, ipl, j, k, &
                                                l, lby, lbz, lo( 2 ), m, n, &
                                                np, ns, nx, ny, nz, rsd( 2 )
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: yz_index
    REAL(KIND=dp)                            :: tfun, tt
    CHARACTER(len=*), PARAMETER :: routineN="pw_grid_distribute",&
         routineP=moduleN//':'//routineN
    logical :: failure

!------------------------------------------------------------------------------

  CALL timeset("pw_grid_distribute","I","",handle)

  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )

  pw_grid % ngpts = PRODUCT ( pw_grid % npts )
  CPPrecondition(ALL(pw_grid%para%rs_dims==0),cp_failure_level,routineP,error,failure)
  IF (PRESENT(rs_dims)) THEN
     pw_grid%para%rs_dims=rs_dims
  END IF

  IF ( pw_grid % para % mode == PW_MODE_LOCAL) THEN

    pw_grid % bounds_local = pw_grid % bounds
    pw_grid % npts_local = pw_grid % npts
    pw_grid % ngpts_cut_local = pw_grid % ngpts_cut
    pw_grid % ngpts_local = PRODUCT ( pw_grid % npts_local )
    pw_grid % para % rs_dims=1
    CALL mp_cart_create ( MPI_COMM_SELF, 2, &
                          pw_grid % para % rs_dims, &
                          pw_grid % para % rs_pos, &
                          pw_grid % para % rs_group )
    CALL mp_cart_rank ( pw_grid % para % rs_group, &
                        pw_grid % para % rs_pos, &
                        pw_grid % para % rs_mpo )

  ELSE

!..find the real space distribution
    nx = pw_grid % npts ( 1 )
    ny = pw_grid % npts ( 2 )
    nz = pw_grid % npts ( 3 )

    np = pw_grid % para % group_size

    IF ( PRODUCT ( pw_grid % para % rs_dims ) == 0 ) THEN

      ns = INT ( SQRT ( REAL ( np,KIND=dp) ) )

      tfun = 1.e20_dp
      DO i = ns, 2, -1
        IF ( MOD ( np, i ) == 0 ) THEN
          j = np/i
          k = nx*ny - i*(nx/i) * j*(ny/j)
          IF ( k == 0 ) THEN
            tfun = 0.0_dp
            m = i
          ELSE
            tt = 1.0_dp - REAL ( k,KIND=dp) / REAL ( np,KIND=dp)
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
        tt = 1.0_dp - REAL ( k,KIND=dp) / REAL ( np,KIND=dp)
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
      pw_grid % para % bo ( 1:2, 1, ip, 1 ) = get_limit (nx,rsd(1),coor(1))
      pw_grid % para % bo ( 1:2, 2, ip, 1 ) = get_limit (ny,rsd(2),coor(2))
      pw_grid % para % bo ( 1, 3, ip, 1 ) = 1
      pw_grid % para % bo ( 2, 3, ip, 1 ) = nz
      ! distribution xYz
      pw_grid % para % bo ( 1:2, 1, ip, 2 ) = get_limit (nx,rsd(1),coor(1))
      pw_grid % para % bo ( 1, 2, ip, 2 ) = 1
      pw_grid % para % bo ( 2, 2, ip, 2 ) = ny
      pw_grid % para % bo ( 1:2, 3, ip, 2 ) = get_limit (nz,rsd(2),coor(2))
      ! distribution Xyz
      pw_grid % para % bo ( 1, 1, ip, 3 ) = 1
      pw_grid % para % bo ( 2, 1, ip, 3 ) = nx
      pw_grid % para % bo ( 1:2, 2, ip, 3 ) = get_limit (ny,rsd(1),coor(1))
      pw_grid % para % bo ( 1:2, 3, ip, 3 ) = get_limit (nz,rsd(2),coor(2))
    END DO

!..find the g space distribution
    pw_grid % ngpts_cut_local = 0

    ALLOCATE ( pw_grid % para % nyzray ( 0: np-1 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
          "pw_grid % para % nyzray", np )

    ALLOCATE ( pw_grid % para % yzq ( ny, nz ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
       "pw_grid % para % yzq", ny*nz )

    IF ( pw_grid % spherical .OR. pw_grid % grid_span == HALFSPACE &
         .OR. .NOT. blocking ) THEN

      pw_grid % para % ray_distribution = .TRUE.

      pw_grid % para % yzq = -1
      IF ( PRESENT ( ref_grid ) ) THEN
        ! tag all vectors from the reference grid
        CALL pre_tag ( pw_grid, yz_mask, ref_grid )
      END IF

      ! Round Robin distribution 
      ! Processors 0 .. NP-1, NP-1 .. 0  get the largest remaining batch
      ! of g vectors in turn
  
      i1 = SIZE ( yz_mask, 1 )
      i2 = SIZE ( yz_mask, 2 )
      ALLOCATE ( yz_index(2,i1*i2), STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "pw_grid_distribute", &
                                          "yz_index", 2*i1*i2 )
      CALL order_mask ( yz_mask, yz_index )
      DO i = 1, i1*i2
        lo(1) = yz_index(1,i)
        lo(2) = yz_index(2,i)
        IF ( lo(1)*lo(2) == 0 ) CYCLE
        gmax = yz_mask ( lo(1), lo(2) )
        IF ( gmax == 0 ) CYCLE
        yz_mask ( lo(1), lo(2) ) = 0
        ip = MOD ( i-1, 2*np )
        IF ( ip > np - 1 ) ip = 2*np - ip - 1
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
      DEALLOCATE ( yz_index, STAT = ierr )
      IF ( ierr /= 0 ) CALL stop_memory ( "", "yz_index" )

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

  ! pos_of_x(i) tells on which cpu pw%cr3d(i,:,:) is located
  ! should be computable in principle, without the need for communication
  IF (pw_grid % para % mode .EQ. PW_MODE_DISTRIBUTED) THEN
    ALLOCATE(pw_grid % para % pos_of_x( pw_grid % bounds(1,1): pw_grid % bounds(2,1) ))
    pw_grid % para % pos_of_x = 0
    pw_grid % para % pos_of_x(pw_grid % bounds_local(1,1) : pw_grid % bounds_local(2,1))=pw_grid % para % my_pos
    CALL mp_sum( pw_grid % para % pos_of_x, pw_grid % para % group )
  ELSE
    ! this should not be needed
    ALLOCATE(pw_grid % para % pos_of_x( pw_grid % bounds(1,1): pw_grid % bounds(2,1) ))
    pw_grid % para % pos_of_x = 0
  ENDIF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_distribute

!******************************************************************************

SUBROUTINE pre_tag ( pw_grid, yz_mask, ref_grid )


    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, DIMENSION(:, :), INTENT(INOUT)  :: yz_mask
    TYPE(pw_grid_type), INTENT(IN)           :: ref_grid

    INTEGER                                  :: gmax, ig, ip, lby, lbz, my, &
                                                mz, ny, nz, uby, ubz, y, yp, &
                                                z, zp

!------------------------------------------------------------------------------

  ny = ref_grid % npts ( 2 )
  nz = ref_grid % npts ( 3 )
  lby = pw_grid % bounds ( 1, 2 )
  lbz = pw_grid % bounds ( 1, 3 )
  uby = pw_grid % bounds ( 2, 2 )
  ubz = pw_grid % bounds ( 2, 3 )
  my = SIZE ( yz_mask, 1 )
  mz = SIZE ( yz_mask, 2 )
 
  ! loop over all processors and all g vectors yz lines on this processor
  DO ip = 0, ref_grid % para % group_size - 1
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
      ! this tag is outside the cutoff range of the new grid
      IF ( pw_grid % grid_span == HALFSPACE ) THEN
        yp = -y - 2*lby + 2
        zp = -z - 2*lbz + 2
        ! if the referenz grid is larger than the mirror point may be
        ! outside the new grid even if the original point is inside
        IF ( yp < 1 .OR. yp > my .OR. zp < 1 .OR. zp > mz ) CYCLE
        gmax = MAX ( yz_mask ( y, z ), yz_mask ( yp, zp ) )
        IF ( gmax == 0 ) CYCLE
        yz_mask ( y, z ) = 0
        yz_mask ( yp, zp ) = 0
        pw_grid % para % yzq ( y, z ) = ip
        pw_grid % para % yzq ( yp, zp ) = ip
      ELSE 
        gmax = yz_mask ( y, z )
        IF ( gmax == 0 ) CYCLE
        yz_mask ( y, z ) = 0
        pw_grid % para % yzq ( y, z ) = ip
      END IF
      IF ( ip == pw_grid % para % my_pos ) THEN
         pw_grid % ngpts_cut_local = pw_grid % ngpts_cut_local + gmax
      END IF
    END DO
  END DO

END SUBROUTINE pre_tag

!------------------------------------------------------------------------------

SUBROUTINE order_mask ( yz_mask, yz_index )


    INTEGER, DIMENSION(:, :), INTENT(IN)     :: yz_mask
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: yz_index

    INTEGER                                  :: i, i1, i2, icount, ierr, j, &
                                                j1, j2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cindex, icol, irow, rindex

!------------------------------------------------------------------------------

  i1 = SIZE ( yz_mask, 1 )
  i2 = SIZE ( yz_mask, 2 )
  ALLOCATE ( irow(i1), rindex(i1), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "order_mask", "irow", 2*i1 )
  ALLOCATE ( icol(i2), cindex(i2), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "order_mask", "icol", 2*i2 )

  yz_index = 0
  DO i = 1, i1
    irow(i) = SUM(yz_mask(i,:))
  END DO
  CALL sort ( irow, i1, rindex )
  icount = 0
  DO i = i1, 1, -1
    j1 = rindex ( i )
    icol = yz_mask(i,:)
    CALL sort ( icol, i2, cindex )
    DO j = i2, 1, -1
      j2 = cindex ( j )
      IF ( yz_mask(j1,j2) /= 0 ) THEN
        icount = icount + 1
        yz_index(1,icount) = j1
        yz_index(2,icount) = j2
      END IF
    END DO
  END DO

  j=HUGE(j)
  DO i=1,icount
    j1=yz_index(1,icount)
    j2=yz_index(2,icount)
    IF(j < yz_mask(j1,j2)) STOP
    j = yz_mask(j1,j2)
  ENDDO

  DEALLOCATE ( irow, rindex, icol, cindex, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "order_mask", "irow" )

END SUBROUTINE order_mask

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
!!*** *************************************************************************

SUBROUTINE pw_grid_count ( h_inv, pw_grid, cutoff, yz_mask )


    REAL(KIND=dp), DIMENSION(3, 3)           :: h_inv
    TYPE(pw_grid_type),   POINTER            :: pw_grid
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: yz_mask

    INTEGER                                  :: gpt, l, m, mm, n, &
                                                n_upperlimit, nlim( 2 ), nn
    INTEGER, DIMENSION(:, :), POINTER        :: bounds
    REAL(KIND=dp)                            :: ggi( 3), gi( 3 ), &
                                                gmat( 3, 3 ), length

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
  IF ( pw_grid % para % mode == PW_MODE_LOCAL ) THEN
    nlim ( 1 ) = bounds ( 1, 3 )
    nlim ( 2 ) = n_upperlimit
  ELSE IF ( pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
    n = n_upperlimit - bounds ( 1, 3 ) + 1
    nlim = get_limit ( n, pw_grid % para % group_size, pw_grid % para % my_pos )
    nlim = nlim + bounds ( 1, 3 ) - 1
  ELSE
     CALL stop_program ( "pw_grid_count", "para % mode not specified" )
  END IF

  yz_mask = 0
  DO n = nlim ( 1 ), nlim ( 2 )
     gi ( 3 ) = REAL(n,dp)
     nn = n - bounds ( 1, 3) + 1
     DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
        gi ( 2 ) = REAL(m,dp)
        mm = m - bounds ( 1, 2) + 1
        DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
           IF ( pw_grid % grid_span == HALFSPACE .AND. n == 0 ) THEN
              IF ( ( m == 0 .AND. l > 0 ) .OR. ( m > 0 ) ) CYCLE
           END IF

           gi ( 1 ) = REAL(l,dp)
           ggi = MATMUL ( gmat, gi )
           length = twopi * twopi * DOT_PRODUCT ( ggi , gi )
           IF ( 0.5_dp * length <= cutoff ) THEN
             gpt = gpt + 1
             yz_mask ( mm, nn ) = yz_mask ( mm, nn ) + 1
           END IF

        END DO
     END DO
  END DO

! number of g-vectors for grid
  IF ( pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
    CALL mp_sum ( gpt, pw_grid % para % group )
    CALL mp_sum ( yz_mask, pw_grid % para % group )
  ENDIF
  pw_grid % ngpts_cut = gpt

END SUBROUTINE pw_grid_count

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
!!*** *************************************************************************

SUBROUTINE pw_grid_assign ( h_inv, pw_grid, cutoff )


    REAL(KIND=dp), DIMENSION(3, 3)           :: h_inv
    TYPE(pw_grid_type), POINTER              :: pw_grid
    REAL(KIND=dp), INTENT(IN)                :: cutoff

    INTEGER                                  :: gpt, handle, i, ip, l, lby, &
                                                lbz, ll, m, mm, n, &
                                                n_upperlimit, nn
    INTEGER, DIMENSION(2, 3)                 :: bol, bounds
    REAL(KIND=dp)                            :: length, length_x, length_y, &
                                                length_z

!------------------------------------------------------------------------------

  CALL timeset("pw_grid_assign","I","",handle)

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
  IF ( pw_grid % para % mode == PW_MODE_LOCAL ) THEN
    gpt = 0
    DO n = bounds ( 1, 3 ), n_upperlimit
       DO m = bounds ( 1, 2 ), bounds ( 2, 2 )
          DO l = bounds ( 1, 1 ), bounds ( 2, 1 )
             IF ( pw_grid % grid_span == HALFSPACE .AND. n == 0 ) THEN
                IF ( ( m == 0 .AND. l > 0 ) .OR. ( m > 0 ) ) CYCLE
             END IF

             length_x &
                  = REAL(l,dp) * h_inv(1,1) &
                  + REAL(m,dp) * h_inv(2,1) &
                  + REAL(n,dp) * h_inv(3,1)
             length_y &
                  = REAL(l,dp) * h_inv(1,2) &
                  + REAL(m,dp) * h_inv(2,2) &
                  + REAL(n,dp) * h_inv(3,2)
             length_z &
                  = REAL(l,dp) * h_inv(1,3) &
                  + REAL(m,dp) * h_inv(2,3) &
                  + REAL(n,dp) * h_inv(3,3)

             length = length_x ** 2 + length_y ** 2 + length_z ** 2
             length = twopi * twopi * length

             IF ( 0.5_dp * length <= cutoff ) THEN
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
                 = REAL(l,dp) * h_inv(1,1) &
                 + REAL(m,dp) * h_inv(2,1) &
                 + REAL(n,dp) * h_inv(3,1)
            length_y &
                 = REAL(l,dp) * h_inv(1,2) &
                 + REAL(m,dp) * h_inv(2,2) &
                 + REAL(n,dp) * h_inv(3,2)
            length_z &
                 = REAL(l,dp) * h_inv(1,3) &
                 + REAL(m,dp) * h_inv(2,3) &
                 + REAL(n,dp) * h_inv(3,3)
            
            length = length_x ** 2 + length_y ** 2 + length_z ** 2
            length = twopi * twopi * length
          
            IF ( 0.5_dp * length <= cutoff ) THEN
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
                   = REAL(l,dp) * h_inv(1,1) &
                   + REAL(m,dp) * h_inv(2,1) &
                   + REAL(n,dp) * h_inv(3,1)
               length_y &
                   = REAL(l,dp) * h_inv(1,2) &
                   + REAL(m,dp) * h_inv(2,2) &
                   + REAL(n,dp) * h_inv(3,2)
               length_z &
                   = REAL(l,dp) * h_inv(1,3) &
                   + REAL(m,dp) * h_inv(2,3) &
                   + REAL(n,dp) * h_inv(3,3)
 
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
  IF ( pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
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

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_assign

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
!!*** *************************************************************************

SUBROUTINE pw_grid_set_maps ( grid_span, g_hat, mapl, mapm, mapn, npts )


    INTEGER, INTENT(IN)                      :: grid_span
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: g_hat
    TYPE(map_pn), INTENT(INOUT)              :: mapl, mapm, mapn
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts

    INTEGER                                  :: gpt, l, m, n, ng

!------------------------------------------------------------------------------

  ng = SIZE ( g_hat ,2 )

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
!!     JGH (21-05-2002) : Optimise code, remove orthorhombic special case
!!
!!*** *************************************************************************

SUBROUTINE pw_grid_find_bounds ( bounds, h_inv, cutoff, fft )


    INTEGER, DIMENSION(2, 3), INTENT(OUT)    :: bounds
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: h_inv
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN)                      :: fft

    INTEGER, PARAMETER                       :: lwork = 20

    INTEGER                                  :: handle, info, rin(3), rout(3)
    REAL(KIND=dp)                            :: eig( 3 ), gmat( 3, 3 ), &
                                                work( lwork )

!------------------------------------------------------------------------------

  CALL timeset("pw_grid_find_bounds","I","",handle)
  bounds ( :, : ) = 0

  gmat = MATMUL ( h_inv, TRANSPOSE ( h_inv ) )
  CALL dsyev("V","U",3,gmat,3,eig,work,lwork,info)
  IF ( info /= 0 ) CALL stop_program ( "pw_grid_find_bounds", "dsyev:info" )

  eig = SQRT ( 2.0_dp * cutoff / eig ) / twopi
  bounds ( 2, : ) = FLOOR ( matvec_3x3 ( gmat, eig ) )

  rin(1) = 2 * bounds ( 2, 1 ) + 1
  rin(2) = 2 * bounds ( 2, 2 ) + 1
  rin(3) = 2 * bounds ( 2, 3 ) + 1
  IF ( fft ) THEN
    CALL fft_radix_operations ( rin(1), rout(1), FFT_RADIX_NEXT )
    bounds ( 1, 1 ) = -rout(1)/2
    CALL fft_radix_operations ( rin(2), rout(2), FFT_RADIX_NEXT )
    bounds ( 1, 2 ) = -rout(2)/2
    CALL fft_radix_operations ( rin(3), rout(3), FFT_RADIX_NEXT )
    bounds ( 1, 3 ) = -rout(3)/2
  ELSE
    rout = rin
    bounds ( 1, : ) = -rout(:)/2
  END IF

  bounds ( 2, : ) = bounds ( 1, : ) + rout - 1

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_find_bounds

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
!!*** *************************************************************************

SUBROUTINE pw_grid_allocate ( pw_grid, ng, bounds )

! Argument
    TYPE(pw_grid_type), INTENT(INOUT)        :: pw_grid
    INTEGER, INTENT(IN)                      :: ng
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: bounds

    INTEGER                                  :: allocstat

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
!!     JGH (19-02-2003) : Order all grids, this makes subgrids also for
!!                        non-spherical cutoffs possible
!!     JGH (21-02-2003) : Introduce gather array for general reference grids
!!
!!*** *************************************************************************

SUBROUTINE pw_grid_sort ( pw_grid, ref_grid )


! Argument
    TYPE(pw_grid_type), INTENT(INOUT)        :: pw_grid
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid

    INTEGER                                  :: allocstat, handle, ig, ih, &
                                                ip, is, it, ng, ngr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: idx, int_tmp
    LOGICAL                                  :: g_found
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: real_tmp

!------------------------------------------------------------------------------

  CALL timeset("pw_grid_sort","I","",handle)

  ng = SIZE ( pw_grid % gsq  )
  ALLOCATE ( idx ( ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "idx", ng )
  ALLOCATE ( int_tmp ( ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "int_tmp", ng )
  ALLOCATE ( real_tmp ( ng ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "real_tmp", ng )

  ! grids are (locally) ordered by length of G-vectors
  CALL sort ( pw_grid % gsq, ng, idx )
  ! within shells order wrt x,y,z
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

  ! check if ordering is compatible to reference grid
  IF ( PRESENT ( ref_grid ) ) THEN
    ngr = SIZE ( ref_grid % gsq  )
    ngr = MIN ( ng, ngr )
    IF ( pw_grid % spherical ) THEN
      IF ( .NOT. ALL ( pw_grid % g_hat ( 1:3, 1:ngr ) &
                  == ref_grid % g_hat ( 1:3, 1:ngr ) ) ) THEN
        CALL stop_program ( "pw_grid_sort", "G space sorting not compatible" )
      END IF
    ELSE
      ALLOCATE ( pw_grid%gidx(1:ngr), STAT = allocstat )
      IF ( allocstat /= 0 ) CALL stop_memory ( "pw_grid_sort", "pw_grid%gidx", ngr )
      pw_grid%gidx = 0
      ! first try as many trivial associations as possible
      it = 0
      DO ig = 1, ngr
        IF ( .NOT. ALL ( pw_grid % g_hat ( 1:3, ig ) &
                  == ref_grid % g_hat ( 1:3, ig ) ) ) EXIT
        pw_grid%gidx(ig) = ig
        it = ig
      END DO
      ! now for the ones that could not be done
      IF ( ng == ngr ) THEN
        ! for the case pw_grid <= ref_grid
        is = it
        DO ig = it + 1, ngr
          g_found=.FALSE.
          DO ih = is + 1, SIZE ( ref_grid % gsq  )
            IF ( ABS ( pw_grid % gsq(ig) - ref_grid % gsq(ih) ) > 1.e-12_dp ) CYCLE
            g_found=.TRUE.
            EXIT
          END DO
          IF ( .NOT. g_found ) THEN
             WRITE(*,"(A,I10,F20.10)") "G-vector", ig, pw_grid % gsq(ig)
             CALL stop_program ( "pw_grid_sort", "G vector not found" )
          END IF
          ip = ih - 1
          DO ih = ip + 1, SIZE ( ref_grid % gsq  )
            IF ( ABS ( pw_grid % gsq(ig) - ref_grid % gsq(ih) ) > 1.e-12_dp ) CYCLE
            IF ( pw_grid % g_hat(1,ig) /= ref_grid % g_hat(1,ih) ) CYCLE
            IF ( pw_grid % g_hat(2,ig) /= ref_grid % g_hat(2,ih) ) CYCLE
            IF ( pw_grid % g_hat(3,ig) /= ref_grid % g_hat(3,ih) ) CYCLE
            pw_grid%gidx(ig) = ih
            EXIT
          END DO
          IF ( pw_grid%gidx(ig) == 0 ) THEN
             WRITE ( *, "(A,2I10)" ) " G-Shell ",is+1,ip+1
             WRITE ( *, "(A,I10,3I6,F20.10)" ) &
                " G-vector", ig, pw_grid % g_hat(1:3,ig), pw_grid % gsq(ig)
             DO ih = 1, SIZE ( ref_grid % gsq  )
               IF ( pw_grid % g_hat(1,ig) /= ref_grid % g_hat(1,ih) ) CYCLE
               IF ( pw_grid % g_hat(2,ig) /= ref_grid % g_hat(2,ih) ) CYCLE
               IF ( pw_grid % g_hat(3,ig) /= ref_grid % g_hat(3,ih) ) CYCLE
               WRITE ( *, "(A,I10,3I6,F20.10)" ) &
                " G-vector", ih, ref_grid % g_hat(1:3,ih),ref_grid % gsq(ih)
             END DO
             CALL stop_program ( "pw_grid_sort", "G vector not found" )
          END IF
          is = ip
        END DO
      ELSE
        ! for the case pw_grid > ref_grid
        is = it
        DO ig = it + 1, ngr
          g_found=.FALSE.
          DO ih = is + 1, ng
            IF ( ABS ( pw_grid % gsq(ih) - ref_grid % gsq(ig) ) > 1.e-12_dp ) CYCLE
            g_found=.TRUE.
            EXIT
          END DO
          IF ( .NOT. g_found ) THEN
             WRITE(*,"(A,I10,F20.10)") "G-vector", ig, ref_grid % gsq(ig)
             CALL stop_program ( "pw_grid_sort", "G vector not found" )
          END IF
          ip = ih - 1
          DO ih = ip + 1, ng
            IF ( ABS ( pw_grid % gsq(ih) - ref_grid % gsq(ig) ) > 1.e-12_dp ) CYCLE
            IF ( pw_grid % g_hat(1,ih) /= ref_grid % g_hat(1,ig) ) CYCLE
            IF ( pw_grid % g_hat(2,ih) /= ref_grid % g_hat(2,ig) ) CYCLE
            IF ( pw_grid % g_hat(3,ih) /= ref_grid % g_hat(3,ig) ) CYCLE
            pw_grid%gidx(ig) = ih
            EXIT
          END DO
          IF ( pw_grid%gidx(ig) == 0 ) THEN
             WRITE ( *, "(A,2I10)" ) " G-Shell ",is+1,ip+1
             WRITE ( *, "(A,I10,3I6,F20.10)" ) &
                " G-vector", ig, ref_grid % g_hat(1:3,ig), ref_grid % gsq(ig)
             DO ih = 1, ng
               IF ( pw_grid % g_hat(1,ih) /= ref_grid % g_hat(1,ig) ) CYCLE
               IF ( pw_grid % g_hat(2,ih) /= ref_grid % g_hat(2,ig) ) CYCLE
               IF ( pw_grid % g_hat(3,ih) /= ref_grid % g_hat(3,ig) ) CYCLE
               WRITE ( *, "(A,I10,3I6,F20.10)" ) &
                " G-vector", ih, pw_grid % g_hat(1:3,ih),pw_grid % gsq(ih)
             END DO
             CALL stop_program ( "pw_grid_sort", "G vector not found" )
          END IF
          is = ip
        END DO
      END IF
      ! test if all g-vectors are associated
      IF ( ANY ( pw_grid%gidx == 0 ) ) THEN
        CALL stop_program ( "pw_grid_sort", "G space sorting not compatible" )
      END IF
    END IF
  END IF

  !check if G=0 is at first position
  IF ( pw_grid % have_g0 ) THEN
    IF ( pw_grid % gsq ( 1 ) /= 0.0_dp ) THEN
      CALL stop_program ( "pw_grid_sort", "G=0 not in first position" )
    END IF
  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_sort

SUBROUTINE sort_shells ( gsq, g_hat, idx )


! Argument
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: gsq
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: g_hat
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: idx

    REAL(KIND=dp), PARAMETER                 :: small = 5.e-16_dp 

    INTEGER                                  :: ig, ng, s1, s2
    REAL(KIND=dp)                            :: s_begin

! Juergs temporary hack to get the grids sorted for large (4000Ry) cutoffs.
! might need to call lapack for machine precision.

  ng = SIZE ( gsq )
  s_begin = -1.0_dp
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


! Argument
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: g_hat
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: idx
    INTEGER, INTENT(IN)                      :: s1, s2

    INTEGER                                  :: i, ii, info, n1, n2, n3, ns
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indl
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: slen

  IF ( s2 <= s1 ) RETURN
  ns = s2 - s1 + 1
  ALLOCATE ( indl ( ns ), STAT = info )
  IF ( info /= 0 ) CALL stop_memory ( "redist", "indl", ns )
  ALLOCATE ( slen ( ns ), STAT = info )
  IF ( info /= 0 ) CALL stop_memory ( "redist", "slen", ns )

  DO i = s1, s2
    ii = idx ( i )
    n1 = g_hat(1,ii)
    n2 = g_hat(2,ii)
    n3 = g_hat(3,ii)
    slen ( i - s1 + 1 ) = 1000.0_dp * REAL(n1,dp) + &
                          REAL(n2,dp) + 0.001_dp * REAL(n3,dp)
  END DO
  CALL sort ( slen, ns, indl )
  DO i = 1, ns
    ii = indl ( i ) + s1 - 1
    indl ( i ) = idx ( ii )
  END DO
  idx ( s1:s2 ) = indl ( 1:ns )

  DEALLOCATE ( indl, STAT = info )
  IF ( info /= 0 ) CALL stop_memory ( "redist", "indl" )
  DEALLOCATE ( slen, STAT = info )
  IF ( info /= 0 ) CALL stop_memory ( "redist", "slen" )

END SUBROUTINE redist

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
!!*** *************************************************************************

SUBROUTINE pw_grid_remap ( pw_grid, yz )

! Argument
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: yz

    INTEGER                                  :: gpt, handle, i, ip, is, j, m, &
                                                n, ny, nz
    INTEGER, DIMENSION(:), POINTER           :: mapm, mapn

!------------------------------------------------------------------------------

  IF ( pw_grid % para % mode == PW_MODE_LOCAL ) RETURN
  IF ( .NOT. pw_grid % para % ray_distribution ) RETURN

  CALL timeset("pw_grid_remap","I","",handle)

  ny = pw_grid % npts ( 2 )
  nz = pw_grid % npts ( 3 )

  yz = 0
  pw_grid % para % yzp = 0
  pw_grid % para % yzq = 0

  mapm => pw_grid % mapm % pos
  mapn => pw_grid % mapn % pos

  DO gpt = 1, SIZE ( pw_grid % gsq  )
    m = mapm ( pw_grid % g_hat ( 2, gpt ) ) + 1
    n = mapn ( pw_grid % g_hat ( 3, gpt ) ) + 1
    yz ( m, n ) = yz ( m, n ) + 1
  END DO
  IF ( pw_grid % grid_span == HALFSPACE ) THEN
     mapm => pw_grid % mapm % neg
     mapn => pw_grid % mapn % neg
     DO gpt = 1, SIZE ( pw_grid % gsq  )
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

  CALL timestop(0.0_dp,handle)

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
!!*** *************************************************************************

SUBROUTINE pw_grid_change ( cell, pw_grid )
  IMPLICIT NONE
  ! Argument
  TYPE(cell_type), POINTER                 :: cell
  TYPE(pw_grid_type), POINTER              :: pw_grid
  ! Local 
  INTEGER                                  :: gpt
  REAL(KIND=dp)                            :: l, m, n
  REAL(KIND=dp), DIMENSION(:, :), POINTER  :: g, h_inv
  
  h_inv => cell % h_inv
  g => pw_grid % g

  pw_grid % vol = ABS ( cell % deth )
  pw_grid % dvol = pw_grid % vol / REAL ( pw_grid % ngpts,KIND=dp)
  pw_grid % dr ( 1 ) = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 1 ),KIND=dp)
  pw_grid % dr ( 2 ) = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 2 ),KIND=dp)
  pw_grid % dr ( 3 ) = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) ) &
       / REAL ( pw_grid % npts ( 3 ),KIND=dp)

  DO gpt = 1, SIZE ( g ,2 )

     l = twopi * REAL ( pw_grid % g_hat ( 1, gpt ),KIND=dp)
     m = twopi * REAL ( pw_grid % g_hat ( 2, gpt ),KIND=dp)
     n = twopi * REAL ( pw_grid % g_hat ( 3, gpt ),KIND=dp)
     
     g ( 1, gpt ) = l * h_inv(1,1) + m * h_inv(2,1) + n * h_inv(3,1)
     g ( 2, gpt ) = l * h_inv(1,2) + m * h_inv(2,2) + n * h_inv(3,2)
     g ( 3, gpt ) = l * h_inv(1,3) + m * h_inv(2,3) + n * h_inv(3,3)
     
     pw_grid % gsq ( gpt ) = g ( 1, gpt ) * g ( 1, gpt ) &
                           + g ( 2, gpt ) * g ( 2, gpt ) &
                           + g ( 3, gpt ) * g ( 3, gpt )

  END DO

END SUBROUTINE pw_grid_change

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
!!*** *************************************************************************

SUBROUTINE pw_find_cutoff ( npts, box, cutoff )


    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(cell_type), INTENT(IN)              :: box
    REAL(KIND=dp), INTENT(OUT)               :: cutoff

    REAL(KIND=dp)                            :: gcut, gdum( 3 ), length

!------------------------------------------------------------------------------
! compute 2*pi*h_inv^t*g  where g = (npts[1],0,0)

  gdum ( : ) = twopi * box % h_inv ( 1, : ) &
       * REAL ( ( npts ( 1 ) - 1 ) / 2,KIND=dp)
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = length

! compute 2*pi*h_inv^t*g  where g = (0,npts[2],0)
  gdum ( : ) = twopi * box % h_inv ( 2, : ) &
       * REAL ( ( npts ( 2 ) - 1 ) / 2,KIND=dp)
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )

! compute 2*pi*h_inv^t*g  where g = (0,0,npts[3])
  gdum ( : ) = twopi * box % h_inv ( 3, : ) &
       * REAL ( ( npts ( 3 ) - 1 ) / 2,KIND=dp)
  length = SQRT ( dotprod_3d ( gdum, gdum ) )
  gcut = MIN ( gcut, length )

  cutoff = gcut - 1.0E-5_dp

END SUBROUTINE pw_find_cutoff

!******************************************************************************

!!****f* pw_grids/pw_grid_retain [1.0] *
!!
!!   NAME
!!     pw_grid_retain
!!
!!   FUNCTION
!!     retains the given pw grid
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - pw_grid: the pw grid to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_grid_retain(pw_grid, error)
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(pw_grid%ref_count>0,cp_failure_level,routineP,error)
     pw_grid%ref_count=pw_grid%ref_count+1
  END IF
END SUBROUTINE pw_grid_retain
!***************************************************************************

!!****f* pw_grids/pw_grid_release [1.0] *
!!
!!   NAME
!!     pw_grid_release
!!
!!   FUNCTION
!!     releases the given pw grid
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - pw_grid: the pw grid to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_grid_release(pw_grid, error)
  TYPE(pw_grid_type), POINTER              :: pw_grid
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_release', &
       routineP = moduleN//':'//routineN

  INTEGER                                  :: stat
  LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(pw_grid)) THEN
     CPPreconditionNoFail(pw_grid%ref_count>0,cp_failure_level,routineP,error)
     pw_grid%ref_count=pw_grid%ref_count-1
     IF (pw_grid%ref_count==0) THEN
        IF ( ASSOCIATED ( pw_grid % gidx ) ) THEN
           DEALLOCATE ( pw_grid % gidx, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % g ) ) THEN
           DEALLOCATE ( pw_grid % g, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % gsq ) ) THEN
           DEALLOCATE ( pw_grid % gsq, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % g_hat ) ) THEN
           DEALLOCATE ( pw_grid % g_hat, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapl % pos ) ) THEN
           DEALLOCATE ( pw_grid % mapl % pos, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapm % pos ) ) THEN
           DEALLOCATE ( pw_grid % mapm % pos, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapn % pos ) ) THEN
           DEALLOCATE ( pw_grid % mapn % pos, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapl % neg ) ) THEN
           DEALLOCATE ( pw_grid % mapl % neg, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapm % neg ) ) THEN
           DEALLOCATE ( pw_grid % mapm % neg, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( pw_grid % mapn % neg ) ) THEN
           DEALLOCATE ( pw_grid % mapn % neg, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( pw_grid % para % mode == PW_MODE_DISTRIBUTED ) THEN
           IF ( ASSOCIATED ( pw_grid % para % yzp ) ) THEN
              DEALLOCATE ( pw_grid % para % yzp, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           IF ( ASSOCIATED ( pw_grid % para % yzq ) ) THEN
              DEALLOCATE ( pw_grid % para % yzq, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           IF ( ASSOCIATED ( pw_grid % para % nyzray ) ) THEN
              DEALLOCATE ( pw_grid % para % nyzray, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           IF ( ASSOCIATED ( pw_grid % para % bo ) ) THEN
              DEALLOCATE ( pw_grid % para % bo, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        END IF
        ! also release groups
        CALL mp_comm_free ( pw_grid % para % group )
        IF (PRODUCT(pw_grid % para % rs_dims) /= 0 ) &
             CALL mp_comm_free ( pw_grid % para % rs_group )
        IF ( ASSOCIATED ( pw_grid % para % pos_of_x ) ) THEN
           DEALLOCATE ( pw_grid % para % pos_of_x, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(pw_grid, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(pw_grid)
END SUBROUTINE pw_grid_release
!***************************************************************************
!!****f* pw_grids/pw_grid_create_copy_no_pbc [1.0] *
!!
!!   NAME
!!     pw_grid_create_copy_no_pbc
!!
!!   FUNCTION
!!     creates a copy of pw_grid_in in which the pbc have been removed
!!     (by adding a point for the upper boundary)
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - pw_grid_in: the pw grid to duplicate
!!     - pw_grid_out: the output pw_grid_type
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi, Teo
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!     04.2005 completly rewritten the duplicate routine, fixed parallel
!!             behaviour, narrowed scope to copy to non pbc and renamed
!!             accordingly [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_grid_create_copy_no_pbc(pw_grid_in, pw_grid_out, cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(pw_grid_type), POINTER              :: pw_grid_in, pw_grid_out
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: cell
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_create_copy_no_pbc', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    INTEGER, DIMENSION(3)                    :: n
    INTEGER, DIMENSION(:),pointer            :: pos_of_x

    failure = .FALSE.
    
    CPPrecondition(pw_grid_in%ngpts_cut>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(pw_grid_out),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL pw_grid_create(pw_grid_out, pw_grid_in%para%group)
       grid_tag = grid_tag + 1
       pw_grid_out %id_nr = grid_tag
       pw_grid_out % ref_count  = 1
       pw_grid_out % reference  = 0

       pw_grid_out%bounds_local=pw_grid_in%bounds_local
       IF (pw_grid_in%bounds_local(2,1)==pw_grid_in%bounds(2,1).AND.&
            pw_grid_in%bounds_local(1,1)<=pw_grid_in%bounds(2,1)) THEN
          pw_grid_out%bounds_local(2,1)=pw_grid_out%bounds_local(2,1)+1
       END IF
       pw_grid_out%bounds_local(2,2)=pw_grid_out%bounds_local(2,2)+1
       pw_grid_out%bounds_local(2,3)=pw_grid_out%bounds_local(2,3)+1

       pw_grid_out%bounds=pw_grid_in%bounds
       pw_grid_out%bounds(2,:)     =  pw_grid_out%bounds(2,:) + 1

       pw_grid_out%npts            =  pw_grid_in%npts  + 1
       pw_grid_out%ngpts           =  PRODUCT ( pw_grid_out% npts )
       pw_grid_out%ngpts_cut=0
       pw_grid_out%npts_local=pw_grid_out%bounds_local(2,:)-pw_grid_out%bounds_local(1,:)+1
       pw_grid_out%ngpts_local=product(pw_grid_out%npts_local)
       pw_grid_out%ngpts_cut_local=0
       pw_grid_out%dr              =  pw_grid_in%dr               
       pw_grid_out%dvol            =  pw_grid_in%dvol             
       pw_grid_out%vol             =  pw_grid_in%vol*REAL(pw_grid_out%ngpts,dp)&
            /REAL(pw_grid_in%ngpts,dp) !FM do not modify?
       pw_grid_out%cutoff          =  pw_grid_in%cutoff
       NULLIFY(pw_grid_out%mapl%pos, pw_grid_out%mapl%neg,&
            pw_grid_out%mapm%pos,pw_grid_out%mapm%neg,&
            pw_grid_out%mapn%pos,pw_grid_out%mapn%neg)
       
       !para
       CALL mp_comm_dup ( pw_grid_in % para % group, pw_grid_out % para % group )
       CALL mp_environ ( pw_grid_out % para % group_size, &
            pw_grid_out % para % my_pos, &
            pw_grid_out % para % group )
       pw_grid_out % para % group_head_id = pw_grid_in % para % group_head_id
       pw_grid_out % para % group_head = &
            ( pw_grid_out % para % group_head_id == pw_grid_out % para % my_pos )
       pw_grid_out % para % mode = pw_grid_in % para % mode
       pw_grid_out % para % ray_distribution = pw_grid_in % para % ray_distribution
       NULLIFY(pw_grid_out % para % yzp, pw_grid_out % para % yzq, &
            pw_grid_out % para % nyzray, pw_grid_out % para % bo)
       ALLOCATE(pos_of_x(pw_grid_out%bounds(1,1):pw_grid_out%bounds(2,1)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       pos_of_x(:pw_grid_out%bounds(2,1)-1)=pw_grid_in%para%pos_of_x
       pos_of_x(pw_grid_out%bounds(2,1))=pos_of_x(pw_grid_out%bounds(2,1)-1)
       pw_grid_out%para%pos_of_x => pos_of_x
       pw_grid_out % para % rs_dims = pw_grid_in % para % rs_dims
       IF (PRODUCT(pw_grid_in % para % rs_dims)/=0) THEN
          CALL mp_comm_dup ( pw_grid_in % para % rs_group, pw_grid_out % para % rs_group )
          CALL mp_cart_rank ( pw_grid_out % para % rs_group, &
               pw_grid_out % para % rs_pos, &
               pw_grid_out % para % rs_mpo )
       ELSE
          pw_grid_out % para % rs_mpo=HUGE(0) !FM copy also from pw_grid_in?
          pw_grid_out % para % rs_pos=HUGE(0)
       END IF
       
       NULLIFY(pw_grid_out%g,pw_grid_out%gsq,pw_grid_out%g_hat)
       CPPrecondition(pw_grid_in%grid_span==FULLSPACE,cp_failure_level,routineP,error,failure)
       pw_grid_out%grid_span=pw_grid_in%grid_span
       pw_grid_out%have_g0=.false.
       pw_grid_out%first_gne0=huge(0)
       pw_grid_out%nglengths=huge(0)
       nullify(pw_grid_out%gidx)
       pw_grid_out%spherical    =  .false.

    END IF
  END SUBROUTINE pw_grid_create_copy_no_pbc

END MODULE pw_grids

!******************************************************************************
