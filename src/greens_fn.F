!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/greens_fn [1.0] *
!!
!!   NAME
!!     greens_fn
!!
!!   FUNCTION
!!     Routines to calculate the Green's function for the Poisson equation '''
!!
!!   AUTHOR
!!     JGH (19-Dec-2000) 
!! 
!!   MODIFICATION HISTORY
!!     JGH (9-Mar-2001) : include influence_function into greens_function_type
!!                        add cell volume as indicator for updates
!!
!!   SOURCE
!******************************************************************************

MODULE greens_fn
  
  USE bessel_lib, ONLY : bessj0, bessj1, bessk0, bessk1
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, coeff_deallocate, &
                                PLANEWAVES
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : fourpi
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE
  USE pw_types, ONLY : pw_type, pw_allocate, pw_deallocate, &
       REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D, &
       NOSPACE, REALSPACE, RECIPROCALSPACE
  USE termination, ONLY : stop_program
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_green_fn, greens_function_type
  PUBLIC :: PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, &
            MT2D, MT1D, MT0D, HOCKNEY2D, HOCKNEY1D, HOCKNEY0D
  
  TYPE greens_function_type
    INTEGER :: method
    INTEGER :: special_dimension
    REAL ( dbl ) :: radius
    REAL ( dbl ) :: slab_size
    REAL ( dbl ) :: volume
    TYPE ( pw_type ), POINTER :: influence_function
  END TYPE greens_function_type

  INTEGER, PARAMETER :: PERIODIC3D = 1000, ANALYTIC2D = 1001,  &
                        ANALYTIC1D = 1002, ANALYTIC0D = 1003,  &
                        MT2D = 1101, MT1D = 1102, MT0D = 1103, &
                        HOCKNEY2D = 1201, HOCKNEY1D = 1202, HOCKNEY0D = 1203

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** greens_fn/pw_green_fn [1.0] *
!!
!!   NAME
!!     pw_green_fn
!!
!!   FUNCTION
!!     Calculates the Green's function in reciprocal space  '''
!!
!!   AUTHOR
!!     JGH (19-12-2000)
!!
!!   MODIFICATION HISTORY
!!     JGH (9-Mar-2001) : update check through voulme
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_green_fn ( gftype )

  IMPLICIT NONE

! Arguments
  TYPE ( greens_function_type ), INTENT ( INOUT ) :: gftype

! Locals
  TYPE ( pw_type ), POINTER :: gf
  TYPE ( pw_grid_type ), POINTER :: grid
  INTEGER :: method, ig, iz, nz
  REAL ( dbl ) :: rlength, zlength, g3d, gxy, g2, j0g, j1g, k0g, k1g, gg

!------------------------------------------------------------------------------

  method = gftype % method
  gf => gftype % influence_function
  grid => gftype % influence_function % pw_grid

  ! some tests
  SELECT CASE ( method )
  CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
    IF ( gf % in_space /= RECIPROCALSPACE .OR. &
         gf % in_use /= COMPLEXDATA1D ) &
      CALL stop_program ( "greens_fn", "gf wrong space or data type" )
  CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
     CALL stop_program ( "pw_green_fn", " Method not implemented ")
  CASE DEFAULT
     CALL stop_program ( "pw_green_fn", " Method not implemented ")
  END SELECT

  ! check if influence function is consisten with current box volume
  IF ( gftype % volume /= grid % vol ) THEN

    gftype % volume = grid % vol

    SELECT CASE ( method )

    CASE ( PERIODIC3D )

      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % glsq ( ig )
        gf % cc ( ig ) = fourpi / g2
      END DO
      IF ( grid % have_g0 ) gf % cc ( 1 ) = 0._dbl

    CASE ( ANALYTIC2D )

      iz = gftype % special_dimension ! iz is the direction with NO PBC
      zlength = gftype % slab_size    ! zlength is the thickness of the cell
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        nz = grid % g_hat ( iz, ig )
        g2 = grid % glsq ( ig )
        g3d = fourpi / g2
        gg = 0.5_dbl * SQRT ( g2 )
        gf % cc ( ig ) = g3d * ( 1._dbl - (-1._dbl)**nz * EXP ( - gg * zlength ) )
      END DO
      IF ( grid % have_g0 ) gf % cc ( 1 ) = 0._dbl

    CASE ( ANALYTIC1D )

      ! iz is the direction of the PBC ( can be 1,2,3 -> x,y,z )
      iz = gftype % special_dimension
      ! rlength is the radius of the tube
      rlength = gftype % radius
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % glsq ( ig )
        g3d = fourpi / g2
        gxy = SQRT ( g2 - grid % g(iz,ig) * grid % g(iz,ig) )
        j0g = bessj0 ( rlength * gxy )
        j1g = bessj1 ( rlength * gxy )
        k0g = bessk0 ( rlength * grid % g(iz,ig) )
        k1g = bessk1 ( rlength * grid % g(iz,ig) )
        gf % cc ( ig ) = g3d * ( 1._dbl - rlength * &
          ( gxy * j1g * k0g - grid % g(iz,ig) * j0g * k1g ) )
      END DO
      IF ( grid % have_g0 ) gf % cc ( 1 ) = 0._dbl

    CASE ( ANALYTIC0D )

      rlength = gftype % radius   ! rlength is the radius of the sphere
      DO ig = grid % first_gne0, grid % ngpts_cut_local
        g2 = grid % glsq ( ig )
        gg = 0.5_dbl * SQRT ( g2 )
        g3d = fourpi / g2
        gf % cc ( ig ) = g3d * ( 1._dbl - COS ( rlength * gg ) )
      END DO
      IF ( grid % have_g0 ) &
        gf % cc ( 1 ) = 0.5_dbl * fourpi * rlength * rlength

    CASE ( MT2D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( MT1D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( MT0D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY2D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY1D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE ( HOCKNEY0D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE DEFAULT
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    END SELECT

  END IF

END SUBROUTINE pw_green_fn

!!*****
!******************************************************************************

END MODULE greens_fn

!******************************************************************************
