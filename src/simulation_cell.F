!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/simulation_cell [1.0] *
!!
!!   NAME
!!     simulation_cell
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE simulation_cell

  USE kinds, ONLY : dbl

  USE mathconstants, ONLY : pi
  USE termination, ONLY : stop_program
  USE greens_fn, ONLY : greens_function_type, &
            PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, &
            MT2D, MT1D, MT0D, HOCKNEY2D, HOCKNEY1D, HOCKNEY0D

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_hinv, cell_type, pbc, get_cell_param, init_cell

  TYPE cell_type
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: hmat
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: h_inv
     REAL ( dbl ) :: deth
     INTEGER :: perd ( 3 )
     TYPE ( greens_function_type ) :: green
  END TYPE cell_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** simulation_cell/get_hinv [1.0] *
!!
!!   NAME
!!     get_hinv
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_hinv ( box )

  IMPLICIT NONE

! Arguments
  TYPE ( cell_type ), INTENT ( INOUT ) :: box

! Locals
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: hmat, hmati
  REAL ( dbl ) :: odet

!------------------------------------------------------------------------------

  hmat = box % hmat
  box % deth = &
       hmat(1,1) * ( hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2) ) + &
       hmat(1,2) * ( hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3) ) + &
       hmat(1,3) * ( hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1) )
  IF ( box % deth < 1.0E-10_dbl ) &
       CALL stop_program ( 'get_hinv', 'box determinant too small')
  odet = 1.0_dbl / box % deth
  hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
  hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
  hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
  hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
  hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
  hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
  hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
  hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
  hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet

  box % h_inv = hmati

END SUBROUTINE get_hinv

!!*****
!******************************************************************************
!!****** simulation_cell/pbc [1.0] *
!!
!!   NAME
!!     pbc
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

FUNCTION pbc ( rin, box, nl ) RESULT ( rout )

  IMPLICIT NONE

! Return value
  REAL ( dbl ) :: rout ( 3 )

! Arguments
  TYPE ( cell_type ), INTENT ( IN ) :: box
  REAL ( dbl ), INTENT ( IN ) :: rin ( 3 )
  INTEGER, INTENT ( IN ), OPTIONAL :: nl ( 3 )

! Locals
  REAL ( dbl ) :: s ( 3 )

!------------------------------------------------------------------------------

  s = MATMUL ( box % h_inv, rin )
  s = s - box % perd * NINT ( s )
  rout = MATMUL ( box % hmat, s )
  IF ( PRESENT ( nl ) ) THEN
     s = REAL ( nl, dbl )
     rout = rout + MATMUL ( box % hmat, s )
  END IF

END FUNCTION pbc

!!*****
!******************************************************************************
!!****** simulation_cell/get_cell_param [1.0] *
!!
!!   NAME
!!     get_cell_param
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_cell_param ( box, cell_length, cell_angle )

  IMPLICIT NONE

! Arguments
  TYPE ( cell_type ), INTENT ( IN ), TARGET :: box
  REAL ( dbl ), DIMENSION ( 3 ), INTENT ( OUT ) :: cell_length
  REAL ( dbl ), DIMENSION ( 3 ), INTENT ( OUT ), OPTIONAL :: cell_angle

! Locals
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: hmat

!------------------------------------------------------------------------------

  hmat => box % hmat

! This code gets the cell parameters given the h-matrix:
! a
  cell_length(1) = &
       SQRT ( hmat(1,1)*hmat(1,1) + hmat(2,1)*hmat(2,1) + hmat(3,1)*hmat(3,1) )
! b
  cell_length(2) = &
       SQRT ( hmat(1,2)*hmat(1,2) + hmat(2,2)*hmat(2,2) + hmat(3,2)*hmat(3,2) )
! c
  cell_length(3) = &
       SQRT ( hmat(1,3)*hmat(1,3) + hmat(2,3)*hmat(2,3) + hmat(3,3)*hmat(3,3) )

  IF ( PRESENT ( cell_angle ) ) THEN
! gamma
     cell_angle(1) = ACOS ( &
          ( hmat(1,1)*hmat(1,2) + hmat(2,1)*hmat(2,2) + hmat(3,1)*hmat(3,2) ) &
          / ( cell_length(1)*cell_length(2) ) )
! beta
     cell_angle(2) = ACOS ( &
          ( hmat(1,1)*hmat(1,3) + hmat(2,1)*hmat(2,3) + hmat(3,1)*hmat(3,3) ) &
          / ( cell_length(1)*cell_length(3) ) )
! alpha
     cell_angle(3) = ACOS ( &
          ( hmat(1,2)*hmat(1,3) + hmat(2,2)*hmat(2,3) + hmat(3,2)*hmat(3,3) ) &
          / ( cell_length(2)*cell_length(3) ) )
     cell_angle = cell_angle * 180.0_dbl / pi
  END IF

END SUBROUTINE get_cell_param

!!*****
!******************************************************************************
!!****** simulation_cell/init_cell [1.0] *
!!
!!   NAME
!!     init_cell
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (12-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_cell ( box, hmat, periodic, poisson_solver )

  IMPLICIT NONE

! Arguments
  TYPE ( cell_type ), INTENT ( OUT ) :: box
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( IN ) :: hmat
  INTEGER , DIMENSION ( 3 ), INTENT ( IN ) :: periodic
  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: poisson_solver

! Locals
  INTEGER :: dim

!------------------------------------------------------------------------------

  box % hmat = hmat
  box % perd = periodic

  box % green % special_dimension = 0
  box % green % radius = 0._dbl
  box % green % slab_size = 0._dbl
  box % green % alpha = 0._dbl
  box % green % volume = 0._dbl

  dim = COUNT ( box % perd == 1 )
  SELECT CASE ( dim )
  CASE ( 0 )
    box % green % method = ANALYTIC0D
    box % green % radius = 0._dbl
  CASE ( 1 )
    box % green % method = ANALYTIC1D
    box % green % special_dimension = MAXLOC ( box % perd, 1 )
    box % green % radius = 0._dbl
  CASE ( 2 )
    box % green % method = ANALYTIC2D
    box % green % special_dimension = MINLOC ( box % perd, 1 )
    box % green % slab_size = 0._dbl
  CASE ( 3 )
    box % green % method = PERIODIC3D
  CASE DEFAULT
    CALL stop_program ( "init_cell", "illegal perd values" )
  END SELECT

  NULLIFY ( box % green % influence_function )

  IF ( PRESENT ( poisson_solver ) ) THEN
    SELECT CASE ( poisson_solver )
    CASE ( "PERIODIC" )
      box % green % method = PERIODIC3D
      IF ( dim /= 3 ) THEN
        CALL stop_program ( "init_cell", &
             "illegal combination of periodicity and Poisson solver (PERIODIC)" )
      END IF
    CASE ( "ANALYTIC" )

      SELECT CASE ( dim )
      CASE ( 0 )
        box % green % method = ANALYTIC0D
      CASE ( 1 )
        box % green % method = ANALYTIC1D
      CASE ( 2 )
        box % green % method = ANALYTIC2D
      CASE ( 3 )
        box % green % method = PERIODIC3D
      CASE DEFAULT
        CALL stop_program ( "init_cell", "illegal perd values" )
      END SELECT

    CASE ( "HOCKNEY" )

      SELECT CASE ( dim )
      CASE ( 0 )
        box % green % method = HOCKNEY0D
      CASE ( 1 )
        box % green % method = HOCKNEY1D
      CASE ( 2 )
        box % green % method = HOCKNEY2D
      CASE ( 3 )
        CALL stop_program ( "init_cell", &
             "illegal combination of periodicity and Poisson solver (HOCKNEY)" )
      CASE DEFAULT
        CALL stop_program ( "init_cell", "illegal perd values" )
      END SELECT

    CASE ( "MT" )

      SELECT CASE ( dim )
      CASE ( 0 )
        box % green % method = MT0D
      CASE ( 1 )
        box % green % method = MT1D
      CASE ( 2 )
        box % green % method = MT2D
      CASE ( 3 )
        CALL stop_program ( "init_cell", &
             "illegal combination of periodicity and Poisson solver (MT)" )
      CASE DEFAULT
        CALL stop_program ( "init_cell", "illegal perd values" )
      END SELECT

    CASE DEFAULT

      CALL stop_program ( "init_cell", "Unknown Poisson solver" )

    END SELECT

  END IF

END SUBROUTINE init_cell

!!*****
!******************************************************************************

END MODULE simulation_cell

!******************************************************************************
