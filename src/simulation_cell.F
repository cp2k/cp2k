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
  USE greens_fn, ONLY : greens_function_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_hinv, cell_type, pbc, get_cell_param

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

!******************************************************************************

END MODULE simulation_cell
