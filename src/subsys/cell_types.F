!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!>      10.2014 Moved many routines from cell_types.F here.
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! *****************************************************************************
MODULE cell_types
  USE cp_units,                        ONLY: cp_unit_to_cp2k,&
                                             cp_units_rad
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: degree,&
                                             sqrt3
  USE mathlib,                         ONLY: angle,&
                                             det_3x3,&
                                             inv_3x3
  USE termination,                     ONLY: stop_program
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_types'

  INTEGER, SAVE, PRIVATE :: last_cell_id = 0

  ! Impose cell symmetry
  INTEGER, PARAMETER, PUBLIC               :: cell_sym_none=0,&
                                              cell_sym_triclinic=1,&
                                              cell_sym_monoclinic=2,&
                                              cell_sym_orthorhombic=3,&
                                              cell_sym_tetragonal_ab=4,&
                                              cell_sym_tetragonal_ac=5,&
                                              cell_sym_tetragonal_bc=6,&
                                              cell_sym_rhombohedral=7,&
                                              cell_sym_hexagonal=8,&
                                              cell_sym_cubic=9

  INTEGER, PARAMETER, PUBLIC               :: use_perd_x=0,&
                                              use_perd_y=1,&
                                              use_perd_z=2,&
                                              use_perd_xy=3,&
                                              use_perd_xz=4,&
                                              use_perd_yz=5,&
                                              use_perd_xyz=6,&
                                              use_perd_none=7

! *****************************************************************************
!> \brief   Type defining parameters related to the simulation cell
!> \version 1.0
! *****************************************************************************
  TYPE cell_type
     INTEGER                           :: ref_count,id_nr
     INTEGER                           :: symmetry_id
     LOGICAL                           :: orthorhombic ! actually means a diagonal hmat
     REAL(KIND = dp)                   :: deth
     INTEGER, DIMENSION(3)             :: perd
     REAL(KIND = dp), DIMENSION(3,3)   :: hmat,h_inv
  END TYPE cell_type

  TYPE cell_p_type
     TYPE(cell_type),POINTER :: cell
  END TYPE cell_p_type

  ! Public data types
  PUBLIC :: cell_type, cell_p_type

  ! Public subroutines
  PUBLIC :: get_cell, get_cell_param, init_cell,&
            cell_create, cell_retain, cell_release,&
            cell_clone, cell_copy, parse_cell_line, set_cell_param,&
            pbc_cp2k_plumed_getset_cell, pbc_cp2k_plumed  !for plumed

  ! Public functions
  PUBLIC :: plane_distance, pbc, real_to_scaled, scaled_to_real

  INTERFACE pbc
     MODULE PROCEDURE pbc1,pbc2,pbc3,pbc4
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cell_in ...
!> \param cell_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cell_clone (cell_in, cell_out, error)

    TYPE(cell_type), POINTER                 :: cell_in, cell_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_clone', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(cell_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell_out),cp_failure_level,routineP,error,failure)

    cell_out%deth         = cell_in%deth
    cell_out%perd         = cell_in%perd
    cell_out%hmat         = cell_in%hmat
    cell_out%h_inv        = cell_in%h_inv
    cell_out%orthorhombic = cell_in%orthorhombic
    cell_out%symmetry_id  = cell_in%symmetry_id
    cell_out%ref_count    = 1
    last_cell_id          = last_cell_id + 1
    cell_out%id_nr        = last_cell_id

  END SUBROUTINE cell_clone

! *****************************************************************************
!> \brief ...
!> \param cell_in ...
!> \param cell_out ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cell_copy (cell_in, cell_out, error)

    TYPE(cell_type), POINTER                 :: cell_in, cell_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_copy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(cell_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell_out),cp_failure_level,routineP,error,failure)

    cell_out%deth         = cell_in%deth
    cell_out%perd         = cell_in%perd
    cell_out%hmat         = cell_in%hmat
    cell_out%h_inv        = cell_in%h_inv
    cell_out%orthorhombic = cell_in%orthorhombic
    cell_out%symmetry_id  = cell_in%symmetry_id

  END SUBROUTINE cell_copy

! *****************************************************************************
!> \brief   Read cell info from a line (parsed from a file)
!> \param input_line ...
!> \param cell_itimes ...
!> \param cell_time ...
!> \param h ...
!> \param vol ...
!> \param error ...
!> \date    19.02.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \version 1.0
! *****************************************************************************
  SUBROUTINE parse_cell_line(input_line, cell_itimes, cell_time, h, vol, error)
    CHARACTER(LEN=*), INTENT(IN)             :: input_line
    INTEGER, INTENT(OUT)                     :: cell_itimes
    REAL(KIND=dp), INTENT(OUT)               :: cell_time
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: h
    REAL(KIND=dp), INTENT(OUT)               :: vol
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'parse_cell_line', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure

    failure = .FALSE.
    READ (input_line,*) cell_itimes, cell_time,&
         h(1,1),h(2,1),h(3,1), h(1,2),h(2,2),h(3,2), h(1,3),h(2,3),h(3,3),vol
    DO i = 1, 3
       DO j = 1, 3
          h(j,i) = cp_unit_to_cp2k(h(j,i), "angstrom", error=error)
       END DO
    END DO

  END SUBROUTINE parse_cell_line

! *****************************************************************************
!> \brief   Get informations about a simulation cell.
!> \param cell ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param deth ...
!> \param orthorhombic ...
!> \param abc ...
!> \param periodic ...
!> \param h ...
!> \param h_inv ...
!> \param id_nr ...
!> \param symmetry_id ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_cell(cell,alpha,beta,gamma,deth,orthorhombic,abc,periodic,&
                      h,h_inv,id_nr,symmetry_id)

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha, beta, gamma, deth
    LOGICAL, INTENT(OUT), OPTIONAL           :: orthorhombic
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: abc
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: periodic
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT), OPTIONAL                  :: h, h_inv
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, symmetry_id

    IF (PRESENT(deth)) deth = cell%deth ! the volume
    IF (PRESENT(orthorhombic)) orthorhombic = cell%orthorhombic
    IF (PRESENT(periodic)) periodic(:) = cell%perd(:)
    IF (PRESENT(h)) h(:,:) = cell%hmat(:,:)
    IF (PRESENT(h_inv)) h_inv(:,:) = cell%h_inv(:,:)

    ! Calculate the lengths of the cell vectors a, b, and c
    IF (PRESENT(abc)) THEN
       abc(1) = SQRT(cell%hmat(1,1)*cell%hmat(1,1) +&
                     cell%hmat(2,1)*cell%hmat(2,1) +&
                     cell%hmat(3,1)*cell%hmat(3,1))
       abc(2) = SQRT(cell%hmat(1,2)*cell%hmat(1,2) +&
                     cell%hmat(2,2)*cell%hmat(2,2) +&
                     cell%hmat(3,2)*cell%hmat(3,2))
       abc(3) = SQRT(cell%hmat(1,3)*cell%hmat(1,3) +&
                     cell%hmat(2,3)*cell%hmat(2,3) +&
                     cell%hmat(3,3)*cell%hmat(3,3))
    END IF

    ! Angles between the cell vectors a, b, and c
    ! alpha = <(b,c)
    IF (PRESENT(alpha)) alpha = angle(cell%hmat(:,2),cell%hmat(:,3))*degree
    ! beta = <(a,c)
    IF (PRESENT(beta))   beta = angle(cell%hmat(:,1),cell%hmat(:,3))*degree
    ! gamma = <(a,b)
    IF (PRESENT(gamma)) gamma = angle(cell%hmat(:,1),cell%hmat(:,2))*degree
    IF (PRESENT(id_nr)) id_nr = cell%id_nr
    IF (PRESENT(symmetry_id)) symmetry_id = cell%symmetry_id
  END SUBROUTINE get_cell

! *****************************************************************************
!> \brief   Access internal type variables
!> \param cell ...
!> \param cell_length ...
!> \param cell_angle ...
!> \param units_angle ...
!> \param periodic ...
!> \param error ...
!> \date    04.04.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_cell_param(cell,cell_length,cell_angle,units_angle,periodic,&
                            error)

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: cell_length
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: cell_angle
    INTEGER, INTENT(IN), OPTIONAL            :: units_angle
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: periodic
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_cell_param', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: alpha, beta, gamma

    CALL get_cell(cell=cell,abc=cell_length)

    IF (PRESENT(cell_angle)) THEN
      CALL get_cell(cell=cell,alpha=alpha,beta=beta,gamma=gamma)
      cell_angle(:) = (/alpha,beta,gamma/)
      IF (PRESENT(units_angle)) THEN
         IF (units_angle == cp_units_rad) cell_angle = cell_angle/degree
      END IF
    END IF

    IF (PRESENT(periodic)) CALL get_cell(cell=cell,periodic=periodic)

  END SUBROUTINE get_cell_param

! *****************************************************************************
!> \brief   Sets the cell using the internal parameters (a,b,c) (alpha,beta,gamma)
!>          using the convention: a parallel to the x axis, b in the x-y plane and
!>          and c univoquely determined; gamma is the angle between a and b; beta
!>          is the angle between c and a and alpha is the angle between c and b
!> \param cell ...
!> \param cell_length ...
!> \param cell_angle ...
!> \param periodic ...
!> \param do_init_cell ...
!> \param error ...
!> \date    03.2008
!> \author  Teodoro Laino
! *****************************************************************************
  SUBROUTINE set_cell_param(cell,cell_length,cell_angle,periodic,do_init_cell,&
                            error)

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: cell_length, cell_angle
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    LOGICAL, INTENT(IN)                      :: do_init_cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_cell_param', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cos_alpha, cos_beta, &
                                                cos_gamma, eps, sin_gamma

    failure = .FALSE.
    CPPrecondition(ALL(cell_angle/=0.0_dp),cp_failure_level,routineP,error,failure)
    eps = EPSILON(0.0_dp)
    cos_gamma = COS(cell_angle(3)); IF (ABS(cos_gamma)<eps) cos_gamma = 0.0_dp
    IF (ABS(ABS(cos_gamma)-1.0_dp)<eps) cos_gamma = SIGN(1.0_dp,cos_gamma)
    sin_gamma = SIN(cell_angle(3)); IF (ABS(sin_gamma)<eps) sin_gamma = 0.0_dp
    IF (ABS(ABS(sin_gamma)-1.0_dp)<eps) sin_gamma = SIGN(1.0_dp,sin_gamma)
    cos_beta  = COS(cell_angle(2)); IF (ABS(cos_beta )<eps) cos_beta  = 0.0_dp
    IF (ABS(ABS(cos_beta)-1.0_dp)<eps) cos_beta = SIGN(1.0_dp,cos_beta)
    cos_alpha = COS(cell_angle(1)); IF (ABS(cos_alpha)<eps) cos_alpha = 0.0_dp
    IF (ABS(ABS(cos_alpha)-1.0_dp)<eps) cos_alpha = SIGN(1.0_dp,cos_alpha)

    cell%hmat(:,1) = (/   1.0_dp,                                   0.0_dp,  0.0_dp/)
    cell%hmat(:,2) = (/cos_gamma,                                sin_gamma,  0.0_dp/)
    cell%hmat(:,3) = (/ cos_beta, (cos_alpha-cos_gamma*cos_beta)/sin_gamma,  0.0_dp/)
    cell%hmat(3,3) = SQRT(1.0_dp-cell%hmat(1,3)**2-cell%hmat(2,3)**2)

    cell%hmat(:,1) = cell%hmat(:,1)*cell_length(1)
    cell%hmat(:,2) = cell%hmat(:,2)*cell_length(2)
    cell%hmat(:,3) = cell%hmat(:,3)*cell_length(3)

    IF (do_init_cell) THEN
       IF (PRESENT(periodic)) THEN
          CALL init_cell(cell=cell,periodic=periodic)
       ELSE
          CALL init_cell(cell=cell)
       END IF
    END IF

  END SUBROUTINE set_cell_param

! *****************************************************************************
!> \brief   Initialise/readjust a simulation cell after hmat has been changed
!> \param cell ...
!> \param hmat ...
!> \param periodic ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_cell(cell,hmat,periodic)

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic

    CHARACTER(len=*), PARAMETER :: routineN = 'init_cell', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: eps_hmat = 1.0E-14_dp

    INTEGER                                  :: dim
    REAL(KIND=dp)                            :: a, acosa, acosah, alpha, &
                                                asina, asinah, beta, norm, &
                                                norm_c
    REAL(KIND=dp), DIMENSION(3)              :: abc

    IF (PRESENT(hmat)) cell%hmat(:,:) = hmat(:,:)
    IF (PRESENT(periodic)) cell%perd(:) = periodic(:)

    cell%deth = ABS(det_3x3(cell%hmat))

    IF (cell%deth < 1.0E-10_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "An invalid set of cell vectors was specified. "//&
                         "The determinant det(h) is too small")
    END IF
    cell%h_inv = inv_3x3(cell%hmat)

    SELECT CASE (cell%symmetry_id)
    CASE (cell_sym_cubic,&
          cell_sym_tetragonal_ab,&
          cell_sym_tetragonal_ac,&
          cell_sym_tetragonal_bc,&
          cell_sym_orthorhombic)
       CALL get_cell(cell=cell,abc=abc)
       abc(2) = plane_distance(0,1,0,cell=cell)
       abc(3) = plane_distance(0,0,1,cell=cell)
       SELECT CASE (cell%symmetry_id)
       CASE (cell_sym_cubic)
          abc(1:3) = SUM(abc(1:3))/3.0_dp
       CASE (cell_sym_tetragonal_ab,&
             cell_sym_tetragonal_ac,&
             cell_sym_tetragonal_bc)
          SELECT CASE (cell%symmetry_id)
          CASE (cell_sym_tetragonal_ab)
             a = 0.5_dp*(abc(1) + abc(2))
             abc(1) = a
             abc(2) = a
          CASE (cell_sym_tetragonal_ac)
             a = 0.5_dp*(abc(1) + abc(3))
             abc(1) = a
             abc(3) = a
          CASE (cell_sym_tetragonal_bc)
             a = 0.5_dp*(abc(2) + abc(3))
             abc(2) = a
             abc(3) = a
          END SELECT
       END SELECT
       cell%hmat(1,1) = abc(1); cell%hmat(1,2) = 0.0_dp; cell%hmat(1,3) = 0.0_dp
       cell%hmat(2,1) = 0.0_dp; cell%hmat(2,2) = abc(2); cell%hmat(2,3) = 0.0_dp
       cell%hmat(3,1) = 0.0_dp; cell%hmat(3,2) = 0.0_dp; cell%hmat(3,3) = abc(3)
    CASE (cell_sym_hexagonal)
       CALL get_cell(cell=cell,abc=abc)
       a = 0.5_dp*(abc(1) + abc(2))
       acosa = 0.5_dp*a
       asina = sqrt3*acosa
       cell%hmat(1,1) = a     ; cell%hmat(1,2) = acosa ; cell%hmat(1,3) = 0.0_dp
       cell%hmat(2,1) = 0.0_dp; cell%hmat(2,2) = asina ; cell%hmat(2,3) = 0.0_dp
       cell%hmat(3,1) = 0.0_dp; cell%hmat(3,2) = 0.0_dp; cell%hmat(3,3) = abc(3)
    CASE (cell_sym_rhombohedral)
       CALL get_cell(cell=cell,abc=abc)
       a = SUM(abc(1:3))/3.0_dp
       alpha = (angle(cell%hmat(:,3),cell%hmat(:,2)) +&
                angle(cell%hmat(:,1),cell%hmat(:,3)) +&
                angle(cell%hmat(:,1),cell%hmat(:,2)))/3.0_dp
       acosa = a*COS(alpha)
       asina = a*SIN(alpha)
       acosah = a*COS(0.5_dp*alpha)
       asinah = a*SIN(0.5_dp*alpha)
       norm = acosa/acosah
       norm_c = SQRT(1.0_dp - norm*norm)
       cell%hmat(1,1) = a     ; cell%hmat(1,2) = acosa ; cell%hmat(1,3) = acosah*norm
       cell%hmat(2,1) = 0.0_dp; cell%hmat(2,2) = asina ; cell%hmat(2,3) = asinah*norm
       cell%hmat(3,1) = 0.0_dp; cell%hmat(3,2) = 0.0_dp; cell%hmat(3,3) = a*norm_c
    CASE (cell_sym_monoclinic)
       CALL get_cell(cell=cell,abc=abc)
       beta = angle(cell%hmat(:,1),cell%hmat(:,3))
       cell%hmat(1,1) = abc(1); cell%hmat(1,2) = 0.0_dp; cell%hmat(1,3) = abc(3)*COS(beta)
       cell%hmat(2,1) = 0.0_dp; cell%hmat(2,2) = abc(2); cell%hmat(2,3) = 0.0_dp
       cell%hmat(3,1) = 0.0_dp; cell%hmat(3,2) = 0.0_dp; cell%hmat(3,3) = abc(3)*SIN(beta)
    CASE (cell_sym_triclinic)
       ! Nothing to do
    END SELECT

    ! Do we have an (almost) orthorhombic cell?
    IF ((ABS(cell%hmat(1,2)) < eps_hmat).AND.(ABS(cell%hmat(1,3)) < eps_hmat).AND.&
        (ABS(cell%hmat(2,1)) < eps_hmat).AND.(ABS(cell%hmat(2,3)) < eps_hmat).AND.&
        (ABS(cell%hmat(3,1)) < eps_hmat).AND.(ABS(cell%hmat(3,2)) < eps_hmat)) THEN
       cell%orthorhombic = .TRUE.
    ELSE
       cell%orthorhombic = .FALSE.
    END IF

    ! Retain an exact orthorhombic cell
    ! (off-diagonal elements must remain zero identically to keep QS fast)
    IF (cell%orthorhombic) THEN
       cell%hmat(1,2) = 0.0_dp
       cell%hmat(1,3) = 0.0_dp
       cell%hmat(2,1) = 0.0_dp
       cell%hmat(2,3) = 0.0_dp
       cell%hmat(3,1) = 0.0_dp
       cell%hmat(3,2) = 0.0_dp
    END IF

    dim = COUNT(cell%perd == 1)
    IF ((dim == 1).AND.(.NOT.cell%orthorhombic)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Non-orthorhombic and not periodic")
    END IF

  END SUBROUTINE init_cell

! *****************************************************************************
!> \brief   Calculate the distance between two lattice planes as defined by
!>          a triple of Miller indices (hkl).
!> \param h ...
!> \param k ...
!> \param l ...
!> \param cell ...
!> \retval distance ...
!> \date    18.11.2004
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION plane_distance(h,k,l,cell) RESULT(distance)

    INTEGER, INTENT(IN)                      :: h, k, l
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp)                            :: distance

    REAL(KIND=dp)                            :: a, alpha, b, beta, c, cosa, &
                                                cosb, cosg, d, gamma, x, y, z
    REAL(KIND=dp), DIMENSION(3)              :: abc

    x = REAL(h,KIND=dp)
    y = REAL(k,KIND=dp)
    z = REAL(l,KIND=dp)

    CALL get_cell(cell=cell,abc=abc)

    a = abc(1)
    b = abc(2)
    c = abc(3)

    IF (cell%orthorhombic) THEN

      d = (x/a)**2 + (y/b)**2 + (z/c)**2

    ELSE

      CALL get_cell(cell=cell,&
                    alpha=alpha,&
                    beta=beta,&
                    gamma=gamma)

      alpha = alpha/degree
      beta = beta/degree
      gamma = gamma/degree

      cosa = COS(alpha)
      cosb = COS(beta)
      cosg = COS(gamma)

      d = ((x*b*c*SIN(alpha))**2 +&
           (y*c*a*SIN(beta))**2 +&
           (z*a*b*SIN(gamma))**2 +&
           2.0_dp*a*b*c*(x*y*c*(cosa*cosb - cosg) +&
                         z*x*b*(cosg*cosa - cosb) +&
                         y*z*a*(cosb*cosg - cosa)))/&
          ((a*b*c)**2*(1.0_dp - cosa**2 - cosb**2 - cosg**2 +&
           2.0_dp*cosa*cosb*cosg))

    END IF

    distance = 1.0_dp/SQRT(d)

  END FUNCTION plane_distance

! *****************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r.
!> \param r ...
!> \param cell ...
!> \retval r_pbc ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION pbc1(r,cell) RESULT(r_pbc)

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc

    REAL(KIND=dp), DIMENSION(3)              :: s

    IF (cell%orthorhombic) THEN
      r_pbc(1) = r(1) - cell%hmat(1,1)*cell%perd(1)*ANINT(cell%h_inv(1,1)*r(1))
      r_pbc(2) = r(2) - cell%hmat(2,2)*cell%perd(2)*ANINT(cell%h_inv(2,2)*r(2))
      r_pbc(3) = r(3) - cell%hmat(3,3)*cell%perd(3)*ANINT(cell%h_inv(3,3)*r(3))
    ELSE
      s(1) = cell%h_inv(1,1)*r(1) + cell%h_inv(1,2)*r(2) + cell%h_inv(1,3)*r(3)
      s(2) = cell%h_inv(2,1)*r(1) + cell%h_inv(2,2)*r(2) + cell%h_inv(2,3)*r(3)
      s(3) = cell%h_inv(3,1)*r(1) + cell%h_inv(3,2)*r(2) + cell%h_inv(3,3)*r(3)
      s(1) = s(1) - cell%perd(1)*ANINT(s(1))
      s(2) = s(2) - cell%perd(2)*ANINT(s(2))
      s(3) = s(3) - cell%perd(3)*ANINT(s(3))
      r_pbc(1) = cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
      r_pbc(2) = cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
      r_pbc(3) = cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
    END IF

  END FUNCTION pbc1

! *****************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r subtracting nl from the periodic images
!> \param r ...
!> \param cell ...
!> \param nl ...
!> \retval r_pbc ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION pbc2(r,cell,nl) RESULT(r_pbc)

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(3), INTENT(IN)        :: nl
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc

    REAL(KIND=dp), DIMENSION(3)              :: s

    IF (cell%orthorhombic) THEN
      r_pbc(1) = r(1) - cell%hmat(1,1)*cell%perd(1)*&
                        REAL(NINT(cell%h_inv(1,1)*r(1)) - nl(1),dp)
      r_pbc(2) = r(2) - cell%hmat(2,2)*cell%perd(2)*&
                        REAL(NINT(cell%h_inv(2,2)*r(2)) - nl(2),dp)
      r_pbc(3) = r(3) - cell%hmat(3,3)*cell%perd(3)*&
                        REAL(NINT(cell%h_inv(3,3)*r(3)) - nl(3),dp)
    ELSE
      s(1) = cell%h_inv(1,1)*r(1) + cell%h_inv(1,2)*r(2) + cell%h_inv(1,3)*r(3)
      s(2) = cell%h_inv(2,1)*r(1) + cell%h_inv(2,2)*r(2) + cell%h_inv(2,3)*r(3)
      s(3) = cell%h_inv(3,1)*r(1) + cell%h_inv(3,2)*r(2) + cell%h_inv(3,3)*r(3)
      s(1) = s(1) - cell%perd(1)*REAL(NINT(s(1)) - nl(1),dp)
      s(2) = s(2) - cell%perd(2)*REAL(NINT(s(2)) - nl(2),dp)
      s(3) = s(3) - cell%perd(3)*REAL(NINT(s(3)) - nl(3),dp)
      r_pbc(1) = cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
      r_pbc(2) = cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
      r_pbc(3) = cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
    END IF

  END FUNCTION pbc2

! *****************************************************************************
!> \brief   Apply the periodic boundary conditions defined by the simulation
!>          cell cell to the vector pointing from atom a to atom b.
!> \param ra ...
!> \param rb ...
!> \param cell ...
!> \retval rab_pbc ...
!> \date    11.03.2004
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  FUNCTION pbc3(ra,rb,cell) RESULT(rab_pbc)

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: rab_pbc

    INTEGER                                  :: icell, jcell, kcell
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(KIND=dp)                            :: rab2, rab2_pbc
    REAL(KIND=dp), DIMENSION(3)              :: r, ra_pbc, rab, rb_image, &
                                                rb_pbc, s2r

    CALL get_cell(cell=cell,periodic=periodic)

    ra_pbc(:) = pbc(ra(:),cell)
    rb_pbc(:) = pbc(rb(:),cell)

    rab2_pbc = HUGE(1.0_dp)

    DO icell=-periodic(1),periodic(1)
      DO jcell=-periodic(2),periodic(2)
        DO kcell=-periodic(3),periodic(3)
          r = REAL((/icell,jcell,kcell/),dp)
          CALL scaled_to_real(s2r,r,cell)
          rb_image(:) = rb_pbc(:) + s2r
          rab(:) = rb_image(:) - ra_pbc(:)
          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
          IF (rab2 < rab2_pbc) THEN
            rab2_pbc = rab2
            rab_pbc(:) = rab(:)
          END IF
        END DO
      END DO
    END DO

  END FUNCTION pbc3
  
  !if positive_range == true, r(i) (or s(i)) in range [0, hmat(i,i)],
  !else, r(i) (s(i)) in range [-hmat(i,i)/2, hmat(i,i)/2]
! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param cell ...
!> \param positive_range ...
!> \retval r_pbc ...
! *****************************************************************************
  FUNCTION pbc4(r,cell,positive_range) RESULT(r_pbc)

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL                                  :: positive_range
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc

    REAL(KIND=dp), DIMENSION(3)              :: s

    IF (positive_range) THEN
      IF (cell%orthorhombic) THEN
        r_pbc(1) = r(1) - cell%hmat(1,1)*cell%perd(1)*FLOOR(cell%h_inv(1,1)*r(1))
        r_pbc(2) = r(2) - cell%hmat(2,2)*cell%perd(2)*FLOOR(cell%h_inv(2,2)*r(2))
        r_pbc(3) = r(3) - cell%hmat(3,3)*cell%perd(3)*FLOOR(cell%h_inv(3,3)*r(3))
      ELSE
        s(1) = cell%h_inv(1,1)*r(1) + cell%h_inv(1,2)*r(2) + cell%h_inv(1,3)*r(3)
        s(2) = cell%h_inv(2,1)*r(1) + cell%h_inv(2,2)*r(2) + cell%h_inv(2,3)*r(3)
        s(3) = cell%h_inv(3,1)*r(1) + cell%h_inv(3,2)*r(2) + cell%h_inv(3,3)*r(3)
        s(1) = s(1) - cell%perd(1)*FLOOR(s(1))
        s(2) = s(2) - cell%perd(2)*FLOOR(s(2))
        s(3) = s(3) - cell%perd(3)*FLOOR(s(3))
        r_pbc(1) = cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
        r_pbc(2) = cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
        r_pbc(3) = cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
      END IF
    ELSE
      r_pbc = pbc1(r,cell)
    END IF

  END FUNCTION pbc4

! *****************************************************************************
!> \brief   Transform real to scaled cell coordinates.
!>          s=h_inv*r
!> \param s ...
!> \param r ...
!> \param cell ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE real_to_scaled(s,r,cell)

    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: s
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell

    IF (cell%orthorhombic) THEN
      s(1) = cell%h_inv(1,1)*r(1)
      s(2) = cell%h_inv(2,2)*r(2)
      s(3) = cell%h_inv(3,3)*r(3)
    ELSE
      s(1) = cell%h_inv(1,1)*r(1) + cell%h_inv(1,2)*r(2) + cell%h_inv(1,3)*r(3)
      s(2) = cell%h_inv(2,1)*r(1) + cell%h_inv(2,2)*r(2) + cell%h_inv(2,3)*r(3)
      s(3) = cell%h_inv(3,1)*r(1) + cell%h_inv(3,2)*r(2) + cell%h_inv(3,3)*r(3)
    END IF

  END SUBROUTINE real_to_scaled

! *****************************************************************************
!> \brief   Transform scaled cell coordinates real coordinates.
!>          r=h*s
!> \param r ...
!> \param s ...
!> \param cell ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE scaled_to_real(r,s,cell)

    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: r
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: s
    TYPE(cell_type), POINTER                 :: cell

    IF (cell%orthorhombic) THEN
      r(1) = cell%hmat(1,1)*s(1)
      r(2) = cell%hmat(2,2)*s(2)
      r(3) = cell%hmat(3,3)*s(3)
    ELSE
      r(1) = cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
      r(2) = cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
      r(3) = cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
    END IF

  END SUBROUTINE scaled_to_real

! *****************************************************************************
!> \brief allocates and initializes a cell
!> \param cell the cell to initialize
!> \param hmat the h matrix that defines the cell
!> \param periodic periodicity of the cell
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cell_create(cell,hmat,periodic,error)

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    ALLOCATE (cell,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_cell_id = last_cell_id + 1
    cell%id_nr = last_cell_id
    cell%ref_count = 1
    IF (PRESENT(periodic)) THEN
       cell%perd = periodic
    ELSE
       cell%perd = 1
    END IF
    cell%orthorhombic = .FALSE.
    cell%symmetry_id = cell_sym_none
    IF (.NOT.PRESENT(hmat)) RETURN
    IF (.NOT.failure) THEN
       CALL init_cell(cell,hmat,periodic)
    END IF

  END SUBROUTINE cell_create

! *****************************************************************************
!> \brief retains the given cell (see doc/ReferenceCounting.html)
!> \param cell the cell to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cell_retain(cell,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(cell%ref_count>0,cp_failure_level,routineP,error)
       cell%ref_count=cell%ref_count+1
    END IF

  END SUBROUTINE cell_retain

! *****************************************************************************
!> \brief releases the given cell (see doc/ReferenceCounting.html)
!> \param cell the cell to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cell_release(cell,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(cell)) THEN
       CPPreconditionNoFail(cell%ref_count>0,cp_failure_level,routineP,error)
       cell%ref_count=cell%ref_count-1
       IF (cell%ref_count==0) THEN
          DEALLOCATE(cell,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       NULLIFY(cell)
    END IF

  END SUBROUTINE cell_release

! *****************************************************************************
!> \brief   Interface with plumed, get a coordinate 3-tuple, return the minimum 
!>          image.
!> \param x ...
!> \param y ...
!> \param z ...
!> \date    12.12.2012
!> \author  RK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE pbc_cp2k_plumed(x,y,z) 
    REAL(KIND=dp), INTENT(INOUT)             :: x, y, z

    REAL(KIND=dp), DIMENSION(3)              :: r, r2
    TYPE(cell_type), POINTER                 :: cell

    CALL pbc_cp2k_plumed_getset_cell(cell, set=.FALSE.) !retrieve the cell pointer

    r(1)=x
    r(2)=y
    r(3)=z

    r2 = pbc(r,cell)

    x=r2(1)
    y=r2(2)
    z=r2(3)

  END SUBROUTINE pbc_cp2k_plumed

! *****************************************************************************
!> \brief   For the interface with plumed, pass a cell pointer and retrieve it 
!>          later. It's a hack, but avoids passing the cell back and forth
!>          across the Fortran/C++ interface
!> \param cell ...
!> \param set ...
!> \date    28.02.2013
!> \author  RK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE pbc_cp2k_plumed_getset_cell(cell, set) 
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL                                  :: set

    TYPE(cell_type), POINTER, SAVE           :: stored_cell

    IF (set .EQV. .TRUE.) THEN
        stored_cell => cell
    ELSE
        cell => stored_cell
    END IF

  END SUBROUTINE pbc_cp2k_plumed_getset_cell

END MODULE cell_types
