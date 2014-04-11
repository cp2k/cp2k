!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! *****************************************************************************
MODULE cell_types
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k,&
                                             cp_units_rad
  USE f77_blas
  USE input_constants,                 ONLY: &
       cell_sym_cubic, cell_sym_hexagonal, cell_sym_monoclinic, &
       cell_sym_none, cell_sym_orthorhombic, cell_sym_rhombohedral, &
       cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
       cell_sym_tetragonal_bc, cell_sym_triclinic, do_cell_cp2k, do_cell_xsc, &
       use_perd_none, use_perd_x, use_perd_xy, use_perd_xyz, use_perd_xz, &
       use_perd_y, use_perd_yz, use_perd_z
  USE input_cp2k_subsys,               ONLY: create_cell_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: &
       section_get_keyword, section_release, section_type, section_vals_get, &
       section_vals_get_subs_vals, section_vals_type, section_vals_val_get, &
       section_vals_val_set, section_vals_val_unset
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: degree,&
                                             sqrt3
  USE mathlib,                         ONLY: angle,&
                                             det_3x3,&
                                             inv_3x3
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_types'

  INTEGER, SAVE, PRIVATE :: last_cell_id = 0

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

  ! Public subroutines
  PUBLIC :: get_cell, get_cell_param, init_cell, read_cell,&
            write_cell, cell_create, cell_retain, cell_release,&
            cell_clone, cell_copy, compare_cells, parse_cell_line, set_cell_param, &
            pbc_cp2k_plumed, pbc_cp2k_plumed_getset_cell !for plumed

  ! Public functions
  PUBLIC :: plane_distance, pbc, real_to_scaled, scaled_to_real

  ! Public data types
  PUBLIC :: cell_type, cell_p_type

  INTERFACE pbc
     MODULE PROCEDURE pbc1,pbc2,pbc3
  END INTERFACE

CONTAINS

! *****************************************************************************
  FUNCTION compare_cells (cell_1, cell_2, error) RESULT(compare)

    TYPE(cell_type), POINTER                 :: cell_1, cell_2
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL                                  :: compare

    CHARACTER(len=*), PARAMETER :: routineN = 'compare_cells', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(cell_1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell_2),cp_failure_level,routineP,error,failure)

    compare = ((cell_2%deth == cell_1%deth).AND.&
               (ALL(cell_2%perd == cell_1%perd)).AND.&
               (ALL(cell_2%hmat == cell_1%hmat)).AND.&
               (ALL(cell_2%h_inv == cell_1%h_inv)).AND.&
               (cell_2%orthorhombic.EQV.cell_1%orthorhombic).AND.&
               (cell_2%symmetry_id == cell_1%symmetry_id))

  END FUNCTION compare_cells

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
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    19.02.2008
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
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \author  Matthias Krack
!> \date    04.04.2002
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
!> \author  Teodoro Laino
!> \date    03.2008
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
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \author  Matthias Krack
!> \date    18.11.2004
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
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \author  Matthias Krack
!> \date    11.03.2004
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

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2005 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE read_cell( cell, cell_ref, use_ref_cell, cell_section,&

    check_for_ref, para_env, error)
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(OUT), OPTIONAL           :: use_ref_cell
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: cell_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_for_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_per
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL :: cell_read_a, cell_read_abc, cell_read_b, cell_read_c, &
      cell_read_file, check, failure, my_check
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_angles, cell_par
    TYPE(section_vals_type), POINTER         :: cell_ref_section

    failure  = .FALSE.
    my_check = .TRUE.
    NULLIFY(cell_ref_section, cell_par, multiple_unit_cell)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(cell)) CALL cell_create(cell,error=error)
       IF (.NOT.ASSOCIATED(cell_ref)) CALL cell_create(cell_ref,error=error)
       IF (PRESENT(check_for_ref)) my_check = check_for_ref

       cell%deth         = 0.0_dp
       cell%orthorhombic = .FALSE.
       cell%perd(:)      = 1
       cell%symmetry_id  = cell_sym_none
       cell%hmat(:,:)    = 0.0_dp
       cell%h_inv(:,:)   = 0.0_dp
       cell_read_file    = .FALSE.
       cell_read_a       = .FALSE.
       cell_read_b       = .FALSE.
       cell_read_c       = .FALSE.
       ! Trying to read cell info from file
       CALL section_vals_val_get(cell_section,"CELL_FILE_NAME",explicit=cell_read_file,error=error)
       IF (cell_read_file) CALL read_cell_from_external_file(cell_section, para_env, error)

       ! Trying to read cell info from the separate A, B, C vectors
       ! If cell information is provided through file A,B,C contain the file information..
       ! a print warning is shown on screen..
       CALL section_vals_val_get(cell_section,"A",explicit=cell_read_a,error=error)
       IF (cell_read_a) THEN
          CALL section_vals_val_get(cell_section,"A",r_vals=cell_par,error=error)
          cell%hmat(:,1) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"B",explicit=cell_read_b,error=error)
       IF (cell_read_b) THEN
          CALL section_vals_val_get(cell_section,"B",r_vals=cell_par,error=error)
          cell%hmat(:,2) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"C",explicit=cell_read_c,error=error)
       IF (cell_read_c) THEN
          CALL section_vals_val_get(cell_section,"C",r_vals=cell_par,error=error)
          cell%hmat(:,3) = cell_par(:)
       END IF
       check = ((cell_read_a.EQV.cell_read_b).AND.(cell_read_b.EQV.cell_read_c))
       CALL cp_assert(check,cp_warning_level,cp_assertion_failed,routineP,&
            "Cell Information provided through vectors A, B and C. Not all three "//&
            "vectors were provided! Cell setup may be incomplete!"//&
CPSourceFileRef,&
            only_ionode=.TRUE.)

       ! Very last option.. Trying to read cell info from ABC keyword
       CALL section_vals_val_get(cell_section,"ABC",explicit=cell_read_abc,error=error)
       IF (cell_read_abc) THEN
          check = (cell_read_a.OR.cell_read_b.OR.cell_read_c)
          CALL cp_assert(.NOT.check,cp_warning_level,cp_assertion_failed,routineP,&
               "Cell Information provided through vectors A, B and C in conjunction with ABC."//&
               " The definition of the ABC keyword will override the one provided by A,B and C."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          cell%hmat = 0.0_dp
          CALL section_vals_val_get(cell_section,"ABC",r_vals=cell_par,error=error)
          CALL section_vals_val_get(cell_section,"ALPHA_BETA_GAMMA",r_vals=cell_angles,error=error)
          CALL set_cell_param(cell,cell_par,cell_angles,do_init_cell=.FALSE.,error=error)
       END IF

       ! Multiple unit cell
       CALL section_vals_val_get(cell_section,"MULTIPLE_UNIT_CELL",i_vals=multiple_unit_cell,error=error)
       IF (ANY(multiple_unit_cell/=1)) CALL set_multiple_unit_cell(cell, multiple_unit_cell, error)

       CALL section_vals_val_get(cell_section,"PERIODIC",i_val=my_per,error=error)
       SELECT CASE(my_per)
       CASE(use_perd_x)
          cell%perd = (/1,0,0/)
       CASE(use_perd_y)
          cell%perd = (/0,1,0/)
       CASE(use_perd_z)
          cell%perd = (/0,0,1/)
       CASE(use_perd_xy)
          cell%perd = (/1,1,0/)
       CASE(use_perd_xz)
          cell%perd = (/1,0,1/)
       CASE(use_perd_yz)
          cell%perd = (/0,1,1/)
       CASE(use_perd_xyz)
          cell%perd = (/1,1,1/)
       CASE(use_perd_none)
          cell%perd = (/0,0,0/)
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       ! Load requested cell symmetry
       CALL section_vals_val_get(cell_section,"SYMMETRY",i_val=cell%symmetry_id,error=error)

       ! Initialize cell
       CALL init_cell(cell)

       IF (.NOT.my_check) RETURN
       cell_ref_section => section_vals_get_subs_vals(cell_section,&
            "CELL_REF",error=error)
       IF (parsed_cp2k_input(cell_ref_section,check_this_section=.TRUE.,error=error)) THEN
          IF(PRESENT(use_ref_cell) ) use_ref_cell = .TRUE.
          CALL read_cell(cell_ref, cell_ref, use_ref_cell, cell_section=cell_ref_section,&
               check_for_ref=.FALSE., para_env=para_env, error=error)
       ELSE
          CALL cell_clone (cell, cell_ref, error)
          IF ( PRESENT ( use_ref_cell ) ) use_ref_cell = .FALSE.
       END IF
    END IF

  END SUBROUTINE read_cell

! *****************************************************************************
!> \brief utility function to ease the transition to the new input.
!>      returns true if the new input was parsed
!> \param input_file the parsed input file
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION parsed_cp2k_input(input_file,check_this_section,error) RESULT(res)
    TYPE(section_vals_type), POINTER         :: input_file
    LOGICAL, INTENT(IN), OPTIONAL            :: check_this_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'parsed_cp2k_input', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_check
    TYPE(section_vals_type), POINTER         :: glob_section

    failure=.FALSE.
    my_check = .FALSE.
    IF (PRESENT(check_this_section)) my_check = check_this_section
    res=ASSOCIATED(input_file)
    IF (res) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          IF (.NOT.my_check) THEN
             glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",&
                  error=error)
             CALL section_vals_get(glob_section,explicit=res,error=error)
          ELSE
             CALL section_vals_get(input_file,explicit=res,error=error)
          END IF
       END IF
    END IF
  END FUNCTION parsed_cp2k_input

! *****************************************************************************
!> \brief   Setup of the multiple unit_cell
!> \author  Teodoro Laino [tlaino]
!> \date    05.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_multiple_unit_cell(cell, multiple_unit_cell, error)

    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_multiple_unit_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    ! Fail is one of the value is set to zero..
    CALL cp_assert(ALL(multiple_unit_cell>0),cp_fatal_level,cp_assertion_failed,routineP,&
         "CELL%MULTIPLE_UNIT_CELL accepts only integer values larger than 0! "//&
         "A value of 0 or negative is meaningless!"//&
CPSourceFileRef)

    ! scale abc accordingly user request
    cell%hmat(:,1) = cell%hmat(:,1)*multiple_unit_cell(1)
    cell%hmat(:,2) = cell%hmat(:,2)*multiple_unit_cell(2)
    cell%hmat(:,3) = cell%hmat(:,3)*multiple_unit_cell(3)

  END SUBROUTINE set_multiple_unit_cell

! *****************************************************************************
!> \brief   Read cell information from an external file
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \date    02.2008
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_cell_from_external_file(cell_section, para_env, error)

    TYPE(section_vals_type), POINTER         :: cell_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cell_from_external_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: cell_file_name
    INTEGER                                  :: i, idum, j, my_format, n_rep, &
                                                stat
    LOGICAL                                  :: explicit, failure, my_end
    REAL(KIND=dp)                            :: xdum
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    TYPE(cp_parser_type), POINTER            :: parser

    failure = .FALSE.
    NULLIFY(parser)
    CALL section_vals_val_get(cell_section,"CELL_FILE_NAME",c_val=cell_file_name,error=error)
    CALL section_vals_val_get(cell_section,"CELL_FILE_FORMAT",i_val=my_format, error=error)
    CALL parser_create(parser,cell_file_name, para_env=para_env,error=error)
    CALL parser_get_next_line(parser,1,error=error)
    SELECT CASE(my_format)
    CASE (do_cell_cp2k)
       my_end = .FALSE.
       DO WHILE (.NOT.my_end)
          READ(parser%input_line,*)idum,xdum,hmat(:,1),hmat(:,2),hmat(:,3)
          CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
       END DO
    CASE (do_cell_xsc)
       READ(parser%input_line,*)idum,hmat(:,1),hmat(:,2),hmat(:,3)
    END SELECT
    CALL parser_release(parser,error=error)
    CALL section_vals_val_unset(cell_section,"CELL_FILE_NAME",error=error)
    CALL section_vals_val_unset(cell_section,"CELL_FILE_FORMAT",error=error)
    ! Conver to CP2K units
    DO i = 1, 3
       DO j = 1, 3
          hmat(j,i) = cp_unit_to_cp2k(hmat(j,i), "angstrom", error=error)
       END DO
    END DO
    ! Check if the cell was already defined
    explicit = .FALSE.
    CALL section_vals_val_get(cell_section,"A",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"B",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"C",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"ABC",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    ! Possibly print a warning
    CALL cp_assert(.NOT.explicit,cp_warning_level,cp_assertion_failed,routineP,&
         "Cell specification (A,B,C or ABC) provided together with the external "//&
         "cell setup! Ignoring (A,B,C or ABC) and proceeding with info read from the "//&
         "external file! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    ! Copy cell information in the A, B, C fields..(we may need them later on..)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,1)
    CALL section_vals_val_set(cell_section,"A",r_vals_ptr=cell_par,error=error)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,2)
    CALL section_vals_val_set(cell_section,"B",r_vals_ptr=cell_par,error=error)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,3)
    CALL section_vals_val_set(cell_section,"C",r_vals_ptr=cell_par,error=error)
    ! Unset possible keywords
    CALL section_vals_val_unset(cell_section,"ABC",error=error)
    CALL section_vals_val_unset(cell_section,"ALPHA_BETA_GAMMA",error=error)

  END SUBROUTINE read_cell_from_external_file

! *****************************************************************************
!> \brief   Transform real to scaled cell coordinates.
!>          s=h_inv*r
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \author  Matthias Krack
!> \date    16.01.2002
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
!> \brief   Write the cell parameters to the output unit.
!> \author  Matthias Krack
!> \date    02.06.2000
!> \version 1.0
!> \par     History
!>    - 11.2008 Teodoro Laino [tlaino] - rewrite and enabling user driven units
! *****************************************************************************
  RECURSIVE SUBROUTINE write_cell(cell,subsys_section,cell_ref,label,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_ref
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_cell', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_label, unit_str
    INTEGER                                  :: output_unit
    REAL(KIND=dp)                            :: alpha, beta, gamma, val
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    NULLIFY (enum)
    NULLIFY (keyword)
    NULLIFY (logger)
    NULLIFY (section)
    logger => cp_error_get_logger(error)
    my_label = "CELL|"
    IF (PRESENT(label)) my_label = TRIM(label)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%CELL",&
                                       extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%CELL%UNIT",c_val=unit_str,error=error)
    IF (output_unit > 0) THEN
       CALL get_cell(cell=cell,abc=abc,alpha=alpha,beta=beta,gamma=gamma)
       WRITE (UNIT=output_unit, FMT='( )')
       val = cp_unit_from_cp2k(cell%deth,TRIM(unit_str)//"^3",error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.3)")&
        TRIM(my_label)//" Volume ["//TRIM(unit_str)//"^3]:",val
       val = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T30,3F10.3,4X,A6,F11.3)")&
        TRIM(my_label)//" Vector a ["//TRIM(unit_str)//"]:",cell%hmat(:,1)*val,&
        "|a| = ",abc(1)*val,&
        TRIM(my_label)//" Vector b ["//TRIM(unit_str)//"]:",cell%hmat(:,2)*val,&
        "|b| = ",abc(2)*val,&
        TRIM(my_label)//" Vector c ["//TRIM(unit_str)//"]:",cell%hmat(:,3)*val,&
        "|c| = ",abc(3)*val
       WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
        TRIM(my_label)//" Angle (b,c), alpha [degree]: ",alpha,&
        TRIM(my_label)//" Angle (a,c), beta  [degree]: ",beta,&
        TRIM(my_label)//" Angle (a,b), gamma [degree]: ",gamma
        IF (cell%symmetry_id /= cell_sym_none) THEN
           CALL create_cell_section(section,error=error)
           keyword => section_get_keyword(section,"SYMMETRY",error=error)
           CALL keyword_get(keyword,enum=enum,error=error)
           WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
            TRIM(my_label)//" Requested initial symmetry: ",&
            ADJUSTR(TRIM(enum_i2c(enum,cell%symmetry_id,error=error)))
           CALL section_release(section,error=error)
        END IF
        IF (cell%orthorhombic) THEN
           WRITE (UNIT=output_unit,FMT="(T2,A,T78,A3)")&
            TRIM(my_label)//" Numerically orthorhombic: ","YES"
        ELSE
           WRITE (UNIT=output_unit,FMT="(T2,A,T78,A3)")&
            TRIM(my_label)//" Numerically orthorhombic: "," NO"
        END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%CELL",error=error)

    IF (PRESENT(cell_ref)) THEN
       CALL write_cell(cell_ref, subsys_section, label="CELL_REF|", error=error)
    END IF

  END SUBROUTINE write_cell

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
!> \author  RK
!> \date    12.12.2012
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
!> \author  RK
!> \date    28.02.2013
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
