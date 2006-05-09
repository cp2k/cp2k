!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/cell_types [1.0] *
!!
!!   NAME
!!     cell_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (16.01.2002, based on a earlier version of CJM and JGH)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE cell_types
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       use_perd_none, use_perd_x, use_perd_xy, use_perd_xyz, use_perd_xz, &
       use_perd_y, use_perd_yz, use_perd_z, use_unit_angstrom, use_unit_bohr, &
       use_unit_scaled_angstrom, use_unit_scaled_bohr
  USE input_cp2k,                      ONLY: parsed_cp2k_input
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: degree
  USE mathlib,                         ONLY: angle,&
                                             det_3x3,&
                                             inv_3x3
  USE physcon,                         ONLY: bohr
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_types'

  INTEGER, SAVE, PRIVATE :: last_cell_id=0

  TYPE cell_type
     !FM to rm
     CHARACTER(LEN=8)           :: unit_of_length_name
     REAL(KIND = dp)            :: unit_of_length
     LOGICAL                    :: scaled_coordinates

     !FM probably to rm
     REAL(KIND = dp)            :: subcells

     !FM should stay
     REAL(KIND = dp)            :: deth
     LOGICAL                    :: orthorhombic
     INTEGER, DIMENSION(3)      :: perd
     REAL(KIND = dp), DIMENSION(3,3)   :: hmat,h_inv
     INTEGER                    :: ref_count, id_nr
   END TYPE cell_type

   TYPE cell_p_type
      TYPE(cell_type),POINTER :: cell
   END TYPE cell_p_type

! *** Public subroutines ***

  PUBLIC :: get_cell,&
            get_cell_param,&
            get_hinv,&
            init_cell,&
            read_cell,&
            write_cell,&
            cell_create,&
            cell_retain,&
            cell_release,&
            cell_clone

! *** Public functions ***

  PUBLIC :: plane_distance,&
            pbc,&
            real_to_scaled,&
            scaled_to_real

! *** Public data types ***

  PUBLIC :: cell_type,&
            cell_p_type

! *****************************************************************************

  INTERFACE pbc
    MODULE PROCEDURE pbc1,pbc2,pbc3
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE cell_clone (cell_in, cell_out, error)
    TYPE(cell_type), POINTER                 :: cell_in, cell_out
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_clone', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(cell_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(cell_out),cp_failure_level,routineP,error,failure)

    cell_out%unit_of_length_name = cell_in%unit_of_length_name
    cell_out%deth =                cell_in%deth
    cell_out%subcells =            cell_in%subcells
    cell_out%unit_of_length =      cell_in%unit_of_length
    cell_out%perd =                cell_in%perd
    cell_out%hmat =                cell_in%hmat
    cell_out%h_inv =               cell_in%h_inv
    cell_out%orthorhombic =        cell_in%orthorhombic
    cell_out%ref_count = 1
    last_cell_id = last_cell_id + 1
    cell_out%id_nr = last_cell_id
  END SUBROUTINE cell_clone
 
! *****************************************************************************

  SUBROUTINE get_cell(cell,unit_of_length_name,alpha,beta,deth,gamma,&
                      subcells,unit_of_length,orthorhombic,scaled_coordinates,&
                      abc,periodic,h,h_inv, id_nr)

!   Purpose: Get informations about a simulation cell.

!   History: - Creation (16.01.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    CHARACTER(LEN=8), INTENT(OUT), OPTIONAL  :: unit_of_length_name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha, beta, deth, gamma, &
                                                subcells, unit_of_length
    LOGICAL, INTENT(OUT), OPTIONAL           :: orthorhombic, &
                                                scaled_coordinates
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: abc
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: periodic
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT), OPTIONAL                  :: h, h_inv
    INTEGER, INTENT(out), OPTIONAL           :: id_nr

    REAL(KIND=dp)                            :: a, b, c

!   ---------------------------------------------------------------------------

    IF (PRESENT(unit_of_length_name))&
      unit_of_length_name = cell%unit_of_length_name
    IF (PRESENT(unit_of_length)) unit_of_length = cell%unit_of_length
    IF (PRESENT(deth)) deth = cell%deth
    IF (PRESENT(orthorhombic)) orthorhombic = cell%orthorhombic
    IF (PRESENT(scaled_coordinates))&
      scaled_coordinates = cell%scaled_coordinates
    IF (PRESENT(periodic)) periodic(:) = cell%perd(:)
    IF (PRESENT(subcells)) subcells = cell%subcells
    IF (PRESENT(h)) h(:,:) = cell%hmat(:,:)
    IF (PRESENT(h_inv)) h_inv(:,:) = cell%h_inv(:,:)

!   *** Calculate the lengths of the cell vectors a, b, and c ***

    IF (PRESENT(abc)) THEN
      a = SQRT(cell%hmat(1,1)*cell%hmat(1,1) +&
               cell%hmat(2,1)*cell%hmat(2,1) +&
               cell%hmat(3,1)*cell%hmat(3,1))
      b = SQRT(cell%hmat(1,2)*cell%hmat(1,2) +&
               cell%hmat(2,2)*cell%hmat(2,2) +&
               cell%hmat(3,2)*cell%hmat(3,2))
      c = SQRT(cell%hmat(1,3)*cell%hmat(1,3) +&
               cell%hmat(2,3)*cell%hmat(2,3) +&
               cell%hmat(3,3)*cell%hmat(3,3))
      abc(:) = (/a,b,c/)
    END IF

!   *** Angles between the cell vectors a, b, and c   ***

!   *** alpha = <(b,c) ***

    IF (PRESENT(alpha)) alpha = angle(cell%hmat(:,2),cell%hmat(:,3))*degree

!   *** beta = <(a,c) ***

    IF (PRESENT(beta)) beta = angle(cell%hmat(:,1),cell%hmat(:,3))*degree

!   *** gamma = <(a,b) ***

    IF (PRESENT(gamma)) gamma = angle(cell%hmat(:,1),cell%hmat(:,2))*degree

    IF (PRESENT(id_nr)) id_nr=cell%id_nr

  END SUBROUTINE get_cell

! *****************************************************************************

  SUBROUTINE get_cell_param(box,cell_length,cell_angle)

!   Purpose: Compatibility routine.

!   History: - Creation (04.04.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: box
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: cell_length
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: cell_angle

    REAL(KIND=dp)                            :: alpha, beta, gamma
   TYPE(cell_type), POINTER                 :: cell

!   ---------------------------------------------------------------------------

    cell => box

    CALL get_cell(cell=cell,abc=cell_length)

    IF (PRESENT(cell_angle)) THEN
      CALL get_cell(cell=cell,alpha=alpha,beta=beta,gamma=gamma)
      cell_angle(:) = (/alpha,beta,gamma/)
    END IF

  END SUBROUTINE get_cell_param

! *****************************************************************************

  SUBROUTINE get_hinv(box)
!   Purpose: Calculate the inverse of the h matrix.
!   History: - Creation (16.01.2002, CJM,JGH,MK)
    TYPE(cell_type), POINTER                 :: box

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'get_hinv'

    box%deth = ABS ( det_3x3(box%hmat) )

    IF (box%deth < 1.0E-10_dp) THEN
      CALL stop_program(moduleN,routineN,__LINE__,&
                        "An invalid set of box vectors was specified. "//&
                        "The determinant det(h) is too small")
    END IF

    box%h_inv = inv_3x3(box%hmat)

  END SUBROUTINE get_hinv

! *****************************************************************************

  SUBROUTINE init_cell(cell,hmat,periodic,orthorhombic)

!   Purpose: Initialise a simulation cell.

!   History: - Creation (16.01.2002, CJM,JGH,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    LOGICAL, INTENT(IN), OPTIONAL            :: orthorhombic

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_cell (MODULE cell_types)"

    INTEGER                                  :: dim

!   ---------------------------------------------------------------------------

    IF (PRESENT(hmat)) cell%hmat(:,:) = hmat(:,:)
    IF (PRESENT(periodic)) cell%perd(:) = periodic(:)

!   *** Calculate the absolute value of the determinant of the matrix h ***

    cell%deth = ABS ( det_3x3(cell%hmat) )

    IF (cell%deth < 1.0E-10_dp) THEN
      CALL stop_program(routine,&
                        "An invalid set of cell vectors was specified. "//&
                        "The determinant det(h) is too small")
    END IF

    cell%h_inv = inv_3x3(cell%hmat)

    cell%orthorhombic = .FALSE.
    IF ((cell%hmat(1,2) == 0.0_dp).AND.(cell%hmat(1,3) == 0.0_dp).AND.&
        (cell%hmat(2,1) == 0.0_dp).AND.(cell%hmat(2,3) == 0.0_dp).AND.&
        (cell%hmat(3,1) == 0.0_dp).AND.(cell%hmat(3,2) == 0.0_dp)) cell%orthorhombic = .TRUE.
    IF ( PRESENT ( orthorhombic ) ) cell % orthorhombic = orthorhombic
    dim = COUNT(cell%perd == 1)
    IF ((dim/=3).AND.(.NOT.cell%orthorhombic)) THEN
       CALL stop_program(routine,"Not ortho and not periodic")
    END IF

  END SUBROUTINE init_cell

! *****************************************************************************

  FUNCTION plane_distance(h,k,l,cell) RESULT(distance)

!   Purpose: Calculate the distance between two lattice planes as defined by
!            a triple of Miller indices (hkl).

!   History: - Creation (18.11.2004,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: h, k, l
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp)                            :: distance

    REAL(KIND=dp)                            :: a, alpha, b, beta, c, cosa, &
                                                cosb, cosg, d, gamma, x, y, z
    REAL(KIND=dp), DIMENSION(3)              :: abc

!   ---------------------------------------------------------------------------

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

  FUNCTION pbc1(r,cell) RESULT(r_pbc)

!   Purpose: Apply the periodic boundary conditions defined by a simulation
!            cell to a position vector r.

!   History: - Creation (16.01.2002, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc

    REAL(KIND=dp), DIMENSION(3)              :: s

!   ---------------------------------------------------------------------------

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

  FUNCTION pbc2(r,cell,nl) RESULT(r_pbc)

!   Purpose: Apply the periodic boundary conditions defined by a simulation
!            cell to a position vector r.

!   History: - Creation (16.01.2002, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(3), INTENT(IN)        :: nl
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc

    REAL(KIND=dp), DIMENSION(3)              :: s

!   ---------------------------------------------------------------------------

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

  FUNCTION pbc3(ra,rb,cell) RESULT(rab_pbc)

!   Purpose: Apply the periodic boundary conditions defined by the simulation
!            cell cell to the vector pointing from atom a to atom b.

!   History: - Creation (11.03.2004,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: rab_pbc

    INTEGER                                  :: icell, jcell, kcell
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(KIND=dp)                            :: rab2, rab2_pbc
    REAL(KIND=dp), DIMENSION(3)              :: r, ra_pbc, rab, rb_image, &
                                                rb_pbc

!   ---------------------------------------------------------------------------

    CALL get_cell(cell=cell,periodic=periodic)

    ra_pbc(:) = pbc(ra(:),cell)
    rb_pbc(:) = pbc(rb(:),cell)

    rab2_pbc = HUGE(1.0_dp)

    DO icell=-periodic(1),periodic(1)
      DO jcell=-periodic(2),periodic(2)
        DO kcell=-periodic(3),periodic(3)
          r = REAL((/icell,jcell,kcell/),dp)
          rb_image(:) = rb_pbc(:) + scaled_to_real(r,cell)
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

!!****f* cell_types/read_cell_new [1.0] *
!!
!!   NAME
!!     read_cell_new
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     03.2005 created [teo]
!!
!!*** **********************************************************************

  RECURSIVE SUBROUTINE read_cell( cell, cell_ref, use_ref_cell, error, cell_section,&
       check_for_ref)
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(OUT), OPTIONAL           :: use_ref_cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: cell_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_for_ref

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_per, my_val, n_rep
    LOGICAL                                  :: failure, my_check
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    TYPE(section_vals_type), POINTER         :: cell_ref_section

    failure  = .FALSE.
    my_check = .TRUE.
    NULLIFY(cell_ref_section, cell_par)
    IF (.NOT.failure) THEN

       IF (.NOT.ASSOCIATED(cell)) THEN
          CALL cell_create(cell)
       END IF
       
       IF (.NOT.ASSOCIATED(cell_ref)) THEN
          CALL cell_create(cell_ref)
       END IF
       
       IF (PRESENT(check_for_ref)) my_check     =  check_for_ref
       
       cell%unit_of_length_name = "BOHR"
       cell%deth = 0.0_dp
       cell%unit_of_length = 1.0_dp
       cell%orthorhombic = .TRUE.
       cell%scaled_coordinates = .FALSE.
       cell%perd(:) = 1
       cell%hmat(:,:) = 0.0_dp
       cell%h_inv(:,:) = 0.0_dp
       
       CALL section_vals_val_get(cell_section,"A",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(cell_section,"A",r_vals=cell_par,error=error)
          cell%hmat(:,1) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"B",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(cell_section,"B",r_vals=cell_par,error=error)
          cell%hmat(:,2) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"C",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(cell_section,"C",r_vals=cell_par,error=error)
          cell%hmat(:,3) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"ABC",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(cell_section,"ABC",r_vals=cell_par,error=error)
          cell%hmat(1,1) = cell_par(1)
          cell%hmat(2,2) = cell_par(2)
          cell%hmat(3,3) = cell_par(3)
       END IF
       CALL section_vals_val_get(cell_section,"SUBCELLS",r_val=cell%subcells,error=error)
       
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

       CALL section_vals_val_get(cell_section,"UNITS",i_val=my_val,error=error)
       SELECT CASE (my_val)
       CASE (use_unit_angstrom)
          cell%unit_of_length = bohr
          cell%unit_of_length_name = "ANGSTROM"
          cell%scaled_coordinates = .FALSE.
       CASE (use_unit_bohr)
          cell%unit_of_length = 1.0_dp
          cell%unit_of_length_name = "BOHR"
          cell%scaled_coordinates = .FALSE.
       CASE (use_unit_scaled_angstrom)
          cell%unit_of_length = bohr
          cell%unit_of_length_name = "ANGSTROM"
          cell%scaled_coordinates = .TRUE.
       CASE (use_unit_scaled_bohr)
          cell%unit_of_length = 1.0_dp
          cell%unit_of_length_name = "BOHR"
          cell%scaled_coordinates = .TRUE.
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       cell%hmat(:,:) = cell%hmat(:,:)*cell%unit_of_length      
       CALL init_cell(cell)
       
       IF (.NOT.my_check) RETURN
       cell_ref_section => section_vals_get_subs_vals(cell_section,&
            "CELL_REF",error=error)
       IF (parsed_cp2k_input(cell_ref_section,check_this_section=.TRUE.,error=error)) THEN
          IF(PRESENT(use_ref_cell) ) use_ref_cell = .TRUE.
          CALL read_cell(cell_ref, cell_ref, use_ref_cell, error=error,&
               cell_section=cell_ref_section, check_for_ref=.FALSE.)
       ELSE
          CALL cell_clone (cell, cell_ref, error)
          IF ( PRESENT ( use_ref_cell ) ) use_ref_cell = .FALSE.
       END IF
    END IF
  END SUBROUTINE read_cell

! *****************************************************************************

  FUNCTION real_to_scaled(r,cell) RESULT(s)

!   Purpose: Transform real to scaled cell coordinates.

!   History: - Creation (16.01.2002, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: s

!   ---------------------------------------------------------------------------

    IF (cell%orthorhombic) THEN
      s(1) = cell%h_inv(1,1)*r(1)
      s(2) = cell%h_inv(2,2)*r(2)
      s(3) = cell%h_inv(3,3)*r(3)
    ELSE
      s(1) = cell%h_inv(1,1)*r(1) + cell%h_inv(1,2)*r(2) + cell%h_inv(1,3)*r(3)
      s(2) = cell%h_inv(2,1)*r(1) + cell%h_inv(2,2)*r(2) + cell%h_inv(2,3)*r(3)
      s(3) = cell%h_inv(3,1)*r(1) + cell%h_inv(3,2)*r(2) + cell%h_inv(3,3)*r(3)
    END IF

  END FUNCTION real_to_scaled

! *****************************************************************************

  FUNCTION scaled_to_real(s,cell) RESULT(r)

!   Purpose: Transform scaled cell coordinates real coordinates.

!   History: - Creation (16.01.2002, Matthias Krack)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: s
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: r

!   ---------------------------------------------------------------------------

    IF (cell%orthorhombic) THEN
      r(1) = cell%hmat(1,1)*s(1)
      r(2) = cell%hmat(2,2)*s(2)
      r(3) = cell%hmat(3,3)*s(3)
    ELSE
      r(1) = cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
      r(2) = cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
      r(3) = cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
    END IF

  END FUNCTION scaled_to_real

! *****************************************************************************

  SUBROUTINE write_cell(cell,force_env_section,cell_ref,error)

!   Purpose: Write the cell parameters to the output unit.

!   History: - Creation (02.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_ref
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: output_unit
    REAL(KIND=dp)                            :: alpha, beta, gamma
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    output_unit = cp_print_key_unit_nr(logger,force_env_section,"SUBSYS%PRINT%CELL",&
         extension=".Log",error=error)
    IF (output_unit>0) THEN

       CALL get_cell(cell=cell,abc=abc,alpha=alpha,beta=beta,gamma=gamma)

     WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.3)")&
          "CELL| Volume ["//TRIM(cell%unit_of_length_name)//"**3]:",&
          cell%deth/cell%unit_of_length**3
     WRITE (UNIT=output_unit,FMT="(T2,A,T30,3F10.3,4X,A6,F11.3)")&
          "CELL| Vector a ["//TRIM(cell%unit_of_length_name)//"]:",&
          cell%hmat(:,1)/cell%unit_of_length,&
          "|a| = ",abc(1)/cell%unit_of_length,&
          "CELL| Vector b ["//TRIM(cell%unit_of_length_name)//"]:",&
          cell%hmat(:,2)/cell%unit_of_length,&
          "|b| = ",abc(2)/cell%unit_of_length,&
          "CELL| Vector c ["//TRIM(cell%unit_of_length_name)//"]:",&
          cell%hmat(:,3)/cell%unit_of_length,&
          "|c| = ",abc(3)/cell%unit_of_length
     WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
          "CELL| Angle (b,c) [degree]: ",alpha,&
          "CELL| Angle (a,c) [degree]: ",beta,&
          "CELL| Angle (a,b) [degree]: ",gamma
     WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
          "CELL| Grid size for subcell generation",cell%subcells
  
     IF ( PRESENT ( cell_ref ) ) THEN
        WRITE (UNIT=output_unit, FMT='( )' ) 
        CALL get_cell(cell=cell_ref,abc=abc,alpha=alpha,beta=beta,gamma=gamma)

        WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.3)")&
          "CELL_REF| Volume ["//TRIM(cell_ref%unit_of_length_name)//"**3]:",&
          cell_ref%deth/cell%unit_of_length**3
        WRITE (UNIT=output_unit,FMT="(T2,A,T30,3F10.3,4X,A6,F11.3)")&
          "CELL_REF| Vector a ["//TRIM(cell_ref%unit_of_length_name)//"]:",&
          cell_ref%hmat(:,1)/cell_ref%unit_of_length,&
          "|a| = ",abc(1)/cell_ref%unit_of_length,&
          "CELL_REF| Vector b ["//TRIM(cell_ref%unit_of_length_name)//"]:",&
          cell_ref%hmat(:,2)/cell_ref%unit_of_length,&
          "|b| = ",abc(2)/cell_ref%unit_of_length,&
          "CELL_REF| Vector c ["//TRIM(cell_ref%unit_of_length_name)//"]:",&
          cell_ref%hmat(:,3)/cell_ref%unit_of_length,&
          "|c| = ",abc(3)/cell_ref%unit_of_length
        WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
          "CELL_REF| Angle (b,c) [degree]: ",alpha,&
          "CELL_REF| Angle (a,c) [degree]: ",beta,&
          "CELL_REF| Angle (a,b) [degree]: ",gamma
        WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
          "CELL_REF| Grid size for subcell generation",cell_ref%subcells
      ENDIF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
       "SUBSYS%PRINT%CELL",error=error)

  END SUBROUTINE write_cell

! *****************************************************************************

!!****f* cell_types/cell_create [1.0] *
!!
!!   NAME
!!     cell_create
!!
!!   FUNCTION
!!     allocates and initializes a cell
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - cell the cell to initialize
!!     - hmat: the h matrix that defines the cell
!!     - periodic: periodicity of the cell
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cell_create(cell,hmat,periodic, orthorhombic, error)
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    LOGICAL, INTENT(IN), OPTIONAL            :: orthorhombic
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    ALLOCATE(cell,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_cell_id=last_cell_id+1
    cell%id_nr=last_cell_id
    cell%ref_count=1
    cell % perd = 1
    IF (.NOT. PRESENT(hmat)) RETURN
    IF (.NOT. failure) THEN
       CALL init_cell(cell,hmat,periodic,orthorhombic)
    END IF
  END SUBROUTINE cell_create
!***************************************************************************

!!****f* cell_types/cell_retain [1.0] *
!!
!!   NAME
!!     cell_retain
!!
!!   FUNCTION
!!     retains the given cell (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - cell: the cell to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cell_retain(cell,error)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(cell%ref_count>0,cp_failure_level,routineP,error)
     cell%ref_count=cell%ref_count+1
  END IF
END SUBROUTINE cell_retain
!***************************************************************************

!!****f* cell_types/cell_release [1.0] *
!!
!!   NAME
!!     cell_release
!!
!!   FUNCTION
!!     releases the given cell (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - cell: the cell to release
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cell_release(cell,error)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(cell)) THEN
     CPPreconditionNoFail(cell%ref_count>0,cp_failure_level,routineP,error)
     cell%ref_count=cell%ref_count-1
     IF (cell%ref_count==0) THEN
        DEALLOCATE(cell,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE cell_release

!***************************************************************************

END MODULE cell_types
