!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/simulation_cell [1.0] *
!!
!!   NAME
!!     simulation_cell
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

MODULE simulation_cell
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE greens_fn,                       ONLY: ANALYTIC0D,&
                                             ANALYTIC1D,&
                                             ANALYTIC2D,&
                                             HOCKNEY0D,&
                                             HOCKNEY1D,&
                                             HOCKNEY2D,&
                                             MT0D,&
                                             MT1D,&
                                             MT2D,&
                                             PERIODIC3D,&
                                             greens_function_type,&
                                             pw_green_fn_destroy
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: degree
  USE physcon,                         ONLY: bohr
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='simulation_cell'
  INTEGER, SAVE, PRIVATE :: last_cell_id=0

  TYPE cell_type
!MK  PRIVATE works currently only for the QS modules
     TYPE(greens_function_type) :: green
     CHARACTER(LEN=8)           :: unit_of_length_name
     REAL(KIND = dp)                   :: deth,subcells,unit_of_length
     LOGICAL                    :: orthorhombic,scaled_coordinates
     INTEGER, DIMENSION(3)      :: perd
     REAL(KIND = dp), DIMENSION(3,3)   :: hmat,h_inv
     INTEGER                    :: ref_count, id_nr
   END TYPE cell_type

! *** Public subroutines ***

  PUBLIC :: get_cell,&
            get_cell_param,&
            get_hinv,&
            init_cell,&
            read_cell,&
            write_cell,&
            cell_create,&
            cell_retain,&
            cell_release

! *** Public functions ***

  PUBLIC :: pbc,&
            real_to_scaled,&
            scaled_to_real

! *** Public data types ***

  PUBLIC :: cell_type

! *****************************************************************************

  INTERFACE pbc
    MODULE PROCEDURE pbc1,pbc2
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_cell(cell,unit_of_length_name,alpha,beta,deth,gamma,&
                      subcells,unit_of_length,orthorhombic,scaled_coordinates,&
                      abc,periodic,h,h_inv, id_nr)

!   Purpose: Get informations about a simulation cell.

!   History: - Creation (16.01.2002, Matthias Krack)

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

    IF (PRESENT(abc).OR.PRESENT(alpha).OR.PRESENT(beta).OR.PRESENT(gamma)) THEN
      a = SQRT(cell%hmat(1,1)*cell%hmat(1,1) +&
               cell%hmat(2,1)*cell%hmat(2,1) +&
               cell%hmat(3,1)*cell%hmat(3,1))
      b = SQRT(cell%hmat(1,2)*cell%hmat(1,2) +&
               cell%hmat(2,2)*cell%hmat(2,2) +&
               cell%hmat(3,2)*cell%hmat(3,2))
      c = SQRT(cell%hmat(1,3)*cell%hmat(1,3) +&
               cell%hmat(2,3)*cell%hmat(2,3) +&
               cell%hmat(3,3)*cell%hmat(3,3))
    END IF

    IF (PRESENT(abc)) abc(:) = (/a,b,c/)

!   *** Angles between the cell vectors a, b, and c   ***

!   *** alpha = <(b,c) ***

    IF (PRESENT(alpha)) THEN
      alpha = ACOS((cell%hmat(1,2)*cell%hmat(1,3) +&
                    cell%hmat(2,2)*cell%hmat(2,3) +&
                    cell%hmat(3,2)*cell%hmat(3,3))/(b*c))*degree
    END IF

!   *** beta = <(a,c) ***

    IF (PRESENT(beta)) THEN
      beta = ACOS((cell%hmat(1,1)*cell%hmat(1,3) +&
                   cell%hmat(2,1)*cell%hmat(2,3) +&
                   cell%hmat(3,1)*cell%hmat(3,3))/(a*c))*degree
    END IF

!   *** gamma = <(a,b) ***

    IF (PRESENT(gamma)) THEN
      gamma = ACOS((cell%hmat(1,1)*cell%hmat(1,2) +&
                    cell%hmat(2,1)*cell%hmat(2,2) +&
                    cell%hmat(3,1)*cell%hmat(3,2))/(a*b))*degree
    END IF

    IF (PRESENT(id_nr)) id_nr=cell%id_nr

  END SUBROUTINE get_cell

! *****************************************************************************

  SUBROUTINE get_cell_param(box,cell_length,cell_angle)

!   Purpose: Compatibility routine.

!   History: - Creation (04.04.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), INTENT(IN), TARGET      :: box
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: cell_length
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: cell_angle

    REAL(KIND=dp)                            :: alpha, beta, gamm
    TYPE(cell_type), POINTER                 :: cell

!   ---------------------------------------------------------------------------

    cell => box

    CALL get_cell(cell=cell,abc=cell_length)

    IF (PRESENT(cell_angle)) THEN
      CALL get_cell(cell=cell,alpha=alpha,beta=beta,gamma=gamm)
      cell_angle(:) = (/alpha,beta,gamm/)
    END IF

  END SUBROUTINE get_cell_param

! *****************************************************************************

  SUBROUTINE get_hinv(box)

!   Purpose: Calculate the inverse of the h matrix.

!   History: - Creation (16.01.2002, CJM,JGH,MK)

!   ***************************************************************************

    TYPE(cell_type), INTENT(INOUT)           :: box

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE get_hinv (MODULE simulation_cell)"

    REAL(KIND=dp)                            :: dethinv
    REAL(KIND=dp), DIMENSION(3, 3)           :: h

!   ---------------------------------------------------------------------------

    h(:,:) = box%hmat(:,:)

!   *** Calculate the determinant of the matrix h ***

    box%deth = h(1,1)*(h(2,2)*h(3,3) - h(2,3)*h(3,2)) + &
               h(1,2)*(h(2,3)*h(3,1) - h(2,1)*h(3,3)) + &
               h(1,3)*(h(2,1)*h(3,2) - h(2,2)*h(3,1))

    IF (box%deth < 1.0E-10_dp) THEN
      CALL stop_program(routine,&
                        "An invalid set of box vectors was specified. "//&
                        "The determinant det(h) is too small")
    END IF

    dethinv = 1.0_dp/box%deth

    box%h_inv(1,1) = (h(2,2)*h(3,3) - h(3,2)*h(2,3))*dethinv
    box%h_inv(2,1) = (h(2,3)*h(3,1) - h(3,3)*h(2,1))*dethinv
    box%h_inv(3,1) = (h(2,1)*h(3,2) - h(3,1)*h(2,2))*dethinv
    box%h_inv(1,2) = (h(1,3)*h(3,2) - h(3,3)*h(1,2))*dethinv
    box%h_inv(2,2) = (h(1,1)*h(3,3) - h(3,1)*h(1,3))*dethinv
    box%h_inv(3,2) = (h(1,2)*h(3,1) - h(3,2)*h(1,1))*dethinv
    box%h_inv(1,3) = (h(1,2)*h(2,3) - h(2,2)*h(1,3))*dethinv
    box%h_inv(2,3) = (h(1,3)*h(2,1) - h(2,3)*h(1,1))*dethinv
    box%h_inv(3,3) = (h(1,1)*h(2,2) - h(2,1)*h(1,2))*dethinv

  END SUBROUTINE get_hinv

! *****************************************************************************

  SUBROUTINE init_cell(box,hmat,periodic,poisson_solver)

!   Purpose: Compatibility routine.

!   History: - Creation (04.04.2002,MK)

!   ***************************************************************************

    TYPE(cell_type), INTENT(OUT), TARGET     :: box
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    INTEGER, DIMENSION(3), INTENT(IN)        :: periodic
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: poisson_solver

    TYPE(cell_type), POINTER                 :: cell

!   ---------------------------------------------------------------------------

    cell => box

    IF (PRESENT(poisson_solver)) THEN
      CALL init_cell_new(cell=cell,&
                         hmat=hmat,&
                         periodic=periodic,&
                         poisson_solver=poisson_solver)
    ELSE
      CALL init_cell_new(cell=cell,&
                         hmat=hmat,&
                         periodic=periodic)
    END IF

  END SUBROUTINE init_cell

! *****************************************************************************

  SUBROUTINE init_cell_new(cell,hmat,periodic,poisson_solver)

!   Purpose: Initialise a simulation cell.

!   History: - Creation (16.01.2002, CJM,JGH,MK)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: poisson_solver

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_cell_new (MODULE simulation_cell)"

    INTEGER                                  :: dim, i
    REAL(KIND=dp)                            :: dethinv
    REAL(KIND=dp), DIMENSION(3)              :: abc
    REAL(KIND=dp), DIMENSION(3, 3)           :: h

!   ---------------------------------------------------------------------------

    IF (PRESENT(hmat)) cell%hmat(:,:) = hmat(:,:)
    IF (PRESENT(periodic)) cell%perd(:) = periodic(:)

    last_cell_id=last_cell_id+1
    cell%id_nr=last_cell_id
    cell%ref_count=1

    h(:,:) = cell%hmat(:,:)

!   *** Calculate the determinant of the matrix h ***

    cell%deth = h(1,1)*(h(2,2)*h(3,3) - h(2,3)*h(3,2)) + &
                h(1,2)*(h(2,3)*h(3,1) - h(2,1)*h(3,3)) + &
                h(1,3)*(h(2,1)*h(3,2) - h(2,2)*h(3,1))

    IF (cell%deth < 1.0E-10_dp) THEN
      CALL stop_program(routine,&
                        "An invalid set of cell vectors was specified. "//&
                        "The determinant det(h) is too small")
    END IF

    dethinv = 1.0_dp/cell%deth

    cell%h_inv(1,1) = (h(2,2)*h(3,3) - h(3,2)*h(2,3))*dethinv
    cell%h_inv(2,1) = (h(2,3)*h(3,1) - h(3,3)*h(2,1))*dethinv
    cell%h_inv(3,1) = (h(2,1)*h(3,2) - h(3,1)*h(2,2))*dethinv
    cell%h_inv(1,2) = (h(1,3)*h(3,2) - h(3,3)*h(1,2))*dethinv
    cell%h_inv(2,2) = (h(1,1)*h(3,3) - h(3,1)*h(1,3))*dethinv
    cell%h_inv(3,2) = (h(1,2)*h(3,1) - h(3,2)*h(1,1))*dethinv
    cell%h_inv(1,3) = (h(1,2)*h(2,3) - h(2,2)*h(1,3))*dethinv
    cell%h_inv(2,3) = (h(1,3)*h(2,1) - h(2,3)*h(1,1))*dethinv
    cell%h_inv(3,3) = (h(1,1)*h(2,2) - h(2,1)*h(1,2))*dethinv

    IF ((cell%hmat(1,2) == 0.0_dp).AND.(cell%hmat(1,3) == 0.0_dp).AND.&
        (cell%hmat(2,1) == 0.0_dp).AND.(cell%hmat(2,3) == 0.0_dp).AND.&
        (cell%hmat(3,1) == 0.0_dp).AND.(cell%hmat(3,2) == 0.0_dp)) THEN
      cell%orthorhombic = .TRUE.
    ELSE
      cell%orthorhombic = .FALSE.
    END IF

    CALL get_cell(cell=cell,abc=abc)

    cell%green%special_dimension = 0
    cell%green%radius = 0.0_dp
    cell%green%slab_size = 0.0_dp
    cell%green%alpha = 0.0_dp
    cell%green%volume = 0.0_dp
    NULLIFY( cell%green%influence_function, cell%green%p3m_charge )
    NULLIFY( cell%green%p3m_coeff, cell%green%p3m_bm2 )

    dim = COUNT(cell%perd == 1)

    SELECT CASE (dim)
    CASE (0)
      IF (.NOT.cell%orthorhombic) THEN
        CALL stop_program(routine,"Not ortho and not periodic")
      END IF
      cell%green%method = ANALYTIC0D
      cell%green%radius = 0.5_dp*MINVAL(abc)
    CASE (1)
      IF (.NOT.cell%orthorhombic) THEN
        CALL stop_program(routine,"Not ortho and not periodic")
      END IF
      cell%green%method = ANALYTIC1D
      cell%green%special_dimension = MAXLOC(cell%perd,1)
      cell%green%radius = MAXVAL(abc)
      DO i=1,3
        IF (i == cell%green%special_dimension) CYCLE
        cell%green%radius = MIN(cell%green%radius,0.5_dp*abc(i))
      END DO
    CASE (2)
      IF (.NOT.cell%orthorhombic) THEN
        CALL stop_program(routine,"Not ortho and not periodic")
      END IF
      cell%green%method = ANALYTIC2D
      i = MINLOC(cell%perd,1)
      cell%green%special_dimension = i
      cell%green%slab_size = abc(i)
     CASE (3)
      cell%green%method = periodic3D
    CASE DEFAULT
      CALL stop_program(routine,"Illegal value for perd found")
    END SELECT

    IF (PRESENT(poisson_solver)) THEN

      SELECT CASE (TRIM(poisson_solver))
      CASE ("PERIODIC")
        cell%green%method = PERIODIC3D
        IF (dim /= 3) THEN
          CALL stop_program(routine,&
                            "Illegal combination of periodicity and "//&
                            "Poisson solver (perd)")
        END IF
      CASE ("ANALYTIC")
        SELECT CASE (dim)
        CASE (0)
          cell%green%method = ANALYTIC0D
        CASE (1)
          cell%green%method = ANALYTIC1D
        CASE (2)
          cell%green%method = ANALYTIC2D
        CASE (3)
          cell%green%method = PERIODIC3D
        CASE DEFAULT
          CALL stop_program(routine,"Illegal value for perd found")
        END SELECT
      CASE ("HOCKNEY")
        SELECT CASE (dim )
        CASE (0)
          cell%green%method = HOCKNEY0D
        CASE (1)
          cell%green%method = HOCKNEY1D
        CASE (2)
          cell%green%method = HOCKNEY2D
        CASE (3)
          CALL stop_program(routine,&
                            "Illegal combination of periodicity and "//&
                            "Poisson solver (HOCKNEY)")
        CASE DEFAULT
          CALL stop_program(routine,"Illegal value for perd found")
        END SELECT
      CASE ("MT")
        SELECT CASE (dim)
        CASE (0)
          cell%green%method = MT0D
        CASE (1)
          cell%green%method = MT1D
        CASE (2)
          cell%green%method = MT2D
        CASE (3)
          CALL stop_program(routine,&
                            "Illegal combination of periodicity and "//&
                            "Poisson solver (MT)")
        CASE DEFAULT
          CALL stop_program(routine,"Illegal value for perd found")
        END SELECT
      CASE DEFAULT
        CALL stop_program(routine,"An unknown Poisson solver was specified")
      END SELECT

    END IF

  END SUBROUTINE init_cell_new

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

  SUBROUTINE read_cell(cell,globenv,new_h_matrix)

!   Purpose: Read the cell parameters from the input file.

!   History: - Creation (19.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    REAL(dp), DIMENSION(1:3, 1:3), &
      INTENT(IN), OPTIONAL                   :: new_h_matrix

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE read_cell (MODULE simulation_cell)"

    CHARACTER(LEN=40)                        :: keyword, string
    INTEGER                                  :: istat
    LOGICAL                                  :: section_found

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(cell)) THEN
      ALLOCATE (cell,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cell",0)
    END IF

!   *** Load the default values ***

    cell%unit_of_length_name = "BOHR"
    cell%deth = 0.0_dp
    cell%unit_of_length = 1.0_dp
    cell%orthorhombic = .TRUE.
    cell%scaled_coordinates = .FALSE.
    cell%subcells = 2.0_dp
    cell%perd(:) = 1
    cell%hmat(:,:) = 0.0_dp
    cell%h_inv(:,:) = 0.0_dp

!   *** Read the input section ***

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="CELL",&
                      section_found=section_found)

    IF (.NOT.section_found) CALL stop_program(routine,"No CELL section found.")

    DO WHILE (test_object(newline=.TRUE.) /= "EOS")
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("A")
        CALL read_object(cell%hmat(1,1))
        CALL read_object(cell%hmat(2,1))
        CALL read_object(cell%hmat(3,1))
      CASE ("ABC")
        CALL read_object(cell%hmat(1,1))
        CALL read_object(cell%hmat(2,2))
        CALL read_object(cell%hmat(3,3))
      CASE ("B")
        CALL read_object(cell%hmat(1,2))
        CALL read_object(cell%hmat(2,2))
        CALL read_object(cell%hmat(3,2))
      CASE ("C")
        CALL read_object(cell%hmat(1,3))
        CALL read_object(cell%hmat(2,3))
        CALL read_object(cell%hmat(3,3))
      CASE ("SUBCELLS")
        CALL read_object(cell%subcells)
      CASE ("UNIT","UNITS")
        CALL read_object(string,lower_to_upper=.TRUE.)
        SELECT CASE (TRIM(string))
        CASE ("ANGSTROM")
          cell%unit_of_length = bohr
          cell%unit_of_length_name = "ANGSTROM"
          cell%scaled_coordinates = .FALSE.
        CASE ("BOHR")
          cell%unit_of_length = 1.0_dp
          cell%unit_of_length_name = "BOHR"
          cell%scaled_coordinates = .FALSE.
        CASE ("SCALED_ANGSTROM")
          cell%unit_of_length = bohr
          cell%unit_of_length_name = "ANGSTROM"
          cell%scaled_coordinates = .TRUE.
        CASE ("SCALED_BOHR")
          cell%unit_of_length = 1.0_dp
          cell%unit_of_length_name = "BOHR"
          cell%scaled_coordinates = .TRUE.
        CASE DEFAULT
          CALL stop_parser(routine,&
                           "Invalid name for the cell unit of length <"//&
                           TRIM(string)//"> found")
        END SELECT
      CASE DEFAULT
          CALL stop_parser(routine,"Invalid keyword found")
      END SELECT
    END DO

    CALL finish_parser()

    cell%hmat(:,:) = cell%hmat(:,:)*cell%unit_of_length

    IF(PRESENT(new_h_matrix)) cell%hmat(:,:) = new_h_matrix(:,:)

!   *** Initialise the remaining components of the cell_type ***

    CALL init_cell_new(cell)

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

  SUBROUTINE write_cell(cell,globenv)

!   Purpose: Write the cell parameters to the output unit.

!   History: - Creation (02.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: output_unit
    REAL(KIND=dp)                            :: alpha, beta, gamm
    REAL(KIND=dp), DIMENSION(3)              :: abc

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%cell_parameters) THEN

      output_unit = globenv%scr

      CALL get_cell(cell=cell,abc=abc,alpha=alpha,beta=beta,gamma=gamm)

      WRITE (UNIT=output_unit,FMT="(/,T2,A,T61,F20.3)")&
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
        "CELL| Angle (a,b) [degree]: ",gamm
      WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
        "CELL| Grid size for subcell generation",cell%subcells

    END IF

  END SUBROUTINE write_cell

! *****************************************************************************

!!****f* simulation_cell/cell_create [1.0] *
!!
!!   NAME
!!     cell_create
!!
!!   SYNOPSIS
!!     Subroutine cell_create(cell, hmat, periodic, poisson_solver, error)
!!       Type(cell_type), Pointer:: cell
!!       Real(KIND = dp), Dimension(3,3), Optional, Intent (IN):: hmat
!!       Integer, Dimension(3), Optional, Intent (IN):: periodic
!!       Character(Len=*), Optional, Intent (IN):: poisson_solver
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cell_create
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
!!     - poisson_solver: the poisson solver to use
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
  SUBROUTINE cell_create(cell,hmat,periodic,poisson_solver, error)
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN), OPTIONAL                   :: hmat
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: periodic
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: poisson_solver
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
    IF (.NOT. failure) THEN
       CALL init_cell_new(cell,hmat,periodic,poisson_solver)
    END IF
  END SUBROUTINE cell_create
!***************************************************************************

!!****f* simulation_cell/cell_retain [1.0] *
!!
!!   NAME
!!     cell_retain
!!
!!   SYNOPSIS
!!     Subroutine cell_retain(cell, error)
!!       Type(cell_type), Pointer:: cell
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cell_retain
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

!!****f* simulation_cell/cell_release [1.0] *
!!
!!   NAME
!!     cell_release
!!
!!   SYNOPSIS
!!     Subroutine cell_release(cell, error)
!!       Type(cell_type), Pointer:: cell
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cell_release
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
        CALL pw_green_fn_destroy( cell % green )
        DEALLOCATE(cell,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE cell_release
!***************************************************************************

END MODULE simulation_cell
