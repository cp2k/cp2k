!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!>      10.2014 Moved many routines from cell_types.F here.
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! **************************************************************************************************
MODULE cell_types
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: degree
   USE mathlib,                         ONLY: angle
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_types'

   ! Impose cell symmetry
   INTEGER, PARAMETER, PUBLIC               :: cell_sym_none = 0, &
                                               cell_sym_triclinic = 1, &
                                               cell_sym_monoclinic = 2, &
                                               cell_sym_monoclinic_gamma_ab = 3, &
                                               cell_sym_orthorhombic = 4, &
                                               cell_sym_tetragonal_ab = 5, &
                                               cell_sym_tetragonal_ac = 6, &
                                               cell_sym_tetragonal_bc = 7, &
                                               cell_sym_rhombohedral = 8, &
                                               cell_sym_hexagonal_gamma_60 = 9, &
                                               cell_sym_hexagonal_gamma_120 = 10, &
                                               cell_sym_cubic = 11

   INTEGER, PARAMETER, PUBLIC               :: use_perd_x = 0, &
                                               use_perd_y = 1, &
                                               use_perd_z = 2, &
                                               use_perd_xy = 3, &
                                               use_perd_xz = 4, &
                                               use_perd_yz = 5, &
                                               use_perd_xyz = 6, &
                                               use_perd_none = 7

! **************************************************************************************************
!> \brief   Type defining parameters related to the simulation cell
!> \version 1.0
! **************************************************************************************************
   TYPE cell_type
      CHARACTER(LEN=12)                 :: tag = "CELL"
      INTEGER                           :: ref_count = -1, &
                                           symmetry_id = use_perd_none
      LOGICAL                           :: orthorhombic = .FALSE. ! actually means a diagonal hmat
      REAL(KIND=dp)                     :: deth = 0.0_dp
      INTEGER, DIMENSION(3)             :: perd = -1
      REAL(KIND=dp), DIMENSION(3, 3)    :: hmat = 0.0_dp, &
                                           h_inv = 0.0_dp
   END TYPE cell_type

   TYPE cell_p_type
      TYPE(cell_type), POINTER :: cell => NULL()
   END TYPE cell_p_type

   ! Public data types
   PUBLIC :: cell_type, &
             cell_p_type

   ! Public subroutines
   PUBLIC :: cell_clone, &
             cell_copy, &
             cell_release, &
             cell_retain, &
             get_cell, &
             parse_cell_line

#if defined (__PLUMED2)
   PUBLIC :: pbc_cp2k_plumed_getset_cell
#endif

   ! Public functions
   PUBLIC :: plane_distance, &
             pbc, &
             real_to_scaled, &
             scaled_to_real

   INTERFACE pbc
      MODULE PROCEDURE pbc1, pbc2, pbc3, pbc4
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Clone cell variable
!> \param cell_in Cell variable to be clone
!> \param cell_out Cloned cell variable
!> \param tag Optional new tag for cloned cell variable
!> \par History
!>      - Optional tag added (17.05.2023, MK)
! **************************************************************************************************
   SUBROUTINE cell_clone(cell_in, cell_out, tag)

      TYPE(cell_type), POINTER                           :: cell_in, cell_out
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: tag

      cell_out = cell_in
      cell_out%ref_count = 1
      IF (PRESENT(tag)) cell_out%tag = tag

   END SUBROUTINE cell_clone

! **************************************************************************************************
!> \brief Copy cell variable
!> \param cell_in Cell variable to be copied
!> \param cell_out Copy of cell variable
!> \param tag Optional new tag
!> \par History
!>      - Optional tag added (17.05.2023, MK)
! **************************************************************************************************
   SUBROUTINE cell_copy(cell_in, cell_out, tag)

      TYPE(cell_type), POINTER                           :: cell_in, cell_out
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: tag

      cell_out%deth = cell_in%deth
      cell_out%perd = cell_in%perd
      cell_out%hmat = cell_in%hmat
      cell_out%h_inv = cell_in%h_inv
      cell_out%orthorhombic = cell_in%orthorhombic
      cell_out%symmetry_id = cell_in%symmetry_id
      IF (PRESENT(tag)) THEN
         cell_out%tag = tag
      ELSE
         cell_out%tag = cell_in%tag
      END IF

   END SUBROUTINE cell_copy

! **************************************************************************************************
!> \brief   Read cell info from a line (parsed from a file)
!> \param input_line ...
!> \param cell_itimes ...
!> \param cell_time ...
!> \param h ...
!> \param vol ...
!> \date    19.02.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE parse_cell_line(input_line, cell_itimes, cell_time, h, vol)

      CHARACTER(LEN=*), INTENT(IN)                       :: input_line
      INTEGER, INTENT(OUT)                               :: cell_itimes
      REAL(KIND=dp), INTENT(OUT)                         :: cell_time
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: h
      REAL(KIND=dp), INTENT(OUT)                         :: vol

      INTEGER                                            :: i, j

      READ (input_line, *) cell_itimes, cell_time, &
         h(1, 1), h(2, 1), h(3, 1), h(1, 2), h(2, 2), h(3, 2), h(1, 3), h(2, 3), h(3, 3), vol
      DO i = 1, 3
         DO j = 1, 3
            h(j, i) = cp_unit_to_cp2k(h(j, i), "angstrom")
         END DO
      END DO

   END SUBROUTINE parse_cell_line

! **************************************************************************************************
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
!> \param symmetry_id ...
!> \param tag ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_cell(cell, alpha, beta, gamma, deth, orthorhombic, abc, periodic, &
                       h, h_inv, symmetry_id, tag)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: alpha, beta, gamma, deth
      LOGICAL, INTENT(OUT), OPTIONAL                     :: orthorhombic
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: abc
      INTEGER, DIMENSION(3), INTENT(OUT), OPTIONAL       :: periodic
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h, h_inv
      INTEGER, INTENT(OUT), OPTIONAL                     :: symmetry_id
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: tag

      CPASSERT(ASSOCIATED(cell))

      IF (PRESENT(deth)) deth = cell%deth ! the volume
      IF (PRESENT(orthorhombic)) orthorhombic = cell%orthorhombic
      IF (PRESENT(periodic)) periodic(:) = cell%perd(:)
      IF (PRESENT(h)) h(:, :) = cell%hmat(:, :)
      IF (PRESENT(h_inv)) h_inv(:, :) = cell%h_inv(:, :)

      ! Calculate the lengths of the cell vectors a, b, and c
      IF (PRESENT(abc)) THEN
         abc(1) = SQRT(cell%hmat(1, 1)*cell%hmat(1, 1) + &
                       cell%hmat(2, 1)*cell%hmat(2, 1) + &
                       cell%hmat(3, 1)*cell%hmat(3, 1))
         abc(2) = SQRT(cell%hmat(1, 2)*cell%hmat(1, 2) + &
                       cell%hmat(2, 2)*cell%hmat(2, 2) + &
                       cell%hmat(3, 2)*cell%hmat(3, 2))
         abc(3) = SQRT(cell%hmat(1, 3)*cell%hmat(1, 3) + &
                       cell%hmat(2, 3)*cell%hmat(2, 3) + &
                       cell%hmat(3, 3)*cell%hmat(3, 3))
      END IF

      ! Angles between the cell vectors a, b, and c
      ! alpha = <(b,c)
      IF (PRESENT(alpha)) alpha = angle(cell%hmat(:, 2), cell%hmat(:, 3))*degree
      ! beta = <(a,c)
      IF (PRESENT(beta)) beta = angle(cell%hmat(:, 1), cell%hmat(:, 3))*degree
      ! gamma = <(a,b)
      IF (PRESENT(gamma)) gamma = angle(cell%hmat(:, 1), cell%hmat(:, 2))*degree
      IF (PRESENT(symmetry_id)) symmetry_id = cell%symmetry_id
      IF (PRESENT(tag)) tag = cell%tag

   END SUBROUTINE get_cell

! **************************************************************************************************
!> \brief   Calculate the distance between two lattice planes as defined by
!>          a triple of Miller indices (hkl).
!> \param h ...
!> \param k ...
!> \param l ...
!> \param cell ...
!> \return ...
!> \date    18.11.2004
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION plane_distance(h, k, l, cell) RESULT(distance)

      INTEGER, INTENT(IN)                                :: h, k, l
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp)                                      :: distance

      REAL(KIND=dp)                                      :: a, alpha, b, beta, c, cosa, cosb, cosg, &
                                                            d, gamma, x, y, z
      REAL(KIND=dp), DIMENSION(3)                        :: abc

      x = REAL(h, KIND=dp)
      y = REAL(k, KIND=dp)
      z = REAL(l, KIND=dp)

      CALL get_cell(cell=cell, abc=abc)

      a = abc(1)
      b = abc(2)
      c = abc(3)

      IF (cell%orthorhombic) THEN

         d = (x/a)**2 + (y/b)**2 + (z/c)**2

      ELSE

         CALL get_cell(cell=cell, &
                       alpha=alpha, &
                       beta=beta, &
                       gamma=gamma)

         alpha = alpha/degree
         beta = beta/degree
         gamma = gamma/degree

         cosa = COS(alpha)
         cosb = COS(beta)
         cosg = COS(gamma)

         d = ((x*b*c*SIN(alpha))**2 + &
              (y*c*a*SIN(beta))**2 + &
              (z*a*b*SIN(gamma))**2 + &
              2.0_dp*a*b*c*(x*y*c*(cosa*cosb - cosg) + &
                            z*x*b*(cosg*cosa - cosb) + &
                            y*z*a*(cosb*cosg - cosa)))/ &
             ((a*b*c)**2*(1.0_dp - cosa**2 - cosb**2 - cosg**2 + &
                          2.0_dp*cosa*cosb*cosg))

      END IF

      distance = 1.0_dp/SQRT(d)

   END FUNCTION plane_distance

! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r.
!> \param r ...
!> \param cell ...
!> \return ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc1(r, cell) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      CPASSERT(ASSOCIATED(cell))

      IF (cell%orthorhombic) THEN
         r_pbc(1) = r(1) - cell%hmat(1, 1)*cell%perd(1)*ANINT(cell%h_inv(1, 1)*r(1))
         r_pbc(2) = r(2) - cell%hmat(2, 2)*cell%perd(2)*ANINT(cell%h_inv(2, 2)*r(2))
         r_pbc(3) = r(3) - cell%hmat(3, 3)*cell%perd(3)*ANINT(cell%h_inv(3, 3)*r(3))
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1) + cell%h_inv(1, 2)*r(2) + cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1) + cell%h_inv(2, 2)*r(2) + cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1) + cell%h_inv(3, 2)*r(2) + cell%h_inv(3, 3)*r(3)
         s(1) = s(1) - cell%perd(1)*ANINT(s(1))
         s(2) = s(2) - cell%perd(2)*ANINT(s(2))
         s(3) = s(3) - cell%perd(3)*ANINT(s(3))
         r_pbc(1) = cell%hmat(1, 1)*s(1) + cell%hmat(1, 2)*s(2) + cell%hmat(1, 3)*s(3)
         r_pbc(2) = cell%hmat(2, 1)*s(1) + cell%hmat(2, 2)*s(2) + cell%hmat(2, 3)*s(3)
         r_pbc(3) = cell%hmat(3, 1)*s(1) + cell%hmat(3, 2)*s(2) + cell%hmat(3, 3)*s(3)
      END IF

   END FUNCTION pbc1

! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r subtracting nl from the periodic images
!> \param r ...
!> \param cell ...
!> \param nl ...
!> \return ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc2(r, cell, nl) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nl
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      CPASSERT(ASSOCIATED(cell))

      IF (cell%orthorhombic) THEN
         r_pbc(1) = r(1) - cell%hmat(1, 1)*cell%perd(1)* &
                    REAL(NINT(cell%h_inv(1, 1)*r(1)) - nl(1), dp)
         r_pbc(2) = r(2) - cell%hmat(2, 2)*cell%perd(2)* &
                    REAL(NINT(cell%h_inv(2, 2)*r(2)) - nl(2), dp)
         r_pbc(3) = r(3) - cell%hmat(3, 3)*cell%perd(3)* &
                    REAL(NINT(cell%h_inv(3, 3)*r(3)) - nl(3), dp)
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1) + cell%h_inv(1, 2)*r(2) + cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1) + cell%h_inv(2, 2)*r(2) + cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1) + cell%h_inv(3, 2)*r(2) + cell%h_inv(3, 3)*r(3)
         s(1) = s(1) - cell%perd(1)*REAL(NINT(s(1)) - nl(1), dp)
         s(2) = s(2) - cell%perd(2)*REAL(NINT(s(2)) - nl(2), dp)
         s(3) = s(3) - cell%perd(3)*REAL(NINT(s(3)) - nl(3), dp)
         r_pbc(1) = cell%hmat(1, 1)*s(1) + cell%hmat(1, 2)*s(2) + cell%hmat(1, 3)*s(3)
         r_pbc(2) = cell%hmat(2, 1)*s(1) + cell%hmat(2, 2)*s(2) + cell%hmat(2, 3)*s(3)
         r_pbc(3) = cell%hmat(3, 1)*s(1) + cell%hmat(3, 2)*s(2) + cell%hmat(3, 3)*s(3)
      END IF

   END FUNCTION pbc2

! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by the simulation
!>          cell cell to the vector pointing from atom a to atom b.
!> \param ra ...
!> \param rb ...
!> \param cell ...
!> \return ...
!> \date    11.03.2004
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc3(ra, rb, cell) RESULT(rab_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: rab_pbc

      INTEGER                                            :: icell, jcell, kcell
      INTEGER, DIMENSION(3)                              :: periodic
      REAL(KIND=dp)                                      :: rab2, rab2_pbc
      REAL(KIND=dp), DIMENSION(3)                        :: r, ra_pbc, rab, rb_image, rb_pbc, s2r

      CALL get_cell(cell=cell, periodic=periodic)

      ra_pbc(:) = pbc(ra(:), cell)
      rb_pbc(:) = pbc(rb(:), cell)

      rab2_pbc = HUGE(1.0_dp)

      DO icell = -periodic(1), periodic(1)
         DO jcell = -periodic(2), periodic(2)
            DO kcell = -periodic(3), periodic(3)
               r = REAL((/icell, jcell, kcell/), dp)
               CALL scaled_to_real(s2r, r, cell)
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
! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param cell ...
!> \param positive_range ...
!> \return ...
! **************************************************************************************************
   FUNCTION pbc4(r, cell, positive_range) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: positive_range
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      CPASSERT(ASSOCIATED(cell))

      IF (positive_range) THEN
         IF (cell%orthorhombic) THEN
            r_pbc(1) = r(1) - cell%hmat(1, 1)*cell%perd(1)*FLOOR(cell%h_inv(1, 1)*r(1))
            r_pbc(2) = r(2) - cell%hmat(2, 2)*cell%perd(2)*FLOOR(cell%h_inv(2, 2)*r(2))
            r_pbc(3) = r(3) - cell%hmat(3, 3)*cell%perd(3)*FLOOR(cell%h_inv(3, 3)*r(3))
         ELSE
            s(1) = cell%h_inv(1, 1)*r(1) + cell%h_inv(1, 2)*r(2) + cell%h_inv(1, 3)*r(3)
            s(2) = cell%h_inv(2, 1)*r(1) + cell%h_inv(2, 2)*r(2) + cell%h_inv(2, 3)*r(3)
            s(3) = cell%h_inv(3, 1)*r(1) + cell%h_inv(3, 2)*r(2) + cell%h_inv(3, 3)*r(3)
            s(1) = s(1) - cell%perd(1)*FLOOR(s(1))
            s(2) = s(2) - cell%perd(2)*FLOOR(s(2))
            s(3) = s(3) - cell%perd(3)*FLOOR(s(3))
            r_pbc(1) = cell%hmat(1, 1)*s(1) + cell%hmat(1, 2)*s(2) + cell%hmat(1, 3)*s(3)
            r_pbc(2) = cell%hmat(2, 1)*s(1) + cell%hmat(2, 2)*s(2) + cell%hmat(2, 3)*s(3)
            r_pbc(3) = cell%hmat(3, 1)*s(1) + cell%hmat(3, 2)*s(2) + cell%hmat(3, 3)*s(3)
         END IF
      ELSE
         r_pbc = pbc1(r, cell)
      END IF

   END FUNCTION pbc4

! **************************************************************************************************
!> \brief   Transform real to scaled cell coordinates.
!>          s=h_inv*r
!> \param s ...
!> \param r ...
!> \param cell ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE real_to_scaled(s, r, cell)

      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell

      CPASSERT(ASSOCIATED(cell))

      IF (cell%orthorhombic) THEN
         s(1) = cell%h_inv(1, 1)*r(1)
         s(2) = cell%h_inv(2, 2)*r(2)
         s(3) = cell%h_inv(3, 3)*r(3)
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1) + cell%h_inv(1, 2)*r(2) + cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1) + cell%h_inv(2, 2)*r(2) + cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1) + cell%h_inv(3, 2)*r(2) + cell%h_inv(3, 3)*r(3)
      END IF

   END SUBROUTINE real_to_scaled

! **************************************************************************************************
!> \brief   Transform scaled cell coordinates real coordinates.
!>          r=h*s
!> \param r ...
!> \param s ...
!> \param cell ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE scaled_to_real(r, s, cell)

      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: r
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: s
      TYPE(cell_type), POINTER                           :: cell

      CPASSERT(ASSOCIATED(cell))

      IF (cell%orthorhombic) THEN
         r(1) = cell%hmat(1, 1)*s(1)
         r(2) = cell%hmat(2, 2)*s(2)
         r(3) = cell%hmat(3, 3)*s(3)
      ELSE
         r(1) = cell%hmat(1, 1)*s(1) + cell%hmat(1, 2)*s(2) + cell%hmat(1, 3)*s(3)
         r(2) = cell%hmat(2, 1)*s(1) + cell%hmat(2, 2)*s(2) + cell%hmat(2, 3)*s(3)
         r(3) = cell%hmat(3, 1)*s(1) + cell%hmat(3, 2)*s(2) + cell%hmat(3, 3)*s(3)
      END IF

   END SUBROUTINE scaled_to_real
! **************************************************************************************************
!> \brief retains the given cell (see doc/ReferenceCounting.html)
!> \param cell the cell to retain
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cell_retain(cell)

      TYPE(cell_type), POINTER                           :: cell

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(cell%ref_count > 0)
      cell%ref_count = cell%ref_count + 1

   END SUBROUTINE cell_retain

! **************************************************************************************************
!> \brief releases the given cell (see doc/ReferenceCounting.html)
!> \param cell the cell to release
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cell_release(cell)

      TYPE(cell_type), POINTER                           :: cell

      IF (ASSOCIATED(cell)) THEN
         CPASSERT(cell%ref_count > 0)
         cell%ref_count = cell%ref_count - 1
         IF (cell%ref_count == 0) THEN
            DEALLOCATE (cell)
         END IF
         NULLIFY (cell)
      END IF

   END SUBROUTINE cell_release

#if defined (__PLUMED2)
! **************************************************************************************************
!> \brief   For the interface with plumed, pass a cell pointer and retrieve it
!>          later. It's a hack, but avoids passing the cell back and forth
!>          across the Fortran/C++ interface
!> \param cell ...
!> \param set ...
!> \date    28.02.2013
!> \author  RK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE pbc_cp2k_plumed_getset_cell(cell, set)

      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: set

      TYPE(cell_type), POINTER, SAVE                     :: stored_cell

      IF (set) THEN
         stored_cell => cell
      ELSE
         cell => stored_cell
      END IF

   END SUBROUTINE pbc_cp2k_plumed_getset_cell
#endif

END MODULE cell_types
