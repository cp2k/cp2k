!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: MIT                                                                   !
!--------------------------------------------------------------------------------------------------!

!
!  libgrpp - a library for the evaluation of integrals over
!            generalized relativistic pseudopotentials.
!
!  Copyright (C) 2021-2023 Alexander Oleynichenko
!

MODULE libgrpp
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT32_T

   INTEGER(4), PARAMETER :: LIBGRPP_CART_ORDER_DIRAC = 0
   INTEGER(4), PARAMETER :: LIBGRPP_CART_ORDER_TURBOMOLE = 1

   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE = 0
   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_CHARGED_BALL = 1
   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_GAUSSIAN = 2
   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_FERMI = 3
   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_FERMI_BUBBLE = 4
   INTEGER(4), PARAMETER :: LIBGRPP_NUCLEAR_MODEL_POINT_CHARGE_NUMERICAL = 5

   INTERFACE

      SUBROUTINE libgrpp_init()
         ! no arguments
      END SUBROUTINE libgrpp_init

      SUBROUTINE libgrpp_finalize()
         ! no arguments
      END SUBROUTINE libgrpp_finalize

      SUBROUTINE libgrpp_set_default_parameters()
         ! no arguments
      END SUBROUTINE libgrpp_set_default_parameters

      SUBROUTINE libgrpp_set_radial_tolerance(tolerance)
      REAL(8), INTENT(in)                                :: tolerance

      END SUBROUTINE libgrpp_set_radial_tolerance

      SUBROUTINE libgrpp_set_angular_screening_tolerance(tolerance)
      REAL(8), INTENT(in)                                :: tolerance

      END SUBROUTINE libgrpp_set_angular_screening_tolerance

      SUBROUTINE libgrpp_set_modified_bessel_tolerance(tolerance)
      REAL(8), INTENT(in)                                :: tolerance

      END SUBROUTINE libgrpp_set_modified_bessel_tolerance

      SUBROUTINE libgrpp_set_cartesian_order(order)
      INTEGER(4), INTENT(in)                             :: order

      END SUBROUTINE libgrpp_set_cartesian_order

      SUBROUTINE libgrpp_type1_integrals( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         matrix &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha
      REAL(8), DIMENSION(*), INTENT(out)                 :: matrix

! shell centered on atom B
! pseudopotential expansion
! output: matrix with PP integrals

      END SUBROUTINE libgrpp_type1_integrals

      SUBROUTINE libgrpp_type2_integrals( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_ang_momentum, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         matrix &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), INTENT(in)                             :: rpp_ang_momentum
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha
      REAL(8), DIMENSION(*), INTENT(out)                 :: matrix

! shell centered on atom B
! pseudopotential expansion
! output: matrix with PP integrals

      END SUBROUTINE libgrpp_type2_integrals

      SUBROUTINE libgrpp_spin_orbit_integrals( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_ang_momentum, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         so_x_matrix, so_y_matrix, so_z_matrix &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), INTENT(in)                             :: rpp_ang_momentum
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha
      REAL(8), DIMENSION(*), INTENT(out)                 :: so_x_matrix, so_y_matrix, so_z_matrix

! shell centered on atom B
! pseudopotential expansion
! output: matrices with PP integrals

      END SUBROUTINE libgrpp_spin_orbit_integrals

      SUBROUTINE libgrpp_type1_integrals_gradient( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         point_3d, grad_arep_x, grad_arep_y, grad_arep_z &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha, point_3d
      REAL(8), DIMENSION(*), INTENT(out)                 :: grad_arep_x, grad_arep_y, grad_arep_z

! shell centered on atom B
! pseudopotential expansion
! differentiation wrt the 3d point (x,y,z)
! output: matrices d<Int>/dx, d<Int>/dy, d<Int>/dZ

      END SUBROUTINE libgrpp_type1_integrals_gradient

      SUBROUTINE libgrpp_type2_integrals_gradient( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_ang_momentum, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         point_3d, grad_arep_x, grad_arep_y, grad_arep_z &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), INTENT(in)                             :: rpp_ang_momentum
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha, point_3d
      REAL(8), DIMENSION(*), INTENT(out)                 :: grad_arep_x, grad_arep_y, grad_arep_z

! shell centered on atom B
! pseudopotential expansion
! differentiation wrt the 3d point (x,y,z)
! output: matrices d<Int>/dx, d<Int>/dy, d<Int>/dZ

      END SUBROUTINE libgrpp_type2_integrals_gradient

      SUBROUTINE libgrpp_spin_orbit_integrals_gradient( &
         origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
         origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
         rpp_origin, rpp_ang_momentum, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
         point_3d, grad_sox_x, grad_sox_y, grad_sox_z, &
         grad_soy_x, grad_soy_y, grad_soy_z, &
         grad_soz_x, grad_soz_y, grad_soz_z &
         )
         ! shell centered on atom A
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_A
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: origin_B
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*)
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_origin
      INTEGER(4), INTENT(in)                             :: rpp_ang_momentum
      INTEGER(4), DIMENSION(*), INTENT(in)               :: rpp_num_primitives, rpp_powers
      REAL(8), DIMENSION(*), INTENT(in)                  :: rpp_coeffs, rpp_alpha, point_3d
      REAL(8), DIMENSION(*), INTENT(out)                 :: grad_sox_x, grad_sox_y, grad_sox_z, &
                                                            grad_soy_x, grad_soy_y, grad_soy_z, &
                                                            grad_soz_x, grad_soz_y, grad_soz_z

! shell centered on atom B
! pseudopotential expansion
! differentiation wrt the 3d point (x,y,z)
! output: matrices d<SO_x>/dx, d<SO_x>/dy, d<SO_x>/dZ
! output: matrices d<SO_y>/dx, d<SO_y>/dy, d<SO_y>/dZ
! output: matrices d<SO_z>/dx, d<SO_z>/dy, d<SO_z>/dZ

      END SUBROUTINE libgrpp_spin_orbit_integrals_gradient

   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param origin_A ...
!> \param L_A ...
!> \param num_primitives_A ...
!> \param coeffs_A ...
!> \param alpha_A ...
!> \param origin_B ...
!> \param L_B ...
!> \param num_primitives_B ...
!> \param coeffs_B ...
!> \param alpha_B ...
!> \param rpp_origin ...
!> \param num_oc_shells ...
!> \param oc_shells_L ...
!> \param oc_shells_J ...
!> \param rpp_num_primitives ...
!> \param rpp_powers ...
!> \param rpp_coeffs ...
!> \param rpp_alpha ...
!> \param oc_shells_num_primitives ...
!> \param oc_shells_coeffs ...
!> \param oc_shells_alpha ...
!> \param arep_matrix ...
!> \param so_x_matrix ...
!> \param so_y_matrix ...
!> \param so_z_matrix ...
! **************************************************************************************************
   SUBROUTINE libgrpp_outercore_potential_integrals( &
      origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
      origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
      rpp_origin, num_oc_shells, &
      oc_shells_L, oc_shells_J, rpp_num_primitives, rpp_powers, rpp_coeffs, rpp_alpha, &
      oc_shells_num_primitives, oc_shells_coeffs, oc_shells_alpha, &
      arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
      )

      ! shell centered on atom A
      REAL(8), INTENT(in)                                :: origin_A(*)
      INTEGER(4), INTENT(in)                             :: L_A, num_primitives_A
      REAL(8), INTENT(in)                                :: coeffs_A(*), alpha_A(*), origin_B(*)
      INTEGER(4), INTENT(in)                             :: L_B, num_primitives_B
      REAL(8), INTENT(in)                                :: coeffs_B(*), alpha_B(*), rpp_origin(*)
      INTEGER(4)                                         :: num_oc_shells
      INTEGER(4), INTENT(in)                             :: oc_shells_L(:), oc_shells_J(:), &
                                                            rpp_num_primitives(:), rpp_powers(:, :)
      REAL(8), INTENT(in)                                :: rpp_coeffs(:, :), rpp_alpha(:, :)
      INTEGER(4)                                         :: oc_shells_num_primitives(:)
      REAL(8)                                            :: oc_shells_coeffs(:, :), &
                                                            oc_shells_alpha(:, :)
      REAL(8), INTENT(out)                               :: arep_matrix(*), so_x_matrix(*), &
                                                            so_y_matrix(*), so_z_matrix(*)

      INTEGER                                            :: i, j, ncart1, ncart2
      INTERFACE
         SUBROUTINE libgrpp_outercore_potential_integrals_part_1( &
            origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
            origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
            pot_origin, pot_L, pot_J, pot_num_primitives, pot_powers, pot_coeffs, pot_alpha, &
            oc_shell_num_primitives, oc_shell_coeffs, oc_shell_alpha, &
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix) &
            BIND(C, name="libgrpp_outercore_potential_integrals_part_1_")
            IMPORT :: C_INT32_T, C_DOUBLE
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: origin_A
            INTEGER(kind=C_INT32_T)                  :: L_A
            INTEGER(kind=C_INT32_T)                  :: num_primitives_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: coeffs_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: alpha_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: origin_B
            INTEGER(kind=C_INT32_T)                  :: L_B
            INTEGER(kind=C_INT32_T)                  :: num_primitives_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: coeffs_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: alpha_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot_origin
            INTEGER(kind=C_INT32_T)                  :: pot_L
            INTEGER(kind=C_INT32_T)                  :: pot_J
            INTEGER(kind=C_INT32_T)                  :: pot_num_primitives
            INTEGER(kind=C_INT32_T), DIMENSION(*)    :: pot_powers
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot_alpha
            INTEGER(kind=C_INT32_T)                  :: oc_shell_num_primitives
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_alpha
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: arep_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_x_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_y_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_z_matrix
         END SUBROUTINE libgrpp_outercore_potential_integrals_part_1
      END INTERFACE
      INTERFACE
         SUBROUTINE libgrpp_outercore_potential_integrals_part_2( &
            origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
            origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
            pot_origin, oc_shell_1_L, oc_shell_1_J, &
            pot1_num_primitives, pot1_powers, pot1_coeffs, pot1_alpha, &
            oc_shell_1_num_primitives, oc_shell_1_coeffs, oc_shell_1_alpha, &
            oc_shell_2_L, oc_shell_2_J, pot2_num_primitives, pot2_powers, pot2_coeffs, &
            pot2_alpha, oc_shell_2_num_primitives, oc_shell_2_coeffs, oc_shell_2_alpha, &
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix) &
            BIND(C, name="libgrpp_outercore_potential_integrals_part_2_")
            IMPORT :: C_INT32_T, C_DOUBLE
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: origin_A
            INTEGER(kind=C_INT32_T)                  :: L_A
            INTEGER(kind=C_INT32_T)                  :: num_primitives_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: coeffs_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: alpha_A
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: origin_B
            INTEGER(kind=C_INT32_T)                  :: L_B
            INTEGER(kind=C_INT32_T)                  :: num_primitives_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: coeffs_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: alpha_B
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot_origin
            INTEGER(kind=C_INT32_T)                  :: oc_shell_1_L
            INTEGER(kind=C_INT32_T)                  :: oc_shell_1_J
            INTEGER(kind=C_INT32_T)                  :: pot1_num_primitives
            INTEGER(kind=C_INT32_T), DIMENSION(*)    :: pot1_powers
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot1_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot1_alpha
            INTEGER(kind=C_INT32_T)                  :: oc_shell_1_num_primitives
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_1_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_1_alpha
            INTEGER(kind=C_INT32_T)                  :: oc_shell_2_L
            INTEGER(kind=C_INT32_T)                  :: oc_shell_2_J
            INTEGER(kind=C_INT32_T)                  :: pot2_num_primitives
            INTEGER(kind=C_INT32_T), DIMENSION(*)    :: pot2_powers
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot2_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: pot2_alpha
            INTEGER(kind=C_INT32_T)                  :: oc_shell_2_num_primitives
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_2_coeffs
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: oc_shell_2_alpha
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: arep_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_x_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_y_matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)        :: so_z_matrix
         END SUBROUTINE libgrpp_outercore_potential_integrals_part_2
      END INTERFACE

! shell centered on atom B
! pseudopotential expansion
! outercore shells
! output: matrices with PP integrals
! local variables

      ncart1 = (L_A + 1)*(L_A + 2)/2
      ncart2 = (L_B + 1)*(L_B + 2)/2

      arep_matrix(1:ncart1*ncart2) = 0.0d0
      so_x_matrix(1:ncart1*ncart2) = 0.0d0
      so_y_matrix(1:ncart1*ncart2) = 0.0d0
      so_z_matrix(1:ncart1*ncart2) = 0.0d0

      ! the first non-local term:
      ! \sum_{nlj} U*|nlj><nlj| + |nlj><nlj|*U
      DO i = 1, num_oc_shells
         CALL libgrpp_outercore_potential_integrals_part_1( &
            origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
            origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
            rpp_origin, oc_shells_L(i), oc_shells_J(i), &
            rpp_num_primitives(i), rpp_powers(i, :), rpp_coeffs(i, :), rpp_alpha(i, :), &
            oc_shells_num_primitives(i), oc_shells_coeffs(i, :), oc_shells_alpha(i, :), &
            arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
            )
      END DO

      ! the second non-local term:
      ! \sum_{nlj,n'lj} |nlj><nlj| U |n'lj><n'lj|
      DO i = 1, num_oc_shells
         DO j = 1, num_oc_shells

            CALL libgrpp_outercore_potential_integrals_part_2( &
               origin_A, L_A, num_primitives_A, coeffs_A, alpha_A, &
               origin_B, L_B, num_primitives_B, coeffs_B, alpha_B, &
               rpp_origin, &
               oc_shells_L(i), oc_shells_J(i), &
               rpp_num_primitives(i), rpp_powers(i, :), rpp_coeffs(i, :), rpp_alpha(i, :), &
               oc_shells_num_primitives(i), oc_shells_coeffs(i, :), oc_shells_alpha(i, :), &
               oc_shells_L(j), oc_shells_J(j), &
               rpp_num_primitives(j), rpp_powers(j, :), rpp_coeffs(j, :), rpp_alpha(j, :), &
               oc_shells_num_primitives(j), oc_shells_coeffs(j, :), oc_shells_alpha(j, :), &
               arep_matrix, so_x_matrix, so_y_matrix, so_z_matrix &
               )

         END DO
      END DO

   END SUBROUTINE libgrpp_outercore_potential_integrals
END MODULE libgrpp

