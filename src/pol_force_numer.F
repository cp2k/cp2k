!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_force_numer [1.0] *
!!
!!   NAME
!!     pol_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE pol_force_numer
!------------------------------------------------------------------------------!
!
   USE coefficient_types, ONLY : coeff_type
   USE dg_types, ONLY : dg_type
   USE ewald_parameters_types, ONLY : ewald_parameters_type
   USE kinds, ONLY : dbl
   USE particle_types, ONLY : particle_type
   USE pol_electrostatics, ONLY : electrostatics
   USE pw_grid_types, ONLY : pw_grid_type
   USE simulation_cell, ONLY : cell_type
   USE termination, ONLY : stop_memory, stop_program

   PRIVATE
   PUBLIC :: part_electrostatics_numer, coef_electrostatics_numer
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!
SUBROUTINE part_electrostatics_numer ( delta, dg_part, dg_coef, part, &
                                      coeff, box, pw_small, pw_big,  &
                                      ewald_param, iref, energy_numer, f_part )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_coef
   TYPE ( coeff_type ), INTENT ( IN ) :: coeff
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_small
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_big
   TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( IN ) :: delta

! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER :: id,  natoms
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus

   energy_numer = 0.0_dbl
   f_part = 0.0_dbl
!
! starting the force loop
!
   natoms = SIZE ( part )
   DO i = 1, natoms
     DO id = 1,3
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
        CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2._dbl * delta
        CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2._dbl / delta
  CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                        pw_big, ewald_param, iref, energy_numer )
  RETURN
END SUBROUTINE part_electrostatics_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_electrostatics_numer ( delta, dg_part, dg_coef, part, &
                                      coeff, box, pw_small, pw_big,  &
                                      ewald_param, iref, energy_numer, f_coef )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_coef
   TYPE ( coeff_type ), INTENT ( INOUT ) :: coeff
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_small
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_big
   TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
   INTEGER, INTENT ( IN ) :: iref
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer
   REAL ( dbl ), INTENT ( IN ) :: delta

! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER ::  ncoeff
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus

   energy_numer = 0.0_dbl
   f_coef = 0.0_dbl
!
! starting the force loop
!
   ncoeff = SIZE ( coeff % ao % cr )
   DO i = 1, ncoeff
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
      CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_plus )
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) - 2._dbl * delta
      CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_minus )
      f_coef ( i ) = energy_minus - energy_plus
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
   END DO 
   f_coef = f_coef / 2._dbl / delta
   CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                        pw_big, ewald_param, iref, energy_numer )
   RETURN
END SUBROUTINE coef_electrostatics_numer


!------------------------------------------------------------------------------!
END MODULE pol_force_numer
!------------------------------------------------------------------------------!
