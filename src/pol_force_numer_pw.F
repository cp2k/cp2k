!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_force_numer_pw [1.0] *
!!
!!   NAME
!!     pol_force_numer_pw
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

MODULE pol_force_numer_pw
!------------------------------------------------------------------------------!
!
   USE ao_types, ONLY : ao_type
   USE atomic_kinds, ONLY: kind_info_type
   USE coefficient_types, ONLY : coeff_type
   USE dg_types, ONLY : dg_type
   USE ewald_parameters_types, ONLY : ewald_parameters_type
   USE empirical_parameters, ONLY : empirical_parameter_type
   USE kinds, ONLY : dbl
   USE molecule_types, ONLY : particle_node_type
   USE particle_types, ONLY : particle_type
   USE pol_electrostatics_pw, ONLY : electrostatics
   USE pol_overlap_pw, ONLY : force_overlap
   USE pol_projectors_pw, ONLY : projectors
   USE pol_thk_pw, ONLY : force_thk
   USE pw_grid_types, ONLY : pw_grid_type
   USE pw_types, ONLY : pw_type
   USE simulation_cell, ONLY : cell_type
   USE termination, ONLY : stop_memory, stop_program

   PRIVATE
   PUBLIC :: part_electrostatics_numer, coef_electrostatics_numer, &
             part_overlap_numer, coef_overlap_numer, part_thk_numer, &
             coef_thk_numer
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!
SUBROUTINE part_electrostatics_numer ( delta, dg_part, part, pw, box, pw_small, &
                                       pw_big, ewald_param, energy_numer, f_part )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_small
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_big
   TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
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
        CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                              pw_big, ewald_param, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2._dbl * delta
        CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                              pw_big, ewald_param,  energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2._dbl / delta
  CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                        pw_big, ewald_param, energy_numer )
  RETURN
END SUBROUTINE part_electrostatics_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_electrostatics_numer ( delta, dg_part, part, &
                                      pw, box, pw_small, pw_big,  &
                                      ewald_param,  energy_numer, f_coef )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg_part
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_small
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: pw_big
   TYPE ( ewald_parameters_type ), INTENT ( IN ) :: ewald_param
   COMPLEX ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef
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
   ncoeff = SIZE ( pw % cc )
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ), dbl )
      CALL electrostatics ( dg_part, pw, part,  box, pw_small,  &
                            pw_big, ewald_param, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) - 2._dbl * delta, &
                              AIMAG ( pw % cc ( i ) ), dbl )
      CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                              pw_big, ewald_param, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0._dbl, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ), dbl )
!bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + delta, dbl )
      CALL electrostatics ( dg_part, pw, part,  box, pw_small,  &
                            pw_big, ewald_param, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) - 2._dbl * delta, dbl )
      CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                              pw_big, ewald_param, energy_minus )
      f_coef ( i ) =  CMPLX ( REAL ( f_coef ( i ), dbl ), energy_minus - energy_plus, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + delta, dbl )
   END DO 
   f_coef = f_coef / 2._dbl / delta
   CALL electrostatics ( dg_part, pw, part, box, pw_small,  &
                        pw_big, ewald_param,  energy_numer )
   RETURN
END SUBROUTINE coef_electrostatics_numer

!------------------------------------------------------------------------------!

SUBROUTINE part_overlap_numer ( delta, dg, pw, ao, drho_basis_info, part, pnode, &
                                box, empparm, grid_s, grid_b, energy_numer, f_part )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   REAL ( dbl ), INTENT ( IN ) :: delta
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg
   TYPE ( pw_type ), INTENT ( IN ) :: pw
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s, grid_b
   TYPE ( ao_type ), INTENT ( INOUT ) :: ao
   TYPE ( kind_info_type ),  INTENT ( IN ), DIMENSION ( : ) :: drho_basis_info
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer

! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER :: ios
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
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                             empparm, grid_s, grid_b, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2._dbl * delta
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                             empparm, grid_s, grid_b, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2._dbl / delta
  CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
  CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                       empparm, grid_s, grid_b, energy_numer )
 RETURN
END SUBROUTINE part_overlap_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_overlap_numer ( delta, dg, pw, ao, drho_basis_info, part, pnode, &
                                box, empparm, grid_s, grid_b, energy_numer, f_coef )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   REAL ( dbl ), INTENT ( IN ) :: delta
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : )  :: grid_s, grid_b
   TYPE ( ao_type ), INTENT ( INOUT ) :: ao
   TYPE ( kind_info_type ),  INTENT (IN), DIMENSION (:) :: drho_basis_info
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
   COMPLEX ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer


! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER :: ios
   INTEGER ::  ncoeff
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus

   ncoeff = SIZE ( pw % cc )
   energy_numer = 0.0_dbl
   f_coef = 0.0_dbl

!
! starting the force loop
!
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) - 2._dbl * delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0._dbl, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ) )

! bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + delta, dbl )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) - 2._dbl * delta, dbl )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( REAL ( f_coef ( i ), dbl ), energy_minus - energy_plus, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + 2._dbl * delta, dbl )
   END DO 
   f_coef = f_coef / 2._dbl / delta
   CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
   CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                        empparm, grid_s, grid_b, energy_numer )
   RETURN
END SUBROUTINE coef_overlap_numer

!------------------------------------------------------------------------------!

SUBROUTINE part_thk_numer ( delta, dg, pw, ao, part, box, empparm, grid_s,  &
                            grid_b, energy_numer, f_part )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   REAL ( dbl ), INTENT ( IN ) :: delta
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg
   TYPE ( pw_type ), INTENT ( IN ) :: pw
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s, grid_b
   TYPE ( ao_type ), INTENT ( INOUT ) :: ao
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
   REAL ( dbl ), INTENT ( OUT ), DIMENSION ( :, : ) :: f_part
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer

! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER :: ios
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
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2._dbl * delta
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2._dbl / delta
  CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
  CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_numer )
 RETURN
END SUBROUTINE part_thk_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_thk_numer ( delta, dg, pw, ao, part, box, empparm, grid_s,  &
                            grid_b, energy_numer, f_coef )
   IMPLICIT NONE
!
! begin passed variable declaration:
!
   REAL ( dbl ), INTENT ( IN ) :: delta
   TYPE ( dg_type ), INTENT ( IN ), DIMENSION ( : ) :: dg
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw
   TYPE ( pw_grid_type ), INTENT ( IN ), DIMENSION ( : ) :: grid_s, grid_b
   TYPE ( ao_type ), INTENT ( INOUT ) :: ao
   TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
   TYPE ( cell_type ), INTENT ( IN ) :: box
   TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
   COMPLEX ( dbl ), INTENT ( OUT ), DIMENSION ( : ) :: f_coef
   REAL ( dbl ), INTENT ( OUT ) :: energy_numer


! begin local variable declaration:
!
   INTEGER :: i, j 
   INTEGER :: ios
   INTEGER ::  ncoeff
   REAL (dbl) :: energy
   REAL (dbl) :: energy_plus, energy_minus

   ncoeff = SIZE ( pw % cc )
   energy_numer = 0.0_dbl
   f_coef = 0.0_dbl

!
! starting the force loop
!
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) - 2._dbl * delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0._dbl, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ) + delta, &
                              AIMAG ( pw % cc ( i ) ) )

! bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + delta, dbl )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) - 2._dbl * delta, dbl )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( REAL ( f_coef ( i ), dbl ), energy_minus - energy_plus, dbl )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ), dbl ), &
                              AIMAG ( pw % cc ( i ) ) + 2._dbl * delta, dbl )
   END DO 
   f_coef = f_coef / 2._dbl / delta
   CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
   CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_numer )
   RETURN
END SUBROUTINE coef_thk_numer

!------------------------------------------------------------------------------!
END MODULE pol_force_numer_pw
!------------------------------------------------------------------------------!
