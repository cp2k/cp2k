!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_force_numer_ao [1.0] *
!!
!!   NAME
!!     pol_force_numer_ao
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

MODULE pol_force_numer_ao
!------------------------------------------------------------------------------!
!
  USE atomic_kinds,                    ONLY: kind_info_type
  USE coefficient_types,               ONLY: coeff_type
  USE dg_types,                        ONLY: dg_type
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE kinds,                           ONLY: dp
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE pol_electrostatics_ao,           ONLY: electrostatics
  USE pol_overlap_ao,                  ONLY: force_overlap
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory

   PRIVATE
   PUBLIC :: part_electrostatics_numer, coef_electrostatics_numer, &
             part_overlap_numer, coef_overlap_numer
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!
SUBROUTINE part_electrostatics_numer ( delta, dg_part, dg_coef, part, &
                                      coeff, box, pw_small, pw_big,  &
                                      ewald_param, iref, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part, dg_coef
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(coeff_type), INTENT(IN)             :: coeff
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: pw_small, pw_big
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    INTEGER, INTENT(IN)                      :: iref
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: f_part

    INTEGER                                  :: i, id, natoms
    REAL(KIND=dp)                                :: energy_minus, energy_plus

! begin local variable declaration:
!

   energy_numer = 0.0_dp
   f_part = 0.0_dp
!
! starting the force loop
!
   natoms = SIZE ( part )
   DO i = 1, natoms
     DO id = 1,3
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
        CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                        pw_big, ewald_param, iref, energy_numer )
  RETURN
END SUBROUTINE part_electrostatics_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_electrostatics_numer ( delta, dg_part, dg_coef, part, &
                                      coeff, box, pw_small, pw_big,  &
                                      ewald_param, iref, energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part, dg_coef
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(coeff_type), INTENT(INOUT)          :: coeff
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: pw_small, pw_big
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    INTEGER, INTENT(IN)                      :: iref
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: f_coef

    INTEGER                                  :: i, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus

! begin local variable declaration:
!

   energy_numer = 0.0_dp
   f_coef = 0.0_dp
!
! starting the force loop
!
   ncoeff = SIZE ( coeff % ao % cr )
   DO i = 1, ncoeff
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
      CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_plus )
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) - 2.0_dp * delta
      CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                              pw_big, ewald_param, iref, energy_minus )
      f_coef ( i ) = energy_minus - energy_plus
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL electrostatics ( dg_part, dg_coef, part, coeff, box, pw_small,  &
                        pw_big, ewald_param, iref, energy_numer )
   RETURN
END SUBROUTINE coef_electrostatics_numer

!------------------------------------------------------------------------------!

SUBROUTINE part_overlap_numer ( delta, coeff, drho_basis_info, part, pnode, &
                                box, empparm, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(coeff_type), INTENT(INOUT)          :: coeff
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: drho_basis_info
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: f_part

    INTEGER                                  :: i, id, ios, natoms, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: fc

! begin local variable declaration:
!

   ncoeff = SIZE ( coeff % ao % cr )
   energy_numer = 0.0_dp
   f_part = 0.0_dp

   IF (.NOT. ALLOCATED (fc)) THEN
      ALLOCATE (fc (ncoeff), stat = ios)
      IF (ios /= 0) THEN
         CALL stop_memory("pol_force_numer","fc", ncoeff )      
      END IF
   END IF
!
! starting the force loop
!
   natoms = SIZE ( part )
   DO i = 1, natoms
     DO id = 1,3
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
        CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode,  &
                             box, empparm, energy_plus, fc )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode, &
                             box, empparm, energy_minus, fc )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode, box, &
                        empparm, energy_numer, fc )             

  IF (ALLOCATED (fc)) THEN
     DEALLOCATE (fc, stat = ios)
     IF (ios /= 0) THEN
        CALL stop_memory("pol_force_numer","fc")      
     END IF
  END IF

 RETURN
END SUBROUTINE part_overlap_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_overlap_numer ( delta, coeff, drho_basis_info, part, pnode, &
                                box, empparm, energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(coeff_type), INTENT(INOUT)          :: coeff
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: drho_basis_info
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)     :: f_coef

    INTEGER                                  :: i, ios, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: fc

! begin local variable declaration:
!

   ncoeff = SIZE ( coeff % ao % cr )
   energy_numer = 0.0_dp
   f_coef = 0.0_dp

   IF (.NOT. ALLOCATED (fc)) THEN
      ALLOCATE (fc (ncoeff), stat = ios)
      IF (ios /= 0) THEN
         CALL stop_memory("pol_force_numer","fc", ncoeff )      
      END IF
   END IF
!
! starting the force loop
!
   DO i = 1, ncoeff
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
      CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode,  &
                             box, empparm, energy_plus, fc )
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) - 2.0_dp * delta
      CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode,  &
                             box, empparm, energy_minus, fc )
      f_coef ( i ) = energy_minus - energy_plus
      coeff % ao % cr ( i ) = coeff % ao % cr ( i ) + delta
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL force_overlap ( coeff % ao, drho_basis_info, part, pnode, box, &
                        empparm, energy_numer, fc )             

   IF (ALLOCATED (fc)) THEN
      DEALLOCATE (fc, stat = ios)
      IF (ios /= 0) THEN
         CALL stop_memory("pol_force_numer","fc")      
      END IF
   END IF

   RETURN
END SUBROUTINE coef_overlap_numer

!------------------------------------------------------------------------------!
END MODULE pol_force_numer_ao
!------------------------------------------------------------------------------!
