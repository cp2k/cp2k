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
  USE ao_types,                        ONLY: ao_type
  USE atomic_kinds,                    ONLY: kind_info_type
  USE coefficient_types,               ONLY: PW_COMPLEXDATA1D,&
                                             coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_type,&
                                             coeff_zero
  USE dg_types,                        ONLY: dg_type
  USE empirical_parameters,            ONLY: empirical_parameter_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE kinds,                           ONLY: dp
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE pol_electrostatics_pw,           ONLY: electrostatics
  USE pol_fo_kernel_pw,                ONLY: force_fo_kernel
  USE pol_overlap_pw,                  ONLY: force_overlap
  USE pol_projectors_pw,               ONLY: projectors
  USE pol_thk_pw,                      ONLY: force_thk
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: RECIPROCALSPACE,&
                                             pw_type
  USE simulation_cell,                 ONLY: cell_type

   PRIVATE
   PUBLIC :: part_electrostatics_numer, coef_electrostatics_numer, &
             part_overlap_numer, coef_overlap_numer, part_thk_numer, &
             coef_thk_numer, coef_fo_kernel_numer, part_fo_kernel_numer
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!
SUBROUTINE part_electrostatics_numer ( delta, dg_part, part, pw, box, pw_small, &
                                       pw_big, ewald_param, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: pw_small, pw_big
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: f_part

    INTEGER                                  :: i, id, natoms
    REAL(KIND=dp)                                :: energy_minus, energy_plus
    TYPE(coeff_type)                         :: rho0

! begin local variable declaration:
!

   CALL coeff_allocate ( rho0 , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
   rho0 % pw % in_space = RECIPROCALSPACE
   CALL coeff_zero ( rho0 )

   energy_numer = 0.0_dp
   f_part = 0.0_dp
!
! starting the force loop
!
   natoms = SIZE ( part )
   DO i = 1, natoms
     DO id = 1,3
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
        CALL coeff_zero ( rho0 )
        CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                              pw_big, ewald_param, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL coeff_zero ( rho0 )
        CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                              pw_big, ewald_param,  energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL coeff_zero ( rho0 )
  CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                        pw_big, ewald_param, energy_numer )
  CALL coeff_deallocate ( rho0 )
  RETURN
END SUBROUTINE part_electrostatics_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_electrostatics_numer ( delta, dg_part, part, &
                                      pw, box, pw_small, pw_big,  &
                                      ewald_param,  energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: pw_small, pw_big
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)  :: f_coef

    INTEGER                                  :: i, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus
    TYPE(coeff_type)                         :: rho0

! begin local variable declaration:
!

   CALL coeff_allocate ( rho0 , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
   rho0 % pw % in_space = RECIPROCALSPACE
   CALL coeff_zero ( rho0 )

   energy_numer = 0.0_dp
   f_coef = 0.0_dp
!
! starting the force loop
!
   ncoeff = SIZE ( pw % cc )
!dbg
   ncoeff = ncoeff/100
!dbg
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ),KIND=dp)
      CALL coeff_zero ( rho0 )
      CALL electrostatics ( dg_part, pw, rho0, part,  box, pw_small,  &
                            pw_big, ewald_param, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) - 2.0_dp * delta, &
                              AIMAG ( pw % cc ( i ) ),KIND=dp)
      CALL coeff_zero ( rho0 )
      CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                              pw_big, ewald_param, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0.0_dp,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ),KIND=dp)
!bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
      CALL coeff_zero ( rho0 )
      CALL electrostatics ( dg_part, pw, rho0, part,  box, pw_small,  &
                            pw_big, ewald_param, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) - 2.0_dp * delta,KIND=dp)
      CALL coeff_zero ( rho0 )
      CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                              pw_big, ewald_param, energy_minus )
      f_coef ( i ) =  CMPLX ( REAL ( f_coef ( i ),KIND=dp), energy_minus - energy_plus,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL coeff_zero ( rho0 )
   CALL electrostatics ( dg_part, pw, rho0, part, box, pw_small,  &
                        pw_big, ewald_param,  energy_numer )
   CALL coeff_deallocate ( rho0 )
   RETURN
END SUBROUTINE coef_electrostatics_numer

!------------------------------------------------------------------------------!

SUBROUTINE part_overlap_numer ( delta, dg, pw, ao, drho_basis_info, part, pnode, &
                                box, empparm, grid_s, grid_b, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: drho_basis_info
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
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
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                             empparm, grid_s, grid_b, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                             empparm, grid_s, grid_b, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
  CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                       empparm, grid_s, grid_b, energy_numer )
 RETURN
END SUBROUTINE part_overlap_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_overlap_numer ( delta, dg, pw, ao, drho_basis_info, part, pnode, &
                                box, empparm, grid_s, grid_b, energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: drho_basis_info
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)  :: f_coef

    INTEGER                                  :: i, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus

! begin local variable declaration:
!

   ncoeff = SIZE ( pw % cc )
   energy_numer = 0.0_dp
   f_coef = 0.0_dp
!dbg
   ncoeff = ncoeff/100 
!dbg

!
! starting the force loop
!
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) - 2.0_dp * delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0.0_dp,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )

! bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) - 2.0_dp * delta,KIND=dp)
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                           empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( REAL ( f_coef ( i ),KIND=dp), energy_minus - energy_plus,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
   CALL force_overlap ( dg, pw % pw_grid, ao, drho_basis_info, part, pnode, box, &
                        empparm, grid_s, grid_b, energy_numer )
   RETURN
END SUBROUTINE coef_overlap_numer

!------------------------------------------------------------------------------!

SUBROUTINE part_thk_numer ( delta, dg, pw, ao, part, box, empparm, grid_s,  &
                            grid_b, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
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
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
        CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
  CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_numer )
 RETURN
END SUBROUTINE part_thk_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_thk_numer ( delta, dg, pw, ao, part, box, empparm, grid_s,  &
                            grid_b, energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(empirical_parameter_type), &
      DIMENSION(:), INTENT(IN)               :: empparm
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)  :: f_coef

    INTEGER                                  :: i, ncoeff
    REAL(KIND=dp)                                :: energy_minus, energy_plus

! begin local variable declaration:
!

   ncoeff = SIZE ( pw % cc )
   energy_numer = 0.0_dp
   f_coef = 0.0_dp
!dbg
   ncoeff = ncoeff/100 
!dbg

!
! starting the force loop
!
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) - 2.0_dp * delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0.0_dp,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )

! bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) - 2.0_dp * delta,KIND=dp)
      CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
      CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( REAL ( f_coef ( i ),KIND=dp), energy_minus - energy_plus,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL projectors ( dg, ao, pw, part, box, grid_s, grid_b )
   CALL force_thk ( dg, ao, pw % pw_grid, part, box, empparm, grid_s, grid_b, energy_numer )
   RETURN
END SUBROUTINE coef_thk_numer
!------------------------------------------------------------------------------!

SUBROUTINE part_fo_kernel_numer ( delta, dg_part, pw, part, box, grid_s,  &
                            grid_b, ewald_param, energy_numer, f_part )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg_part
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: f_part

    INTEGER                                  :: i, id, natoms
    REAL(KIND=dp)                                :: ener, energy_minus, &
                                                energy_plus
    TYPE(coeff_type)                         :: rho0

! begin local variable declaration:
!

   CALL coeff_allocate ( rho0 , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
   rho0 % pw % in_space = RECIPROCALSPACE
   CALL coeff_zero ( rho0 )


   energy_numer = 0.0_dp
   f_part = 0.0_dp

!
! starting the force loop
!
   natoms = SIZE ( part )
   DO i = 1, natoms
     DO id = 1,3
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
        CALL coeff_zero ( rho0 )
        CALL electrostatics ( dg_part, pw, rho0, part,  box, grid_s,  &
                            grid_b, ewald_param, ener )
        CALL  force_fo_kernel ( dg_part, pw, rho0, part, box, grid_s,  &
                              grid_b, energy_plus )
        part ( i ) % r ( id ) = part ( i ) % r ( id ) - 2.0_dp * delta
        CALL coeff_zero ( rho0 )
        CALL electrostatics ( dg_part, pw, rho0, part,  box, grid_s,  &
                            grid_b, ewald_param, ener )
        CALL force_fo_kernel ( dg_part, pw, rho0, part, box, grid_s,  &
                              grid_b, energy_minus )
        f_part ( id, i ) = energy_minus - energy_plus
        part ( i ) % r ( id ) = part ( i ) % r ( id ) + delta
     END DO
  END DO 
  f_part = f_part / 2.0_dp / delta
  CALL coeff_zero ( rho0 )
  CALL electrostatics ( dg_part, pw, rho0, part, box, grid_s, grid_b, ewald_param, ener )
  CALL force_fo_kernel ( dg_part, pw, rho0, part, box, grid_s, grid_b, energy_numer )
  CALL coeff_deallocate (rho0 )
 RETURN
END SUBROUTINE part_fo_kernel_numer

!------------------------------------------------------------------------------!
SUBROUTINE coef_fo_kernel_numer ( delta, dg, pw, part, box, grid_s,  &
                            grid_b, ewald_param, energy_numer, f_coef )
!
! begin passed variable declaration:
!
    REAL(KIND=dp), INTENT(IN)                    :: delta
    TYPE(dg_type), DIMENSION(:), INTENT(IN)  :: dg
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_s, grid_b
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), INTENT(OUT)                   :: energy_numer
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)  :: f_coef

    INTEGER                                  :: i, ncoeff
    REAL(KIND=dp)                                :: ener, energy_minus, &
                                                energy_plus
    TYPE(coeff_type)                         :: rho0

! begin local variable declaration:
!

   CALL coeff_allocate ( rho0 , pw % pw_grid, use_data = PW_COMPLEXDATA1D )
   rho0 % pw % in_space = RECIPROCALSPACE
   CALL coeff_zero ( rho0 )
   CALL electrostatics ( dg, pw, rho0, part, box, grid_s, grid_b, ewald_param, ener )

   ncoeff = SIZE ( pw % cc )
   energy_numer = 0.0_dp
   f_coef = 0.0_dp
!dbg
   ncoeff = ncoeff/100 
!dbg

!
! starting the force loop
!
   DO i = 1, ncoeff
! ags
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL force_fo_kernel ( dg, pw, rho0, part, box, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) - 2.0_dp * delta, &
                              AIMAG ( pw % cc ( i ) ) )
      CALL force_fo_kernel ( dg, pw, rho0, part, box, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( energy_minus - energy_plus, 0.0_dp,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp) + delta, &
                              AIMAG ( pw % cc ( i ) ) )

! bgs
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
      CALL force_fo_kernel ( dg, pw, rho0, part, box, grid_s, grid_b, energy_plus )
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) - 2.0_dp * delta,KIND=dp)
      CALL force_fo_kernel ( dg, pw, rho0, part, box, grid_s, grid_b, energy_minus )
      f_coef ( i ) = CMPLX ( REAL ( f_coef ( i ),KIND=dp), energy_minus - energy_plus,KIND=dp)
      pw % cc ( i ) = CMPLX ( REAL ( pw % cc ( i ),KIND=dp), &
                              AIMAG ( pw % cc ( i ) ) + delta,KIND=dp)
   END DO 
   f_coef = f_coef / 2.0_dp / delta
   CALL force_fo_kernel ( dg, pw, rho0, part, box, grid_s, grid_b, energy_numer )
   CALL coeff_deallocate (rho0)
   RETURN
END SUBROUTINE coef_fo_kernel_numer


!------------------------------------------------------------------------------!
END MODULE pol_force_numer_pw
!------------------------------------------------------------------------------!
