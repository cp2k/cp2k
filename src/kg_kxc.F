!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_force [1.0] *
!!
!!   NAME
!!     kg_kxc 
!!
!!   FUNCTION
!!     calculation of the Kinetic + Exchange correlation
!!
!!   AUTHOR
!!     gloria,30.09.2002 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_kxc

  USE coefficient_types, ONLY: coeff_type, coeff_transform_space
  USE dft_types, ONLY : dft_control_type
  USE functionals
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE message_passing, ONLY : mp_sum, mp_max, mp_sync,&
                  mp_range, mp_allgather
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop
  USE xc,      ONLY: xc_calculate

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_kxc_potential

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kgkxc/calculate_kxc_potential [1.1] *
!!
!!   NAME
!!     calculate_kxc_potential
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE calculate_kxc_potential ( rho_r, rho_g, v_r, e_xc, dft_control )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT (in) :: rho_r
  TYPE ( coeff_type ), INTENT (in) :: rho_g
  TYPE ( coeff_type ), INTENT (out) :: v_r
  REAL ( dbl ), INTENT (out) :: e_xc
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_kxc_potential (MODULE kg_kxc)"
! Locals
  REAL (dbl), DIMENSION(:,:,:), POINTER :: rho
  REAL (dbl), DIMENSION(:,:,:), POINTER :: vxc
  COMPLEX (dbl), DIMENSION(:), POINTER :: rhog
!--------------------------------------------------------

  e_xc = 0.0_dbl
  rho => rho_r%pw%cr3d(:,:,:)
  vxc => v_r%pw%cr3d(:,:,:)
  rhog => rho_g%pw%cc(:)
  rho (:,:,:) = -1.0_dbl * rho(:,:,:)
  rhog (:) = -1.0_dbl * rhog(:)
  call xc_calculate (rho_r%pw,rho_g%pw,v_r%pw, &
                            e_xc,dft_control)
 
  e_xc =  rho_r%pw%pw_grid%dvol*e_xc      
  rho (:,:,:) = -1.0_dbl * rho(:,:,:)
  rhog (:) = -1.0_dbl * rhog(:)
  vxc (:,:,:) = -1.0_dbl * vxc(:,:,:) *rho_r%pw%pw_grid%dvol

 END SUBROUTINE calculate_kxc_potential

!******************************************************************************

END MODULE kg_kxc  

!******************************************************************************
