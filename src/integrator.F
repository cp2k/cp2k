!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/integrator [1.0] *
!!
!!   NAME
!!     integrator
!!
!!   FUNCTION
!!     Provides integrator routines (velocity verlet) for all the
!!     ensemble types
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Mar-2001) : Pass logical for box change to force routine
!!     Harald Forbert (Apr-2001): added path integral routine nvt_pimd
!!     CJM (15-Apr-2001) : added coef integrators and energy routines
!!
!!   SOURCE
!******************************************************************************

MODULE integrator
  
  USE constraint, ONLY: shake_control, rattle_control, pv_constraint, getold, &
      shake_roll_control, rattle_roll_control
  USE eigenvalueproblems, ONLY : diagonalise
  USE force_control, ONLY : force
  USE global_types, ONLY: global_environment_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE md, ONLY : simulation_parameters_type, virial_type, &
                 mdio_parameters_type
  USE md_environment_types, ONLY : md_environment_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
       linklist_atoms
  USE message_passing, ONLY : mp_sum, mp_sync, mp_comm_dup
  USE nose, ONLY : lnhc_parameters_type, npt_info_type, lnhc, yoshida_coef
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  USE structure_types, ONLY : structure_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop
  USE transformations, ONLY : u2r, fr2fu
  USE util, ONLY : get_unit, matvec_3x3, matmul_3x3, transpose_3d

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_integrator, nve, nvt, nvt_pimd, npt_i, npt_f, &
            nve_ao, nvt_ao, nve_pw, nvt_pw, pressure

  TYPE old_variables_type
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: v
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: r
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: eps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: veps
     REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: h
  END TYPE old_variables_type

  TYPE ( old_variables_type ) :: old

  LOGICAL :: ionode
  INTEGER :: int_group

  TYPE ( global_environment_type ) :: glob_env

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** integrator/set_integrator [1.0] *
!!
!!   NAME
!!     set_integrator
!!
!!   FUNCTION
!!     setup of parameters for integrator
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (29-Mar-2001) : duplicate communicator, instead of copy
!!
!!   NOTES
!!     Call this subroutine before the first call to energy or velocity_verlet
!!     or if you want to change ionode and/or output files
!!
!!   SOURCE
!******************************************************************************


SUBROUTINE set_integrator ( globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT (IN) :: globenv

!------------------------------------------------------------------------------

  glob_env = globenv
  CALL mp_comm_dup ( globenv % group, glob_env % group )

  ionode = glob_env % ionode
  int_group = glob_env % group

END SUBROUTINE set_integrator

!******************************************************************************
!!****** integrator/nve_pw [1.0] *
!!
!!   NAME
!!     nve_pw
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta and AO coefficient type
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve_pw ( md_env )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Local
  LOGICAL, PARAMETER :: box_change = .FALSE.
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm, cmass
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc

!------------------------------------------------------------------------------
! Associating local pointers
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos ( 1 ) % pw % cc )
  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
                     + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
     struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
             + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
  END DO
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve_pw

!!*****
!******************************************************************************
!!****** integrator/nve_ao [1.0] *
!!
!!   NAME
!!     nve_ao
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve_ao ( md_env )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Local
  LOGICAL, PARAMETER :: box_change = .FALSE.
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef

!------------------------------------------------------------------------------
! Associating local structures
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc_coef => md_env % nhc_coef ( 1 )
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos ( 1 ) % ao % cr )
  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
                     + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
             + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force (  glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
  END DO
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve_ao
!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_pw
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_pw ( md_env )

  IMPLICIT NONE
  
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env

! Local
  LOGICAL, PARAMETER :: box_change = .FALSE.
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom,dtocm,cmass
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef

!------------------------------------------------------------------------------
! Associating local pointers
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc_part => md_env % nhc_part ( 1 )
  nhc_coef => md_env % nhc_coef ( 1 )
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos (1) % pw % cc )
  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
     struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
                      + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
     dtocm = 0.5_dbl * simpar % dt / cmass
     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
  END DO
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt_pw

!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_ao
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_ao ( md_env )

  IMPLICIT NONE
  
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Local
  LOGICAL, PARAMETER :: box_change = .FALSE.
  INTEGER :: i, nnodes, ncoef
  REAL ( dbl ) :: dtom, dtocm
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef

!------------------------------------------------------------------------------
! Associating local pointers
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc_part => md_env % nhc_part ( 1 )
  nhc_coef => md_env % nhc_coef ( 1 )
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

  ncoef = SIZE ( struc % coef_pos (1) % ao % cr )
  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
                      + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
  END DO

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  DO i = 1, ncoef
     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
  END DO
!
!  IF ( simpar % pol_ensemble == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )

  CALL lnhc ( nhc_part, struc % pnode, int_group )

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt_ao

!!*****
!******************************************************************************
!!****** integrator/nve [1.0] *
!!
!!   NAME
!!     nve
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nve ( md_env )
  IMPLICIT NONE
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env

! locals
  LOGICAL, PARAMETER :: box_change = .FALSE.
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc

! assigning the locals
  simpar => md_env % simpar
  struc => md_env % struc ( 1 )
  thermo => md_env % virial


!------------------------------------------------------------------------------

  IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
  nnodes = size(struc%pnode)
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nve

!!*****
!******************************************************************************
!!****** integrator/nvt [1.0] *
!!
!!   NAME
!!     nvt
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt ( md_env )
  
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Locals
  INTEGER :: i, nnodes
  REAL ( dbl ) :: dtom
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc
  TYPE ( lnhc_parameters_type ), POINTER :: nhc
  LOGICAL, PARAMETER :: box_change = .FALSE.

! assigning the locals
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc => md_env % nhc_part ( 1 )

!------------------------------------------------------------------------------

  CALL lnhc ( nhc, struc % pnode, int_group )

  IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
  nnodes = SIZE ( struc % pnode )
  DO i = 1, nnodes
     dtom = 0.5_dbl * simpar%dt / struc%pnode(i)%p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
          struc%pnode(i) %p%v(1)*simpar%dt
     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
          struc%pnode(i) %p%v(2)*simpar%dt
     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
          struc%pnode(i) %p%v(3)*simpar%dt
  END DO

  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
#if defined(__parallel)
  CALL update_structure ( struc, 'POS' )
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
          dtom*struc%pnode(i) %p%f(1)
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
          dtom*struc%pnode(i) %p%f(2)
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
          dtom*struc%pnode(i) %p%f(3)
  END DO

  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )

  CALL lnhc ( nhc, struc % pnode, int_group )

#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! updating the virial
  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
                                            'FULL',int_group )
  CALL pressure(struc%pnode,thermo)

END SUBROUTINE nvt
!!*****
!******************************************************************************
!!****** integrator/npt_i [1.0] *
!!
!!   NAME
!!     npt_i
!!
!!   FUNCTION
!!     npt_i integrator for particle positions & momenta
!!     isotropic box changes
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE npt_i ( md_env )
  
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Locals
  INTEGER :: i, ii,  nnodes, iroll
  TYPE (cell_type ), POINTER :: box
  TYPE (simulation_parameters_type ), POINTER :: simpar
  TYPE (virial_type ), POINTER :: thermo
  TYPE (structure_type ), POINTER :: struc
  TYPE (lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE (npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: arg_r, poly_r, scale_r, roll_tol
  REAL ( dbl ) :: arg_v, poly_v, scale_v, kin, fdotr, infree
  REAL ( dbl ), DIMENSION (3) :: vector_r, vector_v
  LOGICAL, PARAMETER :: box_change = .TRUE.
  LOGICAL :: first

! assigning the locals
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc_part => md_env % nhc_part ( 1 )
  nhc_baro => md_env % nhc_baro ( 1 )
  npt => md_env % npt
  box => md_env % cell

!------------------------------------------------------------------------------
  infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 



! initialize eps_0 the first time through
  IF ( md_env % itimes == 1 ) eps_0 = npt (1,1) % eps


! Nose-Hoover for the Barostat
     CALL lnhc ( nhc_baro, npt, int_group )
! Nose-Hoover for the Particles
     CALL lnhc ( nhc_part, struc%pnode, int_group )

! setting up for ROLL: saving old variables
  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     iroll = 1
     CALL set ( struc, box, npt, 'F' )
     CALL getold ( struc % molecule )
  END IF

  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP

     IF (simpar%constraint) THEN
        CALL set( struc, box, npt, 'B' )
     END IF
     CALL get_veps ( box, npt, simpar, struc, thermo, infree )

!
! first half of velocity verlet
!
     arg_r = ( 0.5_dbl * npt ( 1, 1 ) % v * simpar % dt ) * &
             ( 0.5_dbl * npt ( 1, 1 ) % v * simpar % dt)
     poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

     arg_v = ( 0.25_dbl * npt ( 1, 1 ) % v * simpar % dt *  &
             (1._dbl + 3._dbl * infree ) )*( 0.25_dbl * npt ( 1, 1 ) % v * &
              simpar % dt * ( 1._dbl + 3._dbl * infree ) )
     poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

     scale_r = EXP( 0.5_dbl * simpar % dt * npt ( 1, 1 ) %v )
     scale_v = EXP( -0.25_dbl * simpar % dt * npt ( 1, 1 ) % v * &
                  (1._dbl + 3._dbl * infree ) )
     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass

        struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(1)* dtom
        struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(2)* dtom
        struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3)*scale_v*scale_v + &
             scale_v*poly_v*struc%pnode(i) %p%f(3)* dtom

        struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(1)*simpar%dt
        struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(2)*simpar%dt
        struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3)*scale_r*scale_r + &
             scale_r*poly_r*struc%pnode(i) %p%v(3)*simpar%dt
     END DO
     roll_tol = 0._dbl
     vector_r ( : ) = scale_r * poly_r
     vector_v ( : ) = scale_v * poly_v

     IF ( simpar%constraint ) CALL shake_roll_control(struc%molecule,simpar, &
          thermo,roll_tol,iroll,vector_r,vector_v,int_group )
  END DO SR

! Update eps:
  npt ( :, : ) % eps =  npt ( :, : ) % eps &
                       + simpar % dt * npt ( :, : ) % v

! Update h_mat
  box % hmat ( :, : ) = box % hmat ( :, : ) &
       * EXP ( npt ( 1, 1 ) % eps - eps_0 )

  eps_0 =  npt ( 1, 1 ) % eps

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(1)* dtom
     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(2)* dtom
     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3)*scale_v*scale_v + &
          scale_v*poly_v*struc%pnode(i) %p%f(3)* dtom
  END DO


  roll_tol = -1.E10_dbl
  IF (simpar%constraint) THEN
     first = .TRUE.
     iroll = 1
     CALL set( struc, box, npt, 'F' )
  END IF

  RR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! RATTLE-ROLL LOOP
     roll_tol = 0._dbl
     IF (simpar%constraint) CALL rattle_roll_setup ( box, npt, struc, &
            thermo, simpar, vector_v, roll_tol, iroll, infree, first )

     CALL get_veps ( box, npt, simpar, struc, thermo, infree )


  END DO RR


! Nose-Hoover for the Particles
  CALL lnhc ( nhc_part, struc%pnode, int_group )
! Nose-Hoover for the Barostat
  CALL lnhc ( nhc_baro, npt, int_group )


#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE npt_i

!!*****
!******************************************************************************
!!****** integrator/pressure [1.0] *
!!
!!   NAME
!!     pressure
!!
!!   FUNCTION
!!     Computes the kinetic part of the pressure tensor and updates
!!     the full VIRIAL (PV)
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pressure ( pnode, thermo )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo

! Locals
  INTEGER :: i, j, iatom, nnodes, handle
  REAL ( dbl ) :: mfl

!------------------------------------------------------------------------------

  CALL timeset ( 'PRESSURE', 'E', 'Mflops', handle )

  thermo % pv_kinetic = 0._dbl
  nnodes = size(pnode)
  DO i = 1, 3
     DO j = 1, 3
        DO iatom = 1, nnodes
           thermo % pv_kinetic ( i, j ) = thermo % pv_kinetic ( i, j ) + &
                pnode ( iatom ) % p % prop % mass * &
                pnode ( iatom ) % p % v ( i )* &
                pnode ( iatom ) % p % v ( j )
        END DO
     END DO
  END DO
  mfl = REAL( 9 * nnodes, dbl ) * 2._dbl * 1.e-6_dbl

#if defined(__parallel)
  CALL mp_sum(thermo%pv_kinetic,int_group)
#endif

! total virial
  thermo%pv_total = thermo%pv_virial + thermo%pv_kinetic + thermo%pv_constraint

  CALL timestop ( mfl, handle )

END SUBROUTINE pressure

!!*****
!******************************************************************************
!!****** integrator/update_structure [1.0] *
!!
!!   NAME
!!     update_structure
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE update_structure(struc,task)

  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  CHARACTER ( LEN = * ), INTENT ( IN ) :: task

! Locals
  REAL ( dbl ), ALLOCATABLE, DIMENSION ( :, : ) :: atot
  INTEGER :: natoms, ios, imol, atombase, iat, i, handle
  TYPE (linklist_atoms), POINTER :: current_atom

!------------------------------------------------------------------------------

  CALL timeset ( 'UPDATE_STRUCTURE', 'E', ' ', handle )

  natoms = SIZE ( struc % part )
  ALLOCATE (atot(3,natoms),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'update_structure', 'atot', 3 * natoms )
  atot = zero
  DO imol = 1, size(struc%molecule)
     current_atom => struc%molecule(imol) %ll_atoms
     DO i = 1, struc%molecule(imol) %natoms_mol
        iat = struc%molecule(imol) %atombase + i - 1
        IF (task=='POS') THEN
           atot(1,iat) = current_atom%part%p%r(1)
           atot(2,iat) = current_atom%part%p%r(2)
           atot(3,iat) = current_atom%part%p%r(3)
        ELSE IF (task=='VEL') THEN
           atot(1,iat) = current_atom%part%p%v(1)
           atot(2,iat) = current_atom%part%p%v(2)
           atot(3,iat) = current_atom%part%p%v(3)
        END IF
        current_atom => current_atom%next
     END DO
  END DO

#ifdef __parallel
  CALL mp_sum(atot,int_group)
#endif

  IF (task=='POS') THEN
     DO iat = 1, natoms
        struc%part(iat) %r(1) = atot(1,iat)
        struc%part(iat) %r(2) = atot(2,iat)
        struc%part(iat) %r(3) = atot(3,iat)
     END DO
  ELSE IF (task=='VEL') THEN
     DO iat = 1, natoms
        struc%part(iat) %v(1) = atot(1,iat)
        struc%part(iat) %v(2) = atot(2,iat)
        struc%part(iat) %v(3) = atot(3,iat)
     END DO
  END IF

! DEALLOCATE ( atot, STAT = ios )
! IF ( ios /= 0 ) CALL stop_memory ( 'update_structure', 'atot' )

  CALL timestop ( zero, handle )

END SUBROUTINE update_structure

!!*****
!******************************************************************************
!!****** integrator/set [1.0] *
!!
!!   NAME
!!     set
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set(struc,box,npt,char)

  IMPLICIT NONE

! Arguments
  TYPE (structure_type ), INTENT ( INOUT ) :: struc
  TYPE (cell_type ), INTENT ( INOUT ) :: box
  TYPE (npt_info_type ), DIMENSION ( :, : ), INTENT ( INOUT ) :: npt
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: idim
  INTEGER :: isos

!------------------------------------------------------------------------------

  isos = 0
  IF ( .NOT. ASSOCIATED ( old % v ) ) &
    ALLOCATE ( old % v ( SIZE ( struc % part ), 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  IF ( .NOT. ASSOCIATED ( old % r ) ) &
    ALLOCATE ( old % r ( SIZE ( struc % part ), 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  IF ( .NOT. ASSOCIATED ( old % eps ) )  &
    ALLOCATE ( old % eps ( SIZE ( npt, 1), SIZE ( npt, 2 ) ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  IF ( .NOT. ASSOCIATED ( old % veps ) ) &
    ALLOCATE ( old % veps ( SIZE ( npt, 1), SIZE ( npt, 2 ) ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  IF ( .NOT. ASSOCIATED ( old % h ) ) &
    ALLOCATE (old%h(SIZE(box%hmat,1), SIZE(box%hmat,2)),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )

  SELECT CASE (char)
  CASE ( 'F') ! forward assigning the old
     DO idim = 1, 3
        old % v ( :, idim ) = struc % part ( : ) % v ( idim )
        old % r ( :, idim ) = struc % part ( : ) % r ( idim )
     END DO
     old % eps ( :, : ) = npt ( :, : ) % eps
     old % veps ( :, : ) = npt ( :, : ) % v
     old % h ( :, : ) = box%hmat ( :, : )
  CASE ( 'B') ! back assigning the original variables
     DO idim = 1, 3
        struc % part ( : ) % r ( idim ) = old % r ( :, idim )
        struc % part ( : ) % v ( idim ) = old % v ( :, idim )
     END DO
     npt ( :, : ) % eps = old % eps ( :, : )
     npt ( :, : ) % v = old % veps ( :, : )
     box % hmat ( :, : )=old % h ( :, : )
  END SELECT

END SUBROUTINE set

!!*****
!******************************************************************************
!!****** integrator/npt_f [1.0] *
!!
!!   NAME
!!     npt_f
!!
!!   FUNCTION
!!     Velocity Verlet integrator for the NPT ensemble with fully flexible cell
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE npt_f ( md_env )
  
  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Locals
  INTEGER :: i, j, nnodes, iroll
  REAL ( dbl ) :: dtom
  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
  REAL ( dbl ), SAVE :: eps_0
  REAL ( dbl ) :: roll_tol, infree, trvg
  REAL ( dbl ), DIMENSION (3,3) :: pv_kinetic, u, uh
  REAL ( dbl ), DIMENSION (3) :: arg_r, poly_r, scale_r, e_val, ur, uv, uf
  REAL ( dbl ), DIMENSION (3) :: arg_v, poly_v, scale_v
  TYPE ( cell_type ), POINTER :: box
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  LOGICAL, PARAMETER :: box_change = .TRUE.
  LOGICAL :: first

! assigning the locals
  simpar => md_env % simpar
  thermo => md_env % virial
  struc => md_env % struc ( 1 )
  nhc_part => md_env % nhc_part ( 1 )
  nhc_baro => md_env % nhc_baro ( 1 )
  npt => md_env % npt
  box => md_env % cell

!------------------------------------------------------------------------------
  infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 
  


! Nose-Hoover for the Barostat
     CALL lnhc ( nhc_baro, npt, int_group )
! Nose-Hoover for the Particles
     CALL lnhc ( nhc_part, struc%pnode, int_group )

! setting up for ROLL: saving old variables
  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     iroll = 1
     CALL set(struc, box, npt, 'F' )
     CALL getold ( struc % molecule )
  END IF
  SR: DO WHILE (abs(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP

     IF (simpar%constraint) THEN
        CALL set ( struc, box, npt, 'B' )
     END IF

     CALL get_veps ( box, npt, simpar, struc, thermo, infree )

     trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v 
     
! find eigenvalues and eigenvectors of npt ( :, : ) % v
!

     CALL diagonalise ( matrix =  npt % v, mysize = 3, &
          storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )
!
! first half of velocity verlet
!
     arg_r ( : ) =  0.5_dbl * e_val ( : ) * simpar % dt * &
                    0.5_dbl * e_val ( : ) * simpar % dt
     poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
     scale_r ( : ) = EXP ( 0.5_dbl * simpar % dt * e_val ( : ) )

     arg_v ( : ) =  0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree ) *&
                    0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree )
     poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
     scale_v ( : ) = EXP ( -0.25_dbl * simpar % dt * ( &
                           e_val ( : ) + trvg * infree ) )

     nnodes = size(struc%pnode)
     DO i = 1, nnodes
        dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
! Transform positions and velocities and forces
        ur = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % r )
        uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
        uf = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )

        uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
        uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
        uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom

        ur(1) = ur(1)*scale_r(1)*scale_r(1) +  &
                uv(1)*scale_r(1)*poly_r(1)*simpar%dt
        ur(2) = ur(2)*scale_r(2)*scale_r(2) +  &
                uv(2)*scale_r(2)*poly_r(2)*simpar%dt
        ur(3) = ur(3)*scale_r(3)*scale_r(3) +  &
                uv(3)*scale_r(3)*poly_r(3)*simpar%dt

        struc % pnode ( i ) % p % r ( : ) = MATVEC_3x3 ( u, ur )
        struc % pnode ( i ) % p % v ( : ) = MATVEC_3x3 ( u, uv )
     END DO

     roll_tol = 0.0_dbl
     IF ( simpar % constraint ) CALL shake_roll_control ( struc % molecule, &
          simpar, thermo,roll_tol,iroll,scale_r*poly_r,scale_v*poly_v,  &
          int_group ,u)
  END DO SR

! Update h_mat

  uh = MATMUL_3X3 ( TRANSPOSE_3D ( u ), box % hmat )

  DO i = 1, 3
     DO j = 1, 3
       uh ( i, j ) = uh ( i, j ) * scale_r ( i ) * scale_r ( i )
     END DO
  END DO

  box % hmat = matmul_3x3 ( u, uh )

#if defined(__parallel)
  CALL update_structure(struc,'POS')
#endif

!
! get new forces
!
  CALL force ( glob_env, box_change )
  
!
! second half of velocity verlet
!
  DO i = 1, nnodes
     dtom = 0.5_dbl*simpar%dt/struc%pnode(i) %p%prop%mass
     uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
     uf = MATVEC_3x3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )

     uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
     uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
     uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom

     struc%pnode(i) %p%v ( : ) = MATVEC_3x3(u,uv)
  END DO

  roll_tol = -1.0E+10_dbl
  IF (simpar%constraint) THEN
     first = .TRUE.
     iroll = 1
     CALL set ( struc, box, npt, 'F' )
  END IF

  RR: DO WHILE ( ABS ( roll_tol ) >= 1.0E-10_dbl ) ! RATTLE-ROLL LOOP
     roll_tol = 0.0_dbl
     IF ( simpar%constraint ) CALL rattle_roll_setup ( box, npt, struc, &
        thermo, simpar, scale_v * poly_v, roll_tol, iroll, infree, first )

     CALL get_veps ( box, npt, simpar, struc, thermo, infree )

  END DO RR

#if defined(__parallel)
  CALL update_structure(struc,'VEL')
#endif

! Nose-Hoover for the Particles
  CALL lnhc ( nhc_part, struc%pnode, int_group )
! Nose-Hoover for the Barostat
  CALL lnhc ( nhc_baro, npt, int_group )

  CALL pressure(struc%pnode,thermo)

END SUBROUTINE npt_f

!!*****
!******************************************************************************
!!****** integrator/rattle_roll_setup [1.0] *
!!
!!   NAME
!!     rattle_roll_setup
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE rattle_roll_setup ( box, npt, struc, thermo, simpar, vector_v,  &
                               roll_tol, iroll, infree, first )

  IMPLICIT NONE

! Arguments
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo
  LOGICAL, INTENT ( INOUT ) :: first
  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
  REAL ( dbl ), INTENT ( IN ) :: infree
  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector_v
  INTEGER, INTENT ( INOUT ) :: iroll

! Locals
  TYPE ( npt_info_type ), DIMENSION ( 3, 3 ) :: npt_loc

!------------------------------------------------------------------------------

  IF ( first ) THEN
    CALL get_veps ( box, npt, simpar, struc, thermo, infree )
  END IF
  first = .FALSE.

! assigning local variable
  SELECT CASE ( simpar % ensemble )
  CASE ( 'NPT_I' )
    npt_loc ( :, : ) % v = 0._dbl
    npt_loc ( :, : ) % mass = 0._dbl
    npt_loc ( 1, 1 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 2, 2 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 3, 3 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 1, 1 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 2, 2 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 3, 3 ) % mass = npt ( 1, 1 ) % mass
  CASE ( 'NPT_F' )
    npt_loc = npt
  END SELECT

! resetting

  CALL set ( struc, box, npt, 'B' )
  CALL rattle_roll_control( struc%molecule, simpar, thermo, npt_loc % v, &
                            npt_loc % mass, roll_tol, iroll, vector_v,  &
                            int_group )

END SUBROUTINE rattle_roll_setup

!!*****
!******************************************************************************
!!****** integrator/get_veps [1.0] *
!!
!!   NAME
!!     get_veps
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_veps ( box, npt, simpar, struc, thermo, infree )

! Arguments
  TYPE (cell_type ), INTENT ( IN ) :: box
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), INTENT ( INOUT ) :: thermo
  TYPE (structure_type ), INTENT ( IN ) :: struc
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: infree

! locals
  INTEGER :: i, j, iatom, ii, jj, nnodes
  REAL ( dbl ) :: kin, fdotr
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_kin, unit
! dbg ISOTROPIC LIMIT
!  REAL ( dbl ) :: trace
! dbg

! initializing locals
  unit = 0.0_dbl
  pv_kin = 0.0_dbl
  unit(1,1) = 1.0_dbl
  unit(2,2) = 1.0_dbl
  unit(3,3) = 1.0_dbl

  nnodes = SIZE ( struc % pnode )
  DO i = 1, 3
     DO j = 1, 3
        DO iatom = 1, nnodes
           pv_kin ( i, j ) = pv_kin ( i, j ) + &
                struc % pnode ( iatom ) % p % prop % mass * &
                struc % pnode ( iatom ) % p % v ( i )* &
                struc % pnode ( iatom ) % p % v ( j )
        END DO
     END DO
  END DO

  kin = 0.0_dbl
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO

#if defined(__parallel)
    CALL mp_sum( pv_kin, int_group )
#endif


! updating the constraint virial
  IF ( simpar % constraint )  &
    CALL pv_constraint ( struc % molecule, thermo, 'FULL', int_group )

  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dbl
     DO ii = 1, 3
       fdotr = fdotr + thermo % pv_virial ( ii, ii ) +  &
                       thermo % pv_constraint ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dbl + ( 3.0_dbl * infree ) ) * kin + fdotr - &
       3.0_dbl * simpar % p_ext * box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = thermo % pv_virial ( :, : ) +                   &
       pv_kin ( :, : ) +  thermo % pv_constraint ( :, : ) -  &
       unit ( :, : ) * simpar % p_ext * box % deth +       &
       infree * kin * unit ( :, : )
! dbg ISOTROPIC LIMIT
!    trace = npt ( 1, 1 ) % f + npt ( 2, 2 ) % f + npt ( 3, 3 ) % f
!    trace = trace / 3.0_dbl
!    npt ( :, : ) % f = trace * unit ( :, : )
! dbg
  ENDIF

! update barostat velocities
     npt ( :, : ) % v = npt ( :, : ) % v + &
     0.5_dbl * simpar % dt * npt ( :, : ) % f / npt ( :, : ) % mass

END SUBROUTINE get_veps
!!*****
!******************************************************************************
!!****** integrator/nvt_fist_pimd [1.0] *
!!
!!   NAME
!!     nvt_pimd
!!
!!   FUNCTION
!!     nvt_pimd integrator for particle positions & momenta
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE nvt_pimd ( md_env )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env  
  
! Locals
  TYPE ( simulation_parameters_type ), POINTER :: simpar
  TYPE ( virial_type ), POINTER :: thermo
  TYPE ( structure_type ), POINTER :: struc ( : )
  TYPE ( lnhc_parameters_type ), POINTER :: nhc ( : )
  INTEGER :: i, j, k, nnodes, beads
  REAL ( dbl ) :: dtom, transdt, old_u, fac
  LOGICAL, PARAMETER :: box_change=.false.
  
!------------------------------------------------------------------------------
! Assigning local pointers
  thermo => md_env % virial
  simpar => md_env % simpar
  struc => md_env % struc
  nhc => md_env % nhc_part

  beads = simpar % pimd_params % beads
  
  DO i = 1, beads
     CALL lnhc ( nhc ( i ), struc ( i ) % pnode, int_group )
  END DO
  
!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
  IF ( simpar % constraint ) CALL stop_program('pimd_nvt', &
                                     'contraints not supported atm with PIMD')
#if defined(__parallel)
  CALL stop_program('pimd_nvt','parallel computation not supported with PIMD atm')
#endif

!
! first half of velocity verlet
!
  fac = 0.5_dbl/beads
  nnodes = size(struc(1)%pnode)
  DO j=1, beads
     transdt = simpar%dt*struc(j)%pimd_env%ismeff
     DO i = 1, nnodes
        dtom = fac*transdt/struc(j)%pnode(i) %p%prop%mass
        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
             dtom*struc(j)%pnode(i) %p%f(1)
        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
             dtom*struc(j)%pnode(i) %p%f(2)
        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
             dtom*struc(j)%pnode(i) %p%f(3)
     END DO
  END DO
  
!  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
#if defined(__parallel)
!  CALL update_structure ( struc, 'POS' )
#endif

! Harmonic bead-bead forces, assuming diagonalizing transformation for now

  DO j=1, beads
    DO i=1, nnodes
      old_u = struc(j)%pimd_env%u(1,i)
      struc(j)%pimd_env%u(1,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_env%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(1) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_env%harmonic_prop(2,2)

      old_u = struc(j)%pimd_env%u(2,i)
      struc(j)%pimd_env%u(2,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_env%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(2) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_env%harmonic_prop(2,2)

      old_u = struc(j)%pimd_env%u(3,i)
      struc(j)%pimd_env%u(3,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_env%harmonic_prop(1,2)
      struc(j)%pnode(i)%p%v(3) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_env%harmonic_prop(2,2)
    END DO
  END DO
  IF (simpar%pimd_params%ndelta > 1) THEN
     ! double time step for thermostats
     DO i=1, beads
        CALL yoshida_coef(nhc(i),2.0_dbl*simpar % dt / &
                          simpar % pimd_params % ndelta)
     END DO
     DO k= 2, simpar%pimd_params%ndelta
        DO i = 1, beads
           CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
           DO j=1, nnodes
              old_u = struc(i)%pimd_env%u(1,j)
              struc(i)%pimd_env%u(1,j)= old_u * &
                         struc(i)%pimd_env%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(1) * &
                         struc(i)%pimd_env%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(1) = old_u * &
                         struc(i)%pimd_env%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(1) * &
                         struc(i)%pimd_env%harmonic_prop(2,2)

              old_u = struc(i)%pimd_env%u(2,j)
              struc(i)%pimd_env%u(2,j)= old_u * &
                         struc(i)%pimd_env%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(2) * &
                         struc(i)%pimd_env%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(2) = old_u * &
                         struc(i)%pimd_env%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(2) * &
                         struc(i)%pimd_env%harmonic_prop(2,2)

              old_u = struc(i)%pimd_env%u(3,j)
              struc(i)%pimd_env%u(3,j)= old_u * &
                         struc(i)%pimd_env%harmonic_prop(1,1) &
                       + struc(i)%pnode(j)%p%v(3) * &
                         struc(i)%pimd_env%harmonic_prop(1,2)
              struc(i)%pnode(j)%p%v(3) = old_u * &
                         struc(i)%pimd_env%harmonic_prop(2,1) &
                       + struc(i)%pnode(j)%p%v(3) * &
                         struc(i)%pimd_env%harmonic_prop(2,2)
           END DO
        END DO
     END DO
     ! set time step for thermostats back to normal
     DO i=1, beads
        CALL yoshida_coef(nhc(i),simpar % dt / &
                          simpar % pimd_params % ndelta)
     END DO
  END IF

  CALL u2r (struc, simpar)
!
! get new forces
!
  CALL force ( glob_env, box_change )

  CALL fr2fu ( struc, simpar )
  
!
! second half of velocity verlet
!

  DO j=1, beads
     transdt = 0.5_dbl*simpar%dt*struc(j)%pimd_env%ismeff / beads
     DO i = 1, nnodes
        dtom = transdt/struc(j)%pnode(i) %p%prop%mass
        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
             dtom*struc(j)%pnode(i) %p%f(1)
        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
             dtom*struc(j)%pnode(i) %p%f(2)
        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
             dtom*struc(j)%pnode(i) %p%f(3)
     END DO
  END DO

!  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )
  
  DO i = 1, beads
     CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
  END DO
  
#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
#endif
  
! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!+++ that is not right yet either:
  CALL pressure ( struc ( 1 ) % pnode, thermo )
  
END SUBROUTINE nvt_pimd

!!*****
!******************************************************************************

END MODULE integrator

!******************************************************************************
