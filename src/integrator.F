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
!!     Joost VandeVondele (Juli-2003): simple version of isokinetic ensemble
!!
!!   SOURCE
!******************************************************************************

MODULE integrator

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE constraint,                      ONLY: getold,&
                                             rattle_control,&
                                             shake_control
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_mass_atom,&
                                             dyn_coeff_mass_coeff,&
                                             dyn_coeff_mass_global,&
                                             dyn_coeff_mass_kind,&
                                             dyn_coeff_p_type,&
                                             dyn_coeff_set_p_type,&
                                             dyn_coeff_type,&
                                             get_dyn_coeff_set
  USE extended_system_dynamics,        ONLY: lnhc
  USE extended_system_types,           ONLY: lnhc_parameters_type, npt_info_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             force_env_p_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp,&
                                             dp_size
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE message_passing,                 ONLY: mp_comm_dup,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type, get_hinv
  USE termination,                     ONLY: stop_memory, stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
!  USE transformations,                 ONLY: fr2fu,&
!                                             u2r

  USE virial_methods,                   ONLY : virial_evaluate
  USE virial_types,                     ONLY : virial_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "integrator"

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

  PUBLIC :: set_integrator, nve, nvt, isokin, npt_i
  INTERFACE get_veps
     MODULE PROCEDURE get_veps_particle_set, get_veps_velocity
  END INTERFACE

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


    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

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
!!   SYNOPSIS
!!     Subroutine nve_pw(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine nve_pw
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

!SUBROUTINE nve_pw ( md_env )
!  
!  IMPLICIT NONE
!  
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
! 
! Local
! LOGICAL, PARAMETER :: box_change = .FALSE.
! INTEGER :: i, nnodes, ncoef
! REAL ( dbl ) :: dtom, dtocm, cmass, mass
! TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef
! TYPE ( simulation_parameters_type ), POINTER :: simpar
! TYPE ( virial_type ), POINTER :: thermo
! TYPE ( structure_type ), POINTER :: struc
! TYPE ( atomic_kind_type ), POINTER :: atomic_kind

!------------------------------------------------------------------------------
! Associating local pointers
! simpar => md_env % simpar
! thermo => md_env % virial
! struc => md_env % struc ( 1 )
! nhc_coef => md_env % nhc_coef
!
! IF ( simpar % ensemble_coef == "NVT" ) & 
! CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
! IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
! nnodes = SIZE(struc%pnode)
! DO i = 1, nnodes
!    atomic_kind => struc%pnode(i)%p%atomic_kind
!    CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!    dtom = 0.5_dbl*simpar%dt/mass
!    struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!         dtom*struc%pnode(i) %p%f(1)
!    struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!         dtom*struc%pnode(i) %p%f(2)
!    struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!         dtom*struc%pnode(i) %p%f(3)
!    struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!         struc%pnode(i) %p%v(1)*simpar%dt
!    struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!         struc%pnode(i) %p%v(2)*simpar%dt
!    struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!         struc%pnode(i) %p%v(3)*simpar%dt
! END DO

! IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

! ncoef = SIZE ( struc % coef_pos ( 1 ) % pw % cc )
! DO i = 1, ncoef
!    cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
!    dtocm = 0.5_dbl * simpar % dt / cmass
!    struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
!                    + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
!    struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
!            + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
! END DO

!if defined(__parallel)
! CALL update_structure(struc,'POS')
!endif

!
! get new forces
!
! CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
  
!
! second half of velocity verlet
!
! DO i = 1, nnodes
!    atomic_kind => struc%pnode(i)%p%atomic_kind
!    CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!    dtom = 0.5_dbl*simpar%dt/mass
!    struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!         dtom*struc%pnode(i) %p%f(1)
!    struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!         dtom*struc%pnode(i) %p%f(2)
!    struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!         dtom*struc%pnode(i) %p%f(3)
! END DO

! DO i = 1, ncoef
!    cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl ) 
!    dtocm = 0.5_dbl * simpar % dt / cmass
!    struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
!                     + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
! END DO
!
! IF ( simpar % ensemble_coef == "NVT" ) & 
! CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
! IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
!if defined(__parallel)
! CALL update_structure(struc,'VEL')
!endif

! updating the virial
! IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                           'FULL',int_group )

! CALL pressure(struc%pnode,thermo)

!END SUBROUTINE nve_pw

!!*****
!******************************************************************************
!!****** integrator/nve_ao [1.0] *
!!
!!   NAME
!!     nve_ao
!!
!!   SYNOPSIS
!!     Subroutine nve_ao(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine nve_ao
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

!SUBROUTINE nve_ao ( md_env )
!  
!  IMPLICIT NONE
  
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Local
!  LOGICAL, PARAMETER :: box_change = .FALSE.
!  INTEGER :: i, nnodes, ncoef
!  REAL ( dbl ) :: dtom, dtocm, mass
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind

!------------------------------------------------------------------------------
! Associating local structures
!  simpar => md_env % simpar
!  thermo => md_env % virial
!  struc => md_env % struc ( 1 )
!  nhc_coef => md_env % nhc_coef 
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
!  nnodes = SIZE(struc%pnode)
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!          struc%pnode(i) %p%v(1)*simpar%dt
!     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!          struc%pnode(i) %p%v(2)*simpar%dt
!     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!          struc%pnode(i) %p%v(3)*simpar%dt
!  END DO

!  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

!  ncoef = SIZE ( struc % coef_pos ( 1 ) % ao % cr )
!  DO i = 1, ncoef
!     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
!     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
!                     + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
!     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
!             + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
!  END DO

!#if defined(__parallel)
!  CALL update_structure(struc,'POS')
!#endif

!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
  
!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!  END DO

!  DO i = 1, ncoef
!     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
!     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
!                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
!  END DO
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
!  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif

! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )

!  CALL pressure(struc%pnode,thermo)

!END SUBROUTINE nve_ao
!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_pw
!!
!!   SYNOPSIS
!!     Subroutine nvt_pw(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine nvt_pw
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

!SUBROUTINE nvt_pw ( md_env )

!  IMPLICIT NONE
  
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env

! Local
!  LOGICAL, PARAMETER :: box_change = .FALSE.
!  INTEGER :: i, nnodes, ncoef
!  REAL ( dbl ) :: dtom,dtocm,cmass,mass
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
!
!------------------------------------------------------------------------------
! Associating local pointers
!  simpar => md_env % simpar
!  thermo => md_env % virial
!  struc => md_env % struc ( 1 )
!  nhc_part => md_env % nhc_part ( 1 )
!  nhc_coef => md_env % nhc_coef 
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
!  CALL lnhc ( nhc_part, struc % pnode, int_group )

!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
!  nnodes = SIZE(struc%pnode)
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!          struc%pnode(i) %p%v(1)*simpar%dt
!     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!          struc%pnode(i) %p%v(2)*simpar%dt
!     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!          struc%pnode(i) %p%v(3)*simpar%dt
!  END DO

!  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

!  ncoef = SIZE ( struc % coef_pos (1) % pw % cc )
!  DO i = 1, ncoef
!     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
!     dtocm = 0.5_dbl * simpar % dt / cmass
!     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw  % cc ( i ) &
!                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
!     struc % coef_pos ( 1 ) % pw % cc ( i ) = struc % coef_pos ( 1 ) % pw % cc ( i ) &
!                      + struc % coef_vel ( 1 ) % pw % cc ( i ) * simpar % dt
!  END DO

!#if defined(__parallel)
!  CALL update_structure(struc,'POS')
!#endif

!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
  
!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!  END DO

!  DO i = 1, ncoef
!     cmass = REAL ( struc % coef_vel ( 1 ) % pw % mass_cc ( i ), dbl )
!     dtocm = 0.5_dbl * simpar % dt / cmass
!     struc % coef_vel ( 1 ) % pw % cc ( i ) = struc % coef_vel ( 1 ) % pw % cc ( i ) &
!                      + dtocm * struc % coef_force ( 1 ) % pw % cc ( i )
!  END DO
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % pw, int_group )
!
!  CALL lnhc ( nhc_part, struc % pnode, int_group )

!  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif

! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )

!  CALL pressure(struc%pnode,thermo)

!END SUBROUTINE nvt_pw

!!*****
!******************************************************************************
!!****** integrator/nvt_ao [1.0] *
!!
!!   NAME
!!     nvt_ao
!!
!!   SYNOPSIS
!!     Subroutine nvt_ao(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine nvt_ao
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

!SUBROUTINE nvt_ao ( md_env )

!  IMPLICIT NONE
  
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Local
!  LOGICAL, PARAMETER :: box_change = .FALSE.
!  INTEGER :: i, nnodes, ncoef
!  REAL ( dbl ) :: dtom, dtocm, mass
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind

!------------------------------------------------------------------------------
! Associating local pointers
!  simpar => md_env % simpar
!  thermo => md_env % virial
!  struc => md_env % struc ( 1 )
!  nhc_part => md_env % nhc_part ( 1 )
!  nhc_coef => md_env % nhc_coef
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
!  CALL lnhc ( nhc_part, struc % pnode, int_group )

!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!
! first half of velocity verlet
!
!  nnodes = SIZE(struc%pnode)
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!          struc%pnode(i) %p%v(1)*simpar%dt
!     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!          struc%pnode(i) %p%v(2)*simpar%dt
!     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!          struc%pnode(i) %p%v(3)*simpar%dt
!  END DO

!  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)

!  ncoef = SIZE ( struc % coef_pos (1) % ao % cr )
!  DO i = 1, ncoef
!     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
!     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao  % cr ( i ) &
!                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
!     struc % coef_pos ( 1 ) % ao % cr ( i ) = struc % coef_pos ( 1 ) % ao % cr ( i ) &
!                      + struc % coef_vel ( 1 ) % ao % cr ( i ) * simpar % dt
!  END DO

!#if defined(__parallel)
!  CALL update_structure(struc,'POS')
!#endif

!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
  
!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!  END DO

!  DO i = 1, ncoef
!     dtocm = 0.5_dbl * simpar % dt / struc % coef_vel ( 1 ) % ao % mass ( i )
!     struc % coef_vel ( 1 ) % ao % cr ( i ) = struc % coef_vel ( 1 ) % ao % cr ( i ) &
!                      + dtocm * struc % coef_force ( 1 ) % ao % cr ( i )
!  END DO
!
!  IF ( simpar % ensemble_coef == "NVT" ) & 
!  CALL lnhc ( nhc_coef, struc % coef_vel ( 1 ) % ao, int_group )
!
!  CALL lnhc ( nhc_part, struc % pnode, int_group )
!
!  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif

! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!
!  CALL pressure(struc%pnode,thermo)
!
!END SUBROUTINE nvt_ao

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

!SUBROUTINE nve ( md_env )
!  IMPLICIT NONE
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env

! locals
!  LOGICAL, PARAMETER :: box_change = .FALSE.
!  INTEGER :: i, nnodes
!  REAL ( dbl ) :: dtom, mass
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
!
! assigning the locals
!  simpar => md_env % simpar
!  struc => md_env % struc ( 1 )
!  thermo => md_env % virial
!

!------------------------------------------------------------------------------

! IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
!  nnodes = SIZE(struc%pnode)
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!          struc%pnode(i) %p%v(1)*simpar%dt
!     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!          struc%pnode(i) %p%v(2)*simpar%dt
!     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!          struc%pnode(i) %p%v(3)*simpar%dt
!  END DO
!
!  IF (simpar%constraint) CALL shake_control(struc%molecule,simpar)
!#if defined(__parallel)
!  CALL update_structure(struc,'POS')
!#endif

!
! get new forces
!

!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )

!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!  END DO
!
!  IF (simpar%constraint) CALL rattle_control(struc%molecule,simpar)
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif
!
! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!
!  CALL pressure(struc%pnode,thermo)
!
!END SUBROUTINE nve

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
!!     - the local particle lists are used instead of pnode (Sep. 2003,MK)
!!     - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE nve ( md_env )

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "nve_new"

    INTEGER :: icoef, iparticle, iparticle_kind, iparticle_local, istat, &
      isubsys, mass_storage, nmolecule, nmolecule_kind, nparticle, &
      nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dmcoef, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(virial_type), POINTER               :: virial

!   ---------------------------------------------------------------------------

    NULLIFY (dyn_coeff_set,coeffs_of_kind,coeffs,nhc_coef)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_coef=nhc_coef)
    dt = simpar%dt
    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("nve","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set, virial = virial )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els

!     *** Allocate work storage for positions and velocities ***

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

      IF (simpar%constraint) CALL getold(local_molecules, molecule_set,  &
                                         molecule_kind_set, particle_set)

!     *** Velocity Verlet (first part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel (:,iparticle) = particle_set(iparticle)%v(:) +&
                              dm * particle_set(iparticle)%f(:)
          pos (:,iparticle) = particle_set(iparticle)%r(:) +&
                              dt * vel ( :, iparticle )
        END DO
      END DO


      IF (simpar%constraint) CALL shake_control(atomic_kind_set, local_particles, &
                                                local_molecules, molecule_set,   &
                                                molecule_kind_set, particle_set, &
                                                pos, vel, simpar )

!     *** Broadcast the new particle positions ***

      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"pos")

!     *** coefficients
    
      IF(ASSOCIATED(dyn_coeff_set)) THEN 
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                      int_group)
        END IF
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)     
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                        coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
                coeffs%pos(iparticle_local,icoef)=coeffs%pos(iparticle_local,icoef)+ &
                       dt*coeffs%vel(iparticle_local,icoef)
              END DO
            END DO
          END IF
        END DO
      END IF

!     *** Update forces ***

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle) = vel ( :, iparticle) +&
                                dm*particle_set(iparticle)%f(:)
        END DO
      END DO

      IF (simpar%constraint) CALL rattle_control(atomic_kind_set, local_particles, &
                                                 local_molecules, molecule_set, &
                                                 molecule_kind_set, particle_set, &
                                                 vel, simpar )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind) 
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                         coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
              END DO
            END DO
          END IF
        END DO
        IF(ASSOCIATED(nhc_coef))THEN 
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                     int_group)
        END IF
      END IF


!     *** Update virial ***

!MK   IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo,"FULL",int_group)

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

    END DO ! next subsys or sub-system

  END SUBROUTINE nve
!!*****
!******************************************************************************
!!****** integrator/isokin [1.0] *
!!
!!   NAME
!!     isokin
!!
!!   FUNCTION
!!     simplest version of the isokinetic gaussian thermostat
!!
!!   NOTES
!!     - time reversible, and conserves the kinetic energy to machine precision
!!     - is not yet supposed to work for e.g. constraints, our the extended version
!!       of this thermostat
!!       see:
!!        - Zhang F. , JCP 106, 6102 (1997)
!!        - Minary P. et al, JCP 118, 2510 (2003)
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     - Created [2004-07]
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE isokin ( md_env )

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "isokin"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, isubsys, &
      nparticle, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: a, b, dm, ds, dt, K, mass, &
                                                rb, s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!   ---------------------------------------------------------------------------

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p)
                    
    dt = simpar%dt

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("isokin","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys)

    IF (SIZE(subsys)/=1) THEN
       ! presumably one just needs to get e.g. kinetic energy of the full system by looping 
       ! over the subsystems
       CALL stop_program("isokin","integrator",__LINE__,&
       "force_env with more than (or less) than one subsys not supported atm!")
    END IF

    IF (simpar%constraint) THEN
        CALL stop_program("isokin","integrator",__LINE__,&
             "Constraints not yet implemented")
    ENDIF

    do_subsys: DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nparticle = particles%n_els
      particle_set => particles%els

!     *** Allocate work storage for positions and velocities ***
      ALLOCATE (work(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "work",3*nparticle*dp_size)
      work(:,:) = 0.0_dp

      ! compute K,a,b,s,ds
      CALL get_s_ds()

!     *** 'Velocity Verlet' (first part) ***
      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = s/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          particle_set(iparticle)%v(:) = (1.0_dp/ds)*( particle_set(iparticle)%v(:) +&
                                                       dm*particle_set(iparticle)%f(:) )
          work(:,iparticle) = particle_set(iparticle)%r(:) +&
                              dt*particle_set(iparticle)%v(:)
        END DO
      END DO

!     *** Broadcast the new particle positions ***

      CALL mp_sum(work,int_group)

      DO iparticle=1,nparticle
        particle_set(iparticle)%r(:) = work(:,iparticle)
        work(:,iparticle) = 0.0_dp
      END DO

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

      ! compute K,a,b,s,ds
      CALL get_s_ds()

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = s/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          work ( :, iparticle ) = (1.0_dp/ds)*( particle_set(iparticle)%v(:) +&
                                                dm*particle_set(iparticle)%f(:))
        END DO
      END DO

!     *** Broadcast the new particle velocities ***

      CALL mp_sum(work,int_group)

      DO iparticle=1,nparticle
        particle_set(iparticle)%v(:) = work(:,iparticle)
      END DO

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

    END DO do_subsys
  CONTAINS
      SUBROUTINE get_s_ds()
         K=0.0_dp
         a=0.0_dp
         b=0.0_dp
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             K=K+ 0.5_dp*mass *DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%v(:))
             a=a+              DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%f(:))
             b=b+(1.0_dp/mass)*DOT_PRODUCT(particle_set(iparticle)%f(:),particle_set(iparticle)%f(:))
           END DO
         END DO
         CALL mp_sum(K,int_group)
         CALL mp_sum(a,int_group)
         CALL mp_sum(b,int_group)
         a=a/(2.0_dp*K)
         b=b/(2.0_dp*K)
         rb=SQRT(b)
         s =(a/b)*(COSH(dt*rb/2.0_dp)-1) +SINH(dt*rb/2.0_dp)/rb
         ds=(a/b)*(SINH(dt*rb/2.0_dp)*rb)+COSH(dt*rb/2.0_dp)
      END SUBROUTINE get_s_ds
  END SUBROUTINE isokin

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
!!     - the local particle lists are used instead of pnode (Sep. 2003,MK)
!!     - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE nvt ( md_env )
    IMPLICIT NONE

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "nvt"

    INTEGER :: icoef, iparticle, iparticle_kind, iparticle_local, istat, &
      isubsys, mass_storage, nmolecule, nmolecule_kind, nparticle, &
      nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dmcoef, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef
    TYPE(lnhc_parameters_type), POINTER      :: nhc( : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE ( virial_type )                     :: virial

!   ---------------------------------------------------------------------------

    NULLIFY (dyn_coeff_set,coeffs_of_kind,coeffs,nhc_coef)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_part=nhc,nhc_coef=nhc_coef)
    dt = simpar%dt
    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("nvt","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set)

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

! perform Nose'-Hoover chain dynamimcs
      CALL lnhc ( nhc ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, int_group )

      IF (simpar%constraint) CALL getold(local_molecules, molecule_set,  &
                                         molecule_kind_set, particle_set)

!     *** Velocity Verlet (first part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = particle_set(iparticle)%v(:) + dm*particle_set(iparticle)%f(:)
          pos ( :, iparticle ) = particle_set(iparticle)%r(:) + dt*vel(:,iparticle)
        END DO
      END DO


      IF (simpar%constraint) CALL shake_control( atomic_kind_set, local_particles, &
                                                 local_molecules, molecule_set,   &
                                                 molecule_kind_set, particle_set, &
                                                 pos, vel, simpar )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                      int_group)
        END IF
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                        coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
                coeffs%pos(iparticle_local,icoef)=coeffs%pos(iparticle_local,icoef)+ &
                       dt*coeffs%vel(iparticle_local,icoef)
              END DO
            END DO
          END IF
        END DO
      END IF

!     *** Update forces ***

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = vel ( :, iparticle ) + &
                                 dm*particle_set(iparticle)%f(:)
        END DO
      END DO

      IF (simpar%constraint) CALL rattle_control( atomic_kind_set, local_particles, &
                                                  local_molecules, molecule_set, &
                                                  molecule_kind_set, particle_set, &  
                                                  vel, simpar )

! perform Nose'-Hoover chain dynamimcs
      CALL lnhc ( nhc ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, vel, int_group )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
!MK   IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo,"FULL",int_group)

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                         coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
              END DO
            END DO
          END IF
        END DO
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                     int_group)
        END IF
      END IF

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

    END DO ! next subsys or sub-system

  END SUBROUTINE nvt

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

!SUBROUTINE nvt ( md_env )
  
!  IMPLICIT NONE

! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Locals
!  INTEGER :: i, nnodes
!  REAL ( dbl ) :: dtom, mass
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc
!  LOGICAL, PARAMETER :: box_change = .FALSE.
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind

! assigning the locals
!  simpar => md_env % simpar
!  thermo => md_env % virial
!  struc => md_env % struc ( 1 )
!  nhc => md_env % nhc_part ( 1 )
!
!------------------------------------------------------------------------------

!  CALL lnhc ( nhc, struc % pnode, int_group )

!  IF ( simpar % constraint ) CALL getold ( struc % molecule )

!
! first half of velocity verlet
!
!  nnodes = SIZE ( struc % pnode )
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl * simpar%dt / mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!     struc%pnode(i) %p%r(1) = struc%pnode(i) %p%r(1) + &
!          struc%pnode(i) %p%v(1)*simpar%dt
!     struc%pnode(i) %p%r(2) = struc%pnode(i) %p%r(2) + &
!          struc%pnode(i) %p%v(2)*simpar%dt
!     struc%pnode(i) %p%r(3) = struc%pnode(i) %p%r(3) + &
!          struc%pnode(i) %p%v(3)*simpar%dt
!  END DO
!
!  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
!#if defined(__parallel)
!  CALL update_structure ( struc, 'POS' )
!#endif
!
!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind => struc%pnode(i)%p%atomic_kind
!     CALL get_atomic_kind (atomic_kind=atomic_kind,mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     struc%pnode(i) %p%v(1) = struc%pnode(i) %p%v(1) + &
!          dtom*struc%pnode(i) %p%f(1)
!     struc%pnode(i) %p%v(2) = struc%pnode(i) %p%v(2) + &
!          dtom*struc%pnode(i) %p%f(2)
!     struc%pnode(i) %p%v(3) = struc%pnode(i) %p%v(3) + &
!          dtom*struc%pnode(i) %p%f(3)
!  END DO
!
!  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )
!
!  CALL lnhc ( nhc, struc % pnode, int_group )
!
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif

! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!  CALL pressure(struc%pnode,thermo)

!END SUBROUTINE nvt
!!*****
!******************************************************************************
!!****** integrator/npt_i [1.0] *
!!
!!   NAME
!!     npt_i
!!
!!   SYNOPSIS
!!     Subroutine npt_i(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_i
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
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "npt_i"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, &
               isubsys,  nmolecule, nmolecule_kind, nparticle, &
               nparticle_kind, nparticle_local
    INTEGER, POINTER :: itimes
    REAL(KIND=dp)                            :: dm, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part( : ), nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER          :: npt ( :, : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cell_type), POINTER                 :: cell
    TYPE ( virial_type ), POINTER            :: virial
    REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
       e6 = e4/42._dbl, e8 = e6/72._dbl
    REAL ( dbl ), SAVE :: eps_0
    REAL ( dbl ) :: arg_r, poly_r, scale_r, roll_tol
    REAL ( dbl ) :: arg_v, poly_v, scale_v, kin, fdotr, infree
    REAL ( dbl ), DIMENSION (3) :: vector_r, vector_v
    LOGICAL, PARAMETER :: box_change = .TRUE.
    LOGICAL :: first


    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_part=nhc_part, nhc_baro=nhc_baro, npt = npt, &
                    itimes = itimes )
    dt = simpar%dt
    infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("npt_i","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       virial = virial, cell = cell )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      CALL allocate_old ( particle_set, npt )

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

!------------------------------------------------------------------------------

! initialize eps_0 the first time through
      IF ( itimes == 1 ) eps_0 = npt (1,1) % eps

! perform Nose'-Hoover chains for barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

! perform Nose'-Hoover chains for particles
      CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, int_group )

! setting up for ROLL: saving old variables
      roll_tol = -1.0E+10_dbl
!      IF (simpar%constraint) THEN
!        iroll = 1
!        CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
!        CALL getold ( local_molecules, molecule_set,  &
!                      molecule_kind_set, particle_set)
!      ENDIF


      SR: DO WHILE (ABS(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP

!          IF (simpar%constraint) THEN
!             CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
!          END IF
          CALL get_veps ( cell, npt, simpar, atomic_kind_set, particle_set, &
                          local_particles, virial, infree )

!
! first half of velocity verlet
!
          arg_r = ( 0.5_dbl * npt ( 1, 1 ) % v * dt ) * &
                  ( 0.5_dbl * npt ( 1, 1 ) % v * dt )
          poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

          arg_v = ( 0.25_dbl * npt ( 1, 1 ) % v * dt *  &
                  (1._dbl + 3._dbl * infree ) )*( 0.25_dbl * npt ( 1, 1 ) % v * &
                   dt * ( 1._dbl + 3._dbl * infree ) )
          poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

          scale_r = EXP( 0.5_dbl * dt * npt ( 1, 1 ) % v )
          scale_v = EXP( -0.25_dbl * dt * npt ( 1, 1 ) % v * &
                       (1._dbl + 3._dbl * infree ) )
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           dm = 0.5_dp*dt/mass
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel ( :, iparticle ) = particle_set(iparticle)%v(:) * scale_v * scale_v + &
                                    scale_v * poly_v * dm * particle_set(iparticle)%f(:)
             pos ( :, iparticle ) = particle_set(iparticle)%r(:) * scale_r * scale_r + &
                                    scale_r * poly_r * dt *vel(:,iparticle)
           END DO
         END DO



        roll_tol = 0._dbl
        vector_r ( : ) = scale_r * poly_r
        vector_v ( : ) = scale_v * poly_v

!     IF ( simpar%constraint ) CALL shake_roll_control(struc%molecule,simpar, &
!          thermo,roll_tol,iroll,vector_r,vector_v,int_group )
  END DO SR

! Update eps:
      npt ( :, : ) % eps =  npt ( :, : ) % eps + dt * npt ( :, : ) % v

! Update h_mat
      cell % hmat ( :, : ) = cell % hmat ( :, : ) &
      * EXP ( npt ( 1, 1 ) % eps - eps_0 )

      eps_0 =  npt ( 1, 1 ) % eps

! Update the inverse
      CALL get_hinv ( cell )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")


!     *** Update forces ***

      CALL force_env_set(force_env_p(1)%force_env,box_changed=.TRUE.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = vel ( :, iparticle ) * scale_v * scale_v + &
                                 scale_v * poly_v * dm * particle_set(iparticle)%f(:)
        END DO
      END DO

      roll_tol = -1.E10_dbl
!      IF (simpar%constraint) THEN
!         first = .TRUE.
!         iroll = 1
!         CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
!      END IF

      RR: DO WHILE (ABS(roll_tol)>=1.E-10_dbl) ! RATTLE-ROLL LOOP
         roll_tol = 0._dbl
!         IF (simpar%constraint) CALL rattle_roll_setup ( cell, npt, struc, &
!                thermo, simpar, vector_v, roll_tol, iroll, infree, first )
!
         CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, &
                         local_particles, virial, infree )

      END DO RR


! Nose-Hoover for the Particles
      CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, vel, int_group )

! Nose-Hoover for the Barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
!MK   IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo,"FULL",int_group)

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

! Deallocate old variables
    CALL deallocate_old ( )
  END DO

  END SUBROUTINE npt_i

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

SUBROUTINE set ( atomic_kind_set, particle_set, local_particles, cell , npt, char )


! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE ( cell_type ), POINTER :: cell
  TYPE (npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: idim, iparticle_kind, nparticle_kind, nparticle_local, iparticle
  INTEGER :: iparticle_local

!------------------------------------------------------------------------------
  nparticle_kind = SIZE ( atomic_kind_set )
  SELECT CASE (char)
  CASE ( 'F') ! forward assigning the old
     DO iparticle_kind = 1, nparticle_kind
       nparticle_local = local_particles%n_el ( iparticle_kind )
       DO iparticle_local = 1, nparticle_local
         iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
           DO idim = 1, 3
              old % v ( iparticle, idim ) =  particle_set ( iparticle ) % v ( idim )
              old % r ( iparticle, idim ) =  particle_set ( iparticle ) % r ( idim )
           END DO
       END DO
     END DO
     old % eps ( :, : ) = npt ( :, : ) % eps
     old % veps ( :, : ) = npt ( :, : ) % v
     old % h ( :, : ) = cell%hmat ( :, : )
  CASE ( 'B') ! back assigning the original variables
     DO iparticle_kind = 1, nparticle_kind
       nparticle_local = local_particles%n_el ( iparticle_kind )
       DO iparticle_local = 1, nparticle_local
         iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
           DO idim = 1, 3
              particle_set ( iparticle ) % v ( idim ) = old % v ( iparticle, idim )
              particle_set ( iparticle ) % r ( idim ) = old % r ( iparticle, idim ) 
           END DO
       END DO
     END DO
     npt ( :, : ) % eps = old % eps ( :, : )
     npt ( :, : ) % v = old % veps ( :, : )
     cell % hmat ( :, : )=old % h ( :, : )
  END SELECT

END SUBROUTINE set
!!*****
!******************************************************************************
!!****** integrator/allocate_old [1.0] *
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

SUBROUTINE allocate_old ( particle_set, npt )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( npt_info_type ), POINTER :: npt ( :, : )

! Locals
  INTEGER :: isos, idim, jdim, natoms

!------------------------------------------------------------------------------

  natoms = SIZE ( particle_set)
  idim = SIZE ( npt, 1)
  jdim = SIZE ( npt, 2)

  isos = 0
  IF ( .NOT. ASSOCIATED ( old % v ) ) &
    ALLOCATE ( old % v ( natoms, 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  old % v = 0._dbl
  IF ( .NOT. ASSOCIATED ( old % r ) ) &
    ALLOCATE ( old % r ( natoms, 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  old % r = 0._dbl
  IF ( .NOT. ASSOCIATED ( old % eps ) )  &
    ALLOCATE ( old % eps ( idim, jdim ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  old % eps = 0._dbl
  IF ( .NOT. ASSOCIATED ( old % veps ) ) &
    ALLOCATE ( old % veps ( idim, jdim ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  old % veps = 0._dbl
  IF ( .NOT. ASSOCIATED ( old % h ) ) &
    ALLOCATE ( old%h(3, 3),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )
  old % h = 0._dbl


END SUBROUTINE allocate_old

!!*****
!******************************************************************************
!!****** integrator/deallocate_old [1.0] *
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

SUBROUTINE deallocate_old 

  IMPLICIT NONE

! Arguments

! Locals
  INTEGER :: isos 

!------------------------------------------------------------------------------


  isos = 0
  DEALLOCATE ( old % v, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  DEALLOCATE ( old % r, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  DEALLOCATE ( old % eps, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  DEALLOCATE ( old % veps, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  DEALLOCATE ( old % h, STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )
  NULLIFY ( old % v )
  NULLIFY ( old % r )
  NULLIFY ( old % eps )
  NULLIFY ( old % veps )
  NULLIFY ( old % h )

END SUBROUTINE deallocate_old

!!*****
!******************************************************************************
!!****** integrator/npt_f [1.0] *
!!
!!   NAME
!!     npt_f
!!
!!   SYNOPSIS
!!     Subroutine npt_f(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_f
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

!SUBROUTINE npt_f ( md_env )
  
!  IMPLICIT NONE

! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env
  
! Locals
!  INTEGER :: i, j, nnodes, iroll
!  REAL ( dbl ) :: dtom, mass
!  REAL ( dbl ), PARAMETER :: e2 = 1._dbl/6._dbl, e4 = e2/20._dbl, &
!       e6 = e4/42._dbl, e8 = e6/72._dbl
!  REAL ( dbl ), SAVE :: eps_0
!  REAL ( dbl ) :: roll_tol, infree, trvg
!  REAL ( dbl ), DIMENSION (3,3) :: pv_kinetic, u, uh
!  REAL ( dbl ), DIMENSION (3) :: arg_r, poly_r, scale_r, e_val, ur, uv, uf
!  REAL ( dbl ), DIMENSION (3) :: arg_v, poly_v, scale_v
!  TYPE ( cell_type ), POINTER :: box
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
!  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
!  LOGICAL, PARAMETER :: box_change = .TRUE.
!  LOGICAL :: first
!
! assigning the locals
!  simpar => md_env % simpar
!  thermo => md_env % virial
!  struc => md_env % struc ( 1 )
!  nhc_part => md_env % nhc_part ( 1 )
!  nhc_baro => md_env % nhc_baro ( 1 )
!  npt => md_env % npt
!  box => md_env % cell

!------------------------------------------------------------------------------
!  infree = 1.0_dbl / REAL ( simpar % nfree, dbl ) 

! Nose-Hoover for the Barostat
!     CALL lnhc ( nhc_baro, npt, int_group )
! Nose-Hoover for the Particles
!     CALL lnhc ( nhc_part, struc%pnode, int_group )

! setting up for ROLL: saving old variables
!  roll_tol = -1.0E+10_dbl
!  IF (simpar%constraint) THEN
!     iroll = 1
!     CALL set(struc, box, npt, 'F' )
!     CALL getold ( struc % molecule )
!  END IF
!  SR: DO WHILE (ABS(roll_tol)>=1.E-10_dbl) ! SHAKE-ROLL LOOP
!
!     IF (simpar%constraint) THEN
!        CALL set ( struc, box, npt, 'B' )
!     END IF
!
!     CALL get_veps ( box, npt, simpar, struc, thermo, infree )
!
!     trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v 
!     
! find eigenvalues and eigenvectors of npt ( :, : ) % v
!

!     CALL diagonalise ( matrix = npt(:,:)%v, mysize = 3, &
!          storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )
!
! first half of velocity verlet
!
!     arg_r ( : ) =  0.5_dbl * e_val ( : ) * simpar % dt * &
!                    0.5_dbl * e_val ( : ) * simpar % dt
!     poly_r = 1._dbl + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
!     scale_r ( : ) = EXP ( 0.5_dbl * simpar % dt * e_val ( : ) )
!
!     arg_v ( : ) =  0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree ) *&
!                    0.25_dbl * simpar % dt * ( e_val ( : ) + trvg * infree )
!     poly_v = 1._dbl + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
!     scale_v ( : ) = EXP ( -0.25_dbl * simpar % dt * ( &
!                           e_val ( : ) + trvg * infree ) )
!
!     nnodes = SIZE(struc%pnode)
!     DO i = 1, nnodes
!
!        atomic_kind=> struc% pnode(i) %p% atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
!        dtom = 0.5_dbl*simpar%dt/mass
! Transform positions and velocities and forces
!        ur = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % r )
!        uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
!        uf = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )
!
!        uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
!        uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
!        uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom
!
!        ur(1) = ur(1)*scale_r(1)*scale_r(1) +  &
!                uv(1)*scale_r(1)*poly_r(1)*simpar%dt
!        ur(2) = ur(2)*scale_r(2)*scale_r(2) +  &
!                uv(2)*scale_r(2)*poly_r(2)*simpar%dt
!        ur(3) = ur(3)*scale_r(3)*scale_r(3) +  &
!                uv(3)*scale_r(3)*poly_r(3)*simpar%dt
!
!        struc % pnode ( i ) % p % r ( : ) = MATVEC_3x3 ( u, ur )
!        struc % pnode ( i ) % p % v ( : ) = MATVEC_3x3 ( u, uv )
!     END DO

!     roll_tol = 0.0_dbl
!     IF ( simpar % constraint ) CALL shake_roll_control ( struc % molecule, &
!          simpar, thermo,roll_tol,iroll,scale_r*poly_r,scale_v*poly_v,  &
!          int_group ,u)
!  END DO SR

! Update h_mat

!  uh = MATMUL_3X3 ( TRANSPOSE_3D ( u ), box % hmat )

!  DO i = 1, 3
!     DO j = 1, 3
!       uh ( i, j ) = uh ( i, j ) * scale_r ( i ) * scale_r ( i )
!     END DO
!  END DO

!  box % hmat = matmul_3x3 ( u, uh )

!#if defined(__parallel)
!  CALL update_structure(struc,'POS')
!#endif

!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )
  
!
! second half of velocity verlet
!
!  DO i = 1, nnodes
!     atomic_kind=> struc % pnode(i) %p% atomic_kind
!     CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
!     dtom = 0.5_dbl*simpar%dt/mass
!     uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % v )
!     uf = MATVEC_3x3 ( TRANSPOSE_3D ( u ), struc % pnode ( i ) % p % f )
!
!     uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dtom
!     uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dtom
!     uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dtom
!
!     struc%pnode(i) %p%v ( : ) = MATVEC_3x3(u,uv)
!  END DO

!  roll_tol = -1.0E+10_dbl
!  IF (simpar%constraint) THEN
!     first = .TRUE.
!     iroll = 1
!     CALL set ( struc, box, npt, 'F' )
!  END IF

!  RR: DO WHILE ( ABS ( roll_tol ) >= 1.0E-10_dbl ) ! RATTLE-ROLL LOOP
!     roll_tol = 0.0_dbl
!     IF ( simpar%constraint ) CALL rattle_roll_setup ( box, npt, struc, &
!        thermo, simpar, scale_v * poly_v, roll_tol, iroll, infree, first )
!
!     CALL get_veps ( box, npt, simpar, struc, thermo, infree )
!
!  END DO RR

!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif

! Nose-Hoover for the Particles
!  CALL lnhc ( nhc_part, struc%pnode, int_group )
! Nose-Hoover for the Barostat
!  CALL lnhc ( nhc_baro, npt, int_group )

!  CALL pressure(struc%pnode,thermo)

!END SUBROUTINE npt_f

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

!SUBROUTINE rattle_roll_setup ( box, npt, struc, thermo, simpar, vector_v,  &
!                               roll_tol, iroll, infree, first )

!  IMPLICIT NONE

! Arguments
!  TYPE ( cell_type ), INTENT ( INOUT ) :: box
!  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
!  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
!  TYPE ( npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
!  TYPE ( virial_type ), INTENT ( INOUT ) :: thermo
!  LOGICAL, INTENT ( INOUT ) :: first
!  REAL ( dbl ), INTENT ( OUT ) :: roll_tol
!  REAL ( dbl ), INTENT ( IN ) :: infree
!  REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: vector_v
!  INTEGER, INTENT ( INOUT ) :: iroll
!
! Locals
!  TYPE ( npt_info_type ), DIMENSION ( 3, 3 ) :: npt_loc
!
!------------------------------------------------------------------------------
!
!  IF ( first ) THEN
!    CALL get_veps ( box, npt, simpar, struc, thermo, infree )
!  END IF
!  first = .FALSE.

! assigning local variable
!  SELECT CASE ( simpar % ensemble )
!  CASE ( 'NPT_I' )
!    npt_loc ( :, : ) % v = 0._dbl
!    npt_loc ( :, : ) % mass = 0._dbl
!    npt_loc ( 1, 1 ) % v = npt ( 1, 1 ) % v
!    npt_loc ( 2, 2 ) % v = npt ( 1, 1 ) % v
!    npt_loc ( 3, 3 ) % v = npt ( 1, 1 ) % v
!    npt_loc ( 1, 1 ) % mass = npt ( 1, 1 ) % mass
!    npt_loc ( 2, 2 ) % mass = npt ( 1, 1 ) % mass
!    npt_loc ( 3, 3 ) % mass = npt ( 1, 1 ) % mass
!  CASE ( 'NPT_F' )
!    npt_loc = npt
!  END SELECT

! resetting

!  CALL set ( struc, box, npt, 'B' )
!  CALL rattle_roll_control( struc%molecule, simpar, thermo, npt_loc % v, &
!                            npt_loc % mass, roll_tol, iroll, vector_v,  &
!                            int_group )

!END SUBROUTINE rattle_roll_setup

!!*****
!******************************************************************************
!!****** integrator/get_veps_particle_set [1.0] *
!!
!!   NAME
!!     get_veps_particle_set, get_veps_velocity
!!
!!   FUNCTION
!!     Overloaded routine to compute veps given the particles
!!     structure or a local copy of the velocity array
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_veps_particle_set ( box, npt, simpar, atomic_kind_set, particle_set, &
                      local_particles, virial, infree )

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE (cell_type ), INTENT ( IN ) :: box
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), INTENT ( INOUT ) :: virial
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: infree

! locals
  INTEGER :: i, j, iatom, ii, jj
  INTEGER :: iparticle_kind, nparticle_local, iparticle_local, iparticle
  REAL ( dbl ) :: kin, fdotr,mass
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_kin, unit
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
! dbg ISOTROPIC LIMIT
!  REAL ( dbl ) :: trace
! dbg

! initializing locals
  unit = 0.0_dbl
  pv_kin = 0.0_dbl
  unit(1,1) = 1.0_dbl
  unit(2,2) = 1.0_dbl
  unit(3,3) = 1.0_dbl
   
  DO iparticle_kind = 1, SIZE ( atomic_kind_set )
    atomic_kind => atomic_kind_set ( iparticle_kind )
    CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
    nparticle_local = local_particles % n_el ( iparticle_kind )
    DO iparticle_local = 1, nparticle_local
      iparticle = local_particles % list ( iparticle_kind ) % array ( iparticle_local )
      DO i = 1, 3
        DO j = 1, 3
           pv_kin ( i, j ) = pv_kin ( i, j ) + &
               mass* particle_set ( iparticle ) % v ( i )* &
                     particle_set ( iparticle ) % v ( j )
        END DO
      END DO
    END DO 
  END DO

  kin = 0.0_dbl
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO

  CALL mp_sum( pv_kin, int_group )


! updating the constraint virial
!  IF ( simpar % constraint )  &
!    CALL pv_constraint ( struc % molecule, thermo, 'FULL', int_group )

  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dbl
     DO ii = 1, 3
       fdotr = fdotr + virial % pv_virial ( ii, ii ) +  &
                       virial % pv_constraint ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dbl + ( 3.0_dbl * infree ) ) * kin + fdotr - &
       3.0_dbl * simpar % p_ext * box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = virial % pv_virial ( :, : ) +                   &
       pv_kin ( :, : ) +  virial % pv_constraint ( :, : ) -  &
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

END SUBROUTINE get_veps_particle_set

!!*****
!******************************************************************************
!!****** integrator/get_veps_particle_set [1.0] *
!!
!!   NAME
!!     get_veps_particle_set, get_veps_velocity
!!
!!   FUNCTION
!!     Overloaded routine to compute veps given the particles
!!     structure or a local copy of the velocity array
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_veps_velocity ( box, npt, simpar, atomic_kind_set, vel, &
                               local_particles, virial, infree )

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  REAL ( dbl ), INTENT ( INOUT ) :: vel ( :, : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE (cell_type ), INTENT ( IN ) :: box
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), INTENT ( INOUT ) :: virial
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: infree

! locals
  INTEGER :: i, j, iatom, ii, jj
  INTEGER :: iparticle_kind, nparticle_local, iparticle_local, iparticle
  REAL ( dbl ) :: kin, fdotr,mass
  REAL ( dbl ), DIMENSION ( 3, 3 ) :: pv_kin, unit
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
! dbg ISOTROPIC LIMIT
!  REAL ( dbl ) :: trace
! dbg

! initializing locals
  unit = 0.0_dbl
  pv_kin = 0.0_dbl
  unit(1,1) = 1.0_dbl
  unit(2,2) = 1.0_dbl
  unit(3,3) = 1.0_dbl
   
  DO iparticle_kind = 1, SIZE ( atomic_kind_set )
    atomic_kind => atomic_kind_set ( iparticle_kind )
    CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
    nparticle_local = local_particles % n_el ( iparticle_kind )
    DO iparticle_local = 1, nparticle_local
      iparticle = local_particles % list ( iparticle_kind ) % array ( iparticle_local )
      DO i = 1, 3
        DO j = 1, 3
           pv_kin ( i, j ) = pv_kin ( i, j ) + &
               mass* vel ( i, iparticle ) * vel ( j, iparticle )
        END DO
      END DO
    END DO 
  END DO

  kin = 0.0_dbl
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO

  CALL mp_sum( pv_kin, int_group )


! updating the constraint virial
!  IF ( simpar % constraint )  &
!    CALL pv_constraint ( struc % molecule, thermo, 'FULL', int_group )

  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dbl
     DO ii = 1, 3
       fdotr = fdotr + virial % pv_virial ( ii, ii ) +  &
                       virial % pv_constraint ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dbl + ( 3.0_dbl * infree ) ) * kin + fdotr - &
       3.0_dbl * simpar % p_ext * box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = virial % pv_virial ( :, : ) +                   &
       pv_kin ( :, : ) +  virial % pv_constraint ( :, : ) -  &
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

END SUBROUTINE get_veps_velocity

!!*****
!******************************************************************************
!!****** integrator/nvt_pimd [1.0] *
!!
!!   NAME
!!     nvt_pimd
!!
!!   SYNOPSIS
!!     Subroutine nvt_pimd(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine nvt_pimd
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

!SUBROUTINE nvt_pimd ( md_env )
  
!  IMPLICIT NONE
  
! Arguments
!  TYPE ( md_environment_type ), INTENT ( INOUT ) :: md_env  
  
! Locals
!  TYPE ( simulation_parameters_type ), POINTER :: simpar
!  TYPE ( virial_type ), POINTER :: thermo
!  TYPE ( structure_type ), POINTER :: struc ( : )
!  TYPE ( lnhc_parameters_type ), POINTER :: nhc ( : )
!  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
!  INTEGER :: i, j, k, nnodes, beads
!  REAL ( dbl ) :: dtom, transdt, old_u, fac, mass
!  LOGICAL, PARAMETER :: box_change=.FALSE.
  
!------------------------------------------------------------------------------
! Assigning local pointers
!  thermo => md_env % virial
!  simpar => md_env % simpar
!  struc => md_env % struc
!  nhc => md_env % nhc_part
!
!  beads = simpar % pimd_params % beads
  
!  DO i = 1, beads
!     CALL lnhc ( nhc ( i ), struc ( i ) % pnode, int_group )
!  END DO
  
!  IF ( simpar % constraint ) CALL getold ( struc % molecule )
!  IF ( simpar % constraint ) CALL stop_program('pimd_nvt', &
!                                     'contraints not supported atm with PIMD')
!#if defined(__parallel)
!  CALL stop_program('pimd_nvt','parallel computation not supported with PIMD atm')
!#endif

!
! first half of velocity verlet
!
!  fac = 0.5_dbl/beads
!  nnodes = SIZE(struc(1)%pnode)
!  DO j=1, beads
!     transdt = simpar%dt*struc(j)%pimd_env%ismeff
!     DO i = 1, nnodes
!        atomic_kind=> struc(j)%pnode(i) %p% atomic_kind
!        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
!        dtom = fac*transdt/mass
!        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
!             dtom*struc(j)%pnode(i) %p%f(1)
!        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
!             dtom*struc(j)%pnode(i) %p%f(2)
!        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
!             dtom*struc(j)%pnode(i) %p%f(3)
!     END DO
!  END DO
  
!  IF ( simpar % constraint ) CALL shake_control ( struc % molecule, simpar )
!#if defined(__parallel)
!  CALL update_structure ( struc, 'POS' )
!#endif

! Harmonic bead-bead forces, assuming diagonalizing transformation for now

!  DO j=1, beads
!    DO i=1, nnodes
!      old_u = struc(j)%pimd_env%u(1,i)
!      struc(j)%pimd_env%u(1,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
!            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_env%harmonic_prop(1,2)
!      struc(j)%pnode(i)%p%v(1) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
!            + struc(j)%pnode(i)%p%v(1) * struc(j)%pimd_env%harmonic_prop(2,2)
!
!      old_u = struc(j)%pimd_env%u(2,i)
!      struc(j)%pimd_env%u(2,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
!            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_env%harmonic_prop(1,2)
!      struc(j)%pnode(i)%p%v(2) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
!            + struc(j)%pnode(i)%p%v(2) * struc(j)%pimd_env%harmonic_prop(2,2)
!
!      old_u = struc(j)%pimd_env%u(3,i)
!      struc(j)%pimd_env%u(3,i)= old_u * struc(j)%pimd_env%harmonic_prop(1,1)&
!            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_env%harmonic_prop(1,2)
!      struc(j)%pnode(i)%p%v(3) = old_u * struc(j)%pimd_env%harmonic_prop(2,1)&
!            + struc(j)%pnode(i)%p%v(3) * struc(j)%pimd_env%harmonic_prop(2,2)
!    END DO
!  END DO
!  IF (simpar%pimd_params%ndelta > 1) THEN
     ! double time step for thermostats
!     DO i=1, beads
!        CALL yoshida_coef(nhc(i),2.0_dbl*simpar % dt / &
!                          simpar % pimd_params % ndelta)
!     END DO
!     DO k= 2, simpar%pimd_params%ndelta
!        DO i = 1, beads
!           CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
!           DO j=1, nnodes
!              old_u = struc(i)%pimd_env%u(1,j)
!              struc(i)%pimd_env%u(1,j)= old_u * &
!                         struc(i)%pimd_env%harmonic_prop(1,1) &
!                       + struc(i)%pnode(j)%p%v(1) * &
!                         struc(i)%pimd_env%harmonic_prop(1,2)
!              struc(i)%pnode(j)%p%v(1) = old_u * &
!                         struc(i)%pimd_env%harmonic_prop(2,1) &
!                       + struc(i)%pnode(j)%p%v(1) * &
!                         struc(i)%pimd_env%harmonic_prop(2,2)
!
!              old_u = struc(i)%pimd_env%u(2,j)
!              struc(i)%pimd_env%u(2,j)= old_u * &
!                         struc(i)%pimd_env%harmonic_prop(1,1) &
!                       + struc(i)%pnode(j)%p%v(2) * &
!                         struc(i)%pimd_env%harmonic_prop(1,2)
!              struc(i)%pnode(j)%p%v(2) = old_u * &
!                         struc(i)%pimd_env%harmonic_prop(2,1) &
!                       + struc(i)%pnode(j)%p%v(2) * &
!                         struc(i)%pimd_env%harmonic_prop(2,2)
!
!              old_u = struc(i)%pimd_env%u(3,j)
!              struc(i)%pimd_env%u(3,j)= old_u * &
!                         struc(i)%pimd_env%harmonic_prop(1,1) &
!                       + struc(i)%pnode(j)%p%v(3) * &
!                         struc(i)%pimd_env%harmonic_prop(1,2)
!              struc(i)%pnode(j)%p%v(3) = old_u * &
!                         struc(i)%pimd_env%harmonic_prop(2,1) &
!                       + struc(i)%pnode(j)%p%v(3) * &
!                         struc(i)%pimd_env%harmonic_prop(2,2)
!           END DO
!        END DO
!     END DO
     ! set time step for thermostats back to normal
!     DO i=1, beads
!        CALL yoshida_coef(nhc(i),simpar % dt / &
!                          simpar % pimd_params % ndelta)
!     END DO
!  END IF

!  CALL u2r (struc, simpar)
!
! get new forces
!
!  CALL force_env_calculate_force ( md_env%force_env, box_change=box_change )

!  CALL fr2fu ( struc, simpar )
  
!
! second half of velocity verlet
!

!  DO j=1, beads
!     transdt = 0.5_dbl*simpar%dt*struc(j)%pimd_env%ismeff / beads
!     DO i = 1, nnodes
!       atomic_kind=> struc(j) % pnode(i) %p% atomic_kind
!       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
!
!        dtom = transdt/mass
!        struc(j)%pnode(i) %p%v(1) = struc(j)%pnode(i) %p%v(1) + &
!             dtom*struc(j)%pnode(i) %p%f(1)
!        struc(j)%pnode(i) %p%v(2) = struc(j)%pnode(i) %p%v(2) + &
!             dtom*struc(j)%pnode(i) %p%f(2)
!        struc(j)%pnode(i) %p%v(3) = struc(j)%pnode(i) %p%v(3) + &
!             dtom*struc(j)%pnode(i) %p%f(3)
!     END DO
!  END DO

!  IF ( simpar % constraint ) CALL rattle_control ( struc % molecule, simpar )
  
!  DO i = 1, beads
!     CALL lnhc ( nhc(i), struc(i) % pnode, int_group )
!  END DO
  
!#if defined(__parallel)
!  CALL update_structure(struc,'VEL')
!#endif
  
! updating the virial
!  IF (simpar%constraint) CALL pv_constraint(struc%molecule,thermo, &
!                                            'FULL',int_group )
!+++ that is not right yet either:
!  CALL pressure ( struc ( 1 ) % pnode, thermo )
  
!END SUBROUTINE nvt_pimd
!******************************************************************************
SUBROUTINE update_particle_set ( particle_set, int_group,  pos, vel )
                                
    TYPE(particle_type), POINTER             :: particle_set( : )
    INTEGER, INTENT(IN)                      :: int_group
    REAL(dbl), INTENT(INOUT), OPTIONAL       :: pos( :, : ), vel( :, : )

    INTEGER                                  :: iparticle, nparticle

  nparticle = SIZE ( particle_set )

  IF ( ( PRESENT ( pos ) ).AND. ( .NOT. PRESENT ( vel ) ) ) THEN
    CALL mp_sum( pos,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%r(:) = pos(:,iparticle)
    END DO
  ELSEIF ( ( PRESENT ( vel ) ) .AND. ( .NOT. PRESENT ( pos ) ) ) THEN
    CALL mp_sum( vel,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%v(:) = vel(:,iparticle)
    END DO
  ELSEIF ( ( PRESENT ( vel ) ) .AND. ( PRESENT ( pos ) ) ) THEN
    CALL mp_sum( pos,int_group)
    CALL mp_sum( vel,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%v(:) = vel(:,iparticle)
      particle_set(iparticle)%r(:) = pos(:,iparticle)
    END DO
  ENDIF

END SUBROUTINE update_particle_set 
!!*****
!******************************************************************************

END MODULE integrator

!******************************************************************************
