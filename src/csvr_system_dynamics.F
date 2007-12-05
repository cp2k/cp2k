!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/csvr_system_dynamics [1.0] *
!!
!!   NAME
!!     csvr_system_dynamics
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE csvr_system_dynamics
 
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE csvr_system_types,               ONLY: csvr_system_type
  USE csvr_system_utils,               ONLY: rescaling_factor
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: map_info_type,&
                                             npt_info_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE parallel_rng_types,              ONLY: rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE thermostat_utils,                ONLY: ke_region_baro,&
                                             ke_region_coeffs,&
                                             ke_region_particles,&
                                             ke_region_shells,&
                                             vel_rescale_baro,&
                                             vel_rescale_coeffs,&
                                             vel_rescale_particles,&
                                             vel_rescale_shells
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PARAMETER :: debug_this_module=.FALSE.
  PUBLIC :: csvr_particles,&
            csvr_barostat,&
            csvr_coeffs,&
            csvr_shells

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_dynamics'

!!*****
!******************************************************************************
CONTAINS

!******************************************************************************
!!****** csvr_system_dynamics/csvr_coeffs [1.0] *
!!
!!   NAME
!!     csvr_coeffs
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!    
!!   MODIFICATION HISTORY
!!    
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_coeffs ( csvr, dyn_coeff_set, group, error )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csvr_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: flops
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineP,'I','Mflops',handle)
    flops = 0.0_dp
    map_info => csvr%map_info

    ! Evaluate the kinetic energy of the region to thermostat
    CALL ke_region_coeffs(map_info, dyn_coeff_set, coeffs, group, flops, error)

    ! Apply the Canonical Sampling through Velocity Rescaling
    CALL do_csvr(csvr, map_info, error)

    ! Rescale coeffs velocities
    CALL vel_rescale_coeffs(map_info, dyn_coeff_set, coeffs, flops, error)

    ! Evaluate the kinetic energy of the region to thermostat
    CALL ke_region_coeffs(map_info, dyn_coeff_set, coeffs, group, flops, error)

    ! Compute thermostat energy
    CALL do_csvr_eval_energy(csvr, map_info, error)
    
    flops = flops * 1.E-6_dp
    CALL timestop(flops,handle)

  END SUBROUTINE csvr_coeffs

!******************************************************************************
!!****** csvr_system_dynamics/csvr_barostat [1.0] *
!!
!!   NAME
!!     csvr_barostat
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!    
!!   MODIFICATION HISTORY
!!     
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_barostat ( csvr, npt, group, error )

    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csvr_barostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: flops
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineP,'I','Mflops',handle)
    flops = 0.0_dp
    map_info => csvr%map_info

    ! Compute the kinetic energy of the barostat
    CALL ke_region_baro(map_info, npt, group, flops, error)

    ! Apply the Canonical Sampling through Velocity Rescaling
    CALL do_csvr(csvr, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_baro(map_info, npt, flops, error)

    ! Re-Compute the kinetic energy of the barostat
    CALL ke_region_baro(map_info, npt, group, flops, error)

    ! Compute thermostat energy
    CALL do_csvr_eval_energy(csvr, map_info, error)
    
    flops = flops * 1.E-6_dp
    CALL timestop(flops,handle)
  END SUBROUTINE csvr_barostat

!******************************************************************************
!!****** csvr_system_dynamics/csvr_particles [1.0] *
!!
!!   NAME
!!     csvr_particles
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!     
!!   MODIFICATION HISTORY
!!     
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_particles( csvr, molecule_kind_set, molecule_set,  &
       particle_set, local_molecules, group, shell_adiabatic,&
       shell_particle_set, core_particle_set, vel, shell_vel, core_vel,&
       error)

    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csvr_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: my_shell_adiabatic
    REAL(KIND=dp)                            :: flops
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineP,'I','Mflops',handle)
    flops = 0.0_dp
    my_shell_adiabatic = .FALSE.
    IF(PRESENT(shell_adiabatic)) my_shell_adiabatic = shell_adiabatic
    map_info => csvr%map_info
    
    ! Compute the kinetic energy for the region to thermostat
    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, flops, vel, error)
    
    ! Apply the Canonical Sampling through Velocity Rescaling
    CALL do_csvr(csvr, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_particles(map_info, molecule_kind_set, molecule_set, particle_set,&
         local_molecules, my_shell_adiabatic, flops, shell_particle_set,core_particle_set,&
         vel, shell_vel, core_vel, error)
       
    ! Re-Compute the kinetic energy for the region to thermostat
    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, flops, vel, error)

    ! Compute thermostat energy
    CALL do_csvr_eval_energy(csvr, map_info, error)

    flops = flops * 1.E-6_dp
    CALL timestop(flops,handle)
  END SUBROUTINE csvr_particles

!******************************************************************************
!!****** csvr_system_dynamics/csvr_shell_set [1.0] *
!!
!!   NAME
!!    csvr_shell_set 
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!     
!!   MODIFICATION HISTORY
!!     
!!   SOURCE
!******************************************************************************
  SUBROUTINE csvr_shells(csvr, atomic_kind_set, particle_set, local_particles, &
       group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel,&
       error)

    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: group
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csvr_shells', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: flops
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineP,'I','Mflops',handle)
    flops = 0.0_dp
    map_info => csvr%map_info

    ! Compute the kinetic energy of the region to thermostat
    CALL ke_region_shells(map_info, particle_set, atomic_kind_set, local_particles,&
    group, flops, core_particle_set, shell_particle_set, core_vel, shell_vel, error)
    
    ! Apply the Canonical Sampling through Velocity Rescaling
    CALL do_csvr(csvr, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_shells(map_info, atomic_kind_set, particle_set, local_particles, &
         flops, shell_particle_set, core_particle_set, shell_vel, core_vel, vel, error)

    ! Re-Compute the kinetic energy of the region to thermostat
    CALL ke_region_shells(map_info, particle_set, atomic_kind_set, local_particles,&
    group, flops, core_particle_set, shell_particle_set, core_vel, shell_vel, error)

    ! Compute thermostat energy
    CALL do_csvr_eval_energy(csvr, map_info, error)

    flops = flops * 1.E-6_dp
    CALL timestop(flops,handle)
  END SUBROUTINE csvr_shells

!******************************************************************************
!!****** csvr_system_dynamics/do_csvr [1.0] *
!!
!!   NAME
!!    do_csvr
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!     
!!   MODIFICATION HISTORY
!!     
!!   SOURCE
!******************************************************************************
  SUBROUTINE do_csvr(csvr, map_info, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_csvr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap, ndeg
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kin_energy, kin_target, taut
    TYPE(rng_stream_type), POINTER           :: rng_stream

    failure = .FALSE.
    DO i = 1, csvr%loc_num_csvr
       imap = map_info%map_index(i)
       kin_energy =  map_info%s_kin(imap)
       csvr%nvt(i)%region_kin_energy = kin_energy
       kin_energy =  kin_energy * 0.5_dp
       kin_target =  csvr%nvt(i)%nkt * 0.5_dp
       ndeg       =  csvr%nvt(i)%degrees_of_freedom
       taut       =  csvr%tau_csvr
       rng_stream => csvr%nvt(i)%gaussian_rng_stream
       map_info%v_scale(imap) = rescaling_factor(kin_energy, kin_target, ndeg, taut,&
            rng_stream, error)
    END DO

  END SUBROUTINE do_csvr

!******************************************************************************
!!****** csvr_system_dynamics/do_csvr_eval_energy [1.0] *
!!
!!   NAME
!!    do_csvr_eval_energy
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   CREATION DATE
!!     
!!   MODIFICATION HISTORY
!!     
!!   SOURCE
!******************************************************************************
  SUBROUTINE do_csvr_eval_energy(csvr, map_info, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_csvr_eval_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kin_energy_ar, kin_energy_br

    failure = .FALSE.
    DO i = 1, csvr%loc_num_csvr
       imap = map_info%map_index(i)
       kin_energy_br = csvr%nvt(i)%region_kin_energy
       kin_energy_ar = map_info%s_kin(imap)
       csvr%nvt(i)%thermostat_energy =  csvr%nvt(i)%thermostat_energy + &
                                        0.5_dp * (kin_energy_br - kin_energy_ar)
    END DO

  END SUBROUTINE do_csvr_eval_energy

END MODULE csvr_system_dynamics
