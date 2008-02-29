!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Molecular Dynamics
!> \author Teodoro Laino [tlaino] - University of Zurich - 09.2007
! *****************************************************************************
MODULE md_util
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_rattle,&
                                             force_env_shake
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             npe_i_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             reftraj_ensemble
  USE input_cp2k_restarts,             ONLY: update_subsys,&
                                             write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: initialize_velocities,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: md_write_output
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_ener_type,&
                                             md_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: kelvin
  USE shell_potential_types,           ONLY: shell_kind_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_util'

  PUBLIC :: md_output,&
            temperature_control,&
            comvel_control,&
            setup_velocities

CONTAINS

! *****************************************************************************
!> \brief collects the part of the MD that, basically, does the output
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE md_output(md_env,force_env_section,root_section,forced_io,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                root_section
    LOGICAL, INTENT(IN)                      :: forced_io
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_print
    TYPE(section_vals_type), POINTER         :: print_section

    CALL timeset(routineN,handle)
    do_print = .TRUE.
    IF (forced_io) THEN
       print_section => section_vals_get_subs_vals(root_section,"MOTION%MD%PRINT",error=error)
       CALL section_vals_val_get(print_section,"FORCE_LAST",l_val=do_print,error=error)
    END IF
    IF (do_print) THEN
       ! Dumps all files related to the MD run
       CALL md_write_output (md_env, error)
       CALL write_restart(md_env=md_env,root_section=root_section,&
            force_env_section=force_env_section,error=error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE md_output

! *****************************************************************************
!> \brief Computes Ekin, VCOM and Temp for particles
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE reset_vcom(subsys, md_ener, vsubtract, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vsubtract
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, handle, iatom, ikind, &
                                                natom, shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_old, imass_c, imass_s, &
                                                mass, v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(particles, atomic_kind, atomic_kinds, atom_list, shell)
    CALL timeset(routineN,"I","",handle)

    CALL cp_subsys_get(subsys,&
                       atomic_kinds=atomic_kinds,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       core_particles=core_particles,&
                       error=error)

    ekin_old = md_ener%ekin
    ! Possibly subtract a quantity from all velocities
       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list,&
             natom=natom, mass=mass, shell_active=is_shell, shell=shell)
          IF (is_shell) THEN
             tmp = 0.5_dp*vsubtract*mass
             imass_s = 1.0_dp/shell%mass_shell
             imass_c = 1.0_dp/shell%mass_core
             DO iatom = 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                shell_particles%els(shell_index)%v = shell_particles%els(shell_index)%v - tmp*imass_s
                core_particles%els(shell_index)%v  = core_particles%els(shell_index)%v  - tmp*imass_c
                particles%els(atom)%v = particles%els(atom)%v - vsubtract
             END DO
          ELSE
             DO iatom = 1, natom
                atom = atom_list(iatom)
                particles%els(atom)%v = particles%els(atom)%v - vsubtract
             END DO
          END IF
       END DO
    ! Compute Kinetic Energy and COM Velocity
    md_ener%vcom       = 0.0_dp
    md_ener%total_mass = 0.0_dp
    md_ener%ekin       = 0.0_dp
    DO ikind=1,atomic_kinds%n_els
       atomic_kind => atomic_kinds%els(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass, natom=natom)
          v2   = 0.0_dp
          v    = 0.0_dp
          DO iatom = 1, natom
             atom = atom_list(iatom)
             v2 = v2 + SUM(particles%els(atom)%v**2)
             v(1) = v(1) + particles%els(atom)%v(1)
             v(2) = v(2) + particles%els(atom)%v(2)
             v(3) = v(3) + particles%els(atom)%v(3)
          END DO
          md_ener%ekin  = md_ener%ekin + 0.5_dp*mass*v2
          md_ener%vcom(1)    = md_ener%vcom(1) + mass*v(1)
          md_ener%vcom(2)    = md_ener%vcom(2) + mass*v(2)
          md_ener%vcom(3)    = md_ener%vcom(3) + mass*v(3)
          md_ener%total_mass = md_ener%total_mass + REAL(natom,KIND=dp)*mass
    END DO
    md_ener%vcom = md_ener%vcom / md_ener%total_mass
    md_ener%constant = md_ener%constant - ekin_old + md_ener%ekin
    IF (md_ener%nfree /=0) THEN
         md_ener%temp_part = 2.0_dp*md_ener%ekin/REAL(md_ener%nfree,KIND=dp)*kelvin
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE reset_vcom

! *****************************************************************************
!> \brief Scale velocities to get the correct temperature
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE scale_velocity(subsys, md_ener, temp_expected, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_old, imass, mass, scale, &
                                                temp_old, v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ABS(temp_expected - md_ener%temp_part/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_part>0.0_dp) scale = SQRT((temp_expected/md_ener%temp_part)*kelvin)
       ekin_old = md_ener%ekin
       temp_old = md_ener%temp_part
       md_ener%ekin       = 0.0_dp
       md_ener%temp_part  = 0.0_dp
       md_ener%vcom       = 0.0_dp
       md_ener%total_mass = 0.0_dp

       CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles,&
          shell_particles=shell_particles, core_particles=core_particles, error=error)

       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass,&
             natom=natom, shell_active=is_shell, shell=shell)
          IF(is_shell) THEN
             imass = 1.0_dp / mass
             v2 = 0.0_dp
             v  = 0.0_dp
             DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = scale*particles%els(atom)%v
                shell_index = particles%els(atom)%shell_index
                vs = shell_particles%els(shell_index)%v
                vc = core_particles %els(shell_index)%v
                tmp(1) = imass*(vs(1)-vc(1))
                tmp(2) = imass*(vs(2)-vc(2))
                tmp(3) = imass*(vs(3)-vc(3))

                shell_particles%els(shell_index)%v(1) = particles%els(atom)%v(1) + tmp(1)*shell%mass_core 
                shell_particles%els(shell_index)%v(2) = particles%els(atom)%v(2) + tmp(2)*shell%mass_core
                shell_particles%els(shell_index)%v(3) = particles%els(atom)%v(3) + tmp(3)*shell%mass_core

                core_particles%els(shell_index)%v(1) =  particles%els(atom)%v(1) - tmp(1)*shell%mass_shell
                core_particles%els(shell_index)%v(2) =  particles%els(atom)%v(2) - tmp(2)*shell%mass_shell
                core_particles%els(shell_index)%v(3) =  particles%els(atom)%v(3) - tmp(3)*shell%mass_shell

               ! kinetic energy and velocity of COM
                v2 = v2 + SUM(particles%els(atom)%v**2)
                v(1) = v(1) + particles%els(atom)%v(1)
                v(2) = v(2) + particles%els(atom)%v(2)
                v(3) = v(3) + particles%els(atom)%v(3)     
             END DO
          ELSE
             v2 = 0.0_dp
             v = 0.0_dp
             DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(1) = scale*particles%els(atom)%v(1)
                particles%els(atom)%v(2) = scale*particles%els(atom)%v(2)
                particles%els(atom)%v(3) = scale*particles%els(atom)%v(3)
               ! kinetic energy and velocity of COM
                v2 = v2 + SUM(particles%els(atom)%v**2)
                v(1) = v(1) + particles%els(atom)%v(1)
                v(2) = v(2) + particles%els(atom)%v(2)
                v(3) = v(3) + particles%els(atom)%v(3)     
             END DO 
          END IF
          md_ener%ekin  = md_ener%ekin + 0.5_dp*mass*v2
          md_ener%vcom(1)    = md_ener%vcom(1) + mass*v(1)
          md_ener%vcom(2)    = md_ener%vcom(2) + mass*v(2)
          md_ener%vcom(3)    = md_ener%vcom(3) + mass*v(3)
          md_ener%total_mass = md_ener%total_mass + REAL(natom,KIND=dp)*mass
       END DO 
       md_ener%vcom = md_ener%vcom / md_ener%total_mass
       IF(md_ener%nfree /=0) THEN
           md_ener%temp_part = 2.0_dp*md_ener%ekin/REAL(md_ener%nfree,KIND=dp)*kelvin
       END IF
       md_ener%constant = md_ener%constant - ekin_old + md_ener%ekin
       
       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")"Temperature scaled to requested temperature:",&
             temp_old," K  ->",md_ener%temp_part," K"
       END IF
    END IF
  END SUBROUTINE scale_velocity

! *****************************************************************************
!> \brief Scale internal motion of CORE-SHELL model to the correct temperature
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE scale_velocity_internal(subsys, md_ener,  temp_expected, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_internal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_shell_old, fac_mass, &
                                                mass, scale, temp_shell_old, &
                                                v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ABS(temp_expected - md_ener%temp_shell/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_shell>EPSILON(0.0_dp)) scale = SQRT((temp_expected/md_ener%temp_shell)*kelvin)
       ekin_shell_old = md_ener%ekin_shell
       temp_shell_old = md_ener%temp_shell
       md_ener%ekin_shell = 0.0_dp
       md_ener%temp_shell = 0.0_dp

       CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles, shell_particles=shell_particles,&
          core_particles=core_particles, error=error)

       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass, natom=natom,&
             shell_active=is_shell, shell=shell)
          IF(is_shell) THEN
             fac_mass = 1.0_dp/mass 
             v2 = 0.0_dp
             DO iatom= 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                vs  = shell_particles%els(shell_index)%v
                vc  = core_particles%els(shell_index)%v
                v   = particles%els(atom)%v
                tmp(1) = fac_mass*(vc(1)-vs(1))
                tmp(2) = fac_mass*(vc(2)-vs(2))
                tmp(3) = fac_mass*(vc(3)-vs(3))
 
                shell_particles%els(shell_index)%v(1) = v(1) - shell%mass_core*scale*tmp(1)
                shell_particles%els(shell_index)%v(2) = v(2) - shell%mass_core*scale*tmp(2)
                shell_particles%els(shell_index)%v(3) = v(3) - shell%mass_core*scale*tmp(3)

                core_particles%els(shell_index)%v(1) =  v(1) + shell%mass_shell*scale*tmp(1)
                core_particles%els(shell_index)%v(2) =  v(2) + shell%mass_shell*scale*tmp(2)
                core_particles%els(shell_index)%v(3) =  v(3) + shell%mass_shell*scale*tmp(3)

                vs  = shell_particles%els(shell_index)%v
                vc  = core_particles%els(shell_index)%v
                tmp(1) = vc(1) - vs(1)
                tmp(2) = vc(2) - vs(2)
                tmp(3) = vc(3) - vs(3)
                v2 = v2 + SUM(tmp**2)
             END DO
             md_ener%ekin_shell = md_ener%ekin_shell + 0.5_dp*shell%mass_core*shell%mass_shell*fac_mass*v2
          END IF
       END DO
       IF(md_ener%nfree_shell>0)THEN
           md_ener%temp_shell = 2.0_dp*md_ener%ekin_shell/REAL(md_ener%nfree_shell,KIND=dp)*kelvin
       END IF
       md_ener%constant = md_ener%constant - ekin_shell_old + md_ener%ekin_shell

       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature  shell internal motion scaled to requested temperature:",&
               temp_shell_old," K  ->",md_ener%temp_shell," K"
       END IF
    ENDIF

  END SUBROUTINE scale_velocity_internal


! *****************************************************************************
!> \brief Scale barostat velocities to get the desired temperature
!> \par History
!>     MI 02.200s87
! *****************************************************************************
  SUBROUTINE scale_velocity_baro(md_env, md_ener, temp_expected, temp_tol, iw, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nfree
    REAL(KIND=dp)                            :: ekin_old, scale, temp_old
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    NULLIFY( npt, simpar)
    CALL get_md_env ( md_env, simpar = simpar, npt = npt,  error=error)
    IF (ABS(temp_expected - md_ener%temp_baro/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_baro>0.0_dp) scale = SQRT((temp_expected/md_ener%temp_baro)*kelvin)
       ekin_old = md_ener%baro_kin
       temp_old = md_ener%temp_baro
       md_ener%baro_kin       = 0.0_dp
       md_ener%temp_baro  = 0.0_dp
      IF ( simpar%ensemble==npt_i_ensemble .OR. simpar%ensemble==npe_i_ensemble) THEN
         npt ( 1, 1 )%v = npt ( 1, 1 )%v*scale
         md_ener%baro_kin = 0.5_dp * npt ( 1, 1 )%v**2 * npt ( 1, 1 )%mass
       ELSE IF (simpar%ensemble==npt_f_ensemble .OR. simpar%ensemble==npe_f_ensemble) THEN
         md_ener%baro_kin = 0.0_dp
         DO i = 1, 3
            DO j = 1, 3
             npt(i,j)%v = npt(i,j)%v*scale
             md_ener%baro_kin = md_ener%baro_kin + 0.5_dp * npt(i,j)%v**2  * npt ( i, j )%mass 
            END DO
         END DO
      END IF

      nfree = SIZE ( npt, 1 ) * SIZE ( npt, 2 )
      md_ener%temp_baro = 2.0_dp * md_ener%baro_kin / REAL(nfree,dp)*kelvin
      IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature  of barostat motion scaled to requested temperature:",&
               temp_old," K  ->",md_ener%temp_baro," K"
      END IF

   END IF

  END SUBROUTINE scale_velocity_baro

  ! *****************************************************************************
!> \brief Perform all temperature manipulations during a QS MD run. 
!> \par History
!>     Creation (15.09.2003,MK)
!>     adapted to force_env (05.10.2003,fawzi)
!>     Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE temperature_control(simpar, md_env, md_ener,force_env,logger, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'temperature_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,"I","",handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)
    iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
         extension=".mdLog",error=error)

    ! Control the particle motion
    IF (simpar%temp_tol > 0.0_dp ) THEN
       CALL scale_velocity(subsys, md_ener, simpar%temp_ext, simpar%temp_tol, iw, error)
    END IF    
    ! Control the internal core-shell motion
    IF(simpar%temp_sh_tol > 0.0_dp) THEN
       CALL scale_velocity_internal(subsys, md_ener, simpar%temp_sh_ext, simpar%temp_sh_tol, iw, error)
    END IF
    ! Control cell motion
    IF(simpar%temp_baro_tol > 0.0_dp) THEN
       CALL scale_velocity_baro(md_env, md_ener, simpar%temp_baro_ext, simpar%temp_baro_tol, iw, error)
    END IF


    CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
         "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE temperature_control

! *****************************************************************************
!> \brief Set to 0 the velocity of the COM along MD runs, if required. 
!> \par History
!>      Creation (29.04.2007,MI)
!>      Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE comvel_control(simpar,md_ener,force_env,logger, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: temp_old, vel_com
    REAL(KIND=dp), DIMENSION(3)              :: vcom_old
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,"I","",handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)

    ! Print COMVEL and COM Position
    iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%CENTER_OF_MASS",&
         extension=".mdLog",error=error)
    IF (iw>0) THEN
       WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
            "Centre of mass motion (COM):","x =",md_ener%vcom(1),"y =",md_ener%vcom(2),"z =",md_ener%vcom(3)
    END IF
    CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
         "MOTION%MD%PRINT%CENTER_OF_MASS", error=error)

    ! If requested rescale COMVEL
    IF (simpar%comvel_tol > 0.0_dp ) THEN
       iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
            extension=".mdLog",error=error)
       vel_com = SQRT(md_ener%vcom(1)**2+md_ener%vcom(2)**2+md_ener%vcom(3)**2)

       ! Subtract the velocity of the COM, if requested
       IF (vel_com > simpar%comvel_tol) THEN
          temp_old = md_ener%temp_part/kelvin
          vcom_old = md_ener%vcom

          CALL reset_vcom( subsys, md_ener, vsubtract=vcom_old, error=error)

          CALL scale_velocity(subsys, md_ener, temp_old, 0.0_dp, iw, error)

          IF (iw>0) THEN
             WRITE (UNIT=iw,FMT="(T2,'MD| ',A,3F16.10,A)") &
                  "Old VCOM = ",vcom_old(1:3)," a.u.",&
                  "New VCOM =", md_ener%vcom(1:3)," a.u"
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
            "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE comvel_control

! *****************************************************************************
!> \brief Initialize Velocities for MD runs
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE setup_velocities(force_env, simpar, globenv, motion_section, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: apply_cns0, failure, &
                                                shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                print_section, subsys_section

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)
    apply_cns0 = .FALSE.
    IF (simpar%constraint) THEN
       CALL section_vals_val_get(motion_section,"CONSTRAINT%CONSTRAINT_INIT",l_val=apply_cns0,error=error)
    END IF
    ! Always initialize velocities and possibly restart them
    CALL force_env_get(force_env, subsys=subsys, cell=cell, cell_ref=cell_ref, para_env=para_env,&
         force_env_section=force_env_section, error=error )
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)

    CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,molecule_kinds_new=molecule_kinds_new,&
         particles=particles,shell_particles=shell_particles,core_particles=core_particles,&
         error=error)
    
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present,&
         shell_adiabatic=shell_adiabatic)

    NULLIFY(particle_set, shell_particle_set, core_particle_set)
    particle_set => particles%els
    IF(shell_present .AND. shell_adiabatic) THEN
       ! shell-model with constraint NYI
       CPPostcondition(.NOT.simpar%constraint,cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_particles),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(core_particles),cp_failure_level,routineP,error,failure)
       shell_particle_set => shell_particles%els
       core_particle_set  =>  core_particles%els
    END IF

    CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds_new,&
       force_env=force_env, para_env=para_env, globenv=globenv, print_section=print_section,&
       subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
       shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.FALSE., &
       error=error)
    
    ! Apply constraints if required and rescale velocities..
    IF (simpar%ensemble /= reftraj_ensemble) THEN
       IF (apply_cns0) THEN
          CALL force_env_calc_energy_force ( force_env, calc_force=.TRUE.,error=error)          
          CALL force_env_shake(force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
          CALL force_env_rattle(force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
          IF (simpar%do_respa)THEN
             CALL force_env_calc_energy_force (force_env%sub_force_env(1)%force_env,&
                  calc_force=.TRUE.,error=error)
             CALL force_env_shake(force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
             CALL force_env_rattle(force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
          END IF
          ! Reinitialize velocities rescaling properly after rattle
          subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
          CALL update_subsys(subsys_section,force_env,.FALSE.,error)

          CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds_new,&
             force_env=force_env, para_env=para_env, globenv=globenv, print_section=print_section, &
             subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
             shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.TRUE., &
             error=error)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE setup_velocities

END MODULE md_util
