!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_util [1.0] *
!!
!!   NAME
!!     md_util
!!
!!   FUNCTION
!!     Utilities for Molecular Dynamics
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 09.2007
!!
!!   MODIFICATION HISTORY
!!  
!!
!!   SOURCE
!******************************************************************************

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
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_rattle,&
                                             force_env_shake
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: reftraj_ensemble
  USE input_cp2k_restarts,             ONLY: update_subsys,&
                                             write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: initialize_velocities,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: md_energy
  USE md_environment_types,            ONLY: md_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE particle_list_types,             ONLY: particle_list_type
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
            compute_ekin_temp_part,&
            setup_velocities

CONTAINS

!!****f* md_util/md_output *
!!
!!   NAME
!!     md_output
!!
!!   FUNCTION
!!     collects the part of the MD that, basically, does the output
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE md_output(md_env,force_env_section,root_section,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)


    CALL md_energy ( md_env, error ) ! also writes trajectories

    CALL write_restart(md_env=md_env,root_section=root_section,&
         force_env_section=force_env_section,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE md_output

!!****f* md_util/compute_ekin_temp_part *
!!
!!   NAME
!!    compute_ekin_temp_part 
!!
!!   FUNCTION
!!    Computes Ekin, VCOM and Temp for particles
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!    Teodoro Laino - University of Zurich - 09.2007 [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE compute_ekin_temp_part(simpar, subsys, vcom, e_kinetic, temperature, do_shell,&
       vsubtract, error)
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: vcom
    REAL(KIND=dp), INTENT(OUT)               :: e_kinetic, temperature
    LOGICAL, INTENT(IN), OPTIONAL            :: do_shell
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(IN), OPTIONAL                   :: vsubtract
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_ekin_temp_part', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, handle, iatom, ikind, &
                                                natom, nfree, shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: do_shell_loc, is_shell, &
                                                shell_adiabatic
    REAL(KIND=dp)                            :: ekin_c, ekin_com, ekin_s, &
                                                fac_massc, fac_masss, mass, &
                                                total_mass, v2
    REAL(KIND=dp), DIMENSION(3)              :: v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(particles, atomic_kind, atomic_kinds, atom_list, shell)
    CALL timeset(routineN,"I","",handle)
    do_shell_loc = .FALSE.
    nfree        = 0
    IF (PRESENT(do_shell))  do_shell_loc = do_shell
    IF (do_shell_loc) THEN
       CALL cp_subsys_get(subsys,&
                          atomic_kinds=atomic_kinds,&
                          particles=particles,&
                          shell_particles=shell_particles,&
                          core_particles=core_particles,&
                          error=error)

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
                                shell_adiabatic=shell_adiabatic)
       IF (shell_adiabatic) nfree = 3*SIZE(shell_particles%els)
    ELSE
       CALL cp_subsys_get(subsys,&
                          atomic_kinds=atomic_kinds,&
                          particles=particles,&
                          error=error)
       nfree = simpar%nfree
    END IF
    ! Possibly subtract a quantity from all velocities
    IF (PRESENT(vsubtract)) THEN
       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               natom=natom,&
                               shell_active=is_shell)
          IF (do_shell_loc.AND.is_shell) THEN
             DO iatom = 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                shell_particles%els(shell_index)%v = shell_particles%els(shell_index)%v - vsubtract
                core_particles%els(shell_index)%v  = core_particles%els(shell_index)%v - vsubtract
                particles%els(atom)%v              = particles%els(atom)%v - vsubtract
             END DO
          ELSE
             DO iatom = 1, natom
                atom = atom_list(iatom)
                particles%els(atom)%v              = particles%els(atom)%v - vsubtract
             END DO
          END IF
       END DO
    END IF
    ! Compute Kinetic Energy and COM Velocity
    vcom       = 0.0_dp
    total_mass = 0.0_dp
    e_kinetic  = 0.0_dp
    DO ikind=1,atomic_kinds%n_els
       atomic_kind => atomic_kinds%els(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            atom_list=atom_list,&
                            mass=mass,&
                            natom=natom,&
                            shell_active=is_shell,&
                            shell=shell)
       IF (do_shell_loc.AND.is_shell) THEN
          fac_masss = shell%mass_shell/mass 
          fac_massc = shell%mass_core/mass
          DO iatom = 1, natom
             atom = atom_list(iatom)
             shell_index = particles%els(atom)%shell_index
             vs     = shell_particles%els(shell_index)%v
             vc     =  core_particles%els(shell_index)%v
             ekin_c = 0.5_dp * shell%mass_core  * (vc(1)**2 + vc(2)**2 + vc(3)**2)
             ekin_s = 0.5_dp * shell%mass_shell * (vs(1)**2 + vs(2)**2 + vs(3)**2)
              
             v =  particles%els(atom)%v
             ekin_com = 0.5_dp*mass * (v(1)**2 + v(2)**2 + v(3)**2)
             e_kinetic = e_kinetic + ekin_c + ekin_s - ekin_com             
          END DO
       ELSE
          v2   = 0.0_dp
          v    = 0.0_dp
          DO iatom = 1, natom
             atom = atom_list(iatom)
             v2 = v2 + particles%els(atom)%v(1)**2 +&
                       particles%els(atom)%v(2)**2 +&
                       particles%els(atom)%v(3)**2
             v(1) = v(1) + particles%els(atom)%v(1)
             v(2) = v(2) + particles%els(atom)%v(2)
             v(3) = v(3) + particles%els(atom)%v(3)
          END DO
          e_kinetic  = e_kinetic + 0.5_dp*mass*v2
       END IF
       vcom(:)    = vcom(:) + mass*v(:)
       total_mass = total_mass + REAL(natom,KIND=dp)*mass
    END DO
    vcom = vcom / total_mass
    temperature = 0.0_dp
    IF (nfree /=0) temperature = 2.0_dp*e_kinetic/REAL(nfree,KIND=dp)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE compute_ekin_temp_part

!!****f* md_util/scale_velocity *
!!
!!   NAME
!!    scale_velocity 
!!
!!   FUNCTION
!!    Scale velocities to get the correct temperature
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!    Teodoro Laino - University of Zurich - 09.2007 [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE scale_velocity(subsys, temp_expected, temperature, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temperature, &
                                                temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: fac_massc, fac_masss, mass, &
                                                scale
    REAL(KIND=dp), DIMENSION(3)              :: vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ABS(temp_expected - temperature) > temp_tol) THEN
       scale = 0.0_dp
       IF (temperature>0.0_dp) scale = SQRT(temp_expected/temperature)
       CALL cp_subsys_get(subsys,&
                          atomic_kinds=atomic_kinds,&
                          particles=particles,&
                          shell_particles=shell_particles,&
                          core_particles=core_particles,&
                          error=error)

       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass,&
                               natom=natom,&
                               shell_active=is_shell,&
                               shell=shell)
          IF(is_shell) THEN
             fac_masss = shell%mass_shell/mass 
             fac_massc = shell%mass_core/mass 
             DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = scale*particles%els(atom)%v
                shell_index = particles%els(atom)%shell_index
                vs = shell_particles%els(shell_index)%v
                vc = core_particles %els(shell_index)%v

                shell_particles%els(shell_index)%v(1) = &
                     particles%els(atom)%v(1) + fac_massc*(vs(1)-vc(1))
                shell_particles%els(shell_index)%v(2) = &
                     particles%els(atom)%v(2) + fac_massc*(vs(2)-vc(2))
                shell_particles%els(shell_index)%v(3) = &
                     particles%els(atom)%v(3) + fac_massc*(vs(3)-vc(3))

                core_particles%els(shell_index)%v(1) = &
                     particles%els(atom)%v(1) + fac_masss*(vc(1)-vs(1))
                core_particles%els(shell_index)%v(2) = &
                     particles%els(atom)%v(2) + fac_masss*(vc(2)-vs(2))
                core_particles%els(shell_index)%v(3) = &
                     particles%els(atom)%v(3) + fac_masss*(vc(3)-vs(3))
             END DO
          ELSE
             DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
             END DO 
          END IF
       END DO 
       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature scaled to requested temperature:",&
               temperature*kelvin," K  ->",temp_expected*kelvin," K"
       END IF
    END IF
  END SUBROUTINE scale_velocity

!!****f* md_util/scale_velocity_internal *
!!
!!   NAME
!!    scale_velocity_internal
!!
!!   FUNCTION
!!    Scale internal motion of CORE-SHELL model to the correct temperature
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!    Teodoro Laino - University of Zurich - 09.2007 [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE scale_velocity_internal(subsys, temp_expected, temperature, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temperature, &
                                                temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_internal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: fac_massc, mass, scale
    REAL(KIND=dp), DIMENSION(3)              :: v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ABS(temp_expected - temperature) > temp_tol) THEN
       scale = 0.0_dp
       IF (temperature>EPSILON(0.0_dp)) scale = SQRT(temp_expected/temperature)

       CALL cp_subsys_get(subsys,&
                          atomic_kinds=atomic_kinds,&
                          particles=particles,&
                          shell_particles=shell_particles,&
                          core_particles=core_particles,&
                          error=error)

       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass,&
                               natom=natom,&
                               shell_active=is_shell,&
                               shell=shell)
          IF(is_shell) THEN
             fac_massc = 1.0_dp/mass 
             DO iatom= 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                vs  = shell_particles%els(shell_index)%v
                vc  = core_particles%els(shell_index)%v
                v   = particles%els(atom)%v
                    
                shell_particles%els(shell_index)%v(1) = v(1) -&
                     shell%mass_core*scale*fac_massc*(vc(1)-vs(1))  
                shell_particles%els(shell_index)%v(2) = v(2) -&
                     shell%mass_core*scale*fac_massc*(vc(2)-vs(2))  
                shell_particles%els(shell_index)%v(3) = v(3) -&
                     shell%mass_core*scale*fac_massc*(vc(3)-vs(3))  

                core_particles%els(shell_index)%v(1) =  v(1) + &
                     shell%mass_shell*scale*fac_massc*(vc(1)-vs(1))
                core_particles%els(shell_index)%v(2) =  v(2) + &
                     shell%mass_shell*scale*fac_massc*(vc(2)-vs(2))
                core_particles%els(shell_index)%v(3) =  v(3) + &
                     shell%mass_shell*scale*fac_massc*(vc(3)-vs(3))

             END DO
          END IF
       END DO
       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature  shell internal motion scaled to requested temperature:",&
               temperature*kelvin," K  ->",temp_expected*kelvin," K"
       END IF
    ENDIF

  END SUBROUTINE scale_velocity_internal

!!****f* md_util/temperature_control *
!!
!!   NAME
!!    temperature_control
!!
!!   FUNCTION
!!    Perform all temperature manipulations during a QS MD run. 
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!    Creation (15.09.2003,MK)
!!    adapted to force_env (05.10.2003,fawzi)
!!    Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
!!
!!*** **********************************************************************
  SUBROUTINE temperature_control(simpar,force_env,logger, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'temperature_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: e_kinetic, temperature
    REAL(KIND=dp), DIMENSION(3)              :: vcom
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,"I","",handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)

    iw = cp_print_key_unit_nr(logger,force_env%root_section,&
         "MOTION%MD%PRINT%PROGRAM_RUN_INFO",extension=".mdLog",error=error)
    ! Control the particle motion
    IF (simpar%temp_tol > 0.0_dp ) THEN
      CALL compute_ekin_temp_part(simpar, subsys, vcom, e_kinetic, temperature, error=error)
      IF (iw>0) THEN
        WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
          "Centre of mass motion (COM):","x =",vcom(1),"y =",vcom(2),"z =",vcom(3)
      END IF
      CALL scale_velocity(subsys, simpar%temp_ext, temperature, simpar%temp_tol, iw, error)
   END IF
   
   ! Control the internal core-shell motion
   IF(simpar%temp_sh_tol > 0.0_dp) THEN
      CALL compute_ekin_temp_part(simpar, subsys, vcom, e_kinetic, temperature, do_shell=.TRUE., error=error)
      IF (iw>0) THEN
         WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
              "Centre of mass motion (COM):","x =",vcom(1),"y =",vcom(2),"z =",vcom(3)
      END IF
      CALL scale_velocity_internal(subsys, simpar%temp_sh_ext, temperature, simpar%temp_sh_tol, iw, error)
    END IF
    CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
         "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE temperature_control

!!****f* md_util/comvel_control *
!!
!!   NAME
!!     comvel_control
!!
!!   FUNCTION
!!     Set to 0 the velocity of the COM along MD runs, if required. 
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     Creation (29.04.2007,MI)
!!     Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
!!
!!*** **********************************************************************
  SUBROUTINE comvel_control(simpar,force_env,logger, error)

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: e_kinetic, temp_new, &
                                                temperature, vel_com
    REAL(KIND=dp), DIMENSION(3)              :: vcom, vcom_new
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,"I","",handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    IF (simpar%comvel_tol > 0.0_dp ) THEN
      CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)

      iw = cp_print_key_unit_nr(logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
           extension=".mdLog",error=error)
      CALL compute_ekin_temp_part(simpar, subsys, vcom, e_kinetic, temperature, error=error)
      vel_com = SQRT(vcom(1)**2+vcom(2)**2+vcom(3)**2)

      ! Subtract the velocity of the COM, if requested
      IF (vel_com > simpar%comvel_tol) THEN
         CALL compute_ekin_temp_part(simpar, subsys, vcom_new, e_kinetic, temp_new, vsubtract= vcom, error=error)
         IF (iw>0) THEN
            WRITE (UNIT=iw,FMT="(T2,'MD| ',A,3F12.6,A)") &
                 "Old VCOM = ",vcom(1:3)," a.u.",&
                 "New VCOM =", vcom_new(1:3)," a.u"
         END IF
         CALL scale_velocity(subsys, temperature, temp_new, 0.0_dp, iw, error)
      END IF
      CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE comvel_control

!!****f* md_util/setup_velocities *
!!
!!   NAME
!!    setup_velocities
!!
!!   FUNCTION
!!    Initialize Velocities for MD runs
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!    Teodoro Laino - University of Zurich - 09.2007 [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE setup_velocities(force_env_p, simpar, globenv, motion_section, error)

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibead
    LOGICAL                                  :: apply_cns0, failure, &
                                                shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER :: force_env_section, print_section, &
      subsys_section, work2_section, work3_section, work_section

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)
    apply_cns0 = .FALSE.
    IF (simpar%constraint) THEN
       CALL section_vals_val_get(motion_section,"CONSTRAINT%CONSTRAIN_INIT",l_val=apply_cns0,error=error)
    END IF
    ! Always initialize velocities and possibly restart them
    DO ibead = 1, SIZE(force_env_p)
       CALL force_env_get(force_env_p(ibead)%force_env, subsys=subsys, &
                          cell=cell, cell_ref=cell_ref, para_env=para_env,&
                          force_env_section=force_env_section, error=error )
       subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)

       CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,molecule_kinds_new=molecule_kinds_new,&
            particles=particles,shell_particles=shell_particles,core_particles=core_particles,&
            error=error)

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present,&
            shell_adiabatic=shell_adiabatic)
       NULLIFY (work_section,work2_section, work3_section)
       work_section  => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)

       IF(shell_present .AND. shell_adiabatic) THEN
          ! shell-model with constraint NYI
          CPPostcondition(.NOT.simpar%constraint,cp_failure_level,routineP,error,failure)
          work2_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
          work3_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
          CALL initialize_velocities(simpar,particles%els,&
                                     molecule_kinds=molecule_kinds_new,&
                                     para_env=para_env,&
                                     globenv=globenv,&
                                     print_section=print_section,&
                                     work_section=work_section,&
                                     work2_section=work2_section,&
                                     work3_section=work3_section,&
                                     shell_part=shell_particles%els,&
                                     core_part=core_particles%els,&
                                     error=error)
       ELSE
          CALL initialize_velocities(simpar,particles%els,&
                                     molecule_kinds=molecule_kinds_new,&
                                     para_env=para_env,&
                                     globenv=globenv,&
                                     print_section=print_section,&
                                     work_section=work_section,&
                                     error=error)
       END IF

       ! Apply constraints if required and rescale velocities..
       IF (simpar%ensemble /= reftraj_ensemble) THEN
          IF (apply_cns0) THEN
             CALL force_env_calc_energy_force ( force_env_p(ibead)%force_env,&
                  calc_force=.TRUE.,error=error)          
             CALL force_env_shake(force_env_p(ibead)%force_env,shake_tol=simpar%shake_tol,&
                  log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
                  dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
             CALL force_env_rattle(force_env_p(ibead)%force_env,shake_tol=simpar%shake_tol,&
                  log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
                  dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
             IF (simpar%do_respa)THEN
                CALL force_env_calc_energy_force ( &
                     force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                     calc_force=.TRUE.,error=error)
                CALL force_env_shake(force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                     shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                     lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
                CALL force_env_rattle(force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                     shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                     lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
             END IF
             ! Reinitialize velocities rescaling properly after rattle
             subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
             CALL update_subsys(subsys_section,force_env_p(ibead)%force_env,.FALSE.,error)
             CALL cp_subsys_get(subsys,molecule_kinds_new=molecule_kinds_new,&
                  particles=particles,error=error)
             NULLIFY (work_section)
             work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
             CALL initialize_velocities(simpar,particles%els,&
                                        molecule_kinds=molecule_kinds_new,&
                                        para_env=para_env,&
                                        globenv=globenv,&
                                        print_section=print_section,&
                                        work_section=work_section,&
                                        force_rescaling=.TRUE.,&
                                        error=error)
          END IF
       END IF
    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE setup_velocities

END MODULE md_util
