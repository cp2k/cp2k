!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_ensembles [1.0] *
!!
!!   NAME
!!     mc_ensembles
!!
!!   FUNCTION
!!     Perform a Monte Carlo (MC) run using QUICKSTEP and Siepmann group
!!     code
!!
!!   AUTHOR
!!     Matthew J. McGrath  (09.26.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_ensembles
  USE cp_blacs_env,                    ONLY: start_blacs
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_calculate_energy,&
                                             force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type,&
                                             read_mc_section
  USE mc_environment_types,            ONLY: mc_environment_type,&
                                             get_mc_env
  USE mc_misc,                         ONLY: get_water_data,&
                                             final_mc_write
  USE mc_move_control,                 ONLY: mc_moves_type,&
                                             init_moves
  USE mc_moves,                        ONLY: transion,&
                                             trans,&
                                             rot
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ensembles"

  PUBLIC :: nvt

!!***
! *****************************************************************************

CONTAINS

  SUBROUTINE nvt ( mc_env, globenv, error )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(mc_environment_type), POINTER :: mc_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "nvt"

!   *** Local variables ***

    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(force_env_type), POINTER :: force_env, test_env
    TYPE(mc_moves_type),POINTER  :: moves
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys
    TYPE(particle_list_type), POINTER :: particles_old
    REAL(KIND = dp)                             :: rand,old_energy,&
                               ave_energy,bond1,bond2,angle
    INTEGER :: nnstep,iparticle

!   ---------------------------------------------------------------------------

! make it so we don't have to type mc_env for every structure    
  CALL get_mc_env ( mc_env, mc_par = mc_par)

   write(6,*)
   write(6,*)
   write(6,*) 'Beginning Monte Carlo calculation in the NVT ensemble'
   write(6,*)
   write(6,*)

   NULLIFY(blacs_env,test_env,force_env,moves,oldsys)
   NULLIFY(particles_old)

! initialize running average variables
   ave_energy = 0.0d0
   bond1 = 0.0d0
   bond2 = 0.0d0
   angle = 0.0d0

! initialize the moves array
  CALL init_moves(moves)

! initialize the blacs enviroment
  CALL start_blacs(blacs_env,globenv)

! create the force environment needed to calculate the energies
  CALL quickstep_create_force_env(force_env, blacs_env, globenv, error=error)

! find the energy of the initial configuration
   CALL force_env_calculate_energy(force_env)
   CALL force_env_get(force_env,&
             potential_energy=old_energy,error=error)
   WRITE(6,*) 'Initial Energy:                        ', old_energy

!    MC STUFF
    DO nnstep=1,mc_par % nstep
      ! CALL _dump_allocated_delta(%val(0))
      CALL random_number(rand)

      IF (rand .lt. mc_par % pmtraion) THEN
! translate a single atom in the system
         CALL transion( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error=error)
      ELSEIF(rand .lt. mc_par % pmtrans) THEN
! translate a whole molecule in the system
         CALL trans( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error)
      ELSE
! rotate a whole molecule in the system
         CALL rot( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error=error)
      ENDIF

! find bond lengths and angle for current water config
      CALL get_water_data( force_env, angle, bond1, bond2, nnstep,&
                           mc_par%nchain)

! once a cycle, write the coordinates and energy to a file
      IF (MOD(nnstep,mc_par%nunits_tot) == 0 ) THEN

         WRITE(78,*) nnstep,old_energy
         CALL force_env_get(force_env,subsys=oldsys,error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
               particles=particles_old, error=error)

         DO iparticle=1,mc_par%nunits_tot
            WRITE(75,*) nnstep,particles_old%els(iparticle)%r(1:3)
          ENDDO

      ENDIF

    ENDDO

! write a bunch of data to the screen
    CALL final_mc_write(ave_energy,bond1,bond2,angle,moves)
 
  END SUBROUTINE nvt

!**********************************************************************

END MODULE mc_ensembles

!**********************************************************************
