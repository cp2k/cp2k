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
  USE cp_blacs_env,                    ONLY: start_blacs,&
                                             end_blacs
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get,&
                                             cp_subsys_release
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE environment,                     ONLY: cp2k_init
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
                                             init_mc_moves
  USE mc_moves,                        ONLY: mc_atom_translation,&
                                             mc_molecule_translation,&
                                             mc_molecule_rotation,&
                                             mc_volume_move
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ensembles"

  PUBLIC :: mc_traditional_ensemble,gemc_nvt

!!***
! *****************************************************************************

CONTAINS

!**********************************************************************
  SUBROUTINE mc_traditional_ensemble ( mc_env, globenv, error )

!   Purpose: Canonical ensemble (NVT) Monte Carlo run using QUICKSTEP.

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(mc_environment_type), POINTER :: mc_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_traditional_ensemble"

!   *** Local variables ***

    TYPE(force_env_type), POINTER :: force_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(cell_type),POINTER  :: cell
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys
    TYPE(particle_list_type), POINTER :: particles_old
    REAL(KIND = dp)                             :: rand,old_energy,&
                               ave_energy,bond1,bond2,angle,ave_density
    INTEGER :: nnstep,iparticle,crd=75,handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! make it so we don't have to type mc_env for every structure    
  CALL get_mc_env ( mc_env, mc_par = mc_par,force_env=force_env,&
                    blacs_env=blacs_env)

   write(6,*)
   write(6,*)
   write(6,*) 'Beginning Monte Carlo calculation',&
                     ' in a traditional ensemble'
   write(6,*)
   write(6,*)

   NULLIFY(moves,move_updates,oldsys)
   NULLIFY(particles_old)

! initialize running average variables
   ave_energy = 0.0d0
   ave_density = 0.0d0
   bond1 = 0.0d0
   bond2 = 0.0d0
   angle = 0.0d0

! initialize the moves array
  CALL init_mc_moves(moves)

! initialize the arrays for updating maximum move displacements
  CALL init_mc_moves(move_updates)

! find the energy of the initial configuration, and write it to a file
   CALL force_env_calculate_energy(force_env)
   CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
             potential_energy=old_energy,error=error)
   CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
   WRITE(78,*) 'Initial Energy:                        ', old_energy
   CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','INITIAL')

! find out how many steps we're doing
   IF (mc_par%lstop .EQ. 'no') THEN
      mc_par%nstep=mc_par%nstep*mc_par%nunits_tot
   ENDIF

!    MC STUFF
    DO nnstep=1,mc_par % nstep

      CALL random_number(rand)

      IF (rand .lt. mc_par % pmvolume) THEN
! translate a single atom in the system
         CALL mc_volume_move( mc_par,force_env, moves,move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,ave_density,error=error)
      ELSEIF(rand .lt. mc_par % pmtraion) THEN
! translate a single atom in the system
         CALL mc_atom_translation( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error=error)
      ELSEIF(rand .lt. mc_par % pmtrans) THEN
! translate a whole molecule in the system
         CALL mc_molecule_translation( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error=error)
      ELSE
! rotate a whole molecule in the system
         CALL mc_molecule_rotation( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error=error)
      ENDIF

! find bond lengths and angle for current water config
      CALL get_water_data( force_env, angle, bond1, bond2, nnstep,&
                           mc_par%nchain)

! once a cycle, write the coordinates and energy to a file
      IF (MOD(nnstep,mc_par%nunits_tot) == 0 ) THEN

         WRITE(78,*) nnstep,old_energy
         CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
                            error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
               particles=particles_old, error=error)

         CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','POSITIONS')

      ENDIF

    ENDDO

   CALL force_env_calculate_energy(force_env)
   CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
             potential_energy=old_energy,error=error)
   CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','FINAL')

   WRITE(78,*) 'Final Energy:                        ', old_energy

! write a bunch of data to the screen
    CALL final_mc_write(mc_par,ave_energy,bond1,bond2,angle,moves,&
                   ave_density)
 
! deallocate some stuff
    CALL end_blacs(blacs_env,globenv)
    CALL force_env_release(force_env,error=error)
    CALL cp_subsys_release(oldsys(1)%subsys)
    NULLIFY(moves,move_updates)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_traditional_ensemble
!*********************************************************************
!**********************************************************************
  SUBROUTINE gemc_nvt ( mc_env, globenv, error )

!   Purpose: Gibbs canonical ensemble Monte Carlo run using QUICKSTEP.

!   History: - Creation (11.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(mc_environment_type), POINTER :: mc_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "gemc_nvt"

!   *** Local variables ***

    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(force_env_type), POINTER :: force_env_1,force_env_2
    TYPE(mc_moves_type),POINTER  :: moves_1,move_updates_1
    TYPE(mc_moves_type),POINTER  :: moves_2,move_updates_2
    TYPE(cell_type),POINTER  :: cell_1,cell_2
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par_1,mc_par_2
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: oldsys_1,&
                                                        oldsys_2
    TYPE(particle_list_type), POINTER :: particles_old_1,particles_old_2
    REAL(KIND = dp)                             :: rand
    REAL(KIND = dp),DIMENSION(1:2)                :: old_energy,&
                               ave_energy,bond1,bond2,angle,ave_density
    INTEGER,DIMENSION(1:2) ::             istep
    INTEGER :: nnstep,iparticle,crd=75,stat,isos,box_flag,handle
    CHARACTER(LEN=200),DIMENSION(1:2) ::       input_file_names

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! make it so we don't have to type mc_env for every structure    
   CALL get_mc_env ( mc_env, mc_par = mc_par_1,force_env=force_env_1,&
                     blacs_env=blacs_env)

! remember the names of the input files
   input_file_names(:) = (/ globenv%input_file_name,&
                           mc_par_1%box2_file /)

   write(6,*)
   write(6,*)
   write(6,*) 'Beginning Gibbs ensemble Monte Carlo calculation',&
                     ' in the NVT ensemble'
   write(6,*)
   write(6,*)

   NULLIFY(moves_1,move_updates_1,oldsys_1)
   NULLIFY(force_env_2,moves_2,move_updates_2,oldsys_2)
   NULLIFY(particles_old_1,particles_old_2,mc_par_2)

! initialize running average variables
   ave_energy(:) = 0.0d0
   ave_density(:) = 0.0d0
   bond1(:) = 0.0d0
   bond2(:) = 0.0d0
   angle(:) = 0.0d0
   istep(:) = 0

! initialize the moves arrays
  CALL init_mc_moves(moves_1)
  CALL init_mc_moves(moves_2)

! initialize the arrays for updating maximum move displacements
  CALL init_mc_moves(move_updates_1)
  CALL init_mc_moves(move_updates_2)

! find the energy of the initial configuration, and write it to a file
   CALL force_env_calculate_energy(force_env_1)
   CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
             potential_energy=old_energy(1),error=error)
   CALL cp_subsys_get(oldsys_1(1)%subsys, &
             particles=particles_old_1, error=error)
   WRITE(78,*) 'Initial Energy(1):                        ',&
                                           old_energy(1)
   CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        crd,'XMOL','POS','INITIAL BOX 1')

! we need to read in the data for box 2, which is in another file
   globenv%input_file_name=input_file_names(2)

! find the Monte Carlo parameters for box 2
  ALLOCATE ( mc_par_2, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par_2', 0 )
  CALL read_mc_section ( mc_par_2, globenv )

! create the force environment for box 2
  CALL quickstep_create_force_env(force_env_2, blacs_env,&
                                  globenv, error=error)

! find the energy of the initial configuration of box 2, and write it to a file
   CALL force_env_calculate_energy(force_env_2)
   CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
             potential_energy=old_energy(2),error=error)
   CALL cp_subsys_get(oldsys_2(1)%subsys, &
             particles=particles_old_2, error=error)
   WRITE(78,*) 'Initial Energy(2):                        ',&
                                          old_energy(2)
   CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        crd,'XMOL','POS','INITIAL BOX 2')


!    MC STUFF
! make sure a cycle includes the number of units in both boxes,
! unless we're doing a specified number of steps
   IF ( mc_par_1%lstop .EQ. 'no') THEN
       mc_par_1%nstep=mc_par_1%nstep*(mc_par_2%nunits_tot+&
                                   mc_par_1%nunits_tot)
   ENDIF


    DO nnstep=1,mc_par_1 % nstep

      CALL random_number(rand)

      IF (rand .lt. mc_par_1 % pmvolume) THEN
! change the volume of the system
!         CALL volume( mc_par,force_env, moves,move_updates,&
!                        nnstep,ave_energy,old_energy,globenv,&
!                        blacs_env,ave_density,error=error)
      ELSEIF(rand .lt. mc_par_1 % pmtraion) THEN
! translate a single atom in the system
         CALL random_number(rand)
         IF (rand .lt. 0.5d0) THEN
! move in box 1
            istep(1)=istep(1)+1
            globenv%input_file_name=input_file_names(1)
            CALL mc_atom_translation( mc_par_1,force_env_1, moves_1,&
                        move_updates_1,&
                        istep(1),ave_energy(1),old_energy(1),globenv,&
                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
            CALL mc_atom_translation( mc_par_2,force_env_2, moves_2,&
                        move_updates_2,&
                        istep(2),ave_energy(2),old_energy(2),globenv,&
                        error=error)
            box_flag=2
         ENDIF
      ELSEIF(rand .lt. mc_par_1 % pmtrans) THEN
! translate a whole molecule in the system
         CALL random_number(rand)
         IF (rand .lt. 0.5d0) THEN
! move in box 1
            istep(1)=istep(1)+1
            globenv%input_file_name=input_file_names(1)
            CALL mc_molecule_translation( mc_par_1,force_env_1,&
                        moves_1,move_updates_1,&
                        istep(1),ave_energy(1),old_energy(1),globenv,&
                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
            CALL mc_molecule_translation( mc_par_2,force_env_2,&
                        moves_2,move_updates_2,&
                        istep(2),ave_energy(2),old_energy(2),globenv,&
                        error=error)
            box_flag=2
         ENDIF
      ELSE
! rotate a whole molecule in the system
         CALL random_number(rand)
         IF (rand .lt. 0.5d0) THEN
! move in box 1
            istep(1)=istep(1)+1
            globenv%input_file_name=input_file_names(1)
            CALL mc_molecule_rotation( mc_par_1,force_env_1,& 
                        moves_1,move_updates_1,&
                        istep(1),ave_energy(1),old_energy(1),globenv,&
                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
            CALL mc_molecule_rotation( mc_par_2,force_env_2,&
                        moves_2,move_updates_2,&
                        istep(2),ave_energy(2),old_energy(2),globenv,&
                        error=error)
            box_flag=2
         ENDIF
      ENDIF

! find bond lengths and angle for current water config for the box
! that was changed
      SELECT CASE (box_flag)
      CASE(1)
         CALL get_water_data( force_env_1, angle(1), bond1(1), bond2(1), &
                              istep(1), mc_par_1%nchain)
      CASE(2)
         CALL get_water_data( force_env_2, angle(2), bond1(2), bond2(2), &
                              istep(2), mc_par_2%nchain)
      END SELECT

! once a cycle, write the coordinates and energy to a file
      IF (MOD(nnstep,mc_par_1%nunits_tot+&
                          mc_par_2%nunits_tot) == 0 ) THEN

         WRITE(78,*) nnstep,old_energy(1)
         WRITE(78,*) nnstep,old_energy(2)
         CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
             potential_energy=old_energy(1),error=error)
         CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        crd,'XMOL','POS','COORDINATES')

         CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
             potential_energy=old_energy(2),error=error)
         CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        crd,'XMOL','POS','COORDINATES')

      ENDIF

    ENDDO

   CALL force_env_calculate_energy(force_env_1)
   CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
             potential_energy=old_energy(1),error=error)
   CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        crd,'XMOL','POS','FINAL BOX 1')

   WRITE(78,*) 'Final Energy:                      ', old_energy(1)

   CALL force_env_calculate_energy(force_env_2)
   CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
             potential_energy=old_energy(2),error=error)
   CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        crd,'XMOL','POS','FINAL BOX 2')

   WRITE(78,*) 'Final Energy:                      ', old_energy(2)


! write a bunch of data to the screen
   WRITE(6,*) '------------------------------------------------'
   WRITE(6,*) '|                   BOX 1                      |'
   WRITE(6,*) '------------------------------------------------'
   CALL final_mc_write(mc_par_1,ave_energy(1),bond1(1),bond2(1),&
                         angle(1),moves_1,ave_density(1))
   WRITE(6,*) '------------------------------------------------'
   WRITE(6,*) '|                   BOX 2                      |'
   WRITE(6,*) '------------------------------------------------'
   CALL final_mc_write(mc_par_2,ave_energy(2),bond1(2),bond2(2),&
                         angle(2),moves_2,ave_density(2))
 
! deallocate some stuff
    CALL end_blacs(blacs_env,globenv)
    CALL force_env_release(force_env_1,error=error)
    CALL force_env_release(force_env_2,error=error)
    CALL cp_subsys_release(oldsys_1(1)%subsys)
    CALL cp_subsys_release(oldsys_2(1)%subsys)
    NULLIFY(moves_1,moves_2,move_updates_1,move_updates_2)
    NULLIFY(mc_par_2)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE gemc_nvt
!*********************************************************************

END MODULE mc_ensembles

!**********************************************************************
