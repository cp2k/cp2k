!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert, &
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_to_string, &
                                             cp_logger_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_calculate_energy,&
                                             force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type,&
                                             read_mc_section,&
                                             write_mc_restart
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             set_mc_env,&
                                             mc_environment_type
  USE mc_misc,                         ONLY: final_mc_write,&
                                             get_water_data,&
                                             init_mc_averages,&
                                             mc_averages_type,&
                                             mc_classical_energy
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_energy_restart,&
                                             mc_moves_type
  USE mc_moves,                        ONLY: mc_Quickstep_move,&
                                             mc_atom_translation,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_ge_moves,                     ONLY: mc_ge_quickstep_move,&
                                             mc_ge_volume_move
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
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

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_traditional_ensemble"

    INTEGER                                  :: crd, ene, handle, i, &
                                                iparticle, iw, nnstep,&
                                                atom,mol
    LOGICAL                                  :: PM0, PM1, PM2, PM3, PM4, &
                                                should_stop
    REAL(KIND=dp)                            :: classical_energy_old, &
                                                energy_check, initial_energy, &
                                                old_energy, rand
    TYPE(mc_averages_type), POINTER          :: averages
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(1:3, 1:mc_env%&
      mc_par%nunits_tot)                     :: r_old
    TYPE(mc_moves_type), POINTER             :: move_updates, moves
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(particle_list_type), POINTER        :: particles_old

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

  NULLIFY(moves,move_updates,oldsys)
  NULLIFY(particles_old)

! deal with the different levels of print statements
  PM0 = .FALSE.
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  CALL write_checkpoint_information("entering "//routineN,globenv)
  iw = globenv%scr
  IF((globenv%ionode).AND.(globenv%print%level .GE. 0)) PM0= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.

! make it so we don't have to type mc_env for every structure    
  CALL get_mc_env ( mc_env, mc_par = mc_par,force_env=force_env,&
                    blacs_env=blacs_env)

  IF (PM0) THEN
      WRITE(iw,*)
      WRITE(iw,*)
      WRITE(iw,*) 'Beginning Monte Carlo calculation',&
                     ' in a traditional ensemble'
      WRITE(iw,*)
      WRITE(iw,*)
  ENDIF

! determine the unit to write coordinates and energies to
   ene=mc_par%ene
   crd=mc_par%crd

! initialize running average variables
   energy_check=0.0d0

! initialize the moves array
  CALL init_mc_moves(moves)

! initialize the arrays for updating maximum move displacements
  CALL init_mc_moves(move_updates)

! initialize the averages array
  CALL init_mc_averages(averages)

! find the energy of the initial configuration, and write it to a file
   CALL force_env_calculate_energy(force_env)
   CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
             potential_energy=old_energy,error=error)
   CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
             error=error)

   CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
   OPEN(UNIT=ene)
      WRITE(ene,*) 'Initial Energy:                        ', old_energy
   CLOSE(ene)
   initial_energy=old_energy
   CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','INITIAL')

! record the old coordinates, in case a move is rejected
   DO iparticle=1,mc_par%nunits_tot
      r_old(1:3,iparticle)=particles_old%els(iparticle)%r(1:3)
   ENDDO

! open the file to write cell length to, if we're doing volume moves
   IF(mc_par%pmvolume .GT. 0.0d0) THEN
      OPEN(UNIT=mc_par%cl)
      WRITE(mc_par%cl,*) 'Initial Cell Length:       ',&
               cell%hmat(1,1)*0.529177
      CLOSE(mc_par%cl)
   ENDIF

! find the classical energy if need be
   IF (mc_par % lclassical) THEN
      CALL mc_classical_energy(force_env,moves%classical_energy,&
                               mc_par)
   ENDIF
   classical_energy_old=moves%classical_energy

! find out how many steps we're doing
   IF (mc_par%lstop .EQ. 'no') THEN
      mc_par%nstep=mc_par%nstep*mc_par%nchain
   ENDIF

!    MC STUFF
    DO nnstep=mc_par % nstart+1,mc_par % nstart + mc_par % nstep

      CALL random_number(rand)
      IF (rand .lt. mc_par % pmvolume) THEN
              CALL mc_volume_move( mc_par,force_env,moves,move_updates,&
                        nnstep,averages,old_energy,&
                        classical_energy_old,globenv,&
                        blacs_env,energy_check,r_old,&
                        error=error)
      ELSE
         DO i=1,mc_par%nmoves

            CALL random_number(rand)
            IF(rand .lt. mc_par % pmtraion) THEN
! translate a single atom in the system
               CALL random_number(rand)
               atom=CEILING(rand*mc_par%nunits_tot)
               CALL mc_atom_translation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,atom,&
                        classical_energy_old,error=error)
            ELSEIF(rand .lt. mc_par % pmtrans) THEN
! translate a whole molecule in the system
! call a random number to figure out which molecule we're moving...the
! second line is necessary because we access a matrix with individual
! particles
               CALL random_number(rand)
               mol=CEILING(mc_par%nchain*rand)
               mol=3*mol-2
               CALL mc_molecule_translation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,&
                        classical_energy_old,error=error)
            ELSE
! rotate a whole molecule in the system
               CALL random_number(rand)
               mol=CEILING(mc_par%nchain*rand)
               mol=3*mol-2
               CALL mc_molecule_rotation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,&
                        classical_energy_old,error=error)
            ENDIF
         ENDDO

! now do a Quickstep calculation to see if we accept the sequence
         CALL mc_Quickstep_move(mc_par,force_env,moves,&
                        move_updates,energy_check,r_old,&
                        nnstep,old_energy,averages,&
                        error=error)

      ENDIF

! find bond lengths and angle for current water config
      CALL get_water_data( force_env,averages,nnstep-mc_par%nstart,&
                           mc_par%nchain,mc_par%data)

! write the coordinates and energy to a file
      OPEN(UNIT=ene)
      DO i=1,nnstep-mc_par%nstart
         READ(ene,*)
      ENDDO
      WRITE(ene,*) nnstep,old_energy
      CLOSE(ENE)

      CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
            error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles_old, error=error)

      CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','POSITIONS')

! check to see if EXIT_MC file exists...if so, end the calculation
      CALL external_control(force_env,should_stop,error=error)
      IF (should_stop) EXIT

! write a restart file
      CALL write_mc_restart(nnstep,mc_par,force_env)

    ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this could cause problems on various platforms!!!!!!!!!!!!!!!!
      CALL mc_energy_restart(force_env,old_energy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
             error=error)
   CALL cp_subsys_get(oldsys(1)%subsys, &
            particles=particles_old, error=error)

   CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','FINAL')
                        
   OPEN(UNIT=ene)
   DO i=1,nnstep-mc_par%nstart
      READ(ene,*)
   ENDDO
   WRITE(ene,*) 'Final Energy:                        ', old_energy
   CLOSE(ENE)

! write a bunch of data to the screen
   IF (PM0) THEN
     CALL final_mc_write(mc_par,moves,iw,energy_check,&
                         initial_energy,old_energy,averages)
   ENDIF


! recombine the mc_env
   CALL set_mc_env( mc_env, mc_par = mc_par, force_env=force_env,&
                   blacs_env=blacs_env)

! deallocate some stuff
   CALL cp_subsys_release(oldsys(1)%subsys)
   NULLIFY(moves,move_updates,averages,particles_old)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_traditional_ensemble
!*********************************************************************
!**********************************************************************
  SUBROUTINE gemc_nvt ( mc_env_1, mc_env_2, globenv, error )

!   Purpose: Gibbs canonical ensemble Monte Carlo run using QUICKSTEP.

!   History: - Creation (11.14.2003)

!   ***************************************************************************

    TYPE(mc_environment_type), POINTER       :: mc_env_1
    TYPE(mc_environment_type), POINTER       :: mc_env_2
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = "gemc_nvt"

    CHARACTER(LEN=200), DIMENSION(1:2)       :: input_file_names
    INTEGER                                  :: handle, atom,mol,i,&
                                                isos, nnstep,iw,&
                                                iparticle,com_crd=15
    INTEGER, DIMENSION(1:2)                  :: istep,ene,crd,box_flag
    REAL(KIND=dp)                            :: rand
    REAL(KIND=dp), DIMENSION(1:2)            :: classical_energy_old, &
                                                old_energy,energy_check,&
                                                initial_energy
    LOGICAL                                  :: PM0, PM1, PM2, PM3, PM4, &
                                                should_stop
    TYPE(cell_type), POINTER                 :: cell_1, cell_2
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_1
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_2
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys_1, oldsys_2
    TYPE(force_env_type), POINTER :: force_env_1,force_env_2
    TYPE(mc_moves_type), POINTER             :: move_updates_1, &
                                                move_updates_2, moves_1, &
                                                moves_2
    TYPE(mc_averages_type), POINTER          :: averages_1,averages_2
    TYPE(mc_simulation_parameters_type), &
      POINTER                 :: mc_par_1,mc_par_2
    TYPE(particle_list_type), POINTER        :: particles_old_1, &
                                                particles_old_2
    REAL(KIND=dp), DIMENSION(1:3, 1:mc_env_1%&
      mc_par%nunits_tot)                     :: r_old_1
    REAL(KIND=dp), DIMENSION(1:3, 1:mc_env_2%&
      mc_par%nunits_tot)                     :: r_old_2
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! make it so we don't have to type mc_env for every structure    
   CALL get_mc_env ( mc_env_1, mc_par = mc_par_1,&
                     force_env=force_env_1,&
                     blacs_env=blacs_env_1)
   CALL get_mc_env ( mc_env_2, mc_par = mc_par_2,&
                     force_env=force_env_2,&
                     blacs_env=blacs_env_2)

! open the coordinate file
  OPEN(UNIT=com_crd,FILE='coordinates')

! deal with the different levels of print statements
  PM0 = .FALSE.
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  CALL write_checkpoint_information("entering "//routineN,globenv)
  iw = globenv%scr
  IF((globenv%ionode).AND.(globenv%print%level .GE. 0)) PM0= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
  IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.

! remember the names of the input files
   input_file_names(:) = (/ mc_par_2%box2_file,&
                           mc_par_1%box2_file /)

! determine where we're sending everything
   ene(1)=mc_par_1%ene
   ene(2)=mc_par_2%ene
   crd(1)=mc_par_1%crd
   crd(2)=mc_par_2%crd

   IF (PM0) THEN
      WRITE(iw,*)
      WRITE(iw,*)
      WRITE(iw,*) 'Beginning Monte Carlo calculation',&
                     ' in the NVT-Gibbs ensemble'
      WRITE(iw,*)
      WRITE(iw,*)
   ENDIF

   NULLIFY(moves_1,move_updates_1,oldsys_1)
   NULLIFY(moves_2,move_updates_2,oldsys_2)
   NULLIFY(particles_old_1,particles_old_2)

! initialize running average variables
   istep(:) = 1
   energy_check(:)=0.0d0
   box_flag(:)=0

! initialize the moves arrays
  CALL init_mc_moves(moves_1)
  CALL init_mc_moves(moves_2)

! initialize the arrays for updating maximum move displacements
  CALL init_mc_moves(move_updates_1)
  CALL init_mc_moves(move_updates_2)

! initialize the averages array
  CALL init_mc_averages(averages_1)
  CALL init_mc_averages(averages_2)

! find the energy of the initial configuration, and write it to a file
   CALL force_env_calculate_energy(force_env_1)
   CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
             potential_energy=old_energy(1),error=error)
   CALL cp_subsys_get(oldsys_1(1)%subsys, &
             particles=particles_old_1, error=error)
   OPEN(UNIT=ene(1))
   WRITE(ene(1),*) 'Initial Energy(1):                        ',&
                                           old_energy(1)
   CLOSE(ene(1))
   initial_energy(1)=old_energy(1)
   CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        com_crd,'XMOL','POS','INITIAL BOX 1')

! find the energy of the initial configuration of box 2, and write it to a file
   CALL force_env_calculate_energy(force_env_2)
   CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
             potential_energy=old_energy(2),error=error)
   CALL cp_subsys_get(oldsys_2(1)%subsys, &
             particles=particles_old_2, error=error)
   OPEN(UNIT=ene(2))
   WRITE(ene(2),*) 'Initial Energy(2):                        ',&
                                          old_energy(2)
   CLOSE(ene(2))
   initial_energy(2)=old_energy(2)
   CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        com_crd,'XMOL','POS','INITIAL BOX 2')

! record the initial total energy
   OPEN(UNIT=44,FILE='energies')
   WRITE(44,*) 'Total initial energy:      ',old_energy(2)+old_energy(1)
   CLOSE(44)

! record the old coordinates, in case a move is rejected
   DO iparticle=1,mc_par_1%nunits_tot
      r_old_1(1:3,iparticle)=particles_old_1%els(iparticle)%r(1:3)
   ENDDO
   DO iparticle=1,mc_par_2%nunits_tot
      r_old_2(1:3,iparticle)=particles_old_2%els(iparticle)%r(1:3)
   ENDDO

! open the file to write cell length to, if we're doing volume moves
   IF(mc_par_1%pmvolume .GT. 0.0d0) THEN
      OPEN(UNIT=mc_par_1%cl)
      WRITE(mc_par_1%cl,*) 'Initial Cell Length:       ',&
               cell_1%hmat(1,1)*0.529177
      CLOSE(mc_par_1%cl)
      OPEN(UNIT=mc_par_2%cl)
      WRITE(mc_par_2%cl,*) 'Initial Cell Length:       ',&
               cell_2%hmat(1,1)*0.529177
      CLOSE(mc_par_2%cl)
   ENDIF

! find the classical energy if need be
   IF (mc_par_1 % lclassical) THEN
      CALL mc_classical_energy(force_env_1,moves_1%classical_energy,&
                               mc_par_1)
      classical_energy_old(1)=moves_1%classical_energy
      CALL mc_classical_energy(force_env_2,moves_2%classical_energy,&
                               mc_par_2)
      classical_energy_old(2)=moves_2%classical_energy
   ENDIF

! find out how many steps we're doing
   IF (mc_par_1%lstop .EQ. 'no') THEN
      mc_par_1%nstep=mc_par_1%nstep*(mc_par_1%nchain+mc_par_2%nchain)
   ENDIF

!    MC STUFF

    DO nnstep=mc_par_1%nstart+1,mc_par_1%nstart + mc_par_1%nstep

      CALL random_number(rand)

      IF (rand .lt. mc_par_1 % pmvolume) THEN
         istep(1)=istep(1)+1
         istep(2)=istep(2)+1
         box_flag(1)=1
         box_flag(2)=1
! change the volume of the system
!         CALL volume( mc_par,force_env, moves,move_updates,&
!                        nnstep,ave_energy,old_energy,globenv,&
!                        blacs_env,ave_density,error=error)
      ELSE
         DO i=1,mc_par_1%nmoves
            CALL random_number(rand)
            IF(rand .lt. mc_par_1 % pmtraion) THEN

! translate a single atom in the system
               CALL random_number(rand)
               atom=CEILING(rand*(mc_par_1%nunits_tot+&
                                               mc_par_2%nunits_tot))
               IF (atom .LE. mc_par_1%nunits_tot) THEN
                  box_flag(1)=1
                  CALL mc_atom_translation(  mc_par_1,force_env_1,&
                        moves_1,move_updates_1,nnstep,i,atom,&
                        classical_energy_old(1),error=error)

               ELSE
                  atom=atom-mc_par_1%nunits_tot
                  box_flag(2)=1
                  CALL mc_atom_translation(  mc_par_2,force_env_2,&
                        moves_2,move_updates_2,nnstep,i,atom,&
                        classical_energy_old(2),error=error)
 
               ENDIF
            ELSEIF(rand .lt. mc_par_1 % pmtrans) THEN
! translate a whole molecule in the system
! call a random number to figure out which molecule we're moving...the
! second line is necessary because we access a matrix with individual
! particles
               CALL random_number(rand)
               mol=CEILING((mc_par_1%nchain+mc_par_2%nchain)*rand)
               mol=3*mol-2
               IF (mol .LT. mc_par_1%nunits_tot) THEN
                  box_flag(1)=1
                  CALL mc_molecule_translation(  mc_par_1,force_env_1,&
                        moves_1,move_updates_1,nnstep,i,mol,&
                        classical_energy_old(1),error=error)
               ELSE
                  mol=mol-mc_par_1%nunits_tot
                  box_flag(2)=1
                  CALL mc_molecule_translation(  mc_par_2,force_env_2,&
                        moves_2,move_updates_2,nnstep,i,mol,&
                        classical_energy_old(2),error=error)
               ENDIF
            ELSE
! rotate a whole molecule in the system
               CALL random_number(rand)
               mol=CEILING((mc_par_1%nchain+mc_par_2%nchain)*rand)
               mol=3*mol-2
               IF (mol .LT. mc_par_1%nunits_tot) THEN
                  box_flag(1)=1
                  CALL mc_molecule_rotation(  mc_par_1,force_env_1,&
                        moves_1,move_updates_1,nnstep,i,mol,&
                        classical_energy_old(1),error=error)
               ELSE
                  mol=mol-mc_par_1%nunits_tot
                  box_flag(2)=1
                  CALL mc_molecule_rotation(  mc_par_2,force_env_2,&
                        moves_2,move_updates_2,nnstep,i,mol,&
                        classical_energy_old(2),error=error)
               ENDIF
            ENDIF
         ENDDO         

! now do a Quickstep calculation to see if we accept the sequence
         CALL mc_ge_quickstep_move(mc_par_1,mc_par_2,force_env_1,&
              force_env_2,moves_1,moves_2,move_updates_1,move_updates_2&
              ,energy_check(:),r_old_1,r_old_2,istep,old_energy(:),&
              averages_1,averages_2,box_flag(:),error=error)
!         IF(box_flag(1) == 1) THEN
!            CALL mc_Quickstep_move(mc_par_1,force_env_1,moves_1,&
!                        move_updates_1,energy_check(1),r_old_1,&
!                        nnstep,old_energy(1),averages_1,&
!                        error=error)
!            box_flag(1)=0
!         ENDIF
!         IF(box_flag(2) == 1) THEN
!            istep(2)=istep(2)+1
!            CALL mc_ge_quickstep_move(mc_par_2,force_env_2,moves_2,&
!                        move_updates_2,energy_check(2),r_old_2,&
!                        nnstep,old_energy(2),averages_2,&
!                        error=error)
!            box_flag(2)=0
!         ENDIF

      ENDIF

! find bond lengths and angle for current water config for the box
! that was changed, and write the energy to a file
      IF (box_flag(1) == 1) THEN
         CALL get_water_data( force_env_1,averages_1,istep(1)&
            ,mc_par_1%nchain,mc_par_1%data)
         OPEN(UNIT=ene(1))
         DO i=1,istep(1)
            READ(ene(1),*)
         ENDDO
         WRITE(ene(1),*) nnstep,old_energy(1)
         CLOSE(ene(1))
         istep(1)=istep(1)+1
      ENDIF
      IF (box_flag(2) == 1) THEN
         CALL get_water_data( force_env_2,averages_2,istep(2)&
            ,mc_par_2%nchain,mc_par_2%data)
         OPEN(UNIT=ene(2))
         DO i=1,istep(2)
            READ(ene(2),*)
         ENDDO
         WRITE(ene(2),*) nnstep,old_energy(2)
         CLOSE(ene(2))
         istep(2)=istep(2)+1
      ENDIF


! write the coordinates and energy to a file
      CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
            error=error)
      CALL cp_subsys_get(oldsys_1(1)%subsys, &
            particles=particles_old_1, error=error)
      CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        com_crd,'XMOL','POS','POSITIONS BOX 1')
      CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
            error=error)
      CALL cp_subsys_get(oldsys_2(1)%subsys, &
            particles=particles_old_2, error=error)
      CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        com_crd,'XMOL','POS','POSITIONS BOX 2')
      OPEN(UNIT=44,FILE='energies')
      DO i=1,nnstep-mc_par_1%nstart
         READ(44,*)
      ENDDO
      WRITE(44,*) nnstep,old_energy(2)+old_energy(1)
      CLOSE(44)

! check to see if EXIT_MC file exists...if so, end the calculation
      CALL external_control(force_env_1,should_stop,error=error)
      IF (should_stop) EXIT

! write a restart file
      CALL write_mc_restart(nnstep,mc_par_1,force_env_1)
      CALL write_mc_restart(nnstep,mc_par_2,force_env_2)

! reset the box flags
      box_flag(1)=0
      box_flag(2)=0

    ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this could cause problems on various platforms!!!!!!!!!!!!!!!!
      CALL mc_energy_restart(force_env_1,old_energy(1))
      CALL mc_energy_restart(force_env_2,old_energy(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL force_env_get(force_env_1,subsys=oldsys_1,cell=cell_1,&
             error=error)
   CALL cp_subsys_get(oldsys_1(1)%subsys, &
            particles=particles_old_1, error=error)
   CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                        com_crd,'XMOL','POS','FINAL BOX 1')
   CALL force_env_get(force_env_2,subsys=oldsys_2,cell=cell_2,&
             error=error)
   CALL cp_subsys_get(oldsys_2(1)%subsys, &
            particles=particles_old_2, error=error)
   CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                        com_crd,'XMOL','POS','FINAL BOX 2')

   OPEN(UNIT=ene(1))
   DO i=1,istep(1)
      READ(ene(1),*)
   ENDDO
   WRITE(ene(1),*) 'Final Energy(1):                    ',old_energy(1)
   CLOSE(ene(1))
   OPEN(UNIT=ene(2))
   DO i=1,istep(2)
      READ(ene(2),*)
   ENDDO
   WRITE(ene(2),*) 'Final Energy(2):                    ',old_energy(2)
   CLOSE(ene(2))
! record the final total energy
   OPEN(UNIT=44,FILE='energies')
   DO i=1,nnstep-mc_par_1%nstart
      READ(44,*)
   ENDDO
   WRITE(44,*) 'Final total energy:      ',old_energy(2)+old_energy(1)
   CLOSE(44)

! write a bunch of data to the screen
   IF (PM0) THEN
      WRITE(iw,*) '------------------------------------------------'
      WRITE(iw,*) '|                   BOX 1                      |'
      WRITE(iw,*) '------------------------------------------------'
      CALL final_mc_write(mc_par_1,moves_1,iw,energy_check(1),&
                         initial_energy(1),old_energy(1),averages_1)
      WRITE(iw,*) '------------------------------------------------'
      WRITE(iw,*) '|                   BOX 2                      |'
      WRITE(iw,*) '------------------------------------------------'
      CALL final_mc_write(mc_par_2,moves_2,iw,energy_check(2),&
                         initial_energy(2),old_energy(2),averages_2)
   ENDIF
 
! deallocate some stuff
   CALL set_mc_env ( mc_env_1, mc_par = mc_par_1,&
                     force_env=force_env_1,&
                     blacs_env=blacs_env_1)
   CALL set_mc_env ( mc_env_2, mc_par = mc_par_2,&
                     force_env=force_env_2,&
                     blacs_env=blacs_env_2)
    CALL cp_subsys_release(oldsys_1(1)%subsys)
    NULLIFY(moves_1,move_updates_1,averages_1,particles_old_1)
    CALL cp_subsys_release(oldsys_2(1)%subsys)
    NULLIFY(moves_2,move_updates_2,averages_2,particles_old_2)

! close the coordinate file
  CLOSE(com_crd)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE gemc_nvt
!*********************************************************************
  SUBROUTINE external_control(force_env,should_stop,error)

!   Purpose: External manipulations during an MD run.

!   History: - Creation (16.09.2003,MK)

!   ***************************************************************************

    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(OUT)                     :: should_stop
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: exit_file_name = "EXIT_MC", &
      routine_name = "external_control", &
      routineP = module_name//':'//routine_name

    INTEGER                                  :: handle, unit_number
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    should_stop = .FALSE.

    IF (logger%para_env%source==logger%para_env%mepos) THEN
       INQUIRE (FILE=exit_file_name,EXIST=should_stop)

       IF (should_stop) THEN
          CALL open_file(file_name=exit_file_name,unit_number=unit_number)
          CALL close_file(unit_number=unit_number,file_status="DELETE")
          WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T2,A,/)")&
               "*** MC run terminated by external request ***"
       END IF
    END IF
    CALL mp_bcast(should_stop,logger%para_env%source,logger%para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE external_control

END MODULE mc_ensembles

!**********************************************************************
