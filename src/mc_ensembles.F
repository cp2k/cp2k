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
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_to_string
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
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
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

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_traditional_ensemble"

    INTEGER                                  :: crd, ene, handle, i, &
                                                iparticle, iw, nnstep
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
   CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
   OPEN(UNIT=ene)
      WRITE(ene,*) 'Initial Energy:                        ', old_energy
   CLOSE(ene)
   OPEN(UNIT=68)
      WRITE(68,*) 'Initial Energy:                        ', old_energy
   CLOSE(68)
   initial_energy=old_energy
   CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','INITIAL')

! record the old coordinates, in case a move is rejected
   DO iparticle=1,mc_par%nunits_tot
      r_old(1:3,iparticle)=particles_old%els(iparticle)%r(1:3)
   ENDDO
!   CALL cp_subsys_release(oldsys(1)%subsys)

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
                        nnstep,averages,old_energy,globenv,&
                        blacs_env,energy_check,r_old,&
                        error=error)
      ELSE
         DO i=1,mc_par%nmoves
            CALL random_number(rand)
            IF(rand .lt. mc_par % pmtraion) THEN
! translate a single atom in the system
               CALL mc_atom_translation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,&
                        classical_energy_old,error=error)
            ELSEIF(rand .lt. mc_par % pmtrans) THEN
! translate a whole molecule in the system
               CALL mc_molecule_translation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,&
                        classical_energy_old,error=error)
            ELSE
! rotate a whole molecule in the system
               CALL mc_molecule_rotation(  mc_par,force_env, moves,&
                        move_updates,nnstep,i,&
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

! deallocate some stuff
    CALL force_env_release(force_env,error=error)
    CALL cp_subsys_release(oldsys(1)%subsys)
    NULLIFY(moves,move_updates,averages,particles_old)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_traditional_ensemble
!*********************************************************************
!**********************************************************************
  SUBROUTINE gemc_nvt ( mc_env, globenv, error )

!   Purpose: Gibbs canonical ensemble Monte Carlo run using QUICKSTEP.

!   History: - Creation (11.14.2003)

!   ***************************************************************************

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = "gemc_nvt"

    CHARACTER(LEN=200), DIMENSION(1:2)       :: input_file_names
    INTEGER                                  :: box_flag, crd = 75, handle, &
                                                isos, nnstep, stat
    INTEGER, DIMENSION(1:2)                  :: istep
    REAL(KIND=dp)                            :: rand
    REAL(KIND=dp), DIMENSION(1:2)            :: angle, ave_density, &
                                                ave_energy, bond1, bond2, &
                                                old_energy
    TYPE(cell_type), POINTER                 :: cell_1, cell_2
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys_1, oldsys_2
    TYPE(force_env_type), POINTER            :: force_env_1, force_env_2
    TYPE(mc_moves_type), POINTER             :: move_updates_1, &
                                                move_updates_2, moves_1, &
                                                moves_2
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par_1, mc_par_2
    TYPE(particle_list_type), POINTER        :: particles_old_1, &
                                                particles_old_2

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
!            CALL mc_atom_translation( mc_par_1,force_env_1, moves_1,&
!                        move_updates_1,energy_check,&
!                        istep(1),ave_energy(1),old_energy(1),globenv,&
!                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
!            CALL mc_atom_translation( mc_par_2,force_env_2, moves_2,&
!                        move_updates_2,energy_check,&
!                        istep(2),ave_energy(2),old_energy(2),globenv,&
!                        error=error)
            box_flag=2
         ENDIF
      ELSEIF(rand .lt. mc_par_1 % pmtrans) THEN
! translate a whole molecule in the system
         CALL random_number(rand)
         IF (rand .lt. 0.5d0) THEN
! move in box 1
            istep(1)=istep(1)+1
            globenv%input_file_name=input_file_names(1)
!            CALL mc_molecule_translation( mc_par_1,force_env_1,&
!                        moves_1,move_updates_1,energy_check,&
!                        istep(1),ave_energy(1),old_energy(1),globenv,&
!                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
!            CALL mc_molecule_translation( mc_par_2,force_env_2,&
!                        moves_2,move_updates_2,energy_check,&
!                        istep(2),ave_energy(2),old_energy(2),globenv,&
!                        error=error)
            box_flag=2
         ENDIF
      ELSE
! rotate a whole molecule in the system
         CALL random_number(rand)
         IF (rand .lt. 0.5d0) THEN
! move in box 1
            istep(1)=istep(1)+1
            globenv%input_file_name=input_file_names(1)
!            CALL mc_molecule_rotation( mc_par_1,force_env_1,& 
!                        moves_1,move_updates_1,energy_check,&
!                        istep(1),ave_energy(1),old_energy(1),globenv,&
!                        error=error)
            box_flag=1
         ELSE
! move in box 2
            istep(2)=istep(2)+1
            globenv%input_file_name=input_file_names(2)
!            CALL mc_molecule_rotation( mc_par_2,force_env_2,&
!                        moves_2,move_updates_2,energy_check,&
!                        istep(2),ave_energy(2),old_energy(2),globenv,&
!                        error=error)
            box_flag=2
         ENDIF
      ENDIF

! find bond lengths and angle for current water config for the box
! that was changed
      SELECT CASE (box_flag)
      CASE(1)
!         CALL get_water_data( force_env_1, angle(1), bond1(1), bond2(1), &
!                              istep(1), mc_par_1%nchain,mc_par_1%data)
      CASE(2)
!         CALL get_water_data( force_env_2, angle(2), bond1(2), bond2(2), &
!                              istep(2), mc_par_2%nchain,mc_par_2%data)
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
!   CALL final_mc_write(mc_par_1,ave_energy(1),bond1(1),bond2(1),&
!                         angle(1),moves_1,ave_density(1))
   WRITE(6,*) '------------------------------------------------'
   WRITE(6,*) '|                   BOX 2                      |'
   WRITE(6,*) '------------------------------------------------'
!   CALL final_mc_write(mc_par_2,ave_energy(2),bond1(2),bond2(2),&
!                         angle(2),moves_2,ave_density(2))
 
! deallocate some stuff
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
