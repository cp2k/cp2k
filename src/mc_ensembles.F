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
!!     Contains the details of the individuals ensembles used in an MC run
!!
!!   AUTHOR
!!     Matthew J. McGrath  (09.26.2003)
!!
!******************************************************************************

#include "cp_prep_globals.h"

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
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calculate_energy
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_retain,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE mc_control,                      ONLY: write_mc_restart
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_environment_type
  USE mc_ge_moves,                     ONLY: mc_ge_quickstep_move,&
                                             mc_ge_swap_move,&
                                             mc_ge_volume_move
  USE mc_misc,                         ONLY: final_mc_write,&
                                             get_water_data,&
                                             init_mc_averages,&
                                             mc_classical_energy
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_move_update,&
                                             mc_moves_release,&
                                             write_move_stats
  USE mc_moves,                        ONLY: mc_Quickstep_move,&
                                             mc_atom_translation,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_types,                        ONLY: mc_averages_p_type,&
                                             mc_averages_type,&
                                             mc_moves_p_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_p_type,&
                                             mc_simulation_parameters_type
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
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

! *****************************************************************************
!!****s* mc_ensembles.F/mc_traditional_ensemble
!!
!!   NAME
!!     mc_traditional_ensemble
!!
!!   FUNCTION
!!     directs the program in running a one box traditional MC
!!     simulation (NVT and NpT ensembles)
!!
!!   ARGUMENTS
!!     - mc_env: the mc_environment for the simulation box
!!     - globenv: the global environment for the simulation
!!     - error: the error type
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE mc_traditional_ensemble ( mc_env, globenv, error )

!   ***************************************************************************

    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_traditional_ensemble"

    INTEGER                                  :: crd, ene, handle, i, &
                                                iparticle, iw, nnstep,&
                                                mol,data,istat,j,&
                                                move_unit,cell_unit
    LOGICAL                               :: PM0, PM1, PM2, PM3, PM4, &
                                              should_stop,ionode,exists,&
                                              should_retain
    REAL(KIND=dp)                    :: energy_check, initial_energy, &
                                         old_energy, rand,final_energy,&
                                         classical_energy,&
                                         classical_energy_new
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: particle_positions,r_old
    TYPE(mc_averages_type), POINTER          :: averages
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(mc_moves_type), POINTER             :: move_updates, moves
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(particle_list_type), POINTER        :: particles_old

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      NULLIFY(moves,move_updates,oldsys)
      NULLIFY(particles_old)

! make it so we don't have to type mc_env for every structure    
      CALL get_mc_env ( mc_env, mc_par = mc_par,force_env=force_env)

! allocate the particle positions array for broadcasting
      ALLOCATE (particle_positions(3,mc_par%nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions",3*mc_par%nunits_tot*dp_size)
      ALLOCATE (r_old(3,mc_par%nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old",3*mc_par%nunits_tot*dp_size)

! deal with the different levels of print statements
      PM0 = .FALSE.
      PM1 = .FALSE.
      PM2 = .FALSE.
      PM3 = .FALSE.
      PM4 = .FALSE.
      CALL write_checkpoint_information("entering "//routineN,globenv)
      iw = globenv%scr

      ionode = mc_par%ionode
      IF((mc_par%ionode).AND.(globenv%print%level .GE. 0)) PM0= .TRUE.
      IF((mc_par%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
      IF((mc_par%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
      IF((mc_par%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
      IF((mc_par%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
      

      IF (PM0) THEN
         WRITE(iw,*)
         WRITE(iw,*)
         WRITE(iw,*) 'Beginning Monte Carlo calculation',&
                ' in a traditional ensemble'
         WRITE(iw,*)
         WRITE(iw,*)
      ENDIF

! change the code to serial mode
      IF (mc_par % ionode) THEN
! initialize running average variables
         energy_check=0.0d0

! initialize the moves array
         CALL init_mc_moves(moves)

! initialize the arrays for updating maximum move displacements
         CALL init_mc_moves(move_updates)

! initialize the averages array
         CALL init_mc_averages(averages)
! back to parallel mode
      ENDIF

! find the energy of the initial configuration
!!!!!!!!!!
!      CALL force_env_calculate_energy(force_env)

      CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
         potential_energy=old_energy,error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
         particles=particles_old, error=error)

! back to seriel
      IF (mc_par % ionode) THEN

         CALL mc_classical_energy(force_env,moves%classical_energy,&
                               mc_par,1,mc_par%nchain,6.211D0/angstrom)
!!!!!!!!!!!
         old_energy=moves%classical_energy
         classical_energy=moves%classical_energy

! open the files for writing coordinates, average structural data, and
! maxiumum displacements
         CALL open_file(file_name=mc_par%data_file,unit_number=DATA,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=mc_par%coords_file,unit_number=crd,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=mc_par%energy_file,unit_number=ene,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=mc_par%moves_file,&
            unit_number=move_unit,file_action='WRITE',&
            file_status='REPLACE')
         CALL open_file(file_name=mc_par%displacement_file,&
                 unit_number=mc_par%rm,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='cell_length',&
                 unit_number=cell_unit,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(cell_unit,*) 'Initial Cell Length:       ',&
               cell%hmat(1,1)*angstrom
         CALL open_file(file_name=mc_par%cell_file,&
                 unit_number=mc_par%cl,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='energy_differences',&
                 unit_number=mc_par%diff,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene,*) 'Initial Energy:                ',&
            old_energy,'0.0D0'
         CALL write_particle_coordinates(particles_old%els,cell,&
                   crd,'XMOL','POS','INITIAL')

         initial_energy=old_energy

! record the old coordinates, in case a move is rejected
         DO iparticle=1,mc_par%nunits_tot
            r_old(1:3,iparticle)=particles_old%els(iparticle)%r(1:3)
         ENDDO

! back to parallel mode
      ENDIF

! find out how many steps we're doing
      IF (mc_par%lstop .EQ. 'no') THEN
         mc_par%nstep=mc_par%nstep*mc_par%nchain
      ENDIF

      DO nnstep=mc_par % nstart+1,mc_par % nstart + mc_par % nstep

         CALL RANDOM_NUMBER(rand)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,mc_par%source,mc_par%group)

         IF (rand .LT. mc_par % pmvolume) THEN

! some memory stuff
#if defined(__HMD)
      IF (mc_par%ionode) CALL _dump_allocated(%VAL(0))
#endif
            CALL mc_volume_move( mc_par,force_env,moves,move_updates,&
                        nnstep,old_energy,globenv,&
                        energy_check,r_old,&
                        error=error)
#if defined(__HMD)
      IF (mc_par%ionode) CALL _dump_allocated(%VAL(0))
#endif

        ELSE
            particle_positions(:,:)=0.0d0
! we want to do all this in serial
            IF( mc_par % ionode ) THEN

               DO i=1,mc_par%nmoves
                  
                  CALL RANDOM_NUMBER(rand)
                  IF(rand .LT. mc_par % pmtraion) THEN
! change molecular conformation
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(rand*mc_par%nchain)
                     mol=3*mol-2
                     CALL mc_atom_translation( mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ELSEIF(rand .LT. mc_par % pmtrans) THEN
! translate a whole molecule in the system
! call a random number to figure out which molecule we're moving...the
! second line is necessary because we access a matrix with individual
! particles
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(mc_par%nchain*rand)
                     mol=3*mol-2
                     CALL mc_molecule_translation(mc_par,force_env,&
                      moves,move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ELSE
! rotate a whole molecule in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(mc_par%nchain*rand)
                     mol=3*mol-2
                     CALL mc_molecule_rotation(mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ENDIF

               ENDDO

! get the new coordinates and dump them into an array so we can broadcast it
               DO iparticle=1,mc_par%nunits_tot
                  particle_positions(1:3,iparticle)=&
                      particles_old%els(iparticle)%r(1:3)
               ENDDO


! back to parallel
            ENDIF

! broadcast all the positions, to make sure we're on the same page
            CALL mp_bcast(particle_positions,mc_par%source,mc_par%group)

            DO iparticle=1,mc_par%nunits_tot
               particles_old%els(iparticle)%r(1:3)= &
                 particle_positions(1:3,iparticle)
            ENDDO

! now do a Quickstep calculation to see if we accept the sequence
            CALL mc_Quickstep_move(mc_par,force_env,moves,&
                        move_updates,energy_check,r_old,&
                        nnstep,old_energy,error=error)

         ENDIF

! make sure the pointers are pointing correctly since the subsys may
!    have changed
         CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
            error=error)
         CALL cp_subsys_get(oldsys(1)%subsys,&
            particles=particles_old, error=error)


         IF(ionode) THEN
!            WRITE(ene,*) nnstep,old_energy/3.16682968D-6,energy_check
            WRITE(ene,*) nnstep,old_energy,energy_check
!!!!!!!!!!
!            CALL write_particle_coordinates(particles_old%els,cell,&
!                        crd,'XMOL','POS','POSITIONS')

! compute some averages
            averages%molecules=averages%molecules*&
                     REAL(nnstep-mc_par%nstart-1,dbl)&
                     /REAL(nnstep-mc_par%nstart,dbl)+mc_par%nchain/&
                     REAL(nnstep-mc_par%nstart,dbl)

            averages%ave_energy=averages%ave_energy*&
                     REAL(nnstep-mc_par%nstart-1,dbl)&
                     /REAL(nnstep-mc_par%nstart,dbl)+old_energy/&
                     REAL(nnstep-mc_par%nstart,dbl)

            averages%ave_cell_length=averages%ave_cell_length*&
                     REAL(nnstep-mc_par%nstart-1,dbl)&
                     /REAL(nnstep-mc_par%nstart,dbl)+cell%hmat(1,1)/&
                     REAL(nnstep-mc_par%nstart,dbl)

            WRITE(cell_unit,*) nnstep,cell%hmat(1,1)*angstrom

! find bond lengths and angle for current water config
            CALL get_water_data(force_env,averages,&
               nnstep-mc_par%nstart,mc_par%nchain,DATA)

! write the move statistics to file
            CALL write_move_stats(moves,move_updates,nnstep,move_unit)

! write a restart file
            CALL write_mc_restart(nnstep,mc_par,force_env)

! flush the buffers to the files
            CALL m_flush(DATA)
            CALL m_flush(mc_par%diff)
            CALL m_flush(ene)
            CALL m_flush(move_unit)
            CALL m_flush(crd)
            CALL m_flush(mc_par%cl)
            CALL m_flush(mc_par%rm)
            CALL m_flush(cell_unit)

! check to see if EXIT file exists...if so, end the calculation
            CALL external_control(force_env,should_stop,error=error)

! update the move displacements, if necessary
            IF (MOD(nnstep-mc_par%nstart,mc_par%iuptrans) == 0) THEN
               CALL mc_move_update(mc_par,move_updates,"trans",nnstep)
            ENDIF
            IF (MOD(nnstep-mc_par%nstart,mc_par%iupvolume) == 0  ) THEN
               CALL mc_move_update(mc_par,move_updates,"volume",nnstep)
            ENDIF

         ENDIF

         CALL mp_bcast(should_stop,mc_par%source,mc_par%group)

         IF (should_stop) EXIT

      ENDDO


! calculate the final energy
!!!!!!!!
!      CALL force_env_calculate_energy(force_env)

! get the energy
      CALL force_env_get(force_env,potential_energy=final_energy,&
         error=error)
   
! do some stuff in serial
      IF (mc_par % ionode) THEN

         CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','FINAL')
                        
!!!!!!!!!!!
         CALL mc_classical_energy(force_env,final_energy,&
                               mc_par,1,mc_par%nchain,6.211D0/angstrom)

         WRITE(ene,*) 'Final Energy:                        ',&
                  final_energy,energy_check

! write a bunch of data to the screen
         CALL final_mc_write(mc_par,moves,iw,energy_check,&
                         initial_energy,final_energy,averages,&
                         move_updates)

! close any open files
         CALL close_file(unit_number=mc_par%diff)
         CALL close_file(unit_number=crd)
         CALL close_file(unit_number=DATA)
         CALL close_file(unit_number=move_unit)
         CALL close_file(unit_number=ene)
         CALL close_file(unit_number=cell_unit)
         CALL close_file(unit_number=mc_par%cl)
         CALL close_file(unit_number=mc_par%rm)

! back to parallel
      ENDIF

      should_retain=.FALSE.
      IF (mc_par%ionode)THEN
         IF( moves%volume%successes .GT. 0) &
               should_retain=.TRUE.
      ENDIF
      CALL mp_bcast(should_retain,mc_par%source,mc_par%group)

      IF (should_retain) THEN
         CALL force_env_retain(force_env)
      ENDIF

! deallocate some stuff
      CALL cp_subsys_release(oldsys(1)%subsys)
      NULLIFY(moves,move_updates,averages,particles_old)

      DEALLOCATE(particle_positions,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions")
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old")

! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_traditional_ensemble
!*********************************************************************

! *****************************************************************************
!!****s* mc_ensembles.F/gemc_nvt
!!
!!   NAME
!!     gemc_nvt
!!
!!   FUNCTION
!!     directs the program in running a two box Gibbs ensemble MC
!!     simulation (NVT ensemble)
!!
!!   ARGUMENTS
!!     - mc_env_1: the mc_environment for the first simulation box
!!     - mc_env_2: the mc_environment for the second simulation box
!!     - globenv: the global environment for the simulation
!!     - error: the error type
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE gemc_nvt ( mc_env_1, mc_env_2, globenv, error )

!   ***************************************************************************

    TYPE(mc_environment_type), POINTER       :: mc_env_1, mc_env_2
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routineN = "gemc_nvt"

    INTEGER                                  :: com_crd, com_ene, com_mol, &
                                                handle, i, ibox, iparticle, &
                                                istat, iw, mol, nnstep,&
                                                cell_unit
    INTEGER, DIMENSION(1:2)                  :: box_flag, crd, DATA, ene, &
                                                move_unit
    LOGICAL                                  :: ionode, PM0, PM1, PM2, PM3, &
                                                PM4, should_retain, &
                                                should_stop
    REAL(KIND=dp)                            :: rand
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: particle_positions, r_old
    REAL(KIND=dp), DIMENSION(1:2)            :: classical_energy, &
                                                energy_check, final_energy, &
                                                initial_energy, old_energy
    TYPE(cell_type), POINTER                 :: cell_1, cell_2
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys_1, oldsys_2
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    TYPE(mc_averages_p_type), DIMENSION(:), &
      POINTER                                :: averages
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: move_updates, moves
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_type), POINTER        :: particles_old_1, &
                                                particles_old_2

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

      NULLIFY(particles_old_1)
      NULLIFY(particles_old_2)

! make it so we don't have to type mc_env for every structure    
!      CALL get_mc_env ( mc_env_1,mc_par=mc_par_1,force_env=force_env_1)
!      CALL get_mc_env ( mc_env_2,mc_par=mc_par_2,force_env=force_env_2)

! allocate the particle positions array for broadcasting
      ALLOCATE(averages(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                   "averages",2)
      ALLOCATE(force_env(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                   "force_env",2)
      ALLOCATE(mc_par(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                   "mc_par",2)
      ALLOCATE(moves(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                   "moves",2)
      ALLOCATE(move_updates(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                   "moves_updates",2)

      CALL get_mc_env ( mc_env_1,mc_par=mc_par(1)%mc_par,&
         force_env=force_env(1)%force_env)
      CALL get_mc_env ( mc_env_2,mc_par=mc_par(2)%mc_par,&
         force_env=force_env(2)%force_env)
  
      ALLOCATE(particle_positions(3,mc_par(1)%mc_par%nunits_tot+&
              mc_par(2)%mc_par%nunits_tot,2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
           "particle_positions",3*(mc_par(1)%mc_par%nunits_tot+&
           mc_par(2)%mc_par%nunits_tot)*dp_size*2)
      ALLOCATE(r_old(3,mc_par(1)%mc_par%nunits_tot+&
              mc_par(2)%mc_par%nunits_tot,2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
           "r_old",3*(mc_par(1)%mc_par%nunits_tot+&
           mc_par(2)%mc_par%nunits_tot)*dp_size*2)

! deal with the different levels of print statements
      PM0 = .FALSE.
      PM1 = .FALSE.
      PM2 = .FALSE.
      PM3 = .FALSE.
      PM4 = .FALSE.
      CALL write_checkpoint_information("entering "//routineN,globenv)
      iw = globenv%scr

      ionode = mc_par(1)%mc_par%ionode
      IF((ionode).AND.(globenv%print%level .GE. 0)) PM0= .TRUE.
      IF((ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
      IF((ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
      IF((ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
      IF((ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
      

      IF (PM0) THEN
         WRITE(iw,*)
         WRITE(iw,*)
         WRITE(iw,*) 'Beginning Monte Carlo calculation',&
                ' in the Gibbs ensemble'
         WRITE(iw,*)
         WRITE(iw,*)
      ENDIF

! initialize running average variables
      energy_check(:)=0.0d0
      box_flag(:)=0

! initialize the moves array,the arrays for updating maximum move 
!   displacements, and the averages array
      DO i=1,2
         CALL init_mc_moves(moves(i)%moves)
         CALL init_mc_moves(move_updates(i)%moves)
         CALL init_mc_averages(averages(i)%averages)
      ENDDO

! find the energy of the initial configuration
      DO i=1,2
!!!!!!!!!!!!
!         CALL force_env_calculate_energy(force_env(i)%force_env)
         CALL mc_classical_energy(force_env(i)%force_env,&
            moves(i)%moves%classical_energy,mc_par(i)%mc_par,1,&
            mc_par(i)%mc_par%nchain,6.211D0/angstrom)
         old_energy(i)=moves(i)%moves%classical_energy
         classical_energy(i)=moves(i)%moves%classical_energy


      ENDDO

!!!!!!!!!!!!!
!      CALL force_env_get(force_env(1)%force_env,subsys=oldsys_1,&
!         cell=cell_1,&
!         potential_energy=old_energy(1),error=error)
      CALL force_env_get(force_env(1)%force_env,subsys=oldsys_1,&
         cell=cell_1,&
         error=error)
      CALL cp_subsys_get(oldsys_1(1)%subsys, &
         particles=particles_old_1, error=error)
      
!      CALL force_env_get(force_env(2)%force_env,subsys=oldsys_2,&
!         cell=cell_2,&
!         potential_energy=old_energy(2),error=error)
      CALL force_env_get(force_env(2)%force_env,subsys=oldsys_2,&
         cell=cell_2,&
         error=error)
      CALL cp_subsys_get(oldsys_2(1)%subsys, &
         particles=particles_old_2, error=error)

! back to seriel
      IF (mc_par(1)%mc_par % ionode) THEN

! open the files for writing coordinates, average structural data, and
! maxiumum displacements

! box 1
         DO i=1,2
            CALL open_file(file_name=mc_par(i)%mc_par%data_file,&
               unit_number=DATA(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=mc_par(i)%mc_par%moves_file,&
               unit_number=move_unit(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=mc_par(i)%mc_par%coords_file,&
               unit_number=crd(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=mc_par(i)%mc_par%energy_file,&
               unit_number=ene(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=mc_par(i)%mc_par%displacement_file,&
               unit_number=mc_par(i)%mc_par%rm,&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=mc_par(i)%mc_par%cell_file,&
                  unit_number=mc_par(i)%mc_par%cl,&
                  file_action='WRITE',file_status='REPLACE')
         ENDDO
            
         CALL open_file(file_name='energy_differences_box1',&
                 unit_number=mc_par(1)%mc_par%diff,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene(1),*) 'Initial Energy Box 1:            ',&
            old_energy(1),' 0.0D0'
         CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                   crd(1),'XMOL','POS','INITIAL')

! box 2
         CALL open_file(file_name='energy_differences_box2',&
                 unit_number=mc_par(2)%mc_par%diff,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene(2),*) 'Initial Energy Box 2:            ',&
            old_energy(2),' 0.0D0'
         CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                   crd(2),'XMOL','POS','INITIAL')

! record the combined energies,coordinates, and cell lengths
         CALL open_file(file_name='cell_length',&
                 unit_number=cell_unit,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(cell_unit,*) 'Initial Cell Lengths: ',&
            cell_1%hmat(1,1)*angstrom,cell_2%hmat(1,1)*angstrom
         CALL open_file(file_name='combined_energies',&
                 unit_number=com_ene,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='combined_coordinates',&
                 unit_number=com_crd,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=mc_par(1)%mc_par%molecules_file,&
                 unit_number=com_mol,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(com_ene,*) 'Initial Total Energy:            ',&
            old_energy(1)+old_energy(2),' 0.0D0'
         CALL write_particle_coordinates(particles_old_1%els,cell_1,&
                   com_crd,'XMOL','POS','INITIAL BOX 1')
         CALL write_particle_coordinates(particles_old_2%els,cell_2,&
                   com_crd,'XMOL','POS','INITIAL BOX 2')

         initial_energy(1)=old_energy(1)
         initial_energy(2)=old_energy(2)

! record the old coordinates, in case a move is rejected
         DO iparticle=1,mc_par(1)%mc_par%nunits_tot
            r_old(1:3,iparticle,1)=particles_old_1%els(iparticle)%r(1:3)
         ENDDO
         DO iparticle=1,mc_par(2)%mc_par%nunits_tot
            r_old(1:3,iparticle,2)=particles_old_2%els(iparticle)%r(1:3)
         ENDDO

! back to parallel mode
      ENDIF

! find out how many steps we're doing
      IF (mc_par(1)%mc_par%lstop .EQ. 'no') &
         mc_par(1)%mc_par%nstep=mc_par(1)%mc_par%nstep*&
            (mc_par(1)%mc_par%nchain+mc_par(2)%mc_par%nchain)

!    MC STUFF

      DO nnstep=mc_par(1)%mc_par%nstart+1,mc_par(1)%mc_par%nstart + &
                                    mc_par(1)%mc_par%nstep

         CALL RANDOM_NUMBER(rand)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,mc_par(1)%mc_par%source,&
            mc_par(1)%mc_par%group)


         IF (rand .LT. mc_par(1)%mc_par % pmvolume) THEN
!     change the volume of the system

            CALL mc_ge_volume_move ( mc_par,force_env, moves,&
               move_updates,nnstep,old_energy,globenv,energy_check,&
               r_old,error)

         ELSEIF ( rand .LT. mc_par(1)%mc_par % pmswap) THEN

! try a swap move
            CALL mc_ge_swap_move(mc_par,force_env,moves,&
                energy_check,r_old,old_energy,&
                nnstep,globenv,error)

         ELSE

            particle_positions(:,:,:) = 0.0D0
            
! go into serial mode
            IF (mc_par(1)%mc_par%ionode) THEN

               DO i=1,mc_par(1)%mc_par%nmoves
                  CALL RANDOM_NUMBER(rand)
                  IF(rand .LT. mc_par(1)%mc_par % pmtraion) THEN
                     
!     translate a single atom in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((mc_par(1)%mc_par%nchain+&
                        mc_par(2)%mc_par%nchain)*rand)
                     mol=3*mol-2
                     IF (mol .LE. mc_par(1)%mc_par%nunits_tot) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-mc_par(1)%mc_par%nunits_tot
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_atom_translation( mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,i,mol,&
                        classical_energy(ibox),error=error)
                        
                  ELSEIF(rand .LT. mc_par(1)%mc_par % pmtrans) THEN
!     translate a whole molecule in the system
!     call a random number to figure out which molecule we're moving...the
!     second line is necessary because we access a matrix with individual
!     particles
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((mc_par(1)%mc_par%nchain+&
                        mc_par(2)%mc_par%nchain)*rand)
                     mol=3*mol-2
                     IF (mol .LE. mc_par(1)%mc_par%nunits_tot) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-mc_par(1)%mc_par%nunits_tot
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_molecule_translation(mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,&
                        i,mol,classical_energy(ibox),error=error)

                  ELSE
!     rotate a whole molecule in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((mc_par(1)%mc_par%nchain+&
                        mc_par(2)%mc_par%nchain)*rand)
                     mol=3*mol-2
                     IF (mol .LE. mc_par(1)%mc_par%nunits_tot) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-mc_par(1)%mc_par%nunits_tot
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_molecule_rotation(mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,i,mol,&
                        classical_energy(ibox),error=error)
                     
                  ENDIF
               ENDDO
               
! get the new coordinates and dump them into an array so we can broadcast it
               DO iparticle=1,mc_par(1)%mc_par%nunits_tot
                  particle_positions(1:3,iparticle,1)=&
                  particles_old_1%els(iparticle)%r(1:3)
               ENDDO
               DO iparticle=1,mc_par(2)%mc_par%nunits_tot
                  particle_positions(1:3,iparticle,2)=&
                     particles_old_2%els(iparticle)%r(1:3)
               ENDDO
               
!     go back to parallel      
            ENDIF
            
!     broadcast all the positions, to make sure we're on the same page
            CALL mp_bcast(particle_positions,mc_par(1)%mc_par%source,&
               mc_par(1)%mc_par%group)
            
            DO iparticle=1,mc_par(1)%mc_par%nunits_tot
               particles_old_1%els(iparticle)%r(1:3)= &
                  particle_positions(1:3,iparticle,1)
            ENDDO
            DO iparticle=1,mc_par(2)%mc_par%nunits_tot
               particles_old_2%els(iparticle)%r(1:3)= &
                  particle_positions(1:3,iparticle,2)
            ENDDO
            
!     now do a Quickstep calculation to see if we accept the sequence
            CALL mc_ge_quickstep_move(mc_par(1)%mc_par,&
               mc_par(2)%mc_par,force_env(1)%force_env,&
               force_env(2)%force_env,moves(1)%moves,moves(2)%moves,&
               move_updates(1)%moves,move_updates(2)%moves,&
               energy_check(:),r_old(:,:,1),r_old(:,:,2),&
               old_energy(:),box_flag(:),nnstep,error=error)
            
         ENDIF
         
! make sure the pointers are pointing correctly since the subsys may
!    have changed
         CALL force_env_get(force_env(1)%force_env,cell=cell_1,&
            subsys=oldsys_1,error=error)
         CALL cp_subsys_get(oldsys_1(1)%subsys,&
            particles=particles_old_1, error=error)
         CALL force_env_get(force_env(2)%force_env,cell=cell_2,&
            subsys=oldsys_2,error=error)
         CALL cp_subsys_get(oldsys_2(1)%subsys, &
            particles=particles_old_2, error=error)


!     find bond lengths, angle, energy, and number of molecules
!     for current water config 
         IF (mc_par(1)%mc_par%ionode) THEN
            DO ibox=1,2
               CALL write_move_stats(moves(ibox)%moves,&
                  move_updates(ibox)%moves,nnstep,&
                  move_unit(ibox))
               CALL get_water_data( force_env(ibox)%force_env,&
                  averages(ibox)%averages,nnstep,&
                  mc_par(ibox)%mc_par%nchain,DATA(ibox))
               WRITE(ene(ibox),*) nnstep,old_energy(ibox)
               averages(ibox)%averages%ave_energy=&
                  averages(ibox)%averages%ave_energy*REAL(nnstep-&
                  mc_par(1)%mc_par%nstart-1,dbl)/&
                  REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)+&
                  old_energy(ibox)/REAL(nnstep,dbl)
               averages(ibox)%averages%molecules=&
                  averages(ibox)%averages%molecules*REAL(nnstep-&
                  mc_par(1)%mc_par%nstart-1,dbl)/&
                  REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)+&
                  REAL(mc_par(ibox)%mc_par%nchain,dbl)/REAL(nnstep-&
                  mc_par(1)%mc_par%nstart,dbl)
            ENDDO
            WRITE(com_ene,*) nnstep,old_energy(2)+old_energy(1)
            IF (mc_par(1)%mc_par%pmswap .GT. 0.0D0) &
               WRITE(com_mol,*) nnstep,mc_par(1)%mc_par%nchain,&
               mc_par(2)%mc_par%nchain

! compute the average cell length, energy, and number of particles of
! both boxes, as well as the total average energy...write to files
            averages(1)%averages%ave_cell_length=&
               averages(1)%averages%ave_cell_length*&
               REAL(nnstep-mc_par(1)%mc_par%nstart-1,dbl)/&
               REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)+&
               cell_1%hmat(1,1)/REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)
            averages(2)%averages%ave_cell_length=&
               averages(2)%averages%ave_cell_length*&
               REAL(nnstep-mc_par(1)%mc_par%nstart-1,dbl)/&
               REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)+&
               cell_2%hmat(1,1)/REAL(nnstep-mc_par(1)%mc_par%nstart,dbl)

            WRITE(cell_unit,*) nnstep,cell_1%hmat(1,1)*angstrom,&
               cell_2%hmat(1,1)*angstrom
            
            

!     write the coordinates and energy to a file
            IF (mc_par(1)%mc_par%nchain == 0) THEN
               WRITE(com_crd,*) ' 0 '
               WRITE(com_crd,*) ' POSITIONS BOX 1 '
               WRITE(crd(1),*) ' 0 '
               WRITE(crd(1),*) ' POSITIONS BOX 1 '
            ELSE
!!!!!!!!!!!
!               CALL write_particle_coordinates(&
!                  particles_old_1%els,cell_1,com_crd,'XMOL',&
!                  'POS','POSITIONS BOX 1')
!               CALL write_particle_coordinates(&
!                  particles_old_1%els,cell_1,crd(1),'XMOL',&
!                  'POS','POSITIONS BOX 1')
            ENDIF
            IF (mc_par(2)%mc_par%nchain == 0) THEN
               WRITE(com_crd,*) ' 0 '
               WRITE(com_crd,*) ' POSITIONS BOX 2 '
               WRITE(crd(2),*) ' 0 '
               WRITE(crd(2),*) ' POSITIONS BOX 2 '
            ELSE
!               CALL write_particle_coordinates(&
!                  particles_old_2%els,cell_2,com_crd,'XMOL',&
!                  'POS','POSITIONS BOX 2')
!               CALL write_particle_coordinates(&
!                  particles_old_2%els,cell_2,crd(2),'XMOL',&
!                  'POS','POSITIONS BOX 2')
            ENDIF
           
            
!     check to see if EXIT file exists...if so, end the calculation
            CALL external_control(force_env(1)%force_env,&
               should_stop,error=error)
            
!     write a restart file
            CALL write_mc_restart(nnstep,mc_par(1)%mc_par,&
               force_env(1)%force_env)
            CALL write_mc_restart(nnstep,mc_par(2)%mc_par,&
               force_env(2)%force_env)
            
! flush the buffers to the files
            DO i=1,2
               CALL m_flush(DATA(i))
               CALL m_flush(move_unit(i))
               CALL m_flush(mc_par(i)%mc_par%diff)
               CALL m_flush(ene(i))
               CALL m_flush(crd(i))
               CALL m_flush(mc_par(i)%mc_par%cl)
               CALL m_flush(mc_par(i)%mc_par%rm)
            ENDDO
            CALL m_flush(com_ene)
            CALL m_flush(com_crd)
            CALL m_flush(com_mol)
            CALL m_flush(cell_unit)

!     reset the box flags
            box_flag(1)=0
            box_flag(2)=0
            
! update the move displacements, if necessary
            IF (MOD(nnstep-mc_par(1)%mc_par%nstart,&
                  mc_par(1)%mc_par%iuptrans) == 0) THEN
               CALL mc_move_update(mc_par(1)%mc_par,&
                  move_updates(1)%moves,"trans",nnstep)
               CALL mc_move_update(mc_par(2)%mc_par,&
                  move_updates(2)%moves,"trans",nnstep)
            ENDIF
            IF (MOD(nnstep-mc_par(1)%mc_par%nstart,&
                  mc_par(1)%mc_par%iupvolume) == 0 .AND. &
                  mc_par(1)%mc_par%pmvolume .NE. 0.0D0 ) THEN
               CALL mc_move_update(mc_par(1)%mc_par,&
                  move_updates(1)%moves,"volume",nnstep)
               CALL mc_move_update(mc_par(2)%mc_par,&
                  move_updates(2)%moves,"volume",nnstep)
            ENDIF


!     back to parallel
         ENDIF
         
!     make sure all the procs exit if necessary
         CALL mp_bcast(should_stop,mc_par(1)%mc_par%source,&
            mc_par(1)%mc_par%group)
         IF (should_stop) EXIT
         
      ENDDO
      
      DO i=1,2
         IF (mc_par(i)%mc_par%nchain .NE. 0) THEN
!!!!!!!!!!!!!!!!!
!            CALL force_env_calculate_energy(force_env(i)%force_env)
!            CALL force_env_get(force_env(i)%force_env,&
!               potential_energy=final_energy(i),error=error)
            CALL mc_classical_energy(force_env(i)%force_env,&
               final_energy(i),mc_par(i)%mc_par,1,&
               mc_par(i)%mc_par%nchain,6.211D0/angstrom)
         ELSE
            final_energy(i)=0.0D0
         ENDIF
      ENDDO

      IF (mc_par(1)%mc_par%ionode) THEN
         CALL force_env_get(force_env(1)%force_env,&
            subsys=oldsys_1,cell=cell_1,error=error)
         CALL cp_subsys_get(oldsys_1(1)%subsys,&
            particles=particles_old_1, error=error)
         CALL write_particle_coordinates(particles_old_1%els,cell_1,&
            com_crd,'XMOL','POS','FINAL BOX 1')
         CALL write_particle_coordinates(particles_old_1%els,cell_1,&
            crd(1),'XMOL','POS','FINAL BOX 1')
         CALL force_env_get(force_env(2)%force_env,&
            subsys=oldsys_2,cell=cell_2,error=error)
         CALL cp_subsys_get(oldsys_2(1)%subsys,&
            particles=particles_old_2, error=error)
         CALL write_particle_coordinates(particles_old_2%els,cell_2,&
            com_crd,'XMOL','POS','FINAL BOX 2')
         CALL write_particle_coordinates(particles_old_2%els,cell_2,&
            crd(2),'XMOL','POS','FINAL BOX 2')
         
         WRITE(ene(1),*) 'Final Energy BOX 1:   ',final_energy(1)
         WRITE(ene(2),*) 'Final Energy BOX 2:   ',final_energy(2)

!     record the final total energy
         WRITE(com_ene,*) 'Final total energy:      ',&
            old_energy(2)+old_energy(1),energy_check(1)+energy_check(2)
      
!     write a bunch of data to the screen
         WRITE(iw,*) '------------------------------------------------'
         WRITE(iw,*) '|                   BOX 1                      |'
         WRITE(iw,*) '------------------------------------------------'
         CALL final_mc_write(mc_par(1)%mc_par,moves(1)%moves,iw,energy_check(1),&
            initial_energy(1),final_energy(1),averages(1)%averages,&
            move_updates(1)%moves)
         WRITE(iw,*) '------------------------------------------------'
         WRITE(iw,*) '|                   BOX 2                      |'
         WRITE(iw,*) '------------------------------------------------'
         CALL final_mc_write(mc_par(2)%mc_par,moves(2)%moves,iw,energy_check(2),&
            initial_energy(2),final_energy(2),averages(2)%averages,&
            move_updates(2)%moves)
      

! close any open files
         DO i=1,2
            CALL close_file(unit_number=mc_par(i)%mc_par%diff)
            CALL close_file(unit_number=crd(i))
            CALL close_file(unit_number=DATA(i))
            CALL close_file(unit_number=move_unit(i))
            CALL close_file(unit_number=ene(i))
            CALL close_file(unit_number=mc_par(i)%mc_par%cl)
            CALL close_file(unit_number=mc_par(i)%mc_par%rm)
         ENDDO
         CALL close_file(unit_number=com_ene)
         CALL close_file(unit_number=com_crd)
         CALL close_file(unit_number=com_mol)
         CALL close_file(unit_number=cell_unit)

      ENDIF

! retain the force environments if we changed where they point
      should_retain=.FALSE.
      IF (mc_par(1)%mc_par%ionode)THEN
         IF(moves(1)%moves%swap%successes .GT. 0 .OR. &
               moves(1)%moves%volume%successes .GT. 0) &
               should_retain=.TRUE.
      ENDIF
      CALL mp_bcast(should_retain,mc_par(1)%mc_par%source,&
            mc_par(1)%mc_par%group)
      IF (should_retain) THEN
         CALL force_env_retain(force_env(1)%force_env)
         CALL force_env_retain(force_env(2)%force_env)
      ENDIF

!     deallocate some stuff
!      CALL set_mc_env ( mc_env_1, mc_par = mc_par(1)%mc_par,&
!         force_env=force_env(1)%force_env)
!      CALL set_mc_env ( mc_env_2, mc_par = mc_par(2)%mc_par,&
!         force_env=force_env(2)%force_env)
      CALL cp_subsys_release(oldsys_1(1)%subsys)
      NULLIFY(particles_old_1)
      CALL cp_subsys_release(oldsys_2(1)%subsys)
      NULLIFY(particles_old_2)

      DEALLOCATE(particle_positions,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions")
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old deallocate")
      DEALLOCATE(averages,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "averages deallocate")
      DEALLOCATE(force_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "force_env deallocate")
      DEALLOCATE(mc_par,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mc_par deallocate")

! release some structures
      DO i=1,2
         CALL mc_moves_release(moves(i)%moves)
         CALL mc_moves_release(move_updates(i)%moves)
      ENDDO
      
      DEALLOCATE(moves,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "moves deallocate")
      DEALLOCATE(move_updates,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "move_updates deallocate")

   

!     end the timing
      CALL timestop(0.0d0,handle)
      
  END SUBROUTINE gemc_nvt
!*********************************************************************
  SUBROUTINE external_control(force_env,should_stop,error)

!   Purpose: External manipulations during an MC run.

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
!    CALL mp_bcast(should_stop,logger%para_env%source,logger%para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE external_control

END MODULE mc_ensembles

!**********************************************************************

