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
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE mc_control,                      ONLY: write_mc_restart
  USE mc_coordinates,                  ONLY: get_water_data,&
                                             mc_classical_energy
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_environment_type,&
                                             set_mc_env
  USE mc_ge_moves,                     ONLY: mc_ge_quickstep_move,&
                                             mc_ge_swap_move,&
                                             mc_ge_volume_move
  USE mc_misc,                         ONLY: final_mc_write,&
                                             mc_averages_create,&
                                             mc_averages_release
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_move_update,&
                                             mc_moves_release,&
                                             write_move_stats
  USE mc_moves,                        ONLY: mc_Quickstep_move,&
                                             mc_atom_translation,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_averages_p_type,&
                                             mc_averages_type,&
                                             mc_moves_p_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_p_type,&
                                             mc_simulation_parameters_type,&
                                             set_mc_par
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
                                                move_unit,cell_unit,&
                                                nunits_tot,nchain,&
                                                source,group,diff,cl,rm,&
                                                nstep,nstart,nmoves,&
                                                iuptrans,iupvolume
    LOGICAL                               :: PM0, PM1, PM2, PM3, PM4, &
                                              should_stop,ionode,exists,&
                                              lstop,lclassical_run,ldum,&
                                              loverlap
    REAL(KIND=dp), DIMENSION(1:3)        :: abc
    REAL(KIND=dp)                    :: energy_check, initial_energy, &
                                         old_energy, rand,final_energy,&
                                         classical_energy,pmtraion,&
                                         classical_energy_new,pmvolume,&
                                         pmtrans,cutoff,rdum
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
   CHARACTER(LEN=40)                         :: data_file,moves_file,&
      energy_file,coords_file,displacement_file,cell_file
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      NULLIFY(moves,move_updates,oldsys)
      NULLIFY(particles_old)

! make it so we don't have to type mc_env for every structure    
      CALL get_mc_env ( mc_env, mc_par = mc_par,force_env=force_env)
! also get some data out of mc_par
      CALL get_mc_par ( mc_par, nunits_tot=nunits_tot,nchain=nchain,&
         ionode=ionode,source=source,group=group,data_file=data_file,&
         moves_file=moves_file,energy_file=energy_file,&
         coords_file=coords_file,displacement_file=displacement_file,&
         cell_file=cell_file,lstop=lstop,nstep=nstep,nstart=nstart,&
         pmvolume=pmvolume,nmoves=nmoves,pmtraion=pmtraion,&
         pmtrans=pmtrans,iuptrans=iuptrans,iupvolume=iupvolume,&
         cutoff=cutoff,lclassical_run=lclassical_run)

! allocate the particle positions array for broadcasting
      ALLOCATE (particle_positions(3,nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions",3*nunits_tot*dp_size)
      ALLOCATE (r_old(3,nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old",3*nunits_tot*dp_size)

! deal with the different levels of print statements
      PM0 = .FALSE.
      PM1 = .FALSE.
      PM2 = .FALSE.
      PM3 = .FALSE.
      PM4 = .FALSE.
      CALL write_checkpoint_information("entering "//routineN,globenv)
      iw = globenv%scr

      IF((ionode).AND.(globenv%print%level .GE. 0)) PM0= .TRUE.
      IF((ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
      IF((ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
      IF((ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
      IF((ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
      

      IF (PM0) THEN
         WRITE(iw,*)
         WRITE(iw,*)
         WRITE(iw,*) 'Beginning Monte Carlo calculation',&
                ' in a traditional ensemble'
         WRITE(iw,*)
         WRITE(iw,*)
      ENDIF

! initialize running average variables
      energy_check=0.0d0

! initialize the moves array
      CALL init_mc_moves(moves)

! initialize the arrays for updating maximum move displacements
      CALL init_mc_moves(move_updates)

! initialize the averages array
      CALL mc_averages_create(averages)

! find the energy of the initial configuration
      IF (.NOT. lclassical_run) THEN
         CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.)
         CALL force_env_get(force_env,potential_energy=old_energy,&
            error=error)

         CALL mc_classical_energy(force_env,moves%classical_energy,&
            nchain,1,nchain,cutoff,ldum)
         IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
          'traditional ensemble found an overlap in the initial config')
         classical_energy=moves%classical_energy

      ELSE
         CALL mc_classical_energy(force_env,moves%classical_energy,&
                               nchain,1,nchain,cutoff,ldum)
         IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
          'traditional ensemble found an overlap in the initial config')
         old_energy=moves%classical_energy
         classical_energy=moves%classical_energy

      ENDIF

      CALL force_env_get(force_env,subsys=oldsys,cell=cell,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys, &
         particles=particles_old, error=error)

! back to seriel
      IF (ionode) THEN

! open the files for writing coordinates, average structural data, and
! maxiumum displacements
         CALL open_file(file_name=data_file,unit_number=DATA,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=coords_file,unit_number=crd,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=energy_file,unit_number=ene,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=moves_file,&
            unit_number=move_unit,file_action='WRITE',&
            file_status='REPLACE')
         CALL open_file(file_name=displacement_file,&
                 unit_number=rm,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='cell_length',&
                 unit_number=cell_unit,&
                 file_action='WRITE',file_status='REPLACE')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         WRITE(cell_unit,*) (abc(1)*angstrom)**3,old_energy
         WRITE(cell_unit,*) 'Initial Cell Length:       ',&
               abc(1)*angstrom

         CALL open_file(file_name=cell_file,&
                 unit_number=cl,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='energy_differences',&
                 unit_number=diff,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene,*) 'Initial Energy:                ',&
            old_energy,'0.0D0'
         CALL write_particle_coordinates(particles_old%els,cell,&
                   crd,'XMOL','POS','INITIAL')

         CALL set_mc_par(mc_par,cl=cl,rm=rm,diff=diff)

         initial_energy=old_energy

! record the old coordinates, in case a move is rejected
         DO iparticle=1,nunits_tot
            r_old(1:3,iparticle)=particles_old%els(iparticle)%r(1:3)
         ENDDO

! back to parallel mode
      ENDIF

! find out how many steps we're doing
      IF (.NOT. lstop ) THEN
         nstep=nstep*nchain
      ENDIF

      DO nnstep=nstart+1,nstart + nstep

         IF(PM0) THEN
            WRITE(iw,*) 
            WRITE(iw,*) "------- On Monte Carlo Step ",nnstep
         ENDIF

         CALL RANDOM_NUMBER(rand)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,source,group)

         IF (rand .LT. pmvolume) THEN

            IF(PM0) THEN
               WRITE(iw,*) "Attempting a volume move"
               WRITE(iw,*) 
            ENDIF

            CALL mc_volume_move( mc_par,force_env,moves,move_updates,&
                        nnstep,old_energy,globenv,&
                        energy_check,r_old,&
                        error=error)

        ELSE

            IF(PM0) THEN
               WRITE(iw,*) "Attempting a Quickstep move"
               WRITE(iw,*) 
            ENDIF

            particle_positions(:,:)=0.0d0
! we want to do all this in serial
            IF( ionode ) THEN

               DO i=1,nmoves
                  
                  CALL RANDOM_NUMBER(rand)
                  IF(rand .LT. pmtraion) THEN
! change molecular conformation
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(rand*nchain)
                     mol=3*mol-2
                     CALL mc_atom_translation( mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ELSEIF(rand .LT. pmtrans) THEN
! translate a whole molecule in the system
! call a random number to figure out which molecule we're moving...the
! second line is necessary because we access a matrix with individual
! particles
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(nchain*rand)
                     mol=3*mol-2
                     CALL mc_molecule_translation(mc_par,force_env,&
                      moves,move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ELSE
! rotate a whole molecule in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING(nchain*rand)
                     mol=3*mol-2
                     CALL mc_molecule_rotation(mc_par,force_env, moves,&
                        move_updates,nnstep,i,mol,classical_energy,&
                        error=error)
                  ENDIF

               ENDDO

! get the new coordinates and dump them into an array so we can broadcast it
               DO iparticle=1,nunits_tot
                  particle_positions(1:3,iparticle)=&
                      particles_old%els(iparticle)%r(1:3)
               ENDDO


! back to parallel
            ENDIF

! broadcast all the positions, to make sure we're on the same page
            CALL mp_bcast(particle_positions,source,group)

            DO iparticle=1,nunits_tot
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
         CALL get_cell(cell,abc=abc)
         CALL cp_subsys_get(oldsys(1)%subsys,&
            particles=particles_old, error=error)


         IF(ionode) THEN
!            WRITE(ene,*) nnstep,old_energy/3.16682968D-6,energy_check
            WRITE(ene,*) nnstep,old_energy,energy_check

            IF(.NOT. lclassical_run) &
               CALL write_particle_coordinates(particles_old%els,cell,&
               crd,'XMOL','POS','POSITIONS')

! compute some averages
            averages%molecules=averages%molecules*&
                     REAL(nnstep-nstart-1,dp)&
                     /REAL(nnstep-nstart,dp)+nchain/&
                     REAL(nnstep-nstart,dp)

            averages%ave_energy=averages%ave_energy*&
                     REAL(nnstep-nstart-1,dp)&
                     /REAL(nnstep-nstart,dp)+old_energy/&
                     REAL(nnstep-nstart,dp)

            averages%ave_cell_length=averages%ave_cell_length*&
                     REAL(nnstep-nstart-1,dp)&
                     /REAL(nnstep-nstart,dp)+abc(1)/&
                     REAL(nnstep-nstart,dp)

!!!!!!!!!!!!!!!!!!!!!!!
!            WRITE(cell_unit,*) (abc(1)*angstrom)**3,old_energy
            WRITE(cell_unit,*) nnstep,abc(1)*angstrom

! find bond lengths and angle for current water config
            CALL get_water_data(force_env,averages,&
               nnstep-nstart,nchain,DATA)

! write the move statistics to file
            CALL write_move_stats(moves,move_updates,nnstep,move_unit)

! write a restart file
            CALL write_mc_restart(nnstep,mc_par,force_env)

! flush the buffers to the files
            CALL m_flush(DATA)
            CALL m_flush(diff)
            CALL m_flush(ene)
            CALL m_flush(move_unit)
            CALL m_flush(crd)
            CALL m_flush(cl)
            CALL m_flush(rm)
            CALL m_flush(cell_unit)

! check to see if EXIT file exists...if so, end the calculation
            CALL external_control(force_env,should_stop,error=error)

! update the move displacements, if necessary
            IF (MOD(nnstep-nstart,iuptrans) == 0) THEN
               CALL mc_move_update(mc_par,move_updates,"trans",nnstep)
            ENDIF
            IF (MOD(nnstep-nstart,iupvolume) == 0  ) THEN
               CALL mc_move_update(mc_par,move_updates,"volume",nnstep)
            ENDIF

         ENDIF

         CALL mp_bcast(should_stop,source,group)

         IF (should_stop) EXIT

! check to see if there are any overlaps in the boxe
         CALL mc_classical_energy(force_env,&
            rdum,nchain,1,nchain,cutoff,loverlap)
         IF(loverlap) THEN
            IF(ionode) WRITE(iw,*) nnstep
            CALL stop_program(routineN,module_name,__LINE__,&
            'coordinate overlap at the end of the above step')
         ENDIF

      ENDDO


! calculate the final energy

      IF (.NOT. lclassical_run) THEN
 
         CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.)
         CALL force_env_get(force_env,potential_energy=final_energy,&
            error=error)
   
      ELSE

         CALL mc_classical_energy(force_env,final_energy,&
            nchain,1,nchain,cutoff,ldum)

      ENDIF

! do some stuff in serial
      IF (ionode) THEN

         CALL write_particle_coordinates(particles_old%els,cell,&
                        crd,'XMOL','POS','FINAL')
                        
         WRITE(ene,*) 'Final Energy:                        ',&
                  final_energy,energy_check

! write a bunch of data to the screen
         CALL final_mc_write(mc_par,moves,iw,energy_check,&
                         initial_energy,final_energy,averages,&
                         move_updates)

! close any open files
         CALL close_file(unit_number=diff)
         CALL close_file(unit_number=crd)
         CALL close_file(unit_number=DATA)
         CALL close_file(unit_number=move_unit)
         CALL close_file(unit_number=ene)
         CALL close_file(unit_number=cell_unit)
         CALL close_file(unit_number=cl)
         CALL close_file(unit_number=rm)

! back to parallel
      ENDIF

      CALL set_mc_env ( mc_env, mc_par = mc_par,force_env=force_env)

! deallocate some stuff
      CALL mc_moves_release(move_updates)
      CALL mc_moves_release(moves)
      CALL mc_averages_release(averages)

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
                                                cell_unit,source,group,nstep,&
                                                nstart,nmoves,iuptrans,&
                                                iupvolume
    INTEGER, DIMENSION(1:2)                  :: box_flag, crd, DATA, ene, &
                                                move_unit,istep,nchain,&
                                                nunits_tot,diff,cl,rm
    LOGICAL                                  :: ionode, PM0, PM1, PM2, PM3, &
                                                PM4, should_stop,lstop,&
                                                lclassical_run,loverlap,&
                                                ldum
    REAL(KIND=dp)                            :: rand,pmvolume,pmswap,&
      pmtraion,pmtrans,cutoff,rdum
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: particle_positions, r_old
    REAL(KIND=dp), DIMENSION(1:3,1:2)        :: abc
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
    CHARACTER(LEN=40),DIMENSION(1:2)         :: data_file,moves_file,&
      energy_file,coords_file,displacement_file,cell_file,molecules_file

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

      NULLIFY(particles_old_1)
      NULLIFY(particles_old_2)

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

      CALL get_mc_par(mc_par(1)%mc_par,nunits_tot=nunits_tot(1),&
         nchain=nchain(1),ionode=ionode,source=source,group=group,&
         data_file=data_file(1),moves_file=moves_file(1),&
         cell_file=cell_file(1),coords_file=coords_file(1),&
         energy_file=energy_file(1),displacement_file=displacement_file(1),&
         lstop=lstop,nstep=nstep,nstart=nstart,pmvolume=pmvolume,&
         molecules_file=molecules_file(1),pmswap=pmswap,nmoves=nmoves,&
         pmtraion=pmtraion,pmtrans=pmtrans,iuptrans=iuptrans,&
         iupvolume=iupvolume,lclassical_run=lclassical_run,cutoff=cutoff)
      CALL get_mc_par(mc_par(2)%mc_par,nunits_tot=nunits_tot(2),&
         nchain=nchain(2),data_file=data_file(2),moves_file=moves_file(2),&
         cell_file=cell_file(2),coords_file=coords_file(2),&
         energy_file=energy_file(2),displacement_file=displacement_file(2),&
         molecules_file=molecules_file(2))
  
      ALLOCATE(particle_positions(3,nunits_tot(1)+&
              nunits_tot(2),2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
           "particle_positions",3*(nunits_tot(1)+nunits_tot(2))*dp_size*2)
      ALLOCATE(r_old(3,nunits_tot(1)+nunits_tot(2),2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
           "r_old",3*(nunits_tot(1)+nunits_tot(2))*dp_size*2)

! deal with the different levels of print statements
      PM0 = .FALSE.
      PM1 = .FALSE.
      PM2 = .FALSE.
      PM3 = .FALSE.
      PM4 = .FALSE.
      CALL write_checkpoint_information("entering "//routineN,globenv)
      iw = globenv%scr

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
      istep(:)=0

! initialize the moves array,the arrays for updating maximum move 
!   displacements, and the averages array
      DO i=1,2
         CALL init_mc_moves(moves(i)%moves)
         CALL init_mc_moves(move_updates(i)%moves)
         CALL mc_averages_create(averages(i)%averages)
      ENDDO

! find the energy of the initial configuration
      DO i=1,2

         IF (nchain(i) == 0) THEN
            old_energy(i)=0.0d0
            moves(i)%moves%classical_energy=0.0d0
            classical_energy(i)=moves(i)%moves%classical_energy
         ELSE

            IF (.NOT. lclassical_run) THEN
               CALL force_env_calc_energy_force(force_env(i)%force_env,&
                    calc_force=.FALSE.)
               CALL force_env_get(force_env(i)%force_env,&
               potential_energy=old_energy(i),error=error)
               
               CALL mc_classical_energy(force_env(i)%force_env,&
                  moves(i)%moves%classical_energy,nchain(i),1,&
                  nchain(i),cutoff,ldum)
               IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
                  'GEMC found an overlap in the initial config')
               classical_energy(i)=moves(i)%moves%classical_energy

            ELSE

               CALL mc_classical_energy(force_env(i)%force_env,&
                  moves(i)%moves%classical_energy,nchain(i),1,&
                  nchain(i),cutoff,ldum)
               IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
                  'GEMC found an overlap in the initial config')
               classical_energy(i)=moves(i)%moves%classical_energy
               old_energy(i)=moves(i)%moves%classical_energy

            ENDIF

         ENDIF

      ENDDO

      CALL force_env_get(force_env(1)%force_env,subsys=oldsys_1,&
         cell=cell_1,error=error)
      CALL get_cell(cell_1,abc=abc(:,1))
      CALL cp_subsys_get(oldsys_1(1)%subsys, &
         particles=particles_old_1, error=error)
      
      CALL force_env_get(force_env(2)%force_env,subsys=oldsys_2,&
         cell=cell_2,error=error)
      CALL get_cell(cell_2,abc=abc(:,2))
      CALL cp_subsys_get(oldsys_2(1)%subsys, &
         particles=particles_old_2, error=error)

! back to seriel
      IF (ionode) THEN

! open the files for writing coordinates, average structural data, and
! maxiumum displacements

! box 1
         DO i=1,2
            CALL open_file(file_name=data_file(i),&
               unit_number=DATA(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=moves_file(i),&
               unit_number=move_unit(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=coords_file(i),&
               unit_number=crd(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=energy_file(i),&
               unit_number=ene(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=displacement_file(i),&
               unit_number=rm(i),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=cell_file(i),&
                  unit_number=cl(i),&
                  file_action='WRITE',file_status='REPLACE')
         ENDDO
            
         CALL open_file(file_name='energy_differences_box1',&
                 unit_number=diff(1),&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene(1),*) 'Initial Energy Box 1:            ',&
            old_energy(1),' 0.0D0'
         IF(nchain(1) == 0) THEN
            WRITE(crd(1),*) ' 0'
            WRITE(crd(1),*) 'INITIAL'
         ELSE
            CALL write_particle_coordinates(particles_old_1%els,cell_1,&
               crd(1),'XMOL','POS','INITIAL')
         ENDIF
! box 2
         CALL open_file(file_name='energy_differences_box2',&
                 unit_number=diff(2),&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(ene(2),*) 'Initial Energy Box 2:            ',&
            old_energy(2),' 0.0D0'
         IF(nchain(2) == 0) THEN
            WRITE(crd(2),*) ' 0'
            WRITE(crd(2),*) 'INITIAL'
         ELSE
            CALL write_particle_coordinates(particles_old_2%els,cell_2,&
               crd(2),'XMOL','POS','INITIAL')
         ENDIF

! set all the units numbers that we just opened in the respective mc_par
         CALL set_mc_par(mc_par(1)%mc_par,cl=cl(1),rm=rm(1),diff=diff(1))
         CALL set_mc_par(mc_par(2)%mc_par,cl=cl(2),rm=rm(2),diff=diff(2))

! record the combined energies,coordinates, and cell lengths
         CALL open_file(file_name='cell_length',&
                 unit_number=cell_unit,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(cell_unit,*) 'Initial Cell Lengths: ',&
            abc(1,1)*angstrom,abc(1,2)*angstrom
         CALL open_file(file_name='combined_energies',&
                 unit_number=com_ene,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='combined_coordinates',&
                 unit_number=com_crd,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name=molecules_file(1),&
                 unit_number=com_mol,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(com_ene,*) 'Initial Total Energy:            ',&
            old_energy(1)+old_energy(2),' 0.0D0'
         IF(nchain(1) == 0) THEN
            WRITE(com_crd,*) ' 0'
            WRITE(com_crd,*) 'INITIAL BOX 1'
         ELSE
            CALL write_particle_coordinates(particles_old_1%els,cell_1,&
               com_crd,'XMOL','POS','INITIAL BOX 1')
         ENDIF
         IF(nchain(2) == 0) THEN
            WRITE(com_crd,*) ' 0'
            WRITE(com_crd,*) 'INITIAL BOX 2'
         ELSE
            CALL write_particle_coordinates(particles_old_2%els,cell_2,&
               com_crd,'XMOL','POS','INITIAL BOX 2')
         ENDIF
         initial_energy(1)=old_energy(1)
         initial_energy(2)=old_energy(2)

! record the old coordinates, in case a move is rejected
         DO iparticle=1,nunits_tot(1)
            r_old(1:3,iparticle,1)=particles_old_1%els(iparticle)%r(1:3)
         ENDDO
         DO iparticle=1,nunits_tot(2)
            r_old(1:3,iparticle,2)=particles_old_2%els(iparticle)%r(1:3)
         ENDDO

! back to parallel mode
      ENDIF

! find out how many steps we're doing
      IF (.NOT. lstop) &
         nstep=nstep*(nchain(1)+nchain(2))

!    MC STUFF

      DO nnstep=nstart+1,nstart + nstep

         IF(PM0) THEN
            WRITE(iw,*) 
            WRITE(iw,*) "------- On Monte Carlo Step ",nnstep
         ENDIF
         CALL RANDOM_NUMBER(rand)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,source,group)


         IF (rand .LT. pmvolume) THEN
!     change the volume of the system
            IF(PM0) THEN
               WRITE(iw,*) "Attempting a volume move"
               WRITE(iw,*) 
            ENDIF

            CALL mc_ge_volume_move ( mc_par,force_env, moves,&
               move_updates,nnstep,old_energy,globenv,energy_check,&
               r_old,error)

         ELSEIF ( rand .LT. pmswap) THEN

! try a swap move
            IF(PM0) THEN
               WRITE(iw,*) "Attempting a swap move"
               WRITE(iw,*) 
            ENDIF

            CALL mc_ge_swap_move(mc_par,force_env,moves,&
                energy_check,r_old,old_energy,&
                nnstep,globenv,error)

! the number of molecules may have changed
            CALL get_mc_par(mc_par(1)%mc_par,nunits_tot=nunits_tot(1),&
               nchain=nchain(1))
            CALL get_mc_par(mc_par(2)%mc_par,nunits_tot=nunits_tot(2),&
               nchain=nchain(2))


         ELSE

            IF(PM0) THEN
               WRITE(iw,*) "Attempting a Quickstep move"
               WRITE(iw,*) 
            ENDIF

            particle_positions(:,:,:) = 0.0D0
            
! go into serial mode
            IF (ionode) THEN

               DO i=1,nmoves
                  CALL RANDOM_NUMBER(rand)
                  IF(rand .LT. pmtraion) THEN
                     
!     translate a single atom in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((nchain(1)+nchain(2))*rand)
                     mol=3*mol-2
                     IF (mol .LE. nunits_tot(1)) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-nunits_tot(1)
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_atom_translation( mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,i,mol,&
                        classical_energy(ibox),error=error)
                        
                  ELSEIF(rand .LT. pmtrans) THEN
!     translate a whole molecule in the system
!     call a random number to figure out which molecule we're moving...the
!     second line is necessary because we access a matrix with individual
!     particles
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((nchain(1)+nchain(2))*rand)
                     mol=3*mol-2
                     IF (mol .LE. nunits_tot(1)) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-nunits_tot(1)
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_molecule_translation(mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,&
                        i,mol,classical_energy(ibox),error=error)

                  ELSE
!     rotate a whole molecule in the system
                     CALL RANDOM_NUMBER(rand)
                     mol=CEILING((nchain(1)+nchain(2))*rand)
                     mol=3*mol-2
                     IF (mol .LE. nunits_tot(1)) THEN
                        ibox=1
                     ELSE
                        ibox=2
                        mol=mol-nunits_tot(1)
                     ENDIF

                     box_flag(ibox)=1
                     CALL mc_molecule_rotation(mc_par(ibox)%mc_par,&
                        force_env(ibox)%force_env,moves(ibox)%moves,&
                        move_updates(ibox)%moves,nnstep,i,mol,&
                        classical_energy(ibox),error=error)
                     
                  ENDIF
               ENDDO
               
! get the new coordinates and dump them into an array so we can broadcast it
               DO iparticle=1,nunits_tot(1)
                  particle_positions(1:3,iparticle,1)=&
                  particles_old_1%els(iparticle)%r(1:3)
               ENDDO
               DO iparticle=1,nunits_tot(2)
                  particle_positions(1:3,iparticle,2)=&
                     particles_old_2%els(iparticle)%r(1:3)
               ENDDO
               
!     go back to parallel      
            ENDIF
            
!     broadcast all the positions, to make sure we're on the same page
            CALL mp_bcast(particle_positions,source,group)
            
            DO iparticle=1,nunits_tot(1)
               particles_old_1%els(iparticle)%r(1:3)= &
                  particle_positions(1:3,iparticle,1)
            ENDDO
            DO iparticle=1,nunits_tot(2)
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
         CALL get_cell(cell_1,abc=abc(:,1))
         CALL cp_subsys_get(oldsys_1(1)%subsys,&
            particles=particles_old_1, error=error)
         CALL force_env_get(force_env(2)%force_env,cell=cell_2,&
            subsys=oldsys_2,error=error)
         CALL get_cell(cell_2,abc=abc(:,2))
         CALL cp_subsys_get(oldsys_2(1)%subsys, &
            particles=particles_old_2, error=error)


!     find bond lengths, angle, energy, and number of molecules
!     for current water config 
         IF (ionode) THEN
            DO ibox=1,2
               CALL write_move_stats(moves(ibox)%moves,&
                  move_updates(ibox)%moves,nnstep,&
                  move_unit(ibox))
               IF(nchain(ibox) .NE. 0) THEN
                  istep(ibox)=istep(ibox)+1
                  CALL get_water_data( force_env(ibox)%force_env,&
                     averages(ibox)%averages,istep(ibox),&
                     nchain(ibox),DATA(ibox))
               ENDIF
               WRITE(ene(ibox),*) nnstep,old_energy(ibox)
               averages(ibox)%averages%ave_energy=&
                  averages(ibox)%averages%ave_energy*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  old_energy(ibox)/REAL(nnstep,dp)
               averages(ibox)%averages%molecules=&
                  averages(ibox)%averages%molecules*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  REAL(nchain(ibox),dp)/REAL(nnstep-nstart,dp)
            ENDDO
            WRITE(com_ene,*) nnstep,old_energy(2)+old_energy(1)
            IF (pmswap .GT. 0.0D0) &
               WRITE(com_mol,*) nnstep,nchain(1),nchain(2)

! compute the average cell length, energy, and number of particles of
! both boxes, as well as the total average energy...write to files
            averages(1)%averages%ave_cell_length=&
               averages(1)%averages%ave_cell_length*&
               REAL(nnstep-nstart-1,dp)/REAL(nnstep-nstart,dp)+&
               abc(1,1)/REAL(nnstep-nstart,dp)
            averages(2)%averages%ave_cell_length=&
               averages(2)%averages%ave_cell_length*&
               REAL(nnstep-nstart-1,dp)/&
               REAL(nnstep-nstart,dp)+&
               abc(1,2)/REAL(nnstep-nstart,dp)

            WRITE(cell_unit,*) nnstep,abc(1,1)*angstrom,&
               abc(1,2)*angstrom
            
            

!     write the coordinates and energy to a file
            IF(.NOT. lclassical_run) THEN
               IF (nchain(1) == 0) THEN
                  WRITE(com_crd,*) ' 0 '
                  WRITE(com_crd,*) ' POSITIONS BOX 1 '
                  WRITE(crd(1),*) ' 0 '
                  WRITE(crd(1),*) ' POSITIONS BOX 1 '
               ELSE
                  CALL write_particle_coordinates(&
                     particles_old_1%els,cell_1,com_crd,'XMOL',&
                     'POS','POSITIONS BOX 1')
                  CALL write_particle_coordinates(&
                     particles_old_1%els,cell_1,crd(1),'XMOL',&
                     'POS','POSITIONS BOX 1')
               ENDIF
               IF (nchain(2) == 0) THEN
                  WRITE(com_crd,*) ' 0 '
                  WRITE(com_crd,*) ' POSITIONS BOX 2 '
                  WRITE(crd(2),*) ' 0 '
                  WRITE(crd(2),*) ' POSITIONS BOX 2 '
               ELSE
                  CALL write_particle_coordinates(&
                     particles_old_2%els,cell_2,com_crd,'XMOL',&
                     'POS','POSITIONS BOX 2')
                  CALL write_particle_coordinates(&
                     particles_old_2%els,cell_2,crd(2),'XMOL',&
                     'POS','POSITIONS BOX 2')
               ENDIF
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
               CALL m_flush(diff(i))
               CALL m_flush(ene(i))
               CALL m_flush(crd(i))
               CALL m_flush(cl(i))
               CALL m_flush(rm(i))
            ENDDO
            CALL m_flush(com_ene)
            CALL m_flush(com_crd)
            CALL m_flush(com_mol)
            CALL m_flush(cell_unit)

!     reset the box flags
            box_flag(1)=0
            box_flag(2)=0
            
! update the move displacements, if necessary
            IF (MOD(nnstep-nstart,iuptrans) == 0) THEN
               CALL mc_move_update(mc_par(1)%mc_par,&
                  move_updates(1)%moves,"trans",nnstep)
               CALL mc_move_update(mc_par(2)%mc_par,&
                  move_updates(2)%moves,"trans",nnstep)
            ENDIF
            IF (MOD(nnstep-nstart,iupvolume) == 0 ) THEN
               CALL mc_move_update(mc_par(1)%mc_par,&
                  move_updates(1)%moves,"volume",nnstep)
               CALL mc_move_update(mc_par(2)%mc_par,&
                  move_updates(2)%moves,"volume",nnstep)
            ENDIF


!     back to parallel
         ENDIF
         
!     make sure all the procs exit if necessary
         CALL mp_bcast(should_stop,source,group)
         IF (should_stop) EXIT
         
! check to see if there are any overlaps in either of the boxes
!         CALL check_for_overlap(r_old(:,:,1),nunits_tot(1),abc(1,1),&
!            loverlap)
         CALL mc_classical_energy(force_env(1)%force_env,&
            rdum,nchain(1),1,nchain(1),cutoff,loverlap)
         IF(loverlap) THEN
            IF(ionode) WRITE(iw,*) nnstep
            CALL stop_program(routineN,module_name,__LINE__,&
            'coordinate overlap  in box 1 at the end of the above step')
         ENDIF
!         CALL check_for_overlap(r_old(:,:,2),nunits_tot(2),abc(1,2),&
!            loverlap)
         CALL mc_classical_energy(force_env(2)%force_env,&
            rdum,nchain(2),1,nchain(2),cutoff,loverlap)
         IF(loverlap) THEN
            IF(ionode) WRITE(iw,*) nnstep
            CALL stop_program(routineN,module_name,__LINE__,&
            'coordinate overlap  in box 2 at the end of the above step')
         ENDIF

      ENDDO
      
      DO i=1,2
         IF (nchain(i) .NE. 0) THEN
            IF (.NOT. lclassical_run) THEN
               CALL force_env_calc_energy_force(force_env(i)%force_env,&
                    calc_force=.FALSE.)
               CALL force_env_get(force_env(i)%force_env,&
                  potential_energy=final_energy(i),error=error)
            ELSE
               CALL mc_classical_energy(force_env(i)%force_env,&
                  final_energy(i),nchain(i),1,nchain(i),cutoff,ldum)
            ENDIF
         ELSE
            final_energy(i)=0.0D0
         ENDIF
      ENDDO

      IF (ionode) THEN
         CALL force_env_get(force_env(1)%force_env,&
            subsys=oldsys_1,cell=cell_1,error=error)
         CALL cp_subsys_get(oldsys_1(1)%subsys,&
            particles=particles_old_1, error=error)
         IF(nchain(1) == 0) THEN
            WRITE(crd(1),*) ' 0'
            WRITE(crd(1),*) 'FINAL'
            WRITE(com_crd,*) ' 0'
            WRITE(com_crd,*) 'FINAL BOX 1'
         ELSE
            CALL write_particle_coordinates(particles_old_1%els,cell_1,&
               com_crd,'XMOL','POS','FINAL BOX 1')
            CALL write_particle_coordinates(particles_old_1%els,cell_1,&
               crd(1),'XMOL','POS','FINAL BOX 1')
         ENDIF
         CALL force_env_get(force_env(2)%force_env,&
            subsys=oldsys_2,cell=cell_2,error=error)
         CALL cp_subsys_get(oldsys_2(1)%subsys,&
            particles=particles_old_2, error=error)
         IF(nchain(2) == 0) THEN
            WRITE(crd(2),*) ' 0'
            WRITE(crd(2),*) 'FINAL'
            WRITE(com_crd,*) ' 0'
            WRITE(com_crd,*) 'FINAL BOX 2'
         ELSE
            CALL write_particle_coordinates(particles_old_2%els,cell_2,&
               com_crd,'XMOL','POS','FINAL BOX 2')
            CALL write_particle_coordinates(particles_old_2%els,cell_2,&
               crd(2),'XMOL','POS','FINAL BOX 2')
         ENDIF
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
            CALL close_file(unit_number=diff(i))
            CALL close_file(unit_number=crd(i))
            CALL close_file(unit_number=DATA(i))
            CALL close_file(unit_number=move_unit(i))
            CALL close_file(unit_number=ene(i))
            CALL close_file(unit_number=cl(i))
            CALL close_file(unit_number=rm(i))
         ENDDO
         CALL close_file(unit_number=com_ene)
         CALL close_file(unit_number=com_crd)
         CALL close_file(unit_number=com_mol)
         CALL close_file(unit_number=cell_unit)

      ENDIF

!     deallocate some stuff
      CALL set_mc_env ( mc_env_1, mc_par = mc_par(1)%mc_par,&
         force_env=force_env(1)%force_env)
      CALL set_mc_env ( mc_env_2, mc_par = mc_par(2)%mc_par,&
         force_env=force_env(2)%force_env)
      NULLIFY(particles_old_1)
      NULLIFY(particles_old_2)

      DEALLOCATE(particle_positions,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions")
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old deallocate")
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
         CALL mc_averages_release(averages(i)%averages)
      ENDDO
      
      DEALLOCATE(averages,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "averages deallocate")
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

