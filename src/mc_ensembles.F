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
                                             cp_subsystem_p_type,&
                                             cp_subsystem_pp_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE mc_control,                      ONLY: write_mc_restart,&
                                             mc_create_bias_force_env
  USE mc_coordinates,                  ONLY: get_water_data,&
                                             mc_classical_energy,&
                                             find_mc_test_molecule,&
                                             create_discrete_array
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_environment_type,&
                                             mc_environment_p_type,&
                                             set_mc_env
  USE mc_ge_moves,                     ONLY: mc_ge_swap_move,&
                                             mc_quickstep_move_new,&
                                             mc_ge_volume_move
  USE mc_misc,                         ONLY: final_mc_write,&
                                             mc_averages_create,&
                                             mc_averages_release
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_move_update,&
                                             mc_moves_release,&
                                             write_move_stats
  USE mc_moves,                        ONLY: mc_atom_translation,&
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
  USE particle_list_types,             ONLY: particle_list_type,&
                                             particle_list_p_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             cell_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_constants,                 ONLY: dump_xmol
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ensembles"

  PUBLIC :: mc_run_ensemble

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_ensembles.F/mc_run_ensemble
!!
!!   NAME
!!     mc_traditional_ensemble
!!
!!   FUNCTION
!!     directs the program in running one or two box MC simulations
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
  SUBROUTINE mc_run_ensemble ( mc_env, globenv, nboxes, error )

!   ***************************************************************************

    TYPE(mc_environment_p_type), POINTER,DIMENSION(:) &
                                             :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER,INTENT(IN)                        :: nboxes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_run_ensemble"

    INTEGER, DIMENSION(1:nboxes)             :: box_flag, data_unit,&
                                                move_unit,istep,nchain,&
                                                nunits_tot,diff,cl,rm
    INTEGER,DIMENSION(1:3,1:2)               :: discrete_array
    INTEGER                                  :: handle, imove, &
                                                iparticle, iw, nnstep,&
                                                mol,istat,cell_unit,&
                                                source,group,&
                                                nstep,nstart,nmoves,&
                                                iuptrans,iupvolume,iprint,&
                                                box_number,ibox,&
                                                com_ene,com_crd,com_mol
    LOGICAL                               :: PM0, PM1, PM2, PM3, PM4, &
                                              should_stop,ionode,&
                                              lstop,ldum,&
                                              loverlap,lbias,lnew_bias_env,&
                                              lreject,ldiscrete
    REAL(KIND=dp), DIMENSION(1:3,1:nboxes)    :: abc
    REAL(KIND=dp)                    ::  rand,pmtraion,&
                                         pmvolume,test_energy,&
                                         pmtrans,cutoff,rdum,inner_cutoff,&
                                         pmswap,discrete_step
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: particle_positions,r_old
    REAL(KIND=dp), DIMENSION(1:nboxes)       :: bias_energy,&
                                              energy_check,final_energy,&
                                              initial_energy,old_energy
    TYPE(cell_p_type),POINTER,DIMENSION(:) &
                                             :: cell
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env,bias_env
    TYPE(mc_averages_p_type), DIMENSION(:), &
      POINTER                                :: averages
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: move_updates, moves
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_p_type), DIMENSION(:),POINTER &
                                             :: particles_old
    CHARACTER(LEN=40),DIMENSION(1:nboxes)    :: data_file,moves_file,&
      energy_file,coords_file,displacement_file,cell_file,molecules_file

    CHARACTER(LEN=40)                        :: fft_lib,cbox,move_type
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(moves,move_updates)

! allocate a whole bunch of stuff based on how many boxes we have
      ALLOCATE (force_env(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "force_env",1*dp_size*nboxes)
      ALLOCATE (bias_env(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "bias_env",1*dp_size*nboxes)
      ALLOCATE (cell(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "cell",1*dp_size*nboxes)
      ALLOCATE (particles_old(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old",1*dp_size*nboxes)
      ALLOCATE (oldsys(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys",1*dp_size*nboxes)
      ALLOCATE (averages(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "averages",1*dp_size*nboxes)
      ALLOCATE (moves(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "moves",1*dp_size*nboxes)
      ALLOCATE (move_updates(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "move_updates",1*dp_size*nboxes)
      ALLOCATE (mc_par(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "mc_par",1*dp_size*nboxes)

      DO ibox=1,nboxes
         CALL get_mc_env ( mc_env(ibox)%mc_env, &
            mc_par = mc_par(ibox)%mc_par,&
            force_env=force_env(ibox)%force_env)
      ENDDO

! get some data out of mc_par
      CALL get_mc_par(mc_par(1)%mc_par,nunits_tot=nunits_tot(1),&
         nchain=nchain(1),ionode=ionode,source=source,group=group,&
         data_file=data_file(1),moves_file=moves_file(1),&
         cell_file=cell_file(1),coords_file=coords_file(1),&
         energy_file=energy_file(1),displacement_file=displacement_file(1),&
         lstop=lstop,nstep=nstep,nstart=nstart,pmvolume=pmvolume,&
         molecules_file=molecules_file(1),pmswap=pmswap,nmoves=nmoves,&
         pmtraion=pmtraion,pmtrans=pmtrans,iuptrans=iuptrans,&
         iupvolume=iupvolume,cutoff=cutoff,ldiscrete=ldiscrete,&
         inner_cutoff=inner_cutoff,lbias=lbias,iprint=iprint,&
         discrete_step=discrete_step,fft_lib=fft_lib)

      IF(nboxes .GT. 1) THEN
         DO ibox=2,nboxes
            CALL get_mc_par(mc_par(ibox)%mc_par,nunits_tot=nunits_tot(ibox),&
            nchain=nchain(ibox),data_file=data_file(ibox),&
            moves_file=moves_file(ibox),&
            cell_file=cell_file(ibox),coords_file=coords_file(ibox),&
            energy_file=energy_file(ibox),&
            displacement_file=displacement_file(ibox),&
            molecules_file=molecules_file(ibox))
         ENDDO
      ENDIF

! allocate the particle positions array for broadcasting
      ALLOCATE (particle_positions(3,SUM(nunits_tot),1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions",3*SUM(nunits_tot)*dp_size*nboxes)
      ALLOCATE (r_old(3,SUM(nunits_tot),1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old",3*SUM(nunits_tot)*dp_size*nboxes)


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
         WRITE(iw,*) 'Beginning the Monte Carlo calculation.'
         WRITE(iw,*)
         WRITE(iw,*)
      ENDIF

! initialize running average variables
      energy_check(:)=0.0E0_dp
      box_flag(:)=0
      istep(:)=0

      DO ibox=1,nboxes
! initialize the moves array, the arrays for updating maximum move 
! displacements, and the averages array
         CALL init_mc_moves(moves(ibox)%moves)
         CALL init_mc_moves(move_updates(ibox)%moves)
         CALL mc_averages_create(averages(ibox)%averages)

! find the energy of the initial configuration
         IF(nchain(ibox) .NE. 0) THEN
            CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.)
            CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=old_energy(ibox),&
               error=error)
         ELSE
            old_energy(ibox)=0.0d0
         ENDIF
         initial_energy(ibox)=old_energy(ibox)
      
! find the bias energy of the initial run
         IF(lbias) THEN
            CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
               force_env(ibox)%force_env,&
               globenv,.FALSE.,TRIM(fft_lib),error)
            IF(nchain(ibox) .NE. 0) THEN
               CALL force_env_calc_energy_force(bias_env(ibox)%force_env,&
                  calc_force=.FALSE.)
               CALL force_env_get(bias_env(ibox)%force_env,&
                  potential_energy=moves(ibox)%moves%bias_energy,error=error)
            ELSE
               moves(ibox)%moves%bias_energy=0.0d0
            ENDIF
            bias_energy(ibox)=moves(ibox)%moves%bias_energy
         ENDIF
         lnew_bias_env=.FALSE.
         CALL mc_classical_energy(force_env(ibox)%force_env,rdum,&
         nchain(ibox),1,nchain(ibox),cutoff,ldum)
         IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
            "overlap in an initial configuration")

! get the subsystems and the cell information
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
            error=error)
         CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
         CALL cp_subsys_get(oldsys(ibox)%subsys(1)%subsys, &
            particles=particles_old(ibox)%list, error=error)
! record the old coordinates, in case a move is rejected
         DO iparticle=1,nunits_tot(ibox)
            r_old(1:3,iparticle,ibox)=&
               particles_old(ibox)%list%els(iparticle)%r(1:3)
         ENDDO
      ENDDO

! back to seriel for a bunch of I/O stuff
      IF (ionode) THEN

! record the combined energies,coordinates, and cell lengths
         CALL open_file(file_name='mc_cell_length',&
                 unit_number=cell_unit,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='mc_energies',&
                 unit_number=com_ene,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='mc_coordinates',&
                 unit_number=com_crd,&
                 file_action='WRITE',file_status='REPLACE')
         CALL open_file(file_name='mc_molecules',&
                 unit_number=com_mol,&
                 file_action='WRITE',file_status='REPLACE')
         WRITE(com_ene,*) 'Initial Energies:       ',&
            old_energy(1:nboxes)
         WRITE(com_mol,*) 'Initial Molecules:       ',&
            nchain(1:nboxes)
         DO ibox=1,nboxes
            WRITE(cell_unit,*) 'Initial: ',&
               abc(1:3,ibox)*angstrom
            WRITE(cbox,'(I4)') ibox
            CALL open_file(file_name='energy_differences_box' // &
                 TRIM(ADJUSTL(cbox)),&
                 unit_number=diff(ibox),&
                 file_action='WRITE',file_status='REPLACE')
            IF(nchain(ibox) == 0) THEN
               WRITE(com_crd,*) ' 0'
               WRITE(com_crd,*) 'INITIAL BOX ' // TRIM(ADJUSTL(cbox))
            ELSE
               CALL write_particle_coordinates(particles_old(ibox)%list%els,&
                  cell(ibox)%cell,com_crd,dump_xmol,'POS',&
                  'INITIAL BOX ' // TRIM(ADJUSTL(cbox)))
            ENDIF
            CALL open_file(file_name=data_file(ibox),&
               unit_number=data_unit(ibox),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=moves_file(ibox),&
               unit_number=move_unit(ibox),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=displacement_file(ibox),&
               unit_number=rm(ibox),&
               file_action='WRITE',file_status='REPLACE')
            CALL open_file(file_name=cell_file(ibox),&
               unit_number=cl(ibox),&
               file_action='WRITE',file_status='REPLACE')

! set all the units numbers that we just opened in the respective mc_par
            CALL set_mc_par(mc_par(ibox)%mc_par,cl=cl(ibox),rm=rm(ibox),&
               diff=diff(ibox))
         ENDDO

! back to parallel mode
      ENDIF

! if we're doing a discrete volume move, we need to set up the array
! that keeps track of which direction we can move in...1 means we can
! increase, -1 means we can decrease

! if we're doing a discrete volume move, we need to set up the array
! that keeps track of which direction we can move in
      IF(ldiscrete) THEN
         IF(nboxes .NE. 1) &
            CALL stop_program(routineN,module_name,__LINE__,&
            'ldiscrete=.true. ONLY for systems with 1 box')
         CALL create_discrete_array(abc(:,1),discrete_array(:,:),&
            discrete_step)
      ENDIF

! find out how many steps we're doing...change the updates to be in cycles
! if the total number of steps is measured in cycles
      IF (.NOT. lstop ) THEN
         nstep=nstep*SUM(nchain(:))
         iuptrans=iuptrans*SUM(nchain(:))
         iupvolume=iupvolume*SUM(nchain(:))
      ENDIF

      DO nnstep=nstart+1,nstart + nstep

         IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
            WRITE(iw,*) 
            WRITE(iw,*) "------- On Monte Carlo Step ",nnstep
         ENDIF

         CALL RANDOM_NUMBER(rand)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,source,group)

         IF (rand .LT. pmvolume) THEN

            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting a volume move"
               WRITE(iw,*) 
            ENDIF

            IF( nboxes .EQ. 1) THEN
               CALL mc_volume_move( mc_par(1)%mc_par,&
                  force_env(1)%force_env,&
                  moves(1)%moves,move_updates(1)%moves,&
                  nnstep,old_energy(1),globenv,&
                  energy_check(1),r_old(:,:,1),iw,discrete_array(:,:),&
                  error=error)

            ELSE
               CALL mc_ge_volume_move ( mc_par,force_env, moves,&
                  move_updates,nnstep,old_energy,globenv,energy_check,&
                  r_old,error)
            ENDIF

! we need a new biasing environment now, if we're into that sort of thing
            IF(lbias) THEN
               DO ibox=1,nboxes
                  CALL force_env_release(bias_env(ibox)%force_env,error=error)
                  CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
                     force_env(ibox)%force_env,&
                     globenv,.FALSE.,TRIM(fft_lib),error)
                  IF(nchain(ibox) .NE. 0) THEN
                     CALL force_env_calc_energy_force(&
                        bias_env(ibox)%force_env,&
                        calc_force=.FALSE.)
                     CALL force_env_get(bias_env(ibox)%force_env,&
                        potential_energy=moves(ibox)%moves%bias_energy,&
                        error=error)
                  ELSE
                    moves(ibox)%moves%bias_energy=0.0d0 
                  ENDIF
                  bias_energy(ibox)=moves(ibox)%moves%bias_energy
               ENDDO
            ENDIF
         ELSEIF ( rand .LT. pmswap) THEN

! try a swap move
            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting a swap move"
               WRITE(iw,*) 
            ENDIF

            CALL mc_ge_swap_move(mc_par,force_env,bias_env,moves,&
                energy_check(:),r_old(:,:,:),old_energy(:),&
                nnstep,globenv,bias_energy(:),error)

! the number of molecules may have changed
            DO ibox=1,nboxes
               CALL get_mc_par(mc_par(ibox)%mc_par,&
                  nunits_tot=nunits_tot(ibox),&
                  nchain=nchain(ibox))
            ENDDO
        ELSE

            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting a Quickstep move"
               WRITE(iw,*) 
            ENDIF

            particle_positions(:,:,:)=0.0E0_dp
! we want to do all this in serial
            IF( ionode ) THEN

               DO imove=1,nmoves
                  
                  CALL RANDOM_NUMBER(rand)
                  IF(rand .LT. pmtraion) THEN
! change molecular conformation
! first, pick a box to do it for
                     CALL RANDOM_NUMBER(rand)
                     IF(nboxes .EQ. 2) THEN
                        IF(rand .LT. 0.75E0_dp) THEN
                           ibox=1
                        ELSE
                           ibox=2
                        ENDIF
                     ELSE
                        ibox=1
                     ENDIF

! now pick a molecule, automatically rejecting the move if the
! box is empty
                     IF(nchain(ibox) == 0) THEN
! indicate that we tried a move
                        moves(ibox)%moves%empty_conf=&
                           moves(ibox)%moves%empty_conf+1
                     ELSE
! pick a molecule in the box
                        CALL find_mc_test_molecule(force_env,mol,&
                           box_number,nchain(:),box=ibox)

! choose if we're changing a bond length or an angle
                        CALL RANDOM_NUMBER(rand)

! used for keeping track of counters for lbias=.false. simulations
                        IF (rand .LT. 0.66666667E0_dp) THEN
                           move_type='bond'
                        ELSE
                           move_type='angle'
                        ENDIF
                        box_flag(box_number)=1
                        CALL mc_atom_translation(mc_par(box_number)%mc_par,&
                           force_env(box_number)%force_env,&
                           bias_env(box_number)%force_env,&
                           moves(box_number)%moves,&
                           move_updates(box_number)%moves,nnstep,imove,&
                           mol,bias_energy(box_number),&
                           move_type,lreject,error=error)
                        IF(lreject) EXIT
                     ENDIF
                  ELSEIF(rand .LT. pmtrans) THEN
! translate a whole molecule in the system
! call a random number to figure out which molecule we're moving...the
! second line is necessary because we access a matrix with individual
! particles
                     move_type='trans'
                     CALL find_mc_test_molecule(force_env,mol,box_number,&
                        nchain(:))
                     box_flag(box_number)=1
                     CALL mc_molecule_translation(mc_par(box_number)%mc_par,&
                        force_env(box_number)%force_env,&
                        bias_env(box_number)%force_env,&
                        moves(box_number)%moves,&
                        move_updates(box_number)%moves,nnstep,imove,mol,&
                        bias_energy(box_number),lreject,error=error)
                     IF(lreject) EXIT
                  ELSE
! rotate a whole molecule in the system
                     move_type='rot'
                     CALL find_mc_test_molecule(force_env,mol,box_number,&
                        nchain(:))
                     box_flag(box_number)=1
                     CALL mc_molecule_rotation(mc_par(box_number)%mc_par,&
                        force_env(box_number)%force_env,&
                        bias_env(box_number)%force_env,&
                        moves(box_number)%moves,&
                        move_updates(box_number)%moves,nnstep,imove,&
                        mol,bias_energy(box_number),&
                        lreject,error=error)
                     IF(lreject) EXIT
                  ENDIF

               ENDDO

! get the new coordinates and dump them into an array so we can broadcast it
               DO ibox=1,nboxes
                  DO iparticle=1,nunits_tot(ibox)
                     particle_positions(1:3,iparticle,ibox)=&
                        particles_old(ibox)%list%els(iparticle)%r(1:3)
                  ENDDO
               ENDDO


! back to parallel
            ENDIF

! broadcast the reject flag
            CALL mp_bcast(lreject,source,group)

! broadcast all the positions, to make sure we're on the same page
            CALL mp_bcast(particle_positions,source,group)

            DO ibox=1,nboxes
               DO iparticle=1,nunits_tot(ibox)
                  particles_old(ibox)%list%els(iparticle)%r(1:3)= &
                     particle_positions(1:3,iparticle,ibox)
               ENDDO
            ENDDO

! now do a Quickstep calculation to see if we accept the sequence
            CALL mc_Quickstep_move_new(mc_par,force_env,bias_env,&
               moves,lreject,move_updates,energy_check(:),r_old(:,:,:),&
               nnstep,old_energy(:),bias_energy(:),&
               nboxes,box_flag(:),oldsys,particles_old,cell,&
               error=error)

         ENDIF

! make sure the pointers are pointing correctly since the subsys may
! have changed
         DO ibox=1,nboxes
            CALL force_env_get(force_env(ibox)%force_env,&
               subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
               error=error)
            CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
            CALL cp_subsys_get(oldsys(ibox)%subsys(1)%subsys,&
               particles=particles_old(ibox)%list, error=error)
         ENDDO

         IF(ionode) THEN
            
            IF(MOD(nnstep,iprint) == 0) THEN
               WRITE(com_ene,*) nnstep,old_energy(1:nboxes)
               WRITE(com_mol,*) nnstep,nchain(1:nboxes)

               DO ibox=1,nboxes

! write the move statistics to file
                  CALL write_move_stats(moves(ibox)%moves,&
                     move_updates(ibox)%moves,nnstep,move_unit(ibox))

! write a restart file
                  CALL write_mc_restart(nnstep,mc_par(ibox)%mc_par,&
                     force_env(ibox)%force_env)

! write cell lengths
                  WRITE(cell_unit,*) nnstep,abc(1:3,ibox)*angstrom

! write particle coordinates
                  WRITE(cbox,'(I4)') ibox
                  IF(nchain(ibox) == 0) THEN
                     WRITE(com_crd,*) ' 0'
                     WRITE(com_crd,*) 'BOX ' // TRIM(ADJUSTL(cbox))
                  ELSE
                     CALL write_particle_coordinates(&
                        particles_old(ibox)%list%els,&
                        cell(ibox)%cell,com_crd,dump_xmol,'POS',&
                        'BOX ' // TRIM(ADJUSTL(cbox)))
                  ENDIF

               ENDDO
            ENDIF

            DO ibox=1,nboxes
! compute some averages
               averages(ibox)%averages%ave_energy=&
                  averages(ibox)%averages%ave_energy*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  old_energy(ibox)/REAL(nnstep,dp)
               averages(ibox)%averages%molecules=&
                  averages(ibox)%averages%molecules*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  REAL(nchain(ibox),dp)/REAL(nnstep-nstart,dp)
               averages(ibox)%averages%ave_volume=&
                  averages(ibox)%averages%ave_volume*&
                  REAL(nnstep-nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  abc(1,ibox)*abc(2,ibox)*abc(3,ibox)/&
                  REAL(nnstep-nstart,dp)

! find bond lengths and angle for current water config
               IF(nchain(ibox) .NE. 0) THEN
                  istep(ibox)=istep(ibox)+1
                  CALL get_water_data(force_env(ibox)%force_env,&
                     averages(ibox)%averages,&
                     istep(ibox),nchain(ibox),data_unit(ibox))
               ENDIF

! update the move displacements, if necessary
               IF (MOD(nnstep-nstart,iuptrans) == 0) THEN
                  CALL mc_move_update(mc_par(ibox)%mc_par,&
                     move_updates(ibox)%moves,"trans",nnstep)
               ENDIF

               IF (MOD(nnstep-nstart,iupvolume) == 0  ) THEN
                  CALL mc_move_update(mc_par(ibox)%mc_par,&
                     move_updates(ibox)%moves,"volume",nnstep)
               ENDIF
                  
! flush the buffers to the files
               CALL m_flush(data_unit(ibox))
               CALL m_flush(diff(ibox))
               CALL m_flush(move_unit(ibox))
               CALL m_flush(cl(ibox))
               CALL m_flush(rm(ibox))

            ENDDO

! flush more buffers to the files
            CALL m_flush(cell_unit)
            CALL m_flush(com_ene)
            CALL m_flush(com_crd)
            CALL m_flush(com_mol)

! reset the box flags
            box_flag(:)=0

! check to see if EXIT file exists...if so, end the calculation
            CALL external_control(force_env(1)%force_env,&
               should_stop,error=error)

         ENDIF

         CALL mp_bcast(should_stop,source,group)

         IF (should_stop) EXIT

! check to see if there are any overlaps in the boxe
         DO ibox=1,nboxes
            IF(nchain(ibox) .NE. 0) THEN
               CALL mc_classical_energy(force_env(ibox)%force_env,&
                  rdum,nchain(ibox),1,nchain(ibox),cutoff,loverlap)
               IF(loverlap) THEN
                  IF(ionode) WRITE(iw,*) nnstep
                  CALL stop_program(routineN,module_name,__LINE__,&
                  'coordinate overlap at the end of the above step')
               ENDIF
            ENDIF
         ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug make sure the energies line up
         IF(.false.) THEN
         DO ibox=1,nboxes
            IF(nchain(ibox) .NE. 0) THEN
               CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.)
               CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=test_energy,&
               error=error)
            ELSE
               test_energy=0.0d0
            ENDIF
            IF(ABS(initial_energy(ibox)+energy_check(ibox)-&
               test_energy) .GT. 0.0000001E0_dp) THEN
               WRITE(iw,*) '!!!!!!! We have an energy problem. !!!!!!!!'
               WRITE(iw,'(A,T64,F16.10)') 'Final Energy = ',test_energy
               WRITE(iw,'(A,T64,F16.10)') 'Inital Energy+energy_check=',&
               initial_energy(ibox)+energy_check(ibox)
               WRITE(6,*) 'Box ',ibox
               WRITE(6,*) 'nchain ',nchain(ibox)
               STOP
            ENDIF
         ENDDO
         ENDIF
      ENDDO

! write a restart file
      IF(ionode) THEN
         DO ibox=1,nboxes
            CALL write_mc_restart(nnstep,mc_par(ibox)%mc_par,&
               force_env(ibox)%force_env)
         ENDDO
      ENDIF

! calculate the final energy
      DO ibox=1,nboxes
         IF(nchain(ibox) .NE. 0) THEN
            CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.)
            CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=final_energy(ibox),&
               error=error)
         ELSE
            final_energy(ibox)=0.0d0
         ENDIF
         IF (lbias) THEN
            CALL force_env_release(bias_env(ibox)%force_env,error=error)
         ENDIF
      ENDDO

! do some stuff in serial
      IF (ionode) THEN

         WRITE(com_ene,*) 'Final Energies:                      ',&
                  final_energy(1:nboxes)

         DO ibox=1,nboxes
            WRITE(cbox,'(I4)') ibox
            IF(nchain(ibox) == 0) THEN
               WRITE(com_crd,*) ' 0'
               WRITE(com_crd,*) 'BOX ' // TRIM(ADJUSTL(cbox))
            ELSE
               CALL write_particle_coordinates(&
               particles_old(ibox)%list%els,&
               cell(ibox)%cell,com_crd,dump_xmol,'POS',&
               'FINAL BOX ' // TRIM(ADJUSTL(cbox)))
            ENDIF

! write a bunch of data to the screen
            WRITE(iw,*) &
               '------------------------------------------------'
            WRITE(iw,'(A,I1,A)') &
               '|                   BOX ',ibox,&
               '                      |'
            WRITE(iw,*) &
               '------------------------------------------------'
            CALL final_mc_write(mc_par(ibox)%mc_par,moves(ibox)%moves,&
               iw,energy_check(ibox),&
               initial_energy(ibox),final_energy(ibox),&
               averages(ibox)%averages,&
               move_updates(ibox)%moves)
            
! close any open files
            CALL close_file(unit_number=diff(ibox))
            CALL close_file(unit_number=data_unit(ibox))
            CALL close_file(unit_number=move_unit(ibox))
            CALL close_file(unit_number=cl(ibox))
            CALL close_file(unit_number=rm(ibox))
         ENDDO

! close some more files
         CALL close_file(unit_number=cell_unit)
         CALL close_file(unit_number=com_ene)
         CALL close_file(unit_number=com_crd)
         CALL close_file(unit_number=com_mol)

! back to parallel
      ENDIF

      DO ibox=1,nboxes
         CALL set_mc_env ( mc_env(ibox)%mc_env, &
            mc_par = mc_par(ibox)%mc_par,&
            force_env=force_env(ibox)%force_env)

! deallocate some stuff
         CALL mc_moves_release(move_updates(ibox)%moves)
         CALL mc_moves_release(moves(ibox)%moves)
         CALL mc_averages_release(averages(ibox)%averages)
      ENDDO

      DEALLOCATE(particle_positions,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions")
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old")
      DEALLOCATE(force_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "force_env")
      DEALLOCATE(bias_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "bias_env")
      DEALLOCATE (cell,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "cell")
      DEALLOCATE (particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particles_old")
      DEALLOCATE (oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "oldsys")
      DEALLOCATE (averages,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "averages")
      DEALLOCATE (moves,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "moves")
      DEALLOCATE (move_updates,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "move_updates")
      DEALLOCATE (mc_par,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mc_par")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_run_ensemble
!*********************************************************************

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
          WRITE (UNIT=cp_logger_get_default_unit_nr(logger),&
             FMT="(/,T2,A,/)")&
             "*** MC run terminated by external request ***"
       END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE external_control

END MODULE mc_ensembles

!**********************************************************************

