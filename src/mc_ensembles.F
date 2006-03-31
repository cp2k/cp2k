!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_ensembles [2.0] *
!!
!!   NAME
!!     mc_ensembles
!!
!!   FUNCTION
!!     Used to run the bulk of the MC simulation, doing things like
!!     choosing move types and writing data to files 
!!
!!   AUTHOR
!!     Matthew J. McGrath  (09.26.2003)
!!
!!   REVISIONS
!!     09.10.05  MJM combined the two subroutines in this module into one
!!
!******************************************************************************


MODULE mc_ensembles
  USE cell_types,                      ONLY: cell_p_type,&
                                             get_cell
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_pp_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE mc_control,                      ONLY: mc_create_bias_force_env,&
                                             write_mc_restart
  USE mc_coordinates,                  ONLY: check_for_overlap,&
                                             create_discrete_array,&
                                             find_mc_test_molecule
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_environment_p_type,&
                                             set_mc_env
  USE mc_ge_moves,                     ONLY: mc_ge_swap_move,&
                                             mc_ge_volume_move,&
                                             mc_quickstep_move
  USE mc_misc,                         ONLY: final_mc_write,&
                                             mc_averages_create,&
                                             mc_averages_release
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_move_update,&
                                             mc_moves_release,&
                                             write_move_stats
  USE mc_moves,                        ONLY: mc_atom_translation,&
                                             mc_avbmc_move,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_averages_p_type,&
                                             mc_moves_p_type,&
                                             mc_simulation_parameters_p_type,&
                                             set_mc_par
  USE message_passing,                 ONLY: mp_bcast
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_p_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: external_control,&
                                             stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_ensembles'

  PUBLIC :: mc_run_ensemble

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_ensembles.F/mc_run_ensemble
!!
!!   NAME
!!     mc_run_ensemble
!!
!!   FUNCTION
!!     directs the program in running one or two box MC simulations
!!
!!   ARGUMENTS
!!     - mc_env: a pointer that contains all mc_env for all the simulation
!!         boxes 
!!     - globenv: the global environment for the simulation
!!     - nboxes: the number of simulation boxes 
!!     - error: the error type
!!     - rng_stream: the stream we pull random numbers from
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
  SUBROUTINE mc_run_ensemble ( mc_env, globenv, nboxes, error,rng_stream )

!   ***************************************************************************

    TYPE(mc_environment_p_type), &
      DIMENSION(:), POINTER                  :: mc_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: nboxes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'mc_run_ensemble'

    CHARACTER(LEN=40)                        :: cbox, fft_lib, move_type, &
                                                move_type_avbmc
    INTEGER, DIMENSION(1:nboxes)             :: box_flag, cl, data_unit, &
                                                diff, istep, move_unit, &
                                                nchain, nunits_tot, rm
    INTEGER, DIMENSION(1:3, 1:2)             :: discrete_array
    INTEGER :: box_number, cell_unit, com_crd, com_ene, com_mol, group, &
      handle, ibox, imove, iparticle, iprint, istat, iuptrans, iupvolume, iw, &
      mol, mol_swap, mol_target, nmoves, nnstep, nstart, nstep, source,&
      print_level
    CHARACTER(LEN=40), DIMENSION(1:nboxes) :: cell_file, coords_file, &
      data_file, displacement_file, energy_file, molecules_file, moves_file
    LOGICAL                                  :: ionode, lbias, ldiscrete, &
                                                lnew_bias_env, loverlap, &
                                                lreject, lstop, should_stop
    REAL(dp), DIMENSION(1:3, 1:2)            :: conf_prob
    REAL(KIND=dp) :: cutoff, discrete_step, inner_cutoff, pbias, pmavbmc, &
      pmswap, pmtraion, pmtrans, pmvolume, rand, test_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r_old
    REAL(KIND=dp), DIMENSION(1:3, 1:nboxes)  :: abc
    REAL(KIND=dp), DIMENSION(1:nboxes)       :: bias_energy, energy_check, &
                                                final_energy, initial_energy, &
                                                old_energy
    TYPE(cell_p_type), DIMENSION(:), POINTER :: cell
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: bias_env, force_env
    TYPE(mc_averages_p_type), DIMENSION(:), &
      POINTER                                :: averages
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: move_updates, moves
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(moves,move_updates)

! allocate a whole bunch of stuff based on how many boxes we have
      ALLOCATE (force_env(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "force_env",1*dp_size*nboxes)
      ALLOCATE (bias_env(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "bias_env",1*dp_size*nboxes)
      ALLOCATE (cell(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell",1*dp_size*nboxes)
      ALLOCATE (particles_old(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_old",1*dp_size*nboxes)
      ALLOCATE (oldsys(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "oldsys",1*dp_size*nboxes)
      ALLOCATE (averages(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "averages",1*dp_size*nboxes)
      ALLOCATE (moves(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "moves",1*dp_size*nboxes)
      ALLOCATE (move_updates(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "move_updates",1*dp_size*nboxes)
      ALLOCATE (mc_par(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
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
         discrete_step=discrete_step,fft_lib=fft_lib,print_level=print_level,&
         conf_prob=conf_prob(:,1),pmavbmc=pmavbmc,pbias=pbias)

      IF(nboxes .GT. 1) THEN
         DO ibox=2,nboxes
            CALL get_mc_par(mc_par(ibox)%mc_par,nunits_tot=nunits_tot(ibox),&
            nchain=nchain(ibox),data_file=data_file(ibox),&
            moves_file=moves_file(ibox),&
            cell_file=cell_file(ibox),coords_file=coords_file(ibox),&
            energy_file=energy_file(ibox),&
            displacement_file=displacement_file(ibox),&
            molecules_file=molecules_file(ibox),conf_prob=conf_prob(:,2))
         ENDDO
      ENDIF

! allocate the particle positions array for broadcasting
      ALLOCATE (r_old(3,SUM(nunits_tot),1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_old",3*SUM(nunits_tot)*dp_size*nboxes)


! figure out what the default write unit is
      iw = globenv%scr

      IF (ionode) THEN
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
               calc_force=.FALSE.,error=error)
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
               globenv,.FALSE.,TRIM(fft_lib),print_level,ionode,error)
            IF(nchain(ibox) .NE. 0) THEN
               CALL force_env_calc_energy_force(bias_env(ibox)%force_env,&
                  calc_force=.FALSE.,error=error)
               CALL force_env_get(bias_env(ibox)%force_env,&
                  potential_energy=moves(ibox)%moves%bias_energy,error=error)
                  
            ELSE
               moves(ibox)%moves%bias_energy=0.0d0
            ENDIF
            bias_energy(ibox)=moves(ibox)%moves%bias_energy
         ENDIF
         lnew_bias_env=.FALSE.
         CALL check_for_overlap(force_env(ibox)%force_env,nchain(ibox),1,&
            nchain(ibox),loverlap)
         IF(loverlap) CALL stop_program(routineN,moduleN,__LINE__,&
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

         ENDDO

! back to parallel mode
      ENDIF

      DO ibox=1,nboxes
         CALL mp_bcast(cl(ibox),source,group)
         CALL mp_bcast(rm(ibox),source,group)
         CALL mp_bcast(diff(ibox),source,group)
! set all the units numbers that we just opened in the respective mc_par
         CALL set_mc_par(mc_par(ibox)%mc_par,cl=cl(ibox),rm=rm(ibox),&
         diff=diff(ibox))
      ENDDO

! if we're doing a discrete volume move, we need to set up the array
! that keeps track of which direction we can move in
      IF(ldiscrete) THEN
         IF(nboxes .NE. 1) &
            CALL stop_program(routineN,moduleN,__LINE__,&
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

         IF (ionode) rand=next_random_number(rng_stream)
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
                  rng_stream,error=error)

            ELSE
               CALL mc_ge_volume_move ( mc_par,force_env, moves,&
                  move_updates,nnstep,old_energy,globenv,energy_check,&
                  r_old,rng_stream,error)
            ENDIF

! we need a new biasing environment now, if we're into that sort of thing
            IF(lbias) THEN
               DO ibox=1,nboxes
                  CALL force_env_release(bias_env(ibox)%force_env,error=error)
                  CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
                     force_env(ibox)%force_env,&
                     globenv,.FALSE.,TRIM(fft_lib),print_level,ionode,error)
                  IF(nchain(ibox) .NE. 0) THEN
                     CALL force_env_calc_energy_force(&
                        bias_env(ibox)%force_env,&
                        calc_force=.FALSE.,error=error)
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
                globenv,bias_energy(:),rng_stream,error)

! the number of molecules may have changed
            DO ibox=1,nboxes
               CALL get_mc_par(mc_par(ibox)%mc_par,&
                  nunits_tot=nunits_tot(ibox),&
                  nchain=nchain(ibox))
            ENDDO
         ELSEIF ( rand .LT. pmavbmc) THEN
! try an AVBMC move
            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting an AVBMC1 move"
               WRITE(iw,*) 
            ENDIF

! first, pick a box to do it for
            IF(ionode) rand=next_random_number(rng_stream)
            CALL mp_bcast(rand,source,group)

            IF(nboxes .EQ. 2) THEN
               IF(rand .LT. 0.25E0_dp) THEN
                  ibox=1
               ELSE
                  ibox=2
               ENDIF
            ELSE
               ibox=1
            ENDIF

! now pick a molecule, automatically rejecting the move if the
! box is empty or only has one molecule
            IF(nchain(ibox) .LE. 1) THEN
! indicate that we tried a move
               moves(ibox)%moves%empty_avbmc=&
                  moves(ibox)%moves%empty_avbmc+1
            ELSE

! pick a molecule to be swapped in the box
               IF (ionode) THEN
                  CALL find_mc_test_molecule(force_env,mol_swap,&
                     box_number,nchain(:),rng_stream,box=ibox,error=error)

! pick a molecule to act as the target in the box
                  DO
                     CALL find_mc_test_molecule(force_env,mol_target,&
                     box_number,nchain(:),rng_stream,box=ibox,error=error)
                     IF(mol_target .NE. mol_swap) EXIT
                  ENDDO

! choose if we're swapping into the bonded region of mol_target, or 
! into the nonbonded region
                  rand=next_random_number(rng_stream)

               ENDIF
               CALL mp_bcast(mol_swap,source,group)
               CALL mp_bcast(box_number,source,group)
               CALL mp_bcast(mol_target,source,group)
               CALL mp_bcast(rand,source,group)

               IF (rand .LT. pbias) THEN
                  move_type_avbmc='in'
               ELSE
                  move_type_avbmc='out'
               ENDIF

               CALL mc_avbmc_move(mc_par(box_number)%mc_par,&
                  force_env(box_number)%force_env,&
                  bias_env(box_number)%force_env,&
                  moves(box_number)%moves,energy_check(box_number),&
                  r_old(:,:,box_number),old_energy(box_number),&
                  mol_swap,mol_target,bias_energy(box_number),&
                  move_type_avbmc,rng_stream,error=error)

            ENDIF

         ELSE

            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting an inner move"
               WRITE(iw,*) 
            ENDIF

            DO imove=1,nmoves
                  
               IF(ionode) rand=next_random_number(rng_stream)
               CALL mp_bcast(rand,source,group)

               IF(rand .LT. pmtraion) THEN
! change molecular conformation
! first, pick a box to do it for
                  IF(ionode) rand=next_random_number(rng_stream)
                  CALL mp_bcast(rand,source,group)
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
                     IF(ionode) THEN
                        CALL find_mc_test_molecule(force_env,mol,&
                           box_number,nchain(:),rng_stream,&
                           box=ibox,error=error)
                        
! choose if we're changing a bond length or an angle
                        rand=next_random_number(rng_stream)
                     ENDIF
                     CALL mp_bcast(rand,source,group)
                     CALL mp_bcast(mol,source,group)
                     CALL mp_bcast(box_number,source,group)

! figure out what kind of move we're doing
                     IF (rand .LT. conf_prob(1,ibox)) THEN
                        move_type='bond'
                     ELSEIF(rand .LT. (conf_prob(1,ibox)+&
                        conf_prob(2,ibox))) THEN
                        move_type='angle'
                     ELSE
                        move_type='dihedral'
                     ENDIF
                     box_flag(box_number)=1
                     CALL mc_atom_translation(mc_par(box_number)%mc_par,&
                        force_env(box_number)%force_env,&
                        bias_env(box_number)%force_env,&
                        moves(box_number)%moves,&
                        move_updates(box_number)%moves,&
                        mol,bias_energy(box_number),&
                        move_type,lreject,rng_stream,error=error)
                     IF(lreject) EXIT
                  ENDIF
               ELSEIF(rand .LT. pmtrans) THEN
! translate a whole molecule in the system
                  IF(ionode) &
                     CALL find_mc_test_molecule(force_env,mol,box_number,&
                     nchain(:),rng_stream,error=error)
                  CALL mp_bcast(mol,source,group)
                  CALL mp_bcast(box_number,source,group)                     
                  box_flag(box_number)=1
                  CALL mc_molecule_translation(mc_par(box_number)%mc_par,&
                     force_env(box_number)%force_env,&
                     bias_env(box_number)%force_env,&
                     moves(box_number)%moves,&
                     move_updates(box_number)%moves,mol,&
                     bias_energy(box_number),lreject,rng_stream,error=error)
                  IF(lreject) EXIT
               ELSE
!     rotate a whole molecule in the system
                  IF(ionode) &
                     CALL find_mc_test_molecule(force_env,mol,box_number,&
                     nchain(:),rng_stream,error=error)
                  CALL mp_bcast(mol,source,group)
                  CALL mp_bcast(box_number,source,group)  
                  box_flag(box_number)=1
                  CALL mc_molecule_rotation(mc_par(box_number)%mc_par,&
                     force_env(box_number)%force_env,&
                     bias_env(box_number)%force_env,&
                     moves(box_number)%moves,&
                     move_updates(box_number)%moves,&
                     mol,bias_energy(box_number),&
                     lreject,rng_stream,error=error)
                  IF(lreject) EXIT
               ENDIF
               
            ENDDO

! now do a Quickstep calculation to see if we accept the sequence
            CALL mc_Quickstep_move(mc_par,force_env,bias_env,&
               moves,lreject,move_updates,energy_check(:),r_old(:,:,:),&
               nnstep,old_energy(:),bias_energy(:),&
               nboxes,box_flag(:),oldsys,particles_old,cell,&
               rng_stream,error=error)

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
                     nnstep,move_unit(ibox))

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

! find bond lengths and angle for current water config
!               IF(nchain(ibox) .NE. 0) THEN
!                  istep(ibox)=istep(ibox)+1
!                  CALL get_water_data(force_env(ibox)%force_env,&
!                     averages(ibox)%averages,&
!                     istep(ibox),nchain(ibox),data_unit(ibox))
!               ENDIF

               ENDDO
            ENDIF  ! end the things we only do every iprint moves

            DO ibox=1,nboxes
! compute some averages
               averages(ibox)%averages%ave_energy=&
                  averages(ibox)%averages%ave_energy*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  old_energy(ibox)/REAL(nnstep-nstart,dp)
               averages(ibox)%averages%molecules=&
                  averages(ibox)%averages%molecules*REAL(nnstep-&
                  nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  REAL(nchain(ibox),dp)/REAL(nnstep-nstart,dp)
               averages(ibox)%averages%ave_volume=&
                  averages(ibox)%averages%ave_volume*&
                  REAL(nnstep-nstart-1,dp)/REAL(nnstep-nstart,dp)+&
                  abc(1,ibox)*abc(2,ibox)*abc(3,ibox)/&
                  REAL(nnstep-nstart,dp)

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

         ENDIF

! reset the box flags
         box_flag(:)=0

! check to see if EXIT file exists...if so, end the calculation
         CALL external_control(should_stop,"MC",globenv,error)
         IF (should_stop) EXIT

! update the move displacements, if necessary
         DO ibox=1,nboxes
            IF (MOD(nnstep-nstart,iuptrans) == 0) THEN
               CALL mc_move_update(mc_par(ibox)%mc_par,&
                  move_updates(ibox)%moves,"trans",nnstep,ionode)
            ENDIF
               
            IF (MOD(nnstep-nstart,iupvolume) == 0  ) THEN
               CALL mc_move_update(mc_par(ibox)%mc_par,&
                  move_updates(ibox)%moves,"volume",nnstep,ionode)
            ENDIF
         ENDDO

! check to see if there are any overlaps in the boxe
         DO ibox=1,nboxes
            IF(nchain(ibox) .NE. 0) THEN
               CALL check_for_overlap(force_env(ibox)%force_env,&
                  nchain(ibox),1,nchain(ibox),loverlap)
               IF(loverlap) THEN
                  IF(ionode) WRITE(iw,*) nnstep
                  CALL stop_program(routineN,moduleN,__LINE__,&
                  'coordinate overlap at the end of the above step')
               ENDIF
            ENDIF
         ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug make sure the energies line up
         IF(.FALSE.) THEN
         DO ibox=1,nboxes
            IF(nchain(ibox) .NE. 0) THEN
               CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.,error=error)
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
               STOP 'mc_ensembles.F'
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
               calc_force=.FALSE.,error=error)
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
               averages(ibox)%averages)
            
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

      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_old")
      DEALLOCATE(force_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "force_env")
      DEALLOCATE(bias_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "bias_env")
      DEALLOCATE (cell,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "cell")
      DEALLOCATE (particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "particles_old")
      DEALLOCATE (oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "oldsys")
      DEALLOCATE (averages,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "averages")
      DEALLOCATE (moves,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "moves")
      DEALLOCATE (move_updates,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "move_updates")
      DEALLOCATE (mc_par,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "mc_par")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_run_ensemble
!*********************************************************************

END MODULE mc_ensembles

!**********************************************************************

