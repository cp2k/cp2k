!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_pp_type,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol
  USE input_section_types,             ONLY: section_vals_release
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             default_string_length
  USE machine,                         ONLY: m_flush
  USE mc_control,                      ONLY: mc_create_bias_force_env,&
                                             write_mc_restart
  USE mc_coordinates,                  ONLY: check_for_overlap,&
                                             create_discrete_array,&
                                             find_mc_test_molecule,&
                                             mc_coordinate_fold
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
  USE mc_moves,                        ONLY: mc_conformation_change,&
                                             mc_avbmc_move,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_averages_p_type,&
                                             mc_moves_p_type,&
                                             mc_simulation_parameters_p_type,&
                                             set_mc_par,&
                                             mc_molecule_info_type,&
                                             get_mc_molecule_info,&
                                             mc_input_file_type
  USE message_passing,                 ONLY: mp_bcast
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
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
  SUBROUTINE mc_run_ensemble ( mc_env, para_env, globenv, nboxes, error, rng_stream)

!   ***************************************************************************

    TYPE(mc_environment_p_type), &
      DIMENSION(:), POINTER                  :: mc_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: nboxes
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'mc_run_ensemble'

    CHARACTER(LEN=20)                        :: ensemble
    CHARACTER(LEN=40)                        :: cbox, fft_lib, move_type, &
                                                move_type_avbmc
    INTEGER, DIMENSION(1:nboxes)             :: box_flag, cl, data_unit, &
                                                diff, istep, move_unit, &
                                                rm
    INTEGER, DIMENSION(1:3, 1:2)             :: discrete_array
    INTEGER :: box_number, cell_unit, com_crd, com_ene, com_mol, group, &
      handle, ibox, imove, iparticle, iprint, istat, iuptrans, iupvolume, iw, &
      start_atom, start_atom_swap, start_atom_target, nmoves, nnstep, nstart,&
      nstep, print_level, jbox,start_mol,end_mol,idum,imol_type,&
      source, molecule_type,nmol_types,itype,nchain_total,&
      molecule_type_swap,molecule_type_target,iunit,imolecule,&
      atom_number
    CHARACTER(LEN=40), DIMENSION(1:nboxes) :: cell_file, coords_file, &
      data_file, displacement_file, energy_file, molecules_file, moves_file
    LOGICAL                                  :: ionode, lbias, ldiscrete, &
                                                lnew_bias_env, loverlap, &
                                                lreject, lstop, should_stop
    REAL(dp), DIMENSION(:),POINTER            :: pmtraion_mol,pmavbmc_mol,&
         pbias,pmrot_mol,pmtrans_mol,pmvol_box
    REAL(dp), DIMENSION(:, :),POINTER            :: conf_prob, mass
    INTEGER, DIMENSION(:, :),POINTER            :: nchains
    INTEGER, DIMENSION(:),POINTER            :: mol_type,avbmc_atom,nunits,&
         nunits_tot,nchains_box
    REAL(KIND=dp)                            :: discrete_step, &
                                                pmavbmc, pmswap, pmtraion, &
                                                pmtrans, pmvolume, rand, &
                                                test_energy
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                     :: r_temp
    REAL(KIND=dp), DIMENSION(1:3, 1:nboxes)  :: abc
    REAL(KIND=dp), DIMENSION(1:nboxes)       :: bias_energy, energy_check, &
                                                final_energy, initial_energy, &
                                                old_energy, last_bias_energy
    TYPE(cell_p_type), DIMENSION(:), POINTER :: cell
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: biassys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: bias_env, force_env
    TYPE(mc_averages_p_type), DIMENSION(:), &
      POINTER                                :: averages
    TYPE(mc_moves_p_type), DIMENSION(:,:), &
      POINTER                                :: move_updates, moves
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: test_moves
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_type),POINTER         :: particles_bias
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old
    TYPE(mc_molecule_info_type),POINTER      :: mc_molecule_info
    CHARACTER(default_string_length),DIMENSION(:),ALLOCATABLE :: atom_names_box
    CHARACTER(default_string_length),DIMENSION(:,:),POINTER:: atom_names
    TYPE(mc_input_file_type),POINTER :: mc_bias_file
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(moves,move_updates,test_moves)

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
      ALLOCATE (mc_par(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par",1*dp_size*nboxes)
      ALLOCATE (pmvol_box(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "pmvol_box",dp_size*nboxes)

      DO ibox=1,nboxes
         CALL get_mc_env ( mc_env(ibox)%mc_env, &
            mc_par = mc_par(ibox)%mc_par,&
            force_env=force_env(ibox)%force_env)
      ENDDO

! get some data out of mc_par
      CALL get_mc_par(mc_par(1)%mc_par,&
         ionode=ionode,source=source,group=group,&
         data_file=data_file(1),moves_file=moves_file(1),&
         cell_file=cell_file(1),coords_file=coords_file(1),&
         energy_file=energy_file(1),displacement_file=displacement_file(1),&
         lstop=lstop,nstep=nstep,nstart=nstart,pmvolume=pmvolume,&
         molecules_file=molecules_file(1),pmswap=pmswap,nmoves=nmoves,&
         pmtraion=pmtraion,pmtrans=pmtrans,iuptrans=iuptrans,&
         iupvolume=iupvolume,ldiscrete=ldiscrete,pmtraion_mol=pmtraion_mol,&
         lbias=lbias,iprint=iprint,pmavbmc_mol=pmavbmc_mol,&
         discrete_step=discrete_step,fft_lib=fft_lib,avbmc_atom=avbmc_atom,&
         pmavbmc=pmavbmc,pbias=pbias,mc_molecule_info=mc_molecule_info,&
         pmrot_mol=pmrot_mol,pmtrans_mol=pmtrans_mol,pmvol_box=pmvol_box(1),&
         ensemble=ensemble)

! get some data from the molecule types
      CALL get_mc_molecule_info(mc_molecule_info,conf_prob=conf_prob,&
           nchains=nchains,nmol_types=nmol_types,nunits_tot=nunits_tot,&
           mol_type=mol_type,nchain_total=nchain_total,nunits=nunits,&
           atom_names=atom_names,mass=mass)

! allocate some stuff based on the number of molecule types we have
      ALLOCATE (moves(1:nmol_types,1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "moves",dp_size*nboxes*nmol_types)
      ALLOCATE (move_updates(1:nmol_types,1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "move_updates",nmol_types*dp_size*nboxes)

      print_level = 1 ! hack

      IF(nboxes .GT. 1) THEN
         DO ibox=2,nboxes
            CALL get_mc_par(mc_par(ibox)%mc_par,&
            data_file=data_file(ibox),&
            moves_file=moves_file(ibox),&
            cell_file=cell_file(ibox),coords_file=coords_file(ibox),&
            energy_file=energy_file(ibox),&
            displacement_file=displacement_file(ibox),&
            molecules_file=molecules_file(ibox),pmvol_box=pmvol_box(ibox))
         ENDDO
      ENDIF

! this is a check we can't do in the input checking
    IF(pmvol_box(nboxes) .LT. 0.9999E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMVOL_BOX needs to be 1.0')       
    ENDIF

! allocate the particle positions array for broadcasting
      ALLOCATE (r_old(3,SUM(nunits_tot),1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_old",3*SUM(nunits_tot)*dp_size*nboxes)


! figure out what the default write unit is
      iw = cp_logger_get_default_unit_nr()

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
         DO itype=1,nmol_types
            CALL init_mc_moves(moves(itype,ibox)%moves)
            CALL init_mc_moves(move_updates(itype,ibox)%moves)
         ENDDO
         CALL mc_averages_create(averages(ibox)%averages)

! find the energy of the initial configuration
         IF(SUM(nchains(:,ibox)) .NE. 0) THEN
            CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.,error=error)
            CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=old_energy(ibox),&
               error=error)
         ELSE
            old_energy(ibox)=0.0d0
         ENDIF
         initial_energy(ibox)=old_energy(ibox)

! check for overlaps
         start_mol=1
         DO jbox=1,ibox-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO
         end_mol=start_mol+SUM(nchains(:,ibox))-1
         CALL check_for_overlap(force_env(ibox)%force_env,nchains(:,ibox),&
              nunits,loverlap,mol_type(start_mol:end_mol))
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

! find the bias energy of the initial run
         IF(lbias) THEN
! determine the atom names of every particle
            ALLOCATE (atom_names_box(1:nunits_tot(ibox)),STAT=istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 "atom_names_box",default_string_length*nunits_tot(ibox))
            atom_number=1
            DO imolecule=1,SUM(nchains(:,ibox))
               DO iunit=1,nunits(mol_type(imolecule+start_mol-1))
                  atom_names_box(atom_number)=&
                       atom_names(iunit,mol_type(imolecule+start_mol-1))
                  atom_number=atom_number+1
               ENDDO
            ENDDO

            CALL get_mc_par(mc_par(ibox)%mc_par,mc_bias_file=mc_bias_file)
            nchains_box=> nchains(:,ibox)
            CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
                 r_old(:,:,ibox),atom_names_box(:),nunits_tot(ibox),&
                 para_env,abc(:,ibox),nchains_box,mc_bias_file,&
                 ionode,error)
            IF(SUM(nchains(:,ibox)) .NE. 0) THEN
               CALL force_env_calc_energy_force(bias_env(ibox)%force_env,&
                  calc_force=.FALSE.,error=error)
               CALL force_env_get(bias_env(ibox)%force_env,&
                  potential_energy=last_bias_energy(ibox),error=error)

            ELSE
               last_bias_energy(ibox)=0.0d0
            ENDIF
            bias_energy(ibox)=last_bias_energy(ibox)
            DEALLOCATE(atom_names_box,STAT=istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 "atom_names_box")

         ENDIF
         lnew_bias_env=.FALSE.

      ENDDO

! back to seriel for a bunch of I/O stuff
      IF (ionode) THEN

! record the combined energies,coordinates, and cell lengths
         CALL open_file(file_name='mc_cell_length',&
                 unit_number=cell_unit,file_position='APPEND',&
                 file_action='WRITE',file_status='UNKNOWN')
         CALL open_file(file_name='mc_energies',&
                 unit_number=com_ene,file_position='APPEND',&
                 file_action='WRITE',file_status='UNKNOWN')
         CALL open_file(file_name='mc_coordinates',&
                 unit_number=com_crd,file_position='APPEND',&
                 file_action='WRITE',file_status='UNKNOWN')
         CALL open_file(file_name='mc_molecules',&
                 unit_number=com_mol,file_position='APPEND',&
                 file_action='WRITE',file_status='UNKNOWN')
         WRITE(com_ene,*) 'Initial Energies:       ',&
            old_energy(1:nboxes)
         DO ibox=1,nboxes
            WRITE(com_mol,*) 'Initial Molecules:       ',&
                 nchains(:,ibox)
         ENDDO
         DO ibox=1,nboxes
            WRITE(cell_unit,*) 'Initial: ',&
               abc(1:3,ibox)*angstrom
            WRITE(cbox,'(I4)') ibox
            CALL open_file(file_name='energy_differences_box' // &
                 TRIM(ADJUSTL(cbox)),&
                 unit_number=diff(ibox),file_position='APPEND',&
                 file_action='WRITE',file_status='UNKNOWN')
            IF(SUM(nchains(:,ibox)) == 0) THEN
               WRITE(com_crd,*) ' 0'
               WRITE(com_crd,*) 'INITIAL BOX ' // TRIM(ADJUSTL(cbox))
            ELSE
               CALL write_particle_coordinates(particles_old(ibox)%list%els,&
                  com_crd,dump_xmol,'POS',&
                  'INITIAL BOX ' // TRIM(ADJUSTL(cbox)))
            ENDIF
            CALL open_file(file_name=data_file(ibox),&
               unit_number=data_unit(ibox),file_position='APPEND',&
               file_action='WRITE',file_status='UNKNOWN')
            CALL open_file(file_name=moves_file(ibox),&
               unit_number=move_unit(ibox),file_position='APPEND',&
               file_action='WRITE',file_status='UNKNOWN')
            CALL open_file(file_name=displacement_file(ibox),&
               unit_number=rm(ibox),file_position='APPEND',&
               file_action='WRITE',file_status='UNKNOWN')
            CALL open_file(file_name=cell_file(ibox),&
               unit_number=cl(ibox),file_position='APPEND',&
               file_action='WRITE',file_status='UNKNOWN')

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
         nstep=nstep*nchain_total
         iuptrans=iuptrans*nchain_total
         iupvolume=iupvolume*nchain_total
      ENDIF

      DO nnstep=nstart+1,nstart + nstep

         IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
            WRITE(iw,*)
            WRITE(iw,*) "------- On Monte Carlo Step ",nnstep
         ENDIF

         IF (ionode) rand=next_random_number(rng_stream,error=error)
! broadcast the random number, to make sure we're on the same move
         CALL mp_bcast(rand,source,group)

         IF (rand .LT. pmvolume) THEN

            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting a volume move"
               WRITE(iw,*)
            ENDIF
            
            SELECT CASE(ensemble)
            CASE("TRADITIONAL")
               CALL mc_volume_move( mc_par(1)%mc_par,&
                    force_env(1)%force_env,&
                    moves(1,1)%moves,move_updates(1,1)%moves,&
                    nnstep,old_energy(1),1,&
                    energy_check(1),r_old(:,:,1),iw,discrete_array(:,:),&
                    rng_stream,error=error)
            CASE("GEMC_NVT")
               CALL mc_ge_volume_move ( mc_par,force_env, moves,&
                    move_updates,nnstep,old_energy,energy_check,&
                    r_old,rng_stream,error)
            CASE("GEMC_NPT")
! we need to select a box based on the probability given in the input file
               IF (ionode) rand=next_random_number(rng_stream,error=error)
               CALL mp_bcast(rand,source,group)               

               DO ibox=1,nboxes
                  IF(rand .LE. pmvol_box(ibox)) THEN
                     box_number=ibox
                     EXIT
                  ENDIF
               ENDDO
               
               CALL mc_volume_move( mc_par(box_number)%mc_par,&
                    force_env(box_number)%force_env,&
                    moves(1,box_number)%moves,&
                    move_updates(1,box_number)%moves,&
                    nnstep,old_energy(box_number),box_number,&
                    energy_check(box_number),r_old(:,:,box_number),iw,&
                    discrete_array(:,:),&
                    rng_stream,error=error)
            END SELECT
            

! we need a new biasing environment now, if we're into that sort of thing
            IF(lbias) THEN
               DO ibox=1,nboxes
                  CALL section_vals_release(bias_env(ibox)%force_env%root_section,&
                       error=error)
                  CALL force_env_release(bias_env(ibox)%force_env,error=error)
! determine the atom names of every particle
                  ALLOCATE (atom_names_box(1:nunits_tot(ibox)),STAT=istat)
                  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_names_box",default_string_length*nunits_tot(ibox))
                  start_mol=1
                  DO jbox=1,ibox-1
                     start_mol=start_mol+SUM(nchains(:,jbox))
                  ENDDO
                  end_mol=start_mol+SUM(nchains(:,ibox))-1
                  atom_number=1
                  DO imolecule=1,SUM(nchains(:,ibox))
                     DO iunit=1,nunits(mol_type(imolecule+start_mol-1))
                        atom_names_box(atom_number)=&
                             atom_names(iunit,mol_type(imolecule+start_mol-1))
                        atom_number=atom_number+1
                     ENDDO
                  ENDDO
                  
                  CALL get_mc_par(mc_par(ibox)%mc_par,&
                       mc_bias_file=mc_bias_file)
                  nchains_box=> nchains(:,ibox)
                  CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
                       r_old(:,:,ibox),atom_names_box(:),nunits_tot(ibox),&
                       para_env,abc(:,ibox),nchains_box,mc_bias_file,&
                       ionode,error)
                  
                  IF(SUM(nchains(:,ibox)) .NE. 0) THEN
                     CALL force_env_calc_energy_force(&
                          bias_env(ibox)%force_env,&
                          calc_force=.FALSE.,error=error)
                     CALL force_env_get(bias_env(ibox)%force_env,&
                          potential_energy=last_bias_energy(ibox),&
                          error=error)
                  ELSE
                     last_bias_energy(ibox)=0.0d0
                  ENDIF
                  bias_energy(ibox)=last_bias_energy(ibox)
                  DEALLOCATE(atom_names_box,STAT=istat)
                  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "atom_names_box")
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
                 para_env,bias_energy(:),last_bias_energy(:),&
                 rng_stream,error)
            
! the number of molecules may have changed, which deallocated the whole
! mc_molecule_info structure
            CALL get_mc_par(mc_par(1)%mc_par,mc_molecule_info=mc_molecule_info)
            CALL get_mc_molecule_info(mc_molecule_info,conf_prob=conf_prob,&
                 nchains=nchains,nmol_types=nmol_types,nunits_tot=nunits_tot,&
                 mol_type=mol_type,nchain_total=nchain_total,nunits=nunits,&
                 atom_names=atom_names,mass=mass)

         ELSEIF ( rand .LT. pmavbmc) THEN
! try an AVBMC move
            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting an AVBMC1 move"
               WRITE(iw,*)
            ENDIF

! first, pick a box to do it for
            IF(ionode) rand=next_random_number(rng_stream,error=error)
            CALL mp_bcast(rand,source,group)

            IF(nboxes .EQ. 2) THEN
               IF(rand .LT. 0.1E0_dp) THEN
                  box_number=1
               ELSE
                  box_number=2
               ENDIF
            ELSE
               box_number=1
            ENDIF

! now pick a molecule type to do it for
            IF(ionode) rand=next_random_number(rng_stream,error=error)
            CALL mp_bcast(rand,source,group)
            molecule_type_swap=0
            DO imol_type=1,nmol_types
               IF(rand .LT. pmavbmc_mol(imol_type)) THEN
                  molecule_type_swap=imol_type
                  EXIT
               ENDIF
            ENDDO
            IF(molecule_type_swap == 0) CALL &
                 stop_program(routineN,moduleN,__LINE__,&
                  'Did not choose a molecule type to swap...check AVBMC input')

! now pick a molecule, automatically rejecting the move if the
! box is empty or only has one molecule
            IF(SUM(nchains(:,box_number)) .LE. 1) THEN
! indicate that we tried a move
               moves(molecule_type_swap,box_number)%moves%empty_avbmc=&
                  moves(molecule_type_swap,box_number)%moves%empty_avbmc+1
            ELSE

! pick a molecule to be swapped in the box
               IF (ionode) THEN
                  CALL find_mc_test_molecule(mc_molecule_info,&
                       start_atom_swap,idum,idum,rng_stream,&
                       box=box_number,molecule_type_old=molecule_type_swap,&
                       error=error)

! pick a molecule to act as the target in the box...we don't care what type
                  DO
                     CALL find_mc_test_molecule(mc_molecule_info,&
                          start_atom_target,idum,molecule_type_target,&
                          rng_stream,box=box_number,error=error)
                     IF(start_atom_swap .NE. start_atom_target) THEN
                        start_atom_target=start_atom_target+&
                             avbmc_atom(molecule_type_target)-1
                        EXIT
                     ENDIF
                  ENDDO

! choose if we're swapping into the bonded region of mol_target, or
! into the nonbonded region
                  rand=next_random_number(rng_stream,error=error)

               ENDIF
               CALL mp_bcast(start_atom_swap,source,group)
               CALL mp_bcast(box_number,source,group)
               CALL mp_bcast(start_atom_target,source,group)
               CALL mp_bcast(rand,source,group)

               IF (rand .LT. pbias(molecule_type_swap)) THEN
                  move_type_avbmc='in'
               ELSE
                  move_type_avbmc='out'
               ENDIF

               CALL mc_avbmc_move(mc_par(box_number)%mc_par,&
                  force_env(box_number)%force_env,&
                  bias_env(box_number)%force_env,&
                  moves(molecule_type_swap,box_number)%moves,&
                  energy_check(box_number),&
                  r_old(:,:,box_number),old_energy(box_number),&
                  start_atom_swap,start_atom_target,molecule_type_swap,&
                  box_number,bias_energy(box_number),&
                  last_bias_energy(box_number),&
                  move_type_avbmc,rng_stream,error=error)

            ENDIF

         ELSE

            IF(MOD(nnstep,iprint) == 0 .AND. ionode) THEN
               WRITE(iw,*) "Attempting an inner move"
               WRITE(iw,*)
            ENDIF

            DO imove=1,nmoves

               IF(ionode) rand=next_random_number(rng_stream,error=error)
               CALL mp_bcast(rand,source,group)
               IF(rand .LT. pmtraion) THEN
! change molecular conformation
! first, pick a box to do it for
                  IF(ionode) rand=next_random_number(rng_stream,error=error)
                  CALL mp_bcast(rand,source,group)
                  IF(nboxes .EQ. 2) THEN
                     IF(rand .LT. 0.75E0_dp) THEN
                        box_number=1
                     ELSE
                        box_number=2
                     ENDIF
                  ELSE
                     box_number=1
                  ENDIF

! figure out which molecule type we're looking for
                  IF(ionode) rand=next_random_number(rng_stream,error=error)
                  CALL mp_bcast(rand,source,group)
                  molecule_type=0
                  DO imol_type=1,nmol_types
                     IF(rand .LT. pmtraion_mol(imol_type)) THEN
                        molecule_type=imol_type
                        EXIT
                     ENDIF
                  ENDDO
                  IF(molecule_type == 0) CALL &
                       stop_program(routineN,moduleN,__LINE__,&
    'Did not choose a molecule type to conf change...PMTRAION_MOL should not be all 0.0')

! now pick a molecule, automatically rejecting the move if the
! box is empty
                  IF(nchains(molecule_type,box_number) == 0) THEN
! indicate that we tried a move
                     moves(molecule_type,box_number)%moves%empty_conf=&
                        moves(molecule_type,box_number)%moves%empty_conf+1
                  ELSE
! pick a molecule in the box
                     IF(ionode) THEN
                        CALL find_mc_test_molecule(mc_molecule_info,&
                             start_atom,idum,&
                             idum,rng_stream,&
                             box=box_number,molecule_type_old=molecule_type,&
                             error=error)

! choose if we're changing a bond length or an angle
                        rand=next_random_number(rng_stream,error=error)
                     ENDIF
                     CALL mp_bcast(rand,source,group)
                     CALL mp_bcast(start_atom,source,group)
                     CALL mp_bcast(box_number,source,group)
                     CALL mp_bcast(molecule_type,source,group)

! figure out what kind of move we're doing
                     IF (rand .LT. conf_prob(1,molecule_type)) THEN
                        move_type='bond'
                     ELSEIF(rand .LT. (conf_prob(1,molecule_type)+&
                        conf_prob(2,molecule_type))) THEN
                        move_type='angle'
                     ELSE
                        move_type='dihedral'
                     ENDIF
                     box_flag(box_number)=1
                     CALL mc_conformation_change(mc_par(box_number)%mc_par,&
                        force_env(box_number)%force_env,&
                        bias_env(box_number)%force_env,&
                        moves(molecule_type,box_number)%moves,&
                        move_updates(molecule_type,box_number)%moves,&
                        start_atom,molecule_type,box_number,&
                        bias_energy(box_number),&
                        move_type,lreject,rng_stream,error=error)
                     IF(lreject) EXIT
                  ENDIF
               ELSEIF(rand .LT. pmtrans) THEN
! translate a whole molecule in the system
! pick a molecule type
                  IF(ionode) rand=next_random_number(rng_stream,error=error)
                  CALL mp_bcast(rand,source,group)
                  molecule_type=0
                  DO imol_type=1,nmol_types
                     IF(rand .LT. pmtrans_mol(imol_type)) THEN
                        molecule_type=imol_type
                        EXIT
                     ENDIF
                  ENDDO
                  IF(molecule_type == 0) CALL &
                       stop_program(routineN,moduleN,__LINE__,&
    'Did not choose a molecule type to translate...PMTRANS_MOL should not be all 0.0')

! now pick a molecule of that type                  
                  IF(ionode) &
                     CALL find_mc_test_molecule(mc_molecule_info,&
                          start_atom,box_number,idum,rng_stream,&
                          molecule_type_old=molecule_type,error=error)
                  CALL mp_bcast(start_atom,source,group)
                  CALL mp_bcast(box_number,source,group)
                  box_flag(box_number)=1
                  CALL mc_molecule_translation(mc_par(box_number)%mc_par,&
                     force_env(box_number)%force_env,&
                     bias_env(box_number)%force_env,&
                     moves(molecule_type,box_number)%moves,&
                     move_updates(molecule_type,box_number)%moves,&
                     start_atom,box_number,bias_energy(box_number),&
                     molecule_type,lreject,rng_stream,error=error)
                  IF(lreject) EXIT
               ELSE
!     rotate a whole molecule in the system
! pick a molecule type
                  IF(ionode) rand=next_random_number(rng_stream,error=error)
                  CALL mp_bcast(rand,source,group)
                  molecule_type=0
                  DO imol_type=1,nmol_types
                     IF(rand .LT. pmrot_mol(imol_type)) THEN
                        molecule_type=imol_type
                        EXIT
                     ENDIF
                  ENDDO
                  IF(molecule_type == 0) CALL &
                       stop_program(routineN,moduleN,__LINE__,&
    'Did not choose a molecule type to rotate...PMROT_MOL should not be all 0.0')

                  IF(ionode) &
                     CALL find_mc_test_molecule(mc_molecule_info,&
                          start_atom,box_number,idum,rng_stream,&
                          molecule_type_old=molecule_type,error=error)
                  CALL mp_bcast(start_atom,source,group)
                  CALL mp_bcast(box_number,source,group)
                  box_flag(box_number)=1
                  CALL mc_molecule_rotation(mc_par(box_number)%mc_par,&
                     force_env(box_number)%force_env,&
                     bias_env(box_number)%force_env,&
                     moves(molecule_type,box_number)%moves,&
                     move_updates(molecule_type,box_number)%moves,&
                     box_number,start_atom,&
                     molecule_type,bias_energy(box_number),&
                     lreject,rng_stream,error=error)
                  IF(lreject) EXIT
               ENDIF

            ENDDO

! now do a Quickstep calculation to see if we accept the sequence
            CALL mc_Quickstep_move(mc_par,force_env,bias_env,&
               moves,lreject,move_updates,energy_check(:),r_old(:,:,:),&
               nnstep,old_energy(:),bias_energy(:),last_bias_energy(:),&
               nboxes,box_flag(:),oldsys,particles_old,&
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

               DO ibox=1,nboxes

! write the molecule information
                  WRITE(com_mol,*) nnstep,nchains(:,ibox)

! write the move statistics to file
                  DO itype=1,nmol_types
                     CALL write_move_stats(moves(itype,ibox)%moves,&
                          nnstep,move_unit(ibox))
                  ENDDO

! write a restart file
                  CALL write_mc_restart(nnstep,mc_par(ibox)%mc_par,&
                       nchains(:,ibox),force_env(ibox)%force_env)

! write cell lengths
                  WRITE(cell_unit,*) nnstep,abc(1:3,ibox)*angstrom

! write particle coordinates
                  WRITE(cbox,'(I4)') ibox
                  IF(SUM(nchains(:,ibox)) == 0) THEN
                     WRITE(com_crd,*) ' 0'
                     WRITE(com_crd,*) 'BOX ' // TRIM(ADJUSTL(cbox))
                  ELSE
                     CALL write_particle_coordinates(&
                          particles_old(ibox)%list%els,&
                          com_crd,dump_xmol,'POS',&
                          'BOX ' // TRIM(ADJUSTL(cbox)))
                  ENDIF
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
                  REAL(SUM(nchains(:,ibox)),dp)/REAL(nnstep-nstart,dp)
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
               DO itype=1,nmol_types
                  CALL mc_move_update(mc_par(ibox)%mc_par,&
                       move_updates(itype,ibox)%moves,itype,&
                       "trans",nnstep,ionode)
               ENDDO
            ENDIF

            IF (MOD(nnstep-nstart,iupvolume) == 0  ) THEN
               CALL mc_move_update(mc_par(ibox)%mc_par,&
                    move_updates(1,ibox)%moves,1337,&
                    "volume",nnstep,ionode)
            ENDIF
         ENDDO

! check to see if there are any overlaps in the boxes, and fold coordinates
         DO ibox=1,nboxes
            IF(SUM(nchains(:,ibox)) .NE. 0) THEN
               start_mol=1
               DO jbox=1,ibox-1
                  start_mol=start_mol+SUM(nchains(:,jbox))
               ENDDO
               end_mol=start_mol+SUM(nchains(:,ibox))-1
               CALL check_for_overlap(force_env(ibox)%force_env,&
                    nchains(:,ibox),nunits,loverlap,&
                    mol_type(start_mol:end_mol))
               IF(loverlap) THEN
                  IF(ionode) WRITE(iw,*) nnstep
                  CALL stop_program(routineN,moduleN,__LINE__,&
                  'coordinate overlap at the end of the above step')
! now fold the coordinates...don't do this anywhere but here, because
! we can get screwed up with the mc_molecule_info stuff (like in swap move)...
! this is kind of ugly, with allocated and deallocating every time
                  ALLOCATE (r_temp(1:3,1:nunits_tot(ibox)),STAT=istat)
                  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_temp",dp_size*nunits_tot(ibox)*3)

                  DO iunit=1,nunits_tot(ibox)
                     r_temp(1:3,iunit)=&
                          particles_old(ibox)%list%els(iunit)%r(1:3)
                  ENDDO
                  
                  CALL mc_coordinate_fold(r_temp(:,:),&
                       SUM(nchains(:,ibox)),mol_type(start_mol:end_mol),&
                       mass,nunits,abc(1:3,ibox))

! save the folded coordinates
                  DO iunit=1,nunits_tot(ibox)
                     r_old(1:3,iunit,ibox)=r_temp(1:3,iunit)
                     particles_old(ibox)%list%els(iunit)%r(1:3)=&
                          r_temp(1:3,iunit)
                  ENDDO

! if we're biasing, we need to do the same
                  IF(lbias) THEN
                     CALL force_env_get(bias_env(ibox)%force_env,&
                          subsys=biassys,error=error)
                     CALL cp_subsys_get(biassys(1)%subsys,&
                          particles=particles_bias, error=error)

                     DO iunit=1,nunits_tot(ibox)
                        particles_bias%els(iunit)%r(1:3)=&
                             r_temp(1:3,iunit)
                     ENDDO
                  ENDIF

                  DEALLOCATE (r_temp,STAT=istat)
                  IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r_temp")

               ENDIF
            ENDIF
         ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! debug make sure the energies line up
         IF(.FALSE.) THEN
         DO ibox=1,nboxes
            IF(SUM(nchains(:,ibox)) .NE. 0) THEN
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
               WRITE(iw,*) 'Box ',ibox
               WRITE(iw,*) 'nchains ',nchains(:,ibox)
               STOP 'mc_ensembles.F'
            ENDIF
         ENDDO
         ENDIF
      ENDDO

! write a restart file
      IF(ionode) THEN
         DO ibox=1,nboxes
            CALL write_mc_restart(nnstep,mc_par(ibox)%mc_par,&
               nchains(:,ibox),force_env(ibox)%force_env)
         ENDDO
      ENDIF

! calculate the final energy
      DO ibox=1,nboxes
         IF(SUM(nchains(:,ibox)) .NE. 0) THEN
            CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.,error=error)
            CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=final_energy(ibox),&
               error=error)
         ELSE
            final_energy(ibox)=0.0d0
         ENDIF
         IF (lbias) THEN
            CALL section_vals_release(bias_env(ibox)%force_env%root_section,error=error)
            CALL force_env_release(bias_env(ibox)%force_env,error=error)
         ENDIF
      ENDDO

! do some stuff in serial
      IF (ionode) THEN

         WRITE(com_ene,*) 'Final Energies:                      ',&
                  final_energy(1:nboxes)

         DO ibox=1,nboxes
            WRITE(cbox,'(I4)') ibox
            IF(SUM(nchains(:,ibox)) == 0) THEN
               WRITE(com_crd,*) ' 0'
               WRITE(com_crd,*) 'BOX ' // TRIM(ADJUSTL(cbox))
            ELSE
               CALL write_particle_coordinates(&
               particles_old(ibox)%list%els,&
               com_crd,dump_xmol,'POS',&
               'FINAL BOX ' // TRIM(ADJUSTL(cbox)))
            ENDIF

! write a bunch of data to the screen
            WRITE(iw,'(A)') &
               '------------------------------------------------'
            WRITE(iw,'(A,I1,A)') &
               '|                   BOX ',ibox,&
               '                      |'
            WRITE(iw,'(A)') &
               '------------------------------------------------'
            test_moves => moves(:,ibox)
            CALL final_mc_write(mc_par(ibox)%mc_par,test_moves,&
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
         DO itype=1,nmol_types
            CALL mc_moves_release(move_updates(itype,ibox)%moves)
            CALL mc_moves_release(moves(itype,ibox)%moves)
         ENDDO
         CALL mc_averages_release(averages(ibox)%averages)
      ENDDO

      DEALLOCATE(pmvol_box,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "pmvol_box")
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

END MODULE mc_ensembles

!**********************************************************************

