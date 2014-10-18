!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Used to run the bulk of the MC simulation, doing things like
!>      choosing move types and writing data to files
!> \author Matthew J. McGrath  (09.26.2003)
!>
!>    REVISIONS
!>      09.10.05  MJM combined the two subroutines in this module into one
! *****************************************************************************
MODULE mc_ensembles
  USE cell_types,                      ONLY: cell_p_type,&
                                             get_cell
  USE cp_external_control,             ONLY: external_control
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_release,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: pi
  USE mc_control,                      ONLY: mc_create_bias_force_env,&
                                             write_mc_restart
  USE mc_coordinates,                  ONLY: check_for_overlap,&
                                             create_discrete_array,&
                                             find_mc_test_molecule,&
                                             get_center_of_mass,&
                                             mc_coordinate_fold,&
                                             rotate_molecule
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
  USE mc_moves,                        ONLY: mc_avbmc_move,&
                                             mc_conformation_change,&
                                             mc_hmc_move,&
                                             mc_molecule_rotation,&
                                             mc_molecule_translation,&
                                             mc_volume_move
  USE mc_types,                        ONLY: get_mc_molecule_info,&
                                             get_mc_par,&
                                             mc_averages_p_type,&
                                             mc_input_file_type,&
                                             mc_molecule_info_type,&
                                             mc_moves_p_type,&
                                             mc_simulation_parameters_p_type,&
                                             set_mc_par
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_bcast
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_methods,                ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule,&
                                             n_avogadro
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_ensembles'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

  PUBLIC :: mc_run_ensemble,mc_compute_virial

CONTAINS

! *****************************************************************************
!> \brief directs the program in running one or two box MC simulations
!> \param mc_env a pointer that contains all mc_env for all the simulation
!>          boxes
!> \param para_env ...
!> \param globenv the global environment for the simulation
!> \param input_declaration ...
!> \param nboxes the number of simulation boxes
!> \param rng_stream the stream we pull random numbers from
!>
!>    Suitable for parallel.
!> \param error the error type
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_run_ensemble ( mc_env, para_env, globenv, input_declaration, nboxes, rng_stream, error)

    TYPE(mc_environment_p_type), &
      DIMENSION(:), POINTER                  :: mc_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_type), POINTER              :: input_declaration
    INTEGER, INTENT(IN)                      :: nboxes
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_run_ensemble', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: atom_names_box
    CHARACTER(default_string_length), &
      DIMENSION(:, :), POINTER               :: atom_names
    CHARACTER(LEN=20)                        :: ensemble
    CHARACTER(LEN=40)                        :: cbox, cstep, fft_lib, &
                                                move_type, move_type_avbmc
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    INTEGER, DIMENSION(:), POINTER           :: avbmc_atom, mol_type, &
                                                nchains_box, nunits, &
                                                nunits_tot
    INTEGER, DIMENSION(1:nboxes)             :: box_flag, cl, data_unit, &
                                                diff, istep, move_unit, rm
    INTEGER, DIMENSION(1:3, 1:2)             :: discrete_array
    INTEGER :: atom_number, box_number, cell_unit, com_crd, com_ene, com_mol, &
      end_mol, group, handle, ibox, idum, imol_type, imolecule, imove, &
      iparticle, iprint, itype, iunit, iuptrans, iupvolume, iw, jbox, jdum, &
      molecule_type, molecule_type_swap, molecule_type_target, nchain_total, &
      nmol_types, nmoves, nnstep, nstart, nstep, source, start_atom, &
      start_atom_swap, start_atom_target, start_mol, stat
    CHARACTER(LEN=default_string_length)     :: unit_str
    CHARACTER(LEN=40), DIMENSION(1:nboxes) :: cell_file, coords_file, &
      data_file, displacement_file, energy_file, molecules_file, moves_file
    LOGICAL                                  :: failure, ionode, lbias, &
                                                ldiscrete, lhmc, &
                                                lnew_bias_env, loverlap, &
                                                lreject, lstop, should_stop
    REAL(dp), DIMENSION(:), POINTER          :: pbias, pmavbmc_mol, &
                                                pmhmc_box, pmrot_mol, &
                                                pmtraion_mol, pmtrans_mol, &
                                                pmvol_box
    REAL(dp), DIMENSION(:, :), POINTER       :: conf_prob, mass
    REAL(KIND=dp) :: discrete_step, pmavbmc, pmhmc, pmswap, pmtraion, &
      pmtrans, pmvolume, rand, test_energy, unit_conv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_temp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r_old
    REAL(KIND=dp), DIMENSION(1:3, 1:nboxes)  :: abc
    REAL(KIND=dp), DIMENSION(1:nboxes)       :: bias_energy, energy_check, &
                                                final_energy, initial_energy, &
                                                last_bias_energy, old_energy
    TYPE(cell_p_type), DIMENSION(:), POINTER :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: oldsys
    TYPE(cp_subsys_type), POINTER            :: biassys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: bias_env, force_env
    TYPE(mc_averages_p_type), DIMENSION(:), &
      POINTER                                :: averages
    TYPE(mc_input_file_type), POINTER        :: mc_bias_file
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: test_moves
    TYPE(mc_moves_p_type), DIMENSION(:, :), &
      POINTER                                :: move_updates, moves
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old
    TYPE(particle_list_type), POINTER        :: particles_bias
    TYPE(section_vals_type), POINTER         :: root_section

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)

    ! nullify some pointers
    NULLIFY(moves,move_updates,test_moves,root_section)

    ! allocate a whole bunch of stuff based on how many boxes we have
    ALLOCATE (force_env(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (bias_env(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (cell(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (particles_old(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (oldsys(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (averages(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (mc_par(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (pmvol_box(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (pmhmc_box(1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ibox=1,nboxes
       CALL get_mc_env ( mc_env(ibox)%mc_env, &
            mc_par = mc_par(ibox)%mc_par,&
            force_env=force_env(ibox)%force_env)
    ENDDO

    ! Gather units of measure for output (if available)
    root_section => force_env(1)%force_env%root_section
    CALL section_vals_val_get(root_section,"MOTION%PRINT%TRAJECTORY%UNIT",&
         c_val=unit_str,error=error)
    unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

    ! get some data out of mc_par
    CALL get_mc_par(mc_par(1)%mc_par,&
         ionode=ionode,source=source,group=group,&
         data_file=data_file(1),moves_file=moves_file(1),&
         cell_file=cell_file(1),coords_file=coords_file(1),&
         energy_file=energy_file(1),displacement_file=displacement_file(1),&
         lstop=lstop,nstep=nstep,nstart=nstart,pmvolume=pmvolume,pmhmc=pmhmc,&
         molecules_file=molecules_file(1),pmswap=pmswap,nmoves=nmoves,&
         pmtraion=pmtraion,pmtrans=pmtrans,iuptrans=iuptrans,&
         iupvolume=iupvolume,ldiscrete=ldiscrete,pmtraion_mol=pmtraion_mol,&
         lbias=lbias,iprint=iprint,pmavbmc_mol=pmavbmc_mol,&
         discrete_step=discrete_step,fft_lib=fft_lib,avbmc_atom=avbmc_atom,&
         pmavbmc=pmavbmc,pbias=pbias,mc_molecule_info=mc_molecule_info,&
         pmrot_mol=pmrot_mol,pmtrans_mol=pmtrans_mol,pmvol_box=pmvol_box(1),&
         ensemble=ensemble,pmhmc_box=pmhmc_box(1),lhmc=lhmc)

    ! get some data from the molecule types
    CALL get_mc_molecule_info(mc_molecule_info,conf_prob=conf_prob,&
         nchains=nchains,nmol_types=nmol_types,nunits_tot=nunits_tot,&
         mol_type=mol_type,nchain_total=nchain_total,nunits=nunits,&
         atom_names=atom_names,mass=mass)

    ! allocate some stuff based on the number of molecule types we have
    ALLOCATE (moves(1:nmol_types,1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (move_updates(1:nmol_types,1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(nboxes .GT. 1) THEN
       DO ibox=2,nboxes
          CALL get_mc_par(mc_par(ibox)%mc_par,&
               data_file=data_file(ibox),&
               moves_file=moves_file(ibox),&
               cell_file=cell_file(ibox),coords_file=coords_file(ibox),&
               energy_file=energy_file(ibox),&
               displacement_file=displacement_file(ibox),&
               molecules_file=molecules_file(ibox),pmvol_box=pmvol_box(ibox),&
               pmhmc_box=pmhmc_box(ibox))
       ENDDO
    ENDIF

    ! this is a check we can't do in the input checking
    IF(pmvol_box(nboxes) .LT. 1.0E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMVOL_BOX needs to be 1.0')
    ENDIF
    IF(pmhmc_box(nboxes) .LT. 1.0E0_dp) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The last value of PMHMC_BOX needs to be 1.0')
    ENDIF

    ! allocate the particle positions array for broadcasting
    ALLOCATE (r_old(3,SUM(nunits_tot),1:nboxes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! figure out what the default write unit is
    iw = cp_logger_get_default_io_unit(logger)

    IF (iw>0) THEN
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
          old_energy(ibox)=0.0E0_dp
       ENDIF
       initial_energy(ibox)=old_energy(ibox)

! don't care about overlaps if we're only doing HMC

       IF(.NOT. lhmc)THEN
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
       ENDIF

       ! get the subsystems and the cell information
       CALL force_env_get(force_env(ibox)%force_env,&
            subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
            error=error)
       CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
       CALL cp_subsys_get(oldsys(ibox)%subsys, &
            particles=particles_old(ibox)%list, error=error)
       ! record the old coordinates, in case a move is rejected
       DO iparticle=1,nunits_tot(ibox)
          r_old(1:3,iparticle,ibox)=&
               particles_old(ibox)%list%els(iparticle)%r(1:3)
       ENDDO

       ! find the bias energy of the initial run
       IF(lbias) THEN
          ! determine the atom names of every particle
          ALLOCATE (atom_names_box(1:nunits_tot(ibox)),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

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
               para_env,abc(:,ibox),nchains_box,input_declaration,mc_bias_file,&
               ionode,error)
          IF(SUM(nchains(:,ibox)) .NE. 0) THEN
             CALL force_env_calc_energy_force(bias_env(ibox)%force_env,&
                  calc_force=.FALSE.,error=error)
             CALL force_env_get(bias_env(ibox)%force_env,&
                  potential_energy=last_bias_energy(ibox),error=error)

          ELSE
             last_bias_energy(ibox)=0.0E0_dp
          ENDIF
          bias_energy(ibox)=last_bias_energy(ibox)
          DEALLOCATE(atom_names_box,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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
                  com_crd,dump_xmol,'POS','INITIAL BOX ' // TRIM(ADJUSTL(cbox)),&
                  unit_conv=unit_conv,error=error)
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

       IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
          WRITE(iw,*)
          WRITE(iw,*) "------- On Monte Carlo Step ",nnstep
       ENDIF

       IF (ionode) rand=next_random_number(rng_stream,error=error)
       ! broadcast the random number, to make sure we're on the same move
       CALL mp_bcast(rand,source,group)

       IF (rand .LT. pmvolume) THEN

          IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
             WRITE(iw,*) "Attempting a volume move"
             WRITE(iw,*)
          ENDIF

          SELECT CASE(ensemble)
          CASE("TRADITIONAL")
             CALL mc_volume_move( mc_par(1)%mc_par,&
                  force_env(1)%force_env,&
                  moves(1,1)%moves,move_updates(1,1)%moves,&
                  old_energy(1),1,&
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
                  old_energy(box_number),box_number,&
                  energy_check(box_number),r_old(:,:,box_number),iw,&
                  discrete_array(:,:),&
                  rng_stream,error=error)
          END SELECT

! update all the pointers here, because otherwise we may pass wrong information when we're making a bias environment
          DO ibox=1,nboxes
             CALL force_env_get(force_env(ibox)%force_env,&
                  subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
                  error=error)
             CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
             CALL cp_subsys_get(oldsys(ibox)%subsys,&
                  particles=particles_old(ibox)%list, error=error)
          ENDDO

          ! we need a new biasing environment now, if we're into that sort of thing
          IF(lbias) THEN
             DO ibox=1,nboxes
                CALL section_vals_release(bias_env(ibox)%force_env%root_section,&
                     error=error)
                CALL force_env_release(bias_env(ibox)%force_env,error=error)
                ! determine the atom names of every particle
                ALLOCATE (atom_names_box(1:nunits_tot(ibox)),STAT=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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

! need to find out what the cell lengths are
                CALL force_env_get(force_env(ibox)%force_env,&
                     subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
                     error=error)
                CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))

                CALL get_mc_par(mc_par(ibox)%mc_par,&
                     mc_bias_file=mc_bias_file)
                nchains_box=> nchains(:,ibox)

                CALL mc_create_bias_force_env(bias_env(ibox)%force_env,&
                     r_old(:,:,ibox),atom_names_box(:),nunits_tot(ibox),&
                     para_env,abc(:,ibox),nchains_box,input_declaration,&
                     mc_bias_file,ionode,error)

                IF(SUM(nchains(:,ibox)) .NE. 0) THEN
                   CALL force_env_calc_energy_force(&
                        bias_env(ibox)%force_env,&
                        calc_force=.FALSE.,error=error)
                   CALL force_env_get(bias_env(ibox)%force_env,&
                        potential_energy=last_bias_energy(ibox),&
                        error=error)
                ELSE
                   last_bias_energy(ibox)=0.0E0_dp
                ENDIF
                bias_energy(ibox)=last_bias_energy(ibox)
                DEALLOCATE(atom_names_box,STAT=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDDO
          ENDIF

       ELSEIF ( rand .LT. pmswap) THEN

          ! try a swap move
          IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
             WRITE(iw,*) "Attempting a swap move"
             WRITE(iw,*)
          ENDIF

          CALL mc_ge_swap_move(mc_par,force_env,bias_env,moves,&
               energy_check(:),r_old(:,:,:),old_energy(:),input_declaration,&
               para_env,bias_energy(:),last_bias_energy(:),rng_stream,error)

          ! the number of molecules may have changed, which deallocated the whole
          ! mc_molecule_info structure
          CALL get_mc_par(mc_par(1)%mc_par,mc_molecule_info=mc_molecule_info)
          CALL get_mc_molecule_info(mc_molecule_info,conf_prob=conf_prob,&
               nchains=nchains,nmol_types=nmol_types,nunits_tot=nunits_tot,&
               mol_type=mol_type,nchain_total=nchain_total,nunits=nunits,&
               atom_names=atom_names,mass=mass)

       ELSEIF ( rand .LT. pmhmc) THEN
! try hybrid Monte Carlo
          IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
             WRITE(iw,*) "Attempting a hybrid Monte Carlo move"
             WRITE(iw,*)
          ENDIF

! pick a box at random
          IF (ionode) rand=next_random_number(rng_stream,error=error)
          CALL mp_bcast(rand,source,group)

          DO ibox=1,nboxes
             IF(rand .LE. pmhmc_box(ibox)) THEN
                box_number=ibox
                EXIT
             ENDIF
          ENDDO

          CALL mc_hmc_move( mc_par(box_number)%mc_par,&
               force_env(box_number)%force_env,globenv,&
               moves(1,box_number)%moves,&
               move_updates(1,box_number)%moves,&
               old_energy(box_number),box_number,&
               energy_check(box_number),r_old(:,:,box_number),iw,&
               rng_stream,error=error)

       ELSEIF ( rand .LT. pmavbmc) THEN
          ! try an AVBMC move
          IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
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
                     start_atom_swap,idum,jdum,rng_stream,&
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

          IF(MOD(nnstep,iprint) == 0 .AND. (iw>0)) THEN
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
                           jdum,rng_stream,&
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
               rng_stream,unit_conv,error=error)

       ENDIF

       ! make sure the pointers are pointing correctly since the subsys may
       ! have changed
       DO ibox=1,nboxes
          CALL force_env_get(force_env(ibox)%force_env,&
               subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,&
               error=error)
          CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
          CALL cp_subsys_get(oldsys(ibox)%subsys,&
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
                WRITE(cstep,'(I8)') nnstep
                IF(SUM(nchains(:,ibox)) == 0) THEN
                   WRITE(com_crd,*) ' 0'
                   WRITE(com_crd,*) 'BOX ' // TRIM(ADJUSTL(cbox)) // &
                        ',  STEP ' // TRIM(ADJUSTL(cstep))
                ELSE
                   CALL write_particle_coordinates(&
                        particles_old(ibox)%list%els,&
                        com_crd,dump_xmol,'POS',&
                        'BOX ' // TRIM(ADJUSTL(cbox))// &
                        ',  STEP ' // TRIM(ADJUSTL(cstep)),&
                        unit_conv=unit_conv,error=error)
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
       CALL external_control(should_stop,"MC",globenv=globenv,error=error)
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
! don't care about overlaps if we're only doing HMC
       IF(.NOT. lhmc)THEN
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
                   IF(iw>0) WRITE(iw,*) nnstep
                   CALL stop_program(routineN,moduleN,__LINE__,&
                        'coordinate overlap at the end of the above step')
                ! now fold the coordinates...don't do this anywhere but here, because
                ! we can get screwed up with the mc_molecule_info stuff (like in swap move)...
                ! this is kind of ugly, with allocated and deallocating every time
                   ALLOCATE (r_temp(1:3,1:nunits_tot(ibox)),STAT=stat)
                   CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

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
                      CALL cp_subsys_get(biassys,&
                           particles=particles_bias, error=error)

                      DO iunit=1,nunits_tot(ibox)
                         particles_bias%els(iunit)%r(1:3)=&
                              r_temp(1:3,iunit)
                      ENDDO
                   ENDIF

                   DEALLOCATE (r_temp,STAT=stat)
                   CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                ENDIF
             ENDIF
          ENDDO
       ENDIF

       !debug code
       IF(debug_this_module) THEN
          DO ibox=1,nboxes
             IF(SUM(nchains(:,ibox)) .NE. 0) THEN
                CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
                     calc_force=.FALSE.,error=error)
                CALL force_env_get(force_env(ibox)%force_env,&
                     potential_energy=test_energy,&
                     error=error)
             ELSE
                test_energy=0.0E0_dp
             ENDIF

             IF(ABS(initial_energy(ibox)+energy_check(ibox)-&
                  test_energy) .GT. 0.0000001E0_dp) THEN
                IF (iw>0) THEN
                   WRITE(iw,*) '!!!!!!! We have an energy problem. !!!!!!!!'
                   WRITE(iw,'(A,T64,F16.10)') 'Final Energy = ',test_energy
                   WRITE(iw,'(A,T64,F16.10)') 'Inital Energy+energy_check=',&
                        initial_energy(ibox)+energy_check(ibox)
                   WRITE(iw,*) 'Box ',ibox
                   WRITE(iw,*) 'nchains ',nchains(:,ibox)
                END IF
                CALL mp_abort('!!!!!!! We have an energy problem. !!!!!!!!')
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
          final_energy(ibox)=0.0E0_dp
       ENDIF
       IF (lbias) THEN
          CALL section_vals_release(bias_env(ibox)%force_env%root_section,error=error)
          CALL force_env_release(bias_env(ibox)%force_env,error=error)
       ENDIF
    ENDDO

    ! do some stuff in serial
    IF (ionode.OR.(iw>0)) THEN

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
                  'FINAL BOX ' // TRIM(ADJUSTL(cbox)),unit_conv=unit_conv,&
                  error=error)
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

    DEALLOCATE(pmhmc_box,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pmvol_box,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r_old,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(force_env,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(bias_env,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (cell,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (particles_old,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (oldsys,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (averages,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (moves,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (move_updates,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (mc_par,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE mc_run_ensemble

! *****************************************************************************
!> \brief Computes the second virial coefficient of a molecule by using the integral form
!>      of the second virial coefficient found in McQuarrie "Statistical Thermodynamics",
!>      B2(T) = -2Pi Int 0toInf [ Exp[-beta*u(r)] -1] r^2 dr     Eq. 15-25
!>      I use trapazoidal integration with various step sizes
!>      (the integral is broken up into three parts, currently, but that's easily
!>      changed by the first variables found below).  It generates nvirial configurations,
!>      doing the integration for each one, and then averages all the B2(T) to produce
!>      the final answer.
!> \param mc_env a pointer that contains all mc_env for all the simulation
!>          boxes
!> \param rng_stream the stream we pull random numbers from
!>
!>    Suitable for parallel.
!> \param error the error type
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_compute_virial ( mc_env, rng_stream, error)

    TYPE(mc_environment_p_type), &
      DIMENSION(:), POINTER                  :: mc_env
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_compute_virial', &
      routineP = moduleN//':'//routineN

    INTEGER :: current_division, end_atom, group, ibin, idivision, iparticle, &
      iprint, itemp, iunit, ivirial, iw, nbins, nchain_total, &
      nintegral_divisions, nmol_types, nvirial, nvirial_temps, source, &
      start_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: failure, ionode, loverlap
    REAL(dp), DIMENSION(:), POINTER          :: BETA, virial_cutoffs, &
                                                virial_stepsize, virial_temps
    REAL(dp), DIMENSION(:, :), POINTER       :: mass, mayer, r_old
    REAL(KIND=dp) :: ave_virial, current_value, distance, exp_max_val, &
      exp_min_val, exponent, integral, previous_value, square_value, &
      trial_energy, triangle_value
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, center_of_mass
    TYPE(cell_p_type), DIMENSION(:), POINTER :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    ! these are current magic numbers for how we compute the virial...
    ! we break it up into three parts to integrate the function so provide
    ! better statistics

    nintegral_divisions=3
    ALLOCATE (virial_cutoffs(1:nintegral_divisions),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (virial_stepsize(1:nintegral_divisions),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    virial_cutoffs(1)=8.0 ! first distance, in bohr
    virial_cutoffs(2)=13.0 ! second distance, in bohr
    virial_cutoffs(3)=22.0 ! maximum distance, in bohr
    virial_stepsize(1)=0.04 ! stepsize from 0 to virial_cutoffs(1)
    virial_stepsize(2)=0.1
    virial_stepsize(3)=0.2

    nbins=CEILING(virial_cutoffs(1)/virial_stepsize(1)+(virial_cutoffs(2)-&
         virial_cutoffs(1))/virial_stepsize(2)+(virial_cutoffs(3)-&
         virial_cutoffs(2))/virial_stepsize(3))

    ! figure out what the default write unit is
    iw = cp_logger_get_default_io_unit(logger)

    ! allocate a whole bunch of stuff based on how many boxes we have
    ALLOCATE (force_env(1:1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (cell(1:1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (particles(1:1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (subsys(1:1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (mc_par(1:1),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_mc_env ( mc_env(1)%mc_env, &
         mc_par = mc_par(1)%mc_par,&
         force_env=force_env(1)%force_env)

    ! get some data out of mc_par
    CALL get_mc_par(mc_par(1)%mc_par,&
         exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,nvirial=nvirial,&
         ionode=ionode,source=source,group=group,&
         mc_molecule_info=mc_molecule_info,virial_temps=virial_temps)

    IF (iw>0) THEN
       WRITE(iw,*)
       WRITE(iw,*)
       WRITE(iw,*) 'Beginning the calculation of the second virial coefficient'
       WRITE(iw,*)
       WRITE(iw,*)
    ENDIF

    ! get some data from the molecule types
    CALL get_mc_molecule_info(mc_molecule_info,&
         nchains=nchains,nmol_types=nmol_types,nunits_tot=nunits_tot,&
         mol_type=mol_type,nchain_total=nchain_total,nunits=nunits,&
         mass=mass)

    nvirial_temps=SIZE(virial_temps)
    ALLOCATE (BETA(1:nvirial_temps),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    DO itemp=1,nvirial_temps
       BETA(itemp) = 1 / virial_temps(itemp) / boltzmann * joule
    ENDDO

    ! get the subsystems and the cell information
    CALL force_env_get(force_env(1)%force_env,&
         subsys=subsys(1)%subsys,cell=cell(1)%cell,&
         error=error)
    CALL get_cell(cell(1)%cell,abc=abc(:))
    CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles(1)%list, error=error)

    ! check and make sure the box is big enough
    IF(abc(1) .NE. abc(2) .OR. abc(2) .NE. abc(3)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            'The box needs to be cubic for a virial calculation (it is easiest).')
    END IF
    IF(virial_cutoffs(nintegral_divisions) .GT. abc(1)/2.0E0_dp) THEN
       IF (iw>0) &
            WRITE(iw,*) "Box length ",abc(1)*angstrom," virial cutoff ",&
            virial_cutoffs(nintegral_divisions)*angstrom
       CALL stop_program(routineN,moduleN,__LINE__,&
            'You need a bigger box to deal with this virial cutoff (see above).')
    END IF

    ! store the coordinates of the molecules in an array so we can work with it
    ALLOCATE (r_old(1:3,1:nunits_tot(1)),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iparticle=1,nunits_tot(1)
       r_old(1:3,iparticle)=&
            particles(1)%list%els(iparticle)%r(1:3)
    ENDDO

    ! move the center of mass of molecule 1 to the origin
    start_atom=1
    end_atom=nunits(mol_type(1))
    CALL get_center_of_mass(r_old(:,start_atom:end_atom),nunits(mol_type(1)),&
         center_of_mass(:),mass(1:nunits(mol_type(1)),mol_type(1)))
    DO iunit=start_atom,end_atom
       r_old(:,iunit)=r_old(:,iunit)-center_of_mass(:)
    ENDDO
    ! set them in the force_env, so the first molecule is ready for the energy calc
    DO iparticle=start_atom,end_atom
       particles(1)%list%els(iparticle)%r(1:3)=r_old(1:3,iparticle)
    ENDDO

    ! print out a notice every 1%
    iprint=FLOOR(REAL(nvirial)/100.0)
    IF(iprint == 0) iprint=1

    ! we'll compute the average potential, and then integrate that, as opposed to
    ! integrating every orientation and then averaging
    ALLOCATE (mayer(1:nvirial_temps,1:nbins),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    mayer(:,:)=0.0_dp

    ! loop over all nvirial random configurations
    DO ivirial=1,nvirial

       ! move molecule two back to the origin
       start_atom=nunits(mol_type(1))+1
       end_atom=nunits_tot(1)
       CALL get_center_of_mass(r_old(:,start_atom:end_atom),nunits(mol_type(2)),&
            center_of_mass(:),mass(1:nunits(mol_type(2)),mol_type(2)))
       DO iunit=start_atom,end_atom
          r_old(:,iunit)=r_old(:,iunit)-center_of_mass(:)
       ENDDO

       ! now we need a random orientation for molecule 2...this routine is
       ! only done in serial since it calls a random number
       IF(ionode) THEN
          CALL rotate_molecule(r_old(:,start_atom:end_atom),&
               mass(1:nunits(mol_type(2)),mol_type(2)),&
               nunits(mol_type(2)),rng_stream,error=error)
       ENDIF
       CALL mp_bcast(r_old(:,:),source,group)

       distance=0.0E0_dp
       ibin=1
       DO
          ! find out what our stepsize is
          current_division=0
          DO idivision=1,nintegral_divisions
             IF(distance .LT. virial_cutoffs(idivision)-virial_stepsize(idivision)/2.0E0_dp) THEN
                current_division=idivision
                EXIT
             ENDIF
          ENDDO
          IF(current_division == 0) EXIT
          distance=distance+virial_stepsize(current_division)

          ! move the second molecule only along the x direction
          DO iparticle=start_atom,end_atom
             particles(1)%list%els(iparticle)%r(1)=r_old(1,iparticle)+distance
             particles(1)%list%els(iparticle)%r(2)=r_old(2,iparticle)
             particles(1)%list%els(iparticle)%r(3)=r_old(3,iparticle)
          ENDDO

          ! check for overlaps
          CALL check_for_overlap (force_env(1)%force_env,nchains(:,1),nunits,loverlap,mol_type)

          ! compute the energy if there is no overlap
          ! exponent is exp(-beta*energy)-1, also called the Mayer term
          IF(loverlap) THEN
             DO itemp=1,nvirial_temps
                mayer(itemp,ibin)=mayer(itemp,ibin)-1.0_dp
             ENDDO
          ELSE
             CALL force_env_calc_energy_force(force_env(1)%force_env,&
                  calc_force=.FALSE.,error=error)
             CALL force_env_get(force_env(1)%force_env,&
                  potential_energy=trial_energy,error=error)

             DO itemp=1,nvirial_temps

                exponent=-BETA(itemp)*trial_energy

                IF    (exponent .GT. exp_max_val) THEN
                   exponent=exp_max_val
                ELSEIF(exponent .LT. exp_min_val) THEN
                   exponent=exp_min_val
                ENDIF
                mayer(itemp,ibin)=mayer(itemp,ibin)+EXP(exponent)-1.0_dp
             ENDDO
          ENDIF

          ibin=ibin+1
       ENDDO
       ! write out some info that keeps track of where we are
       IF(iw>0) THEN
          IF(MOD(ivirial,iprint) == 0) &
               WRITE(iw,'(A,I6,A,I6)') ' Done with config ',ivirial,' out of ',nvirial
       ENDIF
    ENDDO

    ! now we integrate this average potential
    mayer(:,:)=mayer(:,:)/REAL(nvirial,dp)

    DO itemp=1,nvirial_temps
       integral=0.0_dp
       previous_value=0.0_dp
       distance=0.0E0_dp
       ibin=1
       DO
          current_division=0
          DO idivision=1,nintegral_divisions
             IF(distance .LT. virial_cutoffs(idivision)-virial_stepsize(idivision)/2.0E0_dp) THEN
                current_division=idivision
                EXIT
             ENDIF
          ENDDO
          IF(current_division == 0) EXIT
          distance=distance+virial_stepsize(current_division)

          ! now we need to integrate, using the trapazoidal method
          ! first, find the value of the square
          current_value=mayer(itemp,ibin)*distance**2
          square_value=previous_value*virial_stepsize(current_division)
          ! now the triangle that sits on top of it, which is half the size of this square...
          ! notice this is negative if the current value is less than the previous value
          triangle_value=0.5E0_dp*((current_value-previous_value)*virial_stepsize(current_division))

          integral=integral+square_value+triangle_value
          previous_value=current_value
          ibin=ibin+1
       ENDDO

       ! now that the integration is done, compute the second virial that results
       ave_virial=-2.0E0_dp*pi*integral

       ! convert from CP2K units to something else
       ave_virial=ave_virial*n_avogadro*angstrom**3/1.0E8_dp**3

       IF(iw>0) THEN
          WRITE(iw,*)
          WRITE(iw,*) '*********************************************************************'
          WRITE(iw,'(A,F12.6,A)') ' ***                Temperature = ',virial_temps(itemp),&
               '                     ***'
          WRITE(iw,*) '***                                                               ***'
          WRITE(iw,'(A,E12.6,A)') ' ***                  B2(T) = ',ave_virial,&
               ' cm**3/mol               ***'
          WRITE(iw,*) '*********************************************************************'
          WRITE(iw,*)
       ENDIF
    ENDDO

    ! deallocate some stuff
    DEALLOCATE(mc_par,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(subsys,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(force_env,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(particles,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cell,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(virial_cutoffs,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(virial_stepsize,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r_old,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mayer,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(BETA,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE mc_compute_virial

END MODULE mc_ensembles

