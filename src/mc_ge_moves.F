!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the Monte Carlo moves that can handle more than one
!>      box, including the Quickstep move, a volume swap between boxes,
!>       and a particle swap between boxes
!> \par History
!>      MJM (07.28.2005): make the Quickstep move general, and changed
!>                        the swap and volume moves to work with the
!>                        CP2K classical routines
!> \author Matthew J. McGrath  (01.25.2004)
! *****************************************************************************
MODULE mc_ge_moves
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_p_type,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release,&
                                             force_env_set_cell
  USE input_constants,                 ONLY: dump_xmol
  USE input_section_types,             ONLY: section_vals_release
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size
  USE mc_control,                      ONLY: mc_create_force_env
  USE mc_coordinates,                  ONLY: check_for_overlap,&
                                             generate_cbmc_swap_config,&
                                             get_center_of_mass
  USE mc_misc,                         ONLY: mc_make_dat_file_new
  USE mc_move_control,                 ONLY: move_q_reinit,&
                                             q_move_accept
  USE mc_types,                        ONLY: &
       get_mc_molecule_info, get_mc_par, mc_determine_molecule_info, &
       mc_input_file_type, mc_molecule_info_destroy, mc_molecule_info_type, &
       mc_moves_p_type, mc_simulation_parameters_p_type, set_mc_par
  USE message_passing,                 ONLY: mp_bcast
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_ge_moves'

  PUBLIC :: mc_ge_volume_move,mc_ge_swap_move,&
      mc_quickstep_move

CONTAINS

! *****************************************************************************
!> \brief computes the acceptance of a series of biased or unbiased moves
!>      (translation, rotation, conformational changes)
!> \param mc_par the mc parameters for the force envs of the boxes
!> \param force_env the force environments for the boxes
!> \param bias_env the force environments with the biasing potential for the boxes
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected for both boxes
!> \param lreject automatically rejects the move (used when an overlap occurs in
!>                the sequence of moves)
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated for both boxes
!> \param energy_check the running total of how much the energy has changed
!>                      since the initial configuration
!> \param r_old the coordinates of the last accepted move before the sequence
!>         whose acceptance is determined by this call
!> \param nnstep the Monte Carlo step we're on
!> \param old_energy the energy of the last accepted move involving the full potential
!> \param bias_energy_new the energy of the current configuration involving the bias potential
!> \param nboxes the number of boxes (force environments) in the system
!> \param box_flag indicates if a move has been tried in a given box..if not, we don't
!>        recompute the energy
!> \param subsys the pointers for the particle subsystems of both boxes
!> \param particles the pointers for the particle sets
!> \param rng_stream the stream we pull random numbers from
!> \param error the cp_error_type in case something goes wrong
!>
!>      Designed for parallel use.
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_Quickstep_move(mc_par,force_env,bias_env,moves,&
                        lreject,move_updates,energy_check,r_old,&
                        nnstep,old_energy,bias_energy_new,last_bias_energy,&
                        nboxes,box_flag,subsys,particles,rng_stream,&
                        unit_conv,error)

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env, bias_env
    TYPE(mc_moves_p_type), DIMENSION(:, :), &
      POINTER                                :: moves
    LOGICAL, INTENT(IN)                      :: lreject
    TYPE(mc_moves_p_type), DIMENSION(:, :), &
      POINTER                                :: move_updates
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: energy_check
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: r_old
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: old_energy, bias_energy_new, &
                                                last_bias_energy
    INTEGER, INTENT(IN)                      :: nboxes
    INTEGER, DIMENSION(:), INTENT(IN)        :: box_flag
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsys
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp), INTENT(IN)                :: unit_conv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_Quickstep_move', &
      routineP = moduleN//':'//routineN

    INTEGER :: end_mol, group, handle, ibox, iparticle, iprint, istat, itype, &
      jbox, nmol_types, source, start_mol
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(1:nboxes)             :: diff
    LOGICAL                                  :: ionode, lbias, loverlap
    REAL(KIND=dp)                            :: BETA, energies, rand, w
    REAL(KIND=dp), DIMENSION(1:nboxes)       :: bias_energy_old, new_energy
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsys_bias
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_bias

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

      NULLIFY(subsys_bias,particles_bias)

! get a bunch of data from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,lbias=lbias,&
         BETA=BETA,diff=diff(1),source=source,group=group,&
         iprint=iprint,&
         mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types,&
           nchains=nchains,nunits_tot=nunits_tot,nunits=nunits,mol_type=mol_type)

      IF(nboxes .GT. 1) THEN
         DO ibox=2,nboxes
            CALL get_mc_par(mc_par(ibox)%mc_par,diff=diff(ibox))
         ENDDO
      ENDIF

! allocate some stuff
      ALLOCATE(subsys_bias(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "subsys_bias",nboxes)
      ALLOCATE(particles_bias(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_bias",nboxes)

! record the attempt...we really only care about molecule type 1 and box
! type 1, since the acceptance will be identical for all boxes and molecules
      moves(1,1)%moves%Quickstep%attempts=&
           moves(1,1)%moves%Quickstep%attempts+1

! grab the coordinates for the force_env
      DO ibox=1,nboxes
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=subsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(subsys(ibox)%subsys, &
            particles=particles(ibox)%list, error=error)
      ENDDO

! calculate the new energy of the system...if we're biasing,
! force_env hasn't changed but bias_env has
      DO ibox=1,nboxes
         IF(box_flag(ibox) == 1) THEN
            IF(lbias) THEN
! grab the coords from bias_env and put them into force_env
               CALL force_env_get(bias_env(ibox)%force_env,&
                  subsys=subsys_bias(ibox)%subsys,error=error)
               CALL cp_subsys_get(subsys_bias(ibox)%subsys, &
                  particles=particles_bias(ibox)%list, error=error)

               DO iparticle=1,nunits_tot(ibox)
                  particles(ibox)%list%els(iparticle)%r(1:3)=&
                     particles_bias(ibox)%list%els(iparticle)%r(1:3)
               ENDDO

               CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
                  calc_force=.FALSE.,error=error)
               CALL force_env_get(force_env(ibox)%force_env,&
                  potential_energy=new_energy(ibox),error=error)
            ELSE
               IF( .NOT. lreject) THEN
                  CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
                     calc_force=.FALSE.,error=error)
                  CALL force_env_get(force_env(ibox)%force_env,&
                     potential_energy=new_energy(ibox),error=error)
               ENDIF
            ENDIF
         ELSE
            new_energy(ibox)=old_energy(ibox)
         ENDIF

      ENDDO

! accept or reject the move based on Metropolis or the Iftimie rule
      IF (ionode) THEN

! write them out in case something bad happens
         IF(MOD(nnstep,iprint) == 0) THEN
            DO ibox=1,nboxes
               IF(SUM(nchains(:,ibox)) == 0) THEN
                  WRITE(diff(ibox),*) nnstep
                  WRITE(diff(ibox),*) nchains(:,ibox)
               ELSE
                  WRITE(diff(ibox),*) nnstep
                  CALL write_particle_coordinates(&
                     particles(ibox)%list%els,&
                     diff(ibox),dump_xmol,'POS','TRIAL',&
                     unit_conv=unit_conv,error=error)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF(.NOT. lreject) THEN
         IF (lbias) THEN

            DO ibox=1,nboxes
! look for overlap
               IF(SUM(nchains(:,ibox)) .NE. 0) THEN
! find the molecule bounds
                  start_mol=1
                  DO jbox=1,ibox-1
                     start_mol=start_mol+SUM(nchains(:,jbox))
                  ENDDO
                  end_mol=start_mol+SUM(nchains(:,ibox))-1
                  CALL check_for_overlap(bias_env(ibox)%force_env,&
                     nchains(:,ibox),nunits(:),loverlap,mol_type(start_mol:end_mol))
                  IF(loverlap) CALL stop_program(routineN,&
                     moduleN,__LINE__,&
                     'Quickstep move found an overlap in the old config')
               ENDIF
               bias_energy_old(ibox)=last_bias_energy(ibox)
            ENDDO

            energies=-BETA*((SUM(new_energy(:))-SUM(bias_energy_new(:)))&
               -(SUM(old_energy(:))-SUM(bias_energy_old(:))))

! used to prevent over and underflows
            IF(energies .GE. -1.0E-8) THEN
               w=1.0_dp
            ELSEIF(energies .LE. -500.0_dp) THEN
               w=0.0_dp
            ELSE
               w=EXP(energies)
            ENDIF

            IF(ionode) THEN
               DO ibox=1,nboxes
                  WRITE(diff(ibox),*) nnstep,new_energy(ibox)-&
                     old_energy(ibox),&
                     bias_energy_new(ibox)-bias_energy_old(ibox)
               ENDDO
            ENDIF
         ELSE
            energies=-BETA*(SUM(new_energy(:))-SUM(old_energy(:)))
! used to prevent over and underflows
            IF(energies .GE. 0.0_dp) THEN
               w=1.0_dp
            ELSEIF(energies .LE. -500.0_dp) THEN
               w=0.0_dp
            ELSE
               w=EXP(energies)
            ENDIF
         ENDIF
      ELSE
         w=0.0E0_dp
      ENDIF
      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

      IF (rand .LT. w) THEN

! accept the move
         moves(1,1)%moves%Quickstep%successes=&
              moves(1,1)%moves%Quickstep%successes+1

         DO ibox=1,nboxes
! remember what kind of move we did for lbias=.false.
            IF(.NOT. lbias) THEN
               DO itype=1,nmol_types
                  CALL q_move_accept(moves(itype,ibox)%moves,.TRUE.)
                  CALL q_move_accept(move_updates(itype,ibox)%moves,.TRUE.)

! reset the counters
                  CALL move_q_reinit(moves(itype,ibox)%moves, .TRUE.)
                  CALL move_q_reinit(move_updates(itype,ibox)%moves, .TRUE.)
               ENDDO
            ENDIF

            DO itype=1,nmol_types
! we need to record all accepted moves since last Quickstep calculation
               CALL q_move_accept(moves(itype,ibox)%moves,.FALSE.)
               CALL q_move_accept(move_updates(itype,ibox)%moves,.FALSE.)

! reset the counters
               CALL move_q_reinit(moves(itype,ibox)%moves, .FALSE.)
               CALL move_q_reinit(move_updates(itype,ibox)%moves, .FALSE.)
            ENDDO

! update energies
            energy_check(ibox)=energy_check(ibox)+&
               (new_energy(ibox)-old_energy(ibox))
            old_energy(ibox)=new_energy(ibox)

         ENDDO

         IF ( lbias) THEN
            DO ibox=1,nboxes
               last_bias_energy(ibox)=bias_energy_new(ibox)
            ENDDO
         ENDIF

! update coordinates
         DO ibox=1,nboxes
            IF(nunits_tot(ibox) .NE. 0) THEN
               DO iparticle=1,nunits_tot(ibox)
                  r_old(1:3,iparticle,ibox)=&
                     particles(ibox)%list%els(iparticle)%r(1:3)
               ENDDO
            ENDIF
         ENDDO

      ELSE

 ! reject the move
         DO ibox=1,nboxes
            DO itype=1,nmol_types
               CALL move_q_reinit(moves(itype,ibox)%moves,.FALSE.)
               CALL move_q_reinit(move_updates(itype,ibox)%moves,.FALSE.)
               IF(.NOT. lbias) THEN
! reset the counters
                  CALL move_q_reinit(moves(itype,ibox)%moves, .TRUE.)
                  CALL move_q_reinit(move_updates(itype,ibox)%moves, .TRUE.)
               ENDIF
            ENDDO

         ENDDO

         IF ( .NOT. ionode) r_old(:,:,:) = 0.0E0_dp

! coodinates changed, so we need to broadcast those, even for the lbias
! case since bias_env needs to have the same coords as force_env
         CALL mp_bcast(r_old,source,group)

         DO ibox=1,nboxes
            DO iparticle=1,nunits_tot(ibox)
               particles(ibox)%list%els(iparticle)%r(1:3)=&
                  r_old(1:3,iparticle,ibox)
               IF(lbias .AND. box_flag(ibox) == 1) &
                  particles_bias(ibox)%list%els(iparticle)%r(1:3)=&
                  r_old(1:3,iparticle,ibox)
            ENDDO
         ENDDO

! need to reset the energies of the biasing potential
         IF ( lbias) THEN
            DO ibox=1,nboxes
               bias_energy_new(ibox)=last_bias_energy(ibox)
            ENDDO
         ENDIF

      ENDIF

! make sure the coordinates are transferred
      DO ibox=1,nboxes
         CALL cp_subsys_set(subsys(ibox)%subsys,&
            particles=particles(ibox)%list,error=error)
         IF(lbias .AND. box_flag(ibox) == 1) &
            CALL cp_subsys_set(subsys_bias(ibox)%subsys,&
            particles=particles_bias(ibox)%list,error=error)
      ENDDO

 ! deallocate some stuff
      DEALLOCATE(subsys_bias,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"subsys_bias")
      DEALLOCATE(particles_bias,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"particles_bias")

! end the timing
      CALL timestop(handle)

  END SUBROUTINE mc_Quickstep_move

! *****************************************************************************
!> \brief attempts a swap move between two simulation boxes
!> \param mc_par the mc parameters for the force envs of the boxes
!> \param force_env the force environments for the boxes
!> \param bias_env the force environments used to bias moves for the boxes
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected for both boxes
!> \param energy_check the running total of how much the energy has changed
!>                      since the initial configuration
!> \param r_old the coordinates of the last accepted move involving a
!>               full potential calculation for both boxes
!> \param old_energy the energy of the last accepted move involving a
!>                    a full potential calculation
!> \param para_env the parallel environment for this simulation
!> \param bias_energy_old the energies of both boxes computed using the biasing
!>        potential
!> \param last_bias_energy the energy for the biased simulations
!> \param rng_stream the stream we pull random numbers from
!> \param error the cp_error_type in case something goes wrong
!>
!>      Suitable for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_ge_swap_move(mc_par,force_env,bias_env,moves,&
                energy_check,r_old,old_energy,&
                para_env,bias_energy_old,last_bias_energy,&
                rng_stream,error)

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env, bias_env
    TYPE(mc_moves_p_type), DIMENSION(:, :), &
      POINTER                                :: moves
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: energy_check
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: r_old
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: old_energy
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: bias_energy_old, &
                                                last_bias_energy
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_ge_swap_move', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: atom_names_insert, &
                                                atom_names_remove
    CHARACTER(default_string_length), &
      DIMENSION(:, :), POINTER               :: atom_names
    CHARACTER(LEN=200)                       :: fft_lib
    CHARACTER(LEN=40), DIMENSION(1:2)        :: dat_file
    INTEGER :: end_mol, group, handle, iatom, ibox, idim, iiatom, imolecule, &
      ins_atoms, insert_box, ipart, istat, itype, jbox, molecule_type, &
      nmol_types, nswapmoves, print_level, rem_atoms, remove_box, source, &
      start_atom_ins, start_atom_rem, start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, mol_type_test, &
                                                nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains, nchains_test
    LOGICAL                                  :: ionode, lbias, loverlap, &
                                                loverlap_ins, loverlap_rem
    REAL(dp), DIMENSION(:), POINTER          :: eta_insert, eta_remove, &
                                                pmswap_mol
    REAL(dp), DIMENSION(:, :), POINTER       :: insert_coords, remove_coords
    REAL(KIND=dp) :: BETA, del_quickstep_energy, exp_max_val, exp_min_val, &
      max_val, min_val, prefactor, rand, rdum, w, weight_new, weight_old
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cbmc_energies, r_cbmc, &
                                                r_insert_mol
    REAL(KIND=dp), DIMENSION(1:2)            :: bias_energy_new, new_energy
    REAL(KIND=dp), DIMENSION(1:3)            :: abc_insert, abc_remove, &
                                                center_of_mass, &
                                                displace_molecule, pos_insert
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mass
    TYPE(cell_type), POINTER                 :: cell_insert, cell_remove
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: oldsys
    TYPE(cp_subsys_type), POINTER            :: insert_sys, remove_sys
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: test_env, test_env_bias
    TYPE(mc_input_file_type), POINTER        :: mc_bias_file, mc_input_file
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info, &
                                                mc_molecule_info_test
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old
    TYPE(particle_list_type), POINTER        :: particles_insert, &
                                                particles_remove

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! reset the overlap flag
      loverlap=.FALSE.

! nullify some pointers
      NULLIFY(particles_old,mol_type,mol_type_test,mc_input_file,mc_bias_file)
      NULLIFY(oldsys,atom_names,pmswap_mol,insert_coords,remove_coords)
      NULLIFY(eta_insert,eta_remove)

! grab some stuff from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,BETA=BETA,&
         max_val=max_val, min_val=min_val, exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,nswapmoves=nswapmoves,group=group,source=source,&
         lbias=lbias,dat_file=dat_file(1),fft_lib=fft_lib,&
         mc_molecule_info=mc_molecule_info,pmswap_mol=pmswap_mol)
      CALL get_mc_molecule_info(mc_molecule_info,nchains=nchains,&
           nunits=nunits,nunits_tot=nunits_tot,nmol_types=nmol_types,&
           atom_names=atom_names,mass=mass,mol_type=mol_type)

      print_level = 1

      CALL get_mc_par(mc_par(2)%mc_par,dat_file=dat_file(2))

! allocate some stuff
      ALLOCATE(oldsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "oldsys",2)
      ALLOCATE(particles_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_old",2)

! get the old coordinates
      DO ibox=1,2
         CALL force_env_get(force_env(ibox)%force_env,&
             subsys=oldsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(oldsys(ibox)%subsys, &
            particles=particles_old(ibox)%list, &
            error=error)
      ENDDO

!     choose a direction to swap
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

      IF ( rand .LE. 0.50E0_dp ) THEN
         remove_box=1
         insert_box=2
      ELSE
         remove_box=2
         insert_box=1
      ENDIF

! now assign the eta values for the insert and remove boxes
      CALL get_mc_par(mc_par(remove_box)%mc_par,eta=eta_remove)
      CALL get_mc_par(mc_par(insert_box)%mc_par,eta=eta_insert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! testing
!      remove_box=2
!      insert_box=1

! now choose a molecule type at random
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)
      DO itype=1,nmol_types
         IF(rand .LT. pmswap_mol(itype)) THEN
            molecule_type=itype
            EXIT
         ENDIF
      ENDDO

! record the attempt for the box the particle is to be inserted into
      moves(molecule_type,insert_box)%moves%swap%attempts=&
         moves(molecule_type,insert_box)%moves%swap%attempts+1

! now choose a random molecule to remove from the removal box, checking
! to make sure the box isn't empty
      IF(nchains(molecule_type,remove_box) == 0) THEN
         loverlap=.TRUE.
         moves(molecule_type,insert_box)%moves%empty=&
            moves(molecule_type,insert_box)%moves%empty+1
      ELSE

         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
         imolecule=CEILING(rand*nchains(molecule_type,remove_box))
! figure out the atom number this molecule starts on
         start_atom_rem=1
         DO itype=1,nmol_types
            IF(itype == molecule_type) THEN
               start_atom_rem=start_atom_rem+(imolecule-1)*nunits(itype)
               EXIT
            ELSE
               start_atom_rem=start_atom_rem+nchains(itype,remove_box)*nunits(itype)
            ENDIF
         ENDDO

! check for overlap
         start_mol=1
         DO jbox=1,remove_box-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO
         end_mol=start_mol+SUM(nchains(:,remove_box))-1
         CALL check_for_overlap(force_env(remove_box)%force_env,&
            nchains(:,remove_box),nunits,loverlap,mol_type(start_mol:end_mol))
         IF(loverlap) CALL stop_program(routineN,moduleN,__LINE__,&
            'CBMC swap move found an overlap in the old remove config')
         start_mol=1
         DO jbox=1,insert_box-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO
         end_mol=start_mol+SUM(nchains(:,insert_box))-1
         CALL check_for_overlap(force_env(insert_box)%force_env,&
            nchains(:,insert_box),nunits,loverlap,mol_type(start_mol:end_mol))
         IF(loverlap) CALL stop_program(routineN,moduleN,__LINE__,&
            'CBMC swap move found an overlap in the old insert config')
      ENDIF

      IF(loverlap) THEN
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"particles_old")
         CALL timestop(handle)
         RETURN
      ENDIF

! figure out how many atoms will be in each box after the move
      ins_atoms=nunits_tot(insert_box)+nunits(molecule_type)
      rem_atoms=nunits_tot(remove_box)-nunits(molecule_type)
! now allocate the arrays that will hold the coordinates and the
! atom name, for writing to the dat file
      IF(rem_atoms == 0) THEN
         ALLOCATE(remove_coords(1:3,1:nunits(1)),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "remove_coords",nunits(1)*3*dp_size)
         ALLOCATE(atom_names_remove(1:nunits(1)),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "atom_names_remove",nunits(1)*dp_size)
      ELSE
         ALLOCATE(remove_coords(1:3,1:rem_atoms),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "remove_coords",rem_atoms*3*dp_size)
         ALLOCATE(atom_names_remove(1:rem_atoms),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "atom_names_remove",rem_atoms*dp_size)
      ENDIF
      ALLOCATE(insert_coords(1:3,1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "insert_coords",ins_atoms*3*dp_size)
      ALLOCATE(atom_names_insert(1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "atom_names_insert",MAX(ins_atoms,rem_atoms)*dp_size)

! grab the cells for later...acceptance and insertion
      IF(lbias) THEN
         CALL force_env_get(bias_env(insert_box)%force_env,&
            cell=cell_insert,error=error)
         CALL force_env_get(bias_env(remove_box)%force_env,&
            cell=cell_remove,error=error)
      ELSE
         CALL force_env_get(force_env(insert_box)%force_env,&
            cell=cell_insert,error=error)
         CALL force_env_get(force_env(remove_box)%force_env,&
            cell=cell_remove,error=error)
      ENDIF
      CALL get_cell(cell_remove,abc=abc_remove)
      CALL get_cell(cell_insert,abc=abc_insert)

      IF(ionode) THEN
! choose an insertion point
         DO idim=1,3
            rand=next_random_number(rng_stream,error=error)
            pos_insert(idim)=rand*abc_insert(1)
         ENDDO
      ENDIF
      CALL mp_bcast(pos_insert,source,group)

! allocate some arrays we'll be using
      ALLOCATE(r_insert_mol(1:3,1:nunits(molecule_type)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r_insert_mol",3*nunits(molecule_type))

      iiatom=1
      DO iatom=start_atom_rem,start_atom_rem+nunits(molecule_type)-1
         r_insert_mol(1:3,iiatom)=&
            particles_old(remove_box)%list%els(iatom)%r(1:3)
         iiatom=iiatom+1
      ENDDO

!     find the center of mass of the molecule
      CALL get_center_of_mass(r_insert_mol(:,:),nunits(molecule_type),&
      center_of_mass(:),mass(:,molecule_type))

!     move the center of mass to the insertion point
      displace_molecule(1:3)=pos_insert(1:3)-center_of_mass(1:3)
      DO iatom=1,nunits(molecule_type)
         r_insert_mol(1:3,iatom)=r_insert_mol(1:3,iatom)+&
            displace_molecule(1:3)
      ENDDO

! prepare the insertion coordinates to be written to the .dat file so
! we can create a new force environment...remember there is still a particle
! in the box even if nchain=0
      IF (SUM(nchains(:,insert_box)) == 0) THEN
         DO iatom=1,nunits(molecule_type)
            insert_coords(1:3,iatom)=r_insert_mol(1:3,iatom)
            atom_names_insert(iatom)=&
                 particles_old(remove_box)%list%els(start_atom_rem+iatom-1)%atomic_kind%name
         ENDDO
         start_atom_ins=1
      ELSE
! the problem is I can't just tack the new molecule on to the end,
! because of reading in the dat_file...the topology stuff will crash
! if the molecules aren't all grouped together, so I have to insert it
! at the end of the section of molecules with the same type, then
! remember the start number for the CBMC stuff
         start_atom_ins=1
         DO itype=1,nmol_types
            start_atom_ins=start_atom_ins+&
                 nchains(itype,insert_box)*nunits(itype)
            IF(itype == molecule_type) EXIT
         ENDDO

         DO iatom=1,start_atom_ins-1
            insert_coords(1:3,iatom)=&
               particles_old(insert_box)%list%els(iatom)%r(1:3)
            atom_names_insert(iatom)=&
                 particles_old(insert_box)%list%els(iatom)%atomic_kind%name
         ENDDO
         iiatom=1
         DO iatom=start_atom_ins,start_atom_ins+nunits(molecule_type)-1
            insert_coords(1:3,iatom)=r_insert_mol(1:3,iiatom)
            atom_names_insert(iatom)=atom_names(iiatom,molecule_type)
            iiatom=iiatom+1
         ENDDO
         DO iatom=start_atom_ins+nunits(molecule_type),ins_atoms
            insert_coords(1:3,iatom)=&
               particles_old(insert_box)%list%els(iatom-nunits(molecule_type))%r(1:3)
            atom_names_insert(iatom)=&
                 particles_old(insert_box)%list%els(iatom-nunits(molecule_type))%atomic_kind%name
         ENDDO
      ENDIF

! fold the coordinates into the box and check for overlaps
      start_mol=1
      DO jbox=1,insert_box-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,insert_box))-1

! make the .dat file
      IF (ionode) THEN

         nchains(molecule_type,insert_box)=nchains(molecule_type,insert_box)+1
         IF(lbias) THEN
            CALL get_mc_par(mc_par(insert_box)%mc_par,mc_bias_file=mc_bias_file)
            CALL mc_make_dat_file_new(insert_coords(:,:),atom_names_insert(:),ins_atoms,&
               abc_insert(:),dat_file(insert_box),nchains(:,insert_box),&
               mc_bias_file)
         ELSE
            CALL get_mc_par(mc_par(insert_box)%mc_par,mc_input_file=mc_input_file)
            CALL mc_make_dat_file_new(insert_coords(:,:),atom_names_insert(:),ins_atoms,&
                abc_insert(:),dat_file(insert_box),nchains(:,insert_box),&
                mc_input_file)
         ENDIF
         nchains(molecule_type,insert_box)=nchains(molecule_type,insert_box)-1

      ENDIF

! now do the same for the removal box...be careful not to make an empty box
      IF (rem_atoms == 0) THEN
         DO iatom=1,nunits(molecule_type)
            remove_coords(1:3,iatom)=r_insert_mol(1:3,iatom)
            atom_names_remove(iatom)=atom_names(iatom,molecule_type)
         ENDDO

! need to adjust nchains, because otherwise if we are removing a molecule type
! that is not the first molecule, the dat file will have two molecules in it but
! only the coordinates for one
         nchains(molecule_type,remove_box)=nchains(molecule_type,remove_box)-1
         IF(ionode) THEN
            IF(lbias) THEN
               CALL get_mc_par(mc_par(remove_box)%mc_par,mc_bias_file=mc_bias_file)
               CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                    abc_remove(:),dat_file(remove_box),nchains(:,remove_box),&
                    mc_bias_file)
            ELSE
               CALL get_mc_par(mc_par(remove_box)%mc_par,mc_input_file=mc_input_file)
               CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                    abc_remove(:),dat_file(remove_box),nchains(:,remove_box),&
                    mc_input_file)
            ENDIF

         ENDIF
         nchains(molecule_type,remove_box)=nchains(molecule_type,remove_box)+1

      ELSE
         DO iatom=1,start_atom_rem-1
            remove_coords(1:3,iatom)=&
                 particles_old(remove_box)%list%els(iatom)%r(1:3)
            atom_names_remove(iatom)=&
                 particles_old(remove_box)%list%els(iatom)%atomic_kind%name
         ENDDO
         DO iatom=start_atom_rem+nunits(molecule_type),nunits_tot(remove_box)
            remove_coords(1:3,iatom-nunits(molecule_type))=&
                 particles_old(remove_box)%list%els(iatom)%r(1:3)
            atom_names_remove(iatom-nunits(molecule_type))=&
                 particles_old(remove_box)%list%els(iatom)%atomic_kind%name
         ENDDO

! make the .dat file
         IF (ionode) THEN
            nchains(molecule_type,remove_box)=nchains(molecule_type,remove_box)-1
            IF(lbias) THEN
               CALL get_mc_par(mc_par(remove_box)%mc_par,mc_bias_file=mc_bias_file)
               CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                    abc_remove(:),dat_file(remove_box),nchains(:,remove_box),&
                    mc_bias_file)
            ELSE
               CALL get_mc_par(mc_par(remove_box)%mc_par,mc_input_file=mc_input_file)
               CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                    abc_remove(:),dat_file(remove_box),nchains(:,remove_box),&
                    mc_input_file)
            ENDIF
            nchains(molecule_type,remove_box)=nchains(molecule_type,remove_box)+1

         ENDIF
      ENDIF

! deallocate r_insert_mol
      DEALLOCATE(r_insert_mol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r_insert_mol")

! now let's create the two new environments with the different number
! of molecules
      ALLOCATE(test_env(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
      "test_env",3*ins_atoms*dp_size)
      CALL mc_create_force_env(test_env(insert_box)%force_env, para_env, &
           dat_file(insert_box),error=error)
      CALL mc_create_force_env(test_env(remove_box)%force_env, para_env, &
           dat_file(remove_box),error=error)

! allocate an array we'll need
      ALLOCATE(r_cbmc(1:3,1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
      "r_cbmc",3*ins_atoms*dp_size)
      ALLOCATE(cbmc_energies(1:nswapmoves,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
      "cbmc_energies",nswapmoves*dp_size)

      loverlap_ins=.FALSE.
      loverlap_rem=.FALSE.

! compute the new molecule information...we need this for the CBMC part
      IF(rem_atoms == 0) THEN
         CALL mc_determine_molecule_info(test_env,mc_molecule_info_test,error,&
              box_number=remove_box)
      ELSE
         CALL mc_determine_molecule_info(test_env,mc_molecule_info_test,error)
      ENDIF
      CALL get_mc_molecule_info(mc_molecule_info_test,nchains=nchains_test,&
           mol_type=mol_type_test)

! figure out the position of the molecule we're inserting, and the
! Rosenbluth weight
      start_mol=1
      DO jbox=1,insert_box-1
         start_mol=start_mol+SUM(nchains_test(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains_test(:,insert_box))-1

      IF(lbias) THEN
         CALL generate_cbmc_swap_config(test_env(insert_box)%force_env,&
            BETA,max_val, min_val, exp_max_val,&
            exp_min_val,nswapmoves,weight_new,start_atom_ins,ins_atoms,nunits(:),&
            nunits(molecule_type),mass(:,molecule_type),loverlap_ins,bias_energy_new(insert_box),&
            bias_energy_old(insert_box),ionode,.FALSE.,mol_type_test(start_mol:end_mol),&
            nchains_test(:,insert_box),source,group,rng_stream,error)

! the energy that comes out of the above routine is the difference...we want
! the real energy for the acceptance rule...we don't do this for the
! lbias=.false. case because it doesn't appear in the acceptance rule, and
! we compensate in case of acceptance
            bias_energy_new(insert_box)=bias_energy_new(insert_box)+&
               bias_energy_old(insert_box)
      ELSE
         CALL generate_cbmc_swap_config(test_env(insert_box)%force_env,&
            BETA,max_val, min_val, exp_max_val,&
            exp_min_val,nswapmoves,weight_new,start_atom_ins,ins_atoms,nunits(:),&
            nunits(molecule_type),mass(:,molecule_type),loverlap_ins,new_energy(insert_box),&
            old_energy(insert_box),ionode,.FALSE.,mol_type_test(start_mol:end_mol),&
            nchains_test(:,insert_box),source,group,rng_stream,error)
      ENDIF

      CALL force_env_get(test_env(insert_box)%force_env,&
         subsys=insert_sys,error=error)
      CALL cp_subsys_get(insert_sys, &
         particles=particles_insert,error=error)

      DO iatom=1,ins_atoms
         r_cbmc(1:3,iatom)=particles_insert%els(iatom)%r(1:3)
      ENDDO

! make sure there is no overlap

      IF(loverlap_ins .OR. loverlap_rem) THEN
! deallocate some stuff
         CALL mc_molecule_info_destroy(mc_molecule_info_test)
         CALL section_vals_release(test_env(insert_box)%force_env%root_section,error=error)
         CALL section_vals_release(test_env(remove_box)%force_env%root_section,error=error)
         CALL force_env_release(test_env(insert_box)%force_env,error=error)
         CALL force_env_release(test_env(remove_box)%force_env,error=error)
         DEALLOCATE(insert_coords,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"insert_coords")
         DEALLOCATE(remove_coords,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"remove_coords")
         DEALLOCATE(r_cbmc,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"r_cbmc")
         DEALLOCATE(cbmc_energies,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"cbmc_energies")
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
            __LINE__,"particles_old")
         DEALLOCATE(test_env,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
              __LINE__,"test_env")
         CALL timestop(handle)
         RETURN
      ENDIF

! broadcast the choosen coordiantes to all processors

      CALL force_env_get(test_env(insert_box)%force_env,&
         subsys=insert_sys,error=error)
      CALL cp_subsys_get(insert_sys, &
         particles=particles_insert,error=error)

      DO iatom=1,ins_atoms
         particles_insert%els(iatom)%r(1:3)= &
             r_cbmc(1:3,iatom)
      ENDDO

! if we made it this far, we have no overlaps
      moves(molecule_type,insert_box)%moves%grown=&
         moves(molecule_type,insert_box)%moves%grown+1

! if we're biasing, we need to make environments with the non-biasing
! potentials, and calculate the energies
      IF(lbias) THEN

         ALLOCATE(test_env_bias(1:2),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "test_env_bias",2)

! first, the environment to which we added a molecule
         CALL get_mc_par(mc_par(insert_box)%mc_par,mc_input_file=mc_input_file)
          IF(ionode) CALL mc_make_dat_file_new(r_cbmc(:,:),atom_names_insert(:),ins_atoms,&
              abc_insert(:),dat_file(insert_box),nchains_test(:,insert_box),&
              mc_input_file)
         test_env_bias(insert_box)%force_env => test_env(insert_box)%force_env
         NULLIFY(test_env(insert_box)%force_env)
         CALL mc_create_force_env(test_env(insert_box)%force_env, para_env, &
              dat_file(insert_box),error=error)
!         CALL section_vals_release(insert_env%root_section)

         CALL force_env_calc_energy_force(test_env(insert_box)%force_env,&
            calc_force=.FALSE.,error=error)
         CALL force_env_get(test_env(insert_box)%force_env,&
            potential_energy=new_energy(insert_box),error=error)

! now the environment that has one less molecule
         IF (SUM(nchains_test(:,remove_box)) == 0) THEN
            CALL get_mc_par(mc_par(remove_box)%mc_par,mc_input_file=mc_input_file)
             IF(ionode) CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                 abc_remove(:),dat_file(remove_box),nchains_test(:,remove_box),&
                 mc_input_file)
             test_env_bias(remove_box)%force_env => test_env(remove_box)%force_env
            NULLIFY(test_env(remove_box)%force_env)
            CALL mc_create_force_env(test_env(remove_box)%force_env, para_env, &
                 dat_file(remove_box),error=error)
            new_energy(remove_box)=0.0E0_dp
            bias_energy_new(remove_box)=0.0E0_dp
         ELSE
            CALL get_mc_par(mc_par(remove_box)%mc_par,mc_input_file=mc_input_file)
            IF(ionode) CALL mc_make_dat_file_new(remove_coords(:,:),atom_names_remove(:),rem_atoms,&
                 abc_remove(:),dat_file(remove_box),nchains_test(:,remove_box),&
                 mc_input_file)
            test_env_bias(remove_box)%force_env => test_env(remove_box)%force_env
            NULLIFY(test_env(remove_box)%force_env)
            CALL mc_create_force_env(test_env(remove_box)%force_env, para_env, &
                 dat_file(remove_box),error=error)
            CALL force_env_calc_energy_force(test_env(remove_box)%force_env,&
                 calc_force=.FALSE.,error=error)
            CALL force_env_get(test_env(remove_box)%force_env,&
               potential_energy=new_energy(remove_box),error=error)
            CALL force_env_calc_energy_force(test_env_bias(remove_box)%force_env,&
               calc_force=.FALSE.,error=error)
            CALL force_env_get(test_env_bias(remove_box)%force_env,&
               potential_energy=bias_energy_new(remove_box),error=error)
         ENDIF
      ELSE
         IF (SUM(nchains_test(:,remove_box)) == 0) THEN
            new_energy(remove_box)=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(test_env(remove_box)%force_env,&
                 calc_force=.FALSE.,error=error)
            CALL force_env_get(test_env(remove_box)%force_env,&
               potential_energy=new_energy(remove_box),error=error)
         ENDIF
      ENDIF

! now we need to figure out the rosenbluth weight for the old configuration...
! we wait until now to do that because we need the energy of the box that
! has had a molecule removed...notice we use the environment that has not
! had a molecule removed for the CBMC configurations, and therefore nchains
! and mol_type instead of nchains_test and mol_type_test
      start_mol=1
      DO jbox=1,remove_box-1
         start_mol=start_mol+SUM(nchains(:,jbox))
      ENDDO
      end_mol=start_mol+SUM(nchains(:,remove_box))-1
      IF(lbias) THEN
          CALL generate_cbmc_swap_config(bias_env(remove_box)%force_env,&
            BETA,max_val, min_val, exp_max_val,&
            exp_min_val,nswapmoves,weight_old,start_atom_rem,nunits_tot(remove_box),&
            nunits(:),nunits(molecule_type),mass(:,molecule_type),loverlap_rem,rdum,&
            bias_energy_new(remove_box),ionode,.TRUE.,mol_type(start_mol:end_mol),&
            nchains(:,remove_box),source,group,rng_stream,error)
      ELSE
         CALL generate_cbmc_swap_config(force_env(remove_box)%force_env,&
            BETA,max_val, min_val, exp_max_val,&
            exp_min_val,nswapmoves,weight_old,start_atom_rem,nunits_tot(remove_box),&
            nunits(:),nunits(molecule_type),mass(:,molecule_type),loverlap_rem,rdum,&
            new_energy(remove_box),ionode,.TRUE.,mol_type(start_mol:end_mol),&
            nchains(:,remove_box),source,group,rng_stream,error)
      ENDIF

! figure out the prefactor to the boltzmann weight in the acceptance
! rule, based on numbers of particles and volumes

      prefactor=REAL(nchains(molecule_type,remove_box),dp)/&
      REAL(nchains(molecule_type,insert_box)+1,dp)*&
      abc_insert(1)**3/abc_remove(1)**3

      IF(lbias) THEN

         del_quickstep_energy=(-BETA)*(new_energy(insert_box)-&
         old_energy(insert_box)+new_energy(remove_box)-&
         old_energy(remove_box)-(bias_energy_new(insert_box)+&
         bias_energy_new(remove_box)-bias_energy_old(insert_box)&
         -bias_energy_old(remove_box)))

         IF    (del_quickstep_energy .GT. exp_max_val) THEN
            del_quickstep_energy=max_val
         ELSEIF(del_quickstep_energy .LT. exp_min_val) THEN
            del_quickstep_energy=min_val
         ELSE
            del_quickstep_energy=EXP(del_quickstep_energy)
         ENDIF
         w=prefactor*del_quickstep_energy*weight_new/weight_old &
                *EXP(BETA*(eta_remove(molecule_type)-eta_insert(molecule_type)))

      ELSE
         w=prefactor*weight_new/weight_old &
                *EXP(BETA*(eta_remove(molecule_type)-eta_insert(molecule_type)))

      ENDIF

! check if the move is accepted
      IF(w .GE. 1.0E0_dp) THEN
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( rand .LT. w ) THEN

! accept the move

!     accept the move
         moves(molecule_type,insert_box)%moves%swap%successes=&
         moves(molecule_type,insert_box)%moves%swap%successes+1

! we need to compensate for the fact that we take the difference in
! generate_cbmc_config to keep the exponetials small
         IF(.NOT. lbias) THEN
            new_energy(insert_box)=new_energy(insert_box)+&
            old_energy(insert_box)
         ENDIF

         DO ibox=1,2
! update energies
            energy_check(ibox)=energy_check(ibox)+(new_energy(ibox)-&
            old_energy(ibox))
            old_energy(ibox)=new_energy(ibox)
! if we're biasing the update the biasing energy
            IF (lbias) THEN
               last_bias_energy(ibox)=bias_energy_new(ibox)
               bias_energy_old(ibox)=bias_energy_new(ibox)
            ENDIF

         ENDDO

! change particle numbers...basically destory the old mc_molecule_info and attach
! the new stuff to the mc_pars
! figure out the molecule information for the new environments
         CALL mc_molecule_info_destroy(mc_molecule_info)
         CALL set_mc_par(mc_par(insert_box)%mc_par,mc_molecule_info=mc_molecule_info_test)
         CALL set_mc_par(mc_par(remove_box)%mc_par,mc_molecule_info=mc_molecule_info_test)

! update coordinates
         CALL force_env_get(test_env(insert_box)%force_env,&
            subsys=insert_sys,error=error)
         CALL cp_subsys_get(insert_sys, &
            particles=particles_insert,error=error)
         DO ipart=1,ins_atoms
            r_old(1:3,ipart,insert_box)=particles_insert%els(ipart)%r(1:3)
         ENDDO
         CALL force_env_get(test_env(remove_box)%force_env,&
            subsys=remove_sys,error=error)
         CALL cp_subsys_get(remove_sys, &
            particles=particles_remove,error=error)
         DO ipart=1,rem_atoms
            r_old(1:3,ipart,remove_box)=particles_remove%els(ipart)%r(1:3)
         ENDDO

         ! Only the root_section of the second subforce_env needs to be released
         ! The other is pointing to the main force_env and must not be released

         ! insertion box
         IF (insert_box/=1) THEN
            CALL section_vals_release(test_env(insert_box)%force_env%root_section,error=error)
         END IF
         CALL section_vals_release(force_env(insert_box)%force_env%root_section,error=error)
         CALL force_env_release(force_env(insert_box)%force_env,error=error)
         force_env(insert_box)%force_env => test_env(insert_box)%force_env

         ! removal box
         IF (remove_box/=1) THEN
            CALL section_vals_release(test_env(remove_box)%force_env%root_section,error=error)
         END IF
         CALL section_vals_release(force_env(remove_box)%force_env%root_section,error=error)
         CALL force_env_release(force_env(remove_box)%force_env,error=error)
         force_env(remove_box)%force_env => test_env(remove_box)%force_env

! if we're biasing, update the bias_env
         IF(lbias) THEN
            CALL section_vals_release(bias_env(insert_box)%force_env%root_section,error=error)
            CALL section_vals_release(bias_env(remove_box)%force_env%root_section,error=error)

            CALL force_env_release(bias_env(insert_box)%force_env,error=error)
            bias_env(insert_box)%force_env =>  test_env_bias(insert_box)%force_env
            CALL force_env_release(bias_env(remove_box)%force_env,error=error)
            bias_env(remove_box)%force_env => test_env_bias(remove_box)%force_env
            DEALLOCATE(test_env_bias,STAT=istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
                 __LINE__,"test_env_bias")
         ENDIF

      ELSE

! reject the move
         CALL mc_molecule_info_destroy(mc_molecule_info_test)
         CALL section_vals_release(test_env(insert_box)%force_env%root_section,error=error)
         CALL section_vals_release(test_env(remove_box)%force_env%root_section,error=error)
         CALL force_env_release(test_env(insert_box)%force_env,error=error)
         CALL force_env_release(test_env(remove_box)%force_env,error=error)
         IF(lbias) THEN
            CALL section_vals_release(test_env_bias(insert_box)%force_env%root_section,error=error)
            CALL section_vals_release(test_env_bias(remove_box)%force_env%root_section,error=error)
            CALL force_env_release( test_env_bias(insert_box)%force_env,error=error)
            CALL force_env_release( test_env_bias(remove_box)%force_env,error=error)
            DEALLOCATE(test_env_bias,STAT=istat)
            IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
                 __LINE__,"test_env_bias")
         ENDIF
      ENDIF

! deallocate some stuff
      DEALLOCATE(insert_coords,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
           __LINE__,"insert_coords")
      DEALLOCATE(remove_coords,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
           __LINE__,"remove_coords")
      DEALLOCATE(test_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
           __LINE__,"test_env")
      DEALLOCATE(cbmc_energies,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"cbmc_energies")
      DEALLOCATE(r_cbmc,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,&
         __LINE__,"r_cbmc")
      DEALLOCATE(oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "oldsys")
      DEALLOCATE(particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_old")

! end the timing
      CALL timestop(handle)

  END SUBROUTINE mc_ge_swap_move

! *****************************************************************************
!> \brief performs a Monte Carlo move that alters the volume of the simulation boxes,
!>      keeping the total volume of the two boxes the same
!> \param mc_par the mc parameters for the force env
!> \param force_env the force environments used in the move
!> \param moves the structure that keeps track of how many moves have been
!>               accepted/rejected
!> \param move_updates the structure that keeps track of how many moves have
!>               been accepted/rejected since the last time the displacements
!>               were updated
!> \param nnstep the total number of Monte Carlo moves already performed
!> \param old_energy the energy of the last accepted move involving an
!>                    unbiased potential calculation
!> \param energy_check the running total of how much the energy has changed
!>                      since the initial configuration
!> \param r_old the coordinates of the last accepted move involving a
!>               Quickstep calculation
!> \param rng_stream the stream we pull random numbers from
!> \param error the cp_error_type in case of problems
!>
!>      Designed for parallel.
!> \author MJM
! *****************************************************************************
  SUBROUTINE mc_ge_volume_move ( mc_par,force_env, moves,move_updates,&
      nnstep,old_energy,energy_check,r_old,rng_stream,error)

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    TYPE(mc_moves_p_type), DIMENSION(:, :), &
      POINTER                                :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: old_energy, energy_check
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: r_old
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mc_ge_volume_move', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: fft_lib
    CHARACTER(LEN=40), DIMENSION(1:2)        :: dat_file
    INTEGER :: cl, end_atom, end_mol, group, handle, iatom, ibox, imolecule, &
      iside, istat, j, jatom, jbox, max_atoms, molecule_index, molecule_type, &
      print_level, source, start_atom, start_mol
    INTEGER, DIMENSION(:), POINTER           :: mol_type, nunits, nunits_tot
    INTEGER, DIMENSION(:, :), POINTER        :: nchains
    LOGICAL                                  :: ionode
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: loverlap
    LOGICAL, DIMENSION(1:2)                  :: lempty
    REAL(dp), DIMENSION(:, :), POINTER       :: mass
    REAL(KIND=dp)                            :: BETA, prefactor, rand, &
                                                rmvolume, vol_dis, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: r
    REAL(KIND=dp), DIMENSION(1:2)            :: new_energy, volume_new, &
                                                volume_old
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                center_of_mass_new, diff
    REAL(KIND=dp), DIMENSION(1:3, 1:2)       :: abc, new_cell_length, &
                                                old_cell_length
    REAL(KIND=dp), DIMENSION(1:3, 1:3, 1:2)  :: hmat_test
    TYPE(cell_p_type), DIMENSION(:), POINTER :: cell, cell_old, cell_test
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: oldsys
    TYPE(mc_molecule_info_type), POINTER     :: mc_molecule_info
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old

! begin the timing of the subroutine

      CALL timeset(routineN,handle)

! nullify some pointers
      NULLIFY(particles_old,cell,oldsys,cell_old,cell_test)

! get some data from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,source=source,&
         group=group,dat_file=dat_file(1),rmvolume=rmvolume,&
         BETA=BETA,cl=cl,fft_lib=fft_lib,&
         mc_molecule_info=mc_molecule_info)
      CALL get_mc_molecule_info(mc_molecule_info,nunits_tot=nunits_tot,&
           mass=mass,nchains=nchains,nunits=nunits,mol_type=mol_type)

      print_level = 1
      CALL get_mc_par(mc_par(2)%mc_par,dat_file=dat_file(2))

! allocate some stuff
      max_atoms=MAX(nunits_tot(1),nunits_tot(2))
      ALLOCATE(r(1:3,max_atoms,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r",2*3*max_atoms*dp_size)
      ALLOCATE(oldsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "oldsys",2)
      ALLOCATE(particles_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_old",2)
      ALLOCATE(cell(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell",2)
      ALLOCATE(cell_test(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell_test",2)
      ALLOCATE(cell_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell_old",2)
      ALLOCATE(loverlap(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "loverlap",2)

! check for empty boxes...need to be careful because we can't build
! a force_env with no particles
      DO ibox=1,2
         lempty(ibox)=.FALSE.
         IF(SUM(nchains(:,ibox))==0) THEN
            lempty(ibox)=.TRUE.
         ENDIF
      ENDDO

! record the attempt
      DO ibox=1,2
         moves(1,ibox)%moves%volume%attempts= &
            moves(1,ibox)%moves%volume%attempts+1
         move_updates(1,ibox)%moves%volume%attempts=&
            move_updates(1,ibox)%moves%volume%attempts+1
      ENDDO

! now let's grab the cell length and particle positions
      DO ibox=1,2
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=oldsys(ibox)%subsys,cell=cell(ibox)%cell,error=error)
         CALL get_cell(cell(ibox)%cell,abc=abc(:,ibox))
         NULLIFY(cell_old(ibox)%cell)
         CALL cell_create(cell_old(ibox)%cell,error=error)
         CALL cell_clone(cell(ibox)%cell,cell_old(ibox)%cell,error=error)
         CALL cp_subsys_get(oldsys(ibox)%subsys, &
            particles=particles_old(ibox)%list, error=error)

! find the old cell length
         old_cell_length(1:3,ibox)=abc(1:3,ibox)

      ENDDO

      DO ibox=1,2

! save the old coordiantes
         DO iatom=1,nunits_tot(ibox)
            r(1:3,iatom,ibox)=particles_old(ibox)%list%els(iatom)%r(1:3)
         ENDDO

      ENDDO

! call a random number to figure out how far we're moving
      IF(ionode) rand=next_random_number(rng_stream,error=error)
      CALL mp_bcast(rand,source,group)

      vol_dis=rmvolume*(rand-0.5E0_dp)*2.0E0_dp

! add to one box, subtract from the other
      IF(old_cell_length(1,1)*old_cell_length(2,1)*&
         old_cell_length(3,1)+vol_dis .LE. (3.0E0_dp/angstrom)**3) &
         CALL stop_program(routineN,moduleN,__LINE__,&
         'GE_volume moves are trying to make box 1 smaller than 3')
      IF(old_cell_length(1,2)*old_cell_length(2,2)*&
         old_cell_length(3,2)+vol_dis .LE. (3.0E0_dp/angstrom)**3) &
         CALL stop_program(routineN,moduleN,__LINE__,&
         'GE_volume moves are trying to make box 2 smaller than 3')

      DO iside=1,3
         new_cell_length(iside,1)=(old_cell_length(1,1)**3+&
            vol_dis)**(1.0E0_dp/3.0E0_dp)
         new_cell_length(iside,2)=(old_cell_length(1,2)**3-&
            vol_dis)**(1.0E0_dp/3.0E0_dp)
      ENDDO

! now we need to make the new cells
      DO ibox=1,2
         hmat_test(:,:,ibox)=0.0e0_dp
         hmat_test(1,1,ibox)=new_cell_length(1,ibox)
         hmat_test(2,2,ibox)=new_cell_length(2,ibox)
         hmat_test(3,3,ibox)=new_cell_length(3,ibox)
         NULLIFY (cell_test(ibox)%cell)
         CALL cell_create(cell_test(ibox)%cell,hmat=hmat_test(:,:,ibox),&
              periodic=cell(ibox)%cell%perd,error=error)
         CALL force_env_set_cell(force_env(ibox)%force_env,&
              cell_test(ibox)%cell,error=error)
      ENDDO

      DO ibox=1,2

! save the coords
         DO iatom=1,nunits_tot(ibox)
            r(1:3,iatom,ibox)=particles_old(ibox)%list%els(iatom)%r(1:3)
         ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass
         start_atom=1
         molecule_index=1
         DO jbox=1,ibox-1
            IF(jbox == ibox) EXIT
            molecule_index=molecule_index+SUM(nchains(:,jbox))
         ENDDO
         DO imolecule=1,SUM(nchains(:,ibox))
            molecule_type=mol_type(imolecule+molecule_index-1)
            IF(imolecule .NE. 1) THEN
               start_atom=start_atom+nunits(mol_type(imolecule+molecule_index-2))
            ENDIF
            end_atom=start_atom+nunits(molecule_type)-1

! now find the center of mass
            CALL get_center_of_mass(r(:,start_atom:end_atom,ibox),&
               nunits(molecule_type),center_of_mass(:),mass(:,molecule_type))

! scale the center of mass and determine the vector that points from the
!    old COM to the new one
            center_of_mass_new(1:3)=center_of_mass(1:3)*&
               new_cell_length(1:3,ibox)/old_cell_length(1:3,ibox)
            DO j=1,3
               diff(j)=center_of_mass_new(j)-center_of_mass(j)
! now change the particle positions
               DO jatom=start_atom,end_atom
                  particles_old(ibox)%list%els(jatom)%r(j)=&
                     particles_old(ibox)%list%els(jatom)%r(j)+diff(j)
               ENDDO

            ENDDO
         ENDDO

! check for any overlaps we might have
         start_mol=1
         DO jbox=1,ibox-1
            start_mol=start_mol+SUM(nchains(:,jbox))
         ENDDO
         end_mol=start_mol+SUM(nchains(:,ibox))-1
         CALL check_for_overlap(force_env(ibox)%force_env,&
            nchains(:,ibox),nunits,loverlap(ibox),mol_type(start_mol:end_mol),&
            cell_length=new_cell_length(:,ibox))

      ENDDO

! determine the overall energy difference

      DO ibox=1,2
         IF(loverlap(ibox)) CYCLE
! remake the force environment and calculate the energy
         IF (lempty(ibox)) THEN
            new_energy(ibox)=0.0E0_dp
         ELSE

            CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
               calc_force=.FALSE.,error=error)
            CALL force_env_get(force_env(ibox)%force_env,&
               potential_energy=new_energy(ibox),error=error)

         ENDIF
      ENDDO

! accept or reject the move
      DO ibox=1,2
         volume_new(ibox)=new_cell_length(1,ibox)*&
            new_cell_length(2,ibox)*new_cell_length(3,ibox)
         volume_old(ibox)=old_cell_length(1,ibox)*&
            old_cell_length(2,ibox)*old_cell_length(3,ibox)
      ENDDO
      prefactor=(volume_new(1)/volume_old(1))**(SUM(nchains(:,1)))*&
         (volume_new(2)/volume_old(2))**(SUM(nchains(:,2)))

      IF( loverlap(1) .OR. loverlap(2)) THEN
         w=0.0E0_dp
      ELSE
         w=prefactor*EXP(-BETA*&
              (new_energy(1)+new_energy(2)-&
              old_energy(1)-old_energy(2)))

      ENDIF

      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream,error=error)
         CALL mp_bcast(rand,source,group)
      ENDIF

      IF (rand .LT. w) THEN

! write cell length, volume, density, and trial displacement to a file
         IF(ionode) THEN

            WRITE(cl,*) nnstep,new_cell_length(1,1)*&
               angstrom,vol_dis*(angstrom)**3,new_cell_length(1,2)*&
               angstrom
            WRITE(cl,*) nnstep,new_energy(1),&
               old_energy(1),new_energy(2),old_energy(2)
            WRITE(cl,*) prefactor,w
         ENDIF

         DO ibox=1,2
! accept the move
            moves(1,ibox)%moves%volume%successes=&
               moves(1,ibox)%moves%volume%successes+1
            move_updates(1,ibox)%moves%volume%successes=&
               move_updates(1,ibox)%moves%volume%successes+1

! update energies
            energy_check(ibox)=energy_check(ibox)+(new_energy(ibox)-&
               old_energy(ibox))
            old_energy(ibox)=new_energy(ibox)

! and the new "old" coordiantes
            DO iatom=1,nunits_tot(ibox)
               r_old(1:3,iatom,ibox)=&
                  particles_old(ibox)%list%els(iatom)%r(1:3)
            ENDDO

         ENDDO
      ELSE

! reject the move
! write cell length, volume, density, and trial displacement to a file
         IF(ionode) THEN

            WRITE(cl,*) nnstep,new_cell_length(1,1)*&
               angstrom,vol_dis*(angstrom)**3,new_cell_length(1,2)*&
               angstrom
            WRITE(cl,*) nnstep,new_energy(1),&
               old_energy(1),new_energy(2),old_energy(2)
            WRITE(cl,*) prefactor,w

         ENDIF

! reset the cell and particle positions
         DO ibox=1,2
            CALL force_env_set_cell(force_env(ibox)%force_env,&
                 cell_old(ibox)%cell,error=error)
            DO iatom=1,nunits_tot(ibox)
               particles_old(ibox)%list%els(iatom)%r(1:3)=r_old(1:3,iatom,ibox)
            ENDDO
         ENDDO

      ENDIF

! free up some memory
      DO ibox=1,2
         CALL cell_release(cell_test(ibox)%cell,error=error)
         CALL cell_release(cell_old(ibox)%cell,error=error)
      ENDDO
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "r")
      DEALLOCATE(oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "oldsys")
      DEALLOCATE(particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "particles_old")
      DEALLOCATE(cell,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell")
      DEALLOCATE(cell_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell_old")
      DEALLOCATE(cell_test,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "cell_test")
      DEALLOCATE(loverlap,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "loverlap")

! end the timing
      CALL timestop(handle)

  END SUBROUTINE mc_ge_volume_move

END MODULE mc_ge_moves

