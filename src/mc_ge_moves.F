!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_ge_moves [1.0] *
!!
!!   NAME
!!     mc_ge_moves
!!
!!   FUNCTION
!!     contains the Monte Carlo moves that can handle more than one
!!     box, including the Quickstep move, a volume swap between boxes, 
!      and a particle swap between boxes
!!
!!   AUTHOR
!!     Matthew J. McGrath  (01.25.2004)
!!
!!   MODIFICATION HISTORY
!!     MJM (07.28.2005): make the Quickstep move general, and changed
!!                       the swap and volume moves to work with the
!!                       CP2K classical routines
!!
!!   SOURCE
!******************************************************************************

MODULE mc_ge_moves
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_pp_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mc_control,                      ONLY: mc_create_force_env
  USE mc_coordinates,                  ONLY: generate_cbmc_swap_config,&
                                             get_center_of_mass,&
!                                             mc_classical_energy,&
                                             check_for_overlap,&
                                             mc_coordinate_fold
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_move_control,                 ONLY: move_q_reinit,&
                                             q_move_accept
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_moves_p_type,&
                                             mc_simulation_parameters_p_type,&
                                             set_mc_par
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_p_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: rng_stream_type,&
                                             next_random_number
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE simulation_cell,                 ONLY: cell_p_type,&
                                             cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ge_moves"

  PUBLIC :: mc_ge_volume_move,mc_ge_swap_move,&
      mc_quickstep_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_ge_moves.F/mc_quickstep_move
!!
!!   NAME
!!     mc_quickstep_move
!!
!!   FUNCTION
!!     computes the acceptance of a series of biased or unbiased moves 
!!     (translation, rotation, conformational changes)
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force envs of the boxes
!!     - force_env: the force environments for the boxes
!!     - bias_env: the force environments with the biasing potential for the boxes 
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected for both boxes
!!     - lreject: automatically rejects the move (used when an overlap occurs in
!!               the sequence of moves)
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated for both boxes
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move before the sequence
!!        whose acceptance is determined by this call
!!     - nnstep: the Monte Carlo step we're on
!!     - old_energy: the energy of the last accepted move involving the full potential
!!     - bias_energy_new: the energy of the current configuration involving the bias potential
!!     - nboxes: the number of boxes (force environments) in the system
!!     - box_flag: indicates if a move has been tried in a given box..if not, we don't
!!       recompute the energy
!!     - subsys: the pointers for the particle subsystems of both boxes
!!     - particles: the pointers for the particle sets
!!     - cell: the pointers for the cells
!!     - rng_stream: the stream we pull random numbers from
!!     - error: the cp_error_type in case something goes wrong
!!
!!     Designed for parallel use.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
  SUBROUTINE mc_Quickstep_move(mc_par,force_env,bias_env,moves,&
                        lreject,move_updates,energy_check,r_old,&
                        nnstep,old_energy,bias_energy_new,&
                        nboxes,box_flag,subsys,particles,cell,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:),POINTER                     :: mc_par
    TYPE(force_env_p_type),DIMENSION(:), POINTER &
                                        :: force_env,bias_env
    TYPE(mc_moves_p_type),DIMENSION(:), POINTER    &
                                         :: moves, move_updates
    TYPE(cell_p_type),POINTER,DIMENSION(:) &
                                             :: cell
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: subsys,subsys_bias
    TYPE(particle_list_p_type), DIMENSION(:),POINTER &
                                             :: particles,particles_bias
    REAL(KIND=dp), DIMENSION(:,:,:), &
      INTENT(INOUT)                          :: r_old
    INTEGER, INTENT(IN)                      :: nnstep,nboxes
    INTEGER, DIMENSION(:),INTENT(IN)         :: box_flag
    REAL(KIND=dp), DIMENSION(:),INTENT(INOUT) &
                                             :: old_energy,&
                                                bias_energy_new,&
                                                energy_check
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    LOGICAL,INTENT(IN)                       :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_Quickstep_move"

!   *** Local Variables ***
    INTEGER,DIMENSION(1:nboxes)              :: diff,nchain,nunits_tot
    INTEGER                                  :: handle, iparticle,&
      source,group,iprint,ibox,istat
    REAL(KIND=dp)                            :: rand, w,&
      BETA,energies
    REAL(KIND=dp),DIMENSION(1:nboxes)        :: bias_energy_old,&
      new_energy
    LOGICAL                                  :: lbias,ionode,&
      loverlap

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      NULLIFY(subsys_bias,particles_bias)

! get a bunch of data from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,lbias=lbias,&
         BETA=BETA,diff=diff(1),source=source,group=group,&
         nunits_tot=nunits_tot(1),nchain=nchain(1),iprint=iprint)

      IF(nboxes .GT. 1) THEN
         DO ibox=2,nboxes
            CALL get_mc_par(mc_par(ibox)%mc_par,nchain=nchain(ibox),&
               diff=diff(ibox),nunits_tot=nunits_tot(ibox)) 
         ENDDO
      ENDIF

! allocate some stuff
      ALLOCATE(subsys_bias(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "subsys_bias",nboxes)
      ALLOCATE(particles_bias(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_bias",nboxes)

! record the attempt
      DO ibox=1,nboxes
         moves(ibox)%moves%Quickstep%attempts=&
            moves(ibox)%moves%Quickstep%attempts+1
      ENDDO

! grab the coordinates for the force_env
      DO ibox=1,nboxes
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=subsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(subsys(ibox)%subsys(1)%subsys, &
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
               CALL cp_subsys_get(subsys_bias(ibox)%subsys(1)%subsys, &
                  particles=particles_bias(ibox)%list, error=error)
               
               DO iparticle=1,nunits_tot(ibox)
                  particles(ibox)%list%els(iparticle)%r(1:3)=&
                     particles_bias(ibox)%list%els(iparticle)%r(1:3)
               ENDDO

               CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
                  calc_force=.FALSE.)
               CALL force_env_get(force_env(ibox)%force_env,&
                  potential_energy=new_energy(ibox),error=error)
            ELSE
               IF( .NOT. lreject) THEN
                  CALL force_env_calc_energy_force(force_env(ibox)%force_env,&
                     calc_force=.FALSE.)
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
               IF(nchain(ibox) == 0) THEN
                  WRITE(diff(ibox),*) nnstep
                  WRITE(diff(ibox),*) nchain(ibox)
               ELSE
                  WRITE(diff(ibox),*) nnstep
                  CALL write_particle_coordinates(&
                     particles(ibox)%list%els,cell(ibox)%cell,&
                     diff(ibox),dump_xmol,'POS','TRIAL')
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF(.NOT. lreject) THEN
         IF (lbias) THEN

            DO ibox=1,nboxes
! look for overlap
               IF(nchain(ibox) .NE. 0) THEN
                  CALL check_for_overlap(bias_env(ibox)%force_env,&
                     nchain(ibox),1,nchain(ibox),loverlap)
                  IF(loverlap) CALL stop_program(routineN,&
                     module_name,__LINE__,&
                     'Quickstep move found an overlap in the old config') 
               ENDIF
               bias_energy_old(ibox)=moves(ibox)%moves%bias_energy
            ENDDO

            energies=-BETA*((SUM(new_energy(:))-SUM(bias_energy_new(:)))&
               -(SUM(old_energy(:))-SUM(bias_energy_old(:))))

            
! used to prevent over and underflows
            IF(energies .GE. 0.0_dp) THEN
               w=1.0_dp
            ELSEIF(energies .LE. -500.0_dp) THEN
               w=0.0_dp
            ELSE
               w=exp(energies)
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
               w=exp(energies)
            ENDIF
!               WRITE(6,*) 'w = ',w
         ENDIF
      ELSE
         w=0.0E0_dp
      ENDIF
      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream)
!     CALL RANDOM_NUMBER(rand)
         CALL mp_bcast(rand,source,group)      
      ENDIF

      IF (rand .LT. w) THEN

         DO ibox=1,nboxes
! accept the move
            moves(ibox)%moves%Quickstep%successes=&
               moves(ibox)%moves%Quickstep%successes+1

! remember what kind of move we did for lbias=.false.
            IF(.NOT. lbias) THEN
               CALL q_move_accept(moves(ibox)%moves,.TRUE.)
               CALL q_move_accept(move_updates(ibox)%moves,.TRUE.)
                  
! reset the counters
               CALL move_q_reinit(moves(ibox)%moves, .TRUE.)
               CALL move_q_reinit(move_updates(ibox)%moves, .TRUE.)

            ENDIF

! we need to record all accepted moves since last Quickstep calculation
            CALL q_move_accept(moves(ibox)%moves,.FALSE.)
            CALL q_move_accept(move_updates(ibox)%moves,.FALSE.)

! reset the counters
            CALL move_q_reinit(moves(ibox)%moves, .FALSE.)
            CALL move_q_reinit(move_updates(ibox)%moves, .FALSE.)

! update energies
            energy_check(ibox)=energy_check(ibox)+&
               (new_energy(ibox)-old_energy(ibox))
            old_energy(ibox)=new_energy(ibox)

         ENDDO

         IF ( lbias) THEN
            DO ibox=1,nboxes
               moves(ibox)%moves%bias_energy=bias_energy_new(ibox)
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
            CALL move_q_reinit(moves(ibox)%moves,.FALSE.)
            CALL move_q_reinit(move_updates(ibox)%moves,.FALSE.)
            IF(.NOT. lbias) THEN
! reset the counters
               CALL move_q_reinit(moves(ibox)%moves, .TRUE.)
               CALL move_q_reinit(move_updates(ibox)%moves, .TRUE.)
            ENDIF

         ENDDO


         IF ( .NOT. ionode) r_old(:,:,:) = 0.0E0_dp

! coodinates changed, so we need to broadcast those, even for the lbias 
! case since bias_env needs to have the same coords as force_env
         CALL mp_bcast(r_old,source,group)         
         
         DO ibox=1,nboxes
            DO iparticle=1,nunits_tot(ibox)
               particles(ibox)%list%els(iparticle)%r(1:3)=&
                  r_old(1:3,iparticle,ibox)
               IF(lbias) &
                  particles_bias(ibox)%list%els(iparticle)%r(1:3)=&
                  r_old(1:3,iparticle,ibox)
            ENDDO
         ENDDO
      ENDIF

! make sure the coordinates are transferred
      DO ibox=1,nboxes
         CALL cp_subsys_set(subsys(ibox)%subsys(1)%subsys,&
            particles=particles(ibox)%list)
         IF(lbias) &
            CALL cp_subsys_set(subsys_bias(ibox)%subsys(1)%subsys,&
            particles=particles_bias(ibox)%list)
      ENDDO

 ! deallocate some stuff
      DEALLOCATE(subsys_bias,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"subsys_bias")
      DEALLOCATE(particles_bias,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"particles_bias")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_Quickstep_move

! *****************************************************************************

! *****************************************************************************
!!****s* mc_ge_moves.F/mc_ge_swap_move
!!
!!   NAME
!!     mc_ge_swap_move
!!
!!   FUNCTION
!!     attempts a swap move between two simulation boxes
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force envs of the boxes
!!     - force_env: the force environments for the boxes
!!     - bias_env: the force environments used to bias moves for the boxes
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected for both boxes
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move involving a
!!              full potential calculation for both boxes
!!     - old_energy: the energy of the last accepted move involving a
!!                   a full potential calculation
!!     - nnstep: the Monte Carlo step we're on
!!     - globenv: the global environment for the whole simulation
!!     - bias_energy_old: the energies of both boxes computed using the biasing
!!       potential
!!     - rng_stream: the stream we pull random numbers from
!!     - error: the cp_error_type in case something goes wrong
!!
!!     Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_ge_swap_move(mc_par,force_env,bias_env,moves,&
                energy_check,r_old,old_energy,&
                globenv,bias_energy_old,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env,bias_env
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: moves
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: energy_check
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: r_old
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: old_energy,bias_energy_old
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_ge_swap_move"

    INTEGER :: handle, i, iatom, ibox, ii, imolecule, ins_atoms, &
      insert_box, istat, natom, rem_atoms, remove_box,nswapmoves,&
      group,source,ipart
    LOGICAL                                  :: loverlap,ionode,lbias,&
      loverlap_ins,loverlap_rem
    REAL(KIND=dp)                            :: prefactor, rand, w,BETA,&
                                del_quickstep_energy,rdum,&
                                weight_old,weight_new
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_cbmc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:,:)                        :: cbmc_energies
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_insert, r_insert_mol, &
                                                r_remove
    INTEGER, DIMENSION(1:2)                  :: nchain,nunits_tot
    REAL(KIND=dp), DIMENSION(1:2)            :: bias_energy_new, &
                                                new_energy
    REAL(KIND=dp), DIMENSION(:),ALLOCATABLE  :: mass
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                displace_molecule, &
                                                pos_insert,abc_insert,&
                                                abc_remove
    TYPE(mol_kind_new_list_p_type), DIMENSION(:),POINTER  &
                                             :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(cell_type), POINTER                 :: cell_insert, cell_remove
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: insert_sys,remove_sys
    TYPE(force_env_type), POINTER            :: insert_env, remove_env,&
      insert_env_bias,remove_env_bias
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_old
    TYPE(particle_list_type),POINTER         :: particles_insert,&
      particles_remove
    CHARACTER(LEN=40),DIMENSION(1:2)         :: dat_file
    CHARACTER(LEN=200)                       :: fft_lib

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! reset the overlap flag
      loverlap=.FALSE.

! nullify some pointers
      NULLIFY(particles_old)
      NULLIFY(oldsys)

! grab some stuff from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,BETA=BETA,&
         nswapmoves=nswapmoves,group=group,source=source,&
         lbias=lbias,nchain=nchain(1),nunits_tot=nunits_tot(1),&
         dat_file=dat_file(1),fft_lib=fft_lib)
      CALL get_mc_par(mc_par(2)%mc_par,nchain=nchain(2),&
         nunits_tot=nunits_tot(2),dat_file=dat_file(2))

! allocate some stuff
      ALLOCATE(oldsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys",2)
      ALLOCATE(particles_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old",2)
      ALLOCATE(molecule_kinds_new(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecule_kinds_new",2)

! get the old coordinates
      DO ibox=1,2
         CALL force_env_get(force_env(ibox)%force_env,&
             subsys=oldsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(oldsys(ibox)%subsys(1)%subsys, &
            molecule_kinds_new=molecule_kinds_new(ibox)%list,&
            particles=particles_old(ibox)%list, &
            error=error)
      ENDDO

!     choose a direction to swap
      IF(ionode) rand=next_random_number(rng_stream)
!         CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)

      IF ( rand .LE. 0.50E0_dp ) THEN
         remove_box=1
         insert_box=2
      ELSE
         remove_box=2
         insert_box=1
      ENDIF

! record the attempt for the box the particle is to be inserted into
      moves(insert_box)%moves%swap%attempts=&
         moves(insert_box)%moves%swap%attempts+1
         
! now choose a random molecule to remove from the removal box, checking
! to make sure the box isn't empty
      IF(nchain(remove_box) == 0) THEN
         loverlap=.TRUE.
         moves(insert_box)%moves%empty=&
            moves(insert_box)%moves%empty+1
      ELSE

         IF(ionode) rand=next_random_number(rng_stream)
!            CALL RANDOM_NUMBER(rand)    
         CALL mp_bcast(rand,source,group)
         imolecule=CEILING(rand*nchain(remove_box))

! check for overlap
         CALL check_for_overlap(force_env(remove_box)%force_env,&
            nchain(remove_box),1,nchain(remove_box),loverlap)
         IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
            'CBMC swap move found an overlap in the old remove config')   
         CALL check_for_overlap(force_env(insert_box)%force_env,&
            nchain(insert_box),1,nchain(insert_box),loverlap)
         IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
            'CBMC swap move found an overlap in the old insert config')   
      ENDIF

      IF(loverlap) THEN
         DEALLOCATE(molecule_kinds_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecule_kinds_new")
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_old")
         CALL timestop(0.0E0_dp,handle)
         RETURN
      ENDIF

! allocate a couple of arrays that depend on the box
      molecule_kind => molecule_kinds_new(1)%list%els(1)
      CALL get_molecule_kind(molecule_kind,natom=natom)   
      ALLOCATE(mass(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass",natom*dp_size)
      DO iatom=1,natom
         CALL get_atomic_kind(particles_old(1)%list%els(iatom)%atomic_kind,&
            mass=mass(iatom))
      ENDDO
      ins_atoms=nunits_tot(insert_box)+natom
      rem_atoms=nunits_tot(remove_box)-natom

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
         DO i=1,3
!            CALL RANDOM_NUMBER(rand)
            rand=next_random_number(rng_stream)
            pos_insert(i)=rand*abc_insert(1)
         ENDDO
      ENDIF
      CALL mp_bcast(pos_insert,source,group)

! allocate some arrays we'll be using
!      CALL stop_program("natom value not defined here. Maybe ?natoms?",&
!                        "mc_ge_moves.F")
      ALLOCATE(r_insert_mol(1:3,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert_mol",3*natom)
      ALLOCATE(r_insert(1:3,1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert",3*ins_atoms*dp_size)

      DO i=imolecule*natom-natom+1,imolecule*natom
         ii=i-imolecule*natom+natom
         r_insert_mol(1:3,ii)=&
            particles_old(remove_box)%list%els(i)%r(1:3)
      ENDDO
         
!     find the center of mass of the molecule
      CALL get_center_of_mass(r_insert_mol(:,:),natom,&
      center_of_mass(:),mass(:))
         
!     move the center of mass to the insertion point
      displace_molecule(1:3)=pos_insert(1:3)-center_of_mass(1:3)
      DO iatom=1,natom
         r_insert_mol(1:3,iatom)=r_insert_mol(1:3,iatom)+&
            displace_molecule(1:3)
      ENDDO
         
! prepare the insertion coordinates to be written to the .dat file so
! we can create a new force environment...remember there is still a particle
! in the box even if nchain=0
      IF (nchain(insert_box) == 0) THEN
         DO iatom=1,natom
            r_insert(1:3,iatom)=r_insert_mol(1:3,iatom)
         ENDDO
      ELSE
         DO iatom=1,ins_atoms-natom
            r_insert(1:3,iatom)=&
               particles_old(insert_box)%list%els(iatom)%r(1:3)
         ENDDO
         DO iatom=ins_atoms-natom+1,ins_atoms
            ii=iatom-(ins_atoms-natom)
            r_insert(1:3,iatom)=r_insert_mol(1:3,ii)
         ENDDO
      ENDIF
      
! fold the coordinates into the box and check for overlaps
      CALL mc_coordinate_fold(r_insert(:,:),ins_atoms,mass,natom,&
         abc_insert(1:3))
      
! allocate r_remove
      IF (rem_atoms == 0) THEN
         ALLOCATE(r_remove(1:3,1:natom),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove",3*natom*dp_size)
      ELSE
         ALLOCATE(r_remove(1:3,1:rem_atoms),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove",3*rem_atoms*dp_size)
      ENDIF
      
! make the .dat file
      IF (ionode) THEN

         IF(lbias) THEN
            CALL mc_make_dat_file(r_insert(:,:),ins_atoms,&
               abc_insert(1:3)*angstrom,dat_file(insert_box),&
               bias_env(insert_box)%force_env,TRIM(fft_lib))
         ELSE
            CALL mc_make_dat_file(r_insert(:,:),ins_atoms,&
               abc_insert(1:3)*angstrom,dat_file(insert_box),&
               force_env(insert_box)%force_env,TRIM(fft_lib))
         ENDIF
      ENDIF
! now do the same for the removal box...be careful not to make an empty box
      IF (rem_atoms == 0) THEN
         DO iatom=1,natom
            r_remove(1:3,iatom)=r_insert_mol(1:3,iatom)
         ENDDO
!     make the .dat file
         IF(ionode) THEN
            IF(lbias) THEN
               CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
                  abc_remove(1:3)*angstrom,dat_file(remove_box),&
                  bias_env(remove_box)%force_env,TRIM(fft_lib))
            ELSE
               CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
                  abc_remove(1:3)*angstrom,dat_file(remove_box),&
                  force_env(remove_box)%force_env,TRIM(fft_lib))
            ENDIF
         ENDIF
      ELSE
         IF(imolecule == 1) THEN
            DO iatom=natom+1,rem_atoms+natom
               ii=iatom-natom
               r_remove(1:3,ii)=&
                  particles_old(remove_box)%list%els(iatom)%r(1:3)
            ENDDO               
         ELSE
            DO iatom=1,imolecule*natom-natom
               r_remove(1:3,iatom)=&
                  particles_old(remove_box)%list%els(iatom)%r(1:3)
            ENDDO
            DO iatom=imolecule*natom+1,rem_atoms+natom
               ii=iatom-natom
               r_remove(1:3,ii)=&
                  particles_old(remove_box)%list%els(iatom)%r(1:3)
            ENDDO
         ENDIF
         
! make the .dat file
         IF (ionode) THEN
            IF(lbias) THEN
               CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
                  abc_remove(1:3)*angstrom,dat_file(remove_box),&
                  bias_env(remove_box)%force_env,TRIM(fft_lib))
            ELSE
               CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
                  abc_remove(1:3)*angstrom,dat_file(remove_box),&
                  force_env(remove_box)%force_env,TRIM(fft_lib))
            ENDIF
         ENDIF
      ENDIF

! deallocate r_remove,r_insert, and r_insert_mol
      DEALLOCATE(r_insert,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r_insert")
      DEALLOCATE(r_insert_mol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert_mol")

! now let's create the two new environments with the different number
! of molecules

      IF(lbias) THEN
         globenv%input_file_name=dat_file(insert_box)
         CALL mc_create_force_env(insert_env, globenv, .FALSE., &
            bias_env(insert_box)%force_env%in_use,error=error)
         globenv%input_file_name=dat_file(remove_box)
         CALL mc_create_force_env(remove_env, globenv, .FALSE., &
            bias_env(remove_box)%force_env%in_use,error=error)
      ELSE
         globenv%input_file_name=dat_file(insert_box)
         CALL mc_create_force_env(insert_env, globenv, .FALSE., &
            force_env(insert_box)%force_env%in_use,error=error)
         globenv%input_file_name=dat_file(remove_box)
         CALL mc_create_force_env(remove_env, globenv, .FALSE., &
            force_env(remove_box)%force_env%in_use,error=error)
      ENDIF

! allocate an array we'll need
      ALLOCATE(r_cbmc(1:3,1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
      "r_cbmc",3*ins_atoms*dp_size)
      ALLOCATE(cbmc_energies(1:nswapmoves,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
      "cbmc_energies",nswapmoves*dp_size)
      
      loverlap_ins=.FALSE.
      loverlap_rem=.FALSE.


! figure out the position of the molecule we're inserting, and the
! Rosenbluth weight
      IF(lbias) THEN
         CALL generate_cbmc_swap_config(insert_env,&
            BETA,nswapmoves,weight_new,nchain(insert_box)+1,loverlap_ins,&
            bias_energy_new(insert_box),&
            bias_energy_old(insert_box),ionode,source,group,rng_stream)

! the energy that comes out of the above routine is the difference...we want
! the real energy for the acceptance rule...we don't do this for the 
! lbias=.false. case because it doesn't appear in the acceptance rule, and
! we compensate in case of acceptance
            bias_energy_new(insert_box)=bias_energy_new(insert_box)+&
               bias_energy_old(insert_box)
      ELSE
         CALL generate_cbmc_swap_config(insert_env,&
            BETA,nswapmoves,weight_new,nchain(insert_box)+1,&
            loverlap_ins,new_energy(insert_box),&
            old_energy(insert_box),ionode,source,group,&
            rng_stream)
      ENDIF
         
      CALL force_env_get(insert_env,&
         subsys=insert_sys,error=error)
      CALL cp_subsys_get(insert_sys(1)%subsys, &
         particles=particles_insert,error=error)
         
      DO iatom=1,ins_atoms
         r_cbmc(1:3,iatom)=particles_insert%els(iatom)%r(1:3)
      ENDDO
         
! make sure there is no overlap

      IF(loverlap_ins .OR. loverlap_rem) THEN
! deallocate some stuff
         CALL force_env_release(insert_env)
         CALL force_env_release(remove_env)
         DEALLOCATE(r_remove,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove")
         DEALLOCATE(mass,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"mass")
         DEALLOCATE(molecule_kinds_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecule_kinds_new")
         DEALLOCATE(r_cbmc,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_cbmc")
         DEALLOCATE(cbmc_energies,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"cbmc_energies")
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_old")
         CALL timestop(0.0E0_dp,handle)
         RETURN
      ENDIF

! broadcast the choosen coordiantes to all processors

      CALL force_env_get(insert_env,&
         subsys=insert_sys,error=error)
      CALL cp_subsys_get(insert_sys(1)%subsys, &
         particles=particles_insert,error=error)
      
      DO iatom=1,ins_atoms
         particles_insert%els(iatom)%r(1:3)= &
             r_cbmc(1:3,iatom)
      ENDDO
      
! if we made it this far, we have no overlaps
      moves(insert_box)%moves%grown=&
         moves(insert_box)%moves%grown+1

! if we're biasing, we need to make environments with the non-biasing
! potentials, and calculate the energies
      IF(lbias) THEN

! first, the environment to which we added a molecule
         IF(ionode) CALL mc_make_dat_file(r_cbmc(:,:),ins_atoms,&
            abc_insert(1:3)*angstrom,dat_file(insert_box),&
            force_env(insert_box)%force_env,TRIM(fft_lib))
         insert_env_bias => insert_env
         NULLIFY(insert_env)
!         CALL force_env_release(insert_env)
         globenv%input_file_name=dat_file(insert_box)
         CALL mc_create_force_env(insert_env, globenv, .FALSE., &
            force_env(insert_box)%force_env%in_use,error=error)
         CALL force_env_calc_energy_force(insert_env,&
            calc_force=.FALSE.)
         CALL force_env_get(insert_env,&
            potential_energy=new_energy(insert_box),error=error)

! now the environment that has one less molecule
         IF (nchain(remove_box)-1 == 0) THEN
            IF(ionode) CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
               abc_remove(1:3)*angstrom,dat_file(remove_box),&
               force_env(remove_box)%force_env,TRIM(fft_lib))
            remove_env_bias => remove_env
            NULLIFY(remove_env)
            globenv%input_file_name=dat_file(remove_box)
            CALL mc_create_force_env(remove_env, globenv, .FALSE., &
               force_env(remove_box)%force_env%in_use,error=error)
            new_energy(remove_box)=0.0E0_dp
            bias_energy_new(remove_box)=0.0E0_dp
         ELSE
            IF(ionode) CALL mc_make_dat_file(r_remove(:,:),rem_atoms,&
               abc_remove(1:3)*angstrom,dat_file(remove_box),&
               force_env(remove_box)%force_env,TRIM(fft_lib))
            remove_env_bias => remove_env
            NULLIFY(remove_env)
            globenv%input_file_name=dat_file(remove_box)
            CALL mc_create_force_env(remove_env, globenv, .FALSE., &
               force_env(remove_box)%force_env%in_use,error=error)
            CALL force_env_calc_energy_force(remove_env,calc_force=.FALSE.)
            CALL force_env_get(remove_env,&
               potential_energy=new_energy(remove_box),error=error)
            CALL force_env_calc_energy_force(remove_env_bias,&
               calc_force=.FALSE.)
            CALL force_env_get(remove_env_bias,&
               potential_energy=bias_energy_new(remove_box),error=error)
         ENDIF
      ELSE
         IF (nchain(remove_box)-1 == 0) THEN
            new_energy(remove_box)=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(remove_env,calc_force=.FALSE.)
            CALL force_env_get(remove_env,&
               potential_energy=new_energy(remove_box),error=error)
         ENDIF
      ENDIF
      

! now we need to figure out the rosenbluth weight for the old configuration...
! we wait until now to do that because we need the energy of the box that
! has had a molecule removed
      IF(lbias) THEN
          CALL generate_cbmc_swap_config(bias_env(remove_box)%force_env,&
            BETA,nswapmoves,weight_old,nchain(remove_box),loverlap_rem,rdum,&
            bias_energy_new(remove_box),&
            ionode,source,group,rng_stream,rem_molecule=imolecule)
      ELSE
         CALL generate_cbmc_swap_config(force_env(remove_box)%force_env,&
            BETA,nswapmoves,weight_old,nchain(remove_box),loverlap_rem,rdum,&
            new_energy(remove_box),&
            ionode,source,group,rng_stream,rem_molecule=imolecule)
      ENDIF

! figure out the prefactor to the boltzmann weight in the acceptance
! rule, based on numbers of particles and volumes

      prefactor=REAL(nchain(remove_box),dp)/&
      REAL(nchain(insert_box)+1,dp)*&
      abc_insert(1)**3/abc_remove(1)**3
      
      IF(lbias) THEN
         
         del_quickstep_energy=(-BETA)*(new_energy(insert_box)-&
         old_energy(insert_box)+new_energy(remove_box)-&
         old_energy(remove_box)-(bias_energy_new(insert_box)+&
         bias_energy_new(remove_box)-bias_energy_old(insert_box)&
         -bias_energy_old(remove_box)))
         
         IF(del_quickstep_energy .GT. 500.0_dp) THEN
            del_quickstep_energy=1.4E217_dp
         ELSEIF(del_quickstep_energy .LT. -500.0_dp) THEN
            del_quickstep_energy=0.0_dp
         ELSE
            del_quickstep_energy=EXP(del_quickstep_energy)
         ENDIF
         
         w=prefactor*del_quickstep_energy*weight_new/weight_old
         
      ELSE

         w=prefactor*weight_new/weight_old

      ENDIF

! check if the move is accepted
      IF(w .GE. 1.0E0_dp) THEN
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream)
!     CALL RANDOM_NUMBER(rand)
         CALL mp_bcast(rand,source,group)
      ENDIF
         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      w=0.0d0
!      w=1.0d0

      IF ( rand .LT. w ) THEN

! accept the move         

!     accept the move
         moves(insert_box)%moves%swap%successes=&
         moves(insert_box)%moves%swap%successes+1

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
               moves(ibox)%moves%bias_energy=&
               bias_energy_new(ibox)
               bias_energy_old(ibox)=bias_energy_new(ibox)
            ENDIF

         ENDDO


! change particle numbers
! now change the total number of units and chains in the mc_par
         nchain(insert_box)=nchain(insert_box)+1
         nchain(remove_box)=nchain(remove_box)-1

! set the numbers in the mc_par structures
         CALL set_mc_par(mc_par(insert_box)%mc_par,&
            nchain=nchain(insert_box),nunits_tot=ins_atoms)
         CALL set_mc_par(mc_par(remove_box)%mc_par,&
            nchain=nchain(remove_box),nunits_tot=rem_atoms)

! update coordinates
         CALL force_env_get(insert_env,&
            subsys=insert_sys,error=error)
         CALL cp_subsys_get(insert_sys(1)%subsys, &
            particles=particles_insert,error=error)
         DO ipart=1,ins_atoms
            r_old(1:3,ipart,insert_box)=particles_insert%els(ipart)%r(1:3)
         ENDDO
         CALL force_env_get(remove_env,&
            subsys=remove_sys,error=error)
         CALL cp_subsys_get(remove_sys(1)%subsys, &
            particles=particles_remove,error=error)
         DO ipart=1,rem_atoms
            r_old(1:3,ipart,remove_box)=particles_remove%els(ipart)%r(1:3)
         ENDDO         

! insertion box
         CALL force_env_release(force_env(insert_box)%force_env)
         force_env(insert_box)%force_env => insert_env

! removal box
         CALL force_env_release(force_env(remove_box)%force_env)
         force_env(remove_box)%force_env => remove_env

! if we're biasing, update the bias_env
         IF(lbias) THEN
            CALL force_env_release(bias_env(insert_box)%force_env)
            bias_env(insert_box)%force_env => insert_env_bias
            CALL force_env_release(bias_env(remove_box)%force_env)
            bias_env(remove_box)%force_env => remove_env_bias
         ENDIF


!         WRITE(6,*) 'iuovnd ',new_energy(:)
!         WRITE(6,*) 'uovnd ',nchain(:)
      ELSE

!         WRITE(6,*) 'REJECTED'

! reject the move
         CALL force_env_release(insert_env)
         CALL force_env_release(remove_env)
         IF(lbias) THEN
            CALL force_env_release(insert_env_bias)
            CALL force_env_release(remove_env_bias)
         ENDIF
      ENDIF

! deallocate some stuff
      DEALLOCATE(cbmc_energies,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"cbmc_energies")
      DEALLOCATE(r_cbmc,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r_cbmc")
      DEALLOCATE(r_remove,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"r_remove")
      DEALLOCATE(mass,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"mass")
      DEALLOCATE(molecule_kinds_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"molecule_kinds_new")
      DEALLOCATE(oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys")
      DEALLOCATE(particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_ge_swap_move

! *****************************************************************************

! *****************************************************************************
!!****s* mc_ge_moves.F/mc_ge_volume_move
!!
!!   NAME
!!     mc_ge_volume_move
!!
!!   FUNCTION
!!     performs a Monte Carlo move that alters the volume of the simulation boxes,
!!     keeping the total volume of the two boxes the same
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environments used in the move
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of Monte Carlo moves already performed
!!     - old_energy: the energy of the last accepted move involving an
!!                   unbiased potential calculation
!!     - globenv: the global simulation parameters
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move involving a
!!              Quickstep calculation
!!     - rng_stream: the stream we pull random numbers from
!!     - error: the cp_error_type in case of problems
!!                                 
!!     Designed for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_ge_volume_move ( mc_par,force_env, moves,move_updates,&
      nnstep,old_energy,globenv,energy_check,r_old,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:),POINTER                   :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), POINTER &
                                             :: force_env
    TYPE(mc_moves_p_type), DIMENSION(:), POINTER &
                                             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), DIMENSION(:),INTENT(INOUT)&
                                             :: old_energy
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), DIMENSION(:),INTENT(INOUT)&
                                             :: energy_check
    REAL(KIND=dp), &
      DIMENSION(:, :,:), &
      INTENT(INOUT)                          :: r_old
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_ge_volume_move"

! local variables
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass 

    INTEGER                                  :: handle, i, iatom, ibox,&
                                                j,all_atoms,jatom,&
                                                istat,source,group,&
                                                nunits,cl,iside,natoms
    INTEGER,DIMENSION(1:2)                   :: nunits_tot,nchain
    REAL(KIND=dp),DIMENSION(1:2)             :: new_energy, &
      volume_new,volume_old
    REAL(KIND=dp), DIMENSION(1:3,1:2)        :: abc,old_cell_length,&
      new_cell_length
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                center_of_mass_new, diff
    TYPE(force_env_p_type), DIMENSION(:),POINTER &
                                             :: test_env
    TYPE(mol_kind_new_list_p_type), DIMENSION(:),POINTER  &
                                             :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys, testsys
    TYPE(cell_p_type), DIMENSION(:),POINTER  &
                                             :: cell_old
    REAL(KIND=dp), DIMENSION(:,:,:),ALLOCATABLE &
                                             :: r
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_test, particles_old
    REAL(KIND=dp)                            :: rand,w,prefactor,vol_dis,&
      rmvolume,BETA
    LOGICAL,DIMENSION(1:2)                   :: lempty
    LOGICAL                                  :: ionode
    LOGICAL,DIMENSION(:),ALLOCATABLE       :: loverlap
    CHARACTER(LEN=40),DIMENSION(1:2)         :: dat_file
    CHARACTER(LEN=200)                       :: fft_lib
!   ---------------------------------------------------------------------------
      
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles_old,particles_test,cell_old,oldsys,&
         testsys,test_env)

! get some data from mc_par
      CALL get_mc_par(mc_par(1)%mc_par,ionode=ionode,source=source,&
         group=group,nchain=nchain(1),nunits_tot=nunits_tot(1),&
         dat_file=dat_file(1),rmvolume=rmvolume,&
         nunits=nunits,BETA=BETA,cl=cl,fft_lib=fft_lib)
      CALL get_mc_par(mc_par(2)%mc_par,nchain=nchain(2),&
         nunits_tot=nunits_tot(2),dat_file=dat_file(2))

! allocate some stuff
      all_atoms=nunits_tot(1)+nunits_tot(2)
      ALLOCATE(r(1:3,all_atoms,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r",3*2*all_atoms)
      ALLOCATE(oldsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys",2)
      ALLOCATE(particles_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old",2)
      ALLOCATE(cell_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "cell_old",2)
      ALLOCATE(test_env(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "test_env",2)
      ALLOCATE(testsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "testsys",2)
      ALLOCATE(particles_test(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_test",2)      
      ALLOCATE(loverlap(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "loverlap",2)      
      ALLOCATE(molecule_kinds_new(2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecule_kinds_new",2)

! check for empty boxes...need to be careful because we can't build
! a force_env with no particles
      DO ibox=1,2
         lempty(ibox)=.FALSE.
         IF(nchain(ibox)==0) THEN
            lempty(ibox)=.TRUE.
         ENDIF
      ENDDO
      
! record the attempt
      DO ibox=1,2
         moves(ibox)%moves%volume%attempts= &
            moves(ibox)%moves%volume%attempts+1
         move_updates(ibox)%moves%volume%attempts=&
            move_updates(ibox)%moves%volume%attempts+1
      ENDDO

! now let's grab the cell length and particle positions
      DO ibox=1,2
         CALL force_env_get(force_env(ibox)%force_env,&
            subsys=oldsys(ibox)%subsys,cell=cell_old(ibox)%cell,error=error)
         CALL get_cell(cell_old(ibox)%cell,abc=abc(:,ibox))
         CALL cp_subsys_get(oldsys(ibox)%subsys(1)%subsys, &
            molecule_kinds_new=molecule_kinds_new(ibox)%list,&
            particles=particles_old(ibox)%list, error=error)

! find the old cell length
         old_cell_length(1:3,ibox)=abc(1:3,ibox)

      ENDDO

! allocate a couple of arrays that depend on the box
      molecule_kind => molecule_kinds_new(1)%list%els(1)
      CALL get_molecule_kind(molecule_kind,natom=natoms)   
      ALLOCATE(mass(1:natoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "mass",natoms*dp_size)
      DO iatom=1,natoms
         CALL get_atomic_kind(particles_old(1)%list%els(iatom)%atomic_kind,&
            mass=mass(iatom))
      ENDDO

      DO ibox=1,2
         
         IF(lempty(ibox)) THEN

! make the new dat file
            IF (ionode) CALL mc_make_dat_file(r(:,:,ibox),0,&
               old_cell_length(1:3,ibox)*angstrom,dat_file(ibox),&
               force_env(ibox)%force_env,TRIM(fft_lib))
            
         ELSE
            
! save the old coordiantes
            DO iatom=1,nunits_tot(ibox)
               r(1:3,iatom,ibox)=particles_old(ibox)%list%els(iatom)%r(1:3)
            ENDDO

! fold the coordinates into the central box
            CALL mc_coordinate_fold(r(:,:,ibox),&
               nunits_tot(ibox),mass,nunits,&
               old_cell_length(1:3,ibox))

! make the new dat file
            IF(ionode) CALL mc_make_dat_file(r(:,:,ibox),nunits_tot(ibox),&
               old_cell_length(1:3,ibox)*angstrom,dat_file(ibox),&
               force_env(ibox)%force_env,TRIM(fft_lib))

         ENDIF
         
      ENDDO

      DO ibox=1,2
! make a test environment that we can alter to calculate the energy
         globenv%input_file_name=dat_file(ibox)
         CALL mc_create_force_env(test_env(ibox)%force_env, globenv, .FALSE., &
               force_env(ibox)%force_env%in_use,error=error)

! grab the coords so we can change them
         CALL force_env_get(test_env(ibox)%force_env,&
            subsys=testsys(ibox)%subsys,error=error)
         CALL cp_subsys_get(testsys(ibox)%subsys(1)%subsys, &
            particles=particles_test(ibox)%list, error=error)

      ENDDO

! call a random number to figure out how far we're moving
!         CALL RANDOM_NUMBER(rand)
      IF(ionode) rand=next_random_number(rng_stream)
      CALL mp_bcast(rand,source,group)

      vol_dis=rmvolume*(rand-0.5E0_dp)*2.0E0_dp

! add to one box, subtract from the other
      IF(old_cell_length(1,1)*old_cell_length(2,1)*&
         old_cell_length(3,1)+vol_dis .LE. (3.0E0_dp/angstrom)**3) &
         CALL stop_program(routineN,module_name,__LINE__,&
         'GE_volume moves are trying to make box 1 smaller than 3')
      IF(old_cell_length(1,2)*old_cell_length(2,2)*&
         old_cell_length(3,2)+vol_dis .LE. (3.0E0_dp/angstrom)**3) &
         CALL stop_program(routineN,module_name,__LINE__,&
         'GE_volume moves are trying to make box 2 smaller than 3')
      
      DO iside=1,3
         new_cell_length(iside,1)=(old_cell_length(1,1)**3+&
            vol_dis)**(1.0E0_dp/3.0E0_dp)
         new_cell_length(iside,2)=(old_cell_length(1,2)**3-&
            vol_dis)**(1.0E0_dp/3.0E0_dp)
      ENDDO
      
      DO i=1,2

! save the coords
         DO iatom=1,nunits_tot(i)
            r(1:3,iatom,i)=particles_test(i)%list%els(iatom)%r(1:3)
         ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass
         DO iatom=1,nunits_tot(i),nunits

! now find the center of mass
            CALL get_center_of_mass(r(:,iatom:iatom+nunits-1,i),&
               nunits,center_of_mass(:),mass(:))

! scale the center of mass and determine the vector that points from the
!    old COM to the new one
            center_of_mass_new(1:3)=center_of_mass(1:3)*&
               new_cell_length(1:3,i)/old_cell_length(1:3,i)
            DO j=1,3
               diff(j)=center_of_mass_new(j)-center_of_mass(j)
! now change the particle positions
               DO jatom=iatom,iatom+nunits-1
                  particles_test(i)%list%els(jatom)%r(j)=&
                     particles_test(i)%list%els(jatom)%r(j)+diff(j)
               ENDDO
               
            ENDDO
         ENDDO
        
! check for any overlaps we might have
         CALL check_for_overlap(test_env(i)%force_env,&
            nchain(i),1,nchain(i),loverlap(i),&
            cell_length=new_cell_length(:,i))
               
      ENDDO

! determine the overall energy difference

      DO ibox=1,2

         IF(lempty(ibox)) THEN
               
! make the new dat file
            IF(ionode) CALL mc_make_dat_file(r(:,:,ibox),0,&
               new_cell_length(1:3,ibox)*angstrom,dat_file(ibox),&
               force_env(ibox)%force_env,TRIM(fft_lib))

         ELSE

! make the new dat file
            DO iatom=1,nunits_tot(ibox)
               r(1:3,iatom,ibox)=particles_test(ibox)%list%els(iatom)%r(1:3)
            ENDDO
               
            IF(ionode) CALL mc_make_dat_file(r(:,:,ibox),&
               nunits_tot(ibox),new_cell_length(1:3,ibox)*angstrom,&
               dat_file(ibox),force_env(ibox)%force_env,TRIM(fft_lib))
         ENDIF
      ENDDO

      DO ibox=1,2
! remake the force environment and calculate the energy

         CALL force_env_release(test_env(ibox)%force_env,error=error)
            
         globenv%input_file_name=dat_file(ibox)
         CALL mc_create_force_env(test_env(ibox)%force_env, globenv, .FALSE., &
               force_env(ibox)%force_env%in_use,error=error)
         
         IF (lempty(ibox)) THEN
            new_energy(ibox)=0.0E0_dp
         ELSE

            CALL force_env_calc_energy_force(test_env(ibox)%force_env,&
               calc_force=.FALSE.)
            CALL force_env_get(test_env(ibox)%force_env,&
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
      prefactor=(volume_new(1)/volume_old(1))**nchain(1)*&
         (volume_new(2)/volume_old(2))**nchain(2)
      w=prefactor*exp(-BETA*&
         (new_energy(1)+new_energy(2)-&
         old_energy(1)-old_energy(2)))
      
      IF(loverlap(1) .OR. loverlap(2)) w=0.0E0_dp
         

      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream)
!            CALL RANDOM_NUMBER(rand)
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
            moves(ibox)%moves%volume%successes=&
               moves(ibox)%moves%volume%successes+1
            move_updates(ibox)%moves%volume%successes=&
               move_updates(ibox)%moves%volume%successes+1

! update energies
            energy_check(ibox)=energy_check(ibox)+(new_energy(ibox)-&
               old_energy(ibox))
            old_energy(ibox)=new_energy(ibox)

! and the coordiantes
            CALL force_env_get(test_env(ibox)%force_env,subsys=&
               testsys(ibox)%subsys,error=error)
            CALL cp_subsys_get(testsys(ibox)%subsys(1)%subsys, &
               particles=particles_test(ibox)%list, error=error)
            
            DO iatom=1,nunits_tot(ibox)
               r_old(1:3,iatom,ibox)=&
                  particles_test(ibox)%list%els(iatom)%r(1:3)
            ENDDO
            
         ENDDO

! switch the force env to the new env
         DO ibox=1,2
            CALL force_env_release(force_env(ibox)%force_env,error=error)
            force_env(ibox)%force_env => test_env(ibox)%force_env
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

! release some memory
         DO i=1,2
            CALL force_env_release(test_env(i)%force_env,error=error)
         ENDDO

      ENDIF

! free up some memory
      DEALLOCATE(mass,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"mass")
      DEALLOCATE(molecule_kinds_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"molecule_kinds_new")
      DEALLOCATE(test_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,&
         __LINE__,"test_env")
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r")
      DEALLOCATE(oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys")
      DEALLOCATE(particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old")
      DEALLOCATE(cell_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "cell_old")
      DEALLOCATE(testsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "testsys")
      DEALLOCATE(particles_test,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_test")      
      DEALLOCATE(loverlap,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "loverlap")      

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_ge_volume_move


END MODULE mc_ge_moves

!**********************************************************************

