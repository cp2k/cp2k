!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_ge_moves [1.0] *
!!
!!   NAME
!!     mc_moves
!!
!!   FUNCTION
!!     the various moves in Gibbs ensemble Monte Carlo (GEMC) simulations,
!!     including the NVT volume move and a special Quickstep move
!!
!!   AUTHOR
!!     Matthew J. McGrath  (01.25.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_ge_moves
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_pp_type
  USE force_env_methods,               ONLY: force_env_calculate_energy
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mc_misc,                         ONLY: get_center_of_mass,&
                                             mc_classical_energy,&
                                             mc_coordinate_fold,&
                                             mc_make_dat_file
  USE mc_move_control,                 ONLY: move_q_reinit,&
                                             q_move_accept
  USE mc_types,                        ONLY: mc_averages_type,&
                                             mc_moves_p_type,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_p_type,&
                                             mc_simulation_parameters_type
  USE message_passing,                 ONLY: mp_bcast
  USE mol_new_list_types,              ONLY: mol_new_list_p_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE physcon,                         ONLY: angstrom
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type,&
                                             cell_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ge_moves"

  PUBLIC :: mc_ge_quickstep_move,mc_ge_volume_move,mc_ge_swap_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_ge_moves.F/mc_ge_quickstep_move
!!
!!   NAME
!!     mc_ge_quickstep_move
!!
!!   FUNCTION
!!     calls the Quickstep energy routines to calculate the energy of
!!     the current configuration and determine whether we accept or reject
!!     all the moves since the last time a Quickstep calculation was made...
!!     for two simulation boxes
!!
!!   ARGUMENTS
!!     - mc_par_1: the mc parameters for the force env of box 1
!!     - mc_par_2: the mc parameters for the force env of box 2
!!     - force_env_1: the force environment for box 1
!!     - force_env_2: the force environment for box 2
!!     - moves_1: the structure that keeps track of how many moves have been
!!              accepted/rejected for box 1
!!     - moves_2: the structure that keeps track of how many moves have been
!!              accepted/rejected for box 2
!!     - move_updates_1: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated for box 1
!!     - move_updates_2: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated for box 2
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old_1: the coordinates of the last accepted move involving a
!!              Quickstep calculation for box 1
!!     - r_old_2: the coordinates of the last accepted move involving a
!!              Quickstep calculation for box 2
!!     - istep: how many times the energy averages have been taken
!!     - old_energy: the energy of the last accepted move involving a
!!                   Quickstep calculation
!!     - averages_1: the structure that keeps track of the running averages
!!                 in the system for box 1
!!     - averages_1: the structure that keeps track of the running averages
!!                 in the system for box 1
!!     - box_flag: box_flag(i)=1 if a move has been tried in box i since the
!!                 last Quickstep move, 0 otherwise
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_ge_quickstep_move(mc_par_1,mc_par_2,force_env_1,&
                force_env_2,moves_1,moves_2,&
                move_updates_1,move_updates_2,&
                energy_check,r_old_1,r_old_2,old_energy,&
                box_flag,nnstep,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par_1, mc_par_2
    TYPE(force_env_type), POINTER            :: force_env_1, force_env_2
    TYPE(mc_moves_type), POINTER             :: moves_1, moves_2, &
                                                move_updates_1, move_updates_2
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old_1, r_old_2
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: old_energy
    INTEGER, DIMENSION(1:2), INTENT(IN)      :: box_flag
    INTEGER, INTENT(IN)                      :: nnstep
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_ge_quickstep_move"

    INTEGER                                  :: handle, i, iparticle
    REAL(KIND=dp)                            :: rand, w
    REAL(KIND=dp), DIMENSION(1:2)            :: classical_energy_new, &
                                                classical_energy_old, &
                                                new_energy
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_1, subsys_2
    TYPE(particle_list_type), POINTER        :: particles_1, particles_2

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles_1,subsys_1)
      NULLIFY(particles_2,subsys_2)

! in serial
      IF (mc_par_1%ionode) THEN

! record the attempt
         moves_1%Quickstep%attempts=moves_1%Quickstep%attempts+1
         moves_2%Quickstep%attempts=moves_2%Quickstep%attempts+1

      ENDIF

! broadcast the box flags so all procs do the same calculations
      CALL mp_bcast(box_flag,mc_par_1%source,mc_par_1%group)

! calculate the new energy of the system and grab the coordinates
      IF(box_flag(1) == 1) THEN
!!!!!!!!!!!!!!!!!!!!!
!         CALL force_env_calculate_energy(force_env_1)
         CALL force_env_get(force_env_1,subsys=subsys_1,&
            potential_energy=new_energy(1),error=error)
         CALL cp_subsys_get(subsys_1(1)%subsys, &
            particles=particles_1, error=error)
      ELSE
         new_energy(1)=old_energy(1)
         CALL force_env_get(force_env_1,subsys=subsys_1,&
            error=error)
         CALL cp_subsys_get(subsys_1(1)%subsys, &
            particles=particles_1,error=error)
      ENDIF
      IF(box_flag(2) == 1) THEN
!         CALL force_env_calculate_energy(force_env_2)
         CALL force_env_get(force_env_2,subsys=subsys_2,&
            potential_energy=new_energy(2),error=error)
         CALL cp_subsys_get(subsys_2(1)%subsys, &
            particles=particles_2, error=error)
      ELSE
         new_energy(2)=old_energy(2)
         CALL force_env_get(force_env_2,subsys=subsys_2,&
            error=error)
         CALL cp_subsys_get(subsys_2(1)%subsys, &
            particles=particles_2, error=error)
      ENDIF

      IF (mc_par_1%ionode) THEN

! accept or reject the move based on Metropolis or the Iftimie rule
         IF (mc_par_1%lclassical) THEN
            CALL mc_classical_energy(force_env_1,&
               classical_energy_new(1),&
               mc_par_1,1,mc_par_1%nchain,6.211D0/angstrom)
            new_energy(1)=classical_energy_new(1)
            CALL mc_classical_energy(force_env_2,&
               classical_energy_new(2),&
               mc_par_2,1,mc_par_2%nchain,6.211D0/angstrom)
            new_energy(2)=classical_energy_new(2)
            
            classical_energy_old(1)=moves_1%classical_energy
            classical_energy_old(2)=moves_2%classical_energy

            w=dexp(-mc_par_1%BETA*((new_energy(1)-&
               classical_energy_new(1))&
               -(old_energy(1)-classical_energy_old(1))&
               +(new_energy(2)-classical_energy_new(2))&
               -(old_energy(2)-classical_energy_old(2))))
            
            WRITE(mc_par_1%diff,*) nnstep,new_energy(1)-old_energy(1),&
                 classical_energy_new(1)-classical_energy_old(1)
            WRITE(mc_par_2%diff,*) nnstep,new_energy(2)-old_energy(2),&
                 classical_energy_new(2)-classical_energy_old(2)

         ELSE
            w=dexp(-mc_par_1%BETA*(new_energy(1)-old_energy(1)+&
               new_energy(2)-old_energy(2)))
         ENDIF

         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

! back to parallel
      ENDIF

! broadcast these numbers, to make sure we're on the same page
      CALL mp_bcast(rand,mc_par_1%source,mc_par_1%group)      
      CALL mp_bcast(w,mc_par_1%source,mc_par_1%group)     

      IF (rand .LT. w) THEN

         IF ( mc_par_1%ionode) THEN
            
! accept the move
            moves_1%Quickstep%successes=moves_1%Quickstep%successes+1
            moves_2%Quickstep%successes=moves_2%Quickstep%successes+1

            CALL q_move_accept(moves_1)
            CALL q_move_accept(moves_2)

! reset the counters
            CALL move_q_reinit(moves_1)
            CALL move_q_reinit(moves_2)

! update energies for the boxes that were displaced
            IF (box_flag(1) == 1) THEN
               energy_check(1)=energy_check(1)+(new_energy(1)-&
                                              old_energy(1))
               old_energy(1)=new_energy(1)

! update coordinates
               DO iparticle=1,mc_par_1%nunits_tot
                  r_old_1(1:3,iparticle)=&
                     particles_1%els(iparticle)%r(1:3)
               ENDDO
! if we're biasing classically, then we update the classical energy
               IF ( mc_par_1 % lclassical) THEN
                  moves_1%classical_energy=classical_energy_new(1)
               ENDIF
            ENDIF
            IF (box_flag(2) == 1) THEN
               energy_check(2)=energy_check(2)+(new_energy(2)-&
                                              old_energy(2))
               old_energy(2)=new_energy(2)

! update coordinates
               DO iparticle=1,mc_par_2%nunits_tot
                  r_old_2(1:3,iparticle)=&
                     particles_2%els(iparticle)%r(1:3)
               ENDDO
! if we're biasing classically, then we update the classical energy
               IF ( mc_par_1 % lclassical) THEN
                  moves_2%classical_energy=classical_energy_new(2)
               ENDIF
            ENDIF

! back to parallel
         ENDIF

      ELSE

! reject the move
         IF (mc_par_1%ionode) THEN

            CALL move_q_reinit(moves_1)
            CALL move_q_reinit(move_updates_1)
            CALL move_q_reinit(moves_1)
            CALL move_q_reinit(move_updates_1)


         ENDIF ! back to parallel

         IF ( .NOT. mc_par_1%ionode) r_old_1(:,:) = 0.0D0
         IF ( .NOT. mc_par_1%ionode) r_old_2(:,:) = 0.0D0

! coodinates changed, so we need to broadcast those
         CALL mp_bcast(r_old_1,mc_par_1%source,mc_par_1%group)         
         CALL mp_bcast(r_old_2,mc_par_1%source,mc_par_1%group)         
          
         IF (box_flag(1) == 1) THEN 
            DO iparticle=1,mc_par_1%nunits_tot
               particles_1%els(iparticle)%r(1:3)=r_old_1(1:3,iparticle)
            ENDDO
         ENDIF
         IF (box_flag(2) == 1) THEN 
            DO iparticle=1,mc_par_2%nunits_tot
               particles_2%els(iparticle)%r(1:3)=r_old_2(1:3,iparticle)
            ENDDO
         ENDIF

      ENDIF

! make sure the coordinates are transferred
      CALL cp_subsys_set(subsys_1(1)%subsys,particles=particles_1)
      CALL cp_subsys_set(subsys_2(1)%subsys,particles=particles_2)

! nullify some stuff
      NULLIFY(particles_1,particles_2)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_ge_quickstep_move

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
!!     - mc_par_1: the mc parameters for the force env of box 1
!!     - mc_par_2: the mc parameters for the force env of box 2
!!     - force_env_1: the force environment for box 1
!!     - force_env_2: the force environment for box 2
!!     - moves_1: the structure that keeps track of how many moves have been
!!              accepted/rejected for box 1
!!     - moves_2: the structure that keeps track of how many moves have been
!!              accepted/rejected for box 2
!!     - move_updates_1: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated for box 1
!!     - move_updates_2: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated for box 2
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old_1: the coordinates of the last accepted move involving a
!!              Quickstep calculation for box 1
!!     - r_old_2: the coordinates of the last accepted move involving a
!!              Quickstep calculation for box 2
!!     - istep: how many times the energy averages have been taken
!!     - old_energy: the energy of the last accepted move involving a
!!                   Quickstep calculation
!!     - averages_1: the structure that keeps track of the running averages
!!                 in the system for box 1
!!     - averages_1: the structure that keeps track of the running averages
!!                 in the system for box 1
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_ge_swap_move(mc_par,force_env,moves,&
                energy_check,r_old,old_energy,&
                nnstep,globenv,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:), POINTER                  :: mc_par
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env
    TYPE(mc_moves_p_type), DIMENSION(:), &
      POINTER                                :: moves
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: energy_check
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: r_old
    REAL(KIND=dp), DIMENSION(1:2), &
      INTENT(INOUT)                          :: old_energy
    INTEGER, INTENT(IN)                      :: nnstep
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_ge_swap_move"

    INTEGER :: all_atoms, handle, i, iatom, ibox, ii, imolecule, ins_atoms, &
      insert_box, istat, natom, rem_atoms, remove_box
    LOGICAL                                  :: loverlap
    REAL(KIND=dp)                            :: prefactor, rand, w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r_insert, r_insert_mol, &
                                                r_remove
    REAL(KIND=dp), DIMENSION(1:2)            :: classical_energy_new, &
                                                new_energy
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                displace_molecule, mass, &
                                                pos_insert
    TYPE(cell_type), POINTER                 :: cell_insert, cell_remove
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: newsys, oldsys
    TYPE(force_env_type), POINTER            :: insert_env, remove_env
    TYPE(mol_new_list_p_type), &
      DIMENSION(:), POINTER                  :: molecules_new, molecules_old
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_new, particles_old

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! reset the overlap flag
      loverlap=.FALSE.

! nullify some pointers
      NULLIFY(molecules_new)
      NULLIFY(particles_new)
      NULLIFY(newsys)
      NULLIFY(particles_old)
      NULLIFY(molecules_old)
      NULLIFY(oldsys)

! allocate said pointers
      all_atoms=mc_par(1)%mc_par%nunits_tot+mc_par(2)%mc_par%nunits_tot
      ALLOCATE(r_insert_mol(1:3,1:mc_par(1)%mc_par%nunits),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert_mol",3*mc_par(1)%mc_par%nunits)
      ALLOCATE(oldsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys",2)
      ALLOCATE(particles_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old",2)
      ALLOCATE(molecules_old(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecules_old",2) 
      ALLOCATE(newsys(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "newsys",2)
      ALLOCATE(particles_new(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_new",2)      
      ALLOCATE(molecules_new(1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecules_new",2) 

! get the old coordinates
      DO i=1,2
         CALL force_env_get(force_env(i)%force_env,&
            subsys=oldsys(i)%subsys,error=error)
         CALL cp_subsys_get(oldsys(i)%subsys(1)%subsys, &
            particles=particles_old(i)%list, &
            molecules_new=molecules_old(i)%list,error=error)
      ENDDO
      
      IF(mc_par(1)%mc_par%ionode) THEN
!     choose a direction to swap
         CALL RANDOM_NUMBER(rand)
         IF ( rand .LE. 0.50D0 ) THEN
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
         IF(mc_par(remove_box)%mc_par%nchain == 0) THEN
            loverlap=.TRUE.
            moves(insert_box)%moves%empty=&
               moves(insert_box)%moves%empty+1
         ELSE
            CALL RANDOM_NUMBER(rand)    
            imolecule=CEILING(rand*(mc_par(remove_box)%mc_par%nchain))
         ENDIF
         
      ENDIF
 
! exit if the box is empty
      CALL mp_bcast(loverlap,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      IF(loverlap) THEN
! deallocate some stuff
         DEALLOCATE(r_insert_mol,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_insert_mol")
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_old")
         DEALLOCATE(molecules_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecules_old") 
         DEALLOCATE(newsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"newsys")
         DEALLOCATE(particles_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_new")      
         DEALLOCATE(molecules_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecules_new") 
         CALL timestop(0.0d0,handle)
         RETURN
      ENDIF

      CALL mp_bcast(insert_box,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      CALL mp_bcast(remove_box,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      
! allocate a couple of arrays that depend on the box
      ins_atoms=mc_par(insert_box)%mc_par%nunits_tot+&
         mc_par(insert_box)%mc_par%nunits
      ALLOCATE(r_insert(1:3,1:ins_atoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert",3*ins_atoms)
      rem_atoms=mc_par(remove_box)%mc_par%nunits_tot-&
         mc_par(remove_box)%mc_par%nunits
      IF (rem_atoms == 0) THEN
         ALLOCATE(r_remove(1:3,1:mc_par(remove_box)%mc_par%nunits),&
            STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove",3*mc_par(remove_box)%mc_par%nunits)
      ELSE
         ALLOCATE(r_remove(1:3,1:rem_atoms),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove",3*rem_atoms)
      ENDIF
      
      IF(mc_par(1)%mc_par%ionode) THEN
! grab the cells for later...acceptance and insertion
         CALL force_env_get(force_env(insert_box)%force_env,&
            cell=cell_insert,error=error)
         CALL force_env_get(force_env(remove_box)%force_env,&
            cell=cell_remove,error=error)
         
! prepare the insertion box
! choose an insertion point
         DO i=1,3
            CALL RANDOM_NUMBER(rand)
            pos_insert(i)=rand*cell_insert%hmat(1,1)
         ENDDO
      ENDIF

      natom=3
      mass(1:3)=(/ 15.999d0,1.0074d0,1.0074d0 /)
      
      IF(mc_par(1)%mc_par%ionode) THEN
         DO i=molecules_old(remove_box)%list%els(imolecule)%first_atom,&
                     molecules_old(remove_box)%list%els(imolecule)%&
                     first_atom+natom-1
            ii=i-molecules_old(remove_box)%list%els(imolecule)%&
               first_atom+1
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
         IF (mc_par(insert_box)%mc_par%nchain == 0) THEN
            DO iatom=1,natom
               r_insert(1:3,iatom)=&
                  particles_old(insert_box)%list%els(iatom)%r(1:3)
            ENDDO
         ELSE
            DO iatom=1,mc_par(insert_box)%mc_par%nunits_tot
               r_insert(1:3,iatom)=&
                  particles_old(insert_box)%list%els(iatom)%r(1:3)
            ENDDO
            DO iatom=1+mc_par(insert_box)%mc_par%nunits_tot,&
                  natom+mc_par(insert_box)%mc_par%nunits_tot
               ii=iatom-mc_par(insert_box)%mc_par%nunits_tot
               r_insert(1:3,iatom)=r_insert_mol(1:3,ii)
            ENDDO
         ENDIF
      
! fold the coordinates into the box and check for overlaps
         CALL mc_coordinate_fold(r_insert(:,:),natom+&
            mc_par(insert_box)%mc_par%nunits_tot,cell_insert%hmat(1,1))
         loverlap=.FALSE.
         CALL check_for_overlap(r_insert(:,:),&
            mc_par(insert_box)%mc_par%nunits_tot+natom,&
            cell_insert%hmat(1,1),loverlap)
      ENDIF
      
      CALL mp_bcast(loverlap,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      
      IF(loverlap) THEN
! deallocate some stuff
         DEALLOCATE(r_remove,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_remove")
         DEALLOCATE(r_insert,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_insert")
         DEALLOCATE(r_insert_mol,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"r_insert_mol")
         DEALLOCATE(oldsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"oldsys")
         DEALLOCATE(particles_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_old")
         DEALLOCATE(molecules_old,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecules_old") 
         DEALLOCATE(newsys,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"newsys")
         DEALLOCATE(particles_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"particles_new")      
         DEALLOCATE(molecules_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"molecules_new") 
         CALL timestop(0.0d0,handle)
         RETURN
      ENDIF

! make the .dat file
      IF (mc_par(insert_box)%mc_par%ionode) THEN

! if we made it this far, we have no overlaps
         moves(insert_box)%moves%grown=&
            moves(insert_box)%moves%grown+1

         CALL mc_make_dat_file(r_insert(:,:),&
            mc_par(insert_box)%mc_par%nunits_tot+3,&
            cell_insert%hmat(1,1)*angstrom,&
            mc_par(insert_box)%mc_par%dat_file)
         
! now do the same for the removal box...be careful not to make an empty box
         IF (mc_par(remove_box)%mc_par%nchain == 1) THEN
            DO iatom=1,natom
               r_remove(1:3,iatom)=r_insert_mol(1:3,iatom)
            ENDDO
!     make the .dat file
            IF (mc_par(remove_box)%mc_par%ionode) &
            CALL mc_make_dat_file(r_remove(:,:),3,&
               cell_remove%hmat(1,1)*angstrom,&
               mc_par(remove_box)%mc_par%dat_file)
         ELSE
            DO iatom=1,molecules_old(remove_box)%list%els(imolecule)%&
                  first_atom-1
               r_remove(1:3,iatom)=&
                  particles_old(remove_box)%list%els(iatom)%r(1:3)
            ENDDO
            DO iatom=molecules_old(remove_box)%list%els(imolecule)%&
                  first_atom+natom,mc_par(remove_box)%mc_par%nunits_tot
               ii=iatom-natom
               r_remove(1:3,ii)=&
                  particles_old(remove_box)%list%els(iatom)%r(1:3)
            ENDDO

! make the .dat file
            IF (mc_par(remove_box)%mc_par%ionode) &
               CALL mc_make_dat_file(r_remove(:,:),&
                mc_par(remove_box)%mc_par%nunits_tot-3,&
                cell_remove%hmat(1,1)*angstrom,&
                mc_par(remove_box)%mc_par%dat_file)
         ENDIF
      ENDIF

! now let's create the two new environments and calculate the energy...the
! mc_par contain the names of the other box's input file

      CALL mp_bcast(natom,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      
      globenv%input_file_name=mc_par(insert_box)%mc_par%dat_file
      CALL quickstep_create_force_env(insert_env, globenv, error=error)
      globenv%input_file_name=mc_par(remove_box)%mc_par%dat_file
      CALL quickstep_create_force_env(remove_env, globenv, error=error)

! calculate the energies
!!!!!!!!!!!!!!!!!!!!
!      CALL force_env_calculate_energy(insert_env)
!      CALL force_env_get(insert_env,&
!         potential_energy=new_energy(insert_box),error=error)
! now change the total number of units and chains in the mc_par
      mc_par(insert_box)%mc_par%nchain=&
         mc_par(insert_box)%mc_par%nchain+1
      mc_par(insert_box)%mc_par%nunits_tot=&
         mc_par(insert_box)%mc_par%nunits_tot+natom
      CALL mc_classical_energy(insert_env,&
         classical_energy_new(insert_box),&
         mc_par(insert_box)%mc_par,1,&
         mc_par(insert_box)%mc_par%nchain,6.211D0/angstrom)

!!!!!!!!!!!!!!!!!
      new_energy(insert_box)=classical_energy_new(insert_box)

      mc_par(insert_box)%mc_par%nchain=&
         mc_par(insert_box)%mc_par%nchain-1
      mc_par(insert_box)%mc_par%nunits_tot=&
         mc_par(insert_box)%mc_par%nunits_tot-natom
      
      IF (mc_par(remove_box)%mc_par%nchain-1 == 0) THEN
         new_energy(remove_box)=0.0D0
         classical_energy_new(remove_box)=0.0D0
      ELSE
!!!!!!!!!!!!!!!!!
!         CALL force_env_calculate_energy(remove_env)
!         CALL force_env_get(remove_env,&
!            potential_energy=new_energy(remove_box),error=error)
         mc_par(remove_box)%mc_par%nchain=&
            mc_par(remove_box)%mc_par%nchain-1
         mc_par(remove_box)%mc_par%nunits_tot=&
            mc_par(remove_box)%mc_par%nunits_tot-natom
         CALL mc_classical_energy(remove_env,&
            classical_energy_new(remove_box),&
            mc_par(remove_box)%mc_par,1,&
            mc_par(remove_box)%mc_par%nchain,6.211D0/angstrom)
         
!!!!!!!!!!!!!!!!!!!!!
         new_energy(remove_box)=classical_energy_new(remove_box)
         
         mc_par(remove_box)%mc_par%nchain=&
            mc_par(remove_box)%mc_par%nchain+1
         mc_par(remove_box)%mc_par%nunits_tot=&
            mc_par(remove_box)%mc_par%nunits_tot+natom
      ENDIF
      
! figure out the prefactor to the boltzmann weight in the acceptance
! rule, based on numbers of particles and volumes

      IF (mc_par(1)%mc_par%ionode) THEN
         prefactor=REAL(mc_par(remove_box)%mc_par%nchain,dbl)*&
            cell_insert%hmat(1,1)/&
            (REAL(mc_par(insert_box)%mc_par%nchain+1,dbl)&
            *cell_remove%hmat(1,1))

! now to calculate the full weight
         w=prefactor*DEXP(-mc_par(1)%mc_par%BETA*&
            (new_energy(insert_box)-&
            old_energy(insert_box)+new_energy(remove_box)-&
            old_energy(remove_box)))

! check if the move is accepted
            
         IF(w .GT. 1.0D0) THEN
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF
         
!         WRITE(6,*) 'rand,w =',rand,w
!         WRITE(6,*) 'insert =',new_energy(insert_box),&
!            old_energy(insert_box),classical_energy_new(insert_box),&
!            mc_par(insert_box)%mc_par%nchain
!         WRITE(6,*) 'remove =',new_energy(remove_box),&
!            old_energy(remove_box),classical_energy_new(remove_box),&
!            mc_par(remove_box)%mc_par%nchain
      ENDIF

      CALL mp_bcast(rand,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)
      CALL mp_bcast(w,mc_par(1)%mc_par%source,&
         mc_par(1)%mc_par%group)

!!!!!!!!!!!!!!!!!!
!            w=1.0D0
!            w=0.0D0


      IF ( rand .LT. w ) THEN

! accept the move         

         IF ( mc_par(1)%mc_par%ionode) THEN
!            WRITE(6,*) 'ACCEPTING'
!     accept the move
            moves(insert_box)%moves%swap%successes=&
                moves(insert_box)%moves%swap%successes+1

            DO ibox=1,2
! update energies 
               energy_check(ibox)=energy_check(ibox)+(new_energy(ibox)-&
                  old_energy(ibox))
               old_energy(ibox)=new_energy(ibox)

! if we're biasing classically, then we update the classical energy
               IF ( mc_par(ibox)%mc_par % lclassical) &
                  moves(ibox)%moves%classical_energy=&
                  classical_energy_new(ibox)

            ENDDO

! back to parallel
         ENDIF

! change particle numbers
! now change the total number of units and chains in the mc_par
         mc_par(insert_box)%mc_par%nchain=&
            mc_par(insert_box)%mc_par%nchain+1
         mc_par(insert_box)%mc_par%nunits_tot=&
            mc_par(insert_box)%mc_par%nunits_tot+natom
         mc_par(remove_box)%mc_par%nchain=&
            mc_par(remove_box)%mc_par%nchain-1
         mc_par(remove_box)%mc_par%nunits_tot=&
            mc_par(remove_box)%mc_par%nunits_tot-natom

! update coordinates
! insertion box
!         CALL force_env_retain(insert_env)
!         CALL force_env_retain(force_env(insert_box)%force_env)
!         force_env(insert_box)%force_env%ref_count=1
         CALL force_env_release(force_env(insert_box)%force_env)
         CALL cp_subsys_release(oldsys(insert_box)%subsys(1)%subsys)
         force_env(insert_box)%force_env => insert_env

! removal box
!         CALL force_env_retain(remove_env)
!         CALL force_env_retain(force_env(remove_box)%force_env)
!         force_env(remove_box)%force_env%ref_count=1
         CALL force_env_release(force_env(remove_box)%force_env)
         CALL cp_subsys_release(oldsys(remove_box)%subsys(1)%subsys)
         force_env(remove_box)%force_env => remove_env

      ELSE

! reject the move

!         IF (mc_par(1)%mc_par%ionode) THEN
!            WRITE(6,*) 'REJECTING'
!         ENDIF                  ! back to parallel

         CALL force_env_release(insert_env)
         CALL force_env_release(remove_env)
      ENDIF

! deallocate some stuff
      DEALLOCATE(r_remove,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_remove")
      DEALLOCATE(r_insert,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert")
      DEALLOCATE(r_insert_mol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r_insert_mol")
      DEALLOCATE(oldsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "oldsys")
      DEALLOCATE(particles_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_old")
      DEALLOCATE(molecules_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecules_old") 
      DEALLOCATE(newsys,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "newsys")
      DEALLOCATE(particles_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "particles_new")      
      DEALLOCATE(molecules_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "molecules_new") 
      
 
! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_ge_swap_move

! *****************************************************************************

! *****************************************************************************
!!****s* mc_ge_moves.F/mc_ge_volume_move
!!
!!   NAME
!!     mc_volume_move
!!
!!   FUNCTION
!!     performs a Monte Carlo move that alters the volume of the simulation box
!!
!!   ARGUMENTS
!!     - globenv: the global simulation parameters
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of "outer moves" already performed
!!     - counter: the total number of "inner moves" already performed
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move involving a
!!              Quickstep calculation
!!     - old_energy: the energy of the last accepted move involving a
!!                   Quickstep calculation
!!                                 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_ge_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,old_energy,&
                        globenv,&
                        energy_check,r_old,error)

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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_ge_volume_move"
    REAL(KIND=dp), DIMENSION(1:3), PARAMETER :: &
      mass = (/ 15.994d0,1.008d0,1.008d0 /)
    REAL(KIND=dp), PARAMETER                 :: mass_tot = 18.01d0

    INTEGER                                  :: handle, i, iatom, &
                                                ivolmove,j,all_atoms,&
                                                istat
    INTEGER,DIMENSION(1:2)                   :: nunits_tot
    REAL(KIND=dp),DIMENSION(1:2)             :: classical_energy_diff, &
      classical_energy_new, classical_energy_original, &
      classical_energy_test, new_cell_length, new_energy, &
      old_cell_length, test_cell_length, volume_new,&
      volume_old
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                center_of_mass_new, diff
    TYPE(force_env_p_type), DIMENSION(:),POINTER &
                                             :: test_env
    TYPE(cp_subsystem_pp_type), &
      DIMENSION(:), POINTER                  :: oldsys, testsys
    TYPE(cell_p_type), DIMENSION(:),POINTER  &
                                             :: cell_old
    REAL(KIND=dp), DIMENSION(:,:,:),ALLOCATABLE &
                                             :: r
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles_test, particles_old
    REAL(KIND=dp)                            :: rand,w,prefactor,vol_dis
    LOGICAL,DIMENSION(1:2)                   :: lempty
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles_old,particles_test,cell_old,oldsys,&
         testsys,test_env)

! allocate some stuff
      all_atoms=mc_par(1)%mc_par%nunits_tot+mc_par(2)%mc_par%nunits_tot
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

! check for empty boxes...need to be careful because we can't build
! a force_env with no particles
      lempty(:)=.FALSE.
      DO i=1,2
         IF(mc_par(i)%mc_par%nchain==0) THEN
            lempty(i)=.TRUE.
            nunits_tot(i)=mc_par(i)%mc_par%nunits
         ELSE
            nunits_tot(i)=mc_par(i)%mc_par%nunits_tot
         ENDIF
      ENDDO
      

      IF (mc_par(1)%mc_par%ionode) THEN
! record the attempt
         DO i=1,2
            moves(i)%moves%volume%attempts= &
               moves(i)%moves%volume%attempts+1
            move_updates(i)%moves%volume%attempts=&
               move_updates(i)%moves%volume%attempts+1
         ENDDO
      ENDIF

! now let's grab the cell length and particle positions
      DO i=1,2
         CALL force_env_get(force_env(i)%force_env,&
            subsys=oldsys(i)%subsys,cell=cell_old(i)%cell,error=error)
         CALL cp_subsys_get(oldsys(i)%subsys(1)%subsys, &
            particles=particles_old(i)%list, error=error)

! find the old cell length
         old_cell_length(i)=cell_old(i)%cell%hmat(1,1)
         test_cell_length(i)=cell_old(i)%cell%hmat(1,1)

      ENDDO

      IF (mc_par(1)%mc_par%ionode) THEN

         DO i=1,2
! find the old classical energy
            IF(lempty(i)) THEN

               classical_energy_test(i)=0.0d0
               classical_energy_original(i)=0.0d0

            ELSE

              CALL mc_classical_energy(force_env(i)%force_env,&
                  classical_energy_original(i),mc_par(1)%mc_par,1,&
                  mc_par(1)%mc_par%nchain,6.211d0/angstrom,&
                  cell_length=old_cell_length(i))

               classical_energy_test(i)=classical_energy_original(i)
            ENDIF
! save the old coordiantes
            DO iatom=1,nunits_tot(i)
               r(1:3,iatom,i)=particles_old(i)%list%els(iatom)%r(1:3)
            ENDDO
!     fold the coordinates into the central box
            CALL mc_coordinate_fold(r(:,:,i),&
               nunits_tot(i),test_cell_length(i))
! make the new dat file
            CALL mc_make_dat_file(r(:,:,i),nunits_tot(i),&
               test_cell_length(i)*angstrom,mc_par(i)%mc_par%dat_file)


         ENDDO
      ENDIF

      DO i=1,2
! make a test environment that we can alter to calculate the energy
         globenv%input_file_name=mc_par(i)%mc_par%dat_file
         CALL quickstep_create_force_env(test_env(i)%force_env, &
            globenv, error=error)

! grab the coords so we can change them
         CALL force_env_get(test_env(i)%force_env,&
            subsys=testsys(i)%subsys,error=error)
         CALL cp_subsys_get(testsys(i)%subsys(1)%subsys, &
            particles=particles_test(i)%list, error=error)

      ENDDO

! now do some classical moves to find a good move
      IF (mc_par(1)%mc_par%ionode) THEN

! now do a series of volume moves, accepting or rejecting each one based on
!    the fake potential
        DO ivolmove=1,mc_par(1)%mc_par%nvolmoves

! call a random number to figure out how far we're moving
            CALL RANDOM_NUMBER(rand)
!!!!!
            vol_dis=mc_par(1)%mc_par%rmvolume*(rand-0.5d0)*2.0d0
!            vol_dis=5.0d0/angstrom**3
!!!!!
! add to one box, subtract from the other
            new_cell_length(1)=(test_cell_length(1)**3+&
               vol_dis)**(1.0d0/3.0d0)
            new_cell_length(2)=(test_cell_length(2)**3-&
                vol_dis)**(1.0d0/3.0d0)

            DO i=1,2
! save the coords
               DO iatom=1,nunits_tot(i)
                  r(1:3,iatom,i)=&
                     particles_test(i)%list%els(iatom)%r(1:3)
               ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass
               DO iatom=1,nunits_tot(i),mc_par(i)%mc_par%nunits

! now find the center of mass
                  CALL get_center_of_mass(r(:,iatom:iatom+2,i),&
                     mc_par(i)%mc_par%nunits,center_of_mass(:),mass(:))

! scale the center of mass and determine the vector that points from the
!    old COM to the new one
                  center_of_mass_new(1:3)=center_of_mass(1:3)*&
                     new_cell_length(i)/test_cell_length(i)
                  DO j=1,3
                     diff(i)=center_of_mass_new(j)-center_of_mass(j)
! now change the particle positions
                     particles_test(i)%list%els(iatom)%r(j)=&
                        particles_test(i)%list%els(iatom)%r(j)+diff(j)
                     particles_test(i)%list%els(iatom+1)%r(j)=&
                        particles_test(i)%list%els(iatom+1)%r(j)+diff(j)
                     particles_test(i)%list%els(iatom+2)%r(j)=&
                        particles_test(i)%list%els(iatom+2)%r(j)+diff(j)
                  ENDDO
                  
               ENDDO
        
! find and store the classical energy
               IF (lempty(i)) THEN
                  classical_energy_new(i)=0.0D0
               ELSE
                  CALL mc_classical_energy(test_env(i)%force_env,&
                     classical_energy_new(i),mc_par(i)%mc_par,1,&
                     mc_par(i)%mc_par%nchain,6.211d0/angstrom,&
                     cell_length=new_cell_length(i))
               ENDIF
               
            ENDDO

            IF (mc_par(1)%mc_par%nvolmoves .GT. 1) THEN

! determine the acceptance of the move
               volume_new(1:2)=new_cell_length(1:2)**3
               volume_old(1:2)=test_cell_length(1:2)**3
               prefactor=(volume_new(1)**mc_par(1)%mc_par%nchain*&
                  volume_new(2)**mc_par(2)%mc_par%nchain)/&
                  (volume_old(1)**mc_par(1)%mc_par%nchain*&
                  volume_old(2)**mc_par(2)%mc_par%nchain)
               w=prefactor*dexp(-mc_par(1)%mc_par%BETA*&
                  (classical_energy_new(1)+classical_energy_new(2)-&
                  classical_energy_test(1)-classical_energy_test(2)))
               
               IF ( w .GE. 1.0d0 ) THEN
                  w=1.0d0
                  rand=0.0D0
               ELSE
                  CALL RANDOM_NUMBER(rand)
               ENDIF
               
               IF (rand .LT. w ) THEN

                  DO i=1,2
                     move_updates(i)%moves%classical_volume%attempts=&
                        move_updates(i)%moves%classical_volume%attempts+1
                     move_updates(i)%moves%classical_volume%successes=&
                        move_updates(i)%moves%classical_volume%successes+1
                     moves(i)%moves%classical_volume%attempts=&
                        moves(i)%moves%classical_volume%attempts+1
                     moves(i)%moves%classical_volume%successes=&
                        moves(i)%moves%classical_volume%successes+1
                     
! update energy and cell length
                     test_cell_length(i)=new_cell_length(i)
                     classical_energy_test(i)=classical_energy_new(i)
                  ENDDO

               ELSE

                  DO i=1,2
                     move_updates(i)%moves%classical_volume%attempts=&
                       move_updates(i)%moves%classical_volume%attempts+1
                     moves(i)%moves%classical_volume%attempts=&
                       moves(i)%moves%classical_volume%attempts+1

! revert back to old coords                      
                     DO iatom=1,nunits_tot(i)
                        particles_test(i)%list%els(iatom)%r(1:3)=&
                           r(1:3,iatom,i)
                     ENDDO
                  ENDDO
               ENDIF
            
            ELSE
               DO i=1,2
                  test_cell_length(i)=new_cell_length(i)
                  classical_energy_test(i)=classical_energy_new(i)
               ENDDO
            ENDIF
            
         ENDDO

! determine the overall energy difference
         IF (mc_par(1)%mc_par%nvolmoves .GT. 1 ) THEN
            DO i=1,2
               classical_energy_diff(i)=classical_energy_test(i)-&
                  classical_energy_original(i)
            ENDDO
         ELSE
            classical_energy_diff(:)=0.0d0
         ENDIF

         DO i=1,2
            new_energy(i)=classical_energy_test(i)

! make the new dat file
            DO iatom=1,nunits_tot(i)
               r(1:3,iatom,i)=particles_test(i)%list%els(iatom)%r(1:3)
            ENDDO

            CALL mc_make_dat_file(r(:,:,i),&
               nunits_tot(i),test_cell_length(i)*angstrom,&
               mc_par(i)%mc_par%dat_file)
            
         ENDDO
! back to parallel
      ENDIF

      DO i=1,2
! remake the force environment and calculate the energy

            CALL force_env_release(test_env(i)%force_env,error=error)
            CALL cp_subsys_release(testsys(i)%subsys(1)%subsys)
            
            globenv%input_file_name=mc_par(i)%mc_par%dat_file
            CALL quickstep_create_force_env(test_env(i)%force_env, &
               globenv, error=error)

         IF (lempty(i)) THEN
            new_energy(i)=0.0d0
         ELSE
!           CALL force_env_calculate_energy(test_env(i)%force_env)

! now grab the energy
!            CALL force_env_get(test_env(i)%force_env,&
!               potential_energy=new_energy(i),error=error)
         ENDIF
      ENDDO

      IF (mc_par(1)%mc_par%ionode) THEN
! accept or reject the move
         IF (mc_par(1)%mc_par%nvolmoves .GT. 1) THEN
            w=dexp(-mc_par(1)%mc_par%BETA*((new_energy(1)+&
               new_energy(2)-old_energy(1)-old_energy(2))-&
               classical_energy_diff(1)-classical_energy_diff(2)))
         ELSE

            volume_new(1:2)=new_cell_length(1:2)**3
            volume_old(1:2)=test_cell_length(1:2)**3
            prefactor=(volume_new(1)**mc_par(1)%mc_par%nchain*&
               volume_new(2)**mc_par(2)%mc_par%nchain)/&
               (volume_old(1)**mc_par(1)%mc_par%nchain*&
               volume_old(2)**mc_par(2)%mc_par%nchain)
            w=prefactor*dexp(-mc_par(1)%mc_par%BETA*&
               (new_energy(1)+new_energy(2)-&
               old_energy(1)-old_energy(2)))

         ENDIF

!!!!!!!!!!!
!         w=1.0D0
!         w=0.0D0
!!!!!!!!!!!
         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

      ENDIF


      CALL mp_bcast(rand,mc_par(1)%mc_par%source,mc_par(1)%mc_par%group)
      CALL mp_bcast(w,mc_par(1)%mc_par%source,mc_par(1)%mc_par%group)

      IF (rand .LT. w) THEN

! write cell length, volume, density, and trial displacement to a file
         IF(mc_par(1)%mc_par%ionode) THEN

            WRITE(mc_par(1)%mc_par%cl,*) nnstep,test_cell_length(1)*&
               angstrom,vol_dis*(angstrom)**3,test_cell_length(2)*&
               angstrom,vol_dis*(angstrom)**3
            WRITE(mc_par(1)%mc_par%cl,*) nnstep,new_energy(1),&
               old_energy(1),new_energy(2),old_energy(2),&
               classical_energy_diff(1:2)

            DO i=1,2
! accept the move
               moves(i)%moves%volume%successes=&
                  moves(i)%moves%volume%successes+1
               move_updates(i)%moves%volume%successes=&
                  move_updates(i)%moves%volume%successes+1

! update energies
               energy_check(i)=energy_check(i)+(new_energy(i)-&
                  old_energy(i))
               old_energy(i)=new_energy(i)

! and classical energy
               moves(i)%moves%classical_energy=&
                  classical_energy_test(i)

! and the coordiantes
               CALL force_env_get(test_env(i)%force_env,subsys=&
                  testsys(i)%subsys,error=error)
               CALL cp_subsys_get(testsys(i)%subsys(1)%subsys, &
                  particles=particles_test(i)%list, error=error)
               
               DO iatom=1,nunits_tot(i)
                  r_old(1:3,iatom,i)=&
                     particles_test(i)%list%els(iatom)%r(1:3)
               ENDDO
               
            ENDDO
         ENDIF

! switch the force env to the new env
         DO i=1,2
            CALL force_env_release(force_env(i)%force_env,error=error)
            CALL cp_subsys_release(oldsys(i)%subsys(1)%subsys)
            force_env(i)%force_env => test_env(i)%force_env
         ENDDO
      ELSE

! reject the move
! write cell length, volume, density, and trial displacement to a file
         IF(mc_par(1)%mc_par%ionode) THEN

            WRITE(mc_par(1)%mc_par%cl,*) nnstep,test_cell_length(1)*&
               angstrom,vol_dis*(angstrom)**3,test_cell_length(2)*&
               angstrom,vol_dis*(angstrom)**3
            WRITE(mc_par(1)%mc_par%cl,*) nnstep,new_energy(1),&
               old_energy(1),new_energy(2),old_energy(2),&
               classical_energy_diff(1:2)

         ENDIF

! release some memory
         DO i=1,2
            CALL force_env_release(test_env(i)%force_env,error=error)
         ENDDO

         DEALLOCATE(test_env,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,&
            __LINE__,"test_env")
      ENDIF

! free up some memory
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

! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_ge_volume_move


!******************************************************************************
!!****s* mc_misc.F/get_center_of_mass
!!
!!   NAME
!!     get_center_of_mass
!!
!!   SYNOPSIS
!!     Subroutine get_center_of_mass (coordinates,natom,center_of_mass,mass)
!!           Real(Kind=dp), Dimension(1:3),Intent(OUT)
!!                                                 :: center_of_mass
!!           Real(Kind=dp), Dimension(:,:),Intent(IN)
!!                                                 :: coordinates
!!           Real(Kind=dp), Dimension(:),Intent(IN)
!!                                                 :: mass
!!           Integer, Intent(IN)                   :: natom
!!     End Subroutine get_center_of_mass
!!
!!   FUNCTION
!!     calculates the center of mass of a given molecule
!!
!!   ARGUMENTS
!!     - natom: the total number of atoms
!!     - coordinates: the coordiantes of the atoms in the molecule
!!     - mass: the mass of the atoms in the molecule
!!     - center_of_mass: the coordinates of the center of mass
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE check_for_overlap(coordinates,natoms,box_length,loverlap)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coordinates
    INTEGER, INTENT(IN)                      :: natoms
    REAL(KIND=dp), INTENT(IN)                :: box_length
    LOGICAL, INTENT(OUT)                     :: loverlap

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "check_for_overlap"

    INTEGER                                  :: handle, i, iatom, imol, j, &
                                                jmol, nmol
    REAL(KIND=dp)                            :: dist
    REAL(KIND=dp), DIMENSION(1:3)            :: RIJ
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:3, natoms/3)          :: r

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

      DO iatom=1,natoms,3
         r(1:3,1,(iatom+2)/3)=coordinates(1:3,iatom)
         r(1:3,2,(iatom+2)/3)=coordinates(1:3,iatom+1)
         r(1:3,3,(iatom+2)/3)=coordinates(1:3,iatom+2)
      ENDDO

      nmol=natoms/3

      DO imol=1,nmol-1
         DO jmol=imol+1,nmol
            
! check the O-O distance
            RIJ(1)=r(1,1,imol)-r(1,1,jmol)-box_length*ANINT(&
               (r(1,1,imol)-r(1,1,jmol))/box_length)
            RIJ(2)=r(2,1,imol)-r(2,1,jmol)-box_length*ANINT(&
               (r(2,1,imol)-r(2,1,jmol))/box_length)
            RIJ(3)=r(3,1,imol)-r(3,1,jmol)-box_length*ANINT(&
               (r(3,1,imol)-r(3,1,jmol))/box_length)
            
            dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
            
            IF(dist .LT. 2.2D0/angstrom) THEN
               loverlap=.TRUE.
               CALL timestop(0.0d0,handle)
               RETURN
            ELSEIF(dist .LT. 2.6D0/angstrom) THEN
!     check the H-H distance
               DO i=2,3
                  DO j=2,3
                        RIJ(1)=r(1,i,imol)-r(1,j,jmol)-box_length&
                           *ANINT((r(1,i,imol)-r(1,j,jmol))/box_length)
                        RIJ(2)=r(2,i,imol)-r(2,j,jmol)-box_length&
                           *ANINT((r(2,i,imol)-r(2,j,jmol))/box_length)
                        RIJ(3)=r(3,i,imol)-r(3,j,jmol)-box_length&
                           *ANINT((r(3,i,imol)-r(3,j,jmol))/box_length)
                        dist=DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        IF (dist .LT. 1.3D0/angstrom) THEN
                           loverlap=.TRUE.
                           CALL timestop(0.0d0,handle)
                           RETURN
                        ENDIF
                  ENDDO
               ENDDO
            ENDIF

         ENDDO
      ENDDO

! end the timing
      CALL timestop(0.0d0,handle)


END SUBROUTINE check_for_overlap

END MODULE mc_ge_moves

!**********************************************************************
