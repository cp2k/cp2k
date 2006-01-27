!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_moves [1.0] *
!!
!!   NAME
!!     mc_moves
!!
!!   FUNCTION
!!     the various moves in Monte Carlo (MC) simulations, including
!!     change of internal conformation, translation of a molecule, rotation
!!     of a molecule, and changing the size of the simulation box
!!
!!   AUTHOR
!!     Matthew J. McGrath  (10.16.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_moves
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type,&
                                             use_fist_force
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_control,                      ONLY: mc_create_force_env
  USE mc_coordinates,                  ONLY: get_center_of_mass,&
                                             check_for_overlap,&
                                             mc_coordinate_fold,&
                                             generate_cbmc_swap_config,&
                                             create_discrete_array
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_move_control,                 ONLY: move_q_reinit,&
                                             q_move_accept
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             bond_type,&
                                             bend_type,&
                                             torsion_type,&
                                             atom_type
  USE parallel_rng_types,              ONLY: rng_stream_type,&
                                             next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE cell_types,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE  :: change_bond_angle,change_bond_length,depth_first_search,&
      change_dihedral

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_moves"

  PUBLIC :: mc_atom_translation,mc_molecule_translation,&
            mc_molecule_rotation,mc_volume_move,mc_avbmc_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_moves.F/depth_first_search
!!
!!   NAME
!!     depth_first_search
!!
!!   FUNCTION
!!     essentially performs a depth-first search of the molecule structure
!!     to find all atoms connected to a specific atom excluding one branch...
!!     for instance, if water is labelled 1-2-3 for O-H-H, calling this
!!     routine with current_atom=1,avoid_atom=2 returns the array
!!     atom=(0,0,1)
!!
!!   ARGUMENTS
!!     - current_atom: the atom whose connections we're looking at
!!     - avoid_atom: the atom whose direction the search is not supposed to go
!!     - connectivity: an array telling us the neigbors of all atoms
!!     - atom: the array that tells us if one can get to a given atom by
!!       starting at current_atom and not going through avoid_atom...0 is no,
!!       1 is yes
!!              
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  RECURSIVE SUBROUTINE depth_first_search ( current_atom,avoid_atom,&
      connectivity,atom) 

!     ***************************************************************************
    INTEGER, INTENT(IN)                      :: current_atom,avoid_atom
    INTEGER,DIMENSION(:,:),INTENT(IN)        :: connectivity
    INTEGER,DIMENSION(:),INTENT(INOUT)        :: atom

    INTEGER                                  :: iatom


!   ---------------------------------------------------------------------------

      DO iatom=1,6
         IF(connectivity(iatom,current_atom) .NE. 0) THEN
            IF(connectivity(iatom,current_atom) .NE. avoid_atom) THEN
               atom(connectivity(iatom,current_atom))=1
               CALL depth_first_search ( connectivity(iatom,current_atom),&
                 current_atom,connectivity,atom)
            ENDIF
         ELSE
            RETURN
         ENDIF
      ENDDO

  END SUBROUTINE depth_first_search

! *****************************************************************************
!!****s* mc_moves.F/mc_atom_translation
!!
!!   NAME
!!     mc_atom_translation
!!
!!   FUNCTION
!!     performs either a bond or angle change move for a given molecule
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - bias_env: the force environment used to bias the move, if any (it may
!!           be null if lbias=.false. in mc_par)
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - mol: the number of the molecule's first atom, assuming the rest of
!!       the atoms follow sequentially
!!     - bias_energy: the biased energy of the system before the move
!!     - move_type: dictates what kind of conformational change we do
!!     - rng_stream: the random number stream that we draw from
!!     - lreject: set to .true. if there is an overlap
!!     - error: the cp_error type                            
!!
!!     This subroutine is designed to run in parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_atom_translation ( mc_par,force_env,bias_env, moves,&
                        move_updates,mol,bias_energy,&
                        move_type,lreject,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    LOGICAL, INTENT(OUT)                     :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_atom_translation"

    INTEGER                                  :: handle, imol, ipart,&
                                                istat,nchain,nunits_tot,&
                                                nunits,group,source
    LOGICAL                                  :: lbias,&
      loverlap,ionode
    REAL(KIND=dp)                            :: bias_energy_new, &
                                                bias_energy_old, &
                                                dis_length, rand, w, BETA
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_old, r_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)
      
! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,&
         nunits=nunits,group=group,source=source,ionode=ionode)

! nullify some pointers
      NULLIFY(particles,subsys,molecule_kinds_new,molecule_kind)

! do some allocation
      ALLOCATE (r_old(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)
      ALLOCATE (r_new(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new",3*nunits_tot*dp_size)

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)

! save the coordinates
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO

! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles,  molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)

! save the coordinates
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO

      ENDIF

! do the move
      IF (move_type == 'bond') THEN

! record the attempt
         moves%bond%attempts=moves%bond%attempts+1
         move_updates%bond%attempts=move_updates%bond%attempts+1
         moves%bias_bond%attempts=moves%bias_bond%attempts+1
         move_updates%bias_bond%attempts=move_updates%bias_bond%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%bond%qsuccesses=moves%bond%qsuccesses+1
            move_updates%bond%qsuccesses=&
                 move_updates%bond%qsuccesses+1
            moves%bias_bond%qsuccesses=moves%bias_bond%qsuccesses+1
            move_updates%bias_bond%qsuccesses=&
                 move_updates%bias_bond%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_length(r_old,r_new,&
            mc_par,mol,molecule_kind,dis_length,particles,rng_stream)

      ELSEIF( move_type == 'angle') THEN

! record the attempt
         moves%angle%attempts=moves%angle%attempts+1
         move_updates%angle%attempts=move_updates%angle%attempts+1
         moves%bias_angle%attempts=moves%bias_angle%attempts+1
         move_updates%bias_angle%attempts=move_updates%bias_angle%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%angle%qsuccesses=moves%angle%qsuccesses+1
            move_updates%angle%qsuccesses=&
                 move_updates%angle%qsuccesses+1
            moves%bias_angle%qsuccesses=moves%bias_angle%qsuccesses+1
            move_updates%bias_angle%qsuccesses=&
                 move_updates%bias_angle%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_angle(r_old,r_new,&
            mc_par,mol,molecule_kind,particles,rng_stream)
!         CALL change_bond_angle(r_old,r_new,mc_par,mol,rng_stream)
         dis_length=1.0E0_dp
      ELSE
! record the attempt
         moves%dihedral%attempts=moves%dihedral%attempts+1
         move_updates%dihedral%attempts=move_updates%dihedral%attempts+1
         moves%bias_dihedral%attempts=moves%bias_dihedral%attempts+1
         move_updates%bias_dihedral%attempts=move_updates%bias_dihedral%attempts+1
         IF ( .NOT. lbias ) THEN
            moves%dihedral%qsuccesses=moves%dihedral%qsuccesses+1
            move_updates%dihedral%qsuccesses=&
                 move_updates%dihedral%qsuccesses+1
            moves%bias_dihedral%qsuccesses=moves%bias_dihedral%qsuccesses+1
            move_updates%bias_dihedral%qsuccesses=&
                 move_updates%bias_dihedral%qsuccesses+1
         ENDIF

! do the move
         CALL change_dihedral(r_old,r_new,&
            mc_par,mol,molecule_kind,particles,rng_stream)
         dis_length=1.0E0_dp

      ENDIF

! set the coordinates
      DO ipart=1,nunits_tot
         particles%els(ipart)%r(1:3)=r_new(1:3,ipart)
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)


      imol = (mol+nunits-1)/nunits
! check for overlap
      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,&
            nchain,imol,imol,loverlap)
      ELSE
         CALL check_for_overlap(force_env,&
            nchain,imol,imol,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
               potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule with a 
! correction factor for the change in phase space...dis_length is
! made unitless in change_bond_length
            w=exp(-BETA*(bias_energy_new-&
               bias_energy_old))*dis_length**2
         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) THEN
               rand=next_random_number(rng_stream)
!            CALL RANDOM_NUMBER(rand)
            ENDIF
            CALL mp_bcast(rand,source,group)            
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            IF (move_type == 'bond') THEN
               moves%bond%qsuccesses=moves%bond%qsuccesses+1
               move_updates%bond%successes=&
                  move_updates%bond%successes+1
               moves%bias_bond%successes=moves%bias_bond%successes+1
               move_updates%bias_bond%successes=&
                  move_updates%bias_bond%successes+1
            ELSEIF(move_type == 'angle') THEN
               moves%angle%qsuccesses=moves%angle%qsuccesses+1
               move_updates%angle%successes=&
                  move_updates%angle%successes+1
               moves%bias_angle%successes=moves%bias_angle%successes+1
               move_updates%bias_angle%successes=&
                  move_updates%bias_angle%successes+1
            ELSE
               moves%dihedral%qsuccesses=moves%dihedral%qsuccesses+1
               move_updates%dihedral%successes=&
                  move_updates%dihedral%successes+1
               moves%bias_dihedral%successes=moves%bias_dihedral%successes+1
               move_updates%bias_dihedral%successes=&
                  move_updates%bias_dihedral%successes+1
            ENDIF

            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
                particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF


! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")
      DEALLOCATE(r_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new")

! end the timing
  CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_atom_translation

! *****************************************************************************
!!****s* mc_moves.F/mc_molecule_translation
!!
!!   NAME
!!     mc_molecule_translation
!!
!!   FUNCTION
!!     translates the given molecule randomly in either the x,y, or z direction
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - bias_env: the force environment used to bias the move, if any (it may
!!           be null if lbias=.false. in mc_par)
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - mol: the number of the molecule's first atom, assuming the rest of
!!       the atoms follow sequentially
!!     - bias_energy: the biased energy of the system before the move
!!     - lreject: set to .true. if there is an overlap
!!     - rng_stream: the random number stream that we draw from
!!     - error: the cp_error type                            
!!
!!     This subroutine designed to be parallel.              
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

 SUBROUTINE mc_molecule_translation(  mc_par,force_env, bias_env,moves,&
                      move_updates,mol,bias_energy,&
                        lreject,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    LOGICAL,INTENT(OUT)                      :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_translation"

    INTEGER                                  :: handle, i, imol, ipart, &
                                                iparticle,istat,nchain,&
                                                nunits_tot,nunits,&
                                                group,source
    REAL(KIND=dp)                            :: bias_energy_new, &
                                                bias_energy_old, &
                                                dis_mol, rand, w, BETA,&
                                          rmtrans
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_old
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: lbias,loverlap,ionode

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,rmtrans=rmtrans,&
         nunits=nunits,ionode=ionode,source=source,&
         group=group)

! nullify some pointers
      NULLIFY(particles,subsys)

! do some allocation
      ALLOCATE (r_old(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, error=error)

! save the coordinates
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
        ENDDO

! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, error=error)
      ENDIF

! record the attempt
      moves%trans%attempts=moves%trans%attempts+1
      move_updates%trans%attempts=move_updates%trans%attempts+1
      moves%bias_trans%attempts=moves%bias_trans%attempts+1
      move_updates%bias_trans%attempts=move_updates%bias_trans%attempts+1
      IF ( .NOT. lbias ) THEN
         moves%trans%qsuccesses=moves%trans%qsuccesses+1
         move_updates%trans%qsuccesses=move_updates%trans%qsuccesses+1
         moves%bias_trans%qsuccesses=moves%bias_trans%qsuccesses+1
         move_updates%bias_trans%qsuccesses=move_updates%bias_trans%qsuccesses+1
      ENDIF

! move one molecule in the system

! call a random number to figure out which direction we're moving
      IF(ionode) rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)
      i=AINT(3*rand)+1

! call a random number to figure out how far we're moving
!      CALL RANDOM_NUMBER(rand)
      IF(ionode) rand=next_random_number(rng_stream)
      CALL mp_bcast(rand,source,group)
      dis_mol=rmtrans*(rand-0.5E0_dp)*2.0E0_dp

! do the move
      DO iparticle=mol,mol+nunits-1
         particles%els(iparticle)%r(i)=&
             particles%els(iparticle)%r(i)+dis_mol
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! figure out if there is any overlap...need the number of the molecule
      imol = (mol+nunits-1)/nunits

      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,&
            nchain,imol,imol,loverlap)
      ELSE
         CALL check_for_overlap(force_env,&
            nchain,imol,imol,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing with a cheaper potential, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 
         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
               potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule
            w=exp(-BETA*(bias_energy_new-&
               bias_energy_old))
         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) rand=next_random_number(rng_stream)
!            CALL RANDOM_NUMBER(rand)
            CALL mp_bcast(rand,source,group)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%bias_trans%successes=moves%bias_trans%successes+1
            move_updates%bias_trans%successes=move_updates%bias_trans%successes+1
            moves%trans%qsuccesses=moves%trans%qsuccesses+1
            move_updates%trans%successes=&
                      move_updates%trans%successes+1
            moves%qtrans_dis=moves%qtrans_dis+ABS(dis_mol)
            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE
            
! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
               particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF


! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")

! end the timing
  CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_molecule_translation

! *****************************************************************************
!!****s* mc_moves.F/mc_molecule_rotation
!!
!!   NAME
!!     mc_molecule_rotation
!!
!!   FUNCTION
!!     rotates the given molecule randomly around the x,y, or z axis...
!!     only works for water at the moment
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - bias_env: the force environment used to bias the move, if any (it may
!!           be null if lbias=.false. in mc_par)
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - mol: the number of the molecule's first atom, assuming the rest of
!!       the atoms follow sequentially
!!     - bias_energy: the biased energy of the system before the move
!!     - lreject: set to .true. if there is an overlap
!!     - rng_stream: the random number stream that we draw from
!!     - error: the cp_error type 
!!                                 
!!     This subroutine has been designed for parallel use.
!!     
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, bias_env,moves,&
                      move_updates,mol,bias_energy,&
                        lreject,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    LOGICAL,INTENT(OUT)                      :: lreject
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
 
    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_rotation"

    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass 
    INTEGER                                  :: dir, handle, ii, imol, ipart, &
                                                iunit,istat,nchain,nunits_tot,&
                                                nunits,iatom,natoms,&
                                                group,source
    LOGICAL                                  :: lx, ly,lbias,loverlap,ionode
    REAL(KIND=dp) :: bias_energy_new, bias_energy_old, cosdg, &
      dgamma, nxcm, nycm, nzcm, rand, rx, rxnew, ry, rynew, rz, rznew, sindg, &
      w,BETA,rmrot,masstot
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_old
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,nunits=nunits,rmrot=rmrot,&
         ionode=ionode,group=group,source=source)

! nullify some pointers
      NULLIFY(particles,subsys)

! do some allocation
      ALLOCATE (r_old(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)

! initialize some stuff
      lx =.FALSE.
      ly =.FALSE.

! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, molecule_kinds_new=molecule_kinds_new,&
            error=error)

! save the coordinates
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO

! save the energy
         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            molecule_kinds_new=molecule_kinds_new,&
            particles=particles, error=error)
      ENDIF

! grab the masses
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,natom=natoms)   
      ALLOCATE(mass(1:natoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass",natoms*dp_size)
      DO iatom=1,natoms
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=mass(iatom))
      ENDDO
      masstot=SUM(mass(:))

! record the attempt
      moves%bias_rot%attempts=moves%bias_rot%attempts+1
      move_updates%bias_rot%attempts=move_updates%bias_rot%attempts+1
      moves%rot%attempts=moves%rot%attempts+1
      move_updates%rot%attempts=move_updates%rot%attempts+1
      IF ( .NOT. lbias ) THEN
         moves%rot%qsuccesses=moves%rot%qsuccesses+1
         move_updates%rot%qsuccesses=move_updates%rot%qsuccesses+1
         moves%bias_rot%qsuccesses=moves%bias_rot%qsuccesses+1
         move_updates%bias_rot%qsuccesses=move_updates%bias_rot%qsuccesses+1
      ENDIF

! rotate one molecule in the system

! call a random number to figure out which direction we're moving
      IF(ionode) rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)
      dir=AINT(3*rand)+1

      IF (dir .EQ. 1) THEN
         lx = .TRUE.
      ELSEIF (dir .EQ. 2) THEN
         ly = .TRUE.
      ENDIF

! Determine new center of mass for chain i by finding the sum
! of m*r for each unit, then dividing by the total mass of the chain
      nxcm = 0.0E0_dp
      nycm = 0.0E0_dp
      nzcm = 0.0E0_dp
      DO ii = 1, nunits
         nxcm = nxcm + particles%els(mol-1+ii)%r(1)* mass(ii)
         nycm = nycm + particles%els(mol-1+ii)%r(2)* mass(ii)
         nzcm = nzcm + particles%els(mol-1+ii)%r(3)* mass(ii)
      ENDDO
      nxcm = nxcm / masstot
      nycm = nycm / masstot
      nzcm = nzcm / masstot

! call a random number to figure out how far we're moving
      IF(ionode) rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)
      dgamma=rmrot*(rand-0.5E0_dp)*2.0E0_dp

! *** set up the rotation marix ***
 
      cosdg = COS( dgamma )
      sindg = SIN( dgamma )
 
      IF (lx) THEN
 
! ***    ROTATE UNITS OF I AROUND X-AXIS ***
 
         DO  iunit = mol,mol+nunits-1
            ry = particles%els(iunit)%r(2) - nycm
            rz = particles%els(iunit)%r(3) - nzcm
            rynew = cosdg * ry - sindg * rz
            rznew = cosdg * rz + sindg * ry
               
            particles%els(iunit)%r(2) = rynew + nycm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO
      ELSEIF (ly) THEN 
 
! ***    ROTATE UNITS OF I AROUND y-AXIS ***
 
         DO  iunit = mol,mol+nunits-1
            rx = particles%els(iunit)%r(1) - nxcm
            rz = particles%els(iunit)%r(3) - nzcm
            rxnew = cosdg * rx + sindg * rz
            rznew = cosdg * rz - sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO

      ELSE
 
! ***    ROTATE UNITS OF I AROUND z-AXIS ***
 
         DO  iunit = mol,mol+nunits-1
            rx = particles%els(iunit)%r(1) - nxcm
            ry = particles%els(iunit)%r(2) - nycm

            rxnew = cosdg * rx - sindg * ry
            rynew = cosdg * ry + sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(2) = rynew + nycm

         ENDDO

      ENDIF 
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! check for overlap
      imol = (mol+nunits-1)/nunits
      lreject=.FALSE.
      IF(lbias) THEN
         CALL check_for_overlap(bias_env,&
            nchain,imol,imol,loverlap)
      ELSE
         CALL check_for_overlap(force_env,&
            nchain,imol,imol,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.,error=error)
            CALL force_env_get(bias_env,&
            potential_energy=bias_energy_new,error=error)
! accept or reject the move based on the Metropolis rule
            w=exp(-BETA*(bias_energy_new-&
               bias_energy_old))
         ENDIF

         IF ( w .GE. 1.0E0_dp ) THEN
            w=1.0E0_dp
            rand=0.0E0_dp
         ELSE
            IF(ionode) rand=next_random_number(rng_stream)
!            CALL RANDOM_NUMBER(rand)
            CALL mp_bcast(rand,source,group)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%bias_rot%successes=moves%bias_rot%successes+1
            move_updates%bias_rot%successes=move_updates%bias_rot%successes+1
            moves%rot%qsuccesses=moves%rot%qsuccesses+1
            move_updates%rot%successes=move_updates%rot%successes+1
            bias_energy=bias_energy+bias_energy_new-&
                                 bias_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_old")
      DEALLOCATE(mass,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass")

! end the timing
  CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_molecule_rotation

! *****************************************************************************
!!****s* mc_moves.F/mc_volume_move
!!
!!   NAME
!!     mc_volume_move
!!
!!   FUNCTION
!!     performs a Monte Carlo move that alters the volume of the simulation box
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment whose cell we're changing
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of "outer moves" already performed
!!     - old_energy: the energy of the last accepted move involving an
!!                   unbiased calculation
!!     - globenv: the global simulation parameters
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move involving an
!!              unbiased calculation
!!     - iw: the unit number that writes to the screen
!!     - discrete_array: tells use which volumes we can do for the discrete 
!!           case
!!     - rng_stream: the random number stream that we draw from
!!     - error: the cp error type
!!                                 
!!     Designed for parallel use.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,old_energy,globenv,&
                        energy_check,r_old,iw,discrete_array,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep,iw
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    INTEGER,DIMENSION(1:3, 1:2), &
      INTENT(INOUT)                          :: discrete_array
    REAL(KIND=dp),DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_volume_move", &
      routineP = module_name//':'//routineN
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass

    INTEGER                                  :: handle, i, iatom, &
                                                istat,nchain,&
                                                source,group,nunits_tot,&
                                                nunits,cl,&
                                                natoms,iunit,icurrent,&
                                                iside,iside_change
    REAL(KIND=dp) :: new_energy, rand, vol_dis, w,discrete_step,energy_term,&
      BETA,rmvolume,pressure,rcut,value,&
      pressure_term,volume_term
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                center_of_mass_new, diff,&
                                                abc,test_cell_length,&
                                                new_cell_length,&
                                              old_cell_length
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(force_env_type), POINTER            :: test_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys, testsys
    TYPE(cell_type), POINTER                 :: cell_old
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r
    TYPE(particle_list_type), POINTER        :: particles_old, particles_test
    LOGICAL                                  :: ionode,&
      ltoo_small,ldiscrete,lincrease
    LOGICAL                                  :: loverlap
    CHARACTER(LEN=40)                        :: dat_file
    CHARACTER(LEN=200)                       :: fft_lib
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,&
         nchain=nchain,BETA=BETA,source=source,group=group,&
         nunits_tot=nunits_tot,dat_file=dat_file,&
         rmvolume=rmvolume,nunits=nunits,pressure=pressure,cl=cl,&
         fft_lib=fft_lib,discrete_step=discrete_step,&
         ldiscrete=ldiscrete)

! nullify some pointers
      NULLIFY(particles_old,particles_test,cell_old,oldsys,&
         testsys,test_env)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r",3*nunits_tot*dp_size)

! record the attempt
      moves%volume%attempts=moves%volume%attempts+1
      move_updates%volume%attempts=move_updates%volume%attempts+1

! now let's grab the cell length and particle positions
      CALL force_env_get(force_env,subsys=oldsys,cell=cell_old,&
         error=error)
      CALL get_cell(cell_old,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys,particles=particles_old, &
         molecule_kinds_new=molecule_kinds_new,error=error)

! grab the masses
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,natom=natoms)   
      ALLOCATE(mass(1:natoms),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass",natoms*dp_size)
      DO iatom=1,natoms
         CALL get_atomic_kind(particles_old%els(iatom)%atomic_kind,&
            mass=mass(iatom))
      ENDDO

! find the old cell length
      old_cell_length(1)=abc(1)
      old_cell_length(2)=abc(2)
      old_cell_length(3)=abc(3)
      test_cell_length(1)=abc(1)
      test_cell_length(2)=abc(2)
      test_cell_length(3)=abc(3)

! save the old coordiantes
      DO iatom=1,nunits_tot
         r(1:3,iatom)=particles_old%els(iatom)%r(1:3)
      ENDDO

! fold the coordinates into the central box
      CALL mc_coordinate_fold(r(:,:),nunits_tot,mass,nunits,&
         test_cell_length(1:3))

! make the new dat file
      IF (ionode) THEN
         CALL mc_make_dat_file(r(:,:),nunits_tot,&
            test_cell_length(1:3)*angstrom,&
            dat_file,force_env,TRIM(fft_lib))
      ENDIF

! make a test environment that we can alter to calculate the energy
      globenv%input_file_name=dat_file
      CALL mc_create_force_env(test_env, globenv, .FALSE., &
         force_env%in_use,error=error)

! grab the coords so we can change them
      CALL force_env_get(test_env,subsys=testsys,error=error)
      CALL cp_subsys_get(testsys(1)%subsys, &
             particles=particles_test, error=error)

! now do the move

! call a random number to figure out how far we're moving
      IF (ionode) rand=next_random_number(rng_stream)
!         CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)

! find the test cell lenghts for the discrete volume move
      IF(ldiscrete) THEN
         IF(rand .LT. 0.5_dp) THEN
            lincrease=.true.
         ELSE
            lincrease=.false.
         ENDIF
         
         new_cell_length(1:3)=test_cell_length(1:3)

! if we're increasing the volume, we need to find a side we can increase
         IF(lincrease) THEN
            DO
               IF(ionode) rand=next_random_number(rng_stream)
!     CALL RANDOM_NUMBER(rand)
               CALL mp_bcast(rand,source,group)
               iside_change=CEILING(3.0_dp*rand)
               IF(discrete_array(iside_change,1) .EQ. 1) THEN
                  new_cell_length(iside_change)=&
                  new_cell_length(iside_change)+discrete_step
                  EXIT
               ENDIF
            ENDDO
         ELSE
            DO
               IF(ionode) rand=next_random_number(rng_stream)
!                  CALL RANDOM_NUMBER(rand)
               CALL mp_bcast(rand,source,group)
               iside_change=CEILING(3.0_dp*rand)
               IF(discrete_array(iside_change,2) .EQ. 1) THEN
                  new_cell_length(iside_change)=&
                     new_cell_length(iside_change)-discrete_step
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         vol_dis=(new_cell_length(1)*new_cell_length(2)*new_cell_length(3))&
         -old_cell_length(1)*old_cell_length(2)*old_cell_length(3)
      ELSE
! now for the not discrete volume move
!!!!!!!!!!!!!!!! for E_V curves
         vol_dis=rmvolume*(rand-0.5E0_dp)*2.0E0_dp
!     vol_dis=5.0E0_dp/angstrom**3

         new_cell_length(1:3)=(vol_dis+&
         test_cell_length(1)*test_cell_length(2)*&
         test_cell_length(3))**(1.0E0_dp/3.0E0_dp)
      ENDIF

! save the molecule coordinates
      DO iatom=1,nunits_tot
         r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
      ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass, using the minimum image (not all molecules are in
! the central box)

! now we need to scale the coordinates of all the molecules by the
! center of mass
      DO iatom=1,nunits_tot,nunits

! now find the center of mass
         CALL get_center_of_mass(r(:,iatom:iatom+nunits-1),nunits,&
            center_of_mass(:),mass(:))
            
! scale the center of mass and determine the vector that points from the
!    old COM to the new one
         DO iside=1,3
            center_of_mass_new(iside)=center_of_mass(iside)*&
               new_cell_length(iside)/test_cell_length(iside)
         ENDDO
            
         DO i=1,3
            diff(i)=center_of_mass_new(i)-center_of_mass(i)
! now change the particle positions
            DO iunit=0,nunits-1
               icurrent=iatom+iunit
               particles_test%els(icurrent)%r(i)=&
                  particles_test%els(icurrent)%r(i)+diff(i)
            ENDDO
         ENDDO
      ENDDO

! check for overlap
      CALL check_for_overlap(test_env,&
         nchain,1,nchain,loverlap,cell_length=new_cell_length)

      test_cell_length(1)=new_cell_length(1)
      test_cell_length(2)=new_cell_length(2)
      test_cell_length(3)=new_cell_length(3)

! make the new dat file
      DO iatom=1,nunits_tot
         r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
      ENDDO

      IF(ionode) CALL mc_make_dat_file(r(:,:),nunits_tot,&
         test_cell_length(1:3)*angstrom,&
         dat_file,force_env,TRIM(fft_lib))

! realease the environment
      CALL force_env_release(test_env,error=error)

! figure out if we have overlap problems
      CALL mp_bcast(loverlap,source,group)
      IF(loverlap) THEN
! deallocate some stuff
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r")
         DEALLOCATE(mass,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "mass")

         IF(ionode) WRITE(globenv%scr,*) &
            "Volume move rejected due to overlap.",vol_dis
!     end the timing
         CALL timestop(0.0E0_dp,handle)
         RETURN
      ENDIF

! stop if we're trying to change a box to a boxlength smaller than rcut
      IF(ionode) THEN
         ltoo_small=.FALSE.
         IF(force_env%in_use .EQ. use_fist_force) THEN
            CALL get_mc_par(mc_par,rcut=rcut)
            IF(test_cell_length(1) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.
            IF(test_cell_length(2) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.
            IF(test_cell_length(3) .LT. 2.0_dp*rcut) ltoo_small=.TRUE.
            
            IF(ltoo_small) THEN
               WRITE(iw,*) 'test_cell_lengths ',&
               test_cell_length(1:3)/angstrom
               WRITE(iw,*) 'rcut ',rcut/angstrom
            ENDIF
         ENDIF
      ENDIF
      CALL mp_bcast(ltoo_small,source,group)
      IF(ltoo_small) &
         CALL stop_program("mc_volume_move",&
         "Attempted a volume move where box size got too small.")

      globenv%input_file_name=dat_file
      CALL mc_create_force_env(test_env, globenv, .TRUE., &
         force_env%in_use,error=error)

! now compute the energy
      CALL force_env_calc_energy_force(test_env,calc_force=.FALSE.,error=error)
      test_env%box_changed=.false.
      CALL force_env_get(test_env,&
         potential_energy=new_energy,error=error)
           
! accept or reject the move
! to prevent overflows
      energy_term=new_energy-old_energy
      volume_term=-REAL(nchain,dp)/BETA*&
      LOG(test_cell_length(1)*test_cell_length(2)*test_cell_length(3)/&
         (old_cell_length(1)*old_cell_length(2)*old_cell_length(3)))
      pressure_term=pressure*vol_dis

      value=-BETA*(energy_term+volume_term+pressure_term)
      IF(value .GT. 500.0_dp) THEN
         w=10.0_dp
      ELSEIF(value .LT. -500.0_dp) THEN
         w=0.0_dp
      ELSE
         w=exp(value)
      ENDIF
 
!!!!!!!!!!!!!!!! for E_V curves
!         w=1.0E0_dp
!         w=0.0E0_dp

      IF ( w .GE. 1.0E0_dp ) THEN
         w=1.0E0_dp
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream)
         CALL mp_bcast(rand,source,group)
      ENDIF



      IF (rand .LT. w) THEN

! write cell length, volume, density, and trial displacement to a file

!            IF(ionode) WRITE(cl,*) nnstep,test_cell_length(1:3)*angstrom,&
!               (new_energy-old_energy)/3.16682968E-6_dp,classical_energy_diff/3.16682968E-6_dp,&
!               -REAL(nchain,dp)/&
!               BETA*3.0E0_dp*LOG(test_cell_length(1)/&
!               old_cell_length)/3.16682968E-6_dp, &
!               pressure*(test_cell_length(1)**3-old_cell_length**3)/3.16682968E-6_dp,&
!               vol_dis*(angstrom)**3,rand,new_energy/3.16682968E-6_dp,old_energy/3.16682968E-6_dp,w
               
! accept the move
         moves%volume%successes=moves%volume%successes+1
         move_updates%volume%successes=move_updates%volume%successes+1

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy

! and the coordiantes
         CALL force_env_get(test_env,subsys=testsys,error=error)
         CALL cp_subsys_get(testsys(1)%subsys, &
            particles=particles_test, error=error)
         
         DO iatom=1,nunits_tot
            r_old(1:3,iatom)=particles_test%els(iatom)%r(1:3)
         ENDDO
         
! switch the force env to the new env
         CALL force_env_release(force_env,error=error)
         force_env => test_env

! update discrete_array if we're doing a discrete volume move
         IF(ldiscrete) THEN
            CALL create_discrete_array(test_cell_length(:),&
               discrete_array(:,:),discrete_step)
         ENDIF

      ELSE

! reject the move
! write cell length, volume, density, and trial displacement to a file
!            IF(ionode) WRITE(cl,*) nnstep,old_cell_length*angstrom,&
!               (new_energy-old_energy)/3.16682968E-6_dp,classical_energy_diff/3.16682968E-6_dp,&
!               -REAL(nchain,dp)/&
!               BETA*3.0E0_dp*LOG(test_cell_length(1)/&
!               old_cell_length)/3.16682968E-6_dp, &
!               pressure*(test_cell_length(1)*test_cell_length(2)*&
!               test_cell_length(3)-old_cell_length**3)/3.16682968E-6_dp,&
!               vol_dis*(angstrom)**3,rand,new_energy/3.16682968E-6_dp,old_energy/3.16682968E-6_dp,w
               


! release some memory
         CALL force_env_release(test_env,error=error)

      ENDIF


! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")
      DEALLOCATE(mass,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "mass")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_volume_move

! *****************************************************************************
!!****s* mc_moves.F/change_bond_length
!!
!!   NAME
!!     change_bond_length
!!
!!   FUNCTION
!!     alters the length of a random bond for the given molecule, using
!!     a mass weighted scheme so the lightest atoms move the most
!!
!!   ARGUMENTS
!!     - r_old: the initial coordinates of all molecules in the system
!!     - r_new: the new coordinates of all molecules in the system
!!     - mc_par: the mc parameters for the force env
!!     - mol: the number of the first atom in the molecule, assuming the
!!            next natom-1 atoms are the rest of the molecule
!!     - molecule_kind: the structure containing the molecule information
!!     - dis_length: the ratio of the new bond length to the old bond length,
!!                   used in the acceptance rule
!!     - particles: the particle_list_type for all particles in the force_env..
!!            used to grab the mass of each atom
!!     - rng_stream: the random number stream that we draw from
!!             
!!   This subroutine is written to be parallel. 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE change_bond_length ( r_old,r_new,mc_par,mol,molecule_kind,&
      dis_length,particles,rng_stream)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp),DIMENSION(:,:), &
      INTENT(OUT)                            :: r_new
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_old
    INTEGER, INTENT(IN)                      :: mol
    REAL(KIND=dp), INTENT(OUT)               :: dis_length
    TYPE(molecule_kind_type),POINTER         :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER :: routineN = "change_bond_length"

    TYPE(bond_type),DIMENSION(:),POINTER     :: bond_list
    INTEGER                                  :: handle, i, ipart,&
      nunits_tot,nbond,bond_number,natom,istat,iatom,ibond,&
      group,source
    REAL(KIND=dp)                            :: old_length, rand,rmbond,&
      mass_a,mass_b,atom_mass,new_length_a,new_length_b
    REAL(KIND=dp), DIMENSION(1:3)            :: bond_a,bond_b
    INTEGER,DIMENSION(:),ALLOCATABLE                     :: counter,&
      atom_a,atom_b
    INTEGER,DIMENSION(:,:),ALLOCATABLE       :: connection,connectivity
    LOGICAL                                  :: ionode
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,nunits_tot=nunits_tot,rmbond=rmbond,source=source,&
         group=group,ionode=ionode)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,nbond=nbond,&
         bond_list=bond_list)
      bond_number=CEILING(rand*REAL(nbond,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_b(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_b",natom)
      connection(:,:)=0.0_dp
      connectivity(:,:)=0.0_dp
      counter(:)=0
      atom_a(:)=0
      atom_b(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom a's
! side and which are on atom b's
      atom_a(:)=0
      atom_a(bond_list(bond_number)%a)=1
      CALL depth_first_search(bond_list(bond_number)%a,bond_list(bond_number)%b,&
           connectivity(:,:),atom_a(:))
      atom_b(:)=0
      atom_b(bond_list(bond_number)%b)=1
      CALL depth_first_search(bond_list(bond_number)%b,bond_list(bond_number)%a,&
           connectivity(:,:),atom_b(:))

! now figure out the masses of the various sides, so we can weight how far we move each
! group of atoms
      mass_a=0.0_dp
      mass_b=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)
         IF(atom_a(iatom) == 1) THEN
            mass_a=mass_a+atom_mass
         ELSE
            mass_b=mass_b+atom_mass
         ENDIF
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      CALL mp_bcast(rand,source,group)
      
      dis_length=rmbond*2.0E0_dp*(rand-0.5E0_dp)

! find the bond vector that atom a will be moving
      DO i=1,3
         bond_a(i)=r_new(i,mol+bond_list(bond_number)%a-1)-&
            r_new(i,mol+bond_list(bond_number)%b-1)
         bond_b(i)=-bond_a(i)
      ENDDO

! notice we weight by the opposite masses...therefore lighter segments
! will move further
      old_length=SQRT(DOT_PRODUCT(bond_a,bond_a))
      new_length_a=dis_length*mass_b/(mass_a+mass_b)
      new_length_b=dis_length*mass_a/(mass_a+mass_b)

      DO i=1,3
         bond_a(i)=bond_a(i)/old_length*new_length_a
         bond_b(i)=bond_b(i)/old_length*new_length_b
      ENDDO
      
      DO iatom=1,natom
         IF(atom_a(iatom) == 1) THEN
            r_new(1,mol+iatom-1)=r_new(1,mol+iatom-1)+bond_a(1)
            r_new(2,mol+iatom-1)=r_new(2,mol+iatom-1)+bond_a(2)
            r_new(3,mol+iatom-1)=r_new(3,mol+iatom-1)+bond_a(3)
         ELSE
            r_new(1,mol+iatom-1)=r_new(1,mol+iatom-1)+bond_b(1)
            r_new(2,mol+iatom-1)=r_new(2,mol+iatom-1)+bond_b(2)
            r_new(3,mol+iatom-1)=r_new(3,mol+iatom-1)+bond_b(3)
         ENDIF
      ENDDO

! correct the value of dis_length for the acceptance rule
      dis_length=(old_length+dis_length)/old_length

      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_b")
! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE change_bond_length

! *****************************************************************************
!!****s* mc_moves.F/change_bond_angle
!!
!!   NAME
!!     change_bond_angle
!!
!!   FUNCTION
!!     Alters the magnitude of a random angle in a molecule centered on atom C 
!!     (connected to atoms A and B).  Atoms A and B are moved amounts related 
!!     to their masses (and masses of all connecting atoms), so that heavier
!!     segments are moved less.
!!
!!   ARGUMENTS
!!     - r_old: the initial coordinates of all molecules in the system
!!     - r_new: the new coordinates of all molecules in the system
!!     - mc_par: the mc parameters for the force env
!!     - mol: the number of the first atom in the molecule, assuming the
!!            next natom-1 atoms are the rest of the molecule
!!     - molecule_kind: the structure containing the molecule information
!!     - particles: the particle_list_type for all particles in the force_env...
!!            used to grab the mass of each atom
!!     - rng_stream: the random number stream that we draw from
!!             
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE change_bond_angle ( r_old,r_new,mc_par,mol,molecule_kind,&
      particles,rng_stream)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    REAL(KIND=dp), DIMENSION(:,:), &
      INTENT(IN)                             :: r_old
    INTEGER, INTENT(IN)                      :: mol
    TYPE(molecule_kind_type),POINTER         :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "change_bond_angle"

! local variables

    TYPE(bend_type),DIMENSION(:),POINTER     :: bend_list
    TYPE(bond_type),DIMENSION(:),POINTER     :: bond_list
    INTEGER                                  :: handle, i, ipart,nunits_tot,&
      group,source,nbend,bend_number,natom,istat,iatom,ibond,nbond
    REAL(KIND=dp)                            :: bis_length, dis_angle, &
      old_angle, mass_a,mass_c, rand,&
      rmangle,atom_mass,old_length_a,old_length_c,temp_length,&
      new_angle_a,new_angle_c,dis_angle_a,dis_angle_c
    REAL(KIND=dp), DIMENSION(1:3)            :: bisector, &
      bond_a,bond_c,cross_prod,cross_prod_plane,temp
    LOGICAL                                  :: ionode
    INTEGER,DIMENSION(:),ALLOCATABLE                     :: counter,&
      atom_a,atom_c
    INTEGER,DIMENSION(:,:),ALLOCATABLE       :: connection,connectivity

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,nunits_tot=nunits_tot,rmangle=rmangle,source=source,&
         group=group,ionode=ionode)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,nbend=nbend,&
         bend_list=bend_list,bond_list=bond_list,nbond=nbond)
      bend_number=CEILING(rand*REAL(nbend,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_c(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_c",natom)
      connection(:,:)=0.0_dp
      connectivity(:,:)=0.0_dp
      counter(:)=0
      atom_a(:)=0
      atom_c(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom a's
! side and which are on atom c's
      atom_a(:)=0
      atom_a(bend_list(bend_number)%a)=1
      CALL depth_first_search(bend_list(bend_number)%a,bend_list(bend_number)%b,&
           connectivity(:,:),atom_a(:))
      atom_c(:)=0
      atom_c(bend_list(bend_number)%c)=1
      CALL depth_first_search(bend_list(bend_number)%c,bend_list(bend_number)%b,&
           connectivity(:,:),atom_c(:))

! now figure out the masses of the various sides, so we can weight how far we move each
! group of atoms
      mass_a=0.0_dp
      mass_c=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)  
         IF(atom_a(iatom) == 1) mass_a=mass_a+atom_mass
         IF(atom_c(iatom) == 1) mass_c=mass_c+atom_mass
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream)
      CALL mp_bcast(rand,source,group)
      
      dis_angle=rmangle*2.0E0_dp*(rand-0.5E0_dp)

! need to find the A-B-C bisector

! this going to be tough...we need to find the plane of the A-B-C bond and only shift 
! that component for all atoms connected to A and C...otherwise we change other 
! internal degrees of freedom

! find the bond vectors
      DO i=1,3
         bond_a(i)=r_new(i,mol-1+bend_list(bend_number)%a)-&
            r_new(i,mol-1+bend_list(bend_number)%b)
         bond_c(i)=r_new(i,mol-1+bend_list(bend_number)%c)-&
            r_new(i,mol-1+bend_list(bend_number)%b)
      ENDDO
      old_length_a=SQRT(DOT_PRODUCT(bond_a,bond_a))
      old_length_c=SQRT(DOT_PRODUCT(bond_c,bond_c))
      old_angle=ACOS(DOT_PRODUCT(bond_a,bond_c)/(old_length_a*old_length_c))

      DO i=1,3
         bisector(i)=bond_a(i)/old_length_a+& ! not yet normalized
                 bond_c(i)/old_length_c
      ENDDO
      bis_length=SQRT(DOT_PRODUCT(bisector,bisector))
      bisector(1:3)=bisector(1:3)/bis_length

! now we need to find the cross product of the B-A and B-C vectors and normalize
! it, so we have a vector that defines the bend plane
      cross_prod(1)=bond_a(2)*bond_c(3)-bond_a(3)*bond_c(2)
      cross_prod(2)=bond_a(3)*bond_c(1)-bond_a(1)*bond_c(3)
      cross_prod(3)=bond_a(1)*bond_c(2)-bond_a(2)*bond_c(1)
      cross_prod(1:3)=cross_prod(1:3)/SQRT(DOT_PRODUCT(cross_prod,cross_prod))

! we have two axis of a coordinate system...let's get the third
      cross_prod_plane(1)=cross_prod(2)*bisector(3)-cross_prod(3)*bisector(2)
      cross_prod_plane(2)=cross_prod(3)*bisector(1)-cross_prod(1)*bisector(3)
      cross_prod_plane(3)=cross_prod(1)*bisector(2)-cross_prod(2)*bisector(1)
      cross_prod_plane(1:3)=cross_prod_plane(1:3)/&
         SQRT(DOT_PRODUCT(cross_prod_plane,cross_prod_plane))

! now bisector is x, cross_prod_plane is the y vector (pointing towards c), 
! and cross_prod is z
! shift the molecule so that atom b is at the origin
      DO iatom=1,natom
         r_new(1:3,mol-1+iatom)=r_new(1:3,mol-1+iatom)-&
            r_old(1:3,mol-1+bend_list(bend_number)%b)
      ENDDO

! figure out how much we move each side, since we're mass-weighting, by the
! opposite masses, so lighter moves farther..this angle is the angle between
! the bond vector BA or BC and the bisector
      dis_angle_a=dis_angle*mass_c/(mass_a+mass_c)
      dis_angle_c=dis_angle*mass_a/(mass_a+mass_c)      

! now loop through all the atoms, moving the ones that are connected to a or c
      DO iatom=1,natom
! subtract out the z component (perpendicular to the angle plane)
         temp(1:3)=r_new(1:3,mol-1+iatom)-&
            DOT_PRODUCT(cross_prod(1:3),r_new(1:3,mol-1+iatom))*&
            cross_prod(1:3)
         temp_length=SQRT(DOT_PRODUCT(temp,temp))

! we can now compute all three components of the new bond vector along the
! axis defined above
         IF(atom_a(iatom) == 1) THEN

! if the y-coordinate is less than zero, we need to switch the sign when we make the vector,
! as the angle computed by the dot product can't distinguish between that
            IF(DOT_PRODUCT(cross_prod_plane(1:3),r_new(1:3,mol-1+iatom)) &
               .LT. 0.0_dp) THEN
               
! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_a=ACOS(DOT_PRODUCT(bisector,temp(1:3))/&
               (temp_length))+dis_angle_a

            r_new(1:3,mol-1+iatom)=cos(new_angle_a)*temp_length*bisector(1:3)-&
               sin(new_angle_a)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,mol-1+iatom))*&
               cross_prod(1:3)
            ELSE

! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_a=ACOS(DOT_PRODUCT(bisector,temp(1:3))/&
               (temp_length))-dis_angle_a

            r_new(1:3,mol-1+iatom)=cos(new_angle_a)*temp_length*bisector(1:3)+&
               sin(new_angle_a)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,mol-1+iatom))*&
               cross_prod(1:3)
            ENDIF

         ELSEIF(atom_c(iatom) == 1) THEN

! if the y-coordinate is less than zero, we need to switch the sign when we make the vector,
! as the angle computed by the dot product can't distinguish between that
            IF(DOT_PRODUCT(cross_prod_plane(1:3),r_new(1:3,mol-1+iatom)) &
               .LT. 0.0_dp) THEN
! need to figure out the current iatom-B-bisector angle, so we know what the new angle is
            new_angle_c=ACOS(DOT_PRODUCT(bisector(1:3),temp(1:3))/&
               (temp_length))-dis_angle_c

            r_new(1:3,mol-1+iatom)=cos(new_angle_c)*temp_length*bisector(1:3)-&
               sin(new_angle_c)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,mol-1+iatom))*&
               cross_prod(1:3)
            ELSE
            new_angle_c=ACOS(DOT_PRODUCT(bisector(1:3),temp(1:3))/&
               (temp_length))+dis_angle_c

            r_new(1:3,mol-1+iatom)=cos(new_angle_c)*temp_length*bisector(1:3)+&
               sin(new_angle_c)*temp_length*cross_prod_plane(1:3)+&
               DOT_PRODUCT(cross_prod(1:3),r_new(1:3,mol-1+iatom))*&
               cross_prod(1:3)
            ENDIF
         ENDIF

      ENDDO
      
      DO iatom=1,natom
         r_new(1:3,mol-1+iatom)=r_new(1:3,mol-1+iatom)+&
            r_old(1:3,mol-1+bend_list(bend_number)%b)
      ENDDO

! deallocate some stuff
      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_c,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_c")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE change_bond_angle

! *****************************************************************************
!!****s* mc_moves.F/change_dihedral
!!
!!   NAME
!!     change_dihedral
!!
!!   FUNCTION
!!     Alters a dihedral (A-B-C-D) in the molecule so that all other internal
!!     degrees of freedom remain the same.  If other dihedrals are centered
!!     on B-C, they rotate as well to keep the relationship between the
!!     dihedrals the same.  Atoms A and D are moved amounts related to their
!!     masses (and masses of all connecting atoms), so that heavier segments
!!     are moved less.  All atoms except B and C are rotated around the 
!!     B-C bond vector (B and C are not moved).
!!
!!   ARGUMENTS
!!     - r_old: the initial coordinates of all molecules in the system
!!     - r_new: the new coordinates of all molecules in the system
!!     - mc_par: the mc parameters for the force env
!!     - mol: the number of the first atom in the molecule, assuming the
!!            next natom-1 atoms are the rest of the molecule
!!     - molecule_kind: the structure containing the molecule information
!!     - particles: the particle_list_type for all particles in the force_env..
!!            used to grab the mass of each atom
!!     - rng_stream: the random number stream that we draw from
!!             
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE change_dihedral ( r_old,r_new,mc_par,mol,molecule_kind,&
      particles,rng_stream)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    REAL(KIND=dp), DIMENSION(:,:), &
      INTENT(IN)                             :: r_old
    INTEGER, INTENT(IN)                      :: mol
    TYPE(molecule_kind_type),POINTER         :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "change_dihedral"

! local variables

    TYPE(bond_type),DIMENSION(:),POINTER     :: bond_list
    TYPE(torsion_type),DIMENSION(:),POINTER     :: torsion_list
    INTEGER                                  :: handle, i, ipart,nunits_tot,&
      group,source,natom,istat,iatom,ibond,nbond,&
      ntorsion,torsion_number
    REAL(KIND=dp)                            :: dis_angle, &
      mass_a,mass_d, rand,rmdihedral,x,y,z,u,v,w,&
      atom_mass,old_length_a,&
      dis_angle_a,dis_angle_d
    REAL(KIND=dp), DIMENSION(1:3)            :: bond_a,temp
    LOGICAL                                  :: ionode
    INTEGER,DIMENSION(:),ALLOCATABLE                     :: counter,&
      atom_a,atom_d
    INTEGER,DIMENSION(:,:),ALLOCATABLE       :: connection,connectivity

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,nunits_tot=nunits_tot,rmdihedral=rmdihedral,&
         source=source,group=group,ionode=ionode)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      IF(ionode) THEN
         rand=next_random_number(rng_stream)
!      CALL RANDOM_NUMBER(rand)
      ENDIF
      CALL mp_bcast(rand,source,group)
      CALL get_molecule_kind(molecule_kind,natom=natom,&
         bond_list=bond_list,nbond=nbond,&
         ntorsion=ntorsion,torsion_list=torsion_list)
      torsion_number=CEILING(rand*REAL(ntorsion,dp))

      ALLOCATE(connection(1:natom,1:2),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection",2*natom*dp_size)
! assume at most six bonds per atom
      ALLOCATE(connectivity(1:6,1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity",6*natom*dp_size)
      ALLOCATE(counter(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter",natom)
      ALLOCATE(atom_a(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a",natom)
      ALLOCATE(atom_d(1:natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_d",natom)
      connection(:,:)=0.0_dp
      connectivity(:,:)=0.0_dp
      counter(:)=0
      atom_a(:)=0
      atom_d(:)=0

! now we need to find a list of atoms that each atom in this bond is connected
! to
      DO iatom=1,natom
         DO ibond=1,nbond
            IF(bond_list(ibond)%a == iatom) THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%b
            ELSEIF(bond_list(ibond)%b == iatom)THEN
               counter(iatom)=counter(iatom)+1
               connectivity(counter(iatom),iatom)=bond_list(ibond)%a
            ENDIF
         ENDDO
      ENDDO

! now I need to do a depth first search to figure out which atoms are on atom
! a's side and which are on atom d's, but remember we're moving all atoms on a's
! side of b, including atoms not in a's branch
      atom_a(:)=0
      atom_a(torsion_list(torsion_number)%a)=1
      CALL depth_first_search(torsion_list(torsion_number)%b,&
         torsion_list(torsion_number)%c,connectivity(:,:),atom_a(:))
      atom_d(:)=0
      atom_d(torsion_list(torsion_number)%d)=1
      CALL depth_first_search(torsion_list(torsion_number)%c,&
         torsion_list(torsion_number)%b,connectivity(:,:),atom_d(:))

! now figure out the masses of the various sides, so we can weight how far we
! move each group of atoms
      mass_a=0.0_dp
      mass_d=0.0_dp
      DO iatom=1,natom
         CALL get_atomic_kind(particles%els(iatom)%atomic_kind,&
            mass=atom_mass)  
         IF(atom_a(iatom) == 1) mass_a=mass_a+atom_mass
         IF(atom_d(iatom) == 1) mass_d=mass_d+atom_mass
      ENDDO

! choose a displacement
      IF(ionode) rand=next_random_number(rng_stream)
      CALL mp_bcast(rand,source,group)
      
      dis_angle=rmdihedral*2.0E0_dp*(rand-0.5E0_dp)

! find the bond vectors, B-C, so we know what to rotate around
      DO i=1,3
         bond_a(i)=r_new(i,mol-1+torsion_list(torsion_number)%c)-&
            r_new(i,mol-1+torsion_list(torsion_number)%b)
      ENDDO
      old_length_a=SQRT(DOT_PRODUCT(bond_a,bond_a))
      bond_a(1:3)=bond_a(1:3)/old_length_a

! figure out how much we move each side, since we're mass-weighting, by the
! opposite masses, so lighter moves farther...we take the opposite sign of d
! so we're not rotating both angles in the same direction
      dis_angle_a=dis_angle*mass_d/(mass_a+mass_d)
      dis_angle_d=-dis_angle*mass_a/(mass_a+mass_d)   

      DO iatom=1,natom
             
         IF(atom_a(iatom) == 1) THEN
! shift the coords so b is at the origin
            r_new(1:3,mol-1+iatom)=r_new(1:3,mol-1+iatom)-&
               r_new(1:3,mol-1+torsion_list(torsion_number)%b)
               

! multiply by the rotation matrix
               u=bond_a(1)
               v=bond_a(2)
               w=bond_a(3)
               x=r_new(1,mol-1+iatom)
               y=r_new(2,mol-1+iatom)
               z=r_new(3,mol-1+iatom)
               temp(1)=(u*(u*x+v*y+w*z)+(x*(v**2+w**2)-u*(v*y+w*z))*cos(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(v*z-w*y)*sin(dis_angle_a))/(u**2+v**2+w**2)
               temp(2)=(v*(u*x+v*y+w*z)+(y*(u**2+w**2)-v*(u*x+w*z))*cos(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(w*x-u*z)*sin(dis_angle_a))/(u**2+v**2+w**2)
               temp(3)=(w*(u*x+v*y+w*z)+(z*(v**2+u**2)-w*(u*x+v*y))*cos(dis_angle_a)+&
                  SQRT(u**2+v**2+w**2)*(u*y-v*x)*sin(dis_angle_a))/(u**2+v**2+w**2)

! shift back to the original position
            temp(1:3)=temp(1:3)+r_new(1:3,mol-1+torsion_list(torsion_number)%b)
            r_new(1:3,mol-1+iatom)=temp(1:3)

         ELSEIF(atom_d(iatom) == 1) THEN

! shift the coords so c is at the origin
            r_new(1:3,mol-1+iatom)=r_new(1:3,mol-1+iatom)-&
               r_new(1:3,mol-1+torsion_list(torsion_number)%c)
               
! multiply by the rotation matrix
               u=bond_a(1)
               v=bond_a(2)
               w=bond_a(3)
               x=r_new(1,mol-1+iatom)
               y=r_new(2,mol-1+iatom)
               z=r_new(3,mol-1+iatom)
               temp(1)=(u*(u*x+v*y+w*z)+(x*(v**2+w**2)-u*(v*y+w*z))*cos(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(v*z-w*y)*sin(dis_angle_d))/(u**2+v**2+w**2)
               temp(2)=(v*(u*x+v*y+w*z)+(y*(u**2+w**2)-v*(u*x+w*z))*cos(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(w*x-u*z)*sin(dis_angle_d))/(u**2+v**2+w**2)
               temp(3)=(w*(u*x+v*y+w*z)+(z*(v**2+u**2)-w*(u*x+v*y))*cos(dis_angle_d)+&
                  SQRT(u**2+v**2+w**2)*(u*y-v*x)*sin(dis_angle_d))/(u**2+v**2+w**2)

! shift back to the original position
            temp(1:3)=temp(1:3)+r_new(1:3,mol-1+torsion_list(torsion_number)%c)
            r_new(1:3,mol-1+iatom)=temp(1:3)
         ENDIF
      ENDDO

! deallocate some stuff
      DEALLOCATE(connection,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connection")
      DEALLOCATE(connectivity,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "connectivity")
      DEALLOCATE(counter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "counter")
      DEALLOCATE(atom_a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_a")
      DEALLOCATE(atom_d,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "atom_d")

! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE change_dihedral
!*********************************************************************

! *****************************************************************************
!!****s* mc_moves.F/mc_avbmc_move
!!
!!   NAME
!!     mc_avbmc_move
!!
!!   FUNCTION
!!     performs either a bond or angle change move for a given molecule
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - bias_env: the force environment used to bias the move, if any (it may
!!           be null if lbias=.false. in mc_par)
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - enegy_check: the running energy difference between now and the inital energy
!!     - r_old: the coordinates of force_env before the move
!!     - old_energy: the energy of the force_env before the move
!!     - mol_swap: the number of the moving molecule's first atom, assuming the rest of
!!       the atoms follow sequentially
!!     - mol_target: the number of the target molecule's first atom
!!     - bias_energy: the biased energy of the system before the move
!!     - move_type: dictates if we're moving to an "in" or "out" region
!!     - rng_stream: the random number stream that we draw from
!!     - error: the cp_error type                            
!!
!!     Designed for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_avbmc_move ( mc_par,force_env,bias_env, moves,&
                        energy_check,r_old,old_energy,mol_swap,&
                        mol_target,bias_energy_old,move_type,rng_stream,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves
    INTEGER, INTENT(IN)                      :: mol_swap,mol_target
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy_old,old_energy,&
      energy_check
    REAL(KIND=dp), DIMENSION(:,:),INTENT(INOUT)     &
      :: r_old
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_avbmc_move"

    INTEGER                                  :: handle, ipart,&
                                                istat,nchain,nunits_tot,&
                                                avbmc_atom,nswapmoves,&
                                                source,group,imolecule,&
                                                natom
    LOGICAL                                  :: lbias,&
      loverlap,ionode,ldum,lin
    REAL(KIND=dp),DIMENSION(1:3)             :: abc,RIJ
    REAL(KIND=dp)                            :: bias_energy_new, &
      rand, w, BETA,&
      avbmc_rmin,avbmc_rmax,rdum,new_energy,weight_new,weight_old,&
      volume_in,volume_out,del_quickstep_energy,prefactor,distance,pbias
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys,subsys_force
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_new
    TYPE(particle_list_type), POINTER        :: particles,particles_force
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
      TYPE(cell_type), POINTER                 :: cell

!   ---------------------------------------------------------------------------

      rdum=1.0_dp

! begin the timing of the subroutine
      CALL timeset(routineN,'I','',handle)
      
! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,avbmc_atom=avbmc_atom,&
         avbmc_rmin=avbmc_rmin,avbmc_rmax=avbmc_rmax,&
         nswapmoves=nswapmoves,ionode=ionode,source=source,&
         group=group,pbias=pbias)

      
! nullify some pointers
      NULLIFY(particles,subsys,molecule_kinds_new,molecule_kind,&
         particles_force,subsys_force)

! do some allocation
      ALLOCATE (r_new(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new",3*nunits_tot*dp_size)

! now we need to grab and save coordinates, in case we reject
! are we biasing this move?
      IF(lbias) THEN

! grab the coordinates
         CALL force_env_get(bias_env,cell=cell,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)
         CALL get_molecule_kind(molecule_kind,natom=natom)
         CALL get_cell(cell,abc=abc)

! save the energy
!         bias_energy_old=bias_energy

      ELSE

! grab the coordinates
         CALL force_env_get(force_env,cell=cell,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles,  molecule_kinds_new=molecule_kinds_new,&
            error=error)
         molecule_kind => molecule_kinds_new%els(1)
         CALL get_molecule_kind(molecule_kind,natom=natom)
         CALL get_cell(cell,abc=abc)

      ENDIF

! let's determine if the molecule to be moved is in the "in" region or the
! "out" region of the target
      RIJ(1)=particles%els(mol_swap+avbmc_atom-1)%r(1)-&
         particles%els(mol_target+avbmc_atom-1)%r(1)-abc(1)*ANINT(&
         (particles%els(mol_swap+avbmc_atom-1)%r(1)-&
         particles%els(mol_target+avbmc_atom-1)%r(1))/abc(1))
      RIJ(2)=particles%els(mol_swap+avbmc_atom-1)%r(2)-&
         particles%els(mol_target+avbmc_atom-1)%r(2)-abc(2)*ANINT(&
         (particles%els(mol_swap+avbmc_atom-1)%r(2)-&
         particles%els(mol_target+avbmc_atom-1)%r(2))/abc(2))
      RIJ(3)=particles%els(mol_swap+avbmc_atom-1)%r(3)-&
         particles%els(mol_target+avbmc_atom-1)%r(3)-abc(3)*ANINT(&
         (particles%els(mol_swap+avbmc_atom-1)%r(3)-&
         particles%els(mol_target+avbmc_atom-1)%r(3))/abc(3))
      distance=SQRT(RIJ(1)**2+RIJ(2)**2+RIJ(3)**2)
      IF(distance .LE. avbmc_rmax .AND. distance .GE. avbmc_rmin) THEN
         lin=.true.
      ELSE
         lin=.false.
      ENDIF

! increment the counter of the particular move we've done
!     swapping into the "in" region of mol_target         
      IF(lin) THEN
         IF(move_type == 'in' ) THEN
            moves%avbmc_inin%attempts=moves%avbmc_inin%attempts+1
         ELSE
            moves%avbmc_inout%attempts=moves%avbmc_inout%attempts+1
         ENDIF
      ELSE
         IF(move_type == 'in' ) THEN
            moves%avbmc_outin%attempts=moves%avbmc_outin%attempts+1
         ELSE
            moves%avbmc_outout%attempts=moves%avbmc_outout%attempts+1
         ENDIF
      ENDIF

! now try the move
      imolecule=(mol_swap+natom-1)/natom
      IF(lbias) THEN

         IF(move_type == 'in') THEN
            
! do CBMC for the old config
            CALL generate_cbmc_swap_config( bias_env, BETA, nswapmoves,&
               weight_old,imolecule,ldum, rdum,&
               bias_energy_old,ionode,source,group,rng_stream,error,&
               avbmc_atom=avbmc_atom,&
               rmin=avbmc_rmin,rmax=avbmc_rmax,move_type='out',&
               mol_target=mol_target,rem_molecule=imolecule)
            
         ELSE
            
! do CBMC for the old config
           CALL generate_cbmc_swap_config( bias_env, BETA, nswapmoves,&
               weight_old,imolecule,ldum, rdum,&
               bias_energy_old,ionode,source,group,rng_stream,error,&
               avbmc_atom=avbmc_atom,&
               rmin=avbmc_rmin,rmax=avbmc_rmax,move_type='in',&
               mol_target=mol_target,rem_molecule=imolecule)
            
         ENDIF

! generate the new config
         CALL generate_cbmc_swap_config( bias_env, BETA, nswapmoves,&
            weight_new,imolecule,loverlap, bias_energy_new,&
            bias_energy_old,ionode,source,group,rng_stream,error,&
            avbmc_atom=avbmc_atom,&
            rmin=avbmc_rmin,rmax=avbmc_rmax,move_type=move_type,&
            mol_target=mol_target)

! the energy that comes out of the above routine is the difference...we want
! the real energy for the acceptance rule...we don't do this for the 
! lbias=.false. case because it doesn't appear in the acceptance rule, and
! we compensate in case of acceptance
         bias_energy_new=bias_energy_new+bias_energy_old


      ELSE

         IF(move_type == 'in') THEN
            
! find the weight of the old config
            CALL generate_cbmc_swap_config( force_env, BETA, nswapmoves,&
               weight_old,imolecule,ldum,rdum,old_energy,&
               ionode,source,group,rng_stream,error,avbmc_atom=avbmc_atom,&
               rmin=avbmc_rmin,rmax=avbmc_rmax,move_type='out',&
               mol_target=mol_target,rem_molecule=imolecule)
            
         ELSE
            
! find the weight of the old config
            CALL generate_cbmc_swap_config( force_env, BETA, nswapmoves,&
               weight_old,imolecule,ldum,rdum,old_energy,&
               ionode,source,group,rng_stream,error,avbmc_atom=avbmc_atom,&
               rmin=avbmc_rmin,rmax=avbmc_rmax,move_type='in',&
               mol_target=mol_target,rem_molecule=imolecule)
            
         ENDIF

 ! generate the new config...do this after, because it changes the force_env
        CALL generate_cbmc_swap_config( force_env, BETA, nswapmoves,&
            weight_new,imolecule,loverlap,new_energy,old_energy,&
            ionode,source,group,rng_stream,error,avbmc_atom=avbmc_atom,&
            rmin=avbmc_rmin,rmax=avbmc_rmax,move_type=move_type,&
            mol_target=mol_target)

      ENDIF

      IF(loverlap) THEN
         DEALLOCATE(r_new,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r_new")

! need to reset the old coordinates
         IF(lbias) THEN
            CALL force_env_get(bias_env,&
            subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles,error=error)
         ELSE
            CALL force_env_get(force_env,&
               subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
               particles=particles,error=error)
         ENDIF
         DO ipart=1,nunits_tot
            particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
         ENDDO         

         CALL timestop(0.0E0_dp,handle)

         RETURN
      ENDIF

! if we're biasing, we need to compute the new energy with the full
! potential
      IF(lbias) THEN
! need to give the force_env the coords from the bias_env
         CALL force_env_get(force_env,&
            subsys=subsys_force,error=error)
         CALL cp_subsys_get(subsys_force(1)%subsys, &
            particles=particles_force,error=error)
         CALL force_env_get(bias_env,&
            subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles,error=error)
         DO ipart=1,nunits_tot
            particles_force%els(ipart)%r(1:3)=particles%els(ipart)%r(1:3)
         ENDDO

         CALL force_env_calc_energy_force(force_env,&
            calc_force=.FALSE.,error=error)
         CALL force_env_get(force_env,&
            potential_energy=new_energy,error=error)

      ENDIF

      volume_in=4.0_dp/3.0_dp*pi*(avbmc_rmax**3-avbmc_rmin**3)
      volume_out=abc(1)*abc(2)*abc(3)-volume_in

      IF(lin .AND. move_type == 'in' .OR. &
         .NOT. lin .AND. move_type == 'out') THEN
! standard Metropolis rule
         prefactor=1.0_dp
      ELSEIF(.NOT. lin .AND. move_type == 'in') THEN
         prefactor=(1.0_dp-pbias)*volume_in/(pbias*volume_out)
      ELSE
         prefactor=pbias*volume_out/((1.0_dp-pbias)*volume_in)
      ENDIF

      IF(lbias) THEN
! AVBMC with CBMC and a biasing potential...notice that if the biasing
! potential equals the quickstep potential, this cancels out to the 
! acceptance below
         del_quickstep_energy=(-BETA)*(new_energy-old_energy-&
            (bias_energy_new-bias_energy_old))
         
         IF(del_quickstep_energy .GT. 500.0_dp) THEN
            del_quickstep_energy=1.4E217_dp
         ELSEIF(del_quickstep_energy .LT. -500.0_dp) THEN
            del_quickstep_energy=0.0_dp
         ELSE
            del_quickstep_energy=EXP(del_quickstep_energy)
         ENDIF
         
         w=prefactor*del_quickstep_energy*weight_new/weight_old
         
      ELSE

! AVBMC with CBMC
         w=prefactor*weight_new/weight_old
!     WRITE(6,*) 'einvodi ',new_energy,old_energy,lin,move_type
!     WRITE(6,*) 'invodi ',weight_new,weight_old,w
      ENDIF
         
! check if the move is accepted
      IF(w .GE. 1.0E0_dp) THEN
         rand=0.0E0_dp
      ELSE
         IF(ionode) rand=next_random_number(rng_stream)
         CALL mp_bcast(rand,source,group)
      ENDIF
      
      IF ( rand .LT. w ) THEN

! accept the move         

         IF(lin) THEN
            IF(move_type == 'in' ) THEN
               moves%avbmc_inin%successes=moves%avbmc_inin%successes+1
            ELSE
               moves%avbmc_inout%successes=moves%avbmc_inout%successes+1
            ENDIF
         ELSE
            IF(move_type == 'in' ) THEN
               moves%avbmc_outin%successes=moves%avbmc_outin%successes+1
            ELSE
               moves%avbmc_outout%successes=moves%avbmc_outout%successes+1
            ENDIF
         ENDIF

! we need to compensate for the fact that we take the difference in
! generate_cbmc_config to keep the exponetials small
         IF(.NOT. lbias) THEN
            new_energy=new_energy+old_energy
         ENDIF

! update energies 
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy

! if we're biasing the update the biasing energy
         IF (lbias) THEN
            moves%bias_energy=bias_energy_new
            bias_energy_old=bias_energy_new
         ENDIF

! update coordinates
         CALL force_env_get(force_env,&
            subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles,error=error)
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO


      ELSE

!         WRITE(6,*) 'REJECTED'

! reject the move...need to restore the old coordinates
         IF(lbias) THEN
            CALL force_env_get(bias_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
               particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)
         ENDIF
         CALL force_env_get(force_env,subsys=subsys,error=error)
         CALL cp_subsys_get(subsys(1)%subsys, &
            particles=particles, error=error)
         DO ipart=1,nunits_tot
            particles%els(ipart)%r(1:3)=r_old(1:3,ipart)
         ENDDO
         CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

      ENDIF      

! deallocate some stuff
      DEALLOCATE(r_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new")
! end the timing
  CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE mc_avbmc_move

END MODULE mc_moves

!**********************************************************************

