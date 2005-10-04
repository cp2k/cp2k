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
!!     translation of an atom, translation of a molecule, rotation
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
#include "cp_prep_globals.h"

MODULE mc_moves
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level
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
  USE mc_control,                      ONLY: mc_create_force_env
  USE mc_coordinates,                  ONLY: get_center_of_mass,&
                                             mc_classical_energy,&
                                             mc_coordinate_fold,&
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
                                             atom_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE  :: change_bond_angle,change_bond_length

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_moves"

  PUBLIC :: mc_atom_translation,mc_molecule_translation,&
            mc_molecule_rotation,mc_volume_move
 
! *****************************************************************************

CONTAINS

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
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of "outer moves" already performed
!!     - counter: the total number of "inner moves" already performed
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!                                 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_atom_translation ( mc_par,force_env,bias_env, moves,&
                        move_updates,nnstep,counter,mol,bias_energy,&
                        move_type,lreject,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    CHARACTER(LEN=*), INTENT(IN)             :: move_type
    LOGICAL, INTENT(OUT)                     :: lreject

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_atom_translation"

    INTEGER                                  :: handle, imol, ipart,&
                                                istat,nchain,nunits_tot,&
                                                nunits
    LOGICAL                                  :: bond_flag,lbias,&
      loverlap
    REAL(KIND=dp)                            :: bias_energy_new, &
                                                bias_energy_old, &
                                                dis_length, rand, w, BETA,&
                                                cutoff,rdum
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_old, r_new
    TYPE(particle_list_type), POINTER        :: particles

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)
      
! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,inner_cutoff=cutoff,&
         nunits=nunits)

! nullify some pointers
      NULLIFY(particles,subsys)
      bond_flag=.FALSE.

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

! save the coordinates
         DO ipart=1,nunits_tot
            r_old(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO

      ENDIF

! do the move
      IF (move_type == 'bond') THEN

! change a bond length
         bond_flag=.TRUE.

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
            mc_par,mol,dis_length)

      ELSE

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
         CALL change_bond_angle(r_old,r_new,mc_par,mol)
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
         CALL mc_classical_energy(bias_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
      ELSE
         CALL mc_classical_energy(force_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.)
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
            CALL RANDOM_NUMBER(rand)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            IF (bond_flag) THEN
               moves%bond%qsuccesses=moves%bond%qsuccesses+1
               move_updates%bond%successes=&
                  move_updates%bond%successes+1
               moves%bias_bond%successes=moves%bias_bond%successes+1
               move_updates%bias_bond%successes=&
                  move_updates%bias_bond%successes+1
            ELSE
               moves%angle%qsuccesses=moves%angle%qsuccesses+1
               move_updates%angle%successes=&
                  move_updates%angle%successes+1
               moves%bias_angle%successes=moves%bias_angle%successes+1
               move_updates%bias_angle%successes=&
                  move_updates%bias_angle%successes+1
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
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of "outer moves" already performed
!!     - counter: the total number of "inner moves" already performed
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!                                 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

 SUBROUTINE mc_molecule_translation(  mc_par,force_env, bias_env,moves,&
                      move_updates,nnstep,counter,mol,bias_energy,&
                        lreject,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    LOGICAL,INTENT(OUT)                      :: lreject

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_translation"

    INTEGER                                  :: handle, i, imol, ipart, &
                                                iparticle,istat,nchain,&
                                                nunits_tot,nunits
    REAL(KIND=dp)                            :: bias_energy_new, &
                                                bias_energy_old, &
                                                dis_mol, rand, w, BETA,&
                                          rdum,rmtrans,cutoff,energy_beta
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r_old
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: lbias,loverlap

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lbias=lbias,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,rmtrans=rmtrans,&
         nunits=nunits,inner_cutoff=cutoff)

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

!         CALL mc_classical_energy(force_env,classical_energy_old,&
!            nchain,imol,imol,cutoff,loverlap)
!         IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
!            'Rotation found an overlap in the old config') 
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
      CALL RANDOM_NUMBER(rand)
      i=AINT(3*rand)+1

! call a random number to figure out how far we're moving
      CALL RANDOM_NUMBER(rand)
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
         CALL mc_classical_energy(bias_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
      ELSE
         CALL mc_classical_energy(force_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing with a cheaper potential, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 
         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.)
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
            CALL RANDOM_NUMBER(rand)
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
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - nnstep: the total number of "outer moves" already performed
!!     - counter: the total number of "inner moves" already performed
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!                                 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, bias_env,moves,&
                      move_updates,nnstep,counter,mol,bias_energy,&
                        lreject,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env,bias_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: bias_energy
    LOGICAL,INTENT(OUT)                      :: lreject
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
 
    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_rotation"

    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass 
    INTEGER                                  :: dir, handle, ii, imol, ipart, &
                                                iunit,istat,nchain,nunits_tot,&
                                                nunits,iatom,natoms
    LOGICAL                                  :: lx, ly,lbias,loverlap
    REAL(KIND=dp) :: bias_energy_new, bias_energy_old, cosdg, &
      dgamma, nxcm, nycm, nzcm, rand, rx, rxnew, ry, rynew, rz, rznew, sindg, &
      w,BETA,rmrot,cutoff,rdum, masstot
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
         inner_cutoff=cutoff)

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

!         CALL mc_classical_energy(force_env,classical_energy_old,&
!            nchain,imol,imol,cutoff,loverlap)
!         IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
!            'Rotation found an overlap in the old config') 
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
      CALL RANDOM_NUMBER(rand)
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
      CALL RANDOM_NUMBER(rand)
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
         CALL mc_classical_energy(bias_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
      ELSE
         CALL mc_classical_energy(force_env,rdum,&
            nchain,imol,imol,cutoff,loverlap)
         IF(loverlap) lreject=.TRUE.
      ENDIF

! if we're biasing classical, check for acceptance
      IF(lbias) THEN

! here's where we bias the moves 

         IF(loverlap) THEN
            w=0.0E0_dp
         ELSE
            CALL force_env_calc_energy_force(bias_env,calc_force=.FALSE.)
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
            CALL RANDOM_NUMBER(rand)
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

  SUBROUTINE mc_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,old_energy,&
                        globenv,&
                        energy_check,r_old,iw,discrete_array,error)

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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_volume_move", &
      routineP = module_name//':'//routineN
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: mass

    INTEGER                                  :: handle, i, iatom, &
                                                ivolmove,istat,nchain,&
                                                source,group,nunits_tot,&
                                                nvolmoves,nunits,cl,&
                                                natoms,iunit,icurrent,&
                                                iside,iside_change
    REAL(KIND=dp) :: classical_energy_diff, classical_energy_new, &
      classical_energy_original, classical_energy_test, &
      new_energy, rand, vol_dis, w,discrete_step,energy_term,&
      BETA,rmvolume,pressure,cutoff,inner_cutoff,rcut,rdum,value,&
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
      ALLOCATABLE                            :: r,r_mol
    TYPE(particle_list_type), POINTER        :: particles_old, particles_test
    LOGICAL                                  :: ldum,ionode,lbias,&
      ltoo_small,ldiscrete,lincrease
    LOGICAL,DIMENSION(:),ALLOCATABLE         :: loverlap
    CHARACTER(LEN=40)                        :: dat_file
    CHARACTER(LEN=200)                       :: fft_lib
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,lbias=lbias,&
         nchain=nchain,BETA=BETA,source=source,group=group,&
         nunits_tot=nunits_tot,nvolmoves=nvolmoves,dat_file=dat_file,&
         rmvolume=rmvolume,nunits=nunits,pressure=pressure,cl=cl,&
         cutoff=cutoff,fft_lib=fft_lib,discrete_step=discrete_step,&
         inner_cutoff=inner_cutoff,ldiscrete=ldiscrete)

! nullify some pointers
      NULLIFY(particles_old,particles_test,cell_old,oldsys,&
         testsys,test_env)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)
      ALLOCATE (r_mol(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_mol",3*nunits_tot*dp_size)
      ALLOCATE (loverlap(1:nvolmoves),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "loverlap",nvolmoves)

      IF (ionode) THEN
! record the attempt
         moves%volume%attempts= &
            moves%volume%attempts+1
         move_updates%volume%attempts=&
            move_updates%volume%attempts+1
      ENDIF

! now let's grab the cell length and particle positions
      CALL force_env_get(force_env,subsys=oldsys,cell=cell_old,&
                            error=error)
      CALL get_cell(cell_old,abc=abc)
      CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, &
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
      IF (ionode) THEN

! find the old classical energy
         CALL mc_classical_energy(force_env,&
            classical_energy_original,nchain,1,nchain,&
            cutoff,ldum,cell_length=old_cell_length(1))
         IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
         'Volume move found an overlap in the old remove config')   
         classical_energy_test=classical_energy_original

! save the old coordiantes
         DO iatom=1,nunits_tot
            r(1:3,iatom)=particles_old%els(iatom)%r(1:3)
         ENDDO

! fold the coordinates into the central box
         CALL mc_coordinate_fold(r(:,:),nunits_tot,&
             test_cell_length(1))

! make the new dat file
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
      IF (ionode) THEN

! call a random number to figure out how far we're moving
         CALL RANDOM_NUMBER(rand)

! find the test cell lenghts for the discrete volume move
         IF(ldiscrete) THEN
            IF(rand .LT. 0.5d0) THEN
               lincrease=.true.
            ELSE
               lincrease=.false.
            ENDIF

            new_cell_length(1:3)=test_cell_length(1:3)

! if we're increasing the volume, we need to find a side we can increase
            IF(lincrease) THEN
               DO
                  CALL RANDOM_NUMBER(rand)
                  iside_change=CEILING(3.0d0*rand)
                  IF(discrete_array(iside_change,1) .EQ. 1) THEN
                     new_cell_length(iside_change)=&
                        new_cell_length(iside_change)+discrete_step
                     EXIT
                  ENDIF
               ENDDO
            ELSE
               DO
                  CALL RANDOM_NUMBER(rand)
                  iside_change=CEILING(3.0d0*rand)
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

         DO iatom=1,nunits_tot
            r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
         ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass, using the minimum image (not all molecules are in
! the central box)

! save the particle coords
         DO iatom=1,nunits_tot
            r_mol(1:3,iatom)=particles_test%els(iatom)%r(1:3)
         ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass
         DO iatom=1,nunits_tot,nunits

! now find the center of mass
            CALL get_center_of_mass(r_mol(:,iatom:iatom+2),nunits,&
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
         CALL mc_classical_energy(test_env,rdum,&
               nchain,1,nchain,inner_cutoff,loverlap(1),&
               cell_length=new_cell_length(1))

         test_cell_length(1)=new_cell_length(1)
         test_cell_length(2)=new_cell_length(2)
         test_cell_length(3)=new_cell_length(3)

! make the new dat file
         DO iatom=1,nunits_tot
            r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
         ENDDO

         CALL mc_make_dat_file(r(:,:),nunits_tot,&
            test_cell_length(1:3)*angstrom,&
            dat_file,force_env,TRIM(fft_lib))

! back to parallel
      ENDIF

! realease the environment
      CALL force_env_release(test_env,error=error)

! figure out if we have overlap problems
      CALL mp_bcast(loverlap(1),source,group)
      IF(loverlap(1)) THEN
! deallocate some stuff
         DEALLOCATE(r,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r")
         DEALLOCATE(mass,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "mass")
         DEALLOCATE(r_mol,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r_mol")
         DEALLOCATE(loverlap,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "loverlap")
            
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
      CALL force_env_calc_energy_force(test_env,calc_force=.FALSE.)
      test_env%box_changed=.false.
      CALL force_env_get(test_env,&
         potential_energy=new_energy,error=error)
           
      IF (ionode) THEN
! accept or reject the move
! to prevent overflows
	 energy_term=new_energy-old_energy
         volume_term=-REAL(nchain,dp)/BETA*&
            LOG(test_cell_length(1)*test_cell_length(2)*test_cell_length(3)/&
            (old_cell_length(1)*old_cell_length(2)*old_cell_length(3)))
         pressure_term=pressure*vol_dis
         value=-BETA*(energy_term+volume_term+pressure_term)
         IF(value .GT. 500.0d0) THEN
            w=10.0d0
         ELSEIF(value .LT. -500.0d0) THEN
            w=0.0d0
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
            CALL RANDOM_NUMBER(rand)
         ENDIF

      ENDIF


      CALL mp_bcast(rand,source,group)
      CALL mp_bcast(w,source,group)

      IF (rand .LT. w) THEN

! write cell length, volume, density, and trial displacement to a file
         IF(ionode) THEN

!            WRITE(cl,*) nnstep,test_cell_length(1:3)*angstrom,&
!               (new_energy-old_energy)/3.16682968E-6_dp,classical_energy_diff/3.16682968E-6_dp,&
!               -REAL(nchain,dp)/&
!               BETA*3.0E0_dp*LOG(test_cell_length(1)/&
!               old_cell_length)/3.16682968E-6_dp, &
!               pressure*(test_cell_length(1)**3-old_cell_length**3)/3.16682968E-6_dp,&
!               vol_dis*(angstrom)**3,rand,new_energy/3.16682968E-6_dp,old_energy/3.16682968E-6_dp,w
               
! accept the move
            moves%volume%successes=&
               moves%volume%successes+1
            move_updates%volume%successes=&
               move_updates%volume%successes+1

! update energies
            energy_check=energy_check+(new_energy-old_energy)
            old_energy=new_energy

! and classical energy
            IF ( lbias) moves%bias_energy=&
                  classical_energy_test

! and the coordiantes
            CALL force_env_get(test_env,subsys=testsys,error=error)
            CALL cp_subsys_get(testsys(1)%subsys, &
               particles=particles_test, error=error)

            DO iatom=1,nunits_tot
               r_old(1:3,iatom)=particles_test%els(iatom)%r(1:3)
            ENDDO

         ENDIF

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
         IF (ionode) THEN
!            WRITE(cl,*) nnstep,old_cell_length*angstrom,&
!               (new_energy-old_energy)/3.16682968E-6_dp,classical_energy_diff/3.16682968E-6_dp,&
!               -REAL(nchain,dp)/&
!               BETA*3.0E0_dp*LOG(test_cell_length(1)/&
!               old_cell_length)/3.16682968E-6_dp, &
!               pressure*(test_cell_length(1)*test_cell_length(2)*&
!               test_cell_length(3)-old_cell_length**3)/3.16682968E-6_dp,&
!               vol_dis*(angstrom)**3,rand,new_energy/3.16682968E-6_dp,old_energy/3.16682968E-6_dp,w
               

         ENDIF


! release some memory
         CALL force_env_release(test_env,error=error)

      ENDIF


! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")
      DEALLOCATE(r_mol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_mol")
      DEALLOCATE(loverlap,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "loverlap")
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
!!     alters the length of an OH bond for the given water molecule, keeping
!!     the oxygen in the same position
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!     - r_old: the initial coordinates of all molecules in the system
!!     - r_new: the new coordinates of all molecules in the system
!!     - dis_length: the ratio of the new bond length to the old bond length,
!!                   used in the acceptance rule
!!              
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE change_bond_length ( r_old,r_new,mc_par,mol,dis_length)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp),DIMENSION(:,:), &
      INTENT(OUT)                            :: r_new
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r_old
    INTEGER, INTENT(IN)                      :: mol
    REAL(KIND=dp), INTENT(OUT)               :: dis_length

    CHARACTER(LEN=*), PARAMETER :: routineN = "change_bond_length"

    INTEGER                                  :: handle, hydrogen, i, ipart,&
      nunits_tot
    REAL(KIND=dp)                            :: old_length, rand,rmbond
    REAL(KIND=dp), DIMENSION(1:3)            :: OH

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,nunits_tot=nunits_tot,rmbond=rmbond)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      CALL RANDOM_NUMBER(rand)
      IF ( rand .LT. 0.50E0_dp) THEN
         hydrogen=1
      ELSE
         hydrogen=2
      ENDIF

! choose a displacement
      CALL RANDOM_NUMBER(rand)
      dis_length=rmbond*2.0E0_dp*(rand-0.5E0_dp)

! find the bond distance and change it
      DO i=1,3
         OH(i)=r_new(i,mol+hydrogen)-r_new(i,mol)
      ENDDO
      old_length=SQRT(DOT_PRODUCT(OH,OH))
      DO i=1,3
         r_new(i,mol+hydrogen)=r_new(i,mol)+&
              OH(i)*(old_length+dis_length)/old_length
      ENDDO

! correct the value of dis_length for the acceptance rule
      dis_length=(old_length+dis_length)/old_length

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
!!     alters the magnitude of the HOH angle for a water molecule, keeping
!!     the angle bisector in the same position
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - mol: the number of the molecule's oxygen atom, assuming the hydrogens
!!            are the following two
!!     - r_old: the initial coordinates of all molecules in the system
!!     - r_new: the new coordinates of all molecules in the system
!!              
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE change_bond_angle ( r_old,r_new,mc_par,mol)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: r_new
    REAL(KIND=dp), DIMENSION(:,:), &
      INTENT(IN)                             :: r_old
    INTEGER, INTENT(IN)                      :: mol

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "change_bond_angle"

    INTEGER                                  :: handle, i, ipart,nunits_tot
    REAL(KIND=dp)                            :: bis_length, dis_angle, &
                                                new_angle, old_angle, &
                                                old_length1, old_length2,rand,&
                                                rmangle
    REAL(KIND=dp), DIMENSION(1:3)            :: bisector, bisector1, &
                                                bisector2, OH1, OH2, &
                                                OH_between, OH_between1, &
                                                OH_between2

!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from mc_par
      CALL get_mc_par(mc_par,nunits_tot=nunits_tot,rmangle=rmangle)

! copy the incoming coordinates so we can change them
      DO ipart=1,nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! choose a displacement
      CALL RANDOM_NUMBER(rand)
      dis_angle=rmangle*2.0E0_dp*(rand-0.5E0_dp)
!      dis_angle=0.017453E0_dp   ! always increase by one degree for testing

! find the bond vectors
      DO i=1,3
         OH1(i)=r_new(i,mol+1)-r_new(i,mol)
         OH2(i)=r_new(i,mol+2)-r_new(i,mol)
      ENDDO
      old_length1=SQRT(DOT_PRODUCT(OH1,OH1))
      old_length2=SQRT(DOT_PRODUCT(OH2,OH2))
      OH_between(1:3)=OH1(1:3)/old_length1-OH2(1:3)/old_length2
      old_angle=ACOS(DOT_PRODUCT(OH1,OH2)/(old_length1*old_length2))
      new_angle=old_angle+dis_angle

      DO i=1,3
         bisector(i)=OH1(i)/old_length1+& ! not yet normalized
                 OH2(i)/old_length2
      ENDDO
      bis_length=SQRT(DOT_PRODUCT(bisector,bisector))

! now scale the bisector and the vector from H1 to H2
      OH_between1(1:3) = OH_between(1:3)/SQRT(DOT_PRODUCT(OH_between,&
             OH_between))*sin(new_angle/2.0E0_dp)*old_length1
      OH_between2(1:3) = -OH_between(1:3)/SQRT(DOT_PRODUCT(OH_between,&
             OH_between))*sin(new_angle/2.0E0_dp)*old_length2
      bisector1(1:3)=bisector(1:3)/bis_length*&
             cos(new_angle/2.0E0_dp)*old_length1
      bisector2(1:3)=bisector(1:3)/bis_length*&
             cos(new_angle/2.0E0_dp)*old_length2

      r_new(:,:)=r_old(:,:)
      DO i=1,3
         r_new(i,mol+1)=r_old(i,mol)+&
              bisector1(i)+OH_between1(i)
         r_new(i,mol+2)=r_old(i,mol)+&
              bisector2(i)+OH_between2(i)
      ENDDO
      
! end the timing
      CALL timestop(0.0E0_dp,handle)

  END SUBROUTINE change_bond_angle

!*********************************************************************

END MODULE mc_moves

!**********************************************************************

