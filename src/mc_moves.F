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

MODULE mc_moves
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp,&
                                             dp_size
  USE mc_coordinates,                  ONLY: get_center_of_mass,&
                                             mc_classical_energy,&
                                             mc_coordinate_fold
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_move_control,                 ONLY: move_q_reinit,&
                                             q_move_accept
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_moves_type,&
                                             mc_simulation_parameters_type
  USE message_passing,                 ONLY: mp_bcast
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom
  USE qs_main,                         ONLY: quickstep_create_force_env
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
            mc_molecule_rotation,mc_volume_move,mc_Quickstep_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_moves.F/mc_Quickstep_move
!!
!!   NAME
!!     mc_Quickstep_move
!!
!!   SYNOPSIS
!!     Subroutine mc_Quickstep_move(mc_par,force_env,moves,move_updates,
!!                                  energy_check,r_old,nnstep,old_energy,
!!                                  error)
!!           Type(mc_simulation_parameters_type),Pointer
!!                                                 :: mc_par
!!           Type(force_env_type), Pointer         :: force_env
!!           Type(mc_moves_type), Pointer          :: moves, move_updates
!!           Real(Kind=dp), Intent(INOUT)          :: energy_check
!!           Real(Kind=dp),Dimension(1:3, 1:nunits_tot),
!!                Intent(INOUT)                    :: r_old
!!           Integer,Intent(IN)                    :: nnstep
!!           Real(Kind=dp), Intent(INOUT)          :: old_energy
!!           Type(cp_error_type), Intent(INOUT),Optional
!!                                                 :: error
!!     End Subroutine mc_Quickstep_move
!!
!!   FUNCTION
!!     calls the Quickstep energy routines to calculate the energy of
!!     the current configuration and determine whether we accept or reject
!!     all the moves since the last time a Quickstep calculation was made
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment used in the move
!!     - moves: the structure that keeps track of how many moves have been
!!              accepted/rejected
!!     - move_updates: the structure that keeps track of how many moves have
!!              been accepted/rejected since the last time the displacements
!!              were updated
!!     - energy_check: the running total of how much the energy has changed
!!                     since the initial configuration
!!     - r_old: the coordinates of the last accepted move involving a
!!              Quickstep calculation
!!     - nnstep: how many times the energy averages have been taken
!!     - old_energy: the energy of the last accepted move involving a
!!                   Quickstep calculation
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_Quickstep_move (mc_par,force_env,moves,&
                        move_updates,energy_check,r_old,&
                        nnstep,old_energy,&
                        error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_Quickstep_move"

    INTEGER                                  :: handle, iparticle, nchain,&
      diff,source,group,nunits_tot,iunit
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                classical_energy_old, &
                                                new_energy, rand, w,&
                                                BETA,cutoff
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: lclassical,ionode,&
      lclassical_run,ldum

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)

! get a bunch of data from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,lclassical=lclassical,&
         nchain=nchain,BETA=BETA,diff=diff,source=source,group=group,&
         nunits_tot=nunits_tot,lclassical_run=lclassical_run,&
         cutoff=cutoff)

! in serial
      IF (ionode) THEN

! record the attempt
         moves%Quickstep%attempts=moves%Quickstep%attempts+1
! back to parallel
      ENDIF

! calculate the new energy of the system
      IF( .NOT. lclassical_run) THEN

         CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.)

         CALL force_env_get(force_env,potential_energy=new_energy,&
         error=error)
      ENDIF

! now grab the coords
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

! accept or reject the move based on Metropolis or the Iftimie rule
      IF (ionode) THEN
         
! write them out in case something bad happens
         WRITE(diff,*) nnstep
         WRITE(diff,*) nunits_tot
         DO iunit=1,nunits_tot,3
            WRITE(diff,'(A,3X,3(F18.12,3X))') '  O  ',&
               particles%els(iunit)%r(1:3)*angstrom
            WRITE(diff,'(A,3X,3(F18.12,3X))') '  H  ',&
               particles%els(iunit+1)%r(1:3)*angstrom
            WRITE(diff,'(A,3X,3(F18.12,3X))') '  H  ',&
               particles%els(iunit+2)%r(1:3)*angstrom
         ENDDO

         IF (lclassical) THEN
            CALL mc_classical_energy(force_env,&
               classical_energy_new,nchain,1,nchain,cutoff,ldum)
            IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
            'Quickstep move found an overlap in the old config') 
            IF (lclassical_run) new_energy=classical_energy_new

            classical_energy_old=moves%classical_energy
            w=dexp(-BETA*((new_energy-classical_energy_new)&
                   -(old_energy-classical_energy_old)))

            WRITE(diff,*) nnstep,new_energy-old_energy,&
                 classical_energy_new-classical_energy_old
                 
!!!!!!!!!!!!!!!!!!!!!!!
!            w=1.0D0
!            w=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!
         ELSE
            w=dexp(-BETA*(new_energy-old_energy))
         ENDIF

         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

! done with serial
      ENDIF

! broadcast these numbers, to make sure we're on the same page
      CALL mp_bcast(rand,source,group)      
      CALL mp_bcast(w,source,group)      

      IF (rand .LT. w) THEN

! go back to serial
         IF(ionode) THEN
! accept the move
            moves%Quickstep%successes=moves%Quickstep%successes+1

! we need to record all accepted moves since last Quickstep calculation
            CALL q_move_accept(moves)

! reset the counters
            CALL move_q_reinit(moves)

! update energies
            energy_check=energy_check+(new_energy-old_energy)
            old_energy=new_energy

! update coordinates
            DO iparticle=1,nunits_tot
               r_old(1:3,iparticle)=particles%els(iparticle)%r(1:3)
            ENDDO
         
! if we're biasing classically, then we update the classical energy
            IF ( lclassical) THEN
               moves%classical_energy=classical_energy_new
            ENDIF
! back to parallel
         ENDIF

      ELSE

! in serial
         IF (ionode) THEN
 ! reject the move
            CALL move_q_reinit(moves)
            CALL move_q_reinit(move_updates)

! back to parallel
         ENDIF

         IF ( .NOT. ionode) r_old(:,:) = 0.0D0

! coodinates changed, so we need to broadcast those
         CALL mp_bcast(r_old,source,group)         
          
         DO iparticle=1,nunits_tot
            particles%els(iparticle)%r(1:3)=r_old(1:3,iparticle)
         ENDDO

      ENDIF

! make sure the coordinates are transferred
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_Quickstep_move

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

  SUBROUTINE mc_atom_translation ( mc_par,force_env, moves,&
                        move_updates,nnstep,counter,mol,classical_energy,&
                        error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_atom_translation"

    INTEGER                                  :: handle, i, imol, ipart,&
                                                istat,nchain,nunits_tot
    LOGICAL                                  :: bond_flag,lclassical,&
      loverlap
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                classical_energy_old, &
                                                dis_length, rand, w, BETA,&
                                                cutoff
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r, r_new
    TYPE(particle_list_type), POINTER        :: particles

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lclassical=lclassical,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,cutoff=cutoff)

! nullify some pointers
      NULLIFY(particles,subsys)
      bond_flag=.FALSE.

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)
      ALLOCATE (r_new(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new",3*nunits_tot*dp_size)

! find the old classical energy of the moving mol if we're biasing classically
      imol = (mol+2)/3
      IF (lclassical ) THEN
         CALL mc_classical_energy(force_env,classical_energy_old,&
            nchain,imol,imol,cutoff,loverlap)
            IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
            'Atom move found an overlap in the old config') 

      ENDIF

! let's grab the coordinates
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords 
      DO ipart=1,nunits_tot
         r(1:3,ipart)=particles%els(ipart)%r(1:3)
      ENDDO

! choose if we're changing a bond length or an angle
      CALL RANDOM_NUMBER(rand)
!      IF (rand .lt. -0.66666667d0) THEN
      IF (rand .LT. 0.66666667d0) THEN

! change a bond length
         bond_flag=.TRUE.

! record the attempt
         moves%bond%attempts=moves%bond%attempts+1
         move_updates%bond%attempts=move_updates%bond%attempts+1
         IF ( .NOT. lclassical ) THEN
            moves%bond%qsuccesses=moves%bond%qsuccesses+1
            move_updates%bond%qsuccesses=&
                 move_updates%bond%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_length(r,r_new,mc_par,mol,dis_length)

      ELSE

! record the attempt
         moves%angle%attempts=moves%angle%attempts+1
         move_updates%angle%attempts=move_updates%angle%attempts+1
         IF ( .NOT. lclassical ) THEN
            moves%angle%qsuccesses=moves%angle%qsuccesses+1
            move_updates%angle%qsuccesses=&
                 move_updates%angle%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_angle(r,r_new,mc_par,mol)
         dis_length=1.0d0
      ENDIF

! set the coordinates
      DO ipart=1,nunits_tot
         particles%els(ipart)%r(1:3)=r_new(1:3,ipart)
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)


! if we're biasing classical, check for acceptance
      IF(lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,classical_energy_new,&
            nchain,imol,imol,cutoff,loverlap)

! accept or reject the move based on the Metropolis rule with a 
! correction factor for the change in phase space...dis_length is
! made unitless in change_bond_length
         w=dexp(-BETA*(classical_energy_new-&
                    classical_energy_old))*dis_length**2

         IF(loverlap) w=0.0d0

         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            IF (bond_flag) THEN
               moves%bond%qsuccesses=moves%bond%qsuccesses+1
               move_updates%bond%successes=&
                  move_updates%bond%successes+1
            ELSE
               moves%angle%qsuccesses=moves%angle%qsuccesses+1
               move_updates%angle%successes=&
                  move_updates%angle%successes+1
            ENDIF
            classical_energy=classical_energy+classical_energy_new-&
                                 classical_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")
      DEALLOCATE(r_new,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r_new")

! end the timing
  CALL timestop(0.0d0,handle)

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

 SUBROUTINE mc_molecule_translation(  mc_par,force_env, moves,&
                      move_updates,nnstep,counter,mol,classical_energy,&
                        error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_translation"

    INTEGER                                  :: handle, i, imol, ipart, &
                                                iparticle,istat,nchain,&
                                                nunits_tot
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                classical_energy_old, &
                                                dis_mol, rand, w, BETA,&
                                                rmtrans,cutoff
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: lclassical,loverlap

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lclassical=lclassical,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,rmtrans=rmtrans,cutoff=cutoff)

! nullify some pointers
      NULLIFY(particles,subsys)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)

! find the old classical energy of the moving mol if we're biasing classically
      imol = (mol+2)/3
      IF (lclassical ) THEN
        CALL mc_classical_energy(force_env,classical_energy_old,&
           nchain,imol,imol,cutoff,loverlap)
        IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
           'Molecule translation found an overlap in the old config') 
      ENDIF

! record the attempt
      moves%trans%attempts=moves%trans%attempts+1
      move_updates%trans%attempts=move_updates%trans%attempts+1
      IF ( .NOT. lclassical ) THEN
         moves%trans%qsuccesses=moves%trans%qsuccesses+1
         move_updates%trans%qsuccesses=move_updates%trans%qsuccesses+1
      ENDIF

! move one molecule in the system

! call a random number to figure out which direction we're moving
      CALL RANDOM_NUMBER(rand)
      i=AINT(3*rand)+1

! call a random number to figure out how far we're moving
      CALL RANDOM_NUMBER(rand)
      dis_mol=rmtrans*(rand-0.5d0)*2.0d0

! now let's grab the coordinates and do the move
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords if we're biasing classically
      IF (lclassical) THEN
         DO ipart=1,nunits_tot
            r(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO
      ENDIF

! do the move

      DO iparticle=mol,mol+2
         particles%els(iparticle)%r(i)=&
             particles%els(iparticle)%r(i)+dis_mol
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! if we're biasing classical, check for acceptance
      IF(lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,classical_energy_new,&
            nchain,imol,imol,cutoff,loverlap)

! accept or reject the move based on the Metropolis rule
         w=dexp(-BETA*(classical_energy_new-&
                                 classical_energy_old))

         IF(loverlap) w=0.0d0

         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%trans%qsuccesses=moves%trans%qsuccesses+1
            move_updates%trans%successes=&
                      move_updates%trans%successes+1
            moves%qtrans_dis=moves%qtrans_dis+ABS(dis_mol)
            classical_energy=classical_energy+classical_energy_new-&
                                 classical_energy_old
         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF


! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")

! end the timing
  CALL timestop(0.0d0,handle)

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

 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, moves,&
                      move_updates,nnstep,counter,mol,classical_energy,&
                        error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter, mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_rotation"
    REAL(KIND=dp), DIMENSION(1:3), PARAMETER :: &
      mass = (/ 15.999d0,1.008d0,1.008d0 /)
    REAL(KIND=dp), PARAMETER                 :: masstot = 18.01D0

    INTEGER                                  :: dir, handle, ii, imol, ipart, &
                                                iunit,istat,nchain,nunits_tot,&
                                                nunits
    LOGICAL                                  :: lx, ly,lclassical,loverlap
    REAL(KIND=dp) :: classical_energy_new, classical_energy_old, cosdg, &
      dgamma, nxcm, nycm, nzcm, rand, rx, rxnew, ry, rynew, rz, rznew, sindg, &
      w,BETA,rmrot,cutoff
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,lclassical=lclassical,nchain=nchain,&
         BETA=BETA,nunits_tot=nunits_tot,nunits=nunits,rmrot=rmrot,&
         cutoff=cutoff)

! nullify some pointers
      NULLIFY(particles,subsys)

! do some allocation
      ALLOCATE (r(1:3,1:nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r",3*nunits_tot*dp_size)

! initialize some stuff
      lx =.FALSE.
      ly =.FALSE.

! find the old classical energy of the moving mol if we're biasing classically
      imol = (mol+2)/3
      IF (lclassical ) THEN
         CALL mc_classical_energy(force_env,classical_energy_old,&
            nchain,imol,imol,cutoff,loverlap)
         IF(loverlap) CALL stop_program(routineN,module_name,__LINE__,&
            'Rotation found an overlap in the old config') 
      ENDIF
!      STOP

! record the attempt
      moves%rot%attempts=moves%rot%attempts+1
      move_updates%rot%attempts=move_updates%rot%attempts+1
      IF ( .NOT. lclassical ) THEN
         moves%rot%qsuccesses=moves%rot%qsuccesses+1
         move_updates%rot%qsuccesses=move_updates%rot%qsuccesses+1
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

! now let's grab the coordinates and do the move
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords if we're biasing classically
      IF (lclassical) THEN
         DO ipart=1,nunits_tot
            r(1:3,ipart)=particles%els(ipart)%r(1:3)
         ENDDO
      ENDIF

! Determine new center of mass for chain i by finding the sum
! of m*r for each unit, then dividing by the total mass of the chain
      nxcm = 0.0d0
      nycm = 0.0d0
      nzcm = 0.0d0
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
      dgamma=rmrot*(rand-0.5d0)*2.0d0

! *** set up the rotation marix ***
 
      cosdg = DCOS( dgamma )
      sindg = DSIN( dgamma )
 
      IF (lx) THEN
 
! ***    ROTATE UNITS OF I AROUND X-AXIS ***
 
         DO  iunit = mol,mol+2
            ry = particles%els(iunit)%r(2) - nycm
            rz = particles%els(iunit)%r(3) - nzcm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            rynew = cosdg * ry + sindg * rz
!            rznew = cosdg * rz - sindg * ry
            rynew = cosdg * ry - sindg * rz
            rznew = cosdg * rz + sindg * ry
               
            particles%els(iunit)%r(2) = rynew + nycm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO
      ELSEIF (ly) THEN 
 
! ***    ROTATE UNITS OF I AROUND y-AXIS ***
 
         DO  iunit = mol,mol+2
            rx = particles%els(iunit)%r(1) - nxcm
            rz = particles%els(iunit)%r(3) - nzcm
            rxnew = cosdg * rx + sindg * rz
            rznew = cosdg * rz - sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO

      ELSE
 
! ***    ROTATE UNITS OF I AROUND z-AXIS ***
 
         DO  iunit = mol,mol+2
            rx = particles%els(iunit)%r(1) - nxcm
            ry = particles%els(iunit)%r(2) - nycm
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            rxnew = cosdg * rx + sindg * ry
!            rynew = cosdg * ry - sindg * rx
            rxnew = cosdg * rx - sindg * ry
            rynew = cosdg * ry + sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(2) = rynew + nycm

         ENDDO

      ENDIF 
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! if we're biasing classical, check for acceptance
      IF(lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,classical_energy_new,&
            nchain,imol,imol,cutoff,loverlap)

! accept or reject the move based on the Metropolis rule
         w=dexp(-BETA*(classical_energy_new-&
                                 classical_energy_old))

         IF(loverlap) w=0.0d0

         IF ( w .GE. 1.0d0 ) THEN
            w=1.0d0
            rand=0.0D0
         ELSE
            CALL RANDOM_NUMBER(rand)
         ENDIF

         IF (rand .LT. w) THEN

! accept the move
            moves%rot%qsuccesses=moves%rot%qsuccesses+1
            move_updates%rot%successes=move_updates%rot%successes+1
            classical_energy=classical_energy+classical_energy_new-&
                                 classical_energy_old

         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO ipart=1,nunits_tot
               particles%els(ipart)%r(1:3)=r(1:3,ipart)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! deallocate some stuff
      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")

! end the timing
  CALL timestop(0.0d0,handle)

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
                        energy_check,r_old,error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp),DIMENSION(:, :), &
      INTENT(INOUT)                          :: r_old
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: &
      routineN = "mc_volume_move"
    REAL(KIND=dp), DIMENSION(1:3), PARAMETER :: &
      mass = (/ 15.999d0,1.008d0,1.008d0 /)
    REAL(KIND=dp), PARAMETER                 :: mass_tot = 18.01d0

    INTEGER                                  :: handle, i, iatom, &
                                                ivolmove,istat,nchain,&
                                                source,group,nunits_tot,&
                                                nvolmoves,nunits,cl
    REAL(KIND=dp) :: classical_energy_diff, classical_energy_new, &
      classical_energy_original, classical_energy_test, new_cell_length, &
      new_energy, old_cell_length, rand, test_cell_length, vol_dis, w,&
      BETA,rmvolume,pressure,cutoff
    REAL(KIND=dp), DIMENSION(1:3)            :: center_of_mass, &
                                                center_of_mass_new, diff,&
                                                abc
    TYPE(force_env_type), POINTER            :: test_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: oldsys, testsys
    TYPE(cell_type), POINTER                 :: cell_old
    REAL(KIND=dp), DIMENSION(:, :),&
      ALLOCATABLE                            :: r,r_mol
    TYPE(particle_list_type), POINTER        :: particles_old, particles_test
    LOGICAL                                  :: ldum,ionode,lclassical,&
      lclassical_run
    LOGICAL,DIMENSION(:),ALLOCATABLE         :: loverlap
    CHARACTER(LEN=40)                        :: dat_file
    CHARACTER(LEN=200)                       :: fft_lib
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,lclassical=lclassical,&
         nchain=nchain,BETA=BETA,source=source,group=group,&
         nunits_tot=nunits_tot,nvolmoves=nvolmoves,dat_file=dat_file,&
         rmvolume=rmvolume,nunits=nunits,pressure=pressure,cl=cl,&
         lclassical_run=lclassical_run,cutoff=cutoff,fft_lib=fft_lib)

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
         IF(nvolmoves .EQ. 1) THEN
            move_updates%classical_volume%attempts=&
               move_updates%classical_volume%attempts+1
            moves%classical_volume%attempts=&
               moves%classical_volume%attempts+1
         ENDIF
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
             particles=particles_old, error=error)

! find the old cell length
      old_cell_length=abc(1)
      test_cell_length=abc(1)
      IF (ionode) THEN

! find the old classical energy
         CALL mc_classical_energy(force_env,&
            classical_energy_original,nchain,1,nchain,&
            cutoff,ldum,cell_length=old_cell_length)
         IF(ldum) CALL stop_program(routineN,module_name,__LINE__,&
         'Volume move found an overlap in the old remove config')   
         classical_energy_test=classical_energy_original

! save the old coordiantes
         DO iatom=1,nunits_tot
            r(1:3,iatom)=particles_old%els(iatom)%r(1:3)
         ENDDO

! fold the coordinates into the central box
         CALL mc_coordinate_fold(r(:,:),nunits_tot,&
             test_cell_length)

! make the new dat file
         CALL mc_make_dat_file(r(:,:),nunits_tot,&
            test_cell_length*angstrom,dat_file,force_env,TRIM(fft_lib))
      ENDIF

! make a test environment that we can alter to calculate the energy
      globenv%input_file_name=dat_file
      CALL quickstep_create_force_env(test_env, globenv, error=error)

! grab the coords so we can change them
      CALL force_env_get(test_env,subsys=testsys,error=error)
      CALL cp_subsys_get(testsys(1)%subsys, &
             particles=particles_test, error=error)

! now do some classical moves to find a good move
      IF (ionode) THEN

! now do a series of volume moves, accepting or rejecting each one based on
!    the fake potential
        DO ivolmove=1,nvolmoves

! call a random number to figure out how far we're moving
            CALL RANDOM_NUMBER(rand)
!!!!!!!!!!!!
            vol_dis=rmvolume*(rand-0.5d0)*2.0d0
!            vol_dis=5.0d0/angstrom**3

            new_cell_length=(vol_dis+&
               test_cell_length**3)**(1.0d0/3.0d0)
            DO iatom=1,nunits_tot
               r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
            ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass, using the minimum image (not all molecules are in
! the central box)

! save the particle coords
            DO i=1,nunits_tot
               r_mol(1:3,i)=particles_test%els(i)%r(1:3)
            ENDDO

! now we need to scale the coordinates of all the molecules by the
! center of mass
            DO iatom=1,nunits_tot,nunits

! now find the center of mass
               CALL get_center_of_mass(r_mol(:,iatom:iatom+2),nunits,&
                  center_of_mass(:),mass(:))

! scale the center of mass and determine the vector that points from the
!    old COM to the new one
               center_of_mass_new(1:3)=center_of_mass(1:3)*&
                  new_cell_length/test_cell_length

               DO i=1,3
                  diff(i)=center_of_mass_new(i)-center_of_mass(i)
! now change the particle positions
                  particles_test%els(iatom)%r(i)=&
                     particles_test%els(iatom)%r(i)+diff(i)
                  particles_test%els(iatom+1)%r(i)=&
                     particles_test%els(iatom+1)%r(i)+diff(i)
                  particles_test%els(iatom+2)%r(i)=&
                     particles_test%els(iatom+2)%r(i)+diff(i)
               ENDDO
            ENDDO

! find and store the classical energy
            CALL mc_classical_energy(test_env,classical_energy_new,&
               nchain,1,nchain,cutoff,loverlap(ivolmove),&
               cell_length=new_cell_length)

            IF (nvolmoves .GT. 1) THEN

! determine the acceptance of the move
               w=dexp(-BETA*(classical_energy_new-&
                  classical_energy_test-REAL(nchain,dbl)/&
                  BETA*3.0d0*&
                  DLOG(new_cell_length/test_cell_length) + &
                  pressure*vol_dis))
               
               IF(loverlap(ivolmove)) w=0.0d0
               IF ( w .GE. 1.0d0 ) THEN
                  w=1.0d0
                  rand=0.0D0
               ELSE
                  CALL RANDOM_NUMBER(rand)
               ENDIF
               
               IF (rand .LT. w ) THEN

                  move_updates%classical_volume%attempts=&
                     move_updates%classical_volume%attempts+1
                  move_updates%classical_volume%successes=&
                     move_updates%classical_volume%successes+1
                  moves%classical_volume%attempts=&
                     moves%classical_volume%attempts+1
                  moves%classical_volume%successes=&
                     moves%classical_volume%successes+1
         
! update energy and cell length
                  test_cell_length=new_cell_length
                  classical_energy_test=classical_energy_new

               ELSE

                  move_updates%classical_volume%attempts=&
                     move_updates%classical_volume%attempts+1
                  moves%classical_volume%attempts=&
                     moves%classical_volume%attempts+1

! revert back to old coords                      
                  DO iatom=1,nunits_tot
                     particles_test%els(iatom)%r(1:3)=r(1:3,iatom)
                  ENDDO

               ENDIF
            
            ELSE
               test_cell_length=new_cell_length
               classical_energy_test=classical_energy_new
            ENDIF
            
         ENDDO

! determine the overall energy difference
         IF (nvolmoves .GT. 1 ) THEN
            classical_energy_diff=classical_energy_test-&
               classical_energy_original
         ELSE
            classical_energy_diff=0.0d0
         ENDIF


! make the new dat file
         DO iatom=1,nunits_tot
            r(1:3,iatom)=particles_test%els(iatom)%r(1:3)
         ENDDO

         CALL mc_make_dat_file(r(:,:),nunits_tot,&
            test_cell_length*angstrom,dat_file,force_env,TRIM(fft_lib))

!         loverlap=.FALSE.
!         CALL check_for_overlap(r(:,:),&
!            nunits_tot,test_cell_length,loverlap)

! back to parallel
      ENDIF

! set the coordinates and the energy
      CALL force_env_release(test_env,error=error)

! figure out if we have overlap problems
      CALL mp_bcast(loverlap(1),source,group)
      IF(loverlap(1)) THEN
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
            
         IF(ionode) WRITE(globenv%scr,*) &
            "Volume move rejected due to overlap.",vol_dis
!     end the timing
         CALL timestop(0.0d0,handle)
         RETURN
      ENDIF

      globenv%input_file_name=dat_file
      CALL quickstep_create_force_env(test_env, globenv, error=error)

      IF (.NOT. lclassical_run) THEN

         CALL force_env_calc_energy_force(test_env,calc_force=.FALSE.)
         CALL force_env_get(test_env,&
            potential_energy=new_energy,error=error)
           
      ELSE
         new_energy=classical_energy_test
      ENDIF

      IF (ionode) THEN
! accept or reject the move
         IF (nvolmoves .GT. 1) THEN
            w=dexp(-BETA*((new_energy-old_energy)-&
               classical_energy_diff))
         ELSE

            w=dexp(-BETA*((new_energy-old_energy)-&
               REAL(nchain,dbl)/BETA*3.0d0*&
               DLOG(test_cell_length/old_cell_length) &
               + pressure*vol_dis))
            
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


      CALL mp_bcast(rand,source,group)
      CALL mp_bcast(w,source,group)

      IF (rand .LT. w) THEN

! write cell length, volume, density, and trial displacement to a file
         IF(ionode) THEN

            WRITE(cl,*) nnstep,test_cell_length*angstrom,&
               (new_energy-old_energy)/3.16682968d-6,classical_energy_diff/3.16682968d-6,&
               -REAL(nchain,dbl)/&
               BETA*3.0d0*DLOG(test_cell_length/&
               old_cell_length)/3.16682968d-6, &
               pressure*(test_cell_length**3-old_cell_length**3)/3.16682968d-6,&
               vol_dis*(angstrom)**3,rand,new_energy/3.16682968d-6,old_energy/3.16682968d-6,w
               
! accept the move
            IF(nvolmoves .EQ. 1) THEN
               move_updates%classical_volume%successes=&
                  move_updates%classical_volume%successes+1
               moves%classical_volume%successes=&
                  moves%classical_volume%successes+1
            ENDIF
            moves%volume%successes=&
               moves%volume%successes+1
            move_updates%volume%successes=&
               move_updates%volume%successes+1

! update energies
            energy_check=energy_check+(new_energy-old_energy)
            old_energy=new_energy

! and classical energy
            IF ( lclassical) moves%classical_energy=&
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

      ELSE

! reject the move
! write cell length, volume, density, and trial displacement to a file
         IF (ionode) THEN
            WRITE(cl,*) nnstep,old_cell_length*angstrom,&
               (new_energy-old_energy)/3.16682968d-6,classical_energy_diff/3.16682968d-6,&
               -REAL(nchain,dbl)/&
               BETA*3.0d0*DLOG(test_cell_length/&
               old_cell_length)/3.16682968d-6, &
               pressure*(test_cell_length**3-old_cell_length**3)/3.16682968d-6,&
               vol_dis*(angstrom)**3,rand,new_energy/3.16682968d-6,old_energy/3.16682968d-6,w
               

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

! end the timing
      CALL timestop(0.0d0,handle)

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
      IF ( rand .LT. 0.50d0) THEN
         hydrogen=1
      ELSE
         hydrogen=2
      ENDIF

! choose a displacement
      CALL RANDOM_NUMBER(rand)
      dis_length=rmbond*2.0d0*(rand-0.5d0)

! find the bond distance and change it
      DO i=1,3
         OH(i)=r_new(i,mol+hydrogen)-r_new(i,mol)
      ENDDO
      old_length=DSQRT(DOT_PRODUCT(OH,OH))
      DO i=1,3
         r_new(i,mol+hydrogen)=r_new(i,mol)+&
              OH(i)*(old_length+dis_length)/old_length
      ENDDO

! correct the value of dis_length for the acceptance rule
      dis_length=(old_length+dis_length)/old_length

! end the timing
      CALL timestop(0.0d0,handle)

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
      dis_angle=rmangle*2.0d0*(rand-0.5d0)
!      dis_angle=0.017453D0   ! always increase by one degree for testing

! find the bond vectors
      DO i=1,3
         OH1(i)=r_new(i,mol+1)-r_new(i,mol)
         OH2(i)=r_new(i,mol+2)-r_new(i,mol)
      ENDDO
      old_length1=DSQRT(DOT_PRODUCT(OH1,OH1))
      old_length2=DSQRT(DOT_PRODUCT(OH2,OH2))
      OH_between(1:3)=OH1(1:3)/old_length1-OH2(1:3)/old_length2
      old_angle=DACOS(DOT_PRODUCT(OH1,OH2)/(old_length1*old_length2))
      new_angle=old_angle+dis_angle

      DO i=1,3
         bisector(i)=OH1(i)/old_length1+& ! not yet normalized
                 OH2(i)/old_length2
      ENDDO
      bis_length=DSQRT(DOT_PRODUCT(bisector,bisector))

! now scale the bisector and the vector from H1 to H2
      OH_between1(1:3) = OH_between(1:3)/DSQRT(DOT_PRODUCT(OH_between,&
             OH_between))*dsin(new_angle/2.0d0)*old_length1
      OH_between2(1:3) = -OH_between(1:3)/DSQRT(DOT_PRODUCT(OH_between,&
             OH_between))*dsin(new_angle/2.0d0)*old_length2
      bisector1(1:3)=bisector(1:3)/bis_length*&
             dcos(new_angle/2.0d0)*old_length1
      bisector2(1:3)=bisector(1:3)/bis_length*&
             dcos(new_angle/2.0d0)*old_length2

      r_new(:,:)=r_old(:,:)
      DO i=1,3
         r_new(i,mol+1)=r_old(i,mol)+&
              bisector1(i)+OH_between1(i)
         r_new(i,mol+2)=r_old(i,mol)+&
              bisector2(i)+OH_between2(i)
      ENDDO
      
! end the timing
      CALL timestop(0.0d0,handle)

  END SUBROUTINE change_bond_angle

!*********************************************************************

END MODULE mc_moves

!**********************************************************************

