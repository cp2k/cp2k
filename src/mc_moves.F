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
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_retain,&
                                             force_env_release,&
                                             force_env_type,&
                                             force_env_set_cell
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_misc,                         ONLY: mc_averages_type,&
                                             mc_classical_energy
  USE mc_move_control,                 ONLY: mc_energy_restart,&
                                             mc_move_update,&
                                             mc_moves_type,&
                                             move_q_reinit,&
                                             q_move_accept
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type
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
!!                                  averages,classical_energy_old_out,error)
!!           Type(mc_simulation_parameters_type),Pointer
!!                                                 :: mc_par
!!           Type(force_env_type), Pointer         :: force_env
!!           Type(mc_moves_type), Pointer          :: moves, move_updates
!!           Real(Kind=dp), Intent(INOUT)          :: energy_check
!!           Real(Kind=dp),Dimension(1:3, 1:mc_par%nunits_tot),
!!                Intent(INOUT)                    :: r_old
!!           Integer,Intent(IN)                    :: nnstep
!!           Real(Kind=dp), Intent(INOUT)          :: old_energy
!!           Type(mc_averages_type), Pointer       :: averages
!!           Type(cp_error_type), Intent(INOUT),Optional
!!                                                 :: error
!!           Real(Kind=dp),Intent(inout)           :: classical_energy_old_out
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
!!     - averages: the structure that keeps track of the running averages
!!                 in the system
!!     - classical_energy_old_out: the classical energy of the configuration
!!                                 that leaves the routine
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

  SUBROUTINE mc_Quickstep_move (mc_par,force_env,moves,&
                        move_updates,energy_check,r_old,&
                        nnstep,old_energy,averages,&
                        classical_energy_old_out,&
                        error)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    REAL(KIND=dp), INTENT(INOUT)             :: energy_check
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot), &
      INTENT(INOUT)                          :: r_old
    INTEGER, INTENT(IN)                      :: nnstep
    REAL(KIND=dp), INTENT(INOUT)             :: old_energy
    TYPE(mc_averages_type), POINTER          :: averages
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    REAL(KIND=dp),INTENT(inout)          :: classical_energy_old_out

    CHARACTER(LEN=*), PARAMETER        :: routineN = "mc_Quickstep_move"

    INTEGER                                  :: handle,iparticle, istep
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                classical_energy_old, &
                                                new_energy, rand, w
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! create a new variable to keep averages correct
      istep=nnstep-mc_par%nstart

! nullify some pointers
      NULLIFY(particles,subsys)

! record the attempt
      moves%Quickstep%attempts=moves%Quickstep%attempts+1

! calculate the new energy of the system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this could cause problems on various platforms!!!!!!!!!!!!!!!!
      CALL mc_energy_restart(force_env,new_energy,error=error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! now let's grab the coordinates
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! accept or reject the move based on Metropolis or the Iftimie rule
      IF (mc_par%lclassical) THEN
         CALL mc_classical_energy(force_env,&
                                  classical_energy_new,&
                                  mc_par)
         classical_energy_old=moves%classical_energy
         w=dexp(-mc_par%BETA*((new_energy-classical_energy_new)&
                   -(old_energy-classical_energy_old)))
         IF ( w .GT. 1.0d0 ) w=1.0d0
      ELSE
         w=dexp(-mc_par%BETA*(new_energy-old_energy))
         IF ( w .GT. 1.0d0 ) w=1.0d0
      ENDIF

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! accept the move
         moves%Quickstep%successes=moves%Quickstep%successes+1

! we need to record all accepted moves since last Quickstep calculation
         CALL q_move_accept(moves,move_updates)

! reset the counters
         CALL move_q_reinit(moves)
         CALL move_q_reinit(move_updates)

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy
         averages%ave_energy=averages%ave_energy*&
                     REAL(istep-1,dbl)&
                     /REAL(istep,dbl)+old_energy/&
                     REAL(istep,dbl)
! and heat capacity
         averages%ave_energy_squared=averages%ave_energy_squared*&
               REAL(istep-1,dbl)/&
               REAL(istep,dbl)+old_energy**2/&
               REAL(istep,dbl)

! update coordinates
         DO iparticle=1,mc_par%nunits_tot
            r_old(1:3,iparticle)=particles%els(iparticle)%r(1:3)
         ENDDO
         
! if we're biasing classically, then we update the classical energy
         IF ( mc_par % lclassical) THEN
            CALL mc_classical_energy(force_env,&
                                         moves%classical_energy,&
                                         mc_par)
            classical_energy_old_out=moves%classical_energy            
         ENDIF

      ELSE

! reject the move
         CALL move_q_reinit(moves)
         CALL move_q_reinit(move_updates)

         averages%ave_energy=averages%ave_energy*&
              REAL(istep-1,dbl)/&
              REAL(istep,dbl)+old_energy/&
              REAL(istep,dbl)

! and heat capacity
         averages%ave_energy_squared=averages%ave_energy_squared*&
                 REAL(istep-1,dbl)/&
                 REAL(istep,dbl)+old_energy**2&
                 /REAL(istep,dbl)

! replace coordinates
         DO iparticle=1,mc_par%nunits_tot
            particles%els(iparticle)%r(1:3)=r_old(1:3,iparticle)
         ENDDO

! if we need to restore the classical energy
         IF ( mc_par % lclassical ) THEN
            classical_energy_old_out=moves%classical_energy
         ENDIF

      ENDIF

! make sure the coordinates are transferred
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_Quickstep_move

! *****************************************************************************

  SUBROUTINE mc_atom_translation ( mc_par,force_env, moves,&
                        move_updates,nnstep,counter,mol,&
                        classical_energy_old,error)

!   Purpose: A move changing the conformation of the molecule

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep,counter,mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy_old
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_atom_translation"

    INTEGER                                  :: handle, i, imol
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                dis_ion, rand, w,&
                                                dis_length
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot)    :: r,r_new
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL :: bond_flag


!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)
      bond_flag=.FALSE.

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*&
               mc_par%nchain-mc_par%nstart) == 0 .AND. counter==1) THEN
         CALL mc_move_update(mc_par,move_updates,"bond")
         WRITE(mc_par%rm,*) nnstep,' rmbond = ',&
                                    mc_par%rmbond*0.529177d0,' angstoms'
         CALL mc_move_update(mc_par,move_updates,"angle")
         WRITE(mc_par%rm,*) nnstep,' rmangle = ',&
               mc_par%rmangle/3.14159265358979d0*180.0d0,' degrees'
      ENDIF

! let's grab the coordinates
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords 
      DO imol=1,mc_par%nunits_tot
         r(1:3,imol)=particles%els(imol)%r(1:3)
      ENDDO

! choose if we're changing a bond length or an angle
      CALL random_number(rand)
      IF (rand .lt. 0.66666667d0) THEN

! change a bond length
         bond_flag=.TRUE.

! record the attempt
         moves%bond%attempts=moves%bond%attempts+1
         move_updates%bond%attempts=move_updates%bond%attempts+1
         IF ( .NOT. mc_par % lclassical ) THEN
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
         IF ( .NOT. mc_par % lclassical ) THEN
            moves%angle%qsuccesses=moves%angle%qsuccesses+1
            move_updates%angle%qsuccesses=&
                 move_updates%angle%qsuccesses+1
         ENDIF

! do the move
         CALL change_bond_angle(r,r_new,mc_par,mol)
         dis_length=1.0d0
      ENDIF

! set the coordinates
      DO imol=1,mc_par%nunits_tot
         particles%els(imol)%r(1:3)=r_new(1:3,imol)
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)


! if we're biasing classical, check for acceptance
      IF(mc_par%lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,&
                              classical_energy_new,mc_par)

! accept or reject the move based on the Metropolis rule with a 
! correction factor for the change in phase space...dis_length is
! made unitless in change_bond_length
         w=dexp(-mc_par%BETA*(classical_energy_new-&
                    classical_energy_old))*dis_length**2

         IF (w .GT. 1.0d0) w=1.0d0

         CALL random_number(rand)

         IF (rand .lt. w) THEN

! accept the move
            IF (bond_flag) THEN
               moves%bond%qsuccesses=moves%bond%qsuccesses+1
               move_updates%bond%qsuccesses=&
                  move_updates%bond%qsuccesses+1
            ELSE
               moves%angle%qsuccesses=moves%angle%qsuccesses+1
               move_updates%angle%qsuccesses=&
                  move_updates%angle%qsuccesses+1
            ENDIF
            classical_energy_old=classical_energy_new
         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO imol=1,mc_par%nunits_tot
               particles%els(imol)%r(1:3)=r(1:3,imol)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_atom_translation

!**********************************************************************
 SUBROUTINE mc_molecule_translation(  mc_par,force_env, moves,&
                        move_updates,nnstep,counter,mol,&
                        classical_energy_old,error)

!   Purpose: A move displacing a single molecule

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter,mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy_old
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_translation"

    INTEGER                                  :: handle, i, imol, &
                                                iparticle
    REAL(KIND=dp)                            :: classical_energy_new, &
                                                dis_mol, rand, w
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot)    :: r
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)

! update the molecule translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*&
               mc_par%nchain-mc_par%nstart) == 0 .AND. counter==1) THEN
         CALL mc_move_update(mc_par,move_updates,"trans")
         WRITE(mc_par%rm,*) nnstep,' rmtrans = ',&
                         mc_par%rmtrans*0.529177d0,' angstroms'
      ENDIF

! record the attempt
      moves%trans%attempts=moves%trans%attempts+1
      move_updates%trans%attempts=move_updates%trans%attempts+1
      IF ( .NOT. mc_par % lclassical ) THEN
         moves%trans%qsuccesses=moves%trans%qsuccesses+1
         move_updates%trans%qsuccesses=move_updates%trans%qsuccesses+1
      ENDIF

! move one molecule in the system

! call a random number to figure out which direction we're moving
      CALL random_number(rand)
      i=aint(3*rand)+1

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dis_mol=mc_par%rmtrans*(rand-0.5d0)*2.0d0

! now let's grab the coordinates and do the move
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords if we're biasing classically
      IF (mc_par%lclassical) THEN
         DO imol=1,mc_par%nunits_tot
            r(1:3,imol)=particles%els(imol)%r(1:3)
         ENDDO
      ENDIF

! do the move
      DO iparticle=mol,mol+2
         particles%els(iparticle)%r(i)=&
             particles%els(iparticle)%r(i)+dis_mol
      ENDDO
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! if we're biasing classical, check for acceptance
      IF(mc_par%lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,&
                              classical_energy_new,mc_par)

! accept or reject the move based on the Metropolis rule
         w=dexp(-mc_par%BETA*(classical_energy_new-&
                                 classical_energy_old))
         IF (w .GT. 1.0d0) w=1.0d0

         CALL random_number(rand)

         IF (rand .lt. w) THEN

! accept the move
            moves%trans%qsuccesses=moves%trans%qsuccesses+1
            move_updates%trans%qsuccesses=&
                      move_updates%trans%qsuccesses+1
            moves%qtrans_dis=moves%qtrans_dis+ABS(dis_mol)
            classical_energy_old=classical_energy_new
         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO imol=1,mc_par%nunits_tot
               particles%els(imol)%r(1:3)=r(1:3,imol)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_molecule_translation

!**********************************************************************
!**********************************************************************

 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, moves,&
                        move_updates,nnstep,counter,mol,&
                        classical_energy_old,error)

!   Purpose: A move rotating a single molecule around it's center
!            of mass...only works for water

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_moves_type), POINTER             :: moves, move_updates
    INTEGER, INTENT(IN)                      :: nnstep, counter,mol
    REAL(KIND=dp), INTENT(INOUT)             :: classical_energy_old
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_rotation"
    REAL(KIND=dp), DIMENSION(1:3), PARAMETER :: &
      mass = (/ 15.994d0,1.008d0,1.008d0 /)

    INTEGER                                  :: dir, handle, ii, imol, &
                                                iunit
    LOGICAL                                  :: lx, ly
    REAL(KIND=dp)                            :: classical_energy_new, cosdg, &
                                                dgamma, nxcm, nycm, nzcm, &
                                                rand, rx, rxnew, ry, rynew, &
                                                rz, rznew, sindg, w
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot)    :: r
    TYPE(particle_list_type), POINTER        :: particles

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)

! initialize some stuff
      lx =.false.
      ly =.false.

! update the molecule rotation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*&
               mc_par%nchain-mc_par%nstart) == 0 .AND. counter==1) THEN
         CALL mc_move_update(mc_par,move_updates,"rot")
         WRITE(mc_par%rm,*) nnstep,' rmrot = ',&
                    mc_par%rmrot*0.529177d0,' angstroms'
      ENDIF

! record the attempt
      moves%rot%attempts=moves%rot%attempts+1
      move_updates%rot%attempts=move_updates%rot%attempts+1
      IF ( .NOT. mc_par % lclassical ) THEN
         moves%rot%qsuccesses=moves%rot%qsuccesses+1
         move_updates%rot%qsuccesses=move_updates%rot%qsuccesses+1
      ENDIF

! rotate one molecule in the system

! call a random number to figure out which direction we're moving
      CALL random_number(rand)
      dir=aint(3*rand)+1

      IF (dir .eq. 1) THEN
         lx = .true.
      ELSEIF (dir .eq. 2) THEN
         ly = .true.
      ENDIF

! now let's grab the coordinates and do the move
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

! save the coords if we're biasing classically
      IF (mc_par%lclassical) THEN
         DO imol=1,mc_par%nunits_tot
            r(1:3,imol)=particles%els(imol)%r(1:3)
         ENDDO
      ENDIF

! Determine new center of mass for chain i by finding the sum
! of m*r for each unit, then dividing by the total mass of the chain
      nxcm = 0.0d0
      nycm = 0.0d0
      nzcm = 0.0d0
      do ii = 1, mc_par%nunits
         nxcm = nxcm + particles%els(mol-1+ii)%r(1)* mass(ii)
         nycm = nycm + particles%els(mol-1+ii)%r(2)* mass(ii)
         nzcm = nzcm + particles%els(mol-1+ii)%r(3)* mass(ii)
      enddo
      nxcm = nxcm / 18.0d0
      nycm = nycm / 18.0d0
      nzcm = nzcm / 18.0d0

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dgamma=mc_par%rmrot*(rand-0.5d0)*2.0d0

! *** set up the rotation marix ***
 
      cosdg = DCOS( dgamma )
      sindg = DSIN( dgamma )
 
      IF (lx) THEN
 
! ***    ROTATE UNITS OF I AROUND X-AXIS ***
 
         DO  iunit = mol,mol+2
            ry = particles%els(iunit)%r(2) - nycm
            rz = particles%els(iunit)%r(3) - nzcm
            rynew = cosdg * ry + sindg * rz
            rznew = cosdg * rz - sindg * ry
               
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
            rxnew = cosdg * rx + sindg * ry
            rynew = cosdg * ry - sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(2) = rynew + nycm

         ENDDO

      ENDIF 
      CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

! if we're biasing classical, check for acceptance
      IF(mc_par%lclassical) THEN

! here's where we bias the moves with a classical force field 
         CALL mc_classical_energy(force_env,&
                              classical_energy_new,mc_par)

! accept or reject the move based on the Metropolis rule
         w=dexp(-mc_par%BETA*(classical_energy_new-&
                                 classical_energy_old))
         IF (w .GT. 1.0d0) w=1.0d0

         CALL random_number(rand)

         IF (rand .lt. w) THEN

! accept the move
            moves%rot%qsuccesses=moves%rot%qsuccesses+1
            move_updates%rot%qsuccesses=move_updates%rot%qsuccesses+1
!            moves%rot%successes=moves%rot%successes+1
            classical_energy_old=classical_energy_new
         ELSE

! reject the move
! restore the coordinates
            CALL force_env_get(force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
            DO imol=1,mc_par%nunits_tot
               particles%els(imol)%r(1:3)=r(1:3,imol)
            ENDDO
            CALL cp_subsys_set(subsys(1)%subsys,particles=particles)

         ENDIF

      ENDIF

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_molecule_rotation
!**********************************************************************
  SUBROUTINE mc_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,averages,old_energy,&
                        classical_energy_old,globenv,&
                        energy_check,r_old,error)

!   Purpose: A move changing the size of the box

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    REAL(KIND = dp),INTENT(INOUT) :: old_energy,energy_check,&
                                   classical_energy_old
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(mc_averages_type),POINTER  :: averages
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par
    REAL(KIND=dp),DIMENSION(1:3,1:mc_par%nunits_tot), &
                               INTENT(INOUT) :: r_old

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_volume_move"
    REAL(KIND=dp),DIMENSION(1:3),PARAMETER :: &
                             mass=(/ 15.994d0,1.008d0,1.008d0 /)
    REAL(KIND=dp),PARAMETER :: mass_tot=18.01d0

!   *** Local variables ***
    REAL(KIND = dp) ::                 rand,dis_length,w,new_energy,&
                                       old_length,max_error,vol,&
                                       density_old,new_cell_length,&
                                       old_cell_length,vol_dis
    TYPE(force_env_type), POINTER :: test_env
    TYPE(particle_list_type), POINTER :: particles_old
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: oldsys,&
                                                      newsys
    TYPE(cell_type), POINTER :: cell,cell_new
    REAL(KIND =dp),DIMENSION(1:3,1:3) ::       h_matrix,new_h,&
                                               new_h_inverse
    REAL(KIND=dp),DIMENSION(1:3)  :: diff,center_of_mass,&
                                     center_of_mass_new
    INTEGER :: handle,istep,imol,i
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! rewrite a variable so the averages aren't thrown off
      istep=nnstep-mc_par%nstart

! nullify some pointers
      NULLIFY(test_env,particles_old,particles,oldsys)

! update the volume move max displacement, if necessary
      IF (MOD(nnstep,mc_par%iupvolume*&
               mc_par%nchain-mc_par%nstart) == 0 ) THEN
         CALL mc_move_update(mc_par,move_updates,"volume")
         WRITE(mc_par%rm,*) nnstep,' rmvolume = ',&
                 mc_par%rmvolume*(0.529177d0)**3,' angstroms^3'
      ENDIF

! record the attempt
      moves%volume%attempts=moves%volume%attempts+1
      move_updates%volume%attempts=move_updates%volume%attempts+1

! now let's grab the cell length and change it 
      CALL force_env_get(force_env,cell=cell,&
                         error=error)

! create the old h matrix
      h_matrix(:,:)=cell%hmat(:,:)

! find the old cell length
      old_cell_length=cell%hmat(1,1)

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      vol_dis=mc_par%rmvolume*(rand-0.5d0)*2.0d0
      new_cell_length=(vol_dis+old_cell_length**3)**(1.0d0/3.0d0)

! create the new h matrix
      new_h(:,:)=0.0d0
      new_h(1,1)=new_cell_length
      new_h(2,2)=new_cell_length
      new_h(3,3)=new_cell_length

! create the new environment to change particle coordinates
      CALL quickstep_create_force_env(test_env, globenv, &
                   new_h_matrix=new_h,error=error)
!      moves%force_create=moves%force_create+1
      CALL force_env_get(force_env,subsys=oldsys,&
                            error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
      CALL force_env_get(test_env,subsys=newsys,&
          cell=cell_new,error=error)
      CALL cp_subsys_get(newsys(1)%subsys, &
             particles=particles, error=error)

! now we need to scale the coordinates of all the molecules by the
! center of mass, using the minimum image (not all molecules are in
! the central box)
      DO imol=1,mc_par%nchain*mc_par%nunits,3
         DO i=1,3 ! messy because of minimum image
            center_of_mass(i)=((particles_old%els(imol)%r(i)-&
              old_cell_length*FLOOR( particles_old%els(imol)%r(i)/&
              old_cell_length))*mass(1)+&
              (particles_old%els(imol+1)%r(i)-&
              old_cell_length*FLOOR( particles_old%els(imol+1)%r(i)/&
              old_cell_length))*mass(2)+(particles_old%els(imol+2)%r(i)&
              -old_cell_length*FLOOR( particles_old%els(imol+2)%r(i)/&
              old_cell_length))*mass(3))/mass_tot
         ENDDO
         center_of_mass_new(1:3)=center_of_mass(1:3)*new_cell_length&
                     /old_cell_length
         DO i=1,3
            diff(i)=center_of_mass_new(i)-center_of_mass(i)
         ENDDO
         DO i=1,3
           particles%els(imol)%r(i)=particles_old%els(imol)%r(i)+diff(i)
           particles%els(imol+1)%r(i)=&
                 particles_old%els(imol+1)%r(i)+diff(i)
           particles%els(imol+2)%r(i)=&
                 particles_old%els(imol+2)%r(i)+diff(i)
         ENDDO
      ENDDO
        
      CALL cp_subsys_set(newsys(1)%subsys,particles=particles)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this could cause problems on various platforms!!!!!!!!!!!!!!!!
      CALL mc_energy_restart(test_env,new_energy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! accept or reject the move
!      vol=cell%hmat(1,1)**3
      w=dexp(-mc_par%BETA*((new_energy-old_energy)-&
        REAL(mc_par%nchain,dbl)/mc_par%BETA*3.0d0*DLOG(new_cell_length/&
             old_cell_length) &
        + mc_par%pressure*(new_cell_length**3-old_cell_length**3)))

      IF (w > 1.0d0 ) w=1.0

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! write cell length, volume, density, and trial displacement to a file
         OPEN(UNIT=mc_par%cl)
         DO i=1,moves%volume%attempts
            READ(mc_par%cl,*)
         ENDDO
         WRITE(mc_par%cl,*) nnstep,new_cell_length*0.529177d0,&
             (new_energy-old_energy)
         CLOSE(mc_par%cl)
         OPEN(UNIT=80)
         IF (moves%volume%attempts .GT. 1) THEN
            DO i=1,moves%volume%attempts-1
               READ(80,*)
            ENDDO
         ENDIF
         WRITE(80,*) nnstep,&
        -REAL(mc_par%nchain,dbl)/mc_par%BETA*3.0d0*DLOG(new_cell_length/&
             old_cell_length), &
             mc_par%pressure*(new_cell_length**3-old_cell_length**3)
         CLOSE(80)
         OPEN(UNIT=81)
         IF (moves%volume%attempts .GT. 1) THEN
            DO i=1,moves%volume%attempts-1
               READ(81,*)
            ENDDO
         ENDIF
         WRITE(81,*) nnstep,vol_dis*(0.529177d0)**3,rand
         CLOSE(81)

! accept the move
         moves%volume%successes=moves%volume%successes+1
         move_updates%volume%successes=move_updates%volume%successes+1

! update energies
         energy_check=energy_check+(new_energy-old_energy)
         old_energy=new_energy
         averages%ave_energy=averages%ave_energy*REAL(istep-1,dbl)&
                    /REAL(istep,dbl)+old_energy/REAL(istep,dbl)
! and heat capacity
         averages%ave_energy_squared=averages%ave_energy_squared*&
               REAL(istep-1,dbl)/&
               REAL(istep,dbl)+old_energy**2/&
               REAL(istep,dbl)

! update coordinates and cell lengths
         CALL force_env_retain(force_env)
         CALL force_env_release(force_env,error=error)
         CALL cp_subsys_release(oldsys(1)%subsys)
         CALL quickstep_create_force_env(force_env, globenv,&
                   new_h_matrix=new_h,error=error)
         CALL force_env_get(force_env,subsys=oldsys,&
                            error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
!              WRITE(6,*) 'KDFLSSDJKf'
!         CALL force_env_set_cell(force_env,cell_new,error=error)
!               WRITE(6,*) 'KDFLSSDJKf'
         DO imol=1,mc_par%nunits_tot
            particles_old%els(imol)%r(1:3)= &
                  particles%els(imol)%r(1:3)
            r_old(1:3,imol)= particles%els(imol)%r(1:3)
         ENDDO
         CALL cp_subsys_set(oldsys(1)%subsys,particles=particles_old)

! update cell length
         averages%ave_cell_length=averages%ave_cell_length*&
                 REAL(moves%volume%attempts-1,dbl)/&
                 REAL(moves%volume%attempts,dbl)+new_cell_length/&
                 REAL(moves%volume%attempts,dbl)

! update the classical energy, if need be
         IF ( mc_par % lclassical) THEN
            CALL mc_classical_energy(force_env,&
                                         moves%classical_energy,&
                                         mc_par)
            classical_energy_old=moves%classical_energy
         ENDIF

      ELSE

! reject the move
! write cell length, volume, density, and trial displacement to a file
         OPEN(UNIT=mc_par%cl)
         DO i=1,moves%volume%attempts
            READ(mc_par%cl,*)
         ENDDO
         WRITE(mc_par%cl,*) nnstep,old_cell_length*0.529177d0,&
             (new_energy-old_energy)
         CLOSE(mc_par%cl)
         OPEN(UNIT=80)
         IF (moves%volume%attempts .GT. 1) THEN
            DO i=1,moves%volume%attempts-1
               READ(80,*)
            ENDDO
         ENDIF
         WRITE(80,*) nnstep,&
        -REAL(mc_par%nchain,dbl)/mc_par%BETA*3.0d0*DLOG(new_cell_length/&
             old_cell_length), &
             mc_par%pressure*(new_cell_length**3-old_cell_length**3)
         CLOSE(80)
         OPEN(UNIT=81)
         IF (moves%volume%attempts .GT. 1) THEN
            DO i=1,moves%volume%attempts-1
               READ(81,*)
            ENDDO
         ENDIF
         WRITE(81,*) nnstep,vol_dis*(0.529177d0)**3,rand
         CLOSE(81)
! update energy
         averages%ave_energy=averages%ave_energy*REAL(istep-1,dbl)/&
                      REAL(istep,dbl)+old_energy/REAL(istep,dbl)
! and heat capacity
         averages%ave_energy_squared=averages%ave_energy_squared*&
               REAL(istep-1,dbl)/&
               REAL(istep,dbl)+old_energy**2/&
               REAL(istep,dbl)

! and cell length
         averages%ave_cell_length=averages%ave_cell_length*&
                 REAL(moves%volume%attempts-1,dbl)/&
                 REAL(moves%volume%attempts,dbl)+old_cell_length/&
                 REAL(moves%volume%attempts,dbl)

      ENDIF

! release some memory
      CALL force_env_release(test_env,error=error)
      CALL cp_subsys_release(newsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)


  END SUBROUTINE mc_volume_move

!*********************************************************************
  SUBROUTINE change_bond_length ( r_old,r_new,mc_par,mol,dis_length)

!   Purpose: A move changing a bond length in a molecule

!   History: - Creation (2.4.2004)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot),INTENT(IN) :: r_old
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot),INTENT(OUT) :: r_new
    INTEGER,INTENT(IN)    :: mol
    REAL(KIND=dp),INTENT(OUT)    :: dis_length
    
    CHARACTER(LEN=*), PARAMETER :: routineN = "change_bond_length"

! ** Local variables **
      REAL(KIND=dp) :: rand,old_length
      REAL(KIND=dp),DIMENSION(1:3) :: OH
      INTEGER :: hydrogen,ipart,i,handle
      
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! copy the incoming coordinates so we can change them
      DO ipart=1,mc_par%nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! pick which bond in the molecule at random
      CALL random_number(rand)
      IF ( rand .lt. 0.50d0) THEN
         hydrogen=1
      ELSE
         hydrogen=2
      ENDIF

! choose a displacement
      CALL random_number(rand)
      dis_length=mc_par%rmbond*2.0d0*(rand-0.5d0)

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

!*********************************************************************
  SUBROUTINE change_bond_angle ( r_old,r_new,mc_par,mol)

!   Purpose: A move changing a bond angle in water

!   History: - Creation (2.4.2004)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot),INTENT(IN) :: r_old
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par%nunits_tot),INTENT(OUT) :: r_new
    INTEGER,INTENT(IN)    :: mol
    
    CHARACTER(LEN=*), PARAMETER :: routineN = "change_bond_angle"

! ** Local variables **
      REAL(KIND=dp) :: rand,bis_length,old_length1,old_length2,&
                       old_angle,new_angle,dis_angle
      INTEGER :: ipart,i,handle
      REAL(KIND=dp),DIMENSION(1:3) :: OH1,OH2,bisector,bisector2,&
                 OH_between,OH_between1,OH_between2,bisector1
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! copy the incoming coordinates so we can change them
      DO ipart=1,mc_par%nunits_tot
         r_new(1:3,ipart)=r_old(1:3,ipart)
      ENDDO

! choose a displacement
      CALL random_number(rand)
      dis_angle=mc_par%rmangle*2.0d0*(rand-0.5d0)

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
