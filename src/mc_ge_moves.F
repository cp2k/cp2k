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
                                             cp_subsystem_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_misc,                         ONLY: mc_averages_type,&
                                             mc_classical_energy
  USE mc_move_control,                 ONLY: mc_move_update,&
                                             mc_moves_type,&
                                             move_q_reinit,&
                                             q_move_accept
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_ge_moves"

  PUBLIC :: mc_ge_quickstep_move,mc_ge_volume_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE mc_ge_quickstep_move(mc_par_1,mc_par_2,force_env_1,&
                force_env_2,moves_1,moves_2,&
                move_updates_1,move_updates_2,&
                energy_check,r_old_1,r_old_2,istep,old_energy,&
                averages_1,averages_2,box_flag,error)

!   Purpose: A move calling the Quickstep energy routines for the Gibbs
!            ensemble

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par_1,mc_par_2
    TYPE(force_env_type), POINTER            :: force_env_1,force_env_2
    TYPE(mc_moves_type), POINTER             :: moves_1,moves_2,&
                                          move_updates_1,move_updates_2
    REAL(KIND=dp), DIMENSION(1:2),INTENT(INOUT)  :: energy_check,&
                                                    old_energy
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par_1%nunits_tot), &
      INTENT(INOUT)                          :: r_old_1
    REAL(KIND=dp), &
      DIMENSION(1:3, 1:mc_par_2%nunits_tot), &
      INTENT(INOUT)                          :: r_old_2
    INTEGER, DIMENSION(1:2),INTENT(IN)       :: box_flag
    INTEGER, DIMENSION(1:2),INTENT(INOUT)       :: istep
    TYPE(mc_averages_type), POINTER          :: averages_1,averages_2
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER     :: routineN = "mc_ge_quickstep_move"

    INTEGER                                  :: handle, i, iparticle
    REAL(KIND=dp),DIMENSION(1:2)             :: classical_energy_new, &
                                                classical_energy_old, &
                                                new_energy
    REAL(KIND=dp)           :: w,rand
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_1,subsys_2
    TYPE(particle_list_type), POINTER        :: particles_1,particles_2

!   *** Local Counters ***
!   ---------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles_1,subsys_1)
      NULLIFY(particles_2,subsys_2)

! record the attempt
      moves_1%Quickstep%attempts=moves_1%Quickstep%attempts+1
      moves_2%Quickstep%attempts=moves_2%Quickstep%attempts+1

! calculate the new energy of the system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this could cause problems on various platforms!!!!!!!!!!!!!!!!
      IF(box_flag(1) ==1) THEN
!         CALL mc_energy_restart(force_env_1,new_energy(1))
      ELSE
         new_energy(1)=old_energy(1)
      ENDIF
      IF(box_flag(2) ==1) THEN
!         CALL mc_energy_restart(force_env_2,new_energy(2))
      ELSE
         new_energy(2)=old_energy(2)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! now let's grab the coordinates
      CALL force_env_get(force_env_1,subsys=subsys_1,error=error)
      CALL cp_subsys_get(subsys_1(1)%subsys, &
             particles=particles_1, error=error)
      CALL force_env_get(force_env_2,subsys=subsys_2,error=error)
      CALL cp_subsys_get(subsys_2(1)%subsys, &
             particles=particles_2, error=error)

! accept or reject the move based on Metropolis or the Iftimie rule
      IF (mc_par_1%lclassical) THEN
!         CALL mc_classical_energy(force_env_1,&
!                                  classical_energy_new(1),&
!                                  mc_par_1)
!         CALL mc_classical_energy(force_env_2,&
!                                  classical_energy_new(2),&
!                                  mc_par_2)
         classical_energy_old(1)=moves_1%classical_energy
         classical_energy_old(2)=moves_2%classical_energy

       w=dexp(-mc_par_1%BETA*((new_energy(1)-classical_energy_new(1))&
                   -(old_energy(1)-classical_energy_old(1))&
                   +(new_energy(2)-classical_energy_new(2))&
                   -(old_energy(2)-classical_energy_old(2))))

         write(6,*) 'w(1) =',w,new_energy(1)-old_energy(1),&
                 classical_energy_new(1),classical_energy_old(1)
         write(6,*) 'w(2) =',w,new_energy(2)-old_energy(2),&
                 classical_energy_new(2),classical_energy_old(2)
      ELSE
         w=dexp(-mc_par_1%BETA*(new_energy(1)-old_energy(1)+&
                              new_energy(2)-old_energy(2)))
      ENDIF

      IF ( w .GE. 1.0d0 ) THEN
         w=1.0d0
         rand=0.0D0
      ELSE
         CALL random_number(rand)
      ENDIF

      IF (rand .lt. w) THEN

! accept the move
         moves_1%Quickstep%successes=moves_1%Quickstep%successes+1
         moves_2%Quickstep%successes=moves_2%Quickstep%successes+1

! if we're not biasing, we need to record all accepted moves since last
! Quickstep calculation
         IF ( .NOT. mc_par_1%lclassical) THEN
!            CALL q_move_accept(moves_1,move_updates_1)
!            CALL q_move_accept(moves_2,move_updates_2)
            CALL q_move_accept(moves_1)
            CALL q_move_accept(moves_2)

! reset the counters
            CALL move_q_reinit(moves_1)
!            CALL move_q_reinit(move_updates_1)
            CALL move_q_reinit(moves_2)
!            CALL move_q_reinit(move_updates_2)

         ENDIF

! update energies for the boxes that were displaced
         IF (box_flag(1) == 1) THEN
            energy_check(1)=energy_check(1)+(new_energy(1)-&
                                              old_energy(1))
            old_energy(1)=new_energy(1)
            averages_1%ave_energy=averages_1%ave_energy*&
                     REAL(istep(1)-1,dbl)&
                     /REAL(istep(1),dbl)+old_energy(1)/&
                     REAL(istep(1),dbl)
! and heat capacity
            averages_1%ave_energy_squared=&
               averages_1%ave_energy_squared*&
               REAL(istep(1)-1,dbl)/&
               REAL(istep(1),dbl)+old_energy(1)**2/&
               REAL(istep(1),dbl)
! update coordinates
            DO iparticle=1,mc_par_1%nunits_tot
               r_old_1(1:3,iparticle)=particles_1%els(iparticle)%r(1:3)
            ENDDO
! if we're biasing classically, then we update the classical energy
            IF ( mc_par_1 % lclassical) THEN
!               CALL mc_classical_energy(force_env_1,&
!                                         moves_1%classical_energy,&
!                                         mc_par_1)
               moves_1%classical_energy=classical_energy_old(1)
            ENDIF
         ENDIF
         IF (box_flag(2) == 1) THEN
            energy_check(2)=energy_check(2)+(new_energy(2)-&
                                              old_energy(2))
            old_energy(2)=new_energy(2)
            averages_2%ave_energy=averages_2%ave_energy*&
                     REAL(istep(2)-1,dbl)&
                     /REAL(istep(2),dbl)+old_energy(2)/&
                     REAL(istep(2),dbl)
! and heat capacity
            averages_2%ave_energy_squared=&
               averages_2%ave_energy_squared*&
               REAL(istep(2)-1,dbl)/&
               REAL(istep(2),dbl)+old_energy(2)**2/&
               REAL(istep(2),dbl)
! update coordinates
            DO iparticle=1,mc_par_2%nunits_tot
               r_old_2(1:3,iparticle)=particles_2%els(iparticle)%r(1:3)
            ENDDO
! if we're biasing classically, then we update the classical energy
            IF ( mc_par_1 % lclassical) THEN
!               CALL mc_classical_energy(force_env_2,&
!                                         moves_2%classical_energy,&
!                                         mc_par_2)
               moves_2%classical_energy=classical_energy_old(2)
            ENDIF
         ENDIF

      ELSE

! reject the move
         IF ( .NOT. mc_par_1 % lclassical ) THEN
            CALL move_q_reinit(moves_1)
            CALL move_q_reinit(move_updates_1)
            CALL move_q_reinit(moves_1)
            CALL move_q_reinit(move_updates_1)
         ENDIF

         IF(box_flag(1) == 1) THEN
            averages_1%ave_energy=averages_1%ave_energy*&
               REAL(istep(1)-1,dbl)/&
               REAL(istep(1),dbl)+old_energy(1)/&
               REAL(istep(1),dbl)
! and heat capacity
            averages_1%ave_energy_squared=&
                 averages_1%ave_energy_squared*&
                 REAL(istep(1)-1,dbl)/&
                 REAL(istep(1),dbl)+old_energy(1)**2&
                 /REAL(istep(1),dbl)
! replace coordinates
            DO iparticle=1,mc_par_1%nunits_tot
               particles_1%els(iparticle)%r(1:3)=r_old_1(1:3,iparticle)
            ENDDO
            CALL cp_subsys_set(subsys_1(1)%subsys,particles=particles_1)
! if we need to restore the classical energy
            IF ( mc_par_1 % lclassical ) THEN
               classical_energy_old(1)=moves_1%classical_energy
            ENDIF
         ENDIF
         IF(box_flag(2) == 1) THEN
            averages_2%ave_energy=averages_2%ave_energy*&
               REAL(istep(2)-1,dbl)/&
               REAL(istep(2),dbl)+old_energy(2)/&
               REAL(istep(2),dbl)
! and heat capacity
               averages_2%ave_energy_squared=&
                 averages_2%ave_energy_squared*&
                 REAL(istep(2)-1,dbl)/&
                 REAL(istep(2),dbl)+old_energy(2)**2&
                 /REAL(istep(2),dbl)
! replace coordinates
            DO iparticle=1,mc_par_2%nunits_tot
               particles_2%els(iparticle)%r(1:3)=r_old_2(1:3,iparticle)
            ENDDO
            CALL cp_subsys_set(subsys_2(1)%subsys,particles=particles_2)
! if we need to restore the classical energy
            IF ( mc_par_1 % lclassical ) THEN
               classical_energy_old(2)=moves_2%classical_energy
            ENDIF
         ENDIF

      ENDIF

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_ge_quickstep_move

! *****************************************************************************

!**********************************************************************
  SUBROUTINE mc_ge_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,averages,old_energy,globenv,&
                        energy_check,r_old,error)

!   Purpose: A move changing the size of the box
! NOT YET SUITIBLE...COPY FROM mc_moves BEFORE CHANGING TO TWO BOXES
!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    REAL(KIND = dp),INTENT(INOUT) :: old_energy,energy_check
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(mc_averages_type),POINTER  :: averages
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par
    REAL(KIND=dp),DIMENSION(1:3,1:mc_par%nunits_tot), &
                               INTENT(INOUT) :: r_old

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_ge_volume_move"
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
    TYPE(cp_subsystem_type), POINTER :: subsys
    TYPE(cell_type), POINTER :: cell
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
         CALL mc_move_update(mc_par,move_updates,"volume",istep)
         WRITE(mc_par%rm,*) nnstep,' rmvolume = ',&
                 mc_par%rmvolume*(angstrom)**3,' angstroms^3'
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
      CALL force_env_get(force_env,subsys=oldsys,&
                            error=error)
      CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
      CALL force_env_get(test_env,subsys=newsys,&
          error=error)
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
!      CALL mc_energy_restart(test_env,new_energy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! accept or reject the move
      vol=cell%hmat(1,1)**3
      w=dexp(-mc_par%BETA*((new_energy-old_energy)-&
        REAL(mc_par%nchain,dbl)/mc_par%BETA*3.0d0*LOG(new_cell_length/&
             old_cell_length) &
        + mc_par%pressure*(new_cell_length**3-old_cell_length**3)))

      IF ( w .GE. 1.0d0 ) THEN
         w=1.0d0
         rand=0.0D0
      ELSE
         CALL random_number(rand)
      ENDIF

      IF (rand .lt. w) THEN

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
         CALL force_env_release(force_env,error=error)
         CALL cp_subsys_release(oldsys(1)%subsys)
         CALL quickstep_create_force_env(force_env, globenv,&
                   new_h_matrix=new_h,error=error)
         CALL force_env_get(force_env,subsys=oldsys,&
                            error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)

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

! write cell length, volume, density, and trial displacement to a file
         OPEN(UNIT=mc_par%cl)
         DO i=1,moves%volume%attempts
            READ(mc_par%cl,*)
         ENDDO
         WRITE(mc_par%cl,*) nnstep,new_cell_length*angstrom,&
           (new_cell_length*angstrom)**3&
           ,REAL(mc_par%nchain,dbl)/(new_cell_length*angstrom&
            *1.0d-8)**3/6.022d23*mass_tot,&
           vol_dis*(angstrom)**3
         CLOSE(mc_par%cl)

      ELSE

! reject the move
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

! write cell length to a file
         OPEN(UNIT=mc_par%cl)
         DO i=1,moves%volume%attempts
            READ(mc_par%cl,*)
         ENDDO
         WRITE(mc_par%cl,*) nnstep,old_cell_length*0.529177d0
         CLOSE(mc_par%cl)

      ENDIF

! release some memory
      CALL force_env_release(test_env,error=error)
      CALL cp_subsys_release(newsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)


  END SUBROUTINE mc_ge_volume_move

!*********************************************************************

END MODULE mc_ge_moves

!**********************************************************************
