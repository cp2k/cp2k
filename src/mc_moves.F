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
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsystem_type,&
                                             cp_subsys_get,&
                                             cp_subsys_release
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_calculate_energy,&
                                             force_env_release,&
                                             force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invert_matrix
  USE mc_environment_types,            ONLY: mc_environment_type,&
                                             get_mc_env
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: mc_move_update,&
                                             mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_moves"

  PUBLIC :: mc_atom_translation,mc_molecule_translation,&
            mc_molecule_rotation,mc_volume_move
 
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE mc_atom_translation ( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error)

!   Purpose: A move displacing a single atom

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_atom_translation"

!   *** Local variables ***
    REAL(KIND = dp) ::                 rand,dis_ion,w,new_energy,part
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: subsys
!   *** Local Counters ***
    INTEGER  ::                  i,iparticle,handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)

! grab some data from the mc_env
!      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*mc_par%nunits_tot) == 0 ) THEN
         CALL mc_move_update(mc_par,move_updates,"transion")
         WRITE(77,*) nnstep,' rmtraion = ',mc_par%rmion,' a.u.'
      ENDIF


      CALL random_number(rand)
 
! record the attempt
      moves%transion%attempts=moves%transion%attempts+1
      move_updates%transion%attempts=move_updates%transion%attempts+1

! move one ion in the system

! call a random number to figure out which particle we're moving
      CALL random_number(rand)
      iparticle=aint(mc_par%nunits*rand)+1

! call a random number to figure out which direction we're moving
      CALL random_number(rand)
      i=aint(3*rand)+1

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dis_ion=mc_par%rmion*(rand-0.5d0)*2.0d0

! now let's grab the coordinates and do the move
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)
      part=particles%els(iparticle)%r(i)
      particles%els(iparticle)%r(i)=particles%els(iparticle)%r(i)+dis_ion

! calculate the new energy of the system
      CALL force_env_calculate_energy(force_env)
      CALL force_env_get(force_env,&
          potential_energy=new_energy,error=error)
 
! accept or reject the move based on the Metropolis rule
      IF (( new_energy-old_energy ) < 0.0d0 ) THEN
         w=1.0
      ELSE
         w=dexp(-mc_par%BETA*(new_energy-old_energy))
      ENDIF

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! accept the move
         moves%transion%successes=moves%transion%successes+1
         move_updates%transion%successes=&
               move_updates%transion%successes+1

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep) &
                      +old_energy/REAL(nnstep)

      ELSE

! reject the move
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep)+ &
                      old_energy/REAL(nnstep)

! replace coordinates
         particles%els(iparticle)%r(i)=part

      ENDIF

! release some memory
!   CALL cp_subsys_release(subsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_atom_translation

!**********************************************************************
 SUBROUTINE mc_molecule_translation ( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error)

!   Purpose: A move displacing a single molecule

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_translation"

!   *** Local variables ***
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    REAL(KIND = dp) ::                 rand,dis_mol,w,new_energy
    REAL(KIND = dp),DIMENSION(1:3,1:3) :: r          

!   *** Local Counters ***
    INTEGER  ::                  i,imolecule,iparticle,handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(particles,subsys)
 
! grab some data from the mc_env
!      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*mc_par%nunits_tot) == 0 ) THEN
         CALL mc_move_update(mc_par,move_updates,"trans")
         WRITE(77,*) nnstep,'rmtrans = ',mc_par%rmtrans,' a.u.'
      ENDIF


      CALL random_number(rand)

 ! record the attempt
      moves%trans%attempts=moves%trans%attempts+1
      move_updates%trans%attempts=move_updates%trans%attempts+1

! move one molecule in the system

! call a random number to figure out which molecule we're moving
      CALL random_number(rand)
      imolecule=aint(mc_par%nchain*rand)+1
      imolecule=3*imolecule-2

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

! save old coordinates
      r(1:3,1)=particles%els(imolecule)%r(1:3)
      r(1:3,2)=particles%els(imolecule+1)%r(1:3)
      r(1:3,3)=particles%els(imolecule+2)%r(1:3)

! do the move
      DO iparticle=imolecule,imolecule+2
         particles%els(iparticle)%r(i)=&
             particles%els(iparticle)%r(i)+dis_mol
      ENDDO

! calculate the new energy of the system
      CALL force_env_calculate_energy(force_env)
      CALL force_env_get(force_env,&
          potential_energy=new_energy,error=error)
 
 
! accept or reject the move based on the Metropolis rule
      IF (( new_energy-old_energy ) < 0.0d0 ) THEN
         w=1.0
      ELSE
         w=dexp(-mc_par%BETA*(new_energy-old_energy))
      ENDIF

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! accept the move
         moves%trans%successes=moves%trans%successes+1
         move_updates%trans%successes=move_updates%trans%successes+1

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep) &
                      +old_energy/REAL(nnstep)


      ELSE

! reject the move
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep)+ &
                      old_energy/REAL(nnstep)

! restore the old coordinates
         particles%els(imolecule)%r(1:3)=r(1:3,1)
         particles%els(imolecule+1)%r(1:3)=r(1:3,2)
         particles%els(imolecule+2)%r(1:3)=r(1:3,3)

      ENDIF

! release some memory
!   CALL cp_subsys_release(subsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_molecule_translation

!**********************************************************************
!**********************************************************************

 SUBROUTINE mc_molecule_rotation ( mc_par,force_env, moves,&
                        move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        error)

!   Purpose: A move rotating a single molecule around it's center
!            of mass...only works for water

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_molecule_rotation"
    REAL(KIND=dp),DIMENSION(1:3),PARAMETER :: &
                             mass=(/ 16.0d0,1.0d0,1.0d0 /)

!   *** Local variables ***
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(particle_list_type), POINTER :: particles
    REAL(KIND = dp) ::                 rand,dis_mol,w,new_energy
    REAL(KIND = dp),DIMENSION(1:3,1:3) :: r          

!   *** Local Counters ***
    INTEGER  ::                  imolecule,iparticle,dir,ii,iunit
    INTEGER  ::                  handle
    LOGICAL  ::                  lx,ly
    REAL(KIND=dp) ::             nxcm,nycm,nzcm,rx,ry,rz,rxnew,&
                                 rynew,rznew,dgamma,cosdg,sindg
!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! nullify some pointers
      NULLIFY(subsys,particles)

! initialize some stuff
      lx =.false.
      ly =.false.

! grab some data from the mc_env
!      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans*mc_par%nunits_tot) == 0 ) THEN
         CALL mc_move_update(mc_par,move_updates,"rot")
         WRITE(77,*) nnstep,'rmrot = ',mc_par%rmrot,' radians'
      ENDIF

      CALL random_number(rand)

 ! record the attempt
      moves%rot%attempts=moves%rot%attempts+1
      move_updates%rot%attempts=move_updates%rot%attempts+1

! rotate one molecule in the system

! call a random number to figure out which molecule we're moving
      CALL random_number(rand)
      imolecule=aint(mc_par%nchain*rand)+1
      imolecule=3*imolecule-2

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

! save the old coordinates of the molecule
      r(1:3,1)=particles%els(imolecule)%r(1:3)
      r(1:3,2)=particles%els(imolecule+1)%r(1:3)
      r(1:3,3)=particles%els(imolecule+2)%r(1:3)


! Determine new center of mass for chain i by finding the sum
! of m*r for each unit, then dividing by the total mass of the chain
      nxcm = 0.0d0
      nycm = 0.0d0
      nzcm = 0.0d0
      do ii = 1, mc_par%nunits
         nxcm = nxcm + particles%els(imolecule-1+ii)%r(1)* mass(ii)
         nycm = nycm + particles%els(imolecule-1+ii)%r(2)* mass(ii)
         nzcm = nzcm + particles%els(imolecule-1+ii)%r(3)* mass(ii)
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
 
         DO  iunit = imolecule,imolecule+2
            ry = particles%els(iunit)%r(2) - nycm
            rz = particles%els(iunit)%r(3) - nzcm
            rynew = cosdg * ry + sindg * rz
            rznew = cosdg * rz - sindg * ry
               
            particles%els(iunit)%r(2) = rynew + nycm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO
      ELSEIF (ly) THEN 
 
! ***    ROTATE UNITS OF I AROUND y-AXIS ***
 
         DO  iunit = imolecule,imolecule+2
            rx = particles%els(iunit)%r(1) - nxcm
            rz = particles%els(iunit)%r(3) - nzcm
            rxnew = cosdg * rx + sindg * rz
            rznew = cosdg * rz - sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(3) = rznew + nzcm

         ENDDO

      ELSE
 
! ***    ROTATE UNITS OF I AROUND z-AXIS ***
 
         DO  iunit = imolecule,imolecule+2
            rx = particles%els(iunit)%r(1) - nxcm
            ry = particles%els(iunit)%r(2) - nycm
            rxnew = cosdg * rx + sindg * ry
            rynew = cosdg * ry - sindg * rx
               
            particles%els(iunit)%r(1) = rxnew + nxcm
            particles%els(iunit)%r(2) = rynew + nycm

         ENDDO

      ENDIF 

! calculate the new energy of the system
      CALL force_env_calculate_energy(force_env)
      CALL force_env_get(force_env,&
          potential_energy=new_energy,error=error)

! accept or reject the move based on the Metropolis rule
      IF (( new_energy-old_energy ) < 0.0d0 ) THEN
         w=1.0
      ELSE
         w=dexp(-mc_par%BETA*(new_energy-old_energy))
      ENDIF

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! accept the move
         moves%rot%successes=moves%rot%successes+1
         move_updates%rot%successes=move_updates%rot%successes+1

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep) &
                      +old_energy/REAL(nnstep)

      ELSE

! reject the move
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep)+ &
                      old_energy/REAL(nnstep)

! return to the old coordinates
         particles%els(imolecule)%r(1:3)=r(1:3,1)
         particles%els(imolecule+1)%r(1:3)=r(1:3,2)
         particles%els(imolecule+2)%r(1:3)=r(1:3,3)

      ENDIF

! release some memory
!   CALL cp_subsys_release(subsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)

  END SUBROUTINE mc_molecule_rotation
!**********************************************************************
  SUBROUTINE mc_volume_move ( mc_par,force_env, moves,move_updates,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,ave_density,error)

!   Purpose: A move changing the size of the box

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy,ave_density
    TYPE(mc_moves_type),POINTER  :: moves,move_updates
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(mc_simulation_parameters_type),POINTER :: mc_par

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "mc_volume_move"

!   *** Local variables ***
    REAL(KIND = dp) ::                 rand,dis_length,w,new_energy,&
                                       old_length,max_error,vol,&
                                       density_old
    TYPE(force_env_type), POINTER :: test_env
    TYPE(particle_list_type), POINTER :: particles_old
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: oldsys,&
                                                      newsys
    TYPE(cell_type), POINTER :: cell,cell_new
    REAL(KIND =dp),DIMENSION(1:3,1:3) ::       h_matrix,new_h,&
                                               new_h_inverse
!   *** Local Counters ***
    INTEGER  ::                  imol,handle

!   ---------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)


! nullify some pointers
      NULLIFY(test_env,particles_old,particles,oldsys)

! grab some data from the mc_env
!      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the volume move max displacement, if necessary
      IF (MOD(nnstep,mc_par%iupvolume*mc_par%nunits_tot) == 0 ) THEN
         CALL mc_move_update(mc_par,move_updates,"volume")
         WRITE(77,*) nnstep,' rmvolume = ',mc_par%rmvolume,' a.u.^3'
      ENDIF


 ! record the attempt
      moves%volume%attempts=moves%volume%attempts+1
      move_updates%volume%attempts=move_updates%volume%attempts+1

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dis_length=(mc_par%rmvolume)**(1.0d0/3.0d0)*(rand-0.5d0)*2.0d0

! now let's grab the cell length and change it 
      CALL force_env_get(force_env,cell=cell,&
                         error=error)

! find the old density
      density_old=mc_par%nchain/cell%hmat(1,1)**3

! create the old h matrix
      h_matrix(:,:)=cell%hmat(:,:)
      
! create the new h matrix
      new_h(:,:)=0.0d0
      new_h(1,1)=cell%hmat(1,1)+dis_length
      new_h(2,2)=cell%hmat(2,2)+dis_length
      new_h(3,3)=cell%hmat(3,3)+dis_length


! invert it
      CALL invert_matrix(new_h,new_h_inverse,max_error)

! create the new environment to change particle coordinates
      CALL quickstep_create_force_env(test_env, blacs_env, globenv, &
                   new_h_matrix=new_h,error=error)
      CALL force_env_get(test_env,subsys=newsys,&
          error=error)
      CALL cp_subsys_get(newsys(1)%subsys, &
             particles=particles, error=error)

! now we need to scale the coordinates of all the molecules...notice
! this procedure scales the bond lengths as well, which means our
! acceptance rule is slightly different

      DO imol=1,mc_par%nchain*mc_par%nunits
         particles%els(imol)%r(:)=MATMUL(MATMUL(h_matrix,&
                particles%els(imol)%r(:)),new_h_inverse)
      ENDDO

      CALL force_env_calculate_energy(test_env)
      CALL force_env_get(test_env,cell=cell_new,&
          potential_energy=new_energy,error=error)
 
! accept or reject the move
      vol=cell_new%hmat(1,1)**3
      w=dexp(-mc_par%BETA*((new_energy-old_energy)-REAL(mc_par%nchain*&
             mc_par%nunits_tot)/mc_par%BETA*LOG((vol+dis_length**3)/vol)&
             +mc_par%pressure*vol))
      
      IF (w > 1.0d0 ) w=1.0

      CALL random_number(rand)
      IF (rand .lt. w) THEN

! accept the move
         moves%volume%successes=moves%volume%successes+1
         move_updates%volume%successes=move_updates%volume%successes+1

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep) &
                      +old_energy/REAL(nnstep)

! update coordinates and cell lengths
         CALL force_env_release(force_env,error=error)
         CALL quickstep_create_force_env(force_env, blacs_env, globenv, &
                   new_h_matrix=new_h,error=error)

         CALL force_env_get(force_env,subsys=oldsys,&
                            error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)

         DO imol=1,mc_par%nchain*mc_par%nunits
            particles_old%els(imol)%r(1:3)= &
                  particles%els(imol)%r(1:3)
         ENDDO


! update density
         ave_density=ave_density*REAL(moves%volume%attempts-1)/&
                 REAL(moves%volume%attempts)+mc_par%nchain/vol/&
                 REAL(moves%volume%attempts)

! write density to a file
         WRITE(80,*) mc_par%nchain/vol


      ELSE

! reject the move
         ave_energy=ave_energy*REAL(nnstep-1)/REAL(nnstep)+ &
                      old_energy/REAL(nnstep)

! update density
         ave_density=ave_density*REAL(moves%volume%attempts-1)/&
                 REAL(moves%volume%attempts)+density_old/&
                 REAL(moves%volume%attempts)

! write density to a file
         WRITE(80,*) density_old

      ENDIF

! release some memory
      CALL force_env_release(test_env,error=error)
!      CALL cp_subsys_release(oldsys(1)%subsys)
!      CALL cp_subsys_release(newsys(1)%subsys)

! end the timing
  CALL timestop(0.0d0,handle)


  END SUBROUTINE mc_volume_move

!*********************************************************************

END MODULE mc_moves

!**********************************************************************
