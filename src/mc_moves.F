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
!!     the various moves in Monte Carlo (MC) simulations
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
  USE cp_blacs_env,                    ONLY: start_blacs
  USE cp_error_handling,               ONLY: cp_error_type,&
                                             cp_error_init
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsystem_type,&
                                             cp_subsys_get
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE force_env_types,                 ONLY: force_env_calculate_energy,&
                                             force_env_release,&
                                             force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mc_environment_types,            ONLY: mc_environment_type,&
                                             get_mc_env
  USE mc_control,                      ONLY: mc_simulation_parameters_type
  USE mc_move_control,                 ONLY: move_update,&
                                             mc_moves_type
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env


  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_moves"

  PUBLIC :: transion,trans,rot
 
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE transion ( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error)

!   Purpose: A move displacing a single atom

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mc_environment_type), POINTER :: mc_env
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "transion"

!   *** Local variables ***
    real(KIND = dp) ::                 rand,dis_ion,w,new_energy
    TYPE(force_env_type), POINTER :: test_env
    TYPE(particle_list_type), POINTER :: particles_old
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: subsys,oldsys
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par
!   *** Local Counters ***
    INTEGER  ::                  i,iparticle

!   ---------------------------------------------------------------------------

! nullify some pointers
      NULLIFY(mc_par,test_env,particles_old,particles,subsys,oldsys)

! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans) == 0 .AND. &
         MOD(nnstep,mc_par%nunits_tot) == 0 ) THEN
         CALL move_update(mc_par,moves,"transion")
         WRITE(77,*) nnstep,' rmtraion = ',mc_par%rmion
      ENDIF


      CALL quickstep_create_force_env(test_env, blacs_env, globenv, error=error)
      CALL random_number(rand)
 
! record the attempt
      moves%transion%attempts=moves%transion%attempts+1

! move one ion in the system

! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

! call a random number to figure out which particle we're moving
      CALL random_number(rand)
      iparticle=aint(mc_par%nunits*rand)+1

! call a random number to figure out which direction we're moving
      CALL random_number(rand)
      i=aint(3*rand)+1

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dis_ion=mc_par%rmion*(rand-0.5d0)

! now let's grab the coordinates and do the move
      CALL force_env_get(test_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

      particles%els(iparticle)%r(i)=particles%els(iparticle)%r(i)+dis_ion

! calculate the new energy of the system
      CALL force_env_calculate_energy(test_env)
      CALL force_env_get(test_env,&
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

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep) &
                      +new_energy/real(nnstep)

! update coordinates
         CALL force_env_get(force_env,subsys=oldsys,error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)

         particles_old%els(iparticle)%r(i)= &
                     particles%els(iparticle)%r(i)

      ELSE

! reject the move
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep)+ &
                      new_energy/real(nnstep)
      ENDIF

! deallocate
      CALL force_env_release(test_env,error=error)

  END SUBROUTINE transion

!**********************************************************************
 SUBROUTINE trans ( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error)

!   Purpose: A move displacing a single molecule

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mc_environment_type), POINTER :: mc_env
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys,&
                                       oldsys,newsys
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(particle_list_type), POINTER :: particles_old
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "trans"

!   *** Local variables ***
    real(KIND = dp) ::                 rand,dis_mol,w,new_energy
    TYPE(force_env_type), POINTER :: test_env

!   *** Local Counters ***
    INTEGER  ::                  i,imolecule,iparticle

!   ---------------------------------------------------------------------------
 
! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

! initalize the error stuff
!      CALL cp_error_init(error) 

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans) == 0 .AND. &
         MOD(nnstep,mc_par%nunits_tot) == 0 ) THEN
         CALL move_update(mc_par,moves,"trans")
         WRITE(77,*) nnstep,'rmtrans = ',mc_par%rmtrans
      ENDIF


      CALL quickstep_create_force_env(test_env, blacs_env, globenv, error=error)
      CALL random_number(rand)

 ! record the attempt
      moves%trans%attempts=moves%trans%attempts+1

! move one molecule in the system

! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

! call a random number to figure out which molecule we're moving
      CALL random_number(rand)
      imolecule=aint(mc_par%nchain*rand)+1
      imolecule=3*imolecule-2

! call a random number to figure out which direction we're moving
      CALL random_number(rand)
      i=aint(3*rand)+1

! call a random number to figure out how far we're moving
      CALL random_number(rand)
      dis_mol=mc_par%rmtrans*(rand-0.5d0)

! now let's grab the coordinates and do the move
      CALL force_env_get(test_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

      DO iparticle=imolecule,imolecule+2
         particles%els(iparticle)%r(i)=&
             particles%els(iparticle)%r(i)+dis_mol
      ENDDO

! calculate the new energy of the system
      CALL force_env_calculate_energy(test_env)
      CALL force_env_get(test_env,&
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

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep) &
                      +new_energy/real(nnstep)

! update coordinates
         CALL force_env_get(force_env,subsys=oldsys,error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)

         DO iparticle=imolecule,imolecule+2
            particles_old%els(iparticle)%r(1:3)= &
                     particles%els(iparticle)%r(1:3)
         ENDDO


      ELSE

! reject the move
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep)+ &
                      new_energy/real(nnstep)
      ENDIF

! deallocate
      CALL force_env_release(test_env, error=error)

  END SUBROUTINE trans
!**********************************************************************
!**********************************************************************
 SUBROUTINE rot ( mc_env,force_env, moves,&
                        nnstep,ave_energy,old_energy,globenv,&
                        blacs_env,error)

!   Purpose: A move rotating a single molecule around it's center
!            of mass...only works for water

!   History: - Creation (10.14.2003)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mc_environment_type), POINTER :: mc_env
    TYPE(mc_simulation_parameters_type), POINTER :: mc_par
    TYPE(particle_list_type), POINTER :: particles
    TYPE(cp_subsystem_p_type), DIMENSION(:), POINTER :: subsys,&
                                               oldsys,newsys
    TYPE(cp_error_type),OPTIONAL,INTENT(INOUT) :: error
    TYPE(particle_list_type), POINTER :: particles_old
    REAL(KIND = dp),INTENT(INOUT) :: ave_energy,old_energy
    TYPE(mc_moves_type),POINTER  :: moves
    TYPE(force_env_type), POINTER :: force_env
    INTEGER,INTENT(IN)      :: nnstep

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "rot"
    REAL(KIND=dp),DIMENSION(1:3),PARAMETER :: &
                             mass=(/ 16.0d0,1.0d0,1.0d0 /)

!   *** Local variables ***
    real(KIND = dp) ::                 rand,dis_mol,w,new_energy
    TYPE(force_env_type), POINTER :: test_env

!   *** Local Counters ***
    INTEGER  ::                  imolecule,iparticle,dir,ii,iunit
    LOGICAL  ::                  lx,ly
    REAL(KIND=dp) ::             nxcm,nycm,nzcm,rx,ry,rz,rxnew,&
                                 rynew,rznew,dgamma,cosdg,sindg
!   ---------------------------------------------------------------------------
    
! initialize some stuff
      lx =.false.
      ly =.false.

! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

! update the ion translation max displacemnt, if necessary
      IF (MOD(nnstep,mc_par%iuptrans) == 0 .AND. &
         MOD(nnstep,mc_par%nunits_tot) == 0 ) THEN
         CALL move_update(mc_par,moves,"rot")
         WRITE(77,*) nnstep,'rmrot = ',mc_par%rmrot
      ENDIF


      CALL quickstep_create_force_env(test_env, blacs_env, globenv, error=error)
      CALL random_number(rand)

 ! record the attempt
      moves%rot%attempts=moves%rot%attempts+1

! move one molecule in the system

! grab some data from the mc_env
      CALL get_mc_env ( mc_env, mc_par = mc_par )

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
      CALL force_env_get(test_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

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
      dgamma=mc_par%rmrot*(rand-0.5d0)

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
      CALL force_env_calculate_energy(test_env)
      CALL force_env_get(test_env,&
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

! update energies
         old_energy=new_energy
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep) &
                      +new_energy/real(nnstep)

! update coordinates
         CALL force_env_get(force_env,subsys=oldsys,error=error)
         CALL cp_subsys_get(oldsys(1)%subsys, &
             particles=particles_old, error=error)
         DO iparticle=imolecule,imolecule+2
            particles_old%els(iparticle)%r(1:3)= &
                     particles%els(iparticle)%r(1:3)
         ENDDO

      ELSE

! reject the move
         ave_energy=ave_energy*real(nnstep-1)/real(nnstep)+ &
                      new_energy/real(nnstep)
      ENDIF

! deallocate
      CALL force_env_release(test_env, error=error)

  END SUBROUTINE rot
!**********************************************************************

END MODULE mc_moves

!**********************************************************************
