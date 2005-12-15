!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc_run [1.0] *
!!
!!   NAME
!!     mc_run
!!
!!   FUNCTION
!!     preps the system for a Monte Carlo run (sets up some environments,
!!     calls the routines to read in the MC parameters)...converted 
!!     from qs_mc.F
!!
!!   AUTHOR
!!     MJM-Oct-15-03
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!!
!!   REFERENCE LIST: a list of papers for the theory behind various MC moves
!!   Books:
!!      D. Frenkel, B. Smit: Understanding Molecular Simulation (1996)
!!      M.P. Allen, D.J. Tildesley: Computer Simulations of Liquids (1987)
!!
!!   Aggregation volume bias Monte Carlo (AVBMC):
!!      Chen, B.; Siepmann, J.I.  J. Phys. Chem. B 2000, 104, 8725.
!!
!!   Biasing with an inexpensive potential:
!!      Iftimie et al.  J. Chem. Phys. 2000, 113, 4852.
!!      Gelb, L. D.  J. Chem. Phys. 2003, 118, 7747.
!!
!!   Configurational bias Monte Carlo (CBMC):
!!      Siepmann, J.I.; Frenkel, D.  Mol. Phys. 1992, 75, 59.
!!      
!!   Gibbs ensemble Monte Carlo (GEMC):
!!      Panagiotopoulos, A.Z.  Mol. Phys. 1987, 61, 813.
!!      Panagiotopoulos et al.  Mol. Phys. 1988, 63, 527.
!!      Smit et al.  Mol. Phys. 1989, 68, 931.
!!
!!   Isobaric-isothermal ensemble:
!!      McDonald, I.R.  Mol. Phys. 1972, 23, 41.
!!
!!   Original Monte Carlo paper:
!!      Metropolis et al.  J. Chem. Phys. 1953, 21, 1087.
!!   
!!   
!******************************************************************************

MODULE mc_run
  USE cp2k_units,                      ONLY: init_cp2k_units
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             force_env_release,&
                                             use_fist_force
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mc_control,                      ONLY: read_mc_restart,&
                                             mc_create_force_env
  USE mc_ensembles,                    ONLY: mc_run_ensemble
  USE mc_environment_types,            ONLY: mc_env_release,&
                                             mc_env_create,&
                                             mc_environment_p_type,&
                                             set_mc_env,&
                                             get_mc_env
  USE mc_types,                        ONLY: mc_simulation_parameters_type,&
                                             read_mc_section,&
                                             determine_conf_prob,&
                                             get_mc_par,&
                                             set_mc_par,&
                                             find_mc_rcut
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: rng_stream_type,&
                                             create_rng_stream,&
                                             delete_rng_stream
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_run"

  PUBLIC :: do_mon_car
  
!!***
!-----------------------------------------------------------------------------!

CONTAINS

! *****************************************************************************
!!****s* mc_run.F/do_mon_car
!!
!!   NAME
!!     do_mon_car
!!
!!   SYNOPSIS
!!     Subroutine do_mon_car ( force_env, globenv, error )
!!       Type(global_environment_type), Pointer   :: globenv
!!       Type(force_env_type), Pointer            :: force_env
!!       Type(cp_error_type), Intent(inout),Optional      :: error
!!     End Subroutine do_mon_car
!!
!!   FUNCTION
!!     starts the Monte Carlo simulation and determines which ensemble we're
!!     running
!!
!!   ARGUMENTS
!!     - globenv: the global environment for the simulation
!!     - force_env: the force environment for the simulation, or
!!                  the force environment for box 1, depending on which
!!                  ensemble we're running
!!     - error: the error type
!!
!!     Designed for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************

  SUBROUTINE do_mon_car ( force_env, globenv, error )


! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

! Local Parameters
   CHARACTER(LEN=*), PARAMETER :: routine_name = "do_mon_car", &
      routineP = module_name//':'//routine_name

! Local Variables
    TYPE(force_env_type), POINTER            :: force_env_2
    INTEGER                                  :: isos, iw,&
      nchain,nunits_tot,nunits,ibox,nboxes,in_use,&
      nstep,iupvolume,iuptrans,nmoves,nswapmoves,iprint
    LOGICAL                                  :: lrestart,ionode,lterminate,&
      lbias,lstop
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys,subsys_2
    TYPE(mc_environment_p_type), POINTER,DIMENSION(:) &
                                             :: mc_env
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par,mc_par_2
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new_2
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_2
    TYPE(rng_stream_type), POINTER           :: rng_stream
    CHARACTER(LEN=20)                        :: ensemble
    CHARACTER(LEN=40)                        :: dat_file,box2_file
    REAL(KIND=dp)                            :: rcut,pmswap,temp,pressure,&
      pmvolume,pmtraion,pmtrans,BETA
    REAL(KIND=dp),DIMENSION(1:3)             :: conf_prob
!------------------------------------------------------------------------------

      NULLIFY(mc_env,mc_par,mc_par_2,force_env_2,rng_stream)

      CALL init_cp2k_units ()
      iw = globenv % scr

! set some values...will use get_globenv if that ever comes around
      ionode=globenv % ionode

! allocating or nullifying local pointers
      ALLOCATE ( mc_par, STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par', 0 )
 
! initialize the random numbers
      CALL create_rng_stream(rng_stream,'first')

! read the input of the Monte Carlo section
      CALL read_mc_section ( mc_par, globenv ,error)

! if we're doing a classical simulation, figure out the largest
! potential cutoff and write it to the screen
      IF(force_env%in_use .EQ. use_fist_force) THEN
         CALL find_mc_rcut(mc_par,force_env,lterminate)
         CALL get_mc_par(mc_par,rcut=rcut)
         IF(ionode) WRITE (iw,'( A,T73,F8.4 )') &
            ' MC| Interaction cutoff [angstroms]',rcut*angstrom
         IF(lterminate) THEN
            CALL stop_program('mc_run.F/do_mon_car',&
               'Cutoff larger than twice the boxlength')
         ENDIF
      ENDIF

! determine the total number of molecules and atoms
      CALL force_env_get(force_env,subsys=subsys,error=error)
      CALL cp_subsys_get(subsys(1)%subsys, &
         molecule_kinds_new=molecule_kinds_new, error=error)
      molecule_kind => molecule_kinds_new%els(1)
      CALL get_molecule_kind(molecule_kind,nmolecule=nchain,natom=nunits)

! find the total number of interaction sites in the system
      nunits_tot=nunits*nchain

! stuff that info into mc_par
      CALL set_mc_par(mc_par,nchain=nchain,nunits=nunits,&
         nunits_tot=nunits_tot)
      
! write out the data
      IF (ionode) THEN
         WRITE (iw,'( A,T71,I10 )') ' MC| number of molecules', &
            nchain
         WRITE (iw,'( A,A,T71,I10 )') ' MC| number of interactions ',&
            'sites per molecule', nunits
         WRITE (iw,'( A,A,T71,I10 )') ' MC| total number of inter',&
            'action sites', nunits_tot
      ENDIF

! figure out the percentages of internal conformation moves we're trying
      CALL determine_conf_prob(mc_par,molecule_kind)

! get some more data
      CALL get_mc_par(mc_par,lrestart=lrestart,dat_file=dat_file,&
         ensemble=ensemble,conf_prob=conf_prob)

! write out move probabilities
      IF (ionode) THEN
         WRITE (iw,'( A,T75,F6.2 )') ' MC| Percent of bond moves tried', &
            conf_prob(1)*100.0d0
         WRITE (iw,'( A,T75,F6.2 )') ' MC| Percent of angle moves tried', &
            conf_prob(2)*100.0d0
         WRITE (iw,'( A,T75,F6.2 )') &
            ' MC| Percent of dihedral moves tried',conf_prob(3)*100.0d0
      ENDIF

! figure out how many boxes we have
      SELECT CASE(ensemble)
      CASE("TRADITIONAL")
         nboxes=1
      CASE("GEMC")
         nboxes=2
      END SELECT

! create the mc_env pointers
      ALLOCATE(mc_env(1:nboxes),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory('allocstruc','mc_env',nboxes)

      DO ibox=1,nboxes
! allocates an mc_env and sets the variables to zero
         CALL mc_env_create ( mc_env(ibox)%mc_env, &
            error = error )
      ENDDO


! check for restart
      IF ( lrestart ) THEN
         CALL read_mc_restart(mc_par,force_env,iw,rng_stream)
! release the old force env and make the new one
         in_use=force_env%in_use
         CALL force_env_release(force_env)
         globenv%input_file_name=dat_file
         CALL mc_create_force_env(force_env, globenv, .TRUE., &
              in_use,error=error)
      ENDIF

! now perform the calculation
      SELECT CASE(ensemble)
      CASE("TRADITIONAL")

! make sure things aren't wonky
         CALL get_mc_par(mc_par,pmswap=pmswap)
         IF(pmswap .GT. 0.0E0_dp) THEN
            CALL stop_program('mc_run.F/do_mon_car',&
               'pmswap is greater than 0 for a single box system!')
         ENDIF

! pull all the information together into one environment
         CALL set_mc_env( mc_env(1)%mc_env, &
            mc_par = mc_par, force_env=force_env)
         CALL mc_run_ensemble(mc_env,globenv,nboxes,error,rng_stream)

      CASE("GEMC")

         CALL get_mc_par(mc_par,box2_file=box2_file)

! prepare the second box

! allocating or nullifying local pointers
         ALLOCATE ( mc_par_2, STAT = isos )
         IF ( isos /= 0 ) CALL stop_memory ( 'allocstruc', 'mc_par_2', 0 )
 
! read the input of the Monte Carlo section
         globenv%input_file_name=box2_file
         CALL mc_create_force_env(force_env_2, globenv, .FALSE., &
            force_env%in_use,error=error)

         IF (ionode) THEN
            WRITE(iw,*)
            WRITE(iw,*)
            WRITE(iw,*) '******************************** Begin',&
            ' BOX 2 ',&
            '**********************************'
            WRITE(iw,*) '******************************************',&
            '*************************************'
         ENDIF
         
         CALL read_mc_section ( mc_par_2, globenv ,error)

! determine the total number of molecules and atoms
         CALL force_env_get(force_env_2,subsys=subsys_2,error=error)
         CALL cp_subsys_get(subsys_2(1)%subsys, &
         molecule_kinds_new=molecule_kinds_new_2, error=error)
         molecule_kind_2 => molecule_kinds_new_2%els(1)
         CALL get_molecule_kind(molecule_kind_2,nmolecule&
            =nchain,natom=nunits)

! find the total number of interaction sites in the system
         nunits_tot=nunits*nchain

         CALL set_mc_par(mc_par_2,nchain=nchain,nunits=nunits,&
         nunits_tot=nunits_tot)
! write out the data
         IF (ionode) THEN
            WRITE(iw,'( A,T71,I10 )' ) ' MC| number of molecules', &
            nchain
            WRITE(iw,'( A,A,T71,I10 )') ' MC| number of interactions ',&
            'sites per molecule', nunits
            WRITE(iw,'( A,A,T71,I10 )' ) ' MC| total number of inter',&
            'action sites', nunits_tot
         ENDIF

! figure out the percentages of internal conformation moves we're trying
      CALL determine_conf_prob(mc_par_2,molecule_kind_2)

! get some more data
      CALL get_mc_par(mc_par_2,conf_prob=conf_prob)

! write out move probabilities
         IF (ionode) THEN
            WRITE (iw,'( A,T75,F6.2 )') &
               ' MC| Percent of bond moves tried', &
               conf_prob(1)*100.0d0
            WRITE (iw,'( A,T75,F6.2 )') &
               ' MC| Percent of angle moves tried', &
               conf_prob(2)*100.0d0
            WRITE (iw,'( A,T75,F6.2 )') &
               ' MC| Percent of dihedral moves tried',conf_prob(3)*100.0d0
! make sure things are consistent between the two environments
            WRITE(iw,* ) 
            WRITE(iw,'( A,A )' ) 'Ignoring some input for box 2, and ',&
              'using the values for box 1 for the following variables:'
            WRITE(iw,'( A,A )' ) 'nstep,iupvolume,iuptrans,nmoves,',&
            'nswapmoves,iprint,lbias,lstop,temp,pressure,pmswap,'
            WRITE(iw,'( A,A )' ) 'pmvolume,pmtraion,pmtrans,BETA,rcut,',&
            'lrestart'
         ENDIF
         CALL get_mc_par(mc_par,nstep=nstep,iupvolume=iupvolume,&
            iuptrans=iuptrans,nmoves=nmoves,nswapmoves=nswapmoves,&
            iprint=iprint,lbias=lbias,lstop=lstop,temp=temp,&
            pressure=pressure,pmswap=pmswap,pmvolume=pmvolume,&
            pmtraion=pmtraion,pmtrans=pmtrans,BETA=BETA,rcut=rcut,&
            lrestart=lrestart)
         CALL set_mc_par(mc_par_2,nstep=nstep,iupvolume=iupvolume,&
            iuptrans=iuptrans,nmoves=nmoves,nswapmoves=nswapmoves,&
            iprint=iprint,lbias=lbias,lstop=lstop,temp=temp,&
            pressure=pressure,pmswap=pmswap,pmvolume=pmvolume,&
            pmtraion=pmtraion,pmtrans=pmtrans,BETA=BETA,rcut=rcut,&
            lrestart=lrestart)

! check for restart
         CALL get_mc_par(mc_par_2,lrestart=lrestart,dat_file=dat_file)
         IF ( lrestart ) THEN
            CALL read_mc_restart(mc_par_2,force_env_2,iw,rng_stream)
            CALL force_env_release(force_env_2)
            globenv%input_file_name=dat_file
            CALL mc_create_force_env(force_env_2, globenv, .TRUE., &
               force_env%in_use,error=error)
         ENDIF

! pull all the information together into one environment
         CALL set_mc_env( mc_env(1)%mc_env,mc_par=mc_par,&
            force_env=force_env)
         CALL set_mc_env( mc_env(2)%mc_env,mc_par=mc_par_2,&
            force_env=force_env_2)

! do the calculation
         CALL mc_run_ensemble(mc_env,globenv,nboxes,error,rng_stream)

! clean up after ourselves
         CALL get_mc_env( mc_env(2)%mc_env,mc_par=mc_par_2,&
            force_env=force_env_2)

         CALL get_mc_par(mc_par_2,box2_file=box2_file)
         globenv%input_file_name=box2_file
         CALL mc_env_release ( mc_env(2)%mc_env, error )
         CALL force_env_release(force_env_2)
         DEALLOCATE ( mc_par_2, STAT = isos )
         CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)

      CASE DEFAULT
         IF (ionode) THEN
            WRITE(6,*) "No idea how we got here...input section ",&
            "should've caught this: ",ensemble
         ENDIF

      END SELECT

      CALL get_mc_env( mc_env(1)%mc_env, &
         mc_par = mc_par,force_env=force_env)
      DEALLOCATE ( mc_par, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)
      CALL mc_env_release ( mc_env(1)%mc_env, error )
      DEALLOCATE ( mc_env, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)

! delete the random numbers
      CALL delete_rng_stream(rng_stream)

END SUBROUTINE do_mon_car


!******************************************************************************


END MODULE mc_run
