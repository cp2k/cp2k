!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_type,&
                                             force_env_p_type,&
                                             use_fist_force
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_release
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_val_get,&
                                             section_vals_get,&
                                             section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mc_control,                      ONLY: mc_create_force_env,&
                                             read_mc_restart
  USE mc_ensembles,                    ONLY: mc_run_ensemble
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_env_create,&
                                             mc_env_release,&
                                             mc_environment_p_type,&
                                             set_mc_env
  USE mc_types,                        ONLY: find_mc_rcut,&
                                             get_mc_par,&
                                             mc_simulation_parameters_p_type,&
                                             read_mc_section,&
                                             set_mc_par,&
                                             get_mc_molecule_info,&
                                             mc_determine_molecule_info,&
                                             mc_sim_par_create,&
                                             mc_sim_par_destroy,&
                                             mc_molecule_info_type,&
                                             mc_molecule_info_destroy,&
                                             mc_input_file_type,&
                                             mc_input_file_create,&
                                             mc_input_file_destroy,&
                                             mc_input_parameters_check
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             rng_stream_type
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_run'

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
!!   FUNCTION
!!     starts the Monte Carlo simulation and determines which ensemble we're
!!     running
!!
!!   ARGUMENTS
!!     - globenv: the global environment for the simulation
!!     - force_env_1: the force environment for the simulation, or
!!                  the force environment for box 1, depending on which
!!                  ensemble we're running
!!     - input_file_name: the name of the input file that force_env_1 was
!!       created from
!!     - error: the error type
!!
!!     Designed for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************

  SUBROUTINE do_mon_car ( force_env_1, globenv, input_file_name, error )


    TYPE(force_env_type), POINTER            :: force_env_1
    TYPE(force_env_p_type), DIMENSION(:),POINTER &
         :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*)                         :: input_file_name
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_mon_car', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: ensemble
    CHARACTER(LEN=40)                        :: box2_file, dat_file
    INTEGER :: ibox, isos, iw, nboxes, &
      istat, nmol_types,box_number
    LOGICAL                                  :: ionode, lrestart, &
                                                lterminate,lbias,&
                                                lskip_box
    REAL(dp)                                 :: rcut
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(force_env_type), POINTER            :: force_env_2
    TYPE(global_environment_type), POINTER   :: globenv_2
    TYPE(mc_environment_p_type), &
      DIMENSION(:), POINTER                  :: mc_env
    TYPE(mc_simulation_parameters_p_type), &
      DIMENSION(:),POINTER                                :: mc_par
    TYPE(mc_molecule_info_type), POINTER        :: mc_molecule_info
    TYPE(mc_input_file_type), POINTER        :: mc_input_file,&
         mc_bias_file
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                root_section,&
                                                mc_section
    INTEGER,DIMENSION(:),POINTER             :: nunits_tot
    REAL(dp),DIMENSION(:,:),POINTER          :: empty_coordinates
!************************************************************************

      NULLIFY(mc_env,mc_par,force_env_2,rng_stream,force_env_section,&
         root_section,mc_molecule_info)
      para_env=>force_env_1%para_env
      logger=>cp_error_get_logger(error)
      iw = cp_logger_get_default_unit_nr(logger)
      force_env_section => force_env_1%force_env_section
      root_section => force_env_1%root_section
      CALL section_vals_get(force_env_section,n_repetition=isos,error=error)
      CPPostconditionNoFail(isos==1,cp_warning_level,routineP,error)
! set some values...will use get_globenv if that ever comes around
      ionode=para_env%ionode

! initialize the random numbers
      CALL create_rng_stream(rng_stream=rng_stream,&
                             name="first",&
                             distribution_type=UNIFORM,error=error)

! need to figure out how many boxes we have, based on the value
! of mc_par % ensemble
      NULLIFY(mc_section)
      mc_section => section_vals_get_subs_vals(root_section,&
         "MOTION%MC",error=error)
      CALL section_vals_val_get(mc_section,"ENSEMBLE",&
           c_val=ensemble,error=error)

! now we read in the second force_env, if we have another box
      SELECT CASE(ensemble)
      CASE("TRADITIONAL")
         nboxes=1
      CASE("GEMC_NVT","GEMC_NPT")
         nboxes=2
         CALL section_vals_val_get(mc_section,"BOX2_FILE_NAME",&
              c_val=box2_file,error=error)
         CALL mc_create_force_env(force_env_2, para_env,&
              box2_file,globenv_new=globenv_2,error=error)
         CALL section_vals_release(force_env_2%root_section,error=error)
      END SELECT

! now we create the various pointers that contain information for all boxes
      ALLOCATE (force_env(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "force_env",nboxes)
      SELECT CASE(ensemble)
      CASE("TRADITIONAL")
         force_env(1)%force_env => force_env_1
      CASE("GEMC_NVT","GEMC_NPT")
         force_env(1)%force_env => force_env_1
         force_env(2)%force_env => force_env_2
      END SELECT      
      ALLOCATE (mc_par(1:nboxes),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "mc_par",nboxes)
      ALLOCATE(mc_env(1:nboxes),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory('allocstruc','mc_env',nboxes)

! now we need the molecule information
! determine the total number of molecules and atoms
      CALL mc_determine_molecule_info(force_env,mc_molecule_info,error,&
           coordinates_empty=empty_coordinates)
      CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types,&
           nunits_tot=nunits_tot)

      DO ibox=1,nboxes

         IF (ionode) THEN
            WRITE(iw,*)
            WRITE(iw,*)
            WRITE(iw,'(A,I2,A)') '******************************** Begin BOX ',ibox,&
            '  *******************************'
         ENDIF

! allocates an mc_env and sets the variables to zero
         CALL mc_env_create ( mc_env(ibox)%mc_env, &
            error = error )

! now read in the values of the mc_pars
! creating the mc_par
         CALL mc_sim_par_create(mc_par(ibox)%mc_par,nmol_types,nboxes)

! attach molecule information to all mc_par structures, so we know what we have to read in
         CALL set_mc_par(mc_par(ibox)%mc_par,mc_molecule_info=mc_molecule_info)

! read the input of the Monte Carlo section
         force_env_section => force_env(ibox)%force_env%force_env_section
         root_section => force_env(ibox)%force_env%root_section
         IF(ibox == 1) THEN
            CALL read_mc_section ( mc_par(ibox)%mc_par, para_env, globenv ,input_file_name,&
                 root_section,force_env_section,error)
         ELSE
            CALL read_mc_section ( mc_par(ibox)%mc_par, para_env, globenv ,box2_file,&
                 root_section,force_env_section,error)
         ENDIF

! get the input file data, in case we need to make a restart...
! this is also used in the swap move, or anytime we need to make a dat file
! always judge based on lbias from box 1...in case someone only changes the value
! for box 1
         CALL get_mc_par(mc_par(ibox)%mc_par,mc_input_file=mc_input_file)
         CALL get_mc_par(mc_par(1)%mc_par,lbias=lbias)
         IF(ibox == 1) THEN
            CALL mc_input_file_create ( mc_input_file,&
                 input_file_name,mc_molecule_info,empty_coordinates,error)
         ELSE
            CALL mc_input_file_create ( mc_input_file,&
                 box2_file,mc_molecule_info,empty_coordinates,error)
         ENDIF
         CALL set_mc_par(mc_par(ibox)%mc_par,mc_input_file=mc_input_file)

         IF(lbias) THEN
            CALL get_mc_par(mc_par(ibox)%mc_par,mc_bias_file=mc_bias_file)
            CALL mc_input_file_create ( mc_bias_file,&
                 "bias_template.inp",mc_molecule_info,empty_coordinates,error)
            CALL set_mc_par(mc_par(ibox)%mc_par,mc_bias_file=mc_bias_file)
         ENDIF

! check for restart
         CALL get_mc_par(mc_par(ibox)%mc_par,lrestart=lrestart,&
              dat_file=dat_file)
         IF ( lrestart ) THEN
            CALL read_mc_restart(mc_par(ibox)%mc_par,force_env(ibox)%force_env,&
                 iw,nunits_tot(ibox),rng_stream)
! release the old force env and make the new one
            CALL force_env_release(force_env(ibox)%force_env,error=error)
            CALL mc_create_force_env(force_env(ibox)%force_env, &
                 para_env, dat_file,error=error)
         ENDIF

      ENDDO

! figure out if we have an empty box
      box_number=0
      lskip_box=.FALSE.
      DO ibox=1,nboxes
         IF(nunits_tot(ibox) == 0) THEN
            IF(lskip_box) THEN
               CALL stop_program('mc_run.F/do_mon_car',&
               'More than one box has no atoms in it!')
            END IF
            box_number=ibox
            lskip_box=.TRUE.
         ENDIF
      ENDDO

! in case there was a restart, we need to do this again
      CALL mc_molecule_info_destroy(mc_molecule_info)
      CALL mc_determine_molecule_info(force_env,mc_molecule_info,error,box_number=box_number)
      CALL get_mc_molecule_info(mc_molecule_info,nmol_types=nmol_types,&
           nunits_tot=nunits_tot)
      DO ibox=1,nboxes
         CALL set_mc_par(mc_par(ibox)%mc_par,mc_molecule_info=mc_molecule_info)
      ENDDO

! if we're doing a classical simulation, figure out the largest
! potential cutoff and write it to the screen
      IF(force_env(1)%force_env%in_use .EQ. use_fist_force) THEN
         CALL find_mc_rcut(mc_par(1)%mc_par,force_env(1)%force_env,lterminate)
         CALL get_mc_par(mc_par(1)%mc_par,rcut=rcut)
         IF(ionode) WRITE (iw,'( A,T73,F8.4 )') &
            ' MC| Interaction cutoff [angstroms]',rcut*angstrom
         IF(lterminate) THEN
            CALL stop_program('mc_run.F/do_mon_car',&
               'Cutoff larger than twice the boxlength')
         ENDIF
      ENDIF

! make sure some values are the same between boxes
      IF(nboxes == 2) THEN
         CALL equilize_mc_sim_parameters(mc_par,iw,ionode)
      ENDIF

! now check the input parameters and run the simulation
      DO ibox=1,nboxes
         
         CALL mc_input_parameters_check(mc_par(ibox)%mc_par)
         
! attach all the structures to one convientent structure
       CALL set_mc_env( mc_env(ibox)%mc_env,mc_par=mc_par(ibox)%mc_par,&
            force_env=force_env(ibox)%force_env)

      ENDDO
      CALL mc_run_ensemble(mc_env,para_env, globenv,nboxes,error,&
           rng_stream)

! get rid of all the MC molecule information
      CALL get_mc_env( mc_env(1)%mc_env, mc_par = mc_par(1)%mc_par)
      CALL get_mc_par(mc_par(1)%mc_par,mc_molecule_info=mc_molecule_info)
      CALL mc_molecule_info_destroy(mc_molecule_info)
      
      DO ibox=1,nboxes
         CALL get_mc_env( mc_env(ibox)%mc_env, &
              mc_par = mc_par(ibox)%mc_par,force_env=force_env(ibox)%force_env)
         CALL get_mc_par(mc_par(ibox)%mc_par,mc_input_file=mc_input_file)
         
         CALL mc_input_file_destroy(mc_input_file)
         IF(lbias)THEN
            CALL get_mc_par(mc_par(ibox)%mc_par,mc_bias_file=mc_bias_file)
            CALL mc_input_file_destroy(mc_bias_file)
         ENDIF

         CALL mc_sim_par_destroy(mc_par(ibox)%mc_par)
         CALL mc_env_release ( mc_env(ibox)%mc_env, error )

         IF (ibox == 2) THEN
            CALL force_env_release(force_env(ibox)%force_env,error=error)           
            CALL globenv_release(globenv_2,para_env, error)
         ELSE
            force_env_1 => force_env(1)%force_env
         ENDIF
      ENDDO

      DEALLOCATE ( empty_coordinates, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)
      DEALLOCATE ( mc_par, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)
      DEALLOCATE ( mc_env, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)
      DEALLOCATE ( force_env, STAT = isos )
      CPPostconditionNoFail(isos==0,cp_warning_level,routineP,error)

! delete the random numbers
      CALL delete_rng_stream(rng_stream,error=error)


    END SUBROUTINE do_mon_car

! *****************************************************************************
!!****s* mc_run.F/equilize_mc_sim_parameters
!!
!!   NAME
!!     equilize_mc_sim_parameters
!!
!!   FUNCTION
!!     takes some parameters from one set of MC simulation parameters
!!     and transfers them to a second set...used so that we're not using
!!     different parameters between two simulation boxes, if they should
!!     be the same (move probabilities, for instance)
!!
!!   ARGUMENTS
!!     - mc_par: the pointer that contains the simulation parameters
!!          for both boxes
!!     - iw: the unit number that prints to screen
!!     - ionode: equal to .TRUE. if this is the main MPI node
!!
!!     Designed for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!**************************************************************************
    SUBROUTINE equilize_mc_sim_parameters(mc_par,iw,ionode)
      INTEGER,INTENT(IN) :: iw
      LOGICAL,INTENT(IN) :: ionode
      TYPE(mc_simulation_parameters_p_type),DIMENSION(:),POINTER :: mc_par

      INTEGER :: nstep,iupvolume,iuptrans,nmoves,nswapmoves,iprint
      REAL(dp) :: temp,pressure,pmswap,pmvolume,pmtraion,pmtrans,BETA,&
           rcut
      LOGICAL :: lbias,lstop,lrestart
      INTEGER,DIMENSION(:),POINTER :: avbmc_atom
      REAL(dp),DIMENSION(:),POINTER :: avbmc_rmin,avbmc_rmax,pbias,pmswap_mol,&
           pmtrans_mol,pmrot_mol,pmtraion_mol,pmavbmc_mol
      CHARACTER(20) :: ensemble
!***********************************************************************

      IF(ionode) THEN
         WRITE(iw,'( A,A )' ) 'Ignoring some input for box 2, and ',&
              'using the values for box 1 for the following variables:'
         WRITE(iw,'( A,A )' ) 'nstep,iupvolume,iuptrans,nmoves,',&
              'nswapmoves,iprint,lbias,lstop,temp,pressure'
         WRITE(iw,'( A,A )' ) 'pmvolume,pmtraion,pmtrans,BETA,rcut,',&
              'lrestart'
         WRITE(iw,'( A,A )' ) 'pmtraion_mol,pmtrans_mol,pmrot_mol,pmswap_mol,',&
              'avbmc_atom'
         WRITE(iw,'( A,A )' ) 'avbmc_rmin,avmbc_rmax,pmavbmc_mol,pbias,pmswap'
      ENDIF
      
      CALL get_mc_par(mc_par(1)%mc_par,nstep=nstep,iupvolume=iupvolume,&
           iuptrans=iuptrans,nmoves=nmoves,nswapmoves=nswapmoves,&
           iprint=iprint,lbias=lbias,lstop=lstop,temp=temp,&
           pressure=pressure,pmswap=pmswap,pmvolume=pmvolume,&
           pmtraion=pmtraion,pmtrans=pmtrans,BETA=BETA,rcut=rcut,&
           lrestart=lrestart,pmtraion_mol=pmtraion_mol,pmtrans_mol=pmtrans_mol,&
           pmrot_mol=pmrot_mol,pmswap_mol=pmswap_mol,avbmc_atom=avbmc_atom,&
           avbmc_rmin=avbmc_rmin,avbmc_rmax=avbmc_rmax,pmavbmc_mol=pmavbmc_mol,&
           pbias=pbias,ensemble=ensemble)
      CALL set_mc_par(mc_par(2)%mc_par,nstep=nstep,iupvolume=iupvolume,&
           iuptrans=iuptrans,nmoves=nmoves,nswapmoves=nswapmoves,&
           iprint=iprint,lbias=lbias,lstop=lstop,temp=temp,&
           pressure=pressure,pmswap=pmswap,pmvolume=pmvolume,&
           pmtraion=pmtraion,pmtrans=pmtrans,BETA=BETA,rcut=rcut,&
           lrestart=lrestart,pmtraion_mol=pmtraion_mol,pmtrans_mol=pmtrans_mol,&
           pmrot_mol=pmrot_mol,pmswap_mol=pmswap_mol,avbmc_atom=avbmc_atom,&
           avbmc_rmin=avbmc_rmin,avbmc_rmax=avbmc_rmax,pmavbmc_mol=pmavbmc_mol,&
           pbias=pbias,ensemble=ensemble)

    END SUBROUTINE equilize_mc_sim_parameters

    
END MODULE mc_run
