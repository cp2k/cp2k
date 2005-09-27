!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** mc_control.F *********************************************************
!!
!!   NAME
!!     mc_control
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
#include "cp_prep_globals.h"

MODULE mc_control
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_files,                        ONLY: open_file,&
                                             close_file
  USE cp_log_handling,                 ONLY: cp_failure_level, cp_to_string
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             default_string_length
  USE mathconstants,                   ONLY: pi
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_types,                        ONLY: mc_simulation_parameters_type,&
                                             get_mc_par,&
                                             set_mc_par
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             atom_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_control"

  PUBLIC :: write_mc_restart,read_mc_restart,mc_create_force_env,&
      mc_create_bias_force_env


!****************************************************************************

CONTAINS

! *****************************************************************************
!!****s* mc_control.F/write_mc_restart
!!
!!   NAME
!!     write_mc_restart
!!
!!   FUNCTION
!!     writes the coordinates of the current step to a file that can
!!     be read in at the start of the next simulation
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment to write the coords from
!!     - nnstep: how many steps the simulation has run
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE write_mc_restart ( nnstep, mc_par, force_env )

    INTEGER, INTENT(IN)                      :: nnstep
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(LEN=*), PARAMETER    &
      :: routineN = "write_mc_restart"

    INTEGER                                  :: handle, unit,&
      nunits_tot,nunits,ichain,iunit,nchain
    REAL(KIND=dp),DIMENSION(1:3)             :: abc
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(atom_type), DIMENSION(:),POINTER    :: atom_list
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    CHARACTER(LEN=40)                        :: restart_file
    CHARACTER(LEN=20)                        :: ensemble
    REAL(KIND=dp)                            :: temp

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! get some data from mc_par
  CALL get_mc_par(mc_par,restart_file=restart_file,temp=temp,&
      nunits_tot=nunits_tot,nchain=nchain,ensemble=ensemble,&
      nunits=nunits)

! open the file and write some simulation parameters
  CALL open_file(file_name=restart_file,unit_number=unit,&
                 file_action='WRITE',file_status='REPLACE')
  WRITE(unit,*) nnstep
  WRITE(unit,*) temp,nunits_tot
  WRITE(unit,*) ensemble

! get the cell length and coordinates
  CALL force_env_get(force_env,cell=cell,subsys=subsys,&
                         error=error)
  CALL get_cell(cell,abc=abc)
  CALL cp_subsys_get(subsys(1)%subsys, &
             molecule_kinds_new=molecule_kinds_new,&
             particles=particles, error=error)
  molecule_kind => molecule_kinds_new%els(1)
  CALL get_molecule_kind(molecule_kind,atom_list=atom_list)

  WRITE(unit,'(3(F10.6,3X))') abc(1:3)*angstrom  ! in angstroms
  WRITE(unit,*)

! write the coordinates out
  DO ichain=1,nchain
      DO iunit=1,nunits
         WRITE(unit,'(1X,A,1X,3(F15.10,3X))') &
                 TRIM(ADJUSTL(atom_list(iunit)%name)),&
                 particles%els((ichain-1)*nunits+iunit)%r(1:3)*angstrom
      ENDDO
  ENDDO
  
  CALL close_file(unit_number=unit)

! end the timing
  CALL timestop(0.0E0_dp,handle)


END SUBROUTINE write_mc_restart

! *****************************************************************************
!!****s* mc_control.F/read_mc_restart
!!
!!   NAME
!!     read_mc_restart
!!
!!   FUNCTION
!!     reads the input coordinates of the simulation from a file written above
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment to write the coords from
!!     - iw: the unit to write an error message to, in case current
!!           simulation parameters don't match what's in the restart file
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE read_mc_restart ( mc_par, force_env ,iw )

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(LEN=*), PARAMETER            &
      :: routineN = "read_mc_restart"

    CHARACTER(LEN=40)                        :: restart_file,dat_file
    CHARACTER(LEN=200)                       :: fft_lib
    CHARACTER(LEN=20)                        :: chardum, ensemble,&
      mc_ensemble
    INTEGER                                  :: handle, ipart, istat, &
                                                nunits_tot, unit,mc_nunits,&
                                                mc_nunits_tot,mc_nchain,&
                                                nstart,i
    REAL(KIND=dp)                            :: temp,mc_temp,rand
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp),DIMENSION(1:3)             :: abc,box_length
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: ionode

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from the mc_par
      CALL get_mc_par(mc_par,restart_file=restart_file,temp=mc_temp,&
      nunits_tot=mc_nunits_tot,ensemble=mc_ensemble,nunits=mc_nunits,&
      ionode=ionode,dat_file=dat_file,nchain=mc_nchain,fft_lib=fft_lib)

! open the file and read some simulation parameters
      CALL open_file(file_name=restart_file,unit_number=unit,&
         file_action='READ',file_status='OLD')

      READ(unit,*) nstart
      READ(unit,*) temp,nunits_tot
      READ(unit,*) ensemble

! do some checking
      IF(ABS(temp - mc_temp) .GT. 0.01E0_dp) THEN
         IF(ionode) THEN
            WRITE(iw,*) 'The temperature in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file temp =',mc_temp
            WRITE(iw,*) 'Restart file temp =',temp
         ENDIF
         STOP
      ENDIF
      IF(nunits_tot .NE. mc_nunits_tot) THEN
         IF(ionode) THEN
            WRITE(iw,*) 'The total number of units in the restart ',&
            'file is not the same as the input file.'
            WRITE(iw,*) 'Input file units =',mc_nunits_tot
            WRITE(iw,*) 'Restart file units =',nunits_tot
         ENDIF
         mc_nunits_tot=nunits_tot
         mc_nchain=mc_nunits_tot/mc_nunits
!     STOP
      ENDIF
      IF(ensemble .NE. mc_ensemble) THEN
         IF(ionode) THEN
            WRITE(iw,*) 'The ensemble in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file ensemble =',mc_ensemble
            WRITE(iw,*) 'Restart file ensemble =',ensemble
         ENDIF
         STOP
      ENDIF

! get the cell length and coordinates
      CALL force_env_get(force_env,cell=cell,subsys=subsys,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

      READ(unit,*)  box_length(1:3)  ! in angstroms
      IF(ABS(box_length(1) - abc(1)*angstrom) .GT. 0.0001E0_dp .OR. &
         ABS(box_length(2) - abc(2)*angstrom) .GT. 0.0001E0_dp .OR. &
         ABS(box_length(3) - abc(3)*angstrom) .GT. 0.0001E0_dp) THEN
         IF (ionode) THEN
            WRITE(iw,*) 'The cell length in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file cell length =',abc(1:3)*angstrom
            WRITE(iw,*) 'Restart file cell length =',box_length(1:3)
         ENDIF
!     STOP
      ENDIF
      
      READ(unit,*)
! allocate the array holding the coordiantes, and read in the coordiantes,
! and write the dat file so we can make a new force_env
      IF(nunits_tot == 0) THEN
         ALLOCATE (r(3,3),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r",3*3*dp_size)
         r(1:3,1) = (/ 1.088738E0_dp, 1.8323472E0_dp,2.0179184E0_dp /)
         r(1:3,2) = (/ 0.337575E0_dp, 3.0485244E0_dp, 0.80660535E0_dp /)
         r(1:3,3) = (/ 0.8812212E0_dp, 2.615469E0_dp, 3.707380E0_dp /)

         CALL close_file(unit_number=unit)

         IF(ionode) &
            CALL mc_make_dat_file(r(:,:),3,box_length(:),dat_file,&
               force_env,TRIM(fft_lib))

      ELSE
         ALLOCATE (r(3,nunits_tot),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r",3*nunits_tot*dp_size)
         DO ipart=1,mc_nunits_tot
            READ(unit,*) chardum,r(1:3,ipart)
            r(1:3,ipart)=r(1:3,ipart)/angstrom
         ENDDO

         CALL close_file(unit_number=unit)

         IF(ionode) &
            CALL mc_make_dat_file(r(:,:),mc_nunits_tot,&
            box_length(:),dat_file,force_env,TRIM(fft_lib))

      ENDIF

      CALL set_mc_par(mc_par,nunits_tot=mc_nunits_tot,nunits=mc_nunits,&
         nstart=nstart,nchain=mc_nchain)
      
! advance the random number sequence based on the restart step
      DO i=1,nstart+1
         CALL random_number(rand)
      ENDDO

! end the timing
      CALL timestop(0.0E0_dp,handle)

! deallcoate
      DEALLOCATE(r)

END SUBROUTINE read_mc_restart

! *****************************************************************************
!!****s* mc_control.F/mc_create_force_env
!!
!!   NAME
!!     mc_create_force_env
!!
!!   FUNCTION
!!     creates a force environment for any of the different kinds of
!!     simulation
!!
!!   ARGUMENTS
!!     - force_env: the force environment to create
!!     - globenv: the global environment parameters
!!     - box_change: is this a volume or swap move?
!!     - error: the CP error structure
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_create_force_env ( force_env, globenv, box_change, in_use, &
      error )

      LOGICAL, INTENT(IN)                      :: box_change
      INTEGER, INTENT(IN)                      :: in_use
      TYPE(force_env_type), POINTER            :: force_env
      TYPE(global_environment_type), POINTER   :: globenv
      TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
      
      CHARACTER(LEN=*), PARAMETER    &
         :: routineN = "mc_create_force_env",&
            routineP = module_name//':'//routineN
      LOGICAL :: failure
      
      SELECT CASE ( in_use )
      CASE ( use_fist_force )
         CALL fist_create_force_env(force_env, globenv, error=error)
         force_env%box_changed=.true.
      CASE ( use_qs_force )
         CALL quickstep_create_force_env(force_env, globenv, error=error)
      CASE ( use_kg_force )
      CASE default
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

END SUBROUTINE mc_create_force_env

! *****************************************************************************
!!****s* mc_control.F/mc_create_force_env
!!
!!   NAME
!!     mc_create_force_env
!!
!!   FUNCTION
!!     creates a force environment for any of the different kinds of
!!     simulation
!!
!!   ARGUMENTS
!!     - force_env: the force environment to create
!!     - globenv: the global environment parameters
!!     - box_change: is this a volume or swap move?
!!     - error: the CP error structure
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_create_bias_force_env (bias_env,force_env,&
            globenv,box_change,fft_lib,error)

      LOGICAL, INTENT(IN)                      :: box_change
      TYPE(force_env_type), POINTER            :: force_env
      TYPE(force_env_type), POINTER            :: bias_env
      TYPE(global_environment_type), POINTER   :: globenv
      TYPE(cp_error_type), INTENT(INOUT)       :: error
      CHARACTER(LEN=*),INTENT(IN)              :: fft_lib
      
      TYPE(cp_subsystem_p_type), &
         DIMENSION(:), POINTER                 :: subsys
      TYPE(cell_type), POINTER                 :: cell
      TYPE(particle_list_type), POINTER        :: particles
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :)                       :: r
      CHARACTER(LEN=*), PARAMETER    &
         :: routineN = "mc_create_force_env",&
            routineP = module_name//':'//routineN
      LOGICAL :: failure
      INTEGER                                  :: nunits_tot,iparticle,&
                istat
      REAL(KIND=dp), DIMENSION(1:3)        :: abc
                   
!-----------------------------------------------------------------------------
      NULLIFY(cell,subsys)

      globenv%input_file_name='bias_template.inp'
!!!!!!!!!!!!!!!!!!!!!! change to be all bias types
      CALL mc_create_force_env(bias_env, globenv, box_change, &
         use_fist_force,error=error)

      CALL force_env_get(force_env,subsys=subsys,cell=cell,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

      nunits_tot=SIZE(particles%els(:))
      ALLOCATE (r(3,nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "particle_positions",3*nunits_tot*dp_size)

      DO iparticle=1,nunits_tot
         r(1:3,iparticle)=particles%els(iparticle)%r(1:3)
      ENDDO

      CALL mc_make_dat_file(r(:,:),nunits_tot,&
         abc(1:3)*angstrom,&
         'bias_temp.inp',bias_env,TRIM(fft_lib))

      CALL force_env_release(bias_env,error=error)
      globenv%input_file_name='bias_temp.inp'
      CALL mc_create_force_env(bias_env, globenv, box_change, &
         use_fist_force,error=error)

      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                       "r")

END SUBROUTINE mc_create_bias_force_env

END MODULE mc_control

!******************************************************************************
