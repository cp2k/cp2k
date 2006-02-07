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
!!     contains some general routines for dealing with the restart
!!     files and creating force_env for MC use
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_control
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force
  USE global_types,                    ONLY: global_environment_type
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: section_release,&
                                             section_type,&
                                             section_vals_create,&
                                             section_vals_release,&
                                             section_vals_type
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_simulation_parameters_type,&
                                             set_mc_par
  USE message_passing,                 ONLY: mp_bcast
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mc_control'

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
!!   Only use in serial.
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

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'write_mc_restart'

    CHARACTER(LEN=20)                        :: ensemble
    CHARACTER(LEN=40)                        :: restart_file_name
    INTEGER                                  :: handle, ichain, iunit, &
                                                nchain, nunits, nunits_tot, &
                                                unit
    REAL(KIND=dp)                            :: temp
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! get some data from mc_par
  CALL get_mc_par(mc_par,restart_file_name=restart_file_name,temp=temp,&
      nunits_tot=nunits_tot,nchain=nchain,ensemble=ensemble,&
      nunits=nunits)

! open the file and write some simulation parameters
  CALL open_file(file_name=restart_file_name,unit_number=unit,&
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
!!     - rng_stream: the stream we pull random numbers from
!!
!!     Used in parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE read_mc_restart ( mc_par, force_env, iw, rng_stream)

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: iw
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'read_mc_restart'

    CHARACTER(LEN=20)                        :: chardum, ensemble, mc_ensemble
    CHARACTER(LEN=200)                       :: fft_lib
    CHARACTER(LEN=40)                        :: dat_file, restart_file_name
    INTEGER :: group, handle, i, ipart, istat, iunit, mc_nchain, mc_nunits, &
      mc_nunits_tot, nstart, nunits_tot, source, unit
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: mc_temp, rand, temp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp), DIMENSION(1:3)            :: abc, box_length
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! get some stuff from the mc_par
      CALL get_mc_par(mc_par,restart_file_name=restart_file_name,temp=mc_temp,&
         nunits_tot=mc_nunits_tot,ensemble=mc_ensemble,nunits=mc_nunits,&
         ionode=ionode,dat_file=dat_file,nchain=mc_nchain,fft_lib=fft_lib,&
         group=group,source=source)

      IF(ionode) THEN
! open the file and read some simulation parameters
         CALL open_file(file_name=restart_file_name,unit_number=unit,&
            file_action='READ',file_status='OLD')

         READ(unit,*) nstart
         READ(unit,*) temp,nunits_tot
         READ(unit,*) ensemble
      ENDIF
      CALL mp_bcast(nstart,source,group)
      CALL mp_bcast(temp,source,group)
      CALL mp_bcast(nunits_tot,source,group)
      CALL mp_bcast(ensemble,source,group)

! do some checking
      IF(ABS(temp - mc_temp) .GT. 0.01E0_dp) THEN
         IF(ionode) THEN
            WRITE(iw,*) 'The temperature in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file temp =',mc_temp
            WRITE(iw,*) 'Restart file temp =',temp
         ENDIF
         CALL stop_program(routineN,moduleN,__LINE__,&
            "Temperature difference between restart and input")
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
      ENDIF
      IF(ensemble .NE. mc_ensemble) THEN
         IF(ionode) THEN
            WRITE(iw,*) 'The ensemble in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file ensemble =',mc_ensemble
            WRITE(iw,*) 'Restart file ensemble =',ensemble
         ENDIF
         CALL stop_program(routineN,moduleN,__LINE__,&
            "Ensembles different between restart and input")
      ENDIF

! get the cell length and coordinates
      CALL force_env_get(force_env,cell=cell,subsys=subsys,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

      IF(ionode) THEN
         READ(unit,*)  box_length(1:3) ! in angstroms
         READ(unit,*)
      ENDIF
      CALL mp_bcast(box_length,source,group)
      IF(ABS(box_length(1) - abc(1)*angstrom) .GT. 0.0001E0_dp .OR. &
         ABS(box_length(2) - abc(2)*angstrom) .GT. 0.0001E0_dp .OR. &
         ABS(box_length(3) - abc(3)*angstrom) .GT. 0.0001E0_dp) THEN
         IF (ionode) THEN
            WRITE(iw,*) 'The cell length in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file cell length =',abc(1:3)*angstrom
            WRITE(iw,*) 'Restart file cell length =',box_length(1:3)
         ENDIF
      ENDIF
      
! allocate the array holding the coordiantes, and read in the coordiantes,
! and write the dat file so we can make a new force_env
      IF(nunits_tot == 0) THEN
         ALLOCATE (r(3,mc_nunits),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "r",3*mc_nunits*dp_size)
         DO iunit=1,mc_nunits
            r(1:3,iunit) = (/ REAL(iunit,dp), REAL(iunit,dp),REAL(iunit,dp) /)
         ENDDO


         IF(ionode) THEN
            CALL mc_make_dat_file(r(:,:),mc_nunits,box_length(:),dat_file,&
               force_env,TRIM(fft_lib))
            CALL close_file(unit_number=unit)
         ENDIF
      ELSE
         ALLOCATE (r(3,nunits_tot),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "r",3*nunits_tot*dp_size)
         
         IF(ionode) THEN
            DO ipart=1,mc_nunits_tot
               READ(unit,*) chardum,r(1:3,ipart)
               r(1:3,ipart)=r(1:3,ipart)/angstrom
            ENDDO

            CALL close_file(unit_number=unit)

            CALL mc_make_dat_file(r(:,:),mc_nunits_tot,&
               box_length(:),dat_file,force_env,TRIM(fft_lib))
         ENDIF
      ENDIF

      CALL set_mc_par(mc_par,nunits_tot=mc_nunits_tot,nunits=mc_nunits,&
         nstart=nstart,nchain=mc_nchain)
      
! advance the random number sequence based on the restart step
      IF(ionode) THEN
         DO i=1,nstart+1
!         CALL random_number(rand)
            rand=next_random_number(rng_stream)
         ENDDO
      ENDIF

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
!!     MC simulations we can do (FIST, QS, KG)
!!
!!   ARGUMENTS
!!     - force_env: the force environment to create
!!     - globenv: the global environment parameters
!!     - box_change: is this a volume or swap move?
!!     - in_use: an integer that tells what kind of environment it is
!!     - error: the CP error structure
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_create_force_env ( force_env, globenv, box_change, in_use, &
      error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN)                      :: box_change
    INTEGER, INTENT(IN)                      :: in_use
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_create_force_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_parser_type), POINTER            :: my_parser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(section_type), POINTER              :: input_structure
    TYPE(section_vals_type), POINTER         :: file_temp, input_file

!***************************************************************************

      NULLIFY(file_temp,input_file,input_structure,my_parser,default_units)

    ! parse the new input file input
      CALL create_cp2k_root_section(input_structure,error=error)
      CALL section_vals_create(input_file,input_structure, error=error)
      CALL section_release(input_structure,error=error)
      CALL parser_create(my_parser,file_name=globenv%input_file_name,&
         error=error)
      CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
      CALL section_vals_parse(input_file,my_parser,root_section=.FALSE.,&
      default_units=default_units,error=error)
      CALL cp_unit_set_release(default_units,error=error)
      CALL parser_release(my_parser,error=error)
      CALL cp_error_check(error,failure)
      file_temp => globenv%input_file 
      globenv%input_file => input_file

      SELECT CASE ( in_use )
      CASE ( use_fist_force )
         CALL fist_create_force_env(force_env, globenv, error=error)
         force_env%box_changed=.TRUE.
      CASE ( use_qs_force )
         CALL quickstep_create_force_env(force_env, globenv, error=error)
      CASE ( use_kg_force )
         CALL kg_create_force_env ( force_env, globenv, error=error)
      CASE default
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

! release some stuff created here
      globenv%input_file => file_temp
      CALL section_vals_release(input_file)     
      NULLIFY(file_temp)


END SUBROUTINE mc_create_force_env

! *****************************************************************************
!!****s* mc_control.F/mc_create_bias_force_env
!!
!!   NAME
!!     mc_create_bias_force_env
!!
!!   FUNCTION
!!     essentially copies the cell size and coordinates of one force env
!!     to another that we will use to bias some moves with
!!
!!   ARGUMENTS
!!     - bias_env: the force environment to create
!!     - force_env: the force environment we're copying cell length
!!           and coordinates from
!!     - globenv: the global environment parameters
!!     - box_change: is this a volume or swap move?
!!     - fft_lib: the FFT library we're using (e.g. FFTW, FFTSG)
!!     - error: the CP error structure
!!
!!   Suitable for parallel.
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE mc_create_bias_force_env (bias_env,force_env,&
            globenv,box_change,fft_lib,ionode,error)

    TYPE(force_env_type), POINTER            :: bias_env, force_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN)                      :: box_change
    CHARACTER(LEN=*), INTENT(IN)             :: fft_lib
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_create_bias_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, istat, nunits_tot
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp), DIMENSION(1:3)            :: abc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!-----------------------------------------------------------------------------

      NULLIFY(cell,subsys)

! create the bias_env from the template file
      globenv%input_file_name='bias_template.inp'
      CALL mc_create_force_env(bias_env, globenv, box_change, &
         use_fist_force,error=error)

! grab the coordinates and cell parameters of the old environment
      CALL force_env_get(force_env,subsys=subsys,cell=cell,&
         error=error)
      CALL get_cell(cell,abc=abc)
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

!  copy the data we need and make a new force_environment
      nunits_tot=SIZE(particles%els(:))
      ALLOCATE (r(3,nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "particle_positions",3*nunits_tot*dp_size)

      DO iparticle=1,nunits_tot
         r(1:3,iparticle)=particles%els(iparticle)%r(1:3)
      ENDDO

      IF(ionode) &
         CALL mc_make_dat_file(r(:,:),nunits_tot,&
         abc(1:3)*angstrom,&
         'bias_temp.inp',bias_env,TRIM(fft_lib))

      CALL force_env_release(bias_env,error=error)
      globenv%input_file_name='bias_temp.inp'
      CALL mc_create_force_env(bias_env, globenv, box_change, &
         use_fist_force,error=error)

      DEALLOCATE(r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "r")

END SUBROUTINE mc_create_bias_force_env

END MODULE mc_control

!******************************************************************************
