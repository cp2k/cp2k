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

MODULE mc_control
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_types,                        ONLY: mc_simulation_parameters_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule
  USE qs_parser,                       ONLY: open_file,&
                                             close_file
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_control"

  PUBLIC :: write_mc_restart,read_mc_restart


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

    INTEGER                                  :: handle, ipart, unit
    REAL(KIND=dp),DIMENSION(1:3)             :: abc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! open the file and write some simulation parameters
  CALL open_file(file_name=mc_par%restart_file,unit_number=unit,&
                 file_action='WRITE',file_status='REPLACE')
  WRITE(unit,*) nnstep
  WRITE(unit,*) mc_par%temp,mc_par%nunits_tot
  WRITE(unit,*) mc_par%ensemble

! get the cell length and coordinates
  CALL force_env_get(force_env,cell=cell,subsys=subsys,&
                         error=error)
  CALL get_cell(cell,abc=abc)
  CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

  WRITE(unit,'(F10.6)') abc(1)*angstrom  ! in angstroms
  WRITE(unit,*)
  DO ipart=1,mc_par%nunits_tot,mc_par%nunits
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'O ',particles%els(ipart)%r(1:3)*angstrom,' H2O'
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'H ',particles%els(ipart+1)%r(1:3)*angstrom,' H2O'
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'H ',particles%els(ipart+2)%r(1:3)*angstrom,' H2O'
  ENDDO
  
  CALL close_file(unit_number=unit)

! end the timing
  CALL timestop(0.0d0,handle)


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

    CHARACTER(LEN=20)                        :: chardum, ensemble
    INTEGER                                  :: handle, ipart, istat, &
                                                nunits_tot, unit
    REAL(KIND=dp)                            :: box_length, temp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp),DIMENSION(1:3)             :: abc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! open the file and read some simulation parameters
      CALL open_file(file_name=mc_par%restart_file,unit_number=unit,&
         file_action='READ',file_status='OLD')

      READ(unit,*) mc_par%nstart
      READ(unit,*) temp,nunits_tot
      READ(unit,*) ensemble

! do some checking
      IF(ABS(temp - mc_par%temp) .GT. 0.01d0) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The temperature in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file temp =',mc_par%temp
            WRITE(iw,*) 'Restart file temp =',temp
         ENDIF
         STOP
      ENDIF
      IF(nunits_tot .NE. mc_par%nunits_tot) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The total number of units in the restart ',&
            'file is not the same as the input file.'
            WRITE(iw,*) 'Input file units =',mc_par%nunits_tot
            WRITE(iw,*) 'Restart file units =',nunits_tot
         ENDIF
         mc_par%nunits_tot=nunits_tot
         mc_par%nchain=mc_par%nunits_tot/mc_par%nunits
!     STOP
      ENDIF
      IF(ensemble .NE. mc_par%ensemble) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The ensemble in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file ensemble =',mc_par%ensemble
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

      READ(unit,*)  box_length  ! in angstroms
      IF(ABS(box_length - abc(1)*angstrom) .GT. 0.0001d0 ) THEN
         IF (mc_par%ionode) THEN
            WRITE(iw,*) 'The cell length in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file cell length =',abc(1)*angstrom
            WRITE(iw,*) 'Restart file cell length =',box_length
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
         r(1:3,1) = (/ 0.576135d0, 0.969636d0,1.0678366d0 /)
         r(1:3,2) = (/ 0.178637d0, 1.613209d0, 0.426837d0 /)
         r(1:3,3) = (/ 0.466322d0, 1.3840462d0, 1.96186d0 /)

         CALL close_file(unit_number=unit)

         IF(mc_par%ionode) &
            CALL mc_make_dat_file(r(:,:),3,box_length,mc_par%dat_file)

      ELSE
         ALLOCATE (r(3,nunits_tot),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
            "r",3*nunits_tot*dp_size)
         DO ipart=1,mc_par%nunits_tot
            READ(unit,*) chardum,r(1:3,ipart),chardum
            r(1:3,ipart)=r(1:3,ipart)/angstrom
         ENDDO

         CALL close_file(unit_number=unit)

         IF(mc_par%ionode) &
            CALL mc_make_dat_file(r(:,:),mc_par%nunits_tot,&
            box_length,mc_par%dat_file)

      ENDIF

      
! end the timing
      CALL timestop(0.0d0,handle)

! deallcoate
      DEALLOCATE(r)

END SUBROUTINE read_mc_restart


END MODULE mc_control

!******************************************************************************
