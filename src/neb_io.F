!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief I/O Module for Nudged Elastic Band Calculation
!> \note
!>      Numerical accuracy for parallel runs:
!>       Each replica starts the SCF run from the one optimized
!>       in a previous run. It may happen then energies and derivatives
!>       of a serial run and a parallel run could be slightly different
!>       'cause of a different starting density matrix.
!>       Exact results are obtained using:
!>          EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!> \author Teodoro Laino 10.2006
! *****************************************************************************
MODULE neb_io
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get
  USE input_constants,                 ONLY: band_md_opt,&
                                             dump_xmol
  USE input_cp2k_motion,               ONLY: create_band_section
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE neb_md_utils,                    ONLY: get_temperatures
  USE neb_types,                       ONLY: neb_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_io' 

  PUBLIC :: read_neb_section,&
            dump_neb_info,&
            dump_replica_coordinates

CONTAINS

! *****************************************************************************
!> \brief Read data from the NEB input section
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE read_neb_section(neb_env, neb_section,error)
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_neb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(neb_env),cp_failure_level,routineP,error,failure)
    neb_env%istep = 0
    CALL section_vals_val_get(neb_section,"BAND_TYPE",i_val=neb_env%id_type,error=error)
    CALL section_vals_val_get(neb_section,"NUMBER_OF_REPLICA",i_val=neb_env%number_of_replica,error=error)
    CALL section_vals_val_get(neb_section,"K_SPRING",r_val=neb_env%K,error=error)
    CALL section_vals_val_get(neb_section,"ROTATE_FRAMES",l_val=neb_env%rotate_frames,error=error)
    CALL section_vals_val_get(neb_section,"OPTIMIZE_BAND%OPTIMIZE_END_POINTS",l_val=neb_env%optimize_end_points,&
         error=error)
    ! Climb Image NEB
    CALL section_vals_val_get(neb_section,"CI_NEB%NSTEPS_IT",i_val=neb_env%nsteps_it,error=error)
    ! Band Optimization Type
    CALL section_vals_val_get(neb_section,"OPTIMIZE_BAND%OPT_TYPE",i_val=neb_env%opt_type,error=error)
  END SUBROUTINE read_neb_section

! *****************************************************************************
!> \brief dump print info of a NEB run
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE dump_neb_info(neb_env, coords, vels, forces, particle_set, logger,&
       istep, energies, distances, output_unit, error)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: vels, forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(IN)                      :: istep
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: energies, distances
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dump_neb_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: mytype
    CHARACTER(LEN=default_string_length)     :: line, title
    INTEGER                                  :: crd, ener, frc, handle, i, &
                                                irep, natom, ndig, ndigl, &
                                                stat, ttst, vel
    LOGICAL                                  :: explicit, failure, lval
    REAL(KIND=dp)                            :: f_ann
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ekin, temperatures
    TYPE(cell_type), POINTER                 :: cell
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(unit)
    CALL get_cp2k_units ( unit )
    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       natom = SIZE(coords,1)
       ndig  = CEILING(LOG10(REAL(SIZE(coords,2)+1,KIND=dp)))
       CALL force_env_get(neb_env%force_env, cell=cell, error=error)
       DO irep = 1, SIZE(coords,2)
          ndigl = CEILING(LOG10(REAL(irep+1,KIND=dp)))
          WRITE(line,'(A,'//cp_to_string(ndig)//'("0"),T'//cp_to_string(11+ndig+1-ndigl)//',I0)')"Replica_nr_",irep
          crd = cp_print_key_unit_nr(logger,neb_env%motion_print_section,"TRAJECTORY",&
               extension=".xyz",file_form="FORMATTED", middle_name="pos-"//TRIM(line), error=error)
          IF (PRESENT(vels)) THEN
             vel = cp_print_key_unit_nr(logger,neb_env%motion_print_section,"VELOCITIES",&
                  extension=".xyz",file_form="FORMATTED", middle_name="vel-"//TRIM(line), error=error)
          END IF
          IF (PRESENT(forces)) THEN
             frc = cp_print_key_unit_nr(logger,neb_env%motion_print_section,"FORCES",&
                  extension=".xyz",file_form="FORMATTED", middle_name="force-"//TRIM(line), error=error)
          END IF
          ! Dump Trajectory
          IF (crd>0) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
             CALL write_particle_coordinates(particle_set,crd,&
                  dump_xmol,"POS",title,cell=cell,array=coords(:,irep))
             CALL m_flush(crd)             
          END IF
          ! Dump Velocities
          IF (vel>0.AND.PRESENT(vels)) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
             CALL write_particle_coordinates(particle_set,vel,&
                  dump_xmol,"VEL",title,cell=cell,array=vels(:,irep))
             CALL m_flush(vel)
          END IF
          ! Dump Forces
          IF (frc>0.AND.PRESENT(forces)) THEN
             !         *** This information can be digested by Molden ***
             WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
             CALL write_particle_coordinates(particle_set,frc,&
                  dump_xmol,"FRC",title,cell=cell,array=forces(:,irep))
             CALL m_flush(frc)
          END IF
          CALL cp_print_key_finished_output(crd,logger,neb_env%motion_print_section,&
               "TRAJECTORY", error=error)
          IF (PRESENT(vels)) THEN
             CALL cp_print_key_finished_output(vel,logger,neb_env%motion_print_section,&
                  "VELOCITIES", error=error)
          END IF
          IF (PRESENT(forces)) THEN
             CALL cp_print_key_finished_output(frc,logger,neb_env%motion_print_section,&
                  "FORCES", error=error)
          END IF
       END DO
       ! NEB summary info on screen
       IF (output_unit>0) THEN
          tc_section => section_vals_get_subs_vals(neb_env%neb_section,"OPTIMIZE_BAND%MD%TEMP_CONTROL",error=error)
          vc_section => section_vals_get_subs_vals(neb_env%neb_section,"OPTIMIZE_BAND%MD%VEL_CONTROL",error=error)
          ALLOCATE(temperatures(SIZE(coords,2)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(ekin(SIZE(coords,2)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_temperatures(vels,particle_set,temperatures,ekin=ekin,error=error)
          WRITE(output_unit,'(/)',ADVANCE="NO")
          WRITE(output_unit,FMT='(A,A)')' **************************************', &
               '*****************************************'
          NULLIFY(section,keyword,enum)
          CALL create_band_section(section,error=error)
          keyword => section_get_keyword(section,"BAND_TYPE",error=error)
          CALL keyword_get(keyword,enum=enum,error=error)
          mytype=TRIM(enum_i2c(enum,neb_env%id_type,error=error))
          WRITE(output_unit,FMT='(A,T61,A)')&
               ' BAND TYPE                     =',ADJUSTR(mytype)
          CALL section_release(section,error=error)
          WRITE(output_unit,FMT='(A,T61,A)')&
               ' BAND TYPE OPTIMIZATION        =',ADJUSTR ( neb_env%opt_type_label(1:20))
          WRITE ( output_unit, '( A,T71,I10 )' )&
               ' STEP NUMBER                   =', istep
          IF (neb_env%rotate_frames) WRITE ( output_unit, '( A,T71,L10 )' )&
               ' RMSD DISTANCE DEFINITION      =',neb_env%rotate_frames
          ! velocity control parameters output
          CALL section_vals_get(vc_section,explicit=explicit,error=error)
          IF (explicit) THEN          
             CALL section_vals_val_get(vc_section,"PROJ_VELOCITY_VERLET",l_val=lval, error=error)
             IF (lval) WRITE ( output_unit, '( A,T71,L10 )' )&
               ' PROJECTED VELOCITY VERLET     =',lval
             CALL section_vals_val_get(vc_section,"SD_LIKE", l_val=lval, error=error)
             IF (lval) WRITE ( output_unit, '( A,T71,L10)' )&
               ' STEEPEST DESCENT LIKE         =',lval
             CALL section_vals_val_get(vc_section,"ANNEALING", r_val=f_ann, error=error)
             IF (f_ann/=1.0_dp) THEN
                WRITE ( output_unit, '( A,T71,F10.5)' )&
               ' ANNEALING FACTOR              = ',f_ann
             END IF
          END IF
          ! temperature control parameters output
          CALL section_vals_get(tc_section,explicit=explicit,error=error)
          IF (explicit) THEN          
             CALL section_vals_val_get(tc_section,"TEMP_TOL_STEPS",i_val=ttst, error=error)
             IF (istep<=ttst) THEN
                CALL section_vals_val_get(tc_section,"TEMPERATURE",r_val=f_ann, error=error)
                WRITE ( output_unit, '( A,T71,F10.5)' )&
               ' TEMPERATURE TARGET            =',f_ann/unit%kelvin
             END IF
          END IF
          WRITE ( output_unit, '( A,T71,I10 )' )&
               ' NUMBER OF NEB REPLICA         =', &
               SIZE(coords,2)
          WRITE ( output_unit, '( A,T17,4F16.6)' )&
               ' DISTANCES REP =',distances(1:MIN(4,SIZE(distances)))
          IF (SIZE(distances)>4) THEN
             WRITE ( output_unit, '( T17,4F16.6)' )distances(5:SIZE(distances))
          END IF
          WRITE ( output_unit, '( A,T17,4F16.6)' )&
               ' ENERGIES [au] =',energies(1:MIN(4,SIZE(energies)))
          IF (SIZE(energies)>4) THEN
             WRITE ( output_unit, '( T17,4F16.6)' )energies(5:SIZE(energies))
          END IF
          IF (neb_env%opt_type==band_md_opt) THEN
             WRITE ( output_unit, '( A,T33,4(1X,F11.5))' )&
               ' REPLICA TEMPERATURES (K)      =',temperatures(1:MIN(4,SIZE(temperatures)))
             DO i = 5, SIZE(temperatures), 4
                WRITE ( output_unit, '( T33,4(1X,F11.5))' )&
                     temperatures(i:MIN(i+3,SIZE(temperatures)))
             END DO
          END IF
          WRITE ( output_unit, '( A,T56,F25.14)' )&
               ' BAND TOTAL ENERGY [au]        =',SUM(energies(:)+ekin(:))+&
                                                  neb_env%spring_energy
          WRITE(output_unit,FMT='(A,A)')' **************************************', &
               '*****************************************'
          DEALLOCATE(ekin,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(temperatures,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       ! Ener file
       ener = cp_print_key_unit_nr(logger,neb_env%neb_section,"ENERGY",&
            extension=".ener",file_form="FORMATTED", error=error)
       IF (ener>0) THEN
          WRITE(line,'(I0)')2*SIZE(coords,2)-1
          WRITE(ener,'(I10,'//TRIM(line)//'(1X,F20.9))')istep,&
               energies, distances
       END IF
       CALL cp_print_key_finished_output(ener,logger,neb_env%neb_section,&
                  "ENERGY", error=error)
       ! Dump Restarts
       CALL write_restart(force_env=neb_env%force_env,&
                          root_section=neb_env%root_section,&
                          force_env_section=neb_env%force_env_section,&
                          coords=coords,&
                          vels=vels,&
                          error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE dump_neb_info

! *****************************************************************************
!> \brief dump coordinates of a replica NEB
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE dump_replica_coordinates(particle_set, coords, i, ienum, iw)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    INTEGER, INTENT(IN)                      :: i, ienum, iw

    INTEGER                                  :: j

    IF (iw>0) THEN
       WRITE(iw,'(/,T2,"NEB|",75("*"))')
       WRITE(iw,'(T2,"NEB|",1X,A,I6,A)')"Geometry for Replica Nr. ",ienum,"."
       DO j = 1, SIZE(particle_set)
          WRITE(iw,'(T2,"NEB|",1X,A10,5X,3F15.9)')TRIM(particle_set(j)%atomic_kind%name),&
               coords((j-1)*3+1:(j-1)*3+3,i)    
       END DO
       WRITE(iw,'(T2,"NEB|",75("*"))')
    END IF
  END SUBROUTINE dump_replica_coordinates

END MODULE neb_io
