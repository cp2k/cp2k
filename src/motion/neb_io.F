!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp2k_info,                       ONLY: get_runtime_info
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             use_mixed_force
  USE header,                          ONLY: cp2k_footer
  USE input_constants,                 ONLY: band_md_opt,&
                                             do_sm,&
                                             dump_xmol,&
                                             pot_neb_fe,&
                                             pot_neb_full,&
                                             pot_neb_me
  USE input_cp2k_neb,                  ONLY: create_band_section
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
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE neb_md_utils,                    ONLY: get_temperatures
  USE neb_types,                       ONLY: neb_type,&
                                             neb_var_type
  USE particle_methods,                ONLY: write_particle_coordinates
  USE particle_types,                  ONLY: get_particle_pos_or_vel,&
                                             particle_type
  USE physcon,                         ONLY: angstrom
  USE replica_types,                   ONLY: replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_io'

  PUBLIC :: read_neb_section,&
            dump_neb_info,&
            dump_replica_coordinates,&
            handle_band_file_names,&
            neb_rep_env_map_info

CONTAINS

! *****************************************************************************
!> \brief Read data from the NEB input section
!> \param neb_env ...
!> \param neb_section ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE read_neb_section(neb_env, neb_section,error)
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_neb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: wrk_section

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(neb_env),cp_failure_level,routineP,error,failure)
    neb_env%istep = 0
    CALL section_vals_val_get(neb_section,"BAND_TYPE",i_val=neb_env%id_type,error=error)
    CALL section_vals_val_get(neb_section,"NUMBER_OF_REPLICA",i_val=neb_env%number_of_replica,error=error)
    CALL section_vals_val_get(neb_section,"K_SPRING",r_val=neb_env%K,error=error)
    CALL section_vals_val_get(neb_section,"ROTATE_FRAMES",l_val=neb_env%rotate_frames,error=error)
    CALL section_vals_val_get(neb_section,"ALIGN_FRAMES",l_val=neb_env%align_frames,error=error)
    CALL section_vals_val_get(neb_section,"OPTIMIZE_BAND%OPTIMIZE_END_POINTS",l_val=neb_env%optimize_end_points,&
         error=error)
    ! Climb Image NEB
    CALL section_vals_val_get(neb_section,"CI_NEB%NSTEPS_IT",i_val=neb_env%nsteps_it,error=error)
    ! Band Optimization Type
    CALL section_vals_val_get(neb_section,"OPTIMIZE_BAND%OPT_TYPE",i_val=neb_env%opt_type,error=error)
    ! Use colvars
    CALL section_vals_val_get(neb_section,"USE_COLVARS",l_val=neb_env%use_colvar,error=error)
    CALL section_vals_val_get(neb_section,"POT_TYPE",i_val=neb_env%pot_type,error=error)
    ! Before continuing let's do some consistency check between keywords
    IF (neb_env%pot_type/=pot_neb_full) THEN
       ! Requires the use of colvars
       CALL cp_assert(neb_env%use_colvar,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "A potential energy function based on free energy or minimum energy"//&
            " was requested without enabling the usage of COLVARS. Both methods"//&
            " are based on COLVARS definition."//&
 CPSourceFileRef,&
            error=error,failure=failure)
       ! Moreover let's check if the proper sections have been defined..
       SELECT CASE(neb_env%pot_type)
       CASE(pot_neb_fe)
          wrk_section => section_vals_get_subs_vals(neb_env%root_section,"MOTION%MD",error=error)
          CALL section_vals_get(wrk_section,explicit=explicit, error=error)
          CALL cp_assert(explicit,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "A free energy BAND (colvars projected) calculation is requested"//&
               " but NONE MD section was defined in the input."//&
 CPSourceFileRef,&
                error=error,failure=failure)
       CASE(pot_neb_me)
          wrk_section => section_vals_get_subs_vals(neb_env%root_section,"MOTION%GEO_OPT",error=error)
          CALL section_vals_get(wrk_section,explicit=explicit, error=error)
          CALL cp_assert(explicit,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "A minimum energy BAND (colvars projected) calculation is requested"//&
               " but NONE GEO_OPT section was defined in the input."//&
 CPSourceFileRef,&
                error=error,failure=failure)
       END SELECT
    ELSE
       CALL cp_assert(.NOT.neb_env%use_colvar,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "A band calculation was requested with a full potential energy. USE_COLVAR cannot"//&
            " be set for this kind of calculation!"//&
 CPSourceFileRef,&
            error=error,failure=failure)
    END IF
    ! String Method
    CALL section_vals_val_get(neb_section,"STRING_METHOD%SMOOTHING",r_val=neb_env%smoothing,error=error)
    CALL section_vals_val_get(neb_section,"STRING_METHOD%SPLINE_ORDER",i_val=neb_env%spline_order,error=error)
    neb_env%reparametrize_frames =.FALSE.
    IF (neb_env%id_type==do_sm) THEN
       neb_env%reparametrize_frames =.TRUE.
    END IF
  END SUBROUTINE read_neb_section

! *****************************************************************************
!> \brief dump print info of a NEB run
!> \param neb_env ...
!> \param coords ...
!> \param vels ...
!> \param forces ...
!> \param particle_set ...
!> \param logger ...
!> \param istep ...
!> \param energies ...
!> \param distances ...
!> \param output_unit ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE dump_neb_info(neb_env, coords, vels, forces, particle_set, logger,&
       istep, energies, distances, output_unit, error)
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(neb_var_type), POINTER              :: coords
    TYPE(neb_var_type), OPTIONAL, POINTER    :: vels, forces
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
    CHARACTER(LEN=default_string_length)     :: line, title, unit_str
    INTEGER                                  :: crd, ener, frc, handle, i, &
                                                irep, ndig, ndigl, stat, &
                                                ttst, vel
    LOGICAL                                  :: explicit, failure, lval
    REAL(KIND=dp)                            :: f_ann, tmp_r1, unit_conv
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ekin, temperatures
    TYPE(cell_type), POINTER                 :: cell
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    ndig  = CEILING(LOG10(REAL(neb_env%number_of_replica+1,KIND=dp)))
    CALL force_env_get(neb_env%force_env, cell=cell, error=error)
    DO irep = 1, neb_env%number_of_replica
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
          ! Gather units of measure for output
          CALL section_vals_val_get(neb_env%motion_print_section,"TRAJECTORY%UNIT",&
               c_val=unit_str,error=error)
          unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          ! This information can be digested by Molden
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
          CALL write_particle_coordinates(particle_set,crd,dump_xmol,"POS",title,&
               cell=cell,array=coords%xyz(:,irep),unit_conv=unit_conv,error=error)
          CALL m_flush(crd)
       END IF
       ! Dump Velocities
       IF (vel>0.AND.PRESENT(vels)) THEN
          ! Gather units of measure for output
          CALL section_vals_val_get(neb_env%motion_print_section,"VELOCITIES%UNIT",&
               c_val=unit_str,error=error)
          unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
          CALL write_particle_coordinates(particle_set,vel,dump_xmol,"VEL",title,&
               cell=cell,array=vels%xyz(:,irep),unit_conv=unit_conv,error=error)
          CALL m_flush(vel)
       END IF
       ! Dump Forces
       IF (frc>0.AND.PRESENT(forces)) THEN
          ! Gather units of measure for output
          CALL section_vals_val_get(neb_env%motion_print_section,"FORCES%UNIT",&
               c_val=unit_str,error=error)
          unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          WRITE (UNIT=title,FMT="(A,I8,A,F20.10)")" i =",istep,", E =",energies(irep)
          CALL write_particle_coordinates(particle_set,frc,dump_xmol,"FRC",title,&
               cell=cell,array=forces%xyz(:,irep),unit_conv=unit_conv,error=error)
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
       ALLOCATE(temperatures(neb_env%number_of_replica),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ekin(neb_env%number_of_replica),stat=stat)
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
             tmp_r1 = cp_unit_from_cp2k(f_ann,"K",error=error)
             WRITE ( output_unit, '( A,T71,F10.5)' )&
            ' TEMPERATURE TARGET            =',tmp_r1
          END IF
       END IF
       WRITE ( output_unit, '( A,T71,I10 )' )&
            ' NUMBER OF NEB REPLICA         =', neb_env%number_of_replica
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
       WRITE(line,'(I0)')2*neb_env%number_of_replica-1
       WRITE(ener,'(I10,'//TRIM(line)//'(1X,F20.9))')istep,&
            energies, distances
    END IF
    CALL cp_print_key_finished_output(ener,logger,neb_env%neb_section,&
               "ENERGY", error=error)

    ! Dump Restarts
    CALL cp_add_default_logger(logger)
    CALL write_restart(force_env=neb_env%force_env,&
                       root_section=neb_env%root_section,&
                       coords=coords,&
                       vels=vels,&
                       error=error)
    CALL cp_rm_default_logger()

    CALL timestop(handle)

  END SUBROUTINE dump_neb_info

! *****************************************************************************
!> \brief dump coordinates of a replica NEB
!> \param particle_set ...
!> \param coords ...
!> \param i_rep ...
!> \param ienum ...
!> \param iw ...
!> \param use_colvar ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE dump_replica_coordinates(particle_set,coords,i_rep,ienum,iw,use_colvar)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: i_rep, ienum, iw
    LOGICAL, INTENT(IN)                      :: use_colvar

    INTEGER                                  :: iatom, j
    REAL(KIND=dp), DIMENSION(3)              :: r

    IF (iw > 0) THEN
       WRITE (iw,'(/,T2,"NEB|",75("*"))')
       WRITE (iw,'(T2,"NEB|",1X,A,I0,A)')&
          "Geometry for Replica Nr. ",ienum," in Angstrom"
       DO iatom=1,SIZE(particle_set)
          r(1:3) = get_particle_pos_or_vel(iatom,particle_set,coords%xyz(:,i_rep))
          WRITE (iw,'(T2,"NEB|",1X,A10,5X,3F15.9)')&
             TRIM(particle_set(iatom)%atomic_kind%name),r(1:3)*angstrom
       END DO
       IF (use_colvar) THEN
          WRITE (iw,'(/,T2,"NEB|",1X,A10)')"COLLECTIVE VARIABLES:"
          WRITE (iw,'(T2,"NEB|",16X,3F15.9)')&
             (coords%int(j,i_rep),j=1,SIZE(coords%int(:,:),1))
       END IF
       WRITE (iw,'(T2,"NEB|",75("*"))')
       CALL m_flush(iw)
    END IF

  END SUBROUTINE dump_replica_coordinates

! *****************************************************************************
!> \brief Handles the correct file names during a band calculation
!> \param rep_env ...
!> \param irep ...
!> \param n_rep ...
!> \param istep ...
!> \param error ...
!> \author Teodoro Laino  06.2009
! *****************************************************************************
  SUBROUTINE handle_band_file_names(rep_env, irep, n_rep, istep, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    INTEGER, INTENT(IN)                      :: irep, n_rep, istep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'handle_band_file_names', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: output_file_path, &
                                                replica_proj_name
    INTEGER                                  :: handle, handle2, i, ierr, j, &
                                                lp, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger, sub_logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(section_vals_type), POINTER         :: root_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CALL f_env_add_defaults(f_env_id=rep_env%f_env_id, f_env=f_env,&
         new_error=new_error, failure=failure, handle=handle2)
    logger => cp_error_get_logger(new_error)
    CALL force_env_get(f_env%force_env,root_section=root_section,error=new_error)
    j    = irep + (rep_env%local_rep_indices(1)-1)
    ! Get replica_project_name
    replica_proj_name = get_replica_project_name(rep_env, n_rep, j)
    lp=LEN_TRIM(replica_proj_name)
    CALL section_vals_val_set(root_section,"GLOBAL%PROJECT_NAME",&
         c_val=TRIM(replica_proj_name),error=error)
    logger%iter_info%project_name=replica_proj_name

    ! We change the file on which is pointing the global logger and error
    output_file_path=replica_proj_name(1:lp)//".out"
    CALL section_vals_val_set(root_section,"GLOBAL%OUTPUT_FILE_NAME",&
         c_val=TRIM(output_file_path),error=error)
    IF (logger%default_global_unit_nr>0) THEN
       CALL close_file(logger%default_global_unit_nr)
       CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
            file_action="WRITE", file_position="APPEND",&
            unit_number=logger%default_global_unit_nr,&
            skip_get_unit_number=.TRUE.)
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(/,(T2,A79))")&
            "*******************************************************************************",&
            "**                 BAND EVALUATION OF ENERGIES AND FORCES                    **",&
            "*******************************************************************************"
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(T2,A,T79,A)")"**","**"
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(T2,A,T79,A)")"**","**"
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(T2,A,I5,T41,A,I5,T79,A)")&
            "** Replica Env Nr. :",rep_env%local_rep_indices(1)-1,"Replica Band Nr. :",j,"**"
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(T2,A,I5,T79,A)")&
            "** Band  Step  Nr. :",istep,"**"
       WRITE (UNIT=logger%default_global_unit_nr,FMT="(T2,A79)")&
            "*******************************************************************************"
    END IF

    ! Handle specific case for mixed_env
    SELECT CASE (f_env%force_env%in_use)
    CASE (use_mixed_force)
       DO i = 1, f_env%force_env%mixed_env%ngroups
          IF (MODULO(i-1,f_env%force_env%mixed_env%ngroups)==&
              f_env%force_env%mixed_env%group_distribution(f_env%force_env%mixed_env%para_env%mepos)) THEN
             sub_logger => f_env%force_env%mixed_env%sub_logger(i)%p
             sub_logger%iter_info%project_name=replica_proj_name(1:lp)//"-r-"//TRIM(ADJUSTL(cp_to_string(i)))

             unit_nr = sub_logger%default_global_unit_nr
             IF (unit_nr>0) THEN
                CALL close_file(unit_nr)

                output_file_path=replica_proj_name(1:lp)//"-r-"//TRIM(ADJUSTL(cp_to_string(i)))//".out"
                CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
                     file_action="WRITE", file_position="APPEND",&
                     unit_number=unit_nr,skip_get_unit_number=.TRUE.)
             END IF
          END IF
       END DO
    END SELECT

    CALL f_env_rm_defaults(f_env=f_env,error=new_error,ierr=ierr,handle=handle2)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE handle_band_file_names

! *****************************************************************************
!> \brief  Constructs project names for BAND replicas
!> \param rep_env ...
!> \param n_rep ...
!> \param j ...
!> \retval replica_proj_name ...
!> \author Teodoro Laino  06.2009
! *****************************************************************************
  FUNCTION get_replica_project_name(rep_env, n_rep, j) RESULT(replica_proj_name)
    TYPE(replica_env_type), POINTER          :: rep_env
    INTEGER, INTENT(IN)                      :: n_rep, j
    CHARACTER(LEN=default_path_length)       :: replica_proj_name

    CHARACTER(len=*), PARAMETER :: routineN = 'get_replica_project_name', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: padding
    INTEGER                                  :: i, lp, ndigits
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Setup new replica project name and output file
    replica_proj_name = rep_env%original_project_name
    ! Find padding
    ndigits = CEILING(LOG10(REAL(n_rep+1,KIND=dp)))-&
         CEILING(LOG10(REAL(j+1,KIND=dp)))
    padding = ""
    DO i = 1, ndigits
       padding(i:i) = "0"
    END DO
    lp=LEN_TRIM(replica_proj_name)
    replica_proj_name(lp+1:LEN(replica_proj_name))="-BAND"//&
         TRIM(padding)//ADJUSTL(cp_to_string(j))
  END FUNCTION get_replica_project_name

! *****************************************************************************
!> \brief  Print some mapping infos in the replica_env setup output files
!>         i.e. prints in which files one can find information for each band
!>         replica
!> \param rep_env ...
!> \param neb_env ...
!> \param error ...
!> \author Teodoro Laino  06.2009
! *****************************************************************************
  SUBROUTINE neb_rep_env_map_info(rep_env, neb_env, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_rep_env_map_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: replica_proj_name
    INTEGER                                  :: handle2, ierr, irep, n_rep, &
                                                n_rep_neb, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    n_rep_neb = neb_env%number_of_replica
    n_rep     = rep_env%nrep
    CALL f_env_add_defaults(f_env_id=rep_env%f_env_id, f_env=f_env,&
         new_error=new_error, failure=failure, handle=handle2)
    logger => cp_error_get_logger(new_error)
    output_unit = logger%default_global_unit_nr
    IF (output_unit > 0) THEN
       WRITE (UNIT=output_unit,FMT='(/,(T2,A79))')&
            "*******************************************************************************",&
            "**                  MAPPING OF BAND REPLICA TO REPLICA ENV                   **",&
            "*******************************************************************************"
       WRITE (UNIT=output_unit,FMT='(T2,A,I6,T32,A,T79,A)')&
            "** Replica Env Nr.: ",rep_env%local_rep_indices(1)-1,&
            "working on the following BAND replicas","**"
       WRITE (UNIT=output_unit,FMT='(T2,A79)')&
            "**                                                                           **"
    END IF
    DO irep = 1, n_rep_neb, n_rep
       replica_proj_name = get_replica_project_name(rep_env, n_rep_neb, irep+rep_env%local_rep_indices(1)-1)
       IF (output_unit > 0) THEN
          WRITE (UNIT=output_unit,FMT='(T2,A,I6,T32,A,T79,A)')&
               "** Band Replica   Nr.: ",irep+rep_env%local_rep_indices(1)-1,&
               "Output available on file: "//TRIM(replica_proj_name)//".out","**"
       END IF
    END DO
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT='(T2,A79)')&
            "**                                                                           **",&
            "*******************************************************************************"
       WRITE (UNIT=output_unit,FMT='(/)')
    END IF
    ! update runtime info before printing the footer
    CALL get_runtime_info()
    ! print footer
    CALL cp2k_footer(output_unit)
    CALL f_env_rm_defaults(f_env=f_env,error=new_error,ierr=ierr,handle=handle2)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE neb_rep_env_map_info

END MODULE neb_io
