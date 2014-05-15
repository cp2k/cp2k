!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module with utility for  Nudged Elastic Band Calculation
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
MODULE neb_utils
  USE bibliography,                    ONLY: E2002,&
                                             Elber1987,&
                                             Jonsson1998,&
                                             Jonsson2000_1,&
                                             Jonsson2000_2,&
                                             Wales2004,&
                                             cite_reference
  USE colvar_utils,                    ONLY: eval_colvar,&
                                             get_clv_force,&
                                             set_colvars_target
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type,&
                                             get_energy,&
                                             get_force,&
                                             get_pos,&
                                             set_pos
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get
  USE geo_opt,                         ONLY: cp_geo_opt
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       band_diis_opt, do_b_neb, do_band_cartesian, do_ci_neb, do_d_neb, &
       do_eb, do_it_neb, do_sm, pot_neb_fe, pot_neb_full, pot_neb_me
  USE input_cp2k_check,                ONLY: remove_restart_info
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathlib,                         ONLY: matvec_3x3
  USE md_run,                          ONLY: qs_mol_dyn
  USE neb_io,                          ONLY: dump_replica_coordinates,&
                                             handle_band_file_names
  USE neb_md_utils,                    ONLY: neb_initialize_velocity
  USE neb_types,                       ONLY: neb_type,&
                                             neb_var_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: bohr
  USE replica_methods,                 ONLY: rep_env_calc_e_f
  USE replica_types,                   ONLY: rep_env_sync,&
                                             replica_env_type
  USE rmsd,                            ONLY: rmsd3
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_utils'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: build_replica_coords,&
            neb_calc_energy_forces,&
            get_tangent,&
            get_neb_force,&
            reorient_images,&
            reparametrize_images,&
            check_convergence

CONTAINS

! *****************************************************************************
!> \brief Computes the distance between two replica
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_replica_distance(particle_set,coords,i0,i,distance,iw,rotate,&
                                  error)
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: i0, i
    REAL(KIND=dp), INTENT(OUT)               :: distance
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: rotate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_replica_distance', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_rotate

    failure = .FALSE.
    my_rotate = .FALSE.
    IF (PRESENT(rotate)) my_rotate = rotate
    ! The rotation of the replica is enabled exclusively when working in
    ! cartesian coordinates
    IF (my_rotate.AND.(coords%in_use == do_band_cartesian)) THEN
       CPPostcondition(PRESENT(particle_set),cp_failure_level,routineP,error,failure)
       CALL rmsd3(particle_set,coords%xyz(:,i),coords%xyz(:,i0),&
                  iw,rotate=my_rotate,error=error)
    END IF
    distance = SQRT(DOT_PRODUCT(coords%wrk(:,i)-coords%wrk(:,i0),&
                                coords%wrk(:,i)-coords%wrk(:,i0)))

  END SUBROUTINE neb_replica_distance

! *****************************************************************************
!> \brief Constructs or Read the coordinates for all replica
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE build_replica_coords(neb_section, particle_set,&
       coords, vels, neb_env, iw, globenv, para_env, error)
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neb_var_type), POINTER              :: coords, vels
    TYPE(neb_type), POINTER                  :: neb_env
    INTEGER, INTENT(IN)                      :: iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_replica_coords', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: dummy_char
    INTEGER :: handle, i_rep, iatom, ic, input_nr_replica, is, ivar, j, &
      jtarg, n_rep, natom, neb_nr_replica, nr_replica_to_interpolate, nval, &
      nvar, shell_index, stat
    LOGICAL                                  :: check, explicit, failure, &
                                                my_end, skip_vel_section
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distance
    REAL(KIND=dp), DIMENSION(3)              :: r
    REAL(KIND=dp), DIMENSION(:), POINTER     :: initial_colvars, rptr
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(section_vals_type), POINTER         :: coord_section, &
                                                replica_section, vel_section

    NULLIFY (parser)
    CALL timeset(routineN,handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(vels),cp_failure_level,routineP,error,failure)
       neb_nr_replica = neb_env%number_of_replica
       replica_section => section_vals_get_subs_vals(neb_section,"REPLICA",error=error)
       CALL section_vals_get(replica_section,n_repetition=input_nr_replica,error=error)
       ! Calculation is aborted if input replicas are more then the requested ones for the BAND..
       CPPostcondition((input_nr_replica<=neb_nr_replica),cp_failure_level,routineP,error,failure)
       ! Read in replicas coordinates
       skip_vel_section = (input_nr_replica /= neb_nr_replica)
       IF ((iw > 0).AND.skip_vel_section) THEN
          WRITE(iw,'(T2,A)')'NEB| The number of replica in the input is different from the number',&
                            'NEB| of replica requested for NEB. More Replica will be interpolated.',&
                            'NEB| Therefore the possibly provided velocities will not be read.'
       END IF
       ! Further check on velocity section...
       DO i_rep=1,input_nr_replica
          vel_section => section_vals_get_subs_vals(replica_section,"VELOCITY",&
                                                    i_rep_section=i_rep,error=error)
          CALL section_vals_get(vel_section,explicit=explicit,error=error)
          skip_vel_section = skip_vel_section.OR.(.NOT.explicit)
       END DO
       ! Setup cartesian coordinates and COLVAR (if requested)
       coords%xyz(:,:) = 0.0_dp
       DO i_rep=1,input_nr_replica
          coord_section => section_vals_get_subs_vals(replica_section,"COORD",&
                                                      i_rep_section=i_rep,error=error)
          CALL section_vals_get(coord_section,explicit=explicit,error=error)
          ! Cartesian Coordinates
          IF (explicit) THEN
             CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                                       n_rep_val=natom,error=error)
             CPPostcondition((natom == SIZE(particle_set)),cp_failure_level,routineP,error,failure)
             DO iatom=1,natom
                CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                                          i_rep_val=iatom,r_vals=rptr,error=error)
                ic = 3*(iatom - 1)
                coords%xyz(ic+1:ic+3,i_rep) = rptr(1:3)*bohr
                ! Initially core and shell positions are set to the atomic positions
                shell_index = particle_set(iatom)%shell_index
                IF (shell_index /= 0) THEN
                   is = 3*(natom + shell_index - 1)
                   coords%xyz(is+1:is+3,i_rep) = coords%xyz(ic+1:ic+3,i_rep)
                END IF
             END DO
          ELSE
             CALL section_vals_val_get(replica_section,"COORD_FILE_NAME",&
                                       i_rep_section=i_rep,c_val=filename,error=error)
             CPPostcondition(TRIM(filename)/="",cp_failure_level,routineP,error,failure)
             CALL parser_create(parser,filename,para_env=para_env,parse_white_lines=.TRUE.,&
                                error=error)
             CALL parser_get_next_line(parser,1,error=error)
             ! Start parser
             CALL parser_get_object(parser,natom,error=error)
             CPPostcondition((natom == SIZE(particle_set)),cp_failure_level,routineP,error,failure)
             CALL parser_get_next_line(parser,1,error=error)
             DO iatom=1,natom
                ! Atom coordinates
                CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
                CALL cp_assert(.NOT.my_end,cp_fatal_level,cp_assertion_failed,routineP,&
                     "Number of lines in XYZ format not equal to the number of atoms."//&
                     " Error in XYZ format for REPLICA coordinates. Very probably the"//&
                     " line with title is missing or is empty. Please check the XYZ file and rerun your job!"//&
CPSourceFileRef,&
                only_ionode=.TRUE.)
                READ (parser%input_line,*) dummy_char,r(1:3)
                ic = 3*(iatom - 1)
                coords%xyz(ic+1:ic+3,i_rep) = r(1:3)*bohr
                ! Initially core and shell positions are set to the atomic positions
                shell_index = particle_set(iatom)%shell_index
                IF (shell_index /= 0) THEN
                   is = 3*(natom + shell_index - 1)
                   coords%xyz(is+1:is+3,i_rep) = coords%xyz(ic+1:ic+3,i_rep)
                END IF
             END DO
             CALL parser_release(parser,error=error)
          END IF
          ! Collective Variables
          IF (neb_env%use_colvar) THEN
             CALL  section_vals_val_get(replica_section,"COLLECTIVE",&
                                        i_rep_section=i_rep,n_rep_val=n_rep,error=error)
             IF (n_rep /= 0) THEN
                ! Read the values of the collective variables
                NULLIFY (initial_colvars)
                CALL  section_vals_val_get(replica_section,"COLLECTIVE",&
                                           i_rep_section=i_rep,r_vals=initial_colvars,error=error)
                check = (neb_env%nsize_int == SIZE(initial_colvars))
                CPPostcondition(check,cp_failure_level,routineP,error,failure)
                coords%int(:,i_rep) = initial_colvars
             ELSE
                ! Compute the values of the collective variables
                CALL eval_colvar(neb_env%force_env,coords%xyz(:,i_rep),coords%int(:,i_rep),error=error)
             END IF
          END IF
          ! Dump cartesian and colvar info..
          CALL dump_replica_coordinates(particle_set,coords,i_rep,i_rep,iw,neb_env%use_colvar)
          ! Setup Velocities
          IF (skip_vel_section) THEN
             CALL neb_initialize_velocity(vels%wrk,neb_section,particle_set,&
                                          i_rep,iw,globenv,neb_env,error)
          ELSE
             vel_section => section_vals_get_subs_vals(replica_section,"VELOCITY",&
                                                       i_rep_section=i_rep,error=error)
             CALL section_vals_val_get(vel_section,"_DEFAULT_KEYWORD_",&
                                       n_rep_val=nval,error=error)
             ! Setup Velocities for collective or cartesian coordinates
             IF (neb_env%use_colvar) THEN
                nvar = SIZE(vels%wrk,1)
                CPPostcondition((nval == nvar),cp_failure_level,routineP,error,failure)
                DO ivar=1,nvar
                   CALL section_vals_val_get(vel_section,"_DEFAULT_KEYWORD_",&
                                             i_rep_val=ivar,r_vals=rptr,error=error)
                   vels%wrk(ivar,i_rep) = rptr(1)
                END DO
             ELSE
                natom = SIZE(particle_set)
                CPPostcondition((nval == natom),cp_failure_level,routineP,error,failure)
                DO iatom=1,natom
                   CALL section_vals_val_get(vel_section,"_DEFAULT_KEYWORD_",&
                                             i_rep_val=iatom,r_vals=rptr,error=error)
                   ic = 3*(iatom - 1)
                   vels%wrk(ic+1:ic+3,i_rep) = rptr(1:3)
                   ! Initially set shell velocities to core velocity
                   shell_index = particle_set(iatom)%shell_index
                   IF (shell_index /= 0) THEN
                      is = 3*(natom + shell_index - 1)
                      vels%wrk(is+1:is+3,i_rep) = vels%wrk(ic+1:ic+3,i_rep)
                   END IF
                END DO
             END IF
          END IF
       END DO ! i_rep
       ALLOCATE (distance(neb_nr_replica-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (input_nr_replica < neb_nr_replica) THEN
          ! Interpolate missing replicas
          nr_replica_to_interpolate = neb_nr_replica - input_nr_replica
          distance = 0.0_dp
          IF (iw > 0) THEN
             WRITE(iw,'(T2,A,I0,A)')'NEB| Interpolating ',nr_replica_to_interpolate,' missing Replica.'
          END IF
          DO WHILE (nr_replica_to_interpolate > 0)
             ! Compute distances between known images to find the interval
             ! where to add a new image
             DO j=1,input_nr_replica-1
                CALL neb_replica_distance(particle_set,coords,j,j+1,distance(j),iw,&
                                          rotate=neb_env%align_frames,error=error)
             END DO
             jtarg = MAXLOC(distance(1:input_nr_replica),1)
             IF (iw > 0) THEN
                WRITE (iw,'(T2,3(A,I0),A)')'NEB| Interpolating Nr. ',&
                  nr_replica_to_interpolate,' missing Replica between Replica Nr. ',&
                  jtarg,' and ',jtarg+1,'.'
             END IF
             input_nr_replica = input_nr_replica + 1
             nr_replica_to_interpolate = nr_replica_to_interpolate - 1
             ! Interpolation is a simple bisection in XYZ
             coords%xyz(:,jtarg+2:input_nr_replica) = coords%xyz(:,jtarg+1:input_nr_replica-1)
             coords%xyz(:,jtarg+1) = (coords%xyz(:,jtarg) + coords%xyz(:,jtarg+2))/2.0_dp
             IF (neb_env%use_colvar) THEN
                ! Interpolation is a simple bisection also in internal coordinates
                ! in this case the XYZ coordinates need only as a starting point for computing
                ! the potential energy function. The reference are the internal coordinates
                ! interpolated here after..
                coords%int(:,jtarg+2:input_nr_replica) = coords%int(:,jtarg+1:input_nr_replica-1)
                coords%int(:,jtarg+1) = (coords%int(:,jtarg)+coords%int(:,jtarg+2))/2.0_dp
             END IF
             vels%wrk(:,jtarg+2:input_nr_replica) = vels%wrk(:,jtarg+1:input_nr_replica-1)
             vels%wrk(:,jtarg+1) = 0.0_dp
             CALL dump_replica_coordinates(particle_set, coords, jtarg+1,&
                                           input_nr_replica, iw, neb_env%use_colvar)
             CALL neb_initialize_velocity(vels%wrk, neb_section, particle_set,&
                                          jtarg+1, iw, globenv, neb_env, error)
          END DO
       END IF
       vels%wrk(:,1)              = 0.0_dp
       vels%wrk(:,neb_nr_replica) = 0.0_dp
       ! If we perform a DIIS optimization we don't need velocities
       IF (neb_env%opt_type == band_diis_opt) vels%wrk = 0.0_dp
       ! Compute distances between replicas and in case of Cartesian Coordinates
       ! Rotate the frames in order to minimize the RMSD
       DO j=1,input_nr_replica-1
          CALL neb_replica_distance(particle_set,coords,j,j+1,distance(j),iw,&
                                    rotate=neb_env%align_frames,error=error)
       END DO
       DEALLOCATE (distance,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_replica_coords

! *****************************************************************************
!> \brief Driver to compute energy and forces within a NEB,
!>      Based on the use of the replica_env
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
       particle_set, output_unit, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    TYPE(neb_var_type), POINTER              :: coords
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: energies
    TYPE(neb_var_type), POINTER              :: forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_calc_energy_forces', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), DIMENSION(3), &
      PARAMETER                              :: lab = (/"X","Y","Z"/)

    INTEGER                                  :: handle, i, irep, j, n_int, &
                                                n_rep, n_rep_neb, nsize_wrk, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tangent, tmp_a, tmp_b
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cvalues, Mmatrix, Mmatrix_tmp

    CALL timeset(routineN,handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       n_int     = neb_env%nsize_int
       n_rep_neb = neb_env%number_of_replica
       n_rep     = rep_env%nrep
       nsize_wrk = coords%size_wrk(1)
       energies  = 0.0_dp
       ALLOCATE(cvalues(n_int,n_rep),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Mmatrix_tmp(n_int*n_int,n_rep),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Mmatrix(n_int*n_int,n_rep_neb),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (output_unit>0) WRITE(output_unit,'(/,T2,A)')"NEB| Computing Energies and Forces"
       DO irep = 1, n_rep_neb, n_rep
          DO j = 0, n_rep-1
             IF (irep+j <= n_rep_neb) THEN
                ! If the number of replica in replica_env and the number of replica
                ! used in the NEB does not match, the other replica in replica_env
                ! just compute energies and forces keeping the fixed coordinates and
                ! forces
                rep_env%r(:,j+1) = coords%xyz(:,irep+j)
             END IF
          END DO
          ! Fix file name for BAND replicas.. Each BAND replica has its own file
          ! independently from the number of replicas in replica_env..
          CALL handle_band_file_names(rep_env, irep, n_rep_neb, neb_env%istep, error)
          ! Let's select the potential we want to use for the band calculation
          SELECT CASE(neb_env%pot_type)
          CASE(pot_neb_full)
             ! Full potential Energy
             CALL rep_env_calc_e_f(rep_env, calc_f=.TRUE., error=error)
          CASE(pot_neb_fe)
             ! Free Energy Case
             CALL perform_replica_md (rep_env, coords, irep, n_rep_neb, cvalues, Mmatrix_tmp, error)
          CASE(pot_neb_me)
             ! Minimum Potential Energy Case
             CALL perform_replica_geo(rep_env, coords, irep, n_rep_neb, cvalues, Mmatrix_tmp, error)
          END SELECT

          DO j = 0, n_rep-1
             IF (irep+j <= n_rep_neb) THEN
                ! Copy back Forces and Energies
                forces%wrk(:,irep+j) = rep_env%f(1:nsize_wrk,j+1)
                energies(irep+j)     = rep_env%f(rep_env%ndim+1,j+1)
                SELECT CASE (neb_env%pot_type)
                CASE(pot_neb_full)
                   ! Dump Info
                   IF (output_unit > 0) THEN
                      WRITE(output_unit,'(T2,A,I5,A,I5,A)')&
                           "NEB| REPLICA Nr.",irep+j,"- Energy and Forces"
                      WRITE(output_unit,'(T2,A,T42,A,9X,F15.9)')&
                           "NEB|"," Total Energy: ",rep_env%f(rep_env%ndim+1,j+1)
                      WRITE(output_unit,'(T2,"NEB|",T10,"ATOM",T33,3(9X,A,7X))')lab(1),lab(2),lab(3)
                      DO i = 1, SIZE(particle_set)
                         WRITE(output_unit,'(T2,"NEB|",T12,A,T30,3(2X,F15.9))')&
                              particle_set(i)%atomic_kind%name,&
                              rep_env%f((i-1)*3+1:(i-1)*3+3,j+1)
                      END DO
                   END IF
                CASE(pot_neb_fe, pot_neb_me)
                   ! Let's update the cartesian coordinates. This will make
                   ! easier the next evaluation of energies and forces...
                   coords%xyz(:,irep+j) = rep_env%r(1:rep_env%ndim,j+1)
                   Mmatrix(:,irep+j)    = Mmatrix_tmp(:,j+1)
                   IF (output_unit>0) THEN
                      WRITE(output_unit,'(/,T2,A,I5,A,I5,A)')&
                           "NEB| REPLICA Nr.",irep+j,"- Energy, Collective Variables,  Forces"
                      WRITE(output_unit,'(T2,A,T42,A,9X,F15.9)')&
                           "NEB|"," Total Energy: ",rep_env%f(rep_env%ndim+1,j+1)
                      WRITE(output_unit,&
                           '(T2,"NEB|",T10,"CV Nr.",12X,"Expected COLVAR",5X,"Present COLVAR",10X,"Forces")')
                      DO i = 1, n_int
                         WRITE(output_unit,'(T2,"NEB|",T12,I2,7X,3(5X,F15.9))')&
                              i,coords%int(i,irep+j),cvalues(i,j+1),rep_env%f(i,j+1)
                      END DO
                   END IF
                END SELECT
             END IF
          END DO
       END DO
       DEALLOCATE(cvalues,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Mmatrix_tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (PRESENT(neb_env)) THEN
          ! First identify the image of the chain with the higher potential energy
          ! First and last point of the band are never considered
          neb_env%nr_HE_image = MAXLOC(energies(2:n_rep_neb-1),1)+1
          ALLOCATE(tangent(nsize_wrk),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Then modify image forces accordingly to the scheme chosen for the
          ! calculation.
          neb_env%spring_energy = 0.0_dp
          IF (neb_env%optimize_end_points) THEN
             ALLOCATE(tmp_a(SIZE(forces%wrk,1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp_b(SIZE(forces%wrk,1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             tmp_a = forces%wrk(:,1)
             tmp_b = forces%wrk(:,SIZE(forces%wrk,2))
          END IF
          DO i = 2, neb_env%number_of_replica
             CALL get_tangent(neb_env,coords,i,tangent,energies,output_unit,error)
             CALL get_neb_force(neb_env,tangent,coords,i,forces,Mmatrix=Mmatrix,&
                  iw=output_unit,error=error)
          END DO
          IF (neb_env%optimize_end_points) THEN
             forces%wrk(:,1)                  = tmp_a ! Image A
             forces%wrk(:,SIZE(forces%wrk,2)) = tmp_b ! Image B
             DEALLOCATE(tmp_a,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(tmp_b,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
             ! Nullify forces on the two end points images
             forces%wrk(:,1)                  = 0.0_dp ! Image A
             forces%wrk(:,SIZE(forces%wrk,2)) = 0.0_dp ! Image B
          END IF
          DEALLOCATE(tangent,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(Mmatrix,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE neb_calc_energy_forces

! *****************************************************************************
!> \brief Driver to perform an MD run on each single replica to
!>      compute specifically Free Energies in a NEB scheme
!> \author Teodoro Laino  01.2007
! *****************************************************************************
  SUBROUTINE perform_replica_md(rep_env, coords, irep, n_rep_neb, cvalues, Mmatrix, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: irep, n_rep_neb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: cvalues, Mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perform_replica_md', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, ierr, j, n_el
    LOGICAL                                  :: explicit, failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: md_section, root_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure, handle=handle2)
       logger => cp_error_get_logger(new_error)
       CALL force_env_get(f_env%force_env,globenv=globenv,&
            root_section=root_section,error=new_error)
       j    = rep_env%local_rep_indices(1)-1
       n_el = 3*rep_env%nparticle
       Mmatrix = 0.0_dp
       ! Syncronize position on the replica procs
       CALL set_pos(rep_env%f_env_id,rep_env%r(:,j+1),n_el,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       !
       IF (irep+j <= n_rep_neb) THEN
          logger%iter_info%iteration(2)=irep+j
          CALL remove_restart_info(root_section,error=new_error)
          md_section => section_vals_get_subs_vals(root_section,"MOTION%MD",error=error)
          CALL section_vals_get(md_section,explicit=explicit,error=error)
          CPAssert(explicit,cp_failure_level,routineP,error,failure)
          ! Let's syncronize the target of Collective Variables for this run
          CALL set_colvars_target(coords%int(:,irep+j), f_env%force_env, error=error)
          ! Do a molecular dynamics and get back the derivative
          ! of the free energy w.r.t. the colvar and the metric tensor
          CALL qs_mol_dyn(f_env%force_env,globenv=globenv,error=new_error)
          ! Collect the equilibrated coordinates
          CALL get_pos(rep_env%f_env_id, rep_env%r(1:n_el,j+1), n_el, ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          ! Write he gradients in the colvar coordinates into the replica_env array
          ! and copy back also the metric tensor..
          ! work in progress..
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          rep_env%f(:,j+1) = 0.0_dp
          Mmatrix          = 0.0_dp
       ELSE
          rep_env%r(:,j+1) = 0.0_dp
          rep_env%f(:,j+1) = 0.0_dp
          cvalues(:,j+1)   = 0.0_dp
          Mmatrix(:,j+1)   = 0.0_dp
       END IF
       IF (.NOT.failure) THEN
          CALL rep_env_sync(rep_env,rep_env%f,error=error)
          CALL rep_env_sync(rep_env,rep_env%r,error=error)
          CALL rep_env_sync(rep_env,cvalues, error=error)
          CALL rep_env_sync(rep_env,Mmatrix, error=error)
       END IF
       CALL f_env_rm_defaults(f_env=f_env,error=new_error,ierr=ierr,handle=handle2)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE perform_replica_md

! *****************************************************************************
!> \brief Driver to perform a GEO_OPT run on each single replica to
!>        compute specifically minimum energies in a collective variable
!>        NEB scheme
!> \author Teodoro Laino 05.2007
! *****************************************************************************
  SUBROUTINE perform_replica_geo(rep_env, coords, irep, n_rep_neb, cvalues, Mmatrix, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: irep, n_rep_neb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: cvalues, Mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perform_replica_geo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, ierr, j, n_el
    LOGICAL                                  :: explicit, failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: geoopt_section, root_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure, handle=handle2)
       logger => cp_error_get_logger(new_error)
       CALL force_env_get(f_env%force_env,globenv=globenv,&
            root_section=root_section,error=new_error)
       j    = rep_env%local_rep_indices(1)-1
       n_el = 3*rep_env%nparticle
       Mmatrix = 0.0_dp
       ! Syncronize position on the replica procs
       CALL set_pos(rep_env%f_env_id,rep_env%r(:,j+1),n_el,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       IF (irep+j <= n_rep_neb) THEN
          logger%iter_info%iteration(2)=irep+j
          CALL remove_restart_info(root_section,error=new_error)
          geoopt_section => section_vals_get_subs_vals(root_section,"MOTION%GEO_OPT",error=error)
          CALL section_vals_get(geoopt_section,explicit=explicit,error=error)
          CPAssert(explicit,cp_failure_level,routineP,error,failure)
          ! Let's syncronize the target of Collective Variables for this run
          CALL set_colvars_target(coords%int(:,irep+j), f_env%force_env, error=error)
          ! Do a geometry optimization..
          CALL cp_geo_opt(f_env%force_env,globenv=globenv,error=new_error)
          ! Once the geometry optimization is ended let's do a single run
          ! without any constraints/restraints
          CALL force_env_calc_energy_force ( f_env%force_env, &
               calc_force=.TRUE., skip_external_control=.TRUE., error=error)
          ! Collect the optimized coordinates
          CALL get_pos(rep_env%f_env_id, rep_env%r(1:n_el,j+1), n_el, ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          ! Collect the gradients in cartesian coordinates
          CALL get_force(rep_env%f_env_id, rep_env%f(1:n_el,j+1), n_el, ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          ! Copy the energy
          CALL get_energy(rep_env%f_env_id, rep_env%f(n_el+1,j+1), ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          ! The gradients in the colvar coordinates
          CALL get_clv_force(f_env%force_env, rep_env%f(1:n_el,j+1), rep_env%r(1:n_el,j+1),&
               SIZE(coords%xyz,1), SIZE(coords%wrk,1), cvalues(:,j+1), Mmatrix(:,j+1),&
               error)
       ELSE
          rep_env%r(:,j+1) = 0.0_dp
          rep_env%f(:,j+1) = 0.0_dp
          cvalues(:,j+1)   = 0.0_dp
          Mmatrix(:,j+1)   = 0.0_dp
       END IF
       IF (.NOT.failure) THEN
          CALL rep_env_sync(rep_env,rep_env%f,error=error)
          CALL rep_env_sync(rep_env,rep_env%r,error=error)
          CALL rep_env_sync(rep_env,cvalues, error=error)
          CALL rep_env_sync(rep_env,Mmatrix, error=error)
       END IF
       CALL f_env_rm_defaults(f_env=f_env,error=new_error,ierr=ierr,handle=handle2)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE perform_replica_geo

! *****************************************************************************
!> \brief Computes the tangent for point i of the NEB chain
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE get_tangent(neb_env, coords, i, tangent, energies, iw, error)
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: tangent
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: energies
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_tangent', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KINd=dp)                            :: distance0, distance1, &
                                                distance2, DVmax, Dvmin

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
    tangent(:) = 0.0_dp
    ! For the last point we don't need any tangent..
    IF (i==neb_env%number_of_replica) RETURN
    ! Several kind of tangents implemented...
    SELECT CASE(neb_env%id_type)
    CASE(do_eb)
       tangent(:) = 0.0_dp
    CASE(do_b_neb)
       CALL neb_replica_distance(coords=coords,i0=i,i=i-1,distance=distance1,iw=iw,&
            rotate=.FALSE.,error=error)
       CALL neb_replica_distance(coords=coords,i0=i+1,i=i,distance=distance2,iw=iw,&
            rotate=.FALSE.,error=error)
       tangent(:) = (coords%wrk(:,i)  -coords%wrk(:,i-1))/distance1 +&
                    (coords%wrk(:,i+1)-coords%wrk(:,i)  )/distance2
    CASE(do_it_neb, do_ci_neb, do_d_neb)
       IF      ((energies(i+1).GT.energies(i)).AND.(energies(i).GT.(energies(i-1)))) THEN
          tangent(:) = coords%wrk(:,i+1) - coords%wrk(:,i)
       ELSE IF ((energies(i+1).LT.energies(i)).AND.(energies(i).LT.(energies(i-1)))) THEN
          tangent(:) = coords%wrk(:,i)   - coords%wrk(:,i-1)
       ELSE
          DVmax = MAX(ABS(energies(i+1)-energies(i)),ABS(energies(i-1)-energies(i)))
          DVmin = MIN(ABS(energies(i+1)-energies(i)),ABS(energies(i-1)-energies(i)))
          IF (energies(i+1) .GE. energies(i-1)) THEN
             tangent(:) = ( coords%wrk(:,i+1) - coords%wrk(:,i))*DVmax+(coords%wrk(:,i) - coords%wrk(:,i-1))*DVmin
          ELSE
             tangent(:) = ( coords%wrk(:,i+1) - coords%wrk(:,i))*DVmin+(coords%wrk(:,i) - coords%wrk(:,i-1))*DVmax
          END IF
       END IF
    CASE(do_sm)
       ! String method..
       tangent(:) = 0.0_dp
    END SELECT
    distance0  = SQRT(DOT_PRODUCT(tangent(:),tangent(:)))
    IF (distance0 /= 0.0_dp) tangent(:) = tangent(:) / distance0
  END SUBROUTINE get_tangent

! *****************************************************************************
!> \brief Computes the forces for point i of the NEB chain
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE get_neb_force(neb_env,tangent,coords,i,forces,tag,Mmatrix,&
       iw,error)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tangent
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: i
    TYPE(neb_var_type), POINTER              :: forces
    INTEGER, INTENT(IN), OPTIONAL            :: tag
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: Mmatrix
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neb_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, my_tag, nsize_wrk, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: distance1, distance2, fac, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dtmp1, wrk

    failure = .FALSE.
    my_tag  = neb_env%id_type
    IF (PRESENT(tag)) my_tag = tag
    CPPostcondition(ASSOCIATED(forces),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
    nsize_wrk = coords%size_wrk(1)
    ! All methods but not the classical elastic band will skip the force
    ! calculation for the last frame of the band
    SELECT CASE(my_tag)
    CASE(do_b_neb, do_it_neb, do_ci_neb, do_d_neb)
       IF (i==neb_env%number_of_replica) RETURN
    CASE(do_sm)
       ! String Method
       ! The forces do not require any projection. Reparametrization required
       ! after the update of the replica.
       CALL cite_reference(E2002)
       RETURN
    END SELECT
    ! otherwise proceeed normally..
    ALLOCATE (wrk(nsize_wrk),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Spring Energy
    CALL neb_replica_distance(coords=coords,i0=i-1,i=i,distance=distance1,iw=iw,&
         rotate=.FALSE.,error=error)
    tmp   = distance1-neb_env%avg_distance
    neb_env%spring_energy = neb_env%spring_energy + 0.5_dp*neb_env%k*tmp**2
    SELECT CASE(my_tag)
    CASE(do_eb)
       CALL cite_reference(Elber1987)
       ! Elastic band - Hamiltonian formulation according the original Karplus/Elber
       !                formulation
       ALLOCATE(dtmp1(nsize_wrk),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! derivatives of the spring
       tmp   = distance1-neb_env%avg_distance
       dtmp1 = 1.0_dp/distance1 * (coords%wrk(:,i)-coords%wrk(:,i-1))
       wrk   = neb_env%k*tmp*dtmp1
       forces%wrk(:,i  ) = forces%wrk(:,i  ) - wrk
       forces%wrk(:,i-1) = forces%wrk(:,i-1) + wrk
       ! derivatives due to the average length of the spring
       fac = 1.0_dp/(neb_env%avg_distance*REAL(neb_env%number_of_replica-1,KIND=dp))
       wrk = neb_env%k*fac*(coords%wrk(:,i)-coords%wrk(:,i-1))
       tmp = 0.0_dp
       DO j = 2, neb_env%number_of_replica
          CALL neb_replica_distance(coords=coords,i0=j-1,i=j,distance=distance1,iw=iw,&
               rotate=.FALSE.,error=error)
          tmp = tmp + distance1-neb_env%avg_distance
       END DO
       forces%wrk(:,i)   = forces%wrk(:,i)   + wrk * tmp
       forces%wrk(:,i-1) = forces%wrk(:,i-1) - wrk * tmp
       DEALLOCATE(dtmp1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE(do_b_neb)
       ! Bisection NEB
       CALL cite_reference(Jonsson1998)
       wrk = (coords%wrk(:,i+1)-2.0_dp*coords%wrk(:,i)+coords%wrk(:,i-1))
       tmp = neb_env%k * DOT_PRODUCT(wrk,tangent)
       wrk = forces%wrk(:,i)-dot_product_band(neb_env,forces%wrk(:,i),tangent,Mmatrix,error)*tangent
       forces%wrk(:,i) = wrk + tmp * tangent
    CASE(do_it_neb)
       ! Improved tangent NEB
       CALL cite_reference(Jonsson2000_1)
       CALL neb_replica_distance(coords=coords,i0=i,i=i+1,distance=distance1,iw=iw,&
            rotate=.FALSE.,error=error)
       CALL neb_replica_distance(coords=coords,i0=i-1,i=i,distance=distance2,iw=iw,&
            rotate=.FALSE.,error=error)
       tmp = neb_env%k * (distance1-distance2)
       wrk = forces%wrk(:,i)-dot_product_band(neb_env,forces%wrk(:,i),tangent,Mmatrix,error)*tangent
       forces%wrk(:,i) = wrk + tmp * tangent
    CASE(do_ci_neb)
       ! Climbing Image NEB
       CALL cite_reference(Jonsson2000_2)
       IF (neb_env%istep<=neb_env%nsteps_it.OR.i/=neb_env%nr_HE_image) THEN
          CALL get_neb_force(neb_env, tangent, coords, i, forces, do_it_neb, Mmatrix, iw, error)
       ELSE
          wrk = forces%wrk(:,i)
          tmp = -2.0_dp*dot_product_band(neb_env,wrk,tangent,Mmatrix,error=error)
          forces%wrk(:,i) = wrk + tmp * tangent
       END IF
    CASE(do_d_neb)
       ! Doubly NEB
       CALL cite_reference(Wales2004)
       ALLOCATE(dtmp1(nsize_wrk),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       dtmp1       = forces%wrk(:,i)-dot_product_band(neb_env,forces%wrk(:,i),tangent,Mmatrix,error)*tangent
       forces%wrk(:,i) = dtmp1
       dtmp1 = dtmp1 / SQRT(DOT_PRODUCT(dtmp1,dtmp1))
       ! Project out only the spring component interfering with the
       ! orthogonal gradient of the band
       wrk   = (coords%wrk(:,i+1)-2.0_dp*coords%wrk(:,i)+coords%wrk(:,i-1))
       dtmp1 = neb_env%k * (wrk - DOT_PRODUCT(wrk,dtmp1)*dtmp1)
       forces%wrk(:,i) = forces%wrk(:,i) + dtmp1
       DEALLOCATE(dtmp1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_neb_force

! *****************************************************************************
!> \brief  Handles the dot_product when using colvar.. in this case
!>         the scalar product needs to take into account the metric
!>         tensor
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  FUNCTION dot_product_band(neb_env, array1, array2, array3, error) RESULT(value)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: array1, array2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array3
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'dot_product_band', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nsize_int
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    IF (neb_env%use_colvar) THEN
       nsize_int = neb_env%nsize_int
       check =((SIZE(array1)/=SIZE(array2)).OR.&
               (SIZE(array1)/=nsize_int).OR.&
               (SIZE(array3)/=nsize_int*nsize_int))
       ! This condition should always be satisfied..
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       value = DOT_PRODUCT(MATMUL(RESHAPE(array3,(/nsize_int,nsize_int/)),array1),array2)
    ELSE
       value = DOT_PRODUCT(array1,array2)
    END IF
  END FUNCTION dot_product_band

! *****************************************************************************
!> \brief Reorient iteratively all images of the NEB chain in order to
!>      have always the smaller RMSD between two following images
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE reorient_images(rotate_frames, particle_set, coords, vels, iw,&
       distances, number_of_replica, error)
    LOGICAL, INTENT(IN)                      :: rotate_frames
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(neb_var_type), POINTER              :: coords, vels
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: distances
    INTEGER, INTENT(IN)                      :: number_of_replica
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reorient_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, kind
    LOGICAL                                  :: check, failure
    REAL(KIND=dp)                            :: xtmp
    REAL(KIND=dp), DIMENSION(3)              :: tmp
    REAL(KIND=dp), DIMENSION(3, 3)           :: rot

    failure = .FALSE.
    IF (.NOT.failure) THEN
       rot      = 0.0_dp
       rot(1,1) = 1.0_dp
       rot(2,2) = 1.0_dp
       rot(3,3) = 1.0_dp
       DO i = 2, number_of_replica
          ! The rotation of the replica is enabled exclusively when working in
          ! cartesian coordinates
          IF (rotate_frames.AND.(coords%in_use==do_band_cartesian)) THEN
             CALL rmsd3(particle_set, coords%xyz(:,i), coords%xyz(:,i-1), iw,&
                  rotate=.TRUE., rot=rot, error=error)
             ! Rotate velocities
             DO k = 1, SIZE(vels%xyz,1)/3
                kind = (k-1)*3
                tmp  = vels%xyz(kind+1:kind+3,i)
                CALL matvec_3x3(vels%xyz(kind+1:kind+3,i),TRANSPOSE(rot),tmp)
             END DO
           END IF
          IF (PRESENT(distances)) THEN
             check = SIZE(distances)==(number_of_replica-1)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)
             xtmp=DOT_PRODUCT(coords%wrk(:,i)-coords%wrk(:,i-1),&
                              coords%wrk(:,i)-coords%wrk(:,i-1))
             distances(i-1)=SQRT(xtmp)
          END IF
       END DO
    END IF
  END SUBROUTINE reorient_images

! *****************************************************************************
!> \brief Reparametrization of the replica for String Method with splines
!> \author Teodoro Laino - Rodolphe Vuilleumier 09.2008
! *****************************************************************************
  SUBROUTINE reparametrize_images(reparametrize_frames, spline_order, smoothing,&
       coords, sline, distances, error)

    LOGICAL, INTENT(IN)                      :: reparametrize_frames
    INTEGER, INTENT(IN)                      :: spline_order
    REAL(KIND=dp), INTENT(IN)                :: smoothing
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, sline
    REAL(KIND=dp), DIMENSION(:)              :: distances
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reparametrize_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: avg_distance, xtmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_coords

    failure = .FALSE.
    IF (reparametrize_frames) THEN
       ALLOCATE(tmp_coords(SIZE(coords,1),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_coords = coords
       ! Smoothing
       DO i = 2, SIZE(coords,2)-1
          coords(:,i) = tmp_coords(:,i)*(1.0_dp-2.0_dp*smoothing) + &
               tmp_coords(:,i-1)*smoothing + tmp_coords(:,i+1)*smoothing
       END DO
       sline = coords - tmp_coords + sline
       tmp_coords =  coords
       ! Reparametrization
       SELECT CASE(spline_order)
       CASE(1)
          ! Compute distances
          DO i = 2, SIZE(coords,2)
             xtmp = DOT_PRODUCT(coords(:,i)-coords(:,i-1),coords(:,i)-coords(:,i-1))
             distances(i-1) = SQRT(xtmp)
          END DO
          avg_distance = SUM(distances)/REAL(SIZE(coords,2)-1,KIND=dp)
          ! Redistribute frames
          DO i = 2, SIZE(coords,2)-1
             xtmp  = 0.0_dp
             DO j = 1, SIZE(coords,2)-1
                xtmp = xtmp + distances(j)
                IF (xtmp>avg_distance*REAL(i-1,KIND=dp)) THEN
                   xtmp = (xtmp - avg_distance*REAL(i-1,KIND=dp))/distances(j)
                   coords(:,i) = (1.0_dp-xtmp)*tmp_coords(:,j+1) + xtmp*tmp_coords(:,j)
                   EXIT
                END IF
             END DO
          END DO
          ! Re-compute distances
          DO i = 2, SIZE(coords,2)
             xtmp = DOT_PRODUCT(coords(:,i)-coords(:,i-1),coords(:,i)-coords(:,i-1))
             distances(i-1) = SQRT(xtmp)
          END DO
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="String Method: Spline order greater than 1 not implemented.", &
               error=error, error_level=cp_warning_level)
       END SELECT
       sline = coords - tmp_coords + sline
       DEALLOCATE(tmp_coords,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE reparametrize_images

! *****************************************************************************
!> \brief Checks for convergence criteria during a NEB run
!> \author Teodoro Laino 10.2006
! *****************************************************************************
  FUNCTION check_convergence(neb_env, Dcoords, forces, error) RESULT(converged)
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(neb_var_type), POINTER              :: Dcoords, forces
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: converged

    CHARACTER(len=*), PARAMETER :: routineN = 'check_convergence', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3), DIMENSION(4)           :: labels
    INTEGER                                  :: iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: max_dr, max_force, my_max_dr, &
                                                my_max_force, my_rms_dr, &
                                                my_rms_force, rms_dr, &
                                                rms_force
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: cc_section

    failure = .FALSE.
    NULLIFY(logger,cc_section)
    logger => cp_error_get_logger(error)
    cc_section => section_vals_get_subs_vals(neb_env%neb_section,"CONVERGENCE_CONTROL",error=error)
    CALL section_vals_val_get(cc_section,"MAX_DR",r_val=max_dr,error=error)
    CALL section_vals_val_get(cc_section,"MAX_FORCE",r_val=max_force,error=error)
    CALL section_vals_val_get(cc_section,"RMS_DR",r_val=rms_dr,error=error)
    CALL section_vals_val_get(cc_section,"RMS_FORCE",r_val=rms_force,error=error)
    converged    = .FALSE.
    labels       = " NO"
    my_max_dr    = MAXVAL(ABS(Dcoords%wrk))
    my_max_force = MAXVAL(ABS(forces%wrk))
    my_rms_dr    = SQRT(SUM(Dcoords%wrk*Dcoords%wrk)/REAL(SIZE(Dcoords%wrk,1)*SIZE(Dcoords%wrk,2),KIND=dp))
    my_rms_force = SQRT(SUM(forces%wrk*forces%wrk)/REAL(SIZE(forces%wrk,1)*SIZE(forces%wrk,2),KIND=dp))
    IF ( my_max_dr    < max_dr    ) labels(1) = "YES"
    IF ( my_max_force < max_force ) labels(2) = "YES"
    IF ( my_rms_dr    < rms_dr    ) labels(3) = "YES"
    IF ( my_rms_force < rms_force ) labels(4) = "YES"
    IF (ALL(labels=="YES")) converged = .TRUE.

    iw=cp_print_key_unit_nr(logger,neb_env%neb_section,"CONVERGENCE_INFO",&
         extension=".nebLog",error=error)
    IF (iw>0) THEN
       ! Print convergence info
       WRITE(iw,FMT='(A,A)')' **************************************', &
               '*****************************************'
       WRITE(iw,FMT='(1X,A,2X,F8.5,5X,"[",F8.5,"]",1X,T76,"(",A,")")')&
            'RMS DISPLACEMENT =',my_rms_dr,rms_dr,labels(3),&
            'MAX DISPLACEMENT =',my_max_dr,max_dr,labels(1),&
            'RMS FORCE        =',my_rms_force,rms_force,labels(4),&
            'MAX FORCE        =',my_max_force,max_force,labels(2)
       WRITE(iw,FMT='(A,A)')' **************************************', &
            '*****************************************'
    END IF
    CALL cp_print_key_finished_output(iw,logger,neb_env%neb_section,&
         "CONVERGENCE_INFO", error=error)
  END FUNCTION check_convergence

END MODULE neb_utils
