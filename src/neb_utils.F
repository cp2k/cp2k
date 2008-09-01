!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
  USE bibliography,                    ONLY: Elber1987,&
                                             Jonsson1998,&
                                             Jonsson2000_1,&
                                             Jonsson2000_2,&
                                             Wales2004,&
                                             cite_reference
  USE convert_units,                   ONLY: convert_to_cp2k_units
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
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: band_diis_opt,&
                                             do_b_neb,&
                                             do_ci_neb,&
                                             do_d_neb,&
                                             do_eb,&
                                             do_it_neb,&
                                             do_sm
  USE input_cp2k_check,                ONLY: remove_restart_info
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathlib,                         ONLY: matvec_3x3
  USE md_run,                          ONLY: qs_mol_dyn
  USE neb_IO,                          ONLY: dump_replica_coordinates
  USE neb_md_utils,                    ONLY: neb_initialize_velocity
  USE neb_types,                       ONLY: neb_type
  USE particle_types,                  ONLY: particle_type
  USE replica_types,                   ONLY: rep_env_calc_e_f,&
                                             replica_env_type
  USE rmsd,                            ONLY: rmsd3
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_utils'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: neb_replica_distance,&
            build_replica_coords,&
            neb_calc_energy_forces,&
            get_tangent,&
            get_neb_force,&
            reorient_images,&
            check_convergence

CONTAINS

! *****************************************************************************
!> \brief Computes the distance between two replica
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_replica_distance(particle_set, coords, i0, i, distance, iw, rotate, error)
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    INTEGER, INTENT(IN)                      :: i0, i
    REAL(KIND=dp), INTENT(OUT)               :: distance
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: rotate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_replica_distance', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_rotate

    failure = .FALSE.
    IF (.NOT.failure) THEN
       my_rotate = .FALSE.
       IF (PRESENT(rotate)) my_rotate=rotate
       IF (my_rotate) THEN
          CPPostcondition(PRESENT(particle_set),cp_failure_level,routineP,error,failure)
          CALL rmsd3(particle_set, coords(:,i), coords(:,i0), iw,&
               rotate=my_rotate, error=error)
       END IF
       distance = SQRT(DOT_PRODUCT(coords(:,i)-coords(:,i0),coords(:,i)-coords(:,i0)))
    END IF
  END SUBROUTINE neb_replica_distance

! *****************************************************************************
!> \brief Constructs or Read the coordinates for all replica
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE build_replica_coords(neb_section, particle_set,&
       coords, vels, neb_env, ndim, iw, globenv, para_env, error)
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, vels
    TYPE(neb_type), POINTER                  :: neb_env
    INTEGER, INTENT(IN)                      :: ndim, iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_replica_coords', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: dummy_char
    INTEGER :: handle, i, input_nr_replica, j, jtarg, natom, neb_nr_replica, &
      nr_replica_to_interpolate, stat
    LOGICAL                                  :: explicit, failure, my_end, &
                                                skip_vel_section
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distance
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rptr
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(section_vals_type), POINTER         :: coord_section, &
                                                replica_section, vel_section

    NULLIFY(parser)
    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
       CPPostcondition(.NOT.ASSOCIATED(vels),cp_failure_level,routineP,error,failure)
       neb_nr_replica = neb_env%number_of_replica
       replica_section => section_vals_get_subs_vals(neb_section,"REPLICA",error=error)
       CALL section_vals_get(replica_section,n_repetition=input_nr_replica,error=error)
       IF (input_nr_replica<=neb_nr_replica) THEN
          ALLOCATE(coords(ndim,neb_nr_replica),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(vels(ndim,neb_nr_replica),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          coords = 0.0_dp
          vels   = 0.0_dp
       ELSE
          ! More input replicas then requested for NEB.. abort
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       ! Read in replicas coordinates
       skip_vel_section = (input_nr_replica /= neb_nr_replica)
       IF (iw>0.AND.skip_vel_section) THEN
          WRITE(iw,'(T2,A)')'NEB| The number of replica in input is different from  the number of',&
                            'NEB| replica requested for NEB. More Replica  will  be  interpolated',&
                            'NEB| Therefore the possibly provided velocities will not be read.'
       END IF
       DO i = 1, input_nr_replica
          vel_section => section_vals_get_subs_vals(replica_section,"VELOCITY",i_rep_section=i,error=error)
          CALL section_vals_get(vel_section,explicit=explicit,error=error)
          skip_vel_section = skip_vel_section.OR.(.NOT.explicit)
       END DO
       DO i = 1, input_nr_replica
          coord_section => section_vals_get_subs_vals(replica_section,"COORD",i_rep_section=i,error=error)
          CALL section_vals_get(coord_section,explicit=explicit,error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                  n_rep_val=natom,error=error)
             CPPostcondition(natom*3==ndim,cp_failure_level,routineP,error,failure)
             DO j = 1, natom
                CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                     i_rep_val=j,r_vals=rptr,error=error)
               coords((j-1)*3+1:(j-1)*3+3,i) = rptr(1:3)
             END DO
          ELSE
             CALL section_vals_val_get(replica_section,"COORD_FILE_NAME",&
                  i_rep_section=i,c_val=filename,error=error)
             CPPostcondition(TRIM(filename)/="",cp_failure_level,routineP,error,failure)
             CALL parser_create(parser,filename,para_env=para_env,error=error)
             CALL parser_get_next_line(parser,1,error=error)
             ! Start parser             
             CALL parser_get_object(parser,natom,error=error)
             CALL parser_get_next_line(parser,1,error=error)
             DO j = 1, natom
                ! Atom coordinates
                CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
                CALL cp_assert(.NOT.my_end,cp_fatal_level,cp_assertion_failed,routineP,&
                     "Number of lines in XYZ format not equal to the number of atoms."//&
                     " Error in XYZ format for REPLICA coordinates. Very probably the"//&
                     " line with title is missing or is empty. Please check the XYZ file and rerun your job!"//&
CPSourceFileRef,&
                only_ionode=.TRUE.)
                READ(parser%input_line,*)dummy_char,coords((j-1)*3+1:(j-1)*3+3,i)
                CALL convert_to_cp2k_units("ANGSTROM",length=coords((j-1)*3+1,i))
                CALL convert_to_cp2k_units("ANGSTROM",length=coords((j-1)*3+2,i))
                CALL convert_to_cp2k_units("ANGSTROM",length=coords((j-1)*3+3,i))
             END DO
             CALL parser_release(parser,error=error)
          END IF
          CALL dump_replica_coordinates(particle_set, coords, i, i, iw)
          IF (.NOT.skip_vel_section) THEN
             vel_section => section_vals_get_subs_vals(replica_section,"VELOCITY",i_rep_section=i,&
                  error=error)
             CALL section_vals_val_get(vel_section,"_DEFAULT_KEYWORD_",&
                  n_rep_val=natom,error=error)
             CPPostcondition(natom*3==ndim,cp_failure_level,routineP,error,failure)
             DO j = 1, natom
                CALL section_vals_val_get(vel_section,"_DEFAULT_KEYWORD_",&
                     i_rep_val=j,r_vals=rptr,error=error)
                vels((j-1)*3+1:(j-1)*3+3,i) = rptr(1:3)
             END DO
          ELSE
             CALL neb_initialize_velocity(vels, neb_section, particle_set, i, iw, globenv,&
                  error)
          END IF
          IF (iw>0) THEN
             WRITE(iw,'(T2,"NEB|",75("*"),/)')
          END IF
       END DO
       ALLOCATE(distance(neb_nr_replica-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (input_nr_replica < neb_nr_replica) THEN
          ! Interpolate missing replicas
          nr_replica_to_interpolate = neb_nr_replica - input_nr_replica
          distance = 0.0_dp
          IF (iw>0) THEN
             WRITE(iw,'(T2,A,I5,A)')'NEB| Interpolating ',&
                  nr_replica_to_interpolate,' missing Replica.'
          END IF
          DO WHILE (nr_replica_to_interpolate >0)
             ! Compute distances between known images to find the interval
             ! where to add a new image
             DO j=1, input_nr_replica-1
                CALL neb_replica_distance(particle_set, coords,j,j+1,distance(j),iw,&
                     rotate=.TRUE.,error=error)
             END DO
             jtarg = MAXLOC(distance(1:input_nr_replica),1)
             IF (iw>0) THEN
                WRITE(iw,'(T2,3(A,I5),A)')'NEB| Interpolating Nr. ',&
                  nr_replica_to_interpolate,' missing Replica between Replica Nr.',&
                  jtarg,' and',jtarg+1,'.'
             END IF
             input_nr_replica = input_nr_replica + 1
             nr_replica_to_interpolate = nr_replica_to_interpolate - 1
             ! Interpolation is a simple bisection
             coords(:,jtarg+2:input_nr_replica) = coords(:,jtarg+1:input_nr_replica-1)
             coords(:,jtarg+1) = (coords(:,jtarg)+coords(:,jtarg+2))/2.0_dp

             vels(:,jtarg+2:input_nr_replica) = vels(:,jtarg+1:input_nr_replica-1)
             vels(:,jtarg+1) = 0.0_dp
             CALL dump_replica_coordinates(particle_set, coords, jtarg+1, &
                  input_nr_replica, iw)
             CALL neb_initialize_velocity(vels, neb_section, particle_set,&
                  jtarg+1, iw, globenv, error)
          END DO
       END IF
       vels(:,1)            = 0.0_dp
       vels(:,SIZE(vels,2)) = 0.0_dp
       ! If we perform a DIIS optimization we don't need velocities
       IF (neb_env%opt_type==band_diis_opt) vels = 0.0_dp
       ! Compute distances between replicas
       DO j=1, input_nr_replica-1
          CALL neb_replica_distance(particle_set, coords,j,j+1,distance(j),iw,&
               rotate=.TRUE.,error=error)
       END DO
       DEALLOCATE(distance,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
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
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: energies
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_calc_energy_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1), DIMENSION(3)           :: lab = (/"X","Y","Z"/)
    INTEGER                                  :: handle, i, irep, j, n_rep, &
                                                n_rep_neb, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tangent, tmp_a, tmp_b

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       n_rep_neb = SIZE(coords,2)
       n_rep     = rep_env%nrep
       energies  = 0.0_dp
       IF (output_unit>0) WRITE(output_unit,'(/,T2,A)')"NEB| Computing Energies and Forces"
       DO irep = 1, n_rep_neb, n_rep
          DO j = 0, n_rep-1
             IF (irep+j <= n_rep_neb) THEN
                ! If the number of replica in replica_env and the number of replica
                ! used in the NEB does not match, the other replica in replica_env
                ! just compute energies and forces keeping the fixed coordinates and
                ! forces
                rep_env%r(:,j+1)=coords(:,irep+j)
             END IF
          END DO
          CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
          ! Free Energy case.. Work in progress - teo 01.2007
          ! CALL perform_replica_md(rep_env, coords, irep, n_rep_neb, error)
          DO j = 0, n_rep-1
             IF (irep+j <= n_rep_neb) THEN
                forces(:,irep+j) = rep_env%f(1:rep_env%ndim,j+1)
                energies(irep+j) = rep_env%f(rep_env%ndim+1,j+1)
                ! Dump Info
                IF (output_unit>0) THEN
                   WRITE(output_unit,'(T2,A,I5,A,I5,A)')&
                        "NEB| REPLICA Nr.",irep+j,"- Energy and Forces"
                   !
                   WRITE(output_unit,'(T2,A,T42,A,9X,F15.9)')&
                        "NEB|"," Total Energy: ",rep_env%f(rep_env%ndim+1,j+1)
                   WRITE(output_unit,'(T2,"NEB|",T10,"ATOM",T33,3(9X,A,7X))')lab(1),lab(2),lab(3)
                   DO i = 1, SIZE(particle_set)
                      WRITE(output_unit,'(T2,"NEB|",T12,A,T30,3(2X,F15.9))')&
                           particle_set(i)%atomic_kind%name,&
                           rep_env%f((i-1)*3+1:(i-1)*3+3,j+1)
                   END DO
                END IF
             END IF
          END DO
       END DO
       IF (PRESENT(neb_env)) THEN
          ! First identify the image of the chain with the higher potential energy
          neb_env%nr_HE_image = MAXLOC(energies,1)
          ALLOCATE(tangent(SIZE(coords,1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Then modify image forces accordingly to the scheme chosen for the
          ! calculation.
          neb_env%spring_energy = 0.0_dp
          IF (neb_env%optimize_end_points) THEN
             ALLOCATE(tmp_a(SIZE(forces,1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp_b(SIZE(forces,1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             tmp_a = forces(:,1)
             tmp_b = forces(:,SIZE(forces,2))
          END IF
          DO i = 2, SIZE(coords,2)
             CALL get_tangent(neb_env,coords,i,tangent,energies,output_unit,error)
             CALL get_neb_force(neb_env,tangent,coords,i,forces,iw=output_unit,error=error)
          END DO
          IF (neb_env%optimize_end_points) THEN
             forces(:,1)              = tmp_a ! Image A
             forces(:,SIZE(forces,2)) = tmp_b ! Image B             
             DEALLOCATE(tmp_a,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(tmp_b,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)             
          ELSE
             ! Nullify forces on the two end points images
             forces(:,1)              = 0.0_dp ! Image A
             forces(:,SIZE(forces,2)) = 0.0_dp ! Image B
          END IF
          DEALLOCATE(tangent,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE neb_calc_energy_forces

! *****************************************************************************
!> \brief Driver to perform an MD run on each single replica to 
!>      compute specifically Free Energies in a NEB scheme
!> \author Teodoro Laino  01.2007
! *****************************************************************************
  SUBROUTINE perform_replica_md(rep_env, irep, n_rep_neb, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    INTEGER, INTENT(IN)                      :: irep, n_rep_neb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perform_replica_md', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, i, ierr, j
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure, handle=handle2)
       logger => cp_error_get_logger(new_error)
       CALL force_env_get(f_env%force_env,globenv=globenv,&
            root_section=root_section,error=new_error)
       j=rep_env%local_rep_indices(1)-1
       ! Syncronize position on the replica procs
       CALL cp_set_pos(rep_env%f_env_id,rep_env%r(:,j+1),&
            3*rep_env%nat,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       !
       IF (irep+j <= n_rep_neb) THEN
          logger%iter_info%iteration(2)=irep+j
          CALL remove_restart_info(root_section,error=new_error)
          CALL section_vals_val_set(root_section,"MOTION%MD%STEPS",i_val=(irep+j)*1,&
               error=new_error)
          CALL qs_mol_dyn(f_env%force_env,globenv=globenv,error=new_error)
       END IF
       CALL f_env_rm_defaults(f_env=f_env,error=new_error,ierr=ierr,handle=handle2)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       ! Copy Back Free Energies forces
       ! Work in progress - teo 01.2007
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE perform_replica_md

! *****************************************************************************
!> \brief Computes the tangent for point i of the NEB chain
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE get_tangent(neb_env, coords, i, tangent, energies, iw, error)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
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
    IF (i==SIZE(coords,2)) RETURN
    ! Several kind of tangents implemented...
    SELECT CASE(neb_env%id_type)
    CASE(do_eb)
       tangent(:) = 0.0_dp
    CASE(do_b_neb)
       CALL neb_replica_distance(coords=coords,i0=i,i=i-1,distance=distance1,iw=iw,&
            rotate=.FALSE.,error=error)
       CALL neb_replica_distance(coords=coords,i0=i+1,i=i,distance=distance2,iw=iw,&
            rotate=.FALSE.,error=error)
       tangent(:) = (coords(:,i)  -coords(:,i-1))/distance1 +&
                    (coords(:,i+1)-coords(:,i)  )/distance2
    CASE(do_it_neb, do_ci_neb, do_d_neb)
       IF      ((energies(i+1).GT.energies(i)).AND.(energies(i).GT.(energies(i-1)))) THEN
          tangent(:) = coords(:,i+1) - coords(:,i)
       ELSE IF ((energies(i+1).LT.energies(i)).AND.(energies(i).LT.(energies(i-1)))) THEN
          tangent(:) = coords(:,i)   - coords(:,i-1)
       ELSE 
          DVmax = MAX(ABS(energies(i+1)-energies(i)),ABS(energies(i-1)-energies(i)))
          DVmin = MIN(ABS(energies(i+1)-energies(i)),ABS(energies(i-1)-energies(i)))
          IF (energies(i+1) .GE. energies(i-1)) THEN
             tangent(:) = ( coords(:,i+1) - coords(:,i))*DVmax+(coords(:,i) - coords(:,i-1))*DVmin
          ELSE
             tangent(:) = ( coords(:,i+1) - coords(:,i))*DVmin+(coords(:,i) - coords(:,i-1))*DVmax
          END IF
       END IF
    CASE(do_sm)
       STOP
       ! String method..
    END SELECT
    distance0  = SQRT(DOT_PRODUCT(tangent(:),tangent(:)))
    IF (distance0 /= 0.0_dp) tangent(:) = tangent(:) / distance0
  END SUBROUTINE get_tangent

! *****************************************************************************
!> \brief Computes the forces for point i of the NEB chain
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE get_neb_force(neb_env,tangent,coords,i,forces,tag,iw,error)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tangent
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    INTEGER, INTENT(IN), OPTIONAL            :: tag
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neb_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, k, my_tag, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: distance1, distance2, fac, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dtmp1, wrk

    failure = .FALSE.
    my_tag  = neb_env%id_type
    IF (PRESENT(tag)) my_tag = tag
    CPPostcondition(ASSOCIATED(forces),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
    ! All methods but not the classical elastic band will skip the force 
    ! calculation for the last frame of the band
    SELECT CASE(my_tag)
    CASE(do_b_neb, do_it_neb, do_ci_neb, do_d_neb)
       IF (i==SIZE(coords,2)) RETURN
    CASE(do_sm)
       ! Not yet implemented
       STOP
    END SELECT
    ! otherwise proceeed normally..
    ALLOCATE(wrk(SIZE(coords,1)),stat=stat)
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
       ALLOCATE(dtmp1(SIZE(coords,1)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! derivatives of the spring
       tmp   = distance1-neb_env%avg_distance
       dtmp1 = 1.0_dp/distance1 * (coords(:,i)-coords(:,i-1))
       wrk           = neb_env%k*tmp*dtmp1
       forces(:,i  ) = forces(:,i  ) - wrk
       forces(:,i-1) = forces(:,i-1) + wrk
       ! derivatives due to the average length of the spring
       fac = 1.0_dp/(neb_env%avg_distance*REAL(SIZE(coords,2)-1,KIND=dp))
       wrk = neb_env%k*fac*(coords(:,i)-coords(:,i-1))
       tmp = 0.0_dp
       DO j = 2, SIZE(coords,2)
          CALL neb_replica_distance(coords=coords,i0=j-1,i=j,distance=distance1,iw=iw,&
               rotate=.FALSE.,error=error)
          tmp = tmp + distance1-neb_env%avg_distance
       END DO
       forces(:,i)   = forces(:,i)   + wrk * tmp
       forces(:,i-1) = forces(:,i-1) - wrk * tmp
       DEALLOCATE(dtmp1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE(do_b_neb)
       ! Bisection NEB
       CALL cite_reference(Jonsson1998)
       wrk = (coords(:,i+1)-2.0_dp*coords(:,i)+coords(:,i-1))
       tmp = neb_env%k * DOT_PRODUCT(wrk,tangent)
       wrk = forces(:,i)-DOT_PRODUCT(forces(:,i),tangent)*tangent
       forces(:,i) = wrk + tmp * tangent
    CASE(do_it_neb)
       ! Improved tangent NEB
       CALL cite_reference(Jonsson2000_1)
       CALL neb_replica_distance(coords=coords,i0=i,i=i+1,distance=distance1,iw=iw,&
            rotate=.FALSE.,error=error)
       CALL neb_replica_distance(coords=coords,i0=i-1,i=i,distance=distance2,iw=iw,&
            rotate=.FALSE.,error=error)
       tmp = neb_env%k * (distance1-distance2)
       wrk = forces(:,i)-DOT_PRODUCT(forces(:,i),tangent)*tangent
       forces(:,i) = wrk + tmp * tangent       
    CASE(do_ci_neb)
       ! Climbing Image NEB
       CALL cite_reference(Jonsson2000_2)
       IF (neb_env%istep<=neb_env%nsteps_it.OR.i/=neb_env%nr_HE_image) THEN
          CALL get_neb_force(neb_env, tangent, coords, i, forces, do_it_neb, iw, error)
       ELSE
          wrk = forces(:,i) 
          tmp = -2.0_dp*DOT_PRODUCT(wrk,tangent)
          forces(:,i) = wrk + tmp * tangent
       END IF
    CASE(do_d_neb)
       ! Doubly NEB
       CALL cite_reference(Wales2004)
       ALLOCATE(dtmp1(SIZE(coords,1)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       dtmp1       = forces(:,i)-DOT_PRODUCT(forces(:,i),tangent)*tangent
       forces(:,i) = dtmp1
       dtmp1 = dtmp1 / SQRT(DOT_PRODUCT(dtmp1,dtmp1))
       ! Project out only the spring component interfering with the 
       ! orthogonal gradient of the band
       wrk   = (coords(:,i+1)-2.0_dp*coords(:,i)+coords(:,i-1))
       dtmp1 = neb_env%k * (wrk - DOT_PRODUCT(wrk,dtmp1)*dtmp1)
       forces(:,i) = forces(:,i) + dtmp1
       DEALLOCATE(dtmp1,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE(do_sm)
       ! Not yet implemented
       STOP
    END SELECT
    DEALLOCATE(wrk,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_neb_force

! *****************************************************************************
!> \brief Reorient iteratively all images of the NEB chain in order to
!>      have always the smaller RMSD between two following images
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE reorient_images(rotate_frames, particle_set, coords, vels, iw,&
       distances, error)
    LOGICAL, INTENT(IN)                      :: rotate_frames
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, vels
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: distances
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
       DO i = 2, SIZE(coords,2)
          IF (rotate_frames) THEN
             ! Rotate coordinates
             CALL rmsd3(particle_set, coords(:,i), coords(:,i-1), iw,&
                  rotate=.TRUE., rot=rot, error=error)
             ! Rotate velocities
             DO k = 1, SIZE(vels,1)/3
                kind = (k-1)*3
                tmp  = vels(kind+1:kind+3,i)
                CALL matvec_3x3(vels(kind+1:kind+3,i),TRANSPOSE(rot),tmp)
             END DO
           END IF
          IF (PRESENT(distances)) THEN
             check = SIZE(distances)==(SIZE(coords,2)-1)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)
             xtmp=DOT_PRODUCT(coords(:,i)-coords(:,i-1),coords(:,i)-coords(:,i-1))
             distances(i-1)=SQRT(xtmp)
          END IF
       END DO
    END IF
  END SUBROUTINE reorient_images

! *****************************************************************************
!> \brief Checks for convergence criteria during a NEB run
!> \author Teodoro Laino 10.2006
! *****************************************************************************
  FUNCTION check_convergence(neb_env, Dcoords, forces, error) RESULT(converged)
    TYPE(neb_type), POINTER                  :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Dcoords, forces
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
    my_max_dr    = MAXVAL(ABS(Dcoords))
    my_max_force = MAXVAL(ABS(forces))
    my_rms_dr    = SQRT(SUM(Dcoords*Dcoords)/REAL(SIZE(Dcoords,1)*SIZE(Dcoords,2),KIND=dp))
    my_rms_force = SQRT(SUM(forces*forces)/REAL(SIZE(forces,1)*SIZE(forces,2),KIND=dp))
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
