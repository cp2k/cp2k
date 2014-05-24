!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module performing a mdoe selective vibrational analysis
!> \note
!>      Numerical accuracy for parallel runs:
!>       Each replica starts the SCF run from the one optimized
!>       in a previous run. It may happen then energies and derivatives
!>       of a serial run and a parallel run could be slightly different
!>       'cause of a different starting density matrix.
!>       Exact results are obtained using:
!>          EXTRAPOLATION USE_GUESS in QS section (Teo 08.2006)
!> \author Florian Schiffmann 08.2006
! *****************************************************************************
MODULE mode_selective
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: get_results
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: ms_guess_atomic,&
                                             ms_guess_bfgs,&
                                             ms_guess_molden,&
                                             ms_guess_restart,&
                                             ms_guess_restart_vec
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             max_line_length
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_bcast
  USE molden_utils,                    ONLY: molden_out
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: massunit,&
                                             vibfac
  USE replica_methods,                 ONLY: rep_env_calc_e_f
  USE replica_types,                   ONLY: replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mode_selective'
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  TYPE ms_vib_type
     INTEGER                                  :: mat_size
     INTEGER                                  :: select_id
     INTEGER, DIMENSION(:), POINTER           :: inv_atoms
     REAL(KIND=dp)                            :: eps(2)
     REAL(KIND=dp)                            :: sel_freq
     REAL(KIND=dp),POINTER,DIMENSION(:, :)    :: b_vec
     REAL(KIND=dp),POINTER,DIMENSION(:, :)    :: delta_vec
     REAL(KIND=dp),POINTER,DIMENSION(:, :)    :: ms_force
     REAL(KIND=dp), DIMENSION(:), POINTER     :: eig_bfgs
     REAL(KIND=dp), DIMENSION(:), POINTER     :: f_range
     REAL(KIND=dp), DIMENSION(:), POINTER     :: inv_range
     REAL(KIND=dp), POINTER, DIMENSION(:)     :: step_b
     REAL(KIND=dp), POINTER, DIMENSION(:)     :: step_r
     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: b_mat
     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dip_deriv
     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hes_bfgs
     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: s_mat
     INTEGER                                  :: initial_guess
  END TYPE ms_vib_type

  PUBLIC :: ms_vb_anal,&
            ms_vib_type

CONTAINS
  ! *****************************************************************************
  !> \brief Module performing a vibrational analysis
  !> \author Teodoro Laino 08.2006
  ! *****************************************************************************
  SUBROUTINE ms_vb_anal(input,rep_env, para_env, globenv, particles,&
                        nrep,calc_intens,dx,output_unit,output_molden, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER                                  :: nrep
    LOGICAL                                  :: calc_intens
    REAL(KIND=dp)                            :: dx
    INTEGER                                  :: output_unit, output_molden
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ms_vb_anal', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: handle, i, ip1, j, natoms, &
                                                ncoord, stat
    LOGICAL                                  :: converged, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mass, pos0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp_deriv
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: tmp_dip
    TYPE(ms_vib_type)                        :: ms_vib

    CALL timeset(routineN,handle)
    failure = .FALSE.
    converged=.FALSE.
    natoms=SIZE(particles)
    ncoord=3*natoms
    ALLOCATE(mass(3*natoms), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, natoms
       DO j=1,3
          mass((i-1)*3+j) = particles(i)%atomic_kind%mass
          mass((i-1)*3+j)= SQRT(mass((i-1)*3+j))
       END DO
    END DO
    ! Allocate working arrays
    ALLOCATE(ms_vib%delta_vec(ncoord,nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%b_vec(ncoord,nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%step_r(nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%step_b(nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(calc_intens)THEN
       description='[DIPOLE]'
       ALLOCATE(tmp_dip(nrep,3,2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ms_vib%dip_deriv(3,nrep),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL MS_initial_moves( para_env, nrep , input, globenv,ms_vib,&
         particles,&
         mass,&
         dx,output_unit,output_molden,&
         calc_intens,error)
    ncoord=3*natoms
    ALLOCATE(pos0(ncoord), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%ms_force(ncoord,nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, natoms
       DO j=1,3
          pos0((i-1)*3+j) = particles((i))%r(j)
       END DO
    END DO
    ncoord=3*natoms
    DO
       ms_vib%ms_force = HUGE(0.0_dp)
       DO i=1,nrep
          DO j=1,ncoord
             rep_env%r(j,i)=pos0(j)+ms_vib%step_r(i)*ms_vib%delta_vec(j,i)
          END DO
       END DO
       CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)

       DO i = 1, nrep
          IF(calc_intens)THEN
             CALL get_results(results=rep_env%results(i)%results,&
                  description=description,&
                  n_rep=ip1,error=error)
             CALL get_results(results=rep_env%results(i)%results,&
                  description=description,&
                  values=tmp_dip(i,:,1),&
                  nval=ip1,error=error)
          END IF
          DO j = 1, ncoord
             ms_vib%ms_force(j,i) = rep_env%f(j,i)
          END DO
       END DO
       DO i=1,nrep
          DO j=1,ncoord
             rep_env%r(j,i)=pos0(j)-ms_vib%step_r(i)*ms_vib%delta_vec(j,i)
          END DO
       END DO
       CALL rep_env_calc_e_f(rep_env,calc_f=.TRUE.,error=error)
       IF(calc_intens)THEN
          DO i = 1, nrep
             CALL get_results(results=rep_env%results(i)%results,&
                  description=description,&
                  n_rep=ip1,error=error)
             CALL get_results(results=rep_env%results(i)%results,&
                  description=description,&
                  values=tmp_dip(i,:,2),&
                  nval=ip1,error=error)
             ms_vib%dip_deriv(:,ms_vib%mat_size+i)=(tmp_dip(i,:,1)-tmp_dip(i,:,2))/(2*ms_vib%step_b(i))
          END DO
       END IF

       CALL evaluate_H_update_b(rep_env,ms_vib,input,nrep,&
            particles,&
            mass,&
            converged,&
            dx,calc_intens,&
            output_unit,output_molden, error)
       IF(converged)EXIT
       IF(calc_intens)THEN
          ALLOCATE(tmp_deriv(3,ms_vib%mat_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          tmp_deriv=ms_vib%dip_deriv
          DEALLOCATE(ms_vib%dip_deriv,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(ms_vib%dip_deriv(3,ms_vib%mat_size+nrep),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ms_vib%dip_deriv(:,1:ms_vib%mat_size)=tmp_deriv(:,1:ms_vib%mat_size)
          DEALLOCATE(tmp_deriv,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    DEALLOCATE(ms_vib%ms_force, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pos0, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%step_r, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%step_b, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%b_vec, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%delta_vec, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mass, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%b_mat, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%s_mat, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(ms_vib%select_id==3)THEN
       DEALLOCATE(ms_vib%inv_atoms, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(ms_vib%eig_bfgs)) THEN
       DEALLOCATE(ms_vib%eig_bfgs, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(ms_vib%hes_bfgs)) THEN
       DEALLOCATE(ms_vib%hes_bfgs, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(calc_intens)THEN
       DEALLOCATE(ms_vib%dip_deriv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_dip,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ms_vb_anal
! *****************************************************************************
!> \brief Generates the first displacement vector for a mode selctive vibrational
!>      analysis. At the moment this is a random number for selected atoms
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE MS_initial_moves(para_env, nrep , input, globenv ,ms_vib,particles ,&
                              mass, dx,output_unit_ms,output_molden,&
                              calc_intens, error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: nrep
    TYPE(section_vals_type), POINTER         :: input
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(ms_vib_type)                        :: ms_vib
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(Kind=dp), DIMENSION(:)              :: mass
    REAL(KIND=dp)                            :: dx
    INTEGER                                  :: output_unit_ms, output_molden
    LOGICAL                                  :: calc_intens
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'MS_initial_moves', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: guess, handle, i, j, jj, k, &
                                                m, n_rep_val, natoms, ncoord, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_atoms
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: do_involved_atoms, failure, &
                                                ionode
    REAL(KIND=dp)                            :: my_val, norm
    TYPE(section_vals_type), POINTER         :: involved_at_section, &
                                                ms_vib_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       NULLIFY(ms_vib%eig_bfgs,ms_vib%f_range,ms_vib%hes_bfgs,ms_vib%inv_range)
       ms_vib_section => section_vals_get_subs_vals(input,"VIBRATIONAL_ANALYSIS%MODE_SELECTIVE", error=error)
       CALL section_vals_val_get(ms_vib_section,"INITIAL_GUESS",i_val=guess,error=error)
       CALL section_vals_val_get(ms_vib_section,"EPS_MAX_VAL",r_val=ms_vib%eps(1), error=error)
       CALL section_vals_val_get(ms_vib_section,"EPS_NORM",r_val=ms_vib%eps(2), error=error)
       CALL section_vals_val_get(ms_vib_section,"RANGE",n_rep_val=n_rep_val,error=error)
       ms_vib%select_id=0
       IF(n_rep_val.NE.0)THEN
          CALL section_vals_val_get(ms_vib_section,"RANGE",r_vals=ms_vib%f_range,error=error)
          IF(ms_vib%f_range(1).GT.ms_vib%f_range(2))THEN
             my_val=ms_vib%f_range(2)
             ms_vib%f_range(2)=ms_vib%f_range(1)
             ms_vib%f_range(1)=my_val
          END IF
          ms_vib%select_id=2
       END IF
       CALL section_vals_val_get(ms_vib_section,"FREQUENCY",r_val=ms_vib%sel_freq, error=error)
       IF(ms_vib%sel_freq.gt.0._dp)ms_vib%select_id=1
       involved_at_section=>section_vals_get_subs_vals(ms_vib_section,"INVOLVED_ATOMS",&
            error=error)
       CALL section_vals_get(involved_at_section,explicit=do_involved_atoms, error=error)
       IF(do_involved_atoms)THEN
          CALL section_vals_val_get(involved_at_section,"INVOLVED_ATOMS",n_rep_val=n_rep_val,error=error)
          jj = 0
          DO k = 1,n_rep_val
             CALL section_vals_val_get(involved_at_section,"INVOLVED_ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
             DO  j = 1,SIZE(tmplist)
                jj = jj+1
             END DO
          END DO
          IF (jj.GE.1)THEN
             natoms=jj
             ALLOCATE(ms_vib%inv_atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             jj = 0
             DO m=1,n_rep_val
                CALL section_vals_val_get(involved_at_section,"INVOLVED_ATOMS",i_rep_val=m,i_vals=tmplist, error=error)
                DO  j = 1,SIZE(tmplist)
                   ms_vib%inv_atoms (j) = tmplist(j)
                END DO
             END DO
             ms_vib%select_id=3
          END IF
          CALL section_vals_val_get(involved_at_section,"RANGE",n_rep_val=n_rep_val,error=error)
          IF(n_rep_val.NE.0)THEN
             CALL section_vals_val_get(involved_at_section,"RANGE",r_vals=ms_vib%inv_range,error=error)
             IF(ms_vib%inv_range(1).GT.ms_vib%inv_range(2))THEN
                ms_vib%inv_range(2)=my_val
                ms_vib%inv_range(2)=ms_vib%inv_range(1)
                ms_vib%inv_range(1)=my_val
             END IF
          END IF
       END IF
       IF(ms_vib%select_id==0)&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"no frequency, range or involved atoms specified ",&
            error=error,failure=failure)
       ionode=para_env%mepos==para_env%source
       SELECT CASE(guess)
       CASE(ms_guess_atomic)
          ms_vib%initial_guess=1
          CALL section_vals_val_get(ms_vib_section,"ATOMS",n_rep_val=n_rep_val,error=error)
          jj = 0
          DO k = 1,n_rep_val
             CALL section_vals_val_get(ms_vib_section,"ATOMS",i_rep_val=k,i_vals=tmplist, error=error)
             DO  j = 1,SIZE(tmplist)
                jj = jj+1
             END DO
          END DO
          IF (jj<1)THEN
             natoms=SIZE(particles)
             ALLOCATE(map_atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO j=1,natoms
                map_atoms(j)=j
             END DO
          ELSE
             natoms=jj
             ALLOCATE(map_atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             jj = 0
             DO m=1,n_rep_val
                CALL section_vals_val_get(ms_vib_section,"ATOMS",i_rep_val=m,i_vals=tmplist, error=error)
                DO  j = 1,SIZE(tmplist)
                   map_atoms (j) = tmplist(j)
                END DO
             END DO
          END IF

          ! apply random displacement along the mass weighted nuclear cartesian coordinates
          ms_vib%b_vec=0._dp
          ms_vib%delta_vec=0._dp
          jj=0

          DO i=1,nrep
             DO j=1,natoms
                DO k=1,3
                   jj=(map_atoms(j)-1)*3+k
                   ms_vib%b_vec(jj,i)=ABS(next_random_number(globenv%gaussian_rng_stream,error=error))
                END DO
             END DO
             norm=SQRT(DOT_PRODUCT(ms_vib%b_vec(:,i),ms_vib%b_vec(:,i)))
             ms_vib%b_vec(:,i)=ms_vib%b_vec(:,i)/norm
          END DO

          IF(nrep.GT.1)THEN
             DO k=1,10
                DO j=1,nrep
                   DO i=1,nrep
                      IF(i.NE.j)THEN
                         ms_vib%b_vec(:,j)=ms_vib%b_vec(:,j)-DOT_PRODUCT(ms_vib%b_vec(:,j),ms_vib%b_vec(:,i))* ms_vib%b_vec(:,i)
                         ms_vib%b_vec(:,j)=ms_vib%b_vec(:,j)/SQRT(DOT_PRODUCT(ms_vib%b_vec(:,j),ms_vib%b_vec(:,j)))
                      END IF
                   END DO
                END DO
             END DO
          END IF

          ms_vib%mat_size=0
          DO i=1,SIZE(ms_vib%b_vec,1)
             ms_vib%delta_vec(i,:)=ms_vib%b_vec(i,:)/mass(i)
          END DO
       CASE(ms_guess_bfgs)

          ms_vib%initial_guess=2
          CALL bfgs_guess(ms_vib_section,ms_vib,particles,mass,para_env,nrep,error)
          ms_vib%mat_size=0

       CASE(ms_guess_restart_vec)

          ms_vib%initial_guess=3
          ncoord=3*SIZE(particles)
          CALL rest_guess(ms_vib_section,para_env,ms_vib,mass,ionode,particles,nrep,calc_intens,error)


          ms_vib%mat_size=0
       CASE(ms_guess_restart)
          ms_vib%initial_guess=4
          ncoord=3*SIZE(particles)
          CALL rest_guess(ms_vib_section,para_env,ms_vib,mass,ionode,particles,nrep,calc_intens,error)

       CASE(ms_guess_molden)
          ms_vib%initial_guess=5
          ncoord=3*SIZE(particles)
          CALL molden_guess(ms_vib_section,input,para_env,ms_vib,mass,ncoord,nrep,error)
          ms_vib%mat_size=0
       END SELECT
       CALL mp_bcast(ms_vib%b_vec,para_env%source,para_env%group)
       CALL mp_bcast(ms_vib%delta_vec,para_env%source,para_env%group)
    END IF
    DO i=1,nrep
       ms_vib%step_r(i)=dx/SQRT(DOT_PRODUCT( ms_vib%delta_vec(:,i), ms_vib%delta_vec(:,i)))
       ms_vib%step_b(i)=SQRT(DOT_PRODUCT(ms_vib%step_r(i)*ms_vib%b_vec(:,i),ms_vib%step_r(i)*ms_vib%b_vec(:,i)))
    END DO
    CALL timestop(handle)

  END SUBROUTINE MS_initial_moves


! *****************************************************************************
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE bfgs_guess(ms_vib_section,ms_vib,particles,mass,para_env,nrep,error)

    TYPE(section_vals_type), POINTER         :: ms_vib_section
    TYPE(ms_vib_type)                        :: ms_vib
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(Kind=dp), DIMENSION(:)              :: mass
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: nrep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'bfgs_guess', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: hes_filename
    INTEGER                                  :: hesunit, i, j, jj, k, natoms, &
                                                ncoord, stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_val, norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp

    natoms=SIZE(particles)
    ncoord=3*natoms

    ALLOCATE(ms_vib%hes_bfgs(ncoord,ncoord),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%eig_bfgs(ncoord),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(para_env%mepos==para_env%source)THEN
       CALL section_vals_val_get(ms_vib_section,"RESTART_FILE_NAME",c_val=hes_filename,&
            error=error)
       IF(hes_filename=="") hes_filename="HESSIAN"
       CALL open_file(file_name=hes_filename,file_status="OLD",&
            file_form="UNFORMATTED", file_action="READ", unit_number=hesunit)
       ALLOCATE(tmp(ncoord),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmplist(ncoord),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! should use the cp_fm_read_unformatted...
       DO i=1,ncoord
          READ(UNIT=hesunit,IOSTAT=stat) ms_vib%hes_bfgs(:,i)
       ENDDO
       CALL close_file(hesunit)
       IF (stat /= 0) THEN
          WRITE(6,FMT="(/,T2,A)")"**  Error while reading HESSIAN **"
       ELSE
          WRITE(6,FMT="(/,T2,A)")&
               "*** Initial Hessian has been read successfully ***"
       END IF
       DO i = 1, ncoord
          DO j = 1, ncoord
             ms_vib%hes_bfgs(i,j) = ms_vib%hes_bfgs(i,j)/(mass(i)*mass(j))
          END DO
       END DO


       CALL diamat_all(ms_vib%hes_bfgs,ms_vib%eig_bfgs,error=error)
       tmp=0._dp
       IF(ms_vib%select_id==1) my_val=(ms_vib%sel_freq/vibfac)**2/massunit
       IF(ms_vib%select_id==2) my_val=(((ms_vib%f_range(2)+ms_vib%f_range(1))*0.5_dp)/vibfac)**2/massunit
       IF(ms_vib%select_id==1.OR.ms_vib%select_id==2)THEN
          DO i=1,ncoord
             tmp(i)=ABS(my_val-ms_vib%eig_bfgs(i))
          END DO
       ELSE IF(ms_vib%select_id==3)THEN
          DO i=1,ncoord
             DO j=1,SIZE(ms_vib%inv_atoms)
                DO k=1,3
                   jj=(ms_vib%inv_atoms(j)-1)*3+k
                   tmp(i)=tmp(i)+SQRT(ms_vib%hes_bfgs(jj,i)**2)
                END DO
             END DO
             IF((SIGN(1._dp,ms_vib%eig_bfgs(i))*SQRT(ABS(ms_vib%eig_bfgs(i))*massunit)*vibfac).LE.400._dp)tmp(i)=0._dp
          END DO
          tmp=-tmp
       END IF
       CALL sort(tmp,ncoord,tmplist)
       DO i=1,nrep
          ms_vib%b_vec(:,i)=ms_vib%hes_bfgs(:,tmplist(i))
          norm=SQRT(DOT_PRODUCT(ms_vib%b_vec(:,i),ms_vib%b_vec(:,i)))
          ms_vib%b_vec(:,i)=ms_vib%b_vec(:,i)/norm
       END DO
       DO i=1,SIZE(ms_vib%b_vec,1)
          ms_vib%delta_vec(i,:)=ms_vib%b_vec(i,:)/mass(i)
       END DO
       DEALLOCATE(tmp,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmplist,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL mp_bcast(ms_vib%b_vec,para_env%source,para_env%group)
    CALL mp_bcast(ms_vib%delta_vec,para_env%source,para_env%group)

    DEALLOCATE(ms_vib%hes_bfgs,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ms_vib%eig_bfgs,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ms_vib%mat_size=0

  END SUBROUTINE bfgs_guess


! *****************************************************************************
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE rest_guess(ms_vib_section,para_env,ms_vib,mass,ionode,particles,nrep,calc_intens,error)

    TYPE(section_vals_type), POINTER         :: ms_vib_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(ms_vib_type)                        :: ms_vib
    REAL(Kind=dp), DIMENSION(:)              :: mass
    LOGICAL                                  :: ionode
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    INTEGER                                  :: nrep
    LOGICAL                                  :: calc_intens
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rest_guess', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: ms_filename
    INTEGER                                  :: hesunit, i, j, mat, natoms, &
                                                ncoord, stat, statint
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ind
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenval
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: approx_H

    natoms=SIZE(particles)
    ncoord=3*natoms
    IF(calc_intens)THEN
       DEALLOCATE(ms_vib%dip_deriv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF(ionode)THEN

       CALL section_vals_val_get(ms_vib_section,"RESTART_FILE_NAME",c_val=ms_filename,&
            error=error)
       IF(ms_filename=="") ms_filename="MS_RESTART"
       CALL open_file(file_name=ms_filename,&
            file_status="UNKNOWN",&
            file_form="UNFORMATTED",&
            file_action="READ",&
            unit_number=hesunit)
       READ(UNIT=hesunit,IOSTAT=stat)mat
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ms_vib%mat_size=mat
    END IF
    CALL mp_bcast(ms_vib%mat_size,para_env%source,para_env%group)
    ALLOCATE(ms_vib%b_mat(ncoord,ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%s_mat(ncoord,ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(calc_intens)THEN
       ALLOCATE(ms_vib%dip_deriv(3,ms_vib%mat_size+nrep),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ionode)THEN
       statint=0
       READ(UNIT=hesunit,IOSTAT=stat)ms_vib%b_mat
       READ(UNIT=hesunit,IOSTAT=stat)ms_vib%s_mat
       IF(calc_intens) READ(UNIT=hesunit,IOSTAT=statint)ms_vib%dip_deriv(:,1:ms_vib%mat_size)
       IF(statint/=0) WRITE(6,FMT="(/,T2,A)")"**  Error while reading MS_RESTART,",&
            "intensities are requested but not present in restart file **"
       CALL close_file(hesunit)
       IF (stat /= 0) THEN
          WRITE(6,FMT="(/,T2,A)")"**  Error while reading MS_RESTART **"
       ELSE
          WRITE(6,FMT="(/,T2,A)")"*** RESTART has been read successfully ***"
       END IF
    END IF
    CALL mp_bcast(ms_vib%b_mat,para_env%source,para_env%group)
    CALL mp_bcast(ms_vib%s_mat,para_env%source,para_env%group)
    IF(calc_intens)CALL mp_bcast( ms_vib%dip_deriv ,para_env%source,para_env%group)
    ALLOCATE(approx_H(ms_vib%mat_size,ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eigenval(ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ind(ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    CALL dgemm('T','N',ms_vib%mat_size,ms_vib%mat_size,SIZE(ms_vib%s_mat,1),1._dp,ms_vib%b_mat,SIZE(ms_vib%b_mat,1),&
         ms_vib%s_mat,SIZE(ms_vib%s_mat,1),0._dp,approx_H,ms_vib%mat_size)
    CALL diamat_all(approx_H,eigenval,error=error)

    CALL select_vector(ms_vib,nrep,mass,ncoord,approx_H,eigenval,ind,ms_vib%b_vec,error=error)
    IF(ms_vib%initial_guess.NE.4)THEN

       ms_vib%b_vec=0._dp
       DO i=1,nrep
          DO j=1,ms_vib%mat_size
             ms_vib% b_vec(:,i)=ms_vib%b_vec(:,i)+approx_H(j,ind(i))*ms_vib%b_mat(:,j)
          END DO
          ms_vib%b_vec(:,i)=ms_vib%b_vec(:,i)/SQRT(DOT_PRODUCT(ms_vib%b_vec(:,i),ms_vib%b_vec(:,i)))
       END DO

       DEALLOCATE(ms_vib%s_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ms_vib%b_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF(calc_intens)THEN
          DEALLOCATE(ms_vib%dip_deriv,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(ms_vib%dip_deriv(3,nrep),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    DEALLOCATE(approx_H,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eigenval,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ind,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nrep
       ms_vib%delta_vec(:,i)= ms_vib%b_vec(:,i)/mass(:)
    END DO

  END SUBROUTINE rest_guess

! *****************************************************************************
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE molden_guess(ms_vib_section,input,para_env,ms_vib,mass,ncoord,nrep,error)
    TYPE(section_vals_type), POINTER         :: ms_vib_section, input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(ms_vib_type)                        :: ms_vib
    REAL(Kind=dp), DIMENSION(:)              :: mass
    INTEGER                                  :: ncoord, nrep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'molden_guess', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: at_name
    CHARACTER(LEN=default_path_length)       :: ms_filename
    CHARACTER(LEN=max_line_length)           :: info
    INTEGER                                  :: i, istat, iw, j, jj, k, &
                                                nvibs, output_molden, stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_val, norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: freq, tmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: modes
    REAL(KIND=dp), DIMENSION(3, ncoord/3)    :: pos
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL section_vals_val_get(ms_vib_section,"RESTART_FILE_NAME",c_val=ms_filename,&
            error=error)
    IF(ms_filename=="") output_molden=cp_print_key_unit_nr(logger,input,&
         "VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB",extension=".mol",file_status='UNKNOWN',&
         file_action="READ",error=error)
    IF(para_env%mepos==para_env%source)THEN

       IF(ms_filename=="") THEN
          iw=output_molden
       ELSE
          CALL open_file(file_name=TRIM(ms_filename),&
               file_status="UNKNOWN",&
               file_form="FORMATTED",&
               file_action="READ",&
               unit_number=iw)
       END IF
       info=""
       READ(iw,*,IOSTAT=stat)info
       READ(iw,*,IOSTAT=stat)info
       istat=0
       nvibs=0
       DO
          READ(iw,*,IOSTAT=stat)info
          istat=istat+stat
          IF(TRIM(ADJUSTL(info))=="[FR-COORD]")EXIT

          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          nvibs=nvibs+1
       END DO
       REWIND(iw)
       istat=0
       READ(iw,*,IOSTAT=stat)info
       istat=istat+stat
       READ(iw,*,IOSTAT=stat)info
       istat=istat+stat

       ALLOCATE(freq(nvibs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(modes(ncoord,nvibs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO i=1,nvibs
          READ(iw,*,IOSTAT=stat)freq(i)
          istat=istat+stat
       END DO
       READ(iw,*)info
       DO i=1,ncoord/3
          READ(iw,*,IOSTAT=stat)at_name,pos(:,i)
          istat=istat+stat
       END DO
       READ(iw,*)info
       DO i=1,nvibs
          READ(iw,*)info
          istat=istat+stat
          DO j=1,ncoord/3
             k=(j-1)*3+1
             READ(iw,*,IOSTAT=stat)modes(k:k+2,i)
             istat=istat+stat
          END DO
       END DO
       IF(ms_filename.NE."")  CALL close_file(iw)
       IF (istat /= 0) THEN
          WRITE(6,FMT="(/,T2,A)")"**  Error while reading MOLDEN file **"
       ELSE
          WRITE(6,FMT="(/,T2,A)")"*** MOLDEN file has been read successfully ***"
       END IF
       !!!!!!!    select modes     !!!!!!
       ALLOCATE(tmp(nvibs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp = 0.0_dp
       ALLOCATE(tmplist(nvibs),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF(ms_vib%select_id==1) my_val=ms_vib%sel_freq
       IF(ms_vib%select_id==2) my_val=(ms_vib%f_range(2)+ms_vib%f_range(1))*0.5_dp
       IF(ms_vib%select_id==1.OR.ms_vib%select_id==2)THEN
          DO i=1,nvibs
             tmp(i)=ABS(my_val-freq(i))
          END DO
       ELSE IF(ms_vib%select_id==3)THEN
          DO i=1,nvibs
             DO j=1,SIZE(ms_vib%inv_atoms)
                DO k=1,3
                   jj=(ms_vib%inv_atoms(j)-1)*3+k
                   tmp(i)=tmp(i)+SQRT(modes(jj,i)**2)
                END DO
             END DO
             IF(freq(i).LE.400._dp)tmp(i)=0._dp
          END DO
          tmp=-tmp
       END IF
       CALL sort(tmp,nvibs,tmplist)
       DO i=1,nrep
          ms_vib%b_vec(:,i)=modes(:,tmplist(i))*mass(:)
          norm=SQRT(DOT_PRODUCT(ms_vib%b_vec(:,i),ms_vib%b_vec(:,i)))
          ms_vib%b_vec(:,i)=ms_vib%b_vec(:,i)/norm
       END DO
       DO i=1,nrep
          ms_vib%delta_vec(:,i)=ms_vib%b_vec(:,i)/mass(:)
       END DO

       DEALLOCATE(freq,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(modes,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmplist,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF
    CALL mp_bcast(ms_vib%b_vec,para_env%source,para_env%group)
    CALL mp_bcast(ms_vib%delta_vec,para_env%source,para_env%group)

    IF(ms_filename=="")CALL cp_print_key_finished_output(output_molden,logger,input,&
         "VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB",error=error)
  END SUBROUTINE molden_guess

! *****************************************************************************
!> \brief Davidson algorithm for to generate a approximate Hessian for mode
!>      selective vibrational analysis
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE  evaluate_H_update_b(rep_env,ms_vib,input,nrep,&
                                  particles,&
                                  mass, &
                                  converged,dx,&
                                  calc_intens,output_unit_ms,output_molden,error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(ms_vib_type)                        :: ms_vib
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(Kind=dp), DIMENSION(:)              :: mass
    LOGICAL                                  :: converged
    REAL(KIND=dp)                            :: dx
    LOGICAL                                  :: calc_intens
    INTEGER                                  :: output_unit_ms, output_molden
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'evaluate_H_update_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, jj, k, natoms, ncoord, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ind
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenval, freq
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: approx_H, H_save, residuum, &
                                                tmp_b, tmp_s
    REAL(KIND=dp), DIMENSION(2, nrep)        :: criteria
    REAL(Kind=dp), DIMENSION(:), POINTER     :: intensities

    natoms=SIZE(particles)
    ncoord=3*natoms
    nrep=SIZE(rep_env%f,2)
    failure=.FALSE.

    !!!!!!!!   reallocate and update the davidson matrices   !!!!!!!!!!
     IF(ms_vib%mat_size.NE.0)THEN

       ALLOCATE(tmp_b(3*natoms,ms_vib%mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_s(3*natoms,ms_vib%mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       tmp_b(:,:)=ms_vib%b_mat
       tmp_s(:,:)=ms_vib%s_mat

       DEALLOCATE(ms_vib%b_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ms_vib%s_mat,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    ALLOCATE(ms_vib%b_mat(3*natoms,ms_vib%mat_size+nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ms_vib%s_mat(3*natoms,ms_vib%mat_size+nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ms_vib%s_mat=0.0_dp

    DO i=1,3*natoms
       IF(ms_vib%mat_size.NE.0)THEN
          DO j=1,ms_vib%mat_size
             ms_vib%b_mat(i,j)=tmp_b(i,j)
             ms_vib%s_mat(i,j)=tmp_s(i,j)
          END DO
       END IF
       DO j=1,nrep
          ms_vib%b_mat(i,ms_vib%mat_size+j)=ms_vib%b_vec(i,j)
       END DO
    END DO

    IF(ms_vib%mat_size.NE.0)THEN
       DEALLOCATE(tmp_s,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ms_vib%mat_size=ms_vib%mat_size+nrep

    ALLOCATE(approx_H(ms_vib%mat_size,ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(H_save(ms_vib%mat_size,ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eigenval(ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !!!!!!!!!!!!  calculate the new derivativ and the approximate hessian

    DO i=1,nrep
       DO j=1,3*natoms
          ms_vib%s_mat(j,ms_vib%mat_size-nrep+i)= -(ms_vib%ms_force(j,i)-rep_env%f(j,i))/(2*ms_vib%step_b(i)*mass(j))
       END DO
    END DO

    CALL dgemm('T','N',ms_vib%mat_size,ms_vib%mat_size,SIZE(ms_vib%s_mat,1),1._dp,ms_vib%b_mat,SIZE(ms_vib%b_mat,1),&
         ms_vib%s_mat,SIZE(ms_vib%s_mat,1),0._dp,approx_H,ms_vib%mat_size)
    H_save(:,:)=approx_H

    CALL diamat_all(approx_H,eigenval,error=error)

    !!!!!!!!!!!! select eigenvalue(s) and vector(s) and calculate the new displacement vector
    ALLOCATE(ind(ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(residuum(SIZE(ms_vib%s_mat,1),nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL select_vector(ms_vib,nrep,mass,ncoord,approx_H,eigenval,ind,residuum,criteria,error)

    DO i=1,nrep
       DO j=1,natoms
          DO k=1,3
             jj=(j-1)*3+k
             ms_vib%delta_vec(jj,i)= ms_vib%b_vec(jj,i)/mass(jj)
          END DO
       END DO
    END DO

    DO i=1,nrep
       ms_vib%step_r(i)=dx/SQRT(DOT_PRODUCT( ms_vib%delta_vec(:,i), ms_vib%delta_vec(:,i)))
       ms_vib%step_b(i)=SQRT(DOT_PRODUCT(ms_vib%step_r(i)*ms_vib%b_vec(:,i),ms_vib%step_r(i)*ms_vib%b_vec(:,i)))
    END DO
    converged=.FALSE.
    IF(MAXVAL(criteria(1,:)).LE.ms_vib%eps(1).AND.MAXVAL(criteria(2,:))&
         .LE.ms_vib%eps(2).OR.ms_vib%mat_size.GE.ncoord)converged=.TRUE.
    ALLOCATE(freq(nrep),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nrep
       freq(i)=SQRT(ABS(eigenval(ind(i)))*massunit)*vibfac
    END DO

    !!!   write information and output   !!!
    IF(converged)THEN
       eigenval(:)=SIGN(1._dp,eigenval(:))*SQRT(ABS(eigenval(:))*massunit)*vibfac
       ALLOCATE(tmp_b(ncoord,ms_vib%mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_b=0._dp
       ALLOCATE(tmp_s(3,ms_vib%mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_s=0._dp
       IF(calc_intens)THEN
          ALLOCATE(intensities(ms_vib%mat_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          intensities=0._dp
       END IF
       DO i=1,ms_vib%mat_size
          DO j=1,ms_vib%mat_size
             tmp_b(:,i)=tmp_b(:,i)+approx_H(j,i)*ms_vib%b_mat(:,j)/mass(:)
          END DO
          tmp_b(:,i)=tmp_b(:,i)/SQRT(DOT_PRODUCT(tmp_b(:,i),tmp_b(:,i)))
       END DO
       IF(calc_intens)THEN
          DO i=1,ms_vib%mat_size
             DO j=1,ms_vib%mat_size
                tmp_s(:,i)=tmp_s(:,i)+ms_vib%dip_deriv(:,j)*approx_H(j,i)
             END DO
             IF(calc_intens) intensities(i)=SQRT(DOT_PRODUCT(tmp_s(:,i),tmp_s(:,i)))
          END DO
       END IF
       IF(calc_intens)THEN
          CALL ms_out(output_unit_ms, converged, freq, criteria, ms_vib, &
               input, nrep,approx_H, eigenval, calc_intens, &
               intensities=intensities, error=error)
       ELSE
          CALL ms_out(output_unit_ms, converged, freq, criteria, ms_vib, &
               input, nrep, approx_H, eigenval, calc_intens, error=error)
       ENDIF
       CALL molden_out(input,particles,eigenval,tmp_b,intensities,calc_intens,&
                       dump_only_positive=.TRUE.,error=error)
       IF (calc_intens) THEN
          DEALLOCATE(intensities,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(tmp_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(tmp_s,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    IF(.NOT.converged) CALL ms_out(output_unit_ms, converged, freq, criteria, &
        ms_vib, input, nrep,approx_H, eigenval, calc_intens, error=error)

    DEALLOCATE(freq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(approx_H,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eigenval,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(residuum,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ind,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


  END SUBROUTINE evaluate_H_update_b

! *****************************************************************************
!> \brief writes the output for a mode tracking calculation
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE select_vector(ms_vib,nrep,mass,ncoord,approx_H,eigenval,ind,residuum,criteria,error)

    TYPE(ms_vib_type)                        :: ms_vib
    INTEGER                                  :: nrep
    REAL(Kind=dp), DIMENSION(:)              :: mass
    INTEGER                                  :: ncoord
    REAL(KIND=dp), DIMENSION(:, :)           :: approx_H
    REAL(Kind=dp), DIMENSION(:)              :: eigenval
    INTEGER, DIMENSION(:)                    :: ind
    REAL(KIND=dp), DIMENSION(:, :)           :: residuum
    REAL(KIND=dp), DIMENSION(2, nrep), &
      OPTIONAL                               :: criteria
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'select_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, jj, k, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_val, norm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_b

    ALLOCATE(tmp(ms_vib%mat_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    SELECT CASE(ms_vib%select_id)
    CASE(1)
       my_val=(ms_vib%sel_freq/(vibfac))**2/massunit
       DO i=1,ms_vib%mat_size
          tmp(i)=ABS(my_val-eigenval(i))
       END DO
       CALL sort(tmp,(ms_vib%mat_size),ind)
       residuum=0._dp
       DO j=1,nrep
          DO i=1,ms_vib%mat_size
             residuum(:,j)=residuum(:,j)+approx_H(i,ind(j))*(ms_vib%s_mat(:,i)-eigenval(ind(j))*ms_vib%b_mat(:,i))
          END DO
       END DO
    CASE(2)
       CALL get_vibs_in_range(ms_vib,approx_H,eigenval,residuum,nrep,ind)
    CASE(3)

       ALLOCATE(tmp_b(ncoord,ms_vib%mat_size),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp_b=0._dp


       DO i=1,ms_vib%mat_size
          DO j=1,ms_vib%mat_size
             tmp_b(:,i)=tmp_b(:,i)+approx_H(j,i)*ms_vib%b_mat(:,j)/mass(:)
          END DO
          tmp_b(:,i)=tmp_b(:,i)/SQRT(DOT_PRODUCT(tmp_b(:,i),tmp_b(:,i)))
       END DO
       tmp=0._dp
       DO i=1,ms_vib%mat_size
          DO j=1,SIZE(ms_vib%inv_atoms)
             DO k=1,3
                jj=(ms_vib%inv_atoms(j)-1)*3+k
                tmp(i)=tmp(i)+SQRT(tmp_b(jj,i)**2)
             END DO
          END DO
          IF(.NOT.ASSOCIATED(ms_vib%inv_range))THEN
             IF((SIGN(1._dp,eigenval(i))*SQRT(ABS(eigenval(i))*massunit)*vibfac).LE.400._dp)tmp(i)=0._dp
          ELSE
             IF((SIGN(1._dp,eigenval(i))*SQRT(ABS(eigenval(i))*massunit)*vibfac).LE.ms_vib%inv_range(1))tmp(i)=0._dp
             IF((SIGN(1._dp,eigenval(i))*SQRT(ABS(eigenval(i))*massunit)*vibfac).GE.ms_vib%inv_range(2))tmp(i)=0._dp
          END IF
       END DO
       tmp=-tmp
       CALL sort(tmp,(ms_vib%mat_size),ind)
       residuum=0._dp

       DO j=1,nrep
          DO i=1,ms_vib%mat_size
             residuum(:,j)=residuum(:,j)+approx_H(i,ind(j))*(ms_vib%s_mat(:,i)-eigenval(ind(j))*ms_vib%b_mat(:,i))
          END DO
       END DO
       DEALLOCATE(tmp_b,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

    DO j=1,nrep
       DO i=1,ms_vib%mat_size
          residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),ms_vib%b_mat(:,i))* ms_vib%b_mat(:,i)
       END DO
    END DO
    IF(PRESENT(criteria))THEN
       DO i=1,nrep
          criteria(1,i)=MAXVAL((residuum(:,i)))
          criteria(2,i)=SQRT(DOT_PRODUCT(residuum(:,i),residuum(:,i)))
       END DO
    END IF


    DO i=1,nrep
       norm=SQRT(DOT_PRODUCT(residuum(:,i),residuum(:,i)))
       residuum(:,i)=residuum(:,i)/norm
    END DO


    DO k=1,10
       DO j=1,nrep
             DO i=1,ms_vib%mat_size
                residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),ms_vib%b_mat(:,i))* ms_vib%b_mat(:,i)
                residuum(:,j)=residuum(:,j)/SQRT(DOT_PRODUCT(residuum(:,j),residuum(:,j)))
             END DO
             IF(nrep.GT.1)THEN
                DO i=1,nrep
                   IF(i.NE.j)THEN
                      residuum(:,j)=residuum(:,j)-DOT_PRODUCT(residuum(:,j),residuum(:,i))* residuum(:,i)
                      residuum(:,j)=residuum(:,j)/SQRT(DOT_PRODUCT(residuum(:,j),residuum(:,j)))
                   END IF
                END DO
             END IF
       END DO
    END DO
    ms_vib%b_vec=residuum
    DEALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE select_vector

! *****************************************************************************
!> \brief writes the output for a mode tracking calculation
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE ms_out(iw, converged, freq, criter, ms_vib, input, nrep, &
    approx_H, eigenval, calc_intens, intensities, error)

    INTEGER                                  :: iw
    LOGICAL                                  :: converged
    REAL(KIND=dp), DIMENSION(:)              :: freq
    REAL(KIND=dp), DIMENSION(:, :)           :: criter
    TYPE(ms_vib_type)                        :: ms_vib
    TYPE(section_vals_type), POINTER         :: input
    INTEGER                                  :: nrep
    REAL(KIND=dp), DIMENSION(:, :)           :: approx_H
    REAL(KIND=dp), DIMENSION(:)              :: eigenval
    LOGICAL                                  :: calc_intens
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: intensities
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ms_out', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, msunit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: crit_a, crit_b
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: residuum
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: ms_vib_section

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ms_vib_section => section_vals_get_subs_vals(input, &
        "VIBRATIONAL_ANALYSIS%MODE_SELECTIVE", error=error)

    IF(converged)THEN
       IF(iw.GT.0)THEN
          WRITE(iw,'(T2,A)')"MS| DAVIDSON ALGORITHM CONVERGED"
          DO i=1,nrep
             WRITE(iw,'(T2,"MS| TRACKED FREQUENCY (",I0,") IS:",F12.6,3X,A)') i,freq(i),'cm-1'
          END DO
          ALLOCATE(residuum(SIZE(ms_vib%b_mat,1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          WRITE ( iw, '( /, 1X, 79("-") )' )
          WRITE ( iw, '( 25X, A)' ) 'FREQUENCY AND CONVERGENCE LIST'
          IF(PRESENT(intensities))THEN
             WRITE ( iw, '(3X,5(4X, A))')'FREQUENCY','INTENSITY','MAXVAL CRITERIA','NORM CRITERIA', 'CONVERGENCE'
          ELSE
             WRITE ( iw, '(3X,5(4X, A))')'FREQUENCY','MAXVAL CRITERIA','NORM CRITERIA', 'CONVERGENCE'
          ENDIF
          DO i=1,SIZE(ms_vib%b_mat,2)
             residuum=0._dp
             DO j=1,SIZE(ms_vib%b_mat,2)
                residuum(:)=residuum(:)+approx_H(j,i)*(ms_vib%s_mat(:,j)-eigenval(i)*ms_vib%b_mat(:,j))
             END DO
             DO j=1,ms_vib%mat_size
                   residuum(:)=residuum(:)-DOT_PRODUCT(residuum(:),ms_vib%b_mat(:,j))* ms_vib%b_mat(:,j)
             END DO
             crit_a=MAXVAL(residuum(:))
             crit_b=SQRT(DOT_PRODUCT(residuum,residuum))
             IF(PRESENT(intensities))THEN
                IF(crit_a.LE.ms_vib%eps(1).AND.crit_b.LE.ms_vib%eps(2))THEN
                   IF(eigenval(i).GT.0._dp) WRITE(iw,'(2X,A,2X,F8.3,1X,F12.6,3X,E12.3,7X,E12.3,11X,A)')&
                       'VIB|',eigenval(i),intensities(i),crit_a,crit_b,'YES'
                ELSE
                   IF(eigenval(i).GT.0._dp) WRITE(iw,'(2X,A,2X,F8.3,1X,F12.6,3X,E12.3,7X,E12.3,11X,A)')&
                       'VIB|',eigenval(i),intensities(i),crit_a,crit_b,'NO'
                ENDIF
             ELSE
                IF(crit_a.LE.ms_vib%eps(1).AND.crit_b.LE.ms_vib%eps(2))THEN
                   IF(eigenval(i).GT.0._dp) WRITE(iw,'(2X,A,2X,F8.3,5X,E12.6,5X,E12.3,11X,A)')&
                       'VIB|',eigenval(i),crit_a,crit_b,'YES'
                ELSE
                   IF(eigenval(i).GT.0._dp) WRITE(iw,'(2X,A,2X,F8.3,5X,E12.6,5X,E12.3,11X,A)')&
                       'VIB|',eigenval(i),crit_a,crit_b,'NO'
                ENDIF
             END IF
          END DO
          DEALLOCATE(residuum,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          msunit = cp_print_key_unit_nr(logger, ms_vib_section, &
               "PRINT%MS_RESTART", extension=".bin", middle_name="MS_RESTART", &
               file_status="REPLACE", file_form="UNFORMATTED", &
               file_action="WRITE", error=error)

          IF (msunit>0) THEN
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%mat_size
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%b_mat
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%s_mat
             IF(calc_intens)WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%dip_deriv
          END IF

          CALL cp_print_key_finished_output(msunit, logger, ms_vib_section, &
               "PRINT%MS_RESTART", error=error)
       END IF
    ELSE
       IF(iw.GT.0)THEN
          msunit = cp_print_key_unit_nr(logger, ms_vib_section, &
               "PRINT%MS_RESTART", extension=".bin", middle_name="MS_RESTART", &
               file_status="REPLACE", file_form="UNFORMATTED", &
               file_action="WRITE", error=error)

          IF (msunit>0) THEN
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%mat_size
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%b_mat
             WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%s_mat
             IF(calc_intens)WRITE(UNIT=msunit,IOSTAT=stat)ms_vib%dip_deriv
          END IF

          CALL cp_print_key_finished_output(msunit, logger, ms_vib_section, &
               "PRINT%MS_RESTART", error=error)

          WRITE(iw,'(T2,A,3X,I6)')"MS| ITERATION STEP", ms_vib%mat_size/nrep
          DO i=1,nrep
             IF(criter(1,i).LE.1E-7.AND.(criter(2,i)).LE.1E-6)THEN
                WRITE(iw,'(T2,A,3X,F12.6,A)')"MS| TRACKED MODE ", freq(i) , "cm-1  IS  CONVERGED"
             ELSE
                WRITE(iw,'(T2,A,3X,F12.6,A)')"MS| TRACKED MODE ", freq(i) , "cm-1  NOT  CONVERGED"
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE ms_out

! *****************************************************************************
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE get_vibs_in_range(ms_vib,approx_H,eigenval,residuum,nrep,ind)

    TYPE(ms_vib_type)                        :: ms_vib
    REAL(KIND=dp), DIMENSION(:, :)           :: approx_H
    REAL(KIND=dp), DIMENSION(:)              :: eigenval
    REAL(KIND=dp), DIMENSION(:, :)           :: residuum
    INTEGER                                  :: nrep
    INTEGER, DIMENSION(:)                    :: ind

    INTEGER                                  :: count1, count2, i, j
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map2
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: map1
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp, tmp1
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_resid
    REAL(KIND=dp), DIMENSION(2)              :: myrange

    myrange(:)=(ms_vib%f_range(:)/(vibfac))**2/massunit
    count1=0
    count2=0
    residuum=0.0_dp
    ms_vib%mat_size=SIZE(ms_vib%b_mat,2)
    ALLOCATE(map1(SIZE(eigenval),2))
    ALLOCATE(tmp(SIZE(eigenval)))
    DO i=1,SIZE(eigenval)
       IF(ABS(eigenval(i)-myrange(1))+ABS(eigenval(i)-myrange(2)).LE.&
            ABS(myrange(1)-myrange(2))+myrange(1)*0.001_dp)THEN
          count1=count1+1
          map1(count1,1)=i
       ELSE
          count2=count2+1
          map1(count2,2)=i
          tmp(count2)=MIN(ABS(eigenval(i)-myrange(1)),ABS(eigenval(i)-myrange(2)))
       END IF
    END DO

    IF(count1.EQ.nrep)THEN
       DO j=1,count1
          DO i=1,ms_vib%mat_size
             residuum(:,j)=residuum(:,j)+approx_H(i,map1(j,1))*(ms_vib%s_mat(:,i)-eigenval(map1(j,1))*ms_vib%b_mat(:,i))
             ind(j)=map1(j,1)
          END DO
       END DO
    ELSE IF(count1.gt.nrep)THEN
       ALLOCATE(tmp_resid(SIZE(ms_vib%b_mat,1),count1))
       ALLOCATE(tmp1(count1))
       ALLOCATE(map2(count1))
       tmp_resid=0._dp
       DO j=1,count1
          DO i=1,ms_vib%mat_size
             tmp_resid(:,j)=tmp_resid(:,j)+approx_H(i,map1(j,1))*(ms_vib%s_mat(:,i)-eigenval(map1(j,1))*ms_vib%b_mat(:,i))
          END DO
       END DO

       DO j=1,count1
          DO i=1,ms_vib%mat_size
             tmp_resid(:,j)=tmp_resid(:,j)-DOT_PRODUCT(tmp_resid(:,j),ms_vib%b_mat(:,i))* ms_vib%b_mat(:,i)
          END DO
          tmp(j)=MAXVAL(tmp_resid(:,j))
       END DO
       CALL sort(tmp,count1,map2)
       DO j=1,nrep
          residuum(:,j)=tmp_resid(:,map2(count1+1-j))
          ind(j)=map1(map2(count1+1-j),1)
       END DO
       DEALLOCATE(tmp_resid)
       DEALLOCATE(tmp1)
       DEALLOCATE(map2)
    ELSE IF(count1.LT.nrep)THEN

       ALLOCATE(map2(count2))
       IF(count1.NE.0)THEN
          DO j=1,count1
             DO i=1,ms_vib%mat_size
                residuum(:,j)=residuum(:,j)+approx_H(i,map1(j,1))*(ms_vib%s_mat(:,i)-eigenval(map1(j,1))*ms_vib%b_mat(:,i))
             END DO
             ind(j)=map1(j,1)
          END DO
       END IF
       CALL sort(tmp,count2,map2)
       DO j=1,nrep-count1
          DO i=1,ms_vib%mat_size
             residuum(:,count1+j)=residuum(:,count1+j)+approx_H(i,map1(map2(j),2))&
                  *(ms_vib%s_mat(:,i)-eigenval(map1(map2(j),2))*ms_vib%b_mat(:,i))
          END DO
          ind(count1+j)=map1(map2(j),2)
       END DO

       DEALLOCATE(map2)
    END IF

    DEALLOCATE(map1)
    DEALLOCATE(tmp)

  END SUBROUTINE get_vibs_in_range
END MODULE mode_selective
