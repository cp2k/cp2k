!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module to perform a counterpoise correction (BSSE)
!> \par History
!>      6.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE bsse
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp2k_info,                       ONLY: write_restart_header
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: create_small_subsys
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_qs
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,&
                                             section_vals_write
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE string_utilities,                ONLY: compress
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bsse'

  PUBLIC :: do_bsse_calculation

CONTAINS

! *****************************************************************************
!> \brief Perform an COUNTERPOISE CORRECTION (BSSE)
!>      For a 2-body system the correction scheme can be represented as:
!>
!>      E_{AB}^{2}        = E_{AB}(AB) - E_A(AB) - E_B(AB)  [BSSE-corrected interaction energy]
!>      E_{AB}^{2,uncorr} = E_{AB}(AB) - E_A(A)  - E_B(B)
!>      E_{AB}^{CP}       = E_{AB}(AB) + [ E_A(A) - E_A(AB) ] + [ E_B(B) - E_B(AB) ]
!>                                                          [CP-corrected total energy of AB]
!> \param force_env ...
!> \param globenv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE do_bsse_calculation(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_bsse_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istart, k, num_of_conf, &
                                                Num_of_Frag, stat
    INTEGER, DIMENSION(:, :), POINTER        :: conf
    LOGICAL                                  :: explicit, failure, should_stop
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_env_type), POINTER            :: force_env_loc
    TYPE(section_vals_type), POINTER         :: bsse_section, &
                                                fragment_energies_section, &
                                                n_frags, root_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(bsse_section, n_frags, Em, conf, force_env_loc)
       logger => cp_error_get_logger(error)
       root_section => force_env%root_section
       bsse_section => section_vals_get_subs_vals(force_env%force_env_section,"BSSE",error=error)
       n_frags => section_vals_get_subs_vals(bsse_section,"FRAGMENT",error=error)
       CALL section_vals_get(n_frags,n_repetition=Num_of_Frag,error=error)

       ! Number of configurations
       num_of_conf = 0
       DO k = 1, Num_of_frag
          num_of_conf = num_of_conf + FACT(Num_of_frag) / (FACT(k)*FACT(Num_of_frag-k))
       END DO
       ALLOCATE(conf(num_of_conf,Num_of_frag),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Em(num_of_conf), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL gen_Nbody_conf(Num_of_frag, num_of_conf, conf)

       should_stop = .FALSE.
       istart      = 0
       fragment_energies_section => section_vals_get_subs_vals(bsse_section,"FRAGMENT_ENERGIES",error=error)
       CALL section_vals_get(fragment_energies_section, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(fragment_energies_section,"_DEFAULT_KEYWORD_",n_rep_val=istart,error=error)
          DO i = 1, istart
             CALL section_vals_val_get(fragment_energies_section,"_DEFAULT_KEYWORD_",r_val=Em(i),&
                  i_rep_val=i, error=error)
          END DO
       END IF
       ! Setup the iteration level for BSSE
       CALL cp_add_iter_level(logger%iter_info,"BSSE",error=error)
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=istart,error=error)

       ! Evaluating the energy of the N-body cluster terms
       DO i = istart+1, num_of_conf
          CALL cp_iterate(logger%iter_info,last=(i==num_of_conf),iter_nr=i,error=error)
          CALL eval_bsse_energy(conf(i,:),Em(i),force_env,force_env_loc, n_frags, &
               root_section,globenv, should_stop, error=error)
          IF (should_stop) EXIT

          ! If no signal was received in the inner loop let's check also at this stage
          CALL external_control(should_stop,"BSSE",globenv=globenv,error=error)
          IF (should_stop) EXIT

          ! Dump Restart info only if the calculation of the energy of a configuration
          ! ended nicely..
          CALL section_vals_val_set(fragment_energies_section,"_DEFAULT_KEYWORD_",r_val=Em(i),&
               i_rep_val=i,error=error)
          CALL write_bsse_restart(bsse_section, root_section, error=error)
       END DO
       IF (.NOT.should_stop) CALL dump_bsse_results(conf, Em, num_of_frag,  bsse_section, error)
       CALL cp_rm_iter_level(logger%iter_info,"BSSE",error=error)
       DEALLOCATE(Em,   stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(conf, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE do_bsse_calculation

! *****************************************************************************
!> \brief Evaluate the N-body energy contribution to the BSSE evaluation
!> \param conf ...
!> \param Em ...
!> \param force_env ...
!> \param force_env_loc ...
!> \param n_frags ...
!> \param root_section ...
!> \param globenv ...
!> \param should_stop ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE eval_bsse_energy(conf, Em, force_env, force_env_loc, n_frags, root_section,&
       globenv, should_stop, error)
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf
    REAL(KIND=dp), INTENT(OUT)               :: Em
    TYPE(force_env_type), POINTER            :: force_env, force_env_loc
    TYPE(section_vals_type), POINTER         :: n_frags, root_section
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(OUT)                     :: should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_bsse_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, Num_of_sub_conf, &
                                                Num_of_sub_frag, stat
    INTEGER, DIMENSION(:, :), POINTER        :: conf_loc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em_loc

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(conf_loc, Em_loc)
       should_stop = .FALSE.
       ! Count the number of subconfiguration to evaluate..
       Num_of_sub_frag = COUNT(conf == 1)
       Num_of_sub_conf = 0
       IF (Num_of_sub_frag == 1) THEN
          CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf, n_frags, root_section, globenv, Em, error)
       ELSE
          my_energy = 0.0_dp
          DO k = 1, Num_of_sub_frag
             Num_of_sub_conf = Num_of_sub_conf + &
                  FACT(Num_of_sub_frag) / (FACT(k)*FACT(Num_of_sub_frag-k))
          END DO
          ALLOCATE(conf_loc(Num_of_sub_conf,Num_of_sub_frag), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(Em_loc(Num_of_sub_conf), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Em_loc = 0.0_dp
          CALL gen_Nbody_conf(Num_of_sub_frag, Num_of_sub_conf, conf_loc)
          CALL make_plan_conf(conf, conf_loc)
          DO i = 1, Num_of_sub_conf
             CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf_loc(i,:), n_frags,&
                  root_section,globenv, Em_loc(i), error)
             CALL external_control(should_stop,"BSSE",globenv=globenv,error=error)
             IF (should_stop) EXIT
          END DO
          ! Energy
          k = COUNT(conf == 1)
          DO i =  1, Num_of_sub_conf
             j = COUNT(conf_loc(i,:) == 1)
             my_energy = my_energy + (-1.0_dp)**(k+j) * Em_loc(i)
          END DO
          Em = my_energy
          DEALLOCATE(Em_loc, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(conf_loc, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE eval_bsse_energy

! *****************************************************************************
!> \brief Evaluate the N-body energy contribution to the BSSE evaluation
!> \param force_env ...
!> \param force_env_loc ...
!> \param conf ...
!> \param conf_loc ...
!> \param n_frags ...
!> \param root_section ...
!> \param globenv ...
!> \param energy ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2005 created [tlaino]
!>      2014/09/17 made atom list to be read from repeated occurance of LIST [LTong]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE eval_bsse_energy_low(force_env, force_env_loc, conf, conf_loc, n_frags, &
       root_section, globenv, energy, error)
    TYPE(force_env_type), POINTER            :: force_env, force_env_loc
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(section_vals_type), POINTER         :: n_frags, root_section
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_bsse_energy_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_type
    INTEGER :: i, ir, isize, j, k, method_name_id, my_targ, n_rep, &
      num_of_frag, old_size, present_charge, present_multpl, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_index, atom_list, &
                                                my_conf, tmplist
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys, subsys_loc
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: bsse_section, dft_section, &
                                                force_env_section, &
                                                subsys_section

    failure = .FALSE.
    CALL section_vals_get(n_frags,n_repetition=num_of_frag,error=error)
    CPPostcondition(SIZE(conf)==num_of_frag,cp_failure_level,routineP,error,failure)
    NULLIFY(subsys_loc, subsys, particles, para_env, cell, atom_index, atom_type, tmplist,&
         force_env_section )
    CALL force_env_get(force_env, force_env_section=force_env_section, error=error)
    CALL section_vals_val_get(force_env_section,"METHOD",i_val=method_name_id,error=error)
    bsse_section   => section_vals_get_subs_vals(force_env_section,"BSSE",error=error)
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    dft_section    => section_vals_get_subs_vals(force_env_section,"DFT",error=error)
    IF (.NOT.failure) THEN
       ALLOCATE(my_conf(SIZE(conf)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_conf = conf
       CALL force_env_get(force_env=force_env, subsys=subsys, para_env=para_env,&
            cell=cell, error=error)
       CALL cp_subsys_get(subsys, particles=particles, error=error)
       isize = 0
       ALLOCATE(atom_index(isize),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, num_of_frag
          IF (conf(i) == 1) THEN
             !
             ! Get the list of atoms creating the present fragment
             !
             old_size=isize
             CALL section_vals_val_get(n_frags, "LIST", i_rep_section=i, n_rep_val=n_rep, error=error)
             IF (n_rep /= 0) THEN
                DO ir = 1, n_rep
                   CALL section_vals_val_get(n_frags, "LIST", i_rep_section=i, i_rep_val=ir, i_vals=tmplist, error=error)
                   CALL reallocate(atom_index,1,isize+SIZE(tmplist))
                   atom_index(isize+1:isize+SIZE(tmplist)) = tmplist
                   isize = SIZE(atom_index)
                END DO
             END IF
             my_conf(i) = isize - old_size
             CPPostcondition(conf(i)/=0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
       CALL conf_info_setup(present_charge, present_multpl, conf, conf_loc, bsse_section,&
            dft_section, error)
       !
       ! Get names and modify the ghost ones
       !
       ALLOCATE(atom_type(isize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO j = 1, isize
          my_targ = atom_index(j)
          DO k = 1, SIZE(particles%els)
             CALL get_atomic_kind(particles%els(k)%atomic_kind, atom_list=atom_list, name=name)
             IF (ANY(atom_list == my_targ)) EXIT
          END DO
          atom_type(j) = name
       END DO
       DO i = 1, SIZE(conf_loc)
          IF (my_conf(i) /= 0 .AND. conf_loc(i) == 0) THEN
             DO j = SUM(my_conf(1:i-1))+1, SUM(my_conf(1:i))
                atom_type(j) = TRIM(atom_type(j))//"_ghost"
             END DO
          END IF
       END DO
       CALL dump_bsse_info(atom_index, atom_type, conf, conf_loc, bsse_section,&
                           present_charge, present_multpl,  error)
       !
       ! Let's start setting up environments and calculations
       !
       energy = 0.0_dp
       IF (method_name_id == do_qs) THEN
          CALL create_small_subsys(subsys_loc, big_subsys=subsys,&
               small_para_env=para_env, small_cell=cell,sub_atom_index=atom_index,&
               sub_atom_kind_name=atom_type,para_env=para_env,&
               force_env_section=force_env_section,subsys_section=subsys_section,error=error)

          CALL quickstep_create_force_env(force_env_loc,root_section,para_env,globenv,subsys=subsys_loc,&
               force_env_section=force_env_section, subsys_section=subsys_section, use_motion_section=.FALSE.,&
               error=error)

          CALL cp_subsys_release(subsys_loc,error=error)
          !
          ! Evaluate Energy
          !
          CALL force_env_calc_energy_force(force_env_loc,calc_force=.FALSE.,error=error)
          CALL force_env_get(force_env_loc,potential_energy=energy,error=error)
          CALL force_env_release(force_env_loc, error=error)
          NULLIFY(force_env_loc)
       ELSE
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(atom_index, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(atom_type, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_conf, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE eval_bsse_energy_low

! *****************************************************************************
!> \brief Dumps bsse information (configuration fragment)
!> \param atom_index ...
!> \param atom_type ...
!> \param conf ...
!> \param conf_loc ...
!> \param bsse_section ...
!> \param present_charge ...
!> \param present_multpl ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  dump_bsse_info(atom_index, atom_type, conf, conf_loc, bsse_section, &
                             present_charge, present_multpl, error)
    INTEGER, DIMENSION(:), POINTER           :: atom_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_type
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(section_vals_type), POINTER         :: bsse_section
    INTEGER, INTENT(IN)                      :: present_charge, present_multpl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_bsse_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat, iw
    CHARACTER(len=20*SIZE(conf))             :: conf_loc_s, conf_s
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    IF (.NOT.failure) THEN
        logger => cp_error_get_logger(error)
        iw=cp_print_key_unit_nr(logger,bsse_section,"PRINT%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)
        IF (iw>0) THEN
          WRITE (conf_s,fmt="(1000I0)",iostat=istat)conf ;
          IF(istat.NE.0) conf_s="exceeded"
          CALL compress(conf_s,full=.TRUE.)
          WRITE (conf_loc_s,fmt="(1000I0)",iostat=istat)conf_loc ;
          IF(istat.NE.0) conf_loc_s="exceeded"
          CALL compress(conf_loc_s,full=.TRUE.)

          WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T5,A,T30,A,T55,A,T80,A)")&
               "-","BSSE CALCULATION","FRAGMENT CONF: "//TRIM(conf_s),"FRAGMENT SUBCONF: "//TRIM(conf_loc_s),"-"
          WRITE (UNIT=iw,FMT="(T2,A,T30,A,I6,T55,A,I6,T80,A)") "-","CHARGE =",present_charge,"MULTIPLICITY =",&
                                                         present_multpl,"-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,T60,A,T80,A)") "-","ATOM INDEX","ATOM NAME","-"
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,T60,A,T80,A)") "-","----------","---------","-"
          DO i = 1, SIZE(atom_index)
             WRITE (UNIT=iw,FMT="(T2,A,T20,I6,T61,A,T80,A)") "-",atom_index(i),TRIM(atom_type(i)),"-"
          END DO
          WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       END IF

       CALL cp_print_key_finished_output(iw,logger,bsse_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)

    END IF
  END SUBROUTINE dump_bsse_info

! *****************************************************************************
!> \brief Read modified parameters for configurations
!> \param present_charge ...
!> \param present_multpl ...
!> \param conf ...
!> \param conf_loc ...
!> \param bsse_section ...
!> \param dft_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE  conf_info_setup(present_charge, present_multpl, conf, conf_loc,&
       bsse_section, dft_section, error)
    INTEGER, INTENT(OUT)                     :: present_charge, present_multpl
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(section_vals_type), POINTER         :: bsse_section, dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'conf_info_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nconf
    INTEGER, DIMENSION(:), POINTER           :: glb_conf, sub_conf
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: configurations

    present_charge = 0
    present_multpl = 0
    failure = .FALSE.
    NULLIFY(configurations, glb_conf, sub_conf)
    IF (.NOT.failure) THEN
       ! Loop over all configurations to pick up the right one
       configurations => section_vals_get_subs_vals(bsse_section,"CONFIGURATION",error=error)
       CALL section_vals_get(configurations,explicit=explicit,n_repetition=nconf,error=error)
       IF (explicit) THEN
          DO i=1, nconf
             CALL section_vals_val_get(configurations,"GLB_CONF", i_rep_section=i, i_vals=glb_conf, error=error)
             CALL cp_assert(SIZE(glb_conf)==SIZE(conf),cp_fatal_level,cp_assertion_failed,routineP,&
                  "GLB_CONF requires a binary description of the configuration. Number of integer "//&
                  "different from the number of fragments defined!"//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             CALL section_vals_val_get(configurations,"SUB_CONF", i_rep_section=i, i_vals=sub_conf, error=error)
             CALL cp_assert(SIZE(sub_conf)==SIZE(conf),cp_fatal_level,cp_assertion_failed,routineP,&
                  "SUB_CONF requires a binary description of the configuration. Number of integer "//&
                  "different from the number of fragments defined!"//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             IF (ALL(conf==glb_conf).AND.ALL(conf_loc==sub_conf)) THEN
                CALL section_vals_val_get(configurations,"CHARGE", i_rep_section=i,&
                                          i_val=present_charge, error=error)
                CALL section_vals_val_get(configurations,"MULTIPLICITY", i_rep_section=i,&
                                          i_val=present_multpl, error=error)
             END IF
          END DO
       END IF
       ! Setup parameter for this configuration
       CALL section_vals_val_set(dft_section,"CHARGE",i_val=present_charge, error=error)
       CALL section_vals_val_set(dft_section,"MULTIPLICITY",i_val=present_multpl, error=error)
    END IF
  END SUBROUTINE conf_info_setup

! *****************************************************************************
!> \brief Dumps results
!> \param conf ...
!> \param Em ...
!> \param num_of_frag ...
!> \param bsse_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE  dump_bsse_results(conf, Em,  num_of_frag, bsse_section, error)
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: conf
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em
    INTEGER, INTENT(IN)                      :: num_of_frag
    TYPE(section_vals_type), POINTER         :: bsse_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_bsse_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iw
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    IF (.NOT.failure) THEN
        logger => cp_error_get_logger(error)
        iw=cp_print_key_unit_nr(logger,bsse_section,"PRINT%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

        IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T36,A,T80,A)")&
               "-","BSSE RESULTS","-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,F16.6,T80,A)") "-","CP-corrected Total energy:",SUM(Em),"-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          DO i = 1, SIZE(conf,1)
             IF (i.GT.1) THEN
                IF (SUM(conf(i-1,:))==1.AND.SUM(conf(i,:))/=1) THEN
                   WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
                END IF
             END IF
             WRITE (UNIT=iw,FMT="(T2,A,T24,I3,A,F16.6,T80,A)") "-",SUM(conf(i,:)),"-body contribution:",Em(i),"-"
          END DO
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,F16.6,T80,A)") "-","BSSE-free interaction energy:",SUM(Em(Num_of_frag+1:)),"-"
          WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
        END IF

        CALL cp_print_key_finished_output(iw,logger,bsse_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)

    END IF
  END SUBROUTINE dump_bsse_results

! *****************************************************************************
!> \brief generate the N-body configuration for the N-body BSSE evaluation
!> \param Num_of_frag ...
!> \param num_of_conf ...
!> \param conf ...
!> \par History
!>      07.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE gen_Nbody_conf(Num_of_frag, num_of_conf, conf)
    INTEGER, INTENT(IN)                      :: Num_of_frag, num_of_conf
    INTEGER, DIMENSION(:, :), POINTER        :: conf

    INTEGER                                  :: k, my_ind

    my_ind      = 0
    !
    ! Set up the N-body configurations
    !
    conf = 0
    DO k= 1, Num_of_frag
       CALL build_Nbody_conf( 1, Num_of_frag, conf, k, my_ind)
    END DO
  END SUBROUTINE gen_Nbody_conf

! *****************************************************************************
!> \brief ...
!> \param ldown ...
!> \param lup ...
!> \param conf ...
!> \param k ...
!> \param my_ind ...
! *****************************************************************************
  RECURSIVE SUBROUTINE build_Nbody_conf(ldown, lup, conf, k, my_ind)
    INTEGER, INTENT(IN)                      :: ldown, lup
    INTEGER, DIMENSION(:, :), POINTER        :: conf
    INTEGER, INTENT(IN)                      :: k
    INTEGER, INTENT(INOUT)                   :: my_ind

    INTEGER                                  :: i, kloc, my_ind0

    kloc = k - 1
    my_ind0 = my_ind
    IF (kloc /= 0) THEN
       DO i = ldown, lup
          CALL build_Nbody_conf( i+1, lup, conf, kloc, my_ind)
          conf(my_ind0+1:my_ind,i) = 1
          my_ind0 = my_ind
       END DO
    ELSE
       DO i = ldown, lup
          my_ind = my_ind + 1
          conf(my_ind,i) = 1
       END DO
    END IF
  END SUBROUTINE build_Nbody_conf

! *****************************************************************************
!> \brief ...
!> \param num ...
!> \retval my_fact ...
! *****************************************************************************
  RECURSIVE FUNCTION FACT(num) RESULT(my_fact)
    INTEGER, INTENT(IN)                      :: num
    INTEGER                                  :: my_fact

    IF (num<=1) THEN
       my_fact = 1
    ELSE
       my_fact = num * FACT(num-1)
    END IF
  END FUNCTION FACT

! *****************************************************************************
!> \brief ...
!> \param main_conf ...
!> \param conf ...
! *****************************************************************************
  SUBROUTINE make_plan_conf (main_conf, conf)
    INTEGER, DIMENSION(:), INTENT(IN)        :: main_conf
    INTEGER, DIMENSION(:, :), POINTER        :: conf

    INTEGER                                  :: i, ind
    INTEGER, DIMENSION(:, :), POINTER        :: tmp_conf

    ALLOCATE(tmp_conf(SIZE(conf,1),SIZE(main_conf)))
    tmp_conf = 0
    ind = 0
    DO i=1,SIZE(main_conf)
       IF (main_conf(i) /= 0) THEN
          ind = ind + 1
          tmp_conf(:,i) = conf(:,ind)
       END IF
    END DO
    DEALLOCATE(conf)
    ALLOCATE(conf(SIZE(tmp_conf,1),SIZE(tmp_conf,2)))
    conf = tmp_conf
    DEALLOCATE(tmp_conf)

  END SUBROUTINE make_plan_conf

! *****************************************************************************
!> \brief Writes restart for BSSE calculations
!> \param bsse_section ...
!> \param root_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2008 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE write_bsse_restart(bsse_section, root_section, error)

    TYPE(section_vals_type), POINTER         :: bsse_section, root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_bsse_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ires
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    ires = cp_print_key_unit_nr(logger,bsse_section,"PRINT%RESTART",&
         extension=".restart",  do_backup=.FALSE., file_position="REWIND",error=error)

    IF (ires>0) THEN
       CALL write_restart_header(ires)
       CALL section_vals_write(root_section, unit_nr=ires, hide_root=.TRUE., error=error)
    ENDIF

    CALL cp_print_key_finished_output(ires,logger,bsse_section,&
         "PRINT%RESTART", error=error)

  END SUBROUTINE write_bsse_restart

END MODULE bsse
