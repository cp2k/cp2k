!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****h* cp2k/bsse [1.0] *
!!
!!   NAME
!!     bsse
!!
!!   FUNCTION
!!     Module to perform a counterpoise correction (BSSE)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     6.2005 created [tlaino]
!!
!!   SOURCE
!****************************************************************************

MODULE bsse
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
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
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE string_utilities,                ONLY: compress
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bsse'

  PUBLIC :: do_bsse_calculation

CONTAINS

!!****f* bsse/do_bsse_calculation [1.0] *
!!
!!   NAME
!!     do_bsse_calculation
!!
!!   FUNCTION
!!     Perform an COUNTERPOISE CORRECTION (BSSE)
!!     For a 2-body system the correction scheme can be represented as:       
!!
!!     E_{AB}^{2}        = E_{AB}(AB) - E_A(AB) - E_B(AB)  [BSSE-corrected interaction energy]
!!     E_{AB}^{2,uncorr} = E_{AB}(AB) - E_A(A)  - E_B(B)
!!     E_{AB}^{CP}       = E_{AB}(AB) + [ E_A(A) - E_A(AB) ] + [ E_B(B) - E_B(AB) ]  
!!                                                         [CP-corrected total energy of AB]
!!      
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE do_bsse_calculation(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_bsse_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, num_of_conf, &
                                                Num_of_Frag, stat
    INTEGER, DIMENSION(:, :), POINTER        :: conf
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em
    TYPE(force_env_type), POINTER            :: force_env_loc
    TYPE(section_vals_type), POINTER         :: bsse_section, n_frags, &
                                                root_section

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(bsse_section, n_frags, Em, conf, force_env_loc)
       root_section => force_env%root_section
       bsse_section => section_vals_get_subs_vals(force_env%force_env_section,"BSSE",error=error)
       n_frags => section_vals_get_subs_vals(bsse_section,"FRAGMENT",error=error)
       CALL section_vals_get(n_frags,explicit=explicit,n_repetition=Num_of_Frag,error=error)
       !
       ! Number of configurations
       !
       num_of_conf = 0
       DO k = 1, Num_of_frag
          num_of_conf = num_of_conf + FACT(Num_of_frag) / (FACT(k)*FACT(Num_of_frag-k))
       END DO
       ALLOCATE(conf(num_of_conf,Num_of_frag),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Em(num_of_conf), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL gen_Nbody_conf(Num_of_frag, num_of_conf, conf)
       !
       ! Evaluating the energy of the N-body cluster terms
       !    
       DO i = 1, SIZE(conf,1)
          CALL eval_bsse_energy(conf(i,:),Em(i),force_env,force_env_loc, n_frags, &
               root_section,globenv, error=error)
       END DO
       CALL dump_bsse_results(conf, Em, num_of_frag,  bsse_section, error)
       DEALLOCATE(Em,   stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(conf, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE do_bsse_calculation
     
!!****f* bsse/eval_bsse_energy [1.0] *
!!
!!   NAME
!!     eval_bsse_energy
!!
!!   FUNCTION
!!     Evaluate the N-body energy contribution to the BSSE evaluation
!!      
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE eval_bsse_energy(conf, Em, force_env, force_env_loc, n_frags, root_section, globenv, error)
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf
    REAL(KIND=dp), INTENT(OUT)               :: Em
    TYPE(force_env_type), POINTER            :: force_env, force_env_loc
    TYPE(section_vals_type), POINTER         :: n_frags, root_section
    TYPE(global_environment_type), POINTER   :: globenv
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
       !
       ! Count the number of subconfiguration to evaluate..
       !
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
          CALL gen_Nbody_conf(Num_of_sub_frag, Num_of_sub_conf, conf_loc)
          CALL make_plan_conf(conf, conf_loc)           
          DO i = 1, Num_of_sub_conf
             CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf_loc(i,:), n_frags,&
                  root_section,globenv, Em_loc(i), error)             
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

!!****f* bsse/eval_bsse_energy_low [1.0] *
!!
!!   NAME
!!     eval_bsse_energy_low
!!
!!   FUNCTION
!!     Evaluate the N-body energy contribution to the BSSE evaluation
!!      
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [tlaino]
!!
!!*** ********************************************************************** 
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
    INTEGER :: i, isize, j, jj, k, method_name_id, my_targ, n_rep, &
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
             CALL section_vals_val_get(n_frags,"LIST", i_rep_section=i, n_rep_val=n_rep, error=error)
             IF (n_rep /= 0) THEN
                CALL section_vals_val_get(n_frags,"LIST", i_rep_section=i, i_vals=tmplist, error=error)
                CALL reallocate(atom_index,1,isize+SIZE(tmplist))
                atom_index(isize+1:isize+SIZE(tmplist)) = tmplist
                isize = SIZE(atom_index)
             END IF             
             CALL section_vals_val_get(n_frags,"RANGE", i_rep_section=i, n_rep_val=n_rep, error=error)
             IF (n_rep /= 0) THEN
                CALL section_vals_val_get(n_frags,"RANGE", i_rep_section=i, i_vals=tmplist, error=error)
                CALL reallocate(atom_index,1,isize+tmplist(2)-tmplist(1)+1)
                jj = 0
                DO j = tmplist(1), tmplist(2)
                   jj = jj + 1
                   atom_index(isize+jj:isize+jj) = j
                END DO
                isize = SIZE(atom_index)
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
               sub_atom_kind_name=atom_type,root_section=root_section,para_env=para_env,&
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

!!****f* bsse/dump_bsse_info [1.0] *
!!
!!   NAME
!!     dump_bsse_info
!!
!!   FUNCTION
!!     Dumps bsse information (configuration fragment)
!!      
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [tlaino]
!!
!!*** ********************************************************************** 
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

    CHARACTER(len=default_string_length)     :: conf_loc_s, conf_s
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
          WRITE (conf_s,*)conf ;         CALL compress(conf_s,full=.TRUE.)
          WRITE (conf_loc_s,*)conf_loc ; CALL compress(conf_loc_s,full=.TRUE.)

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

!!****f* bsse/conf_info_setup [1.0] *
!!
!!   NAME
!!     conf_info_setup
!!
!!   FUNCTION
!!     Read modified parameters for configurations
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE  conf_info_setup(present_charge, present_multpl, conf, conf_loc,&
       bsse_section, dft_section, error)
    INTEGER, INTENT(OUT)                     :: present_charge, present_multpl
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(section_vals_type), POINTER         :: bsse_section, dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'conf_info_setup', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: conf_loc_s, conf_s, &
                                                glb_conf_s, sub_conf_s
    INTEGER                                  :: i, nconf
    INTEGER, DIMENSION(:), POINTER           :: glb_conf, sub_conf
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: configurations

    present_charge = 0
    present_multpl = 0
    failure = .FALSE.
    NULLIFY(configurations, glb_conf, sub_conf)
    IF (.NOT.failure) THEN
       WRITE (conf_s,*)conf ;         CALL compress(conf_s,full=.TRUE.)
       WRITE (conf_loc_s,*)conf_loc ; CALL compress(conf_loc_s,full=.TRUE.)
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
             WRITE (glb_conf_s,*)glb_conf ; CALL compress(glb_conf_s,full=.TRUE.)
             WRITE (sub_conf_s,*)sub_conf ; CALL compress(sub_conf_s,full=.TRUE.)
             IF ((conf_s==glb_conf_s).AND.(conf_loc_s==sub_conf_s)) THEN
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


!!****f* bsse/dump_bsse_results [1.0] *
!!
!!   NAME
!!     dump_bsse_results
!!
!!   FUNCTION
!!     Dumps results
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** ********************************************************************** 
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

!!****f* bsse/gen_Nbody_conf [1.0] *
!!
!!   NAME
!!     gen_Nbody_conf
!!
!!   FUNCTION
!!     generate the N-body configuration for the N-body BSSE evaluation
!!      
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [tlaino]
!!
!!*** ********************************************************************** 
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

  RECURSIVE FUNCTION FACT(num) RESULT(my_fact)
    INTEGER, INTENT(IN)                      :: num
    INTEGER                                  :: my_fact

    IF (num<=1) THEN 
       my_fact = 1
    ELSE
       my_fact = num * FACT(num-1)
    END IF
  END FUNCTION FACT

  SUBROUTINE make_plan_conf (main_conf, conf)
    INTEGER, DIMENSION(:), INTENT(IN)        :: main_conf
    INTEGER, DIMENSION(:, :), POINTER        :: conf

    INTEGER                                  :: i, ind, stat
    INTEGER, DIMENSION(:, :), POINTER        :: tmp_conf

    ALLOCATE(tmp_conf(SIZE(conf,1),SIZE(main_conf)), stat=stat)
    tmp_conf = 0
    ind = 0
    DO i=1,SIZE(main_conf)
       IF (main_conf(i) /= 0) THEN 
          ind = ind + 1
          tmp_conf(:,i) = conf(:,ind)
       END IF
    END DO
    DEALLOCATE(conf, stat=stat)
    ALLOCATE(conf(SIZE(tmp_conf,1),SIZE(tmp_conf,2)), stat=stat)
    conf = tmp_conf
    DEALLOCATE(tmp_conf, stat=stat)

  END SUBROUTINE make_plan_conf

END MODULE bsse
