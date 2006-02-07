!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
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
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_methods,            ONLY: create_small_subsystem
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_qs
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE string_utilities,                ONLY: compress
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_bsse_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, num_of_conf, &
                                                Num_of_Frag, stat
    INTEGER, DIMENSION(:, :), POINTER        :: conf
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em
    TYPE(force_env_type), POINTER            :: force_env_loc
    TYPE(section_vals_type), POINTER         :: bsse_section, n_frags

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(bsse_section, n_frags, Em, conf, force_env_loc)
       bsse_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%BSSE",error=error)
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
               globenv, error=error)
       END DO
       CALL dump_bsse_results(conf, Em, num_of_frag, globenv, error)
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
  SUBROUTINE eval_bsse_energy(conf, Em, force_env, force_env_loc, n_frags, globenv, error)
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf
    REAL(KIND=dp), INTENT(OUT)               :: Em
    TYPE(force_env_type), POINTER            :: force_env, force_env_loc
    TYPE(section_vals_type), POINTER         :: n_frags
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
          CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf, n_frags, globenv, Em, error)
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
                  globenv, Em_loc(i), error)             
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
       globenv, energy, error)
    TYPE(force_env_type), POINTER            :: force_env, force_env_loc
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(section_vals_type), POINTER         :: n_frags
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_bsse_energy_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_type
    INTEGER                                  :: i, isize, j, jj, k, my_targ, &
                                                n_rep, num_of_frag, old_size, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: atom_index, atom_list, &
                                                my_conf, subsys_index, tmplist
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cp_subsystem_type), POINTER         :: subsys_loc
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: force_env_section

    failure = .FALSE.
    CALL section_vals_get(n_frags,n_repetition=num_of_frag,error=error)
    CPPostcondition(SIZE(conf)==num_of_frag,cp_failure_level,routineP,error,failure)
    NULLIFY(subsys_loc, subsys, particles, para_env, cell, atom_index, atom_type, tmplist,&
         force_env_section )
    force_env_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    IF (.NOT.failure) THEN
       ALLOCATE(my_conf(SIZE(conf)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_conf = conf
       CALL force_env_get(force_env=force_env, subsys=subsys, para_env=para_env,&
            cell=cell, error=error)
       CALL cp_subsys_get(subsys(1)%subsys, particles=particles, error=error)
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
       !
       ! Get names and modify the ghost ones
       !
       ALLOCATE(atom_type(isize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(subsys_index(isize), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       subsys_index = 1
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
       CALL dump_bsse_info(atom_index, atom_type, conf, conf_loc, globenv, error)
       !
       ! Let's start setting up environments and calculations
       !
       energy = 0.0_dp
       IF (globenv%prog_name_id == do_qs) THEN          

          CALL create_small_subsystem(subsys_loc, big_subsys=subsys,&
               small_para_env=para_env, small_cell=cell,sub_atom_index=atom_index,&
               sub_subsys_index=subsys_index, sub_atom_kind_name=atom_type,&
               globenv=globenv, force_env_section=force_env_section, error=error)

          CALL quickstep_create_force_env(force_env_loc,globenv,subsystem=subsys_loc,error=error)

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
       DEALLOCATE(subsys_index, stat=stat)
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
  SUBROUTINE  dump_bsse_info(atom_index, atom_type, conf, conf_loc, globenv, error)
    INTEGER, DIMENSION(:), POINTER           :: atom_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: atom_type
    INTEGER, DIMENSION(:), INTENT(IN)        :: conf, conf_loc
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_bsse_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: conf_loc_s, conf_s
    INTEGER                                  :: i, iw
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: bsse_section, print_info

    failure = .FALSE.
    NULLIFY(bsse_section, logger)
    IF (.NOT.failure) THEN
       logger => cp_error_get_logger(error)
       bsse_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%BSSE",error=error)
       print_info =>  section_vals_get_subs_vals(bsse_section,"PRINT%PROGRAM_RUN_INFO",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_info,error=error),cp_p_file).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          iw = globenv%scr
          WRITE (conf_s,*)conf ;         CALL compress(conf_s,full=.TRUE.)
          WRITE (conf_loc_s,*)conf_loc ; CALL compress(conf_loc_s,full=.TRUE.)

          WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T5,A,T30,A,T55,A,T80,A)")&
               "-","BSSE CALCULATION","FRAGMENT CONF: "//TRIM(conf_s),"FRAGMENT SUBCONF: "//TRIM(conf_loc_s),"-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,T60,A,T80,A)") "-","ATOM INDEX","ATOM NAME","-"
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,T60,A,T80,A)") "-","----------","---------","-"
          DO i = 1, SIZE(atom_index)
             WRITE (UNIT=iw,FMT="(T2,A,T20,I6,T61,A,T80,A)") "-",atom_index(i),TRIM(atom_type(i)),"-"
          END DO
          WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       END IF
    END IF
  END SUBROUTINE dump_bsse_info

  SUBROUTINE  dump_bsse_results(conf, Em,  num_of_frag, globenv, error)
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: conf
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Em
    INTEGER, INTENT(IN)                      :: num_of_frag
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_bsse_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iw
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: bsse_section, print_info

    failure = .FALSE.
    NULLIFY(bsse_section, logger)
    IF (.NOT.failure) THEN
       logger => cp_error_get_logger(error)
       bsse_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%BSSE",error=error)
       print_info =>  section_vals_get_subs_vals(bsse_section,"PRINT%PROGRAM_RUN_INFO",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_info,error=error),cp_p_file).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          iw = globenv%scr

          WRITE (UNIT=iw,FMT="(/,T2,A)") REPEAT("-",79)
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T36,A,T80,A)")&
               "-","BSSE RESULTS","-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          WRITE (UNIT=iw,FMT="(T2,A,T21,A,F15.9,T80,A)") "-","CP-corrected Total energy:",SUM(Em),"-"
          WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
          DO i = 1, SIZE(conf,1)
             IF (i.GT.1) THEN
                IF (SUM(conf(i-1,:))==1.AND.SUM(conf(i,:))/=1) THEN
                   WRITE (UNIT=iw,FMT="(T2,A,T80,A)") "-","-"
                END IF
             END IF
             WRITE (UNIT=iw,FMT="(T2,A,T24,I3,A,F16.9,T80,A)") "-",SUM(conf(i,:)),"-body contribution:",Em(i),"-"
          END DO
          WRITE (UNIT=iw,FMT="(T2,A,T20,A,F16.9,T80,A)") "-","BSSE  interaction  energy:",SUM(Em(Num_of_frag+1:)),"-"
          WRITE (UNIT=iw,FMT="(T2,A)") REPEAT("-",79)
       END IF
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
