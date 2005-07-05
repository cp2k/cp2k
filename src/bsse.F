!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp  
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE input_constants,                 ONLY: do_qs,&
                                             do_fist,&
                                             do_kg,&
                                             do_qmmm,&
                                             do_atom,&
                                             do_farming,&
                                             do_test

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='bsse'

  PUBLIC :: do_bsse_calculation

CONTAINS

!!****f* bsse/do_bsse_calculation [1.0] *
!!
!!   NAME
!!     do_bsse_calculation
!!
!!   FUNCTION
!!     Perform an COUNTERPOISE CORRECTION (BSSE)
!!       
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
    IMPLICIT NONE
    ! Arguments
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_bsse_calculation', &
      routineP = moduleN//':'//routineN    
    INTEGER :: Num_of_Frag, i, stat, k
    TYPE(section_vals_type), POINTER         :: bsse_section, n_frags
    LOGICAL :: explicit
    INTEGER,  POINTER, DIMENSION(:,:)        :: conf
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: Em
    INTEGER :: num_of_conf
    TYPE(force_env_type), POINTER            :: force_env_loc
    
    NULLIFY(bsse_section, n_frags, Em, conf, force_env_loc)
    bsse_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%BSSE",error=error)
    n_frags => section_vals_get_subs_vals(bsse_section,"FRAGMENT",error=error)
    CALL section_vals_get(n_frags,explicit=explicit,n_repetition=Num_of_Frag,error=error)
    !
    ! Number of configurations
    !
    WRITE(*,*)"Numbero of fragments ::",Num_of_Frag
    num_of_conf = 0
    DO k = 1, Num_of_frag
       num_of_conf = num_of_conf + FACT(Num_of_frag) / (FACT(k)*FACT(Num_of_frag-k))
    END DO
    ALLOCATE(conf(num_of_conf,Num_of_frag),stat=stat)
    ALLOCATE(Em(num_of_conf), stat=stat)
    CALL gen_Nbody_conf(Num_of_frag, num_of_conf, conf)
    !
    ! Evaluating the energy of the N-body cluster terms
    !    
    DO i = 1, SIZE(conf,1)
       WRITE(*,*)"Evaluating the configuration ::",conf(i,:)
       CALL eval_bsse_energy(conf(i,:),Em(i),force_env,force_env_loc, n_frags, &
            globenv, error=error)
    END DO

    DEALLOCATE(Em,   stat=stat)
    DEALLOCATE(conf, stat=stat)
    
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
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN), DIMENSION(:) :: conf
    REAL(KIND=dp), INTENT(OUT)        :: Em
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(force_env_type), POINTER            :: force_env_loc
    TYPE(section_vals_type), POINTER         :: n_frags
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_bsse_energy', &
      routineP = moduleN//':'//routineN    
    INTEGER,  POINTER, DIMENSION(:,:)        :: conf_loc
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: Em_loc
    INTEGER :: Num_of_sub_conf, Num_of_sub_frag, k, i, stat
    
    NULLIFY(conf_loc, Em_loc)
    !
    ! Count the number of subconfiguration to evaluate..
    !
    Num_of_sub_frag = COUNT(conf == 1)    
    Num_of_sub_conf = 0
    IF (Num_of_sub_frag == 1) THEN
       WRITE(*,*)"      Evaluating the single term configuration energy:",conf
       CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf, n_frags, globenv, Em, error)
    ELSE
       DO k = 1, Num_of_sub_frag
          Num_of_sub_conf = Num_of_sub_conf + FACT(Num_of_sub_frag) / (FACT(k)*FACT(Num_of_sub_frag-k))
       END DO
       ALLOCATE(conf_loc(Num_of_sub_conf,Num_of_sub_frag), stat=stat)
       ALLOCATE(Em_loc(Num_of_sub_conf), stat=stat)
       CALL gen_Nbody_conf(Num_of_sub_frag, Num_of_sub_conf, conf_loc)
       CALL make_plan_conf(conf, conf_loc) 
   
       DO i = 1, Num_of_sub_conf
             WRITE(*,*)"      Evaluating the sub configuration energy::",conf_loc(i,:)
             CALL eval_bsse_energy_low(force_env, force_env_loc, conf, conf_loc(i,:), n_frags,&
                  globenv, Em_loc(i), error)
       END DO
       !
       ! Scrivere la parte di sum up dell'energia...
       !
    END IF

    DEALLOCATE(Em_loc, stat=stat)
    DEALLOCATE(conf_loc, stat=stat)
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
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN), DIMENSION(:) :: conf, conf_loc
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(force_env_type), POINTER            :: force_env_loc
    TYPE(section_vals_type), POINTER         :: n_frags
    TYPE(global_environment_type), POINTER   :: globenv
    REAL(KIND=dp), INTENT(OUT)               :: energy
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local Variables
    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_bsse_energy_low', &
      routineP = moduleN//':'//routineN    

    IF (globenv%prog_name_id /= do_qs) STOP

!!!TL    CALL create_small_subsystem(subsys_loc,&
!!!TL         big_subsys=big_subsys,small_para_env=para_env,&
!!!TL         small_cell=cell,sub_atom_index=atom_index,&
!!!TL         sub_subsys_index=subsys_index,&
!!!TL         sub_atom_kind_name=atom_type, globenv=globenv, error=error)
!!!TL    CALL read_atomic_kind_set(subsys_loc%atomic_kinds%els, globenv)
!!!TL    
!!!TL    CALL qs_env_create(qs_env, para_env=para_env, error=error)
!!!TL    
!!!TL    CALL set_qs_env(qs_env, input=force_env_section, error=error)
!!!TL    
!!!TL    CALL qs_init_subsys(qs_env,subsystem=subsys_loc,cell=cell,&
!!!TL         cell_ref=cell, use_ref_cell=.FALSE.,qmmm=.FALSE..,&
!!!TL         globenv=globenv,green=green, error=error)
!!!TL    
!!!TL    CALL cp_subsys_release(subsys_loc,error=error)
!!!TL    CALL force_env_create(force_env_loc, qs_env=qs_env, globenv=globenv, error=error)
!!!TL    
!!!TL    !
!!!TL    ! Evaluate Energy
!!!TL    !
!!!TL    CALL force_env_calc_energy_force(force_env_loc,calc_force=.FALSE.,error=error)
!!!TL    CALL force_env_get(force_env_loc,potential_energy=energy,error=error)
!!!TL    CALL force_env_release(force_env_loc, error=error)
!!!TL    NULLIFY(force_env_loc)

  END SUBROUTINE eval_bsse_energy_low

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
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: Num_of_frag, num_of_conf
    INTEGER, POINTER, DIMENSION(:,:) :: conf
    ! Local Variables
    INTEGER :: k,  my_ind
   
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
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: ldown, lup, k
    INTEGER, INTENT(INOUT) :: my_ind
    INTEGER, POINTER, DIMENSION(:,:) :: conf
    INTEGER :: i, kloc, my_ind0
    
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
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: num
    INTEGER              :: my_fact

    IF (num<=1) THEN 
       my_fact = 1
    ELSE
       my_fact = num * FACT(num-1)
    END IF
  END FUNCTION FACT

  SUBROUTINE make_plan_conf (main_conf, conf)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN), DIMENSION(:)        :: main_conf
    INTEGER,  POINTER, DIMENSION(:,:)        :: conf
    ! Local Variables
    INTEGER,  POINTER, DIMENSION(:,:)        :: tmp_conf
    INTEGER :: stat, ind, i
    
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
