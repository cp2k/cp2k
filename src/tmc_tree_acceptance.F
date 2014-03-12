!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief tree nodes acceptance
!>        code is seperated in 3 parts, first the acceptance criteria,
!>        second the tree node acceptance handling, searching etc. and
!>        than the acceptance probability handling
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_tree_acceptance
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: boltzmann,&
                                             joule
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_calculations,                ONLY: compute_estimated_prob,&
                                             get_scaled_cell
  USE tmc_dot_tree,                    ONLY: create_dot_color,&
                                             create_global_tree_dot_color
  USE tmc_file_io,                     ONLY: write_result_list_element
  USE tmc_move_handle,                 ONLY: prob_update
  USE tmc_move_types,                  ONLY: mv_type_MD,&
                                             mv_type_volume_move
  USE tmc_stati,                       ONLY: task_type_gaussian_adaptation
  USE tmc_tree_build,                  ONLY: remove_unused_g_tree
  USE tmc_tree_search,                 ONLY: get_subtree_elements_to_check,&
                                             search_canceling_elements,&
                                             search_next_gt_element_to_check,&
                                             search_parent_element
  USE tmc_tree_types,                  ONLY: &
       add_to_list, global_tree_type, gt_elem_list_type, status_accepted, &
       status_accepted_result, status_calc_approx_ener, status_calculate_MD, &
       status_calculate_NMC_steps, status_calculate_energy, &
       status_calculated, status_cancel_ener, status_cancel_nmc, &
       status_canceled_ener, status_canceled_nmc, status_created, &
       status_deleted, status_deleted_result, status_rejected, &
       status_rejected_result, tree_type
  USE tmc_types,                       ONLY: tmc_env_type,&
                                             tmc_param_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_tree_acceptance'

  PUBLIC :: acceptance_check
  PUBLIC :: check_acceptance_of_depending_subtree_nodes, &
            check_elements_for_acc_prob_update
  PUBLIC :: tree_update
  
  INTEGER, PARAMETER :: DEBUG = 0

CONTAINS

  !============================================================================
  ! acceptance criteria calculations
  !============================================================================
! *****************************************************************************
!> \brief standard Monte Carlo and 2 potential acceptance check
!>        acceptance check of move from old(last accepted) to new configuration
!>        the sum of kinetic and potential energy is used
!>        acc(o->n)=min(1,exp( -beta*(H(n)-H(o)) ))
!> \param tree_element new/actual configuration
!> \param parent_element last accepted configuration
!> \param tmc_params TMC global parameters
!> \param temperature actual temperature configuration should be checked with
!> \param diff_pot_check 2potential check or not?
!> \param accept result (configuration accepted of rejected)
!> \param rnd_nr random number for acceptance check
!> \param approx_ener for NMC the approximated energies schould be used
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE acceptance_check(tree_element, parent_element, tmc_params, &
                              temperature, diff_pot_check, accept, rnd_nr, &
                              approx_ener, error)
    TYPE(tree_type), POINTER                 :: tree_element, parent_element
    TYPE(tmc_param_type), POINTER            :: tmc_params
    REAL(KIND=dp)                            :: temperature
    LOGICAL                                  :: diff_pot_check, accept
    REAL(KIND=dp)                            :: rnd_nr
    LOGICAL, OPTIONAL                        :: approx_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'acceptance_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ekin_last_acc, elem_ener, kB, &
                                                parent_ener

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(tree_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(parent_element),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(temperature.GT.0.0_dp,cp_failure_level,routineP,error,failure)
    CPPrecondition(rnd_nr.GE.0.0_dp.AND.rnd_nr.LE.1.0_dp,cp_failure_level,routineP,error,failure)

    kB = boltzmann/joule

    ! start the timing
    CALL timeset(routineN,handle)

    IF(tmc_params%task_type.EQ.task_type_gaussian_adaptation) THEN
CPPrecondition(.FALSE. ,cp_failure_level,routineP,error,failure)
!TODO       CALL acc_check_gauss_adapt(f=tree_element%potential, ct=tree_element%ekin_before_md, acc=accept)
!DO NOT DO       RETURN
    END IF

    !-- using 2 different potentials, the energy difference between these potentials 
    !   and the two configurations have to be regarded. 
    !   The ensamble should be in equilibrium state of the approximate potential
    IF(diff_pot_check.AND.(tmc_params%NMC_inp_file.NE."")) THEN
      IF((tree_element%potential.EQ.HUGE(tree_element%potential)).OR.&
         tree_element%e_pot_approx.EQ.HUGE(tree_element%e_pot_approx)) THEN
        elem_ener = HUGE(elem_ener)
      ELSE
        !for different potentials we have to regard the differences in energy
        ! min(1,e^{-\beta*[(E_{exact}(n)-E_{approx}(n))-(E_{exact}(o)-E_{approx}(o))]})
        elem_ener   = 1.0_dp/(kB*temperature)*tree_element%potential &
                    - 1.0_dp/(kB*tmc_params%Temp(tree_element%temp_created))&
                      *tree_element%e_pot_approx
      END IF
      parent_ener = 1.0_dp/(kB*temperature)*parent_element%potential &
                  - 1.0_dp/(kB*tmc_params%Temp(tree_element%temp_created))&
                    *parent_element%e_pot_approx

      !-- always accepted if new energy is smaller than old energy
      IF(elem_ener.LE.parent_ener)THEN
        accept = .TRUE.
      ELSE
        !-- gaussian distributed acceptance if new energy is greater than old energy
        IF(rnd_nr.LT.&
          EXP(-(elem_ener-parent_ener)))THEN
          accept = .TRUE.
        ELSE
          accept = .FALSE.
        END IF
      END IF
    ELSE
      !-- standard MC check, but also using the kinetic energy for Hybrid Monte Carlo
      IF(tree_element%move_type.EQ.mv_type_MD)THEN
        ekin_last_acc = tree_element%ekin_before_md
                        ! using the kinetic energy before velocity change 
      ELSE
        ekin_last_acc = parent_element%ekin
      END IF
      ! comparing aproximated energies
      IF(PRESENT(approx_ener))THEN
        elem_ener   = tree_element%e_pot_approx &
                    + tree_element%ekin
        parent_ener = parent_element%e_pot_approx &
                    + ekin_last_acc
      ELSE
        elem_ener   = tree_element%potential &
                    + tree_element%ekin
        parent_ener = parent_element%potential &
                    + ekin_last_acc
      END IF

      !-- always accepted if new energy is smaller than old energy
      IF(elem_ener.LE.parent_ener)THEN
        accept = .TRUE.
      ELSE
        !-- gaussian distributed acceptance if new energy is greater than old energy
        IF(rnd_nr.LT.&
          EXP(-1.0_dp/(kB*temperature)*(elem_ener-parent_ener)))THEN
          accept = .TRUE.
        ELSE
           accept = .FALSE.
        END IF
      END IF
    END IF

    ! update the estimated energy acceptance probability distribution
    IF(diff_pot_check) THEN
      CPPrecondition(ASSOCIATED(tmc_params%prior_NMC_acc),cp_failure_level,routineP,error,failure)
      tmc_params%prior_NMC_acc%counter = tmc_params%prior_NMC_acc%counter +1
      tmc_params%prior_NMC_acc%aver    = (tmc_params%prior_NMC_acc%aver*(tmc_params%prior_NMC_acc%counter-1)+&
                                         ((elem_ener-parent_ener)))/REAL(tmc_params%prior_NMC_acc%counter,KIND=dp)
      tmc_params%prior_NMC_acc%aver_2  = (tmc_params%prior_NMC_acc%aver_2*(tmc_params%prior_NMC_acc%counter-1)+&
                                         (elem_ener-parent_ener)**2)/REAL(tmc_params%prior_NMC_acc%counter,KIND=dp)
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE acceptance_check

! *****************************************************************************
!> \brief parallel tempering swap acceptance check, 
!>        check swapped configurations of different temperatures
!>        using sum of potential and kinetic energy
!>        always accepted if energy of configuration with higher temperature 
!>        is smaller than energy of configuration with lower temperature
!>        acc(o->n)=min(1,exp((beta_i-beta_j)*(U_i-U_j)))
!> \param tree_elem actual global tree element
!> \param conf1 sub tree element of lower temperature
!> \param conf2 sub tree element of higher temperature
!> \param tmc_params TMC environment parameters
!> \param accept acceptance of configurational change
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE swap_acceptance_check(tree_elem, conf1, conf2, tmc_params, accept, error)
    TYPE(global_tree_type), POINTER          :: tree_elem
    TYPE(tree_type), POINTER                 :: conf1, conf2
    TYPE(tmc_param_type), POINTER            :: tmc_params
    LOGICAL                                  :: accept
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'swap_acceptance_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: delta, elem1_ener, &
                                                elem2_ener, kB, vol1, vol2

    failure = .FALSE.
    kB = boltzmann/joule

    CPPrecondition(ASSOCIATED(tree_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(tree_elem%rnd_nr.GE.0.0_dp,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(conf1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(conf2),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    IF(tmc_params%pressure.GT.0.0_dp) THEN
      ! pt-NVT
      elem1_ener = conf1%potential &
                 + conf1%ekin
      elem2_ener = conf2%potential &
                 + conf2%ekin
      ! the swap is done with prob: exp((\beta_i-\beta_j)(U_i-U_j)), 
      ! BUT because they are already swaped we exchange the energies.
      IF(tree_elem%rnd_nr.LT.EXP((1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf))-&
                                  1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf+1)))&
                                  *(elem2_ener-elem1_ener)))THEN
         accept = .TRUE.
      ELSE
         accept = .FALSE.
      END IF
    ELSE
      ! pt-NpT (parallel Tempering with constant pressure, temperature and num. particle)
      CALL get_scaled_cell(cell=tmc_params%cell, box_scale=conf1%box_scale, &
                           vol=vol1, error=error)
      CALL get_scaled_cell(cell=tmc_params%cell, box_scale=conf2%box_scale, &
                           vol=vol2, error=error)
      ! delta= (beta_m-beta_n)(U_n-U_m) + (beta_m*P_m-beta_n*P_n)(V_n-V_m)
      delta = ( 1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf)) &
                - 1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf+1)) ) * &
              ( (conf2%potential+conf2%ekin) - (conf1%potential+conf1%ekin) ) + &
              ( 1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf))*tmc_params%pressure &
                - 1.0_dp/(kB*tmc_params%Temp(tree_elem%mv_conf+1))*tmc_params%pressure) * &
              ( vol2 - vol1 )
       IF(tree_elem%rnd_nr.LT.EXP(delta))THEN
          accept = .TRUE.
       ELSE
          accept = .FALSE.
       END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE swap_acceptance_check

! *****************************************************************************
!> \brief volume move acceptance check
!>        sampling NPT, we need to do volume moves. Their acceptance are 
!>        checked regarding the old and new energy, the volume difference 
!>        and the constant pressure
!> \param parent_elem old tree element with old box size
!> \param new_elem actual tree element with new box size
!> \param temperature actual temperature
!> \param rnd_nr random number to check with
!> \param tmc_params TMC environment parameters
!> \param accept Monte Carlo move acceptance (result)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE volume_acceptance_check(parent_elem, new_elem, temperature, &
                                     rnd_nr, tmc_params, accept, error)
    TYPE(tree_type), POINTER                 :: parent_elem, new_elem
    REAL(KIND=dp)                            :: temperature, rnd_nr
    TYPE(tmc_param_type), POINTER            :: tmc_params
    LOGICAL                                  :: accept
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'volume_acceptance_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d_enthalpy, kB, n_vol, p_vol

    failure = .FALSE.
    kB = boltzmann/joule

    CPPrecondition(ASSOCIATED(parent_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(new_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(temperature.GT.0.0_dp,cp_failure_level,routineP,error,failure)
    CPPrecondition(rnd_nr.GT.0.0_dp,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)
    CPPrecondition(tmc_params%pressure.GE.0.0_dp,cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    CALL get_scaled_cell(cell=tmc_params%cell, box_scale=parent_elem%box_scale,&
                         vol=p_vol, error=error)
    CALL get_scaled_cell(cell=tmc_params%cell, box_scale=new_elem%box_scale, &
                         vol=n_vol, error=error)
    ! H=enthalpy, U=energy, P=pressure, V=volume, kB=Boltzmann constant, T=temperature, N=amount particle
    ! delta_H  =      delta_U                               + P*delta_V           - kB*T*N*ln(V_n/V_p) 
    IF(.FALSE.) THEN
      ! the volume move in volume space (dV)
      d_enthalpy = (new_elem%potential-parent_elem%potential) + &
                   tmc_params%pressure*(n_vol-p_vol) - &
                   kB*temperature*(SIZE(new_elem%pos)/&
                     tmc_params%dim_per_elem)*&
                   LOG(n_vol/p_vol)
    ELSE
      IF(tmc_params%v_isotropic) THEN
        d_enthalpy = (new_elem%potential-parent_elem%potential) + &
                     tmc_params%pressure*(n_vol-p_vol) - &
                     kB*temperature*((SIZE(new_elem%pos)/&
                       tmc_params%dim_per_elem)+2/REAL(3,KIND=dp))*&
                     LOG(n_vol/p_vol)
      ELSE
        d_enthalpy = (new_elem%potential-parent_elem%potential) + &
                     tmc_params%pressure*(n_vol-p_vol) - &
                     kB*temperature*(SIZE(new_elem%pos)/&
                       tmc_params%dim_per_elem)*&
                     LOG(n_vol/p_vol)
      END IF
    END IF
    ! acc(o->n) = min(1, exp(-beta*delta_H))
    IF(d_enthalpy .LE. 0.0_dp)THEN
       accept = .TRUE.
    ELSE
       IF(rnd_nr .LT. EXP(-1.0_dp/(kB*temperature)*d_enthalpy))THEN
          accept = .TRUE.
       ELSE
          accept = .FALSE.
       END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE volume_acceptance_check

  !============================================================================
  ! tree node acceptance
  !============================================================================
! *****************************************************************************
!> \brief check acceptance of energy calculated element and related childs,
!>        when ready 
!> \param tree_elem actual tree element with calculated energy
!> \param tmc_params TMC environment parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE check_acceptance_of_depending_subtree_nodes(tree_elem, tmc_env, error)
    TYPE(tree_type), POINTER                 :: tree_elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'check_acceptance_of_depending_subtree_nodes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, parent_ready
    TYPE(tree_type), POINTER                 :: act_elem, parent_elem

    failure = .FALSE.
    NULLIFY(parent_elem, act_elem)

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tree_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(tree_elem%stat.EQ.status_calculated,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tree_elem%parent),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    act_elem => tree_elem

    ! ------------------------------------------------------
    ! check this element 
    parent_elem => search_parent_element(act_elem, error)
    CPPostconditionNoFail(.NOT.ASSOCIATED(act_elem, parent_elem),cp_failure_level,routineP,error)

    ! check status of parent element
    SELECT CASE(parent_elem%stat)
    CASE(status_created, status_calculate_energy, status_calculate_MD, &
         status_calculate_NMC_steps, status_canceled_ener, &
         status_canceled_nmc, status_cancel_nmc, status_cancel_ener)
      parent_ready = .FALSE.
    CASE(status_accepted_result, status_rejected_result, &
         status_accepted, status_rejected, status_calculated)
      parent_ready = .TRUE.
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown parent element status"//cp_to_string(parent_elem%stat),error,failure)
    END SELECT

    ! ready ?
    IF(parent_ready) THEN
    ! acceptance check
      CALL check_and_change_status_of_subtree_elem(act_elem=act_elem, &
               parent_elem=parent_elem, tmc_env=tmc_env,error=error)
    END IF
    !------------------------------------------------------
    ! check acc child
    parent_elem => tree_elem
    IF(ASSOCIATED(tree_elem%acc)) THEN
      act_elem => tree_elem%acc
      IF(act_elem%stat .EQ. status_calculated) THEN
         CALL check_and_change_status_of_subtree_elem(act_elem=act_elem, &
                  parent_elem=parent_elem, tmc_env=tmc_env,error=error)
      END IF
      !------------------------------------------------------
      ! check all nacc childs
      nacc_loop: DO
        IF(.NOT.ASSOCIATED(act_elem%nacc)) EXIT nacc_loop
        act_elem => act_elem%nacc
        IF(act_elem%stat .EQ. status_calculated) THEN
          CALL check_and_change_status_of_subtree_elem(act_elem=act_elem, &
                   parent_elem=parent_elem, tmc_env=tmc_env, error=error)
        END IF
      END DO nacc_loop
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE check_acceptance_of_depending_subtree_nodes

! *****************************************************************************
!> \brief checking the elements and change the status
!> \param act_elem actual tree element (new configuration)
!> \param parent_elem parent tree element (old configuration)
!> \param tmc_params TMC environment parameters
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE check_and_change_status_of_subtree_elem(act_elem, parent_elem, tmc_env, error)
    TYPE(tree_type), POINTER                 :: act_elem, parent_elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'check_and_change_status_of_subtree_elem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, flag, result_acc
    TYPE(global_tree_type), POINTER          :: tmp_gt_ptr
    TYPE(gt_elem_list_type), POINTER         :: tmp_gt_list_ptr

    failure = .FALSE.
    NULLIFY(tmp_gt_list_ptr, tmp_gt_ptr)

    ! start the timing
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%params),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(act_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(parent_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(act_elem%stat.EQ.status_calculated,cp_failure_level,routineP,error,failure)

    flag = .FALSE.

    tmp_gt_list_ptr => act_elem%gt_nodes_references
    ! check for all global tree elements refering to this subtree element
    ! same subtree element can be accepted at a certain temperature and rejected at another one
    DO WHILE(ASSOCIATED(tmp_gt_list_ptr))
      CPPrecondition(tmp_gt_list_ptr%gt_elem%stat.NE.status_accepted_result,cp_failure_level,routineP,error,failure)
      CPPrecondition(tmp_gt_list_ptr%gt_elem%stat.NE.status_rejected_result,cp_failure_level,routineP,error,failure)

      CALL check_elements(gt_act_elem=tmp_gt_list_ptr%gt_elem, tmc_env=tmc_env,&
                          check_done=flag, result_acc=result_acc, error=error)
      IF(flag)THEN
        ! probability update 
        CALL prob_update(move_types=tmc_env%params%move_types, elem=act_elem, &
                 acc=result_acc, prob_opt=tmc_env%params%esimate_acc_prob, error=error)

        ! change status
        ! accepted case: delete (and cancel) not accepted branch
        IF(result_acc .EQV. .TRUE.) THEN
          tmp_gt_list_ptr%gt_elem%stat = status_accepted
          IF(ASSOCIATED(tmp_gt_list_ptr%gt_elem%nacc))THEN
            tmp_gt_ptr => tmp_gt_list_ptr%gt_elem%nacc
            CALL remove_unused_g_tree(begin_ptr=tmp_gt_ptr, &
                     end_ptr=tmp_gt_list_ptr%gt_elem, removed=flag, &
                     tmc_env=tmc_env, error=error)
          END IF
        ELSE
        ! not accepted case: delete (and cancel) accepted branch
          tmp_gt_list_ptr%gt_elem%stat = status_rejected
          IF(ASSOCIATED(tmp_gt_list_ptr%gt_elem%acc))THEN
            tmp_gt_ptr => tmp_gt_list_ptr%gt_elem%acc
            CALL remove_unused_g_tree(begin_ptr=tmp_gt_ptr, &
                     end_ptr=tmp_gt_list_ptr%gt_elem, removed=flag, &
                     tmc_env=tmc_env, error=error)
          END IF
        END IF

        IF(tmc_env%params%DRAW_TREE) &
          CALL create_global_tree_dot_color(gt_tree_element=tmp_gt_list_ptr%gt_elem, &
                                 tmc_params=tmc_env%params, error=error)
      END IF
      tmp_gt_list_ptr => tmp_gt_list_ptr%next
    END DO
    ! end the timing 
    CALL timestop(handle)
  END SUBROUTINE check_and_change_status_of_subtree_elem

! *****************************************************************************
!> \brief change status flag of all subtree elements
!>        most important to draw the right subtrees
!> \param gt_ptr global tree pointer, the related/changed sub tree element 
!>        should change status
!> \param stat the status of the global tree element
!> \param tmc_params TMC emvironment parameters for drawing
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE subtree_configuration_stat_change(gt_ptr, stat, tmc_params, error)
    TYPE(global_tree_type), POINTER          :: gt_ptr
    INTEGER                                  :: stat
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'subtree_configuration_stat_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(tree_type), POINTER                 :: ptr

    failure = .FALSE.
    NULLIFY(ptr)

    CPPrecondition(ASSOCIATED(gt_ptr),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_params),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! don't respect swaped nodes (subtree nodes are already in the list, 
    !   and we dont want them in marked in subtrees)
    IF(.NOT.gt_ptr%swaped) THEN
      ptr => gt_ptr%conf(gt_ptr%mv_conf)%elem
      ! dependent on the status (don't map each status 1:1, 
      !   e.g. not the result ones)
      SELECT CASE(stat)
      CASE(status_accepted_result)
        ptr%stat = status_accepted
      CASE(status_rejected_result)
        ptr%stat = status_rejected
      CASE(status_accepted, status_rejected)
        ptr%stat = stat
      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown global tree status"//&
                 cp_to_string(stat),error,failure)
      END SELECT

      IF(tmc_params%DRAW_TREE) &
        CALL create_dot_color(tree_element=ptr, tmc_params=tmc_params, &
                              error=error)
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE subtree_configuration_stat_change

! *****************************************************************************
!> \brief checks if the element is ready for an acceptance check
!> \param elem sub tree element
!> \param ready return value
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE elem_ready_to_check(elem, ready, error)
    TYPE(tree_type), POINTER                 :: elem
    LOGICAL                                  :: ready
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'elem_ready_to_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ready = .FALSE.

    SELECT CASE(elem%stat)
    CASE(status_created, status_calculate_energy, &
         status_calculate_MD, status_calculate_NMC_steps, status_calc_approx_ener, &
         status_cancel_nmc, status_cancel_ener, status_canceled_nmc, status_canceled_ener)
       ready = .FALSE.
    CASE(status_calculated, status_accepted_result, status_accepted, status_rejected, status_rejected_result)
       ready = .TRUE.
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"status of actual tree node is unknown"&
               //cp_to_string(elem%stat),error,failure)
    END SELECT
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE elem_ready_to_check

! *****************************************************************************
!> \brief checking the elements (find element and do acceptance check)
!> \param gt_act_elem actual global tree element 
!> \param tmc_env TMC environment
!> \param check_done successful checked? (result)
!> \param result_acc checked configuration accepted? (result)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE check_elements(gt_act_elem, tmc_env, check_done, result_acc, error)
    TYPE(global_tree_type), POINTER          :: gt_act_elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    LOGICAL                                  :: check_done, result_acc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_elements', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: act_ready, failure, tmp_ready
    TYPE(tree_type), POINTER                 :: act_element, tmp_element

    failure = .FALSE.
    NULLIFY(act_element, tmp_element)
    check_done = .FALSE.

    CPPrecondition(ASSOCIATED(gt_act_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env%m_env),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    !====================================================================================
    ! check if global tree element is already checked
    SELECT CASE(gt_act_elem%stat)
    CASE(status_accepted, status_rejected)
      check_done = .TRUE.
      IF(gt_act_elem%stat .EQ. status_accepted) THEN
        result_acc = .TRUE.
      ELSE IF(gt_act_elem%stat .EQ. status_rejected) THEN
        result_acc = .FALSE.
      ELSE
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,"undefinite status of already checked elements:"&
                //cp_to_string(gt_act_elem%stat),error,failure)
      END IF
    CASE DEFAULT
    END SELECT

    IF(.NOT.check_done) THEN
      !====================================================================================
      ! get elements related to global tree element to check
      CALL get_subtree_elements_to_check(gt_act_elem=gt_act_elem, elem1=act_element, &
                                         elem2=tmp_element, error=error)
      CPPostcondition(ASSOCIATED(act_element),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(tmp_element),cp_failure_level,routineP,error,failure)
  
      ! check status of both elements (if they are already calculated and hence ready to check)
      CALL elem_ready_to_check(elem=act_element, ready=act_ready, error=error)
      CALL elem_ready_to_check(elem=tmp_element, ready=tmp_ready, error=error)
  
      ! both ready ? check
      IF(act_ready .AND. tmp_ready) THEN
        ! 2 temperature (swap) check
        IF(gt_act_elem%swaped)THEN
          CALL swap_acceptance_check(tree_elem=gt_act_elem, conf1=act_element, &
                                     conf2=tmp_element, tmc_params=tmc_env%params, &
                                     accept=result_acc, error=error)
        ! volume move check
        ELSE IF(act_element%move_type.EQ.mv_type_volume_move) THEN
          CALL volume_acceptance_check(parent_elem=tmp_element, new_elem=act_element, &
                                       temperature=tmc_env%params%Temp(gt_act_elem%mv_conf),&
                                       rnd_nr=gt_act_elem%rnd_nr, &
                                       tmc_params=tmc_env%params, accept=result_acc, &
                                       error=error)
        ELSE
           IF(tmc_env%m_env%temp_decrease.NE.1.0_dp)THEN
             CALL acceptance_check(tree_element=act_element, parent_element=tmp_element, &
                                   tmc_params=tmc_env%params, temperature=gt_act_elem%Temp, &
                                   diff_pot_check=.TRUE., accept=result_acc, &
                                   rnd_nr=gt_act_elem%rnd_nr, error=error)
           ELSE
             CALL acceptance_check(tree_element=act_element, parent_element=tmp_element, &
                                   tmc_params=tmc_env%params, &
                                   temperature=tmc_env%params%Temp(gt_act_elem%mv_conf), &
                                   diff_pot_check=.TRUE., accept=result_acc, &
                                   rnd_nr=gt_act_elem%rnd_nr, error=error)
           END IF
        END IF
        check_done = .TRUE.
      ELSE
        ! if element is not ready, update status of global tree element
        !TODO maybe update this part (how much status we need in global tree from sub tree??)
        SELECT CASE(gt_act_elem%stat)
        !-- check if at least the MD move is already done
        !-- hence new configurations can be created on basis of this configuration
        CASE(status_calculate_MD, status_calculate_NMC_steps, status_calculate_energy, &
             status_created, status_calc_approx_ener)
          IF(gt_act_elem%conf(gt_act_elem%mv_conf)%elem%stat.NE.gt_act_elem%stat) &
            gt_act_elem%stat = gt_act_elem%conf(gt_act_elem%mv_conf)%elem%stat
        CASE(status_calculated)
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"status of actual checked node is unknown"&
               //cp_to_string(gt_act_elem%stat),error,failure)
        END SELECT
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE check_elements

! *****************************************************************************
!> \brief searching tree nodes to check for Markov Chain, elements are marked 
!>        and stored in lists ... (main entry point)
!> \param 
!> \param tmc_env TMC environment
!> \param 
!> \param result_acc checked configuration accepted? (result)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE tree_update(tmc_env, result_acc, something_updated, &
                         error)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    LOGICAL                                  :: result_acc, something_updated
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tree_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, itmp
    LOGICAL                                  :: failure, found
    TYPE(global_tree_type), POINTER          :: gt_act_elem
    TYPE(tree_type), POINTER                 :: act_element, tmp_element

    failure = .FALSE.
    NULLIFY(gt_act_elem, tmp_element, act_element)

    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    result_acc = .FALSE.
    something_updated = .FALSE.

    gt_act_elem => tmc_env%m_env%gt_act
    search_calculated_element_loop: DO
      NULLIFY(act_element, tmp_element)
      ! search for element to check
      CALL search_next_gt_element_to_check(ptr=gt_act_elem, found=found, error=error)
      IF(.NOT.found) EXIT search_calculated_element_loop

      ! check the elements status end, if possible, do acceptance check
      CALL check_elements(gt_act_elem=gt_act_elem, tmc_env=tmc_env, &
                          check_done=found, result_acc=result_acc,error=error)
      ! if check was not possible, update the status of the global tree element and return
      IF(.NOT.found) EXIT search_calculated_element_loop

      ! get elements related to global tree element, which were checked 
      CALL get_subtree_elements_to_check(gt_act_elem=gt_act_elem, &
                                         elem1=act_element, elem2=tmp_element,&
                                         error=error)

      !========================================================================
      !-- set result counters
      ! counter for certain temperature
      tmc_env%m_env%result_count(gt_act_elem%mv_conf) = &
             tmc_env%m_env%result_count(gt_act_elem%mv_conf) +1
      ! in case of swapped also count the result for 
      !  the other swapped temperature
      IF(gt_act_elem%swaped)&
        tmc_env%m_env%result_count(gt_act_elem%mv_conf+1) = &
               tmc_env%m_env%result_count(gt_act_elem%mv_conf+1) +1
      ! count also for global tree Markov Chain
      tmc_env%m_env%result_count(0) = tmc_env%m_env%result_count(0) +1

      ! flag for doing tree cleaning with canceling and deallocation
      something_updated = .TRUE.

      !IF((.NOT.gt_act_elem%swaped) .AND. act_element%temp_created.NE.0 .AND. &
      !   act_element%temp_created.NE.gt_act_elem%mv_conf)&
      !  count_wrong_temp_move(gt_act_elem%mv_conf) = &
      !       count_wrong_temp_move(gt_act_elem%mv_conf)+1

      !========================================================================
      !-- sort element in lists
      !-- case NOT ACCEPTED
      IF(.NOT.result_acc)THEN               !result NOT accepted
        IF(gt_act_elem%prob_acc.GT.0.0_dp) THEN
          IF(gt_act_elem%prob_acc.GE.0.5_dp) THEN
            ! wrong estimate (estimated accepted)
            tmc_env%m_env%estim_corr_wrong(4) = tmc_env%m_env%estim_corr_wrong(4) + 1
            IF(DEBUG.GT.0) WRITE(tmc_env%m_env%io_unit,*)&
              "Wrong guess for NACC (elem/estim acc prob)", gt_act_elem%nr, gt_act_elem%prob_acc
          ELSE
            tmc_env%m_env%estim_corr_wrong(3) = tmc_env%m_env%estim_corr_wrong(3) + 1
          END IF
        END IF
        gt_act_elem%stat = status_rejected_result
        gt_act_elem%prob_acc = 0.0_dp
        itmp = gt_act_elem%conf(gt_act_elem%mv_conf)%elem%sub_tree_nr

      !-- case ACCEPTED
      ELSE
        IF(gt_act_elem%prob_acc.GT.0.0_dp) THEN
          IF(gt_act_elem%prob_acc.LE.0.5_dp) THEN
            ! wrong estimate (estimated NOT accepted)
            tmc_env%m_env%estim_corr_wrong(2) = tmc_env%m_env%estim_corr_wrong(2) + 1
            IF(DEBUG.GT.0) WRITE(tmc_env%m_env%io_unit,*)&
              "wrong guess for ACC  (elem/estim acc prob)", gt_act_elem%nr, gt_act_elem%prob_acc
          ELSE
            tmc_env%m_env%estim_corr_wrong(1) = tmc_env%m_env%estim_corr_wrong(1) + 1
          END IF
        END IF
        gt_act_elem%stat = status_accepted_result
        gt_act_elem%prob_acc = 1.0_dp

        ! save result in the result list
        IF(.NOT.gt_act_elem%swaped) THEN
          !-- saving moved/changed configuration of result node
          tmc_env%m_env%result_list(gt_act_elem%mv_conf)%elem => &
                 gt_act_elem%conf(gt_act_elem%mv_conf)%elem
        ELSE
          ! ATTENTION: act_element != gt_act_elem%conf(gt_act_elem%mv_conf), 
          !   because we take the last accepted conf
          tmc_env%m_env%result_list(gt_act_elem%mv_conf)%elem => act_element
          tmc_env%m_env%result_list(gt_act_elem%mv_conf+1)%elem => tmp_element
        END IF
        tmc_env%m_env%gt_act => gt_act_elem
      END IF ! result acceptance check
      !======================================================================
      !-- set status accepted or rejected of certain subtree elements
      IF(.NOT.gt_act_elem%swaped) &
        CALL subtree_configuration_stat_change(gt_ptr=gt_act_elem, &
                                               stat=gt_act_elem%stat, &
                                               tmc_params=tmc_env%params, &
                                               error=error)

      IF(tmc_env%params%DRAW_TREE) &
        CALL create_global_tree_dot_color(gt_tree_element=gt_act_elem, &
                                          tmc_params=tmc_env%params, &
                                          error=error)

      ! probability update 
      CALL prob_update(move_types=tmc_env%params%move_types,&
                       pt_el=gt_act_elem, &
                       prob_opt=tmc_env%params%esimate_acc_prob, &
                       error=error)

      !writes only different configurations with repetition file if possible
      CALL write_result_list_element(result_list=tmc_env%m_env%result_list, &
             result_count=tmc_env%m_env%result_count,&
             conf_updated=gt_act_elem%mv_conf, accepted=result_acc, &
             tmc_params=tmc_env%params, error=error)
      IF(gt_act_elem%swaped)&
        CALL write_result_list_element(result_list=tmc_env%m_env%result_list, &
               result_count=tmc_env%m_env%result_count,&
               conf_updated=gt_act_elem%mv_conf+1, accepted=result_acc, &
               tmc_params=tmc_env%params, error=error)

      ! save for analysis
      IF(tmc_env%tmc_comp_set%para_env_m_ana%num_pe.GT.1.AND.result_acc) THEN
        CALL add_to_list(elem=tmc_env%m_env%result_list(gt_act_elem%mv_conf)%elem, &
                         list=tmc_env%m_env%analysis_list, &
                         temp_ind=gt_act_elem%mv_conf, &
                         nr=tmc_env%m_env%result_count(gt_act_elem%mv_conf),&
                         error=error)
        IF(gt_act_elem%swaped) THEN
          CALL add_to_list(elem=tmc_env%m_env%result_list(gt_act_elem%mv_conf+1)%elem, &
                           list=tmc_env%m_env%analysis_list, &
                           temp_ind=gt_act_elem%mv_conf+1, &
                           nr=tmc_env%m_env%result_count(gt_act_elem%mv_conf+1),&
                           error=error)
        END IF
      END IF
    END DO search_calculated_element_loop

    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE tree_update


  !============================================================================
  ! tree node acceptance probability
  !============================================================================
! *****************************************************************************
!> \brief checks if element is ready for accaptance probability update
!>        checks status and the amount of already received intermediate energies
!> \param tree_elem sub tree element to update
!> \param ready return value
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  FUNCTION ready_for_update_acc_prob(elem, error) RESULT(ready)
    TYPE(tree_type), POINTER                 :: elem
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: ready

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ready_for_update_acc_prob', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(elem),cp_failure_level,routineP,error,failure)
    ready = .FALSE.
    IF((elem%scf_energies_count.GE.4)&
       .AND.(elem%stat.NE.status_deleted).AND.(elem%stat.NE.status_deleted_result)&
       .AND.(elem%stat.NE.status_canceled_ener)) &
      ready = .TRUE.
  END FUNCTION ready_for_update_acc_prob

! *****************************************************************************
!> \brief updates the subtree acceptance probability
!>        the swap probabilities are handled within the certain checks of the 
!>        global tree elements (pt references)
!> \param tree_elem sub tree element to update
!> \param tmc_env TMC environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE check_elements_for_acc_prob_update(tree_elem, tmc_env, error)
    TYPE(tree_type), POINTER                 :: tree_elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'check_elements_for_acc_prob_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(tree_type), POINTER                 :: act_elem, parent_elem

    failure = .FALSE.
    NULLIFY(parent_elem, act_elem)

    CPPrecondition(ASSOCIATED(tree_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(tmc_env),cp_failure_level,routineP,error,failure)


    ! start the timing
    CALL timeset(routineN,handle)

    act_elem => tree_elem
    !-- nothing to do if option is disabled or element not ready
    IF(tmc_env%params%esimate_acc_prob.AND.&
       ready_for_update_acc_prob(act_elem,error=error))THEN
      ! ------------------------------------------------------
      ! check this element
      !   for usual moves and swapping
      ! get parent subtree elment for case of not swapped configurations
      parent_elem => search_parent_element(act_elem, error=error)
      ! update the prob of accaptance
      CALL update_prob_gt_node_list(reference_list=act_elem%gt_nodes_references, &
              act_elem=act_elem, parent_elem=parent_elem, act_parent=.TRUE., &
              tmc_env=tmc_env, error=error)

      !------------------------------------------------------
      ! check the childs of this element
      parent_elem => tree_elem

      ! check ACC child (parent element is the entered/updated element)
      IF(ASSOCIATED(tree_elem%acc)) THEN
         act_elem => tree_elem%acc
         ! update the prob of accaptance
         CALL update_prob_gt_node_list(reference_list=act_elem%gt_nodes_references, &
                 act_elem=act_elem, parent_elem=parent_elem, act_parent=.FALSE.,&
                 tmc_env=tmc_env, error=error)
      END IF

      ! check all NACC childs of next accepted one
      nacc_loop: DO
         IF(.NOT.ASSOCIATED(act_elem%nacc)) EXIT nacc_loop
         act_elem => act_elem%nacc
         CALL update_prob_gt_node_list(reference_list=act_elem%gt_nodes_references,&
                 act_elem=act_elem, parent_elem=parent_elem, act_parent=.FALSE.,&
                 tmc_env=tmc_env, error=error)
      END DO nacc_loop
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE check_elements_for_acc_prob_update

! *****************************************************************************
!> \brief search all pt nodes in reference list and update the estimated 
!>        acceptance probabilities
!> \param reference_list list of global tree element refering to the updated 
!>        subtree element
!> \param act_elem actual sub tree element updated or child of the updated one
!> \param parent_elem parent or the actual updated subtree element
!> \param act_parent flag if updated element is the actual (TRUE) or 
!>        the parent (FALSE) element
!> \param tmc_env TMC environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE update_prob_gt_node_list(reference_list, act_elem, parent_elem, act_parent, tmc_env, error)
    TYPE(gt_elem_list_type), POINTER         :: reference_list
    TYPE(tree_type), POINTER                 :: act_elem, parent_elem
    LOGICAL                                  :: act_parent
    TYPE(tmc_env_type), POINTER              :: tmc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_prob_gt_node_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, integration_steps
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kB, tmp_prob
    TYPE(gt_elem_list_type), POINTER         :: tmp_pt_ptr
    TYPE(tree_type), POINTER                 :: tmp_elem, tmp_parent_elem

    failure = .FALSE.
    kB = boltzmann/joule

    NULLIFY(tmp_pt_ptr,tmp_elem,tmp_parent_elem)

    IF(.NOT.ASSOCIATED(reference_list)) RETURN ! could be canceled already

    CPPrecondition(ASSOCIATED(reference_list%gt_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(act_elem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(parent_elem),cp_failure_level,routineP,error,failure)

    ! start the timing
    CALL timeset(routineN,handle)

    ! check if element is ready
    IF(ready_for_update_acc_prob(act_elem, error)) THEN
      ! set the number of integration steps used for 3 point approximation 
      ! of the exact energy, using the sub step energies
      ! still fixed value
      integration_steps = 100

      tmp_pt_ptr => reference_list
      ! do for all references using the updated subtree node
      reference_loop: DO WHILE(ASSOCIATED(tmp_pt_ptr))
        tmp_prob = -1.0_dp
        IF(tmp_pt_ptr%gt_elem%swaped)THEN
          NULLIFY(tmp_elem,tmp_parent_elem)
          ! in case of swap use the other swaped configuration as related one
          CALL get_subtree_elements_to_check(gt_act_elem=tmp_pt_ptr%gt_elem, &
                     elem1=tmp_elem, elem2=tmp_parent_elem, error=error)
          ! NOT if parent is the updated one, and check for correct elements ready
          IF(act_parent .AND. ASSOCIATED(act_elem, tmp_elem) .AND. &
             ready_for_update_acc_prob(elem=tmp_parent_elem, error=error)) THEN
            ! using ln(rnd)/-(beta_i-beta_j) < U_j-U_i)
            tmp_prob = compute_estimated_prob(elem_old=tmp_parent_elem, &
                                              elem_new=act_elem, &
                                              E_classical_diff=0.0_dp, &
                                              rnd_nr=tmp_pt_ptr%gt_elem%rnd_nr, &
                                              beta=1.0_dp/(kB*(tmc_env%params%Temp(tmp_pt_ptr%gt_elem%mv_conf)-&
                                                               tmc_env%params%Temp(tmp_pt_ptr%gt_elem%mv_conf+1))), &
                                              N=integration_steps, tmc_params=tmc_env%params,error=error)
          ELSE
            tmp_pt_ptr => tmp_pt_ptr%next
            CYCLE reference_loop
          END IF
        ELSE
          ! if no swap, use the parent configuration as related one
          IF(ready_for_update_acc_prob(parent_elem, error=error)) THEN
            tmp_prob = compute_estimated_prob(elem_old=parent_elem, &
                                              elem_new=act_elem, &
                                              E_classical_diff=act_elem%e_pot_approx-parent_elem%e_pot_approx, &
                                              rnd_nr=tmp_pt_ptr%gt_elem%rnd_nr, &
                                              beta=1.0_dp/(kB*tmc_env%params%Temp(tmp_pt_ptr%gt_elem%mv_conf)), &
                                              N=integration_steps, tmc_params=tmc_env%params,error=error)
          END IF
        END IF
        !successfull probability update
        IF(tmp_prob.GE.0.0_dp)THEN
          tmp_pt_ptr%gt_elem%prob_acc = tmp_prob
          !-- speculative canceling for the related direction
          IF(tmc_env%params%SPECULATIVE_CANCELING) &
            CALL search_canceling_elements(pt_elem_in=tmp_pt_ptr%gt_elem, &
                    prob=tmp_pt_ptr%gt_elem%prob_acc, tmc_env=tmc_env, error=error)
        END IF

        ! get next related global tree pointer
        tmp_pt_ptr => tmp_pt_ptr%next
      END DO reference_loop ! global tree references of actual subtree element
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE update_prob_gt_node_list

END MODULE tmc_tree_acceptance
