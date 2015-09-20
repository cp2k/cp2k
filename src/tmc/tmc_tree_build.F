!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief tree nodes creation, deallocation, references etc. 
!>        - we distinguish two kinds of tree nodes: global and sub tree nodes 
!>          (because we also are able to do parallel tempering)
!>        - global tree nodes consists of pointers to sub tree nodes
!>        - sub tree nodes consists of position arrays, potential energy, etc.
!>        - furthermore the sub tree elements have references the all global 
!>          tree elements refering to them
!>        - for tree element details see tree_types.F
!>
!>        - for creating we always start with the global tree element
!>          (if not already exist)
!>        - for each new global tree element (depending on the move type):
!>           - two sub tree elements are swapped (Parallel Tempering)
!>             (in global tree element creation)
!>           - the volume of a subtree element is changed 
!>             (directly in sub tree element creation)
!>           - positions in one subtree element changes 
!>             (in sub tree elem creation or NMC)
!>           - ...
!>        - sub tree elements will be deleted only if no reference to 
!>          any global tree element exist anymore
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_tree_build
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: get_rng_stream,&
                                             next_random_number,&
                                             reset_to_next_rng_substream,&
                                             set_rng_stream
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_calculations,                ONLY: calc_e_kin,&
                                             init_vel
  USE tmc_dot_tree,                    ONLY: create_dot,&
                                             create_dot_color,&
                                             create_global_tree_dot,&
                                             create_global_tree_dot_color
  USE tmc_file_io,                     ONLY: read_restart_file,&
                                             write_result_list_element
  USE tmc_move_handle,                 ONLY: select_random_move_type
  USE tmc_move_types,                  ONLY: &
       mv_type_MD, mv_type_NMC_moves, mv_type_atom_swap, mv_type_atom_trans, &
       mv_type_gausian_adapt, mv_type_mol_rot, mv_type_mol_trans, &
       mv_type_none, mv_type_proton_reorder, mv_type_swap_conf, &
       mv_type_volume_move
  USE tmc_moves,                       ONLY: change_pos,&
                                             elements_in_new_subbox
  USE tmc_stati,                       ONLY: TMC_STATUS_FAILED,&
                                             TMC_STATUS_WAIT_FOR_NEW_TASK,&
                                             task_type_MC,&
                                             task_type_gaussian_adaptation
  USE tmc_tree_references,             ONLY: add_to_references,&
                                             remove_gt_references,&
                                             remove_subtree_element_of_all_references,&
                                             search_and_remove_reference_in_list
  USE tmc_tree_search,                 ONLY: most_prob_end,&
                                             search_end_of_clean_g_tree,&
                                             search_end_of_clean_tree,&
                                             search_parent_element
  USE tmc_tree_types,                  ONLY: &
       add_to_list, elem_array_type, global_tree_type, status_accepted, &
       status_accepted_result, status_calc_approx_ener, status_calculate_MD, &
       status_calculate_NMC_steps, status_calculate_energy, &
       status_calculated, status_cancel_ener, status_cancel_nmc, &
       status_canceled_ener, status_canceled_nmc, status_created, &
       status_deleted, status_deleted_result, status_ok, status_rejected, &
       status_rejected_result, tree_type
  USE tmc_types,                       ONLY: tmc_env_type,&
                                             tmc_param_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_tree_build'

  PUBLIC :: allocate_new_sub_tree_node, deallocate_sub_tree_node
  PUBLIC :: init_tree_mod, finalize_init
  PUBLIC :: create_new_gt_tree_node
  PUBLIC :: remove_unused_g_tree
  PUBLIC :: remove_all_trees
  PUBLIC :: finalize_trees
CONTAINS


  !********************************************************************************
  ! ALLOCATION - DEALLOCATION
  !********************************************************************************
! *****************************************************************************
!> \brief allocates an elements of the global element structure 
!> \param next_el ...
!> \param nr_temp ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE allocate_new_global_tree_node(next_el, nr_temp)
    TYPE(global_tree_type), POINTER          :: next_el
    INTEGER                                  :: nr_temp

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'allocate_new_global_tree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, itmp
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(next_el))

    ! start the timing
    CALL timeset(routineN,handle)

    ! allocate everything
    ALLOCATE(next_el)
    ALLOCATE(next_el%conf(nr_temp))
    ALLOCATE(next_el%conf_n_acc(nr_temp))
    next_el%rnd_nr = -1.0_dp

    DO itmp = 1, nr_temp
      NULLIFY(next_el%conf(itmp)%elem)
      next_el%conf_n_acc(itmp) = .FALSE.
    END DO

    next_el%swaped = .FALSE.
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE allocate_new_global_tree_node

! *****************************************************************************
!> \brief deallocates an elements of the global element structure 
!> \param gt_elem ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE deallocate_global_tree_node(gt_elem)
    TYPE(global_tree_type), POINTER          :: gt_elem

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_global_tree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(gt_elem))

    ! start the timing
    CALL timeset(routineN,handle)

    ! deallocate everything
    DEALLOCATE(gt_elem%conf_n_acc)
    DEALLOCATE(gt_elem%conf)
    DEALLOCATE(gt_elem)
    NULLIFY(gt_elem)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE deallocate_global_tree_node  

! *****************************************************************************
!> \brief allocates an elements of the subtree element structure 
!> \param tmc_params structure for storing all (global) parameters
!> \param next_el ...
!> \param nr_dim ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE allocate_new_sub_tree_node(tmc_params, next_el, nr_dim)
    TYPE(tmc_param_type), POINTER            :: tmc_params
    TYPE(tree_type), POINTER                 :: next_el
    INTEGER                                  :: nr_dim

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_new_sub_tree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(next_el))

    ! start the timing
    CALL timeset(routineN,handle)

    ALLOCATE(next_el)
    NULLIFY(next_el%subbox_center, next_el%pos, next_el%mol, next_el%vel, &
            next_el%frc, next_el%dipole, next_el%elem_stat, &
            next_el%gt_nodes_references)

    next_el%scf_energies(:) = HUGE(next_el%scf_energies)
    next_el%scf_energies_count = 0
    ALLOCATE(next_el%pos(nr_dim))
    ALLOCATE(next_el%mol(nr_dim/tmc_params%dim_per_elem))
    ALLOCATE(next_el%vel(nr_dim))
    IF(tmc_params%print_dipole) ALLOCATE(next_el%dipole(tmc_params%dim_per_elem))
    ALLOCATE(next_el%elem_stat(nr_dim))
    next_el%elem_stat = status_ok
    ALLOCATE(next_el%subbox_center(tmc_params%dim_per_elem))
    IF(tmc_params%print_forces.OR.tmc_params%task_type.EQ.task_type_gaussian_adaptation)THEN
       IF(tmc_params%task_type.EQ.task_type_gaussian_adaptation) THEN
          ALLOCATE(next_el%frc(nr_dim*nr_dim))
       ELSE
          ALLOCATE(next_el%frc(nr_dim))
       END IF
       next_el%frc = 0.0_dp
    END IF
    ALLOCATE(next_el%box_scale(3))
    next_el%pos(:) = -1.0_dp
    next_el%mol(:) = -1
    next_el%box_scale(:) = 1.0_dp
    next_el%scf_energies(:) = 0.0_dp
    next_el%e_pot_approx = 0.0_dp
    next_el%potential = 76543.0_dp
    next_el%vel = 0.0_dp ! standart MC don"t uses velocities, but it is used at least in acceptance check
    next_el%ekin = 0.0_dp
    next_el%ekin_before_md = 0.0_dp
    next_el%sub_tree_nr = 0
    next_el%nr = -1
    next_el%rng_seed(:,:,:) = -1.0
    next_el%move_type = mv_type_none

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE allocate_new_sub_tree_node

! *****************************************************************************
!> \brief deallocates an elements of the subtree element structure 
!> \param tree_elem ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE deallocate_sub_tree_node(tree_elem)
    TYPE(tree_type), POINTER                 :: tree_elem

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_sub_tree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(tree_elem))

    ! start the timing
    CALL timeset(routineN,handle)

    ! reference handling
    ! should be not necessary, subtree element should be only deallocated,
    !   if no global tree element points to anymore
    CALL remove_subtree_element_of_all_references(ptr=tree_elem)

    IF(ASSOCIATED(tree_elem%box_scale)) DEALLOCATE(tree_elem%box_scale)
    IF(ASSOCIATED(tree_elem%frc)) DEALLOCATE(tree_elem%frc)
    IF(ASSOCIATED(tree_elem%subbox_center)) DEALLOCATE(tree_elem%subbox_center)
    IF(ASSOCIATED(tree_elem%elem_stat)) DEALLOCATE(tree_elem%elem_stat)
    IF(ASSOCIATED(tree_elem%dipole)) DEALLOCATE(tree_elem%dipole)
    IF(ASSOCIATED(tree_elem%vel)) DEALLOCATE(tree_elem%vel)
    IF(ASSOCIATED(tree_elem%mol)) DEALLOCATE(tree_elem%mol)
    IF(ASSOCIATED(tree_elem%pos)) DEALLOCATE(tree_elem%pos)

    DEALLOCATE(tree_elem)
    NULLIFY(tree_elem)
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE deallocate_sub_tree_node

  !********************************************************************************
  ! INITIALIZATION - FINALIZE
  !********************************************************************************

! *****************************************************************************
!> \brief routine initiate the global and subtrees with the first elements
!> \param start_elem ...
!> \param tmc_env structure for storing all (global) parameters
!> \param job_counts ...
!> \param worker_timings ...
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE init_tree_mod(start_elem, tmc_env, job_counts, worker_timings)
    TYPE(tree_type), POINTER                 :: start_elem
    TYPE(tmc_env_type), POINTER              :: tmc_env
    INTEGER, DIMENSION(:)                    :: job_counts
    REAL(KIND=dp), DIMENSION(4)              :: worker_timings

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_tree_mod', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    TYPE(global_tree_type), POINTER          :: global_tree

    NULLIFY(global_tree)

    CPASSERT(ASSOCIATED(start_elem))
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env))

    ! start the timing
    CALL timeset(routineN,handle)

    ! allocate everything
    CALL allocate_new_global_tree_node(next_el=tmc_env%m_env%gt_act, &
                                       nr_temp=tmc_env%params%nr_temp)

    ! use initial/default values
    CALL get_rng_stream(rng_stream=tmc_env%rng_stream, &
                        bg=tmc_env%m_env%gt_act%rng_seed(:,:,1), &
                        cg=tmc_env%m_env%gt_act%rng_seed(:,:,2), &
                        ig=tmc_env%m_env%gt_act%rng_seed(:,:,3))

    global_tree => tmc_env%m_env%gt_act
    tmc_env%m_env%gt_head => tmc_env%m_env%gt_act

    ! set global random seed
    CALL set_rng_stream(rng_stream=tmc_env%rng_stream, bg=global_tree%rng_seed(:,:,1), &
                        cg=global_tree%rng_seed(:,:,2), ig=global_tree%rng_seed(:,:,3))
    global_tree%rnd_nr = next_random_number(tmc_env%rng_stream)

    !-- SUBTREES: set initial values
    DO i=1, SIZE(global_tree%conf)
       CALL allocate_new_sub_tree_node(tmc_env%params, next_el=global_tree%conf(i)%elem, &
                                       nr_dim=SIZE(start_elem%pos))
       global_tree%conf(i)%elem%move_type     = 0
       global_tree%conf(i)%elem%next_elem_nr  => tmc_env%m_env%tree_node_count(i)
       global_tree%conf(i)%elem%parent        => NULL()
       global_tree%conf(i)%elem%nr            = global_tree%conf(i)%elem%next_elem_nr
       global_tree%conf(i)%elem%sub_tree_nr   = i
       global_tree%conf(i)%elem%elem_stat     = status_ok
       global_tree%conf(i)%elem%pos           = start_elem%pos
       global_tree%conf(i)%elem%mol           = start_elem%mol
       global_tree%conf(i)%elem%e_pot_approx  = start_elem%e_pot_approx
       global_tree%conf(i)%elem%temp_created  = i
       global_tree%conf(i)%elem%stat          = status_calculate_energy
       !it is default already: global_tree%conf(i)%elem%box_scale(:)  = 1.0_dp
       IF(tmc_env%params%task_type.EQ.task_type_gaussian_adaptation)THEN
          global_tree%conf(i)%elem%vel(:)         = start_elem%vel(:)
          global_tree%conf(i)%elem%frc(:)         = start_elem%frc(:)
          global_tree%conf(i)%elem%potential      = start_elem%potential
          global_tree%conf(i)%elem%ekin           = start_elem%ekin
          global_tree%conf(i)%elem%ekin_before_md = start_elem%ekin_before_md
       END IF

       !-- different random seeds for every subtree
       CALL reset_to_next_rng_substream(tmc_env%rng_stream)
       CALL get_rng_stream(rng_stream=tmc_env%rng_stream, bg=global_tree%conf(i)%elem%rng_seed(:,:,1), &
                           cg=global_tree%conf(i)%elem%rng_seed(:,:,2), &
                           ig=global_tree%conf(i)%elem%rng_seed(:,:,3))

       !-- gaussian distributed velocities
       !-- calculating the kinetic energy of the initial configuration velocity
       IF(tmc_env%params%task_type.EQ.task_type_MC)THEN
          IF(tmc_env%params%move_types%mv_weight(mv_type_MD).GT.0.0_dp) THEN
             CALL init_vel(vel=global_tree%conf(i)%elem%vel, atoms=tmc_env%params%atoms, &
                           temerature=tmc_env%params%Temp(i), &
                           rng_stream=tmc_env%rng_stream, &
                           rnd_seed=global_tree%conf(i)%elem%rng_seed)
             global_tree%conf(i)%elem%ekin = calc_e_kin(vel=global_tree%conf(i)%elem%vel, &
                                                        atoms=tmc_env%params%atoms)
          END IF
       END IF

       !-- set tree pointer
       !-- set pointer of first global tree element
       tmc_env%m_env%st_heads(i)%elem         => global_tree%conf(i)%elem
       tmc_env%m_env%st_clean_ends(i)%elem    => global_tree%conf(i)%elem
       !-- set initial pointer of result lists
       tmc_env%m_env%result_list(i)%elem      => global_tree%conf(i)%elem
    END DO
    tmc_env%m_env%tree_node_count(:) = 0 ! initializing the tree node numbering

    !-- initial global tree element
    tmc_env%m_env%gt_head      => global_tree
    tmc_env%m_env%gt_clean_end => global_tree
    global_tree%nr = 0
    global_tree%swaped = .FALSE.
    global_tree%mv_conf = 1
    global_tree%mv_next_conf = MODULO(global_tree%mv_conf,SIZE(global_tree%conf))+1
    global_tree%conf_n_acc = .TRUE.

    global_tree%stat = status_created
    global_tree%prob_acc = 1.0_dp

    ! simulated annealing start temperature
    global_tree%Temp = tmc_env%params%Temp(1)
    CALL cp_assert(.NOT.(tmc_env%params%nr_temp.NE.1 .AND. tmc_env%m_env%temp_decrease.NE.1.0_dp),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "there is no parallel tempering implementation for simulated annealing implemented "//&
                     "(just one Temp per global tree element.",&
                     only_ionode=.TRUE.)

    !-- IF program is restarted, read restart file
    IF(tmc_env%m_env%restart_in_file_name.NE."") THEN
      CALL read_restart_file(tmc_env=tmc_env, job_counts=job_counts, &
                             timings=worker_timings, &
                             file_name=tmc_env%m_env%restart_in_file_name)

      tmc_env%m_env%tree_node_count(0) = global_tree%nr

      DO i=1, SIZE(tmc_env%m_env%result_list(:))
        tmc_env%m_env%tree_node_count(i) = tmc_env%m_env%result_list(i)%elem%nr
        global_tree%conf(i)%elem%stat    = status_accepted
      END DO
      global_tree%prob_acc = 1.0_dp ! accepted (re)start configuration
      WRITE(tmc_env%m_env%io_unit,*)"TMC| restarting at Markov Chain element(s): ",&
                                     tmc_env%m_env%result_count
      !TODO enable calculation of the approx energy for case of fitting potential
      !     and changing the potential in between 
      !     BUT check, there is no double counting (of the last/restarted elem) in the trajectory
      !IF(tmc_env%params%NMC_inp_file.NE."") &
      !  global_tree%conf(1)%elem%stat  = status_calc_approx_ener
      global_tree%stat                 = status_accepted_result
    ELSE IF(tmc_env%params%NMC_inp_file.NE."") THEN
      global_tree%conf(1)%elem%stat    = status_calc_approx_ener
    ELSE
      global_tree%conf(1)%elem%stat    = status_created
    END IF

    !-- set reference of global tree node 
    CALL add_to_references(gt_elem=global_tree)

    !-- draw the first global tree node
    IF(tmc_env%params%DRAW_TREE)THEN
      CALL create_global_tree_dot(new_element=global_tree, &
                                  tmc_params=tmc_env%params)
      CALL create_global_tree_dot_color(gt_tree_element=global_tree, &
                                  tmc_params=tmc_env%params)
    END IF

    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE init_tree_mod

! *****************************************************************************
!> \brief distributes the initial energy to all subtree (if no restart) and 
!>        call analysis for this element (write trajectory...)
!> \param gt_tree_ptr global tree head (initial configuration)
!> \param tmc_env master environment for restart 
!>        (if restart the subtree heads are not equal), result counts and lists
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE finalize_init(gt_tree_ptr, tmc_env)
    TYPE(global_tree_type), POINTER          :: gt_tree_ptr
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(gt_tree_ptr))
    CPASSERT(.NOT.ASSOCIATED(gt_tree_ptr%parent))
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env))
    CPASSERT(ASSOCIATED(tmc_env%params))

    ! start the timing
    CALL timeset(routineN,handle)

    gt_tree_ptr%stat = status_accepted_result
    !-- distribute energy of first element to all subtrees
    DO i=1,SIZE(gt_tree_ptr%conf)
       gt_tree_ptr%conf(i)%elem%stat      = status_accepted_result
       IF(ASSOCIATED(gt_tree_ptr%conf(1)%elem%dipole)) &
         gt_tree_ptr%conf(i)%elem%dipole = gt_tree_ptr%conf(1)%elem%dipole
       IF(tmc_env%m_env%restart_in_file_name.EQ."") &
         gt_tree_ptr%conf(i)%elem%potential = gt_tree_ptr%conf(1)%elem%potential
    END DO

    IF(tmc_env%m_env%restart_in_file_name.EQ."") THEN
      tmc_env%m_env%result_count(:) = tmc_env%m_env%result_count(:) +1
      tmc_env%m_env%result_list(:)  = gt_tree_ptr%conf(:)
      !-- write initial elements in result files
      DO i=1, SIZE(tmc_env%m_env%result_list)
        CALL write_result_list_element(result_list=tmc_env%m_env%result_list, &
               result_count=tmc_env%m_env%result_count,&
               conf_updated=i, accepted=.TRUE., &
               tmc_params=tmc_env%params)
        ! save for analysis
        IF(tmc_env%tmc_comp_set%para_env_m_ana%num_pe.GT.1) THEN
          CALL add_to_list(elem=tmc_env%m_env%result_list(i)%elem, &
                           list=tmc_env%m_env%analysis_list, &
                           nr=tmc_env%m_env%result_count(i), & 
                           temp_ind=i)
        END IF
      END DO
      !CALL write_result_list_element(result_list=tmc_env%m_env%result_list, &
      !         result_count=tmc_env%m_env%result_count,&
      !         conf_updated=0, accepted=.TRUE., &
      !         tmc_params=tmc_env%params)
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE finalize_init

  !============================================================================
  ! tree node creation
  !============================================================================
! *****************************************************************************
!> \brief creates new global tree element and if needed new subtree element
!> \param tmc_env TMC environment with parameters and pointers to gt element
!> \param stat return status value
!> \param new_elem return gt element
!> \param reactivation_cc_count counting the reactivation of subtree elements
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_new_gt_tree_node(tmc_env, stat, new_elem, &
                                     reactivation_cc_count)
    TYPE(tmc_env_type), POINTER              :: tmc_env
    INTEGER, INTENT(OUT)                     :: stat
    TYPE(global_tree_type), INTENT(OUT), &
      POINTER                                :: new_elem
    INTEGER                                  :: reactivation_cc_count

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_new_gt_tree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, swap_conf
    LOGICAL                                  :: failure, keep_on, n_acc
    REAL(KIND=dp)                            :: prob, rnd, rnd2
    TYPE(global_tree_type), POINTER          :: tmp_elem
    TYPE(tree_type), POINTER                 :: tree_elem

    failure = .FALSE.
    NULLIFY(tmp_elem, tree_elem, new_elem)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%params))
    CPASSERT(ASSOCIATED(tmc_env%m_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env%gt_act))

    ! start the timing
    CALL timeset(routineN,handle)

    stat = TMC_STATUS_FAILED
    !-- search most probable end in global tree for new element
    tmp_elem => tmc_env%m_env%gt_act
    n_acc = .TRUE.

    !-- search most probable end to create new element
    CALL most_prob_end(global_tree_elem=tmp_elem, prob=prob, n_acc=n_acc)

    keep_on = .TRUE.
    IF(ASSOCIATED(tmp_elem).AND.(EXP(prob).LT.1.0E-10)) THEN
       new_elem => NULL()
       stat = TMC_STATUS_FAILED
       keep_on = .FALSE.
    !-- if not found, do something else 
    !-- (posible if just one end for further calculations 
    !    and there a MD move is still calculated)
    ELSE IF(.NOT.ASSOCIATED(tmp_elem))THEN
       new_elem => NULL()
       stat = TMC_STATUS_FAILED
       keep_on = .FALSE.
    END IF

    IF(keep_on) THEN
      ! if global tree element already exist use that one 
      !   (skip creating new element)
      ! reactivation
      IF((n_acc.AND.ASSOCIATED(tmp_elem%acc)).OR.&
         ((.NOT.n_acc).AND.ASSOCIATED(tmp_elem%nacc)))THEN

        !set pointer to the actual element
        IF(n_acc)&
           new_elem => tmp_elem%acc
        IF(.NOT.n_acc)&
           new_elem => tmp_elem%nacc

        ! check for existing subtree element
        CPASSERT(ASSOCIATED(new_elem%conf(new_elem%mv_conf)%elem))
        SELECT CASE(new_elem%conf(new_elem%mv_conf)%elem%stat)
        CASE(status_cancel_nmc,status_cancel_ener,status_canceled_nmc,&
             status_canceled_ener)
          ! reactivating subtree element 
          !  (but global tree element already exist)
          CALL add_to_references(gt_elem=new_elem)
          reactivation_cc_count = reactivation_cc_count+1
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"global tree node creation using existing sub tree element, "//&
              "but is not a canceled one, gt elem "//&
              cp_to_string(new_elem%nr)//" st elem "//&
              cp_to_string(new_elem%conf(new_elem%mv_conf)%elem%nr)//&
              " with stat "//&
              cp_to_string(new_elem%conf(new_elem%mv_conf)%elem%stat),&
              failure)
        END SELECT
        ! change the status of the reactivated subtree element
        ! move is only done by the master, 
        !  when standard MC moves with single potential are done
        ! the Nested Monte Carlo routine needs to do the configuration 
        !  to have old configuration to see if change is accepted
        SELECT CASE (new_elem%conf(new_elem%mv_conf)%elem%move_type)
        CASE(mv_type_MD)
          new_elem%conf(new_elem%mv_conf)%elem%stat = status_calculate_MD
        CASE(mv_type_NMC_moves)
          CALL cp_assert(new_elem%conf(new_elem%mv_conf)%elem%stat.EQ.&
              status_canceled_nmc, cp_warning_level,cp_assertion_failed,&
              routineP,"reactivating tree element with wrong status"//&
              cp_to_string(new_elem%conf(new_elem%mv_conf)%elem%stat),& 
              failure)
          new_elem%conf(new_elem%mv_conf)%elem%stat = status_calculate_NMC_steps
          
          !IF(DEBUG.GE.1) WRITE(tmc_out_file_nr,*)"ATTENTION: reactivation of canceled subtree ", &
          !  new_elem%conf(new_elem%mv_conf)%elem%sub_tree_nr, "elem", new_elem%conf(new_elem%mv_conf)%elem%nr, &
          !  " of existing gt elem ",new_elem%nr,", again calculate NMC steps"
        CASE(mv_type_atom_trans, mv_type_mol_trans, mv_type_mol_rot, &
             mv_type_proton_reorder)
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"reactivated st element has no NMC or MD move type, "//&
              "but seems to be canceled. Move type"//&
              cp_to_string(new_elem%conf(new_elem%mv_conf)%elem%move_type), &
              failure)
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"Unknown move type while reactivating subtree element.",&
              failure)
        END SELECT
      ELSE
        !-- if end is found (NOT already existing element), create new elem at the end and if nessecarry new subtree element
        ! set initial values
        CALL allocate_new_global_tree_node(next_el=new_elem, &
                                           nr_temp=tmc_env%params%nr_temp)
        tmc_env%m_env%tree_node_count(0) = tmc_env%m_env%tree_node_count(0) +1
        new_elem%nr = tmc_env%m_env%tree_node_count(0)

        !-- set pointers to and from element one level up
        !-- paste new gt tree node element at right end
        IF(n_acc)THEN
          CALL cp_assert(.NOT.ASSOCIATED(tmp_elem%acc),cp_failure_level,&
              cp_assertion_failed, routineP,&
              "creating new subtree element on an occupied acc branch",&
              failure)
          tmp_elem%acc => new_elem
        ELSE
          CALL cp_assert(.NOT.ASSOCIATED(tmp_elem%nacc),cp_failure_level,&
              cp_assertion_failed, routineP,&
              "creating new subtree element on an occupied nacc branch", &
              failure)
          tmp_elem%nacc => new_elem
        END IF
        new_elem%parent   => tmp_elem

        !-- adopt acceptance flags of elements (old)
        new_elem%conf_n_acc(:) = new_elem%parent%conf_n_acc
        !-- set acceptance flag of modified configuration 
        !    depending on the direction of attaching new element
        IF(.NOT.new_elem%parent%swaped)THEN
          ! set the flag for the direction 
          !  (shows if the configuration is assumed to be acc or rej)
          new_elem%conf_n_acc(new_elem%parent%conf(&
             new_elem%parent%mv_conf)%elem%sub_tree_nr) = n_acc
        ELSE
          !-- in case of swapping the subtree element acceptance do not change
          !-- in case of NOT accepted branch and swapping before, 
          !-- search last NOT swaped gt tree node to take configurations
          IF(.NOT.n_acc) THEN
            DO
              IF(.NOT.ASSOCIATED(tmp_elem%parent)) EXIT
              IF(ASSOCIATED(tmp_elem%parent%acc,tmp_elem))THEN
                 tmp_elem => tmp_elem%parent
                 EXIT
              END IF
              tmp_elem => tmp_elem%parent
              IF(.NOT.tmp_elem%swaped) EXIT
            END DO
          END IF
        END IF

        !-- adapt "old" configurations
        new_elem%conf(:) = tmp_elem%conf(:)

        !-- set rnd nr generator and set next conf to change
        CALL set_rng_stream(rng_stream=tmc_env%rng_stream, &
                            bg=new_elem%parent%rng_seed(:,:,1), &
                            cg=new_elem%parent%rng_seed(:,:,2), &
                            ig=new_elem%parent%rng_seed(:,:,3))
        CALL reset_to_next_rng_substream(tmc_env%rng_stream)
        ! the random number for acceptance check
        new_elem%rnd_nr = next_random_number(tmc_env%rng_stream) 

        ! the next configuration index to move
        !rnd = next_random_number(tmc_env%rng_stream)
        !new_elem%mv_conf = 1+INT(size(new_elem%conf)*rnd)
        ! one temperature after each other
        new_elem%mv_conf      = new_elem%parent%mv_next_conf
        new_elem%mv_next_conf = MODULO(new_elem%mv_conf,SIZE(new_elem%conf))+1

        ! simulated annealing temperature decrease
        new_elem%Temp = tmp_elem%Temp
        IF(n_acc) new_elem%Temp = tmp_elem%Temp * (1-tmc_env%m_env%temp_decrease)

        !-- rnd for swap
        rnd = next_random_number(tmc_env%rng_stream)
        rnd2= next_random_number(tmc_env%rng_stream)
        CALL get_rng_stream(rng_stream=tmc_env%rng_stream, &
                            bg=new_elem%rng_seed(:,:,1), &
                            cg=new_elem%rng_seed(:,:,2), &
                            ig=new_elem%rng_seed(:,:,3))

        ! swap moves are not part of the subtree structure, 
        !  because exisiting elements from DIFFERENT subtrees are swaped
        ! -- do swap ?!
        IF(tmc_env%params%move_types%mv_weight(mv_type_swap_conf).GE.rnd)THEN
          ! set the index for the swaping element 
          !  and the conf to move in next move
          new_elem%mv_next_conf = new_elem%mv_conf
          ! do swap with conf swap_conf and swap_conf+1
          swap_conf = 1+INT((tmc_env%params%nr_temp-1)*rnd2)
          new_elem%mv_conf = swap_conf
          !-- swaping pointers to subtree elements
          ! exchange the pointer to the sub tree elements
          tree_elem => new_elem%conf(new_elem%mv_conf)%elem
          new_elem%conf(new_elem%mv_conf)%elem => &
             new_elem%conf(new_elem%mv_conf+1)%elem
          new_elem%conf(new_elem%mv_conf+1)%elem => tree_elem

          new_elem%stat = status_calculated
          new_elem%swaped = .TRUE.
          new_elem%prob_acc = tmc_env%params%move_types%acc_prob(&
                                 mv_type_swap_conf,new_elem%mv_conf)
          CALL add_to_references(gt_elem=new_elem)
          IF(tmc_env%params%DRAW_TREE) &
            CALL create_global_tree_dot(new_element=new_elem, &
                   tmc_params=tmc_env%params)
          ! nothing to do for the workers
          stat = status_calculated
          keep_on = .FALSE.
        ELSE

          !-- considered subtree node can already exist,
          !    calculated somewhere else in the global tree
          !-- so check if new sub tree node exists, if not, create it
          !-- check if considered configuration is assumed to be 
          !    on accepted or rejected branch
          IF(new_elem%conf_n_acc(new_elem%conf(new_elem%mv_conf)%elem%sub_tree_nr))THEN
            !-- check if child element in ACCEPTED direction already exist
            IF(ASSOCIATED(new_elem%conf(new_elem%mv_conf)%elem%acc))THEN
              new_elem%conf(new_elem%mv_conf)%elem => &
                new_elem%conf(new_elem%mv_conf)%elem%acc
              stat = status_calculated
            ELSE
              !-- if not exist create new subtree element
              CALL create_new_subtree_node(act_gt_el=new_elem, &
                                           tmc_env=tmc_env)
              IF(tmc_env%params%DRAW_TREE) &
                CALL create_dot(new_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                conf=new_elem%conf(new_elem%mv_conf)%elem%sub_tree_nr,&
                                tmc_params=tmc_env%params)
            END IF
          ELSE
            !-- check if child element in REJECTED direction already exist
            IF(ASSOCIATED(new_elem%conf(new_elem%mv_conf)%elem%nacc))THEN
              new_elem%conf(new_elem%mv_conf)%elem => &
                new_elem%conf(new_elem%mv_conf)%elem%nacc
              stat = status_calculated
            ELSE
              !-- if not exist create new subtree element
              CALL create_new_subtree_node(act_gt_el=new_elem, &
                                           tmc_env=tmc_env)
              IF(tmc_env%params%DRAW_TREE) &
                CALL create_dot(new_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                conf=new_elem%conf(new_elem%mv_conf)%elem%sub_tree_nr,&
                                tmc_params=tmc_env%params)
            END IF
          END IF
          ! set approximate probability of acceptance 
          !  (initialization with calculated values from 
          !  (#acc elem in traj)/(#elem in traj))
          new_elem%prob_acc = tmc_env%params%move_types%acc_prob(&
             new_elem%conf(new_elem%mv_conf)%elem%move_type, new_elem%mv_conf)
          ! add refence and dot
          CALL add_to_references(gt_elem=new_elem)
          IF(tmc_env%params%DRAW_TREE) &
            CALL create_global_tree_dot(new_element=new_elem, &
                                        tmc_params=tmc_env%params)
        END IF ! swap or no swap
      END IF ! global tree node already exist. Hence the Subtree node also (it is speculative canceled)
    END IF ! keep on (checking and creating)


    IF(keep_on) THEN ! status changes
      CALL cp_assert(new_elem%stat.NE.status_accepted_result .AND. &
                     new_elem%stat.NE.status_accepted .AND. &
                     new_elem%stat.NE.status_rejected .AND. &
                     new_elem%stat.NE.status_rejected_result, &
                     cp_failure_level, cp_assertion_failed, routineP,&
                     "selected existing RESULT gt node",&
                     failure)
      !-- set status of global tree element for decision in master routine
      SELECT CASE(new_elem%conf(new_elem%mv_conf)%elem%stat)
      CASE(status_rejected_result,status_rejected,status_accepted,&
           status_accepted_result,status_calculated)
        ! energy is already calculated 
        new_elem%stat = status_calculated
        stat = new_elem%conf(new_elem%mv_conf)%elem%stat
        IF(tmc_env%params%DRAW_TREE) &
          CALL create_dot_color(tree_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                tmc_params=tmc_env%params)
      CASE(status_calc_approx_ener)
        new_elem%stat = new_elem%conf(new_elem%mv_conf)%elem%stat
        IF(stat.NE.status_calculated)THEN
          stat = new_elem%conf(new_elem%mv_conf)%elem%stat
          IF(tmc_env%params%DRAW_TREE) &
            CALL create_dot_color(tree_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                  tmc_params=tmc_env%params)
        END IF
      CASE(status_calculate_MD,status_calculate_energy,&
           status_calculate_NMC_steps, status_created)
        ! if not already in progress, set status for new task message
        new_elem%stat = new_elem%conf(new_elem%mv_conf)%elem%stat
        IF(stat.NE.status_calculated)THEN
          stat = new_elem%conf(new_elem%mv_conf)%elem%stat
          IF(tmc_env%params%DRAW_TREE) &
            CALL create_dot_color(tree_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                  tmc_params=tmc_env%params)
        END IF
      CASE(status_cancel_ener,status_canceled_ener)
         ! configuration is already created, 
         !  but energy has to be calculated (again)
         new_elem%conf(new_elem%mv_conf)%elem%stat = status_created
         new_elem%stat = status_created
         ! creation complete, handle energy calculation at a different position
         !  (for different worker group)
         stat = status_calculated
         IF(tmc_env%params%DRAW_TREE) &
           CALL create_dot_color(tree_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                 tmc_params=tmc_env%params)
      CASE(status_cancel_nmc, status_canceled_nmc)
         ! reactivation canceled element (but with new global tree element)
         new_elem%conf(new_elem%mv_conf)%elem%stat = &
            status_calculate_NMC_steps
         new_elem%stat = status_calculate_NMC_steps
         stat = new_elem%conf(new_elem%mv_conf)%elem%stat
         reactivation_cc_count = reactivation_cc_count+1
         IF(tmc_env%params%DRAW_TREE) &
           CALL create_dot_color(tree_element=new_elem%conf(new_elem%mv_conf)%elem,&
                                 tmc_params=tmc_env%params)
      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                       routineP, "unknown stat "//&
                       cp_to_string(new_elem%conf(new_elem%mv_conf)%elem%stat)//&
                       "of subtree element "//&
                       "for creating new gt element", &
                       failure)
      END SELECT

      ! set stat TMC_STATUS_WAIT_FOR_NEW_TASK if no new calculation neccessary 
      !   (energy calculation nodes searched by different routine)
      IF(stat.EQ.TMC_STATUS_FAILED) stat = TMC_STATUS_WAIT_FOR_NEW_TASK
      IF(stat.EQ.status_calculated) stat = TMC_STATUS_WAIT_FOR_NEW_TASK
    END IF
    ! end the timing
    CALL timestop(handle)

  END SUBROUTINE create_new_gt_tree_node

! *****************************************************************************
!> \brief create new subtree element using pointer of global tree
!> \param act_gt_el global tree element 
!> \param tmc_env ...
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE create_new_subtree_node(act_gt_el, tmc_env)
    TYPE(global_tree_type), POINTER          :: act_gt_el
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_new_subtree_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: conf, handle, itmp
    LOGICAL                                  :: failure, mv_rejected, &
                                                new_subbox
    REAL(KIND=dp)                            :: rnd
    TYPE(tree_type), POINTER                 :: new_elem, parent_elem

    failure = .FALSE.
    NULLIFY(new_elem, parent_elem)

    CPASSERT(ASSOCIATED(act_gt_el))
    CPASSERT(ASSOCIATED(act_gt_el%conf(act_gt_el%mv_conf)%elem))
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%params))

    ! start the timing
    CALL timeset(routineN,handle)

    conf = act_gt_el%mv_conf
    CALL allocate_new_sub_tree_node(tmc_params=tmc_env%params,&
           next_el=new_elem,nr_dim=SIZE(act_gt_el%parent%conf(conf)%elem%pos))

    !-- node one level up
    parent_elem     => act_gt_el%conf(conf)%elem
    new_elem%parent => parent_elem

    !-- set initial values
    parent_elem%next_elem_nr = parent_elem%next_elem_nr+1
    new_elem%nr = parent_elem%next_elem_nr
    new_elem%rng_seed = parent_elem%rng_seed

    !-- change to real parent element
    IF(act_gt_el%conf_n_acc(act_gt_el%conf(act_gt_el%mv_conf)%elem%sub_tree_nr))THEN
      parent_elem%acc  => new_elem
    ELSE
      parent_elem%nacc => new_elem
    END IF

    !-- real parent node (taking the configuration from)
    ! search parent
    parent_elem => search_parent_element(current=new_elem)
    new_elem%pos(:)       = parent_elem%pos(:)
    new_elem%mol(:)       = parent_elem%mol(:)
    new_elem%vel(:)       = parent_elem%vel(:)
    new_elem%ekin         = parent_elem%ekin
    new_elem%e_pot_approx = parent_elem%e_pot_approx
    new_elem%next_elem_nr => parent_elem%next_elem_nr
    new_elem%sub_tree_nr  = parent_elem%sub_tree_nr
    new_elem%box_scale    = parent_elem%box_scale
    IF(tmc_env%params%task_type.EQ.task_type_gaussian_adaptation)THEN
      new_elem%frc(:)         = parent_elem%frc(:)
      new_elem%potential      = parent_elem%potential
      new_elem%ekin_before_md = parent_elem%ekin_before_md
    ELSE
      new_elem%potential    = 97589.0_dp
    END IF

    ! set new substream of random number generator
    CALL set_rng_stream(rng_stream=tmc_env%rng_stream, &
                        bg=new_elem%rng_seed(:,:,1), &
                        cg=new_elem%rng_seed(:,:,2), &
                        ig=new_elem%rng_seed(:,:,3))
    CALL reset_to_next_rng_substream(tmc_env%rng_stream)

    ! set the temperature for the NMC moves
    rnd =next_random_number(tmc_env%rng_stream)
    IF(tmc_env%params%NMC_inp_file.NE."") THEN
      new_elem%temp_created = INT(tmc_env%params%nr_temp*rnd)+1
    ELSE
      new_elem%temp_created = act_gt_el%mv_conf
    END IF

    ! rnd nr for selecting move
    rnd= next_random_number(tmc_env%rng_stream)
    !-- set move type
    new_elem%move_type = select_random_move_type(&
                           move_types=tmc_env%params%move_types, &
                           rnd=rnd)
    CALL get_rng_stream(rng_stream=tmc_env%rng_stream, &
                        bg=new_elem%rng_seed(:,:,1), &
                        cg=new_elem%rng_seed(:,:,2), &
                        ig=new_elem%rng_seed(:,:,3))

    ! move is only done by the master, 
    !  when standard MC moves with single potential are done
    ! the Nested Monte Carlo routine needs the old configuration 
    !  to see if change is accepted
    SELECT CASE (new_elem%move_type)
    CASE(mv_type_MD)
      ! velocity change have to be done on workers, 
      !  because of velocity change for NMC acceptance check
      new_elem%stat = status_calculate_MD
      ! set the temperature for creating MD 
      new_elem%temp_created = act_gt_el%mv_conf
      !-- set the subbox (elements in subbox)
      CALL elements_in_new_subbox(tmc_params=tmc_env%params, &
             rng_stream=tmc_env%rng_stream, elem=new_elem, &
             nr_of_sub_box_elements=itmp)
      ! the move is performed on a worker group
    CASE(mv_type_NMC_moves)
      new_elem%stat = status_calculate_NMC_steps
      !-- set the subbox (elements in subbox)
      CALL elements_in_new_subbox(tmc_params=tmc_env%params, &
             rng_stream=tmc_env%rng_stream, elem=new_elem, &
             nr_of_sub_box_elements=itmp)
       ! the move is performed on a worker group
    ! the following moves new no force_env and can be performed on the master directly
    CASE(mv_type_atom_trans, mv_type_atom_swap, mv_type_mol_trans, &
         mv_type_mol_rot, mv_type_proton_reorder, &
         mv_type_volume_move)
        new_subbox = .TRUE.
      ! volume move on whole cell
      IF(new_elem%move_type.EQ.mv_type_volume_move) THEN
        new_subbox = .FALSE.
      END IF
      CALL change_pos(tmc_params=tmc_env%params, &
             move_types=tmc_env%params%move_types, &
             rng_stream=tmc_env%rng_stream, elem=new_elem, &
             mv_conf=conf, new_subbox=new_subbox, &
             move_rejected=mv_rejected)
      IF(mv_rejected) THEN
        new_elem%potential = HUGE(new_elem%potential)
        new_elem%e_pot_approx = HUGE(new_elem%e_pot_approx)
        new_elem%stat = status_calculated
      ELSE
        new_elem%stat = status_created
        IF(tmc_env%params%NMC_inp_file.NE."") &
          new_elem%stat = status_calc_approx_ener
      END IF
    CASE(mv_type_gausian_adapt)
      ! still could be implemented
    CASE DEFAULT
      CALL cp_assert(.FALSE.,&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "unknown move type ("//cp_to_string(new_elem%move_type)//&
                  "), while creating subtree element.")
    END SELECT
    act_gt_el%conf(act_gt_el%mv_conf)%elem => new_elem

    ! end the timing
    CALL timestop(handle)
    CPASSERT(ASSOCIATED(act_gt_el%conf(act_gt_el%mv_conf)%elem))
  END SUBROUTINE create_new_subtree_node

  !============================================================================
  ! tree node deallocation
  !============================================================================
! *****************************************************************************
!> \brief prepares for deallocation of global tree element 
!>        (checks status and set pointers of neighboring elements)
!> \param gt_ptr the global tree element 
!> \param draw if present, changes the coleor in the dot file
!> \param tmc_env tmc environment
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE remove_gt_elem(gt_ptr, draw, tmc_env)
    TYPE(global_tree_type), POINTER          :: gt_ptr
    LOGICAL, OPTIONAL                        :: draw
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_gt_elem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPASSERT(ASSOCIATED(gt_ptr))
    CPASSERT(ASSOCIATED(tmc_env))

    ! start the timing
    CALL timeset(routineN,handle)

    CALL remove_gt_references(gt_ptr=gt_ptr, tmc_env=tmc_env)

    ! set status and draw in tree
    IF((gt_ptr%stat.EQ.status_accepted_result).OR.(gt_ptr%stat.EQ.status_rejected_result))THEN
       gt_ptr%stat = status_deleted_result
    ELSE
       gt_ptr%stat = status_deleted
    END IF
    IF(tmc_env%params%DRAW_TREE .AND. PRESENT(draw)) &
      CALL create_global_tree_dot_color(gt_tree_element=gt_ptr, tmc_params=tmc_env%params)

    !remove pointer from tree parent
    IF(ASSOCIATED(gt_ptr%parent)) THEN
       IF(ASSOCIATED(gt_ptr%parent%acc, gt_ptr))THEN
          gt_ptr%parent%acc => NULL()
       END IF
       IF(ASSOCIATED(gt_ptr%parent%nacc, gt_ptr))THEN
          gt_ptr%parent%nacc => NULL()
       END IF
    END IF

    !remove pointer from tree childs
    IF(ASSOCIATED(gt_ptr%acc))THEN
       gt_ptr%acc%parent => NULL()
    END IF

    IF(ASSOCIATED(gt_ptr%nacc))THEN
       gt_ptr%nacc%parent => NULL()
    END IF

    CALL deallocate_global_tree_node(gt_elem=gt_ptr)
    ! end the timing
    CALL timestop(handle)

    CPASSERT(.NOT.ASSOCIATED(gt_ptr))
  END SUBROUTINE remove_gt_elem

! *****************************************************************************
!> \brief prepares for deallocation of sub tree element 
!>        (checks status and set pointers of neighboring elements)
!> \param ptr the sub tree element 
!> \param draw if present, changes the coleor in the dot file
!> \param tmc_env tmc environment
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE remove_st_elem(ptr, draw, tmc_env)
    TYPE(tree_type), POINTER                 :: ptr
    LOGICAL, OPTIONAL                        :: draw
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_st_elem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, ready

    failure = .FALSE.
    ready = .TRUE.
    CPASSERT(ASSOCIATED(ptr))
    CPASSERT(ASSOCIATED(tmc_env))

    ! start the timing
    CALL timeset(routineN,handle)

    ! if there is still e reference to a global tree pointer, do not deallocate element
    IF(ASSOCIATED(ptr%gt_nodes_references)) THEN
      CALL cp_assert(.NOT.ASSOCIATED(ptr%parent),cp_warning_level,cp_assertion_failed,&
                   routineP,"try to deallocate subtree element"//&
                   cp_to_string(ptr%sub_tree_nr)//cp_to_string(ptr%nr)//&
                   " still with global tree element references e.g."//&
                   cp_to_string(ptr%gt_nodes_references%gt_elem%nr),&
                   failure)
      CPASSERT(ASSOCIATED(ptr%gt_nodes_references%gt_elem))
    ELSE
      SELECT CASE(ptr%stat)
      ! if element is still in progress, do not delete, wait for responding 
      CASE(status_calculate_energy, &
           status_calculate_NMC_steps, status_calculate_MD)
        ! in case of speculative canceling: should be already canceled
        !  try to deallocate subtree element (still in progress)
        CPASSERT(tmc_env%params%SPECULATIVE_CANCELING)
      CASE(status_cancel_nmc, status_cancel_ener)
        ! do not return in case of finalizing (do not wait for canceling receipt)
        IF(PRESENT(draw)) ready = .FALSE. 
      CASE DEFAULT
      END SELECT
  
      ! check if real top to bottom or bottom to top deallocation (no middle element deallocation)
      IF(ASSOCIATED(ptr%parent).AND.&
         (ASSOCIATED(ptr%acc).OR.ASSOCIATED(ptr%nacc)))THEN
        CPASSERT(.FALSE.)
      END IF
 
      IF(ready) THEN 
        ! set status and draw in tree
        IF((ptr%stat.EQ.status_accepted_result).OR.&
           (ptr%stat.EQ.status_rejected_result))THEN
           ptr%stat = status_deleted_result
        ELSE
           ptr%stat = status_deleted
        END IF
        IF(tmc_env%params%DRAW_TREE.AND.PRESENT(draw))&
          CALL create_dot_color(tree_element=ptr,tmc_params=tmc_env%params)
  
        !remove pointer from tree parent
        IF(ASSOCIATED(ptr%parent)) THEN
           IF(ASSOCIATED(ptr%parent%acc,  ptr)) ptr%parent%acc  => NULL()
           IF(ASSOCIATED(ptr%parent%nacc, ptr)) ptr%parent%nacc => NULL()
        END IF
  
        !remove pointer from tree childs
        IF(ASSOCIATED(ptr%acc)) ptr%acc%parent => NULL()
        IF(ASSOCIATED(ptr%nacc)) ptr%nacc%parent => NULL()
        
        ! deallocate
        CALL deallocate_sub_tree_node(tree_elem=ptr)
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE remove_st_elem

! *****************************************************************************
!> \brief deletes the no more used global tree nodes beside the result nodes 
!>        from begin_ptr to end_ptr
!> \param begin_ptr start of the tree region to be cleaned
!> \param end_ptr end of the tree region to be cleaned
!> \param removed retun value if brance is clean
!> \param tmc_env tmc environment
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE remove_unused_g_tree(begin_ptr, end_ptr, removed, tmc_env)
    TYPE(global_tree_type), POINTER          :: begin_ptr, end_ptr
    LOGICAL                                  :: removed
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_unused_g_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: acc_removed, failure, &
                                                nacc_removed
    TYPE(global_tree_type), POINTER          :: acc_ptr, nacc_ptr, tmp_ptr

    failure = .FALSE.
    NULLIFY(acc_ptr, nacc_ptr, tmp_ptr)

    CPASSERT(ASSOCIATED(begin_ptr))
    CPASSERT(ASSOCIATED(end_ptr))
    CPASSERT(ASSOCIATED(tmc_env))

    ! start the timing
    CALL timeset(routineN,handle)

    removed=.FALSE.
    acc_removed=.FALSE.
    nacc_removed=.FALSE.

    IF(.NOT.ASSOCIATED(begin_ptr, end_ptr)) THEN
      !-- go until the ends ot he tree, to deallocate revese
      !-- check if child nodes exist and possibly deallocate child node
      IF(ASSOCIATED(begin_ptr%acc))THEN
         acc_ptr => begin_ptr%acc
         CALL remove_unused_g_tree(acc_ptr, end_ptr, acc_removed, tmc_env)
      ELSE
         acc_removed = .TRUE.
      END IF
      IF(ASSOCIATED(begin_ptr%nacc))THEN
         nacc_ptr => begin_ptr%nacc
         CALL remove_unused_g_tree(nacc_ptr, end_ptr, nacc_removed, tmc_env)
      ELSE
         nacc_removed = .TRUE.
      END IF
  
      !-- deallocate node if no child node exist
      IF(acc_removed.AND.nacc_removed)THEN
        SELECT CASE(begin_ptr%stat)
        CASE(status_accepted, status_rejected, status_calculated, status_created, &
             status_calculate_energy, status_calculate_MD, status_calculate_NMC_steps, status_calc_approx_ener, &
             status_cancel_nmc,status_cancel_ener, status_canceled_nmc, status_canceled_ener)
          ! delete references, cancel elements calculation and deallocate global tree element
          tmp_ptr => begin_ptr
  
          CALL remove_gt_elem(gt_ptr=tmp_ptr, draw=.TRUE., tmc_env=tmc_env)
          IF(.NOT.ASSOCIATED(tmp_ptr)) removed = .TRUE.
        CASE(status_accepted_result, status_rejected_result)
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"try to dealloc unused tree element with status of begin element"& 
                  //cp_to_string(begin_ptr%stat),failure)
        END SELECT
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
    CPASSERT(ASSOCIATED(end_ptr))
  END SUBROUTINE remove_unused_g_tree

! *****************************************************************************
!> \brief deletes the no more used sub tree nodes beside the result nodes 
!>        from begin_ptr to end_ptr
!> \param begin_ptr start of the tree region to be cleaned
!> \param end_ptr end of the tree region to be cleaned
!> \param working_elem_list ...
!> \param removed retun value if brance is clean
!> \param tmc_env tmc environment
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE remove_unused_s_tree(begin_ptr, end_ptr,working_elem_list,&
                                          removed, tmc_env)
    TYPE(tree_type), POINTER                 :: begin_ptr
    TYPE(tree_type), INTENT(IN), POINTER     :: end_ptr
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: working_elem_list
    LOGICAL                                  :: removed
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_unused_s_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: acc_removed, failure, &
                                                nacc_removed, remove_this
    TYPE(tree_type), POINTER                 :: acc_ptr, nacc_ptr, tmp_ptr

    failure = .FALSE.
    NULLIFY(acc_ptr, nacc_ptr, tmp_ptr)
    remove_this  = .FALSE.
    removed      = .FALSE.
    acc_removed  = .FALSE.
    nacc_removed = .FALSE.

    ! start the timing
    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(begin_ptr))
    CPASSERT(ASSOCIATED(end_ptr))
    CPASSERT(ASSOCIATED(working_elem_list))
    CPASSERT(ASSOCIATED(tmc_env))

    !-- if element is last checked in trajectory, go back
    IF(.NOT.ASSOCIATED(begin_ptr, end_ptr)) THEN
      !-- go until the ends on the tree, to deallocate revesely
      !-- check if child nodes exist and possibly deallocate child node
      IF(ASSOCIATED(begin_ptr%acc))THEN
         acc_ptr => begin_ptr%acc
         CALL remove_unused_s_tree(acc_ptr, end_ptr, working_elem_list, &
                                 acc_removed, tmc_env)
      ELSE
         acc_removed = .TRUE.
      END IF
      IF(ASSOCIATED(begin_ptr%nacc))THEN
         nacc_ptr => begin_ptr%nacc
         CALL remove_unused_s_tree(nacc_ptr, end_ptr, working_elem_list, &
                                 nacc_removed, tmc_env)
      ELSE
         nacc_removed = .TRUE.
      END IF
  
      !IF (DEBUG.GE.20) WRITE(tmc_out_file_nr,*)"try to dealloc: node", begin_ptr%nr," sides are removed: ", &
      !                                          acc_removed, nacc_removed
  
      !-- deallocate node if NO child node exist
      ! unused trajectory is deleted when cleaned part is updated
      IF(acc_removed.AND.nacc_removed)THEN
        SELECT CASE(begin_ptr%stat)
        CASE(status_canceled_nmc, status_canceled_ener)
          remove_this = .TRUE.
        CASE(status_accepted, status_rejected, status_calculated, &
             status_accepted_result, status_rejected_result,status_created)
          remove_this = .TRUE.
        ! not to cancel, because still in progress
        CASE(status_calculate_energy, status_calculate_NMC_steps, &
             status_calculate_MD, status_cancel_nmc, status_cancel_ener, &
             status_calc_approx_ener)
           remove_this = .FALSE.
           ! -- delete when calculation is finished or aborted
           ! removed should still be .FALSE.
        CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "unknown status "//cp_to_string(begin_ptr%stat)//&
                     "of sub tree element "//&
                     cp_to_string(begin_ptr%sub_tree_nr)//" "//&
                     cp_to_string(begin_ptr%nr))
        END SELECT
  
        ! delete element
        IF(remove_this)THEN
          !-- mark as deleted and draw it in tree
          CALL cp_assert(ASSOCIATED(begin_ptr%parent),cp_failure_level,&
                         cp_assertion_failed,routineP,&
                         "try to remove unused subtree element "//&
                         cp_to_string(begin_ptr%sub_tree_nr)//" "//&
                         cp_to_string(begin_ptr%nr)//&
                         " but parent does not exist")
          tmp_ptr => begin_ptr
          ! check if a working group is still working on this element
          removed = .TRUE.
          DO i=1, SIZE(working_elem_list(:)) 
            IF(ASSOCIATED(working_elem_list(i)%elem)) THEN
              IF(ASSOCIATED(working_elem_list(i)%elem, tmp_ptr)) &
                removed=.FALSE.
            END IF
          END DO
          IF(removed) THEN
            !IF (DEBUG.GE.20) WRITE(tmc_out_file_nr,*)"deallocation of node ", begin_ptr%nr, "with status ", begin_ptr%stat
            ! if all groups are finished with this element, we can deallocate
            CALL remove_st_elem(ptr=tmp_ptr, draw=.TRUE., tmc_env=tmc_env)
            IF(.NOT.ASSOCIATED(tmp_ptr))THEN
              removed = .TRUE.
            ELSE
              removed = .FALSE.
            END IF
          END IF
        END IF
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE remove_unused_s_tree

! *****************************************************************************
!> \brief deallocates all result nodes (remaining Markov Chain) 
!>        from the tree root to the end of clean tree of the global tree
!> \param end_of_clean_tree ...
!> \param actual_ptr ...
!> \param tmc_env TMC environment for deallocation
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE remove_result_g_tree(end_of_clean_tree, actual_ptr, &
                                              tmc_env)
    TYPE(global_tree_type), POINTER          :: end_of_clean_tree, actual_ptr
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_result_g_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(global_tree_type), POINTER          :: tmp_ptr

    failure = .FALSE.
    CPASSERT(ASSOCIATED(end_of_clean_tree))
    CPASSERT(ASSOCIATED(actual_ptr))

    ! start the timing
    CALL timeset(routineN,handle)

    !-- going up to the head ot the subtree
    IF(ASSOCIATED(actual_ptr%parent))&
      CALL remove_result_g_tree(end_of_clean_tree=end_of_clean_tree, &
                                  actual_ptr=actual_ptr%parent, &
                                  tmc_env=tmc_env)
    !-- new tree head has no parent
    IF(.NOT.ASSOCIATED(actual_ptr, end_of_clean_tree)) THEN
      !-- deallocate node
      !IF(DEBUG.GE.20) WRITE(tmc_out_file_nr,*)"dealloc gt result tree element: ",actual_ptr%nr
      tmp_ptr => actual_ptr
      CALL remove_gt_elem(gt_ptr=tmp_ptr, draw=.TRUE., tmc_env=tmc_env)
      actual_ptr => tmp_ptr
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE remove_result_g_tree

! *****************************************************************************
!> \brief deallocates all result nodes (remaining Markov Chain) 
!>        from the tree root to the end of clean tree of one sub tree
!>        top to buttom deallocation
!> \param end_of_clean_tree ...
!> \param actual_ptr ...
!> \param tmc_env TMC environment for deallocation
!> \author Mandes 12.2012
! *****************************************************************************
  RECURSIVE SUBROUTINE remove_result_s_tree(end_of_clean_tree, actual_ptr,&
                                           tmc_env)
    TYPE(tree_type), POINTER                 :: end_of_clean_tree, actual_ptr
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_result_s_tree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(tree_type), POINTER                 :: tmp_ptr

    failure = .FALSE.
    CPASSERT(ASSOCIATED(end_of_clean_tree))
    CPASSERT(ASSOCIATED(actual_ptr))
    CPASSERT(ASSOCIATED(tmc_env))

    ! start the timing
    CALL timeset(routineN,handle)

    !-- going up to the head ot the subtree
    IF(ASSOCIATED(actual_ptr%parent)) &
      CALL remove_result_s_tree(end_of_clean_tree, actual_ptr%parent, &
                               tmc_env)

    !-- new tree head has no parent
    IF(.NOT.ASSOCIATED(actual_ptr, end_of_clean_tree)) THEN
      ! in trajectory just one direction should exist
      IF(ASSOCIATED(actual_ptr%acc).AND.ASSOCIATED(actual_ptr%nacc)) THEN
      CPASSERT(.FALSE.)
      END IF
      ! the parent should be deleted already, but global tree is allocated to the second last accepted, &
      !   hence there could be still a reference to an element...
      IF(.NOT.ASSOCIATED(actual_ptr%parent)) THEN
        !-- deallocate node
        tmp_ptr => actual_ptr
        CALL remove_st_elem(ptr=tmp_ptr,draw=.TRUE.,tmc_env=tmc_env)
        actual_ptr => tmp_ptr
      END IF
    END IF
    ! end the timing
    CALL timestop(handle)
  END SUBROUTINE remove_result_s_tree

! *****************************************************************************
!> \brief deallocates the no more used tree nodes beside the result nodes 
!>        from begin_ptr to end_ptr
!>        in global and subtrees
!> \param working_elem_list list of actual calculating elements for canceling
!> \param tmc_env TMC environment
!> \author Mandes 12.2012
! *****************************************************************************
  SUBROUTINE remove_all_trees(working_elem_list, tmc_env)
    TYPE(elem_array_type), DIMENSION(:), &
      POINTER                                :: working_elem_list
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_all_trees', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, tree
    LOGICAL                                  :: change_trajec, failure, flag
    TYPE(global_tree_type), POINTER          :: tmp_gt_ptr
    TYPE(tree_type), POINTER                 :: last_acc_st_elem, tmp_ptr

    failure = .FALSE.
    NULLIFY(last_acc_st_elem, tmp_ptr, tmp_gt_ptr)

    CPASSERT(ASSOCIATED(working_elem_list))
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env%gt_act))
    CPASSERT(ASSOCIATED(tmc_env%m_env%gt_clean_end))
    CPASSERT(ASSOCIATED(tmc_env%m_env%result_list))
    CPASSERT(ASSOCIATED(tmc_env%m_env%st_clean_ends))

    flag = .FALSE.
    change_trajec = .FALSE.

    ! start the timing
    CALL timeset(routineN,handle)

    !-- deallocate unused pt tree
    CALL remove_unused_g_tree(begin_ptr=tmc_env%m_env%gt_clean_end, &
                              end_ptr=tmc_env%m_env%gt_act, removed=flag, &
                              tmc_env=tmc_env)
    tmp_gt_ptr => tmc_env%m_env%gt_clean_end
    CALL search_end_of_clean_g_tree(last_acc=tmc_env%m_env%gt_clean_end, &
                                    tree_ptr=tmp_gt_ptr)
    !-- deallocate unused pt trajectory tree elements
    IF(tmc_env%params%USE_REDUCED_TREE)THEN
      tmp_gt_ptr => tmc_env%m_env%gt_clean_end
      CALL remove_result_g_tree(end_of_clean_tree=tmc_env%m_env%gt_clean_end,&
                                actual_ptr=tmp_gt_ptr, tmc_env=tmc_env)

      !check if something changed, if not no deallocation of result subtree necessary
      IF(.NOT.ASSOCIATED(tmc_env%m_env%gt_head,tmc_env%m_env%gt_clean_end)) &
        change_trajec=.TRUE.
      tmc_env%m_env%gt_head => tmc_env%m_env%gt_clean_end
      CPASSERT(.NOT.ASSOCIATED(tmc_env%m_env%gt_head%parent))
      !IF (DEBUG.GE.20) WRITE(tmc_out_file_nr,*)"new head of pt tree is ",tmc_env%m_env%gt_head%nr
    END IF

    !-- deallocate the subtrees
    ! do for all temperatures respectively all subtrees
    DO tree=1, tmc_env%params%nr_temp
      ! get last checked element in trajectory related to the subtree (resultlist order is NOT subtree order)
      conf_loop: DO i=1, SIZE(tmc_env%m_env%result_list)
        last_acc_st_elem => tmc_env%m_env%result_list(i)%elem
        IF(last_acc_st_elem%sub_tree_nr .EQ. tree) &
          EXIT conf_loop
      END DO conf_loop
      CPASSERT(last_acc_st_elem%sub_tree_nr.EQ.tree)
      CALL remove_unused_s_tree(begin_ptr=tmc_env%m_env%st_clean_ends(tree)%elem, &
                 end_ptr=last_acc_st_elem, working_elem_list=working_elem_list, &
                 removed=flag, tmc_env=tmc_env)
      CALL search_end_of_clean_tree(tree_ptr=tmc_env%m_env%st_clean_ends(tree)%elem, &
                                    last_acc=last_acc_st_elem)
    END DO
    !-- deallocate the trajectory subtree elements
    IF(tmc_env%params%USE_REDUCED_TREE .AND. change_trajec)THEN
      DO tree=1, tmc_env%params%nr_temp
        tmp_ptr => tmc_env%m_env%st_clean_ends(tree)%elem
        CPASSERT(tmp_ptr%sub_tree_nr.EQ.tree)
        CALL remove_result_s_tree(end_of_clean_tree=tmc_env%m_env%st_clean_ends(tree)%elem,&
                                 actual_ptr=tmp_ptr, tmc_env=tmc_env)
        tmc_env%m_env%st_heads(tree)%elem => tmc_env%m_env%st_clean_ends(tree)%elem
        !IF(DEBUG.GE.20) &
        !  WRITE(tmc_out_file_nr,*)"new head of tree ",tree," is ",&
        !        tmc_env%m_env%st_heads(tree)%elem%nr
      END DO
    END IF

    ! end the timing
    CALL timestop(handle)
    CPASSERT(ASSOCIATED(tmc_env%m_env%gt_act))
    CPASSERT(ASSOCIATED(tmc_env%m_env%gt_clean_end))
  END SUBROUTINE remove_all_trees

! *****************************************************************************
!> \brief deallocates the whole global tree, to clean up 
!> \param begin_ptr pointer to global tree head
!> \param removed flag, if the this element is removed
!> \param tmc_env ...
!> \author Mandes 01.2013
! *****************************************************************************
  RECURSIVE SUBROUTINE dealloc_whole_g_tree(begin_ptr, removed, tmc_env)
    TYPE(global_tree_type), POINTER          :: begin_ptr
    LOGICAL                                  :: removed
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dealloc_whole_g_tree', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: acc_removed, failure, &
                                                nacc_removed
    TYPE(global_tree_type), POINTER          :: acc_ptr, nacc_ptr, tmp_ptr

    failure = .FALSE.
    CPASSERT(ASSOCIATED(begin_ptr))
    CPASSERT(ASSOCIATED(tmc_env))

    IF(ASSOCIATED(begin_ptr%acc))THEN
      acc_ptr => begin_ptr%acc
      CALL dealloc_whole_g_tree(acc_ptr, acc_removed, tmc_env)
    ELSE
      acc_removed = .TRUE.
    END IF
    IF(ASSOCIATED(begin_ptr%nacc))THEN
      nacc_ptr => begin_ptr%nacc
      CALL dealloc_whole_g_tree(nacc_ptr, nacc_removed, tmc_env)
    ELSE
      nacc_removed = .TRUE.
    END IF

    !-- deallocate node if no child node exist
    IF(acc_removed.AND.nacc_removed)THEN
      CALL search_and_remove_reference_in_list(gt_ptr=begin_ptr, &
             elem=begin_ptr%conf(begin_ptr%mv_conf)%elem, tmc_env=tmc_env)
      tmp_ptr => begin_ptr
      CALL remove_gt_elem(gt_ptr=tmp_ptr, draw=.FALSE., tmc_env=tmc_env)
      !CALL deallocate_global_tree_node(gt_elem=tmp_ptr)
      removed = .TRUE.
    END IF
  END SUBROUTINE dealloc_whole_g_tree
! *****************************************************************************
!> \brief deallocates the whole sub tree, to clean up 
!> \param begin_ptr pointer to sub tree head
!> \param removed flag, if the this element is removed
!> \param tmc_params ...
!> \author Mandes 01.2013
! *****************************************************************************
  RECURSIVE SUBROUTINE dealloc_whole_subtree(begin_ptr, removed, tmc_params)
    TYPE(tree_type), POINTER                 :: begin_ptr
    LOGICAL                                  :: removed
    TYPE(tmc_param_type), POINTER            :: tmc_params

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dealloc_whole_subtree', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: acc_removed, failure, &
                                                nacc_removed
    TYPE(tree_type), POINTER                 :: acc_ptr, nacc_ptr, tmp_ptr

    failure = .FALSE.
    CPASSERT(ASSOCIATED(begin_ptr))
    CPASSERT(ASSOCIATED(tmc_params))

    IF(ASSOCIATED(begin_ptr%acc))THEN
      acc_ptr => begin_ptr%acc
      CALL dealloc_whole_subtree(acc_ptr, acc_removed, tmc_params)
    ELSE
       acc_removed = .TRUE.
    END IF
    IF(ASSOCIATED(begin_ptr%nacc))THEN
      nacc_ptr => begin_ptr%nacc
      CALL dealloc_whole_subtree(nacc_ptr, nacc_removed, tmc_params)
    ELSE
      nacc_removed = .TRUE.
    END IF

    !-- deallocate node if no child node exist
    IF(acc_removed.AND.nacc_removed)THEN
      tmp_ptr => begin_ptr
      CALL deallocate_sub_tree_node(tree_elem=begin_ptr)
      removed = .TRUE.
    END IF
  END SUBROUTINE dealloc_whole_subtree

  !============================================================================
  ! finalizing module (deallocating everything)
  !============================================================================
! *****************************************************************************
!> \brief deallocating every tree node of every trees (clean up)
!> \param tmc_env TMC environment structure
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE finalize_trees(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'finalize_trees', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, flag
    TYPE(global_tree_type), POINTER          :: global_tree

    failure = .FALSE.
    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env))

    global_tree => tmc_env%m_env%gt_act
    !-- deallocate pt tree
    ! start with searching the head
    DO WHILE(ASSOCIATED(global_tree%parent))
      global_tree => global_tree%parent
    END DO
    CALL dealloc_whole_g_tree(begin_ptr=global_tree, removed=flag, &
                               tmc_env=tmc_env)

    !-- deallocate subtrees
    trees_loop: DO i=1, SIZE(tmc_env%m_env%st_clean_ends(:))
      DO WHILE(ASSOCIATED(tmc_env%m_env%st_clean_ends(i)%elem%parent))
        tmc_env%m_env%st_clean_ends(i)%elem => &
          tmc_env%m_env%st_clean_ends(i)%elem%parent
      END DO
      CALL dealloc_whole_subtree(begin_ptr=tmc_env%m_env%st_clean_ends(i)%elem,&
             removed=flag, tmc_params=tmc_env%params)
    END DO trees_loop
    DEALLOCATE(tmc_env%params%atoms)
  END SUBROUTINE finalize_trees

END MODULE tmc_tree_build
