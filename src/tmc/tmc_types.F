!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module handels definition of the tree nodes for the global and 
!>      the subtrees binary tree 
!>                   parent element
!>                      /      \
!>      accepted (acc) /        \  not accepted (nacc)
!>                    /          \
!>                  child       child
!>                   / \         / \
!>     
!>      tree creation assuming acceptance (acc) AND rejectance (nacc) 
!>        of configuration
!>      if configuration is accepted: new configuration (child on acc) on basis 
!>        of last configuration (one level up)
!>      if configuration is rejected: child on nacc on basis of last accepted 
!>        element (last element which is on acc brach of its parent element)
!>      The global tree handles all configurations of different subtrees. 
!>      The structure element "conf" is an array related to the temperature 
!>        (sorted) and points to the subtree elements.
!> \par History
!>      11.2012 created [Mandes Schönherr]
!> \author Mandes
! *****************************************************************************

MODULE tmc_types
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE parallel_rng_types,              ONLY: rng_stream_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE tmc_move_types,                  ONLY: tmc_move_type
  USE tmc_stati,                       ONLY: task_type_MC
  USE tmc_tree_types,                  ONLY: clean_list,&
                                             elem_array_type,&
                                             elem_list_type,&
                                             global_tree_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_types'

  PUBLIC :: tmc_env_type, tmc_comp_set_type, tmc_param_type
  PUBLIC :: tmc_env_create, tmc_env_release,&
            tmc_master_env_create, tmc_master_env_release,&
            tmc_worker_env_create, tmc_worker_env_release
  PUBLIC :: tmc_atom_type
  PUBLIC :: allocate_tmc_atom_type

  ! global environment
  TYPE tmc_env_type
    TYPE(tmc_comp_set_type), POINTER              :: tmc_comp_set
    TYPE(tmc_param_type), POINTER                 :: params
    TYPE(rng_stream_type), POINTER                :: rng_stream
    TYPE(master_env_type), POINTER                :: m_env
    TYPE(worker_env_type), POINTER                :: w_env
  END TYPE tmc_env_type  

  ! structure for remembering the main values used for reordering MPI communicators
  ! \param group_nr: the first group_ener_nr groups are for energy calculation, 
  !                  then the Configurational Change groups,
  !                  the global master has group 0 and unused cores have negative group numbering
  TYPE tmc_comp_set_type
    INTEGER                                       :: group_ener_size
    INTEGER                                       :: group_ener_nr
    INTEGER                                       :: group_cc_size
    INTEGER                                       :: group_cc_nr
    INTEGER                                       :: group_nr
    INTEGER                                       :: ana_on_the_fly
    ! the communicators (para_env) 
    TYPE(cp_para_env_type), POINTER               :: para_env_m_w
    TYPE(cp_para_env_type), POINTER               :: para_env_sub_group
    TYPE(cp_para_env_type), POINTER               :: para_env_m_first_w
    TYPE(cp_para_env_type), POINTER               :: para_env_m_ana
    TYPE(cp_para_env_type), POINTER               :: para_env_m_only
  END TYPE tmc_comp_set_type

  ! struct for TMC global variables
  TYPE tmc_param_type
    INTEGER                                       :: task_type
    INTEGER                                       :: dim_per_elem
    INTEGER                                       :: nr_temp
    REAL(KIND=dp), DIMENSION(:), POINTER          :: Temp
    TYPE(cell_type), POINTER                      :: cell
    REAL(KIND=dp), DIMENSION(:), POINTER          :: sub_box_size
    TYPE(tmc_atom_type), DIMENSION(:), POINTER    :: atoms
    
    INTEGER                                       :: nr_elem_mv
    TYPE(tmc_move_type), POINTER                  :: move_types
    TYPE(tmc_move_type), POINTER                  :: nmc_move_types
    REAL(KIND=dp)                                 :: pressure
    LOGICAL                                       :: v_isotropic
    LOGICAL                                       :: mv_cen_of_mass
    LOGICAL                                       :: esimate_acc_prob
    LOGICAL                                       :: SPECULATIVE_CANCELING
    LOGICAL                                       :: use_scf_energy_info
    LOGICAL                                       :: USE_REDUCED_TREE
    CHARACTER(LEN= default_path_length )          :: energy_inp_file
    CHARACTER(LEN= default_path_length )          :: NMC_inp_file
    LOGICAL                                       :: DRAW_TREE
    CHARACTER(LEN= default_path_length )          :: dot_file_name
    CHARACTER(LEN= default_path_length )          :: all_conf_file_name
    LOGICAL                                       :: print_only_diff_conf
    LOGICAL                                       :: print_trajectory
    LOGICAL                                       :: print_dipole
    LOGICAL                                       :: print_forces
    LOGICAL                                       :: print_cell
    LOGICAL                                       :: print_energies
    TYPE(prior_estimate_acceptance_type), POINTER :: prior_NMC_acc
    LOGICAL                                       :: print_test_output
  END TYPE tmc_param_type

  TYPE tmc_atom_type
    CHARACTER(LEN=default_string_length)   :: name
    REAL(KIND = dp)                        :: mass
    !TYPE(atomic_kind_type), POINTER :: atomic_kind
  END TYPE

  ! to estimate the prior acceptance
  TYPE prior_estimate_acceptance_type
    INTEGER :: counter
    REAL(KIND=dp) :: aver, aver_2
  END TYPE prior_estimate_acceptance_type

  ! environments for the master
  TYPE master_env_type
    INTEGER                                       :: num_MC_elem ! the specified number of Markov Chain elements, to be reached
    CHARACTER(LEN= default_path_length )          :: restart_in_file_name
    CHARACTER(LEN= default_path_length )          :: restart_out_file_name
    INTEGER                                       :: restart_out_step
    INTEGER                                       :: io_unit
    INTEGER                                       :: info_out_step_size
    REAL(KIND=dp)                                 :: walltime
    INTEGER                                       :: rnd_init
    REAL(KIND=dp)                                 :: temp_decrease ! for simulated annealing
    TYPE(elem_list_type), POINTER                 :: cancelation_list
    INTEGER                                       :: count_cancel_ener
    INTEGER                                       :: count_cancel_NMC
    ! masters tree stuff
    TYPE(global_tree_type), POINTER               :: gt_head, gt_act
    INTEGER, DIMENSION(:), POINTER                :: tree_node_count
    INTEGER, DIMENSION(:), POINTER                :: result_count
    TYPE (elem_array_type), DIMENSION(:), &
         POINTER                                  :: result_list, &
                                                     st_heads, &
                                                     st_clean_ends
    TYPE (global_tree_type), POINTER              :: gt_clean_end
    INTEGER, DIMENSION(4)                         :: estim_corr_wrong
    TYPE(elem_list_type), POINTER                 :: analysis_list
  END TYPE master_env_type

  ! environment for the worker
  TYPE worker_env_type
    INTEGER                                       :: env_id_ener, env_id_approx
    INTEGER                                       :: io_unit
    REAL(KIND=dp)                                 :: act_temp
  END TYPE worker_env_type

CONTAINS

! *****************************************************************************
!> \brief creates a new structure environment for TMC
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_env_create(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(.NOT.ASSOCIATED(tmc_env))

    ALLOCATE(tmc_env,STAT=stat)
    CPASSERT(stat==0)

    ALLOCATE(tmc_env%tmc_comp_set,STAT=stat)
    CPASSERT(stat==0)

    tmc_env%tmc_comp_set%group_ener_size = 0 
    tmc_env%tmc_comp_set%group_ener_nr   = 0
    tmc_env%tmc_comp_set%group_cc_size   = 0
    tmc_env%tmc_comp_set%group_cc_nr     = 0
    tmc_env%tmc_comp_set%group_nr        = 0 
    tmc_env%tmc_comp_set%ana_on_the_fly  = 0
    NULLIFY(tmc_env%tmc_comp_set%para_env_m_w)
    NULLIFY(tmc_env%tmc_comp_set%para_env_sub_group)
    NULLIFY(tmc_env%tmc_comp_set%para_env_m_first_w)
    NULLIFY(tmc_env%tmc_comp_set%para_env_m_ana)
    NULLIFY(tmc_env%tmc_comp_set%para_env_m_only)


    NULLIFY(tmc_env%rng_stream)
    NULLIFY(tmc_env%m_env, tmc_env%w_env)

    ! initialize the parameter section
    ALLOCATE(tmc_env%params,STAT=stat)
    CPASSERT(stat==0)

    tmc_env%params%task_type = task_type_MC
    ! dim_per_elem can be variable somewhen, 
    !  I used it originally for the n dimension oscillator, 
    !  where each dimension is a single position
    tmc_env%params%dim_per_elem = 3 
    tmc_env%params%nr_elem_mv = -1
    tmc_env%params%nr_temp = -1 
    ALLOCATE(tmc_env%params%sub_box_size(tmc_env%params%dim_per_elem))
    tmc_env%params%sub_box_size(:) = -1.0_dp
    tmc_env%params%pressure = -1
    tmc_env%params%SPECULATIVE_CANCELING = .FALSE.
    tmc_env%params%use_scf_energy_info = .FALSE.
    tmc_env%params%energy_inp_file = ""
    tmc_env%params%NMC_inp_file = ""
    tmc_env%params%DRAW_TREE = .FALSE.
    tmc_env%params%print_dipole = .FALSE.
    tmc_env%params%print_forces = .FALSE.
    tmc_env%params%print_cell = .FALSE.
    NULLIFY(tmc_env%params%Temp, tmc_env%params%nmc_move_types, &
            tmc_env%params%move_types, &
            tmc_env%params%prior_NMC_acc, tmc_env%params%atoms, &
            tmc_env%params%cell)

    CALL timestop(handle)

  END SUBROUTINE tmc_env_create

! *****************************************************************************
!> \brief releases the structure environment for TMC
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_env_release(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%params))

    DEALLOCATE(tmc_env%params%sub_box_size)
    IF(ASSOCIATED(tmc_env%params%Temp)) &
      DEALLOCATE(tmc_env%params%Temp)
    IF(ASSOCIATED(tmc_env%params%cell)) &
      DEALLOCATE(tmc_env%params%cell)
    IF(ASSOCIATED(tmc_env%params%atoms)) &
      CALL deallocate_tmc_atom_type(tmc_env%params%atoms)
    DEALLOCATE(tmc_env%params,STAT=stat)
    CPASSERT(stat==0)

    CALL cp_para_env_release(tmc_env%tmc_comp_set%para_env_sub_group)
    CALL cp_para_env_release(tmc_env%tmc_comp_set%para_env_m_w)
    IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_first_w))&
      CALL cp_para_env_release(tmc_env%tmc_comp_set%para_env_m_first_w)
    IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_ana))&
    CALL cp_para_env_release(tmc_env%tmc_comp_set%para_env_m_ana)
    IF(ASSOCIATED(tmc_env%tmc_comp_set%para_env_m_only))&
    CALL cp_para_env_release(tmc_env%tmc_comp_set%para_env_m_only)

    DEALLOCATE(tmc_env%tmc_comp_set,STAT=stat)
    CPASSERT(stat==0)

    DEALLOCATE(tmc_env,STAT=stat)
    CPASSERT(stat==0)

    CALL timestop(handle)

  END SUBROUTINE tmc_env_release

! *****************************************************************************
!> \brief creates a new structure environment for TMC master
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_master_env_create(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_master_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%params))
    CPASSERT(tmc_env%params%nr_temp.GT.0)
    
    CPASSERT(.NOT.ASSOCIATED(tmc_env%m_env))

    ALLOCATE(tmc_env%m_env)
    NULLIFY(tmc_env%m_env%gt_head, tmc_env%m_env%gt_act, tmc_env%m_env%tree_node_count, &
            tmc_env%m_env%result_count, tmc_env%m_env%result_list, &
            tmc_env%m_env%st_heads, tmc_env%m_env%st_clean_ends, &
            tmc_env%m_env%gt_clean_end, tmc_env%m_env%cancelation_list, tmc_env%m_env%analysis_list)

    tmc_env%m_env%restart_in_file_name = ""
    tmc_env%m_env%restart_out_file_name = ""
    ALLOCATE(tmc_env%m_env%tree_node_count(0:tmc_env%params%nr_temp))
    tmc_env%m_env%tree_node_count(:) = 0
    ALLOCATE(tmc_env%m_env%result_count(0:tmc_env%params%nr_temp))
    tmc_env%m_env%result_count(:) = 0
    ALLOCATE(tmc_env%m_env%st_heads(tmc_env%params%nr_temp))
    ALLOCATE(tmc_env%m_env%st_clean_ends(tmc_env%params%nr_temp))

    IF(tmc_env%params%USE_REDUCED_TREE) ALLOCATE(tmc_env%m_env%result_list(tmc_env%params%nr_temp))

    DO i=1, tmc_env%params%nr_temp
      tmc_env%m_env%st_heads(i)%elem      => NULL()
      tmc_env%m_env%st_clean_ends(i)%elem => NULL()
      IF(tmc_env%params%USE_REDUCED_TREE) &
        tmc_env%m_env%result_list(i)%elem => NULL()
    END DO
    tmc_env%m_env%gt_head           => NULL()
    tmc_env%m_env%gt_clean_end      => NULL()
    tmc_env%m_env%temp_decrease     = 1.0_dp
    tmc_env%m_env%count_cancel_ener = 0
    tmc_env%m_env%count_cancel_NMC  = 0
    tmc_env%m_env%estim_corr_wrong(:) = 0

    ALLOCATE(tmc_env%params%prior_NMC_acc)
    tmc_env%params%prior_NMC_acc%counter = 0
    tmc_env%params%prior_NMC_acc%aver    = 0.0_dp
    tmc_env%params%prior_NMC_acc%aver_2  = 0.0_dp

    CALL timestop(handle)

  END SUBROUTINE tmc_master_env_create

! *****************************************************************************
!> \brief releases the structure environment for TMC master
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_master_env_release(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_master_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%m_env))

    CALL clean_list(tmc_env%m_env%analysis_list)
    CALL clean_list(tmc_env%m_env%cancelation_list)

    DEALLOCATE(tmc_env%m_env%tree_node_count)
    DEALLOCATE(tmc_env%m_env%result_count)
    DEALLOCATE(tmc_env%m_env%st_heads)
    DEALLOCATE(tmc_env%m_env%st_clean_ends)
    IF(tmc_env%params%USE_REDUCED_TREE) DEALLOCATE(tmc_env%m_env%result_list)
    DEALLOCATE(tmc_env%params%prior_NMC_acc)

    DEALLOCATE(tmc_env%m_env)

    CALL timestop(handle)

  END SUBROUTINE tmc_master_env_release

! *****************************************************************************
!> \brief creates a new structure environment for TMC master
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_worker_env_create(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_worker_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(.NOT.ASSOCIATED(tmc_env%w_env))

    ALLOCATE(tmc_env%w_env)

    tmc_env%w_env%env_id_ener   = -1
    tmc_env%w_env%env_id_approx = -1
    tmc_env%w_env%io_unit       = -1
    tmc_env%w_env%act_temp      = -1.0_dp 

    CALL timestop(handle)

  END SUBROUTINE tmc_worker_env_create

! *****************************************************************************
!> \brief releases the structure environment for TMC master
!> \param tmc_env structure with parameters for TMC
!> \author Mandes 11.2012
! *****************************************************************************
  SUBROUTINE tmc_worker_env_release(tmc_env)
    TYPE(tmc_env_type), POINTER              :: tmc_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_worker_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(tmc_env))
    CPASSERT(ASSOCIATED(tmc_env%w_env))

    DEALLOCATE(tmc_env%w_env)

    CALL timestop(handle)

  END SUBROUTINE tmc_worker_env_release

! *****************************************************************************
!> \brief creates a structure for storing the atom informations
!> \param atoms pointer to a list of tmc_atoms_type
!> \param nr_atoms the amount of atoms
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE allocate_tmc_atom_type(atoms, nr_atoms)
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms
    INTEGER, INTENT(IN)                      :: nr_atoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_tmc_atom_type', &
      routineP = moduleN//':'//routineN

    CPASSERT(.NOT.ASSOCIATED(atoms))
    CPASSERT(nr_atoms.GT.0)

    ALLOCATE(atoms(nr_atoms))
    atoms%name = ""
    atoms%mass = 0.0_dp

    CPASSERT(ASSOCIATED(atoms))
    
  END SUBROUTINE allocate_tmc_atom_type

! *****************************************************************************
!> \brief releases the structure for storing the atom informations
!> \param atoms pointer to a list of tmc_atoms_type
!> \author Mandes 01.2013
! *****************************************************************************
  SUBROUTINE deallocate_tmc_atom_type(atoms)
    TYPE(tmc_atom_type), DIMENSION(:), &
      POINTER                                :: atoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_tmc_atom_type', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(atoms))

    DEALLOCATE(atoms)

    CPASSERT(.NOT.ASSOCIATED(atoms))
  END SUBROUTINE deallocate_tmc_atom_type

END MODULE tmc_types
