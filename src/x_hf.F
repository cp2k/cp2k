!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/x_hf [1.0] *
!!
!!   NAME
!!     x_hf
!!
!!   FUNCTION
!!     routines and types for Hartree-Fock-Exchange
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************
MODULE x_hf   
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE x_hf_compression_methods,        ONLY: bits2ints_specific,&
                                             ints2bits_specific
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  x_hf_type, x_hf_create, x_hf_release,& 
          x_hf_set_distr_energy, &
          x_hf_set_distr_forces, x_hf_create_neighbor_cells,&
          x_hf_cell_type, x_hf_distribution, x_hf_box_scheme,&
          x_hf_potential_type, x_hf_screening_type, x_hf_periodic_type,&
          x_hf_memory_type, x_hf_load_balance_type, x_hf_general_type,&
          x_hf_compression_type,&
          x_hf_container_type, x_hf_container_node, x_hf_cache_type, &
          x_hf_add_single_cache_element, x_hf_get_single_cache_element,&
          x_hf_reset_cache_and_container, x_hf_decompress_first_cache, &
          x_hf_flush_last_cache, x_hf_get_memory_info, x_hf_add_mult_cache_elements,&
          x_hf_get_mult_cache_elements, x_hf_init_container

#define CACHE_SIZE 1024

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'x_hf'
  INTEGER(int_8), PARAMETER, PRIVATE :: one = 1_int_8
!***
!****************************************************************************


  TYPE x_hf_potential_type
    INTEGER                                  :: potential_type     !! 1/r/ erfc(wr)/r ...
    REAL(dp)                                 :: omega              !! w 
  END TYPE

  TYPE x_hf_screening_type
    REAL(dp)                                 :: eps_schwarz         !! threshold
    REAL(dp)                                 :: eps_schwarz_forces  !! threshold
    REAL(dp)                                 :: eps_box             !! threshold
    REAL(dp)                                 :: eps_box_forces      !! threshold
    LOGICAL                                  :: do_p_screening_forces  !! screen on P^2 ?
    LOGICAL                                  :: do_initial_p_screening !! screen on inital guess?
    LOGICAL                                  :: do_init_p_screen_corr  !! correct screen on initial guess?
    REAL(dp)                                 :: initial_p_screening_correction !! corection factor
  END TYPE
  
  TYPE x_hf_memory_type
    INTEGER                                  :: max_memory              !! user def max memory MB
    INTEGER(int_8)                           :: max_compression_counter !! corresponding number of reals 
    LOGICAL                                  :: do_all_on_the_fly       !! max mem == 0 ?
    REAL(dp)                                 :: eps_storage_scaling
    INTEGER                                  :: cache_size
  END TYPE

  TYPE x_hf_periodic_type
    INTEGER                                  :: number_of_shells        !! number of periodic image cells 
    LOGICAL                                  :: do_periodic             !! periodic ?
  END TYPE

  TYPE x_hf_load_balance_type
    INTEGER                                  :: nbins                   
    INTEGER                                  :: moves
    INTEGER                                  :: cycles
    INTEGER                                  :: reduction_steps
    REAL(dp)                                 :: reduction_factor
    REAL(dp)                                 :: temperature
    REAL(dp)                                 :: tolerance
    REAL(dp)                                 :: swap_probability
  END TYPE

  TYPE x_hf_general_type
    REAL(dp)                                 :: fraction                 !! for hybrids
  END TYPE

  TYPE x_hf_cell_type
    INTEGER                                  :: cell(3)
  END TYPE

  TYPE x_hf_distribution
    INTEGER                                  :: istart, jstart, kstart, lstart
    INTEGER(int_8)                           :: number_of_integrals
    INTEGER(int_8)                           :: cost
    INTEGER(int_8)                           :: original_cpu_id
  END TYPE
  
  TYPE x_hf_box_scheme
    REAL(dp)                                 :: center(3)
    REAL(dp)                                 :: K
    REAL(dp)                                 :: R
    REAL(dp)                                 :: alpha
  END TYPE
  
  TYPE x_hf_compression_type
    INTEGER(int_8), DIMENSION(:), POINTER    :: compressed_data
    INTEGER(int_8)                           :: element_pointer
    INTEGER                                  :: bit_pointer
  END TYPE
 
  TYPE x_hf_cache_type
    INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
    INTEGER                                  :: element_counter
  END TYPE

  TYPE x_hf_container_node
    TYPE(x_hf_container_node), POINTER       :: next, prev
    INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
  END TYPE

  TYPE x_hf_container_type
    TYPE(x_hf_container_node), POINTER       :: first, current    
    INTEGER                                  :: element_counter
    INTEGER                                  :: list_counter
  END TYPE


!!***** x_hf/x_hf_type [1.0] *
!!
!!   NAME
!!     x_hf_type
!!
!!   FUNCTION
!!     stores some data used in construction of Kohn-Sham matrix
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - potential_parameter: stores information on the potential (1/r, erfc(wr)/r
!!     - screening_parameter: stores screening infos such as epsilon
!!     - memory_parameter: stores infos on memory used for in-core calculations
!!     - periodic_parameter: stores information on how to apply pbc
!!     - load_balance_parameter: contains infos for Monte Carlo simulated annealing
!!     - general_paramter: at the moment stores the fraction of HF amount to be included
!!     - maxval_container: stores the maxvals in compressed form
!!     - maxval_cache: cache for maxvals in decompressed form
!!     - integral_containers: 64 containers for compressed integrals
!!     - integral_caches: 64 caches for decompressed integrals
!!     - x_screening_buffer: contains screening matrices in non-periodic case
!!     - x_screening_buffer_periodic: contains screening matrices in periodic case
!!     - far_field_boxes: contains information on far field screening
!!     - neighbor_cells: manages handling of periodic cells
!!     - distribution_energy: stores information on parallelization of energy
!!     - distribution_forces: stores information on parallelization of forces
!!     - initial_p: stores the initial guess if requested
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!***************************************************************************


  TYPE x_hf_type
    TYPE(x_hf_potential_type)                :: potential_parameter
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(x_hf_memory_type)                   :: memory_parameter
    TYPE(x_hf_periodic_type)                 :: periodic_parameter
    TYPE(x_hf_load_balance_type)             :: load_balance_parameter
    TYPE(x_hf_general_type)                  :: general_parameter
    TYPE(x_hf_container_type), POINTER       :: maxval_container
    TYPE(x_hf_cache_type)                    :: maxval_cache
    TYPE(x_hf_container_type), POINTER, &
                               DIMENSION(:)  :: integral_containers
    TYPE(x_hf_cache_type)                    :: integral_caches(64)

    REAL(dp), DIMENSION(:, :, :, :), POINTER :: x_screening_buffer
    REAL(dp), DIMENSION(:, :, :, :, :), &
              POINTER                        :: x_screening_buffer_periodic
    TYPE(x_hf_box_scheme), DIMENSION(:,:),&
                           POINTER           :: far_field_boxes
    TYPE(x_hf_cell_type), DIMENSION(:), &
               POINTER                       :: neighbor_cells
    TYPE(x_hf_distribution), DIMENSION(:),&
                             POINTER         :: distribution_energy
    TYPE(x_hf_distribution), DIMENSION(:),&
                             POINTER         :: distribution_forces
        REAL(dp), DIMENSION(:,:,:,:), POINTER &
                                             :: initial_p
    INTEGER, DIMENSION(:,:), POINTER         :: is_assoc_atomic_block
    INTEGER                                  :: number_of_p_entries
  END TYPE

  CONTAINS



!!***** x_hf/x_hf_create [1.0] *
!!
!!   NAME
!!     x_hf_create
!!
!!   FUNCTION
!!     - This routine allocates and initializes all types in x_hf_data
!!
!!   NOTES
!!     - All POINTERS and ALLOCATABLES are allocated, even if their size is
!!       unknown at invocation time
!!
!!   INPUTS
!!     - x_data: contains all relevant data structures for hfx runs
!!     - hf_section: input section
!!     - natom, max_set: quantities that define some array sizes
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_create(x_data,hf_section,natom,max_set,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(section_vals_type), POINTER         :: hf_section
    INTEGER, INTENT(IN)                      :: natom, max_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, int_val, max_memory, &
                                                pbc_shells, stat
    LOGICAL                                  :: do_periodic, failure, &
                                                logic_val
    REAL(dp)                                 :: omega, real_val
    TYPE(section_vals_type), POINTER         :: hf_pbc_section, hf_sub_section

    failure = .FALSE.

    ALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !! Get data from input file
    !!
    !! GENERAL params
    CALL section_vals_val_get(hf_section, "FRACTION", r_val=real_val, error=error)
    x_data%general_parameter%fraction = real_val

    !! MEMORY section
    hf_sub_section => section_vals_get_subs_vals(hf_section,"MEMORY",error=error)
    CALL section_vals_val_get(hf_sub_section, "MAX_MEMORY", i_val=int_val, error=error)
    x_data%memory_parameter%max_memory = int_val
    x_data%memory_parameter%max_compression_counter = int_val*1024*128
    CALL section_vals_val_get(hf_sub_section, "EPS_STORAGE_SCALING", r_val=real_val, error=error)
    x_data%memory_parameter%eps_storage_scaling = real_val
    IF(int_val == 0) THEN
      x_data%memory_parameter%do_all_on_the_fly = .TRUE.
    ELSE
      x_data%memory_parameter%do_all_on_the_fly = .FALSE.
    END IF
    x_data%memory_parameter%cache_size = CACHE_SIZE

    !! PERIODIC section
    hf_sub_section => section_vals_get_subs_vals(hf_section,"PERIODIC",error=error)
    CALL section_vals_get(hf_sub_section, explicit=logic_val, error=error)
    x_data%periodic_parameter%do_periodic = logic_val
    CALL section_vals_val_get(hf_sub_section, "NUMBER_OF_SHELLS", i_val=int_val, error=error)
    x_data%periodic_parameter%number_of_shells = int_val
    
    !! SCREENING section
    hf_sub_section => section_vals_get_subs_vals(hf_section,"SCREENING",error=error)
    CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ", r_val=real_val, error=error)
    x_data%screening_parameter%eps_schwarz = real_val
    CALL section_vals_val_get(hf_sub_section, "EPS_BOX", r_val=real_val, error=error)
    x_data%screening_parameter%eps_box = real_val 
    CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ_FORCES", r_val=real_val, error=error)
    x_data%screening_parameter%eps_schwarz_forces = real_val 
    CALL section_vals_val_get(hf_sub_section, "EPS_BOX_FORCES", r_val=real_val, error=error)
    x_data%screening_parameter%eps_box_forces = real_val 
    CALL section_vals_val_get(hf_sub_section, "SCREEN_P_FORCES", l_val=logic_val, error=error)
    x_data%screening_parameter%do_p_screening_forces = logic_val 
    CALL section_vals_val_get(hf_sub_section, "SCREEN_ON_INITIAL_P", l_val=logic_val, error=error)
    x_data%screening_parameter%do_initial_p_screening = logic_val 
    CALL section_vals_val_get(hf_sub_section, "P_SCREEN_CORRECTION_FACTOR", r_val=real_val, error=error)
    IF( real_val == 0.0_dp ) THEN
      x_data%screening_parameter%do_init_p_screen_corr = .FALSE. 
    ELSE
      x_data%screening_parameter%do_init_p_screen_corr = .TRUE. 
    END IF
    x_data%screening_parameter%initial_p_screening_correction = real_val 

    !! INTERACTION_POTENTIAL section
    hf_sub_section => section_vals_get_subs_vals(hf_section,"INTERACTION_POTENTIAL",error=error)
    CALL section_vals_val_get(hf_sub_section, "POTENTIAL_TYPE", i_val=int_val, error=error)
    x_data%potential_parameter%potential_type = int_val
    CALL section_vals_val_get(hf_sub_section, "OMEGA", r_val=real_val, error=error)
    x_data%potential_parameter%omega = real_val 

    !! LOAD_BALANCE section
    hf_sub_section => section_vals_get_subs_vals(hf_section,"LOAD_BALANCE",error=error)
    CALL section_vals_val_get(hf_sub_section, "NBINS", i_val=int_val, error=error)
    x_data%load_balance_parameter%nbins = int_val
    CALL section_vals_val_get(hf_sub_section, "MOVES", i_val=int_val, error=error)
    x_data%load_balance_parameter%moves = int_val
    CALL section_vals_val_get(hf_sub_section, "CYCLES", i_val=int_val, error=error)
    x_data%load_balance_parameter%cycles = int_val
    CALL section_vals_val_get(hf_sub_section, "REDUCTION_STEPS", i_val=int_val, error=error)
    x_data%load_balance_parameter%reduction_steps = int_val
    CALL section_vals_val_get(hf_sub_section, "REDUCTION_FACTOR", r_val=real_val, error=error)
    x_data%load_balance_parameter%reduction_factor = real_val
    CALL section_vals_val_get(hf_sub_section, "TEMPERATURE", r_val=real_val, error=error)
    x_data%load_balance_parameter%temperature = real_val
    CALL section_vals_val_get(hf_sub_section, "TOLERANCE", r_val=real_val, error=error)
    x_data%load_balance_parameter%tolerance = real_val
    CALL section_vals_val_get(hf_sub_section, "SWAP_PROBABILITY", r_val=real_val, error=error)
    x_data%load_balance_parameter%swap_probability = real_val

    IF(x_data%periodic_parameter%do_periodic) THEN
      hf_pbc_section => section_vals_get_subs_vals(hf_section,"PERIODIC",error=error)
      CALL section_vals_val_get(hf_pbc_section,"NUMBER_OF_SHELLS",i_val=pbc_shells,error=error)
      ALLOCATE(x_data%neighbor_cells(1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL x_hf_create_neighbor_cells(x_data,pbc_shells, error=error)
      ALLOCATE(x_data%x_screening_buffer_periodic(SIZE(x_data%neighbor_cells),&
               max_set,max_set,natom,natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      ALLOCATE(x_data%x_screening_buffer(max_set,max_set,natom,natom),STAT=stat)
    END IF
    
    ALLOCATE(x_data%is_assoc_atomic_block(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%distribution_forces(1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(x_data%distribution_energy(1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%far_field_boxes(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF( x_data%screening_parameter%do_initial_p_screening) THEN
      ALLOCATE(x_data%initial_p(max_set,max_set,natom,natom),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      ALLOCATE(x_data%initial_p(1,1,1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    x_data%maxval_cache%element_counter = 1
    ALLOCATE(x_data%maxval_container,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(x_data%maxval_container%first,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    x_data%maxval_container%first%prev => NULL()
    x_data%maxval_container%first%next => NULL()
    x_data%maxval_container%current => x_data%maxval_container%first
    x_data%maxval_container%current%data = 0
    x_data%maxval_container%element_counter = 1
    x_data%maxval_container%list_counter = 1

    ALLOCATE(x_data%integral_containers(64), STAT=stat)
    DO i=1,64
      x_data%integral_caches(i)%element_counter = 1
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(x_data%integral_containers(i)%first,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      x_data%integral_containers(i)%first%prev => NULL()
      x_data%integral_containers(i)%first%next => NULL()
      x_data%integral_containers(i)%current => x_data%integral_containers(i)%first
      x_data%integral_containers(i)%current%data = 0
      x_data%integral_containers(i)%element_counter = 1
      x_data%integral_containers(i)%list_counter = 1
    END DO

  END SUBROUTINE x_hf_create 
!***************************************************************************

!!***** x_hf/x_hf_release [1.0] *
!!
!!   NAME
!!     x_hf_release
!!
!!   FUNCTION
!!     - This routine deallocates all data structures
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - x_data: contains all relevant data structures for hfx runs
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_release(x_data,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF( x_data%periodic_parameter%do_periodic ) THEN
      DEALLOCATE(x_data%x_screening_buffer_periodic,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(x_data%neighbor_cells,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      DEALLOCATE(x_data%x_screening_buffer,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(x_data%distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data%distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(x_data%far_field_boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(x_data%initial_p,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(x_data%is_assoc_atomic_block,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    !! Deallocate containers
    CALL x_hf_init_container(x_data%maxval_container, error)
    DEALLOCATE(x_data%maxval_container%first,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data%maxval_container,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,64   
      CALL x_hf_init_container(x_data%integral_containers(i), error)
      DEALLOCATE(x_data%integral_containers(i)%first,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(x_data%integral_containers,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    DEALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE x_hf_release
!***************************************************************************

!!***** x_hf/x_hf_set_distr_energy [1.0] *
!!
!!   NAME
!!     x_hf_set_distr_energy
!!
!!   FUNCTION
!!     - This routine stores the data obtained from the load balance routine 
!!       for the energy
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - ptr_to_distr_energy: contains data to store
!!     - x_data: contains all relevant data structures for hfx runs
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_set_distr_energy(ptr_to_distr,x_data, error)
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_set_distr_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%distribution_energy(SIZE(ptr_to_distr)))
    x_data%distribution_energy = ptr_to_distr

  END SUBROUTINE x_hf_set_distr_energy
!***************************************************************************

!!***** x_hf/x_hf_set_distr_forces [1.0] *
!!
!!   NAME
!!     x_hf_set_distr_forces
!!
!!   FUNCTION
!!     - This routine stores the data obtained from the load balance routine 
!!       for the forces
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - ptr_to_distr_energy: contains data to store
!!     - x_data: contains all relevant data structures for hfx runs
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_set_distr_forces(ptr_to_distr, x_data, error)
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_set_distr_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   
    ALLOCATE(x_data%distribution_forces(SIZE(ptr_to_distr)))
    x_data%distribution_forces = ptr_to_distr

  END SUBROUTINE x_hf_set_distr_forces
!***************************************************************************

!!***** x_hf/x_hf_create_neighbor_cells [1.0] *
!!
!!   NAME
!!     x_hf_create_neighbor_cells
!!
!!   FUNCTION
!!     - This routine computes the neighbor cells that are taken into account
!!       in periodic runs 
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - x_data: contains all relevant data structures for hfx runs
!!     - pbc_shells: number of shells taken into account
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_create_neighbor_cells(x_data, pbc_shells, error)
    TYPE(x_hf_type), POINTER                 :: x_data
    INTEGER, INTENT(IN)                      :: pbc_shells
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_create_neighbor_cells', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m(3), stat, &
                                                total_number_of_cells
    LOGICAL                                  :: failure

    failure = .FALSE.
    total_number_of_cells = 0
    DO i = 0,pbc_shells
      total_number_of_cells = total_number_of_cells + count_cells(i)
    END DO
    
    DEALLOCATE(x_data%neighbor_cells,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    m = 0
    i = 1
    DO WHILE(SUM(m**2)<=pbc_shells)
       x_data%neighbor_cells(i)%cell = m
       CALL next_image_cell(m)
       i=i+1
    ENDDO
  END SUBROUTINE x_hf_create_neighbor_cells
!***************************************************************************

!!***** x_hf/count_cells [1.0] *
!!
!!   NAME
!!     count_cells
!!
!!   FUNCTION
!!     - Auxiliary function for creating periodic neighbor cells
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - shell: number of shells in each coordinate direction
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

   FUNCTION count_cells(shell)
    INTEGER, INTENT(IN)                      :: shell
    INTEGER                                  :: count_cells

    INTEGER                                  :: i, j, k

    count_cells=0
    DO i=-shell,shell
      DO j=-shell,shell
        DO k=-shell,shell
          IF ((i**2+j**2+k**2==shell)) count_cells=count_cells+1
        ENDDO
      ENDDO
    ENDDO
  END FUNCTION count_cells
!***************************************************************************

!!***** x_hf/next_image_cell [1.0] *
!!
!!   NAME
!!     next_image_cell
!!
!!   FUNCTION
!!     - Auxiliary function for creating periodic neighbor cells
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - shell: shell for which we need the next, stored as number of shifts
!!              in each direction
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE next_image_cell(m)
    INTEGER                                  :: m(3)

    INTEGER                                  :: i, j, k, shell
    LOGICAL                                  :: found

    found=.FALSE.
    shell=SUM(m**2)
    outer: DO
           DO i=-shell,shell
           DO j=-shell,shell
    inner: DO k=-shell,shell
              IF (.NOT.(i**2+j**2+k**2==shell)) CYCLE inner
              IF (found) THEN
                  m=(/i,j,k/)
                  EXIT outer
              ENDIF
              IF (ALL(M.EQ.(/i,j,k/))) found=.TRUE.
           ENDDO inner
           ENDDO
           ENDDO
           shell=shell+1
           ENDDO outer
  END SUBROUTINE next_image_cell
!***************************************************************************

!!***** x_hf/x_hf_add_single_cache_element [1.0] *
!!
!!   NAME
!!     x_hf_add_single_cache_element
!!
!!   FUNCTION
!!     - This routine adds an int_8 value to a cache. If the cache is full
!!       a compression routine is invoked and the cache is cleared
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - value: value to be added to the cache
!!     - nbits: number of bits to be stored
!!     - cache: cache to which we want to add
!!     - container: container that contains the compressed elements
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_add_single_cache_element(value, nbits, cache, container, error)
    INTEGER(int_8)                           :: value
    INTEGER                                  :: nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER(int_8)                           :: int_val

    int_val = value + ISHFT(one,nbits-1)
    
    IF(cache%element_counter /= CACHE_SIZE ) THEN
      cache%data(cache%element_counter) = int_val
      cache%element_counter = cache%element_counter + 1
    ELSE
      cache%data(CACHE_SIZE) = int_val
      CALL x_hf_compress_cache(cache%data(1), container, CACHE_SIZE, nbits, error)
      cache%element_counter = 1
    END IF
  END SUBROUTINE x_hf_add_single_cache_element
!***************************************************************************

!!***** x_hf/x_hf_compress_cache [1.0] *
!!
!!   NAME
!!     x_hf_compress_cache
!!
!!   FUNCTION
!!     - This routine compresses a full cache and stores its values
!!       in a container. If necessary, a new list entry is allocated 
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - full_array: values from the cache
!!     - container: linked list, that stores the compressed values
!!     - nints: number of int_8 values to be compressed
!!     - nbits: number of bits to be stored
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_compress_cache(full_array, container, nints, nbits, error)
    INTEGER(int_8)                           :: full_array(*)
    TYPE(x_hf_container_type)                :: container
    INTEGER                                  :: nints, nbits
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_compress_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: end_idx, increment_counter, &
                                                start_idx, stat, &
                                                tmp_elements, tmp_nints
    LOGICAL                                  :: failure

    failure = .FALSE.

    start_idx = container%element_counter
    increment_counter = (nbits*CACHE_SIZE + 63)/64
    end_idx = start_idx + increment_counter - 1
    IF( end_idx<CACHE_SIZE) THEN
      CALL ints2bits_specific(nbits,CACHE_SIZE,container%current%data(start_idx),full_array(1))
      container%element_counter = container%element_counter + increment_counter
    ELSE
      !! We have to fill the container first with the remaining number of bits
      tmp_elements = CACHE_SIZE - start_idx + 1
      tmp_nints = (tmp_elements*64)/nbits
      CALL ints2bits_specific(nbits,tmp_nints,container%current%data(start_idx),full_array(1))
      !! Allocate new list entry
      ALLOCATE(container%current%next, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      container%list_counter = container%list_counter + 1
      container%current%next%next => NULL()
      container%current => container%current%next
      !! compress remaining ints
      CALL ints2bits_specific(nbits, CACHE_SIZE-tmp_nints, container%current%data(1), full_array(tmp_nints+1))
      container%element_counter = 1 + (nbits*(CACHE_SIZE-tmp_nints) + 63)/64 
   END IF
  END SUBROUTINE x_hf_compress_cache
!***************************************************************************

!!***** x_hf/x_hf_get_single_cache_element [1.0] *
!!
!!   NAME
!!     x_hf_get_single_cache_element
!!
!!   FUNCTION
!!     - This routine returns an int_8 value from a cache. If the cache is empty
!!       a decompression routine is invoked and the cache is refilled with decompressed
!!       values from a container
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - value: value to be retained from the cache
!!     - nbits: number of bits with which the value has been compressed
!!     - cache: cache from which we get the value
!!     - container: container that contains the compressed elements
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_get_single_cache_element(value, nbits, cache, container)
    INTEGER(int_8)                           :: value
    INTEGER                                  :: nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container

    IF(cache%element_counter /= CACHE_SIZE ) THEN
      value = cache%data(cache%element_counter) 
      cache%element_counter = cache%element_counter + 1
    ELSE
      value = cache%data(CACHE_SIZE)
      CALL x_hf_decompress_cache(cache%data(1), container, CACHE_SIZE, nbits)
      cache%element_counter = 1 
    END IF

    value = value - ISHFT(one,nbits-1)

  END SUBROUTINE x_hf_get_single_cache_element
!***************************************************************************

!!***** x_hf/x_hf_decompress_cache [1.0] *
!!
!!   NAME
!!     x_hf_decompress_cache
!!
!!   FUNCTION
!!     - This routine decompresses data from a container in order to fill
!!       a cache.
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - full_array: values to be retained from container
!!     - container: linked list, that stores the compressed values
!!     - nints: number of int_8 values to be decompressed
!!     - nbits: number of bits with which the values have been stored
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_decompress_cache(full_array, container, nints, nbits)
    INTEGER(int_8)                           :: full_array(*)
    TYPE(x_hf_container_type)                :: container
    INTEGER                                  :: nints, nbits

    INTEGER                                  :: end_idx, increment_counter, &
                                                start_idx, tmp_elements, &
                                                tmp_nints

    start_idx = container%element_counter
    increment_counter = (nbits*CACHE_SIZE + 63)/64
    end_idx = start_idx + increment_counter - 1
    IF( end_idx<CACHE_SIZE) THEN
      CALL bits2ints_specific(nbits,CACHE_SIZE,container%current%data(start_idx),full_array(1))
      container%element_counter = container%element_counter + increment_counter
    ELSE
      !! We have to fill the container first with the remaining number of bits
      tmp_elements = CACHE_SIZE - start_idx + 1
      tmp_nints = (tmp_elements*64)/nbits
      CALL bits2ints_specific(nbits,tmp_nints,container%current%data(start_idx),full_array(1))
      container%current => container%current%next
      container%list_counter = container%list_counter + 1
      !! decompress remaining ints
      CALL bits2ints_specific(nbits, CACHE_SIZE-tmp_nints, container%current%data(1), full_array(tmp_nints+1))
      container%element_counter = 1 + (nbits*(CACHE_SIZE-tmp_nints) + 63)/64 
    END IF
  END SUBROUTINE x_hf_decompress_cache
!***************************************************************************

!!***** x_hf/x_hf_reset_cache_and_container[1.0] *
!!
!!   NAME
!!     x_hf_reset_cache_and_container
!!
!!   FUNCTION
!!     - This routine resets the containers list pointer to the first element and
!!       moves the element counters of container and cache to the beginning
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - cache: cache from which we get the value
!!     - container: container that contains the compressed elements
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_reset_cache_and_container(cache, container)
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container

    cache%element_counter = 1
    container%current => container%first
    container%element_counter = 1
    container%list_counter = 1
  END SUBROUTINE x_hf_reset_cache_and_container
!***************************************************************************

!!***** x_hf/x_hf_init_container[1.0] *
!!
!!   NAME
!!     x_hf_reset_init_container
!!
!!   FUNCTION
!!     - This routine deletes all list entries in a container in order to
!!       deallocate the memory.
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - container: container that contains the compressed elements
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_init_container(container, error)
    TYPE(x_hf_container_type)                :: container
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_init_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(x_hf_container_node), POINTER       :: current, next

    failure =.FALSE.

    !! DEALLOCATE memory

    current => container%first
    DO WHILE(ASSOCIATED(current))
      next => current%next
      DEALLOCATE(current)
      current=>next
    END DO

    !! Allocate first list entry, init members
    ALLOCATE(container%first,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    container%first%prev => NULL()
    container%first%next => NULL()
    container%current => container%first
    container%current%data = 0
    container%element_counter = 1
    container%list_counter = 1

  END SUBROUTINE x_hf_init_container
!***************************************************************************

!!***** x_hf/x_hf_decompress_first_cache[1.0] *
!!
!!   NAME
!!     x_hf_decompress_first_cache
!!
!!   FUNCTION
!!     - This routine decompresses the first bunch of data in a container and
!!       copies them into a cache
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - nbits: number of bits with which the data has been stored
!!     - cache: array where we want to decompress the data
!!     - container: container that contains the compressed elements
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_decompress_first_cache(nbits, cache, container)
    INTEGER                                  :: nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container

    CALL x_hf_decompress_cache(cache%data(1), container, CACHE_SIZE, nbits)
    cache%element_counter = 1
  END SUBROUTINE x_hf_decompress_first_cache
!***************************************************************************

!!***** x_hf/x_hf_flush_last_cache[1.0] *
!!
!!   NAME
!!     x_hf_flush_last_cache
!!
!!   FUNCTION
!!     - This routine compresses the last probably not yet compressed cache into
!!       a container
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - nbits: number of bits with which the data has been stored
!!     - cache: array where we want to decompress the data
!!     - container: container that contains the compressed elements
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_flush_last_cache(nbits, cache, container, error)
    INTEGER                                  :: nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL x_hf_compress_cache(cache%data(1), container, CACHE_SIZE, nbits, error)
  END SUBROUTINE x_hf_flush_last_cache
!***************************************************************************

!!***** x_hf/x_hf_add_mult_cache_elements [1.0] *
!!
!!   NAME
!!     x_hf_add_mult_cache_elements
!!
!!   FUNCTION
!!     - This routine adds an a few real values to a cache. If the cache is full
!!       a compression routine is invoked and the cache is cleared
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - values: values to be added to the cache
!!     - nbits: number of bits to be stored
!!     - cache: cache to which we want to add
!!     - container: container that contains the compressed elements
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_add_mult_cache_elements(values, nints, nbits, cache, container, eps_schwarz, pmax_entry, error)
    REAL(dp)                                 :: values(*)
    INTEGER, INTENT(IN)                      :: nints, nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container
    REAL(dp), INTENT(IN)                     :: eps_schwarz, pmax_entry
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: end_idx, i, start_idx, &
                                                tmp_elements
    INTEGER(int_8)                           :: shift
    REAL(dp)                                 :: eps_schwarz_inv, factor

!! ASSUME nints< CACHE_SIZE

    eps_schwarz_inv = 1.0_dp/eps_schwarz
    factor = eps_schwarz/pmax_entry

    shift = ISHFT(one,nbits-1)

    start_idx = cache%element_counter
    end_idx = start_idx + nints-1
    IF(end_idx< CACHE_SIZE) THEN
      
      DO i=1,nints
        values(i) = values(i)*pmax_entry
        IF(ABS(values(i)) > eps_schwarz) THEN
          cache%data(i+start_idx-1)= ANINT(values(i)*eps_schwarz_inv, dp) + shift
          values(i) = ANINT(values(i)*eps_schwarz_inv, dp)*factor
        ELSE
          values(i) = 0.0_dp
          cache%data(i+start_idx-1) = shift
        END IF
      END DO
      cache%element_counter = end_idx + 1
    ELSE
      tmp_elements = CACHE_SIZE - start_idx + 1    
      DO i=1,tmp_elements
        values(i) = values(i)*pmax_entry
        IF(ABS(values(i)) > eps_schwarz) THEN
          cache%data(i+start_idx-1)= ANINT(values(i)*eps_schwarz_inv,dp) + shift
          values(i) = ANINT(values(i)*eps_schwarz_inv,dp)*factor
        ELSE
          values(i) = 0.0_dp
          cache%data(i+start_idx-1) = shift
        END IF
      END DO
      CALL x_hf_compress_cache(cache%data(1), container, CACHE_SIZE, nbits, error)
      DO i=tmp_elements+1, nints
        values(i) = values(i)*pmax_entry
        IF(ABS(values(i)) > eps_schwarz) THEN
          cache%data(i-tmp_elements)= ANINT(values(i)*eps_schwarz_inv,dp) + shift
          values(i) = ANINT(values(i)*eps_schwarz_inv,dp)*factor
        ELSE
          values(i) = 0.0_dp
          cache%data(i-tmp_elements) = shift
        END IF
      END DO
      cache%element_counter = nints-tmp_elements + 1 
    END IF
  END SUBROUTINE x_hf_add_mult_cache_elements
!***************************************************************************

!!***** x_hf/x_hf_get_mult_cache_elements [1.0] *
!!
!!   NAME
!!     x_hf_get_mult_cache_elements
!!
!!   FUNCTION
!!     - This routine returns a bunch real values from a cache. If the cache is empty
!!       a decompression routine is invoked and the cache is refilled with decompressed
!!       values from a container
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - values: value to be retained from the cache
!!     - nints: number of values to be retained
!!     - nbits: number of bits with which the value has been compressed
!!     - cache: cache from which we get the value
!!     - container: container that contains the compressed elements
!!     - eps_schwarz: threshold for storage
!!     - pmax_entry: multiplication factor for values
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE x_hf_get_mult_cache_elements(values, nints, nbits, cache, container, eps_schwarz, pmax_entry, error)
    REAL(dp)                                 :: values(*)
    INTEGER, INTENT(IN)                      :: nints, nbits
    TYPE(x_hf_cache_type)                    :: cache
    TYPE(x_hf_container_type)                :: container
    REAL(dp), INTENT(IN)                     :: eps_schwarz, pmax_entry
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: end_idx, i, start_idx, &
                                                tmp_elements
    INTEGER(int_8)                           :: shift
    REAL(dp)                                 :: eps_schwarz_inv, factor

!! ASSUME nints< CACHE_SIZE

    eps_schwarz_inv = 1.0_dp/eps_schwarz
    factor = eps_schwarz/pmax_entry

    shift = ISHFT(one,nbits-1)

    start_idx = cache%element_counter
    end_idx = start_idx + nints-1

    IF(end_idx< CACHE_SIZE) THEN
      DO i=1,nints
        values(i) = factor*REAL(cache%data(i+start_idx-1) - shift,dp)
      END DO
      cache%element_counter = end_idx + 1
    ELSE
      tmp_elements = CACHE_SIZE - start_idx + 1
      DO i=1,tmp_elements
        values(i) = factor*REAL(cache%data(i+start_idx-1)-shift,dp)
      END DO
      CALL x_hf_decompress_cache(cache%data(1), container, CACHE_SIZE, nbits)
      DO i=tmp_elements+1,nints
        values(i) = factor*REAL(cache%data(i-tmp_elements)-shift,dp)
      END DO
      cache%element_counter = nints - tmp_elements + 1
    END IF
  END SUBROUTINE x_hf_get_mult_cache_elements
!***************************************************************************

!!***** x_hf/x_hf_get_memory_info [1.0] *
!!
!!   NAME
!!     x_hf_get_memory_info
!!
!!   FUNCTION
!!     - Auxiliary function that returns the memory usage of a container
!!
!!   NOTES
!!     - 
!!
!!   INPUTS
!!     - container: container that contains the compressed elements
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  FUNCTION x_hf_get_memory_info(container, error)
    TYPE(x_hf_container_type)                :: container
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER(int_8)                           :: x_hf_get_memory_info

    INTEGER(int_8)                           :: nlist_entries

    nlist_entries = container%list_counter
    x_hf_get_memory_info = nlist_entries * CACHE_SIZE

  END FUNCTION x_hf_get_memory_info
END MODULE x_hf


