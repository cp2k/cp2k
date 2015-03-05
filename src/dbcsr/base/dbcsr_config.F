!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Configuration options for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-05
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-05
! *****************************************************************************
MODULE dbcsr_config
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_fatal_level, dbcsr_wrong_args_error, &
       external_timeset, external_timestop, timeset_interface, &
       timestop_interface
  USE kinds,                           ONLY: dp

!$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_config'


  PUBLIC :: dbcsr_set_default_config

  PUBLIC :: dbcsr_set_conf_subcomm, dbcsr_get_conf_subcomm
  PUBLIC :: dbcsr_set_conf_mm_stacksize, dbcsr_get_conf_mm_stacksize

  PUBLIC :: dbcsr_set_conf_mpi_mem, dbcsr_get_conf_mpi_mem

  PUBLIC :: dbcsr_set_conf_mm_driver, dbcsr_get_conf_mm_driver, dbcsr_get_conf_max_elements_images

  PUBLIC :: dbcsr_set_conf_max_elements_images, dbcsr_set_conf_nstacks, dbcsr_get_conf_nstacks

  PUBLIC :: dbcsr_set_conf_use_comm_thread, dbcsr_get_conf_use_comm_thread
  PUBLIC :: dbcsr_set_conf_comm_thread_load, dbcsr_get_conf_comm_thread_load
  PUBLIC :: dbcsr_set_conf_external_timing
  PUBLIC :: dbcsr_get_conf_max_ele_block, dbcsr_set_conf_max_ele_block
  

  ! These are for DBCSR internal use.
  !
  PUBLIC :: has_MPI
  PUBLIC :: has_acc
  PUBLIC :: use_subcommunicators
  PUBLIC :: mm_stack_size
  PUBLIC :: mm_driver, mm_host_driver,&
            mm_driver_blas,&
            mm_driver_matmul,&
            mm_driver_smm,&
            mm_driver_acc
  PUBLIC :: mm_async
  PUBLIC :: mm_name_blas,&
            mm_name_matmul,&
            mm_name_smm,&
            mm_name_acc
  PUBLIC :: use_comm_thread, comm_thread_load
  PUBLIC :: max_elements_per_block
  !
  PUBLIC :: is_initialized

  PUBLIC :: accdrv_avoid_after_busy, accdrv_min_flop_process,&
            accdrv_min_flop_sort, accdrv_do_inhomogenous,&
            accdrv_binning_nbins, accdrv_binning_binsize,&
            accdrv_posterior_streams, accdrv_posterior_buffers,&
            accdrv_priority_streams, accdrv_priority_buffers

  PUBLIC :: multrec_limit
  PUBLIC :: default_resize_factor

  ! First the constants are declared.
  REAL(KIND=dp), PARAMETER             :: default_resize_factor = 1.2_dp

  ! Possible drivers to use for matrix multiplications
  INTEGER, PARAMETER :: mm_driver_matmul     = 1
  INTEGER, PARAMETER :: mm_driver_blas       = 2
  INTEGER, PARAMETER :: mm_driver_smm        = 3
  INTEGER, PARAMETER :: mm_driver_acc       = 5

  CHARACTER(len=*), PARAMETER :: mm_name_blas   = "BLAS",&
                                 mm_name_matmul = "MATMUL",&
                                 mm_name_smm    = "SMM",&
                                 mm_name_acc   = "ACC"

  ! Then the capabilities are specified.

  LOGICAL, PARAMETER :: has_blas = .TRUE.

#if defined (__HAS_smm_dnn)
  LOGICAL, PARAMETER :: has_smm_gemm = .TRUE.
#else
  LOGICAL, PARAMETER :: has_smm_gemm = .FALSE.
#endif

#if defined(__HAS_smm_vec)
  LOGICAL, PARAMETER :: has_smm_vec = .TRUE.
#else
  LOGICAL, PARAMETER :: has_smm_vec = .FALSE.
#endif

#if defined (__DBCSR_ACC)
  LOGICAL, PARAMETER :: has_acc = .TRUE.
#else
  LOGICAL, PARAMETER :: has_acc = .FALSE.
#endif

#if defined (__parallel)
  LOGICAL, PARAMETER :: has_MPI = .TRUE.
#else
  LOGICAL, PARAMETER :: has_MPI = .FALSE.
#endif

  LOGICAL, PARAMETER :: has_omp_2_5 = .TRUE.


  ! Now the fallback defaults are defined.  The defaults should be set
  ! by calling the dbcsr_init_conf() subroutine.

  ! Allocates subcommunicators for process rows and columns.
  LOGICAL, SAVE :: use_subcommunicators = .TRUE.

  ! Use MPI-allocated memory.
  LOGICAL, SAVE :: use_MPI_memory = has_MPI

  ! Which driver to use for matrix multiplications.
#if defined (__DBCSR_ACC)
  INTEGER, SAVE :: mm_driver = mm_driver_acc
#else
  INTEGER, SAVE :: mm_driver = mm_driver_smm
#endif

  ! Which driver to use for matrix multiplications (on the host).
  INTEGER, SAVE :: mm_host_driver = mm_driver_smm

  ! Whether an accelerator is used.
  LOGICAL, SAVE :: mm_async = .FALSE.

  ! Stack size to use for multiplication parameters
  INTEGER, SAVE :: mm_stack_size = 1000

  ! Default blocking parameter
  INTEGER, SAVE :: max_elements_per_block = 32

  ! Maximum number of elements for each image
  INTEGER, SAVE :: max_elements_images = -1

  ! Number of stacks to use
  INTEGER, SAVE :: nm_stacks = 3
  INTEGER, SAVE :: nn_stacks = 3
  INTEGER, SAVE :: nk_stacks = 3

  ! Configuration of an MPI progress thread
  LOGICAL, SAVE :: use_comm_thread = .TRUE.
  INTEGER, SAVE :: comm_thread_load = 100

  LOGICAL, SAVE :: is_initialized = .FALSE.


  INTEGER, SAVE :: multrec_limit = 512

  INTEGER, SAVE :: accdrv_priority_streams = 4
  INTEGER, SAVE :: accdrv_priority_buffers = 40
  INTEGER, SAVE :: accdrv_posterior_streams  = 4
  INTEGER, SAVE :: accdrv_posterior_buffers  = 80
  LOGICAL, SAVE :: accdrv_avoid_after_busy = .FALSE.
  INTEGER, SAVE :: accdrv_min_flop_process = 0
  INTEGER, SAVE :: accdrv_min_flop_sort = 4000
  LOGICAL, SAVE :: accdrv_do_inhomogenous = .TRUE.
  INTEGER, SAVE :: accdrv_binning_nbins = 4096
  INTEGER, SAVE :: accdrv_binning_binsize = 16

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param driver ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_mm_driver (driver, error)
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------
!defaults

    mm_stack_size = 1000
    nm_stacks = 3
    nn_stacks = 3
    nk_stacks = 3

    IF(has_smm_vec) THEN
      mm_stack_size = 30000
      nm_stacks = 3
      nn_stacks = 3
      nk_stacks = 3
    ENDIF

    CALL dbcsr_error_set(routineN, error_handle, error)
    SELECT CASE (driver)
    CASE (mm_driver_matmul)
       mm_driver = mm_driver_matmul
    CASE (mm_driver_blas)
       CALL dbcsr_assert (has_blas,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for BLAS not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_blas
    CASE (mm_driver_smm)
       mm_driver = mm_driver_smm
    CASE (mm_driver_acc)
       CALL dbcsr_assert (has_acc,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for ACC not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_acc
       mm_stack_size = 30000
       nm_stacks = 3
       nn_stacks = 3
       nk_stacks = 3
    END SELECT
    mm_async = driver_is_async (mm_driver)
    IF (mm_async) THEN
       mm_host_driver = mm_driver_smm
    ELSE
       mm_host_driver = mm_driver
    ENDIF
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_mm_driver

! *****************************************************************************
!> \brief ...
!> \retval driver ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_mm_driver () RESULT (driver)
    INTEGER                                  :: driver

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    driver = mm_driver
  END FUNCTION dbcsr_get_conf_mm_driver


! *****************************************************************************
!> \brief ...
!> \param use_subcomms ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_subcomm (use_subcomms, error)
    LOGICAL, INTENT(IN)                      :: use_subcomms
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_subcomm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_subcommunicators = use_subcomms
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_subcomm

! *****************************************************************************
!> \brief ...
!> \retval use_subcomms ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_subcomm () RESULT (use_subcomms)
    LOGICAL                                  :: use_subcomms

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_subcomm', &
      routineP = moduleN//':'//routineN

    use_subcomms = use_subcommunicators
  END FUNCTION dbcsr_get_conf_subcomm


! *****************************************************************************
!> \brief ...
!> \param use_mpi_mem ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_mpi_mem (use_mpi_mem, error)
    LOGICAL, INTENT(IN)                      :: use_mpi_mem
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_mpi_mem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_MPI_memory = use_mpi_mem
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_mpi_mem

! *****************************************************************************
!> \brief ...
!> \retval use_mpi_mem ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_mpi_mem () RESULT (use_mpi_mem)
    LOGICAL                                  :: use_mpi_mem

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mpi_mem', &
      routineP = moduleN//':'//routineN

    use_mpi_mem = use_MPI_memory
  END FUNCTION dbcsr_get_conf_mpi_mem

! *****************************************************************************
!> \brief ...
!> \param use_mm_stacksize ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_mm_stacksize (use_mm_stacksize, error)
    INTEGER, INTENT(IN)                      :: use_mm_stacksize
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_mm_stacksize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    mm_stack_size = use_mm_stacksize
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_mm_stacksize

! *****************************************************************************
!> \brief ...
!> \retval use_mm_stacksize ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_mm_stacksize () RESULT (use_mm_stacksize)
    INTEGER                                  :: use_mm_stacksize

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mm_stacksize', &
      routineP = moduleN//':'//routineN

    use_mm_stacksize = mm_stack_size
  END FUNCTION dbcsr_get_conf_mm_stacksize

! *****************************************************************************
!> \brief ...
!> \param max_ele_block ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_max_ele_block (max_ele_block, error)
    INTEGER, INTENT(IN)                      :: max_ele_block
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_max_ele_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    max_elements_per_block = max_ele_block
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_max_ele_block

! *****************************************************************************
!> \brief ...
!> \retval max_ele_block ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_max_ele_block () RESULT (max_ele_block)
    INTEGER                                  :: max_ele_block

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_max_ele_block', &
      routineP = moduleN//':'//routineN

    max_ele_block = max_elements_per_block
  END FUNCTION dbcsr_get_conf_max_ele_block

! *****************************************************************************
!> \brief ...
!> \param new_max_elements_images ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_max_elements_images (new_max_elements_images, error)
    INTEGER, INTENT(IN)                      :: new_max_elements_images
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_conf_max_elements_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    max_elements_images = new_max_elements_images
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_max_elements_images

! *****************************************************************************
!> \brief ...
!> \retval dbcsr_get_conf_max_elements_images ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_max_elements_images ()
    INTEGER :: dbcsr_get_conf_max_elements_images

    dbcsr_get_conf_max_elements_images=max_elements_images
  END FUNCTION dbcsr_get_conf_max_elements_images

! *****************************************************************************
!> \brief ...
!> \param nstacks ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_nstacks (nstacks, error)
    INTEGER, INTENT(IN)                      :: nstacks(:)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_nstacks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    CALL dbcsr_assert (SIZE(nstacks) ,"EQ", 3,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "ntacks size must be 3", __LINE__, error=error)
    nm_stacks        = nstacks(1)
    nn_stacks        = nstacks(2)
    nk_stacks        = nstacks(3)
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_nstacks

! *****************************************************************************
!> \brief ...
!> \param nstacks ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_get_conf_nstacks (nstacks, error)
    INTEGER, INTENT(OUT)                     :: nstacks(:)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_nstacks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    CALL dbcsr_assert (SIZE(nstacks) ,"EQ", 3,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "ntacks size must be 3", __LINE__, error=error)
    nstacks(1) = nm_stacks
    nstacks(2) = nn_stacks
    nstacks(3) = nk_stacks
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_conf_nstacks

! *****************************************************************************
!> \brief ...
!> \param comm_thread ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_use_comm_thread (comm_thread, error)
    LOGICAL, INTENT(IN)                      :: comm_thread
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_conf_use_comm_thread', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_comm_thread = comm_thread
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_use_comm_thread


! *****************************************************************************
!> \brief ...
!> \retval comm_thread ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_use_comm_thread () RESULT (comm_thread)
    LOGICAL                                  :: comm_thread

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_get_conf_use_comm_thread', &
      routineP = moduleN//':'//routineN

    comm_thread = use_comm_thread
  END FUNCTION dbcsr_get_conf_use_comm_thread

! *****************************************************************************
!> \brief ...
!> \param load ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_comm_thread_load (load, error)
    INTEGER, INTENT(IN)                      :: load
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_conf_comm_thread_load', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nthreads

    CALL dbcsr_error_set(routineN, error_handle, error)
    nthreads = 1
    !$ nthreads = OMP_GET_MAX_THREADS()
    comm_thread_load = MAX(0,90-(30*nthreads)/8)
    IF(load >= 0)&
       comm_thread_load = load
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_comm_thread_load


! *****************************************************************************
!> \brief ...
!> \retval load ...
! *****************************************************************************
  FUNCTION dbcsr_get_conf_comm_thread_load () RESULT (load)
    INTEGER                                  :: load

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_get_conf_comm_thread_load', &
      routineP = moduleN//':'//routineN

    load = comm_thread_load
  END FUNCTION dbcsr_get_conf_comm_thread_load


! *****************************************************************************
!> \brief ...
!> \param driver ...
!> \retval driver_is_async ...
! *****************************************************************************
  PURE FUNCTION driver_is_async (driver)
    INTEGER, INTENT(IN)                      :: driver
    LOGICAL                                  :: driver_is_async

    driver_is_async = .FALSE.
    SELECT CASE (driver)
       CASE (mm_driver_acc)
          driver_is_async = .TRUE.
    END SELECT
  END FUNCTION driver_is_async

! *****************************************************************************
!> \brief Set the callback routines for the external timing hooks.
!> \param timeset ...
!> \param timestop ...
!> \param error ...
!> \date    2012-12-20
!> \author  Ole Schuett
!>
!> \note <b>Modification history:</b>
!> \note - Created 2012-12-20
! *****************************************************************************
  SUBROUTINE dbcsr_set_conf_external_timing(timeset, timestop, error)
    PROCEDURE (timeset_interface)            :: timeset
    PROCEDURE (timestop_interface)           :: timestop
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

! no call to dbcsr_error_set / dbcsr_error_stop here,
! otherwise thie would lead to mismatched calls to timeset/timestop

    external_timeset  => timeset
    external_timestop => timestop
                
  END SUBROUTINE dbcsr_set_conf_external_timing

! *****************************************************************************
!> \brief Configure DBCSR driver with default values.
!> \param[inout] error        dbcsr error
!>
!> Basically it reads the default values and reset them.
! *****************************************************************************
  SUBROUTINE dbcsr_set_default_config(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_default_config', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm_thread_load, &
                                                error_handler, max_ele_block, &
                                                mm_driver, mm_ss
    INTEGER, DIMENSION(3)                    :: nstacks
    LOGICAL                                  :: use_comm_thread, use_mpi_mem, &
                                                use_subcomms

    CALL dbcsr_error_set (routineN, error_handler, error=error)

    use_subcomms = dbcsr_get_conf_subcomm()
    use_mpi_mem = dbcsr_get_conf_mpi_mem()
    mm_ss = -1
    mm_driver = dbcsr_get_conf_mm_driver()
    max_ele_block = dbcsr_get_conf_max_ele_block()
    nstacks(1:3) = -1
    use_comm_thread = dbcsr_get_conf_use_comm_thread()
    comm_thread_load = -1

    CALL dbcsr_set_conf_mm_driver (mm_driver, error=error)
    CALL dbcsr_set_conf_subcomm (use_subcomms, error=error)
    CALL dbcsr_set_conf_mpi_mem (use_mpi_mem, error=error)
    CALL dbcsr_set_conf_max_ele_block (max_ele_block, error=error)
    IF (mm_ss > 0) CALL dbcsr_set_conf_mm_stacksize (mm_ss, error=error)
    IF (nstacks(1) > 0) CALL dbcsr_set_conf_nstacks(nstacks, error=error)

    CALL dbcsr_set_conf_use_comm_thread(use_comm_thread, error=error)
    CALL dbcsr_set_conf_comm_thread_load(comm_thread_load, error=error)

    CALL dbcsr_error_stop(error_handler, error=error)

  END SUBROUTINE dbcsr_set_default_config

END MODULE dbcsr_config
