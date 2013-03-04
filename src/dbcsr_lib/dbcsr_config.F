!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
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

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_config'


  PUBLIC :: dbcsr_init_conf
  PUBLIC :: dbcsr_set_conf_subcomm, dbcsr_get_conf_subcomm
  PUBLIC :: dbcsr_set_conf_combtypes, dbcsr_get_conf_combtypes
  PUBLIC :: dbcsr_set_conf_mm_stacksize, dbcsr_get_conf_mm_stacksize

  PUBLIC :: dbcsr_set_conf_mpi_mem, dbcsr_get_conf_mpi_mem
  PUBLIC :: dbcsr_set_conf_cuda_mem, dbcsr_get_conf_cuda_mem

  PUBLIC :: dbcsr_set_conf_mm_driver, dbcsr_get_conf_mm_driver

  PUBLIC :: dbcsr_set_conf_nstacks, dbcsr_get_conf_nstacks

  PUBLIC :: dbcsr_set_conf_use_comm_thread, dbcsr_get_conf_use_comm_thread
  PUBLIC :: dbcsr_set_conf_comm_thread_load, dbcsr_get_conf_comm_thread_load
  PUBLIC :: dbcsr_set_conf_external_timing
  
  PUBLIC :: driver_is_async

  ! These are for DBCSR internal use.
  !
  PUBLIC :: has_MPI
  PUBLIC :: has_blas,&
            has_smm_gemm,&
            has_plasma,&
            has_cuda
  PUBLIC :: has_omp_2_5
  PUBLIC :: use_subcommunicators
  PUBLIC :: use_combined_types
  PUBLIC :: use_MPI_memory
  PUBLIC :: use_CUDA_host_pinned_memory
  PUBLIC :: mm_stack_size
  PUBLIC :: mm_driver, mm_host_driver,&
            mm_driver_blas,&
            mm_driver_matmul,&
            mm_driver_smm,&
            mm_driver_plasma,&
            mm_driver_cuda
  PUBLIC :: mm_async
  PUBLIC :: mm_name_blas,&
            mm_name_matmul,&
            mm_name_smm,&
            mm_name_plasma,&
            mm_name_cuda
  PUBLIC :: nm_stacks, nn_stacks, nk_stacks, nstackbuffers
  PUBLIC :: use_comm_thread, comm_thread_load

  !
  PUBLIC :: is_configured


  ! First the constants are declared.

  ! Possible drivers to use for matrix multiplications
  INTEGER, PARAMETER :: mm_driver_matmul     = 1
  INTEGER, PARAMETER :: mm_driver_blas       = 2
  INTEGER, PARAMETER :: mm_driver_smm        = 3
  INTEGER, PARAMETER :: mm_driver_plasma     = 4
  INTEGER, PARAMETER :: mm_driver_cuda       = 5

  CHARACTER(len=*), PARAMETER :: mm_name_blas   = "BLAS",&
                                 mm_name_matmul = "MATMUL",&
                                 mm_name_smm    = "SMM",&
                                 mm_name_plasma = "PLASMA",&
                                 mm_name_cuda   = "CUDA"

  ! Then the capabilities are specified.

  LOGICAL, PARAMETER :: has_blas = .TRUE.

#if defined (__HAS_smm_dnn)
  LOGICAL, PARAMETER :: has_smm_gemm = .TRUE.
#else
  LOGICAL, PARAMETER :: has_smm_gemm = .FALSE.
#endif

#if defined(__PLASMA)
  LOGICAL, PARAMETER :: has_plasma = .TRUE.
#else
  LOGICAL, PARAMETER :: has_plasma = .FALSE.
#endif

#if defined (__DBCSR_CUDA)
  LOGICAL, PARAMETER :: has_cuda = .TRUE.
#else
  LOGICAL, PARAMETER :: has_cuda = .FALSE.
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

  ! Use combined data types for MPI transfers.
  LOGICAL, SAVE :: use_combined_types = .FALSE.

  ! Use MPI-allocated memory.
  LOGICAL, SAVE :: use_MPI_memory = has_MPI

  ! Use CUDA host-pinned memory.
  LOGICAL, SAVE :: use_CUDA_host_pinned_memory = .FALSE.

  ! Which driver to use for matrix multiplications.
  INTEGER, SAVE :: mm_driver = mm_driver_smm

  ! Which driver to use for matrix multiplications (on the host).
  INTEGER, SAVE :: mm_host_driver = mm_driver_smm

  ! Whether an accelerator such as a GPU is used.
  LOGICAL, SAVE :: mm_async = .FALSE.

  ! Stack size to use for multiplication parameters
  INTEGER, SAVE :: mm_stack_size = 1000
  
  ! Number of stacks to use
  INTEGER, SAVE :: nm_stacks = 0
  INTEGER, SAVE :: nn_stacks = 0
  INTEGER, SAVE :: nk_stacks = 0
  INTEGER, SAVE :: nstackbuffers = 1
  
  ! Configuration of an MPI progress thread
  LOGICAL, SAVE :: use_comm_thread = .TRUE.
  INTEGER, SAVE :: comm_thread_load = 100

  LOGICAL, SAVE :: is_configured = .FALSE.


CONTAINS


  SUBROUTINE dbcsr_init_conf (error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_conf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! subcomm are crucial for make_images for large (>4096) number of MPI tasks
    use_subcommunicators = .TRUE.
    use_combined_types = .FALSE.
    use_MPI_memory = has_MPI
    use_CUDA_host_pinned_memory = .FALSE.

    !
    mm_driver = mm_driver_smm
    mm_stack_size = 1000
    IF (has_cuda) THEN
       mm_driver = mm_driver_cuda
       mm_stack_size = 30000
       nm_stacks = 2
       nn_stacks = 2
       nk_stacks = 2
       nstackbuffers = 4
    ENDIF
    !
    use_comm_thread = .TRUE.
    comm_thread_load = 100
    mm_async = driver_is_async (mm_driver)
    IF (mm_async) THEN
       mm_host_driver = mm_driver_smm
    ELSE
       mm_host_driver = mm_driver
    ENDIF
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_init_conf

  SUBROUTINE dbcsr_set_conf_mm_driver (driver, error)
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

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
    CASE (mm_driver_plasma)
       CALL dbcsr_assert (has_plasma,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for PLASMA not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_plasma
    CASE (mm_driver_cuda)
       CALL dbcsr_assert (has_cuda,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for CUDA not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_cuda
    END SELECT
    mm_async = driver_is_async (mm_driver)
    IF (mm_async) THEN
       mm_host_driver = mm_driver_smm
    ELSE
       mm_host_driver = mm_driver
    ENDIF
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_mm_driver

  FUNCTION dbcsr_get_conf_mm_driver () RESULT (driver)
    INTEGER                                  :: driver

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    driver = mm_driver
  END FUNCTION dbcsr_get_conf_mm_driver


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

  FUNCTION dbcsr_get_conf_subcomm () RESULT (use_subcomms)
    LOGICAL                                  :: use_subcomms

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_subcomm', &
      routineP = moduleN//':'//routineN

    use_subcomms = use_subcommunicators
  END FUNCTION dbcsr_get_conf_subcomm


  SUBROUTINE dbcsr_set_conf_combtypes (use_combtypes, error)
    LOGICAL, INTENT(IN)                      :: use_combtypes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_combtypes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_combined_types = use_combtypes
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_combtypes

  FUNCTION dbcsr_get_conf_combtypes () RESULT (use_combtypes)
    LOGICAL                                  :: use_combtypes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_combtypes', &
      routineP = moduleN//':'//routineN

    use_combtypes = use_combined_types
  END FUNCTION dbcsr_get_conf_combtypes


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

  FUNCTION dbcsr_get_conf_mpi_mem () RESULT (use_mpi_mem)
    LOGICAL                                  :: use_mpi_mem

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mpi_mem', &
      routineP = moduleN//':'//routineN

    use_mpi_mem = use_MPI_memory
  END FUNCTION dbcsr_get_conf_mpi_mem


  SUBROUTINE dbcsr_set_conf_cuda_mem (use_cuda_mem, error)
    LOGICAL, INTENT(IN)                      :: use_cuda_mem
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_cuda_mem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    CALL dbcsr_assert (use_cuda_mem, "IMP", has_cuda,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Can not use CUDA memory with no CUDA support!",&
         __LINE__, error=error)
    use_CUDA_host_pinned_memory = use_cuda_mem
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_cuda_mem

  FUNCTION dbcsr_get_conf_cuda_mem () RESULT (use_cuda_mem)
    LOGICAL                                  :: use_cuda_mem

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_cuda_mem', &
      routineP = moduleN//':'//routineN

    use_cuda_mem = use_CUDA_host_pinned_memory
  END FUNCTION dbcsr_get_conf_cuda_mem


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

  FUNCTION dbcsr_get_conf_mm_stacksize () RESULT (use_mm_stacksize)
    INTEGER                                  :: use_mm_stacksize

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mm_stacksize', &
      routineP = moduleN//':'//routineN

    use_mm_stacksize = mm_stack_size
  END FUNCTION dbcsr_get_conf_mm_stacksize


  SUBROUTINE dbcsr_set_conf_nstacks (nstacks, nbuffers, error)
    INTEGER, INTENT(IN)                      :: nstacks(:), nbuffers
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
    nstackbuffers    = nbuffers
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_nstacks

  SUBROUTINE dbcsr_get_conf_nstacks (nstacks, nbuffers, error)
    INTEGER, INTENT(OUT)                     :: nstacks(:), nbuffers
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
    nbuffers   = nstackbuffers
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_get_conf_nstacks

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


  FUNCTION dbcsr_get_conf_use_comm_thread () RESULT (comm_thread)
    LOGICAL                                  :: comm_thread

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_get_conf_use_comm_thread', &
      routineP = moduleN//':'//routineN

    comm_thread = use_comm_thread
  END FUNCTION dbcsr_get_conf_use_comm_thread

  SUBROUTINE dbcsr_set_conf_comm_thread_load (load, error)
    INTEGER, INTENT(IN)                      :: load
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_conf_comm_thread_load', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    comm_thread_load = load
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_comm_thread_load


  FUNCTION dbcsr_get_conf_comm_thread_load () RESULT (load)
    INTEGER                                  :: load

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_get_conf_comm_thread_load', &
      routineP = moduleN//':'//routineN

    load = comm_thread_load
  END FUNCTION dbcsr_get_conf_comm_thread_load


  PURE FUNCTION driver_is_async (driver)
    INTEGER, INTENT(IN)                      :: driver
    LOGICAL                                  :: driver_is_async

    driver_is_async = .FALSE.
    SELECT CASE (driver)
       CASE (mm_driver_cuda)
          driver_is_async = .TRUE.
    END SELECT
  END FUNCTION driver_is_async

! *****************************************************************************
!> \brief Set the callback routines for the external timing hooks.
!> \author  Ole Schuett
!> \date    2012-12-20
!>
!> <b>Modification history:</b>
!> - Created 2012-12-20
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

END MODULE dbcsr_config
