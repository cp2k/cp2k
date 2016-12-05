!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Configuration options for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-05
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-05
! **************************************************************************************************
MODULE dbcsr_config
   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_config'

   REAL(KIND=dp), PARAMETER             :: default_resize_factor = 1.2_dp

   ! Possible drivers to use for matrix multiplications
   INTEGER, PARAMETER :: mm_driver_auto = 0
   INTEGER, PARAMETER :: mm_driver_matmul = 1
   INTEGER, PARAMETER :: mm_driver_blas = 2
   INTEGER, PARAMETER :: mm_driver_smm = 3
   INTEGER, PARAMETER :: mm_driver_xsmm = 4

   CHARACTER(len=*), PARAMETER :: mm_name_auto = "AUTO", &
                                  mm_name_blas = "BLAS", &
                                  mm_name_matmul = "MATMUL", &
                                  mm_name_smm = "SMM", &
                                  mm_name_xsmm = "XSMM"

#if defined (__HAS_smm_dnn)
   LOGICAL, PARAMETER :: has_smm = .TRUE.
#else
   LOGICAL, PARAMETER :: has_smm = .FALSE.
#endif

#if defined(__HAS_smm_vec)
   LOGICAL, PARAMETER :: has_smm_vec = .TRUE.
#else
   LOGICAL, PARAMETER :: has_smm_vec = .FALSE.
#endif

#if defined(__LIBXSMM)
   LOGICAL, PARAMETER :: has_xsmm = .TRUE.
#else
   LOGICAL, PARAMETER :: has_xsmm = .FALSE.
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

#if defined(__HAS_smm_vec) || defined (__DBCSR_ACC)
   INTEGER, PARAMETER :: mm_stack_default_size = 30000
#else
   INTEGER, PARAMETER :: mm_stack_default_size = 1000
#endif

#if defined(__LIBXSMM)
   INTEGER, PARAMETER :: mm_default_driver = mm_driver_xsmm
#elif defined (__HAS_smm_dnn)
   INTEGER, PARAMETER :: mm_default_driver = mm_driver_smm
#else
   INTEGER, PARAMETER :: mm_default_driver = mm_driver_blas ! always available
#endif

   TYPE dbcsr_config_type
      INTEGER  :: mm_driver = mm_default_driver
      LOGICAL  :: use_MPI_memory = has_MPI
      INTEGER  :: mm_stack_size = mm_stack_default_size
      INTEGER  :: max_elements_per_block = 32
      INTEGER  :: avg_elements_images = -1
      INTEGER  :: num_mult_images = 1
      INTEGER  :: nm_stacks = 3
      INTEGER  :: nn_stacks = 3
      INTEGER  :: nk_stacks = 3
      LOGICAL  :: use_mpi_exp = .FALSE.
      LOGICAL  :: use_mpi_filtering = .TRUE.
      INTEGER  :: num_layers_3D = 1
      LOGICAL  :: use_comm_thread = .TRUE.
      INTEGER  :: comm_thread_load = 100
      INTEGER  :: multrec_limit = 512
      INTEGER  :: accdrv_priority_streams = 4
      INTEGER  :: accdrv_priority_buffers = 40
      INTEGER  :: accdrv_posterior_streams = 4
      INTEGER  :: accdrv_posterior_buffers = 80
      LOGICAL  :: accdrv_avoid_after_busy = .FALSE.
      INTEGER  :: accdrv_min_flop_process = 0
      INTEGER  :: accdrv_min_flop_sort = 4000
      LOGICAL  :: accdrv_do_inhomogenous = .TRUE.
      INTEGER  :: accdrv_binning_nbins = 4096
      INTEGER  :: accdrv_binning_binsize = 16
   END TYPE dbcsr_config_type

   TYPE(dbcsr_config_type), PROTECTED, SAVE :: dbcsr_cfg = dbcsr_config_type() ! defaults

   PUBLIC :: dbcsr_cfg, has_MPI, has_acc, default_resize_factor
   PUBLIC :: mm_driver_blas, mm_driver_matmul, mm_driver_smm, mm_driver_xsmm, mm_driver_auto
   PUBLIC :: dbcsr_set_config, dbcsr_get_default_config, dbcsr_print_config

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mm_driver ...
!> \param use_MPI_memory ...
!> \param mm_stack_size ...
!> \param max_elements_per_block Default blocking parameter
!> \param avg_elements_images Maximum number of elements for each image
!> \param num_mult_images Multiplicative factor for number of virtual images
!> \param nstacks Number of stacks to use
!> \param use_mpi_exp RMA algorithm
!> \param use_mpi_filtering RMA MPI filtering
!> \param num_layers_3D 3D layers
!> \param use_comm_thread ...
!> \param comm_thread_load ...
!> \param multrec_limit ...
!> \param accdrv_priority_streams ...
!> \param accdrv_priority_buffers ...
!> \param accdrv_posterior_streams ...
!> \param accdrv_posterior_buffers ...
!> \param accdrv_avoid_after_busy ...
!> \param accdrv_min_flop_process ...
!> \param accdrv_min_flop_sort ...
!> \param accdrv_do_inhomogenous ...
!> \param accdrv_binning_nbins ...
!> \param accdrv_binning_binsize ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set_config( &
      mm_driver, &
      use_MPI_memory, &
      mm_stack_size, &
      max_elements_per_block, &
      avg_elements_images, &
      num_mult_images, &
      nstacks, &
      use_mpi_exp, &
      use_mpi_filtering, &
      num_layers_3D, &
      use_comm_thread, &
      comm_thread_load, &
      multrec_limit, &
      accdrv_priority_streams, &
      accdrv_priority_buffers, &
      accdrv_posterior_streams, &
      accdrv_posterior_buffers, &
      accdrv_avoid_after_busy, &
      accdrv_min_flop_process, &
      accdrv_min_flop_sort, &
      accdrv_do_inhomogenous, &
      accdrv_binning_nbins, &
      accdrv_binning_binsize)
!
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: mm_driver
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_MPI_memory
      INTEGER, INTENT(IN), OPTIONAL                      :: mm_stack_size, max_elements_per_block, &
                                                            avg_elements_images, num_mult_images, &
                                                            nstacks
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_mpi_exp, use_mpi_filtering
      INTEGER, INTENT(IN), OPTIONAL                      :: num_layers_3D
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_comm_thread
      INTEGER, INTENT(IN), OPTIONAL :: comm_thread_load, multrec_limit, accdrv_priority_streams, &
         accdrv_priority_buffers, accdrv_posterior_streams, accdrv_posterior_buffers
      LOGICAL, INTENT(IN), OPTIONAL                      :: accdrv_avoid_after_busy
      INTEGER, INTENT(IN), OPTIONAL                      :: accdrv_min_flop_process, &
                                                            accdrv_min_flop_sort
      LOGICAL, INTENT(IN), OPTIONAL                      :: accdrv_do_inhomogenous
      INTEGER, INTENT(IN), OPTIONAL                      :: accdrv_binning_nbins, &
                                                            accdrv_binning_binsize

      INTEGER                                            :: nthreads

      IF (PRESENT(use_MPI_memory)) dbcsr_cfg%use_MPI_memory = use_MPI_memory
      IF (PRESENT(max_elements_per_block)) dbcsr_cfg%max_elements_per_block = max_elements_per_block
      IF (PRESENT(avg_elements_images)) dbcsr_cfg%avg_elements_images = avg_elements_images
      IF (PRESENT(num_mult_images)) dbcsr_cfg%num_mult_images = num_mult_images
      IF (PRESENT(use_mpi_exp)) dbcsr_cfg%use_mpi_exp = use_mpi_exp
      IF (PRESENT(use_mpi_filtering)) dbcsr_cfg%use_mpi_filtering = use_mpi_filtering
      IF (PRESENT(num_layers_3D)) dbcsr_cfg%num_layers_3D = num_layers_3D
      IF (PRESENT(use_comm_thread)) dbcsr_cfg%use_comm_thread = use_comm_thread
      IF (PRESENT(multrec_limit)) dbcsr_cfg%multrec_limit = multrec_limit
      IF (PRESENT(accdrv_priority_streams)) dbcsr_cfg%accdrv_priority_streams = accdrv_priority_streams
      IF (PRESENT(accdrv_priority_buffers)) dbcsr_cfg%accdrv_priority_buffers = accdrv_priority_buffers
      IF (PRESENT(accdrv_posterior_streams)) dbcsr_cfg%accdrv_posterior_streams = accdrv_posterior_streams
      IF (PRESENT(accdrv_posterior_buffers)) dbcsr_cfg%accdrv_posterior_buffers = accdrv_posterior_buffers
      IF (PRESENT(accdrv_avoid_after_busy)) dbcsr_cfg%accdrv_avoid_after_busy = accdrv_avoid_after_busy
      IF (PRESENT(accdrv_min_flop_process)) dbcsr_cfg%accdrv_min_flop_process = accdrv_min_flop_process
      IF (PRESENT(accdrv_min_flop_sort)) dbcsr_cfg%accdrv_min_flop_sort = accdrv_min_flop_sort
      IF (PRESENT(accdrv_do_inhomogenous)) dbcsr_cfg%accdrv_do_inhomogenous = accdrv_do_inhomogenous
      IF (PRESENT(accdrv_binning_nbins)) dbcsr_cfg%accdrv_binning_nbins = accdrv_binning_nbins
      IF (PRESENT(accdrv_binning_binsize)) dbcsr_cfg%accdrv_binning_binsize = accdrv_binning_binsize

      IF (PRESENT(comm_thread_load)) THEN
         dbcsr_cfg%comm_thread_load = comm_thread_load
         IF (comm_thread_load < 0) THEN
            nthreads = 1
!$          nthreads = OMP_GET_MAX_THREADS()
            dbcsr_cfg%comm_thread_load = MAX(0, 90-(30*nthreads)/8)
         ENDIF
      ENDIF

      IF (PRESENT(nstacks)) THEN
         dbcsr_cfg%nm_stacks = nstacks
         dbcsr_cfg%nn_stacks = nstacks
         dbcsr_cfg%nk_stacks = nstacks
      ENDIF

      IF (PRESENT(mm_stack_size)) THEN
         dbcsr_cfg%mm_stack_size = mm_stack_size
         IF (mm_stack_size < 0) &
            dbcsr_cfg%mm_stack_size = mm_stack_default_size
      ENDIF

      IF (PRESENT(mm_driver)) THEN
         IF (TRIM(mm_driver) .EQ. mm_name_auto) THEN
            dbcsr_cfg%mm_driver = mm_default_driver
         ELSE IF (TRIM(mm_driver) .EQ. mm_name_blas) THEN
            dbcsr_cfg%mm_driver = mm_driver_blas ! always available
         ELSE IF (TRIM(mm_driver) .EQ. mm_name_matmul) THEN
            dbcsr_cfg%mm_driver = mm_driver_matmul ! always available
         ELSE IF (TRIM(mm_driver) .EQ. mm_name_smm) THEN
            IF (.NOT. has_smm) CPABORT("Support for libsmm not compiled in.")
            dbcsr_cfg%mm_driver = mm_driver_smm
         ELSE IF (TRIM(mm_driver) .EQ. mm_name_xsmm) THEN
            IF (.NOT. has_xsmm) CPABORT("Support for libxsmm not compiled in.")
            dbcsr_cfg%mm_driver = mm_driver_xsmm
         ELSE
            CPABORT("Unkown MM driver: "//TRIM(mm_driver))
         ENDIF
      ENDIF

   END SUBROUTINE dbcsr_set_config

! **************************************************************************************************
!> \brief ...
!> \param use_MPI_memory ...
!> \param mm_stack_size ...
!> \param max_elements_per_block ...
!> \param avg_elements_images ...
!> \param num_mult_images ...
!> \param nstacks ...
!> \param use_mpi_exp ...
!> \param use_mpi_filtering ...
!> \param num_layers_3D ...
!> \param use_comm_thread ...
!> \param comm_thread_load ...
!> \param multrec_limit ...
!> \param accdrv_priority_streams ...
!> \param accdrv_priority_buffers ...
!> \param accdrv_posterior_streams ...
!> \param accdrv_posterior_buffers ...
!> \param accdrv_avoid_after_busy ...
!> \param accdrv_min_flop_process ...
!> \param accdrv_min_flop_sort ...
!> \param accdrv_do_inhomogenous ...
!> \param accdrv_binning_nbins ...
!> \param accdrv_binning_binsize ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_default_config( &
      use_MPI_memory, &
      mm_stack_size, &
      max_elements_per_block, &
      avg_elements_images, &
      num_mult_images, &
      nstacks, &
      use_mpi_exp, &
      use_mpi_filtering, &
      num_layers_3D, &
      use_comm_thread, &
      comm_thread_load, &
      multrec_limit, &
      accdrv_priority_streams, &
      accdrv_priority_buffers, &
      accdrv_posterior_streams, &
      accdrv_posterior_buffers, &
      accdrv_avoid_after_busy, &
      accdrv_min_flop_process, &
      accdrv_min_flop_sort, &
      accdrv_do_inhomogenous, &
      accdrv_binning_nbins, &
      accdrv_binning_binsize)
!
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_MPI_memory
      INTEGER, INTENT(OUT), OPTIONAL                     :: mm_stack_size, max_elements_per_block, &
                                                            avg_elements_images, num_mult_images, &
                                                            nstacks
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_mpi_exp, use_mpi_filtering
      INTEGER, INTENT(OUT), OPTIONAL                     :: num_layers_3D
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_comm_thread
      INTEGER, INTENT(OUT), OPTIONAL :: comm_thread_load, multrec_limit, accdrv_priority_streams, &
         accdrv_priority_buffers, accdrv_posterior_streams, accdrv_posterior_buffers
      LOGICAL, INTENT(OUT), OPTIONAL                     :: accdrv_avoid_after_busy
      INTEGER, INTENT(OUT), OPTIONAL                     :: accdrv_min_flop_process, &
                                                            accdrv_min_flop_sort
      LOGICAL, INTENT(OUT), OPTIONAL                     :: accdrv_do_inhomogenous
      INTEGER, INTENT(OUT), OPTIONAL                     :: accdrv_binning_nbins, &
                                                            accdrv_binning_binsize

      TYPE(dbcsr_config_type)                            :: default_cfg

      IF (PRESENT(use_MPI_memory)) use_MPI_memory = default_cfg%use_MPI_memory
      IF (PRESENT(mm_stack_size)) mm_stack_size = default_cfg%mm_stack_size
      IF (PRESENT(max_elements_per_block)) max_elements_per_block = default_cfg%max_elements_per_block
      IF (PRESENT(avg_elements_images)) avg_elements_images = default_cfg%avg_elements_images
      IF (PRESENT(num_mult_images)) num_mult_images = default_cfg%num_mult_images
      IF (PRESENT(use_mpi_exp)) use_mpi_exp = default_cfg%use_mpi_exp
      IF (PRESENT(use_mpi_filtering)) use_mpi_filtering = default_cfg%use_mpi_filtering
      IF (PRESENT(num_layers_3D)) num_layers_3D = default_cfg%num_layers_3D
      IF (PRESENT(use_comm_thread)) use_comm_thread = default_cfg%use_comm_thread
      IF (PRESENT(comm_thread_load)) comm_thread_load = default_cfg%comm_thread_load
      IF (PRESENT(multrec_limit)) multrec_limit = default_cfg%multrec_limit
      IF (PRESENT(accdrv_priority_streams)) accdrv_priority_streams = default_cfg%accdrv_priority_streams
      IF (PRESENT(accdrv_priority_buffers)) accdrv_priority_buffers = default_cfg%accdrv_priority_buffers
      IF (PRESENT(accdrv_posterior_streams)) accdrv_posterior_streams = default_cfg%accdrv_posterior_streams
      IF (PRESENT(accdrv_posterior_buffers)) accdrv_posterior_buffers = default_cfg%accdrv_posterior_buffers
      IF (PRESENT(accdrv_avoid_after_busy)) accdrv_avoid_after_busy = default_cfg%accdrv_avoid_after_busy
      IF (PRESENT(accdrv_min_flop_process)) accdrv_min_flop_process = default_cfg%accdrv_min_flop_process
      IF (PRESENT(accdrv_min_flop_sort)) accdrv_min_flop_sort = default_cfg%accdrv_min_flop_sort
      IF (PRESENT(accdrv_do_inhomogenous)) accdrv_do_inhomogenous = default_cfg%accdrv_do_inhomogenous
      IF (PRESENT(accdrv_binning_nbins)) accdrv_binning_nbins = default_cfg%accdrv_binning_nbins
      IF (PRESENT(accdrv_binning_binsize)) accdrv_binning_binsize = default_cfg%accdrv_binning_binsize

      CPASSERT(default_cfg%nm_stacks == dbcsr_cfg%nn_stacks)
      CPASSERT(default_cfg%nm_stacks == dbcsr_cfg%nk_stacks)
      IF (PRESENT(nstacks)) nstacks = default_cfg%nm_stacks

   END SUBROUTINE dbcsr_get_default_config

! **************************************************************************************************
!> \brief Prints configuration for DBCSR
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print_config(unit_nr)
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_string_length)               :: mm_name

      IF (unit_nr > 0) THEN
         SELECT CASE (dbcsr_cfg%mm_driver)
         CASE (mm_driver_blas); mm_name = mm_name_blas
         CASE (mm_driver_matmul); mm_name = mm_name_matmul
         CASE (mm_driver_smm); mm_name = mm_name_smm
         CASE (mm_driver_xsmm); mm_name = mm_name_xsmm
         CASE DEFAULT
            CPABORT("Unkown MM driver")
         END SELECT
         WRITE (UNIT=unit_nr, FMT='(1X,A,T41,A40)') &
            "DBCSR| Multiplication driver", ADJUSTR(mm_name(1:40))

         WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
            "DBCSR| Multrec recursion limit", dbcsr_cfg%multrec_limit
         WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
            "DBCSR| Multiplication stack size", dbcsr_cfg%mm_stack_size

         IF (dbcsr_cfg%avg_elements_images > 0) THEN
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Average elements for images", dbcsr_cfg%avg_elements_images
         ELSE
            WRITE (UNIT=unit_nr, FMT='(1X,A,T72,A)') &
               "DBCSR| Maximum elements for images", "UNLIMITED"
         ENDIF
         WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
            "DBCSR| Multiplicative factor virtual images", dbcsr_cfg%num_mult_images

         IF (dbcsr_cfg%nm_stacks == dbcsr_cfg%nn_stacks .AND. dbcsr_cfg%nm_stacks == dbcsr_cfg%nk_stacks) THEN
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Multiplication size stacks", dbcsr_cfg%nm_stacks
         ELSE
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Multiplication size m stacks", dbcsr_cfg%nm_stacks
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Multiplication size n stacks", dbcsr_cfg%nn_stacks
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Multiplication size k stacks", dbcsr_cfg%nk_stacks
         ENDIF

         IF (has_mpi) THEN
            IF (dbcsr_cfg%num_layers_3D < 2) THEN
               WRITE (UNIT=unit_nr, FMT='(1X,A,T75,A)') &
                  "DBCSR| Number of 3D layers", "SINGLE"
            ELSE
               WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
                  "DBCSR| Number of 3D layers", dbcsr_cfg%num_layers_3D
            ENDIF
            WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
               "DBCSR| Use MPI memory allocation", dbcsr_cfg%use_mpi_memory
            WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
               "DBCSR| Use RMA algorithm", dbcsr_cfg%use_mpi_exp
            WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
               "DBCSR| Use Communication thread", dbcsr_cfg%use_comm_thread
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| Communication thread load", dbcsr_cfg%comm_thread_load
            IF (dbcsr_cfg%use_mpi_exp) THEN
               WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
                  "DBCSR| Use MPI filtering", dbcsr_cfg%use_mpi_filtering
            ENDIF
         ENDIF

         IF (has_acc) THEN
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Number of priority stack-buffers", dbcsr_cfg%accdrv_priority_buffers
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Number of posterior stack-buffers", dbcsr_cfg%accdrv_posterior_buffers
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Number of priority streams", dbcsr_cfg%accdrv_priority_streams
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Number of posterior streams", dbcsr_cfg%accdrv_posterior_streams
            WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
               "DBCSR| ACC: Avoid driver after busy ", dbcsr_cfg%accdrv_avoid_after_busy
            WRITE (UNIT=unit_nr, FMT='(1X,A,T80,L1)') &
               "DBCSR| ACC: Process inhomogenous stacks", dbcsr_cfg%accdrv_do_inhomogenous
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Min. flop for processing", dbcsr_cfg%accdrv_min_flop_process
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Min. flop for sorting", dbcsr_cfg%accdrv_min_flop_sort
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Number of binning bins", dbcsr_cfg%accdrv_binning_nbins
            WRITE (UNIT=unit_nr, FMT='(1X,A,T70,I11)') &
               "DBCSR| ACC: Size of binning bins", dbcsr_cfg%accdrv_binning_binsize
         END IF
      ENDIF

   END SUBROUTINE dbcsr_print_config

END MODULE dbcsr_config
