!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type to store integrals for semi-empirical calculations
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   05.2008
! *****************************************************************************
MODULE semi_empirical_store_int_types

  
  USE hfx_compression_methods,         ONLY: hfx_decompress_first_cache,&
                                             hfx_flush_last_cache,&
                                             hfx_reset_cache_and_container
  USE hfx_types,                       ONLY: hfx_cache_type,&
                                             hfx_container_type,&
                                             hfx_init_container,&
                                             hfx_memory_type,&
                                             parse_memory_section
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_store_int_types'

! *****************************************************************************
!> \brief Semi-empirical store integrals type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   05.2008
! *****************************************************************************
  TYPE semi_empirical_si_type
     LOGICAL                                     :: filling_containers, compress
     INTEGER                                     :: nbuffer
     REAL(KIND=dp), POINTER, DIMENSION(:)        :: max_val_buffer, uncompressed_container
     TYPE(hfx_memory_type)                       :: memory_parameter
     TYPE(hfx_cache_type), DIMENSION(:), &
          POINTER                                :: integral_caches
     TYPE(hfx_container_type), DIMENSION(:), &
          POINTER                                :: integral_containers
  END TYPE semi_empirical_si_type


  PUBLIC :: semi_empirical_si_type,&
            semi_empirical_si_create,&
            semi_empirical_si_release,&
            semi_empirical_si_finalize,&
            semi_empirical_si_initialize

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical store integrals type
!> \param store_int_env ...
!> \param se_section ...
!> \param compression ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE semi_empirical_si_create(store_int_env, se_section, compression)
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(section_vals_type), POINTER         :: se_section
    LOGICAL, INTENT(in), OPTIONAL            :: compression

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_si_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: se_mem_section

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(store_int_env))
    ALLOCATE (store_int_env)
    store_int_env%filling_containers = .TRUE.
    store_int_env%nbuffer             = 0
    NULLIFY(store_int_env%max_val_buffer, store_int_env%uncompressed_container)

    ! Memory section
    se_mem_section => section_vals_get_subs_vals(se_section,"MEMORY")
    IF ( PRESENT(compression) ) THEN
      store_int_env%compress = compression
    ELSE
      CALL section_vals_val_get(se_mem_section,"COMPRESS",l_val=store_int_env%compress)
    END IF
    CALL parse_memory_section(store_int_env%memory_parameter, se_mem_section, skip_disk=.TRUE.,&
         skip_in_core_forces =.TRUE.)
    store_int_env%memory_parameter%ram_counter = 0
    ! If we don't compress there's no cache
    IF (.NOT.store_int_env%compress) THEN
       store_int_env%memory_parameter%cache_size = 1
    END IF

    ! Disk Storage disabled for semi-empirical methods
    CALL cp_assert(.NOT.store_int_env%memory_parameter%do_disk_storage,cp_failure_level,&
         cp_assertion_failed,routineP,"Disk storage for SEMIEMPIRICAL methods disabled! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    ! Allocate containers/caches for integral storage if requested
    IF(.NOT.store_int_env%memory_parameter%do_all_on_the_fly.AND.store_int_env%compress) THEN
       ALLOCATE(store_int_env%integral_containers(64))
       ALLOCATE(store_int_env%integral_caches(64))
       DO i=1,64
          store_int_env%integral_caches(i)%element_counter = 1
          store_int_env%integral_caches(i)%data = 0
          ALLOCATE(store_int_env%integral_containers(i)%first)
          store_int_env%integral_containers(i)%first%prev => NULL()
          store_int_env%integral_containers(i)%first%next => NULL()
          store_int_env%integral_containers(i)%current => store_int_env%integral_containers(i)%first
          store_int_env%integral_containers(i)%current%data = 0
          store_int_env%integral_containers(i)%element_counter = 1
       END DO
    END IF
  END SUBROUTINE semi_empirical_si_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical store integrals type
!> \param store_int_env ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE semi_empirical_si_release(store_int_env)
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_si_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(store_int_env)) THEN
       ! Deallocate containers/caches
       IF(.NOT.store_int_env%memory_parameter%do_all_on_the_fly) THEN
          IF (store_int_env%compress) THEN
             ! Deallocate containers/caches
             DO i=1,64
                CALL hfx_init_container(store_int_env%integral_containers(i),&
                     store_int_env%memory_parameter%actual_memory_usage,&
                     .FALSE.)
                DEALLOCATE(store_int_env%integral_containers(i)%first)
             END DO
             IF (ASSOCIATED(store_int_env%max_val_buffer)) THEN
                DEALLOCATE(store_int_env%max_val_buffer)
             END IF
             DEALLOCATE(store_int_env%integral_containers)
             DEALLOCATE(store_int_env%integral_caches)
          ELSE
             IF (ASSOCIATED(store_int_env%uncompressed_container)) THEN
                DEALLOCATE(store_int_env%uncompressed_container)
             END IF
          END IF
       END IF
       ! Deallocate the full store_int_env
       DEALLOCATE (store_int_env)
    END IF

  END SUBROUTINE semi_empirical_si_release

! *****************************************************************************
!> \brief Deallocate the semi-empirical store integrals type
!> \param store_int_env ...
!> \param geometry_did_change ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE semi_empirical_si_initialize(store_int_env, geometry_did_change)
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    LOGICAL, INTENT(IN)                      :: geometry_did_change

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_si_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(store_int_env)) THEN
       IF(.NOT.store_int_env%memory_parameter%do_all_on_the_fly) THEN
          IF (geometry_did_change) THEN
             store_int_env%filling_containers = .TRUE.
             store_int_env%nbuffer            = 0
             store_int_env%memory_parameter%ram_counter = HUGE(store_int_env%memory_parameter%ram_counter)
             IF (store_int_env%compress) THEN
                ! Compress integrals
                CALL reallocate(store_int_env%max_val_buffer, 1, store_int_env%nbuffer)
                ! Clean containers
                DO i=1,64
                   CALL hfx_init_container(store_int_env%integral_containers(i),&
                        store_int_env%memory_parameter%actual_memory_usage,&
                        .FALSE.)
                END DO
             ELSE
                ! Skip compression
                CALL reallocate(store_int_env%uncompressed_container, 1, 0)
                store_int_env%memory_parameter%actual_memory_usage = 1
             END IF
          ELSE
             store_int_env%filling_containers = .FALSE.
             store_int_env%nbuffer            = 0
             IF (store_int_env%compress) THEN
                ! Retrieve data into the cache
                DO i=1,64
                   CALL hfx_decompress_first_cache(i,store_int_env%integral_caches(i),&
                        store_int_env%integral_containers(i),&
                        store_int_env%memory_parameter%actual_memory_usage,.FALSE.)
                END DO
             ELSE
                store_int_env%memory_parameter%actual_memory_usage = 1
             END IF
          END IF
       END IF
    END IF

  END SUBROUTINE semi_empirical_si_initialize

! *****************************************************************************
!> \brief Deallocate the semi-empirical store integrals type
!> \param store_int_env ...
!> \param geometry_did_change ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE semi_empirical_si_finalize(store_int_env, geometry_did_change)
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    LOGICAL, INTENT(IN)                      :: geometry_did_change

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_si_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(store_int_env)) THEN
       IF(.NOT.store_int_env%memory_parameter%do_all_on_the_fly) THEN
          IF (geometry_did_change) THEN
             IF (store_int_env%compress) THEN
                ! Flush last cache
                DO i=1,64
                   CALL hfx_flush_last_cache(i, store_int_env%integral_caches(i),&
                        store_int_env%integral_containers(i),&
                        store_int_env%memory_parameter%actual_memory_usage, .FALSE.)
                END DO
                ! Reallocate this array with the proper size
                CALL reallocate(store_int_env%max_val_buffer, 1, store_int_env%nbuffer)
             ELSE
                ! Skip compression
                CALL reallocate(store_int_env%uncompressed_container, 1,&
                     store_int_env%memory_parameter%actual_memory_usage-1)
             END IF
          END IF
          IF (store_int_env%compress) THEN
             ! Reset caches and containers
             DO i=1,64
                CALL hfx_reset_cache_and_container( store_int_env%integral_caches(i),&
                     store_int_env%integral_containers(i), store_int_env%memory_parameter%actual_memory_usage,&
                     .FALSE.)
             END DO
          END IF
       END IF
    END IF

  END SUBROUTINE semi_empirical_si_finalize

END MODULE semi_empirical_store_int_types
