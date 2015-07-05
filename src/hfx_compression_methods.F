!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines and types for Hartree-Fock-Exchange
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_compression_methods
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE hfx_compression_core_methods,    ONLY: bits2ints_specific,&
                                             ints2bits_specific
  USE hfx_types,                       ONLY: hfx_cache_type,&
                                             hfx_container_type
  USE kinds,                           ONLY: dp,&
                                             int_8
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  hfx_add_single_cache_element, hfx_get_single_cache_element,&
          hfx_reset_cache_and_container, hfx_decompress_first_cache, &
          hfx_flush_last_cache, hfx_add_mult_cache_elements,&
          hfx_get_mult_cache_elements


#define CACHE_SIZE 1024

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_compression_methods'

  INTEGER(kind=int_8), PARAMETER :: ugly_duck = ISHFT(1_int_8,63)
  INTEGER(int_8), PARAMETER :: shifts(0:63) = &
    (/ 1_int_8, 2_int_8, 4_int_8, 8_int_8, 16_int_8,32_int_8, 64_int_8, 128_int_8, 256_int_8, &
       512_int_8, 1024_int_8, 2048_int_8, 4096_int_8, 8192_int_8, 16384_int_8, 32768_int_8, &
       65536_int_8, 131072_int_8, 262144_int_8, 524288_int_8, 1048576_int_8, 2097152_int_8, &
       4194304_int_8, 8388608_int_8, 16777216_int_8, 33554432_int_8, 67108864_int_8, &
       134217728_int_8, 268435456_int_8, 536870912_int_8, 1073741824_int_8, 2147483648_int_8, &
       4294967296_int_8, 8589934592_int_8, 17179869184_int_8, 34359738368_int_8, 68719476736_int_8, &
       137438953472_int_8, 274877906944_int_8, 549755813888_int_8, 1099511627776_int_8, 2199023255552_int_8, &
       4398046511104_int_8, 8796093022208_int_8, 17592186044416_int_8, 35184372088832_int_8, 70368744177664_int_8, &
       140737488355328_int_8, 281474976710656_int_8, 562949953421312_int_8, 1125899906842624_int_8, &
       2251799813685248_int_8, 4503599627370496_int_8, 9007199254740992_int_8, 18014398509481984_int_8, &
       36028797018963968_int_8, 72057594037927936_int_8, 144115188075855872_int_8, 288230376151711744_int_8, &
       576460752303423488_int_8, 1152921504606846976_int_8, 2305843009213693952_int_8, &
       4611686018427387904_int_8, ugly_duck /)



!***


  CONTAINS

! *****************************************************************************
!> \brief - This routine adds an int_8 value to a cache. If the cache is full
!>        a compression routine is invoked and the cache is cleared
!> \param value value to be added to the cache
!> \param nbits number of bits to be stored
!> \param cache cache to which we want to add
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \param max_val_memory ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_add_single_cache_element(value, nbits, cache, container, memory_usage, use_disk_storage, &
                                          max_val_memory, error)
    INTEGER(int_8)                           :: value
    INTEGER                                  :: nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage
    INTEGER(int_8), OPTIONAL                 :: max_val_memory
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER(int_8)                           :: int_val

    int_val = value + shifts(nbits-1)

    IF(cache%element_counter /= CACHE_SIZE ) THEN
      cache%data(cache%element_counter) = int_val
      cache%element_counter = cache%element_counter + 1
    ELSE
      cache%data(CACHE_SIZE) = int_val
      CALL hfx_compress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage, &
                              max_val_memory, error)
      cache%element_counter = 1
    END IF
  END SUBROUTINE hfx_add_single_cache_element

! *****************************************************************************
!> \brief - This routine compresses a full cache and stores its values
!>        in a container. If necessary, a new list entry is allocated
!> \param full_array values from the cache
!> \param container linked list, that stores the compressed values
!> \param nbits number of bits to be stored
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \param max_val_memory ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_compress_cache(full_array, container, nbits, memory_usage, use_disk_storage, &
                                max_val_memory, error)
    INTEGER(int_8)                           :: full_array(*)
    TYPE(hfx_container_type)                 :: container
    INTEGER, INTENT(IN)                      :: nbits
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage
    INTEGER(int_8), OPTIONAL                 :: max_val_memory
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_compress_cache', &
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
      IF( use_disk_storage ) THEN
        !! write to file
        WRITE(container%unit) container%current%data
        !$OMP ATOMIC
        memory_usage = memory_usage + 1
        container%file_counter = container%file_counter + 1
      ELSE
        !! Allocate new list entry
        ALLOCATE(container%current%next, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        !$OMP ATOMIC
        memory_usage = memory_usage + 1
        container%current%next%next => NULL()
        container%current => container%current%next
        IF( PRESENT(max_val_memory) )  max_val_memory = max_val_memory + 1
      END IF
      !! compress remaining ints
      CALL ints2bits_specific(nbits, CACHE_SIZE-tmp_nints, container%current%data(1), full_array(tmp_nints+1))
      container%element_counter = 1 + (nbits*(CACHE_SIZE-tmp_nints) + 63)/64
    END IF

  END SUBROUTINE hfx_compress_cache

! *****************************************************************************
!> \brief - This routine returns an int_8 value from a cache. If the cache is empty
!>        a decompression routine is invoked and the cache is refilled with decompressed
!>        values from a container
!> \param value value to be retained from the cache
!> \param nbits number of bits with which the value has been compressed
!> \param cache cache from which we get the value
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_get_single_cache_element(value, nbits, cache, container, memory_usage, use_disk_storage)
    INTEGER(int_8)                           :: value
    INTEGER                                  :: nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage

    IF(cache%element_counter /= CACHE_SIZE ) THEN
      value = cache%data(cache%element_counter)
      cache%element_counter = cache%element_counter + 1
    ELSE
      value = cache%data(CACHE_SIZE)
      CALL hfx_decompress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage)
      cache%element_counter = 1
    END IF

    value = value - shifts(nbits-1)

  END SUBROUTINE hfx_get_single_cache_element

! *****************************************************************************
!> \brief - This routine decompresses data from a container in order to fill
!>        a cache.
!> \param full_array values to be retained from container
!> \param container linked list, that stores the compressed values
!> \param nbits number of bits with which the values have been stored
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_decompress_cache(full_array, container, nbits, memory_usage, use_disk_storage)
    INTEGER(int_8)                           :: full_array(*)
    TYPE(hfx_container_type)                 :: container
    INTEGER, INTENT(IN)                      :: nbits
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage

    INTEGER                                  :: end_idx, increment_counter, &
                                                start_idx, stat, &
                                                tmp_elements, tmp_nints

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
      IF( use_disk_storage ) THEN
          !! it could happen, that we are at the end of a file and we try to read
          !! This happens in case a container has fully been filled in the compression step
          !! but no other was needed for the current bit size
          !! Therefore we can safely igonore an eof error
          READ(container%unit, IOSTAT=stat) container%current%data
          memory_usage = memory_usage + 1
          container%file_counter = container%file_counter + 1
      ELSE
        container%current => container%current%next
        memory_usage = memory_usage + 1
      END IF
      !! decompress remaining ints
      CALL bits2ints_specific(nbits, CACHE_SIZE-tmp_nints, container%current%data(1), full_array(tmp_nints+1))
      container%element_counter = 1 + (nbits*(CACHE_SIZE-tmp_nints) + 63)/64
    END IF
  END SUBROUTINE hfx_decompress_cache

! *****************************************************************************
!> \brief - This routine resets the containers list pointer to the first element and
!>        moves the element counters of container and cache to the beginning
!> \param cache cache from which we get the value
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param do_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_reset_cache_and_container(cache, container, memory_usage, do_disk_storage)
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: do_disk_storage

    cache%element_counter = 1
    container%current => container%first
    container%element_counter = 1
    memory_usage = 1
    container%file_counter = 1
    IF( do_disk_storage ) THEN
      CALL close_file(container%unit)
      CALL open_file(file_name=container%filename,file_status="OLD",file_form="UNFORMATTED",file_action="READ",&
                     unit_number=container%unit)
      READ(container%unit) container%current%data
    END IF
  END SUBROUTINE hfx_reset_cache_and_container


! *****************************************************************************
!> \brief - This routine decompresses the first bunch of data in a container and
!>        copies them into a cache
!> \param nbits number of bits with which the data has been stored
!> \param cache array where we want to decompress the data
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_decompress_first_cache(nbits, cache, container, memory_usage, use_disk_storage)
    INTEGER                                  :: nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage

    CALL hfx_decompress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage)
    cache%element_counter = 1
  END SUBROUTINE hfx_decompress_first_cache

! *****************************************************************************
!> \brief - This routine compresses the last probably not yet compressed cache into
!>        a container
!> \param nbits number of bits with which the data has been stored
!> \param cache array where we want to decompress the data
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_flush_last_cache(nbits, cache, container, memory_usage, use_disk_storage, error)
    INTEGER                                  :: nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL hfx_compress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage, &
                            error=error)

    !!If we store to file, we have to make sure, that the last container is also written to disk
    IF( use_disk_storage ) THEN
      IF(container%element_counter /= 1) THEN
        WRITE(container%unit) container%current%data
        memory_usage = memory_usage + 1
        container%file_counter = container%file_counter + 1
      END IF
    END IF
  END SUBROUTINE hfx_flush_last_cache

! *****************************************************************************
!> \brief - This routine adds an a few real values to a cache. If the cache is full
!>        a compression routine is invoked and the cache is cleared
!> \param values values to be added to the cache
!> \param nints ...
!> \param nbits number of bits to be stored
!> \param cache cache to which we want to add
!> \param container container that contains the compressed elements
!> \param eps_schwarz ...
!> \param pmax_entry ...
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_add_mult_cache_elements(values, nints, nbits, cache, container, eps_schwarz, pmax_entry, memory_usage, &
                                          use_disk_storage, error)
    REAL(dp)                                 :: values(*)
    INTEGER, INTENT(IN)                      :: nints, nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    REAL(dp), INTENT(IN)                     :: eps_schwarz, pmax_entry
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: end_idx, i, start_idx, &
                                                tmp_elements
    INTEGER(int_8)                           :: shift, tmp
    REAL(dp)                                 :: eps_schwarz_inv, factor

    eps_schwarz_inv = 1.0_dp/eps_schwarz
    factor = eps_schwarz/pmax_entry

    shift = shifts(nbits-1)

    start_idx = cache%element_counter
    end_idx = start_idx + nints-1
    IF(end_idx< CACHE_SIZE) THEN
      DO i=1,nints
        values(i) = values(i)*pmax_entry
        IF(ABS(values(i)) > eps_schwarz) THEN
          tmp=NINT(values(i)*eps_schwarz_inv, KIND=int_8)
          cache%data(i+start_idx-1)= tmp + shift
          values(i) = tmp*factor
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
          tmp=NINT(values(i)*eps_schwarz_inv, KIND=int_8)
          cache%data(i+start_idx-1)= tmp + shift
          values(i) = tmp*factor
        ELSE
          values(i) = 0.0_dp
          cache%data(i+start_idx-1) = shift
        END IF
      END DO
      CALL hfx_compress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage, &
                              error=error)
      DO i=tmp_elements+1, nints
        values(i) = values(i)*pmax_entry
        IF(ABS(values(i)) > eps_schwarz) THEN
          tmp=NINT(values(i)*eps_schwarz_inv, KIND=int_8)
          cache%data(i-tmp_elements)= tmp + shift
          values(i) = tmp*factor
        ELSE
          values(i) = 0.0_dp
          cache%data(i-tmp_elements) = shift
        END IF
      END DO
      cache%element_counter = nints-tmp_elements + 1
    END IF
  END SUBROUTINE hfx_add_mult_cache_elements

! *****************************************************************************
!> \brief - This routine returns a bunch real values from a cache. If the cache is empty
!>        a decompression routine is invoked and the cache is refilled with decompressed
!>        values from a container
!> \param values value to be retained from the cache
!> \param nints number of values to be retained
!> \param nbits number of bits with which the value has been compressed
!> \param cache cache from which we get the value
!> \param container container that contains the compressed elements
!> \param eps_schwarz threshold for storage
!> \param pmax_entry multiplication factor for values
!> \param memory_usage ...
!> \param use_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_get_mult_cache_elements(values, nints, nbits, cache, container, eps_schwarz, pmax_entry, memory_usage, &
                                          use_disk_storage)
    REAL(dp)                                 :: values(*)
    INTEGER, INTENT(IN)                      :: nints, nbits
    TYPE(hfx_cache_type)                     :: cache
    TYPE(hfx_container_type)                 :: container
    REAL(dp), INTENT(IN)                     :: eps_schwarz, pmax_entry
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: use_disk_storage

    INTEGER                                  :: end_idx, i, start_idx, &
                                                tmp_elements
    INTEGER(int_8)                           :: shift
    REAL(dp)                                 :: factor

    factor = eps_schwarz/pmax_entry

    shift = shifts(nbits-1)

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
      CALL hfx_decompress_cache(cache%data(1), container, nbits, memory_usage, use_disk_storage)
      DO i=tmp_elements+1,nints
        values(i) = factor*REAL(cache%data(i-tmp_elements)-shift,dp)
      END DO
      cache%element_counter = nints - tmp_elements + 1
    END IF
  END SUBROUTINE hfx_get_mult_cache_elements

END MODULE hfx_compression_methods

