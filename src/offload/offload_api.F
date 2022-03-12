!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fortran API for the offload package, which is written in C.
!> \author Ole Schuett
! **************************************************************************************************
MODULE offload_api
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_CHAR,&
                                              C_F_POINTER,&
                                              C_INT,&
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR,&
                                              C_SIZE_T
   USE kinds,                           ONLY: dp,&
                                              int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'offload_api'

   PUBLIC :: offload_get_device_count
   PUBLIC :: offload_set_chosen_device, offload_get_chosen_device, offload_activate_chosen_device
   PUBLIC :: offload_timeset, offload_timestop, offload_mem_info
   PUBLIC :: offload_buffer_type, offload_create_buffer, offload_free_buffer
   PUBLIC :: offload_malloc_pinned_mem, offload_free_pinned_mem

   TYPE offload_buffer_type
      REAL(KIND=dp), DIMENSION(:), POINTER :: host_buffer => Null()
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE offload_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief allocate pinned memory.
!> \param buffer address of the buffer
!> \param length length of the buffer
!> \return 0
! **************************************************************************************************
   FUNCTION offload_malloc_pinned_mem(buffer, length) RESULT(res)
      TYPE(C_PTR)                                        :: buffer
      INTEGER(C_SIZE_T), VALUE                           :: length
      INTEGER                                            :: res

      INTERFACE
         FUNCTION offload_malloc_pinned_mem_c(buffer, length) &
            BIND(C, name="offload_host_malloc")
            IMPORT C_SIZE_T, C_PTR, C_INT
            TYPE(C_PTR)              :: buffer
            INTEGER(C_SIZE_T), VALUE :: length
            INTEGER(KIND=C_INT)                :: offload_malloc_pinned_mem_c
         END FUNCTION offload_malloc_pinned_mem_c
      END INTERFACE

      res = offload_malloc_pinned_mem_c(buffer, length)
   END FUNCTION offload_malloc_pinned_mem

! **************************************************************************************************
!> \brief free pinned memory
!> \param buffer address of the buffer
!> \return 0
! **************************************************************************************************
   FUNCTION offload_free_pinned_mem(buffer) RESULT(res)
      TYPE(C_PTR), VALUE                                 :: buffer
      INTEGER                                            :: res

      INTERFACE
         FUNCTION offload_free_pinned_mem_c(buffer) &
            BIND(C, name="offload_host_free")
            IMPORT C_PTR, C_INT
            INTEGER(KIND=C_INT)                :: offload_free_pinned_mem_c
            TYPE(C_PTR), VALUE            :: buffer
         END FUNCTION offload_free_pinned_mem_c
      END INTERFACE

      res = offload_free_pinned_mem_c(buffer)
   END FUNCTION offload_free_pinned_mem

! **************************************************************************************************
!> \brief Returns the number of available devices.
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION offload_get_device_count() RESULT(count)
      INTEGER                                            :: count

      INTERFACE
         FUNCTION offload_get_device_count_c() &
            BIND(C, name="offload_get_device_count")
            IMPORT :: C_INT
            INTEGER(KIND=C_INT)                :: offload_get_device_count_c
         END FUNCTION offload_get_device_count_c
      END INTERFACE

      count = offload_get_device_count_c()

   END FUNCTION offload_get_device_count

! **************************************************************************************************
!> \brief Selects the chosen device to be used.
!> \param device_id ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_set_chosen_device(device_id)
      INTEGER, INTENT(IN)                                :: device_id

      INTERFACE
         SUBROUTINE offload_set_chosen_device_c(device_id) &
            BIND(C, name="offload_set_chosen_device")
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE                :: device_id
         END SUBROUTINE offload_set_chosen_device_c
      END INTERFACE

      CALL offload_set_chosen_device_c(device_id=device_id)

   END SUBROUTINE offload_set_chosen_device

! **************************************************************************************************
!> \brief Returns the chosen device.
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION offload_get_chosen_device() RESULT(device_id)
      INTEGER                                            :: device_id

      INTERFACE
         FUNCTION offload_get_chosen_device_c() &
            BIND(C, name="offload_get_chosen_device")
            IMPORT :: C_INT
            INTEGER(KIND=C_INT)                :: offload_get_chosen_device_c
         END FUNCTION offload_get_chosen_device_c
      END INTERFACE

      device_id = offload_get_chosen_device_c()

      IF (device_id < 0) &
         CPABORT("No offload device has been chosen.")

   END FUNCTION offload_get_chosen_device

! **************************************************************************************************
!> \brief Activates the device selected via offload_set_chosen_device()
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_activate_chosen_device()

      INTERFACE
         SUBROUTINE offload_activate_chosen_device_c() &
            BIND(C, name="offload_activate_chosen_device")
         END SUBROUTINE offload_activate_chosen_device_c
      END INTERFACE

      CALL offload_activate_chosen_device_c()

   END SUBROUTINE offload_activate_chosen_device

! **************************************************************************************************
!> \brief Starts a timing range.
!> \param routineN ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_timeset(routineN)
      CHARACTER(LEN=*), INTENT(IN)                       :: routineN

      INTERFACE
         SUBROUTINE offload_timeset_c(message) BIND(C, name="offload_timeset")
            IMPORT :: C_CHAR
            CHARACTER(kind=C_CHAR), DIMENSION(*), INTENT(IN) :: message
         END SUBROUTINE offload_timeset_c
      END INTERFACE

      CALL offload_timeset_c(TRIM(routineN)//C_NULL_CHAR)

   END SUBROUTINE offload_timeset

! **************************************************************************************************
!> \brief  Ends a timing range.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_timestop()

      INTERFACE
         SUBROUTINE offload_timestop_c() BIND(C, name="offload_timestop")
         END SUBROUTINE offload_timestop_c
      END INTERFACE

      CALL offload_timestop_c()

   END SUBROUTINE offload_timestop

! **************************************************************************************************
!> \brief Gets free and total device memory.
!> \param free ...
!> \param total ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_mem_info(free, total)
      INTEGER(KIND=int_8), INTENT(OUT)                   :: free, total

      INTEGER(KIND=C_SIZE_T)                             :: my_free, my_total
      INTERFACE
         SUBROUTINE offload_mem_info_c(free, total) BIND(C, name="offload_mem_info")
            IMPORT :: C_SIZE_T
            INTEGER(KIND=C_SIZE_T)                   :: free, total
         END SUBROUTINE offload_mem_info_c
      END INTERFACE

      CALL offload_mem_info_c(my_free, my_total)

      ! On 32-bit architectures this converts from int_4 to int_8.
      free = my_free
      total = my_total

   END SUBROUTINE offload_mem_info

! **************************************************************************************************
!> \brief Allocates a buffer of given length, ie. number of elements.
!> \param length ...
!> \param buffer ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_create_buffer(length, buffer)
      INTEGER, INTENT(IN)                                :: length
      TYPE(offload_buffer_type), INTENT(INOUT)           :: buffer

      CHARACTER(LEN=*), PARAMETER :: routineN = 'offload_create_buffer'

      INTEGER                                            :: handle
      TYPE(C_PTR)                                        :: host_buffer_c
      INTERFACE
         SUBROUTINE offload_create_buffer_c(length, buffer) &
            BIND(C, name="offload_create_buffer")
            IMPORT :: C_PTR, C_INT
            INTEGER(KIND=C_INT), VALUE                :: length
            TYPE(C_PTR)                               :: buffer
         END SUBROUTINE offload_create_buffer_c
      END INTERFACE
      INTERFACE

         FUNCTION offload_get_buffer_host_pointer_c(buffer) &
            BIND(C, name="offload_get_buffer_host_pointer")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: buffer
            TYPE(C_PTR)                               :: offload_get_buffer_host_pointer_c
         END FUNCTION offload_get_buffer_host_pointer_c
      END INTERFACE

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(buffer%host_buffer)) THEN
         IF (SIZE(buffer%host_buffer) == 0) DEALLOCATE (buffer%host_buffer)
      END IF

      CALL offload_create_buffer_c(length=length, buffer=buffer%c_ptr)
      CPASSERT(C_ASSOCIATED(buffer%c_ptr))

      IF (length == 0) THEN
         ! While C_F_POINTER usually accepts a NULL pointer it's not standard compliant.
         ALLOCATE (buffer%host_buffer(0))
      ELSE
         host_buffer_c = offload_get_buffer_host_pointer_c(buffer%c_ptr)
         CPASSERT(C_ASSOCIATED(host_buffer_c))
         CALL C_F_POINTER(host_buffer_c, buffer%host_buffer, shape=(/length/))
      END IF

      CALL timestop(handle)
   END SUBROUTINE offload_create_buffer

! **************************************************************************************************
!> \brief Deallocates given buffer.
!> \param buffer ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_free_buffer(buffer)
      TYPE(offload_buffer_type), INTENT(INOUT)           :: buffer

      CHARACTER(LEN=*), PARAMETER :: routineN = 'offload_free_buffer'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE offload_free_buffer_c(buffer) &
            BIND(C, name="offload_free_buffer")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: buffer
         END SUBROUTINE offload_free_buffer_c
      END INTERFACE

      CALL timeset(routineN, handle)

      CPASSERT(C_ASSOCIATED(buffer%c_ptr))

      CALL offload_free_buffer_c(buffer%c_ptr)

      buffer%c_ptr = C_NULL_PTR

      IF (SIZE(buffer%host_buffer) == 0) THEN
         DEALLOCATE (buffer%host_buffer)
      ELSE
         NULLIFY (buffer%host_buffer)
      END IF

      CALL timestop(handle)
   END SUBROUTINE offload_free_buffer
END MODULE offload_api
