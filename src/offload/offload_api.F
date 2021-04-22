!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fortran API for the offload package, which is written in C.
!> \author Ole Schuett
! **************************************************************************************************
MODULE offload_api
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_F_POINTER,&
                                              C_INT,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'offload_api'

   PUBLIC :: offload_set_device_id
   PUBLIC :: offload_buffer_type, offload_create_buffer, offload_free_buffer

   TYPE offload_buffer_type
      REAL(KIND=dp), DIMENSION(:), POINTER :: host_buffer => Null()
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE offload_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Selects the device to be used.
!> \param device_id ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE offload_set_device_id(device_id)
      INTEGER, INTENT(IN)                                :: device_id

      INTERFACE
         SUBROUTINE offload_set_device_id_c(device_id) &
            BIND(C, name="offload_set_device_id")
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE                :: device_id
         END SUBROUTINE offload_set_device_id_c
      END INTERFACE

      CALL offload_set_device_id_c(device_id=device_id)

   END SUBROUTINE offload_set_device_id

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
      ENDIF

      CALL offload_create_buffer_c(length=length, buffer=buffer%c_ptr)
      CPASSERT(C_ASSOCIATED(buffer%c_ptr))

      IF (length == 0) THEN
         ! While C_F_POINTER usually accepts a NULL pointer it's not standard compliant.
         ALLOCATE (buffer%host_buffer(0))
      ELSE
         host_buffer_c = offload_get_buffer_host_pointer_c(buffer%c_ptr)
         CPASSERT(C_ASSOCIATED(host_buffer_c))
         CALL C_F_POINTER(host_buffer_c, buffer%host_buffer, shape=(/length/))
      ENDIF

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
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE offload_free_buffer

END MODULE offload_api
