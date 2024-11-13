!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE openPMD_api
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_CHAR,&
                                              C_INT,&
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ENUM, BIND(C)
      ENUMERATOR :: openpmd_access_create, openpmd_access_read_only
   END ENUM

   TYPE openpmd_series_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE openpmd_series_type

   PUBLIC :: openpmd_series_type
   PUBLIC :: openpmd_series_create, openpmd_series_close
   PUBLIC :: openpmd_access_create, openpmd_access_read_only

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param path ...
!> \param access ...
!> \param series ...
! **************************************************************************************************
   SUBROUTINE openpmd_series_create(path, access, series)
      CHARACTER(len=*), INTENT(IN)                       :: path
      INTEGER(kind=C_INT), INTENT(IN)                    :: access
      TYPE(openpmd_series_type), INTENT(INOUT)           :: series

      INTERFACE
         SUBROUTINE openpmd_c_series_create(path, access, series) BIND(C, name="openPMD_Series_create")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR) :: series
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
            INTEGER(kind=C_INT), VALUE :: access
         END SUBROUTINE
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(series%c_ptr))
      CALL openpmd_c_series_create(path=TRIM(path)//C_NULL_CHAR, access=access, series=series%c_ptr)
      CPASSERT(C_ASSOCIATED(series%c_ptr))
   END SUBROUTINE openpmd_series_create

! **************************************************************************************************
!> \brief ...
!> \param series ...
! **************************************************************************************************
   SUBROUTINE openpmd_series_close(series)
      TYPE(openpmd_series_type), INTENT(INOUT)           :: series

      INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param series ...
! **************************************************************************************************
         SUBROUTINE openpmd_c_series_close(series) BIND(C, name="openPMD_Series_close")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: series
         END SUBROUTINE
      END INTERFACE

      CPASSERT(C_ASSOCIATED(series%c_ptr))
      CALL openpmd_c_series_close(series=series%c_ptr)
      series%c_ptr = C_NULL_PTR
   END SUBROUTINE

END MODULE openPMD_api
