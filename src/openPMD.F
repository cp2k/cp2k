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
   USE message_passing,                 ONLY: mp_comm_type
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
!> \param mpi_comm ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_series_create(path, access, mpi_comm) RESULT(series)
      CHARACTER(len=*), INTENT(IN)                       :: path
      INTEGER(kind=C_INT), INTENT(IN)                    :: access
      TYPE(mp_comm_type), INTENT(in), OPTIONAL           :: mpi_comm
      TYPE(openpmd_series_type)                          :: series

      TYPE(c_ptr)                                        :: comm_c = C_NULL_PTR
      INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param comm_f ...
!> \param comm_c ...
!> \param C ...
!> \param name="CP2K_MPI_Comm_f2c" ...
! **************************************************************************************************
         SUBROUTINE cp2k_c_mpi_comm_f2c(comm_f, comm_c) bind(C, name="CP2K_MPI_Comm_f2c")
            import :: c_ptr, mp_comm_type, c_int
            INTEGER(kind=C_INT), VALUE :: comm_f
            TYPE(c_ptr) :: comm_c
         END SUBROUTINE
      END INTERFACE
      INTERFACE
         SUBROUTINE openpmd_c_series_create_mpi(path, access, mpi_comm, series) BIND(C, name="openPMD_Series_create_mpi")
            IMPORT :: C_PTR, C_CHAR, C_INT
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
            INTEGER(kind=C_INT), VALUE :: access
            TYPE(C_PTR), value :: mpi_comm
            TYPE(C_PTR) :: series
         END SUBROUTINE
      END INTERFACE
      INTERFACE
         SUBROUTINE openpmd_c_series_create(path, access, series) BIND(C, name="openPMD_Series_create")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR) :: series
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
            INTEGER(kind=C_INT), VALUE :: access
         END SUBROUTINE
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(series%c_ptr))
      IF (PRESENT(mpi_comm)) THEN
         CALL cp2k_c_mpi_comm_f2c(mpi_comm%get_handle(), comm_c)
         CALL openpmd_c_series_create_mpi(path=TRIM(path)//C_NULL_CHAR, access=access, mpi_comm=comm_c, series=series%c_ptr)
      ELSE
         CALL openpmd_c_series_create(path=TRIM(path)//C_NULL_CHAR, access=access, series=series%c_ptr)
      END IF
      CPASSERT(C_ASSOCIATED(series%c_ptr))
   END FUNCTION openpmd_series_create

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
