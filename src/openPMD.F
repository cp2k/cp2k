!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE openPMD_api
   USE ISO_C_BINDING,                   ONLY: &
        C_ASSOCIATED, C_CHAR, C_INT, C_INT64_T, C_LOC, C_NULL_CHAR, C_NULL_PTR, C_PTR, C_SIZE_T
   USE kinds,                           ONLY: default_string_length
   USE message_passing,                 ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ENUM, BIND(C)
      ENUMERATOR :: openpmd_access_create, openpmd_access_read_only
   END ENUM

   TYPE openpmd_attributable_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS

   END TYPE openpmd_attributable_type

   TYPE openpmd_series_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: as_pointer
      PROCEDURE, PUBLIC :: openpmd_series_as_attributable
   END TYPE openpmd_series_type

   TYPE openpmd_iteration_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: openpmd_iteration_as_attributable
   END TYPE openpmd_iteration_type

   TYPE openpmd_mesh_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC:: openpmd_mesh_as_record_component
   END TYPE openpmd_mesh_type

   TYPE openpmd_particle_species_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE openpmd_particle_species_type

   TYPE openpmd_record_component_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE openpmd_record_component_type

   TYPE openpmd_record_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: openpmd_record_as_record_component
   END TYPE openpmd_record_type

   PUBLIC :: openpmd_series_type, openpmd_iteration_type, openpmd_mesh_type, openpmd_particle_species_type, openpmd_record_component_type, openpmd_record_type
   PUBLIC :: openpmd_series_create, openpmd_series_close
   PUBLIC :: openpmd_series_present
   PUBLIC :: openpmd_series_write_iteration
   PUBLIC :: openpmd_access_create, openpmd_access_read_only
   PUBLIC :: openpmd_get_default_extension
   PUBLIC :: openpmd_attributable_set_attribute_vec_int
   PUBLIC :: openpmd_iteration_get_mesh, openpmd_iteration_get_particle_species
   PUBLIC :: openpmd_record_component_make_empty_int
   PUBLIC :: openpmd_particle_species_get_record
   PUBLIC :: cp2k_log

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION as_pointer(this) RESULT(res)
      CLASS(openpmd_series_type), INTENT(IN) :: this
      TYPE(C_PTR) :: res

      res = this%c_ptr
   END FUNCTION as_pointer

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_series_as_attributable(this) RESULT(res)
      CLASS(openpmd_series_type), INTENT(IN) :: this
      TYPE(openpmd_attributable_type) :: res

      INTERFACE
         SUBROUTINE openpmd_c_series_as_attributable(series, attributable) BIND(C, name="openPMD_Iteration_upcast_to_attributable")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: series
            TYPE(C_PTR) :: attributable
         END SUBROUTINE
      END INTERFACE

      CALL openpmd_c_series_as_attributable(this%c_ptr, res%c_ptr)
   END FUNCTION openpmd_series_as_attributable

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_mesh_as_record_component(this) RESULT(res)
      CLASS(openpmd_mesh_type), INTENT(IN) :: this
      TYPE(openpmd_record_component_type) :: res

      INTERFACE
        SUBROUTINE openpmd_c_mesh_as_record_component(mesh, record_component) BIND(C, name="openPMD_Mesh_upcast_to_RecordComponent")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: mesh
            TYPE(C_PTR) :: record_component
         END SUBROUTINE
      END INTERFACE

      CALL openpmd_c_mesh_as_record_component(this%c_ptr, res%c_ptr)
   END FUNCTION openpmd_mesh_as_record_component

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_record_as_record_component(this) RESULT(res)
      CLASS(openpmd_record_type), INTENT(IN) :: this
      TYPE(openpmd_record_component_type) :: res

      INTERFACE
    SUBROUTINE openpmd_c_record_as_record_component(record, record_component) BIND(C, name="openPMD_Mesh_upcast_to_RecordComponent")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: record
            TYPE(C_PTR) :: record_component
         END SUBROUTINE
      END INTERFACE

      CALL openpmd_c_record_as_record_component(this%c_ptr, res%c_ptr)
   END FUNCTION openpmd_record_as_record_component

! **************************************************************************************************
!> \brief ...
!> \param series ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_series_present(series) RESULT(res)
      CLASS(openpmd_series_type), INTENT(IN) :: series
      LOGICAL :: res
      INTEGER :: res_internal

      INTERFACE
         FUNCTION openpmd_c_series_present(series) RESULT(res) BIND(C, NAME="openPMD_Series_present")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE, INTENT(IN) :: series
            INTEGER(kind=C_INT) :: res
         END FUNCTION openpmd_c_series_present
      END INTERFACE

      IF (.NOT. C_ASSOCIATED(series%c_ptr)) THEN
         res = .FALSE.
         RETURN
      END IF

      res_internal = openpmd_c_series_present(series%c_ptr)
      res = res_internal /= 0
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_iteration_as_attributable(this) RESULT(res)
      CLASS(openpmd_iteration_type), INTENT(IN) :: this
      TYPE(openpmd_attributable_type) :: res

      INTERFACE
    SUBROUTINE openpmd_c_iteration_as_attributable(iteration, attributable) BIND(C, name="openPMD_Iteration_upcast_to_attributable")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE :: iteration
            TYPE(C_PTR) :: attributable
         END SUBROUTINE
      END INTERFACE

      CALL openpmd_c_iteration_as_attributable(this%c_ptr, res%c_ptr)
   END FUNCTION openpmd_iteration_as_attributable

! **************************************************************************************************
!> \brief ...
!> \param message ...
! **************************************************************************************************
   SUBROUTINE cp2k_log(message)
      CHARACTER(len=*), INTENT(IN)                       :: message

      INTERFACE
         SUBROUTINE cp2k_c_log(message) BIND(C, name="CP2K_log")
            IMPORT :: C_PTR, C_CHAR
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: message
         END SUBROUTINE
      END INTERFACE

      CALL cp2k_c_log(message=TRIM(message)//C_NULL_CHAR)
   END SUBROUTINE cp2k_log

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

      CLASS(mp_comm_type), INTENT(in), OPTIONAL          :: mpi_comm
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

! **************************************************************************************************
!> \brief ...
!> \param series ...
!> \param index ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_series_write_iteration(series, index) RESULT(iteration)
      TYPE(openpmd_series_type)                          :: series
      INTEGER, TARGET                                    :: index
      TYPE(openpmd_iteration_type)                       :: iteration

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
   FUNCTION openpmd_c_series_write_iteration(series, index, iteration) RESULT(status) BIND(C, name="openPMD_Series_write_Iteration")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T
            TYPE(C_PTR), VALUE :: series
            TYPE(C_PTR), VALUE :: index
            TYPE(C_PTR) :: iteration
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_series_write_iteration
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(iteration%c_ptr))
      status = openpmd_c_series_write_iteration(series%c_ptr, c_loc(index), iteration%c_ptr)
      CPASSERT(C_ASSOCIATED(iteration%c_ptr))
   END FUNCTION openpmd_series_write_iteration

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_get_default_extension() RESULT(extension)
      CHARACTER(len=default_string_length)               :: extension

      CHARACTER(kind=C_CHAR), POINTER                    :: c_extension
      INTEGER(C_SIZE_T)                                  :: length_of_c_string
      INTERFACE
         FUNCTION openpmd_c_get_default_extension() RESULT(extension) BIND(c, name="openPMD_get_default_extension")
            IMPORT :: C_CHAR
            CHARACTER(kind=C_CHAR), POINTER :: extension
         END FUNCTION
      END INTERFACE
      INTERFACE
         FUNCTION get_strlen(c_string) RESULT(strlen) BIND(C, NAME="CP2K_strlen")
            IMPORT :: C_CHAR, C_SIZE_T
            CHARACTER(kind=C_CHAR), DIMENSION(*), INTENT(IN) :: c_string
            INTEGER(C_SIZE_T) :: strlen
         END FUNCTION get_strlen
      END INTERFACE

      c_extension => openpmd_c_get_default_extension()

      ! Find the length of the C string
      IF (ASSOCIATED(c_extension)) THEN
         !c_string_array => c_extension
         length_of_c_string = get_strlen(c_extension)
      ELSE
         RETURN
      END IF

      ! Ensure we do not exceed the length of extension
      IF (length_of_c_string > default_string_length) THEN
         length_of_c_string = default_string_length
      END IF

      ! Copy the contents from c_extension to extension
      extension = TRIM(ADJUSTL(TRANSFER(c_extension(1:length_of_c_string + 1), extension)))

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param attributable ...
!> \param name ...
!> \param vec ...
! **************************************************************************************************
   SUBROUTINE openpmd_attributable_set_attribute_vec_int(attributable, name, vec)
      TYPE(openpmd_attributable_type)                    :: attributable
      CHARACTER(len=*), INTENT(IN)                       :: name
      INTEGER, DIMENSION(:), TARGET                      :: vec

      INTEGER, TARGET                                    :: size_
      INTERFACE
         SUBROUTINE openpmd_c_attributable_set_attribute_vec_int(attributable, name, vec, length) BIND(C, NAME="openPMD_attributable_set_attribute_vec_int")
            IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
            TYPE(C_PTR), VALUE :: attributable
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
            TYPE(C_PTR), VALUE :: vec
            INTEGER(kind=C_INT), VALUE :: length
         END SUBROUTINE
      END INTERFACE

      size_ = SIZE(vec)
      CALL openpmd_c_attributable_set_attribute_vec_int(attributable%c_ptr, name, C_LOC(vec(1)), size_)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param iteration ...
!> \param name ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_iteration_get_mesh(iteration, name) RESULT(mesh)
      TYPE(openpmd_iteration_type)                       :: iteration
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(openpmd_mesh_type)                            :: mesh

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_iteration_get_mesh(iteration, name, mesh) RESULT(status) BIND(C, name="openPMD_Iteration_get_mesh")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: iteration
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
            TYPE(C_PTR) :: mesh
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_iteration_get_mesh
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(mesh%c_ptr))
      status = openpmd_c_iteration_get_mesh(iteration%c_ptr, name//C_NULL_CHAR, mesh%c_ptr)
      CPASSERT(C_ASSOCIATED(mesh%c_ptr))
   END FUNCTION openpmd_iteration_get_mesh

! **************************************************************************************************
!> \brief ...
!> \param iteration ...
!> \param name ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_iteration_get_particle_species(iteration, name) RESULT(particle_species)
      TYPE(openpmd_iteration_type)                       :: iteration
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(openpmd_particle_species_type)                :: particle_species

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_iteration_get_particle_species(iteration, name, particle_species) RESULT(status) BIND(C, name="openPMD_Iteration_get_particle_species")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: iteration
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
            TYPE(C_PTR) :: particle_species
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_iteration_get_particle_species
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(particle_species%c_ptr))
      status = openpmd_c_iteration_get_particle_species(iteration%c_ptr, name//C_NULL_CHAR, particle_species%c_ptr)
      CPASSERT(C_ASSOCIATED(particle_species%c_ptr))
   END FUNCTION openpmd_iteration_get_particle_species

! **************************************************************************************************
!> \brief ...
!> \param record_component ...
!> \param dimensionality ...
! **************************************************************************************************
   SUBROUTINE openpmd_record_component_make_empty_int(record_component, dimensionality)
      TYPE(openpmd_record_component_type)                :: record_component
      INTEGER                                            :: dimensionality

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_record_component_make_empty_int(record_component, dimensionality) RESULT(status) BIND(C, name="openPMD_RecordComponent_makeEmpty_int")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: record_component
            INTEGER(KIND=C_INT), VALUE :: dimensionality
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_record_component_make_empty_int
      END INTERFACE

      status = openpmd_c_record_component_make_empty_int(record_component%c_ptr, dimensionality)
   END SUBROUTINE openpmd_record_component_make_empty_int

! **************************************************************************************************
!> \brief ...
!> \param particle_species ...
!> \param name ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_particle_species_get_record(particle_species, name) RESULT(record)
      TYPE(openpmd_particle_species_type)                :: particle_species
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(openpmd_record_type)                          :: record

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_particle_species_get_record(particle_species, name, record) RESULT(status) BIND(C, name="openPMD_ParticleSpecies_get_Record")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: particle_species
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
            TYPE(C_PTR) :: record
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_particle_species_get_record
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(record%c_ptr))
      status = openpmd_c_particle_species_get_record(particle_species%c_ptr, name//C_NULL_CHAR, record%c_ptr)
      CPASSERT(C_ASSOCIATED(record%c_ptr))
   END FUNCTION openpmd_particle_species_get_record

END MODULE openPMD_api
