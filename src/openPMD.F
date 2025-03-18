!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE openPMD_api
   USE ISO_C_BINDING, ONLY: &
      C_ASSOCIATED, C_CHAR, C_INT, C_INT64_T, C_LOC, C_NULL_CHAR, C_NULL_PTR, C_PTR, C_SIZE_T, C_F_POINTER
   USE kinds, ONLY: default_string_length
   USE message_passing, ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ENUM, BIND(C)
      ENUMERATOR :: openpmd_access_create, openpmd_access_read_only
   END ENUM

   ENUM, BIND(C)
      ENUMERATOR :: openpmd_type_char, openpmd_type_uchar, openpmd_type_schar, openpmd_type_short, openpmd_type_int, openpmd_type_long, openpmd_type_longlong, openpmd_type_ushort, openpmd_type_uint, openpmd_type_ulong, openpmd_type_ulonglong, openpmd_type_float, openpmd_type_double, openpmd_type_long_double, openpmd_type_cfloat, openpmd_type_cdouble, openpmd_type_clong_double, openpmd_type_string, openpmd_type_vec_char, openpmd_type_vec_short, openpmd_type_vec_int, openpmd_type_vec_long, openpmd_type_vec_longlong, openpmd_type_vec_uchar, openpmd_type_vec_ushort, openpmd_type_vec_uint, openpmd_type_vec_ulong, openpmd_type_vec_ulonglong, openpmd_type_vec_float, openpmd_type_vec_double, openpmd_type_vec_long_double, openpmd_type_vec_cfloat, openpmd_type_vec_cdouble, openpmd_type_vec_clong_double, openpmd_type_vec_schar, openpmd_type_vec_string, openpmd_type_arr_dbl_7, openpmd_type_bool
   END ENUM

   #:set dataset_types_f = ['REAL(kind=dp)']
   #:set dataset_types_enum = ['double']
   #:set dimensions = [1, 2, 3]

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

   #:for dim in dimensions
      TYPE openpmd_dynamic_memory_view_type_${dim}$d
         PRIVATE
         INTEGER, DIMENSION(${dim}$) :: chunk_extent
         TYPE(C_PTR) :: c_ptr = C_NULL_PTR
      END TYPE openpmd_dynamic_memory_view_type_${dim}$d
   #:endfor

   PUBLIC :: openpmd_series_type, openpmd_iteration_type, openpmd_mesh_type, openpmd_particle_species_type, openpmd_record_component_type, openpmd_record_type
   PUBLIC :: openpmd_series_create, openpmd_series_close
   PUBLIC :: openpmd_series_present
   PUBLIC :: openpmd_series_write_iteration
   PUBLIC :: openpmd_access_create, openpmd_access_read_only
   PUBLIC :: openpmd_type_char, openpmd_type_uchar, openpmd_type_schar, openpmd_type_short, openpmd_type_int, openpmd_type_long, openpmd_type_longlong, openpmd_type_ushort, openpmd_type_uint, openpmd_type_ulong, openpmd_type_ulonglong, openpmd_type_float, openpmd_type_double, openpmd_type_long_double, openpmd_type_cfloat, openpmd_type_cdouble, openpmd_type_clong_double, openpmd_type_string, openpmd_type_vec_char, openpmd_type_vec_short, openpmd_type_vec_int, openpmd_type_vec_long, openpmd_type_vec_longlong, openpmd_type_vec_uchar, openpmd_type_vec_ushort, openpmd_type_vec_uint, openpmd_type_vec_ulong, openpmd_type_vec_ulonglong, openpmd_type_vec_float, openpmd_type_vec_double, openpmd_type_vec_long_double, openpmd_type_vec_cfloat, openpmd_type_vec_cdouble, openpmd_type_vec_clong_double, openpmd_type_vec_schar, openpmd_type_vec_string, openpmd_type_arr_dbl_7, openpmd_type_bool
   PUBLIC :: openpmd_get_default_extension
   PUBLIC :: openpmd_attributable_set_attribute_vec_int, openpmd_attributable_series_flush
   PUBLIC :: openpmd_iteration_get_mesh, openpmd_iteration_get_particle_species, openpmd_iteration_as_attributable
   PUBLIC :: openpmd_record_component_make_empty, openpmd_record_component_make_constant_zero, openpmd_record_get_component, openpmd_record_component_as_attributable, openpmd_record_component_reset_dataset
   #:for type_f, type_enum in zip(dataset_types_f, dataset_types_enum)
      #:for dim in dimensions
         PUBLIC :: openpmd_record_component_store_chunk_${dim}$d_${type_enum}$
      #:endfor
   #:endfor
   PUBLIC :: openpmd_particle_species_get_record
   PUBLIC :: cp2k_log

   INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 30)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)

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

   #:for child, parent in [("Series", "Attributable"), ("Iteration", "Attributable"), ("RecordComponent", "Attributable")]
      #:set fortran_names = {"RecordComponent": "record_component"}
      #:set child_f = fortran_names[child] if child in fortran_names else child
      #:set parent_f = fortran_names[parent] if parent in fortran_names else parent
      FUNCTION openpmd_${child_f}$_as_${parent_f}$ (this) RESULT(res)
         CLASS(openpmd_${child_f}$_type), INTENT(IN) :: this
         TYPE(openpmd_${parent_f}$_type) :: res

         INTERFACE
            SUBROUTINE openpmd_c_${child_f}$_as_${parent_f}$ (child, parent) BIND(C, name="openPMD_${child}$_upcast_to_${parent}$")
               IMPORT :: C_PTR
               TYPE(C_PTR), VALUE :: child
               TYPE(C_PTR) :: parent
            END SUBROUTINE
         END INTERFACE

         CALL openpmd_c_${child_f}$_as_${parent_f}$ (this%c_ptr, res%c_ptr)
      END FUNCTION openpmd_${child_f}$_as_${parent_f}$
   #:endfor

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
!> \param config ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_series_create(path, access, mpi_comm, config) RESULT(series)
      CHARACTER(len=*), INTENT(IN)                       :: path
      INTEGER(kind=C_INT), INTENT(IN)                    :: access
      CLASS(mp_comm_type), INTENT(in), OPTIONAL          :: mpi_comm
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: config
      TYPE(openpmd_series_type)                          :: series

      TYPE(c_ptr)                                        :: comm_c = C_NULL_PTR
      CHARACTER(len=default_string_length) :: my_config = "{}"
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
         SUBROUTINE openpmd_c_series_create_mpi(path, access, mpi_comm, config, series) BIND(C, name="openPMD_Series_create_mpi")
            IMPORT :: C_PTR, C_CHAR, C_INT
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
            INTEGER(kind=C_INT), VALUE :: access
            TYPE(C_PTR), value :: mpi_comm
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: config
            TYPE(C_PTR) :: series
         END SUBROUTINE
      END INTERFACE
      INTERFACE
         SUBROUTINE openpmd_c_series_create(path, access, config, series) BIND(C, name="openPMD_Series_create")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR) :: series
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
            INTEGER(kind=C_INT), VALUE :: access
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: config
         END SUBROUTINE
      END INTERFACE

      IF (PRESENT(config)) my_config = config

      CPASSERT(.NOT. C_ASSOCIATED(series%c_ptr))
      IF (PRESENT(mpi_comm)) THEN
         CALL cp2k_c_mpi_comm_f2c(mpi_comm%get_handle(), comm_c)
         CALL openpmd_c_series_create_mpi(path=TRIM(path)//C_NULL_CHAR, access=access, mpi_comm=comm_c, config=my_config//C_NULL_CHAR, series=series%c_ptr)
      ELSE
       CALL openpmd_c_series_create(path=TRIM(path)//C_NULL_CHAR, access=access, config=my_config//C_NULL_CHAR, series=series%c_ptr)
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

      INTEGER(8)                                         :: cast_index
      INTEGER(kind=C_INT)                                :: status
      INTERFACE
   FUNCTION openpmd_c_series_write_iteration(series, index, iteration) RESULT(status) BIND(C, name="openPMD_Series_write_Iteration")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T
            TYPE(C_PTR), VALUE :: series
            INTEGER(kind=C_INT64_T), VALUE :: index
            TYPE(C_PTR) :: iteration
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_series_write_iteration
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(iteration%c_ptr))
      cast_index = index
      status = openpmd_c_series_write_iteration(series%c_ptr, cast_index, iteration%c_ptr)
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

   SUBROUTINE openpmd_attributable_series_flush(attributable)
      TYPE(openpmd_attributable_type)                    :: attributable
      INTERFACE
         SUBROUTINE openpmd_c_attributable_series_flush(attributable) BIND(C, NAME="openPMD_Attributable_series_flush")
            IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
            TYPE(C_PTR), VALUE :: attributable
         END SUBROUTINE
      END INTERFACE
      CALL openpmd_c_attributable_series_flush(attributable%c_ptr)
   END SUBROUTINE openpmd_attributable_series_flush

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
!> \param dtype ...
!> \param dimensionality ...
! **************************************************************************************************
   SUBROUTINE openpmd_record_component_make_empty(record_component, dtype, dimensionality)
      TYPE(openpmd_record_component_type)                :: record_component
      INTEGER(kind=C_INT)                                :: dtype
      INTEGER                                            :: dimensionality

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_record_component_make_empty(record_component, dtype, dimensionality) RESULT(status) BIND(C, name="openPMD_RecordComponent_makeEmpty")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: record_component
            INTEGER(KIND=C_INT), VALUE :: dtype
            INTEGER(KIND=C_INT), VALUE :: dimensionality
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_record_component_make_empty
      END INTERFACE

      status = openpmd_c_record_component_make_empty(record_component%c_ptr, dtype, dimensionality)
   END SUBROUTINE openpmd_record_component_make_empty

! **************************************************************************************************
!> \brief ...
!> \param record_component ...
!> \param dtype ...
!> \param extent ...
! **************************************************************************************************
   SUBROUTINE openpmd_record_component_make_constant_zero(record_component, dtype, extent)
      TYPE(openpmd_record_component_type)                :: record_component
      INTEGER(kind=C_INT)                                :: dtype
      INTEGER, DIMENSION(:), TARGET                      :: extent

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_record_component_make_constant(record_component, dtype, dimensions, extent, value) RESULT(status) BIND(C, name="openPMD_RecordComponent_makeConstant")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: record_component
            INTEGER(KIND=C_INT), VALUE :: dtype
            INTEGER(KIND=C_INT), VALUE :: dimensions
            TYPE(C_PTR), VALUE :: extent
            TYPE(C_PTR), VALUE :: value
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_record_component_make_constant
      END INTERFACE

      ! TODO: flip the extent?

      status = openpmd_c_record_component_make_constant(record_component%c_ptr, dtype, SIZE(extent), C_LOC(extent(1)), C_NULL_PTR)
   END SUBROUTINE openpmd_record_component_make_constant_zero

   SUBROUTINE openpmd_record_component_reset_dataset(record_component, dtype, extent, cfg)
      TYPE(openpmd_record_component_type)                :: record_component
      INTEGER(kind=C_INT)                                :: dtype
      INTEGER, DIMENSION(:), TARGET                      :: extent
      CHARACTER(len=*), OPTIONAL :: cfg

      CHARACTER(len=default_string_length) :: my_config = "{}"
      INTEGER(kind=C_INT)                                :: status
      INTERFACE
         FUNCTION openpmd_c_record_component_reset_dataset(record_component, dtype, dimensions, extent, cfg) RESULT(status) BIND(C, name="openPMD_RecordComponent_resetDataset")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: record_component
            INTEGER(KIND=C_INT), VALUE :: dtype
            INTEGER(KIND=C_INT), VALUE :: dimensions
            TYPE(C_PTR), VALUE :: extent
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: cfg
            INTEGER(KIND=C_INT) :: status
         END FUNCTION
      END INTERFACE

      IF (PRESENT(cfg)) my_config = cfg

      ! TODO: flip the extent?
      status = openpmd_c_record_component_reset_dataset(record_component%c_ptr, dtype, SIZE(extent), C_LOC(extent(1)), my_config//C_NULL_CHAR)
   END SUBROUTINE openpmd_record_component_reset_dataset

   #:for type_f, type_enum in zip(dataset_types_f, dataset_types_enum)
      #:for dim in [1, 2, 3]
         ! e.g. ':,:,:'
         #:set dimensionality_wildcard = ','.join(':' for _ in range(dim))
         SUBROUTINE openpmd_record_component_store_chunk_${dim}$d_${type_enum}$ (record_component, data_, offset)
            TYPE(openpmd_record_component_type)                        :: record_component
            ${type_f}$, DIMENSION(${dimensionality_wildcard}$), TARGET :: data_
            INTEGER, DIMENSION(${dim}$), TARGET                        :: offset, extent

            INTEGER(kind=C_INT)                                        :: status

            INTERFACE
            FUNCTION openpmd_c_record_component_store_chunk(record_component, dtype, dimensions, offset, extent, data_) RESULT(status) BIND(C, name="openPMD_RecordComponent_storeChunk")
                  IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
                  TYPE(C_PTR), VALUE :: record_component
                  INTEGER(KIND=C_INT), VALUE :: dtype
                  INTEGER(KIND=C_INT), VALUE :: dimensions
                  TYPE(C_PTR), VALUE :: offset, extent
                  TYPE(C_PTR), VALUE :: data_
                  INTEGER(KIND=C_INT) :: status
               END FUNCTION openpmd_c_record_component_store_chunk
            END INTERFACE

            ! TODO flip extent, offset
            extent = SIZE(data_)
            status = openpmd_c_record_component_store_chunk(record_component%c_ptr, openpmd_type_${type_enum}$, ${dim}$, C_LOC(offset(1)), C_LOC(extent(1)), C_LOC(data_))

         END SUBROUTINE openpmd_record_component_store_chunk_${dim}$d_${type_enum}$

    FUNCTION openpmd_record_component_store_chunk_span_${dim}$d_${type_enum}$ (record_component, offset, extent) RESULT(memory_view)
            TYPE(openpmd_record_component_type)                        :: record_component
            INTEGER, DIMENSION(${dim}$), TARGET                        :: offset, extent

            TYPE(openpmd_dynamic_memory_view_type_${dim}$d)            :: memory_view

            INTEGER(kind=C_INT)                                        :: status

            INTERFACE
            FUNCTION openpmd_c_record_component_store_chunk_span(record_component, dtype, dimensions, offset, extent, memory_view) RESULT(status) BIND(C, name="openPMD_RecordComponent_storeChunkSpan")
                  IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
                  TYPE(C_PTR), VALUE :: record_component
                  INTEGER(KIND=C_INT), VALUE :: dtype
                  INTEGER(KIND=C_INT), VALUE :: dimensions
                  TYPE(C_PTR), VALUE :: offset, extent
                  TYPE(C_PTR) :: memory_view
                  INTEGER(KIND=C_INT) :: status
               END FUNCTION openpmd_c_record_component_store_chunk_span
            END INTERFACE

            ! TODO flip extent, offset
            memory_view%chunk_extent = extent
            status = openpmd_c_record_component_store_chunk_span(record_component%c_ptr, openpmd_type_${type_enum}$, ${dim}$, C_LOC(offset(1)), C_LOC(extent(1)), memory_view%c_ptr)

         END FUNCTION openpmd_record_component_store_chunk_span_${dim}$d_${type_enum}$

         FUNCTION openpmd_dynamic_memory_view_resolve_${dim}$d(memory_view) RESULT(buffer)
            TYPE(openpmd_dynamic_memory_view_type_${dim}$d)            :: memory_view
            ${type_f}$, POINTER :: buffer(${dimensionality_wildcard}$)

            INTEGER(kind=C_INT)                                        :: status
            TYPE(C_PTR), TARGET                                        :: c_buffer

            INTERFACE
               FUNCTION openpmd_c_dynamic_memory_view_resolve(memory_view, write_buffer) RESULT(status) BIND(C, name="openPMD_DynamicMemoryView_resolve")
                  IMPORT :: C_PTR, C_INT
                  TYPE(C_PTR), VALUE :: memory_view
                  TYPE(C_PTR) :: write_buffer
                  INTEGER(KIND=C_INT) :: status
               END FUNCTION
            END INTERFACE

            status = openpmd_c_dynamic_memory_view_resolve(memory_view%c_ptr, c_buffer)
            CALL c_f_pointer(c_buffer, buffer, memory_view%chunk_extent)
         END FUNCTION
      #:endfor
   #:endfor

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

! **************************************************************************************************
!> \brief ...
!> \param record ...
!> \param name ...
!> \return ...
! **************************************************************************************************
   FUNCTION openpmd_record_get_component(record, name) RESULT(record_component)
      TYPE(openpmd_record_type)                          :: record
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(openpmd_record_component_type)                :: record_component

      INTEGER(kind=C_INT)                                :: status
      INTERFACE
 FUNCTION openpmd_c_record_get_component(record, name, record_component) RESULT(status) BIND(C, name="openPMD_Record_get_Component")
            IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
            TYPE(C_PTR), VALUE :: record
            CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
            TYPE(C_PTR) :: record_component
            INTEGER(KIND=C_INT) :: status
         END FUNCTION openpmd_c_record_get_component
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(record_component%c_ptr))
      status = openpmd_c_record_get_component(record%c_ptr, name//C_NULL_CHAR, record_component%c_ptr)
      CPASSERT(C_ASSOCIATED(record_component%c_ptr))
   END FUNCTION openpmd_record_get_component

END MODULE openPMD_api
