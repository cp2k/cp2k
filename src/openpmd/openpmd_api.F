!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE openpmd_api

#ifdef __OPENPMD

   USE ISO_C_BINDING, ONLY: &
      C_ASSOCIATED, &
      C_CHAR, &
      C_DOUBLE, &
      C_F_POINTER, &
      C_INT, &
      C_INT64_T, &
      C_LOC, &
      C_NULL_CHAR, &
      C_NULL_PTR, &
      C_PTR, &
      C_SIZE_T
   USE kinds, ONLY: default_string_length, dp, sp
   USE message_passing, ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER :: openpmd_access_create = 0
   INTEGER, PARAMETER :: openpmd_access_read_only = 1

   INTEGER, PARAMETER :: openpmd_type_char = 0
   INTEGER, PARAMETER :: openpmd_type_uchar = 1
   INTEGER, PARAMETER :: openpmd_type_schar = 2
   INTEGER, PARAMETER :: openpmd_type_short = 3
   INTEGER, PARAMETER :: openpmd_type_int = 4
   INTEGER, PARAMETER :: openpmd_type_long = 5
   INTEGER, PARAMETER :: openpmd_type_longlong = 6
   INTEGER, PARAMETER :: openpmd_type_ushort = 7
   INTEGER, PARAMETER :: openpmd_type_uint = 8
   INTEGER, PARAMETER :: openpmd_type_ulong = 9
   INTEGER, PARAMETER :: openpmd_type_ulonglong = 10
   INTEGER, PARAMETER :: openpmd_type_float = 11
   INTEGER, PARAMETER :: openpmd_type_double = 12
   INTEGER, PARAMETER :: openpmd_type_long_double = 13
   INTEGER, PARAMETER :: openpmd_type_cfloat = 14
   INTEGER, PARAMETER :: openpmd_type_cdouble = 15
   INTEGER, PARAMETER :: openpmd_type_clong_double = 16
   INTEGER, PARAMETER :: openpmd_type_string = 17
   INTEGER, PARAMETER :: openpmd_type_vec_char = 18
   INTEGER, PARAMETER :: openpmd_type_vec_short = 19
   INTEGER, PARAMETER :: openpmd_type_vec_int = 20
   INTEGER, PARAMETER :: openpmd_type_vec_long = 21
   INTEGER, PARAMETER :: openpmd_type_vec_longlong = 22
   INTEGER, PARAMETER :: openpmd_type_vec_uchar = 23
   INTEGER, PARAMETER :: openpmd_type_vec_ushort = 24
   INTEGER, PARAMETER :: openpmd_type_vec_uint = 25
   INTEGER, PARAMETER :: openpmd_type_vec_ulong = 26
   INTEGER, PARAMETER :: openpmd_type_vec_ulonglong = 27
   INTEGER, PARAMETER :: openpmd_type_vec_float = 28
   INTEGER, PARAMETER :: openpmd_type_vec_double = 29
   INTEGER, PARAMETER :: openpmd_type_vec_long_double = 30
   INTEGER, PARAMETER :: openpmd_type_vec_cfloat = 31
   INTEGER, PARAMETER :: openpmd_type_vec_cdouble = 32
   INTEGER, PARAMETER :: openpmd_type_vec_clong_double = 33
   INTEGER, PARAMETER :: openpmd_type_vec_schar = 34
   INTEGER, PARAMETER :: openpmd_type_vec_string = 35
   INTEGER, PARAMETER :: openpmd_type_arr_dbl_7 = 36
   INTEGER, PARAMETER :: openpmd_type_bool = 37

   #:set dataset_types_f = ['REAL(kind=dp)']
   #:set dataset_types_enum = ['double']
   #:set dimensions = [1, 2, 3]
   #:set upcasts = [("Series", "Attributable"), ("Iteration", "Attributable"), ("RecordComponent", "Attributable"), ("Mesh", "RecordComponent"), ("Record", "RecordComponent"), ("Mesh", "MeshRecordComponent"), ("MeshRecordComponent", "RecordComponent")]
   #:set fortran_names = {"RecordComponent": "record_component", "MeshRecordComponent": "mesh_record_component"}

   TYPE openpmd_attributable_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: set_attribute_vec_int => openpmd_attributable_set_attribute_vec_int
      PROCEDURE, PUBLIC :: series_flush => openpmd_attributable_series_flush
   END TYPE openpmd_attributable_type

   TYPE openpmd_series_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: as_pointer
      PROCEDURE, PUBLIC :: as_attributable => openpmd_series_as_attributable
      PROCEDURE, PUBLIC :: close => openpmd_series_close
      PROCEDURE, PUBLIC :: present => openpmd_series_present
      PROCEDURE, PUBLIC :: write_iteration => openpmd_series_write_iteration
      PROCEDURE, PUBLIC :: get_iteration => openpmd_series_get_iteration
   END TYPE openpmd_series_type

   TYPE openpmd_iteration_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: as_attributable => openpmd_iteration_as_attributable
      PROCEDURE, PUBLIC :: get_mesh => openpmd_iteration_get_mesh
      PROCEDURE, PUBLIC :: get_particle_species => openpmd_iteration_get_particle_species
      PROCEDURE, PUBLIC :: close => openpmd_iteration_close
      PROCEDURE, PUBLIC :: closed => openpmd_iteration_closed
   END TYPE openpmd_iteration_type

   TYPE openpmd_mesh_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: as_record_component => openpmd_mesh_as_record_component
      PROCEDURE, PUBLIC :: set_axis_labels => openpmd_mesh_set_axis_labels
      PROCEDURE, PUBLIC :: set_grid_global_offset => openpmd_mesh_set_grid_global_offset
      PROCEDURE, PUBLIC :: set_grid_spacing => openpmd_mesh_set_grid_spacing
      PROCEDURE, PUBLIC :: set_position => openpmd_mesh_set_position
   END TYPE openpmd_mesh_type

   TYPE openpmd_particle_species_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: get_record => openpmd_particle_species_get_record
   END TYPE openpmd_particle_species_type

   TYPE openpmd_record_component_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   CONTAINS
      PROCEDURE, PUBLIC :: make_empty => openpmd_record_component_make_empty
      PROCEDURE, PUBLIC :: make_constant_zero => openpmd_record_component_make_constant_zero
      PROCEDURE, PUBLIC :: reset_dataset => openpmd_record_component_reset_dataset
      #:for type_f, type_enum in zip(dataset_types_f, dataset_types_enum)
         #:for dim in dimensions
            PROCEDURE, PRIVATE :: store_chunk_${dim}$d_${type_enum}$ &
               => openpmd_record_component_store_chunk_${dim}$d_${type_enum}$
            PROCEDURE, PUBLIC  :: store_chunk_span_${dim}$d_${type_enum}$ &
               => openpmd_record_component_store_chunk_span_${dim}$d_${type_enum}$
         #:endfor
      #:endfor
      GENERIC, PUBLIC :: store_chunk => &
      #:set zipped = zip(dataset_types_f, dataset_types_enum)
      #:set n = len([x for x in zipped])
      #:set m = len(dimensions)
      #:set count_up_to = m * n
      #:for i, type_f, type_enum in zip(range(n), dataset_types_f, dataset_types_enum)
         #:for j, dim in enumerate(dimensions)
            #:if i * n + j == count_up_to - 1
         store_chunk_${dim}$d_${type_enum}$
         #:else
      store_chunk_${dim}$d_${type_enum}$, &
            #:endif
         #:endfor
      #:endfor
      #:del count_up_to
      #:del zipped
      #:del n
      #:del m
         END TYPE openpmd_record_component_type

      TYPE openpmd_mesh_record_component_type
         PRIVATE
         TYPE(C_PTR) :: c_ptr = C_NULL_PTR
      END TYPE openpmd_mesh_record_component_type

      TYPE openpmd_record_type
         PRIVATE
         TYPE(C_PTR) :: c_ptr = C_NULL_PTR
      CONTAINS
         PROCEDURE, PUBLIC :: as_record_component => openpmd_record_as_record_component
         PROCEDURE, PUBLIC :: get_component => openpmd_record_get_component
      END TYPE openpmd_record_type

      #:for dim in dimensions
         TYPE openpmd_dynamic_memory_view_type_${dim}$d
            PRIVATE
            INTEGER, DIMENSION(${dim}$) :: chunk_extent
            TYPE(C_PTR) :: c_ptr = C_NULL_PTR

         CONTAINS

            #:for type_f, type_enum in zip(dataset_types_f, dataset_types_enum)
               PROCEDURE, PUBLIC :: resolve_${type_enum}$ &
                  => openpmd_dynamic_memory_view_resolve_${dim}$d_${type_enum}$
            #:endfor
         END TYPE openpmd_dynamic_memory_view_type_${dim}$d
      #:endfor

      ! Types
      PUBLIC :: openpmd_attributable_type
      PUBLIC :: openpmd_iteration_type
      PUBLIC :: openpmd_mesh_type
      PUBLIC :: openpmd_particle_species_type
      PUBLIC :: openpmd_record_component_type
      PUBLIC :: openpmd_record_type
      PUBLIC :: openpmd_series_type, openpmd_series_create

      ! Helpers
      PUBLIC :: openpmd_access_create, openpmd_access_read_only
      PUBLIC :: openpmd_get_default_extension
      PUBLIC :: openpmd_json_merge

      ! openPMD datatypes
      PUBLIC :: openpmd_type_arr_dbl_7
      PUBLIC :: openpmd_type_bool
      PUBLIC :: openpmd_type_cdouble
      PUBLIC :: openpmd_type_cfloat
      PUBLIC :: openpmd_type_char
      PUBLIC :: openpmd_type_clong_double
      PUBLIC :: openpmd_type_double
      PUBLIC :: openpmd_type_float
      PUBLIC :: openpmd_type_int
      PUBLIC :: openpmd_type_long
      PUBLIC :: openpmd_type_long_double
      PUBLIC :: openpmd_type_longlong
      PUBLIC :: openpmd_type_schar
      PUBLIC :: openpmd_type_short
      PUBLIC :: openpmd_type_string
      PUBLIC :: openpmd_type_uchar
      PUBLIC :: openpmd_type_uint
      PUBLIC :: openpmd_type_ulong
      PUBLIC :: openpmd_type_ulonglong
      PUBLIC :: openpmd_type_ushort
      PUBLIC :: openpmd_type_vec_cdouble
      PUBLIC :: openpmd_type_vec_cfloat
      PUBLIC :: openpmd_type_vec_char
      PUBLIC :: openpmd_type_vec_clong_double
      PUBLIC :: openpmd_type_vec_double
      PUBLIC :: openpmd_type_vec_float
      PUBLIC :: openpmd_type_vec_int
      PUBLIC :: openpmd_type_vec_long
      PUBLIC :: openpmd_type_vec_long_double
      PUBLIC :: openpmd_type_vec_longlong
      PUBLIC :: openpmd_type_vec_schar
      PUBLIC :: openpmd_type_vec_short
      PUBLIC :: openpmd_type_vec_string
      PUBLIC :: openpmd_type_vec_uchar
      PUBLIC :: openpmd_type_vec_uint
      PUBLIC :: openpmd_type_vec_ulong
      PUBLIC :: openpmd_type_vec_ulonglong
      PUBLIC :: openpmd_type_vec_ushort

      #:for child, parent in upcasts
         #:set child_f = fortran_names[child] if child in fortran_names else child
         #:set parent_f = fortran_names[parent] if parent in fortran_names else parent
         PUBLIC :: openpmd_${child_f}$_as_${parent_f}$
      #:endfor
      #:for dim in dimensions
         PUBLIC :: openpmd_dynamic_memory_view_type_${dim}$d
      #:endfor

#endif

   CONTAINS

#ifdef __OPENPMD

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

      #:for child, parent in upcasts
         #:set child_f = fortran_names[child] if child in fortran_names else child
         #:set parent_f = fortran_names[parent] if parent in fortran_names else parent
         FUNCTION openpmd_${child_f}$_as_${parent_f}$ (this) RESULT(res)
            CLASS(openpmd_${child_f}$_type), INTENT(IN) :: this
            TYPE(openpmd_${parent_f}$_type) :: res

            INTERFACE
               SUBROUTINE openpmd_c_${child_f}$_as_${parent_f}$ (child, parent) &
                  BIND(C, name="openPMD_${child}$_upcast_to_${parent}$")
                  IMPORT :: C_PTR
                  TYPE(C_PTR), VALUE :: child
                  TYPE(C_PTR) :: parent
               END SUBROUTINE openpmd_c_${child_f}$_as_${parent_f}$
            END INTERFACE

            CALL openpmd_c_${child_f}$_as_${parent_f}$ (this%c_ptr, res%c_ptr)
         END FUNCTION openpmd_${child_f}$_as_${parent_f}$
      #:endfor

! **************************************************************************************************
!> \brief ...
!> \param series ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_series_present(this) RESULT(res)
         CLASS(openpmd_series_type), INTENT(IN) :: this
         LOGICAL :: res
         INTEGER :: res_internal

         INTERFACE
            FUNCTION openpmd_c_series_present(series) RESULT(res) &
               BIND(C, NAME="openPMD_Series_present")
               IMPORT :: C_PTR, C_INT
               TYPE(C_PTR), VALUE, INTENT(IN) :: series
               INTEGER(kind=C_INT) :: res
            END FUNCTION openpmd_c_series_present
         END INTERFACE

         IF (.NOT. C_ASSOCIATED(this%c_ptr)) THEN
            res = .FALSE.
            RETURN
         END IF

         res_internal = openpmd_c_series_present(this%c_ptr)
         res = res_internal /= 0
      END FUNCTION openpmd_series_present

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
         CHARACTER(len=*), INTENT(IN), OPTIONAL, TARGET     :: config
         TYPE(openpmd_series_type)                          :: series

         TYPE(c_ptr), SAVE                                  :: comm_c = C_NULL_PTR
         CHARACTER(len=default_string_length), TARGET, SAVE :: default_config = "{}"
         CHARACTER(:), POINTER                              :: my_config
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
            END SUBROUTINE cp2k_c_mpi_comm_f2c
         END INTERFACE
         INTERFACE
            SUBROUTINE openpmd_c_series_create_mpi(path, access, mpi_comm, config, series) &
               BIND(C, name="openPMD_Series_create_mpi")
               IMPORT :: C_PTR, C_CHAR, C_INT
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
               INTEGER(kind=C_INT), VALUE :: access
               TYPE(C_PTR), value :: mpi_comm
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: config
               TYPE(C_PTR) :: series
            END SUBROUTINE openpmd_c_series_create_mpi
         END INTERFACE
         INTERFACE
            SUBROUTINE openpmd_c_series_create(path, access, config, series) &
               BIND(C, name="openPMD_Series_create")
               IMPORT :: C_PTR, C_CHAR, C_INT
               TYPE(C_PTR) :: series
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: path
               INTEGER(kind=C_INT), VALUE :: access
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: config
            END SUBROUTINE openpmd_c_series_create
         END INTERFACE

         IF (PRESENT(config)) THEN
            my_config => config
         ELSE
            my_config => default_config
         END IF

         CPASSERT(.NOT. C_ASSOCIATED(series%c_ptr))
         IF (PRESENT(mpi_comm)) THEN
            CALL cp2k_c_mpi_comm_f2c(mpi_comm%get_handle(), comm_c)
            CALL openpmd_c_series_create_mpi( &
               path=TRIM(path)//C_NULL_CHAR, &
               access=access, &
               mpi_comm=comm_c, &
               config=my_config//C_NULL_CHAR, &
               series=series%c_ptr)
         ELSE
            CALL openpmd_c_series_create( &
               path=TRIM(path)//C_NULL_CHAR, &
               access=access, &
               config=my_config//C_NULL_CHAR, &
               series=series%c_ptr)
         END IF
         CPASSERT(C_ASSOCIATED(series%c_ptr))
      END FUNCTION openpmd_series_create

! **************************************************************************************************
!> \brief ...
!> \param series ...
! **************************************************************************************************
      SUBROUTINE openpmd_series_close(series)
         CLASS(openpmd_series_type), INTENT(INOUT)           :: series

         INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param series ...
! **************************************************************************************************
            SUBROUTINE openpmd_c_series_close(series) &
               BIND(C, name="openPMD_Series_close")
               IMPORT :: C_PTR
               TYPE(C_PTR), VALUE :: series
            END SUBROUTINE openpmd_c_series_close
         END INTERFACE

         CPASSERT(C_ASSOCIATED(series%c_ptr))
         CALL openpmd_c_series_close(series=series%c_ptr)
         series%c_ptr = C_NULL_PTR
      END SUBROUTINE openpmd_series_close

! **************************************************************************************************
!> \brief ...
!> \param series ...
!> \param index ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_series_write_iteration(this, index) RESULT(iteration)
         CLASS(openpmd_series_type)                          :: this
         INTEGER, TARGET                                    :: index
         TYPE(openpmd_iteration_type)                       :: iteration

         INTEGER(8)                                         :: cast_index
         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_series_write_iteration(series, index, iteration) RESULT(status) &
               BIND(C, name="openPMD_Series_write_Iteration")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T
               TYPE(C_PTR), VALUE :: series
               INTEGER(kind=C_INT64_T), VALUE :: index
               TYPE(C_PTR) :: iteration
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_series_write_iteration
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(iteration%c_ptr))
         cast_index = index
         status = openpmd_c_series_write_iteration(this%c_ptr, cast_index, iteration%c_ptr)
         CPASSERT(C_ASSOCIATED(iteration%c_ptr))
      END FUNCTION openpmd_series_write_iteration

      FUNCTION openpmd_series_get_iteration(this, index) RESULT(iteration)
         CLASS(openpmd_series_type)                          :: this
         INTEGER, TARGET                                    :: index
         TYPE(openpmd_iteration_type)                       :: iteration

         INTEGER(8)                                         :: cast_index
         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_series_get_iteration(series, index, iteration) RESULT(status) &
               BIND(C, name="openPMD_Series_get_Iteration")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T
               TYPE(C_PTR), VALUE :: series
               INTEGER(kind=C_INT64_T), VALUE :: index
               TYPE(C_PTR) :: iteration
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_series_get_iteration
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(iteration%c_ptr))
         cast_index = index
         status = openpmd_c_series_get_iteration(this%c_ptr, cast_index, iteration%c_ptr)
         CPASSERT(C_ASSOCIATED(iteration%c_ptr))
      END FUNCTION openpmd_series_get_iteration

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_get_default_extension() RESULT(extension)
         CHARACTER(len=default_string_length)               :: extension

         CHARACTER(kind=C_CHAR), POINTER                    :: c_extension
         INTEGER(C_SIZE_T)                                  :: length_of_c_string, i
         INTERFACE
            FUNCTION openpmd_c_get_default_extension() RESULT(extension) &
               BIND(c, name="openPMD_get_default_extension")
               IMPORT :: C_CHAR
               CHARACTER(kind=C_CHAR), POINTER :: extension
            END FUNCTION openpmd_c_get_default_extension
         END INTERFACE
         INTERFACE
            FUNCTION get_strlen(c_string) RESULT(strlen) &
               BIND(C, NAME="CP2K_strlen")
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
         extension = TRIM(ADJUSTL(TRANSFER(c_extension(1:length_of_c_string), extension)))

         DO i = length_of_c_string + 1, LEN(extension)
            extension(i:i) = ' '
         END DO

      END FUNCTION openpmd_get_default_extension

! **************************************************************************************************
!> \brief ...
!> \param attributable ...
!> \param name ...
!> \param vec ...
! **************************************************************************************************
      SUBROUTINE openpmd_attributable_set_attribute_vec_int(this, name, vec)
         CLASS(openpmd_attributable_type)                    :: this
         CHARACTER(len=*), INTENT(IN)                       :: name
         INTEGER, DIMENSION(:), TARGET                      :: vec

         INTEGER, TARGET                                    :: size_
         INTERFACE
            SUBROUTINE openpmd_c_attributable_set_attribute_vec_int(attributable, name, vec, length) &
               BIND(C, NAME="openPMD_attributable_set_attribute_vec_int")
               IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
               TYPE(C_PTR), VALUE :: attributable
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
               TYPE(C_PTR), VALUE :: vec
               INTEGER(kind=C_INT), VALUE :: length
            END SUBROUTINE openpmd_c_attributable_set_attribute_vec_int
         END INTERFACE

         size_ = SIZE(vec)
         CALL openpmd_c_attributable_set_attribute_vec_int( &
            this%c_ptr, name//C_NULL_CHAR, C_LOC(vec(1)), size_)

      END SUBROUTINE openpmd_attributable_set_attribute_vec_int

      SUBROUTINE openpmd_attributable_series_flush(this, backendconfig)
         CLASS(openpmd_attributable_type)                   :: this
         CHARACTER(len=*), INTENT(IN), OPTIONAL             :: backendconfig
         CHARACTER(len=default_string_length), SAVE         :: my_backendconfig = ""

         INTERFACE
            SUBROUTINE openpmd_c_attributable_series_flush(attributable, config) &
               BIND(C, NAME="openPMD_Attributable_series_flush")
               IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
               TYPE(C_PTR), VALUE :: attributable
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: config
            END SUBROUTINE openpmd_c_attributable_series_flush
         END INTERFACE

         IF (PRESENT(backendconfig)) my_backendconfig = backendconfig
         CALL openpmd_c_attributable_series_flush(this%c_ptr, my_backendconfig//C_NULL_CHAR)
      END SUBROUTINE openpmd_attributable_series_flush

      FUNCTION openpmd_json_merge(into, from, mpi_comm) RESULT(merged)
         CHARACTER(len=*), INTENT(IN)                       :: into, from
         CLASS(mp_comm_type), INTENT(in), OPTIONAL          :: mpi_comm
         CHARACTER(:), ALLOCATABLE                          :: merged
         CHARACTER(kind=C_CHAR), POINTER                    :: c_merged
         INTEGER(C_SIZE_T)                                  :: length_of_c_string
         TYPE(c_ptr), SAVE                                  :: comm_c = C_NULL_PTR

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
            END SUBROUTINE cp2k_c_mpi_comm_f2c
         END INTERFACE
         INTERFACE
            FUNCTION openpmd_c_json_merge(into, from, mpi_comm) RESULT(merged) &
               BIND(C, NAME="openPMD_json_merge")
               IMPORT :: C_CHAR, C_INT, C_PTR
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: into, from
               TYPE(C_PTR), value                   :: mpi_comm
               CHARACTER(kind=C_CHAR), POINTER      :: merged
            END FUNCTION openpmd_c_json_merge
         END INTERFACE
         INTERFACE
            FUNCTION get_strlen(c_string) RESULT(strlen) &
               BIND(C, NAME="CP2K_strlen")
               IMPORT :: C_CHAR, C_SIZE_T
               CHARACTER(kind=C_CHAR), DIMENSION(*), INTENT(IN) :: c_string
               INTEGER(C_SIZE_T) :: strlen
            END FUNCTION get_strlen
         END INTERFACE
         INTERFACE
            SUBROUTINE free(c_string) &
               BIND(C, NAME="CP2K_free")
               IMPORT :: C_CHAR
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: c_string
            END SUBROUTINE free
         END INTERFACE

         IF (PRESENT(mpi_comm)) THEN
            CALL cp2k_c_mpi_comm_f2c(mpi_comm%get_handle(), comm_c)
         END IF
         c_merged => openpmd_c_json_merge(into//C_NULL_CHAR, from//C_NULL_CHAR, comm_c)
         length_of_c_string = get_strlen(c_merged)
         ALLOCATE (character(length_of_c_string) :: merged)
         merged = TRIM(ADJUSTL(TRANSFER(c_merged(1:length_of_c_string), MERGED)))
         CALL free(c_merged)
      END FUNCTION openpmd_json_merge

! **************************************************************************************************
!> \brief ...
!> \param iteration ...
!> \param name ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_iteration_get_mesh(this, name) RESULT(mesh)
         CLASS(openpmd_iteration_type)                       :: this
         CHARACTER(len=*), INTENT(IN)                       :: name
         TYPE(openpmd_mesh_type)                            :: mesh

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_iteration_get_mesh(iteration, name, mesh) RESULT(status) &
               BIND(C, name="openPMD_Iteration_get_mesh")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: iteration
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
               TYPE(C_PTR) :: mesh
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_iteration_get_mesh
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(mesh%c_ptr))
         status = openpmd_c_iteration_get_mesh(this%c_ptr, name//C_NULL_CHAR, mesh%c_ptr)
         CPASSERT(C_ASSOCIATED(mesh%c_ptr))
      END FUNCTION openpmd_iteration_get_mesh

! **************************************************************************************************
!> \brief ...
!> \param iteration ...
!> \param name ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_iteration_get_particle_species(this, name) RESULT(particle_species)
         CLASS(openpmd_iteration_type)                       :: this
         CHARACTER(len=*), INTENT(IN)                       :: name
         TYPE(openpmd_particle_species_type)                :: particle_species

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_iteration_get_particle_species(iteration, name, particle_species) RESULT(status) &
               BIND(C, name="openPMD_Iteration_get_particle_species")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: iteration
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
               TYPE(C_PTR) :: particle_species
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_iteration_get_particle_species
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(particle_species%c_ptr))
         status = openpmd_c_iteration_get_particle_species( &
                  this%c_ptr, name//C_NULL_CHAR, particle_species%c_ptr)
         CPASSERT(C_ASSOCIATED(particle_species%c_ptr))
      END FUNCTION openpmd_iteration_get_particle_species

! **************************************************************************************************
!> \brief ...
!> \param iteration ...
! **************************************************************************************************
      SUBROUTINE openpmd_iteration_close(iteration)
         CLASS(openpmd_iteration_type), INTENT(INOUT)           :: iteration

         INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param iteration ...
! **************************************************************************************************
            SUBROUTINE openpmd_c_iteration_close(iteration) &
               BIND(C, name="openPMD_Iteration_close")
               IMPORT :: C_PTR
               TYPE(C_PTR), VALUE :: iteration
            END SUBROUTINE openpmd_c_iteration_close
         END INTERFACE

         CPASSERT(C_ASSOCIATED(iteration%c_ptr))
         CALL openpmd_c_iteration_close(iteration=iteration%c_ptr)
         iteration%c_ptr = C_NULL_PTR
      END SUBROUTINE openpmd_iteration_close

      FUNCTION openpmd_iteration_closed(iteration) RESULT(res)
         CLASS(openpmd_iteration_type), INTENT(IN)           :: iteration
         LOGICAL :: res

         INTEGER(kind=C_INT) :: res_tmp

         INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param iteration ...
! **************************************************************************************************
            FUNCTION openpmd_c_iteration_closed(iteration) RESULT(res) &
               BIND(C, name="openPMD_Iteration_closed")
               IMPORT :: C_PTR, C_INT
               TYPE(C_PTR), VALUE :: iteration
               INTEGER(kind=C_INT) :: res
            END FUNCTION openpmd_c_iteration_closed
         END INTERFACE

         CPASSERT(C_ASSOCIATED(iteration%c_ptr))
         res_tmp = openpmd_c_iteration_closed(iteration=iteration%c_ptr)
         IF (res_tmp == 0) THEN
            res = .FALSE.
         ELSE
            res = .TRUE.
         END IF
      END FUNCTION openpmd_iteration_closed

! **************************************************************************************************
!> \brief ...
!> \param record_component ...
!> \param dtype ...
!> \param dimensionality ...
! **************************************************************************************************
      SUBROUTINE openpmd_record_component_make_empty(this, dtype, dimensionality)
         CLASS(openpmd_record_component_type)                :: this
         INTEGER(kind=C_INT)                                :: dtype
         INTEGER                                            :: dimensionality

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_record_component_make_empty(record_component, dtype, dimensionality) RESULT(status) &
               BIND(C, name="openPMD_RecordComponent_makeEmpty")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: record_component
               INTEGER(KIND=C_INT), VALUE :: dtype
               INTEGER(KIND=C_INT), VALUE :: dimensionality
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_record_component_make_empty
         END INTERFACE

         status = openpmd_c_record_component_make_empty(this%c_ptr, dtype, dimensionality)
      END SUBROUTINE openpmd_record_component_make_empty

! **************************************************************************************************
!> \brief ...
!> \param record_component ...
!> \param dtype ...
!> \param extent ...
! **************************************************************************************************
      SUBROUTINE openpmd_record_component_make_constant_zero(this, dtype, extent)
         CLASS(openpmd_record_component_type)                :: this
         INTEGER(kind=C_INT)                                :: dtype
         INTEGER, DIMENSION(:), TARGET                      :: extent

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_record_component_make_constant( &
               record_component, dtype, dimensions, extent, invert, value) RESULT(status) &
               BIND(C, name="openPMD_RecordComponent_makeConstant")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: record_component
               INTEGER(KIND=C_INT), VALUE :: dtype
               INTEGER(KIND=C_INT), VALUE :: dimensions, invert
               TYPE(C_PTR), VALUE :: extent
               TYPE(C_PTR), VALUE :: value
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_record_component_make_constant
         END INTERFACE

         status = openpmd_c_record_component_make_constant( &
                  this%c_ptr, dtype, SIZE(extent), C_LOC(extent(1)), 1, C_NULL_PTR)
      END SUBROUTINE openpmd_record_component_make_constant_zero

      SUBROUTINE openpmd_mesh_set_axis_labels(this, labels)
         CLASS(openpmd_mesh_type) :: this
         CHARACTER(len=*), DIMENSION(:) :: labels

         INTEGER :: i
         CHARACTER(len=default_string_length), ALLOCATABLE, TARGET :: labels_as_null_terminated_strings(:)
         TYPE(C_PTR), ALLOCATABLE, TARGET :: labels_as_c_strings(:)

         INTERFACE
            SUBROUTINE openpmd_c_mesh_set_axis_labels(mesh, labels, len_labels, invert) &
               BIND(C, NAME="openPMD_Mesh_set_axis_labels")
               IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
               TYPE(C_PTR), VALUE :: mesh
               TYPE(C_PTR), VALUE :: labels
               INTEGER(kind=C_INT), VALUE :: len_labels, invert
            END SUBROUTINE openpmd_c_mesh_set_axis_labels
         END INTERFACE

         ALLOCATE (labels_as_null_terminated_strings(SIZE(labels)))
         ALLOCATE (labels_as_c_strings(SIZE(labels)))

         DO i = 1, SIZE(labels)
            labels_as_null_terminated_strings(i) = TRIM(labels(i))//C_NULL_CHAR
            labels_as_c_strings(i) = C_LOC(labels_as_null_terminated_strings(i))
         END DO

         CALL openpmd_c_mesh_set_axis_labels(this%c_ptr, C_LOC(labels_as_c_strings(1)), SIZE(labels), invert=1)

         DEALLOCATE (labels_as_c_strings)
         DEALLOCATE (labels_as_null_terminated_strings)
      END SUBROUTINE openpmd_mesh_set_axis_labels

      SUBROUTINE openpmd_record_component_reset_dataset(this, dtype, extent, cfg)
         CLASS(openpmd_record_component_type)               :: this
         INTEGER(kind=C_INT)                                :: dtype
         INTEGER, DIMENSION(:), TARGET                      :: extent
         CHARACTER(len=*), OPTIONAL :: cfg

         CHARACTER(len=default_string_length), SAVE         :: my_config = "{}"
         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_record_component_reset_dataset( &
               record_component, dtype, dimensions, extent, invert, cfg) RESULT(status) &
               BIND(C, name="openPMD_RecordComponent_resetDataset")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: record_component
               INTEGER(KIND=C_INT), VALUE :: dtype
               INTEGER(KIND=C_INT), VALUE :: dimensions, invert
               TYPE(C_PTR), VALUE :: extent
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: cfg
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_record_component_reset_dataset
         END INTERFACE

         IF (PRESENT(cfg)) my_config = cfg

         status = openpmd_c_record_component_reset_dataset( &
                  this%c_ptr, dtype, SIZE(extent), C_LOC(extent(1)), 1, my_config//C_NULL_CHAR)
      END SUBROUTINE openpmd_record_component_reset_dataset

      #:for type_f, type_enum in zip(dataset_types_f, dataset_types_enum)
         #:for dim in [1, 2, 3]
            ! e.g. ':,:,:'
            #:set dimensionality_wildcard = ','.join(':' for _ in range(dim))
            SUBROUTINE openpmd_record_component_store_chunk_${dim}$d_${type_enum}$ (record_component, data_, offset)
               CLASS(openpmd_record_component_type)                        :: record_component
               ${type_f}$, TARGET                                         :: data_(${dimensionality_wildcard}$)
               INTEGER, DIMENSION(${dim}$), TARGET                        :: offset, extent

               INTEGER(kind=C_INT)                                        :: status

               INTERFACE
                  FUNCTION openpmd_c_record_component_store_chunk( &
                     record_component, dtype, dimensions, offset, extent, invert, data_) RESULT(status) &
                     BIND(C, name="openPMD_RecordComponent_storeChunk")
                     IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
                     TYPE(C_PTR), VALUE :: record_component
                     INTEGER(KIND=C_INT), VALUE :: dtype
                     INTEGER(KIND=C_INT), VALUE :: dimensions, invert
                     TYPE(C_PTR), VALUE :: offset, extent
                     TYPE(C_PTR), VALUE :: data_
                     INTEGER(KIND=C_INT) :: status
                  END FUNCTION openpmd_c_record_component_store_chunk
               END INTERFACE

               extent = SHAPE(data_)
               status = openpmd_c_record_component_store_chunk( &
                        record_component%c_ptr, &
                        openpmd_type_${type_enum}$, &
                        ${dim}$, &
                        C_LOC(offset(1)), &
                        C_LOC(extent(1)), &
                        1, &
                        C_LOC(data_))

            END SUBROUTINE openpmd_record_component_store_chunk_${dim}$d_${type_enum}$

            FUNCTION openpmd_record_component_store_chunk_span_${dim}$d_${type_enum}$ ( &
               record_component, offset, extent) RESULT(memory_view)
               CLASS(openpmd_record_component_type)                        :: record_component
               INTEGER, DIMENSION(${dim}$), TARGET                        :: offset, extent

               TYPE(openpmd_dynamic_memory_view_type_${dim}$d)            :: memory_view

               INTEGER(kind=C_INT)                                        :: status

               INTERFACE
                  FUNCTION openpmd_c_record_component_store_chunk_span( &
                     record_component, dtype, dimensions, offset, extent, invert, memory_view) RESULT(status) &
                     BIND(C, name="openPMD_RecordComponent_storeChunkSpan")
                     IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
                     TYPE(C_PTR), VALUE :: record_component
                     INTEGER(KIND=C_INT), VALUE :: dtype
                     INTEGER(KIND=C_INT), VALUE :: dimensions, invert
                     TYPE(C_PTR), VALUE :: offset, extent
                     TYPE(C_PTR) :: memory_view
                     INTEGER(KIND=C_INT) :: status
                  END FUNCTION openpmd_c_record_component_store_chunk_span
               END INTERFACE

               memory_view%chunk_extent = extent
               status = openpmd_c_record_component_store_chunk_span( &
                        record_component%c_ptr, &
                        openpmd_type_${type_enum}$, &
                        ${dim}$, &
                        C_LOC(offset(1)), &
                        C_LOC(extent(1)), &
                        1, &
                        memory_view%c_ptr)

            END FUNCTION openpmd_record_component_store_chunk_span_${dim}$d_${type_enum}$

            FUNCTION openpmd_dynamic_memory_view_resolve_${dim}$d_${type_enum}$ (memory_view, deallocate) RESULT(buffer)
               CLASS(openpmd_dynamic_memory_view_type_${dim}$d)            :: memory_view
               LOGICAL                                                    :: deallocate
               ${type_f}$, POINTER                                        :: buffer(${dimensionality_wildcard}$)

               INTEGER(kind=C_INT)                                        :: status, deallocate_c
               TYPE(C_PTR), TARGET                                        :: c_buffer

               INTERFACE
                  FUNCTION openpmd_c_dynamic_memory_view_resolve(memory_view, deallocate, write_buffer) RESULT(status) &
                     BIND(C, name="openPMD_DynamicMemoryView_resolve")
                     IMPORT :: C_PTR, C_INT
                     TYPE(C_PTR), VALUE :: memory_view
                     INTEGER(KIND=C_INT), VALUE :: deallocate
                     TYPE(C_PTR) :: write_buffer
                     INTEGER(KIND=C_INT) :: status
                  END FUNCTION openpmd_c_dynamic_memory_view_resolve
               END INTERFACE

               IF (deallocate) THEN
                  deallocate_c = 1
               ELSE
                  deallocate_c = 0
               END IF

               status = openpmd_c_dynamic_memory_view_resolve(memory_view%c_ptr, deallocate_c, c_buffer)
               CALL c_f_pointer(c_buffer, buffer, memory_view%chunk_extent)
            END FUNCTION openpmd_dynamic_memory_view_resolve_${dim}$d_${type_enum}$
         #:endfor
      #:endfor

! **************************************************************************************************
!> \brief ...
!> \param particle_species ...
!> \param name ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_particle_species_get_record(this, name) RESULT(record)
         CLASS(openpmd_particle_species_type)                :: this
         CHARACTER(len=*), INTENT(IN)                       :: name
         TYPE(openpmd_record_type)                          :: record

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_particle_species_get_record(particle_species, name, record) RESULT(status) &
               BIND(C, name="openPMD_ParticleSpecies_get_Record")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: particle_species
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
               TYPE(C_PTR) :: record
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_particle_species_get_record
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(record%c_ptr))
         status = openpmd_c_particle_species_get_record(this%c_ptr, name//C_NULL_CHAR, record%c_ptr)
         CPASSERT(C_ASSOCIATED(record%c_ptr))
      END FUNCTION openpmd_particle_species_get_record

! **************************************************************************************************
!> \brief ...
!> \param record ...
!> \param name ...
!> \return ...
! **************************************************************************************************
      FUNCTION openpmd_record_get_component(this, name) RESULT(record_component)
         CLASS(openpmd_record_type)                          :: this
         CHARACTER(len=*), INTENT(IN)                       :: name
         TYPE(openpmd_record_component_type)                :: record_component

         INTEGER(kind=C_INT)                                :: status
         INTERFACE
            FUNCTION openpmd_c_record_get_component(record, name, record_component) RESULT(status) &
               BIND(C, name="openPMD_Record_get_Component")
               IMPORT :: C_PTR, C_INT64_T, C_INT, C_SIZE_T, C_CHAR
               TYPE(C_PTR), VALUE :: record
               CHARACTER(kind=C_CHAR), DIMENSION(*) :: name
               TYPE(C_PTR) :: record_component
               INTEGER(KIND=C_INT) :: status
            END FUNCTION openpmd_c_record_get_component
         END INTERFACE

         CPASSERT(.NOT. C_ASSOCIATED(record_component%c_ptr))
         status = openpmd_c_record_get_component(this%c_ptr, name//C_NULL_CHAR, record_component%c_ptr)
         CPASSERT(C_ASSOCIATED(record_component%c_ptr))
      END FUNCTION openpmd_record_get_component

      #:set f_routines = ["openpmd_mesh_set_grid_global_offset", "openpmd_mesh_set_grid_spacing", "openpmd_mesh_set_position"]
      #:set c_functions = ["openPMD_Mesh_setGridGlobalOffset", "openPMD_Mesh_setGridSpacing", "openPMD_Mesh_setPosition"]
      #:set f_vector_types = ["REAL(KIND=dp)" for _ in range(3)]

      #:for f_routine, c_function, f_vector_type in zip(f_routines, c_functions, f_vector_types)
         SUBROUTINE ${f_routine}$ (this, vec)
         CLASS(openpmd_mesh_type)                        :: this
         ${f_vector_type}$, DIMENSION(:), TARGET        :: vec

         INTEGER, TARGET                                :: size_
         INTERFACE
            SUBROUTINE openpmd_c_set_vector_attribute(mesh, vec, length, invert) &
               BIND(C, NAME="${c_function}$")
               IMPORT :: C_CHAR, C_SIZE_T, C_INT, C_PTR
               TYPE(C_PTR), VALUE :: mesh
               TYPE(C_PTR), VALUE :: vec
               INTEGER(kind=C_INT), VALUE :: length, invert
            END SUBROUTINE openpmd_c_set_vector_attribute
         END INTERFACE

         size_ = SIZE(vec)
         CALL openpmd_c_set_vector_attribute(this%c_ptr, C_LOC(vec(1)), size_, 1)

         END SUBROUTINE ${f_routine}$
         #:endfor

#endif
         ! In two lines of this file, `make pretty` applies a wrong indentation and keeps it for the rest of this file.
         ! Hence, we end up a bit further to the right than we should.
      END MODULE openpmd_api
