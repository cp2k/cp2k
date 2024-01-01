!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A wrapper around the HDF5 Fortran API
!> \par History
!>      04.2023 created [SB]
!> \author Stefano Battaglia
! **************************************************************************************************
MODULE hdf5_wrapper

#ifdef __HDF5
   USE hdf5, ONLY: &
      h5aclose_f, h5acreate_f, h5aopen_f, h5aread_f, h5awrite_f, h5close_f, h5dclose_f, &
      h5dcreate_f, h5dget_space_f, h5dopen_f, h5dread_f, h5dwrite_f, h5f_acc_rdonly_f, &
      h5f_acc_trunc_f, h5fclose_f, h5fcreate_f, h5fopen_f, h5gclose_f, h5gcreate_f, h5gopen_f, &
      h5open_f, h5s_scalar_f, h5sclose_f, h5screate_f, h5screate_simple_f, &
      h5sget_simple_extent_npoints_f, h5t_c_s1, h5t_cset_utf8_f, h5t_enum_f, h5t_native_double, &
      h5t_native_integer, h5t_str_nullpad_f, h5t_string, h5tclose_f, h5tcopy_f, h5tcreate_f, &
      h5tenum_insert_f, h5tset_cset_f, h5tset_size_f, h5tset_strpad_f, hid_t, hsize_t, size_t
#endif
   USE iso_c_binding, ONLY: C_LOC, &
                            c_ptr
   USE kinds, ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hdf5_wrapper'
#ifdef __HDF5
   INTEGER, PARAMETER, PUBLIC           :: hdf5_id = hid_t
#endif

CONTAINS

#ifdef __HDF5
! **************************************************************************************************
!> \brief Initialize the HDF5 fortran API
! **************************************************************************************************
   SUBROUTINE h5open()
      INTEGER                                            :: error

      CALL h5open_f(error)
      IF (error < 0) CPABORT('ERROR: failed to initialize HDF5 interface')

   END SUBROUTINE h5open

! **************************************************************************************************
!> \brief Close the HDF5 fortran API
! **************************************************************************************************
   SUBROUTINE h5close()
      INTEGER                                            :: error

      CALL h5close_f(error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 interface')

   END SUBROUTINE h5close

! **************************************************************************************************
!> \brief Create a HDF5 file
!> \param filename the name of the hdf5 file
!> \param file_id the file id of the hdf5 file
! **************************************************************************************************
   SUBROUTINE h5fcreate(filename, file_id)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER(KIND=hid_t), INTENT(OUT)                   :: file_id

      INTEGER                                            :: error

      CALL h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 file')

   END SUBROUTINE h5fcreate

! **************************************************************************************************
!> \brief Open a HDF5 file
!> \param filename the name of the hdf5 file
!> \param file_id the file id of the hdf5 file
! **************************************************************************************************
   SUBROUTINE h5fopen(filename, file_id)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER(KIND=hid_t), INTENT(OUT)                   :: file_id

      INTEGER                                            :: error

      CALL h5fopen_f(TRIM(filename), h5f_acc_rdonly_f, file_id, error)
      IF (error < 0) CPABORT('ERROR: failed to open HDF5 file')

   END SUBROUTINE h5fopen

! **************************************************************************************************
!> \brief Close a HDF5 file
!> \param file_id the file id of the hdf5 file
! **************************************************************************************************
   SUBROUTINE h5fclose(file_id)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: file_id

      INTEGER                                            :: error

      CALL h5fclose_f(file_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 file')

   END SUBROUTINE h5fclose

! **************************************************************************************************
!> \brief Create a HDF5 group
!> \param loc_id file or group identifier
!> \param name name of the group
!> \param grp_id group identifier
! **************************************************************************************************
   SUBROUTINE h5gcreate(loc_id, name, grp_id)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER(KIND=hid_t), INTENT(OUT)                   :: grp_id

      INTEGER                                            :: error

      CALL h5gcreate_f(loc_id, name, grp_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 group')

   END SUBROUTINE h5gcreate

! **************************************************************************************************
!> \brief Open a HDF5 group
!> \param loc_id file or group identifier
!> \param name name of the group
!> \param grp_id group identifier
! **************************************************************************************************
   SUBROUTINE h5gopen(loc_id, name, grp_id)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER(KIND=hid_t), INTENT(OUT)                   :: grp_id

      INTEGER                                            :: error

      CALL h5gopen_f(loc_id, name, grp_id, error)
      IF (error < 0) CPABORT('ERROR: failed to open HDF5 group')

   END SUBROUTINE h5gopen

! **************************************************************************************************
!> \brief Close a HDF5 group
!> \param grp_id group identifier
! **************************************************************************************************
   SUBROUTINE h5gclose(grp_id)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: grp_id

      INTEGER                                            :: error

      CALL h5gclose_f(grp_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 group')

   END SUBROUTINE h5gclose

! **************************************************************************************************
!> \brief Write a variable-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the string to write
! **************************************************************************************************
   SUBROUTINE h5awrite_varlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET               :: attr_data

      INTEGER                                            :: error, output_unit
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      TYPE(c_ptr)                                        :: buffer
      TYPE(c_ptr), TARGET                                :: in_between_ptr

      ! create a scalar dataspace
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 dataspace'
         RETURN
      END IF

      ! create a variable-length string type
      CALL h5tcopy_f(h5t_string, type_id, error)
      CALL h5tset_cset_f(type_id, h5t_cset_utf8_f, error)
      CALL h5tset_strpad_f(type_id, h5t_str_nullpad_f, error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 attribute'
         RETURN
      END IF

      ! weird in-between pointer needed for variable-length
      ! string to a scalar dataspace
      in_between_ptr = C_LOC(attr_data)
      ! the actual pointer to be passed
      buffer = C_LOC(in_between_ptr)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to write HDF5 attribute'
         RETURN
      END IF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 attribute'
         RETURN
      END IF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 dataspace'
         RETURN
      END IF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 datatype'
         RETURN
      END IF

   END SUBROUTINE h5awrite_varlen_string

! **************************************************************************************************
!> \brief Write a fixed-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the string to write
! **************************************************************************************************
   SUBROUTINE h5awrite_fixlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET               :: attr_data

      INTEGER                                            :: error, output_unit
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      TYPE(c_ptr)                                        :: buffer

      ! create a scalar dataspace
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 dataspace'
         RETURN
      END IF

      ! create a fixed-length string datatype
      CALL h5tcopy_f(h5t_c_s1, type_id, error)
      CALL h5tset_cset_f(type_id, h5t_cset_utf8_f, error)
      CALL h5tset_size_f(type_id, LEN(attr_data, size_t), error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 attribute'
         RETURN
      END IF

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to write HDF5 attribute'
         RETURN
      END IF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 attribute'
         RETURN
      END IF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 dataspace'
         RETURN
      END IF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 datatype'
         RETURN
      END IF

   END SUBROUTINE h5awrite_fixlen_string

! **************************************************************************************************
!> \brief Write a boolean attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the logical to write (.true. or .false.)
! **************************************************************************************************
   SUBROUTINE h5awrite_boolean(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      LOGICAL, INTENT(IN)                                :: attr_data

      INTEGER                                            :: error, output_unit
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      INTEGER(KIND=size_t)                               :: enum_size = 1
      INTEGER, TARGET                                    :: attr_data_to_int
      TYPE(c_ptr)                                        :: buffer

      ! 8-bit integers in enum bool_type

      ! create a scalar dataspace
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 dataspace'
         RETURN
      END IF

      ! create the datatype
      CALL h5tcreate_f(h5t_enum_f, enum_size, type_id, error)
      CALL h5tenum_insert_f(type_id, "FALSE", 0, error)
      CALL h5tenum_insert_f(type_id, "TRUE", 1, error)

      IF (attr_data) THEN
         attr_data_to_int = 1
      ELSE
         attr_data_to_int = 0
      END IF
      ! the C pointer to the actual data
      buffer = C_LOC(attr_data_to_int)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 attribute'
         RETURN
      END IF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to write HDF5 attribute'
         RETURN
      END IF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 attribute'
         RETURN
      END IF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 dataspace'
         RETURN
      END IF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 datatype'
         RETURN
      END IF

   END SUBROUTINE h5awrite_boolean

! **************************************************************************************************
!> \brief Write a (scalar) integer attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the integer to write
! **************************************************************************************************
   SUBROUTINE h5awrite_integer_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      INTEGER, INTENT(IN), TARGET                        :: attr_data

      INTEGER                                            :: error, output_unit
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      TYPE(c_ptr)                                        :: buffer

      ! create a scalar dataspace
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 dataspace'
         RETURN
      END IF

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = h5t_native_integer

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to create HDF5 attribute'
         RETURN
      END IF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to write HDF5 attribute'
         RETURN
      END IF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 attribute'
         RETURN
      END IF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
            ' ERROR: failed to close HDF5 dataspace'
         RETURN
      END IF

   END SUBROUTINE h5awrite_integer_scalar

! **************************************************************************************************
!> \brief Write a (scalar) double precision attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the double to write
! **************************************************************************************************
   SUBROUTINE h5awrite_double_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      REAL(KIND=dp), INTENT(IN), TARGET                  :: attr_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      TYPE(c_ptr)                                        :: buffer

      ! create a scalar dataspace
      CALL h5screate_f(h5s_scalar_f, space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataspace')

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = h5t_native_double

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 attribute')

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to write HDF5 attribute')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 attribute')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

   END SUBROUTINE h5awrite_double_scalar

! **************************************************************************************************
!> \brief Write an array of fixed-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of strings
! **************************************************************************************************
   SUBROUTINE h5awrite_string_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN), TARGET :: attr_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      INTEGER(KIND=hsize_t), DIMENSION(2)                :: dims
      TYPE(c_ptr)                                        :: buffer

      dims(1) = LEN(attr_data(1), kind=hsize_t) ! length of a string entry
      dims(2) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! create a fixed-length string datatype
      CALL h5tcopy_f(h5t_c_s1, type_id, error)
      CALL h5tset_cset_f(type_id, h5t_cset_utf8_f, error)
      CALL h5tset_size_f(type_id, INT(dims(1), size_t), error)

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims(2:2), space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataspace')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 attribute')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to write HDF5 attribute')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 attribute')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 datatype')

   END SUBROUTINE h5awrite_string_simple

! **************************************************************************************************
!> \brief Write an array of doubles attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of doubles
! **************************************************************************************************
   SUBROUTINE h5awrite_double_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), TARGET    :: attr_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      INTEGER(KIND=hsize_t), DIMENSION(1)                :: dims
      TYPE(c_ptr)                                        :: buffer

      dims(1) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! set the type of data
      type_id = h5t_native_double

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataspace')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 attribute')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to write HDF5 attribute')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 attribute')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

   END SUBROUTINE h5awrite_double_simple

! **************************************************************************************************
!> \brief Write an array of integers attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of integers
! **************************************************************************************************
   SUBROUTINE h5awrite_integer_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: attr_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: attr_id, space_id, type_id
      INTEGER(KIND=hsize_t), DIMENSION(1)                :: dims
      TYPE(c_ptr)                                        :: buffer

      dims(1) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! set the type of data
      type_id = h5t_native_integer

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataspace')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 attribute')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to write HDF5 attribute')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 attribute')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

   END SUBROUTINE h5awrite_integer_simple

! **************************************************************************************************
!> \brief Write a dataset containing an array of doubles
!> \param loc_id either file id or group id
!> \param dset_name the name of the dataset
!> \param dset_data the dataset data, i.e. the array of doubles
! **************************************************************************************************
   SUBROUTINE h5dwrite_double_simple(loc_id, dset_name, dset_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: dset_name
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), TARGET    :: dset_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: dset_id, space_id, type_id
      INTEGER(KIND=hsize_t), DIMENSION(1)                :: dims
      TYPE(c_ptr)                                        :: buffer

      dims(1) = SIZE(dset_data, kind=hsize_t)   ! length of array

      ! set the type of data
      type_id = h5t_native_double

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataspace')

      ! create the dataset
      CALL h5dcreate_f(loc_id, dset_name, type_id, space_id, dset_id, error)
      IF (error < 0) CPABORT('ERROR: failed to create HDF5 dataset')

      ! the actual pointer to be passed
      buffer = C_LOC(dset_data(1))

      ! write the string array attribute to file
      CALL h5dwrite_f(dset_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to write HDF5 dataset')

      ! close dataset
      CALL h5dclose_f(dset_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataset')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

   END SUBROUTINE h5dwrite_double_simple

! **************************************************************************************************
!> \brief Read a dataset containing an array of doubles
!> \param loc_id either file id or group id
!> \param dset_name the name of the dataset
!> \param dset_data where the read dataset data will be written
! **************************************************************************************************
   SUBROUTINE h5dread_double_simple(loc_id, dset_name, dset_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: dset_name
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: dset_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: dset_id, npoints, space_id, type_id
      INTEGER(KIND=hsize_t), DIMENSION(1)                :: dims

      dims(1) = SIZE(dset_data, kind=hsize_t)   ! length of array

      ! set the type of data
      type_id = h5t_native_double

      ! open the dataset
      CALL h5dopen_f(loc_id, dset_name, dset_id, error)
      IF (error < 0) CPABORT('ERROR: failed to open HDF5 dataset')

      ! get information on the dataspace
      CALL h5dget_space_f(dset_id, space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to fetch HDF5 dataspace info')

      ! get dataspace dims
      CALL h5sget_simple_extent_npoints_f(space_id, npoints, error)
      IF (error < 0) CPABORT('ERROR: failed to fetch HDF5 dataspace dimension')

      ! read the data
      CALL h5dread_f(dset_id, type_id, dset_data, dims, error)
      IF (error < 0) CPABORT('ERROR: failed to read HDF5 dataset')

      ! close dataset
      CALL h5dclose_f(dset_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataset')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 dataspace')

   END SUBROUTINE h5dread_double_simple

! **************************************************************************************************
!> \brief Read an attribute containing a scalar double
!> \param loc_id either file id or group id
!> \param attr_name ...
!> \param attr_data ...
! **************************************************************************************************
   SUBROUTINE h5aread_double_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      REAL(KIND=dp), INTENT(OUT), TARGET                 :: attr_data

      INTEGER                                            :: error
      INTEGER(KIND=hid_t)                                :: attr_id, type_id
      TYPE(c_ptr)                                        :: buffer

      ! set the type of data
      type_id = h5t_native_double

      ! open the attribute
      CALL h5aopen_f(loc_id, attr_name, attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to open HDF5 attribute')

      buffer = C_LOC(attr_data)
      ! read the data
      CALL h5aread_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: failed to read HDF5 attribute')

      ! close the attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: failed to close HDF5 attribute')

   END SUBROUTINE h5aread_double_scalar

#endif

END MODULE hdf5_wrapper
