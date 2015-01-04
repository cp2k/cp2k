!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Types related to DBCSR data area
!> \author  Ole Schuett
!> \date    10-2014
! *****************************************************************************
MODULE dbcsr_data_types
  USE acc_devmem,                      ONLY: acc_devmem_type
  USE acc_event,                       ONLY: acc_event_type
  USE acc_stream,                      ONLY: acc_stream_type
  USE kinds,                           ONLY: &
       dp, int_4, int_4_size, int_8, int_8_size, real_4, real_4_size, real_8, &
       real_8_size

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_types'

  PUBLIC :: dbcsr_data_obj, dbcsr_data_area_type, dbcsr_scalar_type
  PUBLIC :: dbcsr_datatype_sizeof
  PUBLIC :: dbcsr_mempool_type,&
            dbcsr_mempool_entry_type,&
            dbcsr_memtype_type,&
            dbcsr_memtype_default
  PUBLIC :: dbcsr_type_real_4, dbcsr_type_real_8,&
            dbcsr_type_complex_4, dbcsr_type_complex_8,&
            dbcsr_type_real_default, dbcsr_type_complex_default,&
            dbcsr_type_real_4_2d, dbcsr_type_real_8_2d,&
            dbcsr_type_complex_4_2d, dbcsr_type_complex_8_2d,&
            dbcsr_type_int_4, dbcsr_type_int_8


! *****************************************************************************
!> \brief Stores a scalar in any of the supported data types.
!> \par Reasoning
!>      Easier routine interfaces
!> \var r_sp        stores real values in single precision
!> \var r_dp        stores real values in double precision
!> \var c_sp        stores complex values in single precision
!> \var c_dp        stores complex values in double precision
!> \var refcount    reference counter for current structure
!> \var data_type   which of the data types is actually used
! *****************************************************************************
  TYPE dbcsr_scalar_type
     REAL(KIND=real_4)    :: r_sp
     REAL(KIND=real_8)    :: r_dp
     COMPLEX(KIND=real_4) :: c_sp
     COMPLEX(KIND=real_8) :: c_dp
     INTEGER              :: data_type
  END TYPE dbcsr_scalar_type

 TYPE dbcsr_data_obj
     TYPE(dbcsr_data_area_type), POINTER           :: d => Null()
  END TYPE dbcsr_data_obj

! *****************************************************************************
!> \brief Memory related types
! *****************************************************************************
  TYPE dbcsr_mempool_type
    TYPE(dbcsr_mempool_entry_type), POINTER :: root => Null()
    INTEGER                                 :: capacity = 1
    !$ INTEGER(KIND=omp_lock_kind)          :: lock
  END TYPE dbcsr_mempool_type

  TYPE dbcsr_mempool_entry_type
      TYPE(dbcsr_data_obj)                    :: area
      TYPE(dbcsr_mempool_entry_type), POINTER :: next => Null()
  END TYPE dbcsr_mempool_entry_type

  TYPE dbcsr_memtype_type
     LOGICAL                           :: mpi = .FALSE.
     LOGICAL                           :: acc_hostalloc = .FALSE.
     LOGICAL                           :: acc_devalloc = .FALSE.
     TYPE(acc_stream_type)             :: acc_stream = acc_stream_type()
     TYPE(dbcsr_mempool_type), POINTER :: pool => Null()
     REAL(KIND=dp)                     :: oversize_factor = 1.0
  END TYPE dbcsr_memtype_type

  !providing pool=Null() explicitly to circumvent bug in ifort 12.1
  TYPE(dbcsr_memtype_type), PARAMETER :: dbcsr_memtype_default = dbcsr_memtype_type(pool=Null())

! *****************************************************************************
!> \brief Stores actual data
!> \var r_sp        stores real values in single precision
!> \var r_dp        stores real values in double precision
!> \var c_sp        stores complex values in single precision
!> \var c_dp        stores complex values in double precision
!> \var ref_size    last data element in array actually referenced
!> \var refcount    reference counter for current structure
!> \var memory_type   type of memory where data lives
!> \var data_type   which of the data types is actually used
! *****************************************************************************
  TYPE dbcsr_data_area_type
     INTEGER(KIND=int_4), DIMENSION(:), POINTER    :: i4    => Null()
     INTEGER(KIND=int_8), DIMENSION(:), POINTER    :: i8    => Null()
     REAL(KIND=real_4), DIMENSION(:), POINTER      :: r_sp  => Null()
     REAL(KIND=real_8), DIMENSION(:), POINTER      :: r_dp  => Null()
     COMPLEX(KIND=real_4), DIMENSION(:), POINTER   :: c_sp  => Null()
     COMPLEX(KIND=real_8), DIMENSION(:), POINTER   :: c_dp  => Null()
     REAL(KIND=real_4), DIMENSION(:,:), POINTER    :: r2_sp => Null()
     REAL(KIND=real_8), DIMENSION(:,:), POINTER    :: r2_dp => Null()
     COMPLEX(KIND=real_4), DIMENSION(:,:), POINTER :: c2_sp => Null()
     COMPLEX(KIND=real_8), DIMENSION(:,:), POINTER :: c2_dp => Null()
     INTEGER                                  :: ref_size = -1
     INTEGER                                  :: refcount = -1
     INTEGER                                  :: id = -1
     TYPE(dbcsr_memtype_type)                 :: memory_type = dbcsr_memtype_default
     INTEGER                                  :: data_type  = -1
     TYPE(acc_devmem_type)                    :: acc_devmem
     TYPE(acc_event_type)                     :: acc_ready
  END TYPE dbcsr_data_area_type

!> Type definitions:
!> * bit 0: always 1
!> * bit 1: single (0: 4) vs. double (1: 1)
!> * bit 2: real (0) vs. complex (1)
!> * bit 3: dimension (0: 1, 1: 2)
!> * bit 4: floating point (0) or integer type (1)
  INTEGER, PARAMETER          :: dbcsr_type_real_4 = 1         !001
  INTEGER, PARAMETER          :: dbcsr_type_real_8 = 3         !011
  INTEGER, PARAMETER          :: dbcsr_type_complex_4 = 5      !101
  INTEGER, PARAMETER          :: dbcsr_type_complex_8 = 7      !111

  INTEGER, PARAMETER          :: dbcsr_type_real_4_2d = 9      !1001
  INTEGER, PARAMETER          :: dbcsr_type_real_8_2d = 11     !1011
  INTEGER, PARAMETER          :: dbcsr_type_complex_4_2d = 13  !1101
  INTEGER, PARAMETER          :: dbcsr_type_complex_8_2d = 15  !1111

  INTEGER, PARAMETER          :: dbcsr_type_int_4 = 17         !10001
  INTEGER, PARAMETER          :: dbcsr_type_int_8 = 19         !10011

  INTEGER, PARAMETER          :: dbcsr_type_real_default = dbcsr_type_real_8
  INTEGER, PARAMETER          :: dbcsr_type_complex_default = dbcsr_type_complex_8

CONTAINS

! *****************************************************************************
!> \brief Helper-routine, returns size of given datatype in terms of bytes.
!> \param datatype ...
!> \retval size ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_datatype_sizeof(datatype) RESULT(size)
    INTEGER, INTENT(IN)                      :: datatype
    INTEGER                                  :: size

    SELECT CASE (datatype)
    CASE (dbcsr_type_int_4)
       size = int_4_size
    CASE (dbcsr_type_int_8)
       size = int_8_size
    CASE (dbcsr_type_real_4)
       size = real_4_size
    CASE (dbcsr_type_real_8)
       size = real_8_size
    CASE (dbcsr_type_complex_4)
       size = (2*real_4_size)
    CASE (dbcsr_type_complex_8)
       size = (2*real_8_size)
    CASE default
       STOP "dbcsr_datatype_sizeof: Invalid data type"
    END SELECT
 END FUNCTION dbcsr_datatype_sizeof

END MODULE dbcsr_data_types
