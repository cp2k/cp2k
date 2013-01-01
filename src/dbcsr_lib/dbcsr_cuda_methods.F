!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-06
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-06
! *****************************************************************************
MODULE dbcsr_cuda_methods
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type,&
                                             dbcsr_cuda_mem_type_c4,&
                                             dbcsr_cuda_mem_type_c8,&
                                             dbcsr_cuda_mem_type_i4,&
                                             dbcsr_cuda_mem_type_i8,&
                                             dbcsr_cuda_mem_type_r4,&
                                             dbcsr_cuda_mem_type_r8
  USE dbcsr_error_handling
  USE dbcsr_types,                     ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_int_4,&
                                             dbcsr_type_int_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_methods'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.


  PUBLIC :: dbcsr_cuda_dev_mem_setup

  PUBLIC :: dbcsr_cuda_dev_mem_get_type, dbcsr_cuda_dev_mem_get_alloc,&
       dbcsr_cuda_dev_mem_get_size

  PUBLIC :: dbcsr_cuda_dev_mem_get_r, dbcsr_cuda_dev_mem_get_d,&
            dbcsr_cuda_dev_mem_get_c, dbcsr_cuda_dev_mem_get_z

  INTERFACE dbcsr_cuda_dev_mem_get_alloc
     MODULE PROCEDURE dev_mem_get_alloc_any
     MODULE PROCEDURE dev_mem_get_alloc_i, dev_mem_get_alloc_l
     MODULE PROCEDURE dev_mem_get_alloc_r, dev_mem_get_alloc_d
     MODULE PROCEDURE dev_mem_get_alloc_c, dev_mem_get_alloc_z
  END INTERFACE


  INTERFACE dbcsr_cuda_dev_mem_get_size
     MODULE PROCEDURE dev_mem_get_size_any
     MODULE PROCEDURE dev_mem_get_size_i, dev_mem_get_size_l
     MODULE PROCEDURE dev_mem_get_size_r, dev_mem_get_size_d
     MODULE PROCEDURE dev_mem_get_size_c, dev_mem_get_size_z
  END INTERFACE

CONTAINS

  SUBROUTINE dbcsr_cuda_dev_mem_setup (dev_mem, mem_type)
    TYPE(dbcsr_cuda_mem_type), INTENT(OUT)   :: dev_mem
    INTEGER, INTENT(IN)                      :: mem_type

    dev_mem%data_type = mem_type
  END SUBROUTINE dbcsr_cuda_dev_mem_setup

  PURE FUNCTION dbcsr_cuda_dev_mem_get_type (dev_mem) RESULT (data_type)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    INTEGER                                  :: data_type

    data_type = dev_mem%data_type
  END FUNCTION dbcsr_cuda_dev_mem_get_type


  PURE FUNCTION dev_mem_get_alloc_any (dev_mem) RESULT (data_size)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    INTEGER                                  :: data_size

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_i)
    CASE (dbcsr_type_int_8)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_l)
    CASE (dbcsr_type_real_4)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_r)
    CASE (dbcsr_type_real_8)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_d)
    CASE (dbcsr_type_complex_4)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_c)
    CASE (dbcsr_type_complex_8)
       data_size = dbcsr_cuda_dev_mem_get_alloc(dev_mem%d_z)
    END SELECT
  END FUNCTION dev_mem_get_alloc_any

  PURE FUNCTION dev_mem_get_size_any (dev_mem) RESULT (data_size)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    INTEGER                                  :: data_size

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_i)
    CASE (dbcsr_type_int_8)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_l)
    CASE (dbcsr_type_real_4)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_r)
    CASE (dbcsr_type_real_8)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_d)
    CASE (dbcsr_type_complex_4)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_c)
    CASE (dbcsr_type_complex_8)
       data_size = dbcsr_cuda_dev_mem_get_size(dev_mem%d_z)
    END SELECT
  END FUNCTION dev_mem_get_size_any


  PURE FUNCTION dbcsr_cuda_dev_mem_get_r (dev_mem) RESULT (type_r)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    TYPE(dbcsr_cuda_mem_type_r4)             :: type_r

    type_r = dev_mem%d_r
  END FUNCTION dbcsr_cuda_dev_mem_get_r
  PURE FUNCTION dbcsr_cuda_dev_mem_get_d (dev_mem) RESULT (type_d)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    TYPE(dbcsr_cuda_mem_type_r8)             :: type_d

    type_d = dev_mem%d_d
  END FUNCTION dbcsr_cuda_dev_mem_get_d
  PURE FUNCTION dbcsr_cuda_dev_mem_get_c (dev_mem) RESULT (type_c)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    TYPE(dbcsr_cuda_mem_type_c4)             :: type_c

    type_c = dev_mem%d_c
  END FUNCTION dbcsr_cuda_dev_mem_get_c
  PURE FUNCTION dbcsr_cuda_dev_mem_get_z (dev_mem) RESULT (type_z)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    TYPE(dbcsr_cuda_mem_type_c8)             :: type_z

    type_z = dev_mem%d_z
  END FUNCTION dbcsr_cuda_dev_mem_get_z


#include "dbcsr_cuda_methods_i.F"
#include "dbcsr_cuda_methods_l.F"
#include "dbcsr_cuda_methods_r.F"
#include "dbcsr_cuda_methods_d.F"
#include "dbcsr_cuda_methods_c.F"
#include "dbcsr_cuda_methods_z.F"

END MODULE dbcsr_cuda_methods
