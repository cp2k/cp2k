!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Wrapper for various memory allocation methods and bindings
!> \par Purpose
!>      An MPI implementation may require memory allocation throught
!>      its own mpi_alloc_mem subroutine for remote memory
!>      access. Using its memory allocation may also be more efficient
!>      for other, non-RMA routines.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 1.9
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
!> - 2009-06-09 Fixes for the C interface to the C routines.
!> - 2011-05-11 Removes the C binding routine interfaces and the
!>              interface to C MPI.
! *****************************************************************************

MODULE c_mpi_calls
  USE kinds,                           ONLY: dp, sp, real_4, real_8,&
                                             int_4, int_8
  USE ISO_C_BINDING

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mp_alloc_mem, mp_free_mem

#if defined(__parallel)
  INCLUDE "mpif.h"
#endif


  INTERFACE mp_alloc_mem
     MODULE PROCEDURE mp_alloc_mem_i, mp_alloc_mem_l,&
          mp_alloc_mem_d, mp_alloc_mem_z,&
          mp_alloc_mem_s, mp_alloc_mem_c
  END INTERFACE

  INTERFACE mp_free_mem
     MODULE PROCEDURE mp_free_mem_i, mp_free_mem_l,&
          mp_free_mem_d, mp_free_mem_z,&
          mp_free_mem_s, mp_free_mem_c
  END INTERFACE

CONTAINS

#if defined(__parallel)
! *****************************************************************************
!> \brief Allocates an integer array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_i(DATA, len, stat)
     INTEGER(KIND=int_4), DIMENSION(:), POINTER :: DATA
     INTEGER, INTENT(IN)                        :: len
     INTEGER, INTENT(OUT), OPTIONAL             :: stat

     TYPE(C_PTR)              :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
     INTEGER                  :: integer_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_INTEGER, integer_size, ierr)
     mp_size = length * integer_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_alloc_mem_i

! *****************************************************************************
!> \brief Allocates an integer array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_l(DATA, len, stat)
     INTEGER(KIND=int_8), DIMENSION(:), POINTER :: DATA
     INTEGER, INTENT(IN)                        :: len
     INTEGER, INTENT(OUT), OPTIONAL             :: stat

     TYPE(C_PTR)              :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
     INTEGER                  :: integer_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_INTEGER8, integer_size, ierr)
     mp_size = length * integer_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
  END SUBROUTINE mp_alloc_mem_l


! *****************************************************************************
!> \brief Allocates a double real array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_d(DATA, len, stat)
     REAL(KIND=real_8), DIMENSION(:), POINTER  :: DATA
     INTEGER, INTENT(IN)                   :: len
     INTEGER, INTENT(OUT), OPTIONAL        :: stat

     TYPE(C_PTR)                           :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
     INTEGER                               :: double_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, double_size, ierr)
     mp_size = length * double_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_alloc_mem_d

! *****************************************************************************
!> \brief Allocates a double complex array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_z(DATA, len, stat)
     COMPLEX(KIND=real_8), DIMENSION(:), POINTER :: DATA
     INTEGER, INTENT(IN)                     :: len
     INTEGER, INTENT(OUT), OPTIONAL          :: stat
     TYPE(C_PTR)                             :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
     INTEGER                                 :: double_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX, double_size, ierr)
     mp_size = length * double_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
  END SUBROUTINE mp_alloc_mem_z

! *****************************************************************************
!> \brief Allocates a single real array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_s(DATA, len, stat)
     REAL(KIND=real_4), DIMENSION(:), POINTER  :: DATA
     INTEGER, INTENT(IN)                   :: len
     INTEGER, INTENT(OUT), OPTIONAL        :: stat

     TYPE(C_PTR)                           :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
     INTEGER                               :: single_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_REAL, single_size, ierr)
     mp_size = length * single_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_alloc_mem_s

! *****************************************************************************
!> \brief Allocates a double complex array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_alloc_mem_c(DATA, len, stat)
     COMPLEX(KIND=real_4), DIMENSION(:), POINTER :: DATA
     INTEGER, INTENT(IN)                     :: len
     INTEGER, INTENT(OUT), OPTIONAL          :: stat

     TYPE(C_PTR)                             :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
     INTEGER                                 :: single_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_COMPLEX, single_size, ierr)
     mp_size = length * single_size
     mp_info = MPI_INFO_NULL
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_alloc_mem_c


! *****************************************************************************
!> \brief Deallocates an integer array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_i(DATA, stat)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_i


! *****************************************************************************
!> \brief Deallocates an integer array, version with either C bindings or CRAY
!>        pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_l(DATA, stat)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_l


! *****************************************************************************
!> \brief Deallocates a double real array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_d(DATA, stat)
    REAL(KIND=real_8), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_d

! *****************************************************************************
!> \brief Deallocates a double complex array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_z(DATA, stat)
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_z

! *****************************************************************************
!> \brief Deallocates a double real array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_s(DATA, stat)
    REAL(KIND=real_4), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_s

! *****************************************************************************
!> \brief Deallocates a double complex array, version with either C bindings or
!>        CRAY pointers
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_c(DATA, stat)
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    INTEGER                                  :: mp_res

     CALL MPI_FREE_MEM(DATA, mp_res)
     IF (PRESENT (stat)) stat = mp_res
   END SUBROUTINE mp_free_mem_c

#else 

! *****************************************************************************
!> \brief Allocates an integer array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_i(DATA, size, stat)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_i


! *****************************************************************************
!> \brief Allocates an integer array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_l(DATA, size, stat)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_l


! *****************************************************************************
!> \brief Allocates a double real array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_d(DATA, size, stat)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_d

! *****************************************************************************
!> \brief Allocates a double complex array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_z(DATA, size, stat)
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_z

! *****************************************************************************
!> \brief Allocates a double real array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_s(DATA, size, stat)
    REAL(KIND=sp), DIMENSION(:), POINTER     :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_s

! *****************************************************************************
!> \brief Allocates a double complex array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[in] len        length (in data elements) of data array allocation
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_alloc_mem_c(DATA, size, stat)
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: DATA
    INTEGER, INTENT(IN)                      :: size
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        ALLOCATE(DATA(size), stat=stat)
     ELSE
        ALLOCATE(DATA(size))
     ENDIF
   END SUBROUTINE mp_alloc_mem_c

! *****************************************************************************
!> \brief Deallocates an integer array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_free_mem_i(DATA, stat)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_i


! *****************************************************************************
!> \brief Deallocates an integer array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_free_mem_l(DATA, stat)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_l


! *****************************************************************************
!> \brief Deallocates a double real array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_free_mem_d(DATA, stat)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_d

! *****************************************************************************
!> \brief Deallocates a double complex array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_z(DATA, stat)
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_z

! *****************************************************************************
!> \brief Deallocates a single real array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_free_mem_s(DATA, stat)
    REAL(KIND=sp), DIMENSION(:), POINTER     :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_s

! *****************************************************************************
!> \brief Deallocates a single complex array, generic version
!> \author UB
!> \param data           data array to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
   SUBROUTINE mp_free_mem_c(DATA, stat)
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: DATA
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

     IF (PRESENT (stat)) THEN
        DEALLOCATE(DATA, stat=stat)
     ELSE
        DEALLOCATE(DATA)
     ENDIF
   END SUBROUTINE mp_free_mem_c

#endif 

 END MODULE c_mpi_calls
