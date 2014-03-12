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
!> \par Choices
!>      There are three possible memory (de)allocation calling scenarios
!>      selected by compile-time options.
!>      <ul>
!> <li> (No MPI Memory) Memory management is not handled by MPI but
!>           regular ALLOCATE and DEALLOCATE statements.
!>           (__no_mpi_allocation)
!> <li> (C bindings) Memory management is handled by calling
!>           MPI_ALLOC_MEM using a TYPE(C_PTR) type as the formal
!>           argument to the BASEPTR dummy argument.  Deallocation is
!>           also done by calling MPI_FREE_MEM. (__c_bindings)
!> <li> (CRAY pointers) Memory management is handled by calling
!>           MPI_ALLOC_MEM using a CRAY pointer as the formal argument
!>           to the BASEPTR dummy argument.  Deallocation is also done
!>           by calling MPI_FREE_MEM. (__c_bindings)
!>      </ul>
!> \par Precedence of options
!> <ol>
!> <li> If __NO_MPI_MEMORY is defined, then the No MPI Memory mode
!>      is used; otherwise,
!> <li> The C binding interface is used unless either
!>      __HAS_NO_ISO_C_BINDING or __cray_pointers is defined;
!>      otherwise,
!> <lI> CRAY pointers are used unless
!>       __NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION is set.; otherwise the
!>       compilation aborts with an error.
!> </ol>
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

MODULE dbcsr_c_mpi_calls
  USE dbcsr_kinds,                           ONLY: dp, sp, real_4, real_8,&
                                             int_4, int_8

! Reset the settings for code selection
#ifdef __mpi_allocation
#undef __mpi_allocation
#endif

#ifdef __c_bindings
#undef __c_bindings
#endif

! Logic for code selection

#if defined(__NO_MPI_MEMORY) || !defined(__parallel)
! Do nothing if not in parallel or explicitly not wanting MPI memory.
#elif !defined(__HAS_NO_ISO_C_BINDING) && !defined(__cray_pointers)
! Use C binding features.
#define __mpi_allocation
#define __c_bindings
#elif defined (__CRAY_POINTERS) || defined (__cray_pointers)
! Use CRAY pointers
#define __mpi_allocation
#ifndef __cray_pointers
#define __cray_pointers
#endif
#if defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
  ERROR "Unable to do MPI memory allocation.  Use -D__NO_MPI_MEMORY"
#endif

#else
  ! Technically MPI memory is not needed anymore.
  !ERROR "Unable to do MPI memory allocation.  Use -D__NO_MPI_MEMORY"
#endif


! Real code starts here

#ifndef __HAS_NO_ISO_C_BINDING
  USE ISO_C_BINDING
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mp_alloc_mem, mp_free_mem

#if defined(__mpi_allocation)
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

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, integer version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] iarray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_i(iarray, n, fp)
    INTEGER(KIND=int_4), DIMENSION(*), &
      INTENT(IN), TARGET                     :: iarray
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: fp

    fp => iarray(1:n)
  END SUBROUTINE array_to_fp_i

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, integer version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] iarray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_l(iarray, n, fp)
    INTEGER(KIND=int_8), DIMENSION(*), &
      INTENT(IN), TARGET                     :: iarray
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: fp

    fp => iarray(1:n)
  END SUBROUTINE array_to_fp_l

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, double real version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_d(darray, n, fp)
    REAL(KIND=real_8), DIMENSION(*), &
      INTENT(IN), TARGET                     :: darray
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=real_8), DIMENSION(:), POINTER :: fp

    fp => darray(1:n)
  END SUBROUTINE array_to_fp_d

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, double complex
!>        version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_z(darray, n, fp)
    COMPLEX(KIND=real_8), DIMENSION(*), &
      INTENT(IN), TARGET                     :: darray
    INTEGER, INTENT(IN)                      :: n
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: fp

    fp => darray(1:n)
  END SUBROUTINE array_to_fp_z

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, single real version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_r(darray, n, fp)
    REAL(KIND=real_4), DIMENSION(*), &
      INTENT(IN), TARGET                     :: darray
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=real_4), DIMENSION(:), POINTER :: fp

    fp => darray(1:n)
  END SUBROUTINE array_to_fp_r

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, single complex
!>        version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param[in] n          elemements in data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_c(darray, n, fp)
    COMPLEX(KIND=real_4), DIMENSION(*), &
      INTENT(IN), TARGET                     :: darray
    INTEGER, INTENT(IN)                      :: n
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: fp

    fp => darray(1:n)
  END SUBROUTINE array_to_fp_c


#if defined(__c_bindings) || defined(__cray_pointers)
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

#ifdef __c_bindings
     TYPE(C_PTR)              :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
#else /* (__c_bindings) */
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
     INTEGER(KIND=int_4)                :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif /* (__c_bindings) */
     INTEGER                  :: integer_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_INTEGER, integer_size, ierr)
     mp_size = length * integer_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_i(dynmem(1), len, DATA)
#endif
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

#ifdef __c_bindings
     TYPE(C_PTR)              :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
#else /* (__c_bindings) */
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
     INTEGER(KIND=int_8)                :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif /* (__c_bindings) */
     INTEGER                  :: integer_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_INTEGER8, integer_size, ierr)
     mp_size = length * integer_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_l(dynmem(1), len, DATA)
#endif
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

#ifdef __c_bindings
     TYPE(C_PTR)                           :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
     REAL(KIND=real_8)                     :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif
     INTEGER                               :: double_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, double_size, ierr)
     mp_size = length * double_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_d(dynmem(1), len, DATA)
#endif
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

#ifdef __c_bindings
     TYPE(C_PTR)                             :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
     COMPLEX(KIND=real_8)                    :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif
     INTEGER                                 :: double_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX, double_size, ierr)
     mp_size = length * double_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_z(dynmem(1), len, DATA)
#endif
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

#ifdef __c_bindings
     TYPE(C_PTR)                           :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
     REAL(KIND=real_4)                         :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif
     INTEGER                               :: single_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_REAL, single_size, ierr)
     mp_size = length * single_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_r(dynmem(1), len, DATA)
#endif
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

#ifdef __c_bindings
     TYPE(C_PTR)                             :: mp_baseptr
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
     COMPLEX(KIND=real_4)                        :: dynmem(1:*)
     POINTER (icp, dynmem)
#endif
     INTEGER                                 :: single_size, length, ierr

     length = MAX(len,1)
     CALL MPI_TYPE_SIZE(MPI_COMPLEX, single_size, ierr)
     mp_size = length * single_size
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
     CALL MPI_ALLOC_MEM(mp_size, mp_info, mp_baseptr, mp_res)
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_c(dynmem(1), len, DATA)
#endif
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
     !write(6,*)'int free', ASSOCIATED(data), 'res=', mp_res
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


! -----------------------------------------------------------------------------
! Up to now we had the case with either cray pointers or the C interface to
! C or fortran function calls.
!
! What follows now is when the above facilities are unavailable and we just use
! Fortran's built-in memory allocation and deallocation.
! -----------------------------------------------------------------------------

#else /* __c_bindings || __cray_pointers */

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

#endif /* __c_bindings */

END MODULE dbcsr_c_mpi_calls
