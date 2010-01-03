!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Wrapper for various memory allocation methods and bindings
!> \par Purpose
!>       An MPI implementation may require memory allocation throught its own
!>      interface for remote memory access. Using its memory allocation may
!>      also be more efficient for other, non-RMA routines.
!> \par Choices
!>      There are four possible memory (de)allocation calling scenarios
!>      selected by compile-time options.
!>      <ol>
!> <li> By default memory management is not handled by MPI but regular
!>           ALLOCATE and DEALLOCATE statements.
!> <li> If __cray_pointers is declared, then the MPI MPI_ALLOC_MEM and
!>            MPI_FREE_MEM routine is called through its default
!>            binding as present in mpif.h, but the fortran compiler
!>            must support the Cray pointers extension.
!> <li> If __c_bindings and __mpi_f_bindings are declared, then the
!>           allocation is performed by a call to the assumed binding
!>           MPI_ALLOC_MEM_; this name can be overidden with the
!>           __MPI_F_ALLOC_DECL declaration. Deallocation is through
!>           the default binding MPI_FREE_MEM as present in
!>           mpif.h. This option is a type-safe alternative to
!>           CRAY-pointers.
!> <li> If __c_bindings and __mpi_c_bindings are declared, then the
!>           allocation is performed by a direct call to the
!>           mpi_alloc_mem C function. While in this case the function
!>           name should exist in a correctly-linked program, the
!>           parameter types are not guaranteed to match the ones
!>           declared in the C header.
!>      </ol>
!> \par Multiple declarations
!>      In case of multiple declarations, Cray pointers are preferred
!>      to the fortran or C bindings using the C binding interface. In
!>      case the C binding interface is used, the fortran binding is
!>      preferred due to its type safety.
!> \par Suggested choice
!>      The suggested choice is either Cray pointers or __c_bindings
!>      and __mpi_f_bindings. The routines themselves are meant to be
!>      called from the mp_allocate and mp_deallocate routines
!>      declared in message_passing.F
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 1.1
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
!> - 2009-06-09 Fixes for the C interface to the C routines.
! *****************************************************************************

MODULE c_mpi_calls
  USE kinds,                           ONLY: dp, sp, real_4, real_8,&
                                             int_4, int_8
! Do not use C bindings if we can use CRAY pointers.
#ifdef __cray_pointers
#undef __c_bindings
#endif

! Do not use fortran or C bindings if we can not use C bindings at all.
#ifndef __c_bindings
#undef __mpi_f_bindings
#undef __mpi_c_bindings
#endif

! Can only use fortran or C bindings, not both. fortran takes precedence
#if defined(__c_bindings) && !defined(__mpi_f_bindings) && !defined(__mpi_c_bindings)
#define __mpi_f_bindings
#endif

! Must use either fortran or C bindings, not both. fortran takes precedence
#if defined(__mpi_f_bindings) && defined(__mpi_c_bindings)
#undef __mpi_c_bindings
#endif

! The fortran binding of the MPI_Alloc_mem subroutine
#ifndef __MPI_F_ALLOC_DECL
#define __MPI_F_ALLOC_DECL "mpi_alloc_mem_"
#endif

#ifdef __c_bindings
  USE ISO_C_BINDING
#endif

  IMPLICIT NONE


  PUBLIC :: mp_alloc_mem, mp_free_mem

  PRIVATE 

#if defined(__c_bindings) || defined(__cray_pointers)
  INCLUDE "mpif.h"
#endif

#ifdef __c_bindings

#ifndef __mpi_f_bindings
  INTERFACE
     !> We assume that the type for the argument size corresponds to intptr_t.
     !> It is MPI_ADDRESS_TYPE in fortran.
     FUNCTION mpi_alloc_mem_c(size, info, baseptr)&
          RESULT(res) BIND(C, name="MPI_Alloc_mem")
         USE ISO_C_BINDING
    INTEGER(KIND=C_INTPTR_T), INTENT(IN), &
      VALUE                                  :: size
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: info
    TYPE(C_PTR), INTENT(INOUT)               :: baseptr
    INTEGER(KIND=C_INT)                      :: res

     END FUNCTION mpi_alloc_mem_c
  END INTERFACE

#else /* !__mpi_f_bindings */

  INTERFACE
     !> We guess the name of the fortran binding name.
     SUBROUTINE mpi_alloc_mem_f(size, info, baseptr, ierr)&
          BIND(C, name=__MPI_F_ALLOC_DECL)
        USE ISO_C_BINDING
    INTEGER(KIND=MPI_ADDRESS_KIND), &
      INTENT(IN)                             :: size
    INTEGER, INTENT(IN)                      :: info
    TYPE(C_PTR), INTENT(INOUT)               :: baseptr
    INTEGER, INTENT(INOUT)                   :: ierr

     END SUBROUTINE mpi_alloc_mem_f
  END INTERFACE

#endif /* !__mpi_f_bindings */

#endif /* __c_bindings */


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
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_i(iarray, fp)
    INTEGER(KIND=int_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: iarray
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: fp

    fp => iarray
  END SUBROUTINE array_to_fp_i

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, integer version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] iarray     data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_l(iarray, fp)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: iarray
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: fp

    fp => iarray
  END SUBROUTINE array_to_fp_l

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, double real version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_d(darray, fp)
    REAL(KIND=real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: darray
    REAL(KIND=real_8), DIMENSION(:), POINTER :: fp

    fp => darray
  END SUBROUTINE array_to_fp_d

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, double complex
!>        version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_z(darray, fp)
    COMPLEX(KIND=real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: darray
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: fp

    fp => darray
  END SUBROUTINE array_to_fp_z

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, single real version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_r(darray, fp)
    REAL(KIND=real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: darray
    REAL(KIND=real_4), DIMENSION(:), POINTER :: fp

    fp => darray
  END SUBROUTINE array_to_fp_r

! *****************************************************************************
!> \brief Points a FORTRAN 90 pointer to an existing array, single complex
!>        version
!> \par It is intended that the input array is aliased to a Cray pointer.
!> \author UB
!> \param[in] darray     data array
!> \param fp             pointer to assign
! *****************************************************************************
  SUBROUTINE array_to_fp_c(darray, fp)
    COMPLEX(KIND=real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: darray
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: fp

    fp => darray
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T) :: mp_size
     INTEGER(KIND=C_INT)      :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_i(dynmem(1:len), DATA)
#endif 
     !write(6,*)'int allocation', ASSOCIATED(data), (/length/), &
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'int allocation size', SIZE(data)
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T) :: mp_size
     INTEGER(KIND=C_INT)      :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)     :: mp_size
     INTEGER                            :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_l(dynmem(1:len), DATA)
#endif 
     !write(6,*)'int allocation', ASSOCIATED(data), (/length/), &
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'int allocation size', SIZE(data)
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T)              :: mp_size
     INTEGER(KIND=C_INT)                   :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_d(dynmem(1:len), DATA)
#endif
     !write(6,*)'double allocation', ASSOCIATED(data), (/length/),&
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'double allocation size', SIZE(data)
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T)                :: mp_size
     INTEGER(KIND=C_INT)                     :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_z(dynmem(1:len), DATA)
#endif
     !write(6,*)'double allocation', ASSOCIATED(data), (/length/),&
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'double allocation size', SIZE(data)
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T)              :: mp_size
     INTEGER(KIND=C_INT)                   :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)        :: mp_size
     INTEGER                               :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_r(dynmem(1:len), DATA)
#endif
     !write(6,*)'single allocation', ASSOCIATED(data), (/length/),&
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'single allocation size', SIZE(data)
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
#ifndef __mpi_f_bindings
     INTEGER(KIND=C_INTPTR_T)                :: mp_size
     INTEGER(KIND=C_INT)                     :: mp_info, mp_res
#else
     INTEGER(KIND=MPI_ADDRESS_KIND)          :: mp_size
     INTEGER                                 :: mp_info, mp_res
#endif
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
!> The C bindings may use a different value for MPI_INFO_NULL
     mp_info = MPI_INFO_NULL
#ifdef __c_bindings
#ifndef __mpi_f_bindings
     mp_res = mpi_alloc_mem_c(mp_size, mp_info, mp_baseptr)
#else
     CALL mpi_alloc_mem_f(mp_size, mp_info, mp_baseptr, mp_res)
#endif
     CALL C_F_POINTER(mp_baseptr, DATA, (/length/))
#else /* __c_bindings */
     CALL MPI_ALLOC_MEM(mp_size, mp_info, icp, mp_res)
     CALL array_to_fp_c(dynmem(1:len), DATA)
#endif
     !write(6,*)'single allocation', ASSOCIATED(data), (/length/),&
     !     mp_res
     !IF (ASSOCIATED (data)) &
     !     write(6,*)'single allocation size', SIZE(data)
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
     !write(6,*)'int free', ASSOCIATED(data), 'res=', mp_res
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
     !write(6,*)'int free', ASSOCIATED(data), 'res=', mp_res
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
     !write(6,*)'int free', ASSOCIATED(data), 'res=', mp_res
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

END MODULE c_mpi_calls

