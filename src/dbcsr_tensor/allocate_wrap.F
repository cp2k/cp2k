!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper to allocate and optionally assign array from source array.
!>        Will be redundant once compilers implement Fortran 2008 standard :-/
!> \author Patrick Seewald
! **************************************************************************************************
MODULE allocate_wrap

#:include "dbcsr_tensor.fypp"
#:set maxdim = fortran_max_ndim

   USE kinds, ONLY: ${uselist(dtype_float_prec)}$

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: allocate_any

   INTERFACE allocate_any
#:for dparam, dtype, dsuffix in dtype_all_list
#:for dim in range(1, maxdim+1)
      MODULE PROCEDURE allocate_${dim}$d_${dsuffix}$
#:endfor
#:endfor
   END INTERFACE

CONTAINS

#:for dparam, dtype, dsuffix in dtype_all_list
#:for dim in range(1, maxdim+1)
   PURE SUBROUTINE allocate_${dim}$d_${dsuffix}$(sizes, array, source)
      INTEGER, DIMENSION(${dim}$), INTENT(IN)                              :: sizes
      ${dtype}$, DIMENSION(${shape_colon(dim)}$), ALLOCATABLE, INTENT(OUT) :: array
      ${dtype}$, DIMENSION(${shape_colon(dim)}$), INTENT(IN), OPTIONAL     :: source
      ALLOCATE(array(${arrlist("sizes", nmax=dim)}$))
      IF(PRESENT(source)) THEN
         array(${shape_colon(dim)}$) = source(${shape_colon(dim)}$)
      ENDIF
   END SUBROUTINE
#:endfor
#:endfor
END MODULE
