!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for array operations: allocation, assignment and reshape
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
! **************************************************************************************************
!> \brief Magic routine for array operations.
!> \param sizes sizes (shape) of array to be allocated
!> \param array target array to be allocated
!> \param source source array to be copied to target array
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic)
! **************************************************************************************************
   PURE SUBROUTINE allocate_${dim}$d_${dsuffix}$(sizes, array, source, order)
      INTEGER, DIMENSION(${dim}$), INTENT(IN)                              :: sizes
      ${dtype}$, DIMENSION(${shape_colon(dim)}$), ALLOCATABLE, INTENT(OUT) :: array
      ${dtype}$, DIMENSION(${shape_colon(dim)}$), INTENT(IN), OPTIONAL     :: source
      INTEGER, DIMENSION(${dim}$), INTENT(IN), OPTIONAL                    :: order
      INTEGER, DIMENSION(${dim}$)                                          :: sizes_reordered
      IF (PRESENT(order)) THEN
         sizes_reordered(order(:)) = sizes(:)
      ELSE
         sizes_reordered(:) = sizes(:)
      ENDIF
      ALLOCATE(array(${arrlist("sizes_reordered", nmax=dim)}$))
      IF(PRESENT(source)) THEN
         IF (PRESENT(order)) THEN
            array(${shape_colon(dim)}$) = RESHAPE(source, SHAPE(array), order=order)
         ELSE
            array(${shape_colon(dim)}$) = source(${shape_colon(dim)}$)
         ENDIF
      ENDIF
   END SUBROUTINE
#:endfor
#:endfor
END MODULE
