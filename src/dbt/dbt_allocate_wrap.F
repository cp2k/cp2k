!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for allocating, copying and reshaping arrays.
!> \todo with fortran 2008 support, this should be replaced by plain ALLOCATE
!> \note in particular ALLOCATE(..., SOURCE=...) does not work in gcc 5.4.0, see also
!>       https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_allocate_wrap
   #:include "dbt_macros.fypp"
   #:set maxdim = fortran_max_ndim

   USE kinds, ONLY: dp

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: allocate_any

   INTERFACE allocate_any
      #:for dim in range(1, maxdim+1)
         MODULE PROCEDURE allocate_${dim}$d
      #:endfor
   END INTERFACE

CONTAINS

   #:for dim in range(1, maxdim+1)
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_${dim}$d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(${shape_colon(dim)}$), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(${dim}$), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(${shape_colon(dim)}$), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(${dim}$), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(${dim}$)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(${arrlist("shape_prv", nmax=dim)}$))
               array(${shape_colon(dim)}$) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(${arrlist("shape_prv", nmax=dim)}$), source=source)
            END IF
         ELSE
            ALLOCATE (array(${arrlist("shape_prv", nmax=dim)}$))
         END IF

      END SUBROUTINE
   #:endfor
END MODULE
