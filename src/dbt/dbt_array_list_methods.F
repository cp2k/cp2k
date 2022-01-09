!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Representation of arbitrary number of 1d integer arrays with arbitrary sizes.
!>        This is needed for generic handling of dimension-specific tensor quantities
!>        (such as block index).
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_array_list_methods

   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbt_index, ONLY: dbt_inverse_order
   USE dbt_allocate_wrap, ONLY: allocate_any

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_array_list_methods'

   PUBLIC  :: &
      array_eq_i, &
      array_list, &
      array_offsets, &
      array_sublist, &
      create_array_list, &
      destroy_array_list, &
      get_array_elements, &
      get_arrays, &
      get_ith_array, &
      number_of_arrays, &
      reorder_arrays, &
      sizes_of_arrays, &
      sum_of_arrays, &
      check_equal

   TYPE array_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: col_data
      INTEGER, DIMENSION(:), ALLOCATABLE :: ptr
   END TYPE

   INTERFACE get_ith_array
      MODULE PROCEDURE allocate_and_get_ith_array
      MODULE PROCEDURE get_ith_array
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief number of arrays stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION number_of_arrays(list)
      TYPE(array_list), INTENT(IN) :: list
      INTEGER                      :: number_of_arrays

      number_of_arrays = SIZE(list%ptr) - 1

   END FUNCTION number_of_arrays

! **************************************************************************************************
!> \brief Get an element for each array.
!> \param indices element index for each array
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_array_elements(list, indices)
      TYPE(array_list), INTENT(IN)                           :: list
      INTEGER, DIMENSION(number_of_arrays(list)), INTENT(IN) :: indices
      INTEGER, DIMENSION(number_of_arrays(list))             :: get_array_elements

      INTEGER                                                :: i, ind

      DO i = 1, SIZE(indices)
         ind = indices(i) + list%ptr(i) - 1
         get_array_elements(i) = list%col_data(ind)
      END DO

   END FUNCTION get_array_elements

! **************************************************************************************************
!> \brief collects any number of arrays of different sizes into a single array (list%col_data),
!>        storing the indices that start a new array (list%ptr).
!> \param list list of arrays
!> \param ndata number of arrays
!> \param data arrays 1 and 2
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE create_array_list(list, ndata, ${varlist("data")}$)
      TYPE(array_list), INTENT(OUT)               :: list
      INTEGER, INTENT(IN)                         :: ndata
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("data")}$
      INTEGER                                     :: ptr, size_all

      size_all = 0

      #:for dim in range(1, maxdim+1)
         IF (ndata .GE. ${dim}$) THEN
            CPASSERT(PRESENT(data_${dim}$))
            size_all = size_all + SIZE(data_${dim}$)
         END IF
      #:endfor

      ALLOCATE (list%ptr(ndata + 1))
      ALLOCATE (list%col_data(size_all))

      ptr = 1
      list%ptr(1) = ptr

      #:for dim in range(1, maxdim+1)
         IF (ndata .GE. ${dim}$) THEN
            list%col_data(ptr:ptr + SIZE(data_${dim}$) - 1) = data_${dim}$ (:)
            ptr = ptr + SIZE(data_${dim}$)
            list%ptr(${dim+1}$) = ptr
         END IF
      #:endfor

   END SUBROUTINE

! **************************************************************************************************
!> \brief extract a subset of arrays
!> \param list list of arrays
!> \param i_selected array numbers to retrieve
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION array_sublist(list, i_selected)
      TYPE(array_list), INTENT(IN)                           :: list
      INTEGER, DIMENSION(:), INTENT(IN)                      :: i_selected
      TYPE(array_list)                                       :: array_sublist
      INTEGER :: ndata
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ${varlist("data")}$

      ndata = SIZE(i_selected)

      #:for dim in range(1, maxdim+1)
         IF (ndata == ${dim}$) THEN
            CALL get_arrays(list, ${varlist("data", nmax=dim)}$, i_selected=i_selected)
            CALL create_array_list(array_sublist, ndata, ${varlist("data", nmax=dim)}$)
         END IF
      #:endfor
   END FUNCTION

! **************************************************************************************************
!> \brief destroy array list.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE destroy_array_list(list)
      TYPE(array_list), INTENT(INOUT) :: list

      DEALLOCATE (list%ptr, list%col_data)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Get all arrays contained in list
!> \param data arrays 1 and 2
!> \param i_selected array numbers to retrieve (if not present, all arrays are returned)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE get_arrays(list, ${varlist("data")}$, i_selected)
      !! Get all arrays contained in list
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: ${varlist("data")}$
      INTEGER, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: i_selected
      INTEGER                                            :: i, ndata
      INTEGER, DIMENSION(number_of_arrays(list))         :: o

      o(:) = 0
      IF (PRESENT(i_selected)) THEN
         ndata = SIZE(i_selected)
         o(1:ndata) = i_selected(:)
      ELSE
         ndata = number_of_arrays(list)
         o(1:ndata) = (/(i, i=1, ndata)/)
      END IF

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
         #:for dim in range(1, maxdim+1)
            IF (ndata > ${dim-1}$) THEN
               ALLOCATE (data_${dim}$, source=col_data(ptr(o(${dim}$)):ptr(o(${dim}$) + 1) - 1))
            END IF
         #:endfor
      END ASSOCIATE

   END SUBROUTINE get_arrays

! **************************************************************************************************
!> \brief get ith array
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE get_ith_array(list, i, array_size, array)
      TYPE(array_list), INTENT(IN)                    :: list
      INTEGER, INTENT(IN)                             :: i
      INTEGER, INTENT(IN)                             :: array_size
      INTEGER, DIMENSION(array_size), INTENT(OUT)     :: array

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
         CPASSERT(i <= number_of_arrays(list))

         array(:) = col_data(ptr(i):ptr(i + 1) - 1)

      END ASSOCIATE

   END SUBROUTINE

! **************************************************************************************************
!> \brief get ith array
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE allocate_and_get_ith_array(list, i, array)
      TYPE(array_list), INTENT(IN)                    :: list
      INTEGER, INTENT(IN)                             :: i
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: array

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
         CPASSERT(i <= number_of_arrays(list))

         ALLOCATE (array, source=col_data(ptr(i):ptr(i + 1) - 1))
      END ASSOCIATE
   END SUBROUTINE

! **************************************************************************************************
!> \brief sizes of arrays stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION sizes_of_arrays(list)
      TYPE(array_list), INTENT(IN)       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sizes_of_arrays

      INTEGER                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sizes_of_arrays(num_data))
      DO i_data = 1, num_data
         sizes_of_arrays(i_data) = list%ptr(i_data + 1) - list%ptr(i_data)
      END DO
   END FUNCTION sizes_of_arrays

! **************************************************************************************************
!> \brief sum of all elements for each array stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION sum_of_arrays(list)
      TYPE(array_list), INTENT(IN)       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sum_of_arrays

      INTEGER                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sum_of_arrays(num_data))
      DO i_data = 1, num_data
         sum_of_arrays(i_data) = SUM(list%col_data(list%ptr(i_data):list%ptr(i_data + 1) - 1))
      END DO

   END FUNCTION sum_of_arrays

! **************************************************************************************************
!> \brief partial sums of array elements.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE array_offsets(list_in, list_out)
      TYPE(array_list), INTENT(IN)  :: list_in
      TYPE(array_list), INTENT(OUT) :: list_out

      INTEGER                       :: i_data, i_ptr, num_data, partial_sum

      num_data = number_of_arrays(list_in)
      ALLOCATE (list_out%ptr, source=list_in%ptr)
      ALLOCATE (list_out%col_data(SIZE(list_in%col_data)))
      DO i_data = 1, num_data
         partial_sum = 1
         DO i_ptr = list_out%ptr(i_data), list_out%ptr(i_data + 1) - 1
            list_out%col_data(i_ptr) = partial_sum
            partial_sum = partial_sum + list_in%col_data(i_ptr)
         END DO
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief reorder array list.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE reorder_arrays(list_in, list_out, order)
      TYPE(array_list), INTENT(IN)                     :: list_in
      TYPE(array_list), INTENT(OUT)                    :: list_out
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: ${varlist("data")}$
      INTEGER, DIMENSION(number_of_arrays(list_in)), &
         INTENT(IN)                                    :: order

      #:for ndim in ndims
         IF (number_of_arrays(list_in) == ${ndim}$) THEN
            CALL get_arrays(list_in, ${varlist("data", nmax=ndim)}$, i_selected=dbt_inverse_order(order))
            CALL create_array_list(list_out, number_of_arrays(list_in), &
                                   ${varlist("data", nmax=ndim)}$)
         END IF
      #:endfor

   END SUBROUTINE

! **************************************************************************************************
!> \brief check whether two array lists are equal
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION check_equal(list1, list2)
      TYPE(array_list), INTENT(IN)  :: list1, list2
      LOGICAL :: check_equal

      check_equal = array_eq_i(list1%col_data, list2%col_data) .AND. array_eq_i(list1%ptr, list2%ptr)
   END FUNCTION

! **************************************************************************************************
!> \brief check whether two arrays are equal
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION array_eq_i(arr1, arr2)
      INTEGER, INTENT(IN), DIMENSION(:) :: arr1
      INTEGER, INTENT(IN), DIMENSION(:) :: arr2
      LOGICAL                           :: array_eq_i

      array_eq_i = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)

   END FUNCTION

END MODULE dbt_array_list_methods
