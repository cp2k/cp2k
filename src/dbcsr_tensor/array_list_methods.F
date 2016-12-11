!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Generated from array_list_methods.fypp using Fypp.                                               !
! **DO NOT** modify this file, edit array_list_methods.fypp instead.                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Representation of arbitrary number of 1d integer arrays with arbitrary sizes.
!>
!> This is needed for generic handling of dimension-specific tensor quantities (such as block index).
!> \author Patrick Seewald
! **************************************************************************************************
MODULE array_list_methods

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'array_list_methods'

   PUBLIC  :: &
      array_list, &
      array_offsets, &
      create_array_list, &
      destroy_array_list, &
      get_array_elements, &
      get_arrays, &
      number_of_arrays, &
      reorder_arrays, &
      sizes_of_arrays, &
      sum_of_arrays

   TYPE array_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: col_data
      INTEGER, DIMENSION(:), ALLOCATABLE :: ptr
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief number of arrays stored in list
!> \param list ...
!> \retval number_of_arrays ...
! **************************************************************************************************
   PURE FUNCTION number_of_arrays(list)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER                                            :: number_of_arrays

      number_of_arrays = SIZE(list%ptr)-1

   END FUNCTION number_of_arrays

! **************************************************************************************************
!> \brief Get an element for each array.
!> \param list ...
!> \param indices element index for each array
!> \retval get_array_elements ...
! **************************************************************************************************
   FUNCTION get_array_elements(list, indices)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, DIMENSION(number_of_arrays(list)), &
         INTENT(IN)                                      :: indices
      INTEGER, DIMENSION(number_of_arrays(list))         :: get_array_elements

      INTEGER                                            :: i, ind

      DO i = 1, SIZE(indices)
         ind = indices(i)+list%ptr(i)-1
         get_array_elements(i) = list%col_data(ind)
      ENDDO

   END FUNCTION get_array_elements

! **************************************************************************************************
!> \brief collects any number of arrays of different sizes into a single array (list%col_data),
!>        storing the indices that start a new array (list%ptr).
!> \param list ...
!> \param ndata ...
!> \param data_1 ...
!> \param data_2 ...
!> \param data_3 ...
!> \param data_4 ...
! **************************************************************************************************
   SUBROUTINE create_array_list(list, ndata, data_1, data_2, data_3, data_4)
      TYPE(array_list), INTENT(OUT)                      :: list
      INTEGER, INTENT(IN)                                :: ndata
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: data_1, data_2, data_3, data_4

      INTEGER                                            :: ptr, size_all

      size_all = 0

      IF (ndata .GE. 1) THEN
         CPASSERT(PRESENT(data_1))
         size_all = size_all+SIZE(data_1)
      ENDIF
      IF (ndata .GE. 2) THEN
         CPASSERT(PRESENT(data_2))
         size_all = size_all+SIZE(data_2)
      ENDIF
      IF (ndata .GE. 3) THEN
         CPASSERT(PRESENT(data_3))
         size_all = size_all+SIZE(data_3)
      ENDIF
      IF (ndata .GE. 4) THEN
         CPASSERT(PRESENT(data_4))
         size_all = size_all+SIZE(data_4)
      ENDIF

      ALLOCATE (list%ptr(ndata+1))
      ALLOCATE (list%col_data(size_all))

      ptr = 1
      list%ptr(1) = ptr

      IF (ndata .GE. 1) THEN
         list%col_data(ptr:ptr+SIZE(data_1)-1) = data_1(:)
         ptr = ptr+SIZE(data_1)
         list%ptr(2) = ptr
      ENDIF
      IF (ndata .GE. 2) THEN
         list%col_data(ptr:ptr+SIZE(data_2)-1) = data_2(:)
         ptr = ptr+SIZE(data_2)
         list%ptr(3) = ptr
      ENDIF
      IF (ndata .GE. 3) THEN
         list%col_data(ptr:ptr+SIZE(data_3)-1) = data_3(:)
         ptr = ptr+SIZE(data_3)
         list%ptr(4) = ptr
      ENDIF
      IF (ndata .GE. 4) THEN
         list%col_data(ptr:ptr+SIZE(data_4)-1) = data_4(:)
         ptr = ptr+SIZE(data_4)
         list%ptr(5) = ptr
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief destroy array list.
!> \param list ...
! **************************************************************************************************
   SUBROUTINE destroy_array_list(list)
      TYPE(array_list), INTENT(INOUT)                    :: list

      DEALLOCATE (list%ptr, list%col_data)
   END SUBROUTINE

! **************************************************************************************************
!> \brief helper Get all arrays contained in list
!> \param list ...
!> \param data_1 ...
!> \param data_2 ...
!> \param data_3 ...
!> \param data_4 ...
!> \param order order in which to retrieve arrays
! **************************************************************************************************
   SUBROUTINE get_arrays(list, data_1, data_2, data_3, data_4, order)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: data_1, data_2, data_3, data_4
      INTEGER, DIMENSION(number_of_arrays(list)), &
         INTENT(IN), OPTIONAL                            :: order

      INTEGER                                            :: i, ndata
      INTEGER, DIMENSION(number_of_arrays(list))         :: o

      ndata = number_of_arrays(list)

      o(:) = 0
      IF (PRESENT(order)) THEN
         o(1:ndata) = order(:)
      ELSE
         o(1:ndata) = (/(i, i=1, ndata)/)
      ENDIF

      ASSOCIATE(ptr=>list%ptr, col_data=>list%col_data)
      IF (ndata > 0) THEN
         ALLOCATE (data_1(ptr(o(1)+1)-ptr(o(1))))
         data_1(:) = col_data(ptr(o(1)):ptr(o(1)+1)-1)
      ENDIF
      IF (ndata > 1) THEN
         ALLOCATE (data_2(ptr(o(2)+1)-ptr(o(2))))
         data_2(:) = col_data(ptr(o(2)):ptr(o(2)+1)-1)
      ENDIF
      IF (ndata > 2) THEN
         ALLOCATE (data_3(ptr(o(3)+1)-ptr(o(3))))
         data_3(:) = col_data(ptr(o(3)):ptr(o(3)+1)-1)
      ENDIF
      IF (ndata > 3) THEN
         ALLOCATE (data_4(ptr(o(4)+1)-ptr(o(4))))
         data_4(:) = col_data(ptr(o(4)):ptr(o(4)+1)-1)
      ENDIF
      END ASSOCIATE

   END SUBROUTINE get_arrays

! **************************************************************************************************
!> \brief sizes of arrays stored in list
!> \param list ...
!> \retval sizes_of_arrays ...
! **************************************************************************************************
   FUNCTION sizes_of_arrays(list)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_of_arrays

      INTEGER                                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sizes_of_arrays(num_data))
      DO i_data = 1, num_data
         sizes_of_arrays(i_data) = list%ptr(i_data+1)-list%ptr(i_data)
      ENDDO
   END FUNCTION sizes_of_arrays

! **************************************************************************************************
!> \brief sum of all elements for each array stored in list
!> \param list ...
!> \retval sum_of_arrays ...
! **************************************************************************************************
   FUNCTION sum_of_arrays(list)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sum_of_arrays

      INTEGER                                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sum_of_arrays(num_data))
      DO i_data = 1, num_data
         sum_of_arrays(i_data) = SUM(list%col_data(list%ptr(i_data):list%ptr(i_data+1)-1))
      ENDDO

   END FUNCTION sum_of_arrays

! **************************************************************************************************
!> \brief partial sums of array elements.
!> \param list_in ...
!> \param list_out ...
! **************************************************************************************************
   SUBROUTINE array_offsets(list_in, list_out)
      TYPE(array_list), INTENT(IN)                       :: list_in
      TYPE(array_list), INTENT(OUT)                      :: list_out

      INTEGER                                            :: i_data, i_ptr, num_data, partial_sum

      num_data = number_of_arrays(list_in)
      ALLOCATE (list_out%ptr(num_data+1))
      list_out%ptr(:) = list_in%ptr(:)
      ALLOCATE (list_out%col_data(SIZE(list_in%col_data)))
      DO i_data = 1, num_data
         partial_sum = 1
         DO i_ptr = list_out%ptr(i_data), list_out%ptr(i_data+1)-1
            list_out%col_data(i_ptr) = partial_sum
            partial_sum = partial_sum+list_in%col_data(i_ptr)
         ENDDO
      ENDDO
   END SUBROUTINE

! **************************************************************************************************
!> \brief reorder array list.
!> \param list_in ...
!> \param list_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE reorder_arrays(list_in, list_out, order)
      TYPE(array_list), INTENT(IN)                       :: list_in
      TYPE(array_list), INTENT(OUT)                      :: list_out
      INTEGER, DIMENSION(number_of_arrays(list_in)), &
         INTENT(IN)                                      :: order

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data_1, data_2, data_3, data_4

      CALL get_arrays(list_in, data_1, data_2, data_3, data_4, order=order)
      CALL create_array_list(list_out, number_of_arrays(list_in), &
                             data_1, data_2, data_3, data_4)

   END SUBROUTINE

END MODULE array_list_methods

