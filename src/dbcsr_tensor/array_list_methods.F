!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
      sizes_of_arrays, &
      sum_of_arrays

   TYPE array_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: col_data
      INTEGER, DIMENSION(:), ALLOCATABLE :: ptr
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Get an element for each array.
!> \param list ...
!> \param indices element index for each array
!> \retval get_array_elements ...
! **************************************************************************************************
   FUNCTION get_array_elements(list, indices)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, DIMENSION(:), INTENT(IN)                  :: indices
      INTEGER, DIMENSION(SIZE(indices))                  :: get_array_elements

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
!> \param data1 first array
!> \param data2 second array
!> \param data3 ...
!> \param data4 ...
!> \param data5 ...
!> \param data6 ...
!> \param ndata_check ...
! **************************************************************************************************
   SUBROUTINE create_array_list(list, data1, data2, data3, data4, data5, data6, ndata_check)
      TYPE(array_list), INTENT(OUT)                      :: list
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: data1, data2, data3, data4, data5, data6
      INTEGER, INTENT(IN), OPTIONAL                      :: ndata_check

      INTEGER                                            :: ndata, ptr, size_all

      size_all = 0
      ndata = 0
      IF (present_and_allocated_i(data1)) THEN
         ndata = ndata+1
         size_all = size_all+SIZE(data1)
         IF (present_and_allocated_i(data2)) THEN
            ndata = ndata+1
            size_all = size_all+SIZE(data2)
            IF (present_and_allocated_i(data3)) THEN
               ndata = ndata+1
               size_all = size_all+SIZE(data3)
               IF (present_and_allocated_i(data4)) THEN
                  ndata = ndata+1
                  size_all = size_all+SIZE(data4)
                  IF (present_and_allocated_i(data5)) THEN
                     ndata = ndata+1
                     size_all = size_all+SIZE(data5)
                     IF (present_and_allocated_i(data6)) THEN
                        ndata = ndata+1
                        size_all = size_all+SIZE(data6)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (PRESENT(ndata_check)) THEN
         CPASSERT(ndata_check .EQ. ndata)
      ENDIF

      ALLOCATE (list%ptr(ndata+1))
      ALLOCATE (list%col_data(size_all))

      ptr = 1
      list%ptr(1) = ptr
      IF (present_and_allocated_i(data1)) THEN
         list%col_data(ptr:ptr+SIZE(data1)-1) = data1(:)
         ptr = ptr+SIZE(data1)
         list%ptr(2) = ptr
         IF (present_and_allocated_i(data2)) THEN
            list%col_data(ptr:ptr+SIZE(data2)-1) = data2(:)
            ptr = ptr+SIZE(data2)
            list%ptr(3) = ptr
            IF (present_and_allocated_i(data3)) THEN
               list%col_data(ptr:ptr+SIZE(data3)-1) = data3(:)
               ptr = ptr+SIZE(data3)
               list%ptr(4) = ptr
               IF (present_and_allocated_i(data4)) THEN
                  list%col_data(ptr:ptr+SIZE(data4)-1) = data4(:)
                  ptr = ptr+SIZE(data4)
                  list%ptr(5) = ptr
                  IF (present_and_allocated_i(data5)) THEN
                     list%col_data(ptr:ptr+SIZE(data5)-1) = data5(:)
                     ptr = ptr+SIZE(data5)
                     list%ptr(6) = ptr
                     IF (present_and_allocated_i(data6)) THEN
                        list%col_data(ptr:ptr+SIZE(data6)-1) = data6(:)
                        ptr = ptr+SIZE(data6)
                        list%ptr(7) = ptr
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
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
!> \param data1 first array
!> \param data2 second array
!> \param data3 ...
!> \param data4 ...
!> \param data5 ...
!> \param data6 ...
!> \param order order in which to retrieve arrays
! **************************************************************************************************
   SUBROUTINE get_arrays(list, data1, data2, data3, data4, data5, data6, order)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: data1, data2, data3, data4, data5, data6
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: order

      INTEGER                                            :: I, ndata
      INTEGER, DIMENSION(6)                              :: o

      ndata = number_of_arrays(list)

      o(:) = 0
      IF (PRESENT(order)) THEN
         o(1:ndata) = order(:)
      ELSE
         o(1:ndata) = (/(I, I=1, ndata)/)
      ENDIF

      ASSOCIATE(ptr=>list%ptr, col_data=>list%col_data)
      IF (ndata > 0) THEN
         ALLOCATE (data1(ptr(o(1)+1)-ptr(o(1)))); data1(:) = col_data(ptr(o(1)):ptr(o(1)+1)-1)
         IF (ndata > 1) THEN
            ALLOCATE (data2(ptr(o(2)+1)-ptr(o(2)))); data2(:) = col_data(ptr(o(2)):ptr(o(2)+1)-1)
            IF (ndata > 2) THEN
               ALLOCATE (data3(ptr(o(3)+1)-ptr(o(3)))); data3(:) = col_data(ptr(o(3)):ptr(o(3)+1)-1)
               IF (ndata > 3) THEN
                  ALLOCATE (data4(ptr(o(4)+1)-ptr(o(4)))); data4(:) = col_data(ptr(o(4)):ptr(o(4)+1)-1)
                  IF (ndata > 4) THEN
                     ALLOCATE (data5(ptr(o(5)+1)-ptr(o(5)))); data5(:) = col_data(ptr(o(5)):ptr(o(5)+1)-1)
                     IF (ndata > 5) THEN
                        ALLOCATE (data6(ptr(o(6)+1)-ptr(o(6)))); data6(:) = col_data(ptr(o(6)):ptr(o(6)+1)-1)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      END ASSOCIATE

   END SUBROUTINE get_arrays

! **************************************************************************************************
!> \brief number of arrays stored in list
!> \param list ...
!> \retval number_of_arrays ...
! **************************************************************************************************
   FUNCTION number_of_arrays(list)
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER                                            :: number_of_arrays

      number_of_arrays = SIZE(list%ptr)-1

   END FUNCTION number_of_arrays

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
!> \brief check if argument is present and allocated
!> \param arr ...
!> \retval present_and_allocated_i ...
! **************************************************************************************************
   FUNCTION present_and_allocated_i(arr)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: arr
      LOGICAL                                            :: present_and_allocated_i

      present_and_allocated_i = .FALSE.

      IF (PRESENT(arr)) THEN
         IF (ALLOCATED(arr)) THEN
            present_and_allocated_i = .TRUE.
         ENDIF
      ENDIF

   END FUNCTION

END MODULE array_list_methods

