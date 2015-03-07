!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#define CP_I_LESS_Q(el1,el2) ( el1 < el2 )


! *****************************************************************************
!> \brief Routine for sorting an array
!> \note
!>      CP2K: 
!>      Please use the interface definded in util.F for calling sort().
!>      This could not go into cp_array__nametype1__utils.template because 
!>      there some routines use the error-handling from cp_common_uses.f90,
!>      which in turn uses timings.F, which uses sort() :-(
!>
!>      DBCSR:
!>      Please use the interface definded in dbcsr_toollib.F for calling sort().
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify _lib__array__nametype1__sort.template instead
!>
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! *****************************************************************************
MODULE cp_array_sort_i4
   
  USE kinds,                           ONLY: int_4

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_array_sort_i4'

  PUBLIC :: cp_1d_i4_sort
  
  CONTAINS


! *****************************************************************************
!> \brief Sorts an array inplace using a combination of merge- and bubble-sort.
!>        It also returns the indices, which the elements had before the sort.
!> \param arr the array to sort
!> \param n length of array
!> \param indices returns elements-indices before the sort 
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE cp_1d_i4_sort(arr, n, indices)
    INTEGER, INTENT(in)                      :: n
    INTEGER(kind=int_4), DIMENSION(1:n), &
      INTENT(inout)                          :: arr
    INTEGER, DIMENSION(1:n), INTENT(out)     :: indices

    INTEGER                                  :: i, stat
    INTEGER(kind=int_4), POINTER             :: tmp_arr(:)
    INTEGER, POINTER                         :: tmp_idx(:)

   IF(n == 0) RETURN ! for some reason this is a frequent case in cp2k
   
   ! scratch space used during the merge step
   ALLOCATE(tmp_arr((SIZE(arr)+1)/2), tmp_idx((SIZE(arr)+1)/2), stat=stat)
   IF(stat /= 0) STOP "cp_1d_i4_sort: allocation failed"
   
   indices = (/ (i, i = 1, SIZE(arr)) /)
   
   CALL cp_1d_i4_sort_low(arr(1:n), indices, tmp_arr, tmp_idx)
   
   DEALLOCATE(tmp_arr,tmp_idx)
   
END SUBROUTINE cp_1d_i4_sort



! *****************************************************************************
!> \brief The actual sort routing. 
!>            Only cp_1d_i4_sort and itself should call this.
!> \param arr the array to sort
!> \param indices elements-indices before the sort
!> \param tmp_arr scratch space
!> \param tmp_idx scratch space
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
RECURSIVE SUBROUTINE cp_1d_i4_sort_low(arr, indices, tmp_arr, tmp_idx)
    INTEGER(kind=int_4), DIMENSION(:), &
      INTENT(inout)                          :: arr
    INTEGER, DIMENSION(SIZE(arr)), &
      INTENT(inout)                          :: indices
    INTEGER(kind=int_4), &
      DIMENSION((SIZE(arr)+1)/2), &
      INTENT(inout)                          :: tmp_arr
    INTEGER, DIMENSION((SIZE(arr)+1)/2), &
      INTENT(inout)                          :: tmp_idx

    INTEGER                                  :: i, j, k, m, t
    INTEGER(kind=int_4)                      :: a
    LOGICAL                                  :: swapped

! a,t:  used during swaping of elements in arr and indices

   swapped=.TRUE.
   
   ! If only a few elements are left we switch to bubble-sort for efficiency.
   IF (SIZE(arr) <= 7) THEN ! 7 seems to be a good choice for the moment
          DO j = SIZE(arr)-1, 1, -1
            swapped = .FALSE.
            DO i = 1, j
                 IF (CP_I_LESS_Q(arr(i+1), arr(i))) THEN
                         ! swap arr(i) with arr(i+1)
                        a            = arr(i)
                        arr(i)       = arr(i+1)
                        arr(i+1)     = a
                        ! swap indices(i) with indices(i+1)
                        t            = indices(i)
                        indices(i)   = indices(i+1)
                        indices(i+1) = t
                        swapped      = .TRUE.
                 END IF
            END DO
        IF (.NOT. swapped) EXIT
        END DO
      RETURN
   ENDIF
   
   ! split list in half and recusively sort both sublists
   m = (SIZE(arr)+1)/2 ! index where we going to divide the list in two
   CALL cp_1d_i4_sort_low(arr(1:m), indices(1:m), tmp_arr, tmp_idx)
   CALL cp_1d_i4_sort_low(arr(m+1:), indices(m+1:), tmp_arr, tmp_idx)
   
   ! Check for a special case: Can we just concate the two sorted sublists?
   ! This leads to O(n) scaling if the input is already sorted. 
   IF (CP_I_LESS_Q(arr(m+1), arr(m))) THEN
      ! ...no - let's merge the two sorted sublists arr(:m) and arr(m+1:)
      ! Merge will be performed directly in arr. Need backup of first sublist.
      tmp_arr(1:m) = arr(1:m)
      tmp_idx(1:m) = indices(1:m) 
      i = 1; ! number of elemens consumed from 1st sublist 
      j = 1; ! number of elemens consumed from 2nd sublist 
      k = 1; ! number of elemens already merged
      
      DO WHILE(i <= m .AND. j <= SIZE(arr)-m )
      IF (CP_I_LESS_Q(arr(m+j), tmp_arr(i))) THEN
         arr(k)     = arr(m+j)
         indices(k) = indices(m+j)
         j = j+1
      ELSE
         arr(k)     = tmp_arr(i)
         indices(k) = tmp_idx(i)
         i = i+1
      ENDIF
      k = k + 1
      ENDDO
      
      ! One of the two sublist is now empty.
      ! Copy possibly remaining tail of 1st sublist
      DO WHILE (i <= m) 
         arr(k)     = tmp_arr(i)
         indices(k) = tmp_idx(i)
         i = i + 1
         k = k + 1
         ENDDO
         
     ! The possibly remaining tail of 2nd sublist is already at the right spot.

   ENDIF
   
END SUBROUTINE cp_1d_i4_sort_low


  ! template def put here so that line numbers in template and derived
  ! files are almost the same (multi-line use change it a bit)
  ! [template(type1,nametype1,USE,lessQ,defines)]
! ARGS:
!  USE = "USE kinds,       ONLY:int_4"
!  defines = 
!    "#define CP_I_LESS_Q(el1,el2) ( el1 < el2 )
!     "
!  lessQ = "CP_I_LESS_Q"
!  nametype1 = "i4"
!  type1 = "integer(kind=int_4)"


END MODULE cp_array_sort_i4
