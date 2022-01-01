!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE submatrix_types

   USE kinds,                           ONLY: dp
   USE util,                            ONLY: sort

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER                     :: extvec_alloc_factor = 2, extvec_initial_alloc = 32
   INTEGER, PARAMETER                     :: set_modulus = 257 ! determines the number of buckets, should be a prime

   TYPE :: extvec_type
      INTEGER, DIMENSION(:), ALLOCATABLE   :: darr
      INTEGER                              :: elements = 0, allocated = 0
   CONTAINS
      PROCEDURE :: insert => extvec_insert
      PROCEDURE :: reset => extvec_reset
   END TYPE extvec_type

   TYPE, PUBLIC :: set_type
      TYPE(extvec_type), DIMENSION(0:set_modulus - 1) :: data
      INTEGER, DIMENSION(:), ALLOCATABLE       :: sorted
      INTEGER                                  :: elements = 0
      LOGICAL                                  :: sorted_up_to_date = .FALSE.
   CONTAINS
      PROCEDURE :: insert => set_insert
      PROCEDURE :: reset => set_reset
      PROCEDURE :: find => set_find
      PROCEDURE :: get => set_get
      PROCEDURE :: getall => set_getall
      PROCEDURE :: update_sorted => set_update_sorted
   END TYPE set_type

   TYPE, PUBLIC :: intBuffer_type
      INTEGER, DIMENSION(:), POINTER                    :: data
      INTEGER                                           :: size = 0
      LOGICAL                                           :: allocated = .FALSE.
      INTEGER                                           :: mpi_request = 0
   CONTAINS
      PROCEDURE :: alloc => intbuffer_alloc
      PROCEDURE :: dealloc => intbuffer_dealloc
   END TYPE intBuffer_type

   ! TODO: Make data type generic
   TYPE, PUBLIC :: buffer_type
      REAL(KIND=dp), DIMENSION(:), POINTER     :: data
      INTEGER                                           :: size = 0
      LOGICAL                                           :: allocated = .FALSE.
      INTEGER                                           :: mpi_request = 0
   CONTAINS
      PROCEDURE :: alloc => buffer_alloc
      PROCEDURE :: dealloc => buffer_dealloc
   END TYPE buffer_type

   TYPE, PUBLIC :: bufptr_type
      REAL(KIND=dp), DIMENSION(:), POINTER :: target => NULL()
   END TYPE bufptr_type

   TYPE, PUBLIC :: setarray_type
      TYPE(set_type), DIMENSION(:), ALLOCATABLE :: sets
   END TYPE setarray_type

CONTAINS

! **************************************************************************************************
!> \brief insert element into extendable vector
!> \param this - instance of extvec_type
!> \param elem - element to insert
! **************************************************************************************************
   PURE SUBROUTINE extvec_insert(this, elem)
      CLASS(extvec_type), INTENT(INOUT)       :: this
      INTEGER, INTENT(IN)                     :: elem
      INTEGER, DIMENSION(:), ALLOCATABLE      :: tmp

      IF (this%allocated .EQ. 0) THEN
         this%allocated = extvec_initial_alloc
         ALLOCATE (this%darr(this%allocated))
      ELSE
         IF (this%elements .EQ. this%allocated) THEN
            ALLOCATE (tmp(this%allocated))
            tmp(:) = this%darr
            DEALLOCATE (this%darr)
            ALLOCATE (this%darr(this%allocated*extvec_alloc_factor))
            this%darr(1:this%allocated) = tmp
            DEALLOCATE (tmp)
            this%allocated = this%allocated*extvec_alloc_factor
         END IF
      END IF

      this%elements = this%elements + 1
      this%darr(this%elements) = elem
   END SUBROUTINE extvec_insert

! **************************************************************************************************
!> \brief purge extendable vector and free allocated memory
!> \param this - instance of extvec_type
! **************************************************************************************************
   PURE SUBROUTINE extvec_reset(this)
      CLASS(extvec_type), INTENT(INOUT) :: this

      IF (ALLOCATED(this%darr)) DEALLOCATE (this%darr)
      this%allocated = 0
      this%elements = 0
   END SUBROUTINE extvec_reset

! **************************************************************************************************
!> \brief insert element into set
!> \param this - instance of set_type
!> \param elem - element to insert
! **************************************************************************************************
   PURE SUBROUTINE set_insert(this, elem)
      CLASS(set_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)            :: elem

      IF (.NOT. this%find(elem)) THEN
         CALL this%data(MODULO(elem, set_modulus))%insert(elem)
         this%sorted_up_to_date = .FALSE.
         this%elements = this%elements + 1
      END IF

   END SUBROUTINE set_insert

! **************************************************************************************************
!> \brief purse set and free allocated memory
!> \param this - instance of set_type
! **************************************************************************************************
   PURE SUBROUTINE set_reset(this)
      CLASS(set_type), INTENT(INOUT) :: this
      INTEGER                        :: i

      DO i = 0, set_modulus - 1
         CALL this%data(i)%reset
      END DO
      IF (ALLOCATED(this%sorted)) DEALLOCATE (this%sorted)
      this%elements = 0
      this%sorted_up_to_date = .FALSE.
   END SUBROUTINE set_reset

! **************************************************************************************************
!> \brief find element in set
!> \param this - instance of set_type
!> \param elem - element to look for
!> \return .TRUE. if element is contained in set, .FALSE. otherwise
! **************************************************************************************************
   PURE FUNCTION set_find(this, elem) RESULT(found)
      CLASS(set_type), INTENT(IN)   :: this
      INTEGER, INTENT(IN)           :: elem
      LOGICAL                       :: found
      INTEGER                       :: i, idx

      found = .FALSE.
      idx = MODULO(elem, set_modulus)

      DO i = 1, this%data(idx)%elements
         IF (this%data(idx)%darr(i) .EQ. elem) THEN
            found = .TRUE.
            EXIT
         END IF
      END DO

   END FUNCTION set_find

! **************************************************************************************************
!> \brief get element from specific position in set
!> \param this - instance of set_type
!> \param idx - position in set
!> \return element at position idx
! **************************************************************************************************
   FUNCTION set_get(this, idx) RESULT(elem)
      CLASS(set_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)            :: idx
      INTEGER                        :: elem

      IF (.NOT. this%sorted_up_to_date) CALL this%update_sorted

      elem = this%sorted(idx)
   END FUNCTION set_get

! **************************************************************************************************
!> \brief get all elements in set as sorted list
!> \param this - instance of set_type
!> \return sorted array containing set elements
! **************************************************************************************************
   FUNCTION set_getall(this) RESULT(darr)
      CLASS(set_type), INTENT(INOUT)           :: this
      INTEGER, DIMENSION(this%elements)        :: darr

      IF (.NOT. this%sorted_up_to_date) CALL this%update_sorted

      darr = this%sorted
   END FUNCTION set_getall

! **************************************************************************************************
!> \brief update internal list of set elements
!> \param this - instance of extendable vector
! **************************************************************************************************
   SUBROUTINE set_update_sorted(this)
      CLASS(set_type), INTENT(INOUT)     :: this
      INTEGER                            :: i, idx
      INTEGER, DIMENSION(:), ALLOCATABLE :: tmp

      IF (ALLOCATED(this%sorted)) DEALLOCATE (this%sorted)
      ALLOCATE (this%sorted(this%elements))

      idx = 1
      DO i = 0, set_modulus - 1
         IF (this%data(i)%elements .GT. 0) THEN
            this%sorted(idx:idx + this%data(i)%elements - 1) = this%data(i)%darr(1:this%data(i)%elements)
            idx = idx + this%data(i)%elements
         END IF
      END DO

      ALLOCATE (tmp(this%elements))
      CALL sort(this%sorted, this%elements, tmp)
      DEALLOCATE (tmp)

      this%sorted_up_to_date = .TRUE.
   END SUBROUTINE set_update_sorted

! **************************************************************************************************
!> \brief allocate buffer
!> \param this - instance of buffer_type
!> \param elements - number of elements contained in buffer
! **************************************************************************************************
   PURE SUBROUTINE buffer_alloc(this, elements)
      CLASS(buffer_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)               :: elements

      ALLOCATE (this%data(elements))
      this%allocated = .TRUE.
      this%size = elements
   END SUBROUTINE buffer_alloc

! **************************************************************************************************
!> \brief deallocate buffer
!> \param this - instance of buffer_type
! **************************************************************************************************
   PURE SUBROUTINE buffer_dealloc(this)
      CLASS(buffer_type), INTENT(INOUT) :: this

      IF (this%allocated) DEALLOCATE (this%data)
      this%allocated = .FALSE.
      this%size = 0
   END SUBROUTINE buffer_dealloc

! **************************************************************************************************
!> \brief allocate integer buffer
!> \param this - instance of intBuffer_type
!> \param elements - number of elements contained in buffer
! **************************************************************************************************
   PURE SUBROUTINE intbuffer_alloc(this, elements)
      CLASS(intBuffer_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                  :: elements

      ALLOCATE (this%data(elements))
      this%allocated = .TRUE.
      this%size = elements
   END SUBROUTINE intbuffer_alloc

! **************************************************************************************************
!> \brief deallocate integer buffer
!> \param this - instance of intBuffer_type
! **************************************************************************************************
   PURE SUBROUTINE intbuffer_dealloc(this)
      CLASS(intBuffer_type), INTENT(INOUT) :: this

      IF (this%allocated) DEALLOCATE (this%data)
      this%allocated = .FALSE.
      this%size = 0
   END SUBROUTINE intbuffer_dealloc

END MODULE submatrix_types
