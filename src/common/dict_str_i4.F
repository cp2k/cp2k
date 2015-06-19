!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief A dictionary (also known as hashtable or hashmap).
!>        Internally the dictionary uses an array to holds its data.
!>        If this array reaches a load-factor of 75%, a new array with twice the
!>        size will be allocated and the items are then copied over.
!>        This ensures that the dictionary will perform operations in O(1).
!> \note
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify dict__keytype___valuetype_.template instead
!>
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! ******************************************************************************

#define isequal_str(k1,k2) (k1==k2)

MODULE dict_str_i4
  USE kinds,                           ONLY: default_string_length,&
                                             int_4,&
                                             int_8

  IMPLICIT NONE
  PRIVATE



  !this is an internal type
  !Calculating hashes might be expensive, therefore they are stored
  !for use during change_capacity().
  TYPE private_item_type
     PRIVATE
     CHARACTER(LEN=default_string_length)                            :: key
     INTEGER(kind=int_4)                          :: value
     INTEGER(KIND=int_8)                       :: hash
     TYPE(private_item_type), POINTER          :: next => Null()
  END TYPE private_item_type

  !this is an internal type
  TYPE private_item_p_type
      PRIVATE
      TYPE(private_item_type), POINTER :: p => Null()  
  END TYPE private_item_p_type

  ! this is the public type, which holds a dictionary-instance
  TYPE dict_str_i4_type
     PRIVATE
     TYPE(private_item_p_type), DIMENSION(:), POINTER      :: buckets => Null()
     INTEGER                                               :: size = -1
  END TYPE dict_str_i4_type

  ! this is a public type, its returned by dict_items()
  TYPE dict_str_i4_item_type
     CHARACTER(LEN=default_string_length)        :: key
     INTEGER(kind=int_4)      :: value
  END TYPE dict_str_i4_item_type

  PUBLIC :: dict_str_i4_type
  PUBLIC :: dict_str_i4_init
  PUBLIC :: dict_str_i4_isready
  PUBLIC :: dict_str_i4_set
  PUBLIC :: dict_str_i4_get
  PUBLIC :: dict_str_i4_del
  PUBLIC :: dict_str_i4_size
  PUBLIC :: dict_str_i4_haskey
  PUBLIC :: dict_str_i4_items
  PUBLIC :: dict_str_i4_update
  PUBLIC :: dict_str_i4_destroy

  CONTAINS


! *****************************************************************************
! This is joaat_hash from string_table.F
!
!> \brief generates the hash of a given string
!> \param key a string of any length
!> \retval hash ...
!> \par History
!>       09.2006 created [Joost VandeVondele]
!>       12.2012 copied and adopted [ole]
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
! *****************************************************************************
FUNCTION hash_str(key) RESULT(hash)
    CHARACTER(LEN=*), INTENT(IN)             :: key
    INTEGER(KIND=int_8)                      :: hash

    INTEGER(KIND=int_8), PARAMETER           :: b32 = 2_int_8**32-1_int_8

    INTEGER                                  :: i

    hash=0_int_8
    DO i=1,LEN(key)
       hash=IAND(hash+ICHAR(key(i:i))                ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)  ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32)) ,b32)
    ENDDO
    hash=IAND(     hash+IAND(ISHFT(hash,  3),b32)  ,b32)
    hash=IAND(IEOR(hash,IAND(ISHFT(hash,-11),b32)) ,b32)
    hash=IAND(     hash+IAND(ISHFT(hash, 15),b32)  ,b32)
END FUNCTION hash_str


! *****************************************************************************
!> \brief Allocates the internal data-structures of the given dictionary.
!> \param dict ...
!> \param initial_capacity The initial size of the internal array (default=11).
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dict_str_i4_init(dict, initial_capacity)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    INTEGER, INTENT(in), OPTIONAL            :: initial_capacity

    INTEGER                                  :: initial_capacity_

    initial_capacity_ = 11
    IF(PRESENT(initial_capacity)) initial_capacity_ = initial_capacity

    IF(initial_capacity_ < 1) & 
      STOP "dict_str_i4_init: initial_capacity < 1"

    IF(ASSOCIATED(dict%buckets)) &
      STOP "dict_str_i4_init: dictionary is already initialized."

    ALLOCATE(dict%buckets(initial_capacity_))
    dict%size = 0

END SUBROUTINE dict_str_i4_init


! *****************************************************************************
!> \brief Test if the given dictionary has been initialized.
!> \param dict ...
!> \retval res ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION dict_str_i4_isready(dict) RESULT(res)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    LOGICAL                                  :: res

    res = ASSOCIATED(dict%buckets)
END FUNCTION dict_str_i4_isready


! *****************************************************************************
!> \brief Deallocated the internal data-structures if the given dictionary.
!>        Caution: If the stored keys or values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param dict ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dict_str_i4_destroy(dict)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict

    INTEGER                                  :: i
    TYPE(private_item_type), POINTER         :: item, prev_item

    IF(.NOT. ASSOCIATED(dict%buckets)) &
      STOP "dict_str_i4_destroy: dictionary is not initialized."

    DO i = 1, SIZE(dict%buckets)
        item => dict%buckets(i)%p
        DO WHILE(ASSOCIATED(item))
           prev_item => item
           item => item%next
           DEALLOCATE(prev_item)
        END DO
    END DO

    DEALLOCATE(dict%buckets)
    dict%size = -1

END SUBROUTINE dict_str_i4_destroy

! *****************************************************************************
!> \brief Stores, and possibly overwrites, a given value under a given key.
!> \param dict ...
!> \param key ...
!> \param value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dict_str_i4_set(dict, key, value)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    CHARACTER(LEN=default_string_length), &
      INTENT(in)                             :: key
    INTEGER(kind=int_4), INTENT(in)          :: value

    INTEGER(KIND=int_8)                      :: hash

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_set: dictionary is not initialized."

    hash = hash_str(key)
    CALL set_hashed(dict, key, value, hash)
END SUBROUTINE dict_str_i4_set


! *****************************************************************************
!> \brief Common code used internally by dict_set() and change_capacity().
!> \param dict ...
!> \param key ...
!> \param value ...
!> \param hash ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
RECURSIVE SUBROUTINE set_hashed(dict, key, value, hash)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    CHARACTER(LEN=default_string_length), &
      INTENT(in)                             :: key
    INTEGER(kind=int_4), INTENT(in)          :: value
    INTEGER(KIND=int_8), INTENT(in)          :: hash

    INTEGER(KIND=int_8)                      :: idx
    TYPE(private_item_type), POINTER         :: item, new_item

    idx = MOD(hash, INT(SIZE(dict%buckets),KIND=int_8)) + 1

    ! if allready in dict just update its value
    item => dict%buckets(idx)%p
    DO WHILE(ASSOCIATED(item))
        IF(item%hash == hash) THEN
           IF(isequal_str(item%key, key)) THEN
                item%value = value
               RETURN
           END IF
        END IF
        item => item%next
    END DO

    ! check load-factor
    IF(4*dict%size > 3*SIZE(dict%buckets)) THEN ! load-factor > 75%
       CALL change_capacity(dict, 2*SIZE(dict%buckets)) !double capacity
       idx = MOD(hash, INT(SIZE(dict%buckets),KIND=int_8)) + 1
    END IF


    ! create a new item
    ALLOCATE(new_item)
    new_item%hash = hash
    new_item%key = key
    new_item%value = value
    new_item%next => dict%buckets(idx)%p
    dict%buckets(idx)%p => new_item
    dict%size = dict%size + 1

END SUBROUTINE set_hashed

! *****************************************************************************
!> \brief Internal routine for changing the dictionary's capacity.
!> \param dict ...
!> \param new_capacity ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
RECURSIVE SUBROUTINE change_capacity(dict, new_capacity)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    INTEGER, INTENT(in)                      :: new_capacity

    INTEGER                                  :: i, new_cap, old_size
    TYPE(private_item_p_type), &
      DIMENSION(:), POINTER                  :: old_buckets
    TYPE(private_item_type), POINTER         :: item, prev_item

    new_cap = new_capacity
    ! pre checks
    IF(new_cap > HUGE(i)) THEN
       IF(SIZE(dict%buckets) == HUGE(i)) RETURN ! reached maximum - stay there.
       new_cap = HUGE(i) ! grow as far as possible
    END IF
    IF(new_cap < 1) &
        STOP "dict_str_i4_change_capacity: new_capacity < 1."
    IF(4*dict%size > 3*new_cap) &
     STOP "dict_str_i4_change_capacity: new_capacity too small."

    old_size = dict%size
    old_buckets => dict%buckets
    ALLOCATE(dict%buckets(new_capacity))
    dict%size = 0
    DO i = 1, SIZE(old_buckets)
        item => old_buckets(i)%p
        DO WHILE(ASSOCIATED(item))
           CALL set_hashed(dict, item%key, item%value, item%hash)
           prev_item => item
           item => item%next
           DEALLOCATE(prev_item)
        END DO
    END DO

    DEALLOCATE(old_buckets)

    IF(old_size /= dict%size) &
      STOP "dict_str_i4_change_capacity: assertion failed"
END SUBROUTINE change_capacity

! *****************************************************************************
!> \brief Gets a value for a given key from the dictionary.
!>        If the key is not found the default_value will be returned.
!>        If the key is not found and default_value was not provided the program stops.
!> \param dict ...
!> \param key ...
!> \param default_value ...
!> \retval value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION dict_str_i4_get(dict, key, default_value) RESULT(value)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    CHARACTER(LEN=default_string_length)     :: key
    INTEGER(kind=int_4), INTENT(in), &
      OPTIONAL                               :: default_value
    INTEGER(kind=int_4)                      :: value

    INTEGER(KIND=int_8)                      :: hash, idx
    TYPE(private_item_type), POINTER         :: item

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_get: dictionary is not initialized."

    hash = hash_str(key)
    idx = MOD(hash, INT(SIZE(dict%buckets),KIND=int_8)) + 1

    item => dict%buckets(idx)%p
    DO WHILE(ASSOCIATED(item))
        IF(item%hash == hash) THEN
           IF(isequal_str(item%key, key)) THEN
                value = item%value
                RETURN
           END IF
        END IF
        item => item%next
    END DO

    IF(PRESENT(default_value)) THEN
        value = default_value
        RETURN
    END IF

    STOP "dict_str_i4_get: Key not found in dictionary."
END FUNCTION dict_str_i4_get


! *****************************************************************************
!> \brief Remove the value for a given key from the dictionary.
!>        If the key is not found the proigram stops.           
!> \param dict ...
!> \param key ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dict_str_i4_del(dict, key)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    CHARACTER(LEN=default_string_length)     :: key

    INTEGER(KIND=int_8)                      :: hash, idx
    TYPE(private_item_type), POINTER         :: item, prev_item

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_del: dictionary is not initialized."

    hash = hash_str(key)
    idx = MOD(hash, INT(SIZE(dict%buckets),KIND=int_8)) + 1

    item => dict%buckets(idx)%p
    prev_item => Null()
    DO WHILE(ASSOCIATED(item))
        IF(item%hash == hash) THEN
           IF(isequal_str(item%key, key)) THEN
               IF(ASSOCIATED(prev_item)) THEN
                   prev_item%next => item%next
               ELSE
                   dict%buckets(idx)%p => item%next
               END IF
               DEALLOCATE(item)
               dict%size = dict%size - 1 
               RETURN
           END IF
        END IF
        prev_item => item
        item => item%next
    END DO

    STOP "dict_str_i4_del: Key not found in dictionary."
END SUBROUTINE dict_str_i4_del

! *****************************************************************************
!> \brief Returns the number of key/value-items currently stored in the dictionary.
!> \param dict ...
!> \retval size ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION dict_str_i4_size(dict) RESULT(size)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    INTEGER                                  :: size

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_size: dictionary is not initialized."

    size = dict%size
END FUNCTION dict_str_i4_size


! *****************************************************************************
!> \brief Checks whether a given key is currently stored in the dictionary. 
!> \param dict ...
!> \param key ...
!> \retval res ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION dict_str_i4_haskey(dict, key) RESULT(res)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    CHARACTER(LEN=default_string_length)     :: key
    LOGICAL                                  :: res

    INTEGER(KIND=int_8)                      :: hash, idx
    TYPE(private_item_type), POINTER         :: item

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_haskey: dictionary is not initialized."

    res = .FALSE.
    IF(dict%size==0) RETURN

    hash = hash_str(key)
    idx = MOD(hash, INT(SIZE(dict%buckets),KIND=int_8)) + 1

    item => dict%buckets(idx)%p
    DO WHILE(ASSOCIATED(item))
        IF(item%hash == hash) THEN
           IF(isequal_str(item%key, key)) THEN
               res = .TRUE.
               RETURN
           END IF
        END IF
        item => item%next
    END DO


END FUNCTION dict_str_i4_haskey


! *****************************************************************************
!> \brief Returns a pointer to an array of all key/value-items stored in the dictionary.
!>        Caution: The caller is responsible for deallocating targeted array after usage.        
!> \param dict ...
!> \retval items ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
FUNCTION dict_str_i4_items(dict) RESULT(items)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict
    TYPE(dict_str_i4_item_type), &
      DIMENSION(:), POINTER                  :: items

    INTEGER                                  :: i, j
    TYPE(private_item_type), POINTER         :: item

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_items: dictionary is not initialized."

    ALLOCATE(items(dict%size))
    j = 1
    DO i = 1, SIZE(dict%buckets)
        item => dict%buckets(i)%p
        DO WHILE(ASSOCIATED(item))
           items(j)%key = item%key
           items(j)%value = item%value
           j = j + 1
           item => item%next
        END DO
    END DO

    IF(j /= dict%size + 1) &
       STOP "dict_str_i4_items: assertion failed!"
END FUNCTION dict_str_i4_items


! *****************************************************************************
!> \brief Copies all key/values-items from one dictionary to another.
!>        Afterards dict will contain all items from the from_dict and 
!>        additionaly all its previous items, which were not overwritten.
!>        The two dictionaries have to be of the same type.
!> \param dict destination of items
!> \param from_dict source of items - will not be change
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE dict_str_i4_update(dict, from_dict)
    TYPE(dict_str_i4_type), INTENT(inout)    :: dict, from_dict

    INTEGER                                  :: i
    TYPE(dict_str_i4_item_type), &
      DIMENSION(:), POINTER                  :: from_items

    IF(.NOT. ASSOCIATED(dict%buckets)) &
       STOP "dict_str_i4_update: dictionary is not initialized."
    IF(.NOT. ASSOCIATED(from_dict%buckets)) &
       STOP "dict_str_i4_update: from_dict is not initialized."

    from_items => dict_str_i4_items(from_dict)
    DO i=1, SIZE(from_items)
       CALL dict_str_i4_set(dict, &
                  from_items(i)%key, from_items(i)%value)
    END DO
    DEALLOCATE(from_items)
END SUBROUTINE dict_str_i4_update


END MODULE dict_str_i4
