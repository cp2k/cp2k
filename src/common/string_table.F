!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief generates a unique id number for a string (str2id) that can be used
!>      two compare two strings. I.e.
!>      if (str1==str2) => str2id(str1)==str2id(str2)
!>      if (str1.NE.str2) => str2id(str1).NE.str2id(str2)
!>      and the other way around. Given an id, the string can be retrieved.
!> \note
!>      the purpose of this routine is to speed up string handling,
!>      string searching, ... as an operation on an int is much faster than an
!>      operation on a long string.
!> \par History
!>      9.2006 [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE string_table
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             int_8
#include "cp_common_uses.f90"

  IMPLICIT NONE

  ! user functions
  PUBLIC :: str2id, id2str, s2s

  ! setup function
  PUBLIC :: string_table_allocate,string_table_deallocate

  PRIVATE
  ! For good performance, the hash table should be larger than the largest number
  ! of strings that will be saved, but the memory for an empty table is 16*hash_table_size
  ! the string_table should remain functional for up to ~ 2**32 strings
  INTEGER, PARAMETER :: Nbit=16
  INTEGER, PARAMETER :: hash_table_size=2**Nbit

  ! actual elements in the hash table
  INTEGER, SAVE      :: actual_strings
  INTEGER, SAVE      :: inserted_strings

  ! an element of the linked list of hashed strings
! *****************************************************************************
  TYPE hash_element_type
    CHARACTER(LEN=default_string_length), POINTER :: str=>NULL()
    TYPE(hash_element_type), POINTER :: next=>NULL()
  END TYPE

  ! the array of linked lists of hashed strings
  TYPE(hash_element_type), SAVE, ALLOCATABLE, TARGET, DIMENSION(:) :: hash_table

CONTAINS

! *****************************************************************************
!> \brief returns a unique id for a given string, and stores the string for
!>      later retrieval using the id.
!> \param str the string to be stored (default_string_length)
!> \retval id ...
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      pass literal strings using the s2s function,
!>      which converts strings of any length to default_string_length
!>      id=str2id(s2s("my short string"))
! *****************************************************************************
  FUNCTION str2id(str) RESULT(id)
    CHARACTER(LEN=default_string_length)     :: str
    INTEGER                                  :: id

    INTEGER                                  :: index, ipos
    TYPE(hash_element_type), POINTER         :: this

     inserted_strings=inserted_strings+1
     ! index is the index in the array, ipos is the Nth element of the linked list
     index=joaat_hash(str)
     ipos =0
     this=>hash_table(index)
     DO ! walk the list
       IF (.NOT. ASSOCIATED(this%str)) THEN
          ! str was not in the linked list, add it now
          ALLOCATE(this%str)
          this%str=str
          actual_strings=actual_strings+1
          EXIT
       ELSE
          IF (this%str==str) THEN
             ! str is in the list already
             EXIT
          ELSE
             IF (.NOT. ASSOCIATED(this%next)) ALLOCATE(this%next)
             ipos=ipos+1
             this=>this%next
          ENDIF
       ENDIF
     ENDDO
     id=IOR(index,ISHFT(ipos,Nbit))
  END FUNCTION str2id

! *****************************************************************************
!> \brief returns the string associated with a given id
!> \param id the id to be converted into a string
!> \retval str ...
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      only id's of previously 'registered' strings (str2id) should be passed,
!>      otherwise things crash
! *****************************************************************************
  FUNCTION id2str(id) RESULT(str)
    INTEGER                                  :: id
    CHARACTER(LEN=default_string_length)     :: str

    INTEGER                                  :: i, index, ipos
    TYPE(hash_element_type), POINTER         :: this

     index=IAND(id,2**Nbit-1)
     ipos =ISHFT(id,-Nbit)
     this=>hash_table(index)
     DO i=1,ipos
        this=>this%next
     ENDDO
     str=this%str
  END FUNCTION id2str

! *****************************************************************************
!> \brief converts a string in a string of default_string_length
!> \param str ...
!> \retval res ...
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      useful to pass a literal string to str2id
!>      i.e. id=str2id(s2s("X"))
! *****************************************************************************
  FUNCTION s2s(str) RESULT(res)
    CHARACTER(LEN=*)                         :: str
    CHARACTER(LEN=default_string_length)     :: res

     res=str
  END FUNCTION s2s

! *****************************************************************************
!> \brief allocates the string table
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      this needs to be done only once at program startup, before any use
!>      of other procedures of this module. The scope of this table is global
! *****************************************************************************
  SUBROUTINE string_table_allocate()
     ALLOCATE(hash_table(0:hash_table_size-1))
     actual_strings=0
     inserted_strings=0
  END SUBROUTINE string_table_allocate

! *****************************************************************************
!> \brief deallocates the string table
!> \param iw a unit to which some info about the table usage can be printed
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      This should be done before program termination, all associated ids become meaningless
! *****************************************************************************
  SUBROUTINE string_table_deallocate(iw)
    INTEGER, INTENT(IN)                      :: iw

    INTEGER                                  :: i, ilist, ipos, ipos_max
    TYPE(hash_element_type), POINTER         :: next, this

! clean up all the linked lists of entries

     ipos_max=0
     ilist=0
     DO i=0,hash_table_size-1
        ipos=1
        IF (ASSOCIATED(hash_table(i)%str)) THEN
            DEALLOCATE(hash_table(i)%str)
            ilist=ilist+1
        ENDIF
        this=>hash_table(i)%next
        DO WHILE (ASSOCIATED(this))
           ipos=ipos+1
           next=>this%next
           IF (ASSOCIATED(this%str)) DEALLOCATE(this%str)
           DEALLOCATE(this)
           this=>next
        ENDDO
        ipos_max=MAX(ipos_max,ipos)
     ENDDO
     DEALLOCATE(hash_table)
     IF (iw>0) THEN
        WRITE(iw,*) "string table: # inserted str = ",inserted_strings
        WRITE(iw,*) "              # actual       = ",actual_strings
        WRITE(iw,*) "              # lists        = ",ilist," / ",hash_table_size
        WRITE(iw,*) "              longest list   = ",ipos_max
     ENDIF
     actual_strings=0
     inserted_strings=0
  END SUBROUTINE string_table_deallocate

! *****************************************************************************
!> \brief generates the hash of a string and the index in the table
!> \param key a string of any length
!> \retval hash_index ...
!> \par History
!>       09.2006 created [Joost VandeVondele]
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
!>       However, since fortran doesn't have an unsigned 4 byte int
!>       we compute it using an integer with the appropriate range
!>       we return already the index in the table as a final result
! *****************************************************************************
  FUNCTION joaat_hash(key) RESULT(hash_index)
    CHARACTER(LEN=*), INTENT(IN)             :: key
    INTEGER                                  :: hash_index

    INTEGER(KIND=int_8), PARAMETER           :: b32 = 2_int_8**32-1_int_8

    INTEGER                                  :: i
    INTEGER(KIND=int_8)                      :: hash

    hash=0_int_8
    DO i=1,LEN(key)
       hash=IAND(hash+ICHAR(key(i:i))                ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)  ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32)) ,b32)
    ENDDO
    hash=IAND(     hash+IAND(ISHFT(hash,  3),b32)  ,b32)
    hash=IAND(IEOR(hash,IAND(ISHFT(hash,-11),b32)) ,b32)
    hash=IAND(     hash+IAND(ISHFT(hash, 15),b32)  ,b32)
    ! hash is the real 32bit hash value of the string,
    ! hash_index is an index in the hash_table
    hash_index=MOD(hash,INT(hash_table_size,KIND=int_8))
  END FUNCTION joaat_hash
END MODULE string_table
