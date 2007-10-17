!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_str_cache [1.0] *
!!
!!   NAME
!!     cp_str_cache
!!
!!   FUNCTION
!!     Cache for strings.
!!
!!   NOTES
!!     fortran style 1 based indexing
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     10.2007 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_str_cache
  USE cp_files,                        ONLY: get_eof_stat
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast
#include "cp_common_uses.h"

  IMPLICIT NONE

  INTEGER, SAVE :: last_scache_id=0
  INTEGER, PARAMETER :: baseStrLen=80
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_str_cache'
!***
!****************************************************************************

!!****s* cp_str_cache/str_cache_type [1.0] *
!!
!!   NAME
!!     str_cache_type
!!
!!   FUNCTION
!!     Cache for strings.
!!     Stores multiple lines contiguosly, removes trailing spaces.
!!     Can be syncronized between processors.
!!     Grows automatically if needed.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference count
!!     - cache: the actual cache of the string
!!     - linePos: positions of the newlines. Newline is *after* the
!!       given character. Might be oversized
!!     - nLines: the actual number of lines
!!     - posAtt: actual position in the cache. Reading starts at this char.
!!       (inclusive lower bound)
!!     - lastChar: position of last charcater stored (inclusive upper bound)
!!
!!   SOURCE
!***************************************************************************
  TYPE str_cache_type
     INTEGER :: id_nr, ref_count
     CHARACTER(len=baseStrLen), DIMENSION(:), POINTER :: cache
     INTEGER, DIMENSION(:), POINTER :: linePos
     INTEGER :: nLines
     INTEGER :: posAtt
     INTEGER :: lastChar
  END TYPE str_cache_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_str_cache/sc_create [1.0] *
!!
!!   NAME
!!     sc_create
!!
!!   FUNCTION
!!     create a new string cache with the requested initial capacity
!!
!!   INPUTS
!!     - scache: the cache to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     - startSize: initial capacity (in characters)
!!
!!*** *********************************************************************
  SUBROUTINE sc_create(scache,error,start_capacity)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(in), OPTIONAL            :: start_capacity

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nLines, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.not.ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    ALLOCATE(scache,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    last_scache_id=last_scache_id+1
    scache%id_nr=last_scache_id
    scache%ref_count=1
    nLines=5
    IF (PRESENT(start_capacity)) THEN
       CPPrecondition(start_capacity>0,cp_failure_level,routineP,error,failure)
       nLines=(start_capacity+baseStrLen-1)/baseStrLen
    END IF
    ALLOCATE(scache%cache(nLines),stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    ALLOCATE(scache%linePos(nLines),stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    scache%nLines=0
    scache%posAtt=1
    scache%lastChar=0
  END SUBROUTINE sc_create
!****************************************************************************

!!****f* cp_str_cache/sc_retain [1.0] *
!!
!!   NAME
!!     sc_retain
!!
!!   FUNCTION
!!     retains the given scache
!!
!!   INPUTS
!!     - scache: the cache to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_retain(scache,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       scache%ref_count=scache%ref_count+1
    END IF
  END SUBROUTINE sc_retain
!****************************************************************************

!!****f* cp_str_cache/sc_release [1.0] *
!!
!!   NAME
!!     sc_release
!!
!!   function
!!     releases the given scache
!!
!!   INPUTS
!!     - scache: the cache to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_release(scache,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(scache)) THEN
       CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
       scache%ref_count=scache%ref_count-1
       IF (scache%ref_count==0) THEN
          DEALLOCATE(scache%cache,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(scache%linePos,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(scache,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(scache)
  END SUBROUTINE sc_release
!****************************************************************************

!!****f* cp_str_cache/sc_writeout [1.0] *
!!
!!   NAME
!!     sc_writeout
!!
!!   FUNCTION
!!     writes out the content of the cache
!!     for debugging purposes, always safe
!!
!!   INPUTS
!!     - scache: the cache to write out
!!     - out_unit: output unit were to write the info
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_writeout(scache,out_unit,error)
    TYPE(str_cache_type), POINTER            :: scache
    INTEGER, INTENT(in)                      :: out_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_writeout', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure=.FALSE.
    WRITE(out_unit,"('<str_cache_type ')",advance='NO')
    IF (.not.ASSOCIATED(scache)) THEN
       WRITE(out_unit,"('**NULL** >')")
       RETURN
    END IF
    WRITE(out_unit,"('id_nr',i6,' ref_count',i6)") scache%id_nr,scache%ref_count
    WRITE(out_unit,"('posAtt',i6,' lastChar',i6)") scache%posAtt,scache%lastChar
    WRITE(out_unit,"('nLines',i6,' linePos (')",advance='NO')
    IF (ASSOCIATED(scache%linePos)) THEN
       DO i=1,MIN(SIZE(scache%linePos),scache%nLines)
          WRITE(out_unit,"(i6)",advance='NO') scache%linePos(i)
       END DO
       WRITE(out_unit,"(')')")
    ELSE
       WRITE(out_unit,"('**NULL**)')")
    END IF
    IF (ASSOCIATED(scache%cache)) THEN
       WRITE(out_unit,"(' cache:')")
       DO i=1,MIN(SIZE(scache%cache),scache%nLines)
          WRITE(out_unit,"(i6,a,a,a)",advance='NO') i,'"',TRIM(scache%cache(i)),'"'
       END DO
    ELSE
       WRITE(out_unit,"(' cache: **NULL**')")
    END IF
    WRITE(out_unit,"('>')")
  END SUBROUTINE sc_writeout
!****************************************************************************

!!****f* cp_str_cache/sc_add_chars [1.0] *
!!
!!   NAME
!!     sc_add_chars
!!
!!   FUNCTION
!!     Add the given string to the cache.
!!     If needed the cache increases its size.
!!     
!!
!!   INPUTS
!!     - scache: the cache
!!     - str: the string to add to the cache
!!     - newline: if it is true (the default) a newline is assumed at the
!!       end of the string. 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_add_chars(scache,str,error,newline)
    TYPE(str_cache_type), POINTER            :: scache
    CHARACTER(len=*), INTENT(in)             :: str
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: newline

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_add_chars', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=baseStrLen), &
      DIMENSION(:), POINTER                  :: newCache
    INTEGER                                  :: ic_line, ic_shift, iline, &
                                                last_copied, lenStr, &
                                                shiftLines, stat, toAlloc
    INTEGER, DIMENSION(:), POINTER           :: newLinePos
    LOGICAL                                  :: failure, my_newline

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    my_newline=.TRUE.
    IF (PRESENT(newline)) my_newline=newline
    lenStr=LEN(str)
    IF (my_newline) lenStr=MAX(1,LEN_TRIM(str))
    IF (lenStr < sc_free(scache,error)) THEN
       IF (SIZE(scache%cache)*baseStrLen-scache%lastChar<lenStr) THEN
          CALL sc_compress(scache,error)
       END IF
    ELSE
       toAlloc=SIZE(scache%cache)+((lenStr+baseStrLen-1)/baseStrLen)*baseStrLen
       toAlloc=INT(toAlloc*1.2)
       ALLOCATE(newCache(toAlloc),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL sc_copy_buf(scache%cache,fromPos=scache%posAtt,&
            toPos=scache%lastChar,destBuf=newCache,error=error)
       DEALLOCATE(scache%cache,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scache%cache => newCache
       scache%lastChar=scache%lastChar-scache%posAtt+1
       scache%linePos=scache%linePos-scache%posAtt
       scache%posAtt=1
       shiftLines=0
       DO iline=1,scache%nLines
          IF (scache%linePos(iline)>0) EXIT
          shiftLines=iline
       END DO
       IF (shiftLines>0) THEN
          DO iline=1,scache%nLines-shiftLines
             scache%linePos(iline)=scache%linePos(iline+shiftLines)
          END DO
          scache%nLines=scache%nLines-shiftLines
       END IF
    END IF
    ic_line=(scache%lastChar)/baseStrLen
    ic_shift=MODULO(scache%lastChar,baseStrLen)
    last_copied=0
    IF (LEN(str)==0 .AND. my_newline) THEN
       scache%cache(ic_line)(ic_shift+1:ic_shift+1)=' '
    ELSE
       DO WHILE (last_copied < lenStr)
          scache%cache(ic_line)(ic_shift+1:MIN(ic_shift+lenStr,baseStrLen))=&
               str(1+last_copied:MIN(lenStr,last_copied+baseStrLen-ic_shift))
          last_copied=MIN(lenStr,last_copied+baseStrLen-ic_shift)
          ic_shift=0
       END DO
    END IF
    scache%lastChar=scache%lastChar+lenStr
    IF (my_newline) THEN
       scache%nLines=scache%nLines+1
       IF (SIZE(scache%linePos)<scache%nLines) THEN
          ALLOCATE(newLinePos(INT(CEILING(scache%nLines*1.2))),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          newLinePos(:SIZE(scache%linePos))=scache%linePos
          DEALLOCATE(scache%linePos,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          scache%linePos => newLinePos
       END IF
       scache%linePos(scache%nLines)=scache%lastChar
    END IF
  END SUBROUTINE sc_add_chars
!****************************************************************************

!!****f* cp_str_cache/sc_copy_buf [1.0] *
!!
!!   NAME
!!     sc_copy_buf
!!
!!   FUNCTION
!!     utility routine that copies between cache buffers
!!     works also if sourceBuf==destBuf
!!
!!   INPUTS
!!     - sourceBuf: buffer to copy from
!!     - fromPos: start position to copy from (inclusive bound)
!!     - toPos: end character to copy (inclusive bound)
!!     - destBuf: buffer to copy to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_copy_buf(sourceBuf,fromPos,toPos,destBuf,error)
    CHARACTER(len=baseStrLen), &
      DIMENSION(:), POINTER                  :: sourceBuf
    INTEGER, INTENT(in)                      :: fromPos, toPos
    CHARACTER(len=baseStrLen), &
      DIMENSION(:), POINTER                  :: destBuf
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_copy_buf', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*baseStrLen)              :: buf
    INTEGER                                  :: bufEnd, bufShift, ic2_line, &
                                                ic_endline, ic_line, ic_shift
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(sourceBuf),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(destBuf),cp_failure_level,routineP,error,failure)
    CPPrecondition(fromPos>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(destBuf)*baseStrLen>toPos-fromPos,cp_failure_level,routineP,error,failure)
    IF (fromPos>toPos) RETURN
    ic_line=(fromPos-1)/baseStrLen
    ic_shift=MODULO(fromPos-1,baseStrLen)
    ic_endline=(toPos-1)/baseStrLen
    bufShift=0
    bufEnd=0
    buf(1:baseStrLen-ic_shift)=sourceBuf(ic_line)(ic_shift+1:baseStrLen)
    bufEnd=bufEnd+baseStrLen-ic_shift
    ic_line=ic_line+1
    ic2_line=1
    DO WHILE (ic_line<ic_endline)
       buf(bufEnd+1:bufEnd+baseStrLen)=sourceBuf(ic_line)(1:baseStrLen)
       bufEnd=bufEnd+baseStrLen
       destBuf(ic2_line)=buf(1:baseStrLen)
       buf(1:baseStrLen-ic_shift)=buf(baseStrLen+1:2*baseStrLen-ic_shift)
       ic_line=ic_line+1
       ic2_line=ic2_line+1
    END DO
    destBuf(ic2_line)=buf(1:MODULO(toPos-1,baseStrLen)-bufShift+1)
  END SUBROUTINE sc_copy_buf
!****************************************************************************

!!****f* cp_str_cache/sc_compress [1.0] *
!!
!!   NAME
!!     sc_compress
!!
!!   FUNCTION
!!     Compresses the data to the beginning of the cache buffer (removing
!!     already read characters). Does not change the size of the buffer.
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_compress(scache,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_compress', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iline, shiftLines
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (scache%posAtt==1) RETURN
    CALL sc_copy_buf(scache%cache,fromPos=scache%posAtt,toPos=scache%lastChar,&
         destBuf=scache%cache,error=error)
    scache%lastChar=scache%lastChar-scache%posAtt+1
    scache%linePos=scache%linePos-scache%posAtt
    scache%posAtt=1
    shiftLines=0
    DO iline=1,scache%nLines
       IF (scache%linePos(iline)>0) EXIT
       shiftLines=iline
    END DO
    IF (shiftLines>0) THEN
       DO iline=1,scache%nLines-shiftLines
          scache%linePos(iline)=scache%linePos(iline+shiftLines)
       END DO
       scache%nLines=scache%nLines-shiftLines
    END IF
  END SUBROUTINE sc_compress
!****************************************************************************

!!****f* cp_str_cache/sc_get_chars [1.0] *
!!
!!   NAME
!!     sc_get_chars
!!
!!   FUNCTION
!!     returns a line (or part of it) from the cache to 
!!     outStr(fromChar:toChar). 
!!
!!   INPUTS
!!     - scache: the cache
!!     - outStr: the string where to return the characters
!!     - newLine: is set to true if a whole line was returned, to false 
!!       if only part of it was returned.
!!     - copied_chars: number of characters actually copied (excluding the 
!!       possible spaces that are inserted to get to toChar)
!!     - fromChar: starting char to where to copy, defaults to 1
!!     - toChar: upper bound (inclusive) to where to put the chars,
!!       defaults to len(outStr)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_get_chars(scache, outStr, error, newLine, copied_chars, fromChar, toChar)
    TYPE(str_cache_type), POINTER            :: scache
    CHARACTER(len=*), INTENT(out)            :: outStr
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(out), OPTIONAL           :: newLine
    INTEGER, INTENT(out), OPTIONAL           :: copied_chars
    INTEGER, INTENT(in), OPTIONAL            :: fromChar, toChar

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_get_chars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: copy_now, copy_tot, ic_end, &
                                                ic_line, ic_shift, iline, &
                                                my_from, my_to, nextLine
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    nextLine=HUGE(0)
    DO iline=1,scache%nLines
       IF (scache%linePos(iline)>=scache%posAtt) THEN
          nextLine=scache%linePos(iline)
          EXIT
       END IF
    END DO
    ic_end=MIN(nextLine,scache%lastChar)

    my_from=1
    IF (PRESENT(fromChar)) my_from=fromChar
    my_to=LEN(outStr)
    IF (PRESENT(toChar)) THEN
       my_to=toChar
       CPPrecondition(my_to<LEN(outStr),cp_failure_level,routineP,error,failure)
    END IF

    ic_shift=MODULO(scache%posAtt-1,baseStrLen)
    ic_line=(scache%posAtt-1)/baseStrLen
    copy_tot=0
    DO WHILE (my_from<=my_to .AND. scache%posAtt<=ic_end)
       copy_now=MIN(my_to-my_from+1,MIN(baseStrLen-ic_shift,nextLine-ic_end+1))
       outStr(my_from:my_from+copy_now-1)=scache%cache(ic_line)&
            (ic_shift+1:MIN(baseStrLen,ic_shift+copy_now))
       ic_shift=0
       ic_line=ic_line+1
       scache%posAtt=scache%posAtt+copy_now
       my_from=my_from+copy_now
       copy_tot=copy_tot+copy_now
    END DO
    IF (PRESENT(newLine)) THEN
       newLine=scache%posAtt>nextLine
    END IF
    IF (PRESENT(copied_chars)) copied_chars=copy_tot
    IF (my_from<=my_to) THEN
       outStr(my_from:my_to)=' '
    END IF
  END SUBROUTINE sc_get_chars
!****************************************************************************

!!****f* cp_str_cache/sc_free2 [1.0] *
!!
!!   NAME
!!     sc_free2
!!
!!   FUNCTION
!!     Returns the number of free chars in the cache.
!!     If the cache has zero lines returns still one line, to always 
!!     guarantee a non zero capacity.
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  FUNCTION sc_free2(scache,error) RESULT(res)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_free2', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    res=MAX(1,SIZE(scache%cache))*baseStrLen-scache%lastChar+scache%posAtt
  END FUNCTION sc_free2
!****************************************************************************

!!****f* cp_str_cache/sc_free2 [1.0] *
!!
!!   NAME
!!     sc_free2
!!
!!   FUNCTION
!!     Returns the number of free chars in the cache.
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  FUNCTION sc_free(scache,error) RESULT(res)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_free', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    res=SIZE(scache%cache)*baseStrLen-scache%lastChar+scache%posAtt
  END FUNCTION sc_free
!****************************************************************************

!!****f* cp_str_cache/sc_used [1.0] *
!!
!!   NAME
!!     sc_used
!!
!!   FUNCTION
!!     Returns the number of chars stored in the cache
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  FUNCTION sc_used(scache,error) RESULT(res)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_used', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    res=scache%lastChar-scache%posAtt+1
  END FUNCTION sc_used
!****************************************************************************

!!****f* cp_str_cache/sc_trim_cache [1.0] *
!!
!!   NAME
!!     sc_trim_cache
!!
!!   FUNCTION
!!     Reduces cache size to the minimum size possible without discarding
!!     recorded data.
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_trim_cache(scache,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_trim_cache', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=baseStrLen), &
      DIMENSION(:), POINTER                  :: newCache
    INTEGER                                  :: minLines, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)    
    CALL sc_compress(scache,error=error)
    minLines=(scache%lastChar+baseStrLen-1)/baseStrLen
    IF (SIZE(scache%cache)>minLines) THEN
       ALLOCATE(newCache(minLines),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       newCache(:)=scache%cache(1:minLines)
       DEALLOCATE(scache%cache,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scache%cache => newCache
    END IF
  END SUBROUTINE sc_trim_cache
!****************************************************************************

!!****f* cp_str_cache/sc_clear_cache [1.0] *
!!
!!   NAME
!!     sc_clear_cache
!!
!!   FUNCTION
!!     Clears the cache (but keeps its size constant)
!!
!!   INPUTS
!!     - scache: the cache
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_clear_cache(scache,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_clear_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    DO i=1,SIZE(scache%cache)
       scache%cache(i)=' '
    END DO
    scache%posAtt=1
    scache%lastChar=0
    scache%nLines=0
    scache%linePos=0
  END SUBROUTINE sc_clear_cache
!****************************************************************************

!!****f* cp_str_cache/sc_synchronize_cache [1.0] *
!!
!!   NAME
!!     sc_synchronize_cache
!!
!!   FUNCTION
!!     broadcasts the cache of para_env%source to all the processors
!!
!!   INPUTS
!!     - scache: the cache
!!     - para_env: the parallel environment to do the broadcast
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_synchronize_cache(scache,para_env,error)
    TYPE(str_cache_type), POINTER            :: scache
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_synchronize_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: copy_now, i, ibuf, ic, &
                                                ic_endline, ic_line, ipos, &
                                                stat
    INTEGER, DIMENSION(10*baseStrLen)        :: buf
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (para_env%num_pe==1) RETURN
    IF (para_env%source==para_env%mepos) THEN
       buf(1)=scache%nLines
       buf(2)=scache%posAtt
       buf(3)=scache%lastChar
       buf(4)=SIZE(scache%cache)
       buf(5)=SIZE(scache%linePos)
    END IF
    CALL mp_bcast(buf(1:5),para_env%source,para_env%group)
    scache%nLines              =buf(1)
    scache%posAtt              =buf(2)
    scache%lastChar            =buf(3)
    IF (SIZE(scache%cache) /= buf(4)) THEN
       DEALLOCATE(scache%cache,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scache%cache(buf(4)),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (SIZE(scache%linePos) /= buf(5)) THEN
       DEALLOCATE(scache%linePos,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scache%linePos(buf(5)),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF

    ipos=1
    DO WHILE (ipos<=scache%nLines)
       copy_now=MIN(SIZE(buf),scache%nLines-ipos+1)
       IF (para_env%source==para_env%mepos) THEN
          buf(1:copy_now)=scache%linePos(ipos:ipos+copy_now-1)
       END IF
       CALL mp_bcast(buf(1:copy_now),para_env%source,para_env%group)
       IF (para_env%source/=para_env%mepos) THEN
          scache%linePos(ipos:ipos+copy_now-1)=buf(1:copy_now)
       END IF
       ipos=ipos+copy_now
    END DO

    ipos=1
    ic_line=1
    ic_endline=(scache%lastChar-1)/baseStrLen+1
    DO WHILE (ic_line<=ic_endline)
       copy_now=MIN(SIZE(buf)/baseStrLen,ic_endline-ic_line+1)
       IF (para_env%source==para_env%mepos) THEN
          DO i=1,copy_now
             ibuf=(i-1)*baseStrLen
             DO ic=1,baseStrLen
                buf(ibuf+ic)=ICHAR(scache%cache(i)(ic:ic))
             END DO
          END DO
       END IF
       CALL mp_bcast(buf(1:copy_now*baseStrLen),para_env%source,para_env%group)
       IF (para_env%source/=para_env%mepos) THEN
          DO i=1,copy_now
             ibuf=(i-1)*baseStrLen
             DO ic=1,baseStrLen
                scache%cache(i)(ic:ic)=CHAR(buf(ibuf+ic))
             END DO
          END DO
       END IF
       ic_line=ic_line+copy_now
    END DO
  END SUBROUTINE sc_synchronize_cache
!****************************************************************************

!!****f* cp_str_cache/sc_fill_from_unit [1.0] *
!!
!!   NAME
!!     sc_fill_from_unit
!!
!!   FUNCTION
!!     fills the cache from the given unit
!!
!!   INPUTS
!!     - scache: the cache
!!     - unit_nr: the unit to read from
!!     - at_eof: is set to true if the end of file was encountred.
!!       the cache might (probably will) still contain data.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
  SUBROUTINE sc_fill_from_unit(scache,unit_nr,error,at_eof)
    TYPE(str_cache_type), POINTER            :: scache
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(out), OPTIONAL           :: at_eof

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_fill_from_unit', &
      routineP = moduleN//':'//routineN

    CHARACTER(9)                             :: fmtStr
    CHARACTER(baseStrLen)                    :: buf
    INTEGER                                  :: iostat, readSize
    LOGICAL                                  :: failure, newline

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL sc_compress(scache,error)
    readSize=MIN(baseStrLen,sc_free2(scache,error))
    DO WHILE (readSize>0)
       newline=.FALSE.    
       WRITE(fmtStr,"('(a',i6,')')") readSize
       READ(buf,fmtStr,IOSTAT=iostat) buf
       IF (iostat/=0) THEN
          CPAssert(iostat<0,cp_failure_level,routineP,error,failure)
          IF (iostat==get_eof_stat()) THEN
             PRINT *, 'EOF, buf "',TRIM(buf),'"'
             IF (PRESENT(at_eof)) at_eof=.TRUE.
             EXIT
          ELSE
             newline=.TRUE.
          END IF
       END IF
       IF (failure) EXIT
       PRINT *, 'buf "',TRIM(buf),'" newline',newline
       CALL sc_add_chars(scache,buf(1:readSize),error=error,newline=newline)
    END DO
  END SUBROUTINE sc_fill_from_unit
!****************************************************************************

!!****f* cp_str_cache/sc_guarantee_capacity [1.0] *
!!
!!   NAME
!!     sc_guarantee_capacity
!!
!!   FUNCTION
!!     guarantees that the whole cache can hold at least nchar
!!     (but it might be already full)
!!
!!   INPUTS
!!     - scache: the cache
!!     - nchar: the requested minimum capacity
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** *********************************************************************
! 
  SUBROUTINE sc_guarantee_capacity(scache,nchar,error)
    TYPE(str_cache_type), POINTER            :: scache
    INTEGER, INTENT(in)                      :: nchar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_guarantee_capacity', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=baseStrLen), &
      DIMENSION(:), POINTER                  :: newCache
    INTEGER                                  :: nlines, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(scache),cp_failure_level,routineP,error,failure)
    CPPrecondition(scache%ref_count>0,cp_failure_level,routineP,error,failure)
    nlines=(nchar+baseStrLen-1)/baseStrLen
    IF (nlines>SIZE(scache%cache)) THEN
       ALLOCATE(newCache(nlines),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL sc_compress(scache,error=error)
       CALL sc_copy_buf(scache%cache,fromPos=scache%posAtt,toPos=scache%lastChar,&
            destBuf=newCache,error=error)
       DEALLOCATE(scache%cache,stat=stat)
       scache%cache => newCache
    END IF
  END SUBROUTINE sc_guarantee_capacity
!****************************************************************************

END MODULE cp_str_cache
