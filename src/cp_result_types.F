!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author fschiff (12.2007)
! *****************************************************************************
  MODULE cp_result_types

  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_result_types'

! *** Public data types ***

  PUBLIC :: cp_result_type,&
            cp_result_p_type

! *** Public subroutines ***

  PUBLIC :: cp_result_create,&
            cp_result_release,&
            cp_result_retain,&
            put_results,&
            get_results,&
            erase_result,&
            transfer_result_type

  INTEGER, PARAMETER, PUBLIC           ::  description_length=30

  INTERFACE put_results
     MODULE PROCEDURE put_result_r1, put_result_r2,replace_full_type
  END INTERFACE

  INTERFACE get_results
     MODULE PROCEDURE get_result_r1, get_result_r2,get_nreps
  END INTERFACE
  
! *****************************************************************************
!> \brief contains arbitrary information which need to be stored
!> \note
!>      result_list is a character list, in which everthing can be stored
!>      before passing any variable just name the variable like '[NAME]'
!>      brackets will be used to identify the start of a new set 
!> \author fschiff (12.2007)
! *****************************************************************************
  TYPE cp_result_type
     CHARACTER(LEN=description_length),DIMENSION(:),&
          POINTER                                         :: result_list
     INTEGER,DIMENSION(:),POINTER                         :: desc_map
     INTEGER :: ref_count
  END TYPE cp_result_type

! *****************************************************************************
  TYPE cp_result_p_type
     TYPE(cp_result_type), POINTER                        :: results
  END TYPE cp_result_p_type

CONTAINS

! *****************************************************************************
!> \brief allocates and intitializes the cp_result
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_result_create(results,error)
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I",'',handle)
    failure=.FALSE.

    ALLOCATE(results, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       results%ref_count=1
       NULLIFY(results%result_list,results%desc_map)
       ALLOCATE(results%result_list(2), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       results%result_list(1)="RESULTS"
       results%result_list(2)='END'   
       ALLOCATE(results%desc_map(2), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
       results%desc_map(1)=1
       results%desc_map(2)=2
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_result_create

! *****************************************************************************
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_result_release(results,error)
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I",'',handle)    
    failure=.FALSE.
    IF(ASSOCIATED(results))THEN
       CPPrecondition(results%ref_count>0,cp_failure_level,routineP,error,failure)
       results%ref_count=results%ref_count-1
       IF (results%ref_count==0) THEN
          IF(ASSOCIATED(results%result_list))THEN
             DEALLOCATE(results%result_list,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
          END IF
          IF(ASSOCIATED(results%desc_map))THEN
             DEALLOCATE(results%desc_map,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
          END IF
          
          DEALLOCATE(results,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    NULLIFY(results)
  END SUBROUTINE cp_result_release

! *****************************************************************************
  SUBROUTINE cp_result_retain(results,error)
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_result_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(results%ref_count>0,cp_failure_level,routineP,error,failure)
       results%ref_count=results%ref_count+1
    END IF
  END SUBROUTINE cp_result_retain

! *****************************************************************************
!> \brief writes an arbitrary result as a string in result_list
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE put_result_r1(results,description,values,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length)        :: description
    REAL(KIND=dp), DIMENSION(:)              :: values
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_result_r1', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length)        :: line
    CHARACTER(LEN=description_length), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_string
    INTEGER                                  :: i, isize, j, jsize, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_map
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    isize=SIZE(results%result_list)
    jsize=SIZE(values)

    ALLOCATE(tmp_string(isize+jsize+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp_map(SIZE(results%desc_map)+1))
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !!!! don't transfer the [END] statement
    DO i=1,isize-1
       tmp_string(i)=results%result_list(i)
    END DO
    DO i=1,SIZE(results%desc_map)
       tmp_map(i)=results%desc_map(i)
    END DO
    tmp_map(SIZE(results%desc_map)+1)=isize+jsize+1
    
    DEALLOCATE(results%result_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(results%desc_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(results%result_list(isize+jsize+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(results%desc_map(SIZE(tmp_map)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    tmp_string(isize)=description
    results%desc_map=tmp_map

    DO j=1,jsize
       i=isize+j
       WRITE(line,*)values(j)
       tmp_string(i)=line
    END DO
    tmp_string(isize+jsize+1)='END'
    results%result_list=tmp_string
    DEALLOCATE(tmp_string,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
  END SUBROUTINE put_result_r1
  
! *****************************************************************************
  SUBROUTINE put_result_r2(results,description,values,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length)        :: description
    REAL(KIND=dp), DIMENSION(:, :)           :: values
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_result_r2', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length)        :: line
    CHARACTER(LEN=description_length), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_string
    INTEGER                                  :: i, isize, j, jsize, k, ksize, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_map
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    isize=SIZE(results%result_list)
    jsize=SIZE(values,1)
    ksize=SIZE(values,2)

    ALLOCATE(tmp_string(isize+(jsize*ksize)+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tmp_map(SIZE(results%desc_map)+1))
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !!!! don't transfer the [END] statement
    DO i=1,isize-1
       tmp_string(i)=results%result_list(i)
    END DO
    DO i=1,SIZE(results%desc_map)
       tmp_map(i)=results%desc_map(i)
    END DO
    tmp_map(SIZE(results%desc_map)+1)=isize+jsize*ksize+1

    DEALLOCATE(results%result_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(results%desc_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(results%result_list(isize+(jsize*ksize)+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(results%desc_map(SIZE(tmp_map)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    tmp_string(isize+1)=description
    results%desc_map=tmp_map
    
    DO j=1,jsize
       DO k=1,ksize
          i=isize+(j-1)*jsize+ksize
          tmp_string(i)=line
       END DO
    END DO

    tmp_string(isize+(jsize*ksize)+1)='END'
    results%result_list=tmp_string
    DEALLOCATE(tmp_string,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
  END SUBROUTINE put_result_r2

! *****************************************************************************
  SUBROUTINE replace_full_type(results,list,map,size_list,size_map,displ_list,displ_map,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length), &
      DIMENSION(:)                           :: list
    INTEGER, DIMENSION(:)                    :: map
    INTEGER, OPTIONAL                        :: size_list, size_map, &
                                                displ_list, displ_map
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'replace_full_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    
    !!! displacement (displ_map and displ_list) have a start value of 0 for the first element !!! 

    DEALLOCATE(results%result_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(results%desc_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(PRESENT(size_list))THEN
       CPPostcondition(PRESENT(size_map),cp_failure_level,routineP,error,failure) 
       CPPostcondition(PRESENT(displ_list),cp_failure_level,routineP,error,failure) 
       CPPostcondition(PRESENT(displ_map),cp_failure_level,routineP,error,failure) 
       ALLOCATE(results%result_list(size_list),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%desc_map(size_map),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       results%desc_map(1:size_map)=map(displ_map+1:displ_map+size_map)       
       results%result_list(1:size_list)=list(displ_list+1:displ_list+size_list)

    ELSE
       ALLOCATE(results%result_list(SIZE(list)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%desc_map(SIZE(map)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       results%desc_map(:)=map(:)
       results%result_list(:)=list(:)
    END IF
  END SUBROUTINE replace_full_type

! *****************************************************************************
!> \brief gets the required part out of the result_list
!> \param nval :  if more than one entry for a given description is given you may choose 
!>              which entry you want
!> \param n_rep : integer indicating how many times the section exists in result_list
!> \param n_entries : gets the number of lines used for a given description
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE get_result_r1(results,description,values,nval,n_rep,n_entries,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length)        :: description
    REAL(KIND=dp), DIMENSION(:)              :: values
    INTEGER, OPTIONAL                        :: nval, n_rep, n_entries
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_result_r1', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length)        :: line
    INTEGER                                  :: desc_size, i, k, l, m, nrep, &
                                                size_res
    LOGICAL                                  :: failure

    desc_size=SIZE(results%desc_map)
    nrep=0
    DO i=1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))nrep=nrep+1
    END DO

    IF(PRESENT(n_rep))THEN
       n_rep=nrep
    END IF

    CPPostcondition(nrep.GT.0,cp_failure_level,routineP,error,failure)  
    DO i=1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
          size_res=results%desc_map(i+1)-results%desc_map(i)-1
          EXIT
       END IF
    END DO
    IF(PRESENT(n_entries))n_entries=size_res
    IF(PRESENT(nval))THEN
       CPPostcondition(size_res.EQ.SIZE(values),cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(nrep*size_res.EQ.SIZE(values),cp_failure_level,routineP,error,failure)
    END IF
    k=0
    DO i = 1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
          IF(PRESENT(nval))THEN
             IF(k+1==nval)THEN
                DO l=1,size_res
                   m=results%desc_map(i)+l
                   line=results%result_list(m)
                   READ(line,*)values(l+k*size_res)
                END DO
                EXIT
             ELSE
                k=k+1
             END IF
          END IF
          DO l=1,size_res
             m=results%desc_map(i)+l
             line=results%result_list(m)
             READ(line,*)values(l+k*size_res)
          END DO
          k=k+1
       END IF
    END DO

  END SUBROUTINE get_result_r1

! *****************************************************************************
  SUBROUTINE get_result_r2(results,description,values,nval,n_rep,n_entries,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length)        :: description
    REAL(KIND=dp), DIMENSION(:, :)           :: values
    INTEGER, OPTIONAL                        :: nval, n_rep, n_entries
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_result_r2', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length)        :: line
    INTEGER                                  :: desc_size, i, k, l, m, n, &
                                                nrep, size_res, val_size1, &
                                                val_size2
    LOGICAL                                  :: failure

    desc_size=SIZE(results%desc_map)
    nrep=0
    DO i=1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))nrep=nrep+1
    END DO

    IF(PRESENT(n_rep)) n_rep=nrep

    CPPostcondition(nrep.GT.0,cp_failure_level,routineP,error,failure)  
    DO i=1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
          size_res=results%desc_map(i+1)-results%desc_map(i)-1
          EXIT
       END IF
    END DO
    IF(PRESENT(n_entries))n_entries=size_res
    IF(PRESENT(nval))THEN
       CPPostcondition(size_res.EQ.SIZE(values),cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(nrep*size_res.EQ.SIZE(values),cp_failure_level,routineP,error,failure)
    END IF
    val_size2=SIZE(values,2)
    val_size1=size_res/val_size2
    k=0
    DO i = 1,desc_size
       IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
          IF(PRESENT(nval))THEN
             IF(k+1==nval)THEN
                DO l=1,val_size1
                   DO n=1,val_size2
                      m=results%desc_map(i)+(l-1)*val_size2+n
                      line=results%result_list(m)
                      READ(line,*)values(l+k*size_res,n)
                   END DO
                END DO
                EXIT
             ELSE
                k=k+1
             END IF
             DO l=1,val_size1
                DO n=1,val_size2
                   m=results%desc_map(i)+(l-1)*val_size2+n
                   line=results%result_list(m)
                   READ(line,*)values(l+k*size_res,n)
                END DO
             END DO
             k=k+1
          END IF
       END IF
    END DO

  END SUBROUTINE get_result_r2

! *****************************************************************************
  SUBROUTINE get_nreps(results,description,n_rep,n_entries,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length)        :: description
    INTEGER, OPTIONAL                        :: n_rep, n_entries
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: desc_size, I

    desc_size=SIZE(results%desc_map)
    IF(PRESENT(n_rep))THEN
       n_rep=0
       DO i=1,desc_size
          IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))n_rep=n_rep+1
       END DO
    END IF
    IF(PRESENT(n_entries))THEN
       DO i=1,desc_size
          IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
             n_entries=results%desc_map(i+1)-results%desc_map(i)-1
             EXIT
          END IF
       END DO
    END IF
  END SUBROUTINE get_nreps

! *****************************************************************************
!> \brief erase a part of  result_list
!> \param nval : if more than one entry for a given description is given you may choose 
!>             which entry you want
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE erase_result(results,description,nval,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=description_length), &
      OPTIONAL                               :: description
    INTEGER, OPTIONAL                        :: nval
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'erase_result', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_string
    INTEGER                                  :: count, desc_size, i, isize, &
                                                j, k, n, nrep, size_res, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_map
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(PRESENT(description))THEN
       CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
       desc_size=SIZE(results%desc_map)
       isize=SIZE(results%result_list)
       nrep=0
       DO i=1,desc_size
          IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))nrep=nrep+1
       END DO
       IF(nrep.NE.0)THEN

          DO i=1,desc_size
             IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
                size_res=results%desc_map(i+1)-results%desc_map(i)-1
                EXIT
             END IF
          END DO

          IF(PRESENT(nval))THEN
             CPPostcondition(nval.LE.nrep,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp_string(isize-(size_res+1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp_map(desc_size-1),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
          ELSE
             ALLOCATE(tmp_string(isize-nrep*(size_res+1)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(tmp_map(desc_size-nrep),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF

          j=0
          k=0
          count=0
          DO i=1,desc_size-1
             IF(TRIM(results%result_list(results%desc_map(i)))==TRIM(description))THEN
                IF(PRESENT(nval))THEN
                   IF(k+1==nval)THEN
                      k=k+1
                   ELSE
                      count=count+1
                      tmp_map(count)= results%desc_map(i)
                      DO n=1,results%desc_map(i+1)-results%desc_map(i)
                         j=j+1
                         tmp_string(j)= results%result_list( results%desc_map(i)+n-1)
                      END DO
                      k=k+1
                   END IF
                END IF
             ELSE
                count=count+1
                tmp_map(count)= j+1     
                DO n=1,results%desc_map(i+1)-results%desc_map(i)
                   j=j+1
                   tmp_string(j)= results%result_list( results%desc_map(i)+n-1)
                END DO
             END IF
          END DO
          j=j+1
          tmp_string(j)= results%result_list( results%desc_map(desc_size))   
          tmp_map(SIZE(tmp_map))=j
          DEALLOCATE(results%result_list,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(results%desc_map,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(results%result_list(SIZE(tmp_string)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(results%desc_map(SIZE(tmp_map)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          results%result_list=tmp_string
          results%desc_map=tmp_map

          DEALLOCATE(tmp_string,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(tmp_map,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       END IF
    ELSE
       DEALLOCATE(results%result_list,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(results%desc_map,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%result_list(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%desc_map(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
       results%result_list(1)="RESULTS"
       results%result_list(2)='END'   
       results%desc_map(1)=1
       results%desc_map(2)=2
    END IF
  END SUBROUTINE erase_result

! *****************************************************************************
!> \brief erase a part of  result_list
!> \param res_target :           type to which you want to transfer
!> \param res_source :           type from which you want to transfer
!>      nval(OPTIONAL)::       if more than one entry with description_source 
!>                             is present in res_source you may chose which one to transfer
!>                             if not present all values with this description are transfered
!> \param description_target :   new description in res_target
!> \param description_source :   description of the values you want to transfer
!> \par History
!>      12.2007 created
!> \author fschiff
! *****************************************************************************
  SUBROUTINE transfer_result_type(res_target,res_source,description_target,description_source,nval,error)
    TYPE(cp_result_type), POINTER            :: res_target, res_source
    CHARACTER(LEN=description_length), &
      OPTIONAL                               :: description_target, &
                                                description_source
    INTEGER, OPTIONAL                        :: nval
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transfer_result_type', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=description_length), &
      ALLOCATABLE, DIMENSION(:)              :: tmp_string
    INTEGER                                  :: desc_size, desc_size1, &
                                                handle, i, j, k, n, nrep, &
                                                size_res, size_source, &
                                                size_target, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_map
    LOGICAL                                  :: failure

    CALL timeset(routineN,'I','',handle)
    size_target=SIZE(res_target%result_list)
    size_source=SIZE(res_source%result_list)
    desc_size=SIZE(res_source%desc_map)
    desc_size1=SIZE(res_target%desc_map)
    IF(PRESENT(description_source))THEN
       nrep=0
       DO i=1,desc_size
          IF(TRIM(res_source%result_list(res_source%desc_map(i)))==TRIM(description_source))nrep=nrep+1
       END DO

       CPPostcondition(nrep.GT.0,cp_failure_level,routineP,error,failure)  

       DO i=1,desc_size
          IF(TRIM(res_source%result_list(res_source%desc_map(i)))==(description_source))THEN
             size_res=res_source%desc_map(i+1)-res_source%desc_map(i)-1
             EXIT
          END IF
       END DO

       IF(PRESENT(nval))THEN
          k=size_source+size_target+(size_res+1)-2
          ALLOCATE(tmp_string(k))
          ALLOCATE(tmp_map(desc_size+1))
       ELSE
          k=size_source+size_target+nrep*(size_res+1)-2 
          ALLOCATE(tmp_string(k))
          ALLOCATE(tmp_map(desc_size+nrep))
       END IF

       DO i=1,desc_size1-1
          tmp_map=res_target%desc_map(i)
       END DO
       DO i=1,size_target-1       
          tmp_string(i)=res_target%result_list(i)
       END DO
       k=0
       j=size_target-1
       DO i=1,desc_size-1
          IF(res_source%result_list(res_source%desc_map(i))==(description_source))THEN
             IF(PRESENT(nval))THEN
                IF(k+1==nval)THEN

                   tmp_map(desc_size1+k)=j+1
                   DO n=1,res_source%desc_map(i+1)-res_source%desc_map(i)
                      j=j+1
                      tmp_string(j)= res_source%result_list(res_source%desc_map(i)+n-1)
                   END DO
                   tmp_string(tmp_map(desc_size1+k))= description_target            
                ELSE
                   k=k+1
                END IF
             ELSE
                tmp_map(desc_size1+k)=j+1
                DO n=1,res_source%desc_map(i+1)-res_source%desc_map(i)
                   j=j+1
                   tmp_string(j)= res_source%result_list( res_source%desc_map(i)+n-1)
                END DO
                tmp_string(tmp_map(desc_size1+k))= description_target  
                k=k+1
             END IF
          END IF
       END DO
       tmp_string(j)= res_source%result_list( res_source%desc_map(desc_size))   
       tmp_map(SIZE(tmp_map))=j
    ELSE

       ALLOCATE(tmp_string(size_target+size_source-2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(tmp_map(desc_size+desc_size1-2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO i = 1, size_target-1
          tmp_string(i)=res_target%result_list(i)
       END DO
       DO i=1,desc_size1-1
          tmp_map(i)=res_target%desc_map(i)
       END DO
       DO i=1,desc_size-1
          tmp_map(i+desc_size1-1)=res_source%desc_map(i+1)-1&
               +res_target%desc_map(desc_size1-1)
       END DO
       DO i = 2, size_source
          tmp_string(i+size_target-2)=res_source%result_list(i)
       END DO
    END IF
    DEALLOCATE(res_target%result_list,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(res_target%desc_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(res_target%result_list(SIZE(tmp_string)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(res_target%desc_map(SIZE(tmp_map)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
    ! workaround for a NAG compiler bug that appears randomly from time
    ! to time..
    ! this fix may become obsolete in the future.
    DO i = 1, SIZE(tmp_string)
       res_target%result_list(i) = tmp_string(i)
    END DO
    DO i = 1, SIZE(tmp_map)
       res_target%desc_map(i)    = tmp_map(i)
    END DO

    DEALLOCATE(tmp_string,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_map,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE transfer_result_type
END MODULE cp_result_types
