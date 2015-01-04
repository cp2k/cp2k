!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief stores a lists of integer that are local to a processor.
!>      The idea is that these integers represent objects that are distributed
!>      between the different processors.
!>      The ordering is just to make some operation more efficent, logically
!>      these lists are like sets.
!>      Some operations assume that the integers represent a range of values
!>      from 1 to a (not too big) maxval, and that an element is present just
!>      once, and only on a processor (these assumption are marked in the
!>      documentation of such operations).
!>      The concrete task for which this structure was developed was
!>      distributing atoms between the processors.
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE distribution_1d_types

  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE parallel_rng_types,              ONLY: delete_rng_stream,&
                                             rng_stream_p_type
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'distribution_1d_types'
  INTEGER, SAVE, PRIVATE :: distribution_1d_last_id_nr=0

  PUBLIC :: distribution_1d_type
  PUBLIC :: distribution_1d_create, distribution_1d_retain, distribution_1d_release

! *****************************************************************************
  TYPE local_particle_type
    TYPE(rng_stream_p_type), DIMENSION(:), POINTER :: rng
  END TYPE local_particle_type
!***

! *****************************************************************************
!> \brief structure to store local (to a processor) ordered lists of integers.
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param n_el n_el(i) is number of elements of list(i)
!> \param list list(i) contains an ordered list of integer (the array
!>        might be bigger than n_el(i), but the extra elements should be
!>        ignored)
!> \param para_env the parallel environment underlying the distribution
!> \param listbased_distribution true if each list has its own
!>        distribution
!> \par History
!>      06.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE distribution_1d_type
     INTEGER :: ref_count, id_nr
     LOGICAL :: listbased_distribution
     INTEGER, DIMENSION(:), POINTER :: n_el
     TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER :: list
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(local_particle_type), DIMENSION(:), POINTER :: local_particle_set
  END TYPE distribution_1d_type

! *****************************************************************************
!> \brief pointer to a distribution_1d_type
!> \param distribution_1d: the pointer to the distribution_1d
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE distribution_1d_p_type
     TYPE(distribution_1d_type), POINTER :: distribution_1d
  END TYPE distribution_1d_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a local list
!> \param distribution_1d the lists to create
!> \param para_env parallel environement to be used
!> \param listbased_distribution if each list has its own distribution
!>        (defaults to false)
!> \param n_el number of elements in each list (defaults to 0)
!> \param n_lists number of lists to create (defaults to 1, or size(n_el))
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE distribution_1d_create(distribution_1d,para_env,listbased_distribution,&
     n_el,n_lists,error)
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in), OPTIONAL            :: listbased_distribution
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: n_el
    INTEGER, INTENT(in), OPTIONAL            :: n_lists
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_1d_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, my_n_lists, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  my_n_lists=1
  IF(PRESENT(n_el)) my_n_lists=SIZE(n_el)
  IF(PRESENT(n_lists)) my_n_lists=n_lists

  ALLOCATE(distribution_1d,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     distribution_1d_last_id_nr=distribution_1d_last_id_nr+1
     distribution_1d%id_nr=distribution_1d_last_id_nr

     distribution_1d%ref_count=1

     distribution_1d%para_env => para_env
     CALL cp_para_env_retain(para_env,error=error)

     NULLIFY(distribution_1d%list,distribution_1d%n_el)

     distribution_1d%listbased_distribution=.FALSE.
     IF (PRESENT(listbased_distribution)) &
          distribution_1d%listbased_distribution=listbased_distribution

     ALLOCATE(distribution_1d%n_el(my_n_lists), distribution_1d%list(my_n_lists),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (PRESENT(n_el)) THEN
        distribution_1d%n_el(1:my_n_lists)=n_el(1:my_n_lists)
        DO ilist=1,my_n_lists
          ALLOCATE(distribution_1d%list(ilist)%array(n_el(ilist)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          distribution_1d%list(ilist)%array=-1
        END DO
     ELSE
        distribution_1d%n_el(1:my_n_lists)=0
        DO ilist=1,my_n_lists
          NULLIFY(distribution_1d%list(ilist)%array)
        END DO
     END IF

     NULLIFY (distribution_1d%local_particle_set)

  END IF
END SUBROUTINE distribution_1d_create

! *****************************************************************************
!> \brief retains a distribution_1d
!> \param distribution_1d  the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE distribution_1d_retain(distribution_1d,error)
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_1d_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(distribution_1d),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(distribution_1d%ref_count>0,cp_failure_level,routineP,error)
     distribution_1d%ref_count=distribution_1d%ref_count+1
  END IF
END SUBROUTINE distribution_1d_retain

! *****************************************************************************
!> \brief releases the given distribution_1d
!> \param distribution_1d the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE distribution_1d_release(distribution_1d,error)
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_1d_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilist, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local, stat
    LOGICAL                                  :: failure
    TYPE(local_particle_type), &
      DIMENSION(:), POINTER                  :: local_particle_set

  failure=.FALSE.

  IF (ASSOCIATED(distribution_1d)) THEN
     CPPreconditionNoFail(distribution_1d%ref_count>0,cp_failure_level,routineP,error)
     distribution_1d%ref_count=distribution_1d%ref_count-1
     IF (distribution_1d%ref_count==0) THEN
        DEALLOCATE(distribution_1d%n_el,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DO ilist=1,SIZE(distribution_1d%list)
           DEALLOCATE(distribution_1d%list(ilist)%array,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END DO
        DEALLOCATE(distribution_1d%list,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        !MK Delete Wiener process

        local_particle_set => distribution_1d%local_particle_set

        IF (ASSOCIATED(local_particle_set)) THEN
          nparticle_kind = SIZE(local_particle_set)
          DO iparticle_kind=1,nparticle_kind
            IF (ASSOCIATED(local_particle_set(iparticle_kind)%rng)) THEN
              nparticle_local = SIZE(local_particle_set(iparticle_kind)%rng)
              DO iparticle_local=1,nparticle_local
                IF (ASSOCIATED(local_particle_set(iparticle_kind)%&
                                rng(iparticle_local)%stream)) THEN
                  CALL delete_rng_stream(local_particle_set(iparticle_kind)%&
                                          rng(iparticle_local)%stream,error=error)
                END IF
              END DO
              DEALLOCATE (local_particle_set(iparticle_kind)%rng,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            END IF
          END DO
          DEALLOCATE(local_particle_set,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        CALL cp_para_env_release(distribution_1d%para_env,error=error)

        DEALLOCATE(distribution_1d,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF

END SUBROUTINE distribution_1d_release

! *****************************************************************************
!> \brief writes a description of the local distribution to the given unit
!> \param distribution_1d the list to describe 
!> \param unit_nr the unit to write to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE distribution_1d_write(distribution_1d, unit_nr, error)
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_1d_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ilist
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(distribution_1d)) THEN
     WRITE(unit=unit_nr, fmt="(' <distribution_1d>{ id_nr=',i8,',ref_count=',i8,',')")&
          distribution_1d%id_nr, distribution_1d%ref_count
     IF (ASSOCIATED(distribution_1d%n_el)) THEN
        WRITE(unit=unit_nr, fmt="(' n_el=(')", advance="no")
        DO i=1,SIZE(distribution_1d%n_el)
           WRITE(unit=unit_nr, fmt="(i6,',')", advance="no") distribution_1d%n_el(i)
        END DO
        WRITE(unit=unit_nr, fmt="('),')")
     ELSE
        WRITE(unit=unit_nr, fmt="(' n_el=*null*,')")
     END IF
     IF (ASSOCIATED(distribution_1d%list)) THEN
        WRITE(unit=unit_nr, fmt="(' list=(')")
        DO ilist=1,SIZE(distribution_1d%list)
           IF (ASSOCIATED(distribution_1d%list(ilist)%array)) THEN
              WRITE(unit=unit_nr, fmt="('  (')", advance="no")
              DO i=1,SIZE(distribution_1d%list(ilist)%array)
                 WRITE(unit=unit_nr, fmt="(i6,',')", advance="no")&
                      distribution_1d%list(ilist)%array(i)
              END DO
              WRITE(unit=unit_nr, fmt="('),')")
           ELSE
              WRITE(unit=unit_nr, fmt="('  *null*,')")
           END IF
        END DO
        WRITE(unit=unit_nr, fmt="(' )')")
     ELSE
        WRITE(unit=unit_nr, fmt="(' list=*null*,')")
     END IF
     WRITE(unit=unit_nr, fmt="(' }')")
  ELSE
     WRITE(unit=unit_nr, fmt="(' <distribution_1d *null*>,')")
  END IF
END SUBROUTINE distribution_1d_write

END MODULE distribution_1d_types
