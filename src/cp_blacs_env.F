!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_blacs_env [1.0] *
!!
!!   NAME
!!     cp_blacs_env
!!
!!   FUNCTION
!!     methods related to the blacs parallel environment
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!     02.2004 modified to associate a blacs_env with a given para_env
!!
!!   SOURCE
!****************************************************************************
MODULE cp_blacs_env
  USE cp_array_i_utils,                ONLY: cp_2d_i_write
  USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                             cp_blacs_gridinfo,&
                                             cp_blacs_gridinit
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: BLACS_GRID_COL,&
                                             BLACS_GRID_ROW,&
                                             BLACS_GRID_SQUARE
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_env'

  PUBLIC :: cp_blacs_env_type ! make it accessible only through cp_para_types?
  PUBLIC :: cp_blacs_env_create, cp_blacs_env_retain, cp_blacs_env_release
  PUBLIC :: cp_blacs_env_write, get_blacs_info
!***
!****************************************************************************
CONTAINS

! *****************************************************************************

  SUBROUTINE get_blacs_info(blacs_env,my_process_row,my_process_column,&
                            my_process_number,number_of_process_rows,&
                            number_of_process_columns,number_of_processes,&
                            para_env, blacs2mpi, mpi2blacs)

!   Purpose: Return informations about the specified BLACS context.

!   History: - Creation (19.06.2001, Matthias Krack)
!   moved here from qs_blacs (Joost VandeVondele)

!   ***************************************************************************

    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    INTEGER, INTENT(OUT), OPTIONAL :: my_process_row, my_process_column, &
      my_process_number, number_of_process_rows, number_of_process_columns, &
      number_of_processes
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: blacs2mpi, mpi2blacs

    IF (.NOT.ASSOCIATED(blacs_env)) CALL stop_program("get_blacs_info","no blacs environment")

    IF (PRESENT(my_process_row)) my_process_row = blacs_env%mepos(1)
    IF (PRESENT(my_process_column)) my_process_column = blacs_env%mepos(2)
    IF (PRESENT(my_process_number)) my_process_number = blacs_env%my_pid
    IF (PRESENT(number_of_process_rows)) number_of_process_rows = blacs_env%num_pe(1)
    IF (PRESENT(number_of_process_columns)) number_of_process_columns = blacs_env%num_pe(2)
    IF (PRESENT(number_of_processes)) number_of_processes = blacs_env%n_pid
    IF (PRESENT(para_env)) para_env => blacs_env%para_env
    IF (PRESENT(blacs2mpi)) blacs2mpi => blacs_env%blacs2mpi
    IF (PRESENT(mpi2blacs)) mpi2blacs => blacs_env%mpi2blacs

  END SUBROUTINE get_blacs_info

!!****f* cp_blacs_env/cp_blacs_env_create [1.0] *
!!
!!   NAME
!!     cp_blacs_env_create
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_create(blacs_env, group, para_env,&
!!         owns_group, source, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_create
!!
!!   FUNCTION
!!     allocates and initializes a type that represent a blacs context
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the type to initialize
!!     - para_env : the para_env for which a blacs env should be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_blacs_env_create(blacs_env,para_env,error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blacs_grid_layout, ipcol, &
                                                ipe, iprow, npcol, npe, &
                                                nprow, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

#ifdef __parallel
#ifndef __BLACS
  CALL cp_assert(.FALSE.,cp_failure_level, cp_assertion_failed,&
       fromWhere=routineP,message="to USE the blacs environment "//&
       "you need the blacs library : recompile with -D__BLACS (and link blacs) "//&
       CPSourceFileRef,&
       error=error)
#endif
#endif
  CPPrecondition(.NOT.ASSOCIATED(blacs_env),cp_failure_level,routineP,error,failure)

  ALLOCATE(blacs_env,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     blacs_env%group=0
     blacs_env%ref_count=1
     blacs_env%mepos(:)=0
     blacs_env%num_pe(:)=1
     blacs_env%my_pid=0
     blacs_env%n_pid=1
     CALL cp_para_env_retain(para_env,error=error)
     blacs_env%para_env => para_env

#ifdef __BLACS
     ! get the number of cpus for this blacs grid
     blacs_env%my_pid = para_env%mepos
     blacs_env%n_pid  = para_env%num_pe
     nprow=1
     npcol=1
     npe  =blacs_env%n_pid
     ! get the layout of this grid

     ! hard code for the time being the grid layout
     blacs_grid_layout=BLACS_GRID_SQUARE
     ! XXXXXX
     SELECT CASE (blacs_grid_layout)
     CASE (BLACS_GRID_SQUARE)
       IF (nprow*npcol /= npe) THEN
         DO ipe=CEILING(SQRT(REAL(npe,dp))),npe
           IF (MODULO(npe,ipe) == 0) THEN
             nprow = ipe
             npcol = npe/nprow
             EXIT
           END IF
         END DO
       END IF
     CASE (BLACS_GRID_ROW)
        nprow = 1
        npcol = npe
     CASE (BLACS_GRID_COL)
        nprow = npe
        npcol = 1
     END SELECT

     ! blacs_env%group will be set to the blacs context for this blacs env, this is not the same as the MPI context
     blacs_env%group=para_env%group
     CALL cp_blacs_gridinit(blacs_env%group,"Row-major",nprow,npcol)

     ! fill in the info one needs
     CALL cp_blacs_gridinfo(blacs_env%group,blacs_env%num_pe(1),&
          blacs_env%num_pe(2),blacs_env%mepos(1),blacs_env%mepos(2))
#endif

     ! generate the mappings blacs2mpi and mpi2blacs 
     ALLOCATE(blacs_env%blacs2mpi(0:blacs_env%num_pe(1)-1,0:blacs_env%num_pe(2)-1),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     blacs_env%blacs2mpi=0
     blacs_env%blacs2mpi(blacs_env%mepos(1),blacs_env%mepos(2))=para_env%mepos
     CALL mp_sum(blacs_env%blacs2mpi,para_env%group)
     ALLOCATE(blacs_env%mpi2blacs(2,0:para_env%num_pe-1),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     blacs_env%mpi2blacs=-1
     DO ipcol=0,blacs_env%num_pe(2)-1
        DO iprow=0,blacs_env%num_pe(1)-1
           blacs_env%mpi2blacs(1,blacs_env%blacs2mpi(iprow,ipcol))=iprow
           blacs_env%mpi2blacs(2,blacs_env%blacs2mpi(iprow,ipcol))=ipcol
        END DO
     END DO
  END IF
END SUBROUTINE cp_blacs_env_create
!***************************************************************************

!!****f* cp_blacs_env/cp_blacs_env_retain [1.0] *
!!
!!   NAME
!!     cp_blacs_env_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_retain(blacs_env, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_retain
!!
!!   FUNCTION
!!     retains the given blacs env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_blacs_env_retain(blacs_env,error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(blacs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(blacs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     blacs_env%ref_count=blacs_env%ref_count+1
  END IF
END SUBROUTINE cp_blacs_env_retain
!***************************************************************************

!!****f* cp_blacs_env/cp_blacs_env_release [1.0] *
!!
!!   NAME
!!     cp_blacs_env_release
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_release(blacs_env, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_release
!!
!!   FUNCTION
!!     releases the given blacs_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs env to relase
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_blacs_env_release(blacs_env,error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(blacs_env)) THEN
     CPPrecondition(blacs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     blacs_env%ref_count=blacs_env%ref_count-1
     IF (blacs_env%ref_count<1) THEN
        CALL cp_blacs_gridexit(blacs_env%group)
        CALL cp_para_env_release(blacs_env%para_env,error=error)
        DEALLOCATE(blacs_env%mpi2blacs, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(blacs_env%blacs2mpi,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(blacs_env,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(blacs_env)
END SUBROUTINE cp_blacs_env_release
!***************************************************************************

!!****f* cp_blacs_env/cp_blacs_env_write [1.0] *
!!
!!   NAME
!!     cp_blacs_env_write
!!
!!   SYNOPSIS
!!     Subroutine cp_blacs_env_write(blacs_env, unit_nr, error)
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_blacs_env_write
!!
!!   FUNCTION
!!     writes the description of the given blacs env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - blacs_env: the blacs environment to write
!!     - unit_nr: the unit number where to write the description of the
!!       blacs environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_blacs_env_write(blacs_env, unit_nr, error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_blacs_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iostat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(blacs_env)) THEN
     WRITE (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) blacs_env%group, blacs_env%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  mepos=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%mepos(1), blacs_env%mepos(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  num_pe=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%num_pe(1), blacs_env%num_pe(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     IF (ASSOCIATED(blacs_env%blacs2mpi)) THEN
        WRITE (unit=unit_nr,fmt="('  blacs2mpi=')",advance="no",iostat=iostat)
        CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
        CALL cp_2d_i_write(blacs_env%blacs2mpi,unit_nr=unit_nr, error=error)
     ELSE
        WRITE (unit=unit_nr,fmt="('  blacs2mpi=*null*')",iostat=iostat)
        CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(blacs_env%para_env)) THEN
        WRITE (unit=unit_nr,fmt="('  para_env=<cp_para_env id=',i6,'>,')")&
             blacs_env%para_env%group
     ELSE
        WRITE (unit=unit_nr,fmt="('  para_env=*null*')")
     END IF
     WRITE (unit=unit_nr,fmt="('  my_pid=',i10,', n_pid=',i10,' }')",&
          iostat=iostat) blacs_env%my_pid, blacs_env%n_pid
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)        
  ELSE
     WRITE (unit=unit_nr,&
          fmt="(a)", iostat=iostat) ' <cp_blacs_env>:*null* '
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
  END IF
  CALL m_flush(unit_nr)
END SUBROUTINE cp_blacs_env_write
!***************************************************************************

END MODULE cp_blacs_env
