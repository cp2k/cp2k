!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
!!
!!   SOURCE
!****************************************************************************
MODULE cp_blacs_env
  USE cp_array_i_utils,                ONLY: cp_2d_i_write
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: BLACS_GRID_COL,&
                                             BLACS_GRID_ROW,&
                                             BLACS_GRID_SQUARE,&
                                             global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_blacs_env'

  PUBLIC :: cp_blacs_env_type ! make it accessible only through cp_para_types?
  PUBLIC :: cp_blacs_env_create, cp_blacs_env_retain, cp_blacs_env_release
  PUBLIC :: start_blacs, end_blacs, cp_blacs_env_write, get_blacs_info
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

    TYPE(cp_blacs_env_type), POINTER          :: blacs_env
    INTEGER, OPTIONAL, INTENT(OUT)            :: my_process_column,&
                                                 my_process_number,&
                                                 my_process_row,&
                                                 number_of_process_columns,&
                                                 number_of_process_rows,&
                                                 number_of_processes
    TYPE(cp_para_env_type), POINTER, OPTIONAL  :: para_env
    INTEGER, DIMENSION(:,:), POINTER, OPTIONAL :: blacs2mpi, mpi2blacs

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
!!       Integer, Intent (IN):: group
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Logical, Intent (IN), Optional:: owns_group
!!       Integer, Dimension(2), Intent (IN), Optional:: source
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
!!     - group: the id of the blacs context
!!     - owns_group: if this type owns the group and should release it
!!       (via blacs_gridexit) when deallocated
!!     - source: position of a special processor (defaults to (0,0))
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
SUBROUTINE cp_blacs_env_create(blacs_env,group,para_env,owns_group,&
     source,error)
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  INTEGER, INTENT(in) :: group
  TYPE(cp_para_env_type), POINTER :: para_env
  LOGICAL, INTENT(in), OPTIONAL :: owns_group
  INTEGER, DIMENSION(2), INTENT(in), OPTIONAL :: source
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_blacs_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, ipcol,iprow

  failure=.FALSE.

#ifdef __parallel
#ifndef __BLACS
  CALL cp_assert(.FALSE.,cp_failure_level, cp_assertion_failed,&
       fromWhere=routineP,message="to USE the blacs environment "//&
       "(used by full matrices) you need the blacs library "//&
       CPSourceFileRef,&
       error=error)
#endif
#endif

  ALLOCATE(blacs_env,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     blacs_env%group=group
     blacs_env%owns_group=.TRUE.
     blacs_env%ref_count=1
     blacs_env%mepos=0
     blacs_env%num_pe=1
     blacs_env%source=0
     blacs_env%my_pid=0
     blacs_env%n_pid=1
     IF (PRESENT(owns_group)) blacs_env%owns_group=owns_group
     IF (PRESENT(source)) THEN
        blacs_env%source=source
     END IF
     CALL cp_para_env_retain(para_env,error=error)
     blacs_env%para_env => para_env
#ifdef __BLACS
     CALL blacs_pinfo(blacs_env%my_pid,blacs_env%n_pid)
     CALL blacs_gridinfo(blacs_env%group,blacs_env%num_pe(1),&
          blacs_env%num_pe(2),blacs_env%mepos(1),blacs_env%mepos(2))
#endif
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
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_blacs_env_retain',&
        routineP=moduleN//':'//routineN
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
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_blacs_env_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat
  failure=.FALSE.
  
  IF (ASSOCIATED(blacs_env)) THEN
     CPPrecondition(blacs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     blacs_env%ref_count=blacs_env%ref_count-1
     IF (blacs_env%ref_count<1) THEN
        IF (blacs_env%owns_group) THEN
#ifdef __BLACS
           CALL blacs_gridexit(blacs_env%group)
#endif
        END IF
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
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  INTEGER, INTENT(in) :: unit_nr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_blacs_env_write',&
        routineP=moduleN//':'//routineN
  INTEGER :: iostat
  failure=.FALSE.
  
  IF (ASSOCIATED(blacs_env)) THEN
     WRITE (unit=unit_nr,&
          fmt="(' <cp_blacs_env>:{ owns_group=',l1,',')",&
          iostat=iostat) blacs_env%owns_group
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  group=',i10,', ref_count=',i10,',')",&
          iostat=iostat) blacs_env%group, blacs_env%ref_count
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  mepos=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%mepos(1), blacs_env%mepos(2)
     CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
     WRITE (unit=unit_nr,fmt="('  source=(',i8,',',i8,'),')",&
          iostat=iostat) blacs_env%source(1), blacs_env%source(2)
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

!***************************************************************************

  SUBROUTINE start_blacs(blacs_env,globenv)

    ! here we will need improvement, we want a blacs context for the given MPI group, which might be a subset
    ! of all the processes. Here, we get a grid for MPI_WORLD by default (I think).
    ! blacs_gridmap and related might be the way to go.

!   Purpose: Initialize a BLACS process grid. The BLACS context is returned.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE start_blacs (MODULE cp_blacs_env)"

!   *** Local variables ***

    INTEGER :: context,group,ipe,istat,mypcol,mype,myprow,npcol,npe,nprow,&
               output_unit
    LOGICAL :: ionode
    TYPE(cp_error_type) :: error

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol,pcon,prow

!   ---------------------------------------------------------------------------

    CALL cp_error_init(error)
#if defined(__parallel) && (! defined(__BLACS) || ! defined(__SCALAPACK))
    CALL stop_program(routine,&
                      "Quickstep needs BLACS and SCALAPACK for parallel runs",&
                      globenv)
#endif

#ifdef __parallel
#ifndef __BLACS
  CALL cp_assert(.FALSE.,cp_failure_level, cp_assertion_failed,&
       fromWhere=routine,message="to USE the blacs environment "//&
       "(used by full matrices) you need the blacs library "//&
       CPSourceFileRef,&
       error=error)
#endif
#endif

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    context = 0
    npcol = 0
    nprow = 0
#if defined(__BLACS)
    ! here we will need improvement, we want a blacs context for the given MPI group, which might be a subset
    ! of all the processes. Here, we get a grid for MPI_WORLD by default (I think).
    ! blacs_gridmap and related might be the way to go.
    CALL blacs_pinfo(mype,npe)
    CALL blacs_get(-1,0,context)

    SELECT CASE (globenv%blacs_grid_layout)
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

    CALL blacs_gridinit(context,"Row-major",nprow,npcol)
    CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)

    IF (globenv%print%blacs_info) THEN

      ALLOCATE (prow(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow",npe*int_size)
      prow(:) = 0
      ALLOCATE (pcol(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol",npe*int_size)
      pcol(:) = 0
      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0

      prow(mype) = myprow
      pcol(mype) = mypcol
      pcon(mype) = context

      CALL mp_sum(prow,group)
      CALL mp_sum(pcol,group)
      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,(T3,A,T32,I6))")&
          "BLACS INFORMATION (BLACS started)",&
          "Number of processes:         ",nprow*npcol,&
          "Number of process rows:      ",nprow,&
          "Number of process columns:   ",npcol
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          "Process   Process row   Process column   BLACS context"
        WRITE (UNIT=output_unit,FMT="(T4,I6,T18,I6,T35,I6,T45,I12)")&
          (ipe,prow(ipe),pcol(ipe),pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (prow,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"prow")
      DEALLOCATE (pcol,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcol")
      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

#endif
    CALL cp_blacs_env_create(blacs_env,group = context, para_env=globenv%para_env)

  END SUBROUTINE start_blacs

! *****************************************************************************

  SUBROUTINE end_blacs(blacs_env,globenv)

!   Purpose: Release the resources of a BLACS context.

!   History: - Creation (22.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE(cp_blacs_env_type), POINTER :: blacs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE end_blacs (MODULE cp_blacs_env)"

!   *** Local variables ***

    INTEGER :: group,ipe,istat,mype,npe,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcon
    TYPE(cp_blacs_env_type), POINTER :: context

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    context => blacs_env
#if defined(__BLACS)

    IF (globenv%print%blacs_info) THEN

      npe=context%n_pid
      mype=context%my_pid

      ALLOCATE (pcon(0:npe-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon",npe*int_size)
      pcon(:) = 0

      pcon(mype) = context%group

      CALL mp_sum(pcon,group)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "BLACS INFORMATION (BLACS finished)"
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/)")&
          " PE   BLACS context"
        WRITE (UNIT=output_unit,FMT="(I5,T10,I12)")&
          (ipe,pcon(ipe),ipe=0,npe-1)
      END IF

      DEALLOCATE (pcon,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pcon")

    END IF

    CALL mp_sync(group)

    CALL cp_blacs_env_release(blacs_env)
    !    CALL blacs_gridexit(context)

    !globenv%context = 0
    !globenv%nprow = 0
    !globenv%npcol = 0

#endif
  END SUBROUTINE end_blacs

! *****************************************************************************

END MODULE cp_blacs_env
