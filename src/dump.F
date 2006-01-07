!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dump [1.0] *
!!
!!   NAME
!!     dump
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     CJM 2/01:  Now dumps out all extended_system variables
!!     and the box_ref
!!     - Simple dump for QS without extended variables (MK,15.09.2003)
!!
!!   SOURCE
!******************************************************************************

MODULE dump
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_write_variables
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast, mp_allgather
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="dump"

  PUBLIC :: dump_variables

!!*****
!******************************************************************************

CONTAINS

! *****************************************************************************

!!****f* dump/dump_variables [1.0] *
!!
!!   NAME
!!     dump_force_env_variables
!!
!!   FUNCTION
!!     dumps the variables 
!!
!!   NOTES
!!     this call call the right dump_*_variables, but I find it more clean
!!     to call the force_env routine
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************

  SUBROUTINE dump_variables(md_env,dump_file_name,error)

    TYPE(md_environment_type), POINTER       :: md_env
    CHARACTER(LEN=*), INTENT(IN)             :: dump_file_name
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "dump_force_env_variables", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: dump_unit, handle, i, j, size1, size2
    INTEGER                                  :: nhc_len, num_nhc, isos, numneed
    INTEGER                                  :: tot_nhcneed, iproc, counter
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger, new_logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), POINTER   :: npt( :, : )
    INTEGER, POINTER :: nhc_per_proc ( : ),index(:)
    REAL (KIND=dp), POINTER, DIMENSION ( : ) :: eta, veta, fnhc, mnhc, work

    CHARACTER (LEN=default_string_length)    :: cmdstr
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    CALL get_md_env ( md_env, force_env_p = force_env_p, nhc_part = nhc_part, &
                      nhc_baro = nhc_baro, npt = npt )
    CALL timeset(routineN,"I",'',handle)

    NULLIFY (para_env,new_logger)

    CALL get_md_env ( md_env, para_env = para_env )

    logger => cp_error_get_logger(error)
    CALL cp_logger_create(new_logger,para_env=para_env,template_logger=logger)
    CALL cp_error_init(new_error,logger=new_logger,template_error=error)

    IF (para_env%mepos==para_env%source) THEN

      INQUIRE(FILE=dump_file_name,EXIST=found)
      IF (found) THEN
        WRITE(UNIT=cmdstr,FMT="(A)")&
          "mv "//TRIM(dump_file_name)//" "//TRIM(dump_file_name)//".bak"
        CALL system(cmdstr)
      END IF
      
      CALL open_file(file_name=dump_file_name,&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_status="REPLACE",&
                     unit_number=dump_unit)
    END IF

    IF (ASSOCIATED(force_env_p)) THEN
       DO i = 1, SIZE(force_env_p)
          CALL force_env_write_variables(force_env_p(i)%force_env, &
                                           unit_nr=dump_unit,error=new_error)
       END DO
    END IF


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! write out nhc_part variable in a consistent way for restarts
!-----------------------------------------------------------------------------
 nhc_len = SIZE(nhc_part(1)%nvt,1)
 num_nhc = SIZE(nhc_part(1)%nvt,2)
 numneed = num_nhc
 ALLOCATE ( nhc_per_proc( para_env % num_pe ), STAT = isos )
 IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to allocate nhc_per_proc')
 nhc_per_proc ( : ) = 0

 CALL mp_allgather(numneed,nhc_per_proc,para_env%group)

 IF (nhc_part ( 1 ) % dis_type=='INTER_REP') THEN
   tot_nhcneed = nhc_per_proc(1)
 ELSE
   tot_nhcneed = SUM ( nhc_per_proc )
 ENDIF
  
! This is only activated for NVE runs to print so one can use the
! ALL restart option with NVE
 
 IF (para_env%mepos==para_env%source) THEN
     IF ( ( nhc_len == 0 ).AND. ( tot_nhcneed == 0 ) ) THEN
        WRITE (UNIT=dump_unit,FMT=*) nhc_len, tot_nhcneed
     END IF
 ENDIF

  NULLIFY(work,index)
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! nhc_part%eta  
  !-----------------------------------------------------------------------------
  ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to allocate eta')
  DO iproc=1,para_env%num_pe
    CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
    CALL reallocate(index,1,nhc_per_proc(iproc))
    IF(para_env%mepos == (iproc-1)) THEN
      index(:)=0
      counter=0
      DO i=1,nhc_len
        DO j=1,num_nhc
          counter=counter+1
          work(counter)=nhc_part(1)%nvt(i,j)%eta
          index(j) = nhc_part ( 1 ) % index ( j )
        END DO
      END DO
    ELSE
      work(:) = 0.0_dp
    END IF
    CALL mp_bcast(work,iproc-1,para_env%group)
    CALL mp_bcast(index,iproc-1,para_env%group)
    counter=0
    DO i=1,nhc_len
      DO j=1,nhc_per_proc(iproc)
        counter = counter+1
        eta((index(j)-1)*nhc_len+i) = work(counter)
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! nhc_part % veta   
  !-----------------------------------------------------------------------------
  ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to allocate veta')
  DO iproc=1,para_env%num_pe
    CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
    CALL reallocate(index,1,nhc_per_proc(iproc))
    IF(para_env%mepos == (iproc-1)) THEN
      index(:)=0
      counter=0
      DO i=1,nhc_len
        DO j=1,num_nhc
          counter=counter+1
          work(counter)=nhc_part(1)%nvt(i,j)%v
          index(j) = nhc_part ( 1 ) % index ( j )
        END DO
      END DO
    ELSE
      work(:) = 0.0_dp
    END IF
    CALL mp_bcast(work,iproc-1,para_env%group)
    CALL mp_bcast(index,iproc-1,para_env%group)
    counter=0
    DO i=1,nhc_len
      DO j=1,nhc_per_proc(iproc)
        counter = counter+1
        veta((index(j)-1)*nhc_len+i) = work(counter)
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! nhc_part % force
  !-----------------------------------------------------------------------------
  ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to allocate fnhc')
  DO iproc=1,para_env%num_pe
    CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
    CALL reallocate(index,1,nhc_per_proc(iproc))
    IF(para_env%mepos == (iproc-1)) THEN
      index(:)=0
      counter=0
      DO i=1,nhc_len
        DO j=1,num_nhc
          counter=counter+1
          work(counter)=nhc_part(1)%nvt(i,j)%f
          index(j) = nhc_part ( 1 ) % index ( j )
        END DO
      END DO
    ELSE
      work(:) = 0.0_dp
    END IF
    CALL mp_bcast(work,iproc-1,para_env%group)
    CALL mp_bcast(index,iproc-1,para_env%group)
    counter=0
    DO i=1,nhc_len
      DO j=1,nhc_per_proc(iproc)
        counter = counter+1
        fnhc((index(j)-1)*nhc_len+i) = work(counter)
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! nhc_part % mass
  !-----------------------------------------------------------------------------
  ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = isos)
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to allocate mnhc')
  DO iproc=1,para_env%num_pe
    CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
    CALL reallocate(index,1,nhc_per_proc(iproc))
    IF(para_env%mepos == (iproc-1)) THEN
      index(:)=0
      counter=0
      DO i=1,nhc_len
        DO j=1,num_nhc
          counter=counter+1
          work(counter)=nhc_part(1)%nvt(i,j)%mass
          index(j) = nhc_part ( 1 ) % index ( j )
        END DO
      END DO
    ELSE
      work(:) = 0.0_dp
    END IF
    CALL mp_bcast(work,iproc-1,para_env%group)
    CALL mp_bcast(index,iproc-1,para_env%group)
    counter=0
    DO i=1,nhc_len
      DO j=1,nhc_per_proc(iproc)
        counter = counter+1
        mnhc((index(j)-1)*nhc_len+i) = work(counter)
      END DO
    END DO
  END DO



  DEALLOCATE ( work, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate work')
  DEALLOCATE ( index, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate index')
  DEALLOCATE ( nhc_per_proc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate nhc_per_proc')





  IF ( ( nhc_part ( 1 ) % dis_type=='DIS_REP' ) .OR.  &
       ( nhc_part ( 1 ) % dis_type == 'NONE') ) THEN
    IF (para_env%mepos==para_env%source) THEN

       size1 = SIZE(eta)

       WRITE (UNIT=dump_unit,FMT=*) nhc_len, tot_nhcneed

       DO i=1,size1
        WRITE (UNIT=dump_unit,FMT=*) eta ( i )
       END DO

       DO i=1,size1
        WRITE (UNIT=dump_unit,FMT=*) veta ( i )
       END DO

       DO i=1,size1
        WRITE (UNIT=dump_unit,FMT=*) fnhc ( i )
       END DO

       DO i=1,size1
        WRITE (UNIT=dump_unit,FMT=*) mnhc ( i )
       END DO
    END IF
  ELSEIF ( nhc_part ( 1 ) % dis_type=='INTER_REP' ) THEN
    IF (para_env%mepos==para_env%source) THEN

       size1 = SIZE(nhc_part(1)%nvt,1)
       size2 = SIZE(nhc_part(1)%nvt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%eta
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%f
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%mass
         END DO
       END DO
    END IF
  END IF

  DEALLOCATE ( eta, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate eta')
  DEALLOCATE ( veta, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate veta')
  DEALLOCATE ( fnhc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate fnhc')
  DEALLOCATE ( mnhc, STAT = isos )
  IF ( isos /= 0 ) CALL stop_program ( 'initialize_nhc', &
                               'failed to deallocate mnhc')

  IF (para_env%mepos==para_env%source) THEN

!     *** barostat thermostat variables***
       size1 = SIZE(nhc_baro(1)%nvt,1)
       size2 = SIZE(nhc_baro(1)%nvt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%eta
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%f
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%mass
         END DO
       END DO

!     *** barostat  variables***
       size1 = SIZE(npt,1)
       size2 = SIZE(npt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) npt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) npt(i,j)%mass
         END DO
       END DO

!MK!     *** Close dump unit and keep it ***

      CALL close_file (unit_number=dump_unit)

  END IF

  CALL cp_error_dealloc_ref(new_error)
  CALL cp_logger_release(new_logger)

  CALL timestop(0.0_dp,handle)

  END SUBROUTINE dump_variables

!***************************************************************************

END MODULE dump
