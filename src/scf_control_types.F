!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/scf_control_types [1.0] *
!!
!!   NAME
!!     scf_control_types
!!
!!   FUNCTION
!!     parameters that control an scf iteration
!!
!!   NOTES
!!      not in cp_control_types, to separate operator related parameters from
!!      method related parameters (as suggested by Matthias)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE scf_control_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='scf_control_types'
  INTEGER, SAVE, PRIVATE :: last_scf_c_id_nr=0

  PUBLIC :: scf_control_type
  PUBLIC :: scf_c_create, scf_c_retain, scf_c_release, scf_c_read_parameters,&
       scf_c_write_parameters
!***
!****************************************************************************

!!****s* scf_control/scf_control_type [1.0] *
!!
!!   NAME
!!     scf_control_type
!!
!!   FUNCTION
!!     contains the parameters nedeed by a scf run
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - density_guess: how to choose the initial density
!!       (CORE,RANDOM,RESTART,ATOMIC)
!!     - eps_eigval: wanted error on the eigenvalues
!!     - eps_scf: whanted error on the whole scf
!!     - level_shift: amount of level shift
!!     - p_mix: how to mix the new and old densities in non diss iterations
!!     - smear: amount of smeating applied at the orbitals
!!     - work_syevx: real to control the amount of temporary memory needed
!!       by lapack syevx
!!     - eps_lumos: error on the lumos calculated at the end of the scf
!!     - max_iter_lumus: maxumum number of iterations used to calculate
!!       the lumos at the end of the scf
!!     - max_scf: max scf iterations
!!     - write_restart_each: after how many iterations a restart file 
!!       should be written if the convergence is not yet acheived
!!       (<0 never, 0 only upon convergence)
!!     - nrow_block: number of rows of a full matrix block
!!     - ncol_block: number of of columns in a block of a full matrix
!!     - id_nr: unique number to identify an scf control
!!     - ref_count: reference count (see cp2k/doc/ReferenceCounting.html)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE scf_control_type
     CHARACTER(LEN=10) :: density_guess
     REAL(KIND = dp)          :: eps_eigval,eps_scf,level_shift,p_mix,smear,&
          work_syevx,eps_lumos, eps_diis,eps_jacobi,jacobi_threshold
     INTEGER           :: max_iter_lumos, max_diis
     INTEGER           :: max_scf,nrow_block,&
          ncol_block,maxl,nkind, write_restart_each
     LOGICAL           :: gradient_functional,gth_potential_present,&
          use_cholesky,use_ot
     INTEGER           :: ref_count,id_nr
  END TYPE scf_control_type
!!***
!****************************************************************************

CONTAINS

! *****************************************************************************

!!****f* scf_control/scf_c_create [1.0] *
!!
!!   NAME
!!     scf_c_create
!!
!!   SYNOPSIS
!!     Subroutine scf_c_create(scf_control, error)
!!       Type(scf_control_type), Pointer:: scf_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_c_create
!!
!!   FUNCTION
!!     allocates and initializes an scf control object with the default values
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - scf_control: the object to initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_c_create(scf_control,error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_c_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(scf_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
!   *** Load the default values ***
       scf_control%density_guess = "ATOMIC"
       scf_control%eps_eigval = 1.0E-5_dp
       scf_control%eps_scf = 1.0E-5_dp
       scf_control%eps_lumos = 1.0E-5_dp
       scf_control%max_iter_lumos = 2999
       scf_control%eps_diis = 0.1_dp
       scf_control%level_shift = 0.0_dp
       scf_control%max_diis = 4
       scf_control%max_scf = 50
       scf_control%nrow_block = 32
       scf_control%ncol_block = 32
       scf_control%p_mix = 0.4_dp
       scf_control%smear = 0.0_dp
       scf_control%work_syevx = 1.0_dp
       scf_control%use_cholesky = .TRUE.
       scf_control%use_ot = .FALSE.
       scf_control%max_diis = 4
       scf_control%eps_diis = 0.1_dp
       scf_control%ref_count=1
       last_scf_c_id_nr=last_scf_c_id_nr+1
       scf_control%id_nr=last_scf_c_id_nr
       scf_control%eps_jacobi = 0.0_dp
       scf_control%jacobi_threshold = 1.0E-7_dp
       scf_control%write_restart_each=10
    END IF
  END SUBROUTINE scf_c_create
!***************************************************************************

!!****f* scf_control/scf_c_retain [1.0] *
!!
!!   NAME
!!     scf_c_retain
!!
!!   SYNOPSIS
!!     Subroutine scf_c_retain(scf_control, error)
!!       Type(scf_control_type), Pointer:: scf_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_c_retain
!!
!!   FUNCTION
!!     retains the given scf_control (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - scf_control: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_c_retain(scf_control,error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_c_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
       scf_control%ref_count=scf_control%ref_count+1
    END IF
  END SUBROUTINE scf_c_retain
!***************************************************************************

!!****f* scf_control/scf_c_release [1.0] *
!!
!!   NAME
!!     scf_c_release
!!
!!   SYNOPSIS
!!     Subroutine scf_c_release(scf_control, error)
!!       Type(scf_control_type), Pointer:: scf_control
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine scf_c_release
!!
!!   FUNCTION
!!     releases the given scf_control (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     at the moment does nothing
!!
!!   INPUTS
!!     - scf_control: the object to free
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_c_release(scf_control,error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_c_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(scf_control)) THEN
       CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
       scf_control%ref_count=scf_control%ref_count-1
       IF (scf_control%ref_count<1) THEN
          DEALLOCATE(scf_control,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(scf_control)
  END SUBROUTINE scf_c_release
!***************************************************************************

!!****f* scf_control/scf_c_read_parameters [1.0] *
!!
!!   NAME
!!     scf_c_read_parameters
!!
!!   SYNOPSIS
!!     Subroutine scf_c_read_parameters(scf_control, globenv, error)
!!       Type(scf_control_type), Pointer:: scf_control
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine scf_c_read_parameters
!!
!!   FUNCTION
!!     reads the parameters of the scf section into the given scf_control
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - scf_control: the object that wil contain the values read
!!     - globenv: to define the file to read
!!     - error: controls log and error handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     05.2001 created [Matthias]
!!     09.2002 creaded separated scf_control type [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_c_read_parameters(scf_control, globenv,error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_c_read_parameters', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: keyword
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL start_parser(file_name=globenv%input_file_name,&
            globenv=globenv,&
            start_section_label="SCF")

       DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("NROW_BLOCK")
             CALL read_object(scf_control%nrow_block)
          CASE ("NCOL_BLOCK")
             CALL read_object(scf_control%ncol_block)
          CASE ("OT")
             scf_control%use_ot = .TRUE.
          CASE ("DENSITY_GUESS","SCF_GUESS","GUESS")
             CALL read_object(scf_control%density_guess,lower_to_upper=.TRUE.)
          CASE ("DENSITY_MIXING","MIXING")
             CALL read_object(scf_control%p_mix)
          CASE ("EPS_DIIS")
             CALL read_object(scf_control%eps_diis)
          CASE ("EPS_EIGVAL")
             CALL read_object(scf_control%eps_eigval)
             scf_control%use_cholesky = .FALSE.
          CASE ("EPS_JACOBI")
             CALL read_object(scf_control%eps_jacobi)
          CASE ("JACOBI_THRESHOLD")
             CALL read_object(scf_control%jacobi_threshold)
          CASE ("CHOLESKY_ON")
             scf_control%use_cholesky = .TRUE.
          CASE ("CHOLESKY_OFF")
             scf_control%use_cholesky = .FALSE.
          CASE ("EPS_SCF")
             CALL read_object(scf_control%eps_scf)
          CASE ("LEVEL_SHIFT","LSHIFT")
             CALL read_object(scf_control%level_shift)
             IF (scf_control%level_shift /= 0.0_dp) &
                  scf_control%use_cholesky = .FALSE.
          CASE ("MAX_DIIS_BUFFER_SIZE","MAX_DIIS")
             CALL read_object(scf_control%max_diis)
          CASE ("MAX_SCF_STEP_NUMBER","MAX_SCF")
             CALL read_object(scf_control%max_scf)
          CASE ("SMEAR")
             CALL read_object(scf_control%smear)
          CASE ("WORK_SYEVX")
             CALL read_object(scf_control%work_syevx)
             scf_control%work_syevx = MIN(MAX(0.0_dp,scf_control%work_syevx),1.0_dp)
          CASE ("EPS_LUMOS","EPS_LUMO")
             CALL read_object(scf_control%eps_lumos)
          CASE ("MAX_ITER_LUMOS","MAX_ITER_LUMO")
             CALL read_object(scf_control%max_iter_lumos)
          CASE ("WRITE_RESTART_EACH")
             CALL read_object(scf_control%write_restart_each)
          CASE DEFAULT
            CALL stop_parser(routine=routineN,message="INVALID_KEYWORD "&
                 //keyword)
          END SELECT
       END DO

       CALL finish_parser()

       ! check for ATOMIC .OR. CORE .AND. OT
       IF (scf_control%use_ot .AND. &
            (scf_control%density_guess=="ATOMIC" .OR. scf_control%density_guess=="CORE")) THEN
          CALL stop_program("scf_c_read_parameters","Use GUESS {RESTART,RANDOM} with OT")
       ENDIF

    END IF
  END SUBROUTINE scf_c_read_parameters
!***************************************************************************

!!****f* scf_control/scf_c_write_parameters [1.0] *
!!
!!   NAME
!!     scf_c_write_parameters
!!
!!   SYNOPSIS
!!     Subroutine scf_c_write_parameters(scf_control, globenv, error)
!!       Type(scf_control_type), Pointer:: scf_control
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine scf_c_write_parameters
!!
!!   FUNCTION
!!     writes out the scf parameters
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - scf_control: the object you want to print
!!     - globenv: to know where to print
!!     - error: controls log and error handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     05.2001 created [Matthias]
!!     09.2002 created separated scf_control type [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_c_write_parameters(scf_control,globenv, error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_c_write_parameters', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (.NOT.globenv%ionode) RETURN

       IF (scf_control%max_scf < 1) RETURN

       IF (globenv%print%scf) THEN
          output_unit = globenv%scr
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/)") "SCF PARAMETERS"
          WRITE (UNIT=output_unit,&
               FMT="(T3,A,/,3(/,T3,A,I5),5(/,T3,A,ES9.2),4(/,T3,A,F5.2))")&
               "Density guess: "//TRIM(scf_control%density_guess),&
               "max_scf:           ",scf_control%max_scf,&
               "max_diis:          ",scf_control%max_diis,&
               "write_restart_each:",scf_control%write_restart_each,&
               "eps_scf:           ",scf_control%eps_scf,&
               "eps_diis:          ",scf_control%eps_diis,&
               "eps_eigval:        ",scf_control%eps_eigval,&
               "eps_jacobi:        ",scf_control%eps_jacobi,&
               "jacobi_threshold:  ",scf_control%jacobi_threshold,&
               "p_mix:             ",scf_control%p_mix,&
               "work_syevx:        ",scf_control%work_syevx,&
               "level_shift [a.u.]:",scf_control%level_shift,&
               "smear [a.u.]:      ",scf_control%smear
       END IF
    END IF
  END SUBROUTINE scf_c_write_parameters
! **************************************************************************

END MODULE scf_control_types
