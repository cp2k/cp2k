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
  USE cp_keywords,                     ONLY: keyword_type, &
                                             list_keywords, &
                                             find_keyword
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

  !
  ! define all keywords for use in the scf section
  ! define in the subroutine init_keywords
  ! 
  PRIVATE :: init_keywords
  INTEGER, PARAMETER, PRIVATE :: key_nrow_block=1
  INTEGER, PARAMETER, PRIVATE :: key_ncol_block=2
  INTEGER, PARAMETER, PRIVATE :: key_OT=3
  INTEGER, PARAMETER, PRIVATE :: key_GUESS=4
  INTEGER, PARAMETER, PRIVATE :: key_MIXING=5
  INTEGER, PARAMETER, PRIVATE :: key_EPS_DIIS=6
  INTEGER, PARAMETER, PRIVATE :: key_EPS_EIGVAL=7
  INTEGER, PARAMETER, PRIVATE :: key_EPS_JACOBI=8
  INTEGER, PARAMETER, PRIVATE :: key_JACOBI_THRESHOLD=9
  INTEGER, PARAMETER, PRIVATE :: key_CHOLESKY_ON=10
  INTEGER, PARAMETER, PRIVATE :: key_CHOLESKY_OFF=11
  INTEGER, PARAMETER, PRIVATE :: key_EPS_SCF=12
  INTEGER, PARAMETER, PRIVATE :: key_LEVEL_SHIFT=13
  INTEGER, PARAMETER, PRIVATE :: key_MAX_DIIS=14
  INTEGER, PARAMETER, PRIVATE :: key_MAX_SCF=15
  INTEGER, PARAMETER, PRIVATE :: key_SMEAR=16
  INTEGER, PARAMETER, PRIVATE :: key_WORK_SYEVX=17
  INTEGER, PARAMETER, PRIVATE :: key_EPS_LUMO=18
  INTEGER, PARAMETER, PRIVATE :: key_MAX_ITER_LUMO=19
  INTEGER, PARAMETER, PRIVATE :: key_WRITE_RESTART_EACH=20
  INTEGER, PARAMETER, PRIVATE :: key_ADDED_MOS=21
  INTEGER, PARAMETER, PRIVATE :: Nkeywords=21 ! should updated if keywords are added
  TYPE(keyword_type), PRIVATE, DIMENSION(1:Nkeywords) :: keywords 

!***
!****************************************************************************

!!****s* scf_control/scf_control_type [1.0] *
!!
!!   NAME
!!     scf_control_type
!!
!!   FUNCTION
!!     contains the parameters needed by a scf run
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
!!     - added_mos: additional number of MOs that might be used in the SCF
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
     CHARACTER(LEN=10)     :: density_guess
     REAL(KIND = dp)       :: eps_eigval,eps_scf,level_shift,p_mix,smear,&
                              work_syevx,eps_lumos, eps_diis,eps_jacobi,&
                              jacobi_threshold
     INTEGER               :: max_iter_lumos, max_diis
     INTEGER               :: max_scf,nrow_block,&
                              ncol_block,maxl,nkind, write_restart_each
     LOGICAL               :: gradient_functional,gth_potential_present,&
                              use_cholesky,use_ot
     INTEGER               :: ref_count,id_nr
     INTEGER, DIMENSION(2) :: added_mos
  END TYPE scf_control_type
!!***
!****************************************************************************

CONTAINS
!
! init the keywords of the scf section
! define the constants as module parameters, and increase the lenght before adding here
!
SUBROUTINE init_keywords()
  keywords(key_WRITE_RESTART_EACH)=keyword_type("WRITE_RESTART_EACH", &
                                     "WRITE_RESTART_EACH integer",&
                                     "save the RESTART file every so often",&
                                     "","","")
  keywords(key_MAX_ITER_LUMO)=keyword_type("MAX_ITER_LUMO", &
                                     "MAX_ITER_LUMO integer",&
                                     "The maximum number of iteration for the lumo computation",&
                                     "MAX_ITER_LUMOS","","")
  keywords(key_EPS_LUMO)=keyword_type("EPS_LUMO", &
                                     "EPS_LUMO real",&
                                     "target accuracy of the computation of the lumo energy",&
                                     "EPS_LUMOS","","")
  keywords(key_WORK_SYEVX)=keyword_type("WORK_SYEVX", &
                                     "WORK_SYEVX real",&
                                     "limit the memory usage of SYEVX to a fraction of the maximum size",&
                                     "","","")
  keywords(key_SMEAR)=keyword_type("SMEAR", &
                                     "SMEAR real",&
                                     "Use a finite temperature like smearing of the occupation numbers",&
                                     "","","")
  keywords(key_MAX_SCF)=keyword_type("MAX_SCF", &
                                     "MAX_SCF integer",&
                                     "Maximum number of SCF iteration to be performed for one optimization",&
                                     "MAX_SCF_STEP_NUMBER","","")
  keywords(key_MAX_DIIS)=keyword_type("MAX_DIIS", &
                                     "MAX_DIIS integer",&
                                     "Maximum number of DIIS vectors to be used",&
                                     "MAX_DIIS_BUFFER_SIZE","","")
  keywords(key_LEVEL_SHIFT)=keyword_type("LEVEL_SHIFT", &
                                     "LEVEL_SHIFT real",&
                                     "Use level shifting to improve convergence",&
                                     "LSHIFT","","")
  keywords(key_EPS_SCF)=keyword_type("EPS_SCF", &
                                     "EPS_SCF real",&
                                     "target accuracy for the scf convergence",&
                                     "","","")
  keywords(key_CHOLESKY_OFF)=keyword_type("CHOLESKY_OFF", &
                                        "CHOLESKY_OFF",&
                                        "Do not use cholesky method for computing the inverse of S",&
                                        "","","")
  keywords(key_CHOLESKY_ON)=keyword_type("CHOLESKY_ON", &
                                        "CHOLESKY_ON",&
                                        "Use cholesky method for computing the inverse of S",&
                                        "","","")
  keywords(key_JACOBI_THRESHOLD)=keyword_type("JACOBI_THRESHOLD", &
                                        "EPS_JACOBI real",&
                                        "XXXXXXXX switch to jacobi if ...",&
                                        "","","")
  keywords(key_EPS_JACOBI)=keyword_type("EPS_JACOBI", &
                                        "EPS_JACOBI real",&
                                        "XXXXXXXX target jacobi precision",&
                                        "","","")
  keywords(key_EPS_EIGVAL)=keyword_type("EPS_EIGVAL", &
                                        "EPS_DIIS real",&
                                        "Throw away linear combinations of basis functions with a small eigenvalue in S",&
                                        "","","")
  keywords(key_EPS_DIIS)=  keyword_type("EPS_DIIS", &
                                        "EPS_DIIS real",&
                                        "Threshold to start using DIAG/DIIS",&
                                        "","","")
  keywords(key_MIXING)=    keyword_type("MIXING", &
                                        "MIXING real",&
                                        "Fraction of new density matrix to be mixed in",&
                                        "DENSITY_MIXING","","")
  keywords(key_GUESS)=     keyword_type("GUESS", &
                                        "GUESS {ATOMIC|RESTART|RANDOM}",&
                                        "Change the initial guess for the wavefunction",&
                                        "SCF_GUESS","","DENSITY_GUESS")
  keywords(key_OT)=        keyword_type("OT", &
                                        "OT",&
                                        "Use the orbital transformation method",&
                                        "","","")
  keywords(key_ncol_block)=keyword_type("NCOL_BLOCK", &
                                        "NCOL_BLOCK integer",&
                                        "sets the number of cols in a scalapack block",&
                                        "","","")
  keywords(key_nrow_block)=keyword_type("NROW_BLOCK", &
                                        "NROW_BLOCK integer",&
                                        "sets the number of rows in a scalapack block",&
                                        "","","")
  keywords(key_ADDED_MOS)=keyword_type("ADDED_MOS", &
                                        "ADDED_MOS integer integer",&
                                        "number of additional MOS added for each spin",&
                                        "","","")
END SUBROUTINE

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
       scf_control%added_mos=0
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
       CALL init_keywords()

       CALL start_parser(file_name=globenv%input_file_name,&
            globenv=globenv,&
            start_section_label="SCF")
       DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (find_keyword(keywords,keyword))
          CASE (key_NROW_BLOCK)
             CALL read_object(scf_control%nrow_block)
          CASE (key_NCOL_BLOCK)
             CALL read_object(scf_control%ncol_block)
          CASE (key_OT)
             scf_control%use_ot = .TRUE.
          CASE (key_GUESS)
             CALL read_object(scf_control%density_guess,lower_to_upper=.TRUE.)
          CASE (key_MIXING)
             CALL read_object(scf_control%p_mix)
          CASE (key_EPS_DIIS)
             CALL read_object(scf_control%eps_diis)
          CASE (key_EPS_EIGVAL)
             CALL read_object(scf_control%eps_eigval)
             scf_control%use_cholesky = .FALSE.
          CASE (key_EPS_JACOBI)
             CALL read_object(scf_control%eps_jacobi)
          CASE (key_JACOBI_THRESHOLD)
             CALL read_object(scf_control%jacobi_threshold)
          CASE (key_CHOLESKY_ON)
             scf_control%use_cholesky = .TRUE.
          CASE (key_CHOLESKY_OFF)
             scf_control%use_cholesky = .FALSE.
          CASE (key_EPS_SCF)
             CALL read_object(scf_control%eps_scf)
          CASE (key_LEVEL_SHIFT)
             CALL read_object(scf_control%level_shift)
             IF (scf_control%level_shift /= 0.0_dp) &
                  scf_control%use_cholesky = .FALSE.
          CASE (key_MAX_DIIS)
             CALL read_object(scf_control%max_diis)
          CASE (key_MAX_SCF)
             CALL read_object(scf_control%max_scf)
          CASE (key_SMEAR)
             CALL read_object(scf_control%smear)
          CASE (key_WORK_SYEVX)
             CALL read_object(scf_control%work_syevx)
             scf_control%work_syevx = MIN(MAX(0.0_dp,scf_control%work_syevx),1.0_dp)
          CASE (key_EPS_LUMO)
             CALL read_object(scf_control%eps_lumos)
          CASE (key_MAX_ITER_LUMO)
             CALL read_object(scf_control%max_iter_lumos)
          CASE (key_WRITE_RESTART_EACH)
             CALL read_object(scf_control%write_restart_each)
          CASE (key_ADDED_MOS)
             ! specify the number for each spin
             CALL read_object(scf_control%added_mos(1))
             CALL read_object(scf_control%added_mos(2))
          CASE DEFAULT
            CALL list_keywords(keywords,globenv%scr)
            CALL stop_parser(routine=routineN,message="INVALID_KEYWORD "&
                 //keyword)
          END SELECT
       END DO

       CALL finish_parser()

       ! check for CORE .AND. OT
       IF (scf_control%use_ot .AND. &
            (scf_control%density_guess=="CORE")) THEN
          CALL stop_program("scf_c_read_parameters","Use GUESS {ATOMIC,RESTART,RANDOM} with OT")
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
          WRITE (UNIT=output_unit,&
           FMT="(/,/,T2,A,T25,A,T71,A10,/,T25,56('-'),3(/,T25,A,T76,I5),&
           &T25,56('-'),5(/,T25,A,T72,ES9.2),T25,56('-'),4(/,T25,A,T76,F5.2),1(/,T25,A,T71,2I5))")&
               "SCF PARAMETERS",&
               "Density guess:     ",ADJUSTR(TRIM(scf_control%density_guess)),&
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
               "smear [a.u.]:      ",scf_control%smear,&
               "added MOs          ",scf_control%added_mos
       END IF
    END IF
  END SUBROUTINE scf_c_write_parameters

! **************************************************************************

END MODULE scf_control_types

! **************************************************************************
