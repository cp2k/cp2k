!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 - 2006 CP2K developers group                           !
!-----------------------------------------------------------------------------!
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
!******************************************************************************

MODULE scf_control_types

  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: atomic_guess,&
                                             core_guess,&
                                             general_roks,&
                                             high_spin_roks,&
                                             outer_scf_ddapc_constraint,&
                                             outer_scf_none,&
                                             outer_scf_optimizer_diis,&
                                             outer_scf_optimizer_none,&
                                             outer_scf_optimizer_sd,&
                                             outer_scf_s2_constraint
  USE input_cp2k_dft,                  ONLY: create_scf_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "scf_control_types"
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

  INTEGER, SAVE :: last_scf_c_id_nr = 0

  ! Public data types

  PUBLIC :: outer_scf_control_type,&
            scf_control_type

  ! Public subroutines

  PUBLIC :: scf_c_create,&
            scf_c_read_parameters,&
            scf_c_release,&
            scf_c_retain,&
            scf_c_write_parameters

!!*****************************************************************************
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
!!       (CORE,RANDOM,RESTART,ATOMIC,FROZEN)
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
!******************************************************************************

  TYPE outer_scf_control_type
    LOGICAL               :: have_scf
    INTEGER               :: max_scf
    REAL(KIND = dp)       :: eps_scf,step_size
    INTEGER               :: TYPE
    INTEGER               :: optimizer
    INTEGER               :: diis_buffer_length
    INTEGER               :: extrapolation_order
  END TYPE outer_scf_control_type

  TYPE scf_control_type
    TYPE(outer_scf_control_type) :: outer_scf
    INTEGER               :: density_guess
    REAL(KIND = dp)       :: eps_eigval,eps_scf,level_shift,p_mix,smear,&
                             work_syevx,eps_lumos, eps_diis,eps_jacobi,&
                             jacobi_threshold, gspace_mixing_a, gspace_mixing_b
    INTEGER               :: max_iter_lumos, max_diis
    INTEGER               :: max_scf,nrow_block,max_scf_aspc,&
                             ncol_block,maxl,nkind, write_restart_each, &
                             harris_max_scf
    LOGICAL               :: gradient_functional,gth_potential_present,&
                             use_cholesky,use_ot,fill_wf_history
    INTEGER               :: ref_count,id_nr
    INTEGER, DIMENSION(2) :: added_mos
    INTEGER                               :: roks_scheme
    REAL(KIND=dp)                         :: roks_f
    REAL(KIND=dp), DIMENSION(0:2,0:2,1:2) :: roks_parameter
   END TYPE scf_control_type
!!***
!*****************************************************************************

CONTAINS

!******************************************************************************
!!****f* scf_control/scf_c_create [1.0] *
!!
!!   NAME
!!     scf_c_create
!!
!!   FUNCTION
!!     allocates and initializes an scf control object with the default values
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
!!     - Default ROKS parameters added (05.04.06,MK)
!!
!!*****************************************************************************

  SUBROUTINE scf_c_create(scf_control,error)

    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_c_create",&
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    ALLOCATE (scf_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      ! Load the default values

      scf_control%density_guess = atomic_guess
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
      scf_control%id_nr=last_scf_c_id_nr
      scf_control%eps_jacobi = 0.0_dp
      scf_control%jacobi_threshold = 1.0E-7_dp
      scf_control%write_restart_each=10
      scf_control%added_mos=0
      scf_control%max_scf_aspc = 1
      scf_control%gspace_mixing_a = 1.0_dp
      scf_control%gspace_mixing_b = 0.0_dp
      scf_control%fill_wf_history = .FALSE.
      scf_control%harris_max_scf = 1

      ! ROKS

      scf_control%roks_scheme = high_spin_roks
      scf_control%roks_f = 0.5_dp

      ! Initialize the diagonal blocks with the default ROKS parameters
      ! 0 = v)irtual, 1 = o)pen shell, 2 = c)losed shell

      scf_control%roks_parameter(0,0,1) =  1.5_dp ! avv
      scf_control%roks_parameter(0,0,2) = -0.5_dp ! bvv
      scf_control%roks_parameter(1,1,1) =  0.5_dp ! aoo
      scf_control%roks_parameter(1,1,2) =  0.5_dp ! boo
      scf_control%roks_parameter(2,2,1) = -0.5_dp ! acc
      scf_control%roks_parameter(2,2,2) =  1.5_dp ! bcc

      ! Initialize off-diagonal blocks (fixed)

      scf_control%roks_parameter(0,1,1) =  1.0_dp ! avo
      scf_control%roks_parameter(0,1,2) =  0.0_dp ! bvo
      scf_control%roks_parameter(0,2,1) =  0.5_dp ! avc
      scf_control%roks_parameter(0,2,2) =  0.5_dp ! bvc
      scf_control%roks_parameter(1,2,1) =  0.0_dp ! aoc
      scf_control%roks_parameter(1,2,2) =  1.0_dp ! boc

      ! Symmetry enforces

      scf_control%roks_parameter(1,0,1) = scf_control%roks_parameter(0,1,1) ! aov
      scf_control%roks_parameter(1,0,2) = scf_control%roks_parameter(0,1,2) ! bov
      scf_control%roks_parameter(2,0,1) = scf_control%roks_parameter(0,2,1) ! acv
      scf_control%roks_parameter(2,0,2) = scf_control%roks_parameter(0,2,2) ! bcv
      scf_control%roks_parameter(2,1,1) = scf_control%roks_parameter(1,2,1) ! aco
      scf_control%roks_parameter(2,1,2) = scf_control%roks_parameter(1,2,2) ! bco

      ! Outer SCF default settings

      scf_control%outer_scf%have_scf = .FALSE.
      scf_control%outer_scf%max_scf = 0
      scf_control%outer_scf%eps_scf = 0.0_dp
      scf_control%outer_scf%step_size = 0.0_dp
      scf_control%outer_scf%type = -1
      scf_control%outer_scf%optimizer = -1
      scf_control%outer_scf%diis_buffer_length = -1

      last_scf_c_id_nr = last_scf_c_id_nr + 1

    END IF

  END SUBROUTINE scf_c_create

!******************************************************************************
!!****f* scf_control/scf_c_retain [1.0] *
!!
!!   NAME
!!     scf_c_retain
!!
!!   FUNCTION
!!     retains the given scf_control (see cp2k/doc/ReferenceCounting.html)
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
!!*****************************************************************************

  SUBROUTINE scf_c_retain(scf_control,error)

    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "scf_c_retain", &
                                   routineP = moduleN//":"//routineN

    LOGICAL                                  :: failure

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
      CPPrecondition(scf_control%ref_count > 0,cp_failure_level,routineP,error,failure)
      scf_control%ref_count = scf_control%ref_count + 1
    END IF

  END SUBROUTINE scf_c_retain

!******************************************************************************
!!****f* scf_control/scf_c_release [1.0] *
!!
!!   NAME
!!     scf_c_release
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
!!*****************************************************************************

  SUBROUTINE scf_c_release(scf_control,error)

    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "scf_c_release", &
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    IF (ASSOCIATED(scf_control)) THEN
      CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
      scf_control%ref_count = scf_control%ref_count - 1
      IF (scf_control%ref_count < 1) THEN
        DEALLOCATE(scf_control,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
      END IF
    END IF

    NULLIFY (scf_control)

  END SUBROUTINE scf_c_release

!******************************************************************************
!!****f* scf_control/scf_c_read_parameters [1.0] *
!!
!!   NAME
!!     scf_c_read_parameters
!!
!!   FUNCTION
!!     reads the parameters of the scf section into the given scf_control
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
!!*****************************************************************************

  SUBROUTINE scf_c_read_parameters(scf_control,globenv,force_env_section,error)

    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_c_read_parameters", &
                                   routineP = moduleN//":"//routineN

    REAL(KIND=dp), DIMENSION(:), POINTER     :: parameter
    INTEGER, DIMENSION(:), POINTER           :: added_mos
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_vals
    LOGICAL                                  :: explicit, failure
    TYPE(section_vals_type), POINTER         :: harris_section, &
                                                outer_scf_section, scf_section

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY (r_vals)

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(globenv%input_file),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

       scf_section => section_vals_get_subs_vals(force_env_section,&
                                                 "DFT%SCF",error=error)
       CALL section_vals_val_get(scf_section,"nrow_block",&
            i_val=scf_control%nrow_block,error=error)
       CALL section_vals_val_get(scf_section,"ncol_block",&
            i_val=scf_control%ncol_block,error=error)
       CALL section_vals_val_get(scf_section,"OT%_SECTION_PARAMETERS_",&
            l_val=scf_control%use_ot,error=error)
       CALL section_vals_val_get(scf_section,"SCF_GUESS",&
            i_val=scf_control%density_guess,error=error)
       CALL section_vals_val_get(scf_section,"MIXING",&
            r_val=scf_control%p_mix,error=error)
       CALL section_vals_val_get(scf_section,"GSPACE_MIXING",&
                                 r_vals=r_vals,error=error)
       scf_control%gspace_mixing_a = r_vals(1)
       scf_control%gspace_mixing_b = r_vals(2)
       CALL section_vals_val_get(scf_section,"EPS_DIIS",&
            r_val=scf_control%eps_diis,error=error)
       CALL section_vals_val_get(scf_section,"eps_eigval",&
            r_val=scf_control%eps_eigval,error=error)
       CALL section_vals_val_get(scf_section,"eps_jacobi",&
            r_val=scf_control%eps_jacobi,error=error)
       CALL section_vals_val_get(scf_section,"jacobi_threshold",&
            r_val=scf_control%jacobi_threshold,error=error)
       CALL section_vals_val_get(scf_section,"cholesky",&
            l_val=scf_control%use_cholesky,error=error)
       CALL section_vals_val_get(scf_section,"eps_scf",&
            r_val=scf_control%eps_scf,error=error)
       CALL section_vals_val_get(scf_section,"level_shift",&
            r_val=scf_control%level_shift,error=error)
       CALL section_vals_val_get(scf_section,"max_diis",&
            i_val=scf_control%max_diis,error=error)
       CALL section_vals_val_get(scf_section,"max_scf",&
            i_val=scf_control%max_scf,error=error)
       CALL section_vals_val_get(scf_section,"smear",&
            r_val=scf_control%smear,error=error)
       CALL section_vals_val_get(scf_section,"work_syevx",&
            r_val=scf_control%work_syevx,error=error)
       scf_control%work_syevx = MIN(MAX(0.0_dp,scf_control%work_syevx),1.0_dp)
       CALL section_vals_val_get(scf_section,"write_restart_each",&
            i_val=scf_control%write_restart_each,error=error)

       ! Read ROKS parameters

       CALL section_vals_val_get(scf_section,"ROKS_SCHEME",&
                                 i_val=scf_control%roks_scheme,&
                                 error=error)

       IF (scf_control%roks_scheme == general_roks) THEN

         ! Read parameters for the general ROKS scheme

         CALL section_vals_val_get(scf_section,"ROKS_F",&
                                   r_val=scf_control%roks_f,&
                                   error=error)

       ELSE IF (scf_control%roks_scheme == high_spin_roks) THEN

         ! Read high-spin ROKS parameters for the diagonal block
         ! 0 = v)irtual, 1 = o)pen shell, 2 = c)losed shell

         NULLIFY (parameter)
         CALL section_vals_val_get(scf_section,"ROKS_PARAMETERS",&
                                   r_vals=parameter,&
                                   error=error)
         IF (ASSOCIATED(parameter)) THEN
           scf_control%roks_parameter(2,2,1) = parameter(1) ! acc
           scf_control%roks_parameter(2,2,2) = parameter(2) ! bcc
           scf_control%roks_parameter(1,1,1) = parameter(3) ! aoo
           scf_control%roks_parameter(1,1,2) = parameter(4) ! boo
           scf_control%roks_parameter(0,0,1) = parameter(5) ! avv
           scf_control%roks_parameter(0,0,2) = parameter(6) ! bvv
         END IF

       END IF

       ! should be moved to printkey
       CALL section_vals_val_get(scf_section,"eps_lumo",&
            r_val=scf_control%eps_lumos,error=error)
       CALL section_vals_val_get(scf_section,"max_iter_lumo",&
            i_val=scf_control%max_iter_lumos,error=error)
       CALL section_vals_val_get(scf_section,"added_mos",&
            i_vals=added_mos,error=error)
       CPPrecondition(ASSOCIATED(added_mos),cp_failure_level,routineP,error,failure)
       scf_control%added_mos=0
       IF (SIZE(added_mos)>0) THEN
          scf_control%added_mos(:)=added_mos(1)
          IF (SIZE(added_mos)>1) THEN
             scf_control%added_mos(2)=added_mos(2)
          END IF
       END IF
       CALL section_vals_val_get(scf_section,"max_scf_aspc",&
                                 i_val=scf_control%max_scf_aspc,&
                                 error=error)

       IF (scf_control%level_shift /= 0.0_dp) &
            scf_control%use_cholesky = .FALSE.
       IF (scf_control%use_ot .AND. &
            (scf_control%density_guess==core_guess)) THEN
          CALL stop_program("scf_c_read_parameters","Use GUESS {ATOMIC,RESTART,RANDOM} with OT")
       ENDIF

       outer_scf_section => section_vals_get_subs_vals(scf_section,&
                            "OUTER_SCF",error=error)

       CALL section_vals_get(outer_scf_section,explicit=explicit,error=error)
       IF (explicit) THEN
          scf_control%outer_scf%have_scf=.TRUE.
          CALL section_vals_val_get(outer_scf_section,"EPS_SCF",&
                     r_val=scf_control%outer_scf%eps_scf,error=error)
          CALL section_vals_val_get(outer_scf_section,"STEP_SIZE",&
                     r_val=scf_control%outer_scf%step_size,error=error)
          CALL section_vals_val_get(outer_scf_section,"DIIS_BUFFER_LENGTH",&
                     i_val=scf_control%outer_scf%diis_buffer_length,error=error)
          CALL section_vals_val_get(outer_scf_section,"TYPE",&
                     i_val=scf_control%outer_scf%type,error=error)
          CALL section_vals_val_get(outer_scf_section,"OPTIMIZER",&
                     i_val=scf_control%outer_scf%optimizer,error=error)
          CALL section_vals_val_get(outer_scf_section,"MAX_SCF",&
                     i_val=scf_control%outer_scf%max_scf,error=error)
          CALL section_vals_val_get(outer_scf_section,"EXTRAPOLATION_ORDER",&
                     i_val=scf_control%outer_scf%extrapolation_order,error=error)
       ELSE
          scf_control%outer_scf%have_scf=.FALSE.
       ENDIF

       harris_section => section_vals_get_subs_vals(force_env_section, &
                                                    "DFT%QS%HARRIS",&
                                                    error=error)

       CALL section_vals_val_get(harris_section, "fill_wf_history", &
            l_val=scf_control%fill_wf_history, error=error)
       CALL section_vals_val_get(harris_section, "harris_max_scf", &
            i_val=scf_control%harris_max_scf, error=error)

    END IF

  END SUBROUTINE scf_c_read_parameters

!******************************************************************************
!!****f* scf_control/scf_c_write_parameters [1.0] *
!!
!!   NAME
!!     scf_c_write_parameters
!!
!!   FUNCTION
!!     writes out the scf parameters
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
!!*** ************************************************************************

  SUBROUTINE scf_c_write_parameters(scf_control,globenv,force_env_section,error)

    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "scf_c_write_parameters", &
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: output_unit,roks_scheme
    LOGICAL                                  :: failure,roks
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: dft_section,scf_section

    ! -------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY (logger)
    logger => cp_error_get_logger(error)

    NULLIFY (scf_section)
    NULLIFY (section)

    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      scf_section => section_vals_get_subs_vals(force_env_section,"DFT%SCF",error=error)
      output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
                                         extension=".scfLog",error=error)

      IF (output_unit > 0) THEN
        IF (scf_control%max_scf < 1) RETURN
          CALL create_scf_section(section,error=error)
          keyword => section_get_keyword(section,"SCF_GUESS",error=error)
          CALL keyword_get(keyword,enum=enum,error=error)

          WRITE (UNIT=output_unit,&
                 FMT="(/,/,T2,A,T25,A,T71,A10,/,T25,56('-'),4(/,T25,A,T76,I5),"//&
                     "T25,56('-'),5(/,T25,A,T72,ES9.2),T25,56('-'),"//&
                     "6(/,T25,A,T76,F5.2),1(/,T25,A,T71,2I5))")&
            "SCF PARAMETERS",&
            "Density guess:     ",ADJUSTR(TRIM(enum_i2c(enum,scf_control%density_guess))),&
            "max_scf:           ",scf_control%max_scf,&
            "max_scf_aspc:      ",scf_control%max_scf_aspc,&
            "max_diis:          ",scf_control%max_diis,&
            "write_restart_each:",scf_control%write_restart_each,&
            "eps_scf:           ",scf_control%eps_scf,&
            "eps_diis:          ",scf_control%eps_diis,&
            "eps_eigval:        ",scf_control%eps_eigval,&
            "eps_jacobi:        ",scf_control%eps_jacobi,&
            "jacobi_threshold:  ",scf_control%jacobi_threshold,&
            "p_mix:             ",scf_control%p_mix,&
            "G-space mixing a:  ",scf_control%gspace_mixing_a,&
            "G-space mixing b:  ",scf_control%gspace_mixing_b,&
            "work_syevx:        ",scf_control%work_syevx,&
            "level_shift [a.u.]:",scf_control%level_shift,&
            "smear [a.u.]:      ",scf_control%smear,&
            "added MOs          ",scf_control%added_mos

          dft_section => section_vals_get_subs_vals(force_env_section,"DFT",&
                                                    error=error)

          CALL section_vals_val_get(dft_section,"ROKS",l_val=roks,error=error)

          IF (roks) THEN
            CALL section_vals_val_get(scf_section,"ROKS_SCHEME",&
                                      i_val=roks_scheme,error=error)
            keyword => section_get_keyword(section,"ROKS_SCHEME",error=error)
            CALL keyword_get(keyword,enum=enum,error=error)
            WRITE (UNIT=output_unit,FMT="(T25,A,/,T25,A,T71,A10)")&
              REPEAT("-",56),&
              "ROKS scheme:",ADJUSTR(TRIM(enum_i2c(enum,roks_scheme)))
            IF (roks_scheme == general_roks) THEN
              WRITE (UNIT=output_unit,FMT="(T25,A,T71,F10.6)")&
                "ROKS parameter f:",scf_control%roks_f
            ELSE IF (roks_scheme == high_spin_roks) THEN
              WRITE (UNIT=output_unit,&
                     FMT="(T25,A,6(/,T25,A,T71,F10.6))")&
                "ROKS parameters: a)lpha, b)eta; c)losed, o)pen, v)irtual",&
                "acc",scf_control%roks_parameter(2,2,1),&
                "bcc",scf_control%roks_parameter(2,2,2),&
                "aoo",scf_control%roks_parameter(1,1,1),&
                "boo",scf_control%roks_parameter(1,1,2),&
                "avv",scf_control%roks_parameter(0,0,1),&
                "bvv",scf_control%roks_parameter(0,0,2)
            END IF
          END IF

          CALL section_release(section,error=error)

          IF (scf_control%outer_scf%have_scf) THEN
            WRITE (output_unit,"(T25,56('-'),/,T25,A)") "Outer loop SCF in use "
            SELECT CASE(scf_control%outer_scf%type)
            CASE (outer_scf_none)
              WRITE (output_unit,'(T25,A)') "No variables optimised in outer loop"
            CASE (outer_scf_ddapc_constraint)
              WRITE (output_unit,'(T25,A)') "DDAPC constraint enforced"
            CASE (outer_scf_s2_constraint)
              WRITE (output_unit,'(T25,A)') "S2 constraint enforced"
            CASE DEFAULT
              CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
            END SELECT
            WRITE (output_unit,'(T25,A,T72,ES9.2)') "eps_scf",scf_control%outer_scf%eps_scf
            WRITE (output_unit,'(T25,A,T72,I9)') "max_scf",scf_control%outer_scf%max_scf
            SELECT CASE(scf_control%outer_scf%optimizer)
            CASE (outer_scf_optimizer_none)
              WRITE (output_unit,'(T25,A)') "No outer loop optimization"
            CASE (outer_scf_optimizer_sd)
              WRITE (output_unit,'(T25,A)') "Steepest descent optimization"
            CASE (outer_scf_optimizer_diis)
              WRITE (output_unit,'(T25,A)') "DIIS optimization"
              WRITE (output_unit,'(T25,A,T72,I9)') "DIIS buffer length", &
                scf_control%outer_scf%diis_buffer_length
             CASE DEFAULT
               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             WRITE (output_unit,'(T25,A,T72,ES9.2)') "step_size",scf_control%outer_scf%step_size
          ELSE
             WRITE (output_unit,"(T25,56('-'),/,T25,A)") "No outer SCF"
          ENDIF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
                                         "PRINT%PROGRAM_RUN_INFO",error=error)
    END IF

  END SUBROUTINE scf_c_write_parameters

  ! ***************************************************************************

END MODULE scf_control_types
