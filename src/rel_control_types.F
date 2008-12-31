!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief parameters that control a relativistic calculation 
!> \par History
!>      09.2002 created [fawzi] (as scf_control_types.F)
!>      10.2008 modifed for relativistic control types (Jens Thar)
!> \author Fawzi Mohamed 
! *****************************************************************************
MODULE rel_control_types

  USE f77_blas
  USE input_constants,                 ONLY: rel_none,&
                                             rel_pot_full,&
                                             rel_trans_full
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rel_control_types'
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

  INTEGER, SAVE :: last_rel_c_id_nr = 0

  ! Public data types

  PUBLIC :: rel_control_type 
       !     smear_type  outer_scf_control_type,&

  ! Public subroutines

  PUBLIC :: rel_c_create,&
            rel_c_read_parameters,&
            rel_c_release,&
            rel_c_retain
!           rel_c_write_parameters

! *****************************************************************************
!> \brief contains the parameters needed by a relativistic calculation 
!> \param method used relativistic method (NONE, DKH)
!> \param DKH_order order of the DKH transformation (2,3) 
!> \param transformation used blocks of the full matrix (FULL, MOLECULE, ATOM)
!> \param z_cutoff considered range of the Coulomb interaction
!> \param potential nuclear electron Coulomb potential (FULL, ERFC)
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT] 
!> \author Fawzi Mohamed 
! *****************************************************************************
  TYPE rel_control_type
    INTEGER                               :: rel_method 
    INTEGER                               :: rel_DKH_order
    INTEGER                               :: rel_DKH_transformation
    INTEGER                               :: rel_z_cutoff
    INTEGER                               :: rel_potential
    INTEGER                               :: ref_count
  END TYPE rel_control_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes an rel control object with the default values
!> \param rel_control the object to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      - Default ROKS parameters added (05.04.06,MK)
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed 
! *****************************************************************************
  SUBROUTINE rel_c_create(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE (rel_control,STAT=stat)
    CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      ! Load the default values

      rel_control%rel_method = rel_none 
      rel_control%rel_DKH_order = 2
      rel_control%rel_DKH_transformation = rel_trans_full
      rel_control%rel_z_cutoff = 1
      rel_control%rel_potential = rel_pot_full 
      rel_control%ref_count = 1 

    END IF

  END SUBROUTINE rel_c_create

! *****************************************************************************
!> \brief retains the given rel_control (see cp2k/doc/ReferenceCounting.html)
!> \param rel_control the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed 
! *****************************************************************************
  SUBROUTINE rel_c_retain(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rel_control),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
      CPPrecondition(rel_control%ref_count > 0,cp_failure_level,routineP,error,failure)
      rel_control%ref_count = rel_control%ref_count + 1
    END IF

  END SUBROUTINE rel_c_retain

! *****************************************************************************
!> \brief releases the given rel_control (see cp2k/doc/ReferenceCounting.html)
!> \param rel_control the object to free
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      at the moment does nothing
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed 
! *****************************************************************************
  SUBROUTINE rel_c_release(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(rel_control)) THEN
      CPPrecondition(rel_control%ref_count>0,cp_failure_level,routineP,error,failure)
      rel_control%ref_count = rel_control%ref_count - 1
      IF (rel_control%ref_count < 1) THEN
        DEALLOCATE(rel_control,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
      END IF
    END IF

    NULLIFY (rel_control)

  END SUBROUTINE rel_c_release

! *****************************************************************************
!> \brief reads the parameters of the relativistic section into the given rel_control
!> \param rel_control the object that wil contain the values read
!> \param error controls log and error handling
!> \par History
!>      05.2001 created [Matthias] for scf_control_type
!>      09.2002 created separated scf_control type [fawzi]
!>      10.2008 copied to rel_control_type [JT]
!> \author Matthias Krack 
! *****************************************************************************
  SUBROUTINE rel_c_read_parameters(rel_control,dft_section,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_read_parameters', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: rel_section

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rel_control),cp_failure_level,routineP,error,failure)
    CPPrecondition((rel_control%ref_count > 0),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      rel_section => section_vals_get_subs_vals(dft_section,"RELATIVISTIC",&
                                                error=error)
      CALL section_vals_val_get(rel_section,"method",&
                                i_val=rel_control%rel_method,&
                                error=error)
      CALL section_vals_val_get(rel_section,"DKH_order",&
                                i_val=rel_control%rel_DKH_order,&
                                error=error)
      CALL section_vals_val_get(rel_section,"transformation",&
                                i_val=rel_control%rel_DKH_transformation,&
                                error=error)
      CALL section_vals_val_get(rel_section,"z_cutoff",&
                                i_val=rel_control%rel_z_cutoff,&
                                error=error)
      CALL section_vals_val_get(rel_section,"potential",&
                                i_val=rel_control%rel_potential,&
                                error=error)

    END IF 

  END SUBROUTINE rel_c_read_parameters

! *****************************************************************************
!> \brief writes out the rel parameters
!> \param ref_control the object you want to print
!> \param error controls log and error handling
!> \par History
!>      05.2001 created [Matthias] for scf_control_type
!>      09.2002 created separated scf_control type [fawzi]
!>      10.2008 copied to rel_control_type [JT]
!> \author Matthias Krack 
! *****************************************************************************
! SUBROUTINE rel_c_write_parameters(rel_control,error)
!
!   TYPE(rel_control_type), POINTER          :: rel_control
!   TYPE(cp_error_type), INTENT(inout)       :: error
!
!   CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_write_parameters', &
!     routineP = moduleN//':'//routineN
!
!   INTEGER                                  :: output_unit
!   LOGICAL                                  :: failure
!   REAL(KIND=dp)                            :: elec_temp
!   TYPE(cp_logger_type), POINTER            :: logger
!   TYPE(enumeration_type), POINTER          :: enum
!   TYPE(keyword_type), POINTER              :: keyword
!   TYPE(section_type), POINTER              :: section
!   TYPE(section_vals_type), POINTER         :: dft_section
!
!   failure = .FALSE.
!
!   NULLIFY (logger)
!   logger => cp_error_get_logger(error)
!
!   NULLIFY (scf_section)
!   NULLIFY (section)
!
!   CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
!   CPPrecondition(scf_control%ref_count>0,cp_failure_level,routineP,error,failure)
!
!   IF (.NOT.failure) THEN
!
!     scf_section => section_vals_get_subs_vals(force_env_section,"DFT%SCF",error=error)
!     output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
!                                        extension=".scfLog",error=error)
!
!     IF (output_unit > 0) THEN
!
!       IF (scf_control%max_scf > 0) THEN
!
!         CALL create_scf_section(section,error=error)
!
!         keyword => section_get_keyword(section,"SCF_GUESS",error=error)
!         CALL keyword_get(keyword,enum=enum,error=error)
!
!         WRITE (UNIT=output_unit,&
!                FMT="(/,/,T2,A,T25,A,T51,A30,/,T25,56('-'),3(/,T25,A,T76,I5),/,"//&
!                    "T25,56('-'),6(/,T25,A,T72,ES9.2),/,T25,56('-'),"//&
!                    "5(/,T25,A,T76,F5.2),1(/,T25,A,T71,2I5))")&
!           "SCF PARAMETERS",&
!           "Density guess:     ",ADJUSTR(TRIM(enum_i2c(enum,scf_control%density_guess,error=error))),&
!           "max_scf:           ",scf_control%max_scf,&
!           "max_scf_history:   ",scf_control%max_scf_hist,&
!           "max_diis:          ",scf_control%max_diis,&
!           "eps_scf:           ",scf_control%eps_scf,&
!           "eps_scf_history:   ",scf_control%eps_scf_hist,&
!           "eps_diis:          ",scf_control%eps_diis,&
!           "eps_eigval:        ",scf_control%eps_eigval,&
!           "eps_jacobi:        ",scf_control%eps_jacobi,&
!           "jacobi_threshold:  ",scf_control%jacobi_threshold,&
!           "p_mix:             ",scf_control%p_mix,&
!           "G-space mixing a:  ",scf_control%gspace_mixing_a,&
!           "G-space mixing b:  ",scf_control%gspace_mixing_b,&
!           "work_syevx:        ",scf_control%work_syevx,&
!           "level_shift [a.u.]:",scf_control%level_shift,&
!           "added MOs          ",scf_control%added_mos
!
!         IF (scf_control%smear%do_smear) THEN
!           keyword => section_get_keyword(section,"SMEAR%METHOD",error=error)
!           CALL keyword_get(keyword,enum=enum,error=error)
!           WRITE (UNIT=output_unit,FMT="(T25,A,/,T25,A,T51,A30)")&
!             REPEAT("-",56),&
!             "Smear method:      ",ADJUSTR(TRIM(enum_i2c(enum,scf_control%smear%method,error=error)))
!           SELECT CASE (scf_control%smear%method)
!           CASE (smear_fermi_dirac)
!             elec_temp = cp_unit_from_cp2k(scf_control%smear%electronic_temperature,&
!                                           "K",error=error)
!             WRITE (UNIT=output_unit,FMT="(T25,A,T61,F20.1)")&
!               "Electronic temperature [K]:",elec_temp
!             WRITE (UNIT=output_unit,FMT="(T25,A,T71,ES10.2)")&
!               "Electronic temperature [a.u.]:",scf_control%smear%electronic_temperature,&
!               "Accuracy threshold:",scf_control%smear%eps_fermi_dirac
!           CASE (smear_energy_window)
!             WRITE (UNIT=output_unit,FMT="(T25,A,T71,F10.6)")&
!               "Smear window [a.u.]:       ",scf_control%smear%window_size
!           END SELECT
!         END IF
!         dft_section => section_vals_get_subs_vals(force_env_section,"DFT",&
!                                                   error=error)
!
!         CALL section_vals_val_get(dft_section,"ROKS",l_val=roks,error=error)
!
!         IF (roks.AND.(.NOT.scf_control%use_ot)) THEN
!           CALL section_vals_val_get(scf_section,"ROKS_SCHEME",&
!                                     i_val=roks_scheme,error=error)
!           keyword => section_get_keyword(section,"ROKS_SCHEME",error=error)
!           CALL keyword_get(keyword,enum=enum,error=error)
!           WRITE (UNIT=output_unit,FMT="(T25,A,/,T25,A,T51,A30)")&
!             REPEAT("-",56),&
!             "ROKS scheme:",ADJUSTR(TRIM(enum_i2c(enum,roks_scheme,error=error)))
!           SELECT CASE (roks_scheme)
!           CASE (general_roks)
!             WRITE (UNIT=output_unit,FMT="(T25,A,T71,F10.6)")&
!               "ROKS parameter f:",scf_control%roks_f
!           CASE (high_spin_roks)
!             WRITE (UNIT=output_unit,&
!                    FMT="(T25,A,6(/,T25,A,T71,F10.6))")&
!               "ROKS parameters: a)lpha, b)eta; c)losed, o)pen, v)irtual",&
!               "acc",scf_control%roks_parameter(2,2,1),&
!               "bcc",scf_control%roks_parameter(2,2,2),&
!               "aoo",scf_control%roks_parameter(1,1,1),&
!               "boo",scf_control%roks_parameter(1,1,2),&
!               "avv",scf_control%roks_parameter(0,0,1),&
!               "bvv",scf_control%roks_parameter(0,0,2)
!           END SELECT
!         END IF
!
!         CALL section_release(section,error=error)
!
!         IF (scf_control%outer_scf%have_scf) THEN
!           WRITE (output_unit,"(T25,56('-'),/,T25,A)") "Outer loop SCF in use "
!           SELECT CASE(scf_control%outer_scf%type)
!           CASE (outer_scf_none)
!             WRITE (output_unit,'(T25,A)') "No variables optimised in outer loop"
!           CASE (outer_scf_ddapc_constraint)
!             WRITE (output_unit,'(T25,A)') "DDAPC constraint enforced"
!           CASE (outer_scf_s2_constraint)
!             WRITE (output_unit,'(T25,A)') "S2 constraint enforced"
!           CASE (outer_scf_becke_constraint)
!             WRITE (output_unit,'(T25,A)') "Becke weight population constraint enforced"
!           CASE (outer_scf_scp)
!             WRITE (output_unit,'(T25,A)') "SCP optimization with outer loop enforced"
!           CASE DEFAULT
!             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!           END SELECT
!           WRITE (output_unit,'(T25,A,T72,ES9.2)') "eps_scf",scf_control%outer_scf%eps_scf
!           WRITE (output_unit,'(T25,A,T72,I9)') "max_scf",scf_control%outer_scf%max_scf
!           SELECT CASE(scf_control%outer_scf%optimizer)
!           CASE (outer_scf_optimizer_none)
!             WRITE (output_unit,'(T25,A)') "No outer loop optimization"
!           CASE (outer_scf_optimizer_sd)
!             WRITE (output_unit,'(T25,A)') "Steepest descent optimization"
!           CASE (outer_scf_optimizer_bisect)
!             WRITE (output_unit,'(T25,A)') "Gradient bisection"
!             WRITE (output_unit,'(T25,A,T72,I9)') "bisect_trust_count",scf_control%outer_scf%bisect_trust_count
!           CASE (outer_scf_optimizer_diis)
!             WRITE (output_unit,'(T25,A)') "DIIS optimization"
!             WRITE (output_unit,'(T25,A,T72,I9)') "DIIS buffer length", &
!               scf_control%outer_scf%diis_buffer_length
!           CASE DEFAULT
!             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!           END SELECT
!           WRITE (output_unit,'(T25,A,T72,ES9.2)') "step_size",scf_control%outer_scf%step_size
!         ELSE
!           WRITE (output_unit,"(T25,56('-'),/,T25,A)") "No outer SCF"
!         END IF
!
!       END IF ! max_scf > 0
!
!     END IF ! output_unit > 0
!
!     CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
!                                       "PRINT%PROGRAM_RUN_INFO",error=error)
!
!   END IF ! not failure
!
! END SUBROUTINE scf_c_write_parameters

END MODULE rel_control_types
