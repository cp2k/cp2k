!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_dftb_utils [1.0] *
!!
!!   NAME
!!     qs_dftb_utils
!!
!!   FUNCTION
!!     Working with the DFTB parameter types.
!!
!!   AUTHOR
!!     JGH (24.02.2007)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_dftb_utils 

  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: ptable,&
                                             nelem
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
                                             
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_utils'

  PUBLIC :: allocate_dftb_atom_param,&
            deallocate_dftb_atom_param,&
            get_dftb_atom_param,&
            set_dftb_atom_param,&
            write_dftb_atom_param

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_dftb_atom_param(dftb_parameter,error)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(dftb_parameter)) &
            CALL deallocate_dftb_atom_param(dftb_parameter,error)

    ALLOCATE (dftb_parameter,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    dftb_parameter % defined = .FALSE.
    dftb_parameter % name = ""
    dftb_parameter % typ = "NONE"
    dftb_parameter % z = -1
    dftb_parameter % zeff = -1.0_dp
    dftb_parameter % natorb = 0
    dftb_parameter % lmax = -1
    dftb_parameter % skself = 0.0_dp
    dftb_parameter % occupation = 0.0_dp
    dftb_parameter % eta = 0.0_dp
    dftb_parameter % energy = 0.0_dp

  END SUBROUTINE allocate_dftb_atom_param

! *****************************************************************************

  SUBROUTINE deallocate_dftb_atom_param(dftb_parameter,error)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    CPPrecondition(ASSOCIATED(dftb_parameter),cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_dftb_atom_param

! *****************************************************************************

  SUBROUTINE get_dftb_atom_param(dftb_parameter,name,typ,defined,z,zeff,natorb,&
    lmax,skself,occupation,eta,energy,cutoff,error)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name, typ
    LOGICAL, INTENT(OUT), OPTIONAL           :: defined
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff
    INTEGER, INTENT(OUT), OPTIONAL           :: natorb
    INTEGER, INTENT(OUT), OPTIONAL           :: lmax
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: skself
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: occupation
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: eta
    REAL(KIND=dp), OPTIONAL                  :: energy
    REAL(KIND=dp), OPTIONAL                  :: cutoff
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL       :: error


    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    CPPrecondition(ASSOCIATED(dftb_parameter),cp_failure_level,routineP,error,failure)

    IF (PRESENT(name)) name = dftb_parameter%name
    IF (PRESENT(typ)) typ = dftb_parameter%typ
    IF (PRESENT(defined)) defined = dftb_parameter%defined
    IF (PRESENT(z)) z = dftb_parameter%z
    IF (PRESENT(zeff)) zeff = dftb_parameter%zeff
    IF (PRESENT(natorb)) natorb = dftb_parameter%natorb
    IF (PRESENT(lmax)) lmax = dftb_parameter%lmax
    IF (PRESENT(skself)) skself = dftb_parameter%skself
    IF (PRESENT(eta)) eta = dftb_parameter%eta
    IF (PRESENT(energy)) energy = dftb_parameter%energy
    IF (PRESENT(cutoff)) cutoff = dftb_parameter%cutoff
    IF (PRESENT(occupation)) occupation = dftb_parameter%occupation

  END SUBROUTINE get_dftb_atom_param

! *****************************************************************************

  SUBROUTINE set_dftb_atom_param(dftb_parameter,name,typ,defined,z,zeff,natorb,&
    lmax,skself,occupation,eta,energy,cutoff,error)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name, typ
    LOGICAL, INTENT(IN), OPTIONAL            :: defined
    INTEGER, INTENT(IN), OPTIONAL            :: z
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff
    INTEGER, INTENT(IN), OPTIONAL            :: natorb
    INTEGER, INTENT(IN), OPTIONAL            :: lmax
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: skself
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: occupation
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: eta
    REAL(KIND=dp), OPTIONAL                  :: energy
    REAL(KIND=dp), OPTIONAL                  :: cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    CPPrecondition(ASSOCIATED(dftb_parameter),cp_failure_level,routineP,error,failure)

    IF (PRESENT(name)) dftb_parameter%name = name
    IF (PRESENT(typ)) dftb_parameter%typ = typ
    IF (PRESENT(defined)) dftb_parameter%defined = defined
    IF (PRESENT(z)) dftb_parameter%z = z
    IF (PRESENT(zeff)) dftb_parameter%zeff = zeff
    IF (PRESENT(natorb)) dftb_parameter%natorb = natorb
    IF (PRESENT(lmax)) dftb_parameter%lmax = lmax
    IF (PRESENT(skself)) dftb_parameter%skself = skself
    IF (PRESENT(eta)) dftb_parameter%eta = eta
    IF (PRESENT(occupation)) dftb_parameter%occupation = occupation
    IF (PRESENT(energy)) dftb_parameter%energy = energy
    IF (PRESENT(cutoff)) dftb_parameter%cutoff = cutoff

  END SUBROUTINE set_dftb_atom_param

! *****************************************************************************

  SUBROUTINE write_dftb_atom_param(dftb_parameter,force_env_section,error)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name, typ
    INTEGER                                  :: natorb, output_unit, z, lmax
    LOGICAL                                  :: defined, failure
    REAL(dp)                                 :: zeff
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (ASSOCIATED(dftb_parameter).AND.&
        BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
         "SUBSYS%PRINT%KINDS/BASIS_SET",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,force_env_section,"SUBSYS%PRINT%KINDS",&
            extension=".Log",error=error)

       IF (output_unit >0) THEN
          CALL get_dftb_atom_param(dftb_parameter,name=name,typ=typ,defined=defined,&
               z=z,zeff=zeff,natorb=natorb,lmax=lmax,error=error)

          WRITE (UNIT=output_unit,FMT="(/,A,T67,A14)")&
               " DFTB  parameters: ",TRIM(name)
          IF ( defined ) THEN
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                 "Effective core charge:",zeff
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
                 "Number of orbitals:",natorb
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
                 "Parameters are not defined"
         END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
            "SUBSYS%PRINT%KINDS",error=error)
    END IF

  END SUBROUTINE write_dftb_atom_param

! *****************************************************************************

END MODULE qs_dftb_utils

! *****************************************************************************
