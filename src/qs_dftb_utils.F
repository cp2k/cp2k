!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Working with the DFTB parameter types.
!> \author JGH (24.02.2007)
! *****************************************************************************
MODULE qs_dftb_utils

  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_utils'

  PUBLIC :: allocate_dftb_atom_param,&
            deallocate_dftb_atom_param,&
            get_dftb_atom_param,&
            set_dftb_atom_param,&
            write_dftb_atom_param

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param dftb_parameter ...
! *****************************************************************************
  SUBROUTINE allocate_dftb_atom_param(dftb_parameter)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(dftb_parameter)) &
            CALL deallocate_dftb_atom_param(dftb_parameter)

    ALLOCATE (dftb_parameter)

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
    dftb_parameter % xi = 0.0_dp
    dftb_parameter % di = 0.0_dp
    dftb_parameter % rcdisp = 0.0_dp
    dftb_parameter % dudq = 0.0_dp

  END SUBROUTINE allocate_dftb_atom_param

! *****************************************************************************
!> \brief ...
!> \param dftb_parameter ...
! *****************************************************************************
  SUBROUTINE deallocate_dftb_atom_param(dftb_parameter)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(dftb_parameter))
    DEALLOCATE (dftb_parameter)

  END SUBROUTINE deallocate_dftb_atom_param

! *****************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param skself ...
!> \param occupation ...
!> \param eta ...
!> \param energy ...
!> \param cutoff ...
!> \param xi ...
!> \param di ...
!> \param rcdisp ...
!> \param dudq ...
! *****************************************************************************
  SUBROUTINE get_dftb_atom_param(dftb_parameter,name,typ,defined,z,zeff,natorb,&
    lmax,skself,occupation,eta,energy,cutoff,xi,di,rcdisp,dudq)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name, typ
    LOGICAL, INTENT(OUT), OPTIONAL           :: defined
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff
    INTEGER, INTENT(OUT), OPTIONAL           :: natorb, lmax
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: skself, occupation, eta
    REAL(KIND=dp), OPTIONAL                  :: energy, cutoff, xi, di, &
                                                rcdisp, dudq

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(dftb_parameter))

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
    IF (PRESENT(xi)) xi = dftb_parameter%xi
    IF (PRESENT(di)) di = dftb_parameter%di
    IF (PRESENT(rcdisp)) rcdisp = dftb_parameter%rcdisp
    IF (PRESENT(dudq)) dudq = dftb_parameter%dudq

  END SUBROUTINE get_dftb_atom_param

! *****************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param skself ...
!> \param occupation ...
!> \param eta ...
!> \param energy ...
!> \param cutoff ...
!> \param xi ...
!> \param di ...
!> \param rcdisp ...
!> \param dudq ...
! *****************************************************************************
  SUBROUTINE set_dftb_atom_param(dftb_parameter,name,typ,defined,z,zeff,natorb,&
    lmax,skself,occupation,eta,energy,cutoff,xi,di,rcdisp,dudq)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name, typ
    LOGICAL, INTENT(IN), OPTIONAL            :: defined
    INTEGER, INTENT(IN), OPTIONAL            :: z
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff
    INTEGER, INTENT(IN), OPTIONAL            :: natorb, lmax
    REAL(KIND=dp), DIMENSION(0:3), OPTIONAL  :: skself, occupation, eta
    REAL(KIND=dp), OPTIONAL                  :: energy, cutoff, xi, di, &
                                                rcdisp, dudq

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(dftb_parameter))

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
    IF (PRESENT(xi)) dftb_parameter%xi = xi
    IF (PRESENT(di)) dftb_parameter%di = di
    IF (PRESENT(rcdisp)) dftb_parameter%rcdisp = rcdisp
    IF (PRESENT(dudq)) dftb_parameter%dudq = dudq

  END SUBROUTINE set_dftb_atom_param

! *****************************************************************************
!> \brief ...
!> \param dftb_parameter ...
!> \param subsys_section ...
! *****************************************************************************
  SUBROUTINE write_dftb_atom_param(dftb_parameter,subsys_section)

    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_dftb_atom_param', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name, typ
    INTEGER                                  :: lmax, natorb, output_unit, z
    LOGICAL                                  :: defined
    REAL(dp)                                 :: zeff
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()
    IF (ASSOCIATED(dftb_parameter).AND.&
        BTEST(cp_print_key_should_output(logger%iter_info,subsys_section,&
         "PRINT%KINDS/POTENTIAL"),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%KINDS",&
            extension=".Log")

       IF (output_unit >0) THEN
          CALL get_dftb_atom_param(dftb_parameter,name=name,typ=typ,defined=defined,&
               z=z,zeff=zeff,natorb=natorb,lmax=lmax)

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
       CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
            "PRINT%KINDS")
    END IF

  END SUBROUTINE write_dftb_atom_param

END MODULE qs_dftb_utils

