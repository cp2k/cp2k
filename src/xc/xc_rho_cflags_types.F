!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the structure
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc_rho_cflags_types
  USE f77_blas
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_rho_cflags_types'

  PUBLIC :: xc_rho_cflags_type
  PUBLIC :: xc_rho_cflags_setall,&
            xc_rho_cflags_equal

! *****************************************************************************
!> \brief contains a flag for each component of xc_rho_set, so that you can
!>      use it to tell which components you need, which ones you need,....
!> \param rho flags for rho (the total density)
!> \param rho_spin flag for rhoa and rhob (the spin density with LSD)
!> \param drho flag for drho (the gradient of rho)
!> \param drho_spin flag for drhoa and drhob (the gradient of the spin
!>        density)
!> \param norm_drho flag for norm_drho (the norm of the gradient of rho)
!> \param norm_drho_spin flag for norm_drhoa, norm_drhob (the norm of the
!>        gradient of the spin density)
!> \param drhoa_drhob flag for drhoa_drhob (the scalar product of the
!>        gradient of the two spin densities)
!> \param rho_ 1_3: flag for rho**(1.0_dp/3.0_dp)
!> \param rho_spin_ 1_3: flag for rhoa**(1.0_dp/3.0_dp) and rhob**(1.0_dp/3.0_dp)
!> \param tau flags for the kinetic (KS) part of rho
!> \param tau_spin flags for the kinetic (KS) part of rhoa and rhob
!> \note
!>      low_level type without retain/release
!> \par History
!>      11.2003 created [fawzi]
!>      12.2008 added laplace parts [mguidon]
!> \author fawzi
! *****************************************************************************
  TYPE xc_rho_cflags_type
     LOGICAL :: rho, rho_spin, drho, drho_spin,&
          norm_drho, norm_drho_spin, drhoa_drhob,&
          rho_1_3,rho_spin_1_3, tau, tau_spin, laplace_rho, laplace_rho_spin
  END TYPE xc_rho_cflags_type

CONTAINS

! *****************************************************************************
!> \brief sets all the flags to the given value
!> \param cflags the flags to set
!> \param value the value to set
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_rho_cflags_setall(cflags, value, error)
    TYPE(xc_rho_cflags_type), INTENT(out)    :: cflags
    LOGICAL, INTENT(in)                      :: value
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_setall', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    cflags%rho=value
    cflags%rho_spin=value
    cflags%drho=value
    cflags%drho_spin=value
    cflags%norm_drho=value
    cflags%norm_drho_spin=value
    cflags%drhoa_drhob=value
    cflags%rho_1_3=value
    cflags%rho_spin_1_3=value
    cflags%tau=value
    cflags%tau_spin=value
    cflags%laplace_rho=value
    cflags%laplace_rho_spin=value
  END SUBROUTINE xc_rho_cflags_setall

! *****************************************************************************
!> \brief performs an OR on all the flags (cflags1=cflag1.or.cflags2)
!> \param cflags1 ...
!> \param cflags2 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param cflags 1: the flags to change UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param cflags 2: the flags to OR with UNMATCHED_PROCEDURE_ARGUMENT: please check 
! *****************************************************************************
  SUBROUTINE xc_rho_cflags_or(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_or', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    cflags1%rho=cflags1%rho.OR.cflags2%rho
    cflags1%rho_spin=cflags1%rho_spin.OR.cflags2%rho_spin
    cflags1%drho=cflags1%drho.OR.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.OR.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.OR.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.OR.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.OR.cflags2%drhoa_drhob
    cflags1%rho_1_3=cflags1%rho_1_3.OR.cflags2%rho_1_3
    cflags1%rho_spin_1_3=cflags1%rho_spin_1_3.OR.cflags2%rho_spin_1_3
    cflags1%tau=cflags1%tau.OR.cflags2%tau
    cflags1%tau_spin=cflags1%tau_spin.OR.cflags2%tau_spin
    cflags1%laplace_rho=cflags1%laplace_rho.OR.cflags2%laplace_rho
    cflags1%laplace_rho_spin=cflags1%laplace_rho_spin.OR.cflags2%laplace_rho_spin


  END SUBROUTINE xc_rho_cflags_or

! *****************************************************************************
!> \brief performs an AND on all the flags (cflags1=cflag1.and.cflags2)
!> \param cflags1 ...
!> \param cflags2 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param cflags 1: the flags to change UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param cflags 2: the flags to and with UNMATCHED_PROCEDURE_ARGUMENT: please check 
! *****************************************************************************
  SUBROUTINE xc_rho_cflags_and(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_and', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    cflags1%rho=cflags1%rho.AND.cflags2%rho
    cflags1%rho_spin=cflags1%rho_spin.AND.cflags2%rho_spin
    cflags1%drho=cflags1%drho.AND.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.AND.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.AND.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.AND.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.AND.cflags2%drhoa_drhob
    cflags1%rho_1_3=cflags1%rho_1_3.AND.cflags2%rho_1_3
    cflags1%rho_spin_1_3=cflags1%rho_spin_1_3.AND.cflags2%rho_spin_1_3
    cflags1%tau=cflags1%tau.AND.cflags2%tau
    cflags1%tau_spin=cflags1%tau_spin.AND.cflags2%tau_spin
    cflags1%laplace_rho=cflags1%laplace_rho.AND.cflags2%laplace_rho
    cflags1%laplace_rho_spin=cflags1%laplace_rho_spin.AND.cflags2%laplace_rho_spin

  END SUBROUTINE xc_rho_cflags_and

! *****************************************************************************
!> \brief return true if the two cflags are equal
!> \param cflags1 ...
!> \param cflags2 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param cflags 1, cflags2: the flags to compare UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \retval equal ...
! *****************************************************************************
  FUNCTION xc_rho_cflags_equal(cflags1, cflags2, error) RESULT(equal)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: equal

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_equal', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    equal=((cflags1%rho.EQV.cflags2%rho).AND.&
         (cflags1%rho_spin.EQV.cflags2%rho_spin).AND.&
         (cflags1%drho.EQV.cflags2%drho).AND.&
         (cflags1%drho_spin.EQV.cflags2%drho_spin).AND.&
         (cflags1%norm_drho.EQV.cflags2%norm_drho).AND.&
         (cflags1%norm_drho_spin.EQV.cflags2%norm_drho_spin).AND.&
         (cflags1%drhoa_drhob.EQV.cflags2%drhoa_drhob).AND.&
         (cflags1%rho_1_3.EQV.cflags2%rho_1_3).AND.&
         (cflags1%rho_spin_1_3.EQV.cflags2%rho_spin_1_3).AND.&
         (cflags1%tau.EQV.cflags2%tau).AND.&
         (cflags1%tau_spin.EQV.cflags2%tau_spin).AND.&
         (cflags1%laplace_rho.EQV.cflags2%laplace_rho).AND.&
         (cflags1%laplace_rho_spin.EQV.cflags2%laplace_rho_spin))

  END FUNCTION xc_rho_cflags_equal

! *****************************************************************************
!> \brief return true if all the true flags in cflags2 are also true in cflags1
!> \param cflags1 ...
!> \param cflags2 ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param cflags 1, cflags2: the flags to compare UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \retval contain ...
!> \note
!>      defines a partial ordering on the cflags (cflags1>cflags2)
! *****************************************************************************
  FUNCTION xc_rho_cflags_contain(cflags1, cflags2, error) RESULT(contain)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: contain

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_contain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    contain=((cflags1%rho.OR..NOT.cflags2%rho).AND.&
         (cflags1%rho_spin.OR..NOT.cflags2%rho_spin).AND.&
         (cflags1%drho.OR..NOT.cflags2%drho).AND.&
         (cflags1%drho_spin.OR..NOT.cflags2%drho_spin).AND.&
         (cflags1%norm_drho.OR..NOT.cflags2%norm_drho).AND.&
         (cflags1%norm_drho_spin.OR..NOT.cflags2%norm_drho_spin).AND.&
         (cflags1%drhoa_drhob.OR..NOT.cflags2%drhoa_drhob).AND.&
         (cflags1%rho_1_3.OR..NOT.cflags2%rho_1_3).AND.&
         (cflags1%rho_spin_1_3.OR..NOT.cflags2%rho_spin_1_3).AND.&
         (cflags1%tau.OR..NOT.cflags2%tau).AND.&
         (cflags1%tau_spin.OR..NOT.cflags2%tau_spin).AND.&
         (cflags1%laplace_rho.OR..NOT.cflags2%laplace_rho).AND.&
         (cflags1%laplace_rho_spin.OR..NOT.cflags2%laplace_rho_spin))

  END FUNCTION xc_rho_cflags_contain

! *****************************************************************************
!> \brief writes out the values of the cflags
!> \param cflags the flags to be written out
!> \param unit_nr the unit to write to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xc_rho_cflags_write(cflags,unit_nr,error)
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_cflags_write', &
      routineP = moduleN//':'//routineN

    WRITE(unit=unit_nr,fmt="(' <xc_rho_cflags>{')")
    WRITE(unit=unit_nr,fmt="('   rho=',l1,',')") cflags%rho
    WRITE(unit=unit_nr,fmt="('   rho_spin=',l1,',')") cflags%rho_spin
    WRITE(unit=unit_nr,fmt="('   drho=',l1,',')") cflags%drho
    WRITE(unit=unit_nr,fmt="('   drho_spin=',l1,',')") cflags%drho_spin
    WRITE(unit=unit_nr,fmt="('   norm_drho=',l1,',')") cflags%norm_drho
    WRITE(unit=unit_nr,fmt="('   norm_drho_spin=',l1,',')") cflags%norm_drho_spin
    WRITE(unit=unit_nr,fmt="('   drhoa_drhob=',l1,',')") cflags%drhoa_drhob
    WRITE(unit=unit_nr,fmt="('   rho_1_3=',l1,',')") cflags%rho_1_3
    WRITE(unit=unit_nr,fmt="('   rho_spin_1_3=',l1,',')") cflags%rho_spin_1_3
    WRITE(unit=unit_nr,fmt="('   tau=',l1,',')") cflags%tau
    WRITE(unit=unit_nr,fmt="('   tau_spin=',l1)") cflags%tau_spin
    WRITE(unit=unit_nr,fmt="('   laplace_rho=',l1)") cflags%laplace_rho
    WRITE(unit=unit_nr,fmt="('   laplace_rho_spin=',l1)") cflags%laplace_rho_spin
    WRITE(unit=unit_nr,fmt="(' }')")

  END SUBROUTINE xc_rho_cflags_write

END MODULE xc_rho_cflags_types
