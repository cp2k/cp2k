!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/xc_rho_cflags_types [1.0] *
!!
!!   NAME
!!     xc_rho_set_types
!!
!!   FUNCTION
!!     contains the structure 
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_rho_cflags_types
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_rho_cflags_types'

  PUBLIC :: xc_rho_cflags_type
  PUBLIC :: xc_rho_cflags_setall, xc_rho_cflags_or, xc_rho_cflags_and,&
       xc_rho_cflags_equal, xc_rho_cflags_contain, xc_rho_cflags_write

!!***
  !****************************************************************************

!!****s* xc_rho_set_types/xc_rho_cflags_type [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_type
!!
!!   FUNCTION
!!     contains a flag for each component of xc_rho_set, so that you can
!!     use it to tell which components you need, which ones you need,....
!!
!!   NOTES
!!     low_level type without retain/release
!!
!!   ATTRIBUTES
!!     - rho: flags for rho (the total density)
!!     - rho_spin: flag for rhoa and rhob (the spin density with LSD)
!!     - drho: flag for drho (the gradient of rho)
!!     - drho_spin: flag for drhoa and drhob (the gradient of the spin
!!       density)
!!     - norm_drho: flag for norm_drho (the norm of the gradient of rho)
!!     - norm_drho_spin: flag for norm_drhoa, norm_drhob (the norm of the
!!       gradient of the spin density)
!!     - drhoa_drhob: flag for drhoa_drhob (the scalar product of the
!!       gradient of the two spin densities)
!!     - rho_1_3: flag for rho**(1.0_dp/3.0_dp)
!!     - rho_spin_1_3: flag for rhoa**(1.0_dp/3.0_dp) and rhob**(1.0_dp/3.0_dp)
!!     - tau: flags for the kinetic (KS) part of rho 
!!     - tau_spin: flags for the kinetic (KS) part of rhoa and rhob
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE xc_rho_cflags_type
     LOGICAL :: rho, rho_spin, drho, drho_spin,&
          norm_drho, norm_drho_spin, drhoa_drhob,&
          rho_1_3,rho_spin_1_3, tau, tau_spin
  END TYPE xc_rho_cflags_type
!!***
  !****************************************************************************


CONTAINS

!!****f* xc_rho_set_types/xc_rho_cflags_setall [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_setall
!!
!!   FUNCTION
!!     sets all the flags to the given value
!!
!!   ARGUMENTS
!!     - cflags: the flags to set
!!     - value: the value to set
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_setall(cflags, value, error)
    TYPE(xc_rho_cflags_type), INTENT(out)    :: cflags
    LOGICAL, INTENT(in)                      :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
  END SUBROUTINE xc_rho_cflags_setall
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_or [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_or
!!
!!   FUNCTION
!!     performs an OR on all the flags (cflags1=cflag1.or.cflags2)
!!
!!   ARGUMENTS
!!     - cflags1: the flags to change
!!     - cflags2: the flags to OR with
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_or(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

  END SUBROUTINE xc_rho_cflags_or
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_and [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_and
!!
!!   FUNCTION
!!     performs an AND on all the flags (cflags1=cflag1.and.cflags2)
!!
!!   ARGUMENTS
!!     - cflags1: the flags to change
!!     - cflags2: the flags to and with
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_and(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

  END SUBROUTINE xc_rho_cflags_and
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_equal [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_equal
!!
!!   FUNCTION
!!     return true if the two cflags are equal
!!
!!   ARGUMENTS
!!     - cflags1, cflags2: the flags to compare
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  FUNCTION xc_rho_cflags_equal(cflags1, cflags2, error) RESULT(equal)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
         (cflags1%tau_spin.EQV.cflags2%tau_spin))

  END FUNCTION xc_rho_cflags_equal
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_contain [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_contain
!!
!!   FUNCTION
!!     return true if all the true flags in cflags2 are also true in cflags1
!!
!!   NOTES
!!     defines a partial ordering on the cflags (cflags1>cflags2)
!!
!!   ARGUMENTS
!!     - cflags1, cflags2: the flags to compare
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  FUNCTION xc_rho_cflags_contain(cflags1, cflags2, error) RESULT(contain)
    TYPE(xc_rho_cflags_type), INTENT(inout)  :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags2
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
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
         (cflags1%tau_spin.OR..NOT.cflags2%tau_spin))

  END FUNCTION xc_rho_cflags_contain
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_write *
!!
!!   NAME
!!     xc_rho_cflags_write
!!
!!   FUNCTION
!!     writes out the values of the cflags
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags: the flags to be written out
!!     - unit_nr: the unit to write to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_write(cflags,unit_nr,error)
    TYPE(xc_rho_cflags_type), INTENT(in)     :: cflags
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    WRITE(unit=unit_nr,fmt="(' }')")

  END SUBROUTINE xc_rho_cflags_write
  !***************************************************************************

END MODULE xc_rho_cflags_types
