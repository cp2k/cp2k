!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_rho_set_types [1.0] *
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
MODULE xc_rho_set_types
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE dft_types,                       ONLY: dft_control_type,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_rho_nn10,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_create_cr3d,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_cr3d,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: nn10_coeffs,&
                                             pw_nn_smear_g,&
                                             pw_spline2_deriv_g,&
                                             pw_spline2_evaluate_values_g,&
                                             pw_spline2_interpolate_values_g,&
                                             pw_spline3_deriv_g,&
                                             pw_spline3_evaluate_values_g,&
                                             pw_spline3_interpolate_values_g,&
                                             pw_spline_scale_deriv
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_p_type,&
                                             pw_transfer,&
                                             pw_type
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_rho_set_types'
  INTEGER, SAVE :: last_rho_set_id=0

  PUBLIC :: xc_rho_cflags_type
  PUBLIC :: xc_rho_set_create, xc_rho_set_retain, xc_rho_set_release,&
       xc_rho_set_update, xc_rho_set_get
  PUBLIC :: xc_rho_set_type
  PUBLIC :: xc_rho_cflags_setall, xc_rho_cflags_or, xc_rho_cflags_and,&
       xc_rho_cflags_equal, xc_rho_cflags_contain

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
!!     - rho_1_3: flag for rho**(1./3.)
!!     - rho_spin_1_3: flag for rhoa**(1./3.) and rhob**(1./3.)
!!     - kinetic: flags for the kinetic (KS) part of rho 
!!     - kinetic_spin: flags for the kinetic (KS) part of rhoa and rhob
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
          rho_1_3,rho_spin_1_3, kinetic, kinetic_spin
  END TYPE xc_rho_cflags_type
!!***
!****************************************************************************


!!****s* xc_rho_set_types/xc_rho_set_type [1.0] *
!!
!!   NAME
!!     xc_rho_set_type
!!
!!   FUNCTION
!!     represent a density, with all the representation and data needed
!!     to perform a functional evaluation
!!
!!   NOTES
!!     the use of 3d arrays is the result of trying to use only basic
!!     types (to be generic and independent from the method), and
!!     avoiding copies using the actual structure.
!!     only the part defined by local bounds is guaranteed to be present,
!!     and it is the only meaningful part.
!!
!!     add local kinetic energy?
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - id_nr: identification number (unique)
!!     - local_bounds: the part of the 3d array on which the functional is 
!!       computed
!!     - owns: which components are owned by this structure (and should
!!       be deallocated
!!     - has: which components are present and up to date
!!     - rho: the density
!!     - drho: the gradient of the density (x,y and z direction)
!!     - norm_drho: the norm of the gradient of the density
!!     - rhoa, rhob: spin alpha and beta parts of the density in the LSD case
!!     - drhoa, drhob: gradient of the spin alpha and beta parts of the density
!!       in the LSD case (x,y and z direction)
!!     - norm_drhoa, norm_drhob: norm of the gradient of rhoa and rhob
!!     - drhoa_drhob: the scalar product between the gradient of rhoa and the
!!       one of rhob
!!     - rho_1_3: rho^(1./3.)
!!     - rhoa_1_3, rhob_1_3: rhoa^(1./3.), rhob^(1./3.)
!!     - kinetic: the kinetic (KohnSham) part of rho
!!     - kinetic_a: the kinetic (KohnSham) part of rhoa
!!     - kinetic_b: the kinetic (KohnSham) part of rhob
!!
!!   AUTHOR
!!     fawzi & thomas
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi & thomas]
!!
!!   SOURCE
!***************************************************************************
  TYPE xc_rho_set_type
     INTEGER :: ref_count, id_nr
     INTEGER, DIMENSION(2,3) :: local_bounds
     REAL(kind=dp) :: rho_cutoff, drho_cutoff
     TYPE(xc_rho_cflags_type) :: owns, has
! for spin restricted systems
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: rho
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drho
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: norm_drho
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: rho_1_3
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: kinetic
! for UNrestricted systems
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: rhoa, rhob
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drhoa, drhob
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: norm_drhoa, norm_drhob
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: drhoa_drhob
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: rhoa_1_3, rhob_1_3
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: kinetic_a, kinetic_b
  END TYPE xc_rho_set_type
!!***
!****************************************************************************

CONTAINS

!!****f* xc_rho_set_types/xc_rho_cflags_setall [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_setall
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_cflags_setall(cflags, value, error)
!!       Type(xc_rho_cflags_type), Intent (OUT):: cflags
!!       Logical, Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_cflags_setall
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
    TYPE(xc_rho_cflags_type), INTENT(out) :: cflags
    LOGICAL, INTENT(in) :: value
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_cflags_setall',&
         routineP=moduleN//':'//routineN

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
    cflags%kinetic=value
    cflags%kinetic_spin=value
  END SUBROUTINE xc_rho_cflags_setall
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_or [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_or
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_cflags_or(cflags1, cflags2, error)
!!       Type(xc_rho_cflags_type), Intent (INOUT):: cflags1
!!       Type(xc_rho_cflags_type), Intent (IN):: cflags2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_cflags_or
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
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_cflags_or',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    cflags1%rho=cflags1%rho.or.cflags2%rho
    cflags1%rho_spin=cflags1%rho_spin.or.cflags2%rho_spin
    cflags1%drho=cflags1%drho.or.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.or.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.or.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.or.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.or.cflags2%drhoa_drhob
    cflags1%rho_1_3=cflags1%rho_1_3.or.cflags2%rho_1_3
    cflags1%rho_spin_1_3=cflags1%rho_spin_1_3.or.cflags2%rho_spin_1_3
    cflags1%kinetic=cflags1%kinetic.or.cflags2%kinetic
    cflags1%kinetic_spin=cflags1%kinetic_spin.or.cflags2%kinetic_spin

  END SUBROUTINE xc_rho_cflags_or
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_and [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_and
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_cflags_and(cflags1, cflags2, error)
!!       Type(xc_rho_cflags_type), Intent (INOUT):: cflags1
!!       Type(xc_rho_cflags_type), Intent (IN):: cflags2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_cflags_and
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
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_cflags_and',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    cflags1%rho=cflags1%rho.and.cflags2%rho
    cflags1%rho_spin=cflags1%rho_spin.and.cflags2%rho_spin
    cflags1%drho=cflags1%drho.and.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.and.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.and.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.and.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.and.cflags2%drhoa_drhob
    cflags1%rho_1_3=cflags1%rho_1_3.and.cflags2%rho_1_3
    cflags1%rho_spin_1_3=cflags1%rho_spin_1_3.and.cflags2%rho_spin_1_3
    cflags1%kinetic=cflags1%kinetic.and.cflags2%kinetic
    cflags1%kinetic_spin=cflags1%kinetic_spin.and.cflags2%kinetic_spin

  END SUBROUTINE xc_rho_cflags_and
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_equal [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_equal
!!
!!   SYNOPSIS
!!     Function xc_rho_cflags_equal(cflags1, cflags2, error) Result(equal)
!!       Type(xc_rho_cflags_type), Intent (INOUT):: cflags1
!!       Type(xc_rho_cflags_type), Intent (IN):: cflags2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: equal
!!     End Function xc_rho_cflags_equal
!!
!!   FUNCTION
!!     return true if the two cflags are equal
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags1, cflags2: the flags to compare
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
  FUNCTION xc_rho_cflags_equal(cflags1, cflags2, error) RESULT(equal)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    LOGICAL :: equal

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_cflags_equal',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    equal=(cflags1%rho.eqv.cflags2%rho.AND.&
         cflags1%rho_spin.eqv.cflags2%rho_spin.AND.&
         cflags1%drho.eqv.cflags2%drho.AND.&
         cflags1%drho_spin.eqv.cflags2%drho_spin.AND.&
         cflags1%norm_drho.eqv.cflags2%norm_drho.AND.&
         cflags1%norm_drho_spin.eqv.cflags2%norm_drho_spin.AND.&
         cflags1%drhoa_drhob.eqv.cflags2%drhoa_drhob.AND.&
         cflags1%rho_1_3.eqv.cflags2%rho_1_3.AND.&
         cflags1%rho_spin_1_3.EQV.cflags2%rho_spin_1_3.AND.&
         cflags1%kinetic.EQV.cflags2%kinetic.AND.&
         cflags1%kinetic_spin.eqv.cflags2%kinetic_spin)

  END FUNCTION xc_rho_cflags_equal
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_contain [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_contain
!!
!!   SYNOPSIS
!!     Function xc_rho_cflags_contain(cflags1, cflags2, error)&
!!         Result(contain)
!!       Type(xc_rho_cflags_type), Intent (INOUT):: cflags1
!!       Type(xc_rho_cflags_type), Intent (IN):: cflags2
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: contain
!!     End Function xc_rho_cflags_contain
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
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION xc_rho_cflags_contain(cflags1, cflags2, error) RESULT(contain)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    LOGICAL :: contain

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_cflags_contain',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    contain=((cflags1%rho.or..not.cflags2%rho).AND.&
         (cflags1%rho_spin.or..not.cflags2%rho_spin).AND.&
         (cflags1%drho.OR..not.cflags2%drho).AND.&
         (cflags1%drho_spin.or..not.cflags2%drho_spin).AND.&
         (cflags1%norm_drho.or..not.cflags2%norm_drho).AND.&
         (cflags1%norm_drho_spin.or..not.cflags2%norm_drho_spin).AND.&
         (cflags1%drhoa_drhob.OR..NOT.cflags2%drhoa_drhob).AND.&
         (cflags1%rho_1_3.or..not.cflags2%rho_1_3).AND.&
         (cflags1%rho_spin_1_3.or..not.cflags2%rho_spin_1_3).AND.&
         (cflags1%kinetic.or..not.cflags2%kinetic).AND.&
         (cflags1%kinetic_spin.or..not.cflags2%kinetic_spin))

  END FUNCTION xc_rho_cflags_contain
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_create [1.0] *
!!
!!   NAME
!!     xc_rho_set_create
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_set_create(rho_set, local_bounds, rho_cutoff,&
!!         drho_cutoff, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Integer, Dimension(2,3), Intent (IN):: local_bounds
!!       Real(Kind=dp), Intent (IN), Optional:: rho_cutoff, drho_cutoff
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_set_create
!!
!!   FUNCTION
!!     allocates and does (minimal) initialization of a rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the structure to allocate
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
  SUBROUTINE xc_rho_set_create(rho_set,local_bounds,rho_cutoff,drho_cutoff,&
       error)
    TYPE(xc_rho_set_type), POINTER :: rho_set
    INTEGER, DIMENSION(2,3), INTENT(in) :: local_bounds
    REAL(kind=dp), INTENT(in), OPTIONAL :: rho_cutoff, drho_cutoff
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i

    failure=.FALSE.

    ALLOCATE(rho_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       rho_set%ref_count=1
       last_rho_set_id=last_rho_set_id+1
       rho_set%id_nr=last_rho_set_id
       rho_set%rho_cutoff=EPSILON(0.0_dp)
       IF (PRESENT(rho_cutoff)) rho_set%rho_cutoff=rho_cutoff
       rho_set%drho_cutoff=EPSILON(0.0_dp)
       IF (PRESENT(drho_cutoff)) rho_set%drho_cutoff=drho_cutoff
       rho_set%local_bounds(1,:)=-HUGE(0) ! we want to crash...
       rho_set%local_bounds(1,:)=HUGE(0)
       CALL xc_rho_cflags_setall(rho_set%owns,.TRUE.,error=error)
       CALL xc_rho_cflags_setall(rho_set%has,.FALSE.,error=error)
       NULLIFY(rho_set%rho)
       DO i=1,3
          NULLIFY(rho_set%drho(i)%array)
       END DO
       NULLIFY(rho_set%rho_1_3)
       NULLIFY(rho_set%norm_drho,rho_set%rhoa,rho_set%rhob)
       DO i=1,3
          NULLIFY(rho_set%drhoa(i)%array, rho_set%drhob(i)%array)
       END DO
       NULLIFY(rho_set%norm_drhoa, rho_set%norm_drhob, &
            rho_set%drhoa_drhob,rho_set%rhoa_1_3,rho_set%rhob_1_3,&
            rho_set%kinetic,rho_set%kinetic_a,rho_set%kinetic_b)
    END IF
  END SUBROUTINE xc_rho_set_create
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_retain [1.0] *
!!
!!   NAME
!!     xc_rho_set_retain
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_set_retain(rho_set, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_set_retain
!!
!!   FUNCTION
!!     retains the given rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the object to retain
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
  SUBROUTINE xc_rho_set_retain(rho_set, error)
    TYPE(xc_rho_set_type), POINTER :: rho_set
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_retain',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(rho_set%ref_count>0,cp_failure_level,routineP,error)
       rho_set%ref_count=rho_set%ref_count+1
    END IF
  END SUBROUTINE xc_rho_set_retain
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_release [1.0] *
!!
!!   NAME
!!     xc_rho_set_release
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_set_release(rho_set, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_set_release
!!
!!   FUNCTION
!!     releases the given rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the structure to release
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
  SUBROUTINE xc_rho_set_release(rho_set, error)
    TYPE(xc_rho_set_type), POINTER :: rho_set
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat,i

    failure=.FALSE.

    IF (ASSOCIATED(rho_set)) THEN
       CPPreconditionNoFail(rho_set%ref_count>0,cp_failure_level,routineP,error)
       rho_set%ref_count=rho_set%ref_count-1
       IF (rho_set%ref_count==0) THEN
          rho_set%local_bounds(1,:)=-HUGE(0) ! we want to crash...
          rho_set%local_bounds(1,:)=HUGE(0)
          IF (rho_set%owns%rho .AND. ASSOCIATED(rho_set%rho)) THEN
             DEALLOCATE(rho_set%rho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%rho_spin) THEN
             IF (ASSOCIATED(rho_set%rhoa)) THEN
                DEALLOCATE(rho_set%rhoa, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%rhob)) THEN
                DEALLOCATE(rho_set%rhob, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          IF (rho_set%owns%rho_1_3.AND.ASSOCIATED(rho_set%rho_1_3)) THEN
             DEALLOCATE(rho_set%rho_1_3, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%drho) THEN
             DO i=1,3
                IF (ASSOCIATED(rho_set%drho(i)%array)) THEN
                   DEALLOCATE(rho_set%drho(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
             END DO
          END IF
          IF (rho_set%owns%drho_spin) THEN
             DO i=1,3
                IF (ASSOCIATED(rho_set%drhoa(i)%array)) THEN
                   DEALLOCATE(rho_set%drhoa(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
                IF (ASSOCIATED(rho_set%drhob(i)%array)) THEN
                   DEALLOCATE(rho_set%drhob(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
             END DO
          END IF
          IF (rho_set%owns%norm_drho.and.ASSOCIATED(rho_set%norm_drho)) THEN
             DEALLOCATE(rho_set%norm_drho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%norm_drho_spin) THEN
             IF (ASSOCIATED(rho_set%norm_drhoa)) THEN
                DEALLOCATE(rho_set%norm_drhoa, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%norm_drhob)) THEN
                DEALLOCATE(rho_set%norm_drhob, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          IF (rho_set%owns%drhoa_drhob.AND.ASSOCIATED(rho_set%drhoa_drhob)) THEN
             DEALLOCATE(rho_set%drhoa_drhob, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%rho_spin_1_3) THEN
             IF (ASSOCIATED(rho_set%rhoa_1_3)) THEN
                DEALLOCATE(rho_set%rhoa_1_3, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%rhob_1_3)) THEN
                DEALLOCATE(rho_set%rhob_1_3, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          IF (rho_set%owns%kinetic.AND.ASSOCIATED(rho_set%kinetic)) THEN
             DEALLOCATE(rho_set%kinetic, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%kinetic_spin) THEN
             IF (ASSOCIATED(rho_set%kinetic_a)) THEN
                DEALLOCATE(rho_set%kinetic_a, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%kinetic_b)) THEN
                DEALLOCATE(rho_set%kinetic_b, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          DEALLOCATE(rho_set, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(rho_set)
  END SUBROUTINE xc_rho_set_release
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_get [1.0] *
!!
!!   NAME
!!     xc_rho_set_get
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_set_get(rho_set, can_return_null, rho, drho,&
!!         norm_drho, rhoa, rhob, norm_drhoa, norm_drhob, drhoa_drhob,&
!!         rho_1_3, rhoa_1_3, rhob_1_3, drhoa, drhob, rho_cutoff,&
!!         drho_cutoff, kinetic, kinetic_a, kinetic_b, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Logical, Intent (IN), Optional:: can_return_null
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: rho
!!       Type(cp_3d_r_p_type), Dimension(:), Pointer, Optional:: drho
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: norm_drho
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: rhoa, rhob
!!       Type(cp_3d_r_p_type), Dimension(:), Optional, Pointer:: drhoa,&
!!         drhob
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: norm_drhoa,&
!!         norm_drhob
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: drhoa_drhob
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: rho_1_3
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: rhoa_1_3,&
!!         rhob_1_3
!!       Real(Kind=dp), Intent (OUT), Optional:: rho_cutoff, drho_cutoff
!!       Real(Kind=dp), Dimension(:,:,:), Pointer, Optional:: kinetic,&
!!         kinetic_a, kinetic_b
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_set_get
!!
!!   FUNCTION
!!     returns the various attributes of rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the object you whant info about
!!     - can_return_null: if true the object returned can be null,
!!       if false (the default) it stops with an error if a requested
!!       component is not associated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     for the other components see the attributes of xc_rho_set
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_get(rho_set, can_return_null, rho, drho, norm_drho,&
       rhoa, rhob, norm_drhoa, norm_drhob, drhoa_drhob,rho_1_3,rhoa_1_3,&
       rhob_1_3,drhoa,drhob,rho_cutoff,drho_cutoff,kinetic,kinetic_a,&
       kinetic_b,error)
    TYPE(xc_rho_set_type), POINTER :: rho_set
    LOGICAL, INTENT(in), OPTIONAL :: can_return_null
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rho
    TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER, OPTIONAL:: drho
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: norm_drho
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rhoa, rhob
    TYPE(cp_3d_r_p_type), DIMENSION(:), OPTIONAL, POINTER :: drhoa, drhob
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: norm_drhoa, norm_drhob
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: drhoa_drhob
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rho_1_3
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rhoa_1_3,rhob_1_3
    REAL(kind=dp), INTENT(out), OPTIONAL :: rho_cutoff, drho_cutoff
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: kinetic, &
         kinetic_a, kinetic_b
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_can_return_null, valid
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_get',&
         routineP=moduleN//':'//routineN
    INTEGER :: i

    failure=.FALSE.
    my_can_return_null=.FALSE.
    IF (PRESENT(can_return_null)) my_can_return_null=can_return_null

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho)) THEN
          rho => rho_set%rho
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drho)) THEN
          drho => rho_set%drho
          IF (.NOT.my_can_return_null) THEN
             DO i=1,3
                CPPrecondition(ASSOCIATED(rho_set%drho(i)%array),cp_failure_level,routineP,error,failure)
                NULLIFY(drho)
             END DO
          END IF
       END IF
       IF (PRESENT(norm_drho)) THEN
          norm_drho => rho_set%norm_drho
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drho),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhoa)) THEN
          rhoa => rho_set%rhoa
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhob)) THEN
          rhob => rho_set%rhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drhoa)) THEN
          drhoa => rho_set%drho
          IF (.NOT.my_can_return_null) THEN
             DO i=1,3
                CPAssert(ASSOCIATED(rho_set%drhoa(i)%array),cp_failure_level,routineP,error,failure)
                NULLIFY(drhoa)
             END DO
          END IF
       END IF
       IF (PRESENT(drhob)) THEN
          drho => rho_set%drhob
          IF (.NOT.my_can_return_null) THEN
             DO i=1,3
                CPPrecondition(ASSOCIATED(rho_set%drhob(i)%array),cp_failure_level,routineP,error,failure)
                NULLIFY(drhob)
             END DO
          END IF
       END IF
       IF (PRESENT(norm_drhoa)) THEN
          norm_drhoa => rho_set%norm_drhoa
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(norm_drhob)) THEN
          norm_drhob => rho_set%norm_drhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drhoa_drhob)) THEN
          drhoa_drhob => rho_set%drhoa_drhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(drhoa_drhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rho_1_3)) THEN
          rho_1_3 => rho_set%rho_1_3
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rho_1_3),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhoa_1_3)) THEN
          rhoa_1_3 => rho_set%rhoa_1_3
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhoa_1_3),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhob_1_3)) THEN
          rhob_1_3 => rho_set%rhob_1_3
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhob_1_3),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(kinetic)) THEN
          kinetic => rho_set%kinetic
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(kinetic),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(kinetic_a)) THEN
          kinetic_a => rho_set%kinetic_a
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(kinetic_a),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(kinetic_b)) THEN
          kinetic_b => rho_set%kinetic_b
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(kinetic_b),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rho_cutoff)) rho_cutoff=rho_set%rho_cutoff
       IF (PRESENT(drho_cutoff)) drho_cutoff=rho_set%drho_cutoff
    END IF
  END SUBROUTINE xc_rho_set_get
!***************************************************************************


!!****f* xc_rho_set_types/xc_rho_set_update [1.0] *
!!
!!   NAME
!!     xc_rho_set_update
!!
!!   SYNOPSIS
!!     Subroutine xc_rho_set_update(rho_set, rho_r, rho_g, needs,&
!!         dft_control, cell, pw_pool, error)
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(pw_p_type), Dimension(:), Intent (INOUT), Target:: rho_r
!!       Type(pw_p_type), Dimension(:), Pointer:: rho_g
!!       Type(xc_rho_cflags_type), Intent (IN):: needs
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cell_type), Pointer:: cell
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_rho_set_update
!!
!!   FUNCTION
!!     updates the given rho set with the density given by
!!     rho_r (and rho_g). The rho set will contain the components specified
!!     in needs
!!
!!   NOTES
!!     when 
!!
!!   ARGUMENTS
!!     - rho_set: the rho_set to update
!!     - rho_r: the new density (in r space)
!!     - rho_g: the new density (in g space, needed for smoothing and
!!       derivatives)
!!     - needs: the components of rho that are needed
!!     - dft_control: parameters (LSD/LDA, kind of derivative requested,...)
!!     - cell: cell parameters (to rescale the derivative)
!!     - pw_pool: pool for the allocation of pw and cr3d
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
  SUBROUTINE xc_rho_set_update(rho_set, rho_r, rho_g, needs,&
       dft_control,cell,pw_pool,error)
    TYPE(xc_rho_set_type), POINTER :: rho_set
    TYPE(pw_p_type),DIMENSION(:), INTENT(inout), TARGET :: rho_r
    TYPE(pw_p_type),DIMENSION(:), POINTER        :: rho_g
    TYPE(xc_rho_cflags_type), INTENT(in) :: needs
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(cell_type), POINTER :: cell
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_rho_r_local, my_rho_g_local, gradient_f
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_update',&
         routineP=moduleN//':'//routineN
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER :: stat, ispin, idir, nspins, i,j,k
    TYPE(pw_p_type), DIMENSION(2) :: my_rho_r
    TYPE(pw_p_type), DIMENSION(3,2) :: drho_r
    TYPE(pw_p_type), DIMENSION(3) :: drho_r_att
    TYPE(pw_type), POINTER :: tmp_g, my_rho_g

    failure=.FALSE.
    DO ispin=1,2
       NULLIFY(my_rho_r(ispin)%pw)
       DO idir=1,3
          NULLIFY(drho_r(idir,ispin)%pw)
       END DO
    END DO
    DO idir=1,3
       NULLIFY(drho_r_att(idir)%pw)
    END DO
    NULLIFY(tmp_g,my_rho_g)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CALL cp_assert(ALL(rho_set%local_bounds==pw_pool%pw_grid%bounds_local),&
         cp_failure_level,cp_assertion_failed,routineP,&
         "pw_pool cr3d have different size than expected",error,failure)
    nspins=dft_control%nspins
    rho_set%local_bounds=rho_r(1)%pw%pw_grid%bounds_local

! some checks
    SELECT CASE(dft_control%nspins)
    CASE(1)
       CPPrecondition(SIZE(rho_r)==1,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(1)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(1)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%norm_drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin_1_3,cp_failure_level,routineP,error,failure)
    CASE(2)
       CPPrecondition(SIZE(rho_r)==2,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(1)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(2)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(1)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(2)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

! preclean
    IF (rho_set%owns%rho) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rho, &
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%rho)
    END IF
    IF (rho_set%owns%drho) THEN
       DO idir=1,3
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drho(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
       END DO
    ELSE
       DO idir=1,3
          NULLIFY(rho_set%drho(idir)%array)
       END DO
    END IF
    IF (rho_set%owns%norm_drho) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drho,&
            accept_non_compatible=.TRUE.,error=error)
    ELSE
       NULLIFY(rho_set%norm_drho)
    END IF
    IF (rho_set%owns%rho_spin) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhoa,&
            accept_non_compatible=.TRUE., error=error)
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhob,&
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%rhoa,rho_set%rhob)
    END IF
    IF (rho_set%owns%drho_spin) THEN
       DO idir=1,3
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drhoa(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drhob(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
       END DO
    ELSE
       DO idir=1,3
          NULLIFY(rho_set%drhoa(idir)%array,rho_set%drhob(idir)%array)
       END DO
    END IF
    IF (rho_set%owns%norm_drho_spin) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drhoa,&
            accept_non_compatible=.TRUE.,error=error)
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drhob,&
            accept_non_compatible=.TRUE.,error=error)
    ELSE
       NULLIFY(rho_set%norm_drhoa, rho_set%norm_drhob)
    END IF
    IF (rho_set%owns%drhoa_drhob) THEN
       CALL pw_pool_give_back_cr3d(pw_pool, rho_set%drhoa_drhob,&
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%drhoa_drhob)
    END IF

    CALL xc_rho_cflags_setall(rho_set%has,.FALSE.,error=error)

    gradient_f=(needs%drho_spin.OR.needs%norm_drho_spin.OR.&
         needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho)

    DO ispin=1,nspins
! introduce a smoothing kernel on the density
       IF (dft_control%xc_rho_smooth_id==xc_rho_no_smooth) THEN
          my_rho_r_local=.FALSE.
          my_rho_g_local=.FALSE.
          my_rho_r(ispin)%pw => rho_r(ispin)%pw
          IF (gradient_f) THEN
             my_rho_g => rho_g(ispin)%pw
          END IF
       ELSE
          my_rho_r_local=.TRUE.
          my_rho_g_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_r(ispin)%pw,&
               use_data=REALDATA3D, in_space=REALSPACE, &
               error=error)
          CALL pw_pool_create_pw(pw_pool, my_rho_g,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
               error=error)
          CALL pw_copy(rho_g(ispin)%pw,my_rho_g)
          SELECT CASE(dft_control%xc_rho_smooth_id)
          CASE (xc_rho_spline2_smooth)
             CALL pw_spline2_evaluate_values_g(my_rho_g,error=error)
          CASE (xc_rho_spline3_smooth)
             CALL pw_spline3_evaluate_values_g(my_rho_g,error=error)
          CASE (xc_rho_nn10)
             CALL pw_nn_smear_g(my_rho_g,&
                  nn10_coeffs, error=error)
          END SELECT
          CALL pw_transfer(my_rho_g,my_rho_r(ispin)%pw)
       END IF

       IF (gradient_f) THEN ! calculate the grad of rho
! normally when you need the gradient you need the whole gradient
! (for the partial integration)
! deriv rho
          DO idir=1,3
             NULLIFY(drho_r(idir,ispin)%pw)
             CALL pw_pool_create_pw(pw_pool,drho_r(idir,ispin)%pw, &
                  use_data=REALDATA3D, in_space=REALSPACE, &
                  error=error)
          END DO
          CALL pw_pool_create_pw(pw_pool, tmp_g,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
               error=error)
          SELECT CASE(dft_control%xc_deriv_method_id)
          CASE (xc_deriv_pw)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline2)
             IF (.NOT.my_rho_g_local) THEN
                CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                     use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                     error=error)
                my_rho_g_local=.TRUE.
                CALL pw_copy(rho_g(ispin)%pw, my_rho_g)
             END IF
             CALL pw_spline2_interpolate_values_g(my_rho_g,error=error)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline2_smooth)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline3)
             IF (.NOT.my_rho_g_local) THEN
                CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                     use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                     error=error)
                CALL pw_copy(rho_g(ispin)%pw, my_rho_g)
                my_rho_g_local=.TRUE.
             END IF
             CALL pw_spline3_interpolate_values_g(my_rho_g,error=error)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline3_smooth)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

          IF (dft_control%xc_deriv_method_id/=xc_deriv_pw) THEN
             DO idir=1,3
                drho_r_att(idir)%pw => drho_r(idir,ispin)%pw
             END DO
             CALL pw_spline_scale_deriv(drho_r_att, cell=cell,&
                  error=error)
          END IF

! my_rho_g might be either rho_g or the weights of the spline
! my_rho_r contains the function that has been derived
! drgrid contains the function that should have been derived:
!        the (maybe smoothed) rho

          IF (my_rho_g_local) THEN
             my_rho_g_local=.FALSE.
             CALL pw_pool_give_back_pw(pw_pool, my_rho_g ,error=error)
          END IF
          NULLIFY(my_rho_g)
          CALL pw_pool_give_back_pw(pw_pool, tmp_g ,error=error)

       END IF
    END DO

    SELECT CASE(dft_control%nspins)
    CASE(1)
       IF (needs%rho_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%rho_1_3, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rho_1_3(i,j,k)=(my_rho_r(1)%pw%cr3d(i,j,k))**(1./3.)
                END DO
             END DO
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       IF (needs%rho) THEN
          rho_set%rho => my_rho_r(1)%pw%cr3d
          IF (my_rho_r_local) NULLIFY(my_rho_r(1)%pw%cr3d)
          rho_set%owns%rho=my_rho_r_local
          rho_set%has%rho=.TRUE.
       END IF
       IF (needs%norm_drho) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drho, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drho(i,j,k)=SQRT(&
                        drho_r(1,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,1)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
       IF (needs%drho) THEN
          DO idir=1,3
             rho_set%drho(idir)%array => drho_r(1,1)%pw%cr3d
             NULLIFY(drho_r(1,1)%pw%cr3d)
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%drho=.TRUE.
       END IF
    CASE(2)
       IF (needs%rho) THEN
! this should basically never be the case unless you use LDA functionals
! with LSD

          CALL pw_pool_create_cr3d(pw_pool,rho_set%rho,error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rho(i,j,k)=rho_r(1)%pw%cr3d(i,j,k)+&
                        rho_r(2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       IF (needs%rho_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rho_1_3,error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rho_1_3(i,j,k)=(rho_r(1)%pw%cr3d(i,j,k)+&
                        rho_r(2)%pw%cr3d(i,j,k))**(1./3.)
                END DO
             END DO
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       IF (needs%rho_spin_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rhoa_1_3,error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rhoa_1_3(i,j,k)=(rho_r(1)%pw%cr3d(i,j,k))**(1./3.)
                END DO
             END DO
          END DO
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rhob_1_3,error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rhob_1_3(i,j,k)=(rho_r(2)%pw%cr3d(i,j,k))**(1./3.)
                END DO
             END DO
          END DO
          rho_set%owns%rho_spin_1_3=.TRUE.
          rho_set%has%rho_spin_1_3=.TRUE.
       END IF
       IF (needs%rho_spin) THEN

          rho_set%rhoa => my_rho_r(1)%pw%cr3d
          IF (my_rho_r_local) NULLIFY(my_rho_r(1)%pw%cr3d)

          rho_set%rhob => my_rho_r(2)%pw%cr3d
          IF (my_rho_r_local) NULLIFY(my_rho_r(2)%pw%cr3d)

          rho_set%owns%rho_spin=my_rho_r_local
          rho_set%has%rho_spin=.TRUE.
       END IF
       IF (needs%norm_drho_spin) THEN

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhoa, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drhoa(i,j,k)=SQRT(&
                        drho_r(1,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,1)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhob, &
               error=error)
          rho_set%owns%norm_drho_spin=.TRUE.
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drhob(i,j,k)=SQRT(&
                        drho_r(1,2)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,2)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,2)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO

          rho_set%owns%norm_drho_spin=.TRUE.
          rho_set%has%norm_drho_spin=.TRUE.
       END IF
       IF (needs%drhoa_drhob) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%drhoa_drhob, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%drhoa_drhob(i,j,k)=&
                        drho_r(1,1)%pw%cr3d(i,j,k)*drho_r(1,2)%pw%cr3d(i,j,k)+&
                        drho_r(2,1)%pw%cr3d(i,j,k)*drho_r(2,2)%pw%cr3d(i,j,k)+&
                        drho_r(3,1)%pw%cr3d(i,j,k)*drho_r(3,2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%drhoa_drhob=.TRUE.
          rho_set%has%drhoa_drhob=.TRUE.
       END IF
       IF (needs%drho) THEN
! this should basically never be the case unless you use LDA functionals
! with LSD
          DO idir=1,3
             CALL pw_pool_create_cr3d(pw_pool,rho_set%drho(idir)%array,&
                  error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
             DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                   DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                      rho_set%drho(idir)%array(i,j,k)=&
                           drho_r(idir,1)%pw%cr3d(i,j,k)+&
                           drho_r(idir,2)%pw%cr3d(i,j,k)
                   END DO
                END DO
             END DO
          END DO
          rho_set%owns%drho=.TRUE.
          rho_set%has%drho=.TRUE.
       END IF
       IF (needs%drho_spin) THEN
          DO idir=1,3
             rho_set%drhoa(idir)%array => drho_r(idir,1)%pw%cr3d
             NULLIFY(drho_r(idir,1)%pw%cr3d)
             rho_set%drhob(idir)%array => drho_r(idir,2)%pw%cr3d
             NULLIFY(drho_r(idir,2)%pw%cr3d)
          END DO
          rho_set%owns%drho_spin=.TRUE.
          rho_set%has%drho_spin=.TRUE.
       END IF
    END SELECT
    IF (needs%kinetic.OR.needs%kinetic_spin) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="kinetic part not yet implemented",&
            error=error, error_level=cp_failure_level)
    END IF

! post cleanup
    DO ispin=1,nspins
       DO idir=1,3
          CALL pw_pool_give_back_pw(pw_pool, drho_r(idir,ispin)%pw, &
               accept_non_compatible=.TRUE., error=error)
       END DO
    END DO
    IF (my_rho_r_local) THEN
       DO ispin=1,nspins
          CALL pw_pool_give_back_pw(pw_pool, my_rho_r(ispin)%pw,&
               accept_non_compatible=.TRUE., error=error)
       END DO
    END IF

    CPPostcondition(xc_rho_cflags_equal(rho_set%has,needs),cp_failure_level,routineP,error,failure)

  END SUBROUTINE xc_rho_set_update
!***************************************************************************

END MODULE xc_rho_set_types
