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
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_deriv_fd,&
                                             xc_deriv_nn10_smooth,&
                                             xc_deriv_nn50_smooth,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_rho_nn10,&
                                             xc_rho_nn50,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_pool_types,                   ONLY: pw_pool_create_cr3d,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_cr3d,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: nn10_coeffs,&
                                             nn10_deriv_coeffs,&
                                             nn50_coeffs,&
                                             nn50_deriv_coeffs,&
                                             pw_nn_deriv_r,&
                                             pw_nn_smear_r,&
                                             pw_spline2_deriv_g,&
                                             pw_spline2_interpolate_values_g,&
                                             pw_spline3_deriv_g,&
                                             pw_spline3_interpolate_values_g,&
                                             pw_spline_scale_deriv,&
                                             spline2_coeffs,&
                                             spline2_deriv_coeffs,&
                                             spline3_coeffs,&
                                             spline3_deriv_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_derive_fd,&
                                             pw_p_type,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_rho_set_types'
  INTEGER, SAVE :: last_rho_set_id=0

  PUBLIC :: xc_rho_cflags_type, zero_core, zero_core_force
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
!!     - rho_1_3: flag for rho**(1.D0/3.D0)
!!     - rho_spin_1_3: flag for rhoa**(1.D0/3.D0) and rhob**(1.D0/3.D0)
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
!!     - rho_1_3: rho^(1.D0/3.D0)
!!     - rhoa_1_3, rhob_1_3: rhoa^(1.D0/3.D0), rhob^(1.D0/3.D0)
!!     - tau: the kinetic (KohnSham) part of rho
!!     - tau_a: the kinetic (KohnSham) part of rhoa
!!     - tau_b: the kinetic (KohnSham) part of rhob
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
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: tau
     ! for UNrestricted systems
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: rhoa, rhob
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drhoa, drhob
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: norm_drhoa, norm_drhob
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: drhoa_drhob
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: rhoa_1_3, rhob_1_3
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: tau_a, tau_b
  END TYPE xc_rho_set_type
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

!!****f* xc_rho_set_types/xc_rho_set_create [1.0] *
!!
!!   NAME
!!     xc_rho_set_create
!!
!!   FUNCTION
!!     allocates and does (minimal) initialization of a rho_set
!!
!!   ARGUMENTS
!!     - rho_set: the structure to allocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_create(rho_set,local_bounds,rho_cutoff,drho_cutoff,&
       error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, DIMENSION(2, 3), INTENT(in)     :: local_bounds
    REAL(kind=dp), INTENT(in), OPTIONAL      :: rho_cutoff, drho_cutoff
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
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
       rho_set%local_bounds=local_bounds
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
            rho_set%tau,rho_set%tau_a,rho_set%tau_b)
    END IF
  END SUBROUTINE xc_rho_set_create
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_retain [1.0] *
!!
!!   NAME
!!     xc_rho_set_retain
!!
!!   FUNCTION
!!     retains the given rho_set
!!
!!   ARGUMENTS
!!     - rho_set: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_retain(rho_set, error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
!!   FUNCTION
!!     releases the given rho_set
!!
!!   ARGUMENTS
!!     - rho_set: the structure to release
!!     - pw_pool: the plae where to give back the arrays
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_release(rho_set, pw_pool, error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(rho_set)) THEN
       CPPreconditionNoFail(rho_set%ref_count>0,cp_failure_level,routineP,error)
       rho_set%ref_count=rho_set%ref_count-1
       IF (rho_set%ref_count==0) THEN
          IF (PRESENT(pw_pool)) THEN
             IF (ASSOCIATED(pw_pool)) THEN
                rho_set%ref_count=1
                CALL xc_rho_set_clean(rho_set,pw_pool,error)
                rho_set%ref_count=0
             ELSE
                CPPrecondition(.FALSE.,cp_warning_level,routineP,error,failure)
             END IF
          END IF

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
          IF (rho_set%owns%rho_spin) THEN
             IF (ASSOCIATED(rho_set%rhoa_1_3)) THEN
                DEALLOCATE(rho_set%rhoa_1_3, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%rhob_1_3)) THEN
                DEALLOCATE(rho_set%rhob_1_3, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
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
          IF (rho_set%owns%norm_drho.AND.ASSOCIATED(rho_set%norm_drho)) THEN
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
          IF (rho_set%owns%tau.AND.ASSOCIATED(rho_set%tau)) THEN
             DEALLOCATE(rho_set%tau, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%tau_spin) THEN
             IF (ASSOCIATED(rho_set%tau_a)) THEN
                DEALLOCATE(rho_set%tau_a, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%tau_b)) THEN
                DEALLOCATE(rho_set%tau_b, stat=stat)
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
!!   FUNCTION
!!     returns the various attributes of rho_set
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
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_get(rho_set, can_return_null, rho, drho, norm_drho,&
       rhoa, rhob, norm_drhoa, norm_drhob, drhoa_drhob,rho_1_3,rhoa_1_3,&
       rhob_1_3,drhoa,drhob,rho_cutoff,drho_cutoff,tau,tau_a,&
       tau_b,local_bounds,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(in), OPTIONAL            :: can_return_null
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: rho
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: drho
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: norm_drho, rhoa, rhob, &
                                                norm_drhoa, norm_drhob, &
                                                drhoa_drhob, rho_1_3, &
                                                rhoa_1_3, rhob_1_3
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: drhoa, drhob
    REAL(kind=dp), INTENT(out), OPTIONAL     :: rho_cutoff, drho_cutoff
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: tau, tau_a, tau_b
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: local_bounds
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_can_return_null

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
          drhoa => rho_set%drhoa
          IF (.NOT.my_can_return_null) THEN
             DO i=1,3
                CPAssert(ASSOCIATED(rho_set%drhoa(i)%array),cp_failure_level,routineP,error,failure)
             END DO
          END IF
       END IF
       IF (PRESENT(drhob)) THEN
          drhob => rho_set%drhob
          IF (.NOT.my_can_return_null) THEN
             DO i=1,3
                CPPrecondition(ASSOCIATED(rho_set%drhob(i)%array),cp_failure_level,routineP,error,failure)
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
       IF (PRESENT(tau)) THEN
          tau => rho_set%tau
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(tau),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(tau_a)) THEN
          tau_a => rho_set%tau_a
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(tau_a),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(tau_b)) THEN
          tau_b => rho_set%tau_b
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(tau_b),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rho_cutoff)) rho_cutoff=rho_set%rho_cutoff
       IF (PRESENT(drho_cutoff)) drho_cutoff=rho_set%drho_cutoff
       IF (PRESENT(local_bounds)) local_bounds => rho_set%local_bounds
    END IF
  END SUBROUTINE xc_rho_set_get
  !***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_clean [1.0] *
!!
!!   NAME
!!     xc_rho_set_clean
!!
!!   FUNCTION
!!     cleans (releases) most of the data stored in the rho_set giving back
!!     what it can to the pw_pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_set: the rho_set to be cleaned
!!     - pw_pool: place to give back 3d arrays,...
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_clean(rho_set,pw_pool,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_clean', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       IF (rho_set%owns%rho) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rho, &
               accept_non_compatible=.TRUE., error=error)
       ELSE
          NULLIFY(rho_set%rho)
       END IF
       IF (rho_set%owns%rho_1_3) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rho_1_3,&
               accept_non_compatible=.TRUE.,error=error)
       ELSE
          NULLIFY(rho_set%rho_1_3)
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
       IF (rho_set%owns%tau) THEN
          CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau,&
               accept_non_compatible=.TRUE., error=error)
       ELSE
          NULLIFY(rho_set%tau)
       END IF
       IF (rho_set%owns%rho_spin) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhoa,&
               accept_non_compatible=.TRUE., error=error)
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhob,&
               accept_non_compatible=.TRUE., error=error)
       ELSE
          NULLIFY(rho_set%rhoa,rho_set%rhob)
       END IF
       IF (rho_set%owns%rho_spin_1_3) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhoa_1_3,&
               accept_non_compatible=.TRUE.,error=error)
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhob_1_3,&
               accept_non_compatible=.TRUE.,error=error)
       ELSE
          NULLIFY(rho_set%rhoa_1_3,rho_set%rhob_1_3)
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
       IF (rho_set%owns%tau_spin) THEN
          CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau_a,&
               accept_non_compatible=.TRUE., error=error)
          CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau_b,&
               accept_non_compatible=.TRUE., error=error)
       ELSE
          NULLIFY(rho_set%tau_a,rho_set%tau_b)
       END IF

       CALL xc_rho_cflags_setall(rho_set%has,.FALSE.,error=error)
       CALL xc_rho_cflags_setall(rho_set%owns,.FALSE.,error=error)

    END IF
  END SUBROUTINE xc_rho_set_clean
  !***************************************************************************


!!****f* xc_rho_set_types/xc_rho_set_update [1.0] *
!!
!!   NAME
!!     xc_rho_set_update
!!
!!   FUNCTION
!!     updates the given rho set with the density given by
!!     rho_r (and rho_g). The rho set will contain the components specified
!!     in needs
!!
!!   ARGUMENTS
!!     - rho_set: the rho_set to update
!!     - rho_r: the new density (in r space)
!!     - rho_g: the new density (in g space, needed for some
!!       derivatives)
!!     - needs: the components of rho that are needed
!!     - dft_control: parameters (LSD/LDA, kind of derivative requested,...)
!!     - cell: cell parameters (to rescale the derivative)
!!     - pw_pool: pool for the allocation of pw and cr3d
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs,&
       dft_control,cell,pw_pool,atoms,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    TYPE(xc_rho_cflags_type), INTENT(in)     :: needs
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_update', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f13 = (1.0_dp/3.0_dp)

    INTEGER                                  :: i, idir, ispin, j, k, nspins
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: failure, gradient_f, &
                                                my_rho_g_local, &
                                                my_rho_r_local, needs_rho_g
    REAL(kind=dp)                            :: rho_cutoff
    TYPE(pw_p_type), DIMENSION(2)            :: my_rho_r
    TYPE(pw_p_type), DIMENSION(3)            :: drho_r_att
    TYPE(pw_p_type), DIMENSION(3, 2)         :: drho_r
    TYPE(pw_type), POINTER                   :: my_rho_g, tmp_g

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
    rho_cutoff=0.5*rho_set%rho_cutoff

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

    CALL xc_rho_set_clean(rho_set,pw_pool=pw_pool,error=error)

    gradient_f=(needs%drho_spin.OR.needs%norm_drho_spin.OR.&
         needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho)
    needs_rho_g=(dft_control%xc_deriv_method_id==xc_deriv_spline3.OR.&
         dft_control%xc_deriv_method_id==xc_deriv_spline2.OR.&
         dft_control%xc_deriv_method_id==xc_deriv_pw)

    DO ispin=1,nspins
       ! introduce a smoothing kernel on the density
       IF (dft_control%xc_rho_smooth_id==xc_rho_no_smooth) THEN
          my_rho_r_local=.FALSE.
          my_rho_r(ispin)%pw => rho_r(ispin)%pw
          IF (needs_rho_g) THEN
             IF (ASSOCIATED(rho_g)) THEN
                my_rho_g_local=.FALSE.
                my_rho_g => rho_g(ispin)%pw
             END IF
          END IF
       ELSE
          my_rho_r_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_r(ispin)%pw,&
               use_data=REALDATA3D, in_space=REALSPACE, &
               error=error)

          SELECT CASE(dft_control%xc_rho_smooth_id)
          CASE (xc_rho_no_smooth)
             CALL pw_copy(rho_r(ispin)%pw,my_rho_r(ispin)%pw)
          CASE (xc_rho_spline2_smooth)
             CALL pw_zero(my_rho_r(ispin)%pw)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=spline2_coeffs, error=error)
          CASE (xc_rho_spline3_smooth)
             CALL pw_zero(my_rho_r(ispin)%pw)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=spline3_coeffs, error=error)
          CASE (xc_rho_nn10)
             CALL pw_zero(my_rho_r(ispin)%pw)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=nn10_coeffs, error=error)
          CASE (xc_rho_nn50)
             CALL pw_zero(my_rho_r(ispin)%pw)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=nn50_coeffs, error=error)
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       END IF
       IF (needs_rho_g.AND..NOT.ASSOCIATED(my_rho_g)) THEN
          my_rho_g_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_g,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
               error=error)
          CALL pw_transfer(my_rho_r(ispin)%pw,my_rho_g)
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
          IF (needs_rho_g) THEN
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
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             IF (my_rho_g_local) THEN
                my_rho_g_local=.FALSE.
                CALL pw_pool_give_back_pw(pw_pool, my_rho_g ,error=error)
             END IF
             NULLIFY(my_rho_g)
             CALL pw_pool_give_back_pw(pw_pool, tmp_g ,error=error)
          ELSE
             SELECT CASE(dft_control%xc_deriv_method_id)
             CASE (xc_deriv_spline2_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=spline2_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_spline3_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=spline3_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_nn10_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=nn10_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_nn50_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=nn50_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_fd)
                DO idir=1,3
                   CALL pw_derive_fd(my_rho_r(ispin)%pw,&
                        drho_r(idir,ispin)%pw,&
                        n=nd(:,idir))
                END DO
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
          END IF

          IF (dft_control%xc_deriv_method_id/=xc_deriv_pw &
               .AND.dft_control%xc_deriv_method_id/=xc_deriv_fd &
               ) THEN
             DO idir=1,3
                drho_r_att(idir)%pw => drho_r(idir,ispin)%pw
             END DO
             CALL pw_spline_scale_deriv(drho_r_att, cell=cell,&
                  error=error)
          END IF
          DO idir=1,3
             CALL zero_core(drho_r(idir,ispin)%pw,atoms,error=error)
          END DO
       END IF
    END DO

    SELECT CASE(dft_control%nspins)
    CASE(1)
       IF (needs%rho_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%rho_1_3, &
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,my_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%rho_1_3(i,j,k)=MAX(my_rho_r(1)%pw%cr3d(i,j,k),0.0_dp)**f13
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
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
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
             rho_set%drho(idir)%array => drho_r(idir,1)%pw%cr3d
             NULLIFY(drho_r(idir,1)%pw%cr3d)
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
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,my_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%rho(i,j,k)=my_rho_r(1)%pw%cr3d(i,j,k)+&
                        my_rho_r(2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       IF (needs%rho_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rho_1_3,error=error)
          !assume that the bounds are the same?
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,my_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%rho_1_3(i,j,k)=MAX(my_rho_r(1)%pw%cr3d(i,j,k)+&
                        my_rho_r(2)%pw%cr3d(i,j,k),0.0_dp)**f13
                END DO
             END DO
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       IF (needs%rho_spin_1_3) THEN
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rhoa_1_3,error=error)
          !assume that the bounds are the same?
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,my_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%rhoa_1_3(i,j,k)=MAX(my_rho_r(1)%pw%cr3d(i,j,k),0.0_dp)**f13
                END DO
             END DO
          END DO
          CALL pw_pool_create_cr3d(pw_pool,rho_set%rhob_1_3,error=error)
          !assume that the bounds are the same?
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,my_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%rhob_1_3(i,j,k)=MAX(my_rho_r(2)%pw%cr3d(i,j,k),0.0_dp)**f13
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
       IF (needs%norm_drho) THEN

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drho, &
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%norm_drho(i,j,k)=SQRT(&
                        (drho_r(1,1)%pw%cr3d(i,j,k)+drho_r(1,2)%pw%cr3d(i,j,k))**2+&
                        (drho_r(2,1)%pw%cr3d(i,j,k)+drho_r(2,2)%pw%cr3d(i,j,k))**2+&
                        (drho_r(3,1)%pw%cr3d(i,j,k)+drho_r(3,2)%pw%cr3d(i,j,k))**2)
                END DO
             END DO
          END DO

          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
       IF (needs%norm_drho_spin) THEN

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhoa, &
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
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
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
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
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
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
             !$omp parallel do default(none) private(i,j,k) shared(rho_set,drho_r,idir)
             DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                   DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
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

    ! tau part
    IF (needs%tau.OR.needs%tau_spin) THEN
       CPPrecondition(ASSOCIATED(tau),cp_failure_level,routineP,error,failure)
       DO ispin=1,nspins
          CPPrecondition(ASSOCIATED(tau(ispin)%pw),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(tau(ispin)%pw%cr3d),cp_failure_level,routineP,error,failure)
       END DO
    END IF
    IF (needs%tau) THEN
       IF (nspins==2) THEN
          CALL pw_pool_create_cr3d(pw_pool,rho_set%tau,&
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,tau)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)

                   rho_set%tau(i,j,k)=&
                        tau(1)%pw%cr3d(i,j,k)+&
                        tau(2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%tau=.TRUE.
       ELSE
          rho_set%tau => tau(1)%pw%cr3d
          rho_set%owns%tau=.FALSE.
       END IF
       rho_set%has%tau=.TRUE.
    END IF
    IF (needs%tau_spin) THEN
       CPPrecondition(nspins==2,cp_failure_level,routineP,error,failure)
       rho_set%tau_a => tau(1)%pw%cr3d
       rho_set%tau_b => tau(2)%pw%cr3d
       rho_set%owns%tau_spin=.FALSE.
       rho_set%has%tau_spin=.TRUE.
    END IF

    CPPostcondition(xc_rho_cflags_equal(rho_set%has,needs),cp_failure_level,routineP,error,failure)

  END SUBROUTINE xc_rho_set_update
  !***************************************************************************

!!****f* xc_rho_set_types/zero_core [1.0] *
!!
!!   NAME
!!     zero_core
!!
!!   FUNCTION
!!     switches to zero the pw deriv towrds the cores of the atoms
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - deriv: the pw the zero out in the core region
!!     - subsys: subsytems with the positions of the atoms
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** *********************************************************************
  SUBROUTINE zero_core(deriv,subsys,error)
    TYPE(pw_type), POINTER                   :: deriv
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: subsys
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, iRadius, &
                                                isubsys, j, k, myi, myj, myk, &
                                                natom
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: indexPos, npoints
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: atBoundary, failure
    REAL(dp)                                 :: int_radius, r, r_n, r_n2, &
                                                radius, radius_2, range_2
    REAL(dp), DIMENSION(2)                   :: ldaCoreRadius
    REAL(dp), DIMENSION(3)                   :: dr, pos, posCorr, rv
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
!JV    logger => cp_error_get_logger(error)
    dr=deriv%pw_grid%dr
    npoints=deriv%pw_grid%npts
    NULLIFY(atom_list,atomic_kinds,particles)
    DO isubsys=1,SIZE(subsys)
       CPPrecondition(ASSOCIATED(subsys(isubsys)%subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys(isubsys)%subsys, atomic_kinds=atomic_kinds, &
            particles=particles,error=error)
       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind,&
               ldaCoreRadius=ldaCoreRadius, atom_list=atom_list,natom=natom)

          IF (ldaCoreRadius(2)>0) THEN
             CPPrecondition(ASSOCIATED(atom_list),cp_failure_level,routineP,error,failure)
             CPPrecondition(SIZE(atom_list)>=natom,cp_failure_level,routineP,error,failure)
             int_radius=ldaCoreRadius(1)
             radius=MAXVAL(dr)*ldaCoreRadius(2)
             range_2=(radius-int_radius)/2.0_dp
             radius_2=int_radius+range_2
             iRadius=CEILING(ldaCoreRadius(2))
             DO iatom=1,natom
                pos=particles%els(atom_list(iatom))%r
                bo=deriv%pw_grid%bounds_local
                atBoundary=.FALSE.
                DO i=1,3
                   posCorr(i)=MODULO(pos(i),REAL(npoints(i),dp)*dr(i))&
                        -REAL(npoints(i),dp)/2.0_dp*dr(i)
                   indexPos(i)=MODULO(NINT(pos(i)/dr(i)),&
                        npoints(i))-npoints(i)/2
                   IF (indexPos(i)-iRadius<-npoints(i)/2) THEN
                      atBoundary=.TRUE.
                      bo(1,i)=indexPos(i)-iRadius
                   ELSE
                      bo(1,i)=MAX(bo(1,i),indexPos(i)-iRadius)
                   END IF
                   IF (indexPos(i)+iRadius>(npoints(i)-1)/2) THEN
                      atBoundary=.TRUE.
                      bo(2,i)=indexPos(i)+iRadius
                   ELSE
                      bo(2,i)=MIN(bo(2,i),indexPos(i)+iRadius)
                   END IF
                END DO

!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('indexPos=(',i6,',',i6,',',i6,')')")&
!FM         indexPos(1),indexPos(2),indexPos(3)
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('dr=(',e10.3,e10.3,e10.3,')')")&
!FM         dr
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('posCorr=(',f10.3,f10.3,f10.3,')')")&
!FM         posCorr
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('bo=(',i6,':',i6,',',i6,':',i6,',',i6,':',i6,')')")&
!FM         bo(1,1),bo(2,1),bo(1,2),bo(2,2),bo(1,3),bo(2,3)
!FM    WRITE(cp_logger_get_default_unit_nr(logger),*) 'radii:',int_radius,radius_2,&
!FM         radius, 'range_2=',range_2

                IF (atBoundary) THEN
                   DO k=bo(1,3),bo(2,3)
                      myk=MODULO(k+npoints(3)/2,npoints(3))-npoints(3)/2
                      IF (myk< deriv%pw_grid%bounds_local(1,3).OR.&
                           myk> deriv%pw_grid%bounds_local(2,3)) CYCLE
                      rv(3)=REAL(k,dp)*dr(3)-posCorr(3)
                      DO j=bo(1,2),bo(2,2)
                         myj=MODULO(j+npoints(2)/2,npoints(2))-npoints(2)/2
                         IF (myj< deriv%pw_grid%bounds_local(1,2).OR.&
                              myj> deriv%pw_grid%bounds_local(2,2)) CYCLE
                         rv(2)=REAL(j,dp)*dr(2)-posCorr(2)
                         DO i=bo(1,1),bo(2,1)
                            myi=MODULO(i+npoints(1)/2,npoints(1))-npoints(1)/2
                            IF (myi< deriv%pw_grid%bounds_local(1,1).OR.&
                                 myj> deriv%pw_grid%bounds_local(2,1)) CYCLE
                            rv(1)=REAL(i,dp)*dr(1)-posCorr(1)
                            r=SQRT(rv(1)**2+rv(2)**2+rv(3)**2)

                            IF (r<radius) THEN
                               IF (r<int_radius) THEN
                                  deriv%cr3d(myi,myj,myk)=0.0_dp
                               ELSEIF (r<radius_2) THEN
!FM                      WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                           "('(',i6,',',i6,',',i6,') within radius2')") myi,myj,myk
                                  r_n=(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  deriv%cr3d(myi,myj,myk)=deriv%cr3d(myi,myj,myk)*&
                                       r_n2*(r_n-0.5_dp*r_n2)
                               ELSE
                                  r_n=2.0_dp-(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  deriv%cr3d(myi,myj,myk)=deriv%cr3d(myi,myj,myk)*&
                                       (1.0_dp-r_n2*(r_n-0.5_dp*r_n2))
                               END IF
                            END IF
                         END DO
                      END DO
                   END DO
                ELSE
                   DO k=bo(1,3),bo(2,3)
                      rv(3)=REAL(k,dp)*dr(3)-posCorr(3)
                      DO j=bo(1,2),bo(2,2)
                         rv(2)=REAL(j,dp)*dr(2)-posCorr(2)
                         DO i=bo(1,1),bo(2,1)
                            rv(1)=REAL(i,dp)*dr(1)-posCorr(1)
                            r=SQRT(rv(1)**2+rv(2)**2+rv(3)**2)
!FM                WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                     "('r=',f8.3)") r

                            IF (r<radius) THEN
                               IF (r<int_radius) THEN
                                  deriv%cr3d(i,j,k)=0.0_dp
                               ELSEIF (r<radius_2) THEN
!FM                      WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                           "('(',i6,',',i6,',',i6,') within radius2')") i,j,k
                                  r_n=(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  deriv%cr3d(i,j,k)=deriv%cr3d(i,j,k)*&
                                       r_n2*(r_n-0.5_dp*r_n2)
                               ELSE
                                  r_n=2.0_dp-(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  deriv%cr3d(i,j,k)=deriv%cr3d(i,j,k)*&
                                       (1.0_dp-r_n2*(r_n-0.5_dp*r_n2))
                               END IF
                            END IF
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END IF
       END DO
    END DO
  END SUBROUTINE zero_core

!!****f* xc_rho_set_types/zero_core_force [1.0] *
!!
!!   NAME
!!     zero_core_force
!!
!!   FUNCTION
!!      adds the force due to the zero_core the the atoms.
!!
!!   NOTES
!!     - think about the sign
!!
!!   ARGUMENTS
!!     - energy: the energy density
!!     - subsys: subsytems with the atoms that will get the force added
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** *********************************************************************
  SUBROUTINE zero_core_force(energy,subsys,error)
    TYPE(pw_type), POINTER                   :: energy
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: subsys
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_core_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, iRadius, &
                                                isubsys, j, k, myi, myj, myk, &
                                                natom, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: indexPos, npoints
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: atBoundary, failure
    REAL(dp)                                 :: int_radius, r, r_n, r_n2, &
                                                radius, radius_2, range_2, tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: f
    REAL(dp), DIMENSION(2)                   :: ldaCoreRadius
    REAL(dp), DIMENSION(3)                   :: dr, pos, posCorr, rv
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_list_type), POINTER        :: particles

!JV    logger => cp_error_get_logger(error)
    dr=energy%pw_grid%dr
    npoints=energy%pw_grid%npts
    failure=.FALSE.

    NULLIFY(para_env,particles, atom_list, atomic_kinds)
    DO isubsys=1,SIZE(subsys)
       CPPrecondition(ASSOCIATED(subsys(isubsys)%subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys(isubsys)%subsys, atomic_kinds=atomic_kinds, &
            particles=particles,para_env=para_env,error=error)
       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind,&
               ldaCoreRadius=ldaCoreRadius, atom_list=atom_list,natom=natom)

          IF (ldaCoreRadius(2)>0) THEN
             CPPrecondition(ASSOCIATED(atom_list),cp_failure_level,routineP,error,failure)
             CPPrecondition(SIZE(atom_list)>=natom,cp_failure_level,routineP,error,failure)
             int_radius=ldaCoreRadius(1)
             radius=MAXVAL(dr)*ldaCoreRadius(2)
             range_2=(radius-int_radius)/2.0_dp
             radius_2=int_radius+range_2
             iRadius=CEILING(ldaCoreRadius(2))
             ALLOCATE(f(3,natom),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             f=0.0_dp
             DO iatom=1,natom
                pos=particles%els(atom_list(iatom))%r
                bo=energy%pw_grid%bounds_local
                atBoundary=.FALSE.
                DO i=1,3
                   posCorr(i)=MODULO(pos(i),REAL(npoints(i),dp)*dr(i))&
                        -REAL(npoints(i),dp)/2.0_dp*dr(i)
                   indexPos(i)=MODULO(NINT(pos(i)/dr(i)),&
                        npoints(i))-npoints(i)/2
                   IF (indexPos(i)-iRadius<-npoints(i)/2) THEN
                      atBoundary=.TRUE.
                      bo(1,i)=indexPos(i)-iRadius
                   ELSE
                      bo(1,i)=MAX(bo(1,i),indexPos(i)-iRadius)
                   END IF
                   IF (indexPos(i)+iRadius>(npoints(i)-1)/2) THEN
                      atBoundary=.TRUE.
                      bo(2,i)=indexPos(i)+iRadius
                   ELSE
                      bo(2,i)=MIN(bo(2,i),indexPos(i)+iRadius)
                   END IF
                END DO

!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('indexPos=(',i6,',',i6,',',i6,')')")&
!FM         indexPos(1),indexPos(2),indexPos(3)
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('dr=(',e10.3,e10.3,e10.3,')')")&
!FM         dr
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('posCorr=(',f10.3,f10.3,f10.3,')')")&
!FM         posCorr
!FM    WRITE(cp_logger_get_default_unit_nr(logger),"('bo=(',i6,':',i6,',',i6,':',i6,',',i6,':',i6,')')")&
!FM         bo(1,1),bo(2,1),bo(1,2),bo(2,2),bo(1,3),bo(2,3)
!FM    WRITE(cp_logger_get_default_unit_nr(logger),*) 'radii:',int_radius,radius_2,&
!FM         radius, 'range_2=',range_2

                IF (atBoundary) THEN
                   DO k=bo(1,3),bo(2,3)
                      myk=MODULO(k+npoints(3)/2,npoints(3))-npoints(3)/2
                      IF (myk< energy%pw_grid%bounds_local(1,3).OR.&
                           myk> energy%pw_grid%bounds_local(2,3)) CYCLE
                      rv(3)=REAL(k,dp)*dr(3)-posCorr(3)
                      DO j=bo(1,2),bo(2,2)
                         myj=MODULO(j+npoints(2)/2,npoints(2))-npoints(2)/2
                         IF (myj< energy%pw_grid%bounds_local(1,2).OR.&
                              myj> energy%pw_grid%bounds_local(2,2)) CYCLE
                         rv(2)=REAL(j,dp)*dr(2)-posCorr(2)
                         DO i=bo(1,1),bo(2,1)
                            myi=MODULO(i+npoints(1)/2,npoints(1))-npoints(1)/2
                            IF (myi< energy%pw_grid%bounds_local(1,1).OR.&
                                 myj> energy%pw_grid%bounds_local(2,1)) CYCLE
                            rv(1)=REAL(i,dp)*dr(1)-posCorr(1)
                            r=SQRT(rv(1)**2+rv(2)**2+rv(3)**2)

                            IF (r<radius) THEN
                               IF (r<int_radius) THEN
                                  energy%cr3d(myi,myj,myk)=0.0_dp
                               ELSEIF (r<radius_2) THEN
!FM                      WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                           "('(',i6,',',i6,',',i6,') within radius2')") myi,myj,myk
                                  r_n=(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  tmp=energy%cr3d(myi,myj,myk)*&
                                       r_n2*(3.0_dp-2.0_dp*r_n)&
                                       /range_2/r
                                  f(1,iatom)=f(1,iatom)+rv(1)*tmp
                                  f(2,iatom)=f(2,iatom)+rv(2)*tmp
                                  f(3,iatom)=f(3,iatom)+rv(3)*tmp
                               ELSE
                                  r_n=2.0_dp-(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  tmp=energy%cr3d(myi,myj,myk)*&
                                       r_n2*(3.0_dp-2.0_dp*r_n)&
                                       /range_2/r
                                  f(1,iatom)=f(1,iatom)+rv(1)*tmp
                                  f(2,iatom)=f(2,iatom)+rv(2)*tmp
                                  f(3,iatom)=f(3,iatom)+rv(3)*tmp
                               END IF
                            END IF
                         END DO
                      END DO
                   END DO
                ELSE
                   DO k=bo(1,3),bo(2,3)
                      rv(3)=REAL(k,dp)*dr(3)-posCorr(3)
                      DO j=bo(1,2),bo(2,2)
                         rv(2)=REAL(j,dp)*dr(2)-posCorr(2)
                         DO i=bo(1,1),bo(2,1)
                            rv(1)=REAL(i,dp)*dr(1)-posCorr(1)
                            r=SQRT(rv(1)**2+rv(2)**2+rv(3)**2)
!FM                WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                     "('r=',f8.3)") r

                            IF (r<radius) THEN
                               IF (r<int_radius) THEN
                                  energy%cr3d(i,j,k)=0.0_dp
                               ELSEIF (r<radius_2) THEN
!FM                      WRITE(cp_logger_get_default_unit_nr(logger),&
!FM                           "('(',i6,',',i6,',',i6,') within radius2')") i,j,k
                                  r_n=(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  tmp=energy%cr3d(i,j,k)*&
                                       r_n2*(3.0_dp-2.0_dp*r_n)&
                                       /range_2/r
                                  f(1,iatom)=f(1,iatom)+rv(1)*tmp
                                  f(2,iatom)=f(2,iatom)+rv(2)*tmp
                                  f(3,iatom)=f(3,iatom)+rv(3)*tmp
                               ELSE
                                  r_n=2.0_dp-(r-int_radius)/range_2
                                  r_n2=r_n*r_n
                                  tmp=energy%cr3d(i,j,k)*&
                                       r_n2*(3.0_dp-2.0_dp*r_n)&
                                       /range_2/r
                                  f(1,iatom)=f(1,iatom)+rv(1)*tmp
                                  f(2,iatom)=f(2,iatom)+rv(2)*tmp
                                  f(3,iatom)=f(3,iatom)+rv(3)*tmp
                               END IF
                            END IF
                         END DO
                      END DO
                   END DO
                END IF
             END DO
             CALL mp_sum(f,para_env%group)
             DO iatom=1,natom
                particles%els(atom_list(iatom))%f(1:3)=particles%els(atom_list(iatom))%f(1:3)&
                     +f(1:3,iatom)
             END DO
             DEALLOCATE(f,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
       END DO
    END DO
  END SUBROUTINE zero_core_force

END MODULE xc_rho_set_types
