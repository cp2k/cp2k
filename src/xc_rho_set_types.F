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
MODULE xc_rho_set_types
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       xc_deriv_collocate, xc_deriv_nn10_smooth, xc_deriv_nn50_smooth, &
       xc_deriv_pw, xc_deriv_spline2, xc_deriv_spline2_smooth, &
       xc_deriv_spline3, xc_deriv_spline3_smooth, xc_rho_nn10, xc_rho_nn50, &
       xc_rho_no_smooth, xc_rho_spline2_smooth, xc_rho_spline3_smooth
  USE kinds,                           ONLY: dp
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_derive,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_cr3d,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_cr3d,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: &
       nn10_coeffs, nn10_deriv_coeffs, nn50_coeffs, nn50_deriv_coeffs, &
       pw_nn_deriv_r, pw_nn_smear_r, pw_spline2_deriv_g, &
       pw_spline2_interpolate_values_g, pw_spline3_deriv_g, &
       pw_spline3_interpolate_values_g, pw_spline_scale_deriv, &
       spline2_coeffs, spline2_deriv_coeffs, spline3_coeffs, &
       spline3_deriv_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_equal,&
                                             xc_rho_cflags_setall,&
                                             xc_rho_cflags_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_rho_set_types'
  INTEGER, SAVE :: last_rho_set_id=0

  PUBLIC :: xc_rho_set_type
  PUBLIC :: xc_rho_set_create, xc_rho_set_retain, xc_rho_set_release,&
       xc_rho_set_update, xc_rho_set_get

! *****************************************************************************
!> \brief represent a density, with all the representation and data needed
!>      to perform a functional evaluation
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique)
!> \param local_bounds the part of the 3d array on which the functional is
!>        computed
!> \param owns which components are owned by this structure (and should
!>        be deallocated
!> \param has which components are present and up to date
!> \param rho the density
!> \param drho the gradient of the density (x,y and z direction)
!> \param norm_drho the norm of the gradient of the density
!> \param rhoa , rhob: spin alpha and beta parts of the density in the LSD case
!> \param drhoa , drhob: gradient of the spin alpha and beta parts of the density
!>        in the LSD case (x,y and z direction)
!> \param norm_drhoa , norm_drhob: norm of the gradient of rhoa and rhob
!> \param drhoa_drhob the scalar product between the gradient of rhoa and the
!>        one of rhob
!> \param rho_ 1_3: rho^(1.0_dp/3.0_dp)
!> \param rhoa_ 1_3, rhob_1_3: rhoa^(1.0_dp/3.0_dp), rhob^(1.0_dp/3.0_dp)
!> \param tau the kinetic (KohnSham) part of rho
!> \param tau_a the kinetic (KohnSham) part of rhoa
!> \param tau_b the kinetic (KohnSham) part of rhob
!> \note
!>      the use of 3d arrays is the result of trying to use only basic
!>      types (to be generic and independent from the method), and
!>      avoiding copies using the actual structure.
!>      only the part defined by local bounds is guaranteed to be present,
!>      and it is the only meaningful part.
!> \par History
!>      11.2003 created [fawzi & thomas]
!>      12.2008 added laplace parts [mguidon]
!> \author fawzi & thomas
! *****************************************************************************
  TYPE xc_rho_set_type
     INTEGER :: ref_count, id_nr
     INTEGER, DIMENSION(2,3) :: local_bounds
     REAL(kind=dp) :: rho_cutoff, drho_cutoff, tau_cutoff
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
     REAL(kind = dp), DIMENSION(:,:,:), POINTER :: laplace_rho, laplace_rhoa, laplace_rhob
  END TYPE xc_rho_set_type

CONTAINS

! *****************************************************************************
!> \brief allocates and does (minimal) initialization of a rho_set
!> \param rho_set the structure to allocate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_rho_set_create(rho_set,local_bounds,rho_cutoff,drho_cutoff,&
       tau_cutoff,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, DIMENSION(2, 3), INTENT(in)     :: local_bounds
    REAL(kind=dp), INTENT(in), OPTIONAL      :: rho_cutoff, drho_cutoff, &
                                                tau_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       rho_set%tau_cutoff=EPSILON(0.0_dp)
       IF (PRESENT(tau_cutoff)) rho_set%tau_cutoff=tau_cutoff
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
            rho_set%tau,rho_set%tau_a,rho_set%tau_b, rho_set%laplace_rho, rho_set%laplace_rhoa,&
            rho_set%laplace_rhob)
    END IF
  END SUBROUTINE xc_rho_set_create

! *****************************************************************************
!> \brief retains the given rho_set
!> \param rho_set the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_rho_set_retain(rho_set, error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief releases the given rho_set
!> \param rho_set the structure to release
!> \param pw_pool the plae where to give back the arrays
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_rho_set_release(rho_set, pw_pool, error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

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
          IF (rho_set%owns%laplace_rho.AND.ASSOCIATED(rho_set%laplace_rho)) THEN
             DEALLOCATE(rho_set%laplace_rho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF

          IF (rho_set%owns%norm_drho.AND.ASSOCIATED(rho_set%norm_drho)) THEN
             DEALLOCATE(rho_set%norm_drho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%laplace_rho_spin) THEN
             IF (ASSOCIATED(rho_set%laplace_rhoa)) THEN
                DEALLOCATE(rho_set%laplace_rhoa, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%laplace_rhob)) THEN
                DEALLOCATE(rho_set%laplace_rhob, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
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

! *****************************************************************************
!> \brief returns the various attributes of rho_set
!> \param rho_set the object you whant info about
!> \param can_return_null if true the object returned can be null,
!>        if false (the default) it stops with an error if a requested
!>        component is not associated
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      for the other components see the attributes of xc_rho_set
! *****************************************************************************
  SUBROUTINE xc_rho_set_get(rho_set, can_return_null, rho, drho, norm_drho,&
       rhoa, rhob, norm_drhoa, norm_drhob, drhoa_drhob,rho_1_3,rhoa_1_3,&
       rhob_1_3,laplace_rho,laplace_rhoa,laplace_rhob,drhoa,drhob,rho_cutoff,&
       drho_cutoff,tau_cutoff,tau,tau_a,tau_b,local_bounds,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(in), OPTIONAL            :: can_return_null
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: rho
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: drho
    REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, POINTER :: norm_drho, rhoa, &
      rhob, norm_drhoa, norm_drhob, drhoa_drhob, rho_1_3, rhoa_1_3, rhob_1_3, &
      laplace_rho, laplace_rhoa, laplace_rhob
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: drhoa, drhob
    REAL(kind=dp), INTENT(out), OPTIONAL     :: rho_cutoff, drho_cutoff, &
                                                tau_cutoff
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: tau, tau_a, tau_b
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: local_bounds
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       IF (PRESENT(laplace_rho)) THEN
          laplace_rho => rho_set%laplace_rho
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(laplace_rho),cp_failure_level,routineP,error,failure)
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
       IF (PRESENT(laplace_rhoa)) THEN
          laplace_rhoa => rho_set%laplace_rhoa
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(laplace_rhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(laplace_rhob)) THEN
          laplace_rhob => rho_set%laplace_rhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(laplace_rhob),cp_failure_level,routineP,error,failure)
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
       IF (PRESENT(tau_cutoff)) tau_cutoff=rho_set%tau_cutoff
       IF (PRESENT(local_bounds)) local_bounds => rho_set%local_bounds
    END IF
  END SUBROUTINE xc_rho_set_get

! *****************************************************************************
!> \brief cleans (releases) most of the data stored in the rho_set giving back
!>      what it can to the pw_pool
!> \param rho_set the rho_set to be cleaned
!> \param pw_pool place to give back 3d arrays,...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE xc_rho_set_clean(rho_set,pw_pool,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       IF (rho_set%owns%laplace_rho) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%laplace_rho,&
               accept_non_compatible=.TRUE.,error=error)
       ELSE
          NULLIFY(rho_set%laplace_rho)
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
       IF (rho_set%owns%laplace_rho_spin) THEN
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%laplace_rhoa,&
               accept_non_compatible=.TRUE.,error=error)
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%laplace_rhob,&
               accept_non_compatible=.TRUE.,error=error)
       ELSE
          NULLIFY(rho_set%laplace_rhoa, rho_set%laplace_rhob)
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

! *****************************************************************************
!> \brief updates the given rho set with the density given by
!>      rho_r (and rho_g). The rho set will contain the components specified
!>      in needs
!> \param rho_set the rho_set to update
!> \param rho_r the new density (in r space)
!> \param rho_g the new density (in g space, needed for some
!>        derivatives)
!> \param needs the components of rho that are needed
!> \param pw_pool pool for the allocation of pw and cr3d
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs,&
       xc_deriv_method_id,xc_rho_smooth_id,pw_pool,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r, rho_g, tau
    TYPE(xc_rho_cflags_type), INTENT(in)     :: needs
    INTEGER, INTENT(IN)                      :: xc_deriv_method_id, &
                                                xc_rho_smooth_id
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_update', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f13 = (1.0_dp/3.0_dp)

    INTEGER                                  :: i, idir, ispin, j, k, nspins
    INTEGER, DIMENSION(3, 3)                 :: nd, nd_laplace
    LOGICAL                                  :: failure, gradient_f, &
                                                my_rho_g_local, &
                                                my_rho_r_local, needs_rho_g
    REAL(kind=dp)                            :: rho_cutoff
    TYPE(pw_p_type), DIMENSION(2)            :: my_rho_r
    TYPE(pw_p_type), DIMENSION(3)            :: drho_r_att
    TYPE(pw_p_type), DIMENSION(3, 2)         :: drho_r, laplace_rho_r
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
    nd_laplace = RESHAPE((/2,0,0,0,2,0,0,0,2/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL cp_assert(ALL(rho_set%local_bounds==pw_pool%pw_grid%bounds_local),&
         cp_failure_level,cp_assertion_failed,routineP,&
         "pw_pool cr3d have different size than expected",error,failure)
    nspins = SIZE(rho_r)
    rho_set%local_bounds=rho_r(1)%pw%pw_grid%bounds_local
    rho_cutoff=0.5*rho_set%rho_cutoff

    my_rho_g_local=.FALSE.
    ! some checks
    SELECT CASE(nspins)
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
         needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho.OR.needs%laplace_rho &
         .OR.needs%laplace_rho_spin)
    needs_rho_g=(xc_deriv_method_id==xc_deriv_spline3.OR.&
         xc_deriv_method_id==xc_deriv_spline2.OR.&
         xc_deriv_method_id==xc_deriv_pw)
    DO ispin=1,nspins
       ! introduce a smoothing kernel on the density
       IF (xc_rho_smooth_id==xc_rho_no_smooth) THEN
          my_rho_r_local=.FALSE.
          my_rho_r(ispin)%pw => rho_r(ispin)%pw
          IF (needs_rho_g) THEN
             IF (ASSOCIATED(rho_g)) THEN
                my_rho_g_local=.FALSE.
                my_rho_g => rho_g(ispin)%pw
             END IF
          END IF

          my_rho_r_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_r(ispin)%pw,&
               use_data=REALDATA3D, in_space=REALSPACE, &
               error=error)
          CALL pw_copy(rho_r(ispin)%pw,my_rho_r(ispin)%pw,error=error)
       ELSE
          my_rho_r_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_r(ispin)%pw,&
               use_data=REALDATA3D, in_space=REALSPACE, &
               error=error)

          SELECT CASE(xc_rho_smooth_id)
          CASE (xc_rho_no_smooth)
             CALL pw_copy(rho_r(ispin)%pw,my_rho_r(ispin)%pw,error=error)
          CASE (xc_rho_spline2_smooth)
             CALL pw_zero(my_rho_r(ispin)%pw,error=error)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=spline2_coeffs, error=error)
          CASE (xc_rho_spline3_smooth)
             CALL pw_zero(my_rho_r(ispin)%pw,error=error)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=spline3_coeffs, error=error)
          CASE (xc_rho_nn10)
             CALL pw_zero(my_rho_r(ispin)%pw,error=error)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=nn10_coeffs, error=error)
          CASE (xc_rho_nn50)
             CALL pw_zero(my_rho_r(ispin)%pw,error=error)
             CALL pw_nn_smear_r(pw_in=rho_r(ispin)%pw,&
                  pw_out=my_rho_r(ispin)%pw,&
                  coeffs=nn50_coeffs, error=error)
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
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
             IF (.NOT.ASSOCIATED(my_rho_g)) THEN
                my_rho_g_local=.TRUE.
                CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                     use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                     error=error)
                CALL pw_transfer(my_rho_r(ispin)%pw,my_rho_g,error=error)
             END IF
             CALL pw_pool_create_pw(pw_pool, tmp_g,&
                  use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                  error=error)
             SELECT CASE(xc_deriv_method_id)
             CASE (xc_deriv_pw)
                DO idir=1,3
                   CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                   CALL pw_derive ( tmp_g, nd(:,idir) ,error=error)
                   CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw ,error=error)
                END DO
                IF(needs%laplace_rho.OR.needs%laplace_rho_spin) THEN
                  DO idir=1,3
                    NULLIFY(laplace_rho_r(idir,ispin)%pw)
                    CALL pw_pool_create_pw(pw_pool,laplace_rho_r(idir,ispin)%pw, &
                                           use_data=REALDATA3D, in_space=REALSPACE, &
                                           error=error)
                    CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                    CALL pw_derive ( tmp_g, nd_laplace(:,idir) ,error=error)
                    CALL pw_transfer ( tmp_g, laplace_rho_r(idir,ispin)%pw ,error=error)
                  END DO
                END IF
             CASE (xc_deriv_spline2)
                IF (.NOT.my_rho_g_local) THEN
                   CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                        use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                        error=error)
                   my_rho_g_local=.TRUE.
                   CALL pw_copy(rho_g(ispin)%pw, my_rho_g,error=error)
                END IF
                CALL pw_spline2_interpolate_values_g(my_rho_g,error=error)
                DO idir=1,3
                   CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                   CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                   CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw ,error=error)
                END DO
             CASE (xc_deriv_spline3)
                IF (.NOT.my_rho_g_local) THEN
                   CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                        use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                        error=error)
                   CALL pw_copy(rho_g(ispin)%pw, my_rho_g,error=error)
                   my_rho_g_local=.TRUE.
                END IF
                CALL pw_spline3_interpolate_values_g(my_rho_g,error=error)
                DO idir=1,3
                   CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                   CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                   CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw ,error=error)
                END DO
             CASE (xc_deriv_collocate)
                DO idir=1,3
                   CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                   CALL pw_derive ( tmp_g, nd(:,idir) ,error=error)
                   CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw ,error=error)
                END DO
               CALL cp_unimplemented_error(fromWhere=routineP, &
                    message="Drho collocation not implemented", &
                    error=error, error_level=cp_failure_level)
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             CALL pw_pool_give_back_pw(pw_pool, tmp_g ,error=error)
             IF (my_rho_g_local) THEN
                my_rho_g_local=.FALSE.
                CALL pw_pool_give_back_pw(pw_pool, my_rho_g ,error=error)
             END IF
          ELSE
             SELECT CASE(xc_deriv_method_id)
             CASE (xc_deriv_spline2_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw,error=error)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=spline2_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_spline3_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw,error=error)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=spline3_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_nn10_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw,error=error)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=nn10_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_nn50_smooth)
                DO idir=1,3
                   CALL pw_zero(drho_r(idir,ispin)%pw,error=error)
                   CALL pw_nn_deriv_r(pw_in=my_rho_r(ispin)%pw,&
                        pw_out=drho_r(idir,ispin)%pw,&
                        coeffs=nn50_deriv_coeffs, idir=idir, error=error)
                END DO
             CASE (xc_deriv_collocate)
                DO idir=1,3
                   CALL pw_copy ( my_rho_g, tmp_g ,error=error)
                   CALL pw_derive ( tmp_g, nd(:,idir) ,error=error)
                   CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw ,error=error)
                END DO
               CALL cp_unimplemented_error(fromWhere=routineP, &
                    message="Drho collocation not implemented", &
                    error=error, error_level=cp_failure_level)
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
          END IF

          IF (xc_deriv_method_id/=xc_deriv_pw) THEN
             DO idir=1,3
                drho_r_att(idir)%pw => drho_r(idir,ispin)%pw
             END DO
             CALL pw_spline_scale_deriv(drho_r_att, error=error)
          END IF

       END IF

    END DO

    SELECT CASE(nspins)
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
       IF (needs%laplace_rho) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%laplace_rho, &
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,laplace_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%laplace_rho(i,j,k)=&
                        laplace_rho_r(1,1)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(2,1)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(3,1)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%laplace_rho=.TRUE.
          rho_set%has%laplace_rho=.TRUE.
       END IF

       IF (needs%drho) THEN
          DO idir=1,3
             rho_set%drho(idir)%array => drho_r(idir,1)%pw%cr3d
             NULLIFY(drho_r(idir,1)%pw%cr3d)
          END DO
          rho_set%owns%drho=.TRUE.
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
       IF (needs%laplace_rho_spin) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%laplace_rhoa, &
               error=error)
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,laplace_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%laplace_rhoa(i,j,k)=&
                        laplace_rho_r(1,1)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(2,1)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(3,1)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO

          CALL pw_pool_create_cr3d(pw_pool, rho_set%laplace_rhob, &
               error=error)
          rho_set%owns%laplace_rho_spin=.TRUE.
          !$omp parallel do default(none) private(i,j,k) shared(rho_set,laplace_rho_r)
          DO k=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                   rho_set%laplace_rhob(i,j,k)=&
                        laplace_rho_r(1,2)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(2,2)%pw%cr3d(i,j,k)+&
                        laplace_rho_r(3,2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO

          rho_set%owns%laplace_rho_spin=.TRUE.
          rho_set%has%laplace_rho_spin=.TRUE.
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
          IF(needs%laplace_rho.OR.needs%laplace_rho_spin) THEN
            CALL pw_pool_give_back_pw(pw_pool, laplace_rho_r(idir,ispin)%pw, &
               accept_non_compatible=.TRUE., error=error)
          END IF
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

    CPPostcondition(xc_rho_cflags_equal(rho_set%has,needs,error=error),cp_failure_level,routineP,error,failure)

  END SUBROUTINE xc_rho_set_update

END MODULE xc_rho_set_types
