!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2002  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****s* cp2k/xc [1.0] *
!!
!!   NAME
!!     xc
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     02.2003 Many many changes
!!
!!   SOURCE
!******************************************************************************

MODULE xc
  USE coefficient_types,               ONLY: PLANEWAVES,&
                                             coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_array_r_utils,                ONLY: cp_3d_r_output
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_output_handling,              ONLY: cp_finished_output,&
                                             cp_unitnr_for_output,&
                                             cp_would_output
  USE dft_types,                       ONLY: dft_control_type,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth,&
                                             xc_rho_nn10
  USE kinds,                           ONLY: dbl,&
                                             wp => dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_compare
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_spline_utils,                 ONLY: pw_spline2_deriv_g,&
                                             pw_spline2_deriv_r,&
                                             pw_spline2_evaluate_values_g,&
                                             pw_spline2_interpolate_values_g,&
                                             pw_spline2_test_coeff,&
                                             pw_spline3_deriv_g,&
                                             pw_spline3_deriv_r,&
                                             pw_spline3_evaluate_values_g,&
                                             pw_spline3_interpolate_values_g,&
                                             pw_spline3_test_coeff,&
                                             pw_spline_scale_deriv,&
                                             pw_nn_smear_g
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_allocate,&
                                             pw_copy,&
                                             pw_deallocate,&
                                             pw_derive,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_functionals,                  ONLY: xc_calculate_lsd,&
                                             xc_create_2nd_deriv_info,&
                                             xc_get_pot_size,&
                                             xc_lda,&
                                             xc_lsd,&
                                             xc_set
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xc_calculate_pw_new, xc_calculate_pw1, xc_calculate_pw2, &
       xc_prep_2nd_deriv, xc_calc_2nd_deriv, xc_prep_2nd_deriv1, &
       xc_calc_2nd_deriv1

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc'

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************

!!****f* xc/xc_calculate_pw1 [1.0] *
!!
!!   NAME
!!     xc_calculate_pw1
!!
!!   SYNOPSIS
!!     Subroutine xc_calculate_pw1(rho_r, rho_g, vxc, exc, dft_control,&
!!         pw_pool, error)
!!       Type(pw_type), Pointer:: rho_r
!!       Type(pw_type), Pointer:: rho_g
!!       Type(pw_type), Pointer:: vxc
!!       Real(Kind=wp), Intent (OUT):: exc
!!       Type(dft_control_type), Intent (IN), Target:: dft_control
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calculate_pw1
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations.
!!     Facility call of xc_calculate_pw_new with non-arrays.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - vxc: will contain the resulting xc potential, has to be 
!!       already allocated
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     2.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_calculate_pw1(rho_r,rho_g,vxc,exc,dft_control,&
       pw_pool,error)
    TYPE(pw_type), POINTER                       :: rho_r
    TYPE(pw_type), POINTER                       :: rho_g
    TYPE(pw_type), POINTER                       :: vxc
    REAL(kind=wp), INTENT(out)                   :: exc
    TYPE(dft_control_type), INTENT(in),TARGET    :: dft_control
    TYPE(pw_pool_type), POINTER                  :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_calculate_pw1',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    TYPE(pw_p_type), DIMENSION(1) :: my_rho_r, my_vxc
    TYPE(pw_p_type), DIMENSION(:), POINTER :: my_rho_g
    TYPE(dft_control_type),POINTER :: my_dft_control
    TYPE(cell_type), POINTER :: cell

    failure=.FALSE.
    NULLIFY(my_rho_g, cell)
    my_rho_r(1)%pw => rho_r
    my_vxc(1)%pw => vxc
    my_dft_control => dft_control

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(vxc),cp_failure_level,routineP,error,failure)

    IF (ANY(dft_control%ftags(1,:))) THEN
       CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
       ALLOCATE(my_rho_g(1), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_rho_g(1)%pw => rho_g
    END IF

    CALL xc_calculate_pw_new(rho_r=my_rho_r, rho_g=my_rho_g,&
         vxc=my_vxc, exc=exc, dft_control=my_dft_control,&
         pw_pool=pw_pool,cell=cell, error=error)

    IF (ASSOCIATED(my_rho_g)) THEN
       DEALLOCATE(my_rho_g,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE xc_calculate_pw1
!***************************************************************************

!!****f* xc/xc_calculate_pw2 [1.0] *
!!
!!   NAME
!!     xc_calculate_pw2
!!
!!   SYNOPSIS
!!     Subroutine xc_calculate_pw2(rho_r, rho_g, vxc, exc, dft_control,&
!!         cell, pw_pool, error)
!!       Type(coeff_type), Dimension(:), Intent (IN):: rho_r
!!       Type(coeff_type), Dimension(:), Pointer:: rho_g
!!       Type(coeff_type), Dimension(:), Intent (OUT):: vxc
!!       Real(Kind=wp), Intent (OUT):: exc
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cell_type), Pointer:: cell
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calculate_pw2
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations.
!!     Facility call of xc_calculate_pw_new with coeff_arrays
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - vxc: will contain the resulting xc potential, has to be 
!!       already allocated
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_calculate_pw2(rho_r,rho_g,vxc,exc,dft_control,&
       cell,pw_pool,error)
    TYPE(coeff_type),DIMENSION(:), INTENT(in)    :: rho_r
    TYPE(coeff_type),DIMENSION(:), POINTER       :: rho_g
    TYPE(coeff_type),DIMENSION(:), INTENT(out)   :: vxc
    REAL(kind=wp), INTENT(out)                   :: exc
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(cell_type), POINTER                     :: cell
    TYPE(pw_pool_type), POINTER                  :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_calculate_pw2',&
         routineP=moduleN//':'//routineN
    TYPE(pw_p_type), DIMENSION(:), POINTER :: my_rho_r, my_rho_g, my_vxc
    INTEGER :: stat, ispin

    failure=.FALSE.
    NULLIFY(my_rho_r, my_rho_g, my_vxc)

    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)    
    ALLOCATE(my_rho_r(SIZE(rho_r)), my_vxc(SIZE(vxc)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO ispin=1,SIZE(rho_r)
          CPPrecondition(rho_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
          my_rho_r(ispin)%pw => rho_r(ispin)%pw
       END DO
       DO ispin=1,SIZE(vxc)
          CPPrecondition(vxc(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
          my_vxc(ispin)%pw => vxc(ispin)%pw
       END DO
    END IF
    IF (ASSOCIATED(rho_g)) THEN
       ALLOCATE(my_rho_g(SIZE(rho_g)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.not.failure) THEN
          DO ispin=1,SIZE(rho_g)
             CPPrecondition(rho_g(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho_g(ispin)%pw => rho_g(ispin)%pw
          END DO
       END IF
    END IF
    
    IF (.not.failure) THEN
       CALL xc_calculate_pw_new(rho_r=my_rho_r, rho_g=my_rho_g,&
            vxc=my_vxc, exc=exc, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, error=error)
    END IF
  END SUBROUTINE xc_calculate_pw2
!***************************************************************************

!!****f* xc/xc_calculate_pw_new [1.0] *
!!
!!   NAME
!!     xc_calculate_pw_new
!!
!!   SYNOPSIS
!!     Subroutine xc_calculate_pw_new(rho_r, rho_g, vxc, exc, dft_control,&
!!         cell, pw_pool, error)
!!       Type(pw_p_type), Dimension(:), Intent (INOUT), Target:: rho_r
!!       Type(pw_p_type), Dimension(:), Pointer:: rho_g
!!       Type(pw_p_type), Dimension(:), Intent (OUT):: vxc
!!       Real(Kind=wp), Intent (OUT):: exc
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(cell_type), Pointer:: cell
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calculate_pw_new
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   NOTES
!!     functional that need ||grad(rho_a+rho_b)|| (tot_grad=.true.)
!!     not implemented (none available in cp2k when this was written).
!!     Also mixed LDA/LSD functional are not supported
!!
!!   ARGUMENTS
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - vxc: will contain the resulting xc potential, has to be 
!!       already allocated
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi; based LDA version of JGH, based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (13-Jun-2002): adaptation to new functionals
!!     Fawzi (11.2002): drho_g(1:3)->drho_g
!!     Fawzi (1.2003). lsd version
!!
!!*** **********************************************************************
  SUBROUTINE xc_calculate_pw_new(rho_r,rho_g,vxc,exc,dft_control,&
       cell,pw_pool,error)
    TYPE(pw_p_type),DIMENSION(:), INTENT(inout), TARGET :: rho_r
    TYPE(pw_p_type),DIMENSION(:), POINTER        :: rho_g
    TYPE(pw_p_type),DIMENSION(:), INTENT(out)    :: vxc
    REAL(kind=wp), INTENT(out)                   :: exc
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(cell_type), POINTER                     :: cell
    TYPE(pw_pool_type), POINTER                  :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL, PARAMETER :: tot_grad=.FALSE.
    LOGICAL :: failure, gradient_f, crossterms, my_rho_r_local,my_rho_g_local
    INTEGER :: nspins, ii, i, j, k, ntot, ngrad_spins, handle, stat, ispin,&
         ipot, idir, npot, ifun
    CHARACTER(len=*), PARAMETER :: routineN='xc_calculate_pw_new',&
         routineP=moduleN//':'//routineN
    INTEGER, DIMENSION(3) :: n
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(pw_type), POINTER :: tmp_g, tmp_r, vxc_g, my_rho_g, my_rho_r
    TYPE(pw_p_type), DIMENSION(3) :: drho_r_att
    TYPE(pw_p_type), DIMENSION(:,:), POINTER :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), ALLOCATABLE :: rgrid, drgrid, pot
    REAL ( kind=wp ) :: eps_rho,alpha
    REAL ( kind=wp ), PARAMETER, dimension(4) :: nn10_coeffs=&
         (/ 125._wp/216._wp, 25._wp/432._wp, 5._wp/864._wp, 1._wp/1728._wp /)
    ! factor 10 between subsequent values

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(tmp_g, tmp_r, vxc_g, my_rho_g, my_rho_r, drho_r)
    logger => cp_error_get_logger(error)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    DO idir=1,3
       NULLIFY(drho_r_att(idir)%pw)
    END DO

    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
!FM    CPPrecondition(associated(vxc),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       gradient_f=ANY(dft_control%ftags(1,:))
       nspins=dft_control%nspins
       ngrad_spins=nspins
       crossterms=ANY(dft_control%ftags(2,:)).and.nspins>1

       IF (gradient_f) THEN
          IF (tot_grad) ngrad_spins=nspins+1
          CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
       END IF
       CPPrecondition(SIZE(vxc)==nspins,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       n(1:3) = rho_r(1)%pw%pw_grid%npts_local (1:3)
       bo = rho_r(1)%pw%pw_grid%bounds_local
       ntot = PRODUCT ( n )

       ALLOCATE(rgrid(ntot,nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (gradient_f) THEN
          ALLOCATE(drgrid(ntot,ngrad_spins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(drho_r(3,nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF


    IF (.NOT.failure) THEN
       DO ispin=1,nspins
          ! introduce a smoothing kernel on the density
          IF (dft_control%xc_rho_smooth_id==xc_rho_no_smooth) THEN
             my_rho_r_local=.FALSE.
             my_rho_g_local=.FALSE.
             my_rho_r => rho_r(ispin)%pw
             IF (gradient_f) THEN
               my_rho_g => rho_g(ispin)%pw
             END IF
          ELSE
             my_rho_r_local=.TRUE.
             my_rho_g_local=.TRUE.
             CALL pw_pool_create_pw(pw_pool, my_rho_r,&
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
             CALL pw_transfer(my_rho_g,my_rho_r)
          END IF

          ! sets rgrid
          !$omp parallel do private(k,j,i,ii)
          DO k = bo(1,3), bo(2,3)
             ii = (k-bo(1,3))*n(1)*n(2)
             DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                   ii = ii + 1
                   rgrid(ii,ispin) = my_rho_r%cr3d(i,j,k)
                END DO
             END DO
          END DO

          IF (gradient_f) THEN ! calculate the grad of rho

             ! optionally applies a scaling function
             IF (dft_control%derive_function_rho) THEN
                ! make rho local
                IF (.NOT.my_rho_r_local) THEN
                   CALL pw_pool_create_pw(pw_pool, my_rho_r,&
                        use_data=REALDATA3D, in_space=REALSPACE, &
                        error=error)
                   CALL pw_copy(rho_r(ispin)%pw,my_rho_r)
                   my_rho_r_local=.TRUE.
                END IF
                IF (.NOT.my_rho_g_local) THEN
                   CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                        use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                        error=error)
!FM                   CALL pw_copy(rho_g(ispin)%pw,my_rho_g)
                   my_rho_g_local=.TRUE.
                END IF

                eps_rho = dft_control%density_cut
                alpha = 1.00_wp / 2.0_wp               ! XXXXXXXXXXXXXXXXXXXXxx alpha hard-coded here XXXXXXXXXXXXXXX
                ! compute a simple function of rho e.g. log
                !$omp parallel do private(k,j,i)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ! my_rho_r%cr3d(i,j,k)=log(MAX(eps_rho,my_rho_r%cr3d(i,j,k)))
                         ! my_rho_r%cr3d(i,j,k)=my_rho_r%cr3d(i,j,k)
                         ! my_rho_r%cr3d(i,j,k)= sqrt(MAX(eps_rho,my_rho_r%cr3d(i,j,k)))
                         my_rho_r%cr3d(i,j,k)=MAX(eps_rho,my_rho_r%cr3d(i,j,k))**alpha
                      END DO
                   END DO
                END DO
                ! get the simple function in g-space
                CALL pw_transfer(my_rho_r,my_rho_g)
             END IF

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

             ! if necessary recovers the derivative of rho
             IF (dft_control%derive_function_rho) THEN
                ! get the derivative of the dens as e.g exp(log(rho))*(d/dx log(rho))
                ! where the last term was computed just before
                DO idir=1,3
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bo(1,3), bo(2,3)
                      ii = (k-bo(1,3))*n(1)*n(2)
                      DO j = bo(1,2), bo(2,2)
                         DO i = bo(1,1), bo(2,1)
                            ii=ii+1
                            ! drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                            ! EXP(MAX(eps_rho,rgrid(ii)))
                            ! drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                            ! 1.0_wp
                            ! drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                            ! (2.0_wp) * my_rho_r%cr3d(i,j,k)
                            drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                                 (1.0_wp/alpha) *  my_rho_r%cr3d(i,j,k) **(1.0_wp/alpha-1.0_wp)
                         END DO
                      END DO
                   END DO
                END DO
             END IF

             CPAssert(.NOT.tot_grad,cp_failure_level,routineP,error,failure)
             ! calc drgrid(:,ispin)=||grad rho_ispin||
             !$omp parallel do private(k,j,i,ii)
             DO k = bo(1,3), bo(2,3)
                ii = (k-bo(1,3))*n(1)*n(2)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii = ii + 1
                      drgrid(ii,ispin) = SQRT(drho_r(1,ispin)%pw%cr3d(i,j,k)**2+&
                           drho_r(2,ispin)%pw%cr3d(i,j,k)**2+&
                           drho_r(3,ispin)%pw%cr3d(i,j,k)**2)
                   END DO
                END DO
             END DO

          END IF

          IF (my_rho_g_local) THEN
             my_rho_g_local=.FALSE.
             CALL pw_pool_give_back_pw(pw_pool, my_rho_g ,error=error)
          END IF
          IF (my_rho_r_local) THEN
             my_rho_r_local=.FALSE.
             CALL pw_pool_give_back_pw(pw_pool, my_rho_r ,error=error)
          END IF
          NULLIFY(my_rho_r, my_rho_g)
       END DO
    END IF

    ! rho is in rgrid, and its deriv (if needed) in drho_r,
    ! drgrid contains the norm of the gradient.
    ! do NOT use rho_r, rho_g, my_rho_r, my_rho_g 

    IF (.NOT.failure) THEN
       ! allocate potential
       npot = xc_get_pot_size ( nspins>1, gradient_f, crossterms, 1 )
       ALLOCATE ( pot ( ntot, 0:(npot-1) ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN

       CALL dcopy(size(pot,1)*size(pot,2),0.0_dbl,0,pot,1)

       ! calculate functionals
       CALL xc_set ( dft_control%density_cut, dft_control%gradient_cut )
       SELECT CASE ( nspins )
       CASE (1)
          IF (gradient_f) THEN
             CALL xc_lda ( functional=dft_control%functional,&
                  ftags=dft_control%ftags, &
                  rho=rgrid(:,1), drho=drgrid(:,1), pot=pot, &
                  order=1, error=error )
          ELSE
             CALL xc_lda ( functional=dft_control%functional,&
                  ftags=dft_control%ftags, &
                  rho=rgrid(:,1), pot=pot, &
                  order=1, error=error )
          END IF
       CASE (2)
          DO ifun = 1, SIZE(dft_control%functional)
             IF ( dft_control%functional(ifun)/="NONE" ) THEN
                IF (dft_control%ftags(1,ifun)) THEN
                   CALL xc_calculate_lsd(functional=dft_control%functional(ifun),&
                        rhoa=rgrid(:,1),&
                        rhob=rgrid(:,2), pot=pot, &
                        order=1, drhoa=drgrid(:,1), drhob=drgrid(:,2),&
                        error=error )
                ELSE
                   CALL xc_calculate_lsd(functional=dft_control%functional(ifun),&
                        rhoa=rgrid(:,1),&
                        rhob=rgrid(:,2), pot=pot, &
                        order=1, error=error )
                END IF
             END IF
          END DO
       CASE default
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='ONLY 1 or 2 spin supported '//&
               CPSourceFileRef,&
               error=error)
       END SELECT

       DEALLOCATE ( rgrid, STAT = stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    IF (.NOT.failure) THEN
       exc = SUM ( pot(:,0) )*rho_r(1)%pw%pw_grid%dvol
       IF ( rho_r(1)%pw%pw_grid%para%mode == 1 ) THEN
          CALL mp_sum ( exc, rho_r(1)%pw%pw_grid%para%group )
       END IF

       ! compose vxc
       ! IF (spezgrid ) THEN should transform to lower cutoff
       DO ispin=1,nspins
          !$omp parallel do private(k,j,i,ii)
          DO k = bo(1,3), bo(2,3)
             ii = (k-bo(1,3))*n(1)*n(2)
             DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                   ii = ii + 1
                   vxc(ispin)%pw%cr3d(i,j,k) = pot(ii,ispin)
                END DO
             END DO
          END DO
       END DO

       IF ( ANY(dft_control%ftags(1,:)) .OR.&
            rho_r(1)%pw%pw_grid%spherical.OR.&
            dft_control%xc_rho_smooth_id/=xc_rho_no_smooth) THEN
          CALL pw_pool_create_pw ( pw_pool, vxc_g, use_data=COMPLEXDATA1D,&
               in_space=RECIPROCALSPACE, error=error)

          DO ispin=1,nspins
             IF (rho_r(1)%pw%pw_grid%spherical.OR.&
                  dft_control%xc_rho_smooth_id/=xc_rho_no_smooth) THEN
                ! filter vxc
                CALL pw_transfer ( vxc(ispin)%pw, vxc_g )
             ELSE
                CALL pw_zero(vxc_g)
             END IF

             IF (ANY(dft_control%ftags(1,:))) THEN ! partial integration
                CALL pw_pool_create_pw ( pw_pool, tmp_g,&
                     use_data=COMPLEXDATA1D,&
                     in_space=RECIPROCALSPACE, error=error)

                ipot=nspins+ispin

                !$omp parallel do private(i)
                DO i = 1, ntot
                   IF ( drgrid(i,ispin) > dft_control%gradient_cut ) THEN
                      pot(i,ipot) = -pot(i,ipot) / drgrid(i,ispin)
                   ELSE
                      pot(i,ipot) = 0._dbl
                   END IF
                END DO

                DO idir = 1, 3
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bo(1,3), bo(2,3)
                      ii = (k-bo(1,3))*n(1)*n(2)
                      DO j = bo(1,2), bo(2,2)
                         DO i = bo(1,1), bo(2,1)
                            ii = ii + 1
                            drho_r(idir,ispin)%pw%cr3d(i,j,k) = pot(ii,ipot)*&
                                 drho_r(idir,ispin)%pw%cr3d(i,j,k)
                         END DO
                      END DO
                   END DO
                END DO

                IF (dft_control%xc_deriv_method_id/=xc_deriv_pw) THEN
                   DO idir=1,3
                      drho_r_att(idir)%pw => drho_r(idir,ispin)%pw
                   END DO
                   CALL pw_spline_scale_deriv(drho_r_att, cell=cell,&
                        transpose=.TRUE.,&
                        error=error)
                END IF
                DO idir = 1, 3          
                   CALL pw_transfer ( drho_r(idir,ispin)%pw, tmp_g )

                   SELECT CASE(dft_control%xc_deriv_method_id)
                   CASE (xc_deriv_pw)
                      CALL pw_derive ( tmp_g, nd(:,idir) )
                   CASE (xc_deriv_spline2)
                      CALL pw_spline2_interpolate_values_g(tmp_g,error=error)
                      CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                   CASE (xc_deriv_spline2_smooth)
                      CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                   CASE (xc_deriv_spline3)
                      CALL pw_spline3_interpolate_values_g(tmp_g,error=error)
                      CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                   CASE (xc_deriv_spline3_smooth)
                      CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                   CASE default
                      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                   END SELECT

                   CALL pw_sumup ( tmp_g, vxc_g )
                END DO
                CALL pw_pool_give_back_pw(pw_pool, tmp_g, error=error)

             END IF

             ! final smoothing if rho was smoothed
             SELECT CASE(dft_control%xc_rho_smooth_id)
             CASE (xc_rho_no_smooth)
             CASE (xc_rho_spline2_smooth)
                CALL pw_spline2_evaluate_values_g(vxc_g,error=error)
             CASE (xc_rho_spline3_smooth)
                CALL pw_spline3_evaluate_values_g(vxc_g,error=error)
             CASE (xc_rho_nn10)
                CALL pw_nn_smear_g(vxc_g,&
                     nn10_coeffs,error=error)
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT

             ! transform back vxc
             IF (rho_r(1)%pw%pw_grid%spherical.OR.&
                  dft_control%xc_rho_smooth_id/=xc_rho_no_smooth) THEN
                CALL pw_transfer ( vxc_g, vxc(ispin)%pw )
             ELSE
                CALL pw_pool_create_pw(pw_pool, tmp_r,&
                     use_data=REALDATA3D, in_space=REALSPACE,&
                     error=error)
                CALL pw_transfer ( vxc_g, tmp_r )
                CALL pw_sumup ( tmp_r, vxc(ispin)%pw )
                CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
             END IF
          END DO

          CALL pw_pool_give_back_pw(pw_pool, vxc_g, error=error)
       END IF
       IF (gradient_f) THEN
          DO ispin=1,nspins
             DO idir=1,3
                CALL pw_pool_give_back_pw(pw_pool,drho_r(idir,ispin)%pw,&
                     error=error)
             END DO
          END DO
          DEALLOCATE(drho_r,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

          DEALLOCATE ( drgrid, STAT = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       DEALLOCATE ( pot, STAT = stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE xc_calculate_pw_new
!***************************************************************************

!!****f* xc/xc_prep_2nd_deriv1 [1.0] *
!!
!!   NAME
!!     xc_prep_2nd_deriv1
!!
!!   SYNOPSIS
!!     Subroutine xc_prep_2nd_deriv1(rho_r, rho_g, dft_control, pw_pool,&
!!         drho_r, pot, spin_pot, grad_pot, ndiag_term, error)
!!       Type(coeff_type), Dimension(:), Pointer:: rho_r, rho_g
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(pw_p_type), Dimension(:,:), Pointer:: drho_r
!!       Real(Kind=wp), Dimension(:,:), Pointer:: pot
!!       Integer, Dimension(:,:), Pointer:: spin_pot
!!       Logical, Dimension(:,:), Pointer:: grad_pot
!!       Logical, Dimension(:), Pointer:: ndiag_term
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_prep_2nd_deriv1
!!
!!   FUNCTION
!!     facility call of xc_prep_2nd_deriv with coeff arrays
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_prep_2nd_deriv1(rho_r,rho_g,dft_control,pw_pool,&
       drho_r,pot,spin_pot,grad_pot,ndiag_term,error)
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho_g
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_p_type), DIMENSION(:,:),POINTER :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), POINTER :: pot
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='xc_prep_2nd_deriv1',&
        routineP=moduleN//':'//routineN
  TYPE(pw_p_type), DIMENSION(:), POINTER :: my_rho_r, my_rho_g
  INTEGER :: stat, ispin

  failure=.FALSE.
  NULLIFY(my_rho_r, my_rho_g)
  
  IF (ASSOCIATED(rho_r)) THEN
     ALLOCATE(my_rho_r(SIZE(rho_r)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        DO ispin=1,SIZE(rho_r)
           CPPrecondition(rho_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
           my_rho_r(ispin)%pw => rho_r(ispin)%pw
        END DO
     END IF
  END IF
  IF (ASSOCIATED(rho_g)) THEN
     ALLOCATE(my_rho_g(SIZE(rho_g)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        DO ispin=1,SIZE(rho_g)
           CPPrecondition(rho_g(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
           my_rho_g(ispin)%pw => rho_g(ispin)%pw
        END DO
     END IF
  END IF
  IF (.NOT. failure) THEN
     CALL xc_prep_2nd_deriv(rho_r=my_rho_r,rho_g=my_rho_g,&
          dft_control=dft_control,pw_pool=pw_pool,&
          drho_r=drho_r,pot=pot,spin_pot=spin_pot,grad_pot=grad_pot,&
          ndiag_term=ndiag_term,error=error)
  END IF
END SUBROUTINE xc_prep_2nd_deriv1
!***************************************************************************

!!****f* qs_kpp1_env_methods/xc_prep_2nd_deriv [1.0] *
!!
!!   NAME
!!     xc_prep_2nd_deriv
!!
!!   SYNOPSIS
!!     Subroutine xc_prep_2nd_deriv(rho_r, rho_g, dft_control, pw_pool,&
!!         drho_r, pot, spin_pot, grad_pot, ndiag_term, error)
!!       Implicit None
!!       Type(pw_p_type), Dimension(:), Pointer:: rho_r, rho_g
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(pw_p_type), Dimension(:,:), Pointer:: drho_r
!!       Real(Kind=wp), Dimension(:,:), Pointer:: pot
!!       Integer, Dimension(:,:), Pointer:: spin_pot
!!       Logical, Dimension(:,:), Pointer:: grad_pot
!!       Logical, Dimension(:), Pointer:: ndiag_term
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_prep_2nd_deriv
!!
!!   FUNCTION
!!     Does the first prepartions for the calculation of the 2nd deriv
!!     done by xc_calc_2nd_deriv
!!
!!   NOTES
!!     I am not sure that making this ugly interface public is a good idea...
!!
!!   INPUTS
!!     - rho_r,rho_g: point where to evaluate the (functional) derivative
!!       of exc in r and g space
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - drho_r(idir,ispin) :: derivative of rho_r in the direction x,y,z
!!       in the r space (allocated if the xc potential is gradient corrected)
!!     - pot: the potential (allocated)
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (allocated, see xc_create_2nd_deriv_info)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!     12.2002 many bugfixes, split in prep and calc [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_prep_2nd_deriv(rho_r, rho_g, dft_control, pw_pool, &
                               drho_r, pot, spin_pot, grad_pot, ndiag_term, &
                               error)
    
    IMPLICIT NONE

    ! arguments
    TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_r, rho_g
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_p_type), DIMENSION(:,:),POINTER :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), POINTER :: pot
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL :: failure, gradient_f, crossterms, do_lsd, do_restricted_tddfpt
    CHARACTER(len=*), PARAMETER :: routineN='xc_prep_2nd_deriv',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle, nspins, spin, ii, i,j,k, n(3), n_tot, idir, ipot,ntot,&
         stat, order,npot
    TYPE(pw_type), POINTER :: tmp_g
    REAL ( kind=wp ) , DIMENSION (:,:), ALLOCATABLE :: rgrid, drgrid
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger
    REAL(KIND=wp) :: rgrid_sum
    CHARACTER(LEN=40), DIMENSION(3)           :: functional
    LOGICAL, DIMENSION (2,3)                  :: ftags
    
    ! ------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(tmp_g)
    NULLIFY(drho_r, pot, spin_pot, grad_pot, ndiag_term)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       if (dft_control%do_tddfpt_calculation) then
          functional = dft_control%tddfpt_control%functional
          ftags      = dft_control%tddfpt_control%ftags
       else
          functional = dft_control%functional
          ftags      = dft_control%ftags
       endif

       nspins     = dft_control%nspins
       do_restricted_tddfpt  = (dft_control%do_tddfpt_calculation .AND. (nspins==1))
       do_lsd     = (nspins==2) .OR. do_restricted_tddfpt
       crossterms = do_lsd .AND. ANY(ftags(2,:))
       gradient_f = ANY(ftags(1,:))

       IF (gradient_f) THEN
          CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
       END IF       
    END IF

    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       n(1:3) = rho_r(1)%pw%pw_grid%npts_local (1:3)
       bo = rho_r(1)%pw%pw_grid%bounds_local  
       ntot = PRODUCT ( n )

       ALLOCATE(rgrid(ntot,nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (gradient_f) THEN
          ALLOCATE(drgrid(ntot,nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(drho_r(3,nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN

! set rgrid
       DO spin=1,nspins
!$omp parallel do private(k,j,i,ii)
          DO k = bo(1,3), bo(2,3)
             ii = (k-bo(1,3))*n(1)*n(2)
             DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                   ii = ii + 1
                   rgrid(ii,spin) = rho_r(spin)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
       END DO
       

       IF (gradient_f) THEN
! calc drho_r=grad rho
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          DO spin=1,nspins
             DO idir=1,3
                NULLIFY(drho_r(idir,spin)%pw)
                CALL pw_pool_create_pw(pw_pool,drho_r(idir,spin)%pw, &
                     use_data=REALDATA3D, in_space=REALSPACE, &
                     error=error)
                CALL pw_copy ( rho_g(spin)%pw, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, drho_r(idir,spin)%pw )
             END DO
          END DO
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,&
               error=error)

! calc drgrid=||grad rho||
          DO spin=1,nspins
!$omp parallel do private(k,j,i,ii)
             DO k = bo(1,3), bo(2,3)
                ii = (k-bo(1,3))*n(1)*n(2)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii = ii + 1
                      drgrid(ii,spin) = SQRT(drho_r(1,spin)%pw%cr3d(i,j,k)**2+&
                           drho_r(2,spin)%pw%cr3d(i,j,k)**2+&
                           drho_r(3,spin)%pw%cr3d(i,j,k)**2)
                   END DO
                END DO
             END DO
          END DO
          
       END IF
    END IF

    IF (.NOT.failure) THEN
       order = -2
! allocate potential
       npot = xc_get_pot_size ( do_lsd, gradient_f, crossterms, order )
       ALLOCATE ( pot ( ntot, npot ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN

! calculate functionals
       CALL xc_set ( dft_control%density_cut, dft_control%gradient_cut )

       IF (nspins == 1) THEN

          !-----------------------------------------------!
          ! do lsd but with rhoa=1/2*rho and rhob=1/2*rho !
          ! and drhoa=1/2*drho and drhob=1/2*rhob         !
          !-----------------------------------------------!
          IF (do_restricted_tddfpt) THEN

             CALL dscal (ntot, 0.5_wp, rgrid(:,1), 1)
             IF (gradient_f) THEN
                CALL dscal (ntot, 0.5_wp, drgrid(:,1), 1)
                CALL xc_lsd ( functional, ftags, &
                              rgrid(:,1), rgrid(:,1), &
                              pot, order, &
                              drho_a=drgrid(:,1), drho_b=drgrid(:,1), &
                              error=error )
             ELSE
                CALL xc_lsd ( functional, ftags, &
                              rgrid(:,1), rgrid(:,1), &
                              pot, order, &
                              error=error )
            END IF
             
          ELSE
             IF (gradient_f) THEN
                CALL xc_lda ( functional=functional,&
                              ftags=ftags, &
                              rho=rgrid(:,1), drho=drgrid(:,1), pot=pot, &
                              order=order, error=error )
             ELSE
                CALL xc_lda ( functional=functional,&
                              ftags=ftags, &
                              rho=rgrid(:,1), pot=pot, &
                              order=order, error=error )
             END IF
          END IF
       ELSE ! nspins == 2
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='lsd not yet implemented '//&
               CPSourceFileRef,&
               error=error)
       END IF

       CALL xc_create_2nd_deriv_info( lsd=do_lsd, gradient=gradient_f,&
            crossterms=crossterms, spin_pot=spin_pot, grad_pot=grad_pot,&
            ndiag_term=ndiag_term, error=error)

    END IF

! dealloc rgrid
    IF (.not.failure) THEN
       DEALLOCATE(rgrid,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

! divide the potential by the gradient where necessary
    IF (.NOT.failure) THEN
       DO ipot=1,npot
          IF (grad_pot(1,ipot).AND.grad_pot(2,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

!$omp parallel do private(ii)
                DO ii=1,ntot
                   IF (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut.or.&
                        drgrid(ii,spin_pot(2,ipot))<dft_control%gradient_cut) THEN
                      pot(ii,ipot)=0.0_wp
                   ELSE
                      pot(ii,ipot)=pot(ii,ipot)/(drgrid(ii,spin_pot(1,ipot))*&
                           drgrid(ii,spin_pot(2,ipot)))
                   END IF
                END DO
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

!$omp parallel do private(ii)
                DO ii=1,ntot
                   IF (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut) THEN
                      pot(ii,ipot)=0.0_wp
                   ELSE
                      pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(1,ipot))**2
                   END IF
                END DO
              END IF
          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)
!$omp parallel do private(ii)
             DO ii=1,ntot
                IF (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut.or.&
                     drgrid(ii,spin_pot(2,ipot))<dft_control%gradient_cut) THEN
                   pot(ii,ipot)=0.0_wp
                ELSE
                   pot(ii,ipot)=pot(ii,ipot)/(drgrid(ii,spin_pot(1,ipot))*&
                        drgrid(ii,spin_pot(2,ipot)))
                END IF
             END DO
          ELSE IF (grad_pot(2,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

!$omp parallel do private(ii)
             DO ii=1,ntot
                IF (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut) THEN
                   pot(ii,ipot)=0.0_wp
                ELSE
                   pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(2,ipot))**2
                END IF
             END DO
          ELSE
! no grad
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
          END IF
       END DO
    END IF

    IF (gradient_f) THEN
       IF (.NOT.failure) THEN
! dealloc drgrid
          DEALLOCATE(drgrid,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE xc_prep_2nd_deriv
!***************************************************************************

!!****f* xc/xc_calc_2nd_deriv1 [1.0] *
!!
!!   NAME
!!     xc_calc_2nd_deriv1
!!
!!   SYNOPSIS
!!     Subroutine xc_calc_2nd_deriv1(rho_r, drho_r, rho1_r, rho1_g, pot,&
!!         spin_pot, grad_pot, ndiag_term, dft_control, v_rspace_new,&
!!         pw_pool, error)
!!       Type(coeff_type), Dimension(:), Pointer:: rho_r, rho1_r, rho1_g
!!       Type(pw_p_type), Dimension(:,:), Pointer:: drho_r
!!       Real(Kind=wp), Dimension(:,:), Pointer:: pot
!!       Integer, Dimension(:,:), Pointer:: spin_pot
!!       Logical, Dimension(:,:), Pointer:: grad_pot
!!       Logical, Dimension(:), Pointer:: ndiag_term
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(coeff_type), Dimension(:), Pointer:: v_rspace_new
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calc_2nd_deriv1
!!
!!   FUNCTION
!!     Facility call of xc_calc_2nd_deriv with coeff arrays
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_r: point where to evaluate the (functional) derivative
!!       of exc in r space
!!     - drho_r(idir,ispin): derivatives of rho_r in in the x,y,z
!!       directions in r-space (should be allocated only if you use a
!!       gradient corrected xc potential)
!!     - rho1_r, rho1_g: direction of the first derivative in r and g space
!!     - pot: the derivative of the potential, already divised by the gradient
!!       if nedeed (grad_pot true).
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (see xc_create_2nd_deriv_info)
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - v_rspace_new: will contain the new potential. The array iterates
!!       on the spins. The grids must be already present.
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_calc_2nd_deriv1(rho_r,drho_r,rho1_r,rho1_g,&
       pot,spin_pot,grad_pot,ndiag_term,dft_control,&
       v_rspace_new,pw_pool,error)
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho1_r, rho1_g
    TYPE(pw_p_type), DIMENSION(:,:),POINTER :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), POINTER :: pot
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace_new
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='xc_calc_2nd_deriv1',&
         routineP=moduleN//':'//routineN
    TYPE(pw_p_type), DIMENSION(:), POINTER :: my_rho_r, my_rho1_r,&
         my_rho1_g
    INTEGER :: stat,ispin

    failure=.FALSE.
    NULLIFY(my_rho_r, drho_r, my_rho1_r, my_rho1_g)
  
  IF (ASSOCIATED(rho_r)) THEN
     ALLOCATE(my_rho_r(SIZE(rho_r)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        DO ispin=1,SIZE(rho_r)
           CPPrecondition(rho_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
           my_rho_r(ispin)%pw => rho_r(ispin)%pw
        END DO
     END IF
  END IF
  IF (ASSOCIATED(rho1_r)) THEN
     ALLOCATE(my_rho1_r(SIZE(rho1_r)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        DO ispin=1,SIZE(rho1_r)
           CPPrecondition(rho1_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
           my_rho1_r(ispin)%pw => rho1_r(ispin)%pw
        END DO
     END IF
  END IF
  IF (ASSOCIATED(rho1_g)) THEN
     ALLOCATE(my_rho1_g(SIZE(rho1_g)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.not.failure) THEN
        DO ispin=1,SIZE(rho1_g)
           CPPrecondition(rho1_g(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
           my_rho1_g(ispin)%pw => rho1_g(ispin)%pw
        END DO
     END IF
  END IF
  IF (.NOT. failure) THEN
     CALL  xc_calc_2nd_deriv(rho_r=my_rho_r,drho_r=drho_r,&
          rho1_r=my_rho1_r,rho1_g=my_rho1_g,&
          pot=pot,spin_pot=spin_pot,grad_pot=grad_pot,&
          ndiag_term=ndiag_term,dft_control=dft_control,&
          v_rspace_new=v_rspace_new,pw_pool=pw_pool,error=error)
  END IF
END SUBROUTINE xc_calc_2nd_deriv1
!***************************************************************************

!!****f* qs_kpp1_env_methods/xc_calc_2nd_deriv [1.0] *
!!
!!   NAME
!!     xc_calc_2nd_deriv
!!
!!   SYNOPSIS
!!     Subroutine xc_calc_2nd_deriv(rho_r, drho_r, rho1_r, rho1_g, pot,&
!!         spin_pot, grad_pot, ndiag_term, dft_control, v_rspace_new,&
!!         pw_pool, error)
!!       Implicit None
!!       Type(pw_p_type), Dimension(:), Pointer:: rho_r, rho1_r, rho1_g
!!       Type(pw_p_type), Dimension(:,:), Pointer:: drho_r
!!       Real(Kind=wp), Dimension(:,:), Pointer:: pot
!!       Integer, Dimension(:,:), Pointer:: spin_pot
!!       Logical, Dimension(:,:), Pointer:: grad_pot
!!       Logical, Dimension(:), Pointer:: ndiag_term
!!       Type(dft_control_type), Pointer:: dft_control
!!       Type(coeff_type), Dimension(:), Pointer:: v_rspace_new
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_calc_2nd_deriv
!!
!!   FUNCTION
!!     Does the first prepartions for the calculation of the 2nd deriv
!!     Calculates the second derivative of E_xc at rho in the direction
!!     rho1  (if you see the second derivative as bilinear form)
!!     partial_rho|_(rho=rho) partial_rho|_(rho=rho) E_xc drho(rho1)drho
!!     The other direction is still indetermined, thus it returns
!!     a potential (partial integration is performed to reduce it to
!!     function of rho, removing the dependence from its partial derivs)
!!
!!   NOTES
!!     drho_r, pot, spin_pot, grad_pot, ndiag_term, should be as returned by
!!     xc_prep_2nd_deriv.
!!     I am not sure that making this ugly interface public is a good idea...
!!
!!   INPUTS
!!     - rho_r: point where to evaluate the (functional) derivative
!!       of exc in r space
!!     - drho_r(idir,ispin): derivatives of rho_r in in the x,y,z
!!       directions in r-space (should be allocated only if you use a
!!       gradient corrected xc potential)
!!     - rho1_r, rho1_g: direction of the first derivative in r and g space
!!     - pot: the derivative of the potential, already divised by the gradient
!!       if nedeed (grad_pot true).
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (see xc_create_2nd_deriv_info)
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - v_rspace_new: will contain the new potential. The array iterates
!!       on the spins. The grids must be already present.
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!     12.2002 many bugfixes, split in prep and calc [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_calc_2nd_deriv(rho_r, drho_r, rho1_r, rho1_g, &
                               pot, spin_pot, grad_pot, ndiag_term, dft_control, &
                               v_rspace_new, pw_pool, error)

    IMPLICIT NONE
    
    ! arguments
    TYPE(pw_p_type), DIMENSION(:), POINTER       :: rho_r, rho1_r, rho1_g
    TYPE(pw_p_type), DIMENSION(:,:),POINTER      :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), POINTER  :: pot
    INTEGER, DIMENSION(:,:), POINTER             :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER             :: grad_pot
    LOGICAL, DIMENSION(:), POINTER               :: ndiag_term
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(coeff_type), DIMENSION(:), POINTER      :: v_rspace_new
    TYPE(pw_pool_type), POINTER                  :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    ! locals
    LOGICAL                                      :: failure, gradient_f, &
                                                    crossterms, do_restricted_tddfpt, do_lsd
    REAL(kind=wp)                                :: diag_factor
    INTEGER                                      :: handle, &
                                                    nspins, spin, &
                                                    ii, i,j,k, &
                                                    n(3), n_tot, ntot, npot, &
                                                    idir, ipot, &
                                                    stat
    TYPE(pw_type), POINTER                       :: tmp_g, tmp_r, v_deriv_g
    TYPE(pw_p_type), DIMENSION(:), ALLOCATABLE   :: v_r_to_deriv, drho1_r
    INTEGER, DIMENSION (3,3)                     :: nd 
    INTEGER, DIMENSION (2,3)                     :: bounds
    TYPE(cp_logger_type), POINTER                :: logger
    CHARACTER(LEN=40), DIMENSION(3)              :: functional
    LOGICAL, DIMENSION (2,3)                     :: ftags

    CHARACTER(len=*), PARAMETER :: routineN='xc_calc_2nd_deriv',&
         routineP=moduleN//':'//routineN

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(tmp_g, tmp_r, v_deriv_g)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_rspace_new),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spin_pot),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(grad_pot),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(ndiag_term),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

       if (dft_control%do_tddfpt_calculation) then
          functional = dft_control%tddfpt_control%functional
          ftags      = dft_control%tddfpt_control%ftags
       else
          functional = dft_control%functional
          ftags      = dft_control%ftags
       endif

       nspins     = dft_control%nspins
       do_restricted_tddfpt  = (dft_control%do_tddfpt_calculation .AND. (nspins==1))
       gradient_f = ANY(ftags(1,:))
       do_lsd     = (nspins==2) .OR. do_restricted_tddfpt
       crossterms = do_lsd .AND. ANY(ftags(2,:))
       IF (dft_control%tddfpt_control%restricted_excitations_type=="SINGLET") THEN
          diag_factor = 1.0_wp
       ELSE
          diag_factor = -1.0_wp
       END IF
       npot = SIZE(pot,2)

       IF (gradient_f) THEN
          CPPrecondition(ASSOCIATED(drho_r),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(rho1_g),cp_failure_level,routineP,error,failure)
       END IF       
       CPPrecondition(SIZE(v_rspace_new)==nspins,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       n(1:3) = rho_r(1)%pw%pw_grid%npts_local (1:3)
       bounds = rho_r(1)%pw%pw_grid%bounds_local  
       ntot = PRODUCT ( n )
    END IF
    
    IF (gradient_f) THEN

! alloc & calc drho1_r=<grad rho,grad rho1>
       IF (.NOT.failure) THEN
          ALLOCATE(drho1_r(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          CALL pw_pool_create_pw(pw_pool,tmp_r, &
               use_data=REALDATA3D,in_space=REALSPACE,&
               error=error)
          DO spin=1,nspins
             NULLIFY(drho1_r(spin)%pw)
             CALL pw_pool_create_pw(pw_pool,drho1_r(spin)%pw, &
                  use_data=REALDATA3D, in_space=REALSPACE, &
                  error=error)
             DO idir=1,3
                CALL pw_copy ( rho1_g(spin)%pw, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, tmp_r )

                IF (idir==1) THEN
!$omp parallel do private(k,j,i)
                   DO k = bounds(1,3), bounds(2,3)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            drho1_r(spin)%pw%cr3d(i,j,k) = &
                                 tmp_r%cr3d(i,j,k)*&
                                 drho_r(idir,spin)%pw%cr3d(i,j,k)
                         END DO
                      END DO
                   END DO
                ELSE
!$omp parallel do private(k,j,i)
                   DO k = bounds(1,3), bounds(2,3)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            drho1_r(spin)%pw%cr3d(i,j,k) = &
                                 drho1_r(spin)%pw%cr3d(i,j,k)+&
                                 tmp_r%cr3d(i,j,k)*&
                                 drho_r(idir,spin)%pw%cr3d(i,j,k)
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END DO
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,&
               error=error)
          CALL pw_pool_give_back_pw(pw_pool,tmp_r,&
               error=error)
       END IF

! alloc v_r_to_deriv
       IF (.NOT.failure) THEN
          ALLOCATE(v_r_to_deriv(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO spin=1,nspins
             NULLIFY(v_r_to_deriv(spin)%pw)
             CALL pw_pool_create_pw(pw_pool,v_r_to_deriv(spin)%pw,&
                  use_data=REALDATA3D, in_space=REALSPACE,error=error)
             CALL pw_zero(v_r_to_deriv(spin)%pw)
          END DO
       END IF
    END IF

! calc v_rspace_new and v_r_to_deriv
    IF (.NOT.failure) THEN

       DO ipot=1, npot ! over all pots

          IF (grad_pot(1,ipot) .AND. grad_pot(2,ipot)) THEN

             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

!$omp parallel do private(k,j,i,ii)
                DO k = bounds(1,3), bounds(2,3)
                   ii = (k-bounds(1,3))*n(1)*n(2)
                   DO j = bounds(1,2), bounds(2,2)
                      DO i = bounds(1,1), bounds(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO

!$omp parallel do private(k,j,i,ii)
                DO k = bounds(1,3), bounds(2,3)
                   ii = (k-bounds(1,3))*n(1)*n(2)
                   DO j = bounds(1,2), bounds(2,2)
                      DO i = bounds(1,1), bounds(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

!$omp parallel do private(k,j,i,ii)
                DO k = bounds(1,3), bounds(2,3)
                   ii = (k-bounds(1,3))*n(1)*n(2)
                   DO j = bounds(1,2), bounds(2,2)
                      DO i = bounds(1,1), bounds(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             END IF

          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

!$omp parallel do private(k,j,i,ii)
             DO k = bounds(1,3), bounds(2,3)
                ii = (k-bounds(1,3))*n(1)*n(2)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      ii=ii+1
                      v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                           v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                           drho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

!$omp parallel do private(k,j,i,ii)
             DO k = bounds(1,3), bounds(2,3)
                ii = (k-bounds(1,3))*n(1)*n(2)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      ii=ii+1
                      v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                           v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                           rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

          ELSE IF (grad_pot(2,ipot)) THEN
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

!$omp parallel do private(k,j,i,ii)
             DO k = bounds(1,3), bounds(2,3)
                ii = (k-bounds(1,3))*n(1)*n(2)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      ii=ii+1
                      v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                           v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                           rho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

!$omp parallel do private(k,j,i,ii)
             DO k = bounds(1,3), bounds(2,3)
                ii = (k-bounds(1,3))*n(1)*n(2)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      ii=ii+1
                      v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                           v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                           drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

          ELSE
! no grad

             !----------------------------------------------!
             ! Time-Dependent DFT (restricted ground state) !
             !                                              ! 
             !                    &^2 Exc                   !
             ! K(r,s,t) = ---------------------- rho1(r')   !
             !            &rho(r,s) &rho_t(r',t)            !
             !                                              !
             ! singlet : K(a,a,r) + K(a,b,r)                !
             ! triplet : K(a,a,r) - K(a,b,r)                !
             !----------------------------------------------!
             IF (do_restricted_tddfpt) THEN

                ! case a,a
                IF (spin_pot(1,ipot)==1 .AND. spin_pot(2,ipot)==1) THEN
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bounds(1,3), bounds(2,3)
                      ii = (k-bounds(1,3))*n(1)*n(2)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            ii=ii+1
                            v_rspace_new(1)%pw%cr3d(i,j,k) = &
                                 v_rspace_new(1)%pw%cr3d(i,j,k) + &
                                 rho1_r(1)%pw%cr3d(i,j,k) * pot(ii,ipot)
                         END DO
                      END DO
                   END DO

                ! case a,b   
                ELSE IF (spin_pot(1,ipot)==1 .AND. spin_pot(2,ipot)==2) THEN
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bounds(1,3), bounds(2,3)
                      ii = (k-bounds(1,3))*n(1)*n(2)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            ii=ii+1
                            v_rspace_new(1)%pw%cr3d(i,j,k) = &
                                 v_rspace_new(1)%pw%cr3d(i,j,k) + &
                                 diag_factor * rho1_r(1)%pw%cr3d(i,j,k) * pot(ii,ipot)
                         END DO
                      END DO
                   END DO
                
                END IF

             ELSE

                IF (ndiag_term(ipot)) THEN
                   CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                        cp_failure_level,cp_assertion_failed,routineP,&
                        "error inconsistent spin_pot info "//&
                        CPSourceFileRef,&
                        error,failure)
                   
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bounds(1,3), bounds(2,3)
                      ii = (k-bounds(1,3))*n(1)*n(2)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            ii=ii+1
                            v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                                 v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                                 rho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                                 pot(ii,ipot)
                         END DO
                      END DO
                   END DO
                   
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bounds(1,3), bounds(2,3)
                      ii = (k-bounds(1,3))*n(1)*n(2)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            ii=ii+1
                            v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                                 v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                                 rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                                 pot(ii,ipot)
                         END DO
                      END DO
                   END DO
                ELSE ! (.NOT.ndiag_term)
                   CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                        cp_failure_level,cp_assertion_failed,routineP,&
                        "error inconsistent spin_pot info "//&
                        CPSourceFileRef,&
                        error,failure)
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bounds(1,3), bounds(2,3)
                      ii = (k-bounds(1,3))*n(1)*n(2)
                      DO j = bounds(1,2), bounds(2,2)
                         DO i = bounds(1,1), bounds(2,1)
                            ii=ii+1
                            v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                                 v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                                 rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                                 pot(ii,ipot)
                         END DO
                      END DO
                   END DO
                END IF ! (ndiag_term)

             END IF ! (do_restricted_tddfpt)
          END IF
       END DO
    END IF

    IF (gradient_f) THEN
       IF (.NOT.failure) THEN
! give back drho1_r
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(pw_pool,drho1_r(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(drho1_r,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

! perform partial integration and filter v_rspace_new
    IF (.NOT.failure) THEN
       IF (gradient_f) THEN
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
       END IF
       CALL pw_pool_create_pw(pw_pool,v_deriv_g, &
            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
            error=error)

       DO spin=1,nspins
          CALL pw_transfer(v_rspace_new(spin)%pw,v_deriv_g)
          
          IF (gradient_f) THEN
             DO idir=1,3
                
!$omp parallel do private(k,j,i)
                DO k = bounds(1,3), bounds(2,3)
                   DO j = bounds(1,2), bounds(2,2)
                      DO i = bounds(1,1), bounds(2,1)
                         v_rspace_new(spin)%pw%cr3d(i,j,k)= &
                              -v_r_to_deriv(spin)%pw%cr3d(i,j,k)* &
                              drho_r(idir,spin)%pw%cr3d(i,j,k)
                      END DO
                   END DO
                END DO
                
                CALL pw_transfer(v_rspace_new(spin)%pw,tmp_g)
                CALL pw_derive(tmp_g,nd(:,idir))
                CALL pw_sumup(tmp_g,v_deriv_g)
             END DO
          END IF

          CALL pw_transfer(v_deriv_g,v_rspace_new(spin)%pw)
       END DO

       IF (gradient_f) THEN
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,error=error)
       END IF
       CALL pw_pool_give_back_pw(pw_pool,v_deriv_g,error=error)

       IF (gradient_f) THEN
! dealloc v_r_to_deriv
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(pw_pool,v_r_to_deriv(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(v_r_to_deriv, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       
    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE xc_calc_2nd_deriv
!***************************************************************************

SUBROUTINE pw_smooth(pw_in,pw_out)
   TYPE(pw_type), POINTER :: pw_in
   TYPE(pw_type), POINTER :: pw_out
   REAL(wp) :: alpha,beta,sum,dr(3),radius,sigma,dist
   INTEGER  :: n(3),nc(3),bo(2,3),i,j,k,ir,il,jl,jr,kr,kl,p,q,r,method
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: Kernel

   n(1:3) = pw_in%pw_grid%npts_local (1:3)
   dr(:) = pw_in%pw_grid%dr(:)
   bo = pw_in%pw_grid%bounds_local

   method = 1 ! hard coded right now, like everything in here

   SELECT CASE(method)
   CASE(1) ! just some averaging over neighbors, very fast
      alpha=1.0_wp
      beta =0.1_wp
      sum = alpha + 6*beta
      alpha = alpha/sum
      beta  = beta/sum
      DO k = bo(1,3), bo(2,3)
         DO j = bo(1,2), bo(2,2)
            DO i = bo(1,1), bo(2,1)
                 ir = MODULO(( i + 1 ) - bo(1,1),n(1))+bo(1,1)
                 il = MODULO(( i - 1 ) - bo(1,1),n(1))+bo(1,1)
                 jr = MODULO(( j + 1 ) - bo(1,2),n(2))+bo(1,2)
                 jl = MODULO(( j - 1 ) - bo(1,2),n(2))+bo(1,2)
                 kr = MODULO(( k + 1 ) - bo(1,3),n(3))+bo(1,3)
                 kl = MODULO(( k - 1 ) - bo(1,3),n(3))+bo(1,3)
                 pw_out%cr3d(i,j,k) =  alpha*pw_in%cr3d(i,j,k)+beta*( &
                                       pw_in%cr3d(il,j,k)+pw_in%cr3d(ir,j,k)+ &
                                       pw_in%cr3d(i,jl,k)+pw_in%cr3d(i,jr,k)+ &
                                       pw_in%cr3d(i,j,kl)+pw_in%cr3d(i,j,kr))
            END DO
         END DO
      END DO
   CASE(2) ! allowing for a more advanced functional form and wider mesh for averaging
           ! gets *very* slow rapidly. A g-space smoother would be possible
           ! however, this will most likely not be positive definite
      radius=0.5_wp
      sigma =0.1_wp
      nc(:)=CEILING(radius/dr(:))
      WRITE(6,*) radius,sigma,nc(:)
      ALLOCATE(Kernel(-nc(1):nc(1),-nc(2):nc(2),-nc(3):nc(3)))
      sum = 0.0_wp
      DO r=-nc(3),nc(3)
       DO q=-nc(2),nc(2)
        DO p=-nc(1),nc(1)
           dist=SQRT((r*dr(3))**2+(q*dr(2))**2+(p*dr(1))**2)
           Kernel(p,q,r)=EXP(-(dist/sigma)**2)
           sum = sum + Kernel(p,q,r)
        ENDDO
       ENDDO
      ENDDO
      ! normalize to 1 exactly.
      DO r=-nc(3),nc(3)
       DO q=-nc(2),nc(2)
        DO p=-nc(1),nc(1)
           Kernel(p,q,r)=Kernel(p,q,r)/sum
        ENDDO
       ENDDO
      ENDDO
      pw_out%cr3d(:,:,:) = 0.0_wp
               DO r=-nc(3),nc(3)
                  DO q=-nc(2),nc(2)
      DO k = bo(1,3), bo(2,3)
                  kr = MODULO(( k + r )- bo(1,3),n(3))+bo(1,3)
         DO j = bo(1,2), bo(2,2)
                     jr = MODULO(( j + q )- bo(1,2),n(2))+bo(1,2)
            DO i = bo(1,1), bo(2,1)
                     DO p=-nc(1),nc(1)
                        ir = MODULO(( i + p )- bo(1,1),n(1))+bo(1,1)
                        pw_out%cr3d(i,j,k) =  pw_out%cr3d(i,j,k)  + &
                                              Kernel(p,q,r)*pw_in%cr3d(ir,jr,kr)
                     ENDDO
            END DO
         END DO
      END DO
                  ENDDO
               ENDDO

      DEALLOCATE(Kernel)
   END SELECT

  END SUBROUTINE

END MODULE xc

!******************************************************************************
