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
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     (13-Feb-2001) JGH, based on earlier version of apsi
!!     02.2003 Many many changes [fawzi]
!!     03.2004 new xc interface [fawzi]
!!     04.2004 kinetic functionals [fawzi]
!!
!!   SOURCE
!******************************************************************************

MODULE xc
  USE atomic_kind_types,               ONLY: get_atomic_kind_set
  USE coefficient_types,               ONLY: PLANEWAVES,&
                                             coeff_type
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type,&
                                             xc_debug_new_routine,&
                                             xc_deriv_fd,&
                                             xc_deriv_nn10_smooth,&
                                             xc_deriv_nn50_smooth,&
                                             xc_deriv_pw,&
                                             xc_deriv_spline2,&
                                             xc_deriv_spline2_smooth,&
                                             xc_deriv_spline3,&
                                             xc_deriv_spline3_smooth,&
                                             xc_new_f_routine,&
                                             xc_old_f_routine,&
                                             xc_rho_nn10,&
                                             xc_rho_nn50,&
                                             xc_rho_no_smooth,&
                                             xc_rho_spline2_smooth,&
                                             xc_rho_spline3_smooth,&
                                             xc_test_f_routine,&
                                             xc_test_lsd_f_routine
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE kinds,                           ONLY: dbl,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED,&
                                             pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_cr3d,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: nn10_coeffs,&
                                             nn10_deriv_coeffs,&
                                             nn50_coeffs,&
                                             nn50_deriv_coeffs,&
                                             pw_nn_deriv_r,&
                                             pw_nn_smear_g,&
                                             pw_nn_smear_r,&
                                             pw_spline2_deriv_g,&
                                             pw_spline2_evaluate_values_g,&
                                             pw_spline2_interpolate_values_g,&
                                             pw_spline3_deriv_g,&
                                             pw_spline3_evaluate_values_g,&
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
                                             pw_create,&
                                             pw_derive,&
                                             pw_derive_fd,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             MAX_LABEL_LENGTH
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_calc_derivatives,&
                                             xc_functional_eval,&
                                             xc_functional_get_info
  USE xc_functionals,                  ONLY: xc_create_2nd_deriv_info
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type,&
                                             xc_rho_set_create,&
                                             xc_rho_set_get,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update,&
                                             zero_core,&
                                             zero_core_force

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xc_vxc_pw_create1, xc_vxc_pw_create, xc_calculate_pw1, &
       xc_prep_2nd_deriv, xc_calc_2nd_deriv, xc_prep_2nd_deriv1, &
       xc_calc_2nd_deriv1, xc_rho_set_and_dset_create, xc_exc_calc,&
       xc_calc_2nd_deriv_new, xc_prep_2nd_deriv_new

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc'

!!***
!******************************************************************************

CONTAINS

  !******************************************************************************

!!****f* xc/xc_calculate_pw1 [1.0] *
!!
!!   NAME
!!     xc_calculate_pw1
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations.
!!     Facility call of xc_calculate_pw_old with non-arrays.
!!
!!   NOTES
!!      to remove
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
!!*** **********************************************************************
  SUBROUTINE xc_calculate_pw1(rho_r,rho_g,vxc,exc,dft_control,&
       pw_pool,error)
    TYPE(pw_type), POINTER                   :: rho_r, rho_g, vxc
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), INTENT(in), &
      TARGET                                 :: dft_control
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_pw1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: my_dft_control
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_g, my_rho_r, my_vxc

    failure=.FALSE.
    NULLIFY(my_rho_g, cell, my_rho_r, my_vxc)
    ALLOCATE(my_rho_r(1), my_vxc(1), stat=stat)
    my_rho_r(1)%pw => rho_r
    my_vxc(1)%pw => vxc
    my_dft_control => dft_control

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(vxc),cp_failure_level,routineP,error,failure)

    IF (ANY(dft_control%xc_control%gradient_functionals)) THEN
       CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
       ALLOCATE(my_rho_g(1), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_rho_g(1)%pw => rho_g
    END IF

    CALL xc_calculate_pw_old(rho_r=my_rho_r, rho_g=my_rho_g,&
         vxc=my_vxc, exc=exc, dft_control=my_dft_control,&
         pw_pool=pw_pool,cell=cell, error=error)

    IF (ASSOCIATED(my_rho_g)) THEN
       DEALLOCATE(my_rho_g,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE xc_calculate_pw1
  !***************************************************************************

!!****f* xc/xc_vxc_pw_create1 [1.0] *
!!
!!   NAME
!!     xc_vxc_create1
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations.
!!     depending on the selected functional_routine in dft_control calls
!!     the correct routine
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - vxc_rho: will contain the v_xc part that depend on rho
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - vxc_tau: will contain the kinetic (tau) part of v_xcthe functional
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - tau: value of the kinetic density tau on the grid (can be null,
!!       used only with meta functionals)
!!     - vxc: will contain the resulting xc potential, has to be 
!!       already allocated
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - atoms: a subsys with the atoms whose xc you are calculating
!!     - calculate_force: if the forces arsing from a special treatement
!!       of the core of the atoms should be added to the atoms
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE xc_vxc_pw_create1(vxc_rho,vxc_tau,rho_r,rho_g,tau,exc,dft_control,&
       cell,pw_pool,atoms,calculate_forces,error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout)                          :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, OPTIONAL                        :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_create1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)    
    CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN

       SELECT CASE(dft_control%xc_control%functional_routine)
       CASE(xc_old_f_routine)
          ALLOCATE(vxc_rho(dft_control%nspins),stat=stat)
          DO ispin=1,dft_control%nspins
             NULLIFY(vxc_rho(ispin)%pw)
             CALL pw_pool_create_pw(pw_pool,vxc_rho(ispin)%pw,in_space=REALSPACE,&
                  use_data=REALDATA3D, error=error)
             CALL pw_zero(vxc_rho(ispin)%pw)
          END DO
          CALL xc_calculate_pw_old(rho_r=rho_r, rho_g=rho_g,&
               vxc=vxc_rho, exc=exc, dft_control=dft_control,&
               cell=cell, pw_pool=pw_pool, error=error)
       CASE(xc_new_f_routine)
          CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau,tau=tau,&
               rho_r=rho_r, rho_g=rho_g, exc=exc, dft_control=dft_control,&
               cell=cell, pw_pool=pw_pool, atoms=atoms,&
               calculate_forces=calculate_forces,error=error)
       CASE(xc_debug_new_routine)
          CALL xc_vxc_pw_create_debug(vxc_rho=vxc_rho, vxc_tau=vxc_tau,tau=tau,&
               rho_r=rho_r, rho_g=rho_g, exc=exc, dft_control=dft_control,&
               cell=cell, pw_pool=pw_pool, atoms=atoms,&
               calculate_forces=calculate_forces,error=error)
       CASE(xc_test_f_routine)
          ALLOCATE(vxc_rho(dft_control%nspins),stat=stat)
          DO ispin=1,dft_control%nspins
             NULLIFY(vxc_rho(ispin)%pw)
             CALL pw_pool_create_pw(pw_pool,vxc_rho(ispin)%pw,in_space=REALSPACE,&
                  use_data=REALDATA3D, error=error)
             CALL pw_zero(vxc_rho(ispin)%pw)
          END DO
          CALL xc_calculate_pw_test(rho_r=rho_r, rho_g=rho_g,&
               vxc=vxc_rho, exc=exc, dft_control=dft_control,&
               cell=cell, pw_pool=pw_pool, atoms=atoms,&
               calculate_forces=calculate_forces, error=error)
       CASE(xc_test_lsd_f_routine)
          CALL xc_vxc_pw_create_test_lsd(vxc_rho=vxc_rho, vxc_tau=vxc_tau,&
               tau=tau, rho_r=rho_r, rho_g=rho_g, exc=exc, &
               dft_control=dft_control, cell=cell, pw_pool=pw_pool,&
               atoms=atoms, calculate_forces=calculate_forces,error=error)
       CASE default
       END SELECT
    END IF

  END SUBROUTINE xc_vxc_pw_create1
  !***************************************************************************

!!****f* xc/xc_calculate_pw_test [1.0] *
!!
!!   NAME
!!     calculate_pw_test
!!
!!   FUNCTION
!!     calculates both the new and the old version of the functionals and
!!     compares them
!!
!!   NOTES
!!     for debugging only: leaks, and non parallel
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
!!*** **********************************************************************
  SUBROUTINE xc_calculate_pw_test(rho_r,rho_g,vxc,exc,dft_control,&
       cell,pw_pool, atoms,calculate_forces,error)
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout)                          :: vxc
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, OPTIONAL                        :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_pw_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, ispin, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: diff, exc2, maxdiff
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: tau, vxc2, vxc_tau2
    TYPE(xc_derivative_set_type), POINTER    :: dSet1, dSet2
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_set_type), POINTER           :: rho_set2

  failure=.FALSE.
  NULLIFY(vxc2,vxc_tau2,tau,dSet1,dSet2,rho_set2,pot,deriv)
  
  IF (.NOT. failure) THEN
       bo = rho_r(1)%pw%pw_grid%bounds_local

       PRINT *, "about to calculate xc (old)"
       CALL xc_calculate_pw_old_prep(rho_r=rho_r, rho_g=rho_g,&
            vxc=vxc, exc=exc, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, dSet=dSet1,&
            error=error)
       CALL xc_calculate_pw_old(rho_r=rho_r, rho_g=rho_g,&
            vxc=vxc, exc=exc, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, error=error)
       PRINT *, "did calculate xc (old)"
       PRINT *, "about to calculate xc (new)"
       CALL xc_rho_set_and_dset_create(rho_r=rho_r, rho_g=rho_g,&
            tau=tau,dft_control=dft_control,xc_control=dft_control%xc_control,&
            cell=cell, pw_pool=pw_pool, rho_set=rho_set2,&
            deriv_set=dSet2,deriv_order=1,atoms=atoms,&
            needs_basic_components=.FALSE.,error=error)
       CALL xc_vxc_pw_create1(rho_r=rho_r, rho_g=rho_g,tau=tau,&
            vxc_rho=vxc2, vxc_tau=vxc_tau2,exc=exc2, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, atoms=atoms,&
            calculate_forces=calculate_forces, error=error)
       CPPostcondition(.NOT.ASSOCIATED(vxc_tau2),cp_failure_level,routineP,error,failure)
       PRINT *, "did calculate xc (new)"
       
       OPEN(unit=120,file="rho.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_r(1)%pw%cr3d
       WRITE(unit=120) pot(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="ndrho.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set2%norm_drho
       WRITE(unit=120) pot(:,:,bo(1,3))
       CLOSE(unit=120)
       OPEN(unit=120,file="vxc.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => vxc(1)%pw%cr3d
       WRITE(unit=120) pot(:,:,bo(1,3))
       CLOSE(unit=120)
       OPEN(unit=120,file="vxc2.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => vxc2(1)%pw%cr3d
       WRITE(unit=120) pot(:,:,bo(1,3))
       CLOSE(unit=120)

       DO ispin=0,2
          maxDiff=0.0_dp
          IF (ispin==0) THEN
             deriv => xc_dset_get_derivative(dSet2, "", error=error)
       OPEN(unit=120,file="e_0.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       OPEN(unit=121,file="e_0-2.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
          ELSEIF (ispin==1) THEN
             deriv => xc_dset_get_derivative(dSet2, "(rho)", error=error)
       OPEN(unit=120,file="e_rho.bindata",status="unknown",access='sequential',&
            form="unformatted")
       OPEN(unit=121,file="e_rho-2.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
          ELSE
             deriv => xc_dset_get_derivative(dSet2, "(norm_drho)", error=error)
       OPEN(unit=120,file="e_ndrho.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       OPEN(unit=121,file="e_ndrho-2.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
          END IF
          IF (ASSOCIATED(deriv)) THEN
             CALL xc_derivative_get(deriv, deriv_data=pot, error=error)
             
             WRITE(unit=120) pot(:,:,bo(1,3))
             WRITE(unit=121) dSet1%data(1:(bo(2,2)-bo(1,2)+1)*&
                  (bo(2,1)-bo(1,1)+1),ispin)
             ii=0
             DO k=bo(1,3),bo(2,3)
                DO j=bo(1,2),bo(2,2)
                   DO i=bo(1,1),bo(2,1)
                      ii=ii+1
                      diff=ABS(pot(i,j,k)-dSet1%data(ii,ispin))
                      IF (ii==1) THEN
                         PRINT *,"pot",ispin,"=",pot(i,j,k),"vs",dSet1%data(ii,ispin),"diff=",diff
                      END IF
                      IF(diff>maxDiff) maxDiff=diff
!FM                      IF (diff>5.0e-11) THEN
!FM                         PRINT *, "diff=",diff," at ",i,",",j,",",k," spin=",ispin,"p=",pot(i,j,k),"vs",dSet1%data(ii,ispin)
!FM                         PRINT *,"r=",rho_r(1)%pw%cr3d(i,j,k),"r2=",rho_set2%rho(i,j,k),"dr2=",rho_set2%norm_drho(i,j,k)
!FM                      END IF
                   END DO
                END DO
             END DO
             PRINT *,"diff for pot",ispin,"=",maxDiff
          END IF
          CLOSE(unit=120)
          CLOSE(unit=121)
       END DO

       PRINT *,"calc diff on vxc"
       maxDiff=0.0_dp
       DO ispin=1,SIZE(vxc)
          DO k=bo(1,3),bo(2,3)
             DO j=bo(1,2),bo(2,2)
                DO i=bo(1,1),bo(2,1)
                   diff=ABS(vxc(ispin)%pw%cr3d(i,j,k)-&
                        vxc2(ispin)%pw%cr3d(i,j,k))
                      IF (ii==0) THEN
                         PRINT *,"vxc",ispin,"=",vxc(ispin)%pw%cr3d(i,j,k),"vs",vxc2(ispin)%pw%cr3d(i,j,k),"diff=",diff
                      END IF
                   IF (maxDiff<diff) maxDiff=diff
!FM                   IF (diff>5.0e-11) THEN
!FM                      PRINT *, "diff=",diff," at ",i,",",j,",",k," spin=",ispin,"rho=",rho_r(ispin)%pw%cr3d(i,j,k)
!FM                   END IF
                END DO
             END DO
          END DO
       END DO
       PRINT *,"diff exc=",ABS(exc-exc2),"diff vxc=",maxdiff
!       CPPostcondition(maxdiff<5.e-11,cp_failure_level,routineP,error,failure)
!       CPPostcondition(ABS(exc-exc2)<1.e-14,cp_failure_level,routineP,error,failure)

       DO ispin=1,SIZE(vxc2)
          CALL pw_pool_give_back_pw(pw_pool,vxc2(ispin)%pw,&
               error=error)
       END DO
       DEALLOCATE(vxc2, stat=stat)
       CALL xc_rho_set_release(rho_set2,error=error)
       CALL xc_dset_release(dSet2,error=error)
       CALL xc_dset_release(dSet1, error=error)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     
  END IF
END SUBROUTINE xc_calculate_pw_test
!***************************************************************************

!!****f* xc/xc_vxc_pw_create_test_lsd [1.0] *
!!
!!   NAME
!!     xc_vxc_pw_create_test_lsd
!!
!!   FUNCTION
!!     calculates vxc using lsd with rhoa=rhob=0.5*rho and compares
!!     with the lda result
!!
!!   NOTES
!!     for debugging only: leaks, and non parallel
!!
!!   INPUTS
!!     - vxc_rho: will contain the v_xc part that depend on rho
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - vxc_tau: will contain the kinetic (tau) part of v_xc
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - tau: value of the kinetic density tau on the grid (can be null,
!!       used only with meta functionals)
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE xc_vxc_pw_create_test_lsd(vxc_rho,vxc_tau,rho_r,rho_g,tau,&
     exc,dft_control, cell,pw_pool, atoms, calculate_forces, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, INTENT(in), OPTIONAL            :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create_test_lsd', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=80)                        :: filename
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), POINTER                  :: split_desc
    INTEGER                                  :: i, ii, ispin, j, k, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: diff, exc2, maxdiff, tmp
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot, pot2, pot3
    TYPE(cp_sll_xc_deriv_type), POINTER      :: deriv_iter
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho2_g, rho2_r, tau2, &
                                                vxc_rho2, vxc_tau2
    TYPE(xc_derivative_set_type), POINTER    :: dSet1, dSet2
    TYPE(xc_derivative_type), POINTER        :: deriv, deriv2, deriv3
    TYPE(xc_rho_set_type), POINTER           :: rho_set1, rho_set2

  failure=.FALSE.
  NULLIFY(vxc_rho2,vxc_tau2,tau2,dSet1,dSet2,rho_set1,rho_set2,split_desc,pot,pot3,pot3,&
       deriv,deriv2,deriv3,rho2_g)

  IF (.NOT. failure) THEN
     bo = rho_r(1)%pw%pw_grid%bounds_local

     ALLOCATE(rho2_r(2), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ispin=1,2
        NULLIFY(rho2_r(ispin)%pw)
        CALL pw_pool_create_pw(pw_pool,rho2_r(ispin)%pw,in_space=REALSPACE,&
             use_data=REALDATA3D, error=error)
     END DO
     DO k=bo(1,3),bo(2,3)
        DO j=bo(1,2),bo(2,2)
           DO i=bo(1,1),bo(2,1)
              tmp=rho_r(1)%pw%cr3d(i,j,k)*0.5
              rho2_r(1)%pw%cr3d(i,j,k)=tmp
              rho2_r(2)%pw%cr3d(i,j,k)=tmp
           END DO
        END DO
     END DO

     IF (ASSOCIATED(tau)) THEN
        ALLOCATE(tau2(2),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,2
           NULLIFY(tau2(ispin)%pw)
           CALL pw_pool_create_pw(pw_pool,tau2(ispin)%pw,in_space=REALSPACE,&
                use_data=REALDATA3D, error=error)
        END DO

        DO k=bo(1,3),bo(2,3)
           DO j=bo(1,2),bo(2,2)
              DO i=bo(1,1),bo(2,1)
                 tmp=tau(1)%pw%cr3d(i,j,k)*0.5
                 tau2(1)%pw%cr3d(i,j,k)=tmp
                 tau2(2)%pw%cr3d(i,j,k)=tmp
              END DO
           END DO
        END DO
     END IF
       
       PRINT *, "about to calculate xc (lda)"
       CALL xc_rho_set_and_dset_create(rho_r=rho_r, rho_g=rho_g,&
            tau=tau,dft_control=dft_control,xc_control=dft_control%xc_control,&
            cell=cell, pw_pool=pw_pool,rho_set=rho_set1,&
            deriv_set=dSet1, deriv_order=1,atoms=atoms,&
            needs_basic_components=.FALSE.,error=error)
       CALL xc_vxc_pw_create(rho_r=rho_r, rho_g=rho_g,tau=tau,&
            vxc_rho=vxc_rho,vxc_tau=vxc_tau, exc=exc, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, atoms=atoms,&
            calculate_forces=calculate_forces, error=error)
       PRINT *, "did calculate xc (lda)"
       dft_control%nspins=2
       PRINT *, "about to calculate xc (lsd)"
       CALL xc_rho_set_and_dset_create(rho_set=rho_set2,deriv_set=dSet2,&
            rho_r=rho2_r, rho_g=rho2_g,tau=tau2, dft_control=dft_control,&
            xc_control=dft_control%xc_control,&
            cell=cell, pw_pool=pw_pool, deriv_order=1,atoms=atoms,&
            needs_basic_components=.FALSE.,error=error)
       CALL xc_vxc_pw_create(rho_r=rho2_r, rho_g=rho2_g,tau=tau2,&
            vxc_rho=vxc_rho2,vxc_tau=vxc_tau2,exc=exc2, dft_control=dft_control,&
            cell=cell, pw_pool=pw_pool, atoms=atoms, error=error)
       PRINT *, "did calculate xc (new)"
       dft_control%nspins=1
       PRINT *, "at (0,0,0) rho_r=",rho_r(1)%pw%cr3d(0,0,0),&
            "rho2_r(1)=",rho2_r(1)%pw%cr3d(0,0,0),&
            "rho2_r(2)=",rho2_r(2)%pw%cr3d(0,0,0),&
            "rho_r_sm=",rho_set1%rho(0,0,0), "rhoa2_r_sm=",rho_set2%rhoa(0,0,0),&
            "rhob2_r_sm=",rho_set2%rhob(0,0,0)
       OPEN(unit=120,file="rho.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set1%rho
       WRITE(unit=120) pot(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="rhoa.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set2%rhoa
       WRITE(unit=120) pot(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="rhob.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set2%rhob
       WRITE(unit=120) pot(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="ndrho.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set1%norm_drho
       WRITE(unit=120) pot(:,:,bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="ndrhoa.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set2%norm_drhoa
       WRITE(unit=120) pot(:,:,bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="ndrhob.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => rho_set2%norm_drhob
       WRITE(unit=120) pot(:,:,bo(2,3)) 
       CLOSE(unit=120)
       IF (rho_set1%has%tau) THEN
          OPEN(unit=120,file="tau.bindata",status="unknown",access='sequential',&
               form="unformatted",action="write")
          pot => rho_set1%tau
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
       END IF
       IF (rho_set2%has%tau_spin) THEN
          OPEN(unit=120,file="tau_a.bindata",status="unknown",access='sequential',&
               form="unformatted",action="write")
          pot => rho_set2%tau_a
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
          OPEN(unit=120,file="tau_v.bindata",status="unknown",access='sequential',&
               form="unformatted",action="write")
          pot => rho_set2%tau_b
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
       END IF
       OPEN(unit=120,file="vxc.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => vxc_rho(1)%pw%cr3d
       WRITE(unit=120) pot(:,:,bo(2,3)) 
       CLOSE(unit=120)
       OPEN(unit=120,file="vxc2.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
       pot => vxc_rho2(1)%pw%cr3d
       WRITE(unit=120) pot(:,:,bo(2,3)) 
       CLOSE(unit=120)
       IF (ASSOCIATED(vxc_tau)) THEN
          OPEN(unit=120,file="vxc_tau.bindata",status="unknown",access='sequential',&
            form="unformatted",action="write")
          pot => vxc_tau(1)%pw%cr3d
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
       END IF
       IF (ASSOCIATED(vxc_tau2)) THEN
          OPEN(unit=120,file="vxc_tau2_a.bindata",status="unknown",access='sequential',&
               form="unformatted",action="write")
          pot => vxc_tau2(1)%pw%cr3d
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
          OPEN(unit=120,file="vxc_tau2_b.bindata",status="unknown",access='sequential',&
               form="unformatted",action="write")
          pot => vxc_tau2(2)%pw%cr3d
          WRITE(unit=120) pot(:,:,bo(2,3)) 
          CLOSE(unit=120)
       END IF

       PRINT *,"calc diff on vxc"
       maxDiff=0.0_dp
       DO ispin=1,1
          ii=0
          DO k=bo(1,3),bo(2,3)
             DO j=bo(1,2),bo(2,2)
                DO i=bo(1,1),bo(2,1)
                   ii=ii+1
                   diff=ABS(vxc_rho(ispin)%pw%cr3d(i,j,k)-&
                        vxc_rho2(ispin)%pw%cr3d(i,j,k))
                      IF (ii==1) THEN
                         PRINT *,"vxc",ispin,"=",vxc_rho(ispin)%pw%cr3d(i,j,k),"vs",vxc_rho2(ispin)%pw%cr3d(i,j,k),"diff=",diff
                      END IF
                   IF (maxDiff<diff)THEN
                      maxDiff=diff
                      PRINT *, "diff=",diff," at ",i,",",j,",",k,&
                           " spin=",ispin,"rho=",rho_set1%rho(i,j,k),&
                           " ndrho=",rho_set1%norm_drho(i,j,k)
                   END IF
                END DO
             END DO
          END DO
       END DO
       PRINT *,"diff exc=",ABS(exc-exc2),"diff vxc=",maxdiff
!       CPPostcondition(maxdiff<5.e-11,cp_failure_level,routineP,error,failure)
!       CPPostcondition(ABS(exc-exc2)<1.e-14,cp_failure_level,routineP,error,failure)

       IF (ASSOCIATED(vxc_tau)) THEN
       PRINT *,"calc diff on vxc_tau"
       maxDiff=0.0_dp
       DO ispin=1,1
          ii=0
          DO k=bo(1,3),bo(2,3)
             DO j=bo(1,2),bo(2,2)
                DO i=bo(1,1),bo(2,1)
                   ii=ii+1
                   diff=ABS(vxc_tau(ispin)%pw%cr3d(i,j,k)-&
                        vxc_tau2(ispin)%pw%cr3d(i,j,k))
                      IF (ii==1) THEN
                         PRINT *,"vxc_tau",ispin,"=",vxc_tau(ispin)%pw%cr3d(i,j,k),"vs",vxc_tau2(ispin)%pw%cr3d(i,j,k),"diff=",diff
                      END IF
                   IF (maxDiff<diff)THEN
                      maxDiff=diff
                      PRINT *, "diff=",diff," at ",i,",",j,",",k,&
                           " spin=",ispin,"rho=",rho_set1%rho(i,j,k),&
                           " ndrho=",rho_set1%norm_drho(i,j,k)
                   END IF
                END DO
             END DO
          END DO
       END DO
       PRINT *,"diff exc=",ABS(exc-exc2),"diff vxc_tau=",maxdiff
    END IF
!       CPPostcondition(maxdiff<5.e-11,cp_failure_level,routineP,error,failure)
!       CPPostcondition(ABS(exc-exc2)<1.e-14,cp_failure_level,routineP,error,failure)
!------
       deriv_iter => dSet1%derivs
       DO WHILE (cp_sll_xc_deriv_next(deriv_iter,el_att=deriv,error=error))
          CALL xc_derivative_get(deriv,&
               split_desc=split_desc,deriv_data=pot,&
               error=error)
          SELECT CASE (SIZE(split_desc))
          CASE(0)
             filename="e_0.bindata"
             deriv2 => xc_dset_get_derivative(dSet2, "", error=error)
          CASE(1)
             filename="e_"//TRIM(split_desc(1))//".bindata"
             IF (split_desc(1)=="rho") THEN
                deriv2 => xc_dset_get_derivative(dSet2, "(rhoa)", error=error)
             ELSEIF (split_desc(1)=="tau") THEN
                deriv2 => xc_dset_get_derivative(dSet2,"(tau_a)",error=error)
             ELSEIF (split_desc(1)=="norm_drho") THEN
                deriv2 => xc_dset_get_derivative(dSet2, "(norm_drhoa)", error=error)
                deriv3 => xc_dset_get_derivative(dSet2, "(norm_drho)", error=error)
                IF (ASSOCIATED(deriv3)) THEN
                   IF (ASSOCIATED(deriv2)) THEN
                      CALL xc_derivative_get(deriv2,&
                           deriv_data=pot2,&
                           error=error)
                      CALL xc_derivative_get(deriv3,&
                           deriv_data=pot3,&
                           error=error)
                      pot2=pot2+pot3
                   ELSE
                      deriv2 => deriv3
                   END IF
                   NULLIFY(deriv3,pot2,pot3)
                END IF
             ELSE
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          CASE default
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          CALL xc_derivative_get(deriv2,&
               deriv_data=pot2,&
               error=error)
          PRINT *, "checking ",filename
          maxDiff=0.0_dp
          DO k=bo(1,3),bo(2,3)
             DO j=bo(1,2),bo(2,2)
                DO i=bo(1,1),bo(2,1)
                   diff=ABS(pot(i,j,k)-pot2(i,j,k))
                   IF (maxDiff<diff) THEN
                      maxDiff=diff
                      PRINT *, "ediff(",i,j,k,")=",maxDiff,&
                           "rho=",rho_set1%rho(i,j,k),&
                           "ndrho=",rho_set1%norm_drho(i,j,k)
                   END IF
                END DO
             END DO
          END DO
          PRINT *,"maxdiff ",filename,"=",maxDiff
          OPEN (unit=120,file=TRIM(filename),status="unknown",&
               access='sequential',&
               form="unformatted")
          WRITE (unit=120) pot(:,:,bo(2,3))
          CLOSE (unit=120)
       END DO
       deriv_iter => dSet2%derivs
       DO WHILE (cp_sll_xc_deriv_next(deriv_iter,el_att=deriv,error=error))
          CALL xc_derivative_get(deriv,&
               split_desc=split_desc,deriv_data=pot,&
               error=error)
          SELECT CASE (SIZE(split_desc))
          CASE(0)
             filename="e_0-2.bindata"
          CASE(1)
             filename="e_"//TRIM(split_desc(1))//"-2.bindata"
          CASE default
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          OPEN (unit=120,file=TRIM(filename),status="unknown",&
               access='sequential',&
               form="unformatted")
          WRITE (unit=120) pot(:,:,bo(2,3))
          CLOSE (unit=120)
       END DO
!-------
       CALL xc_rho_set_release(rho_set1,error=error)
       CALL xc_rho_set_release(rho_set2,error=error)
       CALL xc_dset_release(dSet2,error=error)
       CALL xc_dset_release(dSet1, error=error)
       DO ispin=1,2
          CALL pw_pool_give_back_pw(pw_pool,rho2_r(ispin)%pw,&
               error=error)
          CALL pw_pool_give_back_pw(pw_pool,vxc_rho2(ispin)%pw,&
               error=error)
          IF (ASSOCIATED(vxc_tau2)) THEN
             CALL pw_pool_give_back_pw(pw_pool,vxc_tau2(ispin)%pw,&
                  error=error)
          END IF
       END DO
       DEALLOCATE(vxc_rho2,rho2_r,rho2_g, stat=stat)
       IF (ASSOCIATED(vxc_tau2)) THEN
          DEALLOCATE(vxc_tau2,stat=stat)
       END IF
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     
  END IF
END SUBROUTINE xc_vxc_pw_create_test_lsd
!***************************************************************************

!!****f* xc/xc_vxc_pw_create_debug [1.0] *
!!
!!   NAME
!!     xc_vxc_pw_create_debug
!!
!!   FUNCTION
!!     calculates vxc outputting the yz plane of rho, and of the various components
!!     of the the derivatives and of vxc
!!
!!   NOTES
!!     for debugging only.
!!
!!   INPUTS
!!     - vxc_rho: will contain the v_xc part that depend on rho
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - vxc_tau: will contain the kinetic (tau) part of v_xc
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - tau: value of the kinetic density tau on the grid (can be null,
!!       used only with meta functionals)
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
SUBROUTINE xc_vxc_pw_create_debug(vxc_rho,vxc_tau,rho_r,rho_g,tau,exc,&
     dft_control, cell,pw_pool,atoms,calculate_forces, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, INTENT(in), OPTIONAL            :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create_debug', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=80)                        :: filename
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), POINTER                  :: split_desc
    INTEGER                                  :: i, ispin, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_xc_deriv_type), POINTER      :: deriv_iter
    TYPE(xc_derivative_set_type), POINTER    :: dSet1
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_set_type), POINTER           :: rho_set1

  failure=.FALSE.
  NULLIFY(dSet1,rho_set1,split_desc,pot,&
       deriv)
  logger => cp_error_get_logger(error)

  IF (.NOT. failure) THEN
     bo = rho_r(1)%pw%pw_grid%bounds_local

     CALL xc_rho_set_and_dset_create(rho_r=rho_r, rho_g=rho_g,&
          tau=tau,dft_control=dft_control,xc_control=dft_control%xc_control,&
          cell=cell, pw_pool=pw_pool,rho_set=rho_set1,&
          deriv_set=dSet1, deriv_order=1,atoms=atoms,&
          needs_basic_components=.FALSE.,error=error)

     ! outputs 0,:,: plane
     IF (bo(1,1)<=0.AND.0<=bo(2,1)) THEN
        IF (rho_set1%has%rho_spin) THEN
           OPEN(unit=120,file="rhoa.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%rhoa
           WRITE(unit=120) pot(0,:,:)
           CLOSE(unit=120)
           OPEN(unit=120,file="rhob.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%rhob
           WRITE(unit=120) pot(0,:,:) 
           CLOSE(unit=120)
        END IF
        IF (rho_set1%has%norm_drho) THEN
           OPEN(unit=120,file="ndrho.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%norm_drho
           WRITE(unit=120) pot(0,:,:) 
           CLOSE(unit=120)
        END IF
        IF (rho_set1%has%norm_drho_spin) THEN
           OPEN(unit=120,file="ndrhoa.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%norm_drhoa
           WRITE(unit=120) pot(0,:,:) 
           CLOSE(unit=120)
           OPEN(unit=120,file="ndrhob.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%norm_drhob
           WRITE(unit=120) pot(0,:,:) 
           CLOSE(unit=120)
        END IF
        IF (rho_set1%has%rho) THEN
           OPEN(unit=120,file="rho.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%rho
           WRITE(unit=120) pot(0,:,:)
           CLOSE(unit=120)
        END IF
        IF (rho_set1%has%tau) THEN
           OPEN(unit=120,file="tau.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%tau
           WRITE(unit=120) pot(0,:,:)
           CLOSE(unit=120)
        END IF
        IF (rho_set1%has%tau_spin) THEN
           OPEN(unit=120,file="tau_a.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%tau_a
           WRITE(unit=120) pot(0,:,:)
           CLOSE(unit=120)
           OPEN(unit=120,file="tau_b.bindata",status="unknown",access='sequential',&
                form="unformatted",action="write")
           pot => rho_set1%tau_b
           WRITE(unit=120) pot(0,:,:)
           CLOSE(unit=120)
        END IF

        deriv_iter => dSet1%derivs
        DO WHILE (cp_sll_xc_deriv_next(deriv_iter,el_att=deriv,error=error))
           CALL xc_derivative_get(deriv,&
                split_desc=split_desc,deriv_data=pot,&
                error=error)
           SELECT CASE (SIZE(split_desc))
           CASE(0)
              filename="e_0.bindata"
           CASE(1)
              filename="e_"//TRIM(split_desc(1))//".bindata"
           CASE default
              CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
           END SELECT
           OPEN (unit=120,file=TRIM(filename),status="unknown",&
                access='sequential',&
                form="unformatted")
           WRITE (unit=120) pot(0,:,:)
           CLOSE (unit=120)
        END DO
     END IF

     CALL xc_vxc_pw_create(vxc_rho=vxc_rho,vxc_tau=vxc_tau,&
          rho_r=rho_r, rho_g=rho_g,tau=tau,&
          exc=exc, dft_control=dft_control,&
          cell=cell, pw_pool=pw_pool,atoms=atoms,&
          calculate_forces=calculate_forces,error=error)

     ! outputs 0,:,: plane
     IF (bo(1,1)<=0.AND.0<=bo(2,1)) THEN
        IF (ASSOCIATED(vxc_rho)) THEN
           DO ispin=1,SIZE(vxc_rho)
              WRITE (filename,"('vxc-',i1,'.bindata')") ispin
              OPEN(unit=120,file=filename,status="unknown",access='sequential',&
                   form="unformatted",action="write")
              pot => vxc_rho(ispin)%pw%cr3d
              WRITE(unit=120) pot(0,:,:)
              CLOSE(unit=120)

              pot => vxc_rho(ispin)%pw%cr3d
              DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                    DO i=bo(1,1),bo(2,1)
                       IF (ABS(pot(i,j,k))>10.0_dbl) THEN
                          WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' vxc_',i1,'(',i6,',',i6,',',i6,')=',e11.4)",&
                                  advance="no") ispin,i,j,k,pot(i,j,k)
                          IF (rho_set1%has%rho_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' rho=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%rhoa(i,j,k), rho_set1%rhob(i,j,k)
                          ELSE IF (rho_set1%has%rho) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' rho=',e11.4)",advance="no") rho_set1%rho(i,j,k)
                          END IF
                          IF (rho_set1%has%norm_drho_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' ndrho=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%norm_drhoa(i,j,k), &
                                  rho_set1%norm_drhob(i,j,k)
                          ELSE IF (rho_set1%has%norm_drho) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' ndrho=',e11.4)",advance="no") rho_set1%norm_drho(i,j,k)
                          END IF
                          IF (rho_set1%has%tau_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' tau=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%tau_a(i,j,k), &
                                  rho_set1%tau_b(i,j,k)
                          ELSE IF (rho_set1%has%tau) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' tau=',e11.4)",advance="no") rho_set1%tau(i,j,k)
                          END IF

                          deriv_iter => dSet1%derivs
                          DO WHILE (cp_sll_xc_deriv_next(deriv_iter,el_att=deriv,error=error))
                             CALL xc_derivative_get(deriv,&
                                  split_desc=split_desc,deriv_data=pot,&
                                  error=error)
                             SELECT CASE (SIZE(split_desc))
                             CASE(0)
                                filename=" e_0"
                             CASE(1)
                                filename=" e_"//TRIM(split_desc(1))
                             CASE default
                                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                             END SELECT
                             WRITE (unit=cp_logger_get_default_unit_nr(logger),&
                                  fmt="(a,'=',e11.4)",advance="no") &
                                  TRIM(filename),pot(i,j,k)
                          END DO

                          WRITE(cp_logger_get_default_unit_nr(logger),&
                               "()")
                       END IF
                    END DO
                 END DO
              END DO
           END DO
        END IF
        IF (ASSOCIATED(vxc_tau)) THEN
           DO ispin=1,SIZE(vxc_tau)
              WRITE (filename,"('vxc_tau_',i1,'.bindata')") ispin
              OPEN(unit=120,file=filename,status="unknown",access='sequential',&
                   form="unformatted",action="write")
              pot => vxc_tau(ispin)%pw%cr3d
              WRITE(unit=120) pot(0,:,:)
              CLOSE(unit=120)

              pot => vxc_tau(ispin)%pw%cr3d
              DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                    DO i=bo(1,1),bo(2,1)
                       IF (ABS(pot(i,j,k))>10.0_dbl) THEN
                          WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "('vxc_tau_',i1,'(',i6,',',i6,',',i6,')=',e11.4)",&
                                  advance="no") ispin,i,j,k,pot(i,j,k)
                          IF (rho_set1%has%rho_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' rho=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%rhoa(i,j,k), rho_set1%rhob(i,j,k)
                          ELSE IF (rho_set1%has%rho) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' rho=',e11.4)",advance="no") rho_set1%rho(i,j,k)
                          END IF
                          IF (rho_set1%has%norm_drho_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' ndrho=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%norm_drhoa(i,j,k), &
                                  rho_set1%norm_drhob(i,j,k)
                          ELSE IF (rho_set1%has%norm_drho) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' ndrho=',e11.4)",advance="no") rho_set1%norm_drho(i,j,k)
                          END IF
                          IF (rho_set1%has%tau_spin) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' tau=(',e11.4,',',e11.4,')')",advance="no")&
                                  rho_set1%tau_a(i,j,k), &
                                  rho_set1%tau_b(i,j,k)
                          ELSE IF (rho_set1%has%tau) THEN
                             WRITE(cp_logger_get_default_unit_nr(logger),&
                                  "(' tau=',e11.4)",advance="no") rho_set1%tau(i,j,k)
                          END IF

                          deriv_iter => dSet1%derivs
                          DO WHILE (cp_sll_xc_deriv_next(deriv_iter,el_att=deriv,error=error))
                             CALL xc_derivative_get(deriv,&
                                  split_desc=split_desc,deriv_data=pot,&
                                  error=error)
                             SELECT CASE (SIZE(split_desc))
                             CASE(0)
                                filename=" e_0"
                             CASE(1)
                                filename=" e_"//TRIM(split_desc(1))
                             CASE default
                                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                             END SELECT
                             WRITE (unit=cp_logger_get_default_unit_nr(logger),&
                                  fmt="(a,'=',e11.4)",advance="no") &
                                  TRIM(filename),pot(i,j,k)
                          END DO

                          WRITE(cp_logger_get_default_unit_nr(logger),"()")
                       END IF
                    END DO
                 END DO
              END DO
           END DO
        END IF
     END IF

     CALL xc_dset_release(dSet1, error=error)
     CALL xc_rho_set_release(rho_set1,error=error)
  END IF
END SUBROUTINE xc_vxc_pw_create_debug
!***************************************************************************

! debugging only, leaks
  SUBROUTINE xc_calculate_pw_old_prep(rho_r,rho_g,vxc,exc,dft_control,&
       cell,pw_pool,dSet,error)
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout)                          :: vxc
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(xc_derivative_set_type), POINTER    :: dSet
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_pw_old_prep', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: tot_grad = .FALSE.

    INTEGER                                  :: handle, i, idir, ii, ispin, &
                                                j, k, ngrad_spins, nspins, &
                                                ntot, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: crossterms, failure, &
                                                gradient_f, my_rho_g_local, &
                                                my_rho_r_local
    REAL(KIND=dp)                            :: alpha, eps_rho
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drgrid, rgrid
    TYPE(pw_p_type), DIMENSION(3)            :: drho_r_att
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    TYPE(pw_type), POINTER                   :: my_rho_g, my_rho_r, tmp_g
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(tmp_g,  my_rho_g, my_rho_r, drho_r, derivative_set)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    DO idir=1,3
       NULLIFY(drho_r_att(idir)%pw)
    END DO

    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    !FM    CPPrecondition(associated(vxc),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       gradient_f=ANY(dft_control%xc_control%gradient_functionals)
       nspins=dft_control%nspins
       ngrad_spins=nspins
       crossterms=ANY(dft_control%xc_control%crossterms).AND.nspins>1

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

       CALL xc_dset_create(derivative_set, pw_pool=pw_pool, error=error)

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
             CASE (xc_rho_nn50)
                CALL pw_nn_smear_g(my_rho_g,&
                     nn50_coeffs, error=error)
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

                eps_rho = dft_control%xc_control%density_cut
                alpha = 1.00_dp / 2.0_dp               ! XXXXXXXXXXXXXXXXXXXXxx alpha hard-coded here XXXXXXXXXXXXXXX
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
                            ! 1.0_dp
                            ! drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                            ! (2.0_dp) * my_rho_r%cr3d(i,j,k)
                            drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                                 (1.0_dp/alpha) *  my_rho_r%cr3d(i,j,k) **(1.0_dp/alpha-1.0_dp)
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

       SELECT CASE (nspins)
       CASE (1)
          IF (gradient_f) THEN
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                  error=error)
          ELSE
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rgrid(:,1), &
                  error=error)
          END IF
       CASE (2)
          IF (gradient_f) THEN
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                  rhob=rgrid(:,2), drhob=drgrid(:,2), &
                  error=error)
          ELSE
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), &
                  rhob=rgrid(:,2), &
                  error=error)
          END IF
       CASE default
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='ONLY 1 or 2 spin supported '//&
               CPSourceFileRef,&
               error=error)
       END SELECT

       DEALLOCATE ( rgrid, STAT = stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    dSet => derivative_set
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xc_calculate_pw_old_prep

!!****f* xc/xc_calculate_pw_old [1.0] *
!!
!!   NAME
!!     xc_calculate_pw_old
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
  SUBROUTINE xc_calculate_pw_old(rho_r,rho_g,vxc,exc,dft_control,&
       cell,pw_pool,error)
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout)                          :: vxc
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_pw_old', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: tot_grad = .FALSE.

    INTEGER                                  :: handle, i, idir, ii, ipot, &
                                                ispin, j, k, ngrad_spins, &
                                                nspins, ntot, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3)                    :: n
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: crossterms, failure, &
                                                gradient_f, my_rho_g_local, &
                                                my_rho_r_local
    REAL(KIND=dp)                            :: alpha, eps_rho
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drgrid, rgrid
    TYPE(pw_p_type), DIMENSION(3)            :: drho_r_att
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    TYPE(pw_type), POINTER                   :: my_rho_g, my_rho_r, tmp_g, &
                                                tmp_r, vxc_g
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(tmp_g, tmp_r, vxc_g, my_rho_g, my_rho_r, drho_r, derivative_set)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    DO idir=1,3
       NULLIFY(drho_r_att(idir)%pw)
    END DO

    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    !FM    CPPrecondition(associated(vxc),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       gradient_f=ANY(dft_control%xc_control%gradient_functionals)
       nspins=dft_control%nspins
       ngrad_spins=nspins
       crossterms=ANY(dft_control%xc_control%crossterms).AND.nspins>1

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

       CALL xc_dset_create(derivative_set, pw_pool=pw_pool, error=error)

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
             CASE (xc_rho_nn50)
                CALL pw_nn_smear_g(my_rho_g,&
                     nn50_coeffs, error=error)
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

                eps_rho = dft_control%xc_control%density_cut
                alpha = 1.00_dp / 2.0_dp               ! XXXXXXXXXXXXXXXXXXXXxx alpha hard-coded here XXXXXXXXXXXXXXX
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
                            ! 1.0_dp
                            ! drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                            ! (2.0_dp) * my_rho_r%cr3d(i,j,k)
                            drho_r(idir,ispin)%pw%cr3d(i,j,k)=drho_r(idir,ispin)%pw%cr3d(i,j,k) * &
                                 (1.0_dp/alpha) *  my_rho_r%cr3d(i,j,k) **(1.0_dp/alpha-1.0_dp)
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

       SELECT CASE (nspins)
       CASE (1)
          IF (gradient_f) THEN
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                  error=error)
          ELSE
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rgrid(:,1), &
                  error=error)
          END IF
       CASE (2)
          IF (gradient_f) THEN
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                  rhob=rgrid(:,2), drhob=drgrid(:,2), &
                  error=error)
          ELSE
             CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, 1, &
                  rhoa=rgrid(:,1), &
                  rhob=rgrid(:,2), &
                  error=error)
          END IF
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
       exc = SUM ( derivative_set%data(:,0) )*rho_r(1)%pw%pw_grid%dvol
       IF ( rho_r(1)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) THEN
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
                   vxc(ispin)%pw%cr3d(i,j,k) = derivative_set%data(ii,ispin)
                END DO
             END DO
          END DO
       END DO

       IF ( ANY(dft_control%xc_control%gradient_functionals) .OR.&
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

             IF (ANY(dft_control%xc_control%gradient_functionals)) THEN ! partial integration
                CALL pw_pool_create_pw ( pw_pool, tmp_g,&
                     use_data=COMPLEXDATA1D,&
                     in_space=RECIPROCALSPACE, error=error)

                ipot=nspins+ispin

                !$omp parallel do private(i)
                DO i = 1, ntot
                   IF ( drgrid(i,ispin) > dft_control%xc_control%gradient_cut ) THEN
                      derivative_set%data(i,ipot) = -derivative_set%data(i,ipot) / drgrid(i,ispin)
                   ELSE
                      derivative_set%data(i,ipot) = 0._dbl
                   END IF
                END DO

                DO idir = 1, 3
                   !$omp parallel do private(k,j,i,ii)
                   DO k = bo(1,3), bo(2,3)
                      ii = (k-bo(1,3))*n(1)*n(2)
                      DO j = bo(1,2), bo(2,2)
                         DO i = bo(1,1), bo(2,1)
                            ii = ii + 1
                            drho_r(idir,ispin)%pw%cr3d(i,j,k) = derivative_set%data(ii,ipot)*&
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
             CASE (xc_rho_nn50)
                CALL pw_nn_smear_g(vxc_g,&
                     nn50_coeffs,error=error)
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
       
       CALL xc_dset_release(derivative_set, error=error)

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xc_calculate_pw_old
  !***************************************************************************

!!****f* xc/xc_rho_set_and_dset_create [1.0] *
!!
!!   NAME
!!     xc_rho_set_and_dset_create
!!
!!   FUNCTION
!!     creates a xc_rho_set and a derivative set containing the derivatives
!!     of the functionals with the given deriv_order.
!!
!!   NOTES
!!     if any of the functionals is gradient corrected the full gradient is
!!     added to the rho set
!!
!!   ARGUMENTS
!!     - rho_set: will contain the rho set
!!     - deriv_set will contain the derivatives
!!     - deriv_order: the order of the requested derivatives. If positive
!!       0:deriv_order are calculated, if negative only -deriv_order is
!!       guaranteed to be valid. Orders not requested might be present,
!!       but might contain garbage.
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (can be null, used only
!!       without smoothing of rho or deriv)
!!     - tau: value of the kinetic density tau on the grid (can be null,
!!       used only with meta functionals)
!!     - dft_control: varius parameter (number of spin...)
!!     - cell: the simulation cell (to scale the derivatives)
!!     - pw_pool: the pool for the grids
!!     - gradient_f: returns true if any of the functionals is gradient 
!!       corrected
!!     - atoms: positions of the atoms to treat the code in a special way
!!     - needs_basic_components: if the basic components of the arguments 
!!       should be kept in rho set (a basic component is for example drho
!!       when with lda a functional needs norm_drho)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_and_dset_create(rho_set,deriv_set,deriv_order,&
       rho_r,rho_g,tau,dft_control,xc_control,cell,pw_pool,gradient_f,atoms,&
       needs_basic_components,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: deriv_order
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(xc_control_type)                    :: xc_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    LOGICAL, INTENT(out), OPTIONAL           :: gradient_f
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, INTENT(in)                      :: needs_basic_components
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_and_dset_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, lsd
    TYPE(xc_rho_cflags_type)                 :: needs

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(gradient_f)) gradient_f=.FALSE.
       !nspins=dft_control%nspins
       nspins=size(rho_r)
       lsd=(nspins/=1)
       IF (lsd) THEN
          CPPrecondition(nspins==2,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       bo = rho_r(1)%pw%pw_grid%bounds_local

       ! create a place where to put the derivatives
       CALL xc_dset_create(deriv_set, pw_pool, error=error)

       ! create the place where to store the argument for the functionals
       CALL xc_rho_set_create(rho_set,bo,rho_cutoff=xc_control%density_cut,&
            drho_cutoff=xc_control%gradient_cut,tau_cutoff=xc_control%tau_cut,&
            error=error)

       ! assesses the argument needed by the functionals
       CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
       DO i=1,SIZE(xc_control%functionals)
          CALL xc_functional_get_info(xc_control%functionals(i),lsd=lsd,&
               needs=needs,&
               error=error)
       END DO

       IF (needs_basic_components) THEN
          IF (lsd) THEN
             needs%rho_spin=.TRUE.
             needs%tau_spin=needs%tau_spin.OR.needs%tau
          ELSE
             needs%rho=.TRUE.
          END IF
          IF (needs%drho.OR.needs%drho_spin.OR.needs%norm_drho.OR.&
               needs%norm_drho_spin) THEN
             IF (PRESENT(gradient_f)) gradient_f=.TRUE.
             IF (lsd) THEN
                needs%drho_spin=.TRUE.
             ELSE
                needs%drho=.TRUE.
             END IF
          END IF
       END IF

       ! calculate the agruments needed by the functionals
       CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs,&
            dft_control,cell,pw_pool, atoms=atoms,error=error)

       ! evaluate the functionals
       DO i=1,SIZE(xc_control%functionals)
          CALL xc_functional_eval(xc_control%functionals(i), &
               lsd=lsd,&
               rho_set=rho_set, &
               deriv_set=deriv_set,&
               deriv_order=deriv_order, &
               error=error)
       END DO
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xc_rho_set_and_dset_create
  !***************************************************************************

!!****f* xc/smooth_cutoff [1.0] *
!!
!!   NAME
!!     smooth_cutoff
!!
!!   FUNCTION
!!     smooths the cutoff on rho with a function smooth(rho) that is 0 
!!     for rho<rho_cutoff and 1 for rho>rho_cutoff*rho_smooth_cutoff_range:
!!     E= integral e_0*smooth(rho) => dE/d...= de/d... * smooth,
!!     dE/drho = de/drho * smooth + e_0 * dsmooth/drho
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pot: the potential to smooth
!!     - rho, rhoa,rhob: the value of the density (used to apply the cutoff)
!!     - rho_cutoff: the vaule at whch the cutoff function must go to 0
!!     - rho_smooth_cutoff_range: range of the smoothing
!!     - e_0: value of e_0, if given it is assumed that pot is the derivative 
!!       wrt. to rho, and needs the dsmooth*e_0 contribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE smooth_cutoff(pot,rho,rhoa,rhob,rho_cutoff,&
       rho_smooth_cutoff_range, e_0, e_0_scale_factor,error)
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot, rho, rhoa, rhob
    REAL(kind=dp), INTENT(in)                :: rho_cutoff, &
                                                rho_smooth_cutoff_range
    REAL(kind=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: e_0
    REAL(kind=dp), INTENT(in), OPTIONAL      :: e_0_scale_factor
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'smooth_cutoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp) :: my_e_0_scale_factor, my_rho, my_rho_n, my_rho_n2, &
      rho_smooth_cutoff, rho_smooth_cutoff_2, rho_smooth_cutoff_range_2

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pot),cp_failure_level,routineP,error,failure)
    bo(1,:)=LBOUND(pot)
    bo(2,:)=UBOUND(pot)
    my_e_0_scale_factor=1.0_dp
    IF (PRESENT(e_0_scale_factor)) my_e_0_scale_factor=e_0_scale_factor
    IF (.NOT. failure) THEN
       rho_smooth_cutoff=rho_cutoff*rho_smooth_cutoff_range
       rho_smooth_cutoff_2=(rho_cutoff+rho_smooth_cutoff)/2
       rho_smooth_cutoff_range_2=rho_smooth_cutoff_2-rho_cutoff

       IF (rho_smooth_cutoff_range>0.0_dp) THEN
          IF (PRESENT(e_0)) THEN
             CPPrecondition(ASSOCIATED(e_0),cp_failure_level,routineP,error,failure)
             IF (ASSOCIATED(rho)) THEN
                !$omp parallel do default(none) shared(bo,e_0,pot,rho,&
                !$omp             rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
                !$omp             rho_smooth_cutoff_range_2,my_e_0_scale_factor)&
                !$omp       private(k,j,i,my_rho,my_rho_n,my_rho_n2)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         my_rho=rho(i,j,k)
                         IF (my_rho<rho_smooth_cutoff) THEN
                            IF (my_rho<rho_cutoff) THEN
                               pot(i,j,k)=0.0_dp
                            ELSEIF (my_rho<rho_smooth_cutoff_2) THEN
                               my_rho_n=(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2)+&
                                    my_e_0_scale_factor*e_0(i,j,k)*&
                                    my_rho_n2*(3.0_dp-2.0_dp*my_rho_n)&
                                    /rho_smooth_cutoff_range_2
                            ELSE
                               my_rho_n=2.0_dp-(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    (1.0_dp-my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2))&
                                    +my_e_0_scale_factor*e_0(i,j,k)*&
                                    my_rho_n2*(3.0_dp-2.0_dp*my_rho_n)&
                                    /rho_smooth_cutoff_range_2
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
             ELSE
                !$omp parallel do default(none) shared(bo,pot,e_0,rhoa,rhob,&
                !$omp             rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
                !$omp             rho_smooth_cutoff_range_2,my_e_0_scale_factor)&
                !$omp       private(k,j,i,my_rho,my_rho_n,my_rho_n2)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         my_rho=rhoa(i,j,k)+rhob(i,j,k)
                         IF (my_rho<rho_smooth_cutoff) THEN
                            IF (my_rho<rho_cutoff) THEN
                               pot(i,j,k)=0.0_dp
                            ELSEIF (my_rho<rho_smooth_cutoff_2) THEN
                               my_rho_n=(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2)+&
                                    my_e_0_scale_factor*e_0(i,j,k)*&
                                    my_rho_n2*(3.0_dp-2.0_dp*my_rho_n)&
                                    /rho_smooth_cutoff_range_2
                            ELSE
                               my_rho_n=2.0_dp-(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    (1.0_dp-my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2))&
                                    +my_e_0_scale_factor*e_0(i,j,k)*&
                                    my_rho_n2*(3.0_dp-2.0_dp*my_rho_n)&
                                    /rho_smooth_cutoff_range_2
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
             END IF
          ELSE
             IF (ASSOCIATED(rho)) THEN
                !$omp parallel do default(none) shared(bo,pot,&
                !$omp             rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
                !$omp             rho_smooth_cutoff_range_2,rho)&
                !$omp       private(k,j,i,my_rho,my_rho_n,my_rho_n2)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         my_rho=rho(i,j,k)
                         IF (my_rho<rho_smooth_cutoff) THEN
                            IF (my_rho<rho_cutoff) THEN
                               pot(i,j,k)=0.0_dp
                            ELSEIF (my_rho<rho_smooth_cutoff_2) THEN
                               my_rho_n=(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2)
                            ELSE
                               my_rho_n=2.0_dp-(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    (1.0_dp-my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2))
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
             ELSE
                CPPrecondition(ASSOCIATED(rhoa),cp_failure_level,routineP,error,failure)
                CPPrecondition(ASSOCIATED(rhob),cp_failure_level,routineP,error,failure)
                !$omp parallel do default(none) shared(bo,pot,&
                !$omp             rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
                !$omp             rho_smooth_cutoff_range_2,rhoa,rhob)&
                !$omp       private(k,j,i,my_rho,my_rho_n,my_rho_n2)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         my_rho=rhoa(i,j,k)+rhob(i,j,k)
                         IF (my_rho<rho_smooth_cutoff) THEN
                            IF (my_rho<rho_cutoff) THEN
                               pot(i,j,k)=0.0_dp
                            ELSEIF (my_rho<rho_smooth_cutoff_2) THEN
                               my_rho_n=(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2)
                            ELSE
                               my_rho_n=2.0_dp-(my_rho-rho_cutoff)/rho_smooth_cutoff_range_2
                               my_rho_n2=my_rho_n*my_rho_n
                               pot(i,j,k)=pot(i,j,k)*&
                                    (1.0_dp-my_rho_n2*(my_rho_n-0.5_dp*my_rho_n2))
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END IF
    END IF
END SUBROUTINE smooth_cutoff
!***************************************************************************

!!****f* xc/xc_vxc_create [1.0] *
!!
!!   NAME
!!     xc_vxc_create
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   NOTES
!!     Beware: some really dirty pointer handling!
!!     energy should be kept consistent with xc_exc_calc
!!
!!   ARGUMENTS
!!     - vxc_rho: will contain the v_xc part that depend on rho
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - vxc_tau: will contain the kinetic (tau) part of v_xc
!!       (if one of the choosen xc functionals has it it is allocated and you
!!       are responsible for it)
!!     - rho_r: the value of the density in the real space
!!     - rho_g: value of the density in the g space (needs to be associated
!!       only for gradient corrections)
!!     - tau: value of the kinetic density tau on the grid (can be null,
!!       used only with meta functionals)
!!     - exc: the xc energy
!!     - dft_control: varius parameter (number of spin...)
!!     - pw_pool: the pool for the grids
!!     - atoms: a subsys with the atoms, to be able of treating the core
!!       region of the atoms in a special way
!!     - calculate_forces: if we are calculating the forces (and when using
!!       tricks for the core, we should add those forces to the atoms)
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
!!     Fawzi (11.2003): version using the new xc interface
!!     Fawzi (03.2004): fft free for smoothed density and derivs, gga lsd
!!     Fawzi (04.2004): metafunctionals
!!
!!*** **********************************************************************
SUBROUTINE xc_vxc_pw_create(vxc_rho,vxc_tau,exc,rho_r,rho_g,tau,dft_control,&
     cell,pw_pool,atoms,calculate_forces,error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    LOGICAL, OPTIONAL                        :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=30), DIMENSION(2), PARAMETER :: &
      norm_drho_spin_name = (/ "(norm_drhoa)", "(norm_drhob)" /)

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    INTEGER                                  :: handle, i, idir, ispin, j, k, &
                                                n_deriv, npoints, nspins, &
                                                order, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL :: dealloc_pw_to_deriv, failure, gradient_f, has_tau, ldaCore, &
      lsd, my_calculate_forces, zero_result
    REAL(KIND=dp)                            :: drho_cutoff, my_rho, ndr, &
                                                rho_cutoff
    REAL(kind=dp), DIMENSION(4)              :: max_ldaCoreRadius
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: deriv_data, norm_drho, &
                                                norm_drho_spin, rho, rhoa, &
                                                rhob, tmp_cr3d
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      POINTER                                :: drho, drho_spin, drhoa, drhob
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type), DIMENSION(2)            :: vxc_to_deriv
    TYPE(pw_p_type), DIMENSION(3)            :: pw_to_deriv, pw_to_deriv_rho
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r, tmp_r2, vxc_g
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv_att
    TYPE(xc_rho_set_type), POINTER           :: rho_set

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(tmp_g, tmp_r, tmp_r2, vxc_g, norm_drho_spin, norm_drho, drho_spin, drhoa,&
       drhob, pos,deriv_set,rho_set)
  nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  DO idir=1,3
     NULLIFY(pw_to_deriv(idir)%pw, pw_to_deriv_rho(idir)%pw)
  END DO
  DO i=1,2
     NULLIFY(vxc_to_deriv(i)%pw)
  END DO
  my_calculate_forces=.FALSE.
  IF (PRESENT(calculate_forces)) my_calculate_forces=calculate_forces
  ldaCore=.FALSE.
  DO i=1,SIZE(atoms)
     CALL get_atomic_kind_set(atoms(i)%subsys%atomic_kinds%els,&
          max_ldaCoreRadius=max_ldaCoreRadius)
     ldaCore=ldaCore.OR.(max_ldaCoreRadius(2)>0.0_dbl).OR.&
          (max_ldaCoreRadius(4)>0.0_dbl)
  END DO
  pw_grid => rho_r(1)%pw%pw_grid

  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     nspins=dft_control%nspins
     lsd=(nspins/=1)
     IF (lsd) THEN
        CPPrecondition(nspins==2,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  CALL cp_assert(dft_control%integration_grid%pw_default,&
       cp_failure_level,cp_unimplemented_error_nr,routineP,&
       "calculation on a special grid unimplemented "//&
       CPSourceFileRef,&
       error,failure)

  IF (.NOT.failure) THEN
     bo = rho_r(1)%pw%pw_grid%bounds_local
     npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

     ! calculate the potential derivatives
     CALL xc_rho_set_and_dset_create(rho_set=rho_set,deriv_set=deriv_set,&
          deriv_order=1, rho_r=rho_r, rho_g=rho_g, tau=tau,&
          dft_control=dft_control, xc_control=dft_control%xc_control,&
          cell=cell, pw_pool=pw_pool,&
          gradient_f=gradient_f, atoms=atoms, needs_basic_components=.TRUE.,&
          error=error)

  END IF

  IF (.NOT.failure) THEN
     has_tau=.FALSE.
     ! check for unknown derivatives
     n_deriv=0
     pos => deriv_set%derivs
     DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
        CALL xc_derivative_get(deriv_att,order=order,&
             desc=desc,&
             error=error)
        IF (order==1) THEN
           IF (lsd) THEN
              SELECT CASE(desc)
              CASE("(rho)","(rhoa)","(rhob)","(norm_drho)","(norm_drhoa)",&
                   "(norm_drhob)")
                 n_deriv=n_deriv+1
              CASE("(tau)","(tau_a)","(tau_b)")
                 has_tau=.TRUE.
                 n_deriv=n_deriv+1
              CASE default
                 !FM if you are looking at this error probably you are missing the
                 !FM cross term (drhoa_drhob), I never got round to implement it,
                 !FM in the functionals that I have implemented I use norm_drho
                 !FM instead, either do the same or implement it, it shouldn't be
                 !FM difficult (I am just lazy ;)
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "unknown functional derivative (LSD): '"//&
                      TRIM(desc)//"' in "//&
                      CPSourceFileRef,&
                      error,failure)
              END SELECT
           ELSE
              SELECT CASE(desc)
              CASE("(z)")
              CASE("(rho)","(norm_drho)")
                 n_deriv=n_deriv+1
              CASE("(tau)")
                 has_tau=.TRUE.
                 n_deriv=n_deriv+1
              CASE default
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "unknown functional derivative (LDA): '"//&
                      TRIM(desc)//"' in "//&
                      CPSourceFileRef,&
                      error,failure)
              END SELECT
           END IF
        END IF
     END DO
  END IF

  IF (.NOT.failure) THEN
     ALLOCATE(vxc_rho(nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ispin=1,nspins
        NULLIFY(vxc_rho(ispin)%pw)
     END DO

     CALL xc_rho_set_get(rho_set,rho=rho,rhoa=rhoa,rhob=rhob,&
          can_return_null=.TRUE., error=error)
     ! recover the vxc arrays
     IF (lsd) THEN
        deriv_att => xc_dset_get_derivative(deriv_set, "(rhoa)")
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(1)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=deriv_att%deriv_data,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(1)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(1)%pw)
        END IF

        deriv_att => xc_dset_get_derivative(deriv_set, "(rhob)")
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(2)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=deriv_att%deriv_data,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(2)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(2)%pw)
        END IF
     ELSE
        deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(1)%pw,pw_grid=pw_grid,&
                cr3d_ptr=deriv_att%deriv_data,&
                use_data=REALDATA3D,in_space=REALSPACE,&
                error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(1)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(1)%pw)
        END IF
     END IF

     deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
     IF (ASSOCIATED(deriv_att)) THEN
        IF (lsd) THEN ! otherwise already taken care in vxc recovery
           CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                error=error)
           !$omp parallel do default(none) shared(bo,vxc_rho,deriv_data)&
           !$omp private(k,j,i)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    vxc_rho(1)%pw%cr3d(i,j,k) = vxc_rho(1)%pw%cr3d(i,j,k)+deriv_data(i,j,k)
                    vxc_rho(2)%pw%cr3d(i,j,k) = vxc_rho(2)%pw%cr3d(i,j,k)+deriv_data(i,j,k)
                 END DO
              END DO
           END DO
           CALL pw_pool_give_back_cr3d(pw_pool,deriv_att%deriv_data,&
                error=error)
           NULLIFY(deriv_att%deriv_data)
        END IF
     END IF

     ! rhoa,rhob already taken care of in vxc recovery

     deriv_att => xc_dset_get_derivative(deriv_set, "(norm_drho)")
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
             error=error)

        CALL xc_rho_set_get(rho_set,norm_drho=norm_drho,&
             drho=drho,drhoa=drhoa,drhob=drhob,rho_cutoff=rho_cutoff,&
             drho_cutoff=drho_cutoff,&
             can_return_null=.TRUE., error=error)
        IF (my_calculate_forces.AND.ldaCore.AND..FALSE.) THEN
           ! forces of zero core, implemented only for lda
           CALL pw_pool_create_pw(pw_pool, tmp_r,&
                use_data=REALDATA3D, in_space=REALSPACE,&
                error=error)
           IF (ASSOCIATED(norm_drho)) THEN
              !$omp parallel do default(none) shared(bo,deriv_data,norm_drho,&
              !$omp                 drho_cutoff, tmp_r)&
              !$omp             private(k,j,i)
              DO k = bo(1,3), bo(2,3)
                 DO j = bo(1,2), bo(2,2)
                    DO i = bo(1,1), bo(2,1)
                       tmp_r%cr3d(i,j,k)=deriv_data(i,j,k)*&
                            norm_drho(i,j,k)
                    END DO
                 END DO
              END DO
           ELSE IF (ASSOCIATED(drho)) THEN
              !$omp parallel do default(none) shared(bo,deriv_data,drho,&
              !$omp                  drho_cutoff,tmp_r)&
              !$omp             private(k,j,i,ndr)
              DO k = bo(1,3), bo(2,3)
                 DO j = bo(1,2), bo(2,2)
                    DO i = bo(1,1), bo(2,1)
                       ndr=SQRT(drho(1)%array(i,j,k)**2+&
                            drho(2)%array(i,j,k)**2+&
                            drho(3)%array(i,j,k)**2)
                       tmp_r%cr3d(i,j,k)=deriv_data(i,j,k)*ndr
                    END DO
                 END DO
              END DO
           END IF
           CALL zero_core_force(tmp_r,atoms,set_shift=2,error=error)
           CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
        END IF
           
        IF (ASSOCIATED(norm_drho)) THEN
           !$omp parallel do default(none) shared(bo,deriv_data,norm_drho,drho_cutoff)&
           !$omp             private(k,j,i)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    deriv_data(i,j,k)=-deriv_data(i,j,k)/&
                         MAX(norm_drho(i,j,k),drho_cutoff)
                 END DO
              END DO
           END DO
        ELSE IF (ASSOCIATED(drho)) THEN
           !$omp parallel do default(none) shared(bo,deriv_data,drho,drho_cutoff)&
           !$omp             private(k,j,i,ndr)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    ndr=SQRT(drho(1)%array(i,j,k)**2+&
                         drho(2)%array(i,j,k)**2+&
                         drho(3)%array(i,j,k)**2)
                    deriv_data(i,j,k)=deriv_data(i,j,k)/MAX(ndr,drho_cutoff)
                 END DO
              END DO
           END DO
        ELSE
           CPPrecondition(ASSOCIATED(drhoa),cp_failure_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(drhob),cp_failure_level,routineP,error,failure)
           !$omp parallel do default(none) shared(bo,deriv_data,drhoa,drhob,drho_cutoff)&
           !$omp             private(k,j,i,ndr)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    ndr=SQRT((drhoa(1)%array(i,j,k)+drhob(1)%array(i,j,k))**2+&
                         (drhoa(2)%array(i,j,k)+drhob(2)%array(i,j,k))**2+&
                         (drhoa(3)%array(i,j,k)+drhob(3)%array(i,j,k))**2)
                    deriv_data(i,j,k)=-deriv_data(i,j,k)/&
                         MAX(ndr,drho_cutoff)
                 END DO
              END DO
           END DO
        END IF

!FM        CALL smooth_cutoff(pot=deriv_data,rho=rho,rhoa=rhoa,rhob=rhob,&
!FM             rho_cutoff=dft_control%xc_control%density_cut,&
!FM             rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
!FM             error=error)

        IF (ASSOCIATED(drho).AND.ASSOCIATED(drho(1)%array)) THEN
           CPPrecondition(ASSOCIATED(deriv_data),cp_failure_level,routineP,error,failure)
           DO idir=1,3
              CALL pw_create(pw_to_deriv_rho(idir)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=drho(idir)%array,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
              CPPrecondition(ASSOCIATED(drho(idir)%array),cp_failure_level,routineP,error,failure)
              !$omp parallel do default(none) shared(bo,deriv_data,drho,idir)&
              !$omp             private(k,j,i)
              DO k = bo(1,3), bo(2,3)
                 DO j = bo(1,2), bo(2,2)
                    DO i = bo(1,1), bo(2,1)
                       drho(idir)%array(i,j,k)=drho(idir)%array(i,j,k)*&
                            deriv_data(i,j,k)
                    END DO
                 END DO
              END DO
              NULLIFY(drho(idir)%array)
           END DO
        ELSE
           DO idir=1,3
              CALL pw_pool_create_pw(pw_pool,pw_to_deriv_rho(idir)%pw,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
              !$omp parallel do default(none) shared(bo,deriv_data,&
              !$omp                pw_to_deriv_rho,drhoa,drhob,idir)&
              !$omp             private(k,j,i,my_rho)
              DO k = bo(1,3), bo(2,3)
                 DO j = bo(1,2), bo(2,2)
                    DO i = bo(1,1), bo(2,1)
                       my_rho=drhoa(idir)%array(i,j,k)+drhob(idir)%array(i,j,k)
                       pw_to_deriv_rho(idir)%pw%cr3d(i,j,k)=my_rho*&
                            deriv_data(i,j,k)
                    END DO
                 END DO
              END DO
           END DO
        END IF

        CALL pw_pool_give_back_cr3d(pw_pool, deriv_att%deriv_data, error=error)
     END IF

     DO ispin=1,nspins
!FM        IF (n_deriv>0) THEN
!FM           ! smooth cutoff
!FM           deriv_att => xc_dset_get_derivative(deriv_set, "", error=error)
!FM           CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
!FM                error=error)
!FM           CALL smooth_cutoff(pot=vxc_rho(ispin)%pw%cr3d,rho=rho,rhoa=rhoa,rhob=rhob,&
!FM                rho_cutoff=dft_control%xc_control%density_cut,&
!FM                rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
!FM                e_0=deriv_data,error=error)
!FM        END IF

        IF (lsd) THEN
           IF (ispin==1) THEN
              CALL xc_rho_set_get(rho_set,drhoa=drho_spin,&
                   can_return_null=.NOT.gradient_f,error=error)
              CALL xc_rho_set_get(rho_set,norm_drhoa=norm_drho_spin,&
                   can_return_null=.TRUE., error=error)
           ELSE
              CALL xc_rho_set_get(rho_set,drhob=drho_spin,&
                   can_return_null=.NOT.gradient_f,error=error)
              CALL xc_rho_set_get(rho_set,norm_drhob=norm_drho_spin,&
                   can_return_null=.TRUE., error=error)
           END IF

           deriv_att => xc_dset_get_derivative(deriv_set, norm_drho_spin_name(ispin))
           IF (ASSOCIATED(deriv_att)) THEN
              CPPrecondition(lsd,cp_failure_level,routineP,error,failure)
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                   error=error)
              CALL pw_create(tmp_r,pw_grid,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   cr3d_ptr=deriv_data, error=error)

              IF (ASSOCIATED(norm_drho_spin)) THEN
                 !$omp parallel do default(none) shared(bo,deriv_data,norm_drho_spin,drho_cutoff)&
                 !$omp             private(k,j,i)
                 DO k = bo(1,3), bo(2,3)
                    DO j = bo(1,2), bo(2,2)
                       DO i = bo(1,1), bo(2,1)
                          deriv_data(i,j,k)=&
                               -deriv_data(i,j,k)/&
                               MAX(norm_drho_spin(i,j,k),drho_cutoff)
                       END DO
                    END DO
                 END DO
              ELSE
                 !$omp parallel do default(none) shared(bo,deriv_data,drho_spin,drho_cutoff)&
                 !$omp             private(k,j,i, ndr)
                 DO k = bo(1,3), bo(2,3)
                    DO j = bo(1,2), bo(2,2)
                       DO i = bo(1,1), bo(2,1)
                          ndr=SQRT(&
                               drho_spin(1)%array(i,j,k)**2+&
                               drho_spin(2)%array(i,j,k)**2+&
                               drho_spin(3)%array(i,j,k)**2)
                          deriv_data(i,j,k)=&
                               -deriv_data(i,j,k)/MAX(ndr,drho_cutoff)
                       END DO
                    END DO
                 END DO
              END IF

              vxc_to_deriv(ispin)%pw => tmp_r
              NULLIFY(tmp_r,deriv_att%deriv_data)

!FM              CALL smooth_cutoff(pot=vxc_to_deriv(ispin)%pw%cr3d,rho=rho,rhoa=rhoa,rhob=rhob,&
!FM                   rho_cutoff=dft_control%xc_control%density_cut,&
!FM                   rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
!FM                   error=error)

              DO idir=1,3
                 CPPrecondition(ASSOCIATED(drho_spin(idir)%array),cp_failure_level,routineP,error,failure)
                 CPPrecondition(ASSOCIATED(vxc_to_deriv(ispin)%pw%cr3d),cp_failure_level,routineP,error,failure)
                 !$omp parallel do default(none) shared(bo,deriv_data,drho_spin,&
                 !$omp             ispin,idir,vxc_to_deriv,drho_cutoff)&
                 !$omp             private(k,j,i)
                 DO k = bo(1,3), bo(2,3)
                    DO j = bo(1,2), bo(2,2)
                       DO i = bo(1,1), bo(2,1)
                          drho_spin(idir)%array(i,j,k)=&
                               vxc_to_deriv(ispin)%pw%cr3d(i,j,k)*&
                               drho_spin(idir)%array(i,j,k)
                       END DO
                    END DO
                 END DO

                 CALL pw_create(pw_to_deriv(idir)%pw,pw_grid=pw_grid,&
                      cr3d_ptr=drho_spin(idir)%array,&
                      use_data=REALDATA3D,in_space=REALSPACE,&
                      error=error)
                 NULLIFY(drho_spin(idir)%array)
              END DO
              dealloc_pw_to_deriv=.TRUE.
              CALL pw_pool_give_back_pw(pw_pool,vxc_to_deriv(ispin)%pw,error=error)
           END IF
        END IF

        IF (ASSOCIATED(pw_to_deriv_rho(1)%pw)) THEN
           IF (.NOT.ASSOCIATED(pw_to_deriv(1)%pw)) THEN
              DO idir=1,3
                 pw_to_deriv(idir)%pw => pw_to_deriv_rho(idir)%pw
              END DO
              dealloc_pw_to_deriv=(.NOT.lsd).OR.(ispin==2)
              IF (dealloc_pw_to_deriv) THEN
                 DO idir=1,3
                    NULLIFY(pw_to_deriv_rho(idir)%pw)
                 END DO
              END IF
           ELSE
              DO idir=1,3
                 CALL pw_sumup(pw_to_deriv_rho(idir)%pw,pw_to_deriv(idir)%pw)
                 IF (ispin==2) THEN
                    CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv_rho(idir)%pw,&
                         error=error)
                 END IF
              END DO
           END IF
        END IF
        
        IF (ASSOCIATED(pw_to_deriv(1)%pw)) THEN
           ! partial integration
           IF (dft_control%xc_deriv_method_id/=xc_deriv_pw&
                .AND. dft_control%xc_deriv_method_id/=xc_deriv_fd&
                ) THEN
              CALL pw_spline_scale_deriv(pw_to_deriv, cell=cell,&
                   transpose=.TRUE.,&
                   error=error)
           END IF

           IF (dft_control%xc_deriv_method_id==xc_deriv_pw.OR.&
                dft_control%xc_deriv_method_id==xc_deriv_spline3.OR.&
                dft_control%xc_deriv_method_id==xc_deriv_spline2) THEN

              IF (.NOT.ASSOCIATED(vxc_g)) THEN
                 CALL pw_pool_create_pw(pw_pool,vxc_g,&
                      use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                      error=error)
                 zero_result=.TRUE.
              ELSE
                 zero_result=.FALSE.
              END IF

              DO idir = 1,3
                 IF (zero_result .AND. idir==1) THEN
                    tmp_g => vxc_g
                 ELSE
                    CALL pw_pool_create_pw(pw_pool,tmp_g,&
                         use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                         error=error)
                 END IF

                 CALL pw_transfer ( pw_to_deriv(idir)%pw, tmp_g )

                 SELECT CASE(dft_control%xc_deriv_method_id)
                 CASE (xc_deriv_pw)
                    CALL pw_derive ( tmp_g, nd(:,idir) )
                 CASE (xc_deriv_spline2)
                    CALL pw_spline2_interpolate_values_g(tmp_g,error=error)
                    CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                 CASE (xc_deriv_spline3)
                    CALL pw_spline3_interpolate_values_g(tmp_g,error=error)
                    CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                 CASE default
                    CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                 END SELECT

                 IF (zero_result .AND. idir==1) THEN
                    NULLIFY(tmp_g)
                 ELSE
                    CALL pw_sumup ( tmp_g, vxc_g )
                    CALL pw_pool_give_back_pw(pw_pool,tmp_g,error=error)
                 END IF
                 IF(dealloc_pw_to_deriv) THEN
                    CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv(idir)%pw,error=error)
                 END IF
              END DO
              ! transfer vxc in real space
              CALL pw_pool_create_pw(pw_pool, tmp_r,&
                   use_data=REALDATA3D, in_space=REALSPACE,&
                   error=error)
              CALL pw_transfer ( vxc_g, tmp_r )
              CALL pw_sumup ( tmp_r, vxc_rho(ispin)%pw )
              CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
              CALL pw_pool_give_back_pw(pw_pool,vxc_g,error=error)
           ELSE
              IF (ldaCore) THEN
                 CALL pw_pool_create_pw(pw_pool, tmp_r,&
                      use_data=REALDATA3D, in_space=REALSPACE,&
                      error=error)
                 CALL pw_zero(tmp_r)
              ELSE
                 tmp_r => vxc_rho(ispin)%pw
              END IF
              DO idir=1,3
                 SELECT CASE(dft_control%xc_deriv_method_id)
                 CASE (xc_deriv_spline2_smooth)
                    CALL pw_nn_deriv_r ( pw_in=pw_to_deriv(idir)%pw,&
                         pw_out=tmp_r,coeffs=spline2_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_spline3_smooth)
                    CALL pw_nn_deriv_r ( pw_in=pw_to_deriv(idir)%pw,&
                         pw_out=tmp_r,coeffs=spline3_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_nn10_smooth)
                    CALL pw_nn_deriv_r ( pw_in=pw_to_deriv(idir)%pw,&
                         pw_out=tmp_r,coeffs=nn10_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_nn50_smooth)
                    CALL pw_nn_deriv_r ( pw_in=pw_to_deriv(idir)%pw,&
                         pw_out=tmp_r,coeffs=nn50_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_fd)
                    CALL pw_pool_create_pw(pw_pool, tmp_r2,&
                         use_data=REALDATA3D, in_space=REALSPACE,&
                         error=error)
                    CALL pw_derive_fd(pw_to_deriv(idir)%pw,&
                         tmp_r2,nd(:,idir))
                    CALL pw_sumup(tmp_r2,vxc_rho(ispin)%pw)
                    CALL pw_pool_give_back_pw(pw_pool,tmp_r2,error=error)
                 CASE default
                    CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                 END SELECT
                 IF (dealloc_pw_to_deriv) THEN
                    CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv(idir)%pw,error=error)
                 END IF
              END DO
              IF (ldaCore) THEN
                 CALL zero_core(tmp_r,atoms,set_shift=2,error=error)
                 CALL pw_sumup(tmp_r,vxc_rho(ispin)%pw)
                 CALL pw_pool_give_back_pw(pw_pool,tmp_r,error=error)
              ELSE
                 NULLIFY(tmp_r)
              END IF
           END IF
        END IF

        CALL zero_core(vxc_rho(ispin)%pw,atoms,set_shift=0,error=error)

        IF (pw_grid%spherical) THEN
           ! filter vxc
           CALL pw_pool_create_pw(pw_pool,vxc_g,&
                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                error=error)
           CALL pw_transfer ( vxc_rho(ispin)%pw, vxc_g )
           CALL pw_transfer ( vxc_g,vxc_rho(ispin)%pw )
           CALL pw_pool_give_back_pw(pw_pool,vxc_g)
        END IF

        CALL smooth_cutoff(pot=vxc_rho(ispin)%pw%cr3d,rho=rho,rhoa=rhoa,rhob=rhob,&
             rho_cutoff=dft_control%xc_control%density_cut*&
             dft_control%xc_control%density_smooth_cut_range,&
             rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
             error=error)

        ! final smoothing if rho was smoothed
        IF (dft_control%xc_rho_smooth_id/=xc_rho_no_smooth) THEN
           CALL pw_pool_create_pw(pw_pool, tmp_r,&
                use_data=REALDATA3D, in_space=REALSPACE,&
                error=error)
           CALL pw_zero(tmp_r)
           SELECT CASE(dft_control%xc_rho_smooth_id)
           CASE (xc_rho_spline2_smooth)
              CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw,pw_out=tmp_r,&
                   coeffs=spline2_coeffs,error=error)
           CASE (xc_rho_spline3_smooth)
              CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw,pw_out=tmp_r,&
                   coeffs=spline3_coeffs,error=error)
           CASE (xc_rho_nn10)
              CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw,pw_out=tmp_r,&
                   coeffs=nn10_coeffs,error=error)
           CASE (xc_rho_nn50)
              CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw,pw_out=tmp_r,&
                   coeffs=nn50_coeffs,error=error)
           CASE default
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
           END SELECT
           deriv_data => vxc_rho(ispin)%pw%cr3d
           vxc_rho(ispin)%pw%cr3d => tmp_r%cr3d
           tmp_r%cr3d => deriv_data
           CALL pw_pool_give_back_pw(pw_pool,tmp_r,error=error)
        END IF
     END DO

     DO idir=1,3
        CPPostcondition(.NOT.ASSOCIATED(pw_to_deriv(idir)%pw),cp_warning_level,routineP,error,failure)
        CPPostcondition(.NOT.ASSOCIATED(pw_to_deriv_rho(idir)%pw),cp_warning_level,routineP,error,failure)
     END DO

     ! 0-deriv -> value of exc
     ! this has to be kept consistent with xc_exc_calc
     IF (n_deriv>0) THEN
        deriv_att => xc_dset_get_derivative(deriv_set, "", error=error)
        CPPrecondition(ASSOCIATED(deriv_att),cp_failure_level,routineP,error,failure)
        CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
             error=error)

        CALL pw_create(tmp_r,pw_grid,&
             use_data=REALDATA3D,in_space=REALSPACE,&
             cr3d_ptr=deriv_data, error=error)
!FM        IF (my_calculate_forces) THEN
!FM           CALL zero_core_force(tmp_r,atoms,error)
!FM        END IF
!FM        CALL zero_core(tmp_r,atoms,error=error)
        NULLIFY(tmp_r%cr3d)
        CALL pw_release(tmp_r)

        CALL smooth_cutoff(pot=deriv_data,rho=rho,rhoa=rhoa,rhob=rhob,& 
             rho_cutoff=dft_control%xc_control%density_cut,& 
             rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
             error=error)

        exc = SUM ( deriv_data )*pw_grid%dvol
        IF ( pw_grid%para%mode == PW_MODE_DISTRIBUTED ) THEN
           CALL mp_sum ( exc, pw_grid%para%group )
        END IF
     ELSE
        exc=0.0_dp
     END IF

     CALL xc_rho_set_release(rho_set, pw_pool=pw_pool, error=error)

     ! tau part
     IF (has_tau) THEN
        ALLOCATE(vxc_tau(nspins), stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin=1,nspins
           NULLIFY(vxc_tau(ispin)%pw)
        END DO
        IF (lsd) THEN
           deriv_att => xc_dset_get_derivative(deriv_set, "(tau_a)", error=error)
           IF (ASSOCIATED(deriv_att)) THEN
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                   error=error)
              
              CALL pw_create(vxc_tau(1)%pw,pw_grid,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   cr3d_ptr=deriv_data, error=error)
              NULLIFY(deriv_att%deriv_data)
           END IF
           deriv_att => xc_dset_get_derivative(deriv_set, "(tau_b)", error=error)
           IF (ASSOCIATED(deriv_att)) THEN
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                   error=error)
              
              CALL pw_create(vxc_tau(2)%pw,pw_grid,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   cr3d_ptr=deriv_data, error=error)
              NULLIFY(deriv_att%deriv_data)
           END IF
           deriv_att => xc_dset_get_derivative(deriv_set, "(tau)", error=error)
           IF (ASSOCIATED(deriv_att)) THEN
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                   error=error)
              IF (ASSOCIATED(vxc_tau(1)%pw)) THEN
                 DO ispin=1,2
                    CPPrecondition(ASSOCIATED(vxc_tau(ispin)%pw),cp_failure_level,routineP,error,failure)
                    tmp_cr3d => vxc_tau(ispin)%pw%cr3d
                    CALL daxpy(npoints,1.0_dp,deriv_data,1,tmp_cr3d,1)
                 END DO
              ELSE
                 CALL pw_create(vxc_tau(1)%pw,pw_grid,&
                      use_data=REALDATA3D,in_space=REALSPACE,&
                      cr3d_ptr=deriv_data, error=error)
                 NULLIFY(deriv_att%deriv_data)
                 CALL pw_pool_create_pw(pw_pool, vxc_tau(2)%pw,&
                      use_data=REALDATA3D, in_space=REALSPACE,&
                      error=error)
                 CALL pw_copy(vxc_tau(1)%pw,vxc_tau(2)%pw)
              END IF
           END IF
        ELSE
           deriv_att => xc_dset_get_derivative(deriv_set, "(tau)", error=error)
           CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                error=error)
           CALL pw_create(vxc_tau(1)%pw,pw_grid,&
                use_data=REALDATA3D,in_space=REALSPACE,&
                cr3d_ptr=deriv_data, error=error)
           NULLIFY(deriv_att%deriv_data)
        END IF
        DO ispin=1,nspins
           CPPostcondition(ASSOCIATED(vxc_tau(ispin)%pw),cp_failure_level,routineP,error,failure)
        END DO
     END IF
     CALL xc_dset_release(deriv_set, error=error)
  END IF

  CALL timestop(0.0_dp,handle)
  
END SUBROUTINE xc_vxc_pw_create
  !***************************************************************************

!!****f* xc/xc_exc_calc [1.0] *
!!
!!   NAME
!!     xc_exc_calc
!!
!!   FUNCTION
!!     calculates just the exchange and correlation energy
!!     (no vxc)
!!
!!   NOTES
!!     has to be kept consistent with xc_vxc_create
!!
!!   ARGUMENTS
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
FUNCTION xc_exc_calc(rho_r,rho_g,tau,dft_control, cell,pw_pool,atoms,&
     error) RESULT(exc)
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(inout), TARGET                  :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, tau
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_subsystem_p_type), DIMENSION(:)  :: atoms
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(kind=dp)                            :: exc

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_exc_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :, :), POINTER    :: e_0
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_set_type), POINTER           :: rho_set

  CALL timeset(routineN,"I","",handle)
  NULLIFY(rho_set, deriv_set, deriv,e_0)
  failure=.FALSE.
  exc=0.0_dp
  
  ! this has to be consistent with what is done in xc_vxc_create
  CALL xc_rho_set_and_dset_create(rho_set=rho_set,&
       deriv_set=deriv_set,deriv_order=0,&
       rho_r=rho_r,rho_g=rho_g,tau=tau,dft_control=dft_control,&
       xc_control=dft_control%xc_control,&
       cell=cell,pw_pool=pw_pool,atoms=atoms,&
       needs_basic_components=.FALSE.,error=error)
  deriv => xc_dset_get_derivative(deriv_set,"",error=error)
  IF (ASSOCIATED(deriv)) THEN
     CALL xc_derivative_get(deriv, deriv_data=e_0, error=error)
     
     CALL smooth_cutoff(pot=e_0,rho=rho_set%rho,&
          rhoa=rho_set%rhoa,rhob=rho_set%rhob,& 
          rho_cutoff=dft_control%xc_control%density_cut,& 
          rho_smooth_cutoff_range=dft_control%xc_control%density_smooth_cut_range,&
          error=error)
     
     exc = SUM ( e_0 )*rho_r(1)%pw%pw_grid%dvol
     IF ( rho_r(1)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) THEN
        CALL mp_sum ( exc, rho_r(1)%pw%pw_grid%para%group )
     END IF

     CALL xc_rho_set_release(rho_set, pw_pool=pw_pool, error=error)
     CALL xc_dset_release(deriv_set, error=error)
  END IF
  CALL timestop(0.0_dp,handle)
END FUNCTION xc_exc_calc
!***************************************************************************

!!****f* xc/xc_prep_2nd_deriv1 [1.0] *
!!
!!   NAME
!!     xc_prep_2nd_deriv1
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
!!*** **********************************************************************
  SUBROUTINE xc_prep_2nd_deriv1(deriv_set,rho_r,rho_g,dft_control,pw_pool,&
       drho_r,spin_pot,grad_pot,ndiag_term,error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho_r, rho_g
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    INTEGER, DIMENSION(:, :), POINTER        :: spin_pot
    LOGICAL, DIMENSION(:, :), POINTER        :: grad_pot
    LOGICAL, DIMENSION(:), POINTER           :: ndiag_term
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_prep_2nd_deriv1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_g, my_rho_r

    failure=.FALSE.
    NULLIFY(my_rho_r, my_rho_g)

    IF (ASSOCIATED(rho_r)) THEN
       ALLOCATE(my_rho_r(SIZE(rho_r)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,SIZE(rho_r)
             CPPrecondition(rho_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho_r(ispin)%pw => rho_r(ispin)%pw
          END DO
       END IF
    END IF
    IF (ASSOCIATED(rho_g)) THEN
       ALLOCATE(my_rho_g(SIZE(rho_g)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,SIZE(rho_g)
             CPPrecondition(rho_g(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho_g(ispin)%pw => rho_g(ispin)%pw
          END DO
       END IF
    END IF
    IF (.NOT. failure) THEN
       CALL xc_prep_2nd_deriv(deriv_set,rho_r=my_rho_r,rho_g=my_rho_g,&
            dft_control=dft_control,pw_pool=pw_pool,&
            drho_r=drho_r,spin_pot=spin_pot,grad_pot=grad_pot,&
            ndiag_term=ndiag_term,error=error)
    END IF
  END SUBROUTINE xc_prep_2nd_deriv1
  !***************************************************************************

!!****f* qs_kpp1_env_methods/xc_prep_2nd_deriv [1.0] *
!!
!!   NAME
!!     xc_prep_2nd_deriv
!!
!!   FUNCTION
!!     Does the first prepartions for the calculation of the 2nd deriv
!!     done by xc_calc_2nd_deriv
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - deriv_set: contains the result
!!     - rho_r,rho_g: point where to evaluate the (functional) derivative
!!       of exc in r and g space
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - drho_r(idir,ispin) :: derivative of rho_r in the direction x,y,z
!!       in the r space (allocated if the xc potential is gradient corrected)
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (allocated, see xc_create_2nd_deriv_info)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_prep_2nd_deriv(deriv_set, &
                               rho_r, rho_g, dft_control, pw_pool, &
                               drho_r, spin_pot, grad_pot, ndiag_term, &
                               error)


    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r, rho_g
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    INTEGER, DIMENSION(:, :), POINTER        :: spin_pot
    LOGICAL, DIMENSION(:, :), POINTER        :: grad_pot
    LOGICAL, DIMENSION(:), POINTER           :: ndiag_term
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_prep_2nd_deriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, ii, ipot, j, &
                                                k, n(3), nspins, ntot, spin, &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: do_lsd, do_restricted_tddfpt, &
                                                failure, use_crossterms
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drgrid, rgrid
    TYPE(pw_type), POINTER                   :: tmp_g
    TYPE(xc_control_type)                    :: xc_control

! ------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(tmp_g)
    NULLIFY(drho_r, spin_pot, grad_pot, ndiag_term)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       IF (dft_control%do_tddfpt_calculation) THEN
          xc_control = dft_control%tddfpt_control%xc_control
       ELSE
          xc_control = dft_control%xc_control
       ENDIF

       nspins     = dft_control%nspins
       do_restricted_tddfpt  = (dft_control%do_tddfpt_calculation .AND. (nspins==1))
       do_lsd     = (nspins==2) .OR. do_restricted_tddfpt
       use_crossterms = do_lsd .AND. ANY(xc_control%crossterms)

       IF (ANY(xc_control%gradient_functionals)) THEN
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

       IF (.NOT.ASSOCIATED(deriv_set)) &
            CALL xc_dset_create(deriv_set, pw_pool, error=error)

       ALLOCATE(rgrid(ntot,nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ANY(xc_control%gradient_functionals)) THEN
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


       IF (ANY(xc_control%gradient_functionals)) THEN
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

       IF (nspins == 1) THEN

          !-----------------------------------------------!
          ! do lsd but with rhoa=1/2*rho and rhob=1/2*rho !
          ! and drhoa=1/2*drho and drhob=1/2*rhob         !
          !-----------------------------------------------!
          IF (do_restricted_tddfpt) THEN

             CALL dscal (ntot, 0.5_dp, rgrid(:,1), 1)
             IF (ANY(xc_control%gradient_functionals)) THEN
                CALL dscal(ntot, 0.5_dp, drgrid(:,1), 1)
                CALL xc_calc_derivatives(deriv_set, dft_control%xc_control, -2, &
                                         rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                                         rhob=rgrid(:,1), drhob=drgrid(:,1), &
                                         error=error)
             ELSE
                CALL xc_calc_derivatives(deriv_set, dft_control%xc_control, -2, &
                                         rhoa=rgrid(:,1), &
                                         rhob=rgrid(:,1), &
                                         error=error)
             END IF

          ELSE
             IF (ANY(xc_control%gradient_functionals)) THEN
                CALL dscal(ntot, 0.5_dp, drgrid(:,1), 1)
                CALL xc_calc_derivatives(deriv_set, dft_control%xc_control, -2, &
                                         rhoa=rgrid(:,1), drhoa=drgrid(:,1), &
                                         error=error)
             ELSE
                CALL xc_calc_derivatives(deriv_set, dft_control%xc_control, -2, &
                                         rhoa=rgrid(:,1), &
                                         error=error)
             END IF
          END IF
       ELSE ! nspins == 2
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='lsd not yet implemented '//&
               CPSourceFileRef,&
               error=error)
       END IF

       CALL xc_create_2nd_deriv_info( lsd=do_lsd, gradient=ANY(xc_control%gradient_functionals),&
            crossterms=use_crossterms, spin_pot=spin_pot, grad_pot=grad_pot,&
            ndiag_term=ndiag_term, error=error)

    END IF

    ! dealloc rgrid
    IF (.NOT.failure) THEN
       DEALLOCATE(rgrid,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    ! divide the potential by the gradient where necessary
    IF (.NOT.failure) THEN
       DO ipot=1,SIZE(grad_pot,2)
          IF (grad_pot(1,ipot).AND.grad_pot(2,ipot)) THEN
             CPPrecondition(ANY(xc_control%gradient_functionals),cp_failure_level,routineP,error,failure)
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                !$omp parallel do private(ii)
                DO ii=1,ntot
                   IF (drgrid(ii,spin_pot(1,ipot))<dft_control%xc_control%gradient_cut.OR.&
                        drgrid(ii,spin_pot(2,ipot))<dft_control%xc_control%gradient_cut) THEN
                      deriv_set%data(ii,ipot) = 0.0_dp
                   ELSE
                      deriv_set%data(ii,ipot) = &
                           deriv_set%data(ii,ipot) / &
                           (drgrid(ii,spin_pot(1,ipot)) * &
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
                   IF (drgrid(ii,spin_pot(1,ipot))<dft_control%xc_control%gradient_cut) THEN
                      deriv_set%data(ii,ipot) = 0.0_dp
                   ELSE
                      deriv_set%data(ii,ipot) = &
                           deriv_set%data(ii,ipot) / &
                           drgrid(ii,spin_pot(1,ipot))**2
                   END IF
                END DO
             END IF
          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(ANY(xc_control%gradient_functionals),cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)
             !$omp parallel do private(ii)
             DO ii=1,ntot
                IF (drgrid(ii,spin_pot(1,ipot))<dft_control%xc_control%gradient_cut.OR.&
                     drgrid(ii,spin_pot(2,ipot))<dft_control%xc_control%gradient_cut) THEN
                   deriv_set%data(ii,ipot) = 0.0_dp
                ELSE
                   deriv_set%data(ii,ipot) = &
                        deriv_set%data(ii,ipot) / &
                        (drgrid(ii,spin_pot(1,ipot)) * &
                        drgrid(ii,spin_pot(2,ipot)))
                END IF
             END DO
          ELSE IF (grad_pot(2,ipot)) THEN
             CPPrecondition(ANY(xc_control%gradient_functionals),cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

             !$omp parallel do private(ii)
             DO ii=1,ntot
                IF (drgrid(ii,spin_pot(1,ipot))<dft_control%xc_control%gradient_cut) THEN
                   deriv_set%data(ii,ipot) = 0.0_dp
                ELSE
                   deriv_set%data(ii,ipot) = &
                        deriv_set%data(ii,ipot) / &
                        drgrid(ii,spin_pot(2,ipot))**2
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

    IF (ANY(xc_control%gradient_functionals)) THEN
       IF (.NOT.failure) THEN
          ! dealloc drgrid
          DEALLOCATE(drgrid,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE xc_prep_2nd_deriv
  !***************************************************************************

!!****f* xc/xc_calc_2nd_deriv1 [1.0] *
!!
!!   NAME
!!     xc_calc_2nd_deriv1
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
  SUBROUTINE xc_calc_2nd_deriv1(deriv_set, &
       rho_r,drho_r,rho1_r,rho1_g,&
       spin_pot,grad_pot,ndiag_term,dft_control,&
       v_rspace_new,pw_pool,error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho_r
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho1_r, rho1_g
    INTEGER, DIMENSION(:, :), POINTER        :: spin_pot
    LOGICAL, DIMENSION(:, :), POINTER        :: grad_pot
    LOGICAL, DIMENSION(:), POINTER           :: ndiag_term
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace_new
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho1_g, my_rho1_r, my_rho_r

    failure=.FALSE.
    NULLIFY(my_rho_r, drho_r, my_rho1_r, my_rho1_g)

    IF (ASSOCIATED(rho_r)) THEN
       ALLOCATE(my_rho_r(SIZE(rho_r)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,SIZE(rho_r)
             CPPrecondition(rho_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho_r(ispin)%pw => rho_r(ispin)%pw
          END DO
       END IF
    END IF
    IF (ASSOCIATED(rho1_r)) THEN
       ALLOCATE(my_rho1_r(SIZE(rho1_r)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,SIZE(rho1_r)
             CPPrecondition(rho1_r(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho1_r(ispin)%pw => rho1_r(ispin)%pw
          END DO
       END IF
    END IF
    IF (ASSOCIATED(rho1_g)) THEN
       ALLOCATE(my_rho1_g(SIZE(rho1_g)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,SIZE(rho1_g)
             CPPrecondition(rho1_g(ispin)%in_use==PLANEWAVES,cp_failure_level,routineP,error,failure)
             my_rho1_g(ispin)%pw => rho1_g(ispin)%pw
          END DO
       END IF
    END IF
    IF (.NOT. failure) THEN
       CALL  xc_calc_2nd_deriv(deriv_set, &
            rho_r=my_rho_r,drho_r=drho_r,&
            rho1_r=my_rho1_r,rho1_g=my_rho1_g,&
            spin_pot=spin_pot,grad_pot=grad_pot,&
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
!!     drho_r, spin_pot, grad_pot, ndiag_term, should be as returned by
!!     xc_prep_2nd_deriv.
!!
!!   INPUTS
!!     - deriv_set: object containing the potentials
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
!!*** **********************************************************************
  SUBROUTINE xc_calc_2nd_deriv(deriv_set, &
       rho_r, drho_r, rho1_r, rho1_g, &
       spin_pot, grad_pot, ndiag_term, dft_control, &
       v_rspace_new, pw_pool, error)


    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: drho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho1_r, rho1_g
    INTEGER, DIMENSION(:, :), POINTER        :: spin_pot
    LOGICAL, DIMENSION(:, :), POINTER        :: grad_pot
    LOGICAL, DIMENSION(:), POINTER           :: ndiag_term
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace_new
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, ii, ipot, j, &
                                                k, n(3), npot, nspins, ntot, &
                                                spin, stat
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: do_lsd, do_restricted_tddfpt, &
                                                failure
    REAL(KIND=dp)                            :: diag_factor
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: drho1_r, v_r_to_deriv
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r, v_deriv_g
    TYPE(xc_control_type)                    :: xc_control

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
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

       xc_control = dft_control%xc_control
       IF (dft_control%do_tddfpt_calculation) &
            xc_control = dft_control%tddfpt_control%xc_control

       nspins     = dft_control%nspins
       do_restricted_tddfpt  = (dft_control%do_tddfpt_calculation .AND. (nspins==1))
       do_lsd     = (nspins==2) .OR. do_restricted_tddfpt
       IF (dft_control%tddfpt_control%res_etype=="SINGLET") THEN
          diag_factor = 1.0_dp
       ELSE
          diag_factor = -1.0_dp
       END IF
       npot = SIZE(deriv_set%data,2)

       IF (ANY(xc_control%gradient_functionals)) THEN
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

    IF (ANY(xc_control%gradient_functionals)) THEN

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

             CPPrecondition(ANY(xc_control%gradient_functionals),cp_failure_level,routineP,error,failure)
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
                              deriv_set%data(ii,ipot)
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
                              deriv_set%data(ii,ipot)
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
                              deriv_set%data(ii,ipot)
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
                           deriv_set%data(ii,ipot)
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
                           deriv_set%data(ii,ipot)
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
                           deriv_set%data(ii,ipot)
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
                           deriv_set%data(ii,ipot)
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
                                 rho1_r(1)%pw%cr3d(i,j,k) * deriv_set%data(ii,ipot)
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
                                 diag_factor * rho1_r(1)%pw%cr3d(i,j,k) * deriv_set%data(ii,ipot)
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
                                 deriv_set%data(ii,ipot)
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
                                 deriv_set%data(ii,ipot)
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
                                 deriv_set%data(ii,ipot)
                         END DO
                      END DO
                   END DO
                END IF ! (ndiag_term)

             END IF ! (do_restricted_tddfpt)
          END IF
       END DO
    END IF

    IF (ANY(xc_control%gradient_functionals)) THEN
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
       IF (ANY(xc_control%gradient_functionals)) THEN
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
       END IF
       CALL pw_pool_create_pw(pw_pool,v_deriv_g, &
            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
            error=error)

       DO spin=1,nspins
          CALL pw_transfer(v_rspace_new(spin)%pw,v_deriv_g)

          IF (ANY(xc_control%gradient_functionals)) THEN
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

       IF (ANY(xc_control%gradient_functionals)) THEN
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,error=error)
       END IF
       CALL pw_pool_give_back_pw(pw_pool,v_deriv_g,error=error)

       IF (ANY(xc_control%gradient_functionals)) THEN
          ! dealloc v_r_to_deriv
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(pw_pool,v_r_to_deriv(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(v_r_to_deriv, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE xc_calc_2nd_deriv
  !***************************************************************************

!!****f* qs_kpp1_env_methods/xc_calc_2nd_deriv_new [1.0] *
!!
!!   NAME
!!     xc_calc_2nd_deriv_new
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
!!   INPUTS
!!     - deriv_set: object containing the potentials
!!     - rho_set: object containing the density at which the derivatives were calculated
!!     - rho1_set : object containing the density with which to fold
!!     - v_rspace_new: will contain the new potential.
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_calc_2nd_deriv_new(v_xc, &
       deriv_set, rho_set, rho1_set, pw_pool, dft_control, error)


    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_xc
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set, rho1_set
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, ispin, j, k, &
                                                nspins, stat
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: failure, gradient_f, lsd
    REAL(KIND=dp)                            :: dr1dr
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: deriv_data, rho1, rho1a, rho1b
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      POINTER                                :: drho, drho1, drho1a, drho1b, &
                                                drhoa, drhob
    TYPE(pw_p_type)                          :: v_drho
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r, tmp_r2
    TYPE(xc_control_type)                    :: xc_control
    TYPE(xc_derivative_type), POINTER        :: deriv_att
    TYPE(xc_rho_cflags_type)                 :: needs

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(tmp_g, tmp_r)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(v_xc),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       gradient_f = .FALSE.
       nspins     = dft_control%nspins
       lsd = (nspins /= 1)
       xc_control = dft_control%xc_control
       IF (dft_control%do_tddfpt_calculation) THEN
          xc_control = dft_control%tddfpt_control%xc_control
          lsd = .TRUE.
       END IF
       ALLOCATE(v_xc(nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ispin=1, nspins
          NULLIFY(v_xc(ispin)%pw)
          CALL pw_pool_create_pw(pw_pool,v_xc(ispin)%pw,&
                                 use_data = REALDATA3D,&
                                 in_space = REALSPACE, error=error)
          CALL pw_zero(v_xc(ispin)%pw)
       END DO
    END IF

    ! assesses the argument needed by the functionals
    CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
    DO i=1,SIZE(xc_control%functionals)
       CALL xc_functional_get_info(xc_control%functionals(i),lsd=lsd,&
            needs=needs,&
            error=error)
    END DO
    IF (needs%drho.OR.needs%drho_spin.OR.needs%norm_drho.OR.&
         needs%norm_drho_spin) gradient_f=.TRUE.

    bounds = rho_set%local_bounds
    

    IF (lsd) THEN

       !-------------------!
       ! UNrestricted case !
       !-------------------!

       CALL xc_rho_set_get(rho1_set, rhoa=rho1a, rhob=rho1b, error=error)

       IF (gradient_f) THEN
          CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob, error=error)
          CALL xc_rho_set_get(rho1_set,drhoa=drho1a, drhob=drho1b, error=error)
          CALL pw_pool_create_pw(pw_pool,v_drho%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_zero(v_drho%pw)
       END IF

       deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   if (nspins /= 1) then
                      v_xc(1)%pw%cr3d(i,j,k)=&
                        deriv_data(i,j,k)*rho1a(i,j,k)
                   else 
                      v_xc(1)%pw%cr3d(i,j,k)=&
                           2.0_dp*deriv_data(i,j,k)*rho1a(i,j,k)
                   end if
                END DO
             END DO
          END DO
       END IF
       
       deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   IF (nspins /= 1) THEN
                      v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                           deriv_data(i,j,k)*rho1b(i,j,k)
                      v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                           deriv_data(i,j,k)*rho1a(i,j,k)
                   ELSE IF (dft_control%tddfpt_control%res_etype == "SINGLET") THEN
                      v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                           2.0_dp*deriv_data(i,j,k)*rho1b(i,j,k)
                   ELSE
                      v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) - &
                           2.0_dp*deriv_data(i,j,k)*rho1b(i,j,k)
                   END IF
                END DO
             END DO
          END DO
       END IF
       
       IF (.NOT. dft_control%do_tddfpt_calculation) THEN
          deriv_att=> xc_dset_get_derivative(deriv_set, "(rhob)(rhob)", &
               error=error)
          IF (ASSOCIATED(deriv_att)) THEN
             CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
                  error=error)
             !$omp parallel do private(k,j,i)
             DO k = bounds(1,3), bounds(2,3)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      v_xc(2)%pw%cr3d(i,j,k)=&
                           deriv_data(i,j,k)*rho1a(i,j,k)
                   END DO
                END DO
             END DO
          END IF
       END IF

    ELSE

       !-----------------!
       ! restricted case !
       !-----------------!

       CALL xc_rho_set_get(rho1_set,rho=rho1, error=error)

       IF (gradient_f) THEN
          CALL xc_rho_set_get(rho_set,drho=drho, error=error)
          CALL xc_rho_set_get(rho1_set,drho=drho1, error=error)
          CALL pw_pool_create_pw(pw_pool,v_drho%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_zero(v_drho%pw)
       END IF

       deriv_att=> xc_dset_get_derivative(deriv_set, "(rho)(rho)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   v_xc(1)%pw%cr3d(i,j,k)=&
                        deriv_data(i,j,k)*rho1(i,j,k)
                END DO
             END DO
          END DO
       END IF
       
       deriv_att=> xc_dset_get_derivative(deriv_set, "(rho)(norm_drho)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i,idir,dr1dr)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   dr1dr=0._dp
                   DO idir=1,3
                      dr1dr=dr1dr+drho(idir)%array(i,j,k)*drho1(idir)%array(i,j,k)
                   END DO
                   v_xc(1)%pw%cr3d(i,j,k)=v_xc(1)%pw%cr3d(i,j,k)+&
                        deriv_data(i,j,k)*dr1dr
                   v_drho%pw%cr3d(i,j,k) = -1._dp * deriv_data(i,j,k)*rho1(i,j,k)
                END DO
             END DO
          END DO
       END IF


       deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drho)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i,idir,dr1dr)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   dr1dr=0._dp
                   DO idir=1,3
                      dr1dr=dr1dr+drho(idir)%array(i,j,k)*drho1(idir)%array(i,j,k)
                   END DO
                   v_drho%pw%cr3d(i,j,k) = v_drho%pw%cr3d(i,j,k) - deriv_data(i,j,k)*rho1(i,j,k)
                END DO
             END DO
          END DO
       END IF

       
       deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drho)", &
            error=error)
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
               error=error)
          !$omp parallel do private(k,j,i,idir,dr1dr)
          DO k = bounds(1,3), bounds(2,3)
             DO j = bounds(1,2), bounds(2,2)
                DO i = bounds(1,1), bounds(2,1)
                   dr1dr=0._dp
                   DO idir=1,3
                      dr1dr=dr1dr+drho(idir)%array(i,j,k)*drho1(idir)%array(i,j,k)
                   END DO
                   v_drho%pw%cr3d(i,j,k) = v_drho%pw%cr3d(i,j,k) + &
                        deriv_data(i,j,k)*rho1(i,j,k)/(rho_set%norm_drho(i,j,k))**2
                END DO
             END DO
          END DO

          ! partial integration

          ! this does not work with non orthorombic cells 
          ! (you will have to use a vector of pw with 3 components)
          CALL pw_pool_create_pw(pw_pool,tmp_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)

          DO idir=1,3
             !$omp parallel do private(k,j,i)
             DO k = bounds(1,3), bounds(2,3)
                DO j = bounds(1,2), bounds(2,2)
                   DO i = bounds(1,1), bounds(2,1)
                      tmp_r%cr3d(i,j,k)=drho(idir)%array(i,j,k)*v_drho%pw%cr3d(i,j,k)-&
                           drho1(idir)%array(i,j,k)*deriv_data(i,j,k)
                   END DO
                END DO
             END DO

             SELECT CASE(dft_control%xc_deriv_method_id)
             CASE (xc_deriv_spline2_smooth)
                CALL pw_nn_deriv_r ( pw_in=tmp_r,&
                     pw_out=v_xc(1)%pw,coeffs=spline2_deriv_coeffs,&
                     idir=idir, error=error )
             CASE (xc_deriv_spline3_smooth)
                CALL pw_nn_deriv_r ( pw_in=tmp_r,&
                     pw_out=v_xc(1)%pw,coeffs=spline3_deriv_coeffs,&
                     idir=idir, error=error )
             CASE (xc_deriv_nn10_smooth)
                CALL pw_nn_deriv_r ( pw_in=tmp_r,&
                     pw_out=v_xc(1)%pw,coeffs=nn10_deriv_coeffs,&
                     idir=idir, error=error )
             CASE (xc_deriv_nn50_smooth)
                CALL pw_nn_deriv_r ( pw_in=tmp_r,&
                     pw_out=v_xc(1)%pw,coeffs=nn50_deriv_coeffs,&
                     idir=idir, error=error )
             CASE (xc_deriv_fd)
                CALL pw_pool_create_pw(pw_pool, tmp_r2,&
                     use_data=REALDATA3D, in_space=REALSPACE,&
                     error=error)
                CALL pw_derive_fd(tmp_r,&
                     tmp_r2,nd(:,idir))
                CALL pw_sumup(tmp_r2,v_xc(1)%pw)
                CALL pw_pool_give_back_pw(pw_pool,tmp_r2,error=error)
             CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT

          END DO
          CALL pw_pool_give_back_pw(pw_pool,tmp_r,error=error)

       END IF

    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE xc_calc_2nd_deriv_new
  !***************************************************************************

  SUBROUTINE xc_prep_2nd_deriv_new(deriv_set, &
       rho_set, rho_r, pw_pool, dft_control, cell, error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_prep_2nd_deriv_new', &
      routineP = moduleN//':'//routineN

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), POINTER                  :: split_desc
    INTEGER                                  :: handle, i, idesc, ispin, j, &
                                                k, nspins, order
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, lsd
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(cp_subsystem_p_type), DIMENSION(0)  :: atoms
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_g_pw, rho_r_pw, tau
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r
    TYPE(xc_control_type)                    :: xc_control
    TYPE(xc_derivative_type), POINTER        :: deriv_att

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(tmp_g, tmp_r)

    CPPrecondition(.not.ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       nspins     = dft_control%nspins
       lsd = (nspins /= 1)
       xc_control = dft_control%xc_control
       IF (dft_control%do_tddfpt_calculation) THEN
          xc_control = dft_control%tddfpt_control%xc_control
          lsd = .TRUE.
       END IF
    END IF

    IF (dft_control%do_tddfpt_calculation .AND. nspins == 1) THEN
       ! In this case we need to do a unrestricted calculation to
       ! differenciate beween singlets and triplets
       ALLOCATE(rho_r_pw(2))
       DO ispin=1, 2
          CALL pw_create(rho_r_pw(ispin)%pw, rho_r(1)%pw%pw_grid, &
               rho_r(1)%pw%in_use, rho_r(1)%pw%in_space)
          CALL pw_transfer(rho_r(1)%pw, rho_r_pw(ispin)%pw)
          CALL dscal(PRODUCT(rho_r_pw(ispin)%pw%pw_grid%npts_local), &
               0.5_dp, rho_r_pw(ispin)%pw%cr3d, 1)
       END DO
    ELSE
       ALLOCATE(rho_r_pw(nspins), rho_g_pw(nspins))
       DO ispin=1, nspins
          rho_r_pw(ispin)%pw => rho_r(ispin)%pw
       END DO
    END IF

    NULLIFY(rho_g, tau)
    call xc_rho_set_and_dset_create(rho_set,deriv_set,2,&
       rho_r_pw,rho_g,tau,dft_control,xc_control,cell,pw_pool,atoms=atoms,&
       needs_basic_components=.TRUE.,error=error)

    bo = rho_r(1)%pw%pw_grid%bounds_local

    ! check for unknown derivatives and divide by norm_drho where necessary
    pos => deriv_set%derivs
    DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
       CALL xc_derivative_get(deriv_att,order=order,&
            desc=desc, split_desc=split_desc,&
            error=error)
       if (order==1 .or.order==2) then
          do idesc=1,size(split_desc)
             select case(split_desc(idesc))
             case("norm_drho")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              rho_set%norm_drho(i,j,k)
                      END DO
                   END DO
                END DO
             case("norm_drhoa")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              rho_set%norm_drhoa(i,j,k)
                      END DO
                   END DO
                END DO
             case("norm_drhob")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              rho_set%norm_drhob(i,j,k)
                      END DO
                   END DO
                END DO
             CASE("rho")
                CALL cp_assert(.not.lsd,cp_failure_level,cp_assertion_failed,routineP,&
                     "rho not handled in lsd: '"//&
                     TRIM(desc)//"' in "//&
                     CPSourceFileRef,&
                     error,failure)
             case("rhoa","rhob")
             case default
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "unhandled derivative: '"//&
                     TRIM(split_desc(idesc))//"' in '"//&
                     TRIM(desc)//"' in "//&
                     CPSourceFileRef,&
                     error,failure)
             end select
          end do
       end if
    END DO

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE xc_prep_2nd_deriv_new
  !***************************************************************************

  SUBROUTINE pw_smooth(pw_in,pw_out)
    TYPE(pw_type), POINTER                   :: pw_in, pw_out

    INTEGER                                  :: bo(2,3), i, il, ir, j, jl, &
                                                jr, k, kl, kr, method, n(3), &
                                                nc(3), p, q, r
    REAL(KIND=dp)                            :: alpha, beta, dist, dr(3), &
                                                radius, sigma, sum
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Kernel

    n(1:3) = pw_in%pw_grid%npts_local (1:3)
    dr(:) = pw_in%pw_grid%dr(:)
    bo = pw_in%pw_grid%bounds_local

    method = 1 ! hard coded right now, like everything in here

    SELECT CASE(method)
    CASE(1) ! just some averaging over neighbors, very fast
       alpha=1.0_dp
       beta =0.1_dp
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
       radius=0.5_dp
       sigma =0.1_dp
       nc(:)=CEILING(radius/dr(:))
       WRITE(6,*) radius,sigma,nc(:)
       ALLOCATE(Kernel(-nc(1):nc(1),-nc(2):nc(2),-nc(3):nc(3)))
       sum = 0.0_dp
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
       pw_out%cr3d(:,:,:) = 0.0_dp
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

  END SUBROUTINE pw_smooth

END MODULE xc

!******************************************************************************
