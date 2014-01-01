!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Exchange and Correlation functional calculations
!> \par History
!>      (13-Feb-2001) JGH, based on earlier version of apsi
!>      02.2003 Many many changes [fawzi]
!>      03.2004 new xc interface [fawzi]
!>      04.2004 kinetic functionals [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc
  USE cell_types,                      ONLY: cell_type
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE input_constants,                 ONLY: &
       xc_debug_new_routine, xc_deriv_nn10_smooth, xc_deriv_nn50_smooth, &
       xc_deriv_pw, xc_deriv_spline2, xc_deriv_spline2_smooth, &
       xc_deriv_spline3, xc_deriv_spline3_smooth, xc_new_f_routine, &
       xc_rho_nn10, xc_rho_nn50, xc_rho_no_smooth, xc_rho_spline2_smooth, &
       xc_rho_spline3_smooth, xc_test_lsd_f_routine
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_rval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED,&
                                             pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
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
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             MAX_LABEL_LENGTH
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                             xc_functionals_get_needs
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_get,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xc_vxc_pw_create1, xc_vxc_pw_create, &
       xc_rho_set_and_dset_create, xc_exc_calc,&
       xc_calc_2nd_deriv, xc_prep_2nd_deriv, divide_by_norm_drho

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc'

CONTAINS

! *****************************************************************************
!> \brief Exchange and Correlation functional calculations.
!>      depending on the selected functional_routine calls
!>      the correct routine
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xcthe functional
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param vxc will contain the resulting xc potential, has to be
!>        already allocated
!> \param exc the xc energy
!> \param xc_section parameters selecting the xc and the method used to
!>        calculate it
!> \param pw_pool the pool for the grids
!> \param virial virial for calculating the GGA part of the stress
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xc_vxc_pw_create1(vxc_rho,vxc_tau,rho_r,rho_g,tau,exc,xc_section,&
       cell,pw_pool,error,virial)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau, rho_r, &
                                                rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(virial_type), OPTIONAL, POINTER     :: virial

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f_routine
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(xc_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN

       CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",&
            i_val=f_routine,error=error)
       SELECT CASE(f_routine)
       CASE(xc_new_f_routine)
          CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau,tau=tau,&
               rho_r=rho_r, rho_g=rho_g, exc=exc, xc_section=xc_section,&
               cell=cell, pw_pool=pw_pool,error=error,virial=virial)
       CASE(xc_debug_new_routine)
          CALL xc_vxc_pw_create_debug(vxc_rho=vxc_rho, vxc_tau=vxc_tau,tau=tau,&
               rho_r=rho_r, rho_g=rho_g, exc=exc, xc_section=xc_section,&
               cell=cell, pw_pool=pw_pool, error=error)
       CASE(xc_test_lsd_f_routine)
          CALL xc_vxc_pw_create_test_lsd(vxc_rho=vxc_rho, vxc_tau=vxc_tau,&
               tau=tau, rho_r=rho_r, rho_g=rho_g, exc=exc, &
               xc_section=xc_section, cell=cell, pw_pool=pw_pool,&
               error=error)
       CASE default
       END SELECT
    END IF

  END SUBROUTINE xc_vxc_pw_create1

! *****************************************************************************
!> \brief calculates vxc using lsd with rhoa=rhob=0.5*rho and compares
!>      with the lda result
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section which functional to calculate, and how
!> \param pw_pool the pool for the grids
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      for debugging only: leaks, and non parallel
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE xc_vxc_pw_create_test_lsd(vxc_rho,vxc_tau,rho_r,rho_g,tau,&
     exc,xc_section, cell,pw_pool,  error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau, rho_r, &
                                                rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create_test_lsd', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename
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
            tau=tau,xc_section=xc_section,&
            cell=cell, pw_pool=pw_pool,rho_set=rho_set1,&
            deriv_set=dSet1, deriv_order=1,&
            needs_basic_components=.FALSE.,error=error)
       CALL xc_vxc_pw_create(rho_r=rho_r, rho_g=rho_g,tau=tau,&
            vxc_rho=vxc_rho,vxc_tau=vxc_tau, exc=exc, xc_section=xc_section,&
            cell=cell, pw_pool=pw_pool, &
            error=error)
       PRINT *, "did calculate xc (lda)"
       PRINT *, "about to calculate xc (lsd)"
       CALL xc_rho_set_and_dset_create(rho_set=rho_set2,deriv_set=dSet2,&
            rho_r=rho2_r, rho_g=rho2_g,tau=tau2, xc_section=xc_section,&
            cell=cell, pw_pool=pw_pool, deriv_order=1,&
            needs_basic_components=.FALSE.,error=error)
       CALL xc_vxc_pw_create(rho_r=rho2_r, rho_g=rho2_g,tau=tau2,&
            vxc_rho=vxc_rho2,vxc_tau=vxc_tau2,exc=exc2, xc_section=xc_section,&
            cell=cell, pw_pool=pw_pool,  error=error)
       PRINT *, "did calculate xc (new)"
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
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       IF (ASSOCIATED(vxc_tau2)) THEN
          DEALLOCATE(vxc_tau2,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF

  END IF
END SUBROUTINE xc_vxc_pw_create_test_lsd

! *****************************************************************************
!> \brief calculates vxc outputting the yz plane of rho, and of the various components
!>      of the the derivatives and of vxc
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section which functional should be used, and how to do it
!> \param pw_pool the pool for the grids
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      for debugging only.
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE xc_vxc_pw_create_debug(vxc_rho,vxc_tau,rho_r,rho_g,tau,exc,&
     xc_section, cell,pw_pool,error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau, rho_r, &
                                                rho_g, tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create_debug', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename
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
          tau=tau,xc_section=xc_section,&
          cell=cell, pw_pool=pw_pool,rho_set=rho_set1,&
          deriv_set=dSet1, deriv_order=1,&
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
          exc=exc, xc_section=xc_section,&
          cell=cell, pw_pool=pw_pool,&
          error=error)

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
                       IF (ABS(pot(i,j,k))>10.0_dp) THEN
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
                       IF (ABS(pot(i,j,k))>10.0_dp) THEN
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

! *****************************************************************************
!> \brief creates a xc_rho_set and a derivative set containing the derivatives
!>      of the functionals with the given deriv_order.
!> \param rho_set will contain the rho set
!> \param deriv_set will contain the derivatives
!> \param deriv_order the order of the requested derivatives. If positive
!>        0:deriv_order are calculated, if negative only -deriv_order is
!>        guaranteed to be valid. Orders not requested might be present,
!>        but might contain garbage.
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (can be null, used only
!>        without smoothing of rho or deriv)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param cell the simulation cell (to scale the derivatives)
!> \param xc_section the section describing the functional to use
!> \param pw_pool the pool for the grids
!> \param gradient_f returns true if any of the functionals is gradient
!>        corrected
!> \param needs_basic_components if the basic components of the arguments
!>        should be kept in rho set (a basic component is for example drho
!>        when with lda a functional needs norm_drho)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      if any of the functionals is gradient corrected the full gradient is
!>      added to the rho set
!> \author fawzi
! *****************************************************************************
SUBROUTINE xc_rho_set_and_dset_create(rho_set,deriv_set,deriv_order,&
     rho_r,rho_g,tau,xc_section,cell,pw_pool,&
     needs_basic_components,error)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: deriv_order
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r, rho_g, tau
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    LOGICAL, INTENT(in)                      :: needs_basic_components
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_and_dset_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspins
    LOGICAL                                  :: failure, lsd
    TYPE(section_vals_type), POINTER         :: xc_fun_sections

  CALL timeset(routineN,handle)
  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     nspins=SIZE(rho_r)
     lsd=(nspins/=1)
  END IF

  IF (.NOT.failure) THEN

     xc_fun_sections => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",&
                                                   error=error)

     CALL xc_dset_create(deriv_set, pw_pool, error=error)

     CALL xc_rho_set_create(rho_set,&
            rho_r(1)%pw%pw_grid%bounds_local,&
            rho_cutoff=section_get_rval(xc_section,"density_cutoff",error),&
            drho_cutoff=section_get_rval(xc_section,"gradient_cutoff",error),&
            tau_cutoff=section_get_rval(xc_section,"tau_cutoff",error),&
            error=error)

     CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, &
            xc_functionals_get_needs(xc_fun_sections,lsd,needs_basic_components,error),&
            section_get_ival(xc_section,"XC_GRID%XC_DERIV",error),&
            section_get_ival(xc_section,"XC_GRID%XC_SMOOTH_RHO",error),&
            cell,pw_pool, error=error)

     CALL xc_functionals_eval(xc_fun_sections, &
            lsd=lsd,&
            rho_set=rho_set, &
            deriv_set=deriv_set,&
            deriv_order=deriv_order, &
            error=error)

  END IF

  CALL timestop(handle)

END SUBROUTINE xc_rho_set_and_dset_create

! *****************************************************************************
!> \brief smooths the cutoff on rho with a function smooth(rho) that is 0
!>      for rho<rho_cutoff and 1 for rho>rho_cutoff*rho_smooth_cutoff_range:
!>      E= integral e_0*smooth(rho) => dE/d...= de/d... * smooth,
!>      dE/drho = de/drho * smooth + e_0 * dsmooth/drho
!> \param pot the potential to smooth
!> \param rho , rhoa,rhob: the value of the density (used to apply the cutoff)
!> \param rho_cutoff the vaule at whch the cutoff function must go to 0
!> \param rho_smooth_cutoff_range range of the smoothing
!> \param e_ 0: value of e_0, if given it is assumed that pot is the derivative
!>        wrt. to rho, and needs the dsmooth*e_0 contribution
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE smooth_cutoff(pot,rho,rhoa,rhob,rho_cutoff,&
       rho_smooth_cutoff_range, e_0, e_0_scale_factor,error)
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: pot, rho, rhoa, rhob
    REAL(kind=dp), INTENT(in)                :: rho_cutoff, &
                                                rho_smooth_cutoff_range
    REAL(kind=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: e_0
    REAL(kind=dp), INTENT(in), OPTIONAL      :: e_0_scale_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief Exchange and Correlation functional calculations
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the choosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section which functional to calculate, and how to do it
!> \param pw_pool the pool for the grids
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      Beware: some really dirty pointer handling!
!>      energy should be kept consistent with xc_exc_calc
!> \par History
!>      JGH (13-Jun-2002): adaptation to new functionals
!>      Fawzi (11.2002): drho_g(1:3)->drho_g
!>      Fawzi (1.2003). lsd version
!>      Fawzi (11.2003): version using the new xc interface
!>      Fawzi (03.2004): fft free for smoothed density and derivs, gga lsd
!>      Fawzi (04.2004): metafunctionals
!>      mguidon (12.2008) : laplace functionals
!> \author fawzi; based LDA version of JGH, based on earlier version of apsi
! *****************************************************************************
SUBROUTINE xc_vxc_pw_create(vxc_rho,vxc_tau,exc,rho_r,rho_g,tau,xc_section,&
     cell,pw_pool,error,virial)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    REAL(KIND=dp), INTENT(out)               :: exc
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r, rho_g, tau
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(virial_type), OPTIONAL, POINTER     :: virial

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_vxc_pw_create', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=30), DIMENSION(2), PARAMETER :: &
      norm_drho_spin_name = (/ "(norm_drhoa)", "(norm_drhob)" /)

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    INTEGER :: handle, i, idir, ispin, j, jdir, k, n_deriv, npoints, nspins, &
      order, stat, xc_deriv_method_id, xc_rho_smooth_id
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3, 3)                 :: nd, nd_laplace
    LOGICAL                                  :: dealloc_pw_to_deriv, failure, &
                                                has_laplace, has_tau, lsd, &
                                                use_virial, zero_result
    REAL(KIND=dp)                            :: density_smooth_cut_range, &
                                                drho_cutoff, my_rho, ndr, &
                                                rho_cutoff
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
    TYPE(pw_type), POINTER                   :: tmp_g, tmp_r, virial_pw, vxc_g
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv_att
    TYPE(xc_rho_set_type), POINTER           :: rho_set

  CALL timeset(routineN,handle)
  failure=.FALSE.
  NULLIFY(tmp_g, tmp_r, vxc_g, norm_drho_spin, norm_drho, drho_spin, drhoa, &
       drhob, pos, deriv_set, rho_set, virial_pw)
  nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  DO idir=1,3
     NULLIFY(pw_to_deriv(idir)%pw, pw_to_deriv_rho(idir)%pw)
  END DO
  DO i=1,2
     NULLIFY(vxc_to_deriv(i)%pw)
  END DO

  pw_grid => rho_r(1)%pw%pw_grid

  CPPrecondition(ASSOCIATED(xc_section),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_rho),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(vxc_tau),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     nspins=SIZE(rho_r)
     CPPrecondition(ASSOCIATED(rho_r(SIZE(rho_r))%pw),cp_failure_level,routineP,error,failure)
     lsd=(nspins/=1)
     IF (lsd) THEN
        CPPrecondition(nspins==2,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(virial)) THEN
    use_virial = virial%pv_calculate.AND.(.NOT.virial%pv_numer)
  ELSE
    use_virial = .FALSE.
  ENDIF

  IF (.NOT.failure) THEN
     bo = rho_r(1)%pw%pw_grid%bounds_local
     npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

     ! calculate the potential derivatives
     CALL xc_rho_set_and_dset_create(rho_set=rho_set,deriv_set=deriv_set,&
          deriv_order=1, rho_r=rho_r, rho_g=rho_g, tau=tau,&
          xc_section=xc_section,&
          cell=cell, pw_pool=pw_pool,&
          needs_basic_components=.TRUE.,&
          error=error)

  END IF

  IF (.NOT.failure) THEN
     CALL section_vals_val_get(xc_section,"XC_GRID%XC_DERIV",&
          i_val=xc_deriv_method_id,error=error)
     CALL section_vals_val_get(xc_section,"XC_GRID%XC_SMOOTH_RHO",&
          i_val=xc_rho_smooth_id,error=error)
     CALL section_vals_val_get(xc_section,"DENSITY_SMOOTH_CUTOFF_RANGE",&
          r_val=density_smooth_cut_range,error=error)

     CALL xc_rho_set_get(rho_set,rho_cutoff=rho_cutoff,&
             drho_cutoff=drho_cutoff, error=error)

     has_tau=.FALSE.
     has_laplace = .FALSE.
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
              CASE("(laplace_rhoa)","(laplace_rhob)")
                 has_laplace = .TRUE.
                 n_deriv = n_deriv + 1
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
              CASE("(rho)","(norm_drho)")
                 n_deriv=n_deriv+1
              CASE("(tau)")
                 has_tau=.TRUE.
                 n_deriv=n_deriv+1
              CASE("(laplace_rho)")
                 has_laplace = .TRUE.
                 n_deriv = n_deriv + 1
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
        deriv_att => xc_dset_get_derivative(deriv_set, "(rhoa)",error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(1)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=deriv_att%deriv_data,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(1)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(1)%pw, error=error)
        END IF

        deriv_att => xc_dset_get_derivative(deriv_set, "(rhob)",error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(2)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=deriv_att%deriv_data,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(2)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(2)%pw, error=error)
        END IF

     ELSE
        deriv_att => xc_dset_get_derivative(deriv_set, "(rho)",error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL pw_create(vxc_rho(1)%pw,pw_grid=pw_grid,&
                cr3d_ptr=deriv_att%deriv_data,&
                use_data=REALDATA3D,in_space=REALSPACE,&
                error=error)
           NULLIFY(deriv_att%deriv_data)
        ELSE
           CALL pw_pool_create_pw(pw_pool,vxc_rho(1)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
           CALL pw_zero(vxc_rho(1)%pw, error=error)
        END IF
     END IF

     deriv_att => xc_dset_get_derivative(deriv_set, "(rho)",error=error)
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

     deriv_att => xc_dset_get_derivative(deriv_set, "(norm_drho)",error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
             error=error)

        CALL xc_rho_set_get(rho_set,norm_drho=norm_drho,&
             drho=drho,drhoa=drhoa,drhob=drhob,rho_cutoff=rho_cutoff,&
             drho_cutoff=drho_cutoff,&
             can_return_null=.TRUE., error=error)

        IF (ASSOCIATED(norm_drho)) THEN
           !$omp parallel do default(none) shared(bo,deriv_data,norm_drho,drho_cutoff)&
           !$omp             private(k,j,i)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    deriv_data(i,j,k) = -deriv_data(i,j,k)/MAX(norm_drho(i,j,k),drho_cutoff)
                 END DO
              END DO
           END DO
        ELSE IF (ASSOCIATED(drho)) THEN
           !$omp parallel do default(none) shared(bo,deriv_data,drho,drho_cutoff)&
           !$omp             private(k,j,i,ndr)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    ndr = SQRT(drho(1)%array(i,j,k)**2 +&
                               drho(2)%array(i,j,k)**2 +&
                               drho(3)%array(i,j,k)**2)
                    deriv_data(i,j,k) = -deriv_data(i,j,k)/MAX(ndr,drho_cutoff)
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
                    ndr = SQRT((drhoa(1)%array(i,j,k) + drhob(1)%array(i,j,k))**2 +&
                               (drhoa(2)%array(i,j,k) + drhob(2)%array(i,j,k))**2 +&
                               (drhoa(3)%array(i,j,k) + drhob(3)%array(i,j,k))**2)
                    deriv_data(i,j,k) = -deriv_data(i,j,k)/MAX(ndr,drho_cutoff)
                 END DO
              END DO
           END DO
        END IF

        IF (ASSOCIATED(drho).AND.ASSOCIATED(drho(1)%array)) THEN
           CPPrecondition(ASSOCIATED(deriv_data),cp_failure_level,routineP,error,failure)
           IF (use_virial) THEN
             CALL pw_pool_create_pw(pw_pool,virial_pw,&
                                    use_data=REALDATA3D,&
                                    in_space=REALSPACE,&
                                    error=error)
             CALL pw_zero(virial_pw, error=error)
             DO idir=1,3
               DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                   DO i=bo(1,1),bo(2,1)
                     virial_pw%cr3d(i,j,k) = drho(idir)%array(i,j,k)*deriv_data(i,j,k)
                   END DO
                 END DO
               END DO
               DO jdir=1,idir
                 virial%pv_xc(idir,jdir) = pw_grid%dvol*&
                                           accurate_sum(virial_pw%cr3d(:,:,:)*&
                                                        drho(jdir)%array(:,:,:))
                 virial%pv_xc(jdir,idir) = virial%pv_xc(idir,jdir)
               END DO
             END DO
             CALL pw_pool_give_back_pw(pw_pool,virial_pw,error=error)
           END IF ! use_virial
           DO idir=1,3
              CALL pw_create(pw_to_deriv_rho(idir)%pw,pw_grid=pw_grid,&
                   cr3d_ptr=drho(idir)%array,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
              CPPrecondition(ASSOCIATED(drho(idir)%array),cp_failure_level,routineP,error,failure)
              !$omp parallel do default(none) shared(bo,deriv_data,drho,idir,use_virial,virial_pw)&
              !$omp             private(k,j,i)
              DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                    DO i=bo(1,1),bo(2,1)
                       drho(idir)%array(i,j,k) = drho(idir)%array(i,j,k)*deriv_data(i,j,k)
                    END DO
                 END DO
              END DO
              NULLIFY (drho(idir)%array)
           END DO
        ELSE
           IF (use_virial) THEN
             CALL pw_pool_create_pw(pw_pool,virial_pw,&
                                    use_data=REALDATA3D,&
                                    in_space=REALSPACE,&
                                    error=error)
             CALL pw_zero(virial_pw, error=error)
             DO idir=1,3
               DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                   DO i=bo(1,1),bo(2,1)
                     my_rho = drhoa(idir)%array(i,j,k) + drhob(idir)%array(i,j,k)
                     virial_pw%cr3d(i,j,k) = my_rho*deriv_data(i,j,k)
                   END DO
                 END DO
               END DO
               DO jdir=1,idir
                 virial%pv_xc(idir,jdir) = pw_grid%dvol*accurate_sum(virial_pw%cr3d(:,:,:)*&
                                           (drhoa(jdir)%array(:,:,:) + drhob(jdir)%array(:,:,:)))
                 virial%pv_xc(jdir,idir) = virial%pv_xc(idir,jdir)
               END DO
             END DO
             CALL pw_pool_give_back_pw(pw_pool,virial_pw,error=error)
           END IF ! use_virial
           DO idir=1,3
              CALL pw_pool_create_pw(pw_pool,pw_to_deriv_rho(idir)%pw,&
                   use_data=REALDATA3D,in_space=REALSPACE,&
                   error=error)
              !$omp parallel do default(none) shared(bo,deriv_data,&
              !$omp                pw_to_deriv_rho,drhoa,drhob,idir)&
              !$omp             private(k,j,i,my_rho)
              DO k=bo(1,3),bo(2,3)
                 DO j=bo(1,2),bo(2,2)
                    DO i=bo(1,1),bo(2,1)
                       my_rho = drhoa(idir)%array(i,j,k) + drhob(idir)%array(i,j,k)
                       pw_to_deriv_rho(idir)%pw%cr3d(i,j,k) = my_rho*deriv_data(i,j,k)
                    END DO
                 END DO
              END DO
           END DO
        END IF

        CALL pw_pool_give_back_cr3d(pw_pool, deriv_att%deriv_data, error=error)

     END IF

     DO ispin=1,nspins

        IF (lsd) THEN
           IF (ispin==1) THEN
              CALL xc_rho_set_get(rho_set,drhoa=drho_spin,&
                   can_return_null=.TRUE.,error=error)
              CALL xc_rho_set_get(rho_set,norm_drhoa=norm_drho_spin,&
                   can_return_null=.TRUE., error=error)
           ELSE
              CALL xc_rho_set_get(rho_set,drhob=drho_spin,&
                   can_return_null=.TRUE.,error=error)
              CALL xc_rho_set_get(rho_set,norm_drhob=norm_drho_spin,&
                   can_return_null=.TRUE., error=error)
           END IF

           deriv_att => xc_dset_get_derivative(deriv_set, norm_drho_spin_name(ispin),error=error)
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
                          deriv_data(i,j,k) = -deriv_data(i,j,k)/&
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
                          ndr = SQRT(drho_spin(1)%array(i,j,k)**2 +&
                                     drho_spin(2)%array(i,j,k)**2 +&
                                     drho_spin(3)%array(i,j,k)**2)
                          deriv_data(i,j,k) = -deriv_data(i,j,k)/MAX(ndr,drho_cutoff)
                       END DO
                    END DO
                 END DO
              END IF

              IF (use_virial) THEN
                CALL pw_pool_create_pw(pw_pool,virial_pw,&
                                       use_data=REALDATA3D,&
                                       in_space=REALSPACE,&
                                       error=error)
                CALL pw_zero(virial_pw, error=error)
                DO idir=1,3
                  DO k=bo(1,3),bo(2,3)
                    DO j=bo(1,2),bo(2,2)
                      DO i=bo(1,1),bo(2,1)
                        virial_pw%cr3d(i,j,k) = drho_spin(idir)%array(i,j,k)*deriv_data(i,j,k)
                      END DO
                    END DO
                  END DO
                  DO jdir=1,idir
                    virial%pv_xc(idir,jdir) = virial%pv_xc(idir,jdir) + pw_grid%dvol*&
                                              accurate_sum(virial_pw%cr3d(:,:,:)*&
                                                           drho_spin(jdir)%array(:,:,:))
                    virial%pv_xc(jdir,idir) = virial%pv_xc(idir,jdir)
                  END DO
                END DO
                CALL pw_pool_give_back_pw(pw_pool,virial_pw,error=error)
              END IF ! use_virial

              vxc_to_deriv(ispin)%pw => tmp_r
              NULLIFY(tmp_r,deriv_att%deriv_data)

              DO idir=1,3
                 CPPrecondition(ASSOCIATED(drho_spin(idir)%array),cp_failure_level,routineP,error,failure)
                 CPPrecondition(ASSOCIATED(vxc_to_deriv(ispin)%pw%cr3d),cp_failure_level,routineP,error,failure)
                 !$omp parallel do default(none) shared(bo,deriv_data,drho_spin,&
                 !$omp             ispin,idir,vxc_to_deriv,drho_cutoff)&
                 !$omp             private(k,j,i)
                 DO k=bo(1,3),bo(2,3)
                    DO j=bo(1,2),bo(2,2)
                       DO i=bo(1,1),bo(2,1)
                          drho_spin(idir)%array(i,j,k)= vxc_to_deriv(ispin)%pw%cr3d(i,j,k)*&
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
           END IF ! deriv_att

        END IF ! LSD

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
                 CALL pw_axpy(pw_to_deriv_rho(idir)%pw,pw_to_deriv(idir)%pw, error=error)
                 IF (ispin==2) THEN
                    CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv_rho(idir)%pw,&
                         error=error)
                 END IF
              END DO
           END IF
        END IF

        IF (ASSOCIATED(pw_to_deriv(1)%pw)) THEN
           ! partial integration
           IF (xc_deriv_method_id/=xc_deriv_pw) THEN
              CALL pw_spline_scale_deriv(pw_to_deriv, cell=cell,&
                   transpose=.TRUE.,&
                   error=error)
           END IF

           IF (xc_deriv_method_id==xc_deriv_pw.OR.&
                xc_deriv_method_id==xc_deriv_spline3.OR.&
                xc_deriv_method_id==xc_deriv_spline2) THEN

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

                 CALL pw_transfer ( pw_to_deriv(idir)%pw, tmp_g , error=error)

                 SELECT CASE(xc_deriv_method_id)
                 CASE (xc_deriv_pw)
                    CALL pw_derive ( tmp_g, nd(:,idir) , error=error)
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
                    CALL pw_axpy ( tmp_g, vxc_g , error=error)
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
              CALL pw_transfer ( vxc_g, tmp_r , error=error)
              CALL pw_axpy ( tmp_r, vxc_rho(ispin)%pw , error=error)
              CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
              CALL pw_pool_give_back_pw(pw_pool,vxc_g,error=error)
           ELSE
              tmp_r => vxc_rho(ispin)%pw
              DO idir=1,3
                 SELECT CASE(xc_deriv_method_id)
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
                 CASE default
                    CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
                 END SELECT
                 IF (dealloc_pw_to_deriv) THEN
                    CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv(idir)%pw,error=error)
                 END IF
              END DO
              NULLIFY(tmp_r)
           END IF
        END IF

        ! ** Add laplace part to vxc_rho
        IF( has_laplace ) THEN
          nd_laplace = RESHAPE((/2,0,0,0,2,0,0,0,2/),(/3,3/))
          IF(lsd) THEN
            IF( ispin == 1) THEN
              deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rhoa)", error=error)
              CPPrecondition(ASSOCIATED(deriv_att),cp_failure_level,routineP,error,failure)
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                                     error=error)
            ELSE
              deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rhob)", error=error)
              CPPrecondition(ASSOCIATED(deriv_att),cp_failure_level,routineP,error,failure)
              CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                                     error=error)
            END IF

          ELSE
            deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rho)", error=error)
            CPPrecondition(ASSOCIATED(deriv_att),cp_failure_level,routineP,error,failure)
            CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,&
                                   error=error)
          END IF
          DO idir=1,3
            CALL pw_pool_create_pw(pw_pool,pw_to_deriv(idir)%pw,&
                                   use_data=REALDATA3D,in_space=REALSPACE,&
                                   error=error)
            CALL pw_zero(pw_to_deriv(idir)%pw, error=error)
            !$omp parallel do default(none) shared(bo,deriv_data,&
            !$omp                pw_to_deriv,idir)&
            !$omp             private(k,j,i)
            DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                  pw_to_deriv(idir)%pw%cr3d(i,j,k)=deriv_data(i,j,k)
                END DO
              END DO
            END DO
            CALL pw_pool_create_pw(pw_pool,tmp_g,&
                                   use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                                   error=error)
            CALL pw_zero(tmp_g, error=error)
            CALL pw_transfer ( pw_to_deriv(idir)%pw, tmp_g , error=error)

            SELECT CASE(xc_deriv_method_id)
              CASE (xc_deriv_pw)
                CALL pw_derive ( tmp_g, nd_laplace(:,idir) , error=error)
              CASE default
                CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
            END SELECT
            ! Add this to the potential
            CALL pw_pool_create_pw(pw_pool, tmp_r,&
                                   use_data=REALDATA3D, in_space=REALSPACE,&
                                   error=error)
            CALL pw_zero(tmp_r, error=error)
            CALL pw_transfer ( tmp_g, tmp_r , error=error)

            CALL pw_axpy ( tmp_r, vxc_rho(ispin)%pw , error=error)
            CALL pw_pool_give_back_pw(pw_pool, tmp_r, error=error)
            CALL pw_pool_give_back_pw(pw_pool,pw_to_deriv(idir)%pw,error=error)
            CALL pw_pool_give_back_pw(pw_pool,tmp_g,error=error)
          END DO
        END IF


        IF (pw_grid%spherical) THEN
           ! filter vxc
           CALL pw_pool_create_pw(pw_pool,vxc_g,&
                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
                error=error)
           CALL pw_transfer ( vxc_rho(ispin)%pw, vxc_g , error=error)
           CALL pw_transfer ( vxc_g,vxc_rho(ispin)%pw , error=error)
           CALL pw_pool_give_back_pw(pw_pool,vxc_g,error=error)
        END IF
        CALL smooth_cutoff(pot=vxc_rho(ispin)%pw%cr3d,rho=rho,rhoa=rhoa,rhob=rhob,&
             rho_cutoff=rho_cutoff*density_smooth_cut_range,&
             rho_smooth_cutoff_range=density_smooth_cut_range,&
             error=error)

        ! final smoothing if rho was smoothed
        IF (xc_rho_smooth_id/=xc_rho_no_smooth) THEN
           CALL pw_pool_create_pw(pw_pool, tmp_r,&
                use_data=REALDATA3D, in_space=REALSPACE,&
                error=error)
           CALL pw_zero(tmp_r, error=error)
           SELECT CASE(xc_rho_smooth_id)
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
        NULLIFY(tmp_r%cr3d)
        CALL pw_release(tmp_r,error=error)

        CALL smooth_cutoff(pot=deriv_data,rho=rho,rhoa=rhoa,rhob=rhob,&
             rho_cutoff=rho_cutoff,&
             rho_smooth_cutoff_range=density_smooth_cut_range,&
             error=error)

        exc = accurate_sum ( deriv_data )*pw_grid%dvol
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
                 CALL pw_copy(vxc_tau(1)%pw,vxc_tau(2)%pw, error=error)
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

  CALL timestop(handle)

END SUBROUTINE xc_vxc_pw_create

! *****************************************************************************
!> \brief calculates just the exchange and correlation energy
!>      (no vxc)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      has to be kept consistent with xc_vxc_create
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
FUNCTION xc_exc_calc(rho_r,rho_g,tau,xc_section, cell,pw_pool,&
     error) RESULT(exc)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r, rho_g, tau
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: exc

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_exc_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(dp)                                 :: density_smooth_cut_range, &
                                                rho_cutoff
    REAL(dp), DIMENSION(:, :, :), POINTER    :: e_0
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_set_type), POINTER           :: rho_set

  CALL timeset(routineN,handle)
  NULLIFY(rho_set, deriv_set, deriv,e_0)
  failure=.FALSE.
  exc=0.0_dp

  ! this has to be consistent with what is done in xc_vxc_create
  CALL xc_rho_set_and_dset_create(rho_set=rho_set,&
       deriv_set=deriv_set,deriv_order=0,&
       rho_r=rho_r,rho_g=rho_g,tau=tau,xc_section=xc_section,&
       cell=cell,pw_pool=pw_pool,&
       needs_basic_components=.FALSE.,error=error)
  deriv => xc_dset_get_derivative(deriv_set,"",error=error)

  IF (ASSOCIATED(deriv)) THEN
     CALL xc_derivative_get(deriv, deriv_data=e_0, error=error)

     CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",&
          r_val=rho_cutoff,error=error)
     CALL section_vals_val_get(xc_section,"DENSITY_SMOOTH_CUTOFF_RANGE",&
          r_val=density_smooth_cut_range,error=error)
     CALL smooth_cutoff(pot=e_0,rho=rho_set%rho,&
          rhoa=rho_set%rhoa,rhob=rho_set%rhob,&
          rho_cutoff=rho_cutoff,&
          rho_smooth_cutoff_range=density_smooth_cut_range,&
          error=error)

     exc = accurate_sum ( e_0 )*rho_r(1)%pw%pw_grid%dvol
     IF ( rho_r(1)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) THEN
        CALL mp_sum ( exc, rho_r(1)%pw%pw_grid%para%group )
     END IF

     CALL xc_rho_set_release(rho_set, pw_pool=pw_pool, error=error)
     CALL xc_dset_release(deriv_set, error=error)
  END IF
  CALL timestop(handle)
END FUNCTION xc_exc_calc

! *****************************************************************************
!> \brief Calculates the second derivative of E_xc at rho in the direction
!>      rho1  (if you see the second derivative as bilinear form)
!>      partial_rho|_(rho=rho) partial_rho|_(rho=rho) E_xc drho(rho1)drho
!>      The other direction is still indetermined, thus it returns
!>      a potential (partial integration is performed to reduce it to
!>      function of rho, removing the dependence from its partial derivs)
!>      Has to be called after the setup by xc_prep_2nd_deriv.
!> \param deriv_set object containing the potentials
!> \param rho_set object containing the density at which the derivatives were calculated
!> \param rho 1_set : object containing the density with which to fold
!> \param v_rspace_new will contain the new potential.
!> \param pw_pool the pool for the grids
!> \param tddfpt_fac factor that multiplies the crossterms (tddfpt triplets
!>        on a closed shell sistem it should be -1, defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      The old version of this routine was smarter: it handled split_desc(1)
!>      and split_desc(2) separatly, thus the code automatically handled all
!>      possible cross terms (you only had to check if it was diagonal to avoid
!>      double counting). I think that is the way to go if you want to add more
!>      terms (tau,rho in LSD,...). The problem with the old code was that it
!>      because of the old functional structure it sometime guessed wrongly
!>      which derivative was where. There were probably still bugs with gradient
!>      corrected functionals (never tested), and it didn't contain first
!>      derivatives with respect to drho (that contibute also to the second
!>      derivative wrt. rho).
!>      The code was a little complex because it really tried to handle any
!>      functional derivative in the most efficent way with the given contents of
!>      rho_set.
!>      Anyway I strongly encourage whoever wants to modify this code to give a
!>      look to the old version. [fawzi]
! *****************************************************************************
SUBROUTINE xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, rho1_set, &
     pw_pool, xc_section, gapw, vxg, tddfpt_fac, error)

    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_xc
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set, rho1_set
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(section_vals_type), POINTER         :: xc_section
    LOGICAL, INTENT(IN), OPTIONAL            :: gapw
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: vxg
    REAL(kind=dp), INTENT(in), OPTIONAL      :: tddfpt_fac
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv', &
      routineP = moduleN//':'//routineN

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), POINTER                  :: split_desc
    INTEGER                                  :: handle, i, ia, idir, ir, &
                                                ispin, j, k, nspins, order, &
                                                stat, xc_deriv_method_id
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, gradient_f, lsd, &
                                                my_gapw, unknown_deriv
    REAL(KIND=dp)                            :: dr1dr, fac, gradient_cut
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: deriv_data, e_drhoa, e_drhob, &
                                                e_norm_drho, rho1, rho1a, &
                                                rho1b
    TYPE(cp_3d_r_p_type), DIMENSION(:), &
      POINTER                                :: drho, drho1, drho1a, drho1b, &
                                                drhoa, drhob
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(pw_p_type)                          :: v_drho
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: tmp_a, tmp_b, tmp_r, tmp_r2
    TYPE(xc_derivative_type), POINTER        :: deriv_att

! PARAMETER

  CALL timeset(routineN,handle)

  failure=.FALSE.
  NULLIFY(tmp_r, tmp_r2, tmp_a, tmp_b, e_drhoa, e_drhob, e_norm_drho, &
       split_desc)

  my_gapw = .FALSE.
  IF (PRESENT(gapw)) my_gapw = gapw

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(rho1_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(v_xc),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(xc_section),cp_failure_level,routineP,error,failure)
  IF (my_gapw) THEN
     CPPrecondition(PRESENT(vxg),cp_failure_level,routineP,error,failure)
  END IF

  IF (.NOT. failure) THEN
     CALL section_vals_val_get(xc_section,"XC_GRID%XC_DERIV",&
          i_val=xc_deriv_method_id,error=error)
     CALL xc_rho_set_get(rho_set,drho_cutoff=gradient_cut,error=error)
     nspins     = SIZE(v_xc)
     lsd = ASSOCIATED(rho_set%rhoa)
     fac = 0.0_dp
     IF (PRESENT(tddfpt_fac)) fac=tddfpt_fac

     ALLOCATE(tmp_r(nspins), tmp_r2(nspins), tmp_a(nspins), tmp_b(nspins),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ispin=1, nspins
        NULLIFY(tmp_r(ispin)%pw, tmp_r2(ispin)%pw, tmp_a(ispin)%pw, tmp_b(nspins)%pw)
     END DO

  END IF

  bo = rho_set%local_bounds

  gradient_f=.FALSE.
  ! check for unknown derivatives
  pos => deriv_set%derivs
  DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
     unknown_deriv=.FALSE.
     CALL xc_derivative_get(deriv_att,order=order,&
          desc=desc, split_desc=split_desc, error=error)
     SELECT CASE(order)
     CASE(1)
        IF (lsd) THEN
           SELECT CASE(split_desc(1))
           CASE("rho","rhoa","rhob")
           CASE("norm_drho","norm_drhoa","norm_drhob")
              gradient_f=.TRUE.
           CASE default
              unknown_deriv=.TRUE.
           END SELECT
        ELSE
           SELECT CASE(split_desc(1))
           CASE("rho")
           CASE("norm_drho")
              gradient_f=.TRUE.
           CASE default
              unknown_deriv=.TRUE.
           END SELECT
        END IF
     CASE(2)
        IF (lsd) THEN
           SELECT CASE(split_desc(1))
           CASE("rhoa","rhob")
              SELECT CASE(split_desc(2))
              CASE("rhoa","rhob")
              CASE("norm_drhoa", "norm_drhob","norm_drho")
                 gradient_f=.TRUE.
              CASE default
                 unknown_deriv=.TRUE.
              END SELECT
           CASE("norm_drho","norm_drhoa", "norm_drhob")
              gradient_f=.TRUE.
              SELECT CASE(split_desc(2))
              CASE("rhoa","rhob","norm_drhoa", "norm_drhob","norm_drho")
              CASE default
                 unknown_deriv=.TRUE.
              END SELECT
           CASE default
              unknown_deriv=.TRUE.
           END SELECT
        ELSE
           SELECT CASE(split_desc(1))
           CASE("rho")
              SELECT CASE(split_desc(2))
              CASE("rho")
              CASE("norm_drho")
                 gradient_f=.TRUE.
              CASE default
                 unknown_deriv=.TRUE.
              END SELECT
           CASE("norm_drho")
              gradient_f=.TRUE.
              SELECT CASE(split_desc(2))
              CASE("rho","norm_drho")
              CASE default
                 unknown_deriv=.TRUE.
              END SELECT
           CASE default
              unknown_deriv=.TRUE.
           END SELECT
        END IF
     END SELECT
     CALL cp_assert(.NOT.unknown_deriv,cp_failure_level,&
          cp_assertion_failed,routineP,&
          "unknown functional derivative (LSD="//TRIM(cp_to_string(lsd))//&
          "): '"//TRIM(desc)//"' in "//&
CPSourceFileRef,&
          error,failure)
  END DO

  IF (lsd) THEN

     !-------------------!
     ! UNrestricted case !
     !-------------------!

     CALL xc_rho_set_get(rho1_set, rhoa=rho1a, rhob=rho1b, error=error)

     IF (gradient_f) THEN
        CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob, error=error)
        CALL xc_rho_set_get(rho1_set,drhoa=drho1a, drhob=drho1b, error=error)
        DO ispin=1, nspins
           IF (ASSOCIATED(pw_pool)) THEN
              CALL pw_pool_create_pw(pw_pool,tmp_a(ispin)%pw,&
                   use_data = REALDATA3D,&
                   in_space = REALSPACE, error=error)
              CALL pw_zero(tmp_a(ispin)%pw, error=error)
              CALL pw_pool_create_pw(pw_pool,tmp_b(ispin)%pw,&
                   use_data = REALDATA3D,&
                   in_space = REALSPACE, error=error)
              CALL pw_zero(tmp_b(ispin)%pw, error=error)
              CALL pw_pool_create_pw(pw_pool,tmp_r(ispin)%pw,&
                   use_data = REALDATA3D,&
                   in_space = REALSPACE, error=error)
              CALL pw_zero(tmp_r(ispin)%pw, error=error)
           ELSE
              ALLOCATE(tmp_a(ispin)%pw, tmp_b(ispin)%pw, tmp_r(ispin)%pw, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              ALLOCATE(tmp_a(ispin)%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)), &
                   tmp_b(ispin)%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)), &
                   tmp_r(ispin)%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)), &
                   stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              tmp_a(ispin)%pw%cr3d = 0.0_dp
              tmp_b(ispin)%pw%cr3d = 0.0_dp
           END IF
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                      deriv_data(i,j,k) * rho1a(i,j,k)
              END DO
           END DO
        END DO
     END IF

     IF (nspins /= 1) THEN
        deriv_att=> xc_dset_get_derivative(deriv_set, "(rhob)(rhob)", &
             error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
                error=error)
           !$omp parallel do private(k,j,i)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * rho1a(i,j,k)
                 ELSE
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         fac * deriv_data(i,j,k) * rho1b(i,j,k)
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drhoa)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                 END DO
                 v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                      deriv_data(i,j,k) * dr1dr
                 tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                      deriv_data(i,j,k) * rho1a(i,j,k)
              END DO
           END DO
        END DO
     END IF

     IF (nspins /= 1) THEN
        deriv_att=> xc_dset_get_derivative(deriv_set, "(rhob)(norm_drhob)", error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
                error=error)
           !$omp parallel do private(k,j,i,dr1dr,idir)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drhob)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                 END DO
                 IF (nspins /= 1) THEN
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1a(i,j,k)
                 ELSE
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         fac * deriv_data(i,j,k) * dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhob)(norm_drhoa)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                 ELSE
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         fac * deriv_data(i,j,k) * rho1b(i,j,k)
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir) 
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1a(i,j,k)
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1a(i,j,k)
                 ELSE
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr = dr1dr + drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k) + &
                            fac * drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    v_xc(1)%pw%cr3d(i,j,k) = v_xc(1)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1a(i,j,k)
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 ELSE
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k) + &
                            fac * drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rhob)(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    v_xc(2)%pw%cr3d(i,j,k) = v_xc(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * rho1b(i,j,k)
                 ELSE
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         fac * deriv_data(i,j,k) * rho1b(i,j,k)
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhob)(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 ELSE
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         fac*deriv_data(i,j,k) * dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drhoa)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                 END DO
                 tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                      deriv_data(i,j,k) * dr1dr
              END DO
           END DO
        END DO
     END IF

     IF (nspins /= 1) THEN
        deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhob)(norm_drhob)", error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
                error=error)
           !$omp parallel do private(k,j,i,dr1dr,idir)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 END DO
              END DO
           END DO
        END IF
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drhob)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                 END DO
                 IF (nspins /= 1) THEN
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 ELSE
                    tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) - &
                         fac * deriv_data(i,j,k) * dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhoa)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        CALL xc_derivative_get(deriv_att, deriv_data=e_drhoa, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                 END DO
                 tmp_a(1)%pw%cr3d(i,j,k) = tmp_a(1)%pw%cr3d(i,j,k) + &
                      deriv_data(i,j,k) * dr1dr
              END DO
           END DO
        END DO
     END IF

     IF (nspins /= 1) THEN
        deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drhob)", &
             error=error)
        IF (ASSOCIATED(deriv_att)) THEN
           CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
                error=error)
           CALL xc_derivative_get(deriv_att, deriv_data=e_drhob, &
                error=error)
           !$omp parallel do private(k,j,i,dr1dr,idir)
           DO k = bo(1,3), bo(2,3)
              DO j = bo(1,2), bo(2,2)
                 DO i = bo(1,1), bo(2,1)
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(2)%pw%cr3d(i,j,k) = tmp_b(2)%pw%cr3d(i,j,k) + &
                         deriv_data(i,j,k) * dr1dr
                 END DO
              END DO
           END DO
        END IF
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i,dr1dr,idir)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 IF (nspins /= 1) THEN
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                    tmp_a(2)%pw%cr3d(i,j,k) = tmp_a(2)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 ELSE
                    dr1dr=0._dp
                    DO idir=1,3
                       dr1dr=dr1dr+drhob(idir)%array(i,j,k)*drho1a(idir)%array(i,j,k) + &
                            fac * drhoa(idir)%array(i,j,k)*drho1b(idir)%array(i,j,k)
                    END DO
                    tmp_b(1)%pw%cr3d(i,j,k) = tmp_b(1)%pw%cr3d(i,j,k) - &
                         deriv_data(i,j,k) * dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=e_norm_drho, &
             error=error)
     END IF

     IF (gradient_f) THEN

        IF (my_gapw) THEN
           !$omp parallel do private(ia,idir,ispin,ir)
           DO ir = bo(1,2), bo(2,2)
              DO ia = bo(1,1), bo(2,1)
                 DO idir = 1,3
                    DO ispin=1, nspins
                       vxg(idir,ia,ir,ispin) = &
                            tmp_a(ispin)%pw%cr3d(ia,ir,1) * drhoa(idir)%array(ia,ir,1) + &
                            tmp_b(ispin)%pw%cr3d(ia,ir,1) * drhob(idir)%array(ia,ir,1)
                    END DO
                    IF (ASSOCIATED(e_drhoa)) THEN
                       vxg(idir,ia,ir,1) = vxg(idir,ia,ir,1) - &
                            e_drhoa(ia,ir,1) * drho1a(idir)%array(ia,ir,1)
                    END IF
                    IF (nspins /= 1 .AND. ASSOCIATED(e_drhob)) THEN
                       vxg(idir,ia,ir,2) = vxg(idir,ia,ir,2) - &
                            e_drhob(ia,ir,1) * drho1b(idir)%array(ia,ir,1)
                    END IF
                    IF (ASSOCIATED(e_norm_drho)) THEN
                       IF (nspins /= 1) THEN
                          vxg(idir,ia,ir,1) = vxg(idir,ia,ir,1) - &
                               e_norm_drho(ia,ir,1) * drho1b(idir)%array(ia,ir,1)
                          vxg(idir,ia,ir,2) = vxg(idir,ia,ir,2) - &
                               e_norm_drho(ia,ir,1) * drho1a(idir)%array(ia,ir,1)
                       ELSE
                          vxg(idir,ia,ir,1) = vxg(idir,ia,ir,1) - &
                               fac * e_norm_drho(ia,ir,1) * drho1b(idir)%array(ia,ir,1)
                       END IF
                    END IF
                 END DO
              END DO
           END DO
        ELSE

           ! partial integration
           DO idir=1, 3

              DO ispin=1, nspins
                 !$omp parallel do private(k,j,i,dr1dr)
                 DO k = bo(1,3), bo(2,3)
                    DO j = bo(1,2), bo(2,2)
                       DO i = bo(1,1), bo(2,1)
                          tmp_r(ispin)%pw%cr3d(i,j,k) = &
                               tmp_a(ispin)%pw%cr3d(i,j,k) * drhoa(idir)%array(i,j,k) + &
                               tmp_b(ispin)%pw%cr3d(i,j,k) * drhob(idir)%array(i,j,k)
                       END DO
                    END DO
                 END DO
                 IF (ASSOCIATED(e_drhoa)) THEN
                    !$omp parallel do private(k,j,i)
                    DO k = bo(1,3), bo(2,3)
                       DO j = bo(1,2), bo(2,2)
                          DO i = bo(1,1), bo(2,1)
                             tmp_r(1)%pw%cr3d(i,j,k) = tmp_r(1)%pw%cr3d(i,j,k) - &
                                  e_drhoa(i,j,k) * drho1a(idir)%array(i,j,k)
                          END DO
                       END DO
                    END DO
                 END IF
                 IF (nspins /= 1 .AND. ASSOCIATED(e_drhob)) THEN
                    !$omp parallel do private(k,j,i)
                    DO k = bo(1,3), bo(2,3)
                       DO j = bo(1,2), bo(2,2)
                          DO i = bo(1,1), bo(2,1)
                             tmp_r(2)%pw%cr3d(i,j,k) = tmp_r(2)%pw%cr3d(i,j,k) - &
                                  e_drhob(i,j,k) * drho1b(idir)%array(i,j,k)
                          END DO
                       END DO
                    END DO
                 END IF
                 IF (ASSOCIATED(e_norm_drho)) THEN
                    !$omp parallel do private(k,j,i)
                    DO k = bo(1,3), bo(2,3)
                       DO j = bo(1,2), bo(2,2)
                          DO i = bo(1,1), bo(2,1)
                             IF (nspins /= 1) THEN
                                tmp_r(1)%pw%cr3d(i,j,k) = tmp_r(1)%pw%cr3d(i,j,k) - &
                                     e_norm_drho(i,j,k) * drho1b(idir)%array(i,j,k)
                                tmp_r(2)%pw%cr3d(i,j,k) = tmp_r(2)%pw%cr3d(i,j,k) - &
                                     e_norm_drho(i,j,k) * drho1a(idir)%array(i,j,k)
                             ELSE
                                tmp_r(1)%pw%cr3d(i,j,k) = tmp_r(1)%pw%cr3d(i,j,k) - &
                                     fac * e_norm_drho(i,j,k) * drho1b(idir)%array(i,j,k)
                             END IF
                          END DO
                       END DO
                    END DO
                 END IF
              END DO

              DO ispin=1, nspins
                 SELECT CASE(xc_deriv_method_id)
                 CASE (xc_deriv_spline2_smooth)
                    CALL pw_nn_deriv_r ( pw_in=tmp_r(ispin)%pw,&
                         pw_out=v_xc(ispin)%pw,coeffs=spline2_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_spline3_smooth)
                    CALL pw_nn_deriv_r ( pw_in=tmp_r(ispin)%pw,&
                         pw_out=v_xc(ispin)%pw,coeffs=spline3_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_nn10_smooth)
                    CALL pw_nn_deriv_r ( pw_in=tmp_r(ispin)%pw,&
                         pw_out=v_xc(ispin)%pw,coeffs=nn10_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE (xc_deriv_nn50_smooth)
                    CALL pw_nn_deriv_r ( pw_in=tmp_r(ispin)%pw,&
                         pw_out=v_xc(ispin)%pw,coeffs=nn50_deriv_coeffs,&
                         idir=idir, error=error )
                 CASE default
                    CALL cp_unimplemented_error(fromWhere=routineP, &
                            message="XC_DERIV method not implemented for GPW-LSD!",&
                            error=error, error_level=cp_failure_level)
                 END SELECT
              END DO ! ispin

           END DO ! idir

        END IF

        DO ispin=1, nspins
           IF (ASSOCIATED(pw_pool)) THEN
              CALL pw_pool_give_back_pw(pw_pool,tmp_a(ispin)%pw,error=error)
              CALL pw_pool_give_back_pw(pw_pool,tmp_b(ispin)%pw,error=error)
              CALL pw_pool_give_back_pw(pw_pool,tmp_r(ispin)%pw,error=error)
           ELSE
              DEALLOCATE(tmp_a(ispin)%pw%cr3d, tmp_b(ispin)%pw%cr3d, tmp_r(ispin)%pw%cr3d, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(tmp_a(ispin)%pw, tmp_b(ispin)%pw, tmp_r(ispin)%pw, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF
        END DO

     END IF ! gradient_f

  ELSE

     !-----------------!
     ! restricted case !
     !-----------------!

     CALL xc_rho_set_get(rho1_set,rho=rho1, error=error)

     IF (gradient_f) THEN
        CALL xc_rho_set_get(rho_set,drho=drho, error=error)
        CALL xc_rho_set_get(rho1_set,drho=drho1, error=error)
        IF (ASSOCIATED(pw_pool)) THEN
           CALL pw_pool_create_pw(pw_pool,v_drho%pw,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)
           CALL pw_zero(v_drho%pw, error=error)
        ELSE
           ALLOCATE(v_drho%pw, stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(v_drho%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)), &
                stat=stat)
           v_drho%pw%cr3d = 0.0_dp
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(rho)(rho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        !$omp parallel do private(k,j,i)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
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
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
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
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drho(idir)%array(i,j,k)*drho1(idir)%array(i,j,k)
                 END DO
                 v_drho%pw%cr3d(i,j,k) = v_drho%pw%cr3d(i,j,k) - deriv_data(i,j,k)*dr1dr
              END DO
           END DO
        END DO
     END IF

     deriv_att=> xc_dset_get_derivative(deriv_set, "(norm_drho)", &
          error=error)
     IF (ASSOCIATED(deriv_att)) THEN
        CALL xc_derivative_get(deriv_att, deriv_data=deriv_data, &
             error=error)
        CALL xc_derivative_get(deriv_att, deriv_data=e_norm_drho, &
             error=error)
        !$omp parallel do private(k,j,i,idir,dr1dr)
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 dr1dr=0._dp
                 DO idir=1,3
                    dr1dr=dr1dr+drho(idir)%array(i,j,k)*drho1(idir)%array(i,j,k)
                 END DO
                 IF (rho_set%norm_drho(i,j,k) > gradient_cut) THEN
                    dr1dr = dr1dr / (rho_set%norm_drho(i,j,k))**2
                    v_drho%pw%cr3d(i,j,k) = v_drho%pw%cr3d(i,j,k) + deriv_data(i,j,k)*dr1dr
                 END IF
              END DO
           END DO
        END DO
     END IF

     IF (gradient_f) THEN

        IF (my_gapw) THEN

           DO idir=1,3
              !$omp parallel do private(ia,ir)
              DO ia = bo(1,1), bo(2,1)
                 DO ir = bo(1,2), bo(2,2)
                    vxg(idir,ia,ir,1) = drho(idir)%array(ia,ir,1)*v_drho%pw%cr3d(ia,ir,1)
                    IF (ASSOCIATED(e_norm_drho)) THEN
                       vxg(idir,ia,ir,1) = vxg(idir,ia,ir,1) - drho1(idir)%array(ia,ir,1)*e_norm_drho(ia,ir,1)
                    END IF
                 END DO
              END DO
           END DO

        ELSE
           ! partial integration

           ! this does not work with non orthorombic cells
           ! (you will have to use a vector of pw with 3 components)
           IF (ASSOCIATED(pw_pool)) THEN
              CALL pw_pool_create_pw(pw_pool,tmp_r(1)%pw,&
                   use_data = REALDATA3D,&
                   in_space = REALSPACE, error=error)
           ELSE
              ALLOCATE(tmp_r(1)%pw, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              ALLOCATE(tmp_r(1)%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)), &
                   stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF

           DO idir=1,3
              !$omp parallel do private(k,j,i)
              DO k = bo(1,3), bo(2,3)
                 DO j = bo(1,2), bo(2,2)
                    DO i = bo(1,1), bo(2,1)
                       tmp_r(1)%pw%cr3d(i,j,k) = drho(idir)%array(i,j,k)*v_drho%pw%cr3d(i,j,k)-&
                            drho1(idir)%array(i,j,k)*deriv_data(i,j,k)
                    END DO
                 END DO
              END DO

              SELECT CASE(xc_deriv_method_id)
              CASE (xc_deriv_spline2_smooth)
                 CALL pw_nn_deriv_r ( pw_in=tmp_r(1)%pw,&
                      pw_out=v_xc(1)%pw,coeffs=spline2_deriv_coeffs,&
                      idir=idir, error=error )
              CASE (xc_deriv_spline3_smooth)
                 CALL pw_nn_deriv_r ( pw_in=tmp_r(1)%pw,&
                      pw_out=v_xc(1)%pw,coeffs=spline3_deriv_coeffs,&
                      idir=idir, error=error )
              CASE (xc_deriv_nn10_smooth)
                 CALL pw_nn_deriv_r ( pw_in=tmp_r(1)%pw,&
                      pw_out=v_xc(1)%pw,coeffs=nn10_deriv_coeffs,&
                      idir=idir, error=error )
              CASE (xc_deriv_nn50_smooth)
                 CALL pw_nn_deriv_r ( pw_in=tmp_r(1)%pw,&
                      pw_out=v_xc(1)%pw,coeffs=nn50_deriv_coeffs,&
                      idir=idir, error=error )
              CASE default
                 CALL cp_unimplemented_error(fromWhere=routineP, &
                         message="XC_DERIV method not implemented for GPW!",&
                         error=error, error_level=cp_failure_level)
              END SELECT

           END DO

           IF (ASSOCIATED(pw_pool)) THEN
              CALL pw_pool_give_back_pw(pw_pool,tmp_r(1)%pw,error=error)
           ELSE
              DEALLOCATE(tmp_r(1)%pw%cr3d, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(tmp_r(1)%pw, stat=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF
        END IF
        IF (ASSOCIATED(pw_pool)) THEN
           CALL pw_pool_give_back_pw(pw_pool,v_drho%pw,error=error)
        ELSE
           DEALLOCATE(v_drho%pw%cr3d, stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DEALLOCATE(v_drho%pw, stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

     END IF

  END IF

  DEALLOCATE(tmp_r, tmp_r2, tmp_a, tmp_b, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  CALL timestop(handle)
END SUBROUTINE xc_calc_2nd_deriv

! *****************************************************************************
!> \brief Does the first prepartions for the calculation of the 2nd deriv
!>      The calculation must then be performed with xc_calc_2nd_deriv
!>      Calculates the second derivative of E_xc at rho in the direction
!>      rho1  (if you see the second derivative as bilinear form)
!>      partial_rho|_(rho=rho) partial_rho|_(rho=rho) E_xc drho(rho1)drho
!>      The other direction is still indetermined, thus it returns
!>      a potential (partial integration is performed to reduce it to
!>      function of rho, removing the dependence from its partial derivs)
!> \param deriv_set object containing the potentials
!> \param rho_set object that will contain the density at which the
!>        derivatives were calculated
!> \param rho_r the place where you evaluate the derivative
!> \param pw_pool the pool for the grids
!> \param xc_section which functional should be used and how to calculate it
!> \param cell the simulation cell (needed to scale the derivatives)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE xc_prep_2nd_deriv(deriv_set, &
       rho_set, rho_r, pw_pool, xc_section, cell, error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_prep_2nd_deriv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, lsd
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r_pw, tau

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(xc_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       nspins     = SIZE(rho_r)
       lsd = (nspins /= 1)
    END IF

    ALLOCATE(rho_r_pw(nspins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1, nspins
       rho_r_pw(ispin)%pw => rho_r(ispin)%pw
    END DO

    NULLIFY(rho_g, tau)
    CALL xc_rho_set_and_dset_create(rho_set,deriv_set,2,&
       rho_r_pw,rho_g,tau,xc_section,cell,pw_pool,&
       needs_basic_components=.TRUE.,error=error)

    DEALLOCATE(rho_r_pw, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    bo = rho_r(1)%pw%pw_grid%bounds_local

    CALL divide_by_norm_drho(deriv_set, rho_set, lsd, error)

    CALL timestop(handle)
  END SUBROUTINE xc_prep_2nd_deriv

! *****************************************************************************
  SUBROUTINE divide_by_norm_drho(deriv_set, rho_set, lsd, error)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    LOGICAL, INTENT(IN)                      :: lsd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'divide_by_norm_drho', &
      routineP = moduleN//':'//routineN

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    CHARACTER(len=MAX_LABEL_LENGTH), &
      DIMENSION(:), POINTER                  :: split_desc
    INTEGER                                  :: i, idesc, j, k, order
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: drho_cutoff
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(xc_derivative_type), POINTER        :: deriv_att

! check for unknown derivatives and divide by norm_drho where necessary

    failure = .FALSE.

    bo = rho_set%local_bounds
    CALL xc_rho_set_get(rho_set, drho_cutoff=drho_cutoff, error=error)

    pos =>  deriv_set%derivs
    DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
       CALL xc_derivative_get(deriv_att,order=order,&
                              desc=desc, split_desc=split_desc,&
                              error=error)
       IF (order==1 .OR. order==2) THEN
          DO idesc=1,SIZE(split_desc)
             SELECT CASE(split_desc(idesc))
             CASE("norm_drho")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              MAX(rho_set%norm_drho(i,j,k), drho_cutoff)
                      END DO
                   END DO
                END DO
             CASE("norm_drhoa")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              MAX(rho_set%norm_drhoa(i,j,k), drho_cutoff)
                      END DO
                   END DO
                END DO
             CASE("norm_drhob")
                !$omp parallel do private(i,j,k)
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         deriv_att%deriv_data(i,j,k) = deriv_att%deriv_data(i,j,k) / &
                              MAX(rho_set%norm_drhob(i,j,k), drho_cutoff)
                      END DO
                   END DO
                END DO
             CASE("rho")
                CALL cp_assert(.NOT.lsd,cp_failure_level,cp_assertion_failed,routineP,&
                     "rho not handled in lsd: '"//&
                     TRIM(desc)//"' in "//&
CPSourceFileRef,&
                     error,failure)
             CASE("rhoa","rhob")
             CASE default
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "unhandled derivative: '"//&
                     TRIM(split_desc(idesc))//"' in '"//&
                     TRIM(desc)//"' in "//&
CPSourceFileRef,&
                     error,failure)
             END SELECT
          END DO
       END IF
    END DO

  END SUBROUTINE divide_by_norm_drho

! *****************************************************************************
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

