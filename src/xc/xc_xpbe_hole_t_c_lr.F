!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the exchange energy for the pbe hole model in a truncated
!>        coulomb potential, considering the long range part only. Can be used
!>        as longrange correction to a truncated Hartree Fock calculation
!> \par History
!>      Manuel Guidon (01.2009)  : initial version
!> \author Manuel Guidon (01.2009)
! *****************************************************************************

MODULE xc_xpbe_hole_t_c_lr
  USE erf_fn,                          ONLY: erf
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: euler,&
                                             pi
  USE mathlib,                         ONLY: expint
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  PUBLIC :: xpbe_hole_t_c_lr_lda_eval, xpbe_hole_t_c_lr_lda_info,&
            xpbe_hole_t_c_lr_lda_calc_1, xpbe_hole_t_c_lr_lda_calc_2,&
            xpbe_hole_t_c_lr_lsd_eval, xpbe_hole_t_c_lr_lsd_info

  REAL(KIND=dp), PARAMETER :: a1 = 0.00979681_dp, &
                              a2 = 0.04108340_dp, &
                              a3 = 0.18744000_dp, &
                              a4 = 0.00120824_dp, &
                              a5 = 0.0347188_dp
  REAL(KIND=dp), PARAMETER :: A  = 1.0161144_dp, &
                              B  = -0.37170836_dp, &
                              C  = -0.077215461_dp, &
                              D  = 0.57786348_dp, &
                              E  = -0.051955731_dp, &
                              F2 = 0.47965830_dp, &
                              F1 = 6.4753871_dp, &
                              clda = -0.73855876638202240588423_dp
  REAL(KIND=dp), PARAMETER :: expcutoff = 700.0_dp
  REAL(KIND=dp), PARAMETER :: smax = 8.572844_dp, &
                              sconst = 18.79622316_dp, &
                              scutoff = 8.3_dp
  REAL(KIND=dp), PARAMETER :: gcutoff = 0.08_dp, &
                              g1 = -0.02628417880_dp/E, &
                              g2 = -0.07117647788_dp/E, &
                              g3 = 0.08534541323_dp/E, &
                              g4 = 0.0_dp
  REAL(KIND=dp), PARAMETER :: wcutoff = 14.0_dp

  REAL(KIND=dp), PARAMETER :: EPS_RHO = 0.00000001_dp

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_xpbe_hole_t_c_lr'

  CONTAINS


! *****************************************************************************
!> \brief returns various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv controls the number of derivatives
!> \par History
!>        01.2009 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lda_info ( reference, shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LDA}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xpbe_hole_t_c_lr_lda_info

! *****************************************************************************
!> \brief returns various information on the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv controls the number of derivatives
!> \par History
!>        01.2009 created [mguidon]
!> \author mguidon
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lsd_info ( reference, shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    IF ( PRESENT ( reference ) ) THEN
      reference = "{LSD version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "{LSD}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho_spin=.TRUE.
       needs%norm_drho_spin=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=1

  END SUBROUTINE xpbe_hole_t_c_lr_lsd_info

! *****************************************************************************
!> \brief evaluates the pbe-hole exchange in a truncated coulomb potential
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params input parameters (scaling,cutoff)
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lda_eval ( rho_set, deriv_set, order, params)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'xpbe_hole_t_c_lr_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, R, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrho, e_rho, &
                                                norm_drho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)

    CALL section_vals_val_get(params,"SCALE_X",r_val=sx)
    CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R)


    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)


    CALL xc_rho_set_get(rho_set,rho=rho,&
         norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPASSERT(stat==0)
    ELSE
      dummy=> rho
    END IF

    e_0 => dummy
    e_rho => dummy
    e_ndrho => dummy

    IF (order>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (order>=1.OR.order==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rho)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrho)
    END IF
    IF (order>1.OR.order<-1) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 2 not implemented", &
            error_level=cp_failure_level)
    END IF

    IF ( R == 0.0_dp ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Cutoff_Radius 0.0 not implemented", &
                     error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(npoints, order, rho, norm_drho, e_0, e_rho) &
    !$omp          shared(e_ndrho, epsilon_rho) &
    !$omp          shared(sx, r)

    CALL xpbe_hole_t_c_lr_lda_calc(npoints,order,rho=rho, norm_drho=norm_drho,&
         e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,&
         epsilon_rho=epsilon_rho,&
         sx=sx,R=R)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy,stat=stat)
      CPASSERT(stat==0)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xpbe_hole_t_c_lr_lda_eval

! *****************************************************************************
!> \brief intermediate level routine in order to decide which branch has to be
!>        taken
!> \param npoints number of points on the grid
!> \param order order of the derivatives
!> \param rho value of density on the grid
!> \param norm_drho value of gradient on the grid
!> \param e_0 derivatives of the energy on the grid
!> \param e_rho derivatives of the energy on the grid 
!> \param e_ndrho derivatives of the energy on the grid
!> \param epsilon_rho cutoffs
!> \param sx functional parameters
!> \param R functional parameters
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note For numerical reasons there are two different branches
!>
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lda_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho,&
                                       epsilon_rho, sx, R)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, norm_drho, e_0, e_rho, &
                                                e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R

    CHARACTER(len=*), PARAMETER :: routineN = 'xpbe_hole_t_c_lr_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: my_ndrho, my_rho
    REAL(KIND=dp)                            :: ss, ss2, sscale, t1, t2, t3, &
                                                t4, t6, t7, t8

    !$omp do

    DO ip =1,npoints
      my_rho = MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho =MAX(norm_drho(ip),EPSILON(0.0_dp)*1.e4_dp)
        ! *** Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = (smax*ss2-sconst)/(ss2*ss)
        END IF
        IF(ss*sscale>gcutoff) THEN
          CALL xpbe_hole_t_c_lr_lda_calc_1(e_0(ip), e_rho(ip), e_ndrho(ip),&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, order)
        ELSE
          CALL xpbe_hole_t_c_lr_lda_calc_2(e_0(ip), e_rho(ip), e_ndrho(ip),&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, order)
        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE xpbe_hole_t_c_lr_lda_calc

! *****************************************************************************
!> \brief evaluates the pbe-hole exchange in a truncated coulomb potential
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param params input parameters (scaling,cutoff)
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lsd_eval ( rho_set, deriv_set, order, params)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'xpbe_hole_t_c_lr_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_rho, R, sx
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_ndrhoa, &
                                                e_ndrhob, e_rhoa, e_rhob, &
                                                norm_drhoa, norm_drhob, rhoa, &
                                                rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(bo)

    CALL section_vals_val_get(params,"SCALE_X",r_val=sx)
    CALL section_vals_val_get(params,"CUTOFF_RADIUS",r_val=R)


    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)


    CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
         norm_drhoa=norm_drhoa,norm_drhob=norm_drhob,local_bounds=bo,rho_cutoff=epsilon_rho)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

    ! meaningful default for the arrays we don't need: let us make compiler
    ! and debugger happy...
    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPASSERT(stat==0)
    ELSE
      dummy=> rhoa
    END IF

    e_0 => dummy
    e_rhoa => dummy
    e_rhob => dummy
    e_ndrhoa => dummy
    e_ndrhob => dummy

    IF (order>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_0)
    END IF
    IF (order>=1.OR.order==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rhoa)
      deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_rhob)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrhoa)
      deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
            allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv,deriv_data=e_ndrhob)
    END IF
    IF (order>1.OR.order<-1) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 2 not implemented", &
            error_level=cp_failure_level)
    END IF

    IF ( R == 0.0_dp ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Cutoff_Radius 0.0 not implemented", &
                     error_level=cp_failure_level)
    END IF

    !$omp parallel default(none) &
    !$omp          shared(npoints, order, rhoa, norm_drhoa, e_0, e_rhoa) &
    !$omp          shared(epsilon_rho, sx, r) &
    !$omp          shared(rhob, norm_drhob, e_rhob, e_ndrhoa, e_ndrhob)

    CALL xpbe_hole_t_c_lr_lsd_calc(npoints,order,rho=rhoa, norm_drho=norm_drhoa,&
         e_0=e_0,e_rho=e_rhoa,e_ndrho=e_ndrhoa,&
         epsilon_rho=epsilon_rho,sx=sx,R=R)
    CALL xpbe_hole_t_c_lr_lsd_calc(npoints,order,rho=rhob, norm_drho=norm_drhob,&
         e_0=e_0,e_rho=e_rhob,e_ndrho=e_ndrhob,&
         epsilon_rho=epsilon_rho,sx=sx,R=R)

    !$omp end parallel

    IF (cp_debug) THEN
      DEALLOCATE(dummy,stat=stat)
      CPASSERT(stat==0)
    ELSE
      NULLIFY(dummy)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xpbe_hole_t_c_lr_lsd_eval

! *****************************************************************************
!> \brief intermediate level routine in order to decide which branch has to be
!>        taken
!> \param npoints number of points on the grid
!> \param order order of the derivatives
!> \param rho density on the grid
!> \param norm_drho gradient on the grid
!> \param e_0 derivatives of the energy on the grid
!> \param e_rho derivatives of the energy on the grid
!> \param e_ndrho derivatives of the energy on the grid
!> \param epsilon_rho cutoffs
!> \param sx functional parameters
!> \param R functional parameters
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note For numerical reasons there are two different branches. This code calls
!>       the lda version applying spin scaling relations
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lsd_calc(npoints, order, rho, norm_drho, e_0, e_rho, e_ndrho,&
                                       epsilon_rho,  sx, R)

    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(inout)                          :: rho, norm_drho, e_0, e_rho, &
                                                e_ndrho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, sx, R

    CHARACTER(len=*), PARAMETER :: routineN = 'xpbe_hole_t_c_lr_lsd_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(dp)                                 :: my_ndrho, my_rho
    REAL(KIND=dp)                            :: e_tmp, ss, ss2, sscale, t1, &
                                                t2, t3, t4, t6, t7, t8

    !$omp do

    DO ip =1,npoints
      my_rho = 2.0_dp*MAX(rho(ip),0.0_dp)
      IF(my_rho > epsilon_rho) THEN
        my_ndrho = 2.0_dp*MAX(norm_drho(ip),0.0_dp)

        ! *** Do some precalculation in order to catch the correct branch afterwards
        sscale = 1.0_dp
        t1 = pi ** 2
        t2 = t1 * my_rho
        t3 = t2 ** (0.1e1_dp / 0.3e1_dp)
        t4 = 0.1e1_dp / t3
        t6 = my_ndrho * t4
        t7 = 0.1e1_dp / my_rho
        t8 = t7 * sscale
        ss = 0.3466806371753173524216762e0_dp * t6 * t8
        IF( ss > scutoff) THEN
          ss2 = ss*ss
          sscale = (smax*ss2-sconst)/(ss2*ss)
        END IF
        e_tmp = 0.0_dp
        IF(ss*sscale>gcutoff) THEN
          CALL xpbe_hole_t_c_lr_lda_calc_1(e_tmp, e_rho(ip), e_ndrho(ip),&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, order)
        ELSE
          CALL xpbe_hole_t_c_lr_lda_calc_2(e_tmp, e_rho(ip), e_ndrho(ip),&
                                           my_rho,&
                                           my_ndrho, sscale, sx, R, order)
        END IF
        e_0(ip) = e_0(ip) + 0.5_dp * e_tmp

      END IF
    END DO

    !$omp end do

  END SUBROUTINE xpbe_hole_t_c_lr_lsd_calc

! *****************************************************************************
!> \brief low level routine that calculates the energy derivatives in one point
!> \param e_0 derivatives of the energy on the grid
!> \param e_rho derivatives of the energy on the grid
!> \param e_ndrho derivatives of the energy on the grid
!> \param rho value of density on the grid
!> \param ndrho value of gradient on the grid
!> \param sscale functional parameters
!> \param sx functional parameters
!> \param R functional parameters
!> \param order order of the derivatives
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lda_calc_1(e_0, e_rho, e_ndrho, &
                                         rho, ndrho, sscale, sx, R, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, sscale, sx, R
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: p, pndrho, prho, q, qndrho, qrho, sndrho, srho, ssval, &
      t1, t103, t104, t105, t106, t107, t108, t109, t11, t113, t115, t116, &
      t117, t118, t119, t12, t121, t123, t125, t126, t129, t13, t133, t136, &
      t14, t140, t142, t143, t144, t147, t150, t152, t155, t156, t159, t162, &
      t163, t164, t167, t17, t170, t173, t175, t176, t177, t178, t181, t183, &
      t187, t19, t190, t194, t195, t196, t199, t2, t202, t209, t21, t214, &
      t216, t217, t218, t22, t222, t223, t224, t227, t228, t229, t232, t234, &
      t235, t236, t237, t240, t244, t248, t25, t251, t255, t256, t257, t261, &
      t262, t265, t269, t27, t273, t276
    REAL(KIND=dp) :: t279, t280, t281, t29, t292, t3, t308, t31, t312, t314, &
      t318, t320, t322, t33, t331, t334, t339, t34, t341, t347, t349, t35, &
      t359, t362, t368, t37, t370, t373, t38, t382, t385, t39, t393, t4, t40, &
      t401, t402, t404, t405, t409, t41, t42, t420, t421, t424, t425, t426, &
      t429, t430, t431, t437, t44, t445, t446, t45, t451, t46, t473, t475, &
      t479, t48, t484, t497, t498, t5, t50, t503, t504, t51, t517, t52, t523, &
      t524, t527, t53, t533, t534, t54, t541, t56, t568, t57, t58, t581, &
      t590, t6, t60, t603, t604, t61, t611, t612, t64, t640, t65, t653, t667, &
      t675, t677, t69, t7, t71, t716
    REAL(KIND=dp) :: t717, t720, t723, t73, t739, t758, t762, t77, t8, t800, &
      t803, t808, t85, t86, t87, t88, t91, t92, t93, t94, t95, t98, t99

    IF (order >=0 ) THEN
      t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = ndrho * t2
      t4 = pi ** 2
      t5 = t4 * rho
      t6 = t5 ** (0.1e1_dp / 0.3e1_dp)
      t7 = 0.1e1_dp / t6
      t8 = 0.1e1_dp / rho
      ssval = t3 * t7 * t8 / 0.6e1_dp
      t11 = red(rho, ndrho)
      t12 = t11 ** 2
      t13 = sscale ** 2
      t14 = t12 * t13
      t17 = t12 ** 2
      t19 = t13 ** 2
      t21 = a1 * t12 * t13 + a2 * t17 * t19
      t22 = t14 * t21
      t25 = t17 * t11
      t27 = t19 * sscale
      t29 = t17 * t12
      t31 = t19 * t13
      t33 = 1 + a3 * t17 * t19 + a4 * t25 * t27 + a5 * t29 * t31
      t34 = 0.1e1_dp / t33
      t35 = R ** 2
      t37 = t6 ** 2
      t38 = t2 * t37
      t39 = t34 * t35 * t38
      q = t22 * t39
      t40 = t21 * t34
      t41 = 0.1e1_dp / A
      t42 = t40 * t41
      p = 0.9e1_dp / 0.4e1_dp * t14 * t42
      t44 = rho ** (0.1e1_dp / 0.3e1_dp)
      t45 = t44 * rho
      t46 = exei(P,Q)
      t48 = expint(1, q)
      t50 = t14 * t40
      t51 = D + t50
      t52 = t51 * t35
      t53 = t52 * t38
      t54 = expint(1, t53)
      t56 = t51 ** 2
      t57 = t56 * t51
      t58 = 0.1e1_dp / t57
      t60 = D * C
      t61 = D ** 2
      t64 = D * t21
      t65 = t34 * B
      t69 = SQRT(pi)
      t71 = F1 * t21
      t73 = t71 * t34 + F2
      t77 = C * (1 + t73 * t12 * t13)
      t85 = t69 * (15 * E + 6 * t77 * t51 + 4 * B * t56 + 8 * A * t57)
      t86 = SQRT(t51)
      t87 = t86 * t57
      t88 = 0.1e1_dp / t87
      t91 = SQRT(A)
      t92 = pi * t91
      t93 = EXP(p)
      t94 = t11 * sscale
      t95 = SQRT(t42)
      t98 = erf(0.3e1_dp / 0.2e1_dp * t94 * t95)
      t99 = 1 - t98
      t103 = 0.3e1_dp / 0.4e1_dp * pi + t85 * t88 / 0.16e2_dp - 0.3e1_dp / 0.4e1_dp * t92  &
      * t93 * t99
      t104 = 0.1e1_dp / t69
      t105 = t103 * t104
      t106 = 0.1e1_dp / t12
      t107 = 0.1e1_dp / t13
      t108 = t106 * t107
      t109 = t108 * t87
      t113 = t40 * C + REAL(2 * t64 * t65,KIND=dp) - 0.32e2_dp / 0.15e2_dp * t105 * t109  &
      + t60 * t73
      t115 = t17 * t19
      t116 = t21 ** 2
      t117 = t33 ** 2
      t118 = 0.1e1_dp / t117
      t119 = t116 * t118
      t121 = C * t73
      t123 = t119 * B + t40 * t121
      t125 = t35 * t2
      t126 = t61 * C
      t129 = E * t21
      t133 = D * t103 * t104
      t136 = t34 * C
      t140 = REAL(2 * t129 * t34,KIND=dp) - 0.32e2_dp / 0.15e2_dp * t133 * t109 + REAL(2  &
      * t64 * t136,KIND=dp) + t126 * t73
      t142 = t105 * t106
      t143 = t107 * t87
      t144 = t143 * t40
      t147 = t136 * t73
      t150 = C * t116
      t152 = -0.32e2_dp / 0.15e2_dp * t142 * t144 + REAL(2 * t64 * t147,KIND=dp) + t150  &
      * t118
      t155 = t29 * t31 * C
      t156 = t73 * t116
      t159 = t126 + 2 * D * E + t14 * t140 + t115 * t152 + t155 * t156 *  &
      t118
      t162 = t35 ** 2
      t163 = t162 * t1
      t164 = t6 * t5
      t167 = t61 * t103 * t104
      t170 = t34 * E
      t173 = -0.16e2_dp / 0.15e2_dp * t167 * t109 + REAL(2 * t64 * t170,KIND=dp)
      t175 = t34 * t103
      t176 = t64 * t175
      t177 = t104 * t106
      t178 = t177 * t143
      t181 = E * t116
      t183 = -0.32e2_dp / 0.15e2_dp * t176 * t178 + t181 * t118
      t187 = t104 * t87 * t119
      t190 = t61 * E + t14 * t173 + t115 * t183 - 0.16e2_dp / 0.15e2_dp * t115  &
      * t103 * t187
      t194 = 2 * E + t60 + t61 * B + t14 * t113 + t115 * t123 + t125 * t37  &
      * t159 + 3 * t163 * t164 * t190
      t195 = t58 * t194
      t196 = D * t35
      t199 = EXP(-t196 * t38 - q)
      t202 = -0.4e1_dp / 0.9e1_dp * A * t46 + 0.4e1_dp / 0.9e1_dp * A * t48 - 0.4e1_dp /  &
      0.9e1_dp * A * t54 - 0.4e1_dp / 0.9e1_dp * t195 * t199
      e_0 = e_0 + (t45 * t202 * Clda) * sx
    END IF
    IF( order>=1 .OR. order>=-1) THEN
      t209 = rho ** 2
      srho = -t3 / t164 * t8 * t4 / 0.18e2_dp - t3 * t7 / t209 / 0.6e1_dp

      t214 = t2 * t7
      sndrho = t214 * t8 / 0.6e1_dp
      t216 = t11 * t13
      t217 = t216 * t40
      t218 = dsdrho(rho, ndrho)
      t222 = 2 * t217 * t125 * t37 * t218
      t223 = a1 * t11
      t224 = t13 * t218
      t227 = t12 * t11
      t228 = a2 * t227
      t229 = t19 * t218
      t232 = 2 * t223 * t224 + 4 * t228 * t229
      t234 = t14 * t232 * t39
      t235 = t21 * t118
      t236 = t14 * t235
      t237 = a3 * t227
      t240 = a4 * t17
      t244 = a5 * t25
      t248 = 4 * t237 * t229 + 5 * t240 * t27 * t218 + 6 * t244 * t31 * t218
      t251 = t236 * t125 * t37 * t248
      t255 = 0.2e1_dp / 0.3e1_dp * t50 * t125 * t7 * t4
      qrho = t222 + t234 - t251 + t255
      t256 = t216 * t21
      t257 = t34 * t41
      t261 = t232 * t34
      t262 = t261 * t41
      t265 = t118 * t41
      prho = 0.9e1_dp / 0.2e1_dp * t256 * t257 * t218 + 0.9e1_dp / 0.4e1_dp * t14 * t262  &
      - 0.9e1_dp / 0.4e1_dp * t22 * t265 * t248
      t269 = dsdndrho(rho)
      t273 = t13 * t269
      t276 = t19 * t269
      t279 = 2 * t223 * t273 + 4 * t228 * t276
      t280 = t279 * t34
      t281 = t280 * t41
      t292 = 4 * t237 * t276 + 5 * t240 * t27 * t269 + 6 * t244 * t31 * t269
      pndrho = 0.9e1_dp / 0.2e1_dp * t256 * t257 * t269 + 0.9e1_dp / 0.4e1_dp * t14 *  &
      t281 - 0.9e1_dp / 0.4e1_dp * t22 * t265 * t292
      qndrho = 2 * t217 * t125 * t37 * t269 + t14 * t279 * t39 - t236 * t125  &
      * t37 * t292
      t308 = dexeirho(P,Q,Prho,Qrho)
      t312 = EXP(-q)
      t314 = t312 * t106 * t107
      t318 = 0.1e1_dp / t35
      t320 = 0.1e1_dp / t37
      t322 = 0.1e1_dp / t21 * t33 * t318 * t1 * t320
      t331 = 2 * t216 * t40 * t218 + t14 * t261 - t14 * t235 * t248
      t334 = t214 * t4
      t339 = EXP(-t53)
      t341 = 0.1e1_dp / t51
      t347 = t56 ** 2
      t349 = 0.1e1_dp / t347 * t194
      t359 = D * t232
      t362 = t118 * B
      t368 = t118 * t248
      t370 = F1 * t232 * t34 - t71 * t368
      t373 = t73 * t11
      t382 = B * t51
      t385 = A * t56
      t393 = 0.1e1_dp / t86 / t347
      t401 = t92 * t93
      t402 = SQRT(0.31415926535897932385e1_dp)
      t404 = EXP(-p)
      t405 = 0.1e1_dp / t402 * t404
      t409 = 0.1e1_dp / t95
      t420 = REAL(t69 * (6 * C * (t370 * t12 * t13 + 2 * t373 * t224) * t51  &
      + 6 * t77 * t331 + 8 * t382 * t331 + 24 * t385 * t331) * t88,KIND=dp) /  &
      0.16e2_dp - 0.7e1_dp / 0.32e2_dp * REAL(t85,KIND=dp) * REAL(t393,KIND=dp) *&
      REAL(t331,KIND=dp) - 0.3e1_dp  &
      / 0.4e1_dp * t92 * prho * t93 * t99 + 0.3e1_dp / 0.2e1_dp * t401 * t405  &
      * (0.3e1_dp / 0.2e1_dp * t218 * sscale * t95 + 0.3e1_dp / 0.4e1_dp * t94 * t409  &
      * (t262 - t235 * t41 * t248))
      t421 = t420 * t104
      t424 = 0.1e1_dp / t227
      t425 = t105 * t424
      t426 = t143 * t218
      t429 = t86 * t56
      t430 = t107 * t429
      t431 = t430 * t331
      t437 = t227 * t19
      t445 = 0.1e1_dp / t117 / t33
      t446 = t116 * t445
      t451 = t121 * t248
      t473 = t424 * t107
      t475 = t473 * t87 * t218
      t479 = t108 * t429 * t331
      t484 = t118 * C
      t497 = t105 * t473
      t498 = t87 * t21
      t503 = t105 * t108
      t504 = t429 * t21
      t517 = t64 * t118
      t523 = C * t21
      t524 = t118 * t232
      t527 = t445 * t248
      t533 = t25 * t31 * C
      t534 = t118 * t218
      t541 = t73 * t21
      t568 = t118 * E
      t581 = t64 * t118 * t103
      t590 = t104 * t424
      t603 = t437 * t105
      t604 = t87 * t116
      t611 = t115 * t105
      t612 = t429 * t116
      e_rho = e_rho + (0.4e1_dp / 0.3e1_dp * t44 * t202 * Clda + t45 * (-0.4e1_dp / 0.9e1_dp *  &
      A * t308 - 0.4e1_dp / 0.27e2_dp * A * qrho * t314 * t322 + 0.4e1_dp / 0.27e2_dp  &
      * A * (t331 * t35 * t38 + 0.2e1_dp / 0.3e1_dp * t52 * t334) * t339 * t341  &
      * t318 * t1 * t320 + 0.4e1_dp / 0.3e1_dp * t349 * t199 * t331 - 0.4e1_dp  &
      / 0.9e1_dp * t58 * (REAL(2 * t216 * t113 * t218,KIND=dp) + t14 * (t261 * C  &
      - t235 * C * REAL(t248,KIND=dp) + REAL(2 * t359 * t65,KIND=dp) - REAL(2 * t64 * t362  &
      * t248,KIND=dp) - 0.32e2_dp / 0.15e2_dp * t421 * t109 + 0.64e2_dp / 0.15e2_dp * t425  &
      * t426 - 0.112e3_dp / 0.15e2_dp * t142 * t431 + t60 * t370) + REAL(4 *  &
      t437 * t123 * t218,KIND=dp) + t115 * (0.2e1_dp * t235 * B * t232 - 0.2e1_dp * t446  &
      * B * REAL(t248,KIND=dp) + t261 * t121 - t235 * t451 + t40 * C * t370)  &
      + 0.2e1_dp / 0.3e1_dp * t125 * t7 * t159 * t4 + t125 * t37 * (REAL(2 *  &
      t216 * t140 * t218,KIND=dp) + t14 * (0.2e1_dp * E * t232 * t34 - REAL(2 * t129  &
      * t368,KIND=dp) - 0.32e2_dp / 0.15e2_dp * D * t420 * t104 * t109 + 0.64e2_dp / 0.15e2_dp  &
      * t133 * t475 - 0.112e3_dp / 0.15e2_dp * t133 * t479 + REAL(2 * t359  &
      * t136,KIND=dp) - REAL(2 * t64 * t484 * t248,KIND=dp) + t126 * t370) + REAL(4 *  &
      t437 * t152 * t218,KIND=dp) + t115 * (-0.32e2_dp / 0.15e2_dp * t421 * t106 * t144  &
      + 0.64e2_dp / 0.15e2_dp * t497 * t498 * t34 * REAL(t218,KIND=dp) - 0.112e3_dp / 0.15e2_dp  &
      * t503 * t504 * t34 * t331 - 0.32e2_dp / 0.15e2_dp * t142 * t143 *  &
      t261 + 0.32e2_dp / 0.15e2_dp * t503 * t498 * REAL(t368,KIND=dp) + REAL(2 * t359  &
      * t147,KIND=dp) - 0.2e1_dp * t517 * t451 + 0.2e1_dp * REAL(t64,KIND=dp) * REAL(t136,KIND=dp) *  &
      t370 + REAL(2 * t523 * t524,KIND=dp) - REAL(2 * t150 * t527,KIND=dp)) + REAL(6 * t533  &
      * t156 * t534,KIND=dp) + t155 * t370 * t116 * t118 + 0.2e1_dp * t155 * t541  &
      * REAL(t524,KIND=dp) - 0.2e1_dp * t155 * REAL(t156,KIND=dp) * REAL(t527,KIND=dp)) + 0.4e1_dp  &
      * t163 * t6 * t190 * t4 + 0.3e1_dp * t163 * t164 * (REAL(2 * t216 * t173  &
      * t218,KIND=dp) + t14 * (-0.16e2_dp / 0.15e2_dp * t61 * t420 * t104 * t109 +  &
      0.32e2_dp / 0.15e2_dp * t167 * t475 - 0.56e2_dp / 0.15e2_dp * t167 * t479 + REAL(2  &
      * t359 * t170,KIND=dp) - REAL(2 * t64 * t568 * t248,KIND=dp)) + REAL(4 * t437  &
      * t183 * t218,KIND=dp) + t115 * (-0.32e2_dp / 0.15e2_dp * REAL(t359,KIND=dp) * REAL(t175,KIND=dp)  &
      * REAL(t178,KIND=dp) + 0.32e2_dp / 0.15e2_dp * t581 * t177 * t143 * REAL(t248,KIND=dp)  &
      - 0.32e2_dp / 0.15e2_dp * REAL(t64,KIND=dp) * t34 * t420 * REAL(t178,KIND=dp) + 0.64e2_dp  &
      / 0.15e2_dp * t176 * t590 * t426 - 0.112e3_dp / 0.15e2_dp * t176 * t177 *  &
      t431 + REAL(2 * t129 * t524,KIND=dp) - REAL(2 * t181 * t527,KIND=dp)) - 0.64e2_dp /  &
      0.15e2_dp * REAL(t603,KIND=dp) * REAL(t604,KIND=dp) * REAL(t534,KIND=dp) - 0.16e2_dp / 0.15e2_dp *  &
      t115 * t420 * t187 - 0.56e2_dp / 0.15e2_dp * t611 * t612 * t118 * t331 -  &
      0.32e2_dp / 0.15e2_dp * t611 * t498 * REAL(t524,KIND=dp) + 0.32e2_dp / 0.15e2_dp * t611  &
      * REAL(t604,KIND=dp) * REAL(t527,KIND=dp))) * t199 - 0.4e1_dp / 0.9e1_dp * t195 * (-0.2e1_dp  &
      / 0.3e1_dp * t196 * t334 - t222 - t234 + t251 - t255) * t199) *  &
      Clda ) * sx
      t640 = dexeindrho(P,Q,Pndrho,Qndrho)
      t653 = 2 * t216 * t40 * t269 + t14 * t280 - t14 * t235 * t292
      t667 = D * t279
      t675 = t118 * t292
      t677 = F1 * t279 * t34 - t71 * t675
      t716 = REAL(t69 * (6 * C * (t677 * t12 * t13 + 2 * t373 * t273) * t51  &
      + 6 * t77 * t653 + 8 * t382 * t653 + 24 * t385 * t653) * t88,KIND=dp) /  &
      0.16e2_dp - 0.7e1_dp / 0.32e2_dp * REAL(t85,KIND=dp) * REAL(t393,KIND=dp) * &
      REAL(t653,KIND=dp) - 0.3e1_dp  &
      / 0.4e1_dp * t92 * pndrho * t93 * t99 + 0.3e1_dp / 0.2e1_dp * t401 * t405  &
      * (0.3e1_dp / 0.2e1_dp * t269 * sscale * t95 + 0.3e1_dp / 0.4e1_dp * t94 *  &
      t409 * (t281 - t235 * t41 * t292))
      t717 = t716 * t104
      t720 = t143 * t269
      t723 = t430 * t653
      t739 = t121 * t292
      t758 = t473 * t87 * t269
      t762 = t108 * t429 * t653
      t800 = t118 * t279
      t803 = t445 * t292
      t808 = t118 * t269
      e_ndrho = e_ndrho + (t45 * (-0.4e1_dp / 0.9e1_dp * A * t640 - 0.4e1_dp / 0.27e2_dp * A * qndrho  &
      * t314 * t322 + 0.4e1_dp / 0.9e1_dp * A * t653 * t339 * t341 + 0.4e1_dp  &
      / 0.3e1_dp * t349 * t199 * t653 - 0.4e1_dp / 0.9e1_dp * t58 * (REAL(2 * t216  &
      * t113 * t269,KIND=dp) + t14 * (t280 * C - t235 * C * REAL(t292,KIND=dp) + REAL(2  &
      * t667 * t65,KIND=dp) - REAL(2 * t64 * t362 * t292,KIND=dp) - 0.32e2_dp / 0.15e2_dp  &
      * t717 * t109 + 0.64e2_dp / 0.15e2_dp * t425 * t720 - 0.112e3_dp / 0.15e2_dp *  &
      t142 * t723 + t60 * t677) + REAL(4 * t437 * t123 * t269,KIND=dp) + t115 *  &
      (0.2e1_dp * t235 * B * t279 - 0.2e1_dp * t446 * B * REAL(t292,KIND=dp) + t280 *  &
      t121 - t235 * t739 + t40 * C * t677) + t125 * t37 * (REAL(2 * t216  &
      * t140 * t269,KIND=dp) + t14 * (0.2e1_dp * E * t279 * t34 - REAL(2 * t129 *  &
      t675,KIND=dp) - 0.32e2_dp / 0.15e2_dp * D * t716 * t104 * t109 + 0.64e2_dp / 0.15e2_dp  &
      * t133 * t758 - 0.112e3_dp / 0.15e2_dp * t133 * t762 + REAL(2 * t667 *  &
      t136,KIND=dp) - REAL(2 * t64 * t484 * t292,KIND=dp) + t126 * t677) + REAL(4 * t437  &
      * t152 * t269,KIND=dp) + t115 * (-0.32e2_dp / 0.15e2_dp * t717 * t106 * t144 +  &
      0.64e2_dp / 0.15e2_dp * t497 * t498 * t34 * REAL(t269,KIND=dp) - 0.112e3_dp / 0.15e2_dp  &
      * t503 * t504 * t34 * t653 - 0.32e2_dp / 0.15e2_dp * t142 * t143 * t280  &
      + 0.32e2_dp / 0.15e2_dp * t503 * t498 * REAL(t675,KIND=dp) + REAL(2 * t667 *  &
      t147,KIND=dp) - 0.2e1_dp * t517 * t739 + 0.2e1_dp * REAL(t64,KIND=dp) * REAL(t136,KIND=dp) * t677  &
      + REAL(2 * t523 * t800,KIND=dp) - REAL(2 * t150 * t803,KIND=dp)) + REAL(6 * t533  &
      * t156 * t808,KIND=dp) + t155 * t677 * t116 * t118 + 0.2e1_dp * t155 * t541  &
      * REAL(t800,KIND=dp) - 0.2e1_dp * t155 * REAL(t156,KIND=dp) * REAL(t803,KIND=dp)) + 0.3e1_dp * t163  &
      * t164 * (REAL(2 * t216 * t173 * t269,KIND=dp) + t14 * (-0.16e2_dp / 0.15e2_dp  &
      * t61 * t716 * t104 * t109 + 0.32e2_dp / 0.15e2_dp * t167 * t758 - 0.56e2_dp  &
      / 0.15e2_dp * t167 * t762 + REAL(2 * t667 * t170,KIND=dp) - REAL(2 * t64  &
      * t568 * t292,KIND=dp)) + REAL(4 * t437 * t183 * t269,KIND=dp) + t115 * (-0.32e2_dp  &
      / 0.15e2_dp * REAL(t667,KIND=dp) * REAL(t175,KIND=dp) * REAL(t178,KIND=dp) + 0.32e2_dp / 0.15e2_dp  &
      * t581 * t177 * t143 * REAL(t292,KIND=dp) - 0.32e2_dp / 0.15e2_dp * REAL(t64,KIND=dp) *  &
      t34 * t716 * REAL(t178,KIND=dp) + 0.64e2_dp / 0.15e2_dp * t176 * t590 * t720 - 0.112e3_dp  &
      / 0.15e2_dp * t176 * t177 * t723 + REAL(2 * t129 * t800,KIND=dp) - REAL(2  &
      * t181 * t803,KIND=dp)) - 0.64e2_dp / 0.15e2_dp * REAL(t603,KIND=dp) * REAL(t604,KIND=dp) *  &
      REAL(t808,KIND=dp) - 0.16e2_dp / 0.15e2_dp * t115 * t716 * t187 - 0.56e2_dp / 0.15e2_dp  &
      * t611 * t612 * t118 * t653 - 0.32e2_dp / 0.15e2_dp * t611 * t498 * REAL(t800,KIND=dp)  &
      + 0.32e2_dp / 0.15e2_dp * t611 * REAL(t604,KIND=dp) * REAL(t803,KIND=dp))) * t199  &
      + 0.4e1_dp / 0.9e1_dp * t195 * qndrho * t199) * Clda ) * sx
    END IF


  END SUBROUTINE xpbe_hole_t_c_lr_lda_calc_1

! *****************************************************************************
!> \brief low level routine that calculates the energy derivatives in one point
!> \param e_0 derivatives of the energy on the grid
!> \param e_rho derivatives of the energy on the grid
!> \param e_ndrho derivatives of the energy on the grid
!> \param rho value of density on the grid
!> \param ndrho value of gradient on the grid
!> \param sscale functional parameters
!> \param sx functional parameters
!> \param R functional parameters
!> \param order order of the derivatives
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE xpbe_hole_t_c_lr_lda_calc_2(e_0, e_rho, e_ndrho, &
                                         rho, ndrho, sscale, sx, R, order)
    REAL(KIND=dp), INTENT(INOUT)             :: e_0, e_rho, e_ndrho
    REAL(KIND=dp), INTENT(IN)                :: rho, ndrho, sscale, sx, R
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp) :: p, pndrho, prho, q, qndrho, qrho, sndrho, srho, ssval, &
      t1, t102, t106, t11, t110, t113, t115, t117, t118, t119, t12, t122, &
      t125, t126, t127, t128, t13, t130, t133, t135, t138, t14, t140, t142, &
      t143, t146, t150, t151, t152, t155, t158, t165, t17, t170, t172, t173, &
      t174, t178, t179, t180, t183, t184, t185, t188, t19, t190, t191, t192, &
      t193, t196, t197, t2, t200, t204, t207, t21, t211, t212, t213, t217, &
      t22, t221, t225, t229, t232, t235, t236, t242, t248, t25, t264, t268, &
      t27, t272, t276, t278, t280, t287, t289, t29, t292, t297, t299, t3, &
      t305, t307, t31, t317, t320, t324, t327
    REAL(KIND=dp) :: t33, t330, t333, t334, t338, t34, t340, t344, t35, t352, &
      t353, t358, t37, t38, t380, t39, t394, t398, t399, t4, t40, t401, t406, &
      t407, t408, t41, t415, t435, t44, t45, t454, t46, t461, t48, t485, &
      t496, t498, t5, t50, t51, t512, t52, t524, t525, t529, t53, t531, t54, &
      t545, t56, t579, t58, t581, t586, t6, t60, t61, t64, t65, t7, t74, t75, &
      t77, t79, t8, t81, t83, t84, t85, t86, t89, t91, t93, t94, t95, t97

    IF( order >= 0 ) THEN
      t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
      t2 = t1 ** 2
      t3 = ndrho * t2
      t4 = pi ** 2
      t5 = t4 * rho
      t6 = t5 ** (0.1e1_dp / 0.3e1_dp)
      t7 = 0.1e1_dp / t6
      t8 = 0.1e1_dp / rho
      ssval = t3 * t7 * t8 / 0.6e1_dp

      t11 = red(rho, ndrho)
      t12 = t11 ** 2
      t13 = sscale ** 2
      t14 = t12 * t13
      t17 = t12 ** 2
      t19 = t13 ** 2
      t21 = a1 * t12 * t13 + a2 * t17 * t19
      t22 = t14 * t21
      t25 = t17 * t11
      t27 = t19 * sscale
      t29 = t17 * t12
      t31 = t19 * t13
      t33 = 1 + a3 * t17 * t19 + a4 * t25 * t27 + a5 * t29 * t31
      t34 = 0.1e1_dp / t33
      t35 = R ** 2
      t37 = t6 ** 2
      t38 = t2 * t37
      t39 = t34 * t35 * t38
      q = t22 * t39
      t40 = t21 * t34
      t41 = 0.1e1_dp / A
      p = 0.9e1_dp / 0.4e1_dp * t14 * t40 * t41
      t44 = rho ** (0.1e1_dp / 0.3e1_dp)
      t45 = t44 * rho
      t46 = exei(P, Q)
      t48 = expint(1, q)
      t50 = t14 * t40
      t51 = D + t50
      t52 = t51 * t35
      t53 = t52 * t38
      t54 = expint(1, t53)
      t56 = t51 ** 2
      t58 = 0.1e1_dp / t56 / t51
      t60 = D * C
      t61 = D ** 2
      t64 = D * t21
      t65 = t34 * B
      t74 = g1 + g2 * t12 * t13 + g3 * t17 * t19 + g4 * t25 * t27
      t75 = E * t74
      t77 = F1 * t21
      t79 = t77 * t34 + F2
      t81 = t40 * C + 2 * t64 * t65 + 2 * t75 + t60 * t79
      t83 = t17 * t19
      t84 = t21 ** 2
      t85 = t33 ** 2
      t86 = 0.1e1_dp / t85
      t89 = C * t79
      t91 = t84 * t86 * B + t40 * t89
      t93 = t35 * t2
      t94 = t61 * C
      t95 = D * E
      t97 = E * t21
      t102 = t34 * C
      t106 = 2 * t97 * t34 + 2 * t95 * t74 + 2 * t64 * t102 + t94 * t79
      t110 = t102 * t79
      t113 = C * t84
      t115 = 2 * t75 * t40 + 2 * t64 * t110 + t113 * t86
      t117 = t29 * t31
      t118 = t117 * C
      t119 = t79 * t84
      t122 = t94 + 2 * t95 + t14 * t106 + t83 * t115 + t118 * t119 * t86
      t125 = t35 ** 2
      t126 = t125 * t1
      t127 = t6 * t5
      t128 = t61 * E
      t130 = t34 * E
      t133 = t128 * t74 + 2 * t64 * t130
      t135 = t130 * t74
      t138 = E * t84
      t140 = 2 * t64 * t135 + t138 * t86
      t142 = t117 * E
      t143 = t74 * t84
      t146 = t128 + t14 * t133 + t83 * t140 + t142 * t143 * t86
      t150 = 2 * E + t60 + t61 * B + t14 * t81 + t83 * t91 + t93 * t37 *  &
      t122 + 3 * t126 * t127 * t146
      t151 = t58 * t150
      t152 = D * t35
      t155 = EXP(-t152 * t38 - q)
      t158 = -0.4e1_dp / 0.9e1_dp * A * t46 + 0.4e1_dp / 0.9e1_dp * A * t48 - 0.4e1_dp /  &
      0.9e1_dp * A * t54 - 0.4e1_dp / 0.9e1_dp * t151 * t155
      e_0 = e_0  + (t45 * t158 * Clda)* sx
    END IF
    IF (order >=1 .OR. order ==-1) THEN
      t165 = rho ** 2
      srho = -t3 / t127 * t8 * t4 / 0.18e2_dp - t3 * t7 / t165 / 0.6e1_dp
      t170 = t2 * t7
      sndrho = t170 * t8 / 0.6e1_dp
      t172 = t11 * t13
      t173 = t172 * t40
      t174 = dsdrho(rho, ndrho)
      t178 = 2 * t173 * t93 * t37 * t174
      t179 = a1 * t11
      t180 = t13 * t174
      t183 = t12 * t11
      t184 = a2 * t183
      t185 = t19 * t174
      t188 = 2 * t179 * t180 + 4 * t184 * t185
      t190 = t14 * t188 * t39
      t191 = t21 * t86
      t192 = t14 * t191
      t193 = a3 * t183
      t196 = a4 * t17
      t197 = t27 * t174
      t200 = a5 * t25
      t204 = 4 * t193 * t185 + 5 * t196 * t197 + 6 * t200 * t31 * t174
      t207 = t192 * t93 * t37 * t204
      t211 = 0.2e1_dp / 0.3e1_dp * t50 * t93 * t7 * t4
      qrho = t178 + t190 - t207 + t211
      t212 = t172 * t21
      t213 = t34 * t41
      t217 = t188 * t34
      t221 = t86 * t41
      prho = 0.9e1_dp / 0.2e1_dp * t212 * t213 * t174 + 0.9e1_dp / 0.4e1_dp * t14 * t217  &
      * t41 - 0.9e1_dp / 0.4e1_dp * t22 * t221 * t204
      t225 = dsdndrho(rho)
      t229 = t13 * t225
      t232 = t19 * t225
      t235 = 2 * t179 * t229 + 4 * t184 * t232
      t236 = t235 * t34
      t242 = t27 * t225
      t248 = 4 * t193 * t232 + 5 * t196 * t242 + 6 * t200 * t31 * t225
      pndrho = 0.9e1_dp / 0.2e1_dp * t212 * t213 * t225 + 0.9e1_dp / 0.4e1_dp * t14 *  &
      t236 * t41 - 0.9e1_dp / 0.4e1_dp * t22 * t221 * t248
      qndrho = 2 * t173 * t93 * t37 * t225 + t14 * t235 * t39 - t192 * t93  &
      * t37 * t248
      t264 = dexeirho(P,Q,Prho,Qrho)
      t268 = EXP(-q)
      t272 = t268 / t12 / t13
      t276 = 0.1e1_dp / t35
      t278 = 0.1e1_dp / t37
      t280 = 0.1e1_dp / t21 * t33 * t276 * t1 * t278
      t287 = t191 * t204
      t289 = 2 * t172 * t40 * t174 + t14 * t217 - t14 * t287
      t292 = t170 * t4
      t297 = EXP(-t53)
      t299 = 0.1e1_dp / t51
      t305 = t56 ** 2
      t307 = 0.1e1_dp / t305 * t150
      t317 = D * t188
      t320 = t86 * B
      t324 = g2 * t11
      t327 = g3 * t183
      t330 = g4 * t17
      t333 = 2 * t324 * t180 + 4 * t327 * t185 + 5 * t330 * t197
      t334 = E * t333
      t338 = t86 * t204
      t340 = F1 * t188 * t34 - t77 * t338
      t344 = t183 * t19
      t352 = 0.1e1_dp / t85 / t33
      t353 = t84 * t352
      t358 = t89 * t204
      t380 = t86 * C
      t394 = t64 * t86
      t398 = C * t21
      t399 = t86 * t188
      t401 = t352 * t204
      t406 = t25 * t31
      t407 = t406 * C
      t408 = t86 * t174
      t415 = t79 * t21
      t435 = t86 * E
      t454 = t406 * E
      t461 = t74 * t21
      e_rho = e_rho + (0.4e1_dp / 0.3e1_dp * t44 * t158 * Clda + t45 * (-0.4e1_dp / 0.9e1_dp *  &
      A * t264 - 0.4e1_dp / 0.27e2_dp * A * qrho * t272 * t280 + 0.4e1_dp / 0.27e2_dp  &
      * A * (t289 * t35 * t38 + 0.2e1_dp / 0.3e1_dp * t52 * t292) * t297 * t299  &
      * t276 * t1 * t278 + 0.4e1_dp / 0.3e1_dp * t307 * t155 * t289 - 0.4e1_dp  &
      / 0.9e1_dp * t58 * (REAL(2 * t172 * t81 * t174,KIND=dp) + REAL(t14 * (t217  &
      * C - t191 * C * t204 + 2 * t317 * t65 - 2 * t64 * t320 * t204 + 2  &
      * t334 + t60 * t340),KIND=dp) + REAL(4 * t344 * t91 * t174,KIND=dp) + REAL(t83 *  &
      (2 * t191 * B * t188 - 2 * t353 * B * t204 + t217 * t89 - t191 * t358  &
      + t40 * C * t340),KIND=dp) + 0.2e1_dp / 0.3e1_dp * t93 * t7 * t122 * t4 + t93  &
      * t37 * REAL(2 * t172 * t106 * t174 + t14 * (2 * E * t188 * t34  &
      - 2 * t97 * t338 + 2 * t95 * t333 + 2 * t317 * t102 - 2 * t64 * t380  &
      * t204 + t94 * t340) + 4 * t344 * t115 * t174 + t83 * (2 * t334  &
      * t40 + 2 * t75 * t217 - 2 * t75 * t287 + 2 * t317 * t110 - 2 * t394  &
      * t358 + 2 * t64 * t102 * t340 + 2 * t398 * t399 - 2 * t113 *  &
      t401) + 6 * t407 * t119 * t408 + t118 * t340 * t84 * t86 + 2 * t118  &
      * t415 * t399 - 2 * t118 * t119 * t401,KIND=dp) + 0.4e1_dp * t126 * t6 * t146  &
      * t4 + 0.3e1_dp * t126 * t127 * REAL(2 * t172 * t133 * t174 + t14  &
      * (t128 * t333 + 2 * t317 * t130 - 2 * t64 * t435 * t204) + 4 * t344  &
      * t140 * t174 + t83 * (2 * t317 * t135 - 2 * t394 * t75 * t204  &
      + 2 * t64 * t130 * t333 + 2 * t97 * t399 - 2 * t138 * t401) + 6 *  &
      t454 * t143 * t408 + t142 * t333 * t84 * t86 + 2 * t142 * t461 * t399  &
      - 2 * t142 * t143 * t401,KIND=dp)) * t155 - 0.4e1_dp / 0.9e1_dp * t151 * (-0.2e1_dp  &
      / 0.3e1_dp * t152 * t292 - t178 - t190 + t207 - t211) * t155) *  &
      Clda ) * sx
      t485 = dexeindrho(P,Q,Pndrho,Qndrho)
      t496 = t191 * t248
      t498 = 2 * t172 * t40 * t225 + t14 * t236 - t14 * t496
      t512 = D * t235
      t524 = 2 * t324 * t229 + 4 * t327 * t232 + 5 * t330 * t242
      t525 = E * t524
      t529 = t86 * t248
      t531 = F1 * t235 * t34 - t77 * t529
      t545 = t89 * t248
      t579 = t86 * t235
      t581 = t352 * t248
      t586 = t86 * t225
      e_ndrho = e_ndrho + (t45 * (-0.4e1_dp / 0.9e1_dp * A * t485 - 0.4e1_dp / 0.27e2_dp * A * qndrho  &
      * t272 * t280 + 0.4e1_dp / 0.9e1_dp * A * t498 * t297 * t299 + 0.4e1_dp  &
      / 0.3e1_dp * t307 * t155 * t498 - 0.4e1_dp / 0.9e1_dp * t58 * REAL(2 * t172  &
      * t81 * t225 + t14 * (t236 * C - t191 * C * t248 + 2 * t512 * t65  &
      - 2 * t64 * t320 * t248 + 2 * t525 + t60 * t531) + 4 * t344 * t91  &
      * t225 + t83 * (2 * t191 * B * t235 - 2 * t353 * B * t248 + t236  &
      * t89 - t191 * t545 + t40 * C * t531) + t93 * t37 * (2 * t172 *  &
      t106 * t225 + t14 * (2 * E * t235 * t34 - 2 * t97 * t529 + 2 * t95  &
      * t524 + 2 * t512 * t102 - 2 * t64 * t380 * t248 + t94 * t531) +  &
      4 * t344 * t115 * t225 + t83 * (2 * t525 * t40 + 2 * t75 * t236 -  &
      2 * t75 * t496 + 2 * t512 * t110 - 2 * t394 * t545 + 2 * t64 * t102  &
      * t531 + 2 * t398 * t579 - 2 * t113 * t581) + 6 * t407 * t119 *  &
      t586 + t118 * t531 * t84 * t86 + 2 * t118 * t415 * t579 - 2 * t118  &
      * t119 * t581) + 3 * t126 * t127 * (2 * t172 * t133 * t225 + t14  &
      * (t128 * t524 + 2 * t512 * t130 - 2 * t64 * t435 * t248) + 4 * t344  &
      * t140 * t225 + t83 * (2 * t512 * t135 - 2 * t394 * t75 * t248  &
      + 2 * t64 * t130 * t524 + 2 * t97 * t579 - 2 * t138 * t581) + 6 *  &
      t454 * t143 * t586 + t142 * t524 * t84 * t86 + 2 * t142 * t461 * t579  &
      - 2 * t142 * t143 * t581),KIND=dp) * t155 + 0.4e1_dp / 0.9e1_dp * t151 * qndrho  &
      * t155) * Clda ) * sx
    END IF

  END SUBROUTINE xpbe_hole_t_c_lr_lda_calc_2

! *****************************************************************************
!> \brief These functions evaluate products exp(x)*Ei(x) and pi*exp(x)*erfc(sqrt(x)),
!>      as well as their derivatives with respect to various combinations of
!>      rho and norm_drho.
!> \param P = 9/4*s^2*H/A
!> \param Q = s^2*H*R^2*kf^2
!> \retval exei ...
!> \par History
!>      01.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>     - In order to avoid numerical instabilities, these routines use Taylor-
!>       expansions for the above core-products for large arguments.
!>     - red calculates the reduced gradient in a save way (i.e. using a fixed
!>       cutoff EPS_RHO)
! *****************************************************************************
  FUNCTION exei(P,Q)
    REAL(dp), INTENT(IN)                     :: P, Q
    REAL(dp)                                 :: exei

    REAL(dp)                                 :: Q2, Q3, Q4, tmp

    exei = 0.0_dp
    IF( P<expcutoff ) THEN
      !Use exact product
      IF( P+Q < 0.5_dp ) THEN
        tmp = -euler -LOG(P+Q)+P+Q
        tmp = tmp - 0.25_dp*(P+Q)**2+1.0_dp/18.0_dp*(P+Q)**3-1.0_dp/96.0_dp*(P+Q)**4
        tmp = tmp +1.0_dp/600.0_dp*(P+Q)**5
        exei = EXP(P)*tmp
      ELSE
        exei = EXP(P)*expint(1,Q+P)
      END IF
    ELSE
      !Use approximation
      tmp = 1.0_dp / P
      ! *** 1st order
      exei = tmp
      ! *** 2nd order
      tmp = tmp / P
      exei = exei - (Q + 1.0_dp) * tmp
      ! *** 3rd order
      tmp = tmp / P
      Q2 = Q * Q
      exei = exei + (2.0_dp * Q + Q2 + 2.0_dp) * tmp
      ! *** 4th order
      tmp = tmp / P
      Q3 = Q2 * Q
      exei = exei - (3.0_dp * Q2 + 6.0_dp * Q + Q3 + 6.0_dp) * tmp
      ! *** 5th order
      tmp = tmp / P
      Q4 = Q3 * Q
      exei = exei + (24.0_dp + 4.0_dp * Q3 + Q4 + 12.0_dp * Q2 + 24.0_dp * Q) * tmp

      ! *** scaling factor
      exei = EXP(-Q) * exei
    END IF
  END FUNCTION exei

! *****************************************************************************
!> \brief ...
!> \param P ...
!> \param Q ...
!> \param dPrho ...
!> \param dQrho ...
!> \retval dexeirho ...
! *****************************************************************************
  FUNCTION dexeirho(P,Q,dPrho,dQrho)
    REAL(dp), INTENT(IN)                     :: P, Q, dPrho, dQrho
    REAL(dp)                                 :: dexeirho

    dexeirho = dPrho * exei(P,Q) - (dPrho + dQrho)/(P+Q)*EXP(-Q)
  END FUNCTION dexeirho

! *****************************************************************************
!> \brief ...
!> \param P ...
!> \param Q ...
!> \param dPndrho ...
!> \param dQndrho ...
!> \retval dexeindrho ...
! *****************************************************************************
  FUNCTION dexeindrho(P,Q,dPndrho,dQndrho)
    REAL(dp), INTENT(IN)                     :: P, Q, dPndrho, dQndrho
    REAL(dp)                                 :: dexeindrho

    dexeindrho = dPndrho * exei(P,Q) - (dPndrho + dQndrho)/(P+Q)*EXP(-Q)
  END FUNCTION dexeindrho

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param ndrho ...
!> \retval red ...
! *****************************************************************************
  FUNCTION red(rho,ndrho)
    REAL(dp), INTENT(IN)                     :: rho, ndrho
    REAL(dp)                                 :: red

    red = 1.0_dp/6.0_dp*ndrho*3.0_dp**(2.0_dp/3.0_dp)
    red = red /(pi**(2.0_dp/3.0_dp))
    red = red *MAX(1.0_dp/rho**(4.0_dp/3.0_dp),EPS_RHO)
  END FUNCTION red

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param ndrho ...
!> \retval dsdrho ...
! *****************************************************************************
  FUNCTION dsdrho(rho,ndrho)
    REAL(dp), INTENT(IN)                     :: rho, ndrho
    REAL(dp)                                 :: dsdrho

    dsdrho = -2.0_dp/9.0_dp*ndrho*3.0**(2.0_dp/3.0_dp)
    dsdrho = dsdrho/(pi**(2.0_dp/3.0_dp))
    dsdrho = dsdrho*MAX(1.0_dp/rho**(7.0_dp/3.0_dp),EPS_RHO)
  END FUNCTION dsdrho

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \retval dsdndrho ...
! *****************************************************************************
  FUNCTION dsdndrho(rho)
    REAL(dp), INTENT(IN)                     :: rho
    REAL(dp)                                 :: dsdndrho

    dsdndrho = 1.0_dp/6.0_dp*3.0_dp**(2.0_dp/3.0_dp)
    dsdndrho = dsdndrho/(pi**(2.0_dp/3.0_dp))
    dsdndrho = dsdndrho*MAX(1.0_dp/rho**(4.0_dp/3.0_dp),EPS_RHO)
  END FUNCTION dsdndrho

END MODULE xc_xpbe_hole_t_c_lr

