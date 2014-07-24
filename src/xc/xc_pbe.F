!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the pbe correlation functional
!> \note
!>      This was generated with the help of a maple worksheet that you can
!>      find under doc/pbe.mw .
!>      I did not add 3. derivatives in the polazied (lsd) case because it
!>      would have added 2500 lines of code. If you need them the worksheet
!>      is already prepared for them, and by uncommenting a couple of lines
!>      you should be able to generate them.
!> \par History
!>      09.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE xc_pbe
  USE bibliography,                    ONLY: Perdew1996,&
                                             Perdew2008,&
                                             Zhang1998,&
                                             cite_reference
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_input_constants,              ONLY: xc_pbe_orig,&
                                             xc_pbe_rev,&
                                             xc_pbe_sol
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_pbe'
  REAL(kind=dp), PARAMETER, PRIVATE :: a=0.04918_dp, b=0.132_dp,&
       c=0.2533_dp,d=0.349_dp

  PUBLIC :: pbe_lda_info, pbe_lsd_info, pbe_lda_eval, pbe_lsd_eval
CONTAINS

! *****************************************************************************
!> \brief return various information on the functional
!> \param pbe_params section selecting the varius parameters for the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pbe_lda_info(pbe_params,reference,shortform, needs, max_deriv,&
       error)
    TYPE(section_vals_type), POINTER         :: pbe_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: param
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: sc, sx

    failure=.FALSE.
    CALL section_vals_val_get(pbe_params,"parametrization",i_val=param,error=error)
    CALL section_vals_val_get(pbe_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(pbe_params,"scale_c",r_val=sc,error=error)

    SELECT CASE(param)
    CASE(xc_pbe_orig)
       CALL cite_reference(Perdew1996)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "J.P.Perdew, K.Burke, M.Ernzerhof, "//&
                  "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)"//&
                  "{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "PBE Perdew-Burke-Ernzerhof xc functional (unpolarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "J.P.Perdew, K.Burke, M.Ernzerhof, ",&
               "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)",&
               sx,sc,"{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LDA)')")&
                  "PBE, Perdew-Burke-Ernzerhof xc functional (unpolarized)",sx,sc
          END IF
       END IF
    CASE(xc_pbe_rev)
       CALL cite_reference(Perdew1996)
       CALL cite_reference(Zhang1998)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
          reference = "revPBE, Yingkay Zhang and Weitao Yang,"//&
               " Phys. Rev. Letter, vol 80,n 4, p. 890, (1998)"//&
               "{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "revPBE, revised PBE xc functional (unpolarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "revPBE, Yingkay Zhang and Weitao Yang,",&
               " Phys. Rev. Letter, vol 80,n 4, p. 890, (1998)",&
               sx,sc,"{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LDA)')")&
                  "revPBE, revised PBE xc functional (unpolarized)",sx,sc
          END IF
       END IF
    CASE(xc_pbe_sol)
       CALL cite_reference(Perdew1996)
       CALL cite_reference(Perdew2008)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
          reference = "PBEsol, J.P. Perdew et al., "//&
               "Phys. Rev. Letter, vol 100,n 13, p. 136406, (2008) "//&
               "{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "PBEsol, PBE for solids and surfaces xc functional (unpolarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "PBEsol, J.P. Perdew et al., ",&
               "Phys. Rev. Letter, vol 100,n 13, p. 136406, (2008) ",&
               sx,sc,"{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LDA)')")&
                  "PBEsol, PBE for solids and surfaces xc functional (unpolarized)",sx,sc
          END IF
       END IF
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
       needs%rho_1_3=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE pbe_lda_info

! *****************************************************************************
!> \brief return various information on the functional
!> \param pbe_params section selecting the varius parameters for the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE pbe_lsd_info(pbe_params,reference,shortform, needs, max_deriv,&
     error)
    TYPE(section_vals_type), POINTER         :: pbe_params
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: param
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: sc, sx

    failure=.FALSE.
    CALL section_vals_val_get(pbe_params,"parametrization",i_val=param,error=error)
    CALL section_vals_val_get(pbe_params,"scale_x",r_val=sx,error=error)
    CALL section_vals_val_get(pbe_params,"scale_c",r_val=sc,error=error)

    SELECT CASE(param)
    CASE(xc_pbe_orig)
       CALL cite_reference(Perdew1996)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "J.P.Perdew, K.Burke, M.Ernzerhof, "//&
                  "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)"//&
                  "{spin polarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "PBE xc functional (spin polarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "J.P.Perdew, K.Burke, M.Ernzerhof, ",&
               "Phys. Rev. Letter, vol. 77, n 18, pp. 3865-3868, (1996)",&
               sx,sc,"{spin polarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LSD)')")&
                  "PBE xc functional (spin polarized)",sx,sc
          END IF
       END IF
    CASE(xc_pbe_rev)
       CALL cite_reference(Perdew1996)
       CALL cite_reference(Zhang1998)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference="revPBE, Yingkay Zhang and Weitao Yang,"//&
               " Phys. Rev. Letter, vol 80,n 4, p. 890, (1998)"//&
               "{spin polarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "revPBE, revised PBE xc functional (spin polarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "revPBE, Yingkay Zhang and Weitao Yang,",&
               " Phys. Rev. Letter, vol 80,n 4, p. 890, (1998)",&
               sx,sc,"{spin polarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LSD)')")&
                  "revPBE, revised PBE xc functional (spin polarized)",sx,sc
          END IF
       END IF
    CASE(xc_pbe_sol)
       CALL cite_reference(Perdew1996)
       CALL cite_reference(Perdew2008)
       IF (sx==1._dp .AND.sc==1._dp) THEN
          IF ( PRESENT ( reference ) ) THEN
             reference = "PBEsol, J.P. Perdew et al., "//&
               "Phys. Rev. Letter, vol 100,n 13, p. 136406, (2008) "//&
               "{spin polarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             shortform = "PBEsol, PBE for solids and surfaces xc functional (spin polarized)"
          END IF
       ELSE
          IF ( PRESENT ( reference ) ) THEN
          WRITE (reference,"(a,a,'sx=',f5.3,'sc=',f5.3,a)")&
               "PBEsol, J.P. Perdew et al., ",&
               "Phys. Rev. Letter, vol 100,n 13, p. 136406, (2008) ",&
               sx,sc,"{spin unpolarized}"
          END IF
          IF ( PRESENT ( shortform ) ) THEN
             WRITE (shortform,"(a,'sx=',f5.3,'sc=',f5.3,'(LSD)')")&
                  "PBEsol, PBE for solids and surfaces xc functional (spin polarized)",sx,sc
          END IF
       END IF
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
  IF (PRESENT(needs)) THEN
     needs%rho_spin=.TRUE.
     needs%norm_drho_spin=.TRUE.
     needs%norm_drho=.TRUE.
  END IF
  IF (PRESENT(max_deriv)) max_deriv=2

END SUBROUTINE pbe_lsd_info

! *****************************************************************************
!> \brief evaluates the pbe correlation functional for lda
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param pbe_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE pbe_lda_eval(rho_set,deriv_set,grad_deriv,pbe_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: pbe_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, param, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_norm_drho, &
                                                epsilon_rho, scale_ec, &
                                                scale_ex
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho, &
      e_ndrho_rho_rho, e_rho, e_rho_rho, e_rho_rho_rho, norm_drho, rho, &
      rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)

  failure=.FALSE.
  NULLIFY(bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,rho_1_3=rho_1_3,rho=rho,&
          norm_drho=norm_drho,local_bounds=bo,rho_cutoff=epsilon_rho,&
          drho_cutoff=epsilon_norm_drho,error=error)
     npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

     ! meaningful default for the arrays we don't need: let us make compiler
     ! and debugger happy...
     IF (cp_debug) THEN
        ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ELSE
        dummy=> rho
     END IF

     e_0 => dummy
     e_rho => dummy
     e_ndrho => dummy
     e_rho_rho => dummy
     e_ndrho_rho => dummy
     e_ndrho_ndrho => dummy
     e_rho_rho_rho => dummy
     e_ndrho_rho_rho => dummy
     e_ndrho_ndrho_rho => dummy
     e_ndrho_ndrho_ndrho => dummy

     IF (grad_deriv>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
             allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
     END IF
     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(rho)(rho)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_rho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(rho)",allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_rho,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,&
             error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)
     END IF
     IF (grad_deriv>3.OR.grad_deriv<-3) THEN
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="derivatives bigger than 3 not implemented", &
             error=error, error_level=cp_failure_level)
     END IF


 CALL section_vals_val_get(pbe_params,"scale_c",r_val=scale_ec,error=error)
 CALL section_vals_val_get(pbe_params,"scale_x",r_val=scale_ex,error=error)
 CALL section_vals_val_get(pbe_params,"parametrization",i_val=param,error=error)

!$omp parallel default(none), &
!$omp shared(rho,rho_1_3,norm_drho,e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho),&
!$omp shared(e_ndrho_ndrho,e_rho_rho_rho,e_ndrho_rho_rho,e_ndrho_ndrho_rho),&
!$omp shared(e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho),&
!$omp shared(epsilon_norm_drho,pbe_params),&
!$omp shared(error,param,scale_ec,scale_ex)
     CALL pbe_lda_calc(rho=rho, rho_1_3=rho_1_3, norm_drho=norm_drho,&
          e_0=e_0,e_rho=e_rho,e_ndrho=e_ndrho,e_rho_rho=e_rho_rho,&
          e_ndrho_rho=e_ndrho_rho, e_ndrho_ndrho=e_ndrho_ndrho, &
          e_rho_rho_rho=e_rho_rho_rho, e_ndrho_rho_rho=e_ndrho_rho_rho,&
          e_ndrho_ndrho_rho=e_ndrho_ndrho_rho,e_ndrho_ndrho_ndrho=e_ndrho_ndrho_ndrho,&
          grad_deriv=grad_deriv,&
          npoints=npoints,epsilon_rho=epsilon_rho,epsilon_norm_drho=epsilon_norm_drho,&
          param=param,scale_ec=scale_ec,scale_ex=scale_ex,error=error)
!$omp end parallel

     IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     ELSE
        NULLIFY(dummy)
     END IF
  END IF

  CALL timestop(handle)
END SUBROUTINE pbe_lda_eval

! *****************************************************************************
!> \brief evaluates the pbe correlation functional for lda
!> \param rho ...
!> \param rho_1_3 ...
!> \param norm_drho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param e_rho_rho ...
!> \param e_ndrho_rho ...
!> \param e_ndrho_ndrho ...
!> \param e_rho_rho_rho ...
!> \param e_ndrho_rho_rho ...
!> \param e_ndrho_ndrho_rho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param grad_deriv degree of the derivative that should be evalated,
!>        if positive all the derivatives up to the given degree are evaluated,
!>        if negative only the given degree is calculated
!> \param npoints ...
!> \param epsilon_rho ...
!> \param epsilon_norm_drho ...
!> \param ! ...
!> \param pbe_params parameters for the pbe functional
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param rho_set the density where you want to evaluate the functional UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
! *****************************************************************************
SUBROUTINE pbe_lda_calc(rho, rho_1_3, norm_drho,&
     e_0,e_rho,e_ndrho,e_rho_rho,e_ndrho_rho,&
     e_ndrho_ndrho, e_rho_rho_rho, e_ndrho_rho_rho, e_ndrho_ndrho_rho,&
     e_ndrho_ndrho_ndrho,grad_deriv,npoints,epsilon_rho, epsilon_norm_drho,&
!     pbe_params,error)
     param,scale_ec,scale_ex,error)
    INTEGER, INTENT(in)                      :: npoints, grad_deriv
    REAL(kind=dp), DIMENSION(1:npoints), INTENT(inout) :: &
      e_ndrho_ndrho_ndrho, e_ndrho_ndrho_rho, e_ndrho_rho_rho, e_rho_rho_rho, &
      e_ndrho_ndrho, e_ndrho_rho, e_rho_rho, e_ndrho, e_rho, e_0
    REAL(kind=dp), DIMENSION(1:npoints), &
      INTENT(in)                             :: norm_drho, rho_1_3, rho
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_norm_drho
    INTEGER, INTENT(in)                      :: param
    REAL(kind=dp), INTENT(in)                :: scale_ec, scale_ex
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lda_calc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: A, A1rho, A1rhorho, A2rho, A_1, alpha_1_1, Arho, &
      Arho1rho, Arhorho, Arhorhorho, beta, beta_1_1, beta_2_1, beta_3_1, &
      beta_4_1, e_c_u_0, e_c_u_01rho, e_c_u_01rhorho, e_c_u_02rho, &
      e_c_u_0rho, e_c_u_0rho1rho, e_c_u_0rhorho, e_c_u_0rhorhorho, &
      epsilon_cGGA, epsilon_cGGArho, epsilon_cGGArhorho, ex_ldarhorhorho, &
      ex_unif, ex_unif1rho, ex_unif1rhorho, ex_unif2rho, ex_unifrho, &
      ex_unifrho1rho, ex_unifrhorho, Fx, Fx1rho, Fx1rhorho, Fx2rho, &
      Fxnorm_drho, Fxnorm_drho1rho, Fxnorm_drhonorm_drho, Fxnorm_drhorho, &
      Fxrho, Fxrho1rho, Fxrhorho, gamma_var, Hnorm_drho, Hnorm_drhonorm_drho
    REAL(kind=dp) :: Hnorm_drhorho, k_f, k_f2rho, k_frho, k_frhorho, &
      k_frhorhorho, k_s, k_s1rho, k_s1rhorho, k_s2rho, k_srho, k_srho1rho, &
      k_srhorho, kappa, kf, kf2rho, kfrho, kfrhorho, kfrhorhorho, mu, &
      my_norm_drho, my_rho, p_1, p_2, rs, rs2rho, rsrho, rsrhorho, &
      rsrhorhorho, s, s1rho, s1rhorho, s2norm_drho, s2rho, snorm_drho, &
      snorm_drho1rho, snorm_drhorho, srho, srho1rho, srhorho, t, t1, t1001, &
      t1004, t1005, t1006, t1008, t101, t1011, t1012, t1013, t1014, t1016, &
      t1017, t1018, t1019, t1022, t1024, t1026, t1028, t103, t1030, t1031, &
      t1035, t1042, t1046, t1048, t105, t1050, t1052, t1054, t1055
    REAL(kind=dp) :: t1060, t1061, t1062, t1067, t108, t1093, t1097, t11, &
      t1103, t1104, t1106, t1109, t111, t1115, t1118, t1121, t1122, t1126, &
      t1129, t113, t114, t1148, t115, t1152, t1157, t1167, t1187, t119, &
      t1196, t1203, t1218, t1226, t123, t124, t1249, t125, t126, t1262, &
      t1263, t127, t1291, t1292, t1295, t13, t131, t1329, t1342, t135, t138, &
      t1380, t1385, t1386, t1387, t1389, t139, t1393, t14, t140, t142, t146, &
      t1468, t148, t1482, t149, t1492, t1493, t150, t1504, t1505, t151, &
      t1511, t1515, t1521, t1525, t1528, t153, t1532, t1545, t155, t1565, &
      t157, t1573, t158, t1584, t159, t1608, t1612, t162
    REAL(kind=dp) :: t1628, t1629, t163, t1633, t164, t1646, t1652, t1660, &
      t167, t1672, t1676, t168, t17, t170, t171, t1715, t1722, t175, t1758, &
      t1759, t176, t1761, t177, t178, t179, t1797, t182, t183, t1838, t1851, &
      t186, t187, t1878, t1889, t189, t19, t190, t1907, t191, t1922, t1927, &
      t1933, t1937, t195, t1952, t196, t1964, t1965, t1968, t1969, t197, &
      t1972, t198, t1990, t1rho, t1rhorho, t2, t20, t200, t202, t2020, t2024, &
      t2028, t2031, t204, t2041, t208, t21, t214, t217, t218, t22, t226, &
      t229, t230, t238, t239, t240, t241, t246, t252, t253, t255, t259, t26, &
      t260, t261, t262, t265, t266, t267, t27
    REAL(kind=dp) :: t271, t272, t273, t277, t278, t280, t281, t286, t290, &
      t291, t293, t294, t295, t296, t297, t299, t2norm_drho, t2rho, t3, t305, &
      t309, t310, t315, t317, t318, t321, t323, t324, t327, t329, t330, t331, &
      t336, t338, t339, t340, t341, t348, t349, t354, t357, t359, t360, t361, &
      t362, t369, t370, t371, t374, t377, t378, t381, t382, t384, t385, t387, &
      t388, t390, t392, t393, t397, t4, t400, t401, t404, t408, t409, t410, &
      t414, t417, t418, t423, t426, t427, t432, t435, t436, t438, t439, t440, &
      t448, t449, t451, t456, t457, t458, t459, t461, t462, t463, t465, t466, &
      t469, t470, t471, t472
    REAL(kind=dp) :: t476, t487, t491, t496, t498, t5, t500, t503, t506, &
      t507, t510, t513, t517, t521, t525, t529, t530, t535, t541, t548, t553, &
      t556, t557, t559, t562, t566, t581, t586, t590, t591, t594, t598, t6, &
      t605, t609, t611, t612, t614, t627, t645, t65, t654, t656, t661, t664, &
      t669, t670, t671, t673, t675, t685, t69, t693, t7, t70, t701, t702, &
      t71, t714, t717, t72, t720, t73, t74, t740, t743, t748, t75, t758, t77, &
      t776, t78, t8, t80, t801, t809, t81, t812, t821, t825, t83, t831, t84, &
      t85, t86, t868, t87, t877, t879, t88, t880, t885, t90, t91, t94, t940, &
      t942, t943, t945, t946, t948, t95
    REAL(kind=dp) :: t950, t951, t954, t967, t976, t98, t980, t982, t984, &
      t989, t99, t990, t994, t995, t998, t999, tnorm_drho, tnorm_drho1rho, &
      tnorm_drhorho, tnorm_drhorhorho, trho, trho1rho, trhorho, trhorhorho

!TYPE(section_vals_type), POINTER         :: pbe_params
!, param
! scale_ec, scale_ex, snorm_drho, snorm_drho1rho, snorm_drhorho, srho, &

  failure=.FALSE.

  ! Parametrization of PBE
  SELECT CASE(param)
  ! Original PBE
  CASE (xc_pbe_orig)
     kappa = 0.804e0_dp
     beta  = 0.66725e-1_dp
     mu    = beta * pi ** 2 / 0.3e1_dp
  ! Revised PBE (revPBE)
  CASE (xc_pbe_rev)
     kappa = 0.1245e1_dp
     beta  = 0.66725e-1_dp
     mu    = beta * pi ** 2 / 0.3e1_dp
  ! PBE for solids (PBEsol)
  CASE (xc_pbe_sol)
     kappa = 0.804e0_dp
     beta  = 0.46e-1_dp
     mu    = 0.1e2_dp / 0.81e2_dp
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END SELECT

  gamma_var = (0.1e1_dp - LOG(0.2e1_dp)) / pi ** 2
  p_1 = 0.10e1_dp
  A_1 = 0.31091e-1_dp
  alpha_1_1 = 0.21370e0_dp
  beta_1_1 = 0.75957e1_dp
  beta_2_1 = 0.35876e1_dp
  beta_3_1 = 0.16382e1_dp
  beta_4_1 = 0.49294e0_dp
  p_2 = 0.10e1_dp
  t1 = 3 ** (0.1e1_dp / 0.3e1_dp)
  t2 = 4 ** (0.1e1_dp / 0.3e1_dp)
  t3 = t2 ** 2
  t4 = t1 * t3
  t5 = 0.1e1_dp / pi
  t69 = pi ** 2
  t627 = 0.1e1_dp / t69 / pi

  SELECT CASE(grad_deriv)
  CASE default
!$omp do
     DO ii=1,npoints
        my_rho = rho(ii)
        IF (my_rho>epsilon_rho) THEN
           my_norm_drho = norm_drho(ii)

           t6 = 0.1e1_dp / my_rho
           t7 = t5 * t6
           t8 = t7 ** (0.1e1_dp / 0.3e1_dp)
           rs = t4 * t8 / 0.4e1_dp
           t11 = 0.1e1_dp + alpha_1_1 * rs
           t13 = 0.1e1_dp / A_1
           t14 = SQRT(rs)
           t17 = t14 * rs
           t19 = p_1 + 0.1e1_dp
           t20 = rs ** t19
           t21 = beta_4_1 * t20
           t22 = beta_1_1 * t14 + beta_2_1 * rs + beta_3_1 * t17 + t21
           t26 = 0.1e1_dp + t13 / t22 / 0.2e1_dp
           t27 = LOG(t26)
           e_c_u_0 = -0.2e1_dp * A_1 * t11 * t27
           t65 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t70 = t69 * my_rho
           t71 = t70 ** (0.1e1_dp / 0.3e1_dp)
           k_f = t1 * t71
           t72 = k_f * t5
           t73 = SQRT(t72)
           k_s = 0.2e1_dp * t73
           t74 = 0.1e1_dp / k_s
           t75 = my_norm_drho * t74
           t = t75 * t6 / 0.2e1_dp
           t77 = 0.1e1_dp / gamma_var
           t78 = beta * t77
           t80 = EXP(-e_c_u_0 * t77)
           t81 = -0.1e1_dp + t80
           A = t78 / t81
           t83 = t ** 2
           t84 = A * t83
           t85 = 0.1e1_dp + t84
           t86 = t83 * t85
           t87 = A ** 2
           t88 = t83 ** 2
           t90 = 0.1e1_dp + t84 + t87 * t88
           t91 = 0.1e1_dp / t90
           t94 = 0.1e1_dp + t78 * t86 * t91
           t95 = LOG(t94)
           epsilon_cGGA = e_c_u_0 + gamma_var * t95
           kf = k_f
           ex_unif = -0.3e1_dp / 0.4e1_dp * t5 * kf
           t98 = 0.1e1_dp / kf
           t99 = my_norm_drho * t98
           s = t99 * t6 / 0.2e1_dp
           t101 = s ** 2
           t103 = 0.1e1_dp / kappa
           t105 = 0.1e1_dp + mu * t101 * t103
           Fx = 0.1e1_dp + kappa - kappa / t105
           t108 = my_rho * ex_unif

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)+&
                   scale_ex * t108 * Fx + scale_ec * my_rho * epsilon_cGGA
           END IF

           t111 = t8 ** 2
           t113 = 0.1e1_dp / t111 * t5
           t114 = my_rho ** 2
           t115 = 0.1e1_dp / t114
           rsrho = -t4 * t113 * t115 / 0.12e2_dp
           t119 = A_1 * alpha_1_1
           t123 = t22 ** 2
           t124 = 0.1e1_dp / t123
           t125 = t11 * t124
           t126 = 0.1e1_dp / t14
           t127 = beta_1_1 * t126
           t131 = beta_3_1 * t14
           t135 = 0.1e1_dp / rs
           t138 = t127 * rsrho / 0.2e1_dp + beta_2_1 * rsrho + 0.3e1_dp / &
                0.2e1_dp * t131 * rsrho + t21 * t19 * rsrho * t135
           t139 = 0.1e1_dp / t26
           t140 = t138 * t139
           e_c_u_0rho = -0.2e1_dp * t119 * rsrho * t27 + t125 * t140
           t142 = t71 ** 2
           k_frho = t1 / t142 * t69 / 0.3e1_dp
           t146 = 0.1e1_dp / t73
           k_srho = t146 * k_frho * t5
           t148 = k_s ** 2
           t149 = 0.1e1_dp / t148
           t150 = my_norm_drho * t149
           t151 = t6 * k_srho
           t153 = t75 * t115
           trho = -t150 * t151 / 0.2e1_dp - t153 / 0.2e1_dp
           t155 = gamma_var ** 2
           t157 = beta / t155
           t158 = t81 ** 2
           t159 = 0.1e1_dp / t158
           Arho = t157 * t159 * e_c_u_0rho * t80
           t162 = t78 * t
           t163 = t85 * t91
           t164 = t163 * trho
           t167 = Arho * t83
           t168 = A * t
           t170 = 0.2e1_dp * t168 * trho
           t171 = t167 + t170
           t175 = t78 * t83
           t176 = t90 ** 2
           t177 = 0.1e1_dp / t176
           t178 = t85 * t177
           t179 = A * t88
           t182 = t83 * t
           t183 = t87 * t182
           t186 = t167 + t170 + 0.2e1_dp * t179 * Arho + 0.4e1_dp * t183 * trho
           t187 = t178 * t186
           t189 = 0.2e1_dp * t162 * t164 + t78 * t83 * t171 * t91 - t175 * t187
           t190 = gamma_var * t189
           t191 = 0.1e1_dp / t94
           epsilon_cGGArho = e_c_u_0rho + t190 * t191
           kfrho = k_frho
           ex_unifrho = -0.3e1_dp / 0.4e1_dp * t5 * kfrho
           t195 = kf ** 2
           t196 = 0.1e1_dp / t195
           t197 = my_norm_drho * t196
           t198 = t6 * kfrho
           t200 = t99 * t115
           srho = -t197 * t198 / 0.2e1_dp - t200 / 0.2e1_dp
           t202 = t105 ** 2
           t204 = 0.1e1_dp / t202 * mu
           Fxrho = 0.2e1_dp * t204 * s * srho
           t208 = my_rho * ex_unifrho

           IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
              e_rho(ii) = e_rho(ii)+&
                   scale_ex * (ex_unif * Fx + t208 * Fx + t108 * Fxrho) + &
                   scale_ec * (epsilon_cGGA + my_rho * epsilon_cGGArho)
           END IF

           tnorm_drho = t74 * t6 / 0.2e1_dp
           t214 = t163 * tnorm_drho
           t217 = t78 * t182
           t218 = A * tnorm_drho
           t226 = 0.2e1_dp * t168 * tnorm_drho + 0.4e1_dp * t183 * tnorm_drho
           t229 = 0.2e1_dp * t162 * t214 + 0.2e1_dp * t217 * t218 * t91 - &
                t175 * t178 * t226
           t230 = gamma_var * t229
           Hnorm_drho = t230 * t191
           snorm_drho = t98 * t6 / 0.2e1_dp
           Fxnorm_drho = 0.2e1_dp * t204 * s * snorm_drho

           IF (grad_deriv>=1 .OR. grad_deriv==-1) THEN
              e_ndrho(ii) = e_ndrho(ii)+&
                   scale_ex * t108 * Fxnorm_drho + scale_ec * my_rho * &
                   Hnorm_drho
           END IF

           t238 = 0.1e1_dp / t69
           t239 = 0.1e1_dp / t111 / t7 * t238
           t240 = t114 ** 2
           t241 = 0.1e1_dp / t240
           t246 = 0.1e1_dp / t114 / my_rho
           rsrhorho = -t4 * t239 * t241 / 0.18e2_dp + t4 * t113 *&
                t246 / 0.6e1_dp
           t252 = 0.2e1_dp * t119 * rsrhorho * t27
           t253 = alpha_1_1 * rsrho
           t255 = t124 * t138 * t139
           t259 = 0.1e1_dp / t123 / t22
           t260 = t11 * t259
           t261 = t138 ** 2
           t262 = t261 * t139
           t265 = 0.1e1_dp / t17
           t266 = beta_1_1 * t265
           t267 = rsrho ** 2
           t271 = t127 * rsrhorho / 0.2e1_dp
           t272 = beta_2_1 * rsrhorho
           t273 = beta_3_1 * t126
           t277 = 0.3e1_dp / 0.2e1_dp * t131 * rsrhorho
           t278 = t19 ** 2
           t280 = rs ** 2
           t281 = 0.1e1_dp / t280
           t286 = t21 * t19 * rsrhorho * t135
           t290 = -t266 * t267 / 0.4e1_dp + t271 + t272 + 0.3e1_dp / 0.4e1_dp&
                * t273 * t267 + t277 + t21 * t278 * t267 * t281 + t286 - t21 * t19 &
                * t267 * t281
           t291 = t290 * t139
           t293 = t123 ** 2
           t294 = 0.1e1_dp / t293
           t295 = t11 * t294
           t296 = t26 ** 2
           t297 = 0.1e1_dp / t296
           t299 = t261 * t297 * t13
           e_c_u_0rhorho = -t252 + 0.2e1_dp * t253 * t255 - 0.2e1_dp * t260 *&
                t262 + t125 * t291 + t295 * t299 / 0.2e1_dp
           e_c_u_01rho = e_c_u_0rho
           t305 = t69 ** 2
           k_frhorho = -0.2e1_dp / 0.9e1_dp * t1 / t142 / t70 * t305
           t309 = 0.1e1_dp / t73 / t72
           t310 = k_frho ** 2
           t315 = t146 * k_frhorho * t5
           k_srhorho = -t309 * t310 * t238 / 0.2e1_dp + t315
           k_s1rho = k_srho
           t317 = 0.1e1_dp / t148 / k_s
           t318 = my_norm_drho * t317
           t321 = t115 * k_srho
           t323 = t150 * t321 / 0.2e1_dp
           t324 = t6 * k_srhorho
           t327 = t115 * k_s1rho
           t329 = t150 * t327 / 0.2e1_dp
           t330 = t75 * t246
           trhorho = t318 * t151 * k_s1rho + t323 - t150 * t324 / 0.2e1_dp + &
                t329 + t330
           t331 = t6 * k_s1rho
           t1rho = -t150 * t331 / 0.2e1_dp - t153 / 0.2e1_dp
           t336 = beta / t155 / gamma_var
           t338 = 0.1e1_dp / t158 / t81
           t339 = t336 * t338
           t340 = t80 ** 2
           t341 = e_c_u_0rho * t340
           t348 = t336 * t159
           t349 = e_c_u_0rho * e_c_u_01rho
           Arhorho = 0.2e1_dp * t339 * t341 * e_c_u_01rho + t157 * t159 * &
                e_c_u_0rhorho * t80 - t348 * t349 * t80
           A1rho = t157 * t159 * e_c_u_01rho * t80
           t354 = t78 * t1rho
           t357 = A1rho * t83
           t359 = 0.2e1_dp * t168 * t1rho
           t360 = t357 + t359
           t361 = t360 * t91
           t362 = t361 * trho
           t369 = t357 + t359 + 0.2e1_dp * t179 * A1rho + 0.4e1_dp * t183 * t1rho
           t370 = trho * t369
           t371 = t178 * t370
           t374 = t163 * trhorho
           t377 = t171 * t91
           t378 = t377 * t1rho
           t381 = Arhorho * t83
           t382 = Arho * t
           t384 = 0.2e1_dp * t382 * t1rho
           t385 = A1rho * t
           t387 = 0.2e1_dp * t385 * trho
           t388 = A * t1rho
           t390 = 0.2e1_dp * t388 * trho
           t392 = 0.2e1_dp * t168 * trhorho
           t393 = t381 + t384 + t387 + t390 + t392
           t397 = t171 * t177
           t400 = t186 * t1rho
           t401 = t178 * t400
           t404 = t360 * t177
           t408 = 0.1e1_dp / t176 / t90
           t409 = t85 * t408
           t410 = t186 * t369
           t414 = A1rho * t88
           t417 = A * t182
           t418 = Arho * t1rho
           t423 = trho * A1rho
           t426 = t87 * t83
           t427 = trho * t1rho
           t432 = t381 + t384 + t387 + t390 + t392 + 0.2e1_dp * t414 * Arho +&
                0.8e1_dp * t417 * t418 + 0.2e1_dp * t179 * Arhorho + 0.8e1_dp * &
                t417 * t423 + 0.12e2_dp * t426 * t427 + 0.4e1_dp * t183 * trhorho
           t435 = 0.2e1_dp * t354 * t164 + 0.2e1_dp * t162 * t362 - 0.2e1_dp &
                * t162 * t371 + 0.2e1_dp * t162 * t374 + 0.2e1_dp * t162 * t378 + &
                t78 * t83 * t393 * t91 - t175 * t397 * t369 - 0.2e1_dp * t162 * t401&
                - t175 * t404 * t186 + 0.2e1_dp * t175 * t409 * t410 - t175 * t178 &
                * t432
           t436 = gamma_var * t435
           t438 = t94 ** 2
           t439 = 0.1e1_dp / t438
           t440 = t163 * t1rho
           t448 = 0.2e1_dp * t162 * t440 + t78 * t83 * t360 * t91 - t175 * &
                t178 * t369
           t449 = t439 * t448
           t451 = gamma_var * t448
           epsilon_cGGArhorho = e_c_u_0rhorho + t436 * t191 - t190 * t449
           kfrhorho = k_frhorho
           ex_unifrhorho = -0.3e1_dp / 0.4e1_dp * t5 * kfrhorho
           ex_unif1rho = ex_unifrho
           t456 = 0.1e1_dp / t195 / kf
           t457 = my_norm_drho * t456
           t458 = kfrho ** 2
           t459 = t6 * t458
           t461 = t115 * kfrho
           t462 = t197 * t461
           t463 = t6 * kfrhorho
           t465 = t197 * t463 / 0.2e1_dp
           t466 = t99 * t246
           srhorho = t457 * t459 + t462 - t465 + t466
           s1rho = srho
           t469 = mu ** 2
           t470 = 0.1e1_dp / t202 / t105 * t469
           t471 = t470 * t101
           t472 = srho * t103
           t476 = s1rho * srho
           Fxrhorho = -0.8e1_dp * t471 * t472 * s1rho + 0.2e1_dp * t204 * &
                t476 + 0.2e1_dp * t204 * s * srhorho
           Fx1rho = 0.2e1_dp * t204 * s * s1rho
           t487 = my_rho * ex_unifrhorho
           t491 = my_rho * ex_unif1rho

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_rho_rho(ii) = e_rho_rho(ii)+&
                   scale_ex * (ex_unif1rho * Fx + ex_unif * Fx1rho + &
                   ex_unifrho * Fx + t487 * Fx + t208 * Fx1rho + ex_unif * Fxrho + t491&
                   * Fxrho + t108 * Fxrhorho) + scale_ec * (e_c_u_01rho + t451 * t191 &
                   + epsilon_cGGArho + my_rho * epsilon_cGGArhorho)
           END IF

           t496 = t149 * t6
           t498 = t74 * t115
           tnorm_drhorho = -t496 * k_srho / 0.2e1_dp - t498 / 0.2e1_dp
           t500 = t78 * trho
           t503 = t377 * tnorm_drho
           t506 = tnorm_drho * t186
           t507 = t178 * t506
           t510 = t163 * tnorm_drhorho
           t513 = t91 * trho
           t517 = Arho * tnorm_drho
           t521 = A * tnorm_drhorho
           t525 = t177 * t186
           t529 = t226 * trho
           t530 = t178 * t529
           t535 = t226 * t186
           t541 = A * trho
           t548 = tnorm_drho * trho
           t553 = 0.2e1_dp * t382 * tnorm_drho + 0.2e1_dp * t541 * tnorm_drho&
                + 0.2e1_dp * t168 * tnorm_drhorho + 0.8e1_dp * t417 * t517 + &
                0.12e2_dp * t426 * t548 + 0.4e1_dp * t183 * tnorm_drhorho
           t556 = 0.2e1_dp * t500 * t214 + 0.2e1_dp * t162 * t503 - 0.2e1_dp &
                * t162 * t507 + 0.2e1_dp * t162 * t510 + 0.6e1_dp * t175 * t218 * &
                t513 + 0.2e1_dp * t217 * t517 * t91 + 0.2e1_dp * t217 * t521 * t91 -&
                0.2e1_dp * t217 * t218 * t525 - 0.2e1_dp * t162 * t530 - t175 * &
                t397 * t226 + 0.2e1_dp * t175 * t409 * t535 - t175 * t178 * t553
           t557 = gamma_var * t556
           t559 = t439 * t189
           Hnorm_drhorho = t557 * t191 - t230 * t559
           t562 = t196 * t6
           snorm_drhorho = -t562 * kfrho / 0.2e1_dp - t98 * t115 / 0.2e1_dp
           t566 = snorm_drho * t103
           Fxnorm_drhorho = -0.8e1_dp * t471 * t566 * srho + 0.2e1_dp * t204 &
                * srho * snorm_drho + 0.2e1_dp * t204 * s * snorm_drhorho

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_ndrho_rho(ii) = e_ndrho_rho(ii)+&
                   scale_ex * (ex_unif * Fxnorm_drho + t208 * &
                   Fxnorm_drho + t108 * Fxnorm_drhorho) + scale_ec * (Hnorm_drho + my_rho &
                   * Hnorm_drhorho)
           END IF

           t581 = tnorm_drho ** 2
           t586 = A * t581
           t590 = tnorm_drho * t226
           t591 = t178 * t590
           t594 = t177 * t226
           t598 = t226 ** 2
           t605 = 0.2e1_dp * t586 + 0.12e2_dp * t426 * t581
           t609 = gamma_var * (0.2e1_dp * t78 * t581 * t85 * t91 + 0.10e2_dp &
                * t175 * t586 * t91 - 0.4e1_dp * t162 * t591 - 0.4e1_dp * t217 * &
                t218 * t594 + 0.2e1_dp * t175 * t409 * t598 - t175 * t178 * t605)
           t611 = t229 ** 2
           t612 = gamma_var * t611
           Hnorm_drhonorm_drho = t609 * t191 - t612 * t439
           t614 = snorm_drho ** 2
           Fxnorm_drhonorm_drho = -0.8e1_dp * t470 * t101 * t614 * t103 + &
                0.2e1_dp * t204 * t614

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii)+&
                   scale_ex * t108 * Fxnorm_drhonorm_drho +&
                   scale_ec * my_rho * Hnorm_drhonorm_drho
           END IF

           IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
              rsrhorhorho = -0.5e1_dp / 0.54e2_dp * t4 / t111 / t238 / &
                   t115 * t627 / t240 / t114 + t4 * t239 / t240 / my_rho / 0.3e1_dp &
                   - t4 * t113 * t241 / 0.2e1_dp
              rs2rho = rsrho
              t645 = alpha_1_1 * rsrhorho
              t654 = t127 * rs2rho / 0.2e1_dp + beta_2_1 * rs2rho + 0.3e1_dp / &
                   0.2e1_dp * t131 * rs2rho + t21 * t19 * rs2rho * t135
              t656 = t124 * t654 * t139
              t661 = t140 * t654
              t664 = rsrho * rs2rho
              t669 = t21 * t278
              t670 = rs2rho * t281
              t671 = t670 * rsrho
              t673 = t21 * t19
              t675 = -t266 * t664 / 0.4e1_dp + t271 + t272 + 0.3e1_dp / 0.4e1_dp&
                   * t273 * t664 + t277 + t669 * t671 + t286 - t673 * t671
              t685 = alpha_1_1 * rs2rho
              t693 = t675 * t139
              t701 = t297 * t13
              t702 = t701 * t654
              t714 = t267 * rs2rho
              t717 = rsrhorho * rsrho
              t720 = rsrhorho * rs2rho
              t740 = rs2rho / t280 / rs * t267
              t743 = rsrhorho * t281 * rsrho
              t748 = t670 * rsrhorho
              t758 = 0.3e1_dp / 0.8e1_dp * beta_1_1 / t14 / t280 * t714 - t266 *&
                   t717 / 0.2e1_dp - t266 * t720 / 0.4e1_dp + t127 * rsrhorhorho / &
                   0.2e1_dp + beta_2_1 * rsrhorhorho - 0.3e1_dp / 0.8e1_dp * beta_3_1 *&
                   t265 * t714 + 0.3e1_dp / 0.2e1_dp * t273 * t717 + 0.3e1_dp / &
                   0.4e1_dp * t273 * t720 + 0.3e1_dp / 0.2e1_dp * t131 * rsrhorhorho + &
                   t21 * t278 * t19 * t740 + 0.2e1_dp * t669 * t743 - 0.3e1_dp * t669 *&
                   t740 + t669 * t748 + t21 * t19 * rsrhorhorho * t135 - t673 * t748 -&
                   0.2e1_dp * t673 * t743 + 0.2e1_dp * t673 * t740
              t776 = A_1 ** 2
              e_c_u_0rhorhorho = -0.2e1_dp * t119 * rsrhorhorho * t27 + t645 * &
                   t656 + 0.2e1_dp * t645 * t255 - 0.4e1_dp * t253 * t259 * t661 + &
                   0.2e1_dp * t253 * t124 * t675 * t139 + t253 * t294 * t138 * t297 * &
                   t13 * t654 - 0.2e1_dp * t685 * t259 * t261 * t139 + 0.6e1_dp * t295 &
                   * t262 * t654 - 0.4e1_dp * t260 * t693 * t138 - 0.3e1_dp * t11 / &
                   t293 / t22 * t261 * t702 + t685 * t124 * t290 * t139 - 0.2e1_dp * &
                   t260 * t291 * t654 + t125 * t758 * t139 + t295 * t290 * t702 / &
                   0.2e1_dp + t685 * t294 * t299 / 0.2e1_dp + t295 * t675 * t701 * t138&
                   + t11 / t293 / t123 * t261 / t296 / t26 / t776 * t654 / 0.2e1_dp
              e_c_u_0rho1rho = -t252 + t253 * t656 + t685 * t255 - 0.2e1_dp * &
                   t260 * t661 + t125 * t693 + t295 * t138 * t702 / 0.2e1_dp
              e_c_u_01rhorho = e_c_u_0rho1rho
              e_c_u_02rho = -0.2e1_dp * t119 * rs2rho * t27 + t125 * t654 * t139
              k_frhorhorho = 0.10e2_dp / 0.27e2_dp * t1 / t142 / t114 * t69
              k_f2rho = kfrho
              t801 = k_f ** 2
              t809 = t309 * k_frhorho
              t812 = t238 * k_f2rho
              k_srho1rho = -t309 * k_frho * t812 / 0.2e1_dp + t315
              k_s1rhorho = k_srho1rho
              k_s2rho = t146 * k_f2rho * t5
              t821 = t148 ** 2
              t825 = k_srho * k_s1rho
              t831 = t6 * k_srho1rho
              trhorhorho = -0.3e1_dp * my_norm_drho / t821 * t6 * t825 * k_s2rho - &
                   t318 * t321 * k_s1rho + t318 * t831 * k_s1rho + t318 * t151 * &
                   k_s1rhorho - t318 * t321 * k_s2rho - t150 * t246 * k_srho + t150 * &
                   t115 * k_srho1rho / 0.2e1_dp + t318 * t324 * k_s2rho + t150 * t115 *&
                   k_srhorho / 0.2e1_dp - t150 * t6 * (0.3e1_dp / 0.4e1_dp / t73 / &
                   t801 / t238 * t310 * t627 * k_f2rho - t809 * t238 * k_frho - t809 * &
                   t812 / 0.2e1_dp + t146 * k_frhorhorho * t5) / 0.2e1_dp - t318 * t327&
                   * k_s2rho - t150 * t246 * k_s1rho + t150 * t115 * k_s1rhorho / &
                   0.2e1_dp - t150 * t246 * k_s2rho - 0.3e1_dp * t75 * t241
              t868 = t150 * t115 * k_s2rho / 0.2e1_dp
              trho1rho = t318 * t151 * k_s2rho + t323 - t150 * t831 / 0.2e1_dp +&
                   t868 + t330
              t1rhorho = t318 * t331 * k_s2rho + t329 - t150 * t6 * k_s1rhorho /&
                   0.2e1_dp + t868 + t330
              t2rho = -t150 * t6 * k_s2rho / 0.2e1_dp - t153 / 0.2e1_dp
              t877 = t155 ** 2
              t879 = beta / t877
              t880 = t158 ** 2
              t885 = e_c_u_01rho * e_c_u_02rho
              Arhorhorho = 0.6e1_dp * t879 / t880 * e_c_u_0rho * t340 * t80 * &
                   t885 + 0.2e1_dp * t339 * e_c_u_0rho1rho * t340 * e_c_u_01rho - &
                   0.6e1_dp * t879 * t338 * t341 * t885 + 0.2e1_dp * t339 * t341 * &
                   e_c_u_01rhorho + 0.2e1_dp * t339 * e_c_u_0rhorho * t340 * &
                   e_c_u_02rho + t157 * t159 * e_c_u_0rhorhorho * t80 - t348 * &
                   e_c_u_0rhorho * e_c_u_02rho * t80 - t348 * e_c_u_0rho1rho * &
                   e_c_u_01rho * t80 - t348 * e_c_u_0rho * e_c_u_01rhorho * t80 + t879 &
                   * t159 * t349 * e_c_u_02rho * t80
              Arho1rho = 0.2e1_dp * t339 * t341 * e_c_u_02rho + t157 * t159 * &
                   e_c_u_0rho1rho * t80 - t348 * e_c_u_0rho * e_c_u_02rho * t80
              A1rhorho = 0.2e1_dp * t339 * e_c_u_01rho * t340 * e_c_u_02rho + &
                   t157 * t159 * e_c_u_01rhorho * t80 - t348 * t885 * t80
              A2rho = t157 * t159 * e_c_u_02rho * t80
              t940 = Arho1rho * t83
              t942 = 0.2e1_dp * t382 * t2rho
              t943 = A2rho * t
              t945 = 0.2e1_dp * t943 * trho
              t946 = A * t2rho
              t948 = 0.2e1_dp * t946 * trho
              t950 = 0.2e1_dp * t168 * trho1rho
              t951 = A2rho * t88
              t954 = Arho * t2rho
              t967 = t940 + t942 + t945 + t948 + t950 + 0.2e1_dp * t951 * Arho +&
                   0.8e1_dp * t417 * t954 + 0.2e1_dp * t179 * Arho1rho + 0.8e1_dp * &
                   t417 * trho * A2rho + 0.12e2_dp * t426 * trho * t2rho + 0.4e1_dp * &
                   t183 * trho1rho
              t976 = t78 * t2rho
              t980 = t78 * t * t85
              t982 = A2rho * t83
              t984 = 0.2e1_dp * t168 * t2rho
              t989 = t982 + t984 + 0.2e1_dp * t179 * A2rho + 0.4e1_dp * t183 * t2rho
              t990 = t369 * t989
              t994 = t982 + t984
              t995 = t994 * t177
              t998 = Arhorhorho * t83
              t999 = Arhorho * t
              t1001 = 0.2e1_dp * t999 * t2rho
              t1004 = 0.2e1_dp * Arho1rho * t * t1rho
              t1005 = t954 * t1rho
              t1006 = 0.2e1_dp * t1005
              t1008 = 0.2e1_dp * t382 * t1rhorho
              t1011 = 0.2e1_dp * A1rhorho * t * trho
              t1012 = A1rho * t2rho
              t1013 = t1012 * trho
              t1014 = 0.2e1_dp * t1013
              t1016 = 0.2e1_dp * t385 * trho1rho
              t1017 = A2rho * t1rho
              t1018 = t1017 * trho
              t1019 = 0.2e1_dp * t1018
              t1022 = 0.2e1_dp * A * t1rhorho * trho
              t1024 = 0.2e1_dp * t388 * trho1rho
              t1026 = 0.2e1_dp * t943 * trhorho
              t1028 = 0.2e1_dp * t946 * trhorho
              t1030 = 0.2e1_dp * t168 * trhorhorho
              t1031 = t998 + t1001 + t1004 + t1006 + t1008 + t1011 + t1014 + &
                   t1016 + t1019 + t1022 + t1024 + t1026 + t1028 + t1030
              t1035 = t940 + t942 + t945 + t948 + t950
              t1042 = t369 * t2rho
              t1046 = A1rhorho * t83
              t1048 = 0.2e1_dp * t385 * t2rho
              t1050 = 0.2e1_dp * t943 * t1rho
              t1052 = 0.2e1_dp * t946 * t1rho
              t1054 = 0.2e1_dp * t168 * t1rhorho
              t1055 = t1046 + t1048 + t1050 + t1052 + t1054
              t1060 = t78 * t86
              t1061 = t176 ** 2
              t1062 = 0.1e1_dp / t1061
              t1067 = -0.2e1_dp * t162 * t178 * t967 * t1rho + 0.2e1_dp * t175 *&
                   t409 * t967 * t369 + 0.2e1_dp * t976 * t374 + 0.4e1_dp * t980 * &
                   t408 * trho * t990 - t175 * t995 * t432 + t78 * t83 * t1031 * t91 - &
                   t175 * t1035 * t177 * t369 - 0.2e1_dp * t162 * t995 * t370 - &
                   0.2e1_dp * t162 * t397 * t1042 + 0.2e1_dp * t162 * t1055 * t91 * &
                   trho - 0.6e1_dp * t1060 * t1062 * t186 * t990
              t1093 = t1046 + t1048 + t1050 + t1052 + t1054 + 0.2e1_dp * t951 * &
                   A1rho + 0.8e1_dp * t417 * t1012 + 0.2e1_dp * t179 * A1rhorho + &
                   0.8e1_dp * t417 * t1017 + 0.12e2_dp * t426 * t1rho * t2rho + &
                   0.4e1_dp * t183 * t1rhorho
              t1097 = t1rho * t989
              t1103 = trho * t989
              t1104 = t178 * t1103
              t1106 = t994 * t91
              t1109 = t171 * t408
              t1115 = t976 * t362 + t162 * t361 * trho1rho - t162 * t178 * t432 &
                   * t2rho + t175 * t994 * t408 * t410 - t162 * t178 * trho1rho * t369 &
                   - t162 * t178 * trho * t1093 - t162 * t397 * t1097 + t175 * t409 * &
                   t186 * t1093 - t354 * t1104 + t162 * t1106 * trhorho + t175 * t1109 &
                   * t990 + t175 * t409 * t432 * t989
              t1118 = t1106 * trho
              t1121 = t360 * t408
              t1122 = t186 * t989
              t1126 = t393 * t177
              t1129 = t408 * t186
              t1148 = t186 * t2rho
              t1152 = 0.2e1_dp * t354 * t1118 + 0.2e1_dp * t175 * t1121 * t1122 &
                   - t175 * t1126 * t989 + 0.4e1_dp * t980 * t1129 * t1042 + 0.2e1_dp *&
                   t976 * t378 - t175 * t404 * t967 + 0.2e1_dp * t162 * t377 * &
                   t1rhorho - 0.2e1_dp * t976 * t401 + 0.2e1_dp * t78 * t1rhorho * t164&
                   - 0.2e1_dp * t162 * t404 * t1103 - 0.2e1_dp * t162 * t404 * t1148
              t1157 = t393 * t91
              t1167 = A2rho * t182
              t1187 = A1rho * t182
              t1196 = t87 * t
              t1203 = t1008 + 0.8e1_dp * t417 * trhorho * A2rho + 0.12e2_dp * &
                   t426 * trhorho * t2rho + 0.8e1_dp * t1167 * t423 + 0.8e1_dp * t417 *&
                   trho * A1rhorho + 0.8e1_dp * t417 * Arho * t1rhorho + 0.8e1_dp * &
                   t1167 * t418 + 0.8e1_dp * t417 * Arho1rho * t1rho + 0.12e2_dp * t426&
                   * trho1rho * t1rho + 0.12e2_dp * t426 * trho * t1rhorho + t998 + &
                   0.8e1_dp * t1187 * t954 + 0.8e1_dp * t417 * trho1rho * A1rho + &
                   0.8e1_dp * t417 * Arhorho * t2rho + 0.24e2_dp * t1196 * t427 * t2rho&
                   + 0.2e1_dp * A1rhorho * t88 * Arho + t1014
              t1218 = t1016 + t1001 + 0.2e1_dp * t951 * Arhorho + t1026 + t1028 &
                   + t1022 + t1011 + t1030 + 0.2e1_dp * t179 * Arhorhorho + 0.4e1_dp * &
                   t183 * trhorhorho + 0.2e1_dp * t414 * Arho1rho + t1004 + t1006 + &
                   t1019 + t1024 + 0.24e2_dp * t84 * t1018 + 0.24e2_dp * t84 * t1005 + &
                   0.24e2_dp * t84 * t1013
              t1226 = t163 * trho1rho
              t1249 = -0.2e1_dp * t162 * t178 * t186 * t1rhorho + 0.2e1_dp * &
                   t162 * t1157 * t2rho - t175 * t178 * (t1203 + t1218) + 0.2e1_dp * &
                   t162 * t1035 * t91 * t1rho + 0.2e1_dp * t354 * t1226 - 0.2e1_dp * &
                   t162 * t995 * t400 + 0.2e1_dp * t162 * t163 * trhorhorho - 0.2e1_dp &
                   * t162 * t178 * trhorho * t989 - t175 * t1055 * t177 * t186 - &
                   0.2e1_dp * t976 * t371 - t175 * t397 * t1093 + 0.4e1_dp * t980 * &
                   t1129 * t1097
              t1262 = 0.2e1_dp * t162 * t163 * t2rho + t78 * t83 * t994 * t91 - &
                   t175 * t178 * t989
              t1263 = t439 * t1262
              t1291 = 0.2e1_dp * t976 * t164 + 0.2e1_dp * t162 * t1118 - &
                   0.2e1_dp * t162 * t1104 + 0.2e1_dp * t162 * t1226 + 0.2e1_dp * t162 &
                   * t377 * t2rho + t78 * t83 * t1035 * t91 - t175 * t397 * t989 - &
                   0.2e1_dp * t162 * t178 * t1148 - t175 * t995 * t186 + 0.2e1_dp * &
                   t175 * t409 * t1122 - t175 * t178 * t967
              t1292 = gamma_var * t1291
              t1295 = 0.1e1_dp / t438 / t94
              t1329 = 0.2e1_dp * t976 * t440 + 0.2e1_dp * t162 * t1106 * t1rho -&
                   0.2e1_dp * t162 * t178 * t1097 + 0.2e1_dp * t162 * t163 * t1rhorho &
                   + 0.2e1_dp * t162 * t361 * t2rho + t78 * t83 * t1055 * t91 - t175 * &
                   t404 * t989 - 0.2e1_dp * t162 * t178 * t1042 - t175 * t995 * t369 + &
                   0.2e1_dp * t175 * t409 * t990 - t175 * t178 * t1093
              kfrhorhorho = k_frhorhorho
              kf2rho = k_f2rho
              ex_unifrho1rho = ex_unifrhorho
              ex_unif1rhorho = ex_unifrho1rho
              ex_unif2rho = -0.3e1_dp / 0.4e1_dp * t5 * kf2rho
              t1342 = t195 ** 2
              srho1rho = t457 * t198 * kf2rho + t462 / 0.2e1_dp - t465 + t197 * &
                   t115 * kf2rho / 0.2e1_dp + t466
              s1rhorho = srho1rho
              s2rho = -t197 * t6 * kf2rho / 0.2e1_dp - t200 / 0.2e1_dp
              t1380 = t202 ** 2
              t1385 = 0.1e1_dp / t1380 * t469 * mu * t101 * s
              t1386 = kappa ** 2
              t1387 = 0.1e1_dp / t1386
              t1389 = s1rho * s2rho
              t1393 = t470 * s
              Fxrho1rho = -0.8e1_dp * t471 * t472 * s2rho + 0.2e1_dp * t204 * &
                   s2rho * srho + 0.2e1_dp * t204 * s * srho1rho
              Fx1rhorho = -0.8e1_dp * t471 * s1rho * t103 * s2rho + 0.2e1_dp * &
                   t204 * t1389 + 0.2e1_dp * t204 * s * s1rhorho
              Fx2rho = 0.2e1_dp * t204 * s * s2rho
              ex_ldarhorhorho = ex_unif1rhorho * Fx + ex_unif1rho * Fx2rho + &
                   ex_unif2rho * Fx1rho + ex_unif * Fx1rhorho + ex_unifrho1rho * Fx + &
                   ex_unifrho * Fx2rho + ex_unifrhorho * Fx - 0.3e1_dp / 0.4e1_dp * my_rho&
                   * t5 * kfrhorhorho * Fx + t487 * Fx2rho + ex_unifrho * Fx1rho + my_rho&
                   * ex_unifrho1rho * Fx1rho + t208 * Fx1rhorho + ex_unif2rho * Fxrho &
                   + ex_unif * Fxrho1rho + ex_unif1rho * Fxrho + my_rho * ex_unif1rhorho *&
                   Fxrho + t491 * Fxrho1rho + ex_unif * Fxrhorho + my_rho * ex_unif2rho *&
                   Fxrhorho + t108 * (0.48e2_dp * t1385 * srho * t1387 * t1389 - &
                   0.24e2_dp * t1393 * t472 * t1389 - 0.8e1_dp * t471 * srho1rho * t103&
                   * s1rho - 0.8e1_dp * t471 * t472 * s1rhorho + 0.2e1_dp * t204 * &
                   s1rhorho * srho + 0.2e1_dp * t204 * s1rho * srho1rho - 0.8e1_dp * &
                   t471 * srhorho * t103 * s2rho + 0.2e1_dp * t204 * s2rho * srhorho + &
                   0.2e1_dp * t204 * s * (-0.3e1_dp * my_norm_drho / t1342 * t459 * kf2rho&
                   - t457 * t115 * t458 + 0.2e1_dp * t457 * t463 * kfrho - 0.2e1_dp * &
                   t457 * t461 * kf2rho - 0.2e1_dp * t197 * t246 * kfrho + 0.3e1_dp / &
                   0.2e1_dp * t197 * t115 * kfrhorho + t457 * t463 * kf2rho - t197 * t6&
                   * kfrhorhorho / 0.2e1_dp - t197 * t246 * kf2rho - 0.3e1_dp * t99 * &
                   t241))

              e_rho_rho_rho(ii) = e_rho_rho_rho(ii)+&
                   scale_ex * ex_ldarhorhorho + scale_ec * (&
                   e_c_u_01rhorho + gamma_var * t1329 * t191 - t451 * t1263 + &
                   e_c_u_0rho1rho + t1292 * t191 - t190 * t1263 + epsilon_cGGArhorho + &
                   my_rho * (e_c_u_0rhorhorho + gamma_var * (t1067 + 0.2e1_dp * t1115 + &
                   t1152 + t1249) * t191 - t436 * t1263 - t1292 * t449 + 0.2e1_dp * &
                   t190 * t1295 * t448 * t1262 - t190 * t439 * t1329))

              t1468 = t149 * t115
              tnorm_drhorhorho = t317 * t6 * t825 + t1468 * k_srho / 0.2e1_dp - &
                   t496 * k_srhorho / 0.2e1_dp + t1468 * k_s1rho / 0.2e1_dp + t74 * &
                   t246
              tnorm_drho1rho = -t496 * k_s1rho / 0.2e1_dp - t498 / 0.2e1_dp
              t1482 = A1rho * tnorm_drhorho
              t1492 = t78 * t84
              t1493 = tnorm_drho * t177
              t1504 = t408 * tnorm_drho
              t1505 = t1504 * t410
              t1511 = A1rho * tnorm_drho
              t1515 = t226 * t1rho
              t1521 = A * tnorm_drho1rho
              t1525 = 0.2e1_dp * t175 * t409 * t553 * t369 + 0.2e1_dp * t217 * &
                   t1482 * t91 - 0.2e1_dp * t162 * t178 * tnorm_drhorho * t369 + &
                   0.2e1_dp * t354 * t510 - 0.6e1_dp * t1492 * t1493 * t400 + 0.6e1_dp &
                   * t175 * t218 * t91 * trhorho + 0.2e1_dp * t162 * t1157 * tnorm_drho&
                   + 0.4e1_dp * t980 * t1505 - 0.6e1_dp * t1492 * t1493 * t370 + &
                   0.6e1_dp * t175 * t1511 * t513 - 0.2e1_dp * t162 * t397 * t1515 - &
                   0.2e1_dp * t354 * t507 + 0.6e1_dp * t175 * t1521 * t513
              t1528 = Arhorho * tnorm_drho
              t1532 = t177 * t369
              t1545 = t78 * t417
              t1565 = 0.2e1_dp * t385 * tnorm_drho + 0.2e1_dp * t388 * &
                   tnorm_drho + 0.2e1_dp * t168 * tnorm_drho1rho + 0.8e1_dp * t417 * &
                   t1511 + 0.12e2_dp * t426 * tnorm_drho * t1rho + 0.4e1_dp * t183 * &
                   tnorm_drho1rho
              t1573 = t91 * t1rho
              t1584 = -0.2e1_dp * t354 * t530 + 0.2e1_dp * t217 * t1528 * t91 - &
                   0.2e1_dp * t217 * t517 * t1532 - 0.2e1_dp * t217 * t1511 * t525 + &
                   0.2e1_dp * t162 * t377 * tnorm_drho1rho + 0.2e1_dp * t162 * t361 * &
                   tnorm_drhorho + 0.4e1_dp * t1545 * t1505 + 0.2e1_dp * t217 * A * &
                   tnorm_drhorhorho * t91 + 0.2e1_dp * t175 * t409 * t1565 * t186 - &
                   0.2e1_dp * t217 * t521 * t1532 + 0.6e1_dp * t175 * t521 * t1573 - &
                   0.6e1_dp * t1060 * t1062 * t226 * t410 + 0.6e1_dp * t175 * t517 * &
                   t1573
              t1608 = t408 * t226
              t1612 = t163 * tnorm_drho1rho
              t1628 = -0.2e1_dp * t162 * t404 * t506 + 0.2e1_dp * t354 * t503 - &
                   0.2e1_dp * t162 * t178 * tnorm_drho * t432 - 0.2e1_dp * t162 * t178 &
                   * t553 * t1rho - 0.2e1_dp * t162 * t404 * t529 - 0.2e1_dp * t162 * &
                   t178 * t1565 * trho - t175 * t1126 * t226 + 0.4e1_dp * t980 * t1608 &
                   * t370 + 0.2e1_dp * t500 * t1612 + 0.12e2_dp * t78 * t168 * &
                   tnorm_drho * t91 * t427 + 0.2e1_dp * t162 * t163 * tnorm_drhorhorho &
                   - t175 * t404 * t553 + 0.2e1_dp * t175 * t1121 * t535
              t1629 = Arho * tnorm_drho1rho
              t1633 = tnorm_drho * t369
              t1646 = t361 * tnorm_drho
              t1652 = t226 * t369
              t1660 = t178 * t1633
              t1672 = t418 * tnorm_drho
              t1676 = t423 * tnorm_drho
              t1715 = 0.2e1_dp * t999 * tnorm_drho + 0.2e1_dp * t1672 + 0.2e1_dp&
                   * t382 * tnorm_drho1rho + 0.2e1_dp * t1676 + 0.2e1_dp * A * trhorho&
                   * tnorm_drho + 0.2e1_dp * t541 * tnorm_drho1rho + 0.2e1_dp * t385 *&
                   tnorm_drhorho + 0.2e1_dp * t388 * tnorm_drhorho + 0.2e1_dp * t168 *&
                   tnorm_drhorhorho + 0.8e1_dp * t1187 * t517 + 0.24e2_dp * t84 * &
                   t1672 + 0.8e1_dp * t417 * t1629 + 0.8e1_dp * t417 * t1528 + &
                   0.24e2_dp * t84 * t1676 + 0.24e2_dp * t1196 * t548 * t1rho + &
                   0.12e2_dp * t426 * tnorm_drho1rho * trho + 0.12e2_dp * t426 * &
                   tnorm_drho * trhorho + 0.8e1_dp * t417 * t1482 + 0.12e2_dp * t426 * &
                   tnorm_drhorho * t1rho + 0.4e1_dp * t183 * tnorm_drhorhorho
              t1722 = 0.2e1_dp * t217 * t1629 * t91 - 0.2e1_dp * t162 * t397 * &
                   t1633 - t175 * t397 * t1565 - 0.2e1_dp * t162 * t178 * t226 * &
                   trhorho + 0.2e1_dp * t78 * trhorho * t214 + 0.2e1_dp * t500 * t1646 &
                   - 0.2e1_dp * t217 * t1521 * t525 + 0.2e1_dp * t175 * t1109 * t1652 -&
                   0.2e1_dp * t162 * t178 * tnorm_drho1rho * t186 - 0.2e1_dp * t500 * &
                   t1660 + 0.4e1_dp * t980 * t1608 * t400 + 0.2e1_dp * t175 * t409 * &
                   t226 * t432 - t175 * t178 * t1715 - 0.2e1_dp * t217 * t218 * t177 * &
                   t432
              t1758 = 0.2e1_dp * t354 * t214 + 0.2e1_dp * t162 * t1646 - &
                   0.2e1_dp * t162 * t1660 + 0.2e1_dp * t162 * t1612 + 0.6e1_dp * t175 &
                   * t218 * t1573 + 0.2e1_dp * t217 * t1511 * t91 + 0.2e1_dp * t217 * &
                   t1521 * t91 - 0.2e1_dp * t217 * t218 * t1532 - 0.2e1_dp * t162 * &
                   t178 * t1515 - t175 * t404 * t226 + 0.2e1_dp * t175 * t409 * t1652 -&
                   t175 * t178 * t1565
              t1759 = gamma_var * t1758
              t1761 = t1295 * t189
              snorm_drho1rho = snorm_drhorho
              t1797 = snorm_drhorho * t103
              Fxnorm_drho1rho = -0.8e1_dp * t471 * t566 * s1rho + 0.2e1_dp * &
                   t204 * s1rho * snorm_drho + 0.2e1_dp * t204 * s * snorm_drho1rho

              e_ndrho_rho_rho(ii) = e_ndrho_rho_rho(ii)+&
                   scale_ex * (ex_unif1rho * Fxnorm_drho + &
                   ex_unif * Fxnorm_drho1rho + ex_unifrho * Fxnorm_drho + t487 * &
                   Fxnorm_drho + t208 * Fxnorm_drho1rho + ex_unif * Fxnorm_drhorho + &
                   t491 * Fxnorm_drhorho + t108 * (0.48e2_dp * t1385 * snorm_drho * &
                   t1387 * t476 - 0.24e2_dp * t1393 * t566 * t476 - 0.8e1_dp * t471 * &
                   snorm_drho1rho * t103 * srho - 0.8e1_dp * t471 * t566 * srhorho + &
                   0.2e1_dp * t204 * srhorho * snorm_drho + 0.2e1_dp * t204 * srho * &
                   snorm_drho1rho - 0.8e1_dp * t471 * t1797 * s1rho + 0.2e1_dp * t204 *&
                   s1rho * snorm_drhorho + 0.2e1_dp * t204 * s * (t456 * t6 * t458 + &
                   t196 * t115 * kfrho - t562 * kfrhorho / 0.2e1_dp + t98 * t246))) + &
                   scale_ec * (t1759 * t191 - t230 * t449 + Hnorm_drhorho + my_rho * (&
                   gamma_var * (t1525 + t1584 + t1628 + t1722) * t191 - t557 * t449 - &
                   t1759 * t559 + 0.2e1_dp * t230 * t1761 * t448 - t230 * t439 * t435))

              t1838 = t1504 * t535
              t1851 = Arho * t581
              t1878 = 0.4e1_dp * t175 * t409 * t226 * t553 - 0.2e1_dp * t162 * &
                   t178 * t605 * trho - 0.4e1_dp * t217 * t218 * t177 * t553 + 0.8e1_dp&
                   * t1545 * t1838 + 0.2e1_dp * t78 * t581 * t171 * t91 - 0.4e1_dp * &
                   t162 * t397 * t590 - 0.10e2_dp * t175 * t586 * t525 - t175 * t178 * &
                   (0.2e1_dp * t1851 + 0.4e1_dp * t521 * tnorm_drho + 0.24e2_dp * t84 *&
                   t1851 + 0.24e2_dp * t1196 * t581 * trho + 0.24e2_dp * t426 * &
                   tnorm_drhorho * tnorm_drho) - 0.12e2_dp * t1492 * t1493 * t529 + &
                   0.8e1_dp * t980 * t1838 - 0.4e1_dp * t162 * t178 * tnorm_drhorho * &
                   t226 + 0.20e2_dp * t162 * t586 * t513
              t1889 = t78 * t581
              t1907 = t85 * t1062
              t1922 = -0.4e1_dp * t500 * t591 + 0.2e1_dp * t175 * t409 * t605 * &
                   t186 + 0.4e1_dp * t162 * t409 * t598 * trho - 0.2e1_dp * t1889 * &
                   t187 + 0.2e1_dp * t175 * t1109 * t598 + 0.20e2_dp * t175 * t218 * &
                   t91 * tnorm_drhorho + 0.10e2_dp * t175 * t1851 * t91 - 0.4e1_dp * &
                   t217 * t517 * t594 - t175 * t397 * t605 - 0.6e1_dp * t175 * t1907 * &
                   t598 * t186 - 0.4e1_dp * t162 * t178 * tnorm_drho * t553 + 0.4e1_dp &
                   * t78 * tnorm_drhorho * t214 - 0.4e1_dp * t217 * t521 * t594
              t1927 = t439 * t229
              t1933 = t614 * t1387
              t1937 = t614 * t103

              e_ndrho_ndrho_rho(ii) = e_ndrho_ndrho_rho(ii)+&
                   scale_ex * (ex_unif * &
                   Fxnorm_drhonorm_drho + t208 * Fxnorm_drhonorm_drho + t108 * (&
                   0.48e2_dp * t1385 * t1933 * srho - 0.24e2_dp * t1393 * t1937 * srho &
                   - 0.16e2_dp * t471 * t1797 * snorm_drho + 0.4e1_dp * t204 * &
                   snorm_drhorho * snorm_drho)) + scale_ec * (Hnorm_drhonorm_drho + my_rho&
                   * (gamma_var * (t1878 + t1922) * t191 - t609 * t559 - 0.2e1_dp * &
                   t557 * t1927 + 0.2e1_dp * t612 * t1761))

              t2norm_drho = tnorm_drho
              t1952 = t226 * t2norm_drho
              t1964 = 0.2e1_dp * t168 * t2norm_drho + 0.4e1_dp * t183 * t2norm_drho
              t1965 = t178 * t1964
              t1968 = t226 * t1964
              t1969 = t1504 * t1968
              t1972 = A * t2norm_drho
              t1990 = 0.2e1_dp * t1972 * tnorm_drho + 0.12e2_dp * t426 * &
                   tnorm_drho * t2norm_drho
              t2020 = t177 * t1964
              t2024 = t91 * t2norm_drho
              t2028 = t78 * t2norm_drho
              t2031 = -0.20e2_dp * t1492 * t1493 * t1952 + 0.4e1_dp * t162 * &
                   t409 * t598 * t2norm_drho - 0.2e1_dp * t1889 * t1965 + 0.8e1_dp * &
                   t1545 * t1969 + 0.4e1_dp * t217 * t1972 * t408 * t598 - 0.6e1_dp * &
                   t175 * t1907 * t598 * t1964 - 0.2e1_dp * t217 * t1972 * t177 * t605 &
                   - 0.4e1_dp * t217 * t218 * t177 * t1990 + 0.4e1_dp * t175 * t409 * &
                   t1990 * t226 - 0.24e2_dp * t78 * t182 * t85 * t177 * t87 * t581 * &
                   t2norm_drho + 0.2e1_dp * t175 * t409 * t605 * t1964 + 0.8e1_dp * &
                   t980 * t1969 - 0.2e1_dp * t162 * t178 * t605 * t2norm_drho - &
                   0.4e1_dp * t162 * t178 * t1990 * tnorm_drho - 0.10e2_dp * t175 * &
                   t586 * t2020 + 0.24e2_dp * t162 * t586 * t2024 - 0.4e1_dp * t2028 * &
                   t591
              t2041 = 0.2e1_dp * t162 * t163 * t2norm_drho + 0.2e1_dp * t217 * &
                   t1972 * t91 - t175 * t1965
              s2norm_drho = snorm_drho

              e_ndrho_ndrho_ndrho(ii) = e_ndrho_ndrho_ndrho(ii)+&
                   scale_ex * t108 * (0.48e2_dp *&
                   t1385 * t1933 * s2norm_drho - 0.24e2_dp * t1393 * t1937 * &
                   s2norm_drho) + scale_ec * my_rho * (gamma_var * t2031 * t191 - t609 * &
                   t439 * t2041 - 0.2e1_dp * gamma_var * (0.2e1_dp * t2028 * t214 + &
                   0.10e2_dp * t175 * t218 * t2024 - 0.2e1_dp * t162 * t178 * &
                   tnorm_drho * t1964 - 0.2e1_dp * t217 * t218 * t2020 - 0.2e1_dp * &
                   t162 * t178 * t1952 - 0.2e1_dp * t217 * t1972 * t594 + 0.2e1_dp * &
                   t175 * t409 * t1968 - t175 * t178 * t1990) * t1927 + 0.2e1_dp * t612&
                   * t1295 * t2041)
           END IF
        END IF
     END DO
!$omp end do
  END SELECT

END SUBROUTINE pbe_lda_calc

! *****************************************************************************
!> \brief evaluates the becke 88 exchange functional for lsd
!> \param rho_set ...
!> \param deriv_set ...
!> \param grad_deriv ...
!> \param pbe_params ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE pbe_lsd_eval(rho_set,deriv_set,grad_deriv,pbe_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: grad_deriv
    TYPE(section_vals_type), POINTER         :: pbe_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, param, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: epsilon_drho, epsilon_rho, &
                                                scale_ec, scale_ex
    REAL(kind=dp), DIMENSION(:, :, :), POINTER :: dummy, e_0, e_ndr, &
      e_ndr_ndr, e_ndr_ra, e_ndr_rb, e_ndra, e_ndra_ndra, e_ndra_ra, e_ndrb, &
      e_ndrb_ndrb, e_ndrb_rb, e_ra, e_ra_ra, e_ra_rb, e_rb, e_rb_rb, &
      norm_drho, norm_drhoa, norm_drhob, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
  failure=.FALSE.
  NULLIFY(deriv, bo)

  CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL xc_rho_set_get(rho_set,&
          rhoa=rhoa, rhob=rhob, norm_drhoa=norm_drhoa, &
          norm_drhob=norm_drhob, norm_drho=norm_drho, &
          rho_cutoff=epsilon_rho,&
          drho_cutoff=epsilon_drho, local_bounds=bo, error=error)
     npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

     ! meaningful default for the arrays we don't need: let us make compiler
     ! and debugger happy...
     IF (cp_debug) THEN
        ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ELSE
        dummy=> rhoa
     END IF
     e_0 => dummy
     e_ra => dummy
     e_rb => dummy
     e_ndra_ra => dummy
     e_ndrb_rb => dummy
     e_ndr_ndr => dummy
     e_ndra_ndra => dummy
     e_ndrb_ndrb => dummy
     e_ndr => dummy
     e_ndra => dummy
     e_ndrb => dummy
     e_ra_ra => dummy
     e_ra_rb => dummy
     e_rb_rb => dummy
     e_ndr_ra => dummy
     e_ndr_rb => dummy

     IF (grad_deriv>=0) THEN
        deriv => xc_dset_get_derivative(deriv_set,"",&
             allocate_deriv=.TRUE., error=error)
        CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
     END IF
     IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb,error=error)
     END IF
     IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ra_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_rb_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)(rhoa)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra_ra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)(rhob)",&
             allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb_rb,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndr_ndr,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhoa)(norm_drhoa)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndra_ndra,error=error)
        deriv => xc_dset_get_derivative(deriv_set,&
             "(norm_drhob)(norm_drhob)", allocate_deriv=.TRUE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=e_ndrb_ndrb,error=error)
     END IF
     IF (grad_deriv>=3.OR.grad_deriv==-3) THEN
        ! to do
     END IF

 CALL section_vals_val_get(pbe_params,"scale_c",r_val=scale_ec,error=error)
 CALL section_vals_val_get(pbe_params,"scale_x",r_val=scale_ex,error=error)
 CALL section_vals_val_get(pbe_params,"parametrization",i_val=param,error=error)

!$omp parallel default(none),&
!$omp shared(rhoa,rhob,norm_drho,norm_drhoa,norm_drhob,e_0,e_ra,e_rb,e_ndra_ra),&
!$omp shared(e_ndrb_rb,e_ndr_ndr,e_ndra_ndra,e_ndrb_ndrb,e_ndr,e_ndra,e_ndrb),&
!$omp shared(e_ra_ra,e_ra_rb,e_rb_rb,e_ndr_ra,e_ndr_rb,grad_deriv,npoints),&
!$omp shared(epsilon_rho,epsilon_drho,error,param,scale_ec,scale_ex)
     CALL pbe_lsd_calc(&
          rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, norm_drhoa=norm_drhoa,&
          norm_drhob=norm_drhob, e_0=e_0, e_ra=e_ra, e_rb=e_rb,&
          e_ra_ndra=e_ndra_ra,&
          e_rb_ndrb=e_ndrb_rb, e_ndr_ndr=e_ndr_ndr,&
          e_ndra_ndra=e_ndra_ndra, e_ndrb_ndrb=e_ndrb_ndrb, e_ndr=e_ndr,&
          e_ndra=e_ndra, e_ndrb=e_ndrb, e_ra_ra=e_ra_ra, &
          e_ra_rb=e_ra_rb, e_rb_rb=e_rb_rb, e_ra_ndr=e_ndr_ra,&
          e_rb_ndr=e_ndr_rb,&
          grad_deriv=grad_deriv, npoints=npoints, &
          epsilon_rho=epsilon_rho,epsilon_drho=epsilon_drho,&
          param=param,scale_ec=scale_ec,scale_ex=scale_ex,error=error)
!$omp end parallel

     IF (cp_debug) THEN
        DEALLOCATE(dummy,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     ELSE
        NULLIFY(dummy)
     END IF
  END IF
  CALL timestop(handle)
END SUBROUTINE pbe_lsd_eval

! *****************************************************************************
!> \brief low level calculation of the pbe exchange-correlation functional for lsd
!> \param rhoa ,rhob: alpha or beta spin density
!> \param rhob ...
!> \param norm_drho ...
!> \param norm_drhoa ,norm_drhob,norm_drho: || grad rhoa |||,| grad rhoa ||,
!>        || grad (rhoa+rhob) ||
!> \param norm_drhob ...
!> \param e_0 ...
!> \param e_ra ...
!> \param e_rb ...
!> \param e_ra_ndra ...
!> \param e_rb_ndrb ...
!> \param e_ndr_ndr ...
!> \param e_ndra_ndra ...
!> \param e_ndrb_ndrb ...
!> \param e_ndr ...
!> \param e_ndra ...
!> \param e_ndrb ...
!> \param e_ra_ra ...
!> \param e_ra_rb ...
!> \param e_rb_rb ...
!> \param e_ra_ndr ...
!> \param e_rb_ndr ...
!> \param grad_deriv ...
!> \param npoints ...
!> \param epsilon_rho ...
!> \param epsilon_drho ...
!> \param param ...
!> \param scale_ec ...
!> \param scale_ex ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param e_ 0: adds to it the local value of the functional UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param e_ *: derivative of the functional wrt. to the variables
!>        named where the * is. UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author fawzi
! *****************************************************************************
SUBROUTINE pbe_lsd_calc(rhoa, rhob, norm_drho, norm_drhoa, norm_drhob,&
     e_0, e_ra, e_rb, e_ra_ndra, e_rb_ndrb, e_ndr_ndr,&
     e_ndra_ndra, e_ndrb_ndrb, e_ndr,&
     e_ndra, e_ndrb, e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr,&
     grad_deriv,npoints,epsilon_rho,epsilon_drho,param,scale_ec,scale_ex,error)
    REAL(kind=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob, norm_drho, &
                                                norm_drhoa, norm_drhob
    REAL(kind=dp), DIMENSION(*), INTENT(inout) :: e_0, e_ra, e_rb, e_ra_ndra, &
      e_rb_ndrb, e_ndr_ndr, e_ndra_ndra, e_ndrb_ndrb, e_ndr, e_ndra, e_ndrb, &
      e_ra_ra, e_ra_rb, e_rb_rb, e_ra_ndr, e_rb_ndr
    INTEGER, INTENT(in)                      :: grad_deriv, npoints
    REAL(kind=dp), INTENT(in)                :: epsilon_rho, epsilon_drho
    INTEGER, INTENT(in)                      :: param
    REAL(kind=dp), INTENT(in)                :: scale_ec, scale_ex
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pbe_lsd_calc', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER :: small = 1.0e-20, t_1_9 = 1.0_dp/9.0_dp, &
      t_2_3 = 2.0_dp/3.0_dp, t_4_27 = 4.0_dp/27.0_dp, t_4_3 = 4.0_dp/3.0_dp, &
      t_4_9 = 4.0_dp/9.0_dp, t_8_27 = 8.0_dp/27.0_dp

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(kind=dp) :: A, A1rhoa, A1rhob, A_1, A_2, A_3, alpha_1_1, alpha_1_2, &
      alpha_1_3, alpha_c, alpha_c1rhoa, alpha_c1rhob, alpha_crhoa, &
      alpha_crhob, Arhoa, Arhoarhoa, Arhoarhob, Arhob, Arhobrhob, beta, &
      beta_1_1, beta_1_2, beta_1_3, beta_2_1, beta_2_2, beta_2_3, beta_3_1, &
      beta_3_2, beta_3_3, beta_4_1, beta_4_2, beta_4_3, chi, chirhoa, &
      chirhoarhoa, chirhoarhob, chirhob, chirhobrhob, e_c_u_0, e_c_u_01rhoa, &
      e_c_u_01rhob, e_c_u_0rhoa, e_c_u_0rhoarhoa, e_c_u_0rhoarhob, &
      e_c_u_0rhob, e_c_u_0rhobrhob, e_c_u_1rhoa, e_c_u_1rhob, epsilon_c_unif, &
      epsilon_c_unif1rhoa, epsilon_c_unif1rhob, epsilon_c_unifrhoa
    REAL(kind=dp) :: epsilon_c_unifrhoarhoa, epsilon_c_unifrhoarhob, &
      epsilon_c_unifrhob, epsilon_c_unifrhobrhob, epsilon_cGGA, &
      epsilon_cGGArhoa, epsilon_cGGArhob, ex_unif_a, ex_unif_a1rhoa, &
      ex_unif_arhoa, ex_unif_b, ex_unif_b1rhob, ex_unif_brhob, f, f1rhoa, &
      f1rhob, frhoa, frhoarhoa, frhoarhob, frhob, frhobrhob, Fx_a, Fx_a1rhoa, &
      Fx_anorm_drhoa, Fx_arhoa, Fx_b, Fx_b1rhob, Fx_bnorm_drhob, Fx_brhob, &
      gamma_var, Hnorm_drho, k_frhoa, k_frhoarhoa, k_frhoarhob, k_frhob, k_s, &
      k_s1rhoa, k_s1rhob, k_srhoa, k_srhob, kappa, kf_a, kf_arhoa, &
      kf_arhoarhoa, kf_b, kf_brhob, kf_brhobrhob, mu, my_norm_drho
    REAL(kind=dp) :: my_norm_drhoa, my_norm_drhob, my_rho, my_rhoa, my_rhob, &
      p_1, p_2, p_3, phi, phi1rhoa, phi1rhob, phirhoa, phirhoarhoa, &
      phirhoarhob, phirhob, phirhobrhob, rs, rsrhoa, rsrhoarhoa, rsrhoarhob, &
      rsrhob, rsrhobrhob, s_a, s_a1rhoa, s_anorm_drhoa, s_arhoa, s_b, &
      s_b1rhob, s_bnorm_drhob, s_brhob, t, t1, t100, t1000, t1001, t101, &
      t102, t103, t1031, t1033, t1038, t104, t105, t1050, t1069, t107, t1071, &
      t108, t110, t1104, t1106, t111, t112, t113, t115, t116, t1164, t118, &
      t119, t1193, t12, t122, t1228, t123, t1231, t124, t125, t126, t1269, &
      t1283, t1285, t1286, t1288, t129, t1291, t1293
    REAL(kind=dp) :: t130, t1304, t131, t1327, t133, t1330, t1342, t1346, &
      t135, t137, t1377, t14, t140, t1411, t142, t143, t1440, t146, t1462, &
      t1465, t147, t148, t1482, t1489, t1492, t15, t150, t1501, t1514, t1520, &
      t1529, t153, t1550, t1552, t1555, t156, t1588, t1590, t1591, t1599, &
      t1602, t1603, t162, t163, t1630, t1632, t1635, t1638, t164, t1640, &
      t165, t167, t1674, t1677, t1680, t1698, t171, t1711, t1712, t1713, &
      t1739, t1741, t1743, t1748, t175, t176, t1765, t1766, t1767, t177, &
      t178, t179, t18, t1801, t1829, t183, t1830, t1831, t1865, t187, t1871, &
      t1876, t1888, t190, t1901, t191, t192, t1922
    REAL(kind=dp) :: t194, t1949, t198, t199, t1rhoa, t1rhob, t2, t20, t200, &
      t201, t205, t21, t211, t212, t213, t215, t219, t22, t220, t221, t222, &
      t226, t23, t232, t233, t234, t240, t242, t244, t245, t246, t248, t249, &
      t250, t252, t254, t256, t257, t259, t262, t266, t268, t269, t27, t272, &
      t273, t274, t277, t278, t28, t280, t281, t282, t284, t285, t286, t289, &
      t293, t294, t297, t298, t299, t3, t302, t303, t305, t306, t310, t311, &
      t312, t313, t314, t317, t318, t32, t321, t324, t325, t326, t329, t335, &
      t336, t337, t34, t340, t341, t344, t346, t350, t368, t38, t382, t39, &
      t396, t4, t40, t403, t405, t407
    REAL(kind=dp) :: t409, t41, t410, t411, t413, t415, t417, t427, t429, &
      t432, t436, t439, t442, t444, t445, t45, t453, t456, t457, t46, t460, &
      t466, t467, t468, t471, t472, t475, t477, t481, t488, t493, t494, t5, &
      t50, t502, t505, t506, t518, t519, t52, t523, t525, t536, t538, t543, &
      t544, t548, t549, t550, t556, t56, t561, t563, t564, t57, t576, t578, &
      t579, t58, t580, t588, t59, t590, t595, t596, t6, t600, t606, t611, &
      t624, t626, t627, t628, t63, t636, t638, t64, t643, t644, t648, t654, &
      t659, t66, t672, t674, t675, t676, t681, t682, t687, t69, t7, t70, &
      t705, t708, t71, t711, t72, t726, t73, t733
    REAL(kind=dp) :: t736, t74, t745, t75, t750, t755, t763, t767, t768, t77, &
      t775, t776, t779, t78, t785, t789, t79, t795, t798, t8, t80, t801, &
      t812, t82, t820, t821, t822, t823, t825, t828, t839, t84, t840, t85, &
      t851, t858, t865, t867, t868, t87, t876, t879, t88, t880, t9, t90, &
      t904, t908, t909, t91, t911, t914, t917, t919, t92, t924, t93, t936, &
      t94, t944, t95, t953, t959, t96, t962, t965, t966, t967, t97, t98, &
      t985, t998, t999, tnorm_drho, trhoa, trhoanorm_drho, trhoarhoa, &
      trhoarhob, trhob, trhobnorm_drho, trhobrhob

  failure=.FALSE.

  ! Parametrization of PBE
  SELECT CASE(param)
  ! Original PBE
  CASE (xc_pbe_orig)
     kappa = 0.804e0_dp
     beta  = 0.66725e-1_dp
     mu    = beta * pi ** 2 / 0.3e1_dp
  ! Revised PBE (revPBE)
  CASE (xc_pbe_rev)
     kappa = 0.1245e1_dp
     beta  = 0.66725e-1_dp
     mu    = beta * pi ** 2 / 0.3e1_dp
  ! PBE for solids (PBEsol)
  CASE (xc_pbe_sol)
     kappa = 0.804e0_dp
     beta  = 0.46e-1_dp
     mu    = 0.1e2_dp / 0.81e2_dp
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
  END SELECT

  gamma_var = (0.1e1_dp - LOG(0.2e1_dp)) / pi ** 2
  p_1 = 0.10e1_dp
  A_1 = 0.31091e-1_dp
  alpha_1_1 = 0.21370e0_dp
  beta_1_1 = 0.75957e1_dp
  beta_2_1 = 0.35876e1_dp
  beta_3_1 = 0.16382e1_dp
  beta_4_1 = 0.49294e0_dp
  p_2 = 0.10e1_dp
  A_2 = 0.15545e-1_dp
  alpha_1_2 = 0.20548e0_dp
  beta_1_2 = 0.141189e2_dp
  beta_2_2 = 0.61977e1_dp
  beta_3_2 = 0.33662e1_dp
  beta_4_2 = 0.62517e0_dp
  p_3 = 0.10e1_dp
  A_3 = 0.16887e-1_dp
  alpha_1_3 = 0.11125e0_dp
  beta_1_3 = 0.10357e2_dp
  beta_2_3 = 0.36231e1_dp
  beta_3_3 = 0.88026e0_dp
  beta_4_3 = 0.49671e0_dp
  t3 = 3 ** (0.1e1_dp / 0.3e1_dp)
  t4 = 4 ** (0.1e1_dp / 0.3e1_dp)
  t5 = t4 ** 2
  t6 = t3 * t5
  t7 = 0.1e1_dp / pi
  t90 = pi ** 2

  SELECT CASE(grad_deriv)
  CASE default
!$omp do
     DO ii=1,npoints
        my_rhoa=MAX(rhoa(ii),0.0_dp)
        my_rhob=MAX(rhob(ii),0.0_dp)
        my_rho=my_rhoa+my_rhob
        IF (my_rho>epsilon_rho) THEN
           my_rhoa=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhoa)
           my_rhob=MAX(EPSILON(0.0_dp)*1.e4_dp,my_rhob)
           my_rho = my_rhoa + my_rhob
           my_norm_drho=norm_drho(ii)
           my_norm_drhoa=norm_drhoa(ii)
           my_norm_drhob=norm_drhob(ii)

           t1 = my_rhoa - my_rhob
           t2 = 0.1e1_dp / my_rho
           chi = t1 * t2
           t8 = t7 * t2
           t9 = t8 ** (0.1e1_dp / 0.3e1_dp)
           rs = t6 * t9 / 0.4e1_dp
           t12 = 0.1e1_dp + alpha_1_1 * rs
           t14 = 0.1e1_dp / A_1
           t15 = SQRT(rs)
           t18 = t15 * rs
           t20 = p_1 + 0.1e1_dp
           t21 = rs ** t20
           t22 = beta_4_1 * t21
           t23 = beta_1_1 * t15 + beta_2_1 * rs + beta_3_1 * t18 + t22
           t27 = 0.1e1_dp + t14 / t23 / 0.2e1_dp
           t28 = LOG(t27)
           e_c_u_0 = -0.2e1_dp * A_1 * t12 * t28
           t32 = 0.1e1_dp + alpha_1_2 * rs
           t34 = 0.1e1_dp / A_2
           t38 = p_2 + 0.1e1_dp
           t39 = rs ** t38
           t40 = beta_4_2 * t39
           t41 = beta_1_2 * t15 + beta_2_2 * rs + beta_3_2 * t18 + t40
           t45 = 0.1e1_dp + t34 / t41 / 0.2e1_dp
           t46 = LOG(t45)
           t50 = 0.1e1_dp + alpha_1_3 * rs
           t52 = 0.1e1_dp / A_3
           t56 = p_3 + 0.1e1_dp
           t57 = rs ** t56
           t58 = beta_4_3 * t57
           t59 = beta_1_3 * t15 + beta_2_3 * rs + beta_3_3 * t18 + t58
           t63 = 0.1e1_dp + t52 / t59 / 0.2e1_dp
           t64 = LOG(t63)
           alpha_c = 0.2e1_dp * A_3 * t50 * t64
           t66 = 2 ** (0.1e1_dp / 0.3e1_dp)
           t69 = 1 / (2 * t66 - 2)
           t70 = 0.1e1_dp + chi
           t71 = t70 ** (0.1e1_dp / 0.3e1_dp)
           t72 = t71 * t70
           t73 = 0.1e1_dp - chi
           t74 = t73 ** (0.1e1_dp / 0.3e1_dp)
           t75 = t74 * t73
           f = (t72 + t75 - 0.2e1_dp) * t69
           t77 = alpha_c * f
           t78 = 0.9e1_dp / 0.8e1_dp / t69
           t79 = chi ** 2
           t80 = t79 ** 2
           t82 = t78 * (0.1e1_dp - t80)
           t84 = -0.2e1_dp * A_2 * t32 * t46 - e_c_u_0
           t85 = t84 * f
           epsilon_c_unif = e_c_u_0 + t77 * t82 + t85 * t80
           t87 = t71 ** 2
           t88 = t74 ** 2
           phi = t87 / 0.2e1_dp + t88 / 0.2e1_dp
           t91 = t90 * my_rho
           t92 = t91 ** (0.1e1_dp / 0.3e1_dp)
           t93 = t3 * t92 * t7
           t94 = SQRT(t93)
           k_s = 0.2e1_dp * t94
           t95 = 0.1e1_dp / phi
           t96 = my_norm_drho * t95
           t97 = 0.1e1_dp / k_s
           t98 = t97 * t2
           t = t96 * t98 / 0.2e1_dp
           t100 = 0.1e1_dp / gamma_var
           t101 = beta * t100
           t102 = epsilon_c_unif * t100
           t103 = phi ** 2
           t104 = t103 * phi
           t105 = 0.1e1_dp / t104
           t107 = EXP(-t102 * t105)
           t108 = t107 - 0.1e1_dp
           A = t101 / t108
           t110 = gamma_var * t104
           t111 = t ** 2
           t112 = A * t111
           t113 = 0.1e1_dp + t112
           t115 = A ** 2
           t116 = t111 ** 2
           t118 = 0.1e1_dp + t112 + t115 * t116
           t119 = 0.1e1_dp / t118
           t122 = 0.1e1_dp + t101 * t111 * t113 * t119
           t123 = LOG(t122)
           epsilon_cGGA = epsilon_c_unif + t110 * t123
           t124 = t3 * t66
           t125 = t90 * my_rhoa
           t126 = t125 ** (0.1e1_dp / 0.3e1_dp)
           kf_a = t124 * t126
           ex_unif_a = -0.3e1_dp / 0.4e1_dp * t7 * kf_a
           t129 = 0.1e1_dp / kf_a
           t130 = my_norm_drhoa * t129
           t131 = 0.1e1_dp / my_rhoa
           s_a = t130 * t131 / 0.2e1_dp
           t133 = s_a ** 2
           t135 = 0.1e1_dp / kappa
           t137 = 0.1e1_dp + mu * t133 * t135
           Fx_a = 0.1e1_dp + kappa - kappa / t137
           t140 = my_rhoa * ex_unif_a
           t142 = t90 * my_rhob
           t143 = t142 ** (0.1e1_dp / 0.3e1_dp)
           kf_b = t124 * t143
           ex_unif_b = -0.3e1_dp / 0.4e1_dp * t7 * kf_b
           t146 = 0.1e1_dp / kf_b
           t147 = my_norm_drhob * t146
           t148 = 0.1e1_dp / my_rhob
           s_b = t147 * t148 / 0.2e1_dp
           t150 = s_b ** 2
           t153 = 0.1e1_dp + mu * t150 * t135
           Fx_b = 0.1e1_dp + kappa - kappa / t153
           t156 = my_rhob * ex_unif_b

           IF (grad_deriv>=0) THEN
              e_0(ii) = e_0(ii)+&
                   scale_ex * (0.2e1_dp * t140 * Fx_a + 0.2e1_dp * t156 * Fx_b)&
                   / 0.2e1_dp + scale_ec * my_rho * epsilon_cGGA
           END IF

           t162 = my_rho ** 2
           t163 = 0.1e1_dp / t162
           t164 = t1 * t163
           chirhoa = t2 - t164
           t165 = t9 ** 2
           t167 = 0.1e1_dp / t165 * t7
           rsrhoa = -t6 * t167 * t163 / 0.12e2_dp
           t171 = A_1 * alpha_1_1
           t175 = t23 ** 2
           t176 = 0.1e1_dp / t175
           t177 = t12 * t176
           t178 = 0.1e1_dp / t15
           t179 = beta_1_1 * t178
           t183 = beta_3_1 * t15
           t187 = 0.1e1_dp / rs
           t190 = t179 * rsrhoa / 0.2e1_dp + beta_2_1 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t183 * rsrhoa + t22 * t20 * rsrhoa * t187
           t191 = 0.1e1_dp / t27
           t192 = t190 * t191
           e_c_u_0rhoa = -0.2e1_dp * t171 * rsrhoa * t28 + t177 * t192
           t194 = A_2 * alpha_1_2
           t198 = t41 ** 2
           t199 = 0.1e1_dp / t198
           t200 = t32 * t199
           t201 = beta_1_2 * t178
           t205 = beta_3_2 * t15
           t211 = t201 * rsrhoa / 0.2e1_dp + beta_2_2 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t205 * rsrhoa + t40 * t38 * rsrhoa * t187
           t212 = 0.1e1_dp / t45
           t213 = t211 * t212
           e_c_u_1rhoa = -0.2e1_dp * t194 * rsrhoa * t46 + t200 * t213
           t215 = A_3 * alpha_1_3
           t219 = t59 ** 2
           t220 = 0.1e1_dp / t219
           t221 = t50 * t220
           t222 = beta_1_3 * t178
           t226 = beta_3_3 * t15
           t232 = t222 * rsrhoa / 0.2e1_dp + beta_2_3 * rsrhoa + 0.3e1_dp / &
                0.2e1_dp * t226 * rsrhoa + t58 * t56 * rsrhoa * t187
           t233 = 0.1e1_dp / t63
           t234 = t232 * t233
           alpha_crhoa = 0.2e1_dp * t215 * rsrhoa * t64 - t221 * t234
           frhoa = (0.4e1_dp / 0.3e1_dp * t71 * chirhoa - 0.4e1_dp / 0.3e1_dp&
                * t74 * chirhoa) * t69
           t240 = alpha_crhoa * f
           t242 = alpha_c * frhoa
           t244 = t79 * chi
           t245 = t78 * t244
           t246 = t245 * chirhoa
           t248 = 0.4e1_dp * t77 * t246
           t249 = e_c_u_1rhoa - e_c_u_0rhoa
           t250 = t249 * f
           t252 = t84 * frhoa
           t254 = t244 * chirhoa
           t256 = 0.4e1_dp * t85 * t254
           epsilon_c_unifrhoa = e_c_u_0rhoa + t240 * t82 + t242 * t82 - t248 &
                + t250 * t80 + t252 * t80 + t256
           t257 = 0.1e1_dp / t71
           t259 = 0.1e1_dp / t74
           phirhoa = t257 * chirhoa / 0.3e1_dp - t259 * chirhoa / 0.3e1_dp
           t262 = t92 ** 2
           k_frhoa = t3 / t262 * t90 / 0.3e1_dp
           t266 = 0.1e1_dp / t94
           k_srhoa = t266 * k_frhoa * t7
           t268 = 0.1e1_dp / t103
           t269 = my_norm_drho * t268
           t272 = k_s ** 2
           t273 = 0.1e1_dp / t272
           t274 = t273 * t2
           t277 = t97 * t163
           t278 = t96 * t277
           trhoa = -t269 * t98 * phirhoa / 0.2e1_dp - t96 * t274 * k_srhoa / &
                0.2e1_dp - t278 / 0.2e1_dp
           t280 = t108 ** 2
           t281 = 0.1e1_dp / t280
           t282 = epsilon_c_unifrhoa * t100
           t284 = t103 ** 2
           t285 = 0.1e1_dp / t284
           t286 = t285 * phirhoa
           t289 = -t282 * t105 + 0.3e1_dp * t102 * t286
           Arhoa = -t101 * t281 * t289 * t107
           t293 = gamma_var * t103
           t294 = t123 * phirhoa
           t297 = t101 * t
           t298 = t113 * t119
           t299 = t298 * trhoa
           t302 = Arhoa * t111
           t303 = A * t
           t305 = 0.2e1_dp * t303 * trhoa
           t306 = t302 + t305
           t310 = t101 * t111
           t311 = t118 ** 2
           t312 = 0.1e1_dp / t311
           t313 = t113 * t312
           t314 = A * t116
           t317 = t111 * t
           t318 = t115 * t317
           t321 = t302 + t305 + 0.2e1_dp * t314 * Arhoa + 0.4e1_dp * t318 * trhoa
           t324 = 0.2e1_dp * t297 * t299 + t101 * t111 * t306 * t119 - t310 *&
                t313 * t321
           t325 = 0.1e1_dp / t122
           t326 = t324 * t325
           epsilon_cGGArhoa = epsilon_c_unifrhoa + 0.3e1_dp * t293 * t294 + &
                t110 * t326
           t329 = t126 ** 2
           kf_arhoa = t124 / t329 * t90 / 0.3e1_dp
           ex_unif_arhoa = -0.3e1_dp / 0.4e1_dp * t7 * kf_arhoa
           t335 = kf_a ** 2
           t336 = 0.1e1_dp / t335
           t337 = my_norm_drhoa * t336
           t340 = my_rhoa ** 2
           t341 = 0.1e1_dp / t340
           s_arhoa = -t337 * t131 * kf_arhoa / 0.2e1_dp - t130 * t341 / 0.2e1_dp
           t344 = t137 ** 2
           t346 = 0.1e1_dp / t344 * mu
           Fx_arhoa = 0.2e1_dp * t346 * s_a * s_arhoa
           t350 = my_rhoa * ex_unif_arhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ra(ii) = e_ra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a * Fx_a + 0.2e1_dp * &
                   t350 * Fx_a + 0.2e1_dp * t140 * Fx_arhoa) / 0.2e1_dp + scale_ec * (&
                   epsilon_cGGA + my_rho * epsilon_cGGArhoa)
           END IF

           chirhob = -t2 - t164
           rsrhob = rsrhoa
           t368 = t179 * rsrhob / 0.2e1_dp + beta_2_1 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t183 * rsrhob + t22 * t20 * rsrhob * t187
           e_c_u_0rhob = -0.2e1_dp * t171 * rsrhob * t28 + t177 * t368 * t191
           t382 = t201 * rsrhob / 0.2e1_dp + beta_2_2 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t205 * rsrhob + t40 * t38 * rsrhob * t187
           e_c_u_1rhob = -0.2e1_dp * t194 * rsrhob * t46 + t200 * t382 * t212
           t396 = t222 * rsrhob / 0.2e1_dp + beta_2_3 * rsrhob + 0.3e1_dp / &
                0.2e1_dp * t226 * rsrhob + t58 * t56 * rsrhob * t187
           alpha_crhob = 0.2e1_dp * t215 * rsrhob * t64 - t221 * t396 * t233
           frhob = (0.4e1_dp / 0.3e1_dp * t71 * chirhob - 0.4e1_dp / 0.3e1_dp&
                * t74 * chirhob) * t69
           t403 = alpha_crhob * f
           t405 = alpha_c * frhob
           t407 = t245 * chirhob
           t409 = 0.4e1_dp * t77 * t407
           t410 = e_c_u_1rhob - e_c_u_0rhob
           t411 = t410 * f
           t413 = t84 * frhob
           t415 = t244 * chirhob
           t417 = 0.4e1_dp * t85 * t415
           epsilon_c_unifrhob = e_c_u_0rhob + t403 * t82 + t405 * t82 - t409 &
                + t411 * t80 + t413 * t80 + t417
           phirhob = t257 * chirhob / 0.3e1_dp - t259 * chirhob / 0.3e1_dp
           k_frhob = k_frhoa
           k_srhob = t266 * k_frhob * t7
           trhob = -t269 * t98 * phirhob / 0.2e1_dp - t96 * t274 * k_srhob / &
                0.2e1_dp - t278 / 0.2e1_dp
           t427 = epsilon_c_unifrhob * t100
           t429 = t285 * phirhob
           t432 = -t427 * t105 + 0.3e1_dp * t102 * t429
           Arhob = -t101 * t281 * t432 * t107
           t436 = t123 * phirhob
           t439 = t298 * trhob
           t442 = Arhob * t111
           t444 = 0.2e1_dp * t303 * trhob
           t445 = t442 + t444
           t453 = t442 + t444 + 0.2e1_dp * t314 * Arhob + 0.4e1_dp * t318 * trhob
           t456 = 0.2e1_dp * t297 * t439 + t101 * t111 * t445 * t119 - t310 *&
                t313 * t453
           t457 = t456 * t325
           epsilon_cGGArhob = epsilon_c_unifrhob + 0.3e1_dp * t293 * t436 + &
                t110 * t457
           t460 = t143 ** 2
           kf_brhob = t124 / t460 * t90 / 0.3e1_dp
           ex_unif_brhob = -0.3e1_dp / 0.4e1_dp * t7 * kf_brhob
           t466 = kf_b ** 2
           t467 = 0.1e1_dp / t466
           t468 = my_norm_drhob * t467
           t471 = my_rhob ** 2
           t472 = 0.1e1_dp / t471
           s_brhob = -t468 * t148 * kf_brhob / 0.2e1_dp - t147 * t472 / 0.2e1_dp
           t475 = t153 ** 2
           t477 = 0.1e1_dp / t475 * mu
           Fx_brhob = 0.2e1_dp * t477 * s_b * s_brhob
           t481 = my_rhob * ex_unif_brhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_rb(ii) = e_rb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b * Fx_b + 0.2e1_dp * &
                   t481 * Fx_b + 0.2e1_dp * t156 * Fx_brhob) / 0.2e1_dp + scale_ec * (&
                   epsilon_cGGA + my_rho * epsilon_cGGArhob)
           END IF

           t488 = t95 * t97
           tnorm_drho = t488 * t2 / 0.2e1_dp
           t493 = t101 * t317
           t494 = A * tnorm_drho
           t502 = 0.2e1_dp * t303 * tnorm_drho + 0.4e1_dp * t318 * tnorm_drho
           t505 = 0.2e1_dp * t297 * t298 * tnorm_drho + 0.2e1_dp * t493 * &
                t494 * t119 - t310 * t313 * t502
           t506 = t505 * t325
           Hnorm_drho = t110 * t506

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndr(ii) = e_ndr(ii)+&
                   scale_ec * my_rho * Hnorm_drho
           END IF

           s_anorm_drhoa = t129 * t131 / 0.2e1_dp
           Fx_anorm_drhoa = 0.2e1_dp * t346 * s_a * s_anorm_drhoa

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndra(ii) = e_ndra(ii)+&
                   scale_ex * t140 * Fx_anorm_drhoa
           END IF

           s_bnorm_drhob = t146 * t148 / 0.2e1_dp
           Fx_bnorm_drhob = 0.2e1_dp * t477 * s_b * s_bnorm_drhob

           IF (grad_deriv>=1.OR.grad_deriv==-1) THEN
              e_ndrb(ii) = e_ndrb(ii)+&
                   scale_ex * t156 * Fx_bnorm_drhob
           END IF

           IF (grad_deriv>=2.OR.grad_deriv==-2) THEN
              t518 = 0.1e1_dp / t162 / my_rho
              t519 = t1 * t518
              chirhoarhoa = -0.2e1_dp * t163 + 0.2e1_dp * t519
              t523 = 0.1e1_dp / t90
              t525 = t162 ** 2
              rsrhoarhoa = -t6 / t165 / t8 * t523 / t525 / 0.18e2_dp + &
                   t6 * t167 * t518 / 0.6e1_dp
              t536 = alpha_1_1 * rsrhoa
              t538 = t176 * t190 * t191
              t543 = t12 / t175 / t23
              t544 = t190 ** 2
              t548 = 0.1e1_dp / t18
              t549 = beta_1_1 * t548
              t550 = rsrhoa ** 2
              t556 = beta_3_1 * t178
              t561 = t20 ** 2
              t563 = rs ** 2
              t564 = 0.1e1_dp / t563
              t576 = t175 ** 2
              t578 = t12 / t576
              t579 = t27 ** 2
              t580 = 0.1e1_dp / t579
              e_c_u_0rhoarhoa = -0.2e1_dp * t171 * rsrhoarhoa * t28 + 0.2e1_dp *&
                   t536 * t538 - 0.2e1_dp * t543 * t544 * t191 + t177 * (-t549 * t550 &
                   / 0.4e1_dp + t179 * rsrhoarhoa / 0.2e1_dp + beta_2_1 * rsrhoarhoa + &
                   0.3e1_dp / 0.4e1_dp * t556 * t550 + 0.3e1_dp / 0.2e1_dp * t183 * &
                   rsrhoarhoa + t22 * t561 * t550 * t564 + t22 * t20 * rsrhoarhoa * &
                   t187 - t22 * t20 * t550 * t564) * t191 + t578 * t544 * t580 * t14 / &
                   0.2e1_dp
              e_c_u_01rhoa = e_c_u_0rhoa
              t588 = alpha_1_2 * rsrhoa
              t590 = t199 * t211 * t212
              t595 = t32 / t198 / t41
              t596 = t211 ** 2
              t600 = beta_1_2 * t548
              t606 = beta_3_2 * t178
              t611 = t38 ** 2
              t624 = t198 ** 2
              t626 = t32 / t624
              t627 = t45 ** 2
              t628 = 0.1e1_dp / t627
              t636 = alpha_1_3 * rsrhoa
              t638 = t220 * t232 * t233
              t643 = t50 / t219 / t59
              t644 = t232 ** 2
              t648 = beta_1_3 * t548
              t654 = beta_3_3 * t178
              t659 = t56 ** 2
              t672 = t219 ** 2
              t674 = t50 / t672
              t675 = t63 ** 2
              t676 = 0.1e1_dp / t675
              alpha_c1rhoa = alpha_crhoa
              t681 = 0.1e1_dp / t87
              t682 = chirhoa ** 2
              t687 = 0.1e1_dp / t88
              frhoarhoa = (0.4e1_dp / 0.9e1_dp * t681 * t682 + 0.4e1_dp / &
                   0.3e1_dp * t71 * chirhoarhoa + 0.4e1_dp / 0.9e1_dp * t687 * t682 - &
                   0.4e1_dp / 0.3e1_dp * t74 * chirhoarhoa) * t69
              f1rhoa = frhoa
              t705 = alpha_c1rhoa * f
              t708 = alpha_c * f1rhoa
              t711 = t78 * t79
              t726 = e_c_u_1rhoa - e_c_u_01rhoa
              t733 = t726 * f
              t736 = t84 * f1rhoa
              t745 = -0.4e1_dp * t77 * t245 * chirhoarhoa + (-0.2e1_dp * t194 * &
                   rsrhoarhoa * t46 + 0.2e1_dp * t588 * t590 - 0.2e1_dp * t595 * t596 *&
                   t212 + t200 * (-t600 * t550 / 0.4e1_dp + t201 * rsrhoarhoa / &
                   0.2e1_dp + beta_2_2 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t606 * t550&
                   + 0.3e1_dp / 0.2e1_dp * t205 * rsrhoarhoa + t40 * t611 * t550 * &
                   t564 + t40 * t38 * rsrhoarhoa * t187 - t40 * t38 * t550 * t564) * &
                   t212 + t626 * t596 * t628 * t34 / 0.2e1_dp - e_c_u_0rhoarhoa) * f * &
                   t80 + t249 * f1rhoa * t80 + 0.4e1_dp * t250 * t254 + t726 * frhoa * &
                   t80 + t84 * frhoarhoa * t80 + 0.4e1_dp * t252 * t254 + 0.4e1_dp * &
                   t733 * t254 + 0.4e1_dp * t736 * t254 + 0.12e2_dp * t85 * t79 * t682 &
                   + 0.4e1_dp * t85 * t244 * chirhoarhoa
              epsilon_c_unifrhoarhoa = e_c_u_0rhoarhoa + (0.2e1_dp * t215 * &
                   rsrhoarhoa * t64 - 0.2e1_dp * t636 * t638 + 0.2e1_dp * t643 * t644 *&
                   t233 - t221 * (-t648 * t550 / 0.4e1_dp + t222 * rsrhoarhoa / &
                   0.2e1_dp + beta_2_3 * rsrhoarhoa + 0.3e1_dp / 0.4e1_dp * t654 * t550&
                   + 0.3e1_dp / 0.2e1_dp * t226 * rsrhoarhoa + t58 * t659 * t550 * &
                   t564 + t58 * t56 * rsrhoarhoa * t187 - t58 * t56 * t550 * t564) * &
                   t233 - t674 * t644 * t676 * t52 / 0.2e1_dp) * f * t82 + alpha_crhoa &
                   * f1rhoa * t82 - 0.4e1_dp * t240 * t246 + alpha_c1rhoa * frhoa * t82&
                   + alpha_c * frhoarhoa * t82 - 0.4e1_dp * t242 * t246 - 0.4e1_dp * &
                   t705 * t246 - 0.4e1_dp * t708 * t246 - 0.12e2_dp * t77 * t711 * t682&
                   + t745
              epsilon_c_unif1rhoa = e_c_u_01rhoa + t705 * t82 + t708 * t82 - &
                   t248 + t733 * t80 + t736 * t80 + t256
              t750 = 0.1e1_dp / t72
              t755 = 0.1e1_dp / t75
              phirhoarhoa = -t750 * t682 / 0.9e1_dp + t257 * chirhoarhoa / &
                   0.3e1_dp - t755 * t682 / 0.9e1_dp - t259 * chirhoarhoa / 0.3e1_dp
              phi1rhoa = phirhoa
              t763 = t90 ** 2
              k_frhoarhoa = -0.2e1_dp / 0.9e1_dp * t3 / t262 / t91 * t763
              t767 = 0.1e1_dp / t94 / t93
              t768 = k_frhoa ** 2
              k_s1rhoa = k_srhoa
              t775 = my_norm_drho * t105 * t97
              t776 = t2 * phirhoa
              t779 = t269 * t273
              t785 = t269 * t277 * phirhoa / 0.2e1_dp
              t789 = t2 * k_srhoa
              t795 = t96 / t272 / k_s
              t798 = t273 * t163
              t801 = t96 * t798 * k_srhoa / 0.2e1_dp
              t812 = t96 * t97 * t518
              trhoarhoa = t775 * t776 * phi1rhoa + t779 * t776 * k_s1rhoa / &
                   0.2e1_dp + t785 - t269 * t98 * phirhoarhoa / 0.2e1_dp + t779 * t789 &
                   * phi1rhoa / 0.2e1_dp + t795 * t789 * k_s1rhoa + t801 - t96 * t274 *&
                   (-t767 * t768 * t523 / 0.2e1_dp + t266 * k_frhoarhoa * t7) / &
                   0.2e1_dp + t269 * t277 * phi1rhoa / 0.2e1_dp + t96 * t798 * k_s1rhoa&
                   / 0.2e1_dp + t812
              t1rhoa = -t269 * t98 * phi1rhoa / 0.2e1_dp - t96 * t274 * k_s1rhoa&
                   / 0.2e1_dp - t278 / 0.2e1_dp
              t820 = t101 / t280 / t108
              t821 = t107 ** 2
              t822 = t289 * t821
              t823 = epsilon_c_unif1rhoa * t100
              t825 = t285 * phi1rhoa
              t828 = -t823 * t105 + 0.3e1_dp * t102 * t825
              t839 = 0.1e1_dp / t284 / phi
              t840 = t839 * phirhoa
              t851 = t101 * t281
              Arhoarhoa = 0.2e1_dp * t820 * t822 * t828 - t101 * t281 * (&
                   -epsilon_c_unifrhoarhoa * t100 * t105 + 0.3e1_dp * t282 * t825 + &
                   0.3e1_dp * t823 * t286 - 0.12e2_dp * t102 * t840 * phi1rhoa + &
                   0.3e1_dp * t102 * t285 * phirhoarhoa) * t107 - t851 * t289 * t828 * &
                   t107
              A1rhoa = -t101 * t281 * t828 * t107
              t858 = gamma_var * phi
              t865 = A1rhoa * t111
              t867 = 0.2e1_dp * t303 * t1rhoa
              t868 = t865 + t867
              t876 = t865 + t867 + 0.2e1_dp * t314 * A1rhoa + 0.4e1_dp * t318 * t1rhoa
              t879 = 0.2e1_dp * t297 * t298 * t1rhoa + t101 * t111 * t868 * t119&
                   - t310 * t313 * t876
              t880 = t879 * t325
              t904 = t306 * t119
              t908 = Arhoarhoa * t111
              t909 = Arhoa * t
              t911 = 0.2e1_dp * t909 * t1rhoa
              t914 = 0.2e1_dp * A1rhoa * t * trhoa
              t917 = 0.2e1_dp * A * t1rhoa * trhoa
              t919 = 0.2e1_dp * t303 * trhoarhoa
              t924 = t306 * t312
              t936 = t113 / t311 / t118
              t944 = A * t317
              t953 = t115 * t111
              t959 = t908 + t911 + t914 + t917 + t919 + 0.2e1_dp * A1rhoa * t116&
                   * Arhoa + 0.8e1_dp * t944 * Arhoa * t1rhoa + 0.2e1_dp * t314 * &
                   Arhoarhoa + 0.8e1_dp * t944 * trhoa * A1rhoa + 0.12e2_dp * t953 * &
                   trhoa * t1rhoa + 0.4e1_dp * t318 * trhoarhoa
              t962 = 0.2e1_dp * t101 * t1rhoa * t299 + 0.2e1_dp * t297 * t868 * &
                   t119 * trhoa - 0.2e1_dp * t297 * t313 * trhoa * t876 + 0.2e1_dp * &
                   t297 * t298 * trhoarhoa + 0.2e1_dp * t297 * t904 * t1rhoa + t101 * &
                   t111 * (t908 + t911 + t914 + t917 + t919) * t119 - t310 * t924 * &
                   t876 - 0.2e1_dp * t297 * t313 * t321 * t1rhoa - t310 * t868 * t312 *&
                   t321 + 0.2e1_dp * t310 * t936 * t321 * t876 - t310 * t313 * t959
              t965 = t122 ** 2
              t966 = 0.1e1_dp / t965
              t967 = t324 * t966
              kf_arhoarhoa = -0.2e1_dp / 0.9e1_dp * t124 / t329 / t125 * t763
              ex_unif_a1rhoa = ex_unif_arhoa
              t985 = kf_arhoa ** 2
              s_a1rhoa = s_arhoa
              t998 = mu ** 2
              t999 = 0.1e1_dp / t344 / t137 * t998
              t1000 = t999 * t133
              t1001 = s_arhoa * t135
              Fx_a1rhoa = 0.2e1_dp * t346 * s_a * s_a1rhoa

              e_ra_ra(ii) = e_ra_ra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a1rhoa * Fx_a + &
                   0.2e1_dp * ex_unif_a * Fx_a1rhoa + 0.2e1_dp * ex_unif_arhoa * Fx_a -&
                   0.3e1_dp / 0.2e1_dp * my_rhoa * t7 * kf_arhoarhoa * Fx_a + 0.2e1_dp * &
                   t350 * Fx_a1rhoa + 0.2e1_dp * ex_unif_a * Fx_arhoa + 0.2e1_dp * my_rhoa&
                   * ex_unif_a1rhoa * Fx_arhoa + 0.2e1_dp * t140 * (-0.8e1_dp * t1000 &
                   * t1001 * s_a1rhoa + 0.2e1_dp * t346 * s_a1rhoa * s_arhoa + 0.2e1_dp&
                   * t346 * s_a * (my_norm_drhoa / t335 / kf_a * t131 * t985 + t337 * &
                   t341 * kf_arhoa - t337 * t131 * kf_arhoarhoa / 0.2e1_dp + t130 / &
                   t340 / my_rhoa))) / 0.2e1_dp + scale_ec * (epsilon_c_unif1rhoa + &
                   0.3e1_dp * t293 * t123 * phi1rhoa + t110 * t880 + epsilon_cGGArhoa +&
                   my_rho * (epsilon_c_unifrhoarhoa + 0.6e1_dp * t858 * t294 * phi1rhoa +&
                   0.3e1_dp * t293 * t880 * phirhoa + 0.3e1_dp * t293 * t123 * &
                   phirhoarhoa + 0.3e1_dp * t293 * t326 * phi1rhoa + t110 * t962 * t325&
                   - t110 * t967 * t879))

              chirhoarhob = 0.2e1_dp * t519
              rsrhoarhob = rsrhoarhoa
              t1031 = t176 * t368 * t191
              t1033 = alpha_1_1 * rsrhob
              t1038 = rsrhoa * rsrhob
              t1050 = rsrhob * t564 * rsrhoa
              e_c_u_0rhoarhob = -0.2e1_dp * t171 * rsrhoarhob * t28 + t536 * &
                   t1031 + t1033 * t538 - 0.2e1_dp * t543 * t192 * t368 + t177 * (-t549&
                   * t1038 / 0.4e1_dp + t179 * rsrhoarhob / 0.2e1_dp + beta_2_1 * &
                   rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t556 * t1038 + 0.3e1_dp / &
                   0.2e1_dp * t183 * rsrhoarhob + t22 * t561 * t1050 + t22 * t20 * &
                   rsrhoarhob * t187 - t22 * t20 * t1050) * t191 + t578 * t190 * t580 *&
                   t14 * t368 / 0.2e1_dp
              t1069 = t199 * t382 * t212
              t1071 = alpha_1_2 * rsrhob
              t1104 = t220 * t396 * t233
              t1106 = alpha_1_3 * rsrhob
              frhoarhob = (0.4e1_dp / 0.9e1_dp * t681 * chirhoa * chirhob + &
                   0.4e1_dp / 0.3e1_dp * t71 * chirhoarhob + 0.4e1_dp / 0.9e1_dp * t687&
                   * chirhoa * chirhob - 0.4e1_dp / 0.3e1_dp * t74 * chirhoarhob) * &
                   t69
              t1164 = t79 * chirhoa * chirhob
              t1193 = -0.4e1_dp * t77 * t245 * chirhoarhob + (-0.2e1_dp * t194 *&
                   rsrhoarhob * t46 + t588 * t1069 + t1071 * t590 - 0.2e1_dp * t595 * &
                   t213 * t382 + t200 * (-t600 * t1038 / 0.4e1_dp + t201 * rsrhoarhob /&
                   0.2e1_dp + beta_2_2 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t606 * &
                   t1038 + 0.3e1_dp / 0.2e1_dp * t205 * rsrhoarhob + t40 * t611 * t1050&
                   + t40 * t38 * rsrhoarhob * t187 - t40 * t38 * t1050) * t212 + t626 &
                   * t211 * t628 * t34 * t382 / 0.2e1_dp - e_c_u_0rhoarhob) * f * t80 +&
                   t249 * frhob * t80 + 0.4e1_dp * t250 * t415 + t410 * frhoa * t80 + &
                   t84 * frhoarhob * t80 + 0.4e1_dp * t252 * t415 + 0.4e1_dp * t411 * &
                   t254 + 0.4e1_dp * t413 * t254 + 0.12e2_dp * t85 * t1164 + 0.4e1_dp *&
                   t85 * t244 * chirhoarhob
              epsilon_c_unifrhoarhob = e_c_u_0rhoarhob + (0.2e1_dp * t215 * &
                   rsrhoarhob * t64 - t636 * t1104 - t1106 * t638 + 0.2e1_dp * t643 * &
                   t234 * t396 - t221 * (-t648 * t1038 / 0.4e1_dp + t222 * rsrhoarhob /&
                   0.2e1_dp + beta_2_3 * rsrhoarhob + 0.3e1_dp / 0.4e1_dp * t654 * &
                   t1038 + 0.3e1_dp / 0.2e1_dp * t226 * rsrhoarhob + t58 * t659 * t1050&
                   + t58 * t56 * rsrhoarhob * t187 - t58 * t56 * t1050) * t233 - t674 &
                   * t232 * t676 * t52 * t396 / 0.2e1_dp) * f * t82 + alpha_crhoa * &
                   frhob * t82 - 0.4e1_dp * t240 * t407 + alpha_crhob * frhoa * t82 + &
                   alpha_c * frhoarhob * t82 - 0.4e1_dp * t242 * t407 - 0.4e1_dp * t403&
                   * t246 - 0.4e1_dp * t405 * t246 - 0.12e2_dp * t77 * t78 * t1164 + &
                   t1193
              phirhoarhob = -t750 * chirhoa * chirhob / 0.9e1_dp + t257 * &
                   chirhoarhob / 0.3e1_dp - t755 * chirhoa * chirhob / 0.9e1_dp - t259 &
                   * chirhoarhob / 0.3e1_dp
              k_frhoarhob = k_frhoarhoa
              t1228 = t269 * t277 * phirhob / 0.2e1_dp
              t1231 = t96 * t798 * k_srhob / 0.2e1_dp
              trhoarhob = t775 * t776 * phirhob + t779 * t776 * k_srhob / &
                   0.2e1_dp + t785 - t269 * t98 * phirhoarhob / 0.2e1_dp + t779 * t789 &
                   * phirhob / 0.2e1_dp + t795 * t789 * k_srhob + t801 - t96 * t274 * (&
                   -t767 * k_frhoa * t523 * k_frhob / 0.2e1_dp + t266 * k_frhoarhob * &
                   t7) / 0.2e1_dp + t1228 + t1231 + t812
              Arhoarhob = 0.2e1_dp * t820 * t822 * t432 - t101 * t281 * (&
                   -epsilon_c_unifrhoarhob * t100 * t105 + 0.3e1_dp * t282 * t429 + &
                   0.3e1_dp * t427 * t286 - 0.12e2_dp * t102 * t840 * phirhob + &
                   0.3e1_dp * t102 * t285 * phirhoarhob) * t107 - t851 * t289 * t432 * &
                   t107
              t1269 = t445 * t119
              t1283 = Arhoarhob * t111
              t1285 = 0.2e1_dp * t909 * trhob
              t1286 = Arhob * t
              t1288 = 0.2e1_dp * t1286 * trhoa
              t1291 = 0.2e1_dp * A * trhob * trhoa
              t1293 = 0.2e1_dp * t303 * trhoarhob
              t1304 = t445 * t312
              t1327 = t1283 + t1285 + t1288 + t1291 + t1293 + 0.2e1_dp * Arhob *&
                   t116 * Arhoa + 0.8e1_dp * t944 * Arhoa * trhob + 0.2e1_dp * t314 * &
                   Arhoarhob + 0.8e1_dp * t944 * trhoa * Arhob + 0.12e2_dp * t953 * &
                   trhoa * trhob + 0.4e1_dp * t318 * trhoarhob
              t1330 = 0.2e1_dp * t101 * trhob * t299 + 0.2e1_dp * t297 * t1269 *&
                   trhoa - 0.2e1_dp * t297 * t313 * trhoa * t453 + 0.2e1_dp * t297 * &
                   t298 * trhoarhob + 0.2e1_dp * t297 * t904 * trhob + t101 * t111 * (&
                   t1283 + t1285 + t1288 + t1291 + t1293) * t119 - t310 * t924 * t453 -&
                   0.2e1_dp * t297 * t313 * t321 * trhob - t310 * t1304 * t321 + &
                   0.2e1_dp * t310 * t936 * t321 * t453 - t310 * t313 * t1327

              e_ra_rb(ii) = e_ra_rb(ii)+&
                   scale_ec * (epsilon_cGGArhob + epsilon_cGGArhoa + &
                   my_rho * (epsilon_c_unifrhoarhob + 0.6e1_dp * t858 * t294 * phirhob + &
                   0.3e1_dp * t293 * t457 * phirhoa + 0.3e1_dp * t293 * t123 * &
                   phirhoarhob + 0.3e1_dp * t293 * t326 * phirhob + t110 * t1330 * t325&
                   - t110 * t967 * t456))

              chirhobrhob = 0.2e1_dp * t163 + 0.2e1_dp * t519
              rsrhobrhob = rsrhoarhob
              t1342 = t368 ** 2
              t1346 = rsrhob ** 2
              e_c_u_0rhobrhob = -0.2e1_dp * t171 * rsrhobrhob * t28 + 0.2e1_dp *&
                   t1033 * t1031 - 0.2e1_dp * t543 * t1342 * t191 + t177 * (-t549 * &
                   t1346 / 0.4e1_dp + t179 * rsrhobrhob / 0.2e1_dp + beta_2_1 * &
                   rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t556 * t1346 + 0.3e1_dp / &
                   0.2e1_dp * t183 * rsrhobrhob + t22 * t561 * t1346 * t564 + t22 * t20&
                   * rsrhobrhob * t187 - t22 * t20 * t1346 * t564) * t191 + t578 * &
                   t1342 * t580 * t14 / 0.2e1_dp
              e_c_u_01rhob = e_c_u_0rhob
              t1377 = t382 ** 2
              t1411 = t396 ** 2
              alpha_c1rhob = alpha_crhob
              t1440 = chirhob ** 2
              frhobrhob = (0.4e1_dp / 0.9e1_dp * t681 * t1440 + 0.4e1_dp / &
                   0.3e1_dp * t71 * chirhobrhob + 0.4e1_dp / 0.9e1_dp * t687 * t1440 - &
                   0.4e1_dp / 0.3e1_dp * t74 * chirhobrhob) * t69
              f1rhob = frhob
              t1462 = alpha_c1rhob * f
              t1465 = alpha_c * f1rhob
              t1482 = e_c_u_1rhob - e_c_u_01rhob
              t1489 = t1482 * f
              t1492 = t84 * f1rhob
              t1501 = -0.4e1_dp * t77 * t245 * chirhobrhob + (-0.2e1_dp * t194 *&
                   rsrhobrhob * t46 + 0.2e1_dp * t1071 * t1069 - 0.2e1_dp * t595 * &
                   t1377 * t212 + t200 * (-t600 * t1346 / 0.4e1_dp + t201 * rsrhobrhob &
                   / 0.2e1_dp + beta_2_2 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t606 * &
                   t1346 + 0.3e1_dp / 0.2e1_dp * t205 * rsrhobrhob + t40 * t611 * t1346&
                   * t564 + t40 * t38 * rsrhobrhob * t187 - t40 * t38 * t1346 * t564) &
                   * t212 + t626 * t1377 * t628 * t34 / 0.2e1_dp - e_c_u_0rhobrhob) * f&
                   * t80 + t410 * f1rhob * t80 + 0.4e1_dp * t411 * t415 + t1482 * &
                   frhob * t80 + t84 * frhobrhob * t80 + 0.4e1_dp * t413 * t415 + &
                   0.4e1_dp * t1489 * t415 + 0.4e1_dp * t1492 * t415 + 0.12e2_dp * t85 &
                   * t79 * t1440 + 0.4e1_dp * t85 * t244 * chirhobrhob
              epsilon_c_unifrhobrhob = e_c_u_0rhobrhob + (0.2e1_dp * t215 * &
                   rsrhobrhob * t64 - 0.2e1_dp * t1106 * t1104 + 0.2e1_dp * t643 * &
                   t1411 * t233 - t221 * (-t648 * t1346 / 0.4e1_dp + t222 * rsrhobrhob &
                   / 0.2e1_dp + beta_2_3 * rsrhobrhob + 0.3e1_dp / 0.4e1_dp * t654 * &
                   t1346 + 0.3e1_dp / 0.2e1_dp * t226 * rsrhobrhob + t58 * t659 * t1346&
                   * t564 + t58 * t56 * rsrhobrhob * t187 - t58 * t56 * t1346 * t564) &
                   * t233 - t674 * t1411 * t676 * t52 / 0.2e1_dp) * f * t82 + &
                   alpha_crhob * f1rhob * t82 - 0.4e1_dp * t403 * t407 + alpha_c1rhob *&
                   frhob * t82 + alpha_c * frhobrhob * t82 - 0.4e1_dp * t405 * t407 - &
                   0.4e1_dp * t1462 * t407 - 0.4e1_dp * t1465 * t407 - 0.12e2_dp * t77 &
                   * t711 * t1440 + t1501
              epsilon_c_unif1rhob = e_c_u_01rhob + t1462 * t82 + t1465 * t82 - &
                   t409 + t1489 * t80 + t1492 * t80 + t417
              phirhobrhob = -t750 * t1440 / 0.9e1_dp + t257 * chirhobrhob / &
                   0.3e1_dp - t755 * t1440 / 0.9e1_dp - t259 * chirhobrhob / 0.3e1_dp
              phi1rhob = phirhob
              t1514 = k_frhob ** 2
              k_s1rhob = k_srhob
              t1520 = t2 * phirhob
              t1529 = t2 * k_srhob
              trhobrhob = t775 * t1520 * phi1rhob + t779 * t1520 * k_s1rhob / &
                   0.2e1_dp + t1228 - t269 * t98 * phirhobrhob / 0.2e1_dp + t779 * &
                   t1529 * phi1rhob / 0.2e1_dp + t795 * t1529 * k_s1rhob + t1231 - t96 &
                   * t274 * (-t767 * t1514 * t523 / 0.2e1_dp + t266 * k_frhoarhob * t7)&
                   / 0.2e1_dp + t269 * t277 * phi1rhob / 0.2e1_dp + t96 * t798 * &
                   k_s1rhob / 0.2e1_dp + t812
              t1rhob = -t269 * t98 * phi1rhob / 0.2e1_dp - t96 * t274 * k_s1rhob&
                   / 0.2e1_dp - t278 / 0.2e1_dp
              t1550 = epsilon_c_unif1rhob * t100
              t1552 = t285 * phi1rhob
              t1555 = -t1550 * t105 + 0.3e1_dp * t102 * t1552
              Arhobrhob = 0.2e1_dp * t820 * t432 * t821 * t1555 - t101 * t281 * &
                   (-epsilon_c_unifrhobrhob * t100 * t105 + 0.3e1_dp * t427 * t1552 + &
                   0.3e1_dp * t1550 * t429 - 0.12e2_dp * t102 * t839 * phirhob * &
                   phi1rhob + 0.3e1_dp * t102 * t285 * phirhobrhob) * t107 - t851 * &
                   t432 * t1555 * t107
              A1rhob = -t101 * t281 * t1555 * t107
              t1588 = A1rhob * t111
              t1590 = 0.2e1_dp * t303 * t1rhob
              t1591 = t1588 + t1590
              t1599 = t1588 + t1590 + 0.2e1_dp * t314 * A1rhob + 0.4e1_dp * t318&
                   * t1rhob
              t1602 = 0.2e1_dp * t297 * t298 * t1rhob + t101 * t111 * t1591 * &
                   t119 - t310 * t313 * t1599
              t1603 = t1602 * t325
              t1630 = Arhobrhob * t111
              t1632 = 0.2e1_dp * t1286 * t1rhob
              t1635 = 0.2e1_dp * A1rhob * t * trhob
              t1638 = 0.2e1_dp * A * t1rhob * trhob
              t1640 = 0.2e1_dp * t303 * trhobrhob
              t1674 = t1630 + t1632 + t1635 + t1638 + t1640 + 0.2e1_dp * A1rhob &
                   * t116 * Arhob + 0.8e1_dp * t944 * Arhob * t1rhob + 0.2e1_dp * t314 &
                   * Arhobrhob + 0.8e1_dp * t944 * trhob * A1rhob + 0.12e2_dp * t953 * &
                   trhob * t1rhob + 0.4e1_dp * t318 * trhobrhob
              t1677 = 0.2e1_dp * t101 * t1rhob * t439 + 0.2e1_dp * t297 * t1591 &
                   * t119 * trhob - 0.2e1_dp * t297 * t313 * trhob * t1599 + 0.2e1_dp *&
                   t297 * t298 * trhobrhob + 0.2e1_dp * t297 * t1269 * t1rhob + t101 *&
                   t111 * (t1630 + t1632 + t1635 + t1638 + t1640) * t119 - t310 * &
                   t1304 * t1599 - 0.2e1_dp * t297 * t313 * t453 * t1rhob - t310 * &
                   t1591 * t312 * t453 + 0.2e1_dp * t310 * t936 * t453 * t1599 - t310 *&
                   t313 * t1674
              t1680 = t456 * t966
              kf_brhobrhob = -0.2e1_dp / 0.9e1_dp * t124 / t460 / t142 * t763
              ex_unif_b1rhob = ex_unif_brhob
              t1698 = kf_brhob ** 2
              s_b1rhob = s_brhob
              t1711 = 0.1e1_dp / t475 / t153 * t998
              t1712 = t1711 * t150
              t1713 = s_brhob * t135
              Fx_b1rhob = 0.2e1_dp * t477 * s_b * s_b1rhob

              e_rb_rb(ii) = e_rb_rb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b1rhob * Fx_b + &
                   0.2e1_dp * ex_unif_b * Fx_b1rhob + 0.2e1_dp * ex_unif_brhob * Fx_b -&
                   0.3e1_dp / 0.2e1_dp * my_rhob * t7 * kf_brhobrhob * Fx_b + 0.2e1_dp * &
                   t481 * Fx_b1rhob + 0.2e1_dp * ex_unif_b * Fx_brhob + 0.2e1_dp * my_rhob&
                   * ex_unif_b1rhob * Fx_brhob + 0.2e1_dp * t156 * (-0.8e1_dp * t1712 &
                   * t1713 * s_b1rhob + 0.2e1_dp * t477 * s_b1rhob * s_brhob + 0.2e1_dp&
                   * t477 * s_b * (my_norm_drhob / t466 / kf_b * t148 * t1698 + t468 * &
                   t472 * kf_brhob - t468 * t148 * kf_brhobrhob / 0.2e1_dp + t147 / &
                   t471 / my_rhob))) / 0.2e1_dp + scale_ec * (epsilon_c_unif1rhob + &
                   0.3e1_dp * t293 * t123 * phi1rhob + t110 * t1603 + epsilon_cGGArhob &
                   + my_rho * (epsilon_c_unifrhobrhob + 0.6e1_dp * t858 * t436 * phi1rhob &
                   + 0.3e1_dp * t293 * t1603 * phirhob + 0.3e1_dp * t293 * t123 * &
                   phirhobrhob + 0.3e1_dp * t293 * t457 * phi1rhob + t110 * t1677 * &
                   t325 - t110 * t1680 * t1602))
              t1739 = t268 * t97
              t1741 = t95 * t273
              t1743 = t488 * t163
              trhoanorm_drho = -t1739 * t776 / 0.2e1_dp - t1741 * t789 / &
                   0.2e1_dp - t1743 / 0.2e1_dp
              t1748 = t101 * tnorm_drho
              t1765 = t909 * tnorm_drho
              t1766 = t494 * trhoa
              t1767 = t303 * trhoanorm_drho
              t1801 = 0.2e1_dp * t1748 * t299 + 0.4e1_dp * t310 * t494 * t119 * &
                   trhoa - 0.2e1_dp * t297 * t313 * trhoa * t502 + 0.2e1_dp * t297 * &
                   t298 * trhoanorm_drho + 0.2e1_dp * t297 * t904 * tnorm_drho + t101 *&
                   t111 * (0.2e1_dp * t1765 + 0.2e1_dp * t1766 + 0.2e1_dp * t1767) * &
                   t119 - t310 * t924 * t502 - 0.2e1_dp * t297 * t313 * t321 * &
                   tnorm_drho - 0.2e1_dp * t493 * t494 * t312 * t321 + 0.2e1_dp * t310 &
                   * t936 * t321 * t502 - t310 * t313 * (0.2e1_dp * t1765 + 0.2e1_dp * &
                   t1766 + 0.2e1_dp * t1767 + 0.8e1_dp * t944 * Arhoa * tnorm_drho + &
                   0.12e2_dp * t953 * trhoa * tnorm_drho + 0.4e1_dp * t318 * &
                   trhoanorm_drho)

              e_ra_ndr(ii) = e_ra_ndr(ii)+&
                   scale_ec * (Hnorm_drho + my_rho * (0.3e1_dp * &
                   t293 * t506 * phirhoa + t110 * t1801 * t325 - t110 * t967 * t505))

              trhobnorm_drho = -t1739 * t1520 / 0.2e1_dp - t1741 * t1529 / &
                   0.2e1_dp - t1743 / 0.2e1_dp
              t1829 = t1286 * tnorm_drho
              t1830 = t494 * trhob
              t1831 = t303 * trhobnorm_drho
              t1865 = 0.2e1_dp * t1748 * t439 + 0.4e1_dp * t310 * t494 * t119 * &
                   trhob - 0.2e1_dp * t297 * t313 * trhob * t502 + 0.2e1_dp * t297 * &
                   t298 * trhobnorm_drho + 0.2e1_dp * t297 * t1269 * tnorm_drho + t101 &
                   * t111 * (0.2e1_dp * t1829 + 0.2e1_dp * t1830 + 0.2e1_dp * t1831) * &
                   t119 - t310 * t1304 * t502 - 0.2e1_dp * t297 * t313 * t453 * &
                   tnorm_drho - 0.2e1_dp * t493 * t494 * t312 * t453 + 0.2e1_dp * t310 &
                   * t936 * t453 * t502 - t310 * t313 * (0.2e1_dp * t1829 + 0.2e1_dp * &
                   t1830 + 0.2e1_dp * t1831 + 0.8e1_dp * t944 * Arhob * tnorm_drho + &
                   0.12e2_dp * t953 * trhob * tnorm_drho + 0.4e1_dp * t318 * &
                   trhobnorm_drho)

              e_rb_ndr(ii) = e_rb_ndr(ii)+&
                   scale_ec * (Hnorm_drho + my_rho * (0.3e1_dp * &
                   t293 * t506 * phirhob + t110 * t1865 * t325 - t110 * t1680 * t505))

              t1871 = tnorm_drho ** 2
              t1876 = A * t1871
              t1888 = t502 ** 2
              t1901 = t505 ** 2

              e_ndr_ndr(ii) = e_ndr_ndr(ii)+&
                   scale_ec * my_rho * (t110 * (0.2e1_dp * &
                   t101 * t1871 * t113 * t119 + 0.10e2_dp * t310 * t1876 * t119 - &
                   0.4e1_dp * t297 * t313 * tnorm_drho * t502 - 0.4e1_dp * t493 * t494 &
                   * t312 * t502 + 0.2e1_dp * t310 * t936 * t1888 - t310 * t313 * (&
                   0.2e1_dp * t1876 + 0.12e2_dp * t953 * t1871)) * t325 - t110 * t1901 &
                   * t966)

              e_ra_ndra(ii) = e_ra_ndra(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_a * &
                   Fx_anorm_drhoa + 0.2e1_dp * t350 * Fx_anorm_drhoa + 0.2e1_dp * t140 &
                   * (-0.8e1_dp * t1000 * t1001 * s_anorm_drhoa + 0.2e1_dp * t346 * &
                   s_anorm_drhoa * s_arhoa + 0.2e1_dp * t346 * s_a * (-t336 * t131 * &
                   kf_arhoa / 0.2e1_dp - t129 * t341 / 0.2e1_dp))) / 0.2e1_dp

              t1922 = s_anorm_drhoa ** 2

              e_ndra_ndra(ii) = e_ndra_ndra(ii)+&
                   scale_ex * t140 * (-0.8e1_dp * t999 * &
                   t133 * t1922 * t135 + 0.2e1_dp * t346 * t1922)

              e_rb_ndrb(ii) = e_rb_ndrb(ii)+&
                   scale_ex * (0.2e1_dp * ex_unif_b * &
                   Fx_bnorm_drhob + 0.2e1_dp * t481 * Fx_bnorm_drhob + 0.2e1_dp * t156 &
                   * (-0.8e1_dp * t1712 * t1713 * s_bnorm_drhob + 0.2e1_dp * t477 * &
                   s_bnorm_drhob * s_brhob + 0.2e1_dp * t477 * s_b * (-t467 * t148 * &
                   kf_brhob / 0.2e1_dp - t146 * t472 / 0.2e1_dp))) / 0.2e1_dp

              t1949 = s_bnorm_drhob ** 2
              e_ndrb_ndrb(ii) = e_ndrb_ndrb(ii)+&
                   scale_ex * t156 * (-0.8e1_dp * t1711 *&
                   t150 * t1949 * t135 + 0.2e1_dp * t477 * t1949)
           END IF
        END IF
     END DO
!$omp end do
  END SELECT
END SUBROUTINE pbe_lsd_calc

END MODULE xc_pbe
