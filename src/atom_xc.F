!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic code
! *****************************************************************************
MODULE atom_xc

  USE atom_types,                      ONLY: atom_type,&
                                             opmat_type
  USE atom_utils,                      ONLY: atom_density,&
                                             integrate_grid,&
                                             numpot_matrix
  USE input_constants,                 ONLY: xc_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_atom,                         ONLY: xc_rho_set_atom_update
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                             xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_xc'

  PUBLIC :: calculate_atom_vxc

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_atom_vxc(xcmat,atom,method_section,error)
    TYPE(opmat_type), POINTER                :: xcmat
    TYPE(atom_type), INTENT(INOUT)           :: atom
    TYPE(section_vals_type), POINTER         :: method_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_vxc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deriv_order, handle, ierr, l, &
                                                myfun, n1, n2, n3, nr, nspins
    INTEGER, DIMENSION(2, 3)                 :: bounds
    LOGICAL                                  :: failure = .FALSE., lsd
    REAL(KIND=dp)                            :: density_cut, gradient_cut, &
                                                tau_cut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: drho, rho, tau
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: taumat, xcpot
    TYPE(section_vals_type), POINTER         :: xc_fun_section, xc_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    NULLIFY(rho_set) 

    xc_section => section_vals_get_subs_vals(method_section,"XC",error=error)
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
    CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)

    IF(myfun == xc_none) THEN
      atom%energy%exc = 0._dp
    ELSE
      CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",r_val=density_cut,error=error)
      CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",r_val=gradient_cut,error=error)
      CALL section_vals_val_get(xc_section,"TAU_CUTOFF",r_val=tau_cut,error=error)

      ! hard coded for now
      lsd = .FALSE.
      nspins = 1
      needs = xc_functionals_get_needs(xc_fun_section,lsd=lsd,add_basic_components=.FALSE.,error=error)

      ! Prepare the structures needed to calculate and store the xc derivatives
  
      ! Array dimension: here anly one dimensional arrays are used,
      ! i.e. only the first column of deriv_data is read.
      ! The other to dimensions  are set to size equal 1
      nr = atom%basis%grid%nr
      bounds(1:2,1:3) = 1
      bounds(2,1) = nr
  
      ! create a place where to put the derivatives
      NULLIFY(deriv_set)
      CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
      ! create the place where to store the argument for the functionals
      CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=density_cut,&
             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,error=error)
      ! allocate the required 3d arrays where to store rho and drho
      CALL xc_rho_set_atom_update(rho_set,needs,nspins,bounds)

      NULLIFY(rho,drho,tau)
      IF ( needs%rho ) THEN
        ALLOCATE(rho(nr),STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        CALL atom_density(rho,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
      END IF
      IF ( needs%norm_drho ) THEN
        ALLOCATE(drho(nr),STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        CALL atom_density(drho,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="DER",error=error)
      END IF
      IF ( needs%tau ) THEN
        ALLOCATE(tau(nr),STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        CALL atom_density(tau,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,&
                          typ="KIN",rr=atom%basis%grid%rad2,error=error)
      END IF

      CALL fill_rho_set(rho_set,nspins,needs,rho,drho,tau,nr,error=error)

      CALL xc_dset_zero_all(deriv_set, error)

      deriv_order = 1
      CALL xc_functionals_eval(xc_fun_section,lsd=lsd,rho_set=rho_set,deriv_set=deriv_set,&
                               deriv_order=deriv_order,error=error)

      ! Integration to get the matrix elements and energy
      deriv => xc_dset_get_derivative(deriv_set,"",allocate_deriv=.FALSE., error=error)
      CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
      atom%energy%exc = fourpi*integrate_grid(xcpot(:,1,1),atom%basis%grid)

      IF ( needs%rho ) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(rho)",allocate_deriv=.FALSE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
        CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,0,error)
        DEALLOCATE(rho,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF
      IF ( needs%norm_drho ) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",allocate_deriv=.FALSE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
        CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,1,error)
        DEALLOCATE(drho,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF
      IF ( needs%tau ) THEN
        deriv => xc_dset_get_derivative(deriv_set,"(tau)",allocate_deriv=.FALSE.,error=error)
        CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
        n1 = SIZE(xcmat%op,1)
        n2 = SIZE(xcmat%op,2)
        n3 = SIZE(xcmat%op,3)
        ALLOCATE(taumat(n1,n2,0:n3-1),STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        taumat = 0._dp

        xcpot(:,1,1) = 0.5_dp * xcpot(:,1,1)
        CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,2,error)
        xcpot(:,1,1) = xcpot(:,1,1)/atom%basis%grid%rad2(:)
        CALL numpot_matrix(taumat,xcpot(:,1,1),atom%basis,0,error)
        DO l=0,3
          xcmat%op(:,:,l) = xcmat%op(:,:,l) + REAL(l*(l+1),dp)*taumat(:,:,l)
        END DO

        DEALLOCATE(tau,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        DEALLOCATE(taumat,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF

      ! Release the xc structure used to store the xc derivatives
      CALL xc_dset_release(deriv_set, error=error)
      CALL xc_rho_set_release(rho_set,error=error)

    END IF !xc_none

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_vxc

! *****************************************************************************
  SUBROUTINE fill_rho_set(rho_set,nspins,needs,rho,drho,tau,na,error)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(xc_rho_cflags_type), INTENT(in)     :: needs
    REAL(dp), DIMENSION(:), POINTER          :: rho, drho, tau
    INTEGER, INTENT(IN)                      :: na
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_rho_set', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f13 = (1.0_dp/3.0_dp)

    INTEGER                                  :: ia
    LOGICAL                                  :: failure

    failure = .FALSE.

    SELECT CASE(nspins)
    CASE(1)
       CPPrecondition(.NOT.needs%rho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%norm_drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin_1_3,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%tau_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho,cp_failure_level,routineP,error,failure)
       ! Give rho to 1/3
       IF (needs%rho_1_3) THEN
          DO ia = 1,na
             rho_set%rho_1_3(ia,1,1) = MAX(rho(ia),0.0_dp)**f13
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       ! Give the density
       IF (needs%rho) THEN
          DO ia = 1,na
             rho_set%rho(ia,1,1) = rho(ia)
          END DO
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       ! Give the norm of the gradient of the density
       IF (needs%norm_drho) THEN
          DO ia = 1,na
             rho_set%norm_drho(ia,1,1) = drho(ia)
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
    CASE(2)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !
    END SELECT

    ! tau part
    IF (needs%tau) THEN
      DO ia = 1,na
         rho_set%tau(ia,1,1) = tau(ia)
      END DO
      rho_set%owns%tau=.TRUE.
      rho_set%has%tau=.TRUE.
    END IF

  END SUBROUTINE fill_rho_set
! *****************************************************************************

END MODULE atom_xc
