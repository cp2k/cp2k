!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic code
! *****************************************************************************
MODULE atom_xc

  USE atom_types,                      ONLY: atom_type
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: xc_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_atom,                         ONLY: xc_rho_set_atom_update
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
  USE xc_derivatives,                  ONLY: xc_functionals_get_needs,&
                                             xc_functionals_eval
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type


#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_xc'

  PUBLIC :: 

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_atom_vxc(atom,atom_input,error)
    TYPE(atom_type), INTENT(INOUT)           :: atom
    TYPE(section_vals_type), POINTER         :: atom_input
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_vxc', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure=.FALSE.
    INTEGER                                  :: handle, myfun, nspins
    TYPE(section_vals_type), POINTER         :: xc_section, xc_fun_section
    TYPE(xc_rho_cflags_type), POINTER        :: needs
    REAL(KIND=dp)                            :: density_cut, gradient_cut, tau_cut
    REAL(KIND=dp), DIMENSION(2,3)            :: bounds
    LOGICAL                                  :: lsd, gradient_f
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(rho_set) 

    xc_section => section_vals_get_subs_vals(atom_input,"METHOD%XC",error=error)
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
    CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)

    IF(myfun == xc_none) THEN
!     energy%exc1 = 0.0_dp
    ELSE
      CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",r_val=density_cut,error=error)
      CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",r_val=gradient_cut,error=error)
      CALL section_vals_val_get(xc_section,"TAU_CUTOFF",r_val=tau_cut,error=error)

      ! hard coded for now
      lsd = .FALSE.
      nspins = 1
      needs = xc_functionals_get_needs(xc_fun_section,lsd=lsd,add_basic_components=.TRUE.,error=error)

      ! Prepare the structures needed to calculate and store the xc derivatives
  
      ! Array dimension: here anly one dimensional arrays are used,
      ! i.e. only the first column of deriv_data is read.
      ! The other to dimensions  are set to size equal 1
      bounds(1:2,1:3) = 1
      bounds(2,1) = atom%basis%grid%nr
  
      ! create a place where to put the derivatives
      CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
      ! create the place where to store the argument for the functionals
      CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=density_cut,&
                             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                             error=error)
      ! allocate the required 3d arrays where to store rho and drho
      CALL xc_rho_set_atom_update(rho_set,needs,nspins,bounds)

!         CALL fill_rho_set(rho_set,lsd,nspins,needs,rho,drho,nr,1,error=error)

      CALL xc_dset_zero_all(deriv_set, error)

      CALL xc_functionals_eval(xc_fun_section,lsd=lsd,rho_set=rho_set,deriv_set=deriv_set,&
                               deriv_order=deriv_order,error=error)

      ! Integration to get the matrix elements relative to the vxc_atom
      ! here the products with the primitives is done: gaVxcgb

!         IF (.NOT.energy_only) THEN
!           IF (gradient_f) THEN
!             CALL gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
!                             rho_atom,drho_h,drho_s,nspins,error=error)
!           ELSE
!             CALL gaVxcgb_noGC(vxc_h,vxc_s,atom_kind,rho_atom,nspins,error=error)
!           END IF
!         END IF  ! energy_only

      ! Release the xc structure used to store the xc derivatives
      CALL xc_dset_release(deriv_set, error=error)
      CALL xc_rho_set_release(rho_set,error=error)

    END IF !xc_none

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_atom_vxc

! *****************************************************************************

END MODULE atom_xc
