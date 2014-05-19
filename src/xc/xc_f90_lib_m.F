!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: libxc.f90 3550 2007-11-19 14:32:49Z marques $
!-------------------------------------------------------------------
MODULE xc_f90_lib_m

  USE libxc_funcs_m
  USE xc_f90_types_m,                  ONLY: xc_f90_kind,&
                                             xc_f90_pointer_t

  IMPLICIT NONE

  PUBLIC

  ! Families of xc functionals
  INTEGER, PARAMETER :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32

  INTEGER, PARAMETER :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized

  INTEGER, PARAMETER :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.

  ! Kinds
  INTEGER, PARAMETER :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2, &
    XC_KINETIC = 3

  INTEGER, PARAMETER :: &
    XC_FLAGS_HAVE_EXC = 1, &
    XC_FLAGS_HAVE_VXC = 2, &
    XC_FLAGS_HAVE_FXC = 4, &
    XC_FLAGS_HAVE_KXC = 8, &
    XC_FLAGS_HAVE_LXC = 16, &
    XC_FLAGS_1D = 32, &
    XC_FLAGS_2D = 64, &
    XC_FLAGS_3D = 128, &
    XC_FLAGS_STABLE = 512, &
    XC_FLAGS_DEVELOPMENT = 1024

  ! This value was redefined as XC_GGA_X_LB, we define it here to keep
  ! compatibility.
  INTEGER, PARAMETER :: XC_GGA_XC_LB = 160

  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_version(major, minor)
    INTEGER, INTENT(OUT)                     :: major, minor

    END SUBROUTINE xc_f90_version
  END INTERFACE

  !----------------------------------------------------------------
  INTERFACE
    INTEGER FUNCTION xc_f90_info_number(info)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
      TYPE(xc_f90_pointer_t), INTENT(in) :: info
    END FUNCTION xc_f90_info_number

    INTEGER FUNCTION xc_f90_info_kind(info)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
      TYPE(xc_f90_pointer_t), INTENT(in) :: info
    END FUNCTION xc_f90_info_kind

    SUBROUTINE xc_f90_info_name(info, s)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: info
    CHARACTER(len=*), INTENT(out)            :: s

    END SUBROUTINE xc_f90_info_name

    INTEGER FUNCTION xc_f90_info_family(info)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
      TYPE(xc_f90_pointer_t), INTENT(in) :: info
    END FUNCTION xc_f90_info_family

    INTEGER FUNCTION xc_f90_info_flags(info)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
      TYPE(xc_f90_pointer_t), INTENT(in) :: info
    END FUNCTION xc_f90_info_flags

    SUBROUTINE xc_f90_info_refs(info, number, str, s)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: info
    INTEGER, INTENT(inout)                   :: number
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: str
    CHARACTER(len=*), INTENT(out)            :: s

! number of the reference. Must be 0 in the first call
! this will hold a (char **) pointer
! the string that is output

    END SUBROUTINE xc_f90_info_refs

    INTEGER FUNCTION xc_f90_family_from_id(id)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
      INTEGER, INTENT(in) :: id
    END FUNCTION xc_f90_family_from_id
  END INTERFACE


  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_func_init(p, info, functional, nspin)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(out)      :: p, info
    INTEGER, INTENT(in)                      :: functional, nspin

    END SUBROUTINE xc_f90_func_init

    SUBROUTINE xc_f90_func_end(p)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p

    END SUBROUTINE xc_f90_func_end
  END INTERFACE


  ! LDAs
  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_lda(p, np, rho, zk, vrho, fxc, kxc)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: zk, vrho, fxc, kxc

! rho(nspin) the density
! the energy per unit particle
! v(nspin) the potential
! v(nspin,nspin) the xc kernel
! v(nspin,nspin,nspin) the derivative of xc kernel

    END SUBROUTINE xc_f90_lda

    SUBROUTINE xc_f90_lda_exc(p, np, rho, zk)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: zk

! rho(nspin) the density
! the energy per unit particle

    END SUBROUTINE xc_f90_lda_exc

    SUBROUTINE xc_f90_lda_exc_vxc(p, np, rho, e, v)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: e, v

! rho(nspin) the density
! the energy per unit particle
! v(nspin) the potential

    END SUBROUTINE xc_f90_lda_exc_vxc

    SUBROUTINE xc_f90_lda_vxc(p, np, rho, v)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: v

! rho(nspin) the density
! v(nspin) the potential

    END SUBROUTINE xc_f90_lda_vxc

    SUBROUTINE xc_f90_lda_fxc(p, np, rho, fxc)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: fxc

! rho(nspin) the density
! v(nspin,nspin) the xc kernel

    END SUBROUTINE xc_f90_lda_fxc

    SUBROUTINE xc_f90_lda_kxc(p, np, rho, kxc)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho
    REAL(xc_f90_kind), INTENT(out)           :: kxc

! rho(nspin) the density

    END SUBROUTINE xc_f90_lda_kxc
  END INTERFACE


  INTERFACE
    SUBROUTINE xc_f90_lda_x_1d_set_par(p, interaction, bb)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    INTEGER, INTENT(in)                      :: interaction
    REAL(xc_f90_kind), INTENT(in)            :: bb

    END SUBROUTINE xc_f90_lda_x_1d_set_par

    SUBROUTINE xc_f90_lda_c_xalpha_set_par(p, alpha)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    REAL(xc_f90_kind), INTENT(in)            :: alpha

    END SUBROUTINE xc_f90_lda_c_xalpha_set_par

    SUBROUTINE xc_f90_lda_x_set_par(p, alpha, relativistic, omega)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    REAL(xc_f90_kind), INTENT(in)            :: alpha
    INTEGER, INTENT(in)                      :: relativistic
    REAL(xc_f90_kind), INTENT(in)            :: omega

! of Xalpha, set to 4/3 to obtain standard LDA

    END SUBROUTINE xc_f90_lda_x_set_par

    SUBROUTINE xc_f90_lda_c_1d_csc_set_par(p, interaction, bb)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    INTEGER, INTENT(in)                      :: interaction
    REAL(xc_f90_kind), INTENT(in)            :: bb

    END SUBROUTINE xc_f90_lda_c_1d_csc_set_par

    SUBROUTINE xc_f90_lda_c_2d_prm_set_par(p, N)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    REAL(xc_f90_kind), INTENT(in)            :: N

    END SUBROUTINE xc_f90_lda_c_2d_prm_set_par
  END INTERFACE

  ! GGAs
  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_gga(p, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma
    REAL(xc_f90_kind), INTENT(out)           :: zk, vrho, vsigma, v2rho2, &
                                                v2rhosigma, v2sigma2

    END SUBROUTINE xc_f90_gga

    SUBROUTINE xc_f90_gga_exc(p, np, rho, sigma, zk)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma
    REAL(xc_f90_kind), INTENT(out)           :: zk

    END SUBROUTINE xc_f90_gga_exc

    SUBROUTINE xc_f90_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma
    REAL(xc_f90_kind), INTENT(out)           :: zk, vrho, vsigma

    END SUBROUTINE xc_f90_gga_exc_vxc

    SUBROUTINE xc_f90_gga_vxc(p, np, rho, sigma, vrho, vsigma)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma
    REAL(xc_f90_kind), INTENT(out)           :: vrho, vsigma

    END SUBROUTINE xc_f90_gga_vxc

    SUBROUTINE xc_f90_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma
    REAL(xc_f90_kind), INTENT(out)           :: v2rho2, v2rhosigma, v2sigma2

    END SUBROUTINE xc_f90_gga_fxc
  END INTERFACE

  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_gga_lb_set_par(p, modified, threshold, ip, qtot)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: modified
    REAL(xc_f90_kind), INTENT(in)            :: threshold, ip, qtot

! should we use the modified version
! if so, the threshold to use the asymtotic version
! ionization potential
! total charge

    END SUBROUTINE xc_f90_gga_lb_set_par
  END INTERFACE


  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_gga_lb_modified(p, np, rho, grho, r, dedd)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, grho, r
    REAL(xc_f90_kind), INTENT(out)           :: dedd

! rho(nspin) the density
! grho(3,nspin) the gradient of the density
! distance from center of finite system

    END SUBROUTINE xc_f90_gga_lb_modified
  END INTERFACE


  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_hyb_exx_coef(p, coef)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    REAL(xc_f90_kind), INTENT(out)           :: coef

    END SUBROUTINE xc_f90_hyb_exx_coef
  END INTERFACE


  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_gga_x_wpbeh_set_par(p, omega)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    REAL(xc_f90_kind), INTENT(in)            :: omega

! range separation

    END SUBROUTINE xc_f90_gga_x_wpbeh_set_par
  END INTERFACE

  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_gga_x_hjs_set_par(p, omega)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    REAL(xc_f90_kind), INTENT(in)            :: omega

! range separation

    END SUBROUTINE xc_f90_gga_x_hjs_set_par
  END INTERFACE

  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_hyb_gga_xc_hse_set_par(p, omega)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    REAL(xc_f90_kind), INTENT(in)            :: omega

! range separation

    END SUBROUTINE xc_f90_hyb_gga_xc_hse_set_par
  END INTERFACE


  ! the meta-GGAs
  !----------------------------------------------------------------
  INTERFACE
    SUBROUTINE xc_f90_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma, lapl, tau
    REAL(xc_f90_kind), INTENT(out) :: zk, vrho, vsigma, vlapl, vtau, v2rho2, &
      v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau

    END SUBROUTINE xc_f90_mgga

    SUBROUTINE xc_f90_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma, lapl, tau
    REAL(xc_f90_kind), INTENT(out)           :: zk

    END SUBROUTINE xc_f90_mgga_exc

    SUBROUTINE xc_f90_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma, lapl, tau
    REAL(xc_f90_kind), INTENT(out)           :: zk, vrho, vsigma, vlapl, vtau

    END SUBROUTINE xc_f90_mgga_exc_vxc

    SUBROUTINE xc_f90_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma, lapl, tau
    REAL(xc_f90_kind), INTENT(out)           :: vrho, vsigma, vlapl, vtau

    END SUBROUTINE xc_f90_mgga_vxc

    SUBROUTINE xc_f90_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(in)       :: p
    INTEGER, INTENT(in)                      :: np
    REAL(xc_f90_kind), INTENT(in)            :: rho, sigma, lapl, tau
    REAL(xc_f90_kind), INTENT(out) :: v2rho2, v2sigma2, v2lapl2, v2tau2, &
      v2rhosigma, v2rholapl, v2rhotau, v2sigmalapl, v2sigmatau, v2lapltau

    END SUBROUTINE xc_f90_mgga_fxc
  END INTERFACE

  INTERFACE
    SUBROUTINE xc_f90_mgga_x_tb09_set_par(p, cc)
      USE xc_f90_types_m, ONLY: xc_f90_pointer_t, xc_f90_kind
    TYPE(xc_f90_pointer_t), INTENT(inout)    :: p
    REAL(xc_f90_kind), INTENT(in)            :: cc

    END SUBROUTINE xc_f90_mgga_x_tb09_set_par
  END INTERFACE


END MODULE xc_f90_lib_m
