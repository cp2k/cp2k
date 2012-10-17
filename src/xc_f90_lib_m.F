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

  use xc_f90_types_m
  use libxc_funcs_m

  implicit none

  public

  ! Families of xc functionals
  integer, parameter :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32

  integer, parameter :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized

  integer, parameter :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.

  ! Kinds
  integer, parameter :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2, &
    XC_KINETIC = 3

  integer, parameter :: &
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
  integer, parameter :: XC_GGA_XC_LB = 160

  !----------------------------------------------------------------
  interface
    integer function xc_f90_info_number(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_number

    integer function xc_f90_info_kind(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_family

    integer function xc_f90_info_flags(info)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
    end function xc_f90_info_flags

    subroutine xc_f90_info_refs(info, number, str, s)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: info
      integer, intent(inout) :: number ! number of the reference. Must be 0 in the first call
      type(xc_f90_pointer_t), intent(inout) :: str ! this will hold a (char **) pointer
      character(len=*), intent(out) :: s ! the string that is output
    end subroutine xc_f90_info_refs

    integer function xc_f90_family_from_id(id)
      use xc_f90_types_m
      integer, intent(in) :: id
    end function xc_f90_family_from_id
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_func_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(out) :: p
      type(xc_f90_pointer_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_func_init

    subroutine xc_f90_func_end(p)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
    end subroutine xc_f90_func_end
  end interface


  ! LDAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lda(p, np, rho, zk, vrho, fxc, kxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: zk ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: vrho ! v(nspin) the potential
      real(xc_f90_kind), intent(out) :: fxc ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind), intent(out) :: kxc ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda

    subroutine xc_f90_lda_exc(p, np, rho, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: zk ! the energy per unit particle
    end subroutine xc_f90_lda_exc

    subroutine xc_f90_lda_exc_vxc(p, np, rho, e, v)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: v ! v(nspin) the potential
    end subroutine xc_f90_lda_exc_vxc

    subroutine xc_f90_lda_vxc(p, np, rho, v)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: v ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc

    subroutine xc_f90_lda_fxc(p, np, rho, fxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: fxc ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc

    subroutine xc_f90_lda_kxc(p, np, rho, kxc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: kxc
    end subroutine xc_f90_lda_kxc
  end interface


  interface
    subroutine xc_f90_lda_x_1d_set_par(p, interaction, bb)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      integer, intent(in) :: interaction
      real(xc_f90_kind), intent(in) :: bb
    end subroutine xc_f90_lda_x_1d_set_par

    subroutine xc_f90_lda_c_xalpha_set_par(p, alpha)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: alpha
    end subroutine xc_f90_lda_c_xalpha_set_par

    subroutine xc_f90_lda_x_set_par(p, relativistic)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      integer, intent(in) :: relativistic
    end subroutine xc_f90_lda_x_set_par

    subroutine xc_f90_lda_c_1d_csc_set_par(p, interaction, bb)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      integer, intent(in) :: interaction
      real(xc_f90_kind), intent(in) :: bb
    end subroutine xc_f90_lda_c_1d_csc_set_par

    subroutine xc_f90_lda_c_2d_prm_set_par(p, N)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: N
    end subroutine xc_f90_lda_c_2d_prm_set_par
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga(p, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
    end subroutine xc_f90_gga

    subroutine xc_f90_gga_exc(p, np, rho, sigma, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
    end subroutine xc_f90_gga_exc

    subroutine xc_f90_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
    end subroutine xc_f90_gga_exc_vxc

    subroutine xc_f90_gga_vxc(p, np, rho, sigma, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
    end subroutine xc_f90_gga_vxc

    subroutine xc_f90_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
    end subroutine xc_f90_gga_fxc
  end interface

  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_lb_set_par(p, modified, threshold, ip, qtot)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: modified ! should we use the modified version
      real(xc_f90_kind), intent(in) :: threshold ! if so, the threshold to use the asymtotic version
      real(xc_f90_kind), intent(in) :: ip ! ionization potential
      real(xc_f90_kind), intent(in) :: qtot ! total charge
    end subroutine xc_f90_gga_lb_set_par
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_lb_modified(p, np, rho, grho, r, dedd)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: grho ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(in) :: r ! distance from center of finite system
      real(xc_f90_kind), intent(out) :: dedd
    end subroutine xc_f90_gga_lb_modified
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_gga_exx_coef(p, coef)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(out) :: coef
    end subroutine xc_f90_hyb_gga_exx_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: vlapl
      real(xc_f90_kind), intent(out) :: vtau
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2sigma2
      real(xc_f90_kind), intent(out) :: v2lapl2
      real(xc_f90_kind), intent(out) :: v2tau2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2rholapl
      real(xc_f90_kind), intent(out) :: v2rhotau
      real(xc_f90_kind), intent(out) :: v2sigmalapl
      real(xc_f90_kind), intent(out) :: v2sigmatau
      real(xc_f90_kind), intent(out) :: v2lapltau
    end subroutine xc_f90_mgga

    subroutine xc_f90_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
    end subroutine xc_f90_mgga_exc

    subroutine xc_f90_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: vlapl
      real(xc_f90_kind), intent(out) :: vtau
    end subroutine xc_f90_mgga_exc_vxc

    subroutine xc_f90_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: vlapl
      real(xc_f90_kind), intent(out) :: vtau
    end subroutine xc_f90_mgga_vxc

    subroutine xc_f90_mgga_fxc(p, np, rho, sigma, lapl, tau, &
      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, &
      v2sigmalapl, v2sigmatau, v2lapltau)

      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      integer, intent(in) :: np
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: lapl
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2sigma2
      real(xc_f90_kind), intent(out) :: v2lapl2
      real(xc_f90_kind), intent(out) :: v2tau2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2rholapl
      real(xc_f90_kind), intent(out) :: v2rhotau
      real(xc_f90_kind), intent(out) :: v2sigmalapl
      real(xc_f90_kind), intent(out) :: v2sigmatau
      real(xc_f90_kind), intent(out) :: v2lapltau
    end subroutine xc_f90_mgga_fxc
  end interface

  interface
    subroutine xc_f90_mgga_x_tb09_set_par(p, cc)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: cc
    end subroutine xc_f90_mgga_x_tb09_set_par
  end interface


  ! the LCAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lca(p, rho, v, e, dedd, dedv)
      use xc_f90_types_m
      type(xc_f90_pointer_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: v ! v(3,nspin) the vorticity
      real(xc_f90_kind), intent(out) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: dedd ! dedd(nspin) the derivative of the energy
                                              ! in terms of the density
      real(xc_f90_kind), intent(out) :: dedv ! and in terms of the vorticity
    end subroutine xc_f90_lca
  end interface

END MODULE xc_f90_lib_m
