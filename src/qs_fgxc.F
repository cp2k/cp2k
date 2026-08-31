!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  https://en.wikipedia.org/wiki/Finite_difference_coefficient
!---------------------------------------------------------------------------------------------------
!Derivative    Accuracy            4       3       2       1       0       1       2      3       4
!---------------------------------------------------------------------------------------------------
!    1             2                                    -1/2       0     1/2
!                  4                            1/12    -2/3       0     2/3   -1/12
!                  6                   -1/60    3/20    -3/4       0     3/4   -3/20   1/60
!                  8           1/280  -4/105     1/5    -4/5       0     4/5    -1/5  4/105  -1/280
!---------------------------------------------------------------------------------------------------
!    2             2                                       1      -2       1
!                  4                           -1/12     4/3    -5/2     4/3   -1/12
!                  6                    1/90   -3/20     3/2  -49/18     3/2   -3/20   1/90
!                  8          -1/560   8/315    -1/5     8/5 -205/72     8/5    -1/5  8/315  -1/560
!---------------------------------------------------------------------------------------------------
!> \par History
!>     init 17.03.2020
!>     complete refactoring 08.2026
!> \author JGH
! **************************************************************************************************
MODULE qs_fgxc

   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc,                          ONLY: qs_fxc_create,&
                                              qs_fxc_fdiff
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_methods,                  ONLY: qs_rho_copy,&
                                              qs_rho_scale_and_add,&
                                              qs_rho_scale_and_add_b,&
                                              qs_rho_transfer
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE qs_vxc_atom,                     ONLY: fgxc_atom_calc,&
                                              fgxc_atom_diff
   USE xc,                              ONLY: xc_calc_2nd_deriv_analytical,&
                                              xc_calc_3rd_deriv_analytical,&
                                              xc_prep_3rd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: qs_fgxc_create, qs_fgxc_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fgxc'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param rho0_atom_set ...
!> \param xc_section ...
!> \param do_onecenter ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
!> \param rhof_atom_set ...
!> \param rhog_atom_set ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_create(qs_env, rho0_struct, rho1_struct, rho0_atom_set, &
                             xc_section, do_onecenter, is_triplet, spinflip, &
                             fxc_rho, fxc_tau, gxc_rho, gxc_tau, &
                             rhof_atom_set, rhog_atom_set, qs_kind_set)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: do_onecenter, is_triplet, spinflip
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau, gxc_rho, gxc_tau
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rhof_atom_set, rhog_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fgxc_create'

      INTEGER                                            :: handle, ispin, nspins, nsteps, order
      LOGICAL                                            :: analytic_2nd_deriv, analytic_3rd_deriv, &
                                                            has_rho, has_tau, uf_grid
      REAL(KIND=dp)                                      :: eps_delta, factor
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, uf_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho_uf, fxc_tau_uf, gxc_rho_uf, &
                                                            gxc_tau_uf, rho_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, weights, weights_uf
      TYPE(qs_rho_type), POINTER                         :: rho0_uf, rho1_uf

      CALL timeset(routineN, handle)

      ! Only rho0 and rho1 should be associated
      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(.NOT. ASSOCIATED(gxc_rho))
      CPASSERT(.NOT. ASSOCIATED(gxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      nsteps = section_get_ival(xc_section, "NSTEPS")
      order = 2*nsteps
      eps_delta = section_get_rval(xc_section, "STEP_SIZE")
      analytic_2nd_deriv = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")
      analytic_3rd_deriv = section_get_lval(xc_section, "3RD_DERIV_ANALYTICAL")
      IF (analytic_3rd_deriv .AND. .NOT. analytic_2nd_deriv) THEN
         CALL cp_warn(__LOCATION__, "XC Analytic 3rd derivatives cannot be combined with "// &
                      "numeric 2nd derivatives.")
         CPABORT("XC derivative options")
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, xc_pw_pool=xc_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

      NULLIFY (fxc_rho_uf, fxc_tau_uf, gxc_rho_uf, gxc_tau_uf)
      IF (analytic_2nd_deriv .AND. analytic_3rd_deriv) THEN

         ! NLCC
         NULLIFY (rho_nlcc, rho_nlcc_g)
         CALL get_qs_env(qs_env, rho_nlcc=rho_nlcc, rho_nlcc_g=rho_nlcc_g)
         IF (ASSOCIATED(rho_nlcc)) THEN
            NULLIFY (rho_r, rho_g)
            CALL qs_rho_get(rho0_struct, rho_r=rho_r, rho_g=rho_g)
            factor = 1.0_dp
            DO ispin = 1, nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF

         ! Integration weights
         NULLIFY (weights)
         CALL get_qs_env(qs_env, xcint_weights=weights)

         NULLIFY (weights_uf)
         IF (uf_grid) THEN
            IF (ASSOCIATED(weights)) THEN
               ALLOCATE (weights_uf)
               CALL xc_pw_pool%create_pw(weights_uf)
               BLOCK
                  TYPE(pw_c1d_gs_type) :: weights_g, weights_g_uf
                  CALL auxbas_pw_pool%create_pw(weights_g)
                  CALL xc_pw_pool%create_pw(weights_g_uf)
                  CALL pw_transfer(weights, weights_g)
                  CALL pw_transfer(weights_g, weights_g_uf)
                  CALL pw_transfer(weights_g_uf, weights_uf)
                  CALL xc_pw_pool%give_back_pw(weights_g_uf)
                  CALL auxbas_pw_pool%give_back_pw(weights_g)
               END BLOCK
            END IF
            !
            ALLOCATE (rho0_uf, rho1_uf)
            CALL qs_rho_create(rho0_uf)
            CALL qs_rho_create(rho1_uf)
            CALL qs_rho_transfer(rho0_struct, rho0_uf, auxbas_pw_pool, xc_pw_pool)
            CALL qs_rho_transfer(rho1_struct, rho1_uf, auxbas_pw_pool, xc_pw_pool)
            uf_pw_pool => xc_pw_pool
         ELSE
            weights_uf => weights
            rho0_uf => rho0_struct
            rho1_uf => rho1_struct
            uf_pw_pool => auxbas_pw_pool
         END IF

         ! only spinflip is activated !?
         IF (.NOT. spinflip) THEN
            CPABORT("Analytic 3rd EXC derivatives not available")
         END IF
         CALL qs_fgxc_analytic(rho0_uf, rho1_uf, xc_section, uf_pw_pool, weights_uf, &
                               is_triplet, spinflip, &
                               fxc_rho_uf, fxc_tau_uf, gxc_rho_uf, gxc_tau_uf)

         ! check for results
         has_rho = .FALSE.
         IF (ASSOCIATED(fxc_rho_uf) .AND. ASSOCIATED(gxc_rho_uf)) has_rho = .TRUE.
         has_tau = .FALSE.
         IF (ASSOCIATED(fxc_tau_uf) .AND. ASSOCIATED(gxc_tau_uf)) has_tau = .TRUE.
         ! transfer to standard grid
         IF (uf_grid) THEN
            ! rho
            IF (has_rho) THEN
               ALLOCATE (fxc_rho(nspins), gxc_rho(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(fxc_rho(ispin))
                  CALL auxbas_pw_pool%create_pw(gxc_rho(ispin))
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: xc_g, xc_g_uf
                     CALL auxbas_pw_pool%create_pw(xc_g)
                     CALL xc_pw_pool%create_pw(xc_g_uf)
                     CALL pw_transfer(fxc_rho_uf(ispin), xc_g_uf)
                     CALL pw_transfer(xc_g_uf, xc_g)
                     CALL pw_transfer(xc_g, fxc_rho(ispin))
                     CALL pw_transfer(gxc_rho_uf(ispin), xc_g_uf)
                     CALL pw_transfer(xc_g_uf, xc_g)
                     CALL pw_transfer(xc_g, gxc_rho(ispin))
                     CALL xc_pw_pool%give_back_pw(xc_g_uf)
                     CALL auxbas_pw_pool%give_back_pw(xc_g)
                  END BLOCK
                  CALL xc_pw_pool%give_back_pw(fxc_rho_uf(ispin))
                  CALL xc_pw_pool%give_back_pw(gxc_rho_uf(ispin))
               END DO
               DEALLOCATE (fxc_rho_uf, gxc_rho_uf)
            ELSE
               NULLIFY (fxc_rho, gxc_rho)
            END IF
            ! tau
            IF (has_tau) THEN
               ALLOCATE (fxc_tau(nspins), gxc_tau(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(fxc_tau(ispin))
                  CALL auxbas_pw_pool%create_pw(gxc_tau(ispin))
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: xc_g, xc_g_uf
                     CALL auxbas_pw_pool%create_pw(xc_g)
                     CALL xc_pw_pool%create_pw(xc_g_uf)
                     CALL pw_transfer(fxc_tau_uf(ispin), xc_g_uf)
                     CALL pw_transfer(xc_g_uf, xc_g)
                     CALL pw_transfer(xc_g, fxc_tau(ispin))
                     CALL pw_transfer(gxc_tau_uf(ispin), xc_g_uf)
                     CALL pw_transfer(xc_g_uf, xc_g)
                     CALL pw_transfer(xc_g, gxc_tau(ispin))
                     CALL xc_pw_pool%give_back_pw(xc_g_uf)
                     CALL auxbas_pw_pool%give_back_pw(xc_g)
                  END BLOCK
                  CALL xc_pw_pool%give_back_pw(fxc_tau_uf(ispin))
                  CALL xc_pw_pool%give_back_pw(gxc_tau_uf(ispin))
               END DO
               DEALLOCATE (fxc_tau_uf, gxc_tau_uf)
            ELSE
               NULLIFY (fxc_tau, gxc_tau)
            END IF
         ELSE
            fxc_rho => fxc_rho_uf
            fxc_tau => fxc_tau_uf
            gxc_rho => gxc_rho_uf
            gxc_tau => gxc_tau_uf
         END IF

         ! de-apply NLCC density
         IF (ASSOCIATED(rho_nlcc)) THEN
            factor = -1.0_dp
            DO ispin = 1, nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF

      ELSE IF (analytic_2nd_deriv .AND. .NOT. analytic_3rd_deriv) THEN

         NULLIFY (weights)
         CALL get_qs_env(qs_env, xcint_weights=weights)
         CALL qs_fgxc_gdiff(qs_env, rho0_struct, rho1_struct, xc_section, &
                            order, eps_delta, is_triplet, spinflip, &
                            fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      ELSE IF (.NOT. (analytic_2nd_deriv .OR. analytic_3rd_deriv)) THEN

         CALL qs_fgxc_fgdiff(qs_env, rho0_struct, rho1_struct, xc_section, &
                             order, is_triplet, spinflip, &
                             fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      ELSE
         CPABORT("EXC derivative option not available")
      END IF

      IF (do_onecenter) THEN
         IF (analytic_2nd_deriv .AND. analytic_3rd_deriv) THEN
            CPABORT("Analytic 3rd EXC derivatives not available")
         ELSE IF (analytic_2nd_deriv) THEN
            CALL fgxc_atom_diff(qs_env, rho0_atom_set, rhof_atom_set, rhog_atom_set, &
                                qs_kind_set, xc_section, is_triplet, order, eps_delta)
         ELSE IF (analytic_3rd_deriv) THEN
            CALL cp_warn(__LOCATION__, "XC Analytic 3rd derivatives cannot be combined with "// &
                         "numeric 2nd derivatives.")
            CPABORT("XC derivative options")
         ELSE
            CALL fgxc_atom_calc(qs_env, rho0_atom_set, rhof_atom_set, rhog_atom_set, &
                                qs_kind_set, xc_section, is_triplet, order)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fgxc_create

! **************************************************************************************************
!> \brief Calculates the values at the grid points  in real space (r_i), of the second and third
!>        functional derivatives of the exchange-correlation energy functional.
!>         fxc_rho(r_i) =  fxc[n](r_i)*n^(1)(r_i)            ! Second functional derivative
!>         gxc_rho(r_i) =  n^(1)(r_i)*gxc[n](r_i)*n^(1)(r_i) ! Third functional derivative
!> \param rho0_struct Ground state density, n(r).
!> \param rho1_struct Density used to fold the functional derivatives, n^(1)(r).
!> \param xc_section ...
!> \param pw_pool ...
!> \param weights ...
!> \param is_triplet Flag used to activate triplet kernel (not available)
!> \param spinflip Flag used to activate the spin-flip noncollinear kernel and kernel derivatives.
!> \param fxc_rho Second functional derivative with respect to the density, n(r).
!> \param fxc_tau mGGA contribution to the second functional derivative with respect to the density.
!> \param gxc_rho Third functional derivative with respect to the density, n(r).
!> \param gxc_tau mGGA contribution to the third functional derivative with respect to the density.
!> \par History
!>    * 07.2024 Created [LHS]
! **************************************************************************************************
   SUBROUTINE qs_fgxc_analytic(rho0_struct, rho1_struct, xc_section, pw_pool, weights, &
                               is_triplet, spinflip, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      LOGICAL, INTENT(IN)                                :: is_triplet, spinflip
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fgxc_analytic'

      INTEGER                                            :: handle, ispin, nspins, spindim
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: do_sf, lsd
      REAL(KIND=dp)                                      :: fac
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho0_g, rho1_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho0_r, rho1_r, tau0_r, tau1_r
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho0_set, rho1_set

      CALL timeset(routineN, handle)

      ! Only rho0 and rho1 should be associated
      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(.NOT. ASSOCIATED(gxc_rho))
      CPASSERT(.NOT. ASSOCIATED(gxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      ! only tested for spinflip
      CPASSERT(.NOT. is_triplet)
      CPASSERT(spinflip)

      ! Initialize parameters
      do_sf = spinflip
      !
      ! Get the values on the gridpoints of the rho0 density
      CALL qs_rho_get(rho0_struct, rho_r=rho0_r, rho_g=rho0_g, tau_r=tau0_r)
      nspins = SIZE(rho0_r)
      lsd = (nspins == 2)
      !
      IF (do_sf) THEN
         spindim = 1
      ELSE
         spindim = nspins
      END IF
      !
      fac = 0._dp
      IF (nspins == 1) THEN
         fac = 1.0_dp
      END IF

      ! Read xc functional section and find out what the functional actually needs
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! Create fields for the kernel and kernel derivative
      ALLOCATE (fxc_rho(spindim), gxc_rho(nspins))
      DO ispin = 1, spindim
         CALL pw_pool%create_pw(fxc_rho(ispin))
         CALL pw_zero(fxc_rho(ispin))
      END DO
      DO ispin = 1, nspins
         CALL pw_pool%create_pw(gxc_rho(ispin))
         CALL pw_zero(gxc_rho(ispin))
      END DO
      ! Create fields for mGGA functionals. This implementation is not ready yet!
      IF (needs%tau .OR. needs%tau_spin) THEN
         IF (.NOT. ASSOCIATED(tau1_r)) THEN
            CPABORT("Tau-dependent functionals requires allocated kinetic energy density grid")
         END IF
         ALLOCATE (fxc_tau(spindim), gxc_tau(nspins))
         DO ispin = 1, spindim
            CALL pw_pool%create_pw(fxc_tau(ispin))
            CALL pw_zero(fxc_tau(ispin))
         END DO
         DO ispin = 1, nspins
            CALL pw_pool%create_pw(gxc_tau(ispin))
            CALL pw_zero(gxc_tau(ispin))
         END DO
      END IF

      ! Build rho0_set
      ! calculate the arguments needed by the functionals
      ! Needs
      !   deriv_set  xc_derivative_set_type just declared
      !   rho0_set   xc_rho_set_type just declared
      !   rho0_r     pw_type calculated by qs_rho_get
      !   pw_pool    given by the calling subroutine
      !   xc_section given by the calling subroutine
      !   tau0_r     pw_type calculated by qs_rho_get
      CALL xc_prep_3rd_deriv(deriv_set, rho0_set, rho0_r, pw_pool, weights, &
                             xc_section, tau_r=tau0_r, do_sf=do_sf)

      ! Build rho1_set
      ! Get the values on the gridpoints of the rho1 density
      CALL qs_rho_get(rho1_struct, rho_r=rho1_r, rho_g=rho1_g, tau_r=tau1_r)
      bo = rho1_r(1)%pw_grid%bounds_local
      ! create the place where to store the argument for the functionals
      ! Needs
      !   rho1_set xc_rho_set_type just declared
      !   bo 2x3   integer matrix should have bounds_local or rho1_r
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))
      ! calculate the arguments needed by the functionals
      ! Needs
      !   rho1_set             object created by xc_rho_set_create
      !   rho1_r,rho1_g,tau1_r pw_type values of rho1 in real space grid
      !   needs                xc_rho_cflags_type defined through xc_functionals_get_needs
      !   pw_pool              Given by the calling subroutine
      CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             pw_pool, spinflip=do_sf)

      ! Calculate exchange correlation kernel
      ! Needs
      !   fxc_rho, fxc_tau pw_type not associated
      !   deriv_set        created and defined by xc_prep_3rd_deriv
      !   rho0_set         xc_rho_set_type build by xc_prep_3rd_deriv
      !   rho1_set         xc_rho_set_type build by xc_rho_set_create/update
      !   pw_pool          needs to be given by the calling subroutine
      !   xc_section       needs to be given by the calling subroutine
      CALL xc_calc_2nd_deriv_analytical(fxc_rho, fxc_tau, deriv_set, rho0_set, rho1_set, pw_pool, &
                                        xc_section, .FALSE., spinflip=do_sf, tddfpt_fac=fac)
      ! Calculate exchange correlation kernel derivative
      CALL xc_calc_3rd_deriv_analytical(gxc_rho, gxc_tau, deriv_set, rho0_set, rho1_set, pw_pool, &
                                        xc_section, spinflip=do_sf)

      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho0_set)
      CALL xc_rho_set_release(rho1_set)

      CALL timestop(handle)

   END SUBROUTINE qs_fgxc_analytic

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param accuracy ...
!> \param epsrho ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_gdiff(qs_env, rho0_struct, rho1_struct, xc_section, accuracy, epsrho, &
                            is_triplet, spinflip, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      INTEGER, INTENT(IN)                                :: accuracy
      REAL(KIND=dp), INTENT(IN)                          :: epsrho
      LOGICAL, INTENT(IN)                                :: is_triplet, spinflip
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fgxc_gdiff'

      INTEGER                                            :: handle, ispin, istep, nspins, nstep
      LOGICAL                                            :: do_sf
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_tau_rspace, vxc00, vxc00b
      TYPE(qs_rho_type), POINTER                         :: rhoin
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set, rho1_atom_set

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(.NOT. ASSOCIATED(gxc_rho))
      CPASSERT(.NOT. ASSOCIATED(gxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      do_sf = spinflip

      ak = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
      END SELECT

      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = dft_control%nspins
      exc = 0.0_dp

      IF (do_sf) THEN
         CALL qs_fxc_create(qs_env, rho0_struct, rho1_struct, rho0_atom_set, xc_section, &
                            .FALSE., fxc_rho, fxc_tau, rho1_atom_set, &
                            is_triplet=is_triplet, spinflip=do_sf)
      ELSE
         CALL qs_fxc_fdiff(qs_env, rho0_struct, rho1_struct, xc_section, &
                           accuracy, fxc_rho, fxc_tau, is_triplet)
      END IF

      DO istep = -nstep, nstep

         IF (ak(istep) /= 0.0_dp) THEN
            alpha = 1.0_dp
            beta = REAL(istep, KIND=dp)*epsrho
            NULLIFY (rhoin)
            ALLOCATE (rhoin)
            CALL qs_rho_create(rhoin)
            NULLIFY (vxc00, vxc00b, v_tau_rspace)
            CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, nspins)
            CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, beta)
            IF (do_sf) THEN
               ! variation in alpha density
               CALL qs_fxc_create(qs_env, rhoin, rho1_struct, rho0_atom_set, xc_section, &
                                  .FALSE., vxc00, v_tau_rspace, rho1_atom_set, &
                                  is_triplet=is_triplet, spinflip=do_sf)
               ! variation in beta density
               CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, nspins)
               CALL qs_rho_scale_and_add_b(rhoin, rho1_struct, alpha, beta)
               CALL qs_fxc_create(qs_env, rhoin, rho1_struct, rho0_atom_set, xc_section, &
                                  .FALSE., vxc00b, v_tau_rspace, rho1_atom_set, &
                                  is_triplet=is_triplet, spinflip=do_sf)
            ELSE
               CALL qs_fxc_fdiff(qs_env, rhoin, rho1_struct, xc_section, accuracy, &
                                 fxc_rho=vxc00, fxc_tau=v_tau_rspace, is_triplet=is_triplet)
            END IF
            CALL qs_rho_release(rhoin)
            DEALLOCATE (rhoin)
            IF (.NOT. ASSOCIATED(gxc_rho)) THEN
               ALLOCATE (gxc_rho(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(gxc_rho(ispin))
                  CALL pw_zero(gxc_rho(ispin))
               END DO
            END IF
            IF (do_sf) THEN
               CALL pw_axpy(vxc00(1), gxc_rho(1), ak(istep))
               CALL pw_axpy(vxc00b(1), gxc_rho(2), ak(istep))
            ELSE
               DO ispin = 1, nspins
                  CALL pw_axpy(vxc00(ispin), gxc_rho(ispin), ak(istep))
               END DO
            END IF
            DO ispin = 1, SIZE(vxc00)
               CALL auxbas_pw_pool%give_back_pw(vxc00(ispin))
            END DO
            DEALLOCATE (vxc00)
            IF (ASSOCIATED(vxc00b)) THEN
               CALL auxbas_pw_pool%give_back_pw(vxc00b(1))
               DEALLOCATE (vxc00b)
            END IF
            IF (ASSOCIATED(v_tau_rspace)) THEN
               IF (.NOT. ASSOCIATED(gxc_tau)) THEN
                  ALLOCATE (gxc_tau(nspins))
                  DO ispin = 1, nspins
                     CALL auxbas_pw_pool%create_pw(gxc_tau(ispin))
                     CALL pw_zero(gxc_tau(ispin))
                  END DO
               END IF
               DO ispin = 1, nspins
                  CALL pw_axpy(v_tau_rspace(ispin), gxc_tau(ispin), ak(istep))
               END DO
               DO ispin = 1, SIZE(v_tau_rspace)
                  CALL auxbas_pw_pool%give_back_pw(v_tau_rspace(ispin))
               END DO
               DEALLOCATE (v_tau_rspace)
            END IF
         END IF

      END DO

      oeps1 = 1.0_dp/epsrho
      DO ispin = 1, nspins
         CALL pw_scale(gxc_rho(ispin), oeps1)
      END DO
      IF (ASSOCIATED(gxc_tau)) THEN
         DO ispin = 1, nspins
            CALL pw_scale(gxc_tau(ispin), oeps1)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fgxc_gdiff

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param accuracy ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_fgdiff(qs_env, rho0_struct, rho1_struct, xc_section, &
                             accuracy, is_triplet, spinflip, &
                             fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      INTEGER, INTENT(IN)                                :: accuracy
      LOGICAL, INTENT(IN)                                :: is_triplet, spinflip
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fgxc_fgdiff'
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: handle, ispin, istep, nspins, nstep
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1, oeps2
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak, bl
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_tau_rspace, vxc00
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rhoin

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(.NOT. ASSOCIATED(gxc_rho))
      CPASSERT(.NOT. ASSOCIATED(gxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      ! no spinflip tested
      CPASSERT(.NOT. spinflip)

      ak = 0.0_dp
      bl = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
         bl(-2:2) = [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
         bl(-3:3) = [2.0_dp, -27.0_dp, 270.0_dp, -490.0_dp, 270.0_dp, -27.0_dp, 2.0_dp]/180.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
         bl(-4:4) = [-1.0_dp, 128.0_dp/9.0_dp, -112.0_dp, 896.0_dp, -14350.0_dp/9.0_dp, &
                     896.0_dp, -112.0_dp, 128.0_dp/9.0_dp, -1.0_dp]/560.0_dp
      END SELECT

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = dft_control%nspins
      exc = 0.0_dp

      DO istep = -nstep, nstep

         alpha = 1.0_dp
         beta = REAL(istep, KIND=dp)*epsrho
         NULLIFY (rhoin)
         ALLOCATE (rhoin)
         CALL qs_rho_create(rhoin)
         NULLIFY (vxc00, v_tau_rspace)
         IF (is_triplet) THEN
            CPASSERT(nspins == 1)
            ! rhoin = (0.5 rho0, 0.5 rho0)
            CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, 2)
            ! rhoin = (0.5 rho0 + 0.5 rho1, 0.5 rho0)
            CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, 0.5_dp*beta)
            CALL qs_vxc_create(ks_env, rhoin, xc_section, vxc00, v_tau_rspace, exc)
            CALL pw_axpy(vxc00(2), vxc00(1), -1.0_dp)
         ELSE
            CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, nspins)
            CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, beta)
            CALL qs_vxc_create(ks_env, rhoin, xc_section, vxc00, v_tau_rspace, exc)
         END IF
         CALL qs_rho_release(rhoin)
         DEALLOCATE (rhoin)
         IF (.NOT. ASSOCIATED(fxc_rho)) THEN
            ALLOCATE (fxc_rho(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(fxc_rho(ispin))
               CALL pw_zero(fxc_rho(ispin))
            END DO
         END IF
         IF (.NOT. ASSOCIATED(gxc_rho)) THEN
            ALLOCATE (gxc_rho(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(gxc_rho(ispin))
               CALL pw_zero(gxc_rho(ispin))
            END DO
         END IF
         CPASSERT(.NOT. ASSOCIATED(v_tau_rspace))
         DO ispin = 1, nspins
            IF (ak(istep) /= 0.0_dp) THEN
               CALL pw_axpy(vxc00(ispin), fxc_rho(ispin), ak(istep))
            END IF
            IF (bl(istep) /= 0.0_dp) THEN
               CALL pw_axpy(vxc00(ispin), gxc_rho(ispin), bl(istep))
            END IF
         END DO
         DO ispin = 1, SIZE(vxc00)
            CALL auxbas_pw_pool%give_back_pw(vxc00(ispin))
         END DO
         DEALLOCATE (vxc00)

      END DO

      oeps1 = 1.0_dp/epsrho
      oeps2 = 1.0_dp/(epsrho**2)
      DO ispin = 1, nspins
         CALL pw_scale(fxc_rho(ispin), oeps1)
         CALL pw_scale(gxc_rho(ispin), oeps2)
      END DO

      CALL timestop(handle)

   END SUBROUTINE qs_fgxc_fgdiff

! **************************************************************************************************
!> \brief ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
!> \param pw_pool ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_release(fxc_rho, fxc_tau, gxc_rho, gxc_tau, pw_pool)

      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau, gxc_rho, gxc_tau
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: ispin

      IF (ASSOCIATED(fxc_rho)) THEN
         DO ispin = 1, SIZE(fxc_rho)
            CALL pw_pool%give_back_pw(fxc_rho(ispin))
         END DO
         DEALLOCATE (fxc_rho)
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         DO ispin = 1, SIZE(fxc_tau)
            CALL pw_pool%give_back_pw(fxc_tau(ispin))
         END DO
         DEALLOCATE (fxc_tau)
      END IF
      IF (ASSOCIATED(gxc_rho)) THEN
         DO ispin = 1, SIZE(gxc_rho)
            CALL pw_pool%give_back_pw(gxc_rho(ispin))
         END DO
         DEALLOCATE (gxc_rho)
      END IF
      IF (ASSOCIATED(gxc_tau)) THEN
         DO ispin = 1, SIZE(gxc_tau)
            CALL pw_pool%give_back_pw(gxc_tau(ispin))
         END DO
         DEALLOCATE (gxc_tau)
      END IF

   END SUBROUTINE qs_fgxc_release

END MODULE qs_fgxc
