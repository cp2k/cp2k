!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_vxc

   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: sic_ad,&
                                              sic_eo,&
                                              sic_mauri_spz,&
                                              sic_mauri_us,&
                                              sic_none,&
                                              xc_none,&
                                              xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_multiply,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_dispersion_nonloc,            ONLY: calculate_dispersion_nonloc
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: calc_xc_density,&
                                              xc_exc_calc,&
                                              xc_vxc_pw_create1
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: qs_vxc_create, qs_xc_density

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc'

CONTAINS

! **************************************************************************************************
!> \brief calculates and allocates the xc potential, already reducing it to
!>      the dependence on rho and the one on tau
!> \param ks_env to get all the needed things
!> \param rho_struct density for which v_xc is calculated
!> \param xc_section ...
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic tau part of v_xc
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param exc ...
!> \param just_energy if true calculates just the energy, and does not
!>        allocate v_*_rspace
!> \param edisp ...
!> \param dispersion_env ...
!> \param adiabatic_rescale_factor ...
!> \param pw_env_external    external plane wave environment
!> \par History
!>      - 05.2002 modified to use the mp_allgather function each pe
!>        computes only part of the grid and this is broadcasted to all
!>        instead of summed.
!>        This scales significantly better (e.g. factor 3 on 12 cpus
!>        32 H2O) [Joost VdV]
!>      - moved to qs_ks_methods [fawzi]
!>      - sic alterations [Joost VandeVondele]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_vxc_create(ks_env, rho_struct, xc_section, vxc_rho, vxc_tau, exc, &
                            just_energy, edisp, dispersion_env, adiabatic_rescale_factor, &
                            pw_env_external)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      LOGICAL, INTENT(in), OPTIONAL                      :: just_energy
      REAL(KIND=dp), INTENT(out), OPTIONAL               :: edisp
      TYPE(qs_dispersion_type), OPTIONAL, POINTER        :: dispersion_env
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: adiabatic_rescale_factor
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_vxc_create'

      INTEGER                                            :: handle, ispin, mspin, myfun, &
                                                            nelec_spin(2), vdw
      LOGICAL :: compute_virial, do_adiabatic_rescaling, my_just_energy, rho_g_valid, &
         sic_scaling_b_zero, tau_r_valid, uf_grid, vdW_nl
      REAL(KIND=dp)                                      :: exc_m, factor, &
                                                            my_adiabatic_rescale_factor, &
                                                            my_scaling, nelec_s_inv
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc_tmp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g, rho_m_gspace, rho_struct_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g, tmp_g, tmp_g2
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, vdw_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: my_vxc_rho, my_vxc_tau, rho_m_rspace, &
                                                            rho_r, rho_struct_r, tau, tau_struct_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, tmp_pw
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(vxc_rho))
      CPASSERT(.NOT. ASSOCIATED(vxc_tau))
      NULLIFY (dft_control, pw_env, auxbas_pw_pool, xc_pw_pool, vdw_pw_pool, cell, my_vxc_rho, &
               tmp_pw, tmp_g, tmp_g2, my_vxc_tau, rho_g, rho_r, tau, rho_m_rspace, &
               rho_m_gspace, rho_nlcc, rho_nlcc_g, rho_struct_r, rho_struct_g, tau_struct_r)

      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      my_adiabatic_rescale_factor = 1.0_dp
      do_adiabatic_rescaling = .FALSE.
      IF (PRESENT(adiabatic_rescale_factor)) THEN
         my_adiabatic_rescale_factor = adiabatic_rescale_factor
         do_adiabatic_rescaling = .TRUE.
      END IF

      CALL get_ks_env(ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      cell=cell, &
                      virial=virial, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)

      CALL qs_rho_get(rho_struct, &
                      tau_r_valid=tau_r_valid, &
                      rho_g_valid=rho_g_valid, &
                      rho_r=rho_struct_r, &
                      rho_g=rho_struct_g, &
                      tau_r=tau_struct_r)

      compute_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      IF (compute_virial) THEN
         virial%pv_xc = 0.0_dp
      END IF

      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      CALL section_vals_val_get(xc_section, "VDW_POTENTIAL%POTENTIAL_TYPE", &
                                i_val=vdw)

      vdW_nl = (vdw == xc_vdw_fun_nonloc)
      ! this combination has not been investigated
      CPASSERT(.NOT. (do_adiabatic_rescaling .AND. vdW_nl))
      ! are the necessary inputs available
      IF (.NOT. (PRESENT(dispersion_env) .AND. PRESENT(edisp))) THEN
         vdW_nl = .FALSE.
      END IF
      IF (PRESENT(edisp)) edisp = 0.0_dp

      IF (myfun /= xc_none .OR. vdW_nl) THEN

         ! test if the real space density is available
         CPASSERT(ASSOCIATED(rho_struct))
         IF (dft_control%nspins /= 1 .AND. dft_control%nspins /= 2) &
            CPABORT("nspins must be 1 or 2")
         mspin = SIZE(rho_struct_r)
         IF (dft_control%nspins == 2 .AND. mspin == 1) &
            CPABORT("Spin count mismatch")

         ! there are some options related to SIC here.
         ! Normal DFT computes E(rho_alpha,rho_beta) (or its variant E(2*rho_alpha) for non-LSD)
         ! SIC can             E(rho_alpha,rho_beta)-b*(E(rho_alpha,rho_beta)-E(rho_beta,rho_beta))
         ! or compute          E(rho_alpha,rho_beta)-b*E(rho_alpha-rho_beta,0)

         ! my_scaling is the scaling needed of the standard E(rho_alpha,rho_beta) term
         my_scaling = 1.0_dp
         SELECT CASE (dft_control%sic_method_id)
         CASE (sic_none)
            ! all fine
         CASE (sic_mauri_spz, sic_ad)
            ! no idea yet what to do here in that case
            CPASSERT(.NOT. tau_r_valid)
         CASE (sic_mauri_us)
            my_scaling = 1.0_dp - dft_control%sic_scaling_b
            ! no idea yet what to do here in that case
            CPASSERT(.NOT. tau_r_valid)
         CASE (sic_eo)
            ! NOTHING TO BE DONE
         CASE DEFAULT
            ! this case has not yet been treated here
            CPABORT("NYI")
         END SELECT

         IF (dft_control%sic_scaling_b .EQ. 0.0_dp) THEN
            sic_scaling_b_zero = .TRUE.
         ELSE
            sic_scaling_b_zero = .FALSE.
         END IF

         IF (PRESENT(pw_env_external)) &
            pw_env => pw_env_external
         CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool, auxbas_pw_pool=auxbas_pw_pool)
         uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

         IF (.NOT. uf_grid) THEN
            rho_r => rho_struct_r

            IF (tau_r_valid) THEN
               tau => tau_struct_r
            END IF

            ! for gradient corrected functional the density in g space might
            ! be useful so if we have it, we pass it in
            IF (rho_g_valid) THEN
               rho_g => rho_struct_g
            END IF
         ELSE
            CPASSERT(rho_g_valid)
            ALLOCATE (rho_r(mspin))
            ALLOCATE (rho_g(mspin))
            DO ispin = 1, mspin
               CALL xc_pw_pool%create_pw(rho_g(ispin))
               CALL pw_transfer(rho_struct_g(ispin), rho_g(ispin))
            END DO
            DO ispin = 1, mspin
               CALL xc_pw_pool%create_pw(rho_r(ispin))
               CALL pw_transfer(rho_g(ispin), rho_r(ispin))
            END DO
            IF (tau_r_valid) THEN
               ! tau with finer grids is not implemented (at least not correctly), which this asserts
               CPABORT("tau with finer grids not implemented")
            END IF
         END IF

         ! add the nlcc densities
         IF (ASSOCIATED(rho_nlcc)) THEN
            factor = 1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF

         !
         ! here the rho_r, rho_g, tau is what it should be
         ! we get back the right my_vxc_rho and my_vxc_tau as required
         !
         IF (my_just_energy) THEN
            exc = xc_exc_calc(rho_r=rho_r, tau=tau, &
                              rho_g=rho_g, xc_section=xc_section, &
                              pw_pool=xc_pw_pool)

         ELSE
            CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_r, &
                                   rho_g=rho_g, tau=tau, exc=exc, &
                                   xc_section=xc_section, &
                                   pw_pool=xc_pw_pool, &
                                   compute_virial=compute_virial, &
                                   virial_xc=virial%pv_xc)
         END IF

         ! remove the nlcc densities (keep stuff in original state)
         IF (ASSOCIATED(rho_nlcc)) THEN
            factor = -1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF

         ! calclulate non-local vdW functional
         ! only if this XC_SECTION has it
         ! if yes, we use the dispersion_env from ks_env
         ! this is dangerous, as it assumes a special connection xc_section -> qs_env
         IF (vdW_nl) THEN
            CALL get_ks_env(ks_env=ks_env, para_env=para_env)
            ! no SIC functionals allowed
            CPASSERT(dft_control%sic_method_id == sic_none)
            !
            CALL pw_env_get(pw_env, vdw_pw_pool=vdw_pw_pool)
            IF (my_just_energy) THEN
               CALL calculate_dispersion_nonloc(my_vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                                my_just_energy, vdw_pw_pool, xc_pw_pool, para_env)
            ELSE
               CALL calculate_dispersion_nonloc(my_vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                                my_just_energy, vdw_pw_pool, xc_pw_pool, para_env, virial=virial)
            END IF
         END IF

         !! Apply rescaling to the potential if requested
         IF (.NOT. my_just_energy) THEN
            IF (do_adiabatic_rescaling) THEN
               IF (ASSOCIATED(my_vxc_rho)) THEN
                  DO ispin = 1, SIZE(my_vxc_rho)
                     CALL pw_scale(my_vxc_rho(ispin), my_adiabatic_rescale_factor)
                  END DO
               END IF
            END IF
         END IF

         IF (my_scaling .NE. 1.0_dp) THEN
            exc = exc*my_scaling
            IF (ASSOCIATED(my_vxc_rho)) THEN
               DO ispin = 1, SIZE(my_vxc_rho)
                  CALL pw_scale(my_vxc_rho(ispin), my_scaling)
               END DO
            END IF
            IF (ASSOCIATED(my_vxc_tau)) THEN
               DO ispin = 1, SIZE(my_vxc_tau)
                  CALL pw_scale(my_vxc_tau(ispin), my_scaling)
               END DO
            END IF
         END IF

         ! we have pw data for the xc, qs_ks requests coeff structure, here we transfer
         ! pw -> coeff
         IF (ASSOCIATED(my_vxc_rho)) THEN
            vxc_rho => my_vxc_rho
            NULLIFY (my_vxc_rho)
         END IF
         IF (ASSOCIATED(my_vxc_tau)) THEN
            vxc_tau => my_vxc_tau
            NULLIFY (my_vxc_tau)
         END IF
         IF (uf_grid) THEN
            DO ispin = 1, SIZE(rho_r)
               CALL xc_pw_pool%give_back_pw(rho_r(ispin))
            END DO
            DEALLOCATE (rho_r)
            IF (ASSOCIATED(rho_g)) THEN
               DO ispin = 1, SIZE(rho_g)
                  CALL xc_pw_pool%give_back_pw(rho_g(ispin))
               END DO
               DEALLOCATE (rho_g)
            END IF
         END IF

         ! compute again the xc but now for Exc(m,o) and the opposite sign
         IF (dft_control%sic_method_id .EQ. sic_mauri_spz .AND. .NOT. sic_scaling_b_zero) THEN
            ALLOCATE (rho_m_rspace(2), rho_m_gspace(2))
            CALL xc_pw_pool%create_pw(rho_m_gspace(1))
            CALL xc_pw_pool%create_pw(rho_m_rspace(1))
            CALL pw_copy(rho_struct_r(1), rho_m_rspace(1))
            CALL pw_axpy(rho_struct_r(2), rho_m_rspace(1), alpha=-1._dp)
            CALL pw_copy(rho_struct_g(1), rho_m_gspace(1))
            CALL pw_axpy(rho_struct_g(2), rho_m_gspace(1), alpha=-1._dp)
            ! bit sad, these will be just zero...
            CALL xc_pw_pool%create_pw(rho_m_gspace(2))
            CALL xc_pw_pool%create_pw(rho_m_rspace(2))
            CALL pw_zero(rho_m_rspace(2))
            CALL pw_zero(rho_m_gspace(2))

            IF (my_just_energy) THEN
               exc_m = xc_exc_calc(rho_r=rho_m_rspace, tau=tau, &
                                   rho_g=rho_m_gspace, xc_section=xc_section, &
                                   pw_pool=xc_pw_pool)
            ELSE
               ! virial untested
               CPASSERT(.NOT. compute_virial)
               CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_m_rspace, &
                                      rho_g=rho_m_gspace, tau=tau, exc=exc_m, &
                                      xc_section=xc_section, &
                                      pw_pool=xc_pw_pool, &
                                      compute_virial=.FALSE., &
                                      virial_xc=virial_xc_tmp)
            END IF

            exc = exc - dft_control%sic_scaling_b*exc_m

            ! and take care of the potential only vxc_rho is taken into account
            IF (.NOT. my_just_energy) THEN
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(1), -dft_control%sic_scaling_b)
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(2), dft_control%sic_scaling_b)
               CALL my_vxc_rho(1)%release()
               CALL my_vxc_rho(2)%release()
               DEALLOCATE (my_vxc_rho)
            END IF

            DO ispin = 1, 2
               CALL xc_pw_pool%give_back_pw(rho_m_rspace(ispin))
               CALL xc_pw_pool%give_back_pw(rho_m_gspace(ispin))
            END DO
            DEALLOCATE (rho_m_rspace)
            DEALLOCATE (rho_m_gspace)

         END IF

         ! now we have - sum_s N_s * Exc(rho_s/N_s,0)
         IF (dft_control%sic_method_id .EQ. sic_ad .AND. .NOT. sic_scaling_b_zero) THEN

            ! find out how many elecs we have
            CALL get_ks_env(ks_env, nelectron_spin=nelec_spin)

            ALLOCATE (rho_m_rspace(2), rho_m_gspace(2))
            DO ispin = 1, 2
               CALL xc_pw_pool%create_pw(rho_m_gspace(ispin))
               CALL xc_pw_pool%create_pw(rho_m_rspace(ispin))
            END DO

            DO ispin = 1, 2
               IF (nelec_spin(ispin) .GT. 0.0_dp) THEN
                  nelec_s_inv = 1.0_dp/nelec_spin(ispin)
               ELSE
                  ! does it matter if there are no electrons with this spin (H) ?
                  nelec_s_inv = 0.0_dp
               END IF
               CALL pw_copy(rho_struct_r(ispin), rho_m_rspace(1))
               CALL pw_copy(rho_struct_g(ispin), rho_m_gspace(1))
               CALL pw_scale(rho_m_rspace(1), nelec_s_inv)
               CALL pw_scale(rho_m_gspace(1), nelec_s_inv)
               CALL pw_zero(rho_m_rspace(2))
               CALL pw_zero(rho_m_gspace(2))

               IF (my_just_energy) THEN
                  exc_m = xc_exc_calc(rho_r=rho_m_rspace, tau=tau, &
                                      rho_g=rho_m_gspace, xc_section=xc_section, &
                                      pw_pool=xc_pw_pool)
               ELSE
                  ! virial untested
                  CPASSERT(.NOT. compute_virial)
                  CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_m_rspace, &
                                         rho_g=rho_m_gspace, tau=tau, exc=exc_m, &
                                         xc_section=xc_section, &
                                         pw_pool=xc_pw_pool, &
                                         compute_virial=.FALSE., &
                                         virial_xc=virial_xc_tmp)
               END IF

               exc = exc - dft_control%sic_scaling_b*nelec_spin(ispin)*exc_m

               ! and take care of the potential only vxc_rho is taken into account
               IF (.NOT. my_just_energy) THEN
                  CALL pw_axpy(my_vxc_rho(1), vxc_rho(ispin), -dft_control%sic_scaling_b)
                  CALL my_vxc_rho(1)%release()
                  CALL my_vxc_rho(2)%release()
                  DEALLOCATE (my_vxc_rho)
               END IF
            END DO

            DO ispin = 1, 2
               CALL xc_pw_pool%give_back_pw(rho_m_rspace(ispin))
               CALL xc_pw_pool%give_back_pw(rho_m_gspace(ispin))
            END DO
            DEALLOCATE (rho_m_rspace)
            DEALLOCATE (rho_m_gspace)

         END IF

         ! compute again the xc but now for Exc(n_down,n_down)
         IF (dft_control%sic_method_id .EQ. sic_mauri_us .AND. .NOT. sic_scaling_b_zero) THEN
            ALLOCATE (rho_r(2))
            rho_r(1) = rho_struct_r(2)
            rho_r(2) = rho_struct_r(2)
            IF (rho_g_valid) THEN
               ALLOCATE (rho_g(2))
               rho_g(1) = rho_struct_g(2)
               rho_g(2) = rho_struct_g(2)
            END IF

            IF (my_just_energy) THEN
               exc_m = xc_exc_calc(rho_r=rho_r, tau=tau, &
                                   rho_g=rho_g, xc_section=xc_section, &
                                   pw_pool=xc_pw_pool)
            ELSE
               ! virial untested
               CPASSERT(.NOT. compute_virial)
               CALL xc_vxc_pw_create1(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_r, &
                                      rho_g=rho_g, tau=tau, exc=exc_m, &
                                      xc_section=xc_section, &
                                      pw_pool=xc_pw_pool, &
                                      compute_virial=.FALSE., &
                                      virial_xc=virial_xc_tmp)
            END IF

            exc = exc + dft_control%sic_scaling_b*exc_m

            ! and take care of the potential
            IF (.NOT. my_just_energy) THEN
               ! both go to minority spin
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(2), 2.0_dp*dft_control%sic_scaling_b)
               CALL my_vxc_rho(1)%release()
               CALL my_vxc_rho(2)%release()
               DEALLOCATE (my_vxc_rho)
            END IF
            DEALLOCATE (rho_r, rho_g)

         END IF

         !
         ! cleanups
         !
         IF (uf_grid .AND. (ASSOCIATED(vxc_rho) .OR. ASSOCIATED(vxc_tau))) THEN
            BLOCK
               TYPE(pw_r3d_rs_type) :: tmp_pw
               TYPE(pw_c1d_gs_type) :: tmp_g, tmp_g2
               CALL xc_pw_pool%create_pw(tmp_g)
               CALL auxbas_pw_pool%create_pw(tmp_g2)
               IF (ASSOCIATED(vxc_rho)) THEN
                  DO ispin = 1, SIZE(vxc_rho)
                     CALL auxbas_pw_pool%create_pw(tmp_pw)
                     CALL pw_transfer(vxc_rho(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, tmp_g2)
                     CALL pw_transfer(tmp_g2, tmp_pw)
                     CALL xc_pw_pool%give_back_pw(vxc_rho(ispin))
                     vxc_rho(ispin) = tmp_pw
                  END DO
               END IF
               IF (ASSOCIATED(vxc_tau)) THEN
                  DO ispin = 1, SIZE(vxc_tau)
                     CALL auxbas_pw_pool%create_pw(tmp_pw)
                     CALL pw_transfer(vxc_tau(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, tmp_g2)
                     CALL pw_transfer(tmp_g2, tmp_pw)
                     CALL xc_pw_pool%give_back_pw(vxc_tau(ispin))
                     vxc_tau(ispin) = tmp_pw
                  END DO
               END IF
               CALL auxbas_pw_pool%give_back_pw(tmp_g2)
               CALL xc_pw_pool%give_back_pw(tmp_g)
            END BLOCK
         END IF
         IF (ASSOCIATED(tau) .AND. uf_grid) THEN
            DO ispin = 1, SIZE(tau)
               CALL xc_pw_pool%give_back_pw(tau(ispin))
            END DO
            DEALLOCATE (tau)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_vxc_create

! **************************************************************************************************
!> \brief calculates the XC density: E_xc(r) - V_xc(r)*rho(r)  or  E_xc(r)/rho(r)
!> \param ks_env to get all the needed things
!> \param rho_struct density
!> \param xc_section ...
!> \param dispersion_env ...
!> \param xc_ener will contain the xc energy density E_xc(r) - V_xc(r)*rho(r)
!> \param xc_den will contain the xc energy density E_xc(r)/rho(r)
!> \param vxc ...
!> \param vtau ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env, &
                            xc_ener, xc_den, vxc, vtau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(qs_dispersion_type), OPTIONAL, POINTER        :: dispersion_env
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: xc_ener, xc_den
      TYPE(pw_r3d_rs_type), DIMENSION(:), OPTIONAL       :: vxc, vtau

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_xc_density'

      INTEGER                                            :: handle, ispin, myfun, nspins, vdw
      LOGICAL                                            :: rho_g_valid, tau_r_valid, uf_grid, vdW_nl
      REAL(KIND=dp)                                      :: edisp, exc, factor, rho_cutoff
      REAL(KIND=dp), DIMENSION(3, 3)                     :: vdum
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, vdw_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type)                               :: exc_r
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau_r, vxc_rho, vxc_tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc

      CALL timeset(routineN, handle)

      CALL get_ks_env(ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      cell=cell, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)

      CALL qs_rho_get(rho_struct, &
                      tau_r_valid=tau_r_valid, &
                      rho_g_valid=rho_g_valid, &
                      rho_r=rho_r, &
                      rho_g=rho_g, &
                      tau_r=tau_r)
      nspins = dft_control%nspins

      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)
      CALL section_vals_val_get(xc_section, "VDW_POTENTIAL%POTENTIAL_TYPE", i_val=vdw)
      vdW_nl = (vdw == xc_vdw_fun_nonloc)
      IF (PRESENT(xc_ener)) THEN
         IF (tau_r_valid) THEN
            CALL cp_warn(__LOCATION__, "Tau contribution will not be correctly handled")
         END IF
      END IF
      IF (vdW_nl) THEN
         CALL cp_warn(__LOCATION__, "vdW functional contribution will be ignored")
      END IF

      CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool, auxbas_pw_pool=auxbas_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)
      IF (uf_grid) THEN
         CALL cp_warn(__LOCATION__, "Fine grid option not possible with local energy")
         CPABORT("Fine Grid in Local Energy")
      END IF

      IF (PRESENT(xc_ener)) THEN
         CALL pw_zero(xc_ener)
      END IF
      IF (PRESENT(xc_den)) THEN
         CALL pw_zero(xc_den)
      END IF
      IF (PRESENT(vxc)) THEN
         DO ispin = 1, nspins
            CALL pw_zero(vxc(ispin))
         END DO
      END IF
      IF (PRESENT(vtau)) THEN
         DO ispin = 1, nspins
            CALL pw_zero(vtau(ispin))
         END DO
      END IF

      IF (myfun /= xc_none) THEN

         CPASSERT(ASSOCIATED(rho_struct))
         CPASSERT(dft_control%sic_method_id == sic_none)

         ! add the nlcc densities
         IF (ASSOCIATED(rho_nlcc)) THEN
            factor = 1.0_dp
            DO ispin = 1, nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF
         NULLIFY (vxc_rho, vxc_tau)
         CALL xc_vxc_pw_create1(vxc_rho=vxc_rho, vxc_tau=vxc_tau, rho_r=rho_r, &
                                rho_g=rho_g, tau=tau_r, exc=exc, &
                                xc_section=xc_section, &
                                pw_pool=xc_pw_pool, &
                                compute_virial=.FALSE., &
                                virial_xc=vdum, &
                                exc_r=exc_r)
         ! calclulate non-local vdW functional
         ! only if this XC_SECTION has it
         ! if yes, we use the dispersion_env from ks_env
         ! this is dangerous, as it assumes a special connection xc_section -> qs_env
         IF (vdW_nl) THEN
            CALL get_ks_env(ks_env=ks_env, para_env=para_env)
            ! no SIC functionals allowed
            CPASSERT(dft_control%sic_method_id == sic_none)
            !
            CALL pw_env_get(pw_env, vdw_pw_pool=vdw_pw_pool)
            CALL calculate_dispersion_nonloc(vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                             .FALSE., vdw_pw_pool, xc_pw_pool, para_env)
         END IF

         ! remove the nlcc densities (keep stuff in original state)
         IF (ASSOCIATED(rho_nlcc)) THEN
            factor = -1.0_dp
            DO ispin = 1, dft_control%nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF
         !
         IF (PRESENT(xc_den)) THEN
            CALL pw_copy(exc_r, xc_den)
            rho_cutoff = 1.E-14_dp
            CALL calc_xc_density(xc_den, rho_r, rho_cutoff)
         END IF
         IF (PRESENT(xc_ener)) THEN
            CALL pw_copy(exc_r, xc_ener)
            DO ispin = 1, nspins
               CALL pw_multiply(xc_ener, vxc_rho(ispin), rho_r(ispin), alpha=-1.0_dp)
            END DO
         END IF
         IF (PRESENT(vxc)) THEN
            DO ispin = 1, nspins
               CALL pw_copy(vxc_rho(ispin), vxc(ispin))
            END DO
         END IF
         IF (PRESENT(vtau)) THEN
            DO ispin = 1, nspins
               CALL pw_copy(vxc_tau(ispin), vtau(ispin))
            END DO
         END IF
         ! remove arrays
         IF (ASSOCIATED(vxc_rho)) THEN
            DO ispin = 1, nspins
               CALL vxc_rho(ispin)%release()
            END DO
            DEALLOCATE (vxc_rho)
         END IF
         IF (ASSOCIATED(vxc_tau)) THEN
            DO ispin = 1, nspins
               CALL vxc_tau(ispin)%release()
            END DO
            DEALLOCATE (vxc_tau)
         END IF
         CALL exc_r%release()
         !
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_xc_density

END MODULE qs_vxc
