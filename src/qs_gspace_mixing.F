!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_gspace_mixing

   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invert_matrix
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_density_mixing_types,         ONLY: broyden_mixing_nr,&
                                              gspace_mixing_nr,&
                                              mixing_storage_type,&
                                              multisecant_mixing_nr,&
                                              pulay_mixing_nr
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gspace_mixing'

   PUBLIC :: gspace_mixing

CONTAINS

! **************************************************************************************************
!> \brief  Driver for the g-space mixing, calls the proper routine given the
!>requested method
!> \param qs_env ...
!> \param mixing_method ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \param iter_count ...
!> \par History
!>      05.2009
!>      02.2015 changed input to make g-space mixing available in linear scaling
!>              SCF [Patrick Seewald]
!> \author MI
! **************************************************************************************************
   SUBROUTINE gspace_mixing(qs_env, mixing_method, mixing_store, rho, para_env, iter_count)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: mixing_method
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: iter_count

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gspace_mixing'

      INTEGER                                            :: handle, iatom, ig, ispin, natom, ng, &
                                                            nimg, nspin
      LOGICAL                                            :: gapw
      REAL(dp)                                           :: alpha
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type)                               :: rho_tmp
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env)
      IF (.NOT. (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb .OR. &
                 dft_control%qs_control%semi_empirical)) THEN

         CPASSERT(ASSOCIATED(mixing_store))
         NULLIFY (auxbas_pw_pool, rho_atom, rho_g, rho_r, tot_rho_r)
         CALL qs_rho_get(rho, rho_g=rho_g, rho_r=rho_r, tot_rho_r=tot_rho_r)

         nspin = SIZE(rho_g, 1)
         nimg = dft_control%nimages
         ng = SIZE(rho_g(1)%pw_grid%gsq)
         CPASSERT(ng == SIZE(mixing_store%rhoin(1)%cc))

         alpha = mixing_store%alpha
         gapw = dft_control%qs_control%gapw

         IF (nspin == 2) THEN
            CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL auxbas_pw_pool%create_pw(rho_tmp)
            CALL pw_zero(rho_tmp)
            CALL pw_copy(rho_g(1), rho_tmp)
            CALL pw_axpy(rho_g(2), rho_g(1), 1.0_dp)
            CALL pw_axpy(rho_tmp, rho_g(2), -1.0_dp)
            CALL pw_scale(rho_g(2), -1.0_dp)
         END IF

         IF (iter_count + 1 <= mixing_store%nskip_mixing) THEN
            ! skip mixing
            DO ispin = 1, nspin
               DO ig = 1, ng
                  mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%array(ig)
               END DO
            END DO
            IF (mixing_store%gmix_p .AND. gapw) THEN
               CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
               natom = SIZE(rho_atom)
               DO ispin = 1, nspin
                  DO iatom = 1, natom
                     IF (mixing_store%paw(iatom)) THEN
                        mixing_store%cpc_h_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                        mixing_store%cpc_s_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
                     END IF
                  END DO
               END DO
            END IF

            mixing_store%iter_method = "NoMix"
            IF (nspin == 2) THEN
               CALL pw_axpy(rho_g(2), rho_g(1), 1.0_dp)
               CALL pw_scale(rho_g(1), 0.5_dp)
               CALL pw_axpy(rho_g(1), rho_g(2), -1.0_dp)
               CALL pw_scale(rho_g(2), -1.0_dp)
               CALL auxbas_pw_pool%give_back_pw(rho_tmp)
            END IF
            CALL timestop(handle)
            RETURN
         END IF

         IF ((iter_count + 1 - mixing_store%nskip_mixing) <= mixing_store%n_simple_mix) THEN
            CALL gmix_potential_only(qs_env, mixing_store, rho)
            mixing_store%iter_method = "Kerker"
         ELSEIF (mixing_method == gspace_mixing_nr) THEN
            CALL gmix_potential_only(qs_env, mixing_store, rho)
            mixing_store%iter_method = "Kerker"
         ELSEIF (mixing_method == pulay_mixing_nr) THEN
            CALL pulay_mixing(qs_env, mixing_store, rho, para_env)
            mixing_store%iter_method = "Pulay"
         ELSEIF (mixing_method == broyden_mixing_nr) THEN
            CALL broyden_mixing(qs_env, mixing_store, rho, para_env)
            mixing_store%iter_method = "Broy."
         ELSEIF (mixing_method == multisecant_mixing_nr) THEN
            CPASSERT(.NOT. gapw)
            CALL multisecant_mixing(mixing_store, rho, para_env)
            mixing_store%iter_method = "MSec."
         END IF

         IF (nspin == 2) THEN
            CALL pw_axpy(rho_g(2), rho_g(1), 1.0_dp)
            CALL pw_scale(rho_g(1), 0.5_dp)
            CALL pw_axpy(rho_g(1), rho_g(2), -1.0_dp)
            CALL pw_scale(rho_g(2), -1.0_dp)
            CALL auxbas_pw_pool%give_back_pw(rho_tmp)
         END IF

         DO ispin = 1, nspin
            CALL pw_transfer(rho_g(ispin), rho_r(ispin))
            tot_rho_r(ispin) = pw_integrate_function(rho_r(ispin), isign=-1)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE gspace_mixing

! **************************************************************************************************
!> \brief G-space mixing performed via the Kerker damping on the density on the grid
!>        thus affecting only the caluclation of the potential, but not the denisity matrix
!> \param qs_env ...
!> \param mixing_store ...
!> \param rho ...
!> \par History
!>      02.2009
!> \author MI
! **************************************************************************************************

   SUBROUTINE gmix_potential_only(qs_env, mixing_store, rho)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(qs_rho_type), POINTER                         :: rho

      CHARACTER(len=*), PARAMETER :: routineN = 'gmix_potential_only'

      COMPLEX(dp), DIMENSION(:), POINTER                 :: cc_new
      INTEGER                                            :: handle, iatom, ig, ispin, natom, ng, &
                                                            nspin
      LOGICAL                                            :: gapw
      REAL(dp)                                           :: alpha, f_mix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CPASSERT(ASSOCIATED(mixing_store%rhoin))
      CPASSERT(ASSOCIATED(mixing_store%kerker_factor))

      CALL timeset(routineN, handle)
      NULLIFY (cc_new, dft_control, rho_atom, rho_g)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL qs_rho_get(rho, rho_g=rho_g)

      nspin = SIZE(rho_g, 1)
      ng = SIZE(rho_g(1)%pw_grid%gsq)

      gapw = dft_control%qs_control%gapw
      alpha = mixing_store%alpha

      DO ispin = 1, nspin
         cc_new => rho_g(ispin)%array
         DO ig = 1, mixing_store%ig_max ! ng
            f_mix = mixing_store%alpha*mixing_store%kerker_factor(ig)
            cc_new(ig) = (1.0_dp - f_mix)*mixing_store%rhoin(ispin)%cc(ig) + f_mix*cc_new(ig)
            mixing_store%rhoin(ispin)%cc(ig) = cc_new(ig)
         END DO
         DO ig = mixing_store%ig_max + 1, ng
            f_mix = mixing_store%alpha
            cc_new(ig) = (1.0_dp - f_mix)*mixing_store%rhoin(ispin)%cc(ig) + f_mix*cc_new(ig)
            mixing_store%rhoin(ispin)%cc(ig) = cc_new(ig)
         END DO

      END DO

      IF (mixing_store%gmix_p .AND. gapw) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         rho_atom_set=rho_atom)
         natom = SIZE(rho_atom)
         DO ispin = 1, nspin
            DO iatom = 1, natom
               IF (mixing_store%paw(iatom)) THEN
                  rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                                                        mixing_store%cpc_h_in(iatom, ispin)%r_coef*(1._dp - alpha)
                  rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                                                        mixing_store%cpc_s_in(iatom, ispin)%r_coef*(1._dp - alpha)
                  mixing_store%cpc_h_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                  mixing_store%cpc_s_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
               END IF
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE gmix_potential_only

! **************************************************************************************************
!> \brief Pulay Mixing using as metrics for the residual the Kerer damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \param qs_env ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      03.2009
!> \author MI
! **************************************************************************************************

   SUBROUTINE pulay_mixing(qs_env, mixing_store, rho, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pulay_mixing'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: cc_mix
      INTEGER                                            :: handle, i, iatom, ib, ibb, ig, ispin, &
                                                            jb, n1, n2, natom, nb, nb1, nbuffer, &
                                                            ng, nspin
      LOGICAL                                            :: gapw
      REAL(dp)                                           :: alpha_kerker, alpha_pulay, beta, f_mix, &
                                                            inv_err, norm, norm_c_inv, res_norm, &
                                                            vol
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: alpha_c, ev
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: b, c, c_inv, cpc_h_mix, cpc_s_mix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CPASSERT(ASSOCIATED(mixing_store%res_buffer))
      CPASSERT(ASSOCIATED(mixing_store%rhoin_buffer))
      NULLIFY (dft_control, rho_atom, rho_g)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL qs_rho_get(rho, rho_g=rho_g)
      nspin = SIZE(rho_g, 1)
      ng = SIZE(mixing_store%res_buffer(1, 1)%cc)
      vol = rho_g(1)%pw_grid%vol

      alpha_kerker = mixing_store%alpha
      beta = mixing_store%pulay_beta
      alpha_pulay = mixing_store%pulay_alpha
      nbuffer = mixing_store%nbuffer
      gapw = dft_control%qs_control%gapw
      IF (gapw) THEN
         CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
         natom = SIZE(rho_atom)
      END IF

      ALLOCATE (cc_mix(ng))

      IF (nspin == 2) THEN
         CALL pw_axpy(rho_g(1), rho_g(2), 1.0_dp, -1.0_dp)
         DO ig = 1, ng
            mixing_store%rhoin(2)%cc(ig) = mixing_store%rhoin(1)%cc(ig) - mixing_store%rhoin(2)%cc(ig)
         END DO
         IF (gapw .AND. mixing_store%gmix_p) THEN
            DO iatom = 1, natom
               IF (mixing_store%paw(iatom)) THEN
                  rho_atom(iatom)%cpc_h(2)%r_coef = rho_atom(iatom)%cpc_h(1)%r_coef - rho_atom(iatom)%cpc_h(2)%r_coef
                  rho_atom(iatom)%cpc_s(2)%r_coef = rho_atom(iatom)%cpc_s(1)%r_coef - rho_atom(iatom)%cpc_s(2)%r_coef
                  mixing_store%cpc_h_in(iatom, 2)%r_coef = mixing_store%cpc_h_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_h_in(iatom, 2)%r_coef
                  mixing_store%cpc_s_in(iatom, 2)%r_coef = mixing_store%cpc_s_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_s_in(iatom, 2)%r_coef
               END IF
            END DO
         END IF
      END IF

      DO ispin = 1, nspin
         ib = MODULO(mixing_store%ncall_p(ispin), nbuffer) + 1

         mixing_store%ncall_p(ispin) = mixing_store%ncall_p(ispin) + 1
         nb = MIN(mixing_store%ncall_p(ispin), nbuffer)
         ibb = MODULO(mixing_store%ncall_p(ispin), nbuffer) + 1

         mixing_store%res_buffer(ib, ispin)%cc(:) = CMPLX(0._dp, 0._dp, KIND=dp)
         norm = 0.0_dp
         IF (nb == 1) mixing_store%rhoin_buffer(1, ispin)%cc = mixing_store%rhoin(ispin)%cc
         res_norm = 0.0_dp
         DO ig = 1, ng
            f_mix = mixing_store%kerker_factor(ig)
            mixing_store%res_buffer(ib, ispin)%cc(ig) = f_mix*(rho_g(ispin)%array(ig) - &
                                                               mixing_store%rhoin_buffer(ib, ispin)%cc(ig))
            res_norm = res_norm + &
                       REAL(mixing_store%res_buffer(ib, ispin)%cc(ig), dp)*REAL(mixing_store%res_buffer(ib, ispin)%cc(ig), dp) + &
                       AIMAG(mixing_store%res_buffer(ib, ispin)%cc(ig))*AIMAG(mixing_store%res_buffer(ib, ispin)%cc(ig))
         END DO

         CALL para_env%sum(res_norm)
         res_norm = SQRT(res_norm)

         IF (res_norm < 1.E-14_dp) THEN
            mixing_store%ncall_p(ispin) = mixing_store%ncall_p(ispin) - 1
         ELSE
            nb1 = nb + 1
            IF (ALLOCATED(b)) DEALLOCATE (b)
            ALLOCATE (b(nb1, nb1))
            b = 0.0_dp
            IF (ALLOCATED(c)) DEALLOCATE (c)
            ALLOCATE (c(nb, nb))
            c = 0.0_dp
            IF (ALLOCATED(c_inv)) DEALLOCATE (c_inv)
            ALLOCATE (c_inv(nb, nb))
            c_inv = 0.0_dp
            IF (ALLOCATED(alpha_c)) DEALLOCATE (alpha_c)
            ALLOCATE (alpha_c(nb))
            alpha_c = 0.0_dp
            norm_c_inv = 0.0_dp
            IF (ALLOCATED(ev)) DEALLOCATE (ev)
            ALLOCATE (ev(nb1))
            ev = 0.0_dp

         END IF

         DO jb = 1, nb
            mixing_store%pulay_matrix(jb, ib, ispin) = 0.0_dp
            DO ig = 1, ng
               f_mix = mixing_store%special_metric(ig)
               mixing_store%pulay_matrix(jb, ib, ispin) = mixing_store%pulay_matrix(jb, ib, ispin) + &
                                                          f_mix*(REAL(mixing_store%res_buffer(jb, ispin)%cc(ig), dp) &
                                                                 *REAL(mixing_store%res_buffer(ib, ispin)%cc(ig), dp) + &
                                                                 AIMAG(mixing_store%res_buffer(jb, ispin)%cc(ig))* &
                                                                 AIMAG(mixing_store%res_buffer(ib, ispin)%cc(ig)))
            END DO
            CALL para_env%sum(mixing_store%pulay_matrix(jb, ib, ispin))
            mixing_store%pulay_matrix(jb, ib, ispin) = mixing_store%pulay_matrix(jb, ib, ispin)*vol*vol
            mixing_store%pulay_matrix(ib, jb, ispin) = mixing_store%pulay_matrix(jb, ib, ispin)
         END DO

         IF (nb == 1 .OR. res_norm < 1.E-14_dp) THEN
            DO ig = 1, ng
               f_mix = alpha_kerker*mixing_store%kerker_factor(ig)
               cc_mix(ig) = rho_g(ispin)%array(ig) - &
                            mixing_store%rhoin_buffer(ib, ispin)%cc(ig)
               rho_g(ispin)%array(ig) = f_mix*cc_mix(ig) + &
                                        mixing_store%rhoin_buffer(ib, ispin)%cc(ig)
               mixing_store%rhoin_buffer(ibb, ispin)%cc(ig) = rho_g(ispin)%array(ig)
            END DO
            IF (mixing_store%gmix_p) THEN
               IF (gapw) THEN
                  DO iatom = 1, natom
                     IF (mixing_store%paw(iatom)) THEN
                        mixing_store%cpc_h_res_buffer(ib, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef - &
                                                                                 mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        mixing_store%cpc_s_res_buffer(ib, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef - &
                                                                                 mixing_store%cpc_s_in(iatom, ispin)%r_coef

                        rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha_kerker*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                                                              (1.0_dp - alpha_kerker)*mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha_kerker*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                                                              (1.0_dp - alpha_kerker)*mixing_store%cpc_s_in(iatom, ispin)%r_coef

                        mixing_store%cpc_h_in_buffer(ib, iatom, ispin)%r_coef = mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        mixing_store%cpc_s_in_buffer(ib, iatom, ispin)%r_coef = mixing_store%cpc_s_in(iatom, ispin)%r_coef

                        mixing_store%cpc_h_in_buffer(ibb, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                        mixing_store%cpc_s_in_buffer(ibb, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
                     END IF
                  END DO
               END IF
            END IF
         ELSE

            b(1:nb, 1:nb) = mixing_store%pulay_matrix(1:nb, 1:nb, ispin)
            c(1:nb, 1:nb) = b(1:nb, 1:nb)
            b(nb1, 1:nb) = -1.0_dp
            b(1:nb, nb1) = -1.0_dp
            b(nb1, nb1) = 0.0_dp

            CALL invert_matrix(c, c_inv, inv_err, improve=.TRUE.)
            alpha_c = 0.0_dp
            DO i = 1, nb
               DO jb = 1, nb
                  alpha_c(i) = alpha_c(i) + c_inv(jb, i)
                  norm_c_inv = norm_c_inv + c_inv(jb, i)
               END DO
            END DO
            alpha_c(1:nb) = alpha_c(1:nb)/norm_c_inv
            cc_mix = CMPLX(0._dp, 0._dp, KIND=dp)
            DO jb = 1, nb
               DO ig = 1, ng
                  cc_mix(ig) = cc_mix(ig) + &
                               alpha_c(jb)*(mixing_store%rhoin_buffer(jb, ispin)%cc(ig) + &
                                            mixing_store%pulay_beta*mixing_store%res_buffer(jb, ispin)%cc(ig))
               END DO
            END DO
            mixing_store%rhoin_buffer(ibb, ispin)%cc = CMPLX(0._dp, 0._dp, KIND=dp)
            IF (alpha_pulay > 0.0_dp) THEN
               DO ig = 1, ng
                  f_mix = alpha_pulay*mixing_store%kerker_factor(ig)
                  rho_g(ispin)%array(ig) = f_mix*rho_g(ispin)%array(ig) + &
                                           (1.0_dp - f_mix)*cc_mix(ig)
                  mixing_store%rhoin_buffer(ibb, ispin)%cc(ig) = rho_g(ispin)%array(ig)
               END DO
            ELSE
               DO ig = 1, ng
                  rho_g(ispin)%array(ig) = cc_mix(ig)
                  mixing_store%rhoin_buffer(ibb, ispin)%cc(ig) = rho_g(ispin)%array(ig)
               END DO
            END IF

            IF (mixing_store%gmix_p .AND. gapw) THEN
               CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
               DO iatom = 1, natom
                  IF (mixing_store%paw(iatom)) THEN
                     n1 = SIZE(rho_atom(iatom)%cpc_h(ispin)%r_coef, 1)
                     n2 = SIZE(rho_atom(iatom)%cpc_h(ispin)%r_coef, 2)
                     ALLOCATE (cpc_h_mix(n1, n2))
                     ALLOCATE (cpc_s_mix(n1, n2))

                     mixing_store%cpc_h_res_buffer(ib, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef - &
                                                                              mixing_store%cpc_h_in_buffer(ib, iatom, ispin)%r_coef
                     mixing_store%cpc_s_res_buffer(ib, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef - &
                                                                              mixing_store%cpc_s_in_buffer(ib, iatom, ispin)%r_coef
                     cpc_h_mix = 0.0_dp
                     cpc_s_mix = 0.0_dp
                     DO jb = 1, nb
                        cpc_h_mix(:, :) = cpc_h_mix(:, :) + &
                                          alpha_c(jb)*mixing_store%cpc_h_in_buffer(jb, iatom, ispin)%r_coef(:, :) + &
                                          alpha_c(jb)*beta*mixing_store%cpc_h_res_buffer(jb, iatom, ispin)%r_coef(:, :)
                        cpc_s_mix(:, :) = cpc_s_mix(:, :) + &
                                          alpha_c(jb)*mixing_store%cpc_s_in_buffer(jb, iatom, ispin)%r_coef(:, :) + &
                                          alpha_c(jb)*beta*mixing_store%cpc_s_res_buffer(jb, iatom, ispin)%r_coef(:, :)
                     END DO
                     rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha_pulay*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                                                           (1.0_dp - alpha_pulay)*cpc_h_mix
                     rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha_pulay*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                                                           (1.0_dp - alpha_pulay)*cpc_s_mix
                     mixing_store%cpc_h_in_buffer(ibb, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                     mixing_store%cpc_s_in_buffer(ibb, iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
                     DEALLOCATE (cpc_h_mix)
                     DEALLOCATE (cpc_s_mix)
                  END IF
               END DO
            END IF
         END IF ! nb==1

         IF (res_norm > 1.E-14_dp) THEN
            DEALLOCATE (b)
            DEALLOCATE (ev)
            DEALLOCATE (c, c_inv, alpha_c)
         END IF

      END DO ! ispin

      DEALLOCATE (cc_mix)

      IF (nspin == 2) THEN
         CALL pw_axpy(rho_g(1), rho_g(2), 1.0_dp, -1.0_dp)
         DO ig = 1, ng
            mixing_store%rhoin(2)%cc(ig) = mixing_store%rhoin(1)%cc(ig) - mixing_store%rhoin(2)%cc(ig)
         END DO
         IF (gapw .AND. mixing_store%gmix_p) THEN
            DO iatom = 1, natom
               IF (mixing_store%paw(iatom)) THEN
                  rho_atom(iatom)%cpc_h(2)%r_coef = rho_atom(iatom)%cpc_h(1)%r_coef - rho_atom(iatom)%cpc_h(2)%r_coef
                  rho_atom(iatom)%cpc_s(2)%r_coef = rho_atom(iatom)%cpc_s(1)%r_coef - rho_atom(iatom)%cpc_s(2)%r_coef
                  mixing_store%cpc_h_in(iatom, 2)%r_coef = mixing_store%cpc_h_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_h_in(iatom, 2)%r_coef
                  mixing_store%cpc_s_in(iatom, 2)%r_coef = mixing_store%cpc_s_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_s_in(iatom, 2)%r_coef
               END IF
            END DO
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE pulay_mixing

! **************************************************************************************************
!> \brief Broyden Mixing using as metrics for the residual the Kerer damping factor
!>        The mixing is applied directly on the density expansion in reciprocal space
!> \param qs_env ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      03.2009
!> \author MI
! **************************************************************************************************

   SUBROUTINE broyden_mixing(qs_env, mixing_store, rho, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'broyden_mixing'

      COMPLEX(dp)                                        :: cc_mix
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: res_rho
      INTEGER                                            :: handle, i, iatom, ib, ig, ispin, j, jb, &
                                                            kb, n1, n2, natom, nb, nbuffer, ng, &
                                                            nspin
      LOGICAL                                            :: gapw, only_kerker
      REAL(dp)                                           :: alpha, dcpc_h_res, dcpc_s_res, &
                                                            delta_norm, f_mix, inv_err, res_norm, &
                                                            rho_norm, valh, vals, w0
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: c, g
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: a, b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CPASSERT(ASSOCIATED(mixing_store%res_buffer))
      CPASSERT(ASSOCIATED(mixing_store%rhoin))
      CPASSERT(ASSOCIATED(mixing_store%rhoin_old))
      CPASSERT(ASSOCIATED(mixing_store%drho_buffer))
      NULLIFY (dft_control, rho_atom, rho_g)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL qs_rho_get(rho, rho_g=rho_g)

      nspin = SIZE(rho_g, 1)
      ng = SIZE(mixing_store%res_buffer(1, 1)%cc)

      alpha = mixing_store%alpha
      w0 = mixing_store%broy_w0
      nbuffer = mixing_store%nbuffer
      gapw = dft_control%qs_control%gapw

      ALLOCATE (res_rho(ng))

      mixing_store%ncall = mixing_store%ncall + 1
      IF (mixing_store%ncall == 1) THEN
         only_kerker = .TRUE.
      ELSE
         only_kerker = .FALSE.
         nb = MIN(mixing_store%ncall - 1, nbuffer)
         ib = MODULO(mixing_store%ncall - 2, nbuffer) + 1
      END IF

      IF (gapw) THEN
         CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
         natom = SIZE(rho_atom)
      ELSE
         natom = 0
      END IF

      IF (nspin == 2) THEN
         CALL pw_axpy(rho_g(1), rho_g(2), 1.0_dp, -1.0_dp)
         DO ig = 1, ng
            mixing_store%rhoin(2)%cc(ig) = mixing_store%rhoin(1)%cc(ig) - mixing_store%rhoin(2)%cc(ig)
         END DO
         IF (gapw .AND. mixing_store%gmix_p) THEN
            DO iatom = 1, natom
               IF (mixing_store%paw(iatom)) THEN
                  rho_atom(iatom)%cpc_h(2)%r_coef = rho_atom(iatom)%cpc_h(1)%r_coef - rho_atom(iatom)%cpc_h(2)%r_coef
                  rho_atom(iatom)%cpc_s(2)%r_coef = rho_atom(iatom)%cpc_s(1)%r_coef - rho_atom(iatom)%cpc_s(2)%r_coef
                  mixing_store%cpc_h_in(iatom, 2)%r_coef = mixing_store%cpc_h_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_h_in(iatom, 2)%r_coef
                  mixing_store%cpc_s_in(iatom, 2)%r_coef = mixing_store%cpc_s_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_s_in(iatom, 2)%r_coef
               END IF
            END DO
         END IF
      END IF

      DO ispin = 1, nspin
         res_rho = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ig = 1, ng
            res_rho(ig) = rho_g(ispin)%array(ig) - mixing_store%rhoin(ispin)%cc(ig)
         END DO

         IF (only_kerker) THEN
            DO ig = 1, ng
               mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
               f_mix = alpha*mixing_store%kerker_factor(ig)
               rho_g(ispin)%array(ig) = mixing_store%rhoin(ispin)%cc(ig) + f_mix*res_rho(ig)
               mixing_store%rhoin_old(ispin)%cc(ig) = mixing_store%rhoin(ispin)%cc(ig)
               mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%array(ig)
            END DO

            IF (mixing_store%gmix_p) THEN
               IF (gapw) THEN
                  DO iatom = 1, natom
                     IF (mixing_store%paw(iatom)) THEN
                        mixing_store%cpc_h_lastres(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef - &
                                                                          mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        mixing_store%cpc_s_lastres(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef - &
                                                                          mixing_store%cpc_s_in(iatom, ispin)%r_coef

                        rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                                                              mixing_store%cpc_h_in(iatom, ispin)%r_coef*(1._dp - alpha)
                        rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                                                              mixing_store%cpc_s_in(iatom, ispin)%r_coef*(1._dp - alpha)

                        mixing_store%cpc_h_old(iatom, ispin)%r_coef = mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        mixing_store%cpc_s_old(iatom, ispin)%r_coef = mixing_store%cpc_s_in(iatom, ispin)%r_coef
                        mixing_store%cpc_h_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                        mixing_store%cpc_s_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
                     END IF
                  END DO
               END IF
            END IF
         ELSE

            ALLOCATE (a(nb, nb))
            a = 0.0_dp
            ALLOCATE (b(nb, nb))
            b = 0.0_dp
            ALLOCATE (c(nb))
            c = 0.0_dp
            ALLOCATE (g(nb))
            g = 0.0_dp

            delta_norm = 0.0_dp
            res_norm = 0.0_dp
            rho_norm = 0.0_dp
            DO ig = 1, ng
               mixing_store%res_buffer(ib, ispin)%cc(ig) = res_rho(ig) - mixing_store%last_res(ispin)%cc(ig)
               mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
               res_norm = res_norm + &
                          REAL(res_rho(ig), dp)*REAL(res_rho(ig), dp) + &
                          AIMAG(res_rho(ig))*AIMAG(res_rho(ig))
               delta_norm = delta_norm + &
                            REAL(mixing_store%res_buffer(ib, ispin)%cc(ig), dp)* &
                            REAL(mixing_store%res_buffer(ib, ispin)%cc(ig), dp) + &
                            AIMAG(mixing_store%res_buffer(ib, ispin)%cc(ig))* &
                            AIMAG(mixing_store%res_buffer(ib, ispin)%cc(ig))
               rho_norm = rho_norm + &
                          REAL(rho_g(ispin)%array(ig), dp)*REAL(rho_g(ispin)%array(ig), dp) + &
                          AIMAG(rho_g(ispin)%array(ig))*AIMAG(rho_g(ispin)%array(ig))
            END DO
            DO ig = 1, ng
               mixing_store%drho_buffer(ib, ispin)%cc(ig) = &
                  mixing_store%rhoin(ispin)%cc(ig) - &
                  mixing_store%rhoin_old(ispin)%cc(ig)
            END DO
            CALL para_env%sum(delta_norm)
            delta_norm = SQRT(delta_norm)
            CALL para_env%sum(res_norm)
            res_norm = SQRT(res_norm)
            CALL para_env%sum(rho_norm)
            rho_norm = SQRT(rho_norm)

            IF (res_norm > 1.E-14_dp) THEN
               mixing_store%res_buffer(ib, ispin)%cc(:) = mixing_store%res_buffer(ib, ispin)%cc(:)/delta_norm
               mixing_store%drho_buffer(ib, ispin)%cc(:) = mixing_store%drho_buffer(ib, ispin)%cc(:)/delta_norm

               IF (mixing_store%gmix_p .AND. gapw) THEN
                  DO iatom = 1, natom
                     IF (mixing_store%paw(iatom)) THEN
                        n1 = SIZE(mixing_store%cpc_s_in(iatom, ispin)%r_coef, 1)
                        n2 = SIZE(mixing_store%cpc_s_in(iatom, ispin)%r_coef, 2)
                        DO i = 1, n2
                           DO j = 1, n1
                              mixing_store%dcpc_h_in(ib, iatom, ispin)%r_coef(j, i) = &
                                 (mixing_store%cpc_h_in(iatom, ispin)%r_coef(j, i) - &
                                  mixing_store%cpc_h_old(iatom, ispin)%r_coef(j, i))/delta_norm
                              dcpc_h_res = ((rho_atom(iatom)%cpc_h(ispin)%r_coef(j, i) - &
                                             mixing_store%cpc_h_in(iatom, ispin)%r_coef(j, i)) - &
                                            mixing_store%cpc_h_lastres(iatom, ispin)%r_coef(j, i))/delta_norm
                              mixing_store%cpc_h_lastres(iatom, ispin)%r_coef(j, i) = rho_atom(iatom)%cpc_h(ispin)%r_coef(j, i) - &
                                                                                    mixing_store%cpc_h_in(iatom, ispin)%r_coef(j, i)

                              mixing_store%dcpc_s_in(ib, iatom, ispin)%r_coef(j, i) = &
                                 (mixing_store%cpc_s_in(iatom, ispin)%r_coef(j, i) - &
                                  mixing_store%cpc_s_old(iatom, ispin)%r_coef(j, i))/delta_norm
                              dcpc_s_res = ((rho_atom(iatom)%cpc_s(ispin)%r_coef(j, i) - &
                                             mixing_store%cpc_s_in(iatom, ispin)%r_coef(j, i)) - &
                                            mixing_store%cpc_s_lastres(iatom, ispin)%r_coef(j, i))/delta_norm
                              mixing_store%cpc_s_lastres(iatom, ispin)%r_coef(j, i) = rho_atom(iatom)%cpc_s(ispin)%r_coef(j, i) - &
                                                                                    mixing_store%cpc_s_in(iatom, ispin)%r_coef(j, i)

                              mixing_store%dcpc_h_in(ib, iatom, ispin)%r_coef(j, i) = &
                                 alpha*dcpc_h_res + &
                                 mixing_store%dcpc_h_in(ib, iatom, ispin)%r_coef(j, i)
                              mixing_store%dcpc_s_in(ib, iatom, ispin)%r_coef(j, i) = &
                                 alpha*dcpc_s_res + &
                                 mixing_store%dcpc_s_in(ib, iatom, ispin)%r_coef(j, i)
                           END DO
                        END DO
                     END IF
                  END DO
               END IF

               a(:, :) = 0.0_dp
               DO ig = 1, ng
                  f_mix = alpha*mixing_store%kerker_factor(ig)
                  mixing_store%drho_buffer(ib, ispin)%cc(ig) = &
                     f_mix*mixing_store%res_buffer(ib, ispin)%cc(ig) + &
                     mixing_store%drho_buffer(ib, ispin)%cc(ig)
               END DO
               DO jb = 1, nb
                  DO kb = jb, nb
                     DO ig = 1, mixing_store%ig_max !ng
                        a(kb, jb) = a(kb, jb) + mixing_store%p_metric(ig)*( &
                                    REAL(mixing_store%res_buffer(jb, ispin)%cc(ig), dp)* &
                                    REAL(mixing_store%res_buffer(kb, ispin)%cc(ig), dp) + &
                                    AIMAG(mixing_store%res_buffer(jb, ispin)%cc(ig))* &
                                    AIMAG(mixing_store%res_buffer(kb, ispin)%cc(ig)))
                     END DO
                     a(jb, kb) = a(kb, jb)
                  END DO
               END DO
               CALL para_env%sum(a)

               C = 0.0_dp
               DO jb = 1, nb
                  a(jb, jb) = w0 + a(jb, jb)
                  DO ig = 1, mixing_store%ig_max !ng
                     c(jb) = c(jb) + mixing_store%p_metric(ig)*( &
                             REAL(mixing_store%res_buffer(jb, ispin)%cc(ig), dp)*REAL(res_rho(ig), dp) + &
                             AIMAG(mixing_store%res_buffer(jb, ispin)%cc(ig))*AIMAG(res_rho(ig)))
                  END DO
               END DO
               CALL para_env%sum(c)
               CALL invert_matrix(a, b, inv_err)

               CALL dgemv('T', nb, nb, 1.0_dp, B, nb, C, 1, 0.0_dp, G, 1)
            ELSE
               G = 0.0_dp
            END IF

            DO ig = 1, ng
               cc_mix = CMPLX(0.0_dp, 0.0_dp, kind=dp)
               DO jb = 1, nb
                  cc_mix = cc_mix - G(jb)*mixing_store%drho_buffer(jb, ispin)%cc(ig)
               END DO
               f_mix = alpha*mixing_store%kerker_factor(ig)

               IF (res_norm > 1.E-14_dp) rho_g(ispin)%array(ig) = mixing_store%rhoin(ispin)%cc(ig) + &
                                                                  f_mix*res_rho(ig) + cc_mix
               mixing_store%rhoin_old(ispin)%cc(ig) = mixing_store%rhoin(ispin)%cc(ig)
               mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%array(ig)
            END DO

            IF (mixing_store%gmix_p) THEN
               IF (gapw) THEN
                  DO iatom = 1, natom
                     IF (mixing_store%paw(iatom)) THEN
                        n1 = SIZE(mixing_store%cpc_s_in(iatom, ispin)%r_coef, 1)
                        n2 = SIZE(mixing_store%cpc_s_in(iatom, ispin)%r_coef, 2)
                        DO i = 1, n2
                           DO j = 1, n1
                              valh = 0.0_dp
                              vals = 0.0_dp
                              DO jb = 1, nb
                                 valh = valh - G(jb)*mixing_store%dcpc_h_in(jb, iatom, ispin)%r_coef(j, i)
                                 vals = vals - G(jb)*mixing_store%dcpc_s_in(jb, iatom, ispin)%r_coef(j, i)
                              END DO
                              IF (res_norm > 1.E-14_dp) THEN
                                 rho_atom(iatom)%cpc_h(ispin)%r_coef(j, i) = &
                                    alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef(j, i) + &
                                    mixing_store%cpc_h_in(iatom, ispin)%r_coef(j, i)*(1._dp - alpha) + valh
                                 rho_atom(iatom)%cpc_s(ispin)%r_coef(j, i) = &
                                    alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef(j, i) + &
                                    mixing_store%cpc_s_in(iatom, ispin)%r_coef(j, i)*(1._dp - alpha) + vals
                              END IF
                           END DO
                        END DO

                        mixing_store%cpc_h_old(iatom, ispin)%r_coef = mixing_store%cpc_h_in(iatom, ispin)%r_coef
                        mixing_store%cpc_s_old(iatom, ispin)%r_coef = mixing_store%cpc_s_in(iatom, ispin)%r_coef
                        mixing_store%cpc_h_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
                        mixing_store%cpc_s_in(iatom, ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
                     END IF
                  END DO
               END IF
            END IF

            DEALLOCATE (a, b, c, g)
         END IF

      END DO ! ispin
      IF (nspin == 2) THEN
         CALL pw_axpy(rho_g(1), rho_g(2), 1.0_dp, -1.0_dp)
         DO ig = 1, ng
            mixing_store%rhoin(2)%cc(ig) = mixing_store%rhoin(1)%cc(ig) - mixing_store%rhoin(2)%cc(ig)
         END DO
         IF (gapw .AND. mixing_store%gmix_p) THEN
            DO iatom = 1, natom
               IF (mixing_store%paw(iatom)) THEN
                  rho_atom(iatom)%cpc_h(2)%r_coef = rho_atom(iatom)%cpc_h(1)%r_coef - rho_atom(iatom)%cpc_h(2)%r_coef
                  rho_atom(iatom)%cpc_s(2)%r_coef = rho_atom(iatom)%cpc_s(1)%r_coef - rho_atom(iatom)%cpc_s(2)%r_coef
                  mixing_store%cpc_h_in(iatom, 2)%r_coef = mixing_store%cpc_h_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_h_in(iatom, 2)%r_coef
                  mixing_store%cpc_s_in(iatom, 2)%r_coef = mixing_store%cpc_s_in(iatom, 1)%r_coef - &
                                                           mixing_store%cpc_s_in(iatom, 2)%r_coef
               END IF
            END DO
         END IF

      END IF

      DEALLOCATE (res_rho)

      CALL timestop(handle)

   END SUBROUTINE broyden_mixing

! **************************************************************************************************
!> \brief Multisecant scheme to perform charge density Mixing
!>        as preconditioner we use the Kerker damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      05.2009
!> \author MI
! **************************************************************************************************
   SUBROUTINE multisecant_mixing(mixing_store, rho, para_env)

      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'multisecant_mixing'
      COMPLEX(KIND=dp), PARAMETER                        :: cmone = (-1.0_dp, 0.0_dp), &
                                                            cone = (1.0_dp, 0.0_dp), &
                                                            czero = (0.0_dp, 0.0_dp)

      COMPLEX(dp)                                        :: saa, yaa
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: gn_global, pgn, res_matrix_global, &
                                                            tmp_vec, ugn
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: a_matrix, res_matrix, sa, step_matrix, ya
      COMPLEX(dp), DIMENSION(:), POINTER                 :: gn
      INTEGER                                            :: handle, ib, ib_next, ib_prev, ig, &
                                                            ig_global, iig, ispin, jb, kb, nb, &
                                                            nbuffer, ng, ng_global, nspin
      LOGICAL                                            :: use_zgemm, use_zgemm_rev
      REAL(dp) :: alpha, f_mix, gn_norm, gn_norm_old, inv_err, n_low, n_up, pgn_norm, prec, &
         r_step, reg_par, sigma_max, sigma_tilde, step_size
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: norm_res, norm_res_low, norm_res_up
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: b_matrix, binv_matrix
      REAL(dp), DIMENSION(:), POINTER                    :: g2
      REAL(dp), SAVE                                     :: sigma_old = 1.0_dp
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g

      CPASSERT(ASSOCIATED(mixing_store))

      CALL timeset(routineN, handle)

      NULLIFY (gn, rho_g)

      use_zgemm = .FALSE.
      use_zgemm_rev = .TRUE.

      ! prepare the parameters
      CALL qs_rho_get(rho, rho_g=rho_g)

      nspin = SIZE(rho_g, 1)
      ! not implemented for large grids.
      CPASSERT(rho_g(1)%pw_grid%ngpts < HUGE(ng_global))
      ng_global = INT(rho_g(1)%pw_grid%ngpts)
      ng = SIZE(mixing_store%rhoin_buffer(1, 1)%cc)
      alpha = mixing_store%alpha

      sigma_max = mixing_store%sigma_max
      reg_par = mixing_store%reg_par
      r_step = mixing_store%r_step
      nbuffer = mixing_store%nbuffer

      ! determine the step number, and multisecant iteration
      nb = MIN(mixing_store%ncall, nbuffer - 1)
      ib = MODULO(mixing_store%ncall, nbuffer) + 1
      IF (mixing_store%ncall > 0) THEN
         ib_prev = MODULO(mixing_store%ncall - 1, nbuffer) + 1
      ELSE
         ib_prev = 0
      END IF
      mixing_store%ncall = mixing_store%ncall + 1
      ib_next = MODULO(mixing_store%ncall, nbuffer) + 1

      ! compute the residual gn and its norm gn_norm
      DO ispin = 1, nspin
         gn => mixing_store%res_buffer(ib, ispin)%cc
         gn_norm = 0.0_dp
         DO ig = 1, ng
            gn(ig) = (rho_g(ispin)%array(ig) - mixing_store%rhoin_buffer(ib, ispin)%cc(ig))
            gn_norm = gn_norm + &
                      REAL(gn(ig), dp)*REAL(gn(ig), dp) + AIMAG(gn(ig))*AIMAG(gn(ig))
         END DO
         CALL para_env%sum(gn_norm)
         gn_norm = SQRT(gn_norm)
         mixing_store%norm_res_buffer(ib, ispin) = gn_norm
      END DO

      IF (nb == 0) THEN
         !simple mixing
         DO ispin = 1, nspin
            DO ig = 1, ng
               f_mix = alpha*mixing_store%kerker_factor(ig)
               rho_g(ispin)%array(ig) = mixing_store%rhoin_buffer(1, ispin)%cc(ig) + &
                                        f_mix*mixing_store%res_buffer(1, ispin)%cc(ig)
               mixing_store%rhoin_buffer(ib_next, ispin)%cc(ig) = rho_g(ispin)%array(ig)
            END DO
         END DO
         CALL timestop(handle)
         RETURN
      END IF

      ! allocate temporary arrays
      ! step_matrix  S ngxnb
      ALLOCATE (step_matrix(ng, nb))
      ! res_matrix Y  ngxnb
      ALLOCATE (res_matrix(ng, nb))
      ! matrix A  nbxnb
      ALLOCATE (a_matrix(nb, ng_global))
      ! PSI nb vector of norms
      ALLOCATE (norm_res(nb))
      ALLOCATE (norm_res_low(nb))
      ALLOCATE (norm_res_up(nb))

      ! matrix B   nbxnb
      ALLOCATE (b_matrix(nb, nb))
      ! matrix B_inv   nbxnb
      ALLOCATE (binv_matrix(nb, nb))

      ALLOCATE (gn_global(ng_global))
      ALLOCATE (res_matrix_global(ng_global))
      IF (use_zgemm) THEN
         ALLOCATE (sa(ng, ng_global))
         ALLOCATE (ya(ng, ng_global))
      END IF
      IF (use_zgemm_rev) THEN
         ALLOCATE (tmp_vec(nb))
      END IF
      ALLOCATE (pgn(ng))
      ALLOCATE (ugn(ng))

      DO ispin = 1, nspin
         ! generate the global vector with the present residual
         gn => mixing_store%res_buffer(ib, ispin)%cc
         gn_global = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ig = 1, ng
            ig_global = mixing_store%ig_global_index(ig)
            gn_global(ig_global) = gn(ig)
         END DO
         CALL para_env%sum(gn_global)

         ! compute steps (matrix S) and residual differences (matrix Y)
         ! with respect to the present
         ! step and the present residual (use stored rho_in and res_buffer)

         ! These quantities are pre-conditioned by means of Kerker factor multipliccation
         g2 => rho_g(1)%pw_grid%gsq

         DO jb = 1, nb + 1
            IF (jb < ib) THEN
               kb = jb
            ELSEIF (jb > ib) THEN
               kb = jb - 1
            ELSE
               CYCLE
            END IF
            norm_res(kb) = 0.0_dp
            norm_res_low(kb) = 0.0_dp
            norm_res_up(kb) = 0.0_dp
            DO ig = 1, ng

               prec = 1.0_dp

               step_matrix(ig, kb) = prec*(mixing_store%rhoin_buffer(jb, ispin)%cc(ig) - &
                                           mixing_store%rhoin_buffer(ib, ispin)%cc(ig))
               res_matrix(ig, kb) = (mixing_store%res_buffer(jb, ispin)%cc(ig) - &
                                     mixing_store%res_buffer(ib, ispin)%cc(ig))
               norm_res(kb) = norm_res(kb) + REAL(res_matrix(ig, kb), dp)*REAL(res_matrix(ig, kb), dp) + &
                              AIMAG(res_matrix(ig, kb))*AIMAG(res_matrix(ig, kb))
               IF (g2(ig) < 4.0_dp) THEN
                  norm_res_low(kb) = norm_res_low(kb) + &
                                     REAL(res_matrix(ig, kb), dp)*REAL(res_matrix(ig, kb), dp) + &
                                     AIMAG(res_matrix(ig, kb))*AIMAG(res_matrix(ig, kb))
               ELSE
                  norm_res_up(kb) = norm_res_up(kb) + &
                                    REAL(res_matrix(ig, kb), dp)*REAL(res_matrix(ig, kb), dp) + &
                                    AIMAG(res_matrix(ig, kb))*AIMAG(res_matrix(ig, kb))
               END IF
               res_matrix(ig, kb) = prec*res_matrix(ig, kb)
            END DO
         END DO !jb

         ! normalize each column of S and Y => Snorm Ynorm
         CALL para_env%sum(norm_res)
         CALL para_env%sum(norm_res_up)
         CALL para_env%sum(norm_res_low)
         norm_res(1:nb) = 1.0_dp/SQRT(norm_res(1:nb))
         n_low = 0.0_dp
         n_up = 0.0_dp
         DO jb = 1, nb
            n_low = n_low + norm_res_low(jb)/norm_res(jb)
            n_up = n_up + norm_res_up(jb)/norm_res(jb)
         END DO
         DO ig = 1, ng
            IF (g2(ig) > 4.0_dp) THEN
               step_matrix(ig, 1:nb) = step_matrix(ig, 1:nb)*SQRT(n_low/n_up)
               res_matrix(ig, 1:nb) = res_matrix(ig, 1:nb)*SQRT(n_low/n_up)
            END IF
         END DO
         DO kb = 1, nb
            step_matrix(1:ng, kb) = step_matrix(1:ng, kb)*norm_res(kb)
            res_matrix(1:ng, kb) = res_matrix(1:ng, kb)*norm_res(kb)
         END DO

         ! compute A as [(Ynorm^t Ynorm) + (alpha I)]^(-1) Ynorm^t
         ! compute B
         DO jb = 1, nb
            DO kb = 1, nb
               b_matrix(kb, jb) = 0.0_dp
               DO ig = 1, ng
                  ! it is assumed that summing over all G vector gives a real, because
                  !  y(-G,kb) = (y(G,kb))*
                  b_matrix(kb, jb) = b_matrix(kb, jb) + REAL(res_matrix(ig, kb)*res_matrix(ig, jb), dp)
               END DO
            END DO
         END DO

         CALL para_env%sum(b_matrix)
         DO jb = 1, nb
            b_matrix(jb, jb) = b_matrix(jb, jb) + reg_par
         END DO
         ! invert B
         CALL invert_matrix(b_matrix, binv_matrix, inv_err)

         ! A = Binv Ynorm^t
         a_matrix = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO kb = 1, nb
            res_matrix_global = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO ig = 1, ng
               ig_global = mixing_store%ig_global_index(ig)
               res_matrix_global(ig_global) = res_matrix(ig, kb)
            END DO
            CALL para_env%sum(res_matrix_global)

            DO jb = 1, nb
               DO ig = 1, ng
                  ig_global = mixing_store%ig_global_index(ig)
                  a_matrix(jb, ig_global) = a_matrix(jb, ig_global) + &
                                            binv_matrix(jb, kb)*res_matrix_global(ig_global)
               END DO
            END DO
         END DO
         CALL para_env%sum(a_matrix)

         ! compute the two components of gn that will be used to update rho
         gn => mixing_store%res_buffer(ib, ispin)%cc
         pgn_norm = 0.0_dp

         IF (use_zgemm) THEN

            CALL zgemm("N", "N", ng, ng_global, nb, cmone, step_matrix(1, 1), ng, &
                       a_matrix(1, 1), nb, czero, sa(1, 1), ng)
            CALL zgemm("N", "N", ng, ng_global, nb, cmone, res_matrix(1, 1), ng, &
                       a_matrix(1, 1), nb, czero, ya(1, 1), ng)
            DO ig = 1, ng
               ig_global = mixing_store%ig_global_index(ig)
               ya(ig, ig_global) = ya(ig, ig_global) + CMPLX(1.0_dp, 0.0_dp, KIND=dp)
            END DO

            CALL zgemv("N", ng, ng_global, cone, sa(1, 1), &
                       ng, gn_global(1), 1, czero, pgn(1), 1)
            CALL zgemv("N", ng, ng_global, cone, ya(1, 1), &
                       ng, gn_global(1), 1, czero, ugn(1), 1)

            DO ig = 1, ng
               pgn_norm = pgn_norm + REAL(pgn(ig), dp)*REAL(pgn(ig), dp) + &
                          AIMAG(pgn(ig))*AIMAG(pgn(ig))
            END DO
            CALL para_env%sum(pgn_norm)
         ELSEIF (use_zgemm_rev) THEN

            CALL zgemv("N", nb, ng_global, cone, a_matrix(1, 1), &
                       nb, gn_global(1), 1, czero, tmp_vec(1), 1)

            CALL zgemv("N", ng, nb, cmone, step_matrix(1, 1), ng, &
                       tmp_vec(1), 1, czero, pgn(1), 1)

            CALL zgemv("N", ng, nb, cmone, res_matrix(1, 1), ng, &
                       tmp_vec(1), 1, czero, ugn(1), 1)

            DO ig = 1, ng
               pgn_norm = pgn_norm + REAL(pgn(ig), dp)*REAL(pgn(ig), dp) + &
                          AIMAG(pgn(ig))*AIMAG(pgn(ig))
               ugn(ig) = ugn(ig) + gn(ig)
            END DO
            CALL para_env%sum(pgn_norm)

         ELSE
            DO ig = 1, ng
               pgn(ig) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               ugn(ig) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
               ig_global = mixing_store%ig_global_index(ig)
               DO iig = 1, ng_global
                  saa = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
                  yaa = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

                  IF (ig_global == iig) yaa = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

                  DO jb = 1, nb
                     saa = saa - step_matrix(ig, jb)*a_matrix(jb, iig)
                     yaa = yaa - res_matrix(ig, jb)*a_matrix(jb, iig)
                  END DO
                  pgn(ig) = pgn(ig) + saa*gn_global(iig)
                  ugn(ig) = ugn(ig) + yaa*gn_global(iig)
               END DO
            END DO
            DO ig = 1, ng
               pgn_norm = pgn_norm + REAL(pgn(ig), dp)*REAL(pgn(ig), dp) + &
                          AIMAG(pgn(ig))*AIMAG(pgn(ig))
            END DO
            CALL para_env%sum(pgn_norm)
         END IF

         gn_norm = mixing_store%norm_res_buffer(ib, ispin)
         gn_norm_old = mixing_store%norm_res_buffer(ib_prev, ispin)
         IF (ib_prev /= 0) THEN
            sigma_tilde = sigma_old*MAX(0.5_dp, MIN(2.0_dp, gn_norm_old/gn_norm))
         ELSE
            sigma_tilde = 0.5_dp
         END IF
         sigma_tilde = 0.1_dp
         ! Step size for the unpredicted component
         step_size = MIN(sigma_tilde, r_step*pgn_norm/gn_norm, sigma_max)
         sigma_old = step_size

         ! update the density
         DO ig = 1, ng
            prec = mixing_store%kerker_factor(ig)
            rho_g(ispin)%array(ig) = mixing_store%rhoin_buffer(ib, ispin)%cc(ig) &
                                     - prec*step_size*ugn(ig) + prec*pgn(ig) ! - 0.1_dp * prec* gn(ig)
            mixing_store%rhoin_buffer(ib_next, ispin)%cc(ig) = rho_g(ispin)%array(ig)
         END DO

      END DO ! ispin

      ! Deallocate  temporary arrays
      DEALLOCATE (step_matrix, res_matrix)
      DEALLOCATE (norm_res)
      DEALLOCATE (norm_res_up)
      DEALLOCATE (norm_res_low)
      DEALLOCATE (a_matrix, b_matrix, binv_matrix)
      DEALLOCATE (ugn, pgn)
      IF (use_zgemm) THEN
         DEALLOCATE (sa, ya)
      END IF
      IF (use_zgemm_rev) THEN
         DEALLOCATE (tmp_vec)
      END IF
      DEALLOCATE (gn_global, res_matrix_global)

      CALL timestop(handle)

   END SUBROUTINE multisecant_mixing

END MODULE qs_gspace_mixing
