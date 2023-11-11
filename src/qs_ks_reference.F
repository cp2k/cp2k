!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the KS reference potentials
!> \par History
!>       07.2022 created
!> \author JGH
! **************************************************************************************************
MODULE qs_ks_reference
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals,&
                                              init_coulomb_local
   USE hartree_local_types,             ONLY: hartree_local_create,&
                                              hartree_local_release,&
                                              hartree_local_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace,&
                                              rho0_s_grid_create
   USE qs_rho0_methods,                 ONLY: init_rho0
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom_coeff
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_reference'

   PUBLIC :: ks_ref_potential, ks_ref_potential_atom

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief calculate the Kohn-Sham reference potential
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param vadmm_rspace ...
!> \param ehartree ...
!> \param exc ...
!> \param h_stress container for the stress tensor of the Hartree term
!> \par History
!>      10.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, &
                               ehartree, exc, h_stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_type), INTENT(INOUT)                       :: vh_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: vxc_rspace, vtau_rspace, vadmm_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: ehartree, exc
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ks_ref_potential'

      INTEGER                                            :: handle, iab, ispin, nspins
      REAL(dp)                                           :: eadmm, eovrl, eself
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: v_admm_rspace, v_admm_tau_rspace, &
                                                            v_rspace, v_tau_rspace
      TYPE(pw_type), POINTER                             :: rho_core
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_xc
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, dft_control=dft_control, &
                      para_env=para_env, ks_env=ks_env, rho_core=rho_core)

      nspins = dft_control%nspins

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      CALL auxbas_pw_pool%create_pw(v_hartree_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(v_hartree_rspace, &
                                    use_data=REALDATA3D, in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(rho_tot_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace, ehartree, &
                            v_hartree_gspace, h_stress=h_stress, rho_core=rho_core)
      CALL pw_transfer(v_hartree_gspace, v_hartree_rspace)
      CALL pw_scale(v_hartree_rspace, v_hartree_rspace%pw_grid%dvol)

      CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
      CALL auxbas_pw_pool%give_back_pw(rho_tot_gspace)
      !
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, PRESENT(h_stress), E_overlap_core=eovrl)
      ehartree = ehartree + eovrl + eself

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      NULLIFY (v_rspace, v_tau_rspace)
      IF (dft_control%qs_control%gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, rho_xc=rho_xc)
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_xc, xc_section=xc_section, &
                            vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
      ELSE
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                            vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
      END IF

      NULLIFY (v_admm_rspace, v_admm_tau_rspace)
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            ! For the virial, we have to save the pv_xc component because it will be reset in qs_vxc_create
            IF (PRESENT(h_stress)) THEN
               CALL get_qs_env(qs_env, virial=virial)
               virial_xc = virial%pv_xc
            END IF
            CALL get_admm_env(admm_env, rho_aux_fit=rho)
            xc_section => admm_env%xc_section_aux
            CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                               vxc_rho=v_admm_rspace, vxc_tau=v_admm_tau_rspace, exc=eadmm, just_energy=.FALSE.)
            IF (PRESENT(h_stress)) virial%pv_xc = virial%pv_xc + virial_xc
         END IF
      END IF

      ! allocate potentials
      IF (ASSOCIATED(vh_rspace%pw_grid)) THEN
         CALL vh_rspace%release()
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL vxc_rspace(iab)%release()
         END DO
      ELSE
         ALLOCATE (vxc_rspace(nspins))
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO iab = 1, SIZE(vtau_rspace)
               CALL vtau_rspace(iab)%release()
            END DO
         ELSE
            ALLOCATE (vtau_rspace(nspins))
         END IF
      ELSE
         NULLIFY (vtau_rspace)
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         IF (ASSOCIATED(vadmm_rspace)) THEN
            DO iab = 1, SIZE(vadmm_rspace)
               CALL vadmm_rspace(iab)%release()
            END DO
         ELSE
            ALLOCATE (vadmm_rspace(nspins))
         END IF
      ELSE
         NULLIFY (vadmm_rspace)
      END IF

      pw_grid => v_hartree_rspace%pw_grid
      CALL vh_rspace%create(pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      DO ispin = 1, nspins
         CALL vxc_rspace(ispin)%create(pw_grid, &
                                       use_data=REALDATA3D, in_space=REALSPACE)
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL vtau_rspace(ispin)%create(pw_grid, &
                                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
         IF (ASSOCIATED(vadmm_rspace)) THEN
            CALL vadmm_rspace(ispin)%create(pw_grid, &
                                            use_data=REALDATA3D, in_space=REALSPACE)
         END IF
      END DO
      !
      CALL pw_transfer(v_hartree_rspace, vh_rspace)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_rspace(ispin), vxc_rspace(ispin))
            CALL pw_scale(vxc_rspace(ispin), v_rspace(ispin)%pw_grid%dvol)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_transfer(v_tau_rspace(ispin), vtau_rspace(ispin))
               CALL pw_scale(vtau_rspace(ispin), v_tau_rspace(ispin)%pw_grid%dvol)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(vxc_rspace(ispin))
         END DO
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_admm_rspace(ispin), vadmm_rspace(ispin))
            CALL pw_scale(vadmm_rspace(ispin), vadmm_rspace(ispin)%pw_grid%dvol)
         END DO
      END IF

      ! return pw grids
      CALL auxbas_pw_pool%give_back_pw(v_hartree_rspace)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_rspace(ispin))
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL auxbas_pw_pool%give_back_pw(v_tau_rspace(ispin))
            END IF
         END DO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_admm_rspace(ispin))
         END DO
         DEALLOCATE (v_admm_rspace)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential

! **************************************************************************************************
!> \brief calculate the Kohn-Sham GAPW reference potentials
!> \param qs_env ...
!> \param local_rho_set ...
!> \param local_rho_set_admm ...
!> \param v_hartree_rspace ...
!> \par History
!>      07.2022 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential_atom(qs_env, local_rho_set, local_rho_set_admm, v_hartree_rspace)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(local_rho_type), POINTER                      :: local_rho_set, local_rho_set_admm
      TYPE(pw_type), INTENT(IN)                          :: v_hartree_rspace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ks_ref_potential_atom'

      INTEGER                                            :: handle, natom, nspins
      LOGICAL                                            :: gapw, gapw_xc
      REAL(KIND=dp)                                      :: eh1c, exc1, exc1_admm
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_aux, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hartree_local_type), POINTER                  :: hartree_local
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      CALL get_qs_env(qs_env=qs_env, rho=rho, pw_env=pw_env, &
                      dft_control=dft_control, para_env=para_env)

      nspins = dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      IF (gapw .OR. gapw_xc) THEN
         NULLIFY (hartree_local, local_rho_set, local_rho_set_admm)
         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set)
         CALL local_rho_set_create(local_rho_set)
         CALL allocate_rho_atom_internals(local_rho_set%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         IF (gapw) THEN
            CALL get_qs_env(qs_env, natom=natom)
            CALL init_rho0(local_rho_set, qs_env, dft_control%qs_control%gapw_control)
            CALL rho0_s_grid_create(pw_env, local_rho_set%rho0_mpole)
            CALL hartree_local_create(hartree_local)
            CALL init_coulomb_local(hartree_local, natom)
         END IF

         CALL get_qs_env(qs_env=qs_env, oce=oce, sab_orb=sab)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         CALL calculate_rho_atom_coeff(qs_env, rho_ao_kp, local_rho_set%rho_atom_set, &
                                       qs_kind_set, oce, sab, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set, do_rho0=gapw)

         IF (gapw) THEN
            CALL Vh_1c_gg_integrals(qs_env, eh1c, hartree_local%ecoul_1c, local_rho_set, para_env, .FALSE.)
            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, calculate_forces=.FALSE., &
                                       local_rho_set=local_rho_set)
         END IF
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            xc_section => admm_env%xc_section_primary
         ELSE
            xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
         END IF
         CALL calculate_vxc_atom(qs_env, .FALSE., exc1=exc1, xc_section_external=xc_section, &
                                 rho_atom_set_external=local_rho_set%rho_atom_set)

         IF (dft_control%do_admm) THEN
            IF (admm_env%do_gapw) THEN
               CALL local_rho_set_create(local_rho_set_admm)
               CALL allocate_rho_atom_internals(local_rho_set_admm%rho_atom_set, atomic_kind_set, &
                                                admm_env%admm_gapw_env%admm_kind_set, dft_control, para_env)
               oce => admm_env%admm_gapw_env%oce
               sab => admm_env%sab_aux_fit
               CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit)
               CALL qs_rho_get(rho, rho_ao_kp=rho_ao_aux)
               CALL calculate_rho_atom_coeff(qs_env, rho_ao_aux, local_rho_set_admm%rho_atom_set, &
                                             admm_env%admm_gapw_env%admm_kind_set, oce, sab, para_env)
               CALL prepare_gapw_den(qs_env, local_rho_set=local_rho_set_admm, &
                                     do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
               !compute the potential due to atomic densities
               xc_section => admm_env%xc_section_aux
               CALL calculate_vxc_atom(qs_env, energy_only=.FALSE., exc1=exc1_admm, &
                                       kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                                       xc_section_external=xc_section, &
                                       rho_atom_set_external=local_rho_set_admm%rho_atom_set)
            END IF
         END IF

         ! clean up
         CALL hartree_local_release(hartree_local)

      ELSE

         NULLIFY (local_rho_set, local_rho_set_admm)

      END IF

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential_atom

! **************************************************************************************************

END MODULE qs_ks_reference
