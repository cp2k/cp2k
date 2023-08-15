!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      none
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE atprop_types,                    ONLY: atprop_array_add,&
                                              atprop_array_init,&
                                              atprop_type
   USE bandstructure_methods,           ONLY: post_SCF_bandstructure
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_control_utils,                ONLY: read_ddapc_section
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set
   USE et_coupling,                     ONLY: calc_et_coupling
   USE et_coupling_proj,                ONLY: calc_et_coupling_proj
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE hartree_local_types,             ONLY: ecoul_1center_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mulliken,                        ONLY: atom_trace
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_module,                ONLY: linres_calculation_low
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_atom_types,               ONLY: rho_atom_type,&
                                              zero_rho_atom_integrals
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf,                          ONLY: scf
   USE qs_tddfpt2_methods,              ONLY: tddfpt
   USE qs_vxc,                          ONLY: qs_xc_density
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
   USE tip_scan_methods,                ONLY: tip_scanning
   USE xas_methods,                     ONLY: xas
   USE xas_tdp_methods,                 ONLY: xas_tdp
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_utils'

   PUBLIC :: qs_energies_properties

CONTAINS

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of properties
!>        into separate subroutine
!> \param qs_env ...
!> \param calc_forces ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_properties(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_properties'

      INTEGER                                            :: handle, natom
      LOGICAL                                            :: do_et, do_et_proj, &
                                                            do_post_SCF_bandstructure, do_tip_scan
      REAL(KIND=dp)                                      :: ekts
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, post_SCF_bandstr_section, &
                                                            proj_section, rest_b_section, &
                                                            tip_scan_section

      NULLIFY (atprop, energy, pw_env)
      CALL timeset(routineN, handle)

      ! atomic energies using Mulliken partition
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      atprop=atprop, &
                      energy=energy, &
                      v_hartree_rspace=v_hartree_rspace, &
                      para_env=para_env, &
                      pw_env=pw_env)
      IF (atprop%energy) THEN
         CALL qs_energies_mulliken(qs_env)
         CALL get_qs_env(qs_env, natom=natom)
         IF (.NOT. dft_control%qs_control%semi_empirical .AND. &
             .NOT. dft_control%qs_control%xtb .AND. &
             .NOT. dft_control%qs_control%dftb) THEN
            ! Nuclear charge correction
            CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)
            IF (.NOT. ASSOCIATED(atprop%ateb)) THEN
               CALL atprop_array_init(atprop%ateb, natom)
            END IF
            ! Kohn-Sham Functional corrections
            CALL ks_xc_correction(qs_env)
         END IF
         CALL atprop_array_add(atprop%atener, atprop%ateb)
         CALL atprop_array_add(atprop%atener, atprop%ateself)
         CALL atprop_array_add(atprop%atener, atprop%atexc)
         CALL atprop_array_add(atprop%atener, atprop%atecoul)
         CALL atprop_array_add(atprop%atener, atprop%atevdw)
         CALL atprop_array_add(atprop%atener, atprop%ategcp)
         CALL atprop_array_add(atprop%atener, atprop%atecc)
         CALL atprop_array_add(atprop%atener, atprop%ate1c)
         ! entropic energy
         ekts = energy%kts/REAL(natom, KIND=dp)/REAL(para_env%num_pe, KIND=dp)
         atprop%atener(:) = atprop%atener(:) + ekts
      END IF

      ! ET coupling - projection-operator approach
      NULLIFY (proj_section)
      proj_section => &
         section_vals_get_subs_vals(input, "PROPERTIES%ET_COUPLING%PROJECTION")
      CALL section_vals_get(proj_section, explicit=do_et_proj)
      IF (do_et_proj) THEN
         CALL calc_et_coupling_proj(qs_env)
      END IF

      ! **********  Calculate the electron transfer coupling elements********
      do_et = .FALSE.
      do_et = dft_control%qs_control%et_coupling_calc
      IF (do_et) THEN
         qs_env%et_coupling%energy = energy%total
         qs_env%et_coupling%keep_matrix = .TRUE.
         qs_env%et_coupling%first_run = .TRUE.
         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
         qs_env%et_coupling%first_run = .FALSE.
         IF (dft_control%qs_control%ddapc_restraint) THEN
            rest_b_section => section_vals_get_subs_vals(input, "PROPERTIES%ET_COUPLING%DDAPC_RESTRAINT_B")
            CALL read_ddapc_section(qs_control=dft_control%qs_control, &
                                    ddapc_restraint_section=rest_b_section)
         END IF
         CALL scf(qs_env=qs_env)
         qs_env%et_coupling%keep_matrix = .TRUE.

         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
         CALL calc_et_coupling(qs_env)
      END IF

      !Properties
      IF (dft_control%do_xas_calculation) THEN
         CALL xas(qs_env, dft_control)
      END IF

      IF (dft_control%do_xas_tdp_calculation) THEN
         CALL xas_tdp(qs_env)
      END IF

      ! Compute Linear Response properties as post-scf
      IF (.NOT. qs_env%linres_run) THEN
         CALL linres_calculation_low(qs_env)
      END IF

      IF (dft_control%tddfpt2_control%enabled) THEN
         CALL tddfpt(qs_env, calc_forces)
      END IF

      ! post-SCF bandstructure calculation from higher level methods
      NULLIFY (post_SCF_bandstr_section)
      post_SCF_bandstr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%BANDSTRUCTURE")
      CALL section_vals_get(post_SCF_bandstr_section, explicit=do_post_SCF_bandstructure)
      IF (do_post_SCF_bandstructure) THEN
         CALL post_SCF_bandstructure(qs_env)
      END IF

      ! tip scan
      NULLIFY (tip_scan_section)
      tip_scan_section => section_vals_get_subs_vals(input, "PROPERTIES%TIP_SCAN")
      CALL section_vals_get(tip_scan_section, explicit=do_tip_scan)
      IF (do_tip_scan) THEN
         CALL tip_scanning(qs_env, tip_scan_section)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies_properties

! **************************************************************************************************
!> \brief   Use a simple Mulliken-like energy decomposition
!> \param qs_env ...
!> \date    07.2011
!> \author  JHU
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE qs_energies_mulliken(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ispin
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_ks, rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (atprop, matrix_h, matrix_ks, rho, rho_ao)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_h=matrix_h, &
                      rho=rho, atprop=atprop)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (atprop%energy) THEN
         ! E = 0.5*Tr(H*P+F*P)
         atprop%atener = 0._dp
         DO ispin = 1, SIZE(rho_ao)
            CALL atom_trace(matrix_h(1)%matrix, rho_ao(ispin)%matrix, &
                            0.5_dp, atprop%atener)
            CALL atom_trace(matrix_ks(ispin)%matrix, rho_ao(ispin)%matrix, &
                            0.5_dp, atprop%atener)
         END DO
      END IF

   END SUBROUTINE qs_energies_mulliken

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE ks_xc_correction(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ks_xc_correction'

      INTEGER                                            :: handle, iatom, ispin, natom, nspins
      LOGICAL                                            :: gapw, gapw_xc
      REAL(KIND=dp)                                      :: eh1, exc1
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao, xcmat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ecoul_1center_type), DIMENSION(:), POINTER    :: ecoul_1c
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: xc_den
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:)           :: vtau, vxc
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_fun_section, xc_section
      TYPE(xc_rho_cflags_type)                           :: needs

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, pw_env=pw_env, atprop=atprop)

      IF (atprop%energy) THEN

         nspins = dft_control%nspins
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
         xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         needs = xc_functionals_get_needs(xc_fun_section, (nspins == 2), .TRUE.)
         gapw = dft_control%qs_control%gapw
         gapw_xc = dft_control%qs_control%gapw_xc

         ! Nuclear charge correction
         CALL get_qs_env(qs_env, v_hartree_rspace=v_hartree_rspace)
         IF (gapw .OR. gapw_xc) THEN
            CALL get_qs_env(qs_env=qs_env, local_rho_set=local_rho_set, &
                            rho_atom_set=rho_atom_set, ecoul_1c=ecoul_1c, &
                            natom=natom, para_env=para_env)
            CALL zero_rho_atom_integrals(rho_atom_set)
            CALL calculate_vxc_atom(qs_env, .FALSE., exc1)
            IF (gapw) THEN
               CALL Vh_1c_gg_integrals(qs_env, eh1, ecoul_1c, local_rho_set, para_env, tddft=.FALSE.)
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
               CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, calculate_forces=.FALSE., &
                                          local_rho_set=local_rho_set, atener=atprop%ateb)
            END IF
         END IF

         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, xc_den, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         ALLOCATE (vxc(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, vxc(ispin), &
                                   use_data=REALDATA3D, in_space=REALSPACE)
         END DO
         IF (needs%tau .OR. needs%tau_spin) THEN
            ALLOCATE (vtau(nspins))
            DO ispin = 1, nspins
               CALL pw_pool_create_pw(auxbas_pw_pool, vtau(ispin), &
                                      use_data=REALDATA3D, in_space=REALSPACE)
            END DO
         END IF

         IF (gapw_xc) THEN
            CALL get_qs_env(qs_env, rho_xc=rho_struct, dispersion_env=dispersion_env)
         ELSE
            CALL get_qs_env(qs_env, rho=rho_struct, dispersion_env=dispersion_env)
         END IF
         IF (needs%tau .OR. needs%tau_spin) THEN
            CALL qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env=dispersion_env, &
                               xc_den=xc_den, vxc=vxc, vtau=vtau)
         ELSE
            CALL qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env=dispersion_env, &
                               xc_den=xc_den, vxc=vxc)
         END IF
         CALL get_qs_env(qs_env, rho=rho_struct)
         CALL qs_rho_get(rho_struct, rho_ao=rho_ao)
         CALL get_qs_env(qs_env, natom=natom, matrix_s=matrix_s)
         CALL atprop_array_init(atprop%atexc, natom)
         ALLOCATE (xcmat(nspins))
         DO ispin = 1, nspins
            ALLOCATE (xcmat(ispin)%matrix)
            CALL dbcsr_create(xcmat(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_copy(xcmat(ispin)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(xcmat(ispin)%matrix, 0.0_dp)
            CALL pw_scale(vxc(ispin), -0.5_dp)
            CALL pw_axpy(xc_den, vxc(ispin))
            CALL pw_scale(vxc(ispin), vxc(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=vxc(ispin), hmat=xcmat(ispin), &
                                    calculate_forces=.FALSE., gapw=(gapw .OR. gapw_xc))
            IF (needs%tau .OR. needs%tau_spin) THEN
               CALL pw_scale(vtau(ispin), -0.5_dp*vtau(ispin)%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=vtau(ispin), &
                                       hmat=xcmat(ispin), calculate_forces=.FALSE., &
                                       gapw=(gapw .OR. gapw_xc), compute_tau=.TRUE.)
            END IF
         END DO
         IF (gapw .OR. gapw_xc) THEN
            ! remove one-center potential matrix part
            CALL qs_rho_get(rho_struct, rho_ao_kp=matrix_p)
            CALL update_ks_atom(qs_env, xcmat, matrix_p, forces=.FALSE., kscale=-0.5_dp)
            CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom_set)
            CALL atprop_array_init(atprop%ate1c, natom)
            atprop%ate1c = 0.0_dp
            DO iatom = 1, natom
               atprop%ate1c(iatom) = atprop%ate1c(iatom) + &
                                     rho_atom_set(iatom)%exc_h - rho_atom_set(iatom)%exc_s
            END DO
            IF (gapw) THEN
               CALL get_qs_env(qs_env=qs_env, ecoul_1c=ecoul_1c)
               DO iatom = 1, natom
                  atprop%ate1c(iatom) = atprop%ate1c(iatom) + &
                                        ecoul_1c(iatom)%ecoul_1_h - ecoul_1c(iatom)%ecoul_1_s + &
                                        ecoul_1c(iatom)%ecoul_1_z - ecoul_1c(iatom)%ecoul_1_0
               END DO
            END IF
         END IF
         DO ispin = 1, nspins
            CALL atom_trace(xcmat(ispin)%matrix, rho_ao(ispin)%matrix, 1.0_dp, atprop%atexc)
            CALL dbcsr_release(xcmat(ispin)%matrix)
            DEALLOCATE (xcmat(ispin)%matrix)
         END DO
         DEALLOCATE (xcmat)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, xc_den)
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc(ispin))
         END DO
         IF (needs%tau .OR. needs%tau_spin) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, vtau(ispin))
            END DO
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE ks_xc_correction

END MODULE qs_energy_utils
