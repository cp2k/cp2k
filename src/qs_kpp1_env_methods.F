!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief module that builds the second order perturbation kernel
!>      kpp1 = delta_rho|_P delta_rho|_P E drho(P1) drho
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_kpp1_env_methods
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_to_string
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none
   USE input_section_types,             ONLY: section_vals_type
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc,                          ONLY: qs_fxc_apply,&
                                              qs_fxc_prep
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_methods'

   PUBLIC :: kpp1_create, kpp1_check_i_alloc, calc_kpp1

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes a kpp1_env
!> \param kpp1_env the environment to initialize
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE kpp1_create(kpp1_env)
      TYPE(qs_kpp1_env_type)                             :: kpp1_env

      NULLIFY (kpp1_env%v_ao)

   END SUBROUTINE kpp1_create

! **************************************************************************************************
!> \brief ...
!> \param rho1_xc ...
!> \param rho1 ...
!> \param xc_section ...
!> \param lrigpw ...
!> \param qs_env ...
!> \param p_env ...
!> \param calc_forces ...
!> \param calc_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calc_kpp1(rho1_xc, rho1, xc_section, lrigpw, qs_env, p_env, &
                        calc_forces, calc_virial, virial)

      TYPE(qs_rho_type), POINTER                         :: rho1_xc, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lrigpw
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_forces, calc_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_kpp1'

      INTEGER                                            :: handle, ikind, ispin, nkind, ns, nspins
      LOGICAL                                            :: do_onecenter, gapw, gapw_xc, &
                                                            my_calc_forces
      REAL(KIND=dp)                                      :: alpha, energy_hartree, energy_hartree_1c
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k1mat, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, psmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type)                               :: rho1_tot_gspace
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho1_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: v_hartree_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_rspace_new, v_xc, v_xc_tau
      TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho0_fxc, rho1_fxc
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho1_atom_set, rho_atom_set

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(p_env%kpp1))
      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      CPASSERT(ASSOCIATED(rho1))

      my_calc_forces = .FALSE.
      IF (PRESENT(calc_forces)) my_calc_forces = calc_forces

      CALL get_qs_env(qs_env, &
                      pw_env=pw_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      rho=rho)

      IF (lrigpw) THEN
         CALL get_qs_env(qs_env, &
                         lri_env=lri_env, &
                         lri_density=lri_density, &
                         atomic_kind_set=atomic_kind_set)
      END IF

      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      IF (gapw_xc) THEN
         CPASSERT(ASSOCIATED(rho1_xc))
      END IF
      do_onecenter = gapw .OR. gapw_xc
      nspins = dft_control%nspins

      kpp1_env => p_env%kpp1_env
      CALL kpp1_check_i_alloc(kpp1_env, qs_env, xc_section)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho1, rho_g=rho1_g)

      ! gets the tmp grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL auxbas_pw_pool%create_pw(v_hartree_rspace)

      IF (gapw .OR. gapw_xc) THEN
         CALL prepare_gapw_den(qs_env, p_env%local_rho_set, do_rho0=(.NOT. gapw_xc))
      END IF

      ! *** calculate the hartree potential on the total density ***
      CALL auxbas_pw_pool%create_pw(rho1_tot_gspace)
      CALL pw_zero(rho1_tot_gspace)
      DO ispin = 1, nspins
         CALL pw_axpy(rho1_g(ispin), rho1_tot_gspace)
      END DO
      IF (gapw) THEN
         CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs, rho1_tot_gspace)
         IF (ASSOCIATED(p_env%local_rho_set%rho0_mpole%rhoz_cneo_s_gs)) THEN
            CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rhoz_cneo_s_gs, rho1_tot_gspace)
         END IF
      END IF
      !
      BLOCK
         TYPE(pw_c1d_gs_type) :: v_hartree_gspace
         CALL auxbas_pw_pool%create_pw(v_hartree_gspace)
         CALL pw_poisson_solve(poisson_env, rho1_tot_gspace, &
                               energy_hartree, &
                               v_hartree_gspace)
         CALL pw_transfer(v_hartree_gspace, v_hartree_rspace)
         CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
      END BLOCK
      !
      CALL pw_scale(v_hartree_rspace, v_hartree_rspace%pw_grid%dvol)

      CALL auxbas_pw_pool%give_back_pw(rho1_tot_gspace)

      ! *** calculate the xc potential ***
      NULLIFY (v_xc, v_xc_tau)
      IF (gapw_xc) THEN
         CALL get_qs_env(qs_env, rho_xc=rho0_fxc)
         rho1_fxc => rho1_xc
      ELSE
         CALL get_qs_env(qs_env, rho=rho0_fxc)
         rho1_fxc => rho1
      END IF

      NULLIFY (rho_atom_set, rho1_atom_set)
      IF (do_onecenter) THEN
         CALL get_qs_env(qs_env, rho_atom_set=rho_atom_set)
         rho1_atom_set => p_env%local_rho_set%rho_atom_set
      END IF
      CALL qs_fxc_apply(qs_env, kpp1_env%deriv_set, kpp1_env%rho_set, &
                        rho1_fxc, rho_atom_set, xc_section, &
                        do_onecenter, v_xc, v_xc_tau, rho1_atom_set, &
                        compute_virial=calc_virial, virial_xc=virial)

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
      END DO
      v_rspace_new => v_xc
      IF (SIZE(v_xc) /= nspins) THEN
         CALL auxbas_pw_pool%give_back_pw(v_xc(2))
      END IF
      NULLIFY (v_xc)
      IF (ASSOCIATED(v_xc_tau)) THEN
         DO ispin = 1, nspins
            CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
         END DO
         IF (SIZE(v_xc_tau) /= nspins) THEN
            CALL auxbas_pw_pool%give_back_pw(v_xc_tau(2))
         END IF
      END IF

      alpha = 1.0_dp
      IF (nspins == 1) alpha = 2.0_dp

      !-------------------------------!
      ! Add both hartree and xc terms !
      !-------------------------------!
      DO ispin = 1, nspins
         CALL dbcsr_set(kpp1_env%v_ao(ispin)%matrix, 0.0_dp)

         IF (gapw_xc) THEN
            ! XC and Hartree are integrated separatedly
            ! XC uses the soft basis set only
            CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                    pmat=rho_ao(ispin), &
                                    hmat=kpp1_env%v_ao(ispin), &
                                    qs_env=qs_env, &
                                    calculate_forces=my_calc_forces, gapw=gapw_xc)
            IF (ASSOCIATED(v_xc_tau)) THEN
               CALL integrate_v_rspace(v_rspace=v_xc_tau(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       compute_tau=.TRUE., &
                                       calculate_forces=my_calc_forces, gapw=gapw_xc)
            END IF
            ! add Hartree for SINGLETS
            CALL pw_copy(v_hartree_rspace, v_rspace_new(1))
            CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                    pmat=rho_ao(ispin), &
                                    hmat=kpp1_env%v_ao(ispin), &
                                    qs_env=qs_env, &
                                    calculate_forces=my_calc_forces, gapw=gapw)
         ELSE
            CALL pw_axpy(v_hartree_rspace, v_rspace_new(ispin))
            IF (lrigpw) THEN
               IF (ASSOCIATED(v_xc_tau)) CPABORT("Meta-GGA functionals not supported with LRI!")

               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, .FALSE., "LRI_AUX")
               DO ikind = 1, nkind
                  CALL para_env%sum(lri_v_int(ikind)%v_int)
               END DO
               ALLOCATE (k1mat(1))
               k1mat(1)%matrix => kpp1_env%v_ao(ispin)%matrix
               IF (lri_env%exact_1c_terms) THEN
                  CALL integrate_v_rspace_diagonal(v_rspace_new(ispin), k1mat(1)%matrix, &
                                                   rho_ao(ispin)%matrix, qs_env, my_calc_forces, "ORB")
               END IF
               CALL calculate_lri_ks_matrix(lri_env, lri_v_int, k1mat, atomic_kind_set)
               DEALLOCATE (k1mat)
            ELSE
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw)
               IF (ASSOCIATED(v_xc_tau)) THEN
                  CALL integrate_v_rspace(v_rspace=v_xc_tau(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          compute_tau=.TRUE., &
                                          calculate_forces=my_calc_forces, gapw=gapw)
               END IF
            END IF
         END IF

         CALL dbcsr_add(p_env%kpp1(ispin)%matrix, kpp1_env%v_ao(ispin)%matrix, 1.0_dp, alpha)
      END DO

      IF (gapw) THEN
         CALL Vh_1c_gg_integrals(qs_env, energy_hartree_1c, &
                                 p_env%hartree_local%ecoul_1c, &
                                 p_env%local_rho_set, &
                                 para_env, tddft=.TRUE., core_2nd=.TRUE.)
         CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, &
                                    calculate_forces=my_calc_forces, &
                                    local_rho_set=p_env%local_rho_set)
         !  ***  Add single atom contributions to the KS matrix ***
         ! remap pointer
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)
         CALL update_ks_atom(qs_env, ksmat, psmat, forces=my_calc_forces, tddft=.TRUE., &
                             rho_atom_external=p_env%local_rho_set%rho_atom_set)
      ELSE IF (gapw_xc) THEN
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)
         CALL update_ks_atom(qs_env, ksmat, psmat, forces=my_calc_forces, tddft=.TRUE., &
                             rho_atom_external=p_env%local_rho_set%rho_atom_set)
      END IF

      CALL auxbas_pw_pool%give_back_pw(v_hartree_rspace)
      DO ispin = 1, SIZE(v_rspace_new)
         CALL auxbas_pw_pool%give_back_pw(v_rspace_new(ispin))
      END DO
      DEALLOCATE (v_rspace_new)
      IF (ASSOCIATED(v_xc_tau)) THEN
         DO ispin = 1, SIZE(v_xc_tau)
            CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
         END DO
         DEALLOCATE (v_xc_tau)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_kpp1

! **************************************************************************************************
!> \brief checks that the intenal storage is allocated, and allocs it if needed
!> \param kpp1_env the environment to check
!> \param qs_env the qs environment this kpp1_env lives in
!> \param xc_section ...
!> \author Fawzi Mohamed
!> \note
!>      private routine
! **************************************************************************************************
   SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, xc_section)

      TYPE(qs_kpp1_env_type)                             :: kpp1_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section

      INTEGER                                            :: ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: admm_xc_section

! ------------------------------------------------------------------

      CALL get_qs_env(qs_env, pw_env=pw_env, matrix_s=matrix_s, &
                      admm_env=admm_env, dft_control=dft_control)

      nspins = dft_control%nspins

      IF (.NOT. ASSOCIATED(kpp1_env%v_ao)) THEN
         CALL dbcsr_allocate_matrix_set(kpp1_env%v_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (kpp1_env%v_ao(ispin)%matrix)
            CALL dbcsr_copy(kpp1_env%v_ao(ispin)%matrix, matrix_s(1)%matrix, &
                            name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%deriv_set)) THEN
         IF (dft_control%qs_control%gapw_xc) THEN
            CALL get_qs_env(qs_env, rho_xc=rho)
         ELSE
            CALL get_qs_env(qs_env, rho=rho)
         END IF
         ALLOCATE (kpp1_env%deriv_set, kpp1_env%rho_set)
         CALL qs_fxc_prep(qs_env, rho, &
                          kpp1_env%rho_set, kpp1_env%deriv_set, &
                          xc_section, pw_env, is_triplet=.FALSE.)
      END IF

      ! ADMM Correction
      IF (dft_control%do_admm) THEN
         IF (admm_env%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            IF (.NOT. ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
               ALLOCATE (kpp1_env%deriv_set_admm, kpp1_env%rho_set_admm)
               admm_xc_section => admm_env%xc_section_aux
               CALL get_admm_env(admm_env, rho_aux_fit=rho)
               CALL qs_fxc_prep(qs_env, rho, kpp1_env%rho_set_admm, kpp1_env%deriv_set_admm, &
                                admm_xc_section, pw_env, is_triplet=.FALSE.)
            END IF
         END IF
      END IF

   END SUBROUTINE kpp1_check_i_alloc
! **************************************************************************************************
!> \brief ...
!> \param rho1 ...
!> \param rho1_tot_gspace ...
!> \param out_unit ...
! **************************************************************************************************
   SUBROUTINE print_densities(rho1, rho1_tot_gspace, out_unit)

      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: rho1_tot_gspace
      INTEGER                                            :: out_unit

      REAL(KIND=dp)                                      :: total_rho_gspace
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho1_r

      NULLIFY (tot_rho1_r)

      total_rho_gspace = pw_integrate_function(rho1_tot_gspace, isign=-1)
      IF (out_unit > 0) THEN
         CALL qs_rho_get(rho1, tot_rho_r=tot_rho1_r)
         WRITE (UNIT=out_unit, FMT="(T3,A,T60,F20.10)") &
            "KPP1 total charge density (r-space):", &
            accurate_sum(tot_rho1_r), &
            "KPP1 total charge density (g-space):", &
            total_rho_gspace
      END IF

   END SUBROUTINE print_densities

END MODULE qs_kpp1_env_methods
