!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the CPKS equation and the resulting forces
!> \par History
!>       03.2014 created
!>       09.2019 Moved from KG to Kohn-Sham
!>       11.2019 Moved from energy_correction
!> \author JGH
! **************************************************************************************************
MODULE response_solver
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_constants,                 ONLY: ot_precond_full_single_inverse
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type,&
                                              total_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_methods,               ONLY: build_dm_response,&
                                              linres_solver
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
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

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'response_solver'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

   PUBLIC ::  response_equation, response_force, ks_ref_potential

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param cpmos RHS of equation as Ax + b = 0 (sign of b)
!> \param lr_section ...
! **************************************************************************************************
   SUBROUTINE response_equation(qs_env, p_env, cpmos, lr_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos
      TYPE(section_vals_type), OPTIONAL, POINTER         :: lr_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_equation', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, nspins
      LOGICAL                                            :: should_stop
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: psi0, psi1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      ENDIF

      ! initialized linres_control
      NULLIFY (linres_control)
      CALL linres_control_create(linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      IF (PRESENT(lr_section)) THEN
         CALL section_vals_val_get(lr_section, "RESTART", l_val=linres_control%linres_restart)
         CALL section_vals_val_get(lr_section, "MAX_ITER", i_val=linres_control%max_iter)
         CALL section_vals_val_get(lr_section, "EPS", r_val=linres_control%eps)
         CALL section_vals_val_get(lr_section, "EPS_FILTER", r_val=linres_control%eps_filter)
         CALL section_vals_val_get(lr_section, "RESTART_EVERY", i_val=linres_control%restart_every)
         CALL section_vals_val_get(lr_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
         CALL section_vals_val_get(lr_section, "ENERGY_GAP", r_val=linres_control%energy_gap)
      ELSE
         linres_control%linres_restart = .FALSE.
         linres_control%max_iter = 100
         linres_control%eps = 1.0e-10_dp
         linres_control%eps_filter = 1.0e-15_dp
         linres_control%restart_every = 50
         linres_control%preconditioner_type = ot_precond_full_single_inverse
         linres_control%energy_gap = 0.02_dp
      END IF

      ! initialized p_env
      CALL p_env_release(p_env)
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE., &
                        linres_control=linres_control)
      CALL set_qs_env(qs_env, linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)
      p_env%os_valid = .FALSE.
      p_env%new_preconditioner = .TRUE.
      !
      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos)
      nspins = dft_control%nspins

      ALLOCATE (psi0(nspins), psi1(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         psi0(ispin)%matrix => mo_coeff
         CALL cp_fm_get_info(mo_coeff, matrix_struct=fm_struct)
         NULLIFY (psi1(ispin)%matrix)
         CALL cp_fm_create(psi1(ispin)%matrix, fm_struct)
         CALL cp_fm_set_all(psi1(ispin)%matrix, 0.0_dp)
      ENDDO

      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      CALL dbcsr_allocate_matrix_set(p_env%p1, nspins)
      CALL dbcsr_allocate_matrix_set(p_env%w1, nspins)
      DO ispin = 1, nspins
         ALLOCATE (p_env%p1(ispin)%matrix, p_env%w1(ispin)%matrix)
         CALL dbcsr_create(matrix=p_env%p1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=p_env%w1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(p_env%p1(ispin)%matrix, sab_orb)
         CALL cp_dbcsr_alloc_block_from_nbl(p_env%w1(ispin)%matrix, sab_orb)
      END DO

      CALL linres_solver(p_env, qs_env, psi1, cpmos, psi0, iounit, should_stop)

      DO ispin = 1, nspins
         CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
      END DO
      CALL build_dm_response(psi0, psi1, p_env%p1)
      DO ispin = 1, nspins
         CALL dbcsr_scale(p_env%p1(ispin)%matrix, 0.5_dp)
      END DO
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      DO ispin = 1, nspins
         CALL calculate_wz_matrix(mos(ispin)%mo_set, psi1(ispin)%matrix, matrix_ks(ispin)%matrix, &
                                  p_env%w1(ispin)%matrix)
      ENDDO
      DO ispin = 1, nspins
         CALL cp_fm_release(psi1(ispin)%matrix)
      ENDDO
      DEALLOCATE (psi0, psi1)
      CALL linres_control_release(linres_control)

      CALL timestop(handle)

   END SUBROUTINE response_equation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param matrix_hz ...
! **************************************************************************************************
   SUBROUTINE response_force(qs_env, p_env, vh_rspace, vxc_rspace, vtau_rspace, matrix_hz)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_force', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, natom, nder, &
                                                            nimages, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: lsd, use_virial
      REAL(KIND=dp)                                      :: eps_ppnl, focc, total_rhoz, zehartree
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2, ftot3
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pz, matrix_wz, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoz_tot_gspace, vhxc_rspace, &
                                                            zv_hartree_gspace, zv_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, rhoz_g, rhoz_r, tau_r, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_fun_section, xc_section
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      CPASSERT(ASSOCIATED(qs_env))

      NULLIFY (ks_env, sab_orb, sac_ppl, sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, &
                      sab_orb=sab_orb, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, force=force)
      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      matrix_pz => p_env%p1
      matrix_wz => p_env%w1

      nspins = SIZE(matrix_pz, 1)
      nder = 1
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_pz(1)%matrix, matrix_pz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      ! Kinetic energy matrix
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_pz(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1)-fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dT      ", fodeb
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_pz(1)%matrix, matrix_pz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, dft_control=dft_control, force=force)
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
      nimages = 1
      ALLOCATE (matrix_p(nspins, 1), matrix_h(nspins, 1))
      DO ispin = 1, nspins
         matrix_p(ispin, 1)%matrix => matrix_pz(ispin)%matrix
         matrix_h(ispin, 1)%matrix => scrm(ispin)%matrix
      END DO
      IF (ASSOCIATED(sac_ppl)) THEN
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         CALL build_core_ppl(matrix_h, matrix_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1)-fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppl   ", fodeb
         END IF
      END IF
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         CALL build_core_ppnl(matrix_h, matrix_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)-fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppnl  ", fodeb
         END IF
      END IF
      DEALLOCATE (matrix_p)
      IF (debug_forces) THEN
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1)-ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pz*dHcore", fodeb
      END IF

      ! Vhxc
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL pw_pool_create_pw(auxbas_pw_pool, vhxc_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_transfer(vh_rspace%pw, vhxc_rspace%pw)
         CALL pw_axpy(vxc_rspace(ispin)%pw, vhxc_rspace%pw)
         CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                 hmat=matrix_h(1, 1), pmat=matrix_pz(ispin), &
                                 qs_env=qs_env, calculate_forces=.TRUE.)
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                    hmat=scrm(1), pmat=matrix_pz(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
         END IF
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1)-fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVhxc   ", fodeb
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vhxc_rspace%pw)

      ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoz_r(ispin)%pw, rhoz_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_pz(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                 total_rho=total_rhoz)
         CALL pw_axpy(rhoz_g(ispin)%pw, rhoz_tot_gspace%pw)
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace%pw, zehartree, &
                            zv_hartree_gspace%pw)
      CALL pw_transfer(zv_hartree_gspace%pw, zv_hartree_rspace%pw)
      CALL pw_scale(zv_hartree_rspace%pw, zv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1)-fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(rhoz)*dncore ", fodeb
      END IF
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = rhoz_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rhoz_r, rhoz_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)
      !
      !
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(zv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=matrix_hz(ispin), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1)-fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz ", fodeb
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoz_r, rhoz_g, v_xc)
      IF (debug_forces) THEN
         ALLOCATE (ftot3(3, natom))
         CALL total_qs_force(ftot3, force, atomic_kind_set)
         fodeb(1:3) = ftot3(1:3, 1)-ftot2(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pin*V(rhoz)", fodeb
      END IF

      ! Overlap matrix
      ! H(drho+dz) + Wz
      focc = 1.0_dp
      IF (nspins == 1) focc = 2.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL calculate_whz_matrix(mos(ispin)%mo_set, matrix_hz(ispin)%matrix, &
                                   matrix_wz(ispin)%matrix, focc)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1)-fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1)-ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         fodeb(1:3) = ftot2(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Total Force ", fodeb
         DEALLOCATE (ftot1, ftot2, ftot3)
      END IF

      CALL dbcsr_deallocate_matrix_set(scrm)
      DEALLOCATE (matrix_h)

      CALL timestop(handle)

   END SUBROUTINE response_force

! **************************************************************************************************
!> \brief Calculate the response W matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param mo_set type containing the full matrix of the MO and the eigenvalues
!> \param psi1 response orbitals
!> \param ks_matrix Kohn-Sham sparse matrix
!> \param w_matrix sparse matrix
!> \par History
!>               adapted from calculate_w_matrix_1
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_wz_matrix(mo_set, psi1, ks_matrix, w_matrix)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_fm_type), POINTER                          :: psi1
      TYPE(dbcsr_type), POINTER                          :: ks_matrix, w_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_wz_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: ksmat, scrv

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=mo_set%mo_coeff, ncol_global=ncol, nrow_global=nrow)
      CALL cp_fm_create(scrv, mo_set%mo_coeff%matrix_struct, "scr vectors")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=mo_set%mo_coeff%matrix_struct%para_env, &
                               context=mo_set%mo_coeff%matrix_struct%context)
      CALL cp_fm_create(ksmat, fm_struct_tmp, name="KS")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_dbcsr_sm_fm_multiply(ks_matrix, mo_set%mo_coeff, scrv, ncol)
      CALL cp_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, mo_set%mo_coeff, scrv, 0.0_dp, ksmat)
      CALL cp_gemm("N", "N", nrow, ncol, ncol, 1.0_dp, mo_set%mo_coeff, ksmat, 0.0_dp, scrv)
      CALL dbcsr_set(w_matrix, 0.0_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=scrv, matrix_g=psi1, &
                                 ncol=mo_set%homo, alpha=0.5_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=psi1, matrix_g=scrv, &
                                 ncol=mo_set%homo, alpha=0.5_dp)
      CALL cp_fm_release(scrv)
      CALL cp_fm_release(ksmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_wz_matrix

! **************************************************************************************************
!> \brief Calculate the Wz matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param mo_set type containing the full matrix of the MO and the eigenvalues
!> \param hzm ...
!> \param w_matrix sparse matrix
!> \param focc ...
!> \par History
!>               adapted from calculate_w_matrix_1
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_whz_matrix(mo_set, hzm, w_matrix, focc)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(dbcsr_type), POINTER                          :: hzm, w_matrix
      REAL(KIND=dp), INTENT(IN)                          :: focc

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_whz_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao, norb
      REAL(KIND=dp)                                      :: falpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_mat
      TYPE(cp_fm_type), POINTER                          :: c0vec, chcmat, hcvec

      CALL timeset(routineN, handle)

      falpha = focc

      c0vec => mo_set%mo_coeff
      CALL cp_fm_create(hcvec, c0vec%matrix_struct, "hcvec")
      CALL cp_fm_get_info(hcvec, matrix_struct=fm_struct, nrow_global=nao, ncol_global=norb)
      CALL cp_fm_struct_create(fm_struct_mat, context=fm_struct%context, nrow_global=norb, &
                               ncol_global=norb, para_env=fm_struct%para_env)
      CALL cp_fm_create(chcmat, fm_struct_mat)
      CALL cp_fm_struct_release(fm_struct_mat)

      CALL cp_dbcsr_sm_fm_multiply(hzm, c0vec, hcvec, norb)
      CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, c0vec, hcvec, 0.0_dp, chcmat)
      CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, c0vec, chcmat, 0.0_dp, hcvec)

      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=hcvec, matrix_g=c0vec, ncol=norb, alpha=falpha)

      CALL cp_fm_release(hcvec)
      CALL cp_fm_release(chcmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_whz_matrix

! **************************************************************************************************
!> \brief calculate the Kohn-Sham reference potential
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param ehartree ...
!> \param exc ...
!> \par History
!>      10.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, ehartree, exc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: ehartree, exc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ks_ref_potential', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iab, ispin, nspins
      REAL(dp)                                           :: eovrl, eself
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env)

      nspins = dft_control%nspins

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      !
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, .FALSE., E_overlap_core=eovrl)
      ehartree = ehartree+eovrl+eself

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)

      ! allocate potentials
      IF (ASSOCIATED(vh_rspace)) THEN
         CALL pw_release(vh_rspace%pw)
      ELSE
         ALLOCATE (vh_rspace)
         NULLIFY (vh_rspace%pw)
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL pw_release(vxc_rspace(iab)%pw)
         END DO
      ELSE
         ALLOCATE (vxc_rspace(nspins))
         DO iab = 1, nspins
            NULLIFY (vxc_rspace(iab)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO iab = 1, SIZE(vtau_rspace)
               CALL pw_release(vtau_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vtau_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vtau_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vtau_rspace)
      END IF

      pw_grid => v_hartree_rspace%pw%pw_grid
      CALL pw_create(vh_rspace%pw, pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      DO ispin = 1, nspins
         NULLIFY (vxc_rspace(ispin)%pw)
         CALL pw_create(vxc_rspace(ispin)%pw, pw_grid, &
                        use_data=REALDATA3D, in_space=REALSPACE)
         IF (ASSOCIATED(vtau_rspace)) THEN
            NULLIFY (vtau_rspace(ispin)%pw)
            CALL pw_create(vtau_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
      END DO
      !
      CALL pw_transfer(v_hartree_rspace%pw, vh_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_rspace(ispin)%pw, vxc_rspace(ispin)%pw)
            CALL pw_scale(vxc_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_transfer(v_tau_rspace(ispin)%pw, vtau_rspace(ispin)%pw)
               CALL pw_scale(vtau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(vxc_rspace(ispin)%pw)
         END DO
      END IF

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
            END IF
         ENDDO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential

! **************************************************************************************************

END MODULE response_solver
