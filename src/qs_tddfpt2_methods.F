!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_methods
   USE admm_types,                      ONLY: admm_type
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_solve
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_fm,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_contracted_trace,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              cp_fm_pool_type,&
                                              fm_pool_create,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_maxabsval, &
        cp_fm_p_type, cp_fm_read_unformatted, cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, &
        cp_fm_to_fm, cp_fm_type, cp_fm_write_unformatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_distribution_type, dbcsr_get_info, &
        dbcsr_init_p, dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE input_constants,                 ONLY: cholesky_dbcsr,&
                                              cholesky_inverse,&
                                              cholesky_off,&
                                              cholesky_restore,&
                                              tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_dot_product,&
                                              accurate_sum
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_min,&
                                              mp_sum,&
                                              mp_wait
   USE moments_utils,                   ONLY: get_reference_point
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
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
                                              pw_p_type,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: rRc_xyz_ao
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_post_gpw,                 ONLY: make_lumo
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_sub_env_init,&
                                              tddfpt_sub_env_release,&
                                              tddfpt_subgroup_env_type
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt

! **************************************************************************************************
!> \brief Ground state molecular orbitals.
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_ground_state_mos
      !> occupied MOs stored in a matrix form [nao x nmo_occ]
      TYPE(cp_fm_type), POINTER :: mos_occ
      !> virtual MOs stored in a matrix form [nao x nmo_virt]
      TYPE(cp_fm_type), POINTER :: mos_virt
      !> occupied orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_occ
      !> virtual orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_virt
      !> phase of occupied MOs; +1.0 -- positive, -1.0 -- negative;
      !> it is mainly needed to make the restart file transferable
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: phases_occ
   END TYPE tddfpt_ground_state_mos

! **************************************************************************************************
!> \brief Set of temporary ("work") matrices.
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_work_matrices
      !
      ! *** globally distributed dense matrices ***
      !
      !> pool of dense [nao x nmo_occ(spin)] matrices;
      !> used mainly to dynamically expand the list of trial vectors
      TYPE(cp_fm_pool_p_type), ALLOCATABLE, DIMENSION(:) :: fm_pool_ao_mo_occ
      !> S * mos_occ(spin)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0
      !> S * \rho_0(spin)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0_C0T
      !> globally distributed dense matrices with shape [nao x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_ao_mo_occ
      !> globally distributed dense matrices with shape [nmo_occ(spin) x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_mo_occ_mo_occ
      !> globally distributed dense matrices with shape [nmo_virt(spin) x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_mo_virt_mo_occ
      !
      ! *** dense matrices distributed across parallel (sub)groups ***
      !
      !> evects_sub(1:nspins, 1:nstates): a copy of the last 'nstates' trial vectors distributed
      !> across parallel (sub)groups. Here 'nstates' is the number of requested excited states which
      !> is typically much smaller than the total number of Krylov's vectors. Allocated only if
      !> the number of parallel groups > 1, otherwise we use the original globally distributed vectors.
      !> evects_sub(spin, state) == null() means that the trial vector is assigned to a different (sub)group
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)    :: evects_sub
      !> action of TDDFPT operator on trial vectors distributed across parallel (sub)groups
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)    :: Aop_evects_sub
      !> electron density expressed in terms of atomic orbitals using primary basis set
      TYPE(cp_fm_type), POINTER                        :: rho_ao_orb_fm_sub
      !
      ! NOTE: we do not need the next 2 matrices in case of a sparse matrix 'tddfpt_subgroup_env_type%admm_A'
      !
      !> electron density expressed in terms of atomic orbitals using auxiliary basis set;
      !> can be seen as a group-specific version of the matrix 'admm_type%work_aux_aux'
      TYPE(cp_fm_type), POINTER                        :: rho_ao_aux_fit_fm_sub
      !> group-specific version of the matrix 'admm_type%work_aux_orb' with shape [nao_aux x nao]
      TYPE(cp_fm_type), POINTER                        :: wfm_aux_orb_sub
      !
      ! *** sparse matrices distributed across parallel (sub)groups ***
      !
      !> sparse matrix with shape [nao x nao] distributed across subgroups;
      !> Aop_evects_sub(spin,:) = A_ia_munu_sub(spin) * mos_occ(spin)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: A_ia_munu_sub
      !
      ! *** structures to store electron densities distributed across parallel (sub)groups ***
      !
      !> electron density in terms of primary basis set
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct_sub
      !> electron density in terms of auxiliary basis set
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit_struct_sub
      !> group-specific copy of a Coulomb/xc-potential on a real-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace_sub
      !> group-specific copy of a reciprocal-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: wpw_gspace_sub
      !> group-specific copy of a real-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: wpw_rspace_sub
      !
      ! *** globally distributed matrices required to compute exact exchange terms ***
      !
      !> globally distributed version of the matrix 'rho_ao_orb_fm_sub' to store the electron density
      TYPE(cp_fm_type), POINTER                          :: hfx_fm_ao_ao
      !> sparse matrix to store the electron density in terms of auxiliary (ADMM calculation)
      !> or primary (regular calculation) basis set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hfx_rho_ao
      !> exact exchange expressed in terms of auxiliary or primary basis set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hfx_hmat
   END TYPE tddfpt_work_matrices

! **************************************************************************************************
!> \brief Collection of variables required to evaluate adiabatic TDDFPT kernel.
!> \par History
!>   * 12.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_kernel_env_type
      ! ground state electron density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      ! response density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set
      !> first and second derivatives of exchange-correlation functional
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      !> XC input section
      TYPE(section_vals_type), POINTER                   :: xc_section
      !> flags which indicate required components of the exchange-correlation functional
      !> (density, gradient, etc)
      TYPE(xc_rho_cflags_type)                           :: xc_rho1_cflags
      !> the method used to compute derivatives
      INTEGER                                            :: deriv_method_id
      !> the density smoothing method
      INTEGER                                            :: rho_smooth_id
      !> scaling coefficients in the linear combination:
      !> K = alpha * K_{\alpha,\alpha} + beta * K_{\alpha,\beta}
      REAL(kind=dp)                                      :: alpha, beta
   END TYPE tddfpt_kernel_env_type

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT calculation.
!> \param qs_env  Quickstep environment
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!>    * 03.2017 cleaned and refactored [Sergey Chulkov]
!> \note Based on the subroutines tddfpt_env_init(), and tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER :: energy_unit, handle, ispin, istate, iter, log_unit, mult, niters, nmo_avail, &
         nmo_occ, nmo_virt, nspins, nstates, nstates_read
      LOGICAL                                            :: do_admm, do_hfx, explicit_xc, &
                                                            is_restarted
      REAL(kind=dp)                                      :: C_hf, conv
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      REAL(kind=dp), DIMENSION(:), POINTER               :: evals_virt_spin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: evals_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dipole_op_mos_occ, evects, S_evects
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_virt
      TYPE(cp_fm_type), POINTER                          :: mos_virt_spin
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, &
                                                            tddfpt_print_section, tddfpt_section, &
                                                            xc_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), ALLOCATABLE, &
         DIMENSION(:)                                    :: gs_mos
      TYPE(tddfpt_kernel_env_type)                       :: kernel_env, kernel_env_admm_aux
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env, cell, dft_control, input, matrix_ks, matrix_s, mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control

      nspins = dft_control%nspins
      CPASSERT(dft_control%nimages <= 1)

      IF (tddfpt_control%nstates <= 0) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL cp_warn(__LOCATION__, "TDDFPT calculation was requested for "// &
                      TRIM(nstates_str)//" excited states: nothing to do.")
         CALL timestop(handle)
         RETURN
      END IF

      tddfpt_section => section_vals_get_subs_vals(input, "PROPERTIES%TDDFPT")
      tddfpt_print_section => section_vals_get_subs_vals(tddfpt_section, "PRINT")

      xc_section => section_vals_get_subs_vals(tddfpt_section, "XC%XC_FUNCTIONAL")
      CALL section_vals_get(xc_section, explicit=explicit_xc)
      IF (explicit_xc) THEN
         xc_section => section_vals_get_subs_vals(tddfpt_section, "XC")
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")

      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
         do_hfx = (C_hf /= 0.0_dp)
      END IF

      do_admm = do_hfx .AND. dft_control%do_admm
      IF (do_admm) THEN
         IF (explicit_xc) THEN
            ! 'admm_env%xc_section_primary' and 'admm_env%xc_section_aux' need to be redefined
            CALL cp_abort(__LOCATION__, &
                          "ADMM is not implemented for a TDDFT kernel XC-functional which is different from "// &
                          "the one used for the ground-state calculation. A ground-state 'admm_env' cannot be reused.")
         END IF

         CALL get_qs_env(qs_env, admm_env=admm_env)
      END IF

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarised calculations")
      END IF

      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      ALLOCATE (gs_mos(nspins))

      ! when the number of unoccupied orbitals is limited and OT has been used for the ground-state DFT calculation,
      ! compute the missing unoccupied orbitals using OT as well.
      NULLIFY (evals_virt, evals_virt_spin, mos_virt, mos_virt_spin)
      IF (ASSOCIATED(scf_env)) THEN
         IF (scf_env%method == ot_method_nr .AND. tddfpt_control%nlumo > 0) THEN
            ! As OT with ADDED_MOS/=0 is currently not implemented, the following block is equivalent to:
            ! nmo_virt = tddfpt_control%nlumo

            ! number of already computed unoccupied orbitals (added_mos) .
            nmo_virt = HUGE(nmo_virt)
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo_avail, homo=nmo_occ)
               nmo_virt = MIN(nmo_virt, nmo_avail - nmo_occ)
            END DO
            ! number of unoccupied orbitals to compute
            nmo_virt = tddfpt_control%nlumo - nmo_virt

            IF (nmo_virt > 0) THEN
               ALLOCATE (evals_virt(nspins), mos_virt(nspins))
               ! the number of actually computed unoccupied orbitals will be stored as 'nmo_avail'
               CALL make_lumo(qs_env, scf_env, mos_virt, evals_virt, nmo_virt, nmo_avail)
            END IF
         END IF
      END IF

      DO ispin = 1, nspins
         IF (ASSOCIATED(evals_virt)) &
            evals_virt_spin => evals_virt(ispin)%array
         IF (ASSOCIATED(mos_virt)) &
            mos_virt_spin => mos_virt(ispin)%matrix

         CALL tddfpt_init_ground_state_mos(gs_mos=gs_mos(ispin), mo_set=mos(ispin)%mo_set, nlumo=tddfpt_control%nlumo, &
                                           blacs_env=blacs_env, cholesky_method=cholesky_restore, &
                                           matrix_ks=matrix_ks(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                                           mos_virt=mos_virt_spin, evals_virt=evals_virt_spin)
      END DO

      IF (ASSOCIATED(evals_virt)) THEN
         DO ispin = 1, SIZE(evals_virt)
            IF (ASSOCIATED(evals_virt(ispin)%array)) DEALLOCATE (evals_virt(ispin)%array)
         END DO
         DEALLOCATE (evals_virt)
      END IF

      IF (ASSOCIATED(mos_virt)) THEN
         DO ispin = 1, SIZE(mos_virt)
            IF (ASSOCIATED(mos_virt(ispin)%matrix)) CALL cp_fm_release(mos_virt(ispin)%matrix)
         END DO
         DEALLOCATE (mos_virt)
      END IF

      ! components of the dipole operator
      CALL tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control, gs_mos, qs_env)

      ! multiplicity of molecular system
      IF (nspins > 1) THEN
         mult = ABS(SIZE(gs_mos(1)%evals_occ) - SIZE(gs_mos(2)%evals_occ)) + 1
         IF (mult > 2) &
            CALL cp_warn(__LOCATION__, "There is a convergence issue for multiplicity >= 3")
      ELSE
         IF (tddfpt_control%rks_triplets) THEN
            mult = 3
         ELSE
            mult = 1
         END IF
      END IF

      ! split mpi communicator
      ALLOCATE (evects(nspins, 1))
      DO ispin = 1, nspins
         evects(ispin, 1)%matrix => gs_mos(ispin)%mos_occ
      END DO
      CALL tddfpt_sub_env_init(sub_env, qs_env, mos_occ=evects(:, 1))
      DEALLOCATE (evects)

      ! allocate pools and work matrices
      nstates = tddfpt_control%nstates
      CALL tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, qs_env, sub_env)

      ! create kernel environment
      CALL tddfpt_construct_ground_state_orb_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                     is_rks_triplets=tddfpt_control%rks_triplets, &
                                                     qs_env=qs_env, sub_env=sub_env, &
                                                     wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub)

      IF (do_admm) THEN
         CALL tddfpt_create_kernel_env(kernel_env=kernel_env, &
                                       rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                       xc_section=admm_env%xc_section_primary, &
                                       is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)

         CALL tddfpt_construct_aux_fit_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                               rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
                                               qs_env=qs_env, sub_env=sub_env, &
                                               wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
                                               wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
                                               wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

         CALL tddfpt_create_kernel_env(kernel_env=kernel_env_admm_aux, &
                                       rho_struct_sub=work_matrices%rho_aux_fit_struct_sub, &
                                       xc_section=admm_env%xc_section_aux, &
                                       is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)
      ELSE
         CALL tddfpt_create_kernel_env(kernel_env=kernel_env, &
                                       rho_struct_sub=work_matrices%rho_orb_struct_sub, &
                                       xc_section=xc_section, &
                                       is_rks_triplets=tddfpt_control%rks_triplets, sub_env=sub_env)
      END IF

      ALLOCATE (evals(nstates))
      ALLOCATE (evects(nspins, nstates), S_evects(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (evects(ispin, istate)%matrix, S_evects(ispin, istate)%matrix)
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_evects(ispin, istate)%matrix)
         END DO
      END DO

      ! reuse Ritz vectors from the previous calculation if available
      IF (tddfpt_control%is_restart) THEN
         nstates_read = tddfpt_read_restart(evects=evects, evals=evals, gs_mos=gs_mos, &
                                            logger=logger, tddfpt_section=tddfpt_section, &
                                            tddfpt_print_section=tddfpt_print_section, &
                                            fm_pool_ao_mo_occ=work_matrices%fm_pool_ao_mo_occ, &
                                            blacs_env_global=blacs_env)
      ELSE
         nstates_read = 0
      END IF

      is_restarted = nstates_read >= nstates

      ! build the list of missed singly excited states and sort them in ascending order according to their excitation energies
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_guess_vectors(evects=evects, evals=evals, gs_mos=gs_mos, log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "GUESS_VECTORS")

      CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)
      CALL tddfpt_orthonormalize_psi1_psi1(evects, nstates, S_evects, matrix_s(1)%matrix)

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "DETAILED_ENERGY", extension=".tddfptLog")

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') 'TDDFPT WAVEFUNCTION OPTIMIZATION'
            WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
            WRITE (log_unit, '(1X,50("-"))')
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         DO
            ! *** perform Davidson iterations ***
            conv = tddfpt_davidson_solver(evects=evects, evals=evals, S_evects=S_evects, gs_mos=gs_mos, &
                                          do_hfx=do_hfx, tddfpt_control=tddfpt_control, qs_env=qs_env, &
                                          kernel_env=kernel_env, kernel_env_admm_aux=kernel_env_admm_aux, &
                                          sub_env=sub_env, logger=logger, &
                                          iter_unit=log_unit, energy_unit=energy_unit, &
                                          tddfpt_print_section=tddfpt_print_section, work_matrices=work_matrices)

            ! at this point at least one of the following conditions are met:
            ! a) convergence criteria has been achieved;
            ! b) maximum number of iterations has been reached;
            ! c) Davidson iterations must be restarted due to lack of Krylov vectors or numerical instability

            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)
            ! terminate the loop if either (a) or (b) is true ...
            IF ((conv <= tddfpt_control%conv .AND. is_restarted) .OR. iter >= niters) EXIT

            ! ... otherwise restart Davidson iterations
            is_restarted = .TRUE.
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,10("-"),1X,A,1X,11("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF
         END DO

         ! write TDDFPT restart file at the last iteration if requested to do so
         CALL cp_iterate(logger%iter_info, increment=0, last=.TRUE.)
         CALL tddfpt_write_restart(evects=evects, evals=evals, gs_mos=gs_mos, &
                                   logger=logger, tddfpt_print_section=tddfpt_print_section)

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! print convergence summary
         IF (log_unit > 0) THEN
            CALL integer_to_string(iter, nstates_str)
            IF (conv <= tddfpt_control%conv) THEN
               WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run converged in "//TRIM(nstates_str)//" iteration(s) ***"
            ELSE
               WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ***"
            END IF
         END IF

         CALL cp_print_key_finished_output(energy_unit, logger, tddfpt_print_section, "DETAILED_ENERGY")
         CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, "Skipping TDDFPT wavefunction optimization")
      END IF

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit()

      CALL tddfpt_print_summary(log_unit, evects, evals, mult, dipole_op_mos_occ)
      CALL tddfpt_print_excitation_analysis(log_unit, evects, gs_mos, matrix_s(1)%matrix, &
                                            min_amplitude=tddfpt_control%min_excitation_amplitude)

      ! -- clean up all useless stuff
      DO istate = SIZE(evects, 2), 1, -1
         DO ispin = nspins, 1, -1
            CALL cp_fm_release(evects(ispin, istate)%matrix)
            CALL cp_fm_release(S_evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (evects, S_evects, evals)

      IF (do_admm) THEN
         CALL tddfpt_release_kernel_env(kernel_env_admm_aux)
      END IF
      CALL tddfpt_release_kernel_env(kernel_env)
      CALL tddfpt_release_work_matrices(work_matrices, sub_env)
      CALL tddfpt_sub_env_release(sub_env)

      IF (ALLOCATED(dipole_op_mos_occ)) THEN
         DO ispin = nspins, 1, -1
            DO istate = SIZE(dipole_op_mos_occ, 1), 1, -1
               CALL cp_fm_release(dipole_op_mos_occ(istate, ispin)%matrix)
            END DO
         END DO
         DEALLOCATE (dipole_op_mos_occ)
      END IF

      DO ispin = nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos(ispin))
      END DO
      DEALLOCATE (gs_mos)

      CALL timestop(handle)
   END SUBROUTINE tddfpt

! **************************************************************************************************
!> \brief Allocate work matrices.
!> \param work_matrices  work matrices (allocated on exit)
!> \param gs_mos         occupied and virtual molecular orbitals optimised for the ground state
!> \param nstates        number of excited states to converge
!> \param do_hfx         flag that requested to allocate work matrices required for computation
!>                       of exact-exchange terms
!> \param qs_env         Quickstep environment
!> \param sub_env        parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, qs_env, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(out)            :: work_matrices
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: nstates
      LOGICAL, INTENT(in)                                :: do_hfx
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_create_work_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igroup, ispin, istate, nao, &
                                                            nao_aux, ngroups, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: do_admm
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit, rho_ia_ao
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template_hfx
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_hfx
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      do_admm = do_hfx .AND. ASSOCIATED(sub_env%admm_A)
      IF (do_admm) THEN
         CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s_aux_fit)
         CALL dbcsr_get_info(matrix_s_aux_fit(1)%matrix, nfullrows_total=nao_aux)
      END IF

      NULLIFY (fm_struct)
      ALLOCATE (work_matrices%fm_pool_ao_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL fm_pool_create(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ALLOCATE (work_matrices%S_C0_C0T(nspins))
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0_C0T(ispin)%matrix)
         CALL cp_fm_create(work_matrices%S_C0_C0T(ispin)%matrix, fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (work_matrices%S_C0(nspins), work_matrices%wfm_ao_mo_occ(nspins))
      ALLOCATE (work_matrices%wfm_mo_occ_mo_occ(nspins), work_matrices%wfm_mo_virt_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0(ispin)%matrix, work_matrices%wfm_ao_mo_occ(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%S_C0(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%wfm_ao_mo_occ(ispin)%matrix)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%mos_occ, work_matrices%S_C0(ispin)%matrix, &
                                      ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, work_matrices%S_C0(ispin)%matrix, &
                      gs_mos(ispin)%mos_occ, 0.0_dp, work_matrices%S_C0_C0T(ispin)%matrix)

         NULLIFY (work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         NULLIFY (work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      IF (sub_env%is_split) THEN
         DO ispin = 1, nspins
            NULLIFY (fm_struct_evects(ispin)%struct)
            CALL cp_fm_struct_create(fm_struct_evects(ispin)%struct, nrow_global=nao, &
                                     ncol_global=nmo_occ(ispin), context=sub_env%blacs_env)
         END DO

         ALLOCATE (work_matrices%evects_sub(nspins, nstates), work_matrices%Aop_evects_sub(nspins, nstates))
         DO istate = 1, nstates
            DO ispin = 1, nspins
               NULLIFY (work_matrices%evects_sub(ispin, istate)%matrix)
               NULLIFY (work_matrices%Aop_evects_sub(ispin, istate)%matrix)
            END DO
         END DO

         CALL get_blacs_info(blacs_env, para_env=para_env)
         igroup = sub_env%group_distribution(para_env%mepos)
         ngroups = sub_env%ngroups

         DO istate = ngroups - igroup, nstates, ngroups
            DO ispin = 1, nspins
               CALL cp_fm_create(work_matrices%evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
               CALL cp_fm_create(work_matrices%Aop_evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
            END DO
         END DO

         DO ispin = nspins, 1, -1
            CALL cp_fm_struct_release(fm_struct_evects(ispin)%struct)
         END DO
      END IF

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=sub_env%blacs_env)
      NULLIFY (work_matrices%rho_ao_orb_fm_sub)
      CALL cp_fm_create(work_matrices%rho_ao_orb_fm_sub, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (work_matrices%rho_ao_aux_fit_fm_sub, work_matrices%wfm_aux_orb_sub)
      IF (do_admm) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao_aux, context=sub_env%blacs_env)
         CALL cp_fm_create(work_matrices%rho_ao_aux_fit_fm_sub, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao, context=sub_env%blacs_env)
         CALL cp_fm_create(work_matrices%wfm_aux_orb_sub, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END IF

      ! group-specific dbcsr matrices
      NULLIFY (work_matrices%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(work_matrices%A_ia_munu_sub, nspins)
      DO ispin = 1, nspins
         CALL dbcsr_init_p(work_matrices%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(work_matrices%A_ia_munu_sub(ispin)%matrix, template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_orb)
      END DO

      ! group-specific response density
      NULLIFY (rho_ia_ao)
      CALL dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
      DO ispin = 1, nspins
         CALL dbcsr_init_p(rho_ia_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(rho_ia_ao(ispin)%matrix, template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_orb)
      END DO

      NULLIFY (work_matrices%rho_orb_struct_sub)
      CALL qs_rho_create(work_matrices%rho_orb_struct_sub)
      CALL qs_rho_set(work_matrices%rho_orb_struct_sub, rho_ao=rho_ia_ao)
      CALL qs_rho_rebuild(work_matrices%rho_orb_struct_sub, qs_env, rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE., pw_env_external=sub_env%pw_env)

      NULLIFY (work_matrices%rho_aux_fit_struct_sub)
      IF (do_admm) THEN
         NULLIFY (rho_ia_ao)
         CALL dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(rho_ia_ao(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(rho_ia_ao(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                                             dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_aux_fit)
         END DO

         CALL qs_rho_create(work_matrices%rho_aux_fit_struct_sub)
         CALL qs_rho_set(work_matrices%rho_aux_fit_struct_sub, rho_ao=rho_ia_ao)
         CALL qs_rho_rebuild(work_matrices%rho_aux_fit_struct_sub, qs_env, rebuild_ao=.FALSE., &
                             rebuild_grids=.TRUE., pw_env_external=sub_env%pw_env)
      END IF

      ! work plain-wave grids
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      ALLOCATE (work_matrices%A_ia_rspace_sub(nspins))
      ALLOCATE (work_matrices%wpw_gspace_sub(nspins), work_matrices%wpw_rspace_sub(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%A_ia_rspace_sub(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%A_ia_rspace_sub(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)

         NULLIFY (work_matrices%wpw_gspace_sub(ispin)%pw, work_matrices%wpw_rspace_sub(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%wpw_gspace_sub(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%wpw_rspace_sub(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
      END DO

      ! HFX-related globally distributed matrices
      NULLIFY (work_matrices%hfx_fm_ao_ao, work_matrices%hfx_rho_ao, work_matrices%hfx_hmat)
      IF (do_hfx) THEN
         IF (do_admm) THEN
            CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, sab_aux_fit=sab_hfx)
            dbcsr_template_hfx => matrix_s_aux_fit(1)%matrix
         ELSE
            CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, sab_orb=sab_hfx)
            dbcsr_template_hfx => matrix_s(1)%matrix
         END IF

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(work_matrices%hfx_fm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL dbcsr_allocate_matrix_set(work_matrices%hfx_rho_ao, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(work_matrices%hfx_rho_ao(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(work_matrices%hfx_rho_ao(ispin)%matrix, &
                                             template=dbcsr_template_hfx, dbcsr_dist=dbcsr_dist, sab=sab_hfx)
         END DO

         CALL dbcsr_allocate_matrix_set(work_matrices%hfx_hmat, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(work_matrices%hfx_hmat(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(work_matrices%hfx_hmat(ispin)%matrix, &
                                             template=dbcsr_template_hfx, dbcsr_dist=dbcsr_dist, sab=sab_hfx)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_create_work_matrices

! **************************************************************************************************
!> \brief Release work matrices.
!> \param work_matrices  work matrices (destroyed on exit)
!> \param sub_env        parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_work_matrices(work_matrices, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_work_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! HFX-ralated matrices
      IF (ASSOCIATED(work_matrices%hfx_hmat)) THEN
         DO ispin = SIZE(work_matrices%hfx_hmat), 1, -1
            CALL dbcsr_deallocate_matrix(work_matrices%hfx_hmat(ispin)%matrix)
         END DO
         DEALLOCATE (work_matrices%hfx_hmat)
      END IF

      IF (ASSOCIATED(work_matrices%hfx_rho_ao)) THEN
         DO ispin = SIZE(work_matrices%hfx_rho_ao), 1, -1
            CALL dbcsr_deallocate_matrix(work_matrices%hfx_rho_ao(ispin)%matrix)
         END DO
         DEALLOCATE (work_matrices%hfx_rho_ao)
      END IF

      IF (ASSOCIATED(work_matrices%hfx_fm_ao_ao)) &
         CALL cp_fm_release(work_matrices%hfx_fm_ao_ao)

      ! real-space and reciprocal-space grids
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      DO ispin = SIZE(work_matrices%wpw_rspace_sub), 1, -1
         CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%wpw_rspace_sub(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%wpw_gspace_sub(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%A_ia_rspace_sub(ispin)%pw)
      END DO
      DEALLOCATE (work_matrices%A_ia_rspace_sub, work_matrices%wpw_gspace_sub, work_matrices%wpw_rspace_sub)

      IF (ASSOCIATED(work_matrices%rho_aux_fit_struct_sub)) &
         CALL qs_rho_release(work_matrices%rho_aux_fit_struct_sub)
      CALL qs_rho_release(work_matrices%rho_orb_struct_sub)

      DO ispin = SIZE(work_matrices%A_ia_munu_sub), 1, -1
         CALL dbcsr_deallocate_matrix(work_matrices%A_ia_munu_sub(ispin)%matrix)
      END DO
      DEALLOCATE (work_matrices%A_ia_munu_sub)

      IF (ASSOCIATED(work_matrices%wfm_aux_orb_sub)) &
         CALL cp_fm_release(work_matrices%wfm_aux_orb_sub)
      IF (ASSOCIATED(work_matrices%rho_ao_aux_fit_fm_sub)) &
         CALL cp_fm_release(work_matrices%rho_ao_aux_fit_fm_sub)
      CALL cp_fm_release(work_matrices%rho_ao_orb_fm_sub)

      IF (ALLOCATED(work_matrices%evects_sub)) THEN
         DO istate = SIZE(work_matrices%evects_sub, 2), 1, -1
            DO ispin = SIZE(work_matrices%evects_sub, 1), 1, -1
               CALL cp_fm_release(work_matrices%Aop_evects_sub(ispin, istate)%matrix)
               CALL cp_fm_release(work_matrices%evects_sub(ispin, istate)%matrix)
            END DO
         END DO
         DEALLOCATE (work_matrices%Aop_evects_sub, work_matrices%evects_sub)
      END IF

      DO ispin = SIZE(work_matrices%fm_pool_ao_mo_occ), 1, -1
         CALL cp_fm_release(work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%wfm_ao_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%S_C0(ispin)%matrix)
         CALL cp_fm_release(work_matrices%S_C0_C0T(ispin)%matrix)
      END DO
      DEALLOCATE (work_matrices%S_C0, work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)
      DEALLOCATE (work_matrices%wfm_mo_occ_mo_occ, work_matrices%wfm_mo_virt_mo_occ)

      DO ispin = SIZE(work_matrices%fm_pool_ao_mo_occ), 1, -1
         CALL fm_pool_release(work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
      END DO
      DEALLOCATE (work_matrices%fm_pool_ao_mo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_work_matrices

! **************************************************************************************************
!> \brief Create kernel environment.
!> \param kernel_env       kernel environment (allocated and initialised on exit)
!> \param rho_struct_sub   ground state charge density
!> \param xc_section       input section which defines an exchange-correlation functional
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param sub_env          parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
!>    * 06.2018 the charge density needs to be provided via a dummy argument [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_create_kernel_env(kernel_env, rho_struct_sub, xc_section, is_rks_triplets, sub_env)
      TYPE(tddfpt_kernel_env_type), INTENT(out)          :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct_sub
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_create_kernel_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      LOGICAL                                            :: lsd
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r, rho_ij_r2, tau_ij_r, tau_ij_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_fun_section

      CALL timeset(routineN, handle)

      nspins = SIZE(sub_env%mos_occ)
      lsd = (nspins > 1) .OR. is_rks_triplets

      DO ispin = 1, nspins
         CALL cp_fm_get_info(sub_env%mos_occ(ispin)%matrix, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL qs_rho_get(rho_struct_sub, rho_r=rho_ij_r, tau_r=tau_ij_r)

      NULLIFY (kernel_env%xc_rho_set, kernel_env%xc_rho1_set, kernel_env%xc_deriv_set)

      IF (is_rks_triplets) THEN
         ! we are about to compute triplet states using spin-restricted reference MOs;
         ! we still need the beta-spin density component in order to compute the TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1)%pw => rho_ij_r(1)%pw
         rho_ij_r2(2)%pw => rho_ij_r(1)%pw

         IF (ASSOCIATED(tau_ij_r)) THEN
            ALLOCATE (tau_ij_r2(2))
            tau_ij_r2(1)%pw => tau_ij_r(1)%pw
            tau_ij_r2(2)%pw => tau_ij_r(1)%pw
         END IF

         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r2, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r2)

         IF (ASSOCIATED(tau_ij_r)) &
            DEALLOCATE (tau_ij_r2)

         DEALLOCATE (rho_ij_r2)
      ELSE
         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r)
      END IF

      ! ++ allocate structure for response density
      kernel_env%xc_section => xc_section
      kernel_env%deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      kernel_env%rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      kernel_env%xc_rho1_cflags = xc_functionals_get_needs(functionals=xc_fun_section, lsd=lsd, add_basic_components=.TRUE.)

      CALL xc_rho_set_create(kernel_env%xc_rho1_set, auxbas_pw_pool%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      kernel_env%alpha = 1.0_dp
      kernel_env%beta = 0.0_dp

      ! kernel_env%beta is taken into account in spin-restricted case only
      IF (nspins == 1) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            kernel_env%beta = -1.0_dp
         ELSE
            !                                                 alpha                 beta
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            kernel_env%alpha = 2.0_dp
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_create_kernel_env

! **************************************************************************************************
!> \brief Release kernel environment.
!> \param kernel_env  kernel environment (destroyed on exit)
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_kernel_env(kernel_env)
      TYPE(tddfpt_kernel_env_type), INTENT(inout)        :: kernel_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_kernel_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL xc_rho_set_release(kernel_env%xc_rho1_set)
      CALL xc_dset_release(kernel_env%xc_deriv_set)
      CALL xc_rho_set_release(kernel_env%xc_rho_set)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_kernel_env

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param gs_mos           structure to store occupied and virtual molecular orbitals
!>                         (allocated and initialised on exit)
!> \param mo_set           ground state molecular orbitals for a given spin
!> \param nlumo            number of unoccupied states to consider (-1 means all states)
!> \param blacs_env        BLACS parallel environment
!> \param cholesky_method  Cholesky method to compute the inverse overlap matrix
!> \param matrix_ks        Kohn-Sham matrix for a given spin
!> \param matrix_s         overlap matrix
!> \param mos_virt         precomputed (OT) expansion coefficients of virtual molecular orbitals
!>                         (in addition to the ADDED_MOS, if present)
!> \param evals_virt       orbital energies of precomputed (OT) virtual molecular orbitals
!> \par History
!>    * 05.2016 created as tddfpt_lumos() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_init_ground_state_mos(gs_mos, mo_set, nlumo, blacs_env, cholesky_method, matrix_ks, matrix_s, &
                                           mos_virt, evals_virt)
      TYPE(tddfpt_ground_state_mos), INTENT(out)         :: gs_mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      INTEGER, INTENT(in)                                :: nlumo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(dbcsr_type), POINTER                          :: matrix_ks, matrix_s
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), DIMENSION(:), POINTER               :: evals_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_ground_state_mos', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER                           :: eps_dp = EPSILON(0.0_dp)

      INTEGER :: cholesky_method_inout, handle, icol_global, icol_local, imo, irow_global, &
         irow_local, nao, ncol_local, nelectrons, nmo_occ, nmo_scf, nmo_virt, nrow_local, sign_int
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: minrow_neg_array, minrow_pos_array, &
                                                            sum_sign_array
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_eigen
      REAL(kind=dp)                                      :: element, maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block
      TYPE(cp_fm_pool_type), POINTER                     :: wfn_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct, wfn_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      CALL get_blacs_info(blacs_env, para_env=para_env)

      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=nmo_occ, maxocc=maxocc, &
                      nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_virt = nao - nmo_occ
      IF (nlumo >= 0) &
         nmo_virt = MIN(nmo_virt, nlumo)

      IF (nmo_virt <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'At least one unoccupied molecular orbital is required to calculate excited states.')

      do_eigen = .FALSE.
      ! diagonalise the Kohn-Sham matrix one more time if the number of available unoccupied states are too small
      IF (ASSOCIATED(evals_virt)) THEN
         CPASSERT(ASSOCIATED(mos_virt))
         IF (nmo_virt > nmo_scf - nmo_occ + SIZE(evals_virt)) do_eigen = .TRUE.
      ELSE
         IF (nmo_virt > nmo_scf - nmo_occ) do_eigen = .TRUE.
      END IF

      ! ++ allocate storage space for gs_mos
      NULLIFY (ao_mo_occ_fm_struct, ao_mo_virt_fm_struct)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt, context=blacs_env)

      NULLIFY (gs_mos%mos_occ, gs_mos%mos_virt)
      CALL cp_fm_create(gs_mos%mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(gs_mos%mos_virt, ao_mo_virt_fm_struct)

      ALLOCATE (gs_mos%evals_occ(nmo_occ))
      ALLOCATE (gs_mos%evals_virt(nmo_virt))
      ALLOCATE (gs_mos%phases_occ(nmo_occ))

      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      ! ++ nullify pointers
      NULLIFY (ao_ao_fm_struct, wfn_fm_struct, wfn_fm_pool)
      NULLIFY (mos_extended, mo_coeff_extended, mo_evals_extended, mo_occ_extended)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)

      IF (do_eigen) THEN
         ! ++ set of molecular orbitals
         CALL cp_fm_struct_create(wfn_fm_struct, nrow_global=nao, ncol_global=nmo_occ + nmo_virt, context=blacs_env)
         CALL fm_pool_create(wfn_fm_pool, wfn_fm_struct)

         CALL allocate_mo_set(mos_extended, nao, nmo_occ + nmo_virt, nelectrons, &
                              REAL(nelectrons, dp), maxocc, flexible_electron_count=0.0_dp)
         CALL init_mo_set(mos_extended, fm_pool=wfn_fm_pool, name="mos-extended")
         CALL fm_pool_release(wfn_fm_pool)
         CALL cp_fm_struct_release(wfn_fm_struct)
         CALL get_mo_set(mos_extended, mo_coeff=mo_coeff_extended, &
                         eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)

         ! use the explicit loop in order to avoid temporary arrays.
         !
         ! The assignment statement : mo_occ_extended(1:nmo_scf) = mo_occ_scf(1:nmo_scf)
         ! implies temporary arrays as a compiler does not know in advance that the pointers
         ! on both sides of the statement point to non-overlapped memory regions
         DO imo = 1, nmo_scf
            mo_occ_extended(imo) = mo_occ_scf(imo)
         END DO
         mo_occ_extended(nmo_scf + 1:) = 0.0_dp

         ! ++ allocate temporary matrices
         NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
         CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
         CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
         CALL cp_fm_create(work_fm, ao_ao_fm_struct)

         CALL cp_fm_struct_release(ao_ao_fm_struct)

         ! some stuff from the subroutine general_eigenproblem()
         CALL copy_dbcsr_to_fm(matrix_s, ortho_fm)
         CALL copy_dbcsr_to_fm(matrix_ks, matrix_ks_fm)

         IF (cholesky_method == cholesky_dbcsr) THEN
            CPABORT('CHOLESKY DBCSR_INVERSE is not implemented in TDDFT.')
         ELSE IF (cholesky_method == cholesky_off) THEN
            CPABORT('CHOLESKY OFF is not implemented in TDDFT.')
         ELSE
            CALL cp_fm_cholesky_decompose(ortho_fm)
            IF (cholesky_method == cholesky_inverse) THEN
               CALL cp_fm_triangular_invert(ortho_fm)
            END IF

            ! need to store 'cholesky_method' in a temporary variable, as the subroutine eigensolver() will update this variable
            cholesky_method_inout = cholesky_method
            CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=mos_extended, ortho=ortho_fm, &
                             work=work_fm, cholesky_method=cholesky_method_inout, &
                             do_level_shift=.FALSE., level_shift=0.0_dp, matrix_u_fm=null(), use_jacobi=.FALSE.)
         END IF

         ! -- clean up needless matrices
         CALL cp_fm_release(work_fm)
         CALL cp_fm_release(ortho_fm)
         CALL cp_fm_release(matrix_ks_fm)
      ELSE
         CALL get_mo_set(mo_set, mo_coeff=mo_coeff_extended, &
                         eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)
      END IF

      ! compute the phase of molecular orbitals;
      ! matrix work_fm holds occupied molecular orbital coefficients distributed among all the processors
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env)
      CALL cp_fm_create(work_fm, ao_mo_occ_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      CALL cp_fm_to_fm(mo_coeff_extended, work_fm, ncol=nmo_occ, source_start=1, target_start=1)
      CALL cp_fm_get_info(work_fm, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=my_block)

      ALLOCATE (minrow_neg_array(nmo_occ), minrow_pos_array(nmo_occ), sum_sign_array(nmo_occ))
      minrow_neg_array(:) = nao
      minrow_pos_array(:) = nao
      sum_sign_array(:) = 0
      DO icol_local = 1, ncol_local
         icol_global = col_indices(icol_local)

         DO irow_local = 1, nrow_local
            element = my_block(irow_local, icol_local)

            sign_int = 0
            IF (element >= eps_dp) THEN
               sign_int = 1
            ELSE IF (element <= -eps_dp) THEN
               sign_int = -1
            END IF

            sum_sign_array(icol_global) = sum_sign_array(icol_global) + sign_int

            irow_global = row_indices(irow_local)
            IF (sign_int > 0) THEN
               IF (minrow_pos_array(icol_global) > irow_global) &
                  minrow_pos_array(icol_global) = irow_global
            ELSE IF (sign_int < 0) THEN
               IF (minrow_neg_array(icol_global) > irow_global) &
                  minrow_neg_array(icol_global) = irow_global
            END IF
         END DO
      END DO

      CALL mp_sum(sum_sign_array, para_env%group)
      CALL mp_min(minrow_neg_array, para_env%group)
      CALL mp_min(minrow_pos_array, para_env%group)

      DO icol_local = 1, nmo_occ
         IF (sum_sign_array(icol_local) > 0) THEN
            ! most of the expansion coefficients are positive => MO's phase = +1
            gs_mos%phases_occ(icol_local) = 1.0_dp
         ELSE IF (sum_sign_array(icol_local) < 0) THEN
            ! most of the expansion coefficients are negative => MO's phase = -1
            gs_mos%phases_occ(icol_local) = -1.0_dp
         ELSE
            ! equal number of positive and negative expansion coefficients
            IF (minrow_pos_array(icol_local) <= minrow_neg_array(icol_local)) THEN
               ! the first positive expansion coefficient has a lower index then
               ! the first negative expansion coefficient; MO's phase = +1
               gs_mos%phases_occ(icol_local) = 1.0_dp
            ELSE
               ! MO's phase = -1
               gs_mos%phases_occ(icol_local) = -1.0_dp
            END IF
         END IF
      END DO

      DEALLOCATE (minrow_neg_array, minrow_pos_array, sum_sign_array)
      CALL cp_fm_release(work_fm)

      ! return the requested occupied and virtual molecular orbitals and corresponding orbital energies
      CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_occ, ncol=nmo_occ, source_start=1, target_start=1)
      gs_mos%evals_occ(1:nmo_occ) = mo_evals_extended(1:nmo_occ)

      IF (ASSOCIATED(evals_virt) .AND. (.NOT. do_eigen) .AND. nmo_virt > nmo_scf - nmo_occ) THEN
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_virt, ncol=nmo_scf - nmo_occ, &
                          source_start=nmo_occ + 1, target_start=1)
         CALL cp_fm_to_fm(mos_virt, gs_mos%mos_virt, ncol=nmo_virt - (nmo_scf - nmo_occ), &
                          source_start=1, target_start=nmo_scf - nmo_occ + 1)
         gs_mos%evals_virt(1:nmo_scf - nmo_occ) = evals_virt(nmo_occ + 1:nmo_occ + nmo_scf)
         gs_mos%evals_virt(nmo_scf - nmo_occ + 1:nmo_virt) = evals_virt(1:nmo_virt - (nmo_scf - nmo_occ))
      ELSE
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_virt, ncol=nmo_virt, source_start=nmo_occ + 1, target_start=1)
         gs_mos%evals_virt(1:nmo_virt) = mo_evals_extended(nmo_occ + 1:nmo_occ + nmo_virt)
      END IF

      IF (do_eigen) &
         CALL deallocate_mo_set(mos_extended)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_ground_state_mos

! **************************************************************************************************
!> \brief Release molecular orbitals.
!> \param gs_mos  structure that holds occupied and virtual molecular orbitals
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_ground_state_mos(gs_mos)
      TYPE(tddfpt_ground_state_mos), INTENT(inout)       :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_ground_state_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(gs_mos%phases_occ)) &
         DEALLOCATE (gs_mos%phases_occ)

      IF (ALLOCATED(gs_mos%evals_virt)) &
         DEALLOCATE (gs_mos%evals_virt)

      IF (ALLOCATED(gs_mos%evals_occ)) &
         DEALLOCATE (gs_mos%evals_occ)

      IF (ASSOCIATED(gs_mos%mos_virt)) &
         CALL cp_fm_release(gs_mos%mos_virt)

      IF (ASSOCIATED(gs_mos%mos_occ)) &
         CALL cp_fm_release(gs_mos%mos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_ground_state_mos

! **************************************************************************************************
!> \brief Compute the number of possible singly excited states (occ -> virt)
!> \param gs_mos          occupied and virtual molecular orbitals optimised for the ground state
!> \return the number of possible single excitations
!> \par History
!>    * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   PURE FUNCTION tddfpt_total_number_of_states(gs_mos) RESULT(nstates_total)
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER(kind=int_8)                                :: nstates_total

      INTEGER                                            :: ispin, nspins

      nstates_total = 0
      nspins = SIZE(gs_mos)

      DO ispin = 1, nspins
         nstates_total = nstates_total + &
                         SIZE(gs_mos(ispin)%evals_occ, kind=int_8)* &
                         SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO
   END FUNCTION tddfpt_total_number_of_states

! **************************************************************************************************
!> \brief Generate missed guess vectors.
!> \param evects   guess vectors distributed across all processors (initialised on exit)
!> \param evals    guessed transition energies (initialised on exit)
!> \param gs_mos   occupied and virtual molecular orbitals optimised for the ground state
!> \param log_unit output unit
!> \par History
!>    * 05.2016 created as tddfpt_guess() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype, supports spin-polarised density [Sergey Chulkov]
!>    * 01.2017 simplified prototype, do not compute all possible singly-excited states
!>              [Sergey Chulkov]
!> \note \parblock
!>       Based on the subroutine co_initial_guess() which was originally created by
!>       Thomas Chassaing on 06.2003.
!>
!>       Only not associated guess vectors 'evects(spin, state)%matrix' are allocated and
!>       initialised; associated vectors assumed to be initialised elsewhere (e.g. using
!>       a restart file).
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_guess_vectors(evects, evals, gs_mos, log_unit)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_guess_vectors', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=5)                                   :: spin_label
      INTEGER                                            :: handle, imo_occ, imo_virt, ind, ispin, &
                                                            istate, jspin, nspins, nstates, &
                                                            nstates_occ_virt_alpha, &
                                                            nstates_selected
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ_avail, nmo_occ_selected, &
                                                            nmo_virt_selected
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nstates > 0)
         CPASSERT(nspins == 1 .OR. nspins == 2)
         CPASSERT(SIZE(gs_mos) == nspins)
      END IF

      DO ispin = 1, nspins
         ! number of occupied orbitals for each spin component
         nmo_occ_avail(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         ! number of occupied and virtual orbitals which can potentially
         ! contribute to the excited states in question.
         nmo_occ_selected(ispin) = MIN(nmo_occ_avail(ispin), nstates)
         nmo_virt_selected(ispin) = MIN(SIZE(gs_mos(ispin)%evals_virt), nstates)

         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct_evects(ispin)%struct)
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = DOT_PRODUCT(nmo_occ_selected(1:nspins), nmo_virt_selected(1:nspins))

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO ispin = 1, nspins
         DO imo_occ = 1, nmo_occ_selected(ispin)
            ! Here imo_occ enumerate Occupied orbitals in inverse order (from the last to the first element)
            e_occ = gs_mos(ispin)%evals_occ(nmo_occ_avail(ispin) - imo_occ + 1)

            DO imo_virt = 1, nmo_virt_selected(ispin)
               istate = istate + 1
               e_virt_minus_occ(istate) = gs_mos(ispin)%evals_virt(imo_virt) - e_occ
            END DO
         END DO
      END DO

      IF (debug_this_module) THEN
         CPASSERT(istate == nstates_selected)
      END IF

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      IF (nspins == 1) THEN
         ispin = 1
         spin_label = '     '
      END IF

      nstates_occ_virt_alpha = nmo_occ_selected(1)*nmo_virt_selected(1)
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,21X,A)') "TDDFPT initial guess"
         WRITE (log_unit, '(1X,60("-"))')
         WRITE (log_unit, '(5X,A)') "State       Occupied    ->    Virtual        Excitation"
         WRITE (log_unit, '(5X,A)') "number       orbital          orbital        energy (eV)"
         WRITE (log_unit, '(1X,60("-"))')
      END IF

      DO istate = 1, nstates
         IF (ASSOCIATED(evects(1, istate)%matrix)) THEN
            IF (log_unit > 0) &
               WRITE (log_unit, '(1X,I8,11X,A19,8X,F14.5)') &
               istate, "***  restarted  ***", evals(istate)*evolt
         ELSE
            ind = inds(istate) - 1
            IF (nspins > 1) THEN
               IF (ind < nstates_occ_virt_alpha) THEN
                  ispin = 1
                  spin_label = '(alp)'
               ELSE
                  ispin = 2
                  ind = ind - nstates_occ_virt_alpha
                  spin_label = '(bet)'
               END IF
            END IF

            imo_occ = nmo_occ_avail(ispin) - ind/nmo_virt_selected(ispin)
            imo_virt = MOD(ind, nmo_virt_selected(ispin)) + 1
            evals(istate) = e_virt_minus_occ(istate)

            IF (log_unit > 0) &
               WRITE (log_unit, '(1X,I8,5X,I8,1X,A5,3X,I8,1X,A5,2X,F14.5)') &
               istate, imo_occ, spin_label, nmo_occ_avail(ispin) + imo_virt, spin_label, e_virt_minus_occ(istate)*evolt

            DO jspin = 1, nspins
               ! .NOT. ASSOCIATED(evects(jspin, istate)%matrix))
               CALL cp_fm_create(evects(jspin, istate)%matrix, fm_struct_evects(jspin)%struct)
               CALL cp_fm_set_all(evects(jspin, istate)%matrix, 0.0_dp)

               IF (jspin == ispin) &
                  CALL cp_fm_to_fm(gs_mos(ispin)%mos_virt, evects(ispin, istate)%matrix, &
                                   ncol=1, source_start=imo_virt, target_start=imo_occ)
            END DO
         END IF
      END DO

      IF (log_unit > 0) &
         WRITE (log_unit, '(/,1X,A,T30,I24)') 'Number of active states:', tddfpt_total_number_of_states(gs_mos)

      DEALLOCATE (e_virt_minus_occ)
      DEALLOCATE (inds)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_guess_vectors

! **************************************************************************************************
!> \brief Make TDDFPT trial vectors orthogonal to all occupied molecular orbitals.
!> \param evects            trial vectors distributed across all processors (modified on exit)
!> \param S_C0_C0T          matrix product S * C_0 * C_0^T, where C_0 is the ground state
!>                          wave function for each spin expressed in atomic basis set,
!>                          and S is the corresponding overlap matrix
!> \param work_fm_ao_mo_occ work matrices with shape [nao x nmo_occ(spin)] to store matrix products
!>                          C_0 * C_0^T * S * evects (modified on exit)
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note  Based on the subroutine p_preortho() which was created by Thomas Chassaing on 09.2002.
!>        Should be useless when ground state MOs are computed with extremely high accuracy,
!>        as all virtual orbitals are already orthogonal to the occupied ones by design.
!>        However, when the norm of residual vectors is relatively small (e.g. less then SCF_EPS),
!>        new Krylov's vectors seem to be random and should be orthogonalised even with respect to
!>        the occupied MOs.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthogonalize_psi1_psi0(evects, S_C0_C0T, work_fm_ao_mo_occ)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0_C0T, work_fm_ao_mo_occ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthogonalize_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0 .AND. nspins <= maxspins)
         CPASSERT(SIZE(S_C0_C0T) == nspins)
         CPASSERT(SIZE(work_fm_ao_mo_occ) == nspins)
      END IF

      CALL cp_fm_get_info(matrix=work_fm_ao_mo_occ(1)%matrix, nrow_global=nao, ncol_global=nmo_occ(1))
      DO ispin = 2, nspins
         CALL cp_fm_get_info(matrix=work_fm_ao_mo_occ(ispin)%matrix, ncol_global=nmo_occ(ispin))
      END DO

      IF (nvects > 0) THEN

         DO ivect = 1, nvects
            DO ispin = 1, nspins
               ! work_fm_ao_mo_occ: C0 * C0^T * S * C1 == (S * C0 * C0^T)^T * C1
               CALL cp_gemm('T', 'N', nao, nmo_occ(ispin), nao, 1.0_dp, S_C0_C0T(ispin)%matrix, &
                            evects(ispin, ivect)%matrix, 0.0_dp, work_fm_ao_mo_occ(ispin)%matrix)

               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, ivect)%matrix, -1.0_dp, work_fm_ao_mo_occ(ispin)%matrix)
            END DO
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthogonalize_psi1_psi0

! **************************************************************************************************
!> \brief Check that orthogonalised TDDFPT trial vectors remain orthogonal to
!>        occupied molecular orbitals.
!> \param evects    trial vectors
!> \param S_C0      matrix product S * C_0, where C_0 is the ground state wave function
!>                  for each spin in atomic basis set, and S is the corresponding overlap matrix
!> \param max_norm  the largest possible overlap between the ground state and
!>                  excited state wave functions
!> \param work_fm_mo_occ_mo_occ work matrices with shape [nmo_occ(spin) x nmo_occ(spin)]
!>                              to store matrix products S_0^T * S * evects (modified on exit)
!> \return true if trial vectors are non-orthogonal to occupied molecular orbitals
!> \par History
!>    * 07.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_is_nonorthogonal_psi1_psi0(evects, S_C0, max_norm, work_fm_mo_occ_mo_occ) RESULT(is_nonortho)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0
      REAL(kind=dp), INTENT(in)                          :: max_norm
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: work_fm_mo_occ_mo_occ
      LOGICAL                                            :: is_nonortho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_is_nonorthogonal_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      REAL(kind=dp)                                      :: maxabs_val

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0 .AND. nspins <= maxspins)
         CPASSERT(SIZE(S_C0) == nspins)
         CPASSERT(SIZE(work_fm_mo_occ_mo_occ) == nspins)
      END IF

      CALL cp_fm_get_info(matrix=S_C0(1)%matrix, nrow_global=nao, ncol_global=nmo_occ(1))
      DO ispin = 2, nspins
         CALL cp_fm_get_info(matrix=S_C0(ispin)%matrix, ncol_global=nmo_occ(ispin))
      END DO

      is_nonortho = .FALSE.

      loop: DO ivect = 1, nvects
         DO ispin = 1, nspins
            ! work_fm_mo_occ_mo_occ = S_0^T * S * C_1
            CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, 1.0_dp, S_C0(ispin)%matrix, &
                         evects(ispin, ivect)%matrix, 0.0_dp, work_fm_mo_occ_mo_occ(ispin)%matrix)

            CALL cp_fm_maxabsval(work_fm_mo_occ_mo_occ(ispin)%matrix, maxabs_val)
            is_nonortho = maxabs_val > max_norm
            IF (is_nonortho) EXIT loop
         END DO
      END DO loop

      CALL timestop(handle)
   END FUNCTION tddfpt_is_nonorthogonal_psi1_psi0

! **************************************************************************************************
!> \brief Make new TDDFPT trial vectors orthonormal to all previous trial vectors.
!> \param evects      trial vectors (modified on exit)
!> \param nvects_new  number of new trial vectors to orthogonalise
!> \param S_evects    set of matrices to store matrix product S * evects (modified on exit)
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 02.2017 caching the matrix product S * evects [Sergey Chulkov]
!> \note \parblock
!>       Based on the subroutines reorthogonalize() and normalize() which were originally created
!>       by Thomas Chassaing on 03.2003.
!>
!>       In order to orthogonalise a trial vector C3 = evects(:,3) with respect to previously
!>       orthogonalised vectors C1 = evects(:,1) and C2 = evects(:,2) we need to compute the
!>       quantity C3'' using the following formulae:
!>          C3'  = C3  - Tr(C3^T  * S * C1) * C1,
!>          C3'' = C3' - Tr(C3'^T * S * C2) * C2,
!>       which can be expanded as:
!>          C3'' = C3 - Tr(C3^T  * S * C1) * C1 - Tr(C3^T * S * C2) * C2 +
!>                 Tr(C3^T * S * C1) * Tr(C2^T * S * C1) * C2 .
!>       In case of unlimited float-point precision, the last term in above expression is exactly 0,
!>       due to orthogonality condition between C1 and C2. In this case the expression could be
!>       simplified as (taking into account the identity: Tr(A * S * B) = Tr(B * S * A)):
!>          C3'' = C3 - Tr(C1^T  * S * C3) * C1 - Tr(C2^T * S * C3) * C2 ,
!>       which means we do not need the variable S_evects to keep the matrix products S * Ci .
!>
!>       In reality, however, we deal with limited float-point precision arithmetic meaning that
!>       the trace Tr(C2^T * S * C1) is close to 0 but does not equal to 0 exactly. The term
!>          Tr(C3^T * S * C1) * Tr(C2^T * S * C1) * C2
!>       can not be ignored anymore. Ignorance of this term will lead to numerical instability
!>       when the trace Tr(C3^T * S * C1) is large enough.
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_orthonormalize_psi1_psi1(evects, nvects_new, S_evects, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      INTEGER, INTENT(in)                                :: nvects_new
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: S_evects
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthonormalize_psi1_psi1', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, jvect, nspins, &
                                                            nvects_old, nvects_total
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), DIMENSION(maxspins)                 :: weights

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects_total = SIZE(evects, 2)
      nvects_old = nvects_total - nvects_new

      IF (debug_this_module) THEN
         CPASSERT(SIZE(S_evects, 1) == nspins)
         CPASSERT(SIZE(S_evects, 2) == nvects_total)
         CPASSERT(nvects_old >= 0)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, ncol_global=nmo_occ(ispin))
      END DO

      DO jvect = nvects_old + 1, nvects_total
         ! <psi1_i | psi1_j>
         DO ivect = 1, jvect - 1
            CALL cp_fm_trace(evects(:, jvect), S_evects(:, ivect), weights(1:nspins))
            norm = accurate_sum(weights(1:nspins))

            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, jvect)%matrix, -norm, evects(ispin, ivect)%matrix)
            END DO
         END DO

         ! <psi1_j | psi1_j>
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, jvect)%matrix, S_evects(ispin, jvect)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         END DO

         CALL cp_fm_trace(evects(:, jvect), S_evects(:, jvect), weights(1:nspins))

         norm = accurate_sum(weights(1:nspins))
         norm = 1.0_dp/SQRT(norm)

         DO ispin = 1, nspins
            CALL cp_fm_scale(norm, evects(ispin, jvect)%matrix)
            CALL cp_fm_scale(norm, S_evects(ispin, jvect)%matrix)
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************
!> \brief Apply orbital energy difference term:
!>        Aop_evects(spin,state) += KS(spin) * evects(spin,state) -
!>                                  S * evects(spin,state) * diag(evals_occ(spin))
!> \param Aop_evects  action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects      trial vectors C_{1,i}
!> \param S_evects    S * C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied orbital
!>                    energies [component %evals_occ] are needed)
!> \param matrix_ks   Kohn-Sham matrix
!> \param work_fm_ao_mo_occ work matrices with shape [nao x nmo_occ(spin)] to store matrix
!>                          products S * evects * diag(evals_occ).
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!>    * 03.2017 renamed from tddfpt_init_energy_diff(), altered prototype [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_energy_diff(Aop_evects, evects, S_evects, gs_mos, matrix_ks, work_fm_ao_mo_occ)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects, S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_ks
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: work_fm_ao_mo_occ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_energy_diff', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, evects(ispin, ivect)%matrix, &
                                         Aop_evects(ispin, ivect)%matrix, ncol=nmo_occ(ispin), &
                                         alpha=1.0_dp, beta=1.0_dp)

            CALL cp_fm_to_fm(S_evects(ispin, ivect)%matrix, work_fm_ao_mo_occ(ispin)%matrix)
            CALL cp_fm_column_scale(work_fm_ao_mo_occ(ispin)%matrix, gs_mos(ispin)%evals_occ)

            ! KS * C1 - S * C1 * occupied_orbital_energies
            CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect)%matrix, -1.0_dp, work_fm_ao_mo_occ(ispin)%matrix)
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_energy_diff

! **************************************************************************************************
!> \brief Update v_rspace by adding coulomb term.
!> \param A_ia_rspace    action of TDDFPT operator on the trial vector expressed in a plane wave
!>                       representation (modified on exit)
!> \param rho_ia_g       response density in reciprocal space for the given trial vector
!> \param pw_env         plain wave environment
!> \param work_v_gspace  work reciprocal-space grid to store Coulomb potential (modified on exit)
!> \param work_v_rspace  work real-space grid to store Coulomb potential (modified on exit)
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!>    * 03.2017 proceed excited states sequentially; minimise the number of conversions between
!>              DBCSR and FM matrices [Sergey Chulkov]
!>    * 06.2018 return the action expressed in the plane wave representation instead of the one
!>              in the atomic basis set representation
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_coulomb(A_ia_rspace, rho_ia_g, pw_env, work_v_gspace, work_v_rspace)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(pw_type), POINTER                             :: rho_ia_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), INTENT(inout)                     :: work_v_gspace, work_v_rspace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(kind=dp)                                      :: alpha, pair_energy
      TYPE(pw_poisson_type), POINTER                     :: poisson_env

      CALL timeset(routineN, handle)

      nspins = SIZE(A_ia_rspace)
      CALL pw_env_get(pw_env, poisson_env=poisson_env)

      IF (nspins > 1) THEN
         alpha = 1.0_dp
      ELSE
         ! spin-restricted case: alpha == 2 due to singlet state.
         ! In case of triplet states alpha == 0, so we should not call this subroutine at all.
         alpha = 2.0_dp
      END IF

      CALL pw_poisson_solve(poisson_env, rho_ia_g, pair_energy, work_v_gspace%pw)
      CALL pw_transfer(work_v_gspace%pw, work_v_rspace%pw)

      ! (i a || j b) = ( i_alpha a_alpha + i_beta a_beta || j_alpha b_alpha + j_beta b_beta) =
      !                tr (Cj_alpha^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_alpha) +
      !                tr (Cj_beta^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_beta)
      DO ispin = 1, nspins
         CALL pw_axpy(work_v_rspace%pw, A_ia_rspace(ispin)%pw, alpha)
      END DO

      ! TO DO: tau component

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_coulomb

! **************************************************************************************************
!> \brief Update A_ia_munu by adding exchange-correlation term.
!> \param A_ia_rspace      action of TDDFPT operator on trial vectors expressed in a plane wave
!>                         representation (modified on exit)
!> \param kernel_env       kernel environment
!> \param rho_ia_struct    response density for the given trial vector
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param pw_env           plain wave environment
!> \param work_v_xc        work real-space grid to store the gradient of the exchange-correlation
!>                         potential with respect to the response density (modified on exit)
!> \par History
!>    * 05.2016 compute all kernel terms in one go [Sergey Chulkov]
!>    * 03.2017 proceed excited states sequentially; minimise the number of conversions between
!>              DBCSR and FM matrices [Sergey Chulkov]
!>    * 06.2018 return the action expressed in the plane wave representation instead of the one
!>              in the atomic basis set representation
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_xc(A_ia_rspace, kernel_env, rho_ia_struct, is_rks_triplets, pw_env, work_v_xc)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace
      TYPE(tddfpt_kernel_env_type), INTENT(in)           :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: work_v_xc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_xc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g2, rho_ia_r, &
                                                            rho_ia_r2, tau_ia_r, tau_ia_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      nspins = SIZE(A_ia_rspace)
      CALL qs_rho_get(rho_ia_struct, rho_g=rho_ia_g, rho_r=rho_ia_r, tau_r=tau_ia_r)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(rho_ia_g) == nspins)
         CPASSERT(SIZE(rho_ia_r) == nspins)
         CPASSERT((.NOT. ASSOCIATED(tau_ia_r)) .OR. SIZE(tau_ia_r) == nspins)
         CPASSERT((.NOT. is_rks_triplets) .OR. nspins == 1)
      END IF

      NULLIFY (tau_ia_r2)
      IF (is_rks_triplets) THEN
         ALLOCATE (rho_ia_r2(2))
         ALLOCATE (rho_ia_g2(2))
         rho_ia_r2(1)%pw => rho_ia_r(1)%pw
         rho_ia_r2(2)%pw => rho_ia_r(1)%pw
         rho_ia_g2(1)%pw => rho_ia_g(1)%pw
         rho_ia_g2(2)%pw => rho_ia_g(1)%pw

         IF (ASSOCIATED(tau_ia_r)) THEN
            ALLOCATE (tau_ia_r2(2))
            tau_ia_r2(1)%pw => tau_ia_r(1)%pw
            tau_ia_r2(2)%pw => tau_ia_r(1)%pw
         END IF
      ELSE
         ALLOCATE (rho_ia_r2(nspins))
         ALLOCATE (rho_ia_g2(nspins))
         DO ispin = 1, nspins
            rho_ia_r2(ispin)%pw => rho_ia_r(ispin)%pw
            rho_ia_g2(ispin)%pw => rho_ia_g(ispin)%pw
         END DO

         IF (ASSOCIATED(tau_ia_r)) THEN
            ALLOCATE (tau_ia_r2(nspins))
            DO ispin = 1, nspins
               tau_ia_r2(ispin)%pw => tau_ia_r(ispin)%pw
            END DO
         END IF
      END IF

      DO ispin = 1, nspins
         CALL pw_zero(work_v_xc(ispin)%pw)
      END DO

      CALL xc_rho_set_update(rho_set=kernel_env%xc_rho1_set, rho_r=rho_ia_r2, rho_g=rho_ia_g2, tau=tau_ia_r2, &
                             needs=kernel_env%xc_rho1_cflags, xc_deriv_method_id=kernel_env%deriv_method_id, &
                             xc_rho_smooth_id=kernel_env%rho_smooth_id, pw_pool=auxbas_pw_pool)

      CALL xc_calc_2nd_deriv(v_xc=work_v_xc, deriv_set=kernel_env%xc_deriv_set, rho_set=kernel_env%xc_rho_set, &
                             rho1_set=kernel_env%xc_rho1_set, pw_pool=auxbas_pw_pool, &
                             xc_section=kernel_env%xc_section, gapw=.FALSE., tddfpt_fac=kernel_env%beta)

      DO ispin = 1, nspins
         ! pw2 = pw2 + alpha * pw1
         CALL pw_axpy(work_v_xc(ispin)%pw, A_ia_rspace(ispin)%pw, kernel_env%alpha)
      END DO

      DEALLOCATE (rho_ia_g2, rho_ia_r2)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_xc

! **************************************************************************************************
!> \brief Compute the ground-state charge density expressed in primary basis set.
!> \param rho_orb_struct      ground-state density in primary basis set
!> \param is_rks_triplets     indicates that the triplet excited states calculation using
!>                            spin-unpolarised molecular orbitals has been requested
!> \param qs_env              Quickstep environment
!> \param sub_env             parallel (sub)group environment
!> \param wfm_rho_orb         work dense matrix with shape [nao x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \par History
!>    * 06.2018 created by splitting the subroutine tddfpt_apply_admm_correction() in two
!>              subroutines tddfpt_construct_ground_state_orb_density() and
!>              tddfpt_construct_aux_fit_density [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_construct_ground_state_orb_density(rho_orb_struct, is_rks_triplets, qs_env, sub_env, wfm_rho_orb)
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_fm_type), POINTER                          :: wfm_rho_orb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_construct_ground_state_orb_density', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ij_ao

      CALL timeset(routineN, handle)

      nspins = SIZE(sub_env%mos_occ)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(sub_env%mos_occ(ispin)%matrix, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      CALL qs_rho_get(rho_orb_struct, rho_ao=rho_ij_ao)
      DO ispin = 1, nspins
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, &
                      sub_env%mos_occ(ispin)%matrix, sub_env%mos_occ(ispin)%matrix, &
                      0.0_dp, wfm_rho_orb)

         CALL copy_fm_to_dbcsr(wfm_rho_orb, rho_ij_ao(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      ! take into account that all MOs are doubly occupied in spin-restricted case
      IF (nspins == 1 .AND. (.NOT. is_rks_triplets)) &
         CALL dbcsr_scale(rho_ij_ao(1)%matrix, 2.0_dp)

      CALL qs_rho_update_rho(rho_orb_struct, qs_env, &
                             pw_env_external=sub_env%pw_env, &
                             task_list_external=sub_env%task_list_orb)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_construct_ground_state_orb_density

! **************************************************************************************************
!> \brief Project a charge density expressed in primary basis set into the auxiliary basis set.
!> \param rho_orb_struct      response density in primary basis set
!> \param rho_aux_fit_struct  response density in auxiliary basis set (modified on exit)
!> \param qs_env              Quickstep environment
!> \param sub_env             parallel (sub)group environment
!> \param wfm_rho_orb         work dense matrix with shape [nao x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \param wfm_rho_aux_fit     work dense matrix with shape [nao_aux x nao_aux] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \param wfm_aux_orb         work dense matrix with shape [nao_aux x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \par History
!>    * 03.2017 the subroutine tddfpt_apply_admm_correction() was originally created by splitting
!>              the subroutine tddfpt_apply_hfx() in two parts [Sergey Chulkov]
!>    * 06.2018 created by splitting the subroutine tddfpt_apply_admm_correction() in two subroutines
!>              tddfpt_construct_ground_state_orb_density() and tddfpt_construct_aux_fit_density()
!>              in order to avoid code duplication [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_construct_aux_fit_density(rho_orb_struct, rho_aux_fit_struct, qs_env, sub_env, &
                                               wfm_rho_orb, wfm_rho_aux_fit, wfm_aux_orb)
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct, rho_aux_fit_struct
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_fm_type), POINTER                          :: wfm_rho_orb, wfm_rho_aux_fit, wfm_aux_orb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_construct_aux_fit_density', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      REAL(kind=dp), DIMENSION(:), POINTER               :: tot_rho_aux_fit_r
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_aux_fit, rho_ao_orb
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_aux_fit_g, rho_aux_fit_r
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(sub_env%admm_A))

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL qs_rho_get(rho_orb_struct, rho_ao=rho_ao_orb)
      CALL qs_rho_get(rho_aux_fit_struct, rho_ao=rho_ao_aux_fit, rho_g=rho_aux_fit_g, &
                      rho_r=rho_aux_fit_r, tot_rho_r=tot_rho_aux_fit_r)

      nspins = SIZE(rho_ao_orb)

      CALL cp_fm_get_info(sub_env%admm_A, nrow_global=nao_aux, ncol_global=nao)
      DO ispin = 1, nspins
         ! TO DO: consider sub_env%admm_A to be a DBCSR matrix
         CALL copy_dbcsr_to_fm(rho_ao_orb(ispin)%matrix, wfm_rho_orb)
         CALL cp_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, sub_env%admm_A, &
                      wfm_rho_orb, 0.0_dp, wfm_aux_orb)
         CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, sub_env%admm_A, wfm_aux_orb, &
                      0.0_dp, wfm_rho_aux_fit)
         CALL copy_fm_to_dbcsr(wfm_rho_aux_fit, rho_ao_aux_fit(ispin)%matrix, keep_sparsity=.TRUE.)

         CALL calculate_rho_elec(matrix_p=rho_ao_aux_fit(ispin)%matrix, &
                                 rho=rho_aux_fit_r(ispin), rho_gspace=rho_aux_fit_g(ispin), &
                                 total_rho=tot_rho_aux_fit_r(ispin), ks_env=ks_env, &
                                 soft_valid=.FALSE., basis_type="AUX_FIT", &
                                 pw_env_external=sub_env%pw_env, task_list_external=sub_env%task_list_aux_fit)
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_construct_aux_fit_density

! **************************************************************************************************
!> \brief Update action of TDDFPT operator on trial vectors by adding exact-exchange term.
!> \param Aop_evects      action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects          trial vectors
!> \param gs_mos          molecular orbitals optimised for the ground state (only occupied
!>                        molecular orbitals [component %mos_occ] are needed)
!> \param do_admm         perform auxiliary density matrix method calculations
!> \param qs_env          Quickstep environment
!> \param work_rho_ia_ao  work sparse matrix with shape [nao x nao] distributed globally
!>                        to store response density (modified on exit)
!> \param work_hmat       work sparse matrix with shape [nao x nao] distributed globally
!>                        (modified on exit)
!> \param wfm_rho_orb     work dense matrix with shape [nao x nao] distributed globally
!>                        (modified on exit)
!> \par History
!>    * 05.2016 compute all exact-exchange terms in one go [Sergey Chulkov]
!>    * 03.2017 code related to ADMM correction is now moved to tddfpt_apply_admm_correction()
!>              in order to compute this correction within parallel groups [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_hfx(Aop_evects, evects, gs_mos, do_admm, qs_env, &
                               work_rho_ia_ao, work_hmat, wfm_rho_orb)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: do_admm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: work_rho_ia_ao, work_hmat
      TYPE(cp_fm_type), POINTER                          :: wfm_rho_orb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_hfx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nao_aux, &
                                                            nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      REAL(kind=dp)                                      :: alpha
      TYPE(admm_type), POINTER                           :: admm_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (nspins > 1) THEN
         alpha = 2.0_dp
      ELSE
         alpha = 4.0_dp
      END IF

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      IF (do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
      END IF

      ! some stuff from qs_ks_build_kohn_sham_matrix
      ! TO DO: add SIC support
      DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, gs_mos(ispin)%mos_occ, &
                         evects(ispin, ivect)%matrix, 0.0_dp, wfm_rho_orb)
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, evects(ispin, ivect)%matrix, &
                         gs_mos(ispin)%mos_occ, 1.0_dp, wfm_rho_orb)

            CALL dbcsr_set(work_hmat(ispin)%matrix, 0.0_dp)
            IF (do_admm) THEN
               CALL cp_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, admm_env%A, &
                            wfm_rho_orb, 0.0_dp, admm_env%work_aux_orb)
               CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%A, admm_env%work_aux_orb, &
                            0.0_dp, admm_env%work_aux_aux)
               CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, work_rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
            ELSE
               CALL copy_fm_to_dbcsr(wfm_rho_orb, work_rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
            END IF
         END DO

         CALL tddft_hfx_matrix(work_hmat, work_rho_ia_ao, qs_env)

         IF (do_admm) THEN
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(work_hmat(ispin)%matrix, admm_env%A, admm_env%work_aux_orb, &
                                            ncol=nao, alpha=1.0_dp, beta=0.0_dp)

               CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, &
                            admm_env%work_aux_orb, 0.0_dp, wfm_rho_orb)

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, alpha, wfm_rho_orb, &
                            gs_mos(ispin)%mos_occ, 1.0_dp, Aop_evects(ispin, ivect)%matrix)
            END DO
         ELSE
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(work_hmat(ispin)%matrix, gs_mos(ispin)%mos_occ, &
                                            Aop_evects(ispin, ivect)%matrix, ncol=nmo_occ(ispin), &
                                            alpha=alpha, beta=1.0_dp)
            END DO
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_hfx

! **************************************************************************************************
!> \brief Compute action matrix-vector products.
!> \param Aop_evects            action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects                TDDFPT trial vectors
!> \param S_evects              cached matrix product S * evects where S is the overlap matrix
!>                              in primary basis set
!> \param gs_mos                molecular orbitals optimised for the ground state
!> \param is_rks_triplets       indicates that a triplet excited states calculation using
!>                              spin-unpolarised molecular orbitals has been requested
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param matrix_ks             Kohn-Sham matrix
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param kernel_env_admm_aux   kernel environment for ADMM correction
!> \param sub_env               parallel (sub)group environment
!> \param work_matrices         collection of work matrices (modified on exit)
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_Aop_evects(Aop_evects, evects, S_evects, gs_mos, is_rks_triplets, &
                                        do_hfx, matrix_ks, qs_env, kernel_env, kernel_env_admm_aux, &
                                        sub_env, work_matrices)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects, S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_rks_triplets, do_hfx
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_ks
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kernel_env_type), INTENT(in)           :: kernel_env, kernel_env_admm_aux
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_Aop_evects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      LOGICAL                                            :: do_admm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ia_ao, rho_ia_ao_aux_fit
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g_aux_fit, rho_ia_r, &
                                                            rho_ia_r_aux_fit, tau_ia_r, &
                                                            tau_ia_r_aux_fit

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      do_admm = ASSOCIATED(sub_env%admm_A)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0)
         CPASSERT(SIZE(Aop_evects, 1) == nspins)
         CPASSERT(SIZE(Aop_evects, 2) == nvects)
         CPASSERT(SIZE(S_evects, 1) == nspins)
         CPASSERT(SIZE(S_evects, 2) == nvects)
         CPASSERT(SIZE(gs_mos) == nspins)
      END IF

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      IF (nvects > 0) THEN
         CALL cp_fm_get_info(evects(1, 1)%matrix, para_env=para_env)
         CALL qs_rho_get(work_matrices%rho_orb_struct_sub, rho_ao=rho_ia_ao, rho_g=rho_ia_g, rho_r=rho_ia_r, tau_r=tau_ia_r)
         IF (do_hfx .AND. do_admm) THEN
            CALL qs_rho_get(work_matrices%rho_aux_fit_struct_sub, rho_ao=rho_ia_ao_aux_fit, rho_g=rho_ia_g_aux_fit, &
                            rho_r=rho_ia_r_aux_fit, tau_r=tau_ia_r_aux_fit)
         END IF

         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            DO ivect = 1, nvects
               DO ispin = 1, nspins
                  CALL cp_fm_copy_general(evects(ispin, ivect)%matrix, work_matrices%evects_sub(ispin, ivect)%matrix, para_env)
               END DO
            END DO
         END IF

         DO ivect = 1, nvects
            IF (ALLOCATED(work_matrices%evects_sub)) THEN
               IF (ASSOCIATED(work_matrices%evects_sub(1, ivect)%matrix)) THEN
                  DO ispin = 1, nspins
                     CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), &
                                  0.5_dp, sub_env%mos_occ(ispin)%matrix, &
                                  work_matrices%evects_sub(ispin, ivect)%matrix, &
                                  0.0_dp, work_matrices%rho_ao_orb_fm_sub)
                     CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), &
                                  0.5_dp, work_matrices%evects_sub(ispin, ivect)%matrix, &
                                  sub_env%mos_occ(ispin)%matrix, &
                                  1.0_dp, work_matrices%rho_ao_orb_fm_sub)

                     CALL copy_fm_to_dbcsr(work_matrices%rho_ao_orb_fm_sub, &
                                           rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
                  END DO
               ELSE
                  ! skip trial vectors which are assigned to different parallel groups
                  CYCLE
               END IF
            ELSE
               DO ispin = 1, nspins
                  CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, sub_env%mos_occ(ispin)%matrix, &
                               evects(ispin, ivect)%matrix, 0.0_dp, work_matrices%rho_ao_orb_fm_sub)
                  CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, evects(ispin, ivect)%matrix, &
                               sub_env%mos_occ(ispin)%matrix, 1.0_dp, work_matrices%rho_ao_orb_fm_sub)

                  CALL copy_fm_to_dbcsr(work_matrices%rho_ao_orb_fm_sub, &
                                        rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
               END DO
            END IF

            CALL qs_rho_update_rho(work_matrices%rho_orb_struct_sub, qs_env, &
                                   pw_env_external=sub_env%pw_env, task_list_external=sub_env%task_list_orb)

            DO ispin = 1, nspins
               CALL dbcsr_set(work_matrices%A_ia_munu_sub(ispin)%matrix, 0.0_dp)
            END DO

            ! electron-hole exchange-correlation interaction
            DO ispin = 1, nspins
               CALL pw_zero(work_matrices%A_ia_rspace_sub(ispin)%pw)
            END DO

            ! C_x d^{2}E_{x}^{DFT}[\rho] / d\rho^2
            ! + C_{HF} d^{2}E_{x, ADMM}^{DFT}[\rho] / d\rho^2 in case of ADMM calculation
            CALL tddfpt_apply_xc(A_ia_rspace=work_matrices%A_ia_rspace_sub, kernel_env=kernel_env, &
                                 rho_ia_struct=work_matrices%rho_orb_struct_sub, &
                                 is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                                 work_v_xc=work_matrices%wpw_rspace_sub)

            ! ADMM correction
            IF (do_admm) THEN
               CALL tddfpt_construct_aux_fit_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                     rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
                                                     qs_env=qs_env, sub_env=sub_env, &
                                                     wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
                                                     wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
                                                     wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

               ! - C_{HF} d^{2}E_{x, ADMM}^{DFT}[\hat{\rho}] / d\hat{\rho}^2
               CALL tddfpt_apply_xc(A_ia_rspace=work_matrices%A_ia_rspace_sub, &
                                    kernel_env=kernel_env_admm_aux, &
                                    rho_ia_struct=work_matrices%rho_aux_fit_struct_sub, &
                                    is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                                    work_v_xc=work_matrices%wpw_rspace_sub)
            END IF

            ! electron-hole Coulomb interaction
            IF (.NOT. is_rks_triplets) THEN
               ! a sum J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu can be computed by solving
               ! the Poisson equation for combined density (rho_{ia,alpha} + rho_{ia,beta}) .
               ! The following action will destroy reciprocal-space grid in spin-unrestricted case.
               DO ispin = 2, nspins
                  CALL pw_axpy(rho_ia_g(ispin)%pw, rho_ia_g(1)%pw)
               END DO

               CALL tddfpt_apply_coulomb(A_ia_rspace=work_matrices%A_ia_rspace_sub, &
                                         rho_ia_g=rho_ia_g(1)%pw, pw_env=sub_env%pw_env, &
                                         work_v_gspace=work_matrices%wpw_gspace_sub(1), &
                                         work_v_rspace=work_matrices%wpw_rspace_sub(1))
            END IF

            ! convert from the plane-wave representation into the Gaussian basis set representation
            DO ispin = 1, nspins
               CALL pw_scale(work_matrices%A_ia_rspace_sub(ispin)%pw, work_matrices%A_ia_rspace_sub(ispin)%pw%pw_grid%dvol)

               CALL integrate_v_rspace(v_rspace=work_matrices%A_ia_rspace_sub(ispin), &
                                       hmat=work_matrices%A_ia_munu_sub(ispin), &
                                       qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE., &
                                       pw_env_external=sub_env%pw_env, task_list_external=sub_env%task_list_orb)
            END DO

            IF (ALLOCATED(work_matrices%evects_sub)) THEN
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(work_matrices%A_ia_munu_sub(ispin)%matrix, sub_env%mos_occ(ispin)%matrix, &
                                               work_matrices%Aop_evects_sub(ispin, ivect)%matrix, &
                                               ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               END DO
            ELSE
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(work_matrices%A_ia_munu_sub(ispin)%matrix, sub_env%mos_occ(ispin)%matrix, &
                                               Aop_evects(ispin, ivect)%matrix, ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               END DO
            END IF
         END DO

         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            DO ivect = 1, nvects
               DO ispin = 1, nspins
                  CALL cp_fm_copy_general(work_matrices%Aop_evects_sub(ispin, ivect)%matrix, &
                                          Aop_evects(ispin, ivect)%matrix, para_env)
               END DO
            END DO
         END IF

         ! orbital energy difference term
         CALL tddfpt_apply_energy_diff(Aop_evects=Aop_evects, evects=evects, S_evects=S_evects, gs_mos=gs_mos, &
                                       matrix_ks=matrix_ks, work_fm_ao_mo_occ=work_matrices%wfm_ao_mo_occ)

         IF (do_hfx) THEN
            CALL tddfpt_apply_hfx(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, do_admm=do_admm, &
                                  qs_env=qs_env, work_rho_ia_ao=work_matrices%hfx_rho_ao, &
                                  work_hmat=work_matrices%hfx_hmat, wfm_rho_orb=work_matrices%hfx_fm_ao_ao)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_Aop_evects

! **************************************************************************************************
!> \brief Solve eigenproblem for the reduced action matrix and find new Ritz eigenvectors and
!>        eigenvalues.
!> \param ritz_vects       Ritz eigenvectors (initialised on exit)
!> \param Aop_ritz         approximate action of TDDFPT operator on Ritz vectors (initialised on exit)
!> \param evals            Ritz eigenvalues (initialised on exit)
!> \param krylov_vects     Krylov's vectors
!> \param Aop_krylov       action of TDDFPT operator on Krylov's vectors
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 altered prototype, OpenMP parallelisation [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_ritz_vects(ritz_vects, Aop_ritz, evals, krylov_vects, Aop_krylov)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: krylov_vects, Aop_krylov

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_ritz_vects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikv, irv, ispin, nkvs, nrvs, &
                                                            nspins
      REAL(kind=dp)                                      :: act
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: Atilde, evects_Atilde
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: Atilde_fm, evects_Atilde_fm

      CALL timeset(routineN, handle)

      nspins = SIZE(krylov_vects, 1)
      nkvs = SIZE(krylov_vects, 2)
      nrvs = SIZE(ritz_vects, 2)

      CALL cp_fm_get_info(krylov_vects(1, 1)%matrix, context=blacs_env_global)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nkvs, ncol_global=nkvs, context=blacs_env_global)
      CALL cp_fm_create(Atilde_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! *** compute upper-diagonal reduced action matrix ***
      ALLOCATE (Atilde(nkvs, nkvs))
      ! TO DO: the subroutine 'cp_fm_contracted_trace' will compute all elements of
      ! the matrix 'Atilde', however only upper-triangular elements are actually needed
      CALL cp_fm_contracted_trace(Aop_krylov, krylov_vects, Atilde)
      CALL cp_fm_set_submatrix(Atilde_fm, Atilde)
      DEALLOCATE (Atilde)

      ! *** solve an eigenproblem for the reduced matrix ***
      CALL choose_eigv_solver(Atilde_fm, evects_Atilde_fm, evals(1:nkvs))

      ALLOCATE (evects_Atilde(nkvs, nrvs))
      CALL cp_fm_get_submatrix(evects_Atilde_fm, evects_Atilde, start_row=1, start_col=1, n_rows=nkvs, n_cols=nrvs)
      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(Atilde_fm)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(act, ikv, irv, ispin), &
!$OMP             SHARED(Aop_krylov, Aop_ritz, krylov_vects, evects_Atilde, nkvs, nrvs, nspins, ritz_vects)
      DO irv = 1, nrvs
         DO ispin = 1, nspins
            CALL cp_fm_set_all(ritz_vects(ispin, irv)%matrix, 0.0_dp)
            CALL cp_fm_set_all(Aop_ritz(ispin, irv)%matrix, 0.0_dp)
         END DO

         DO ikv = 1, nkvs
            act = evects_Atilde(ikv, irv)
            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, ritz_vects(ispin, irv)%matrix, &
                                        act, krylov_vects(ispin, ikv)%matrix)
               CALL cp_fm_scale_and_add(1.0_dp, Aop_ritz(ispin, irv)%matrix, &
                                        act, Aop_krylov(ispin, ikv)%matrix)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      DEALLOCATE (evects_Atilde)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_ritz_vects

! **************************************************************************************************
!> \brief Expand Krylov space by computing residual vectors.
!> \param residual_vects          residual vectors (modified on exit)
!> \param evals                   Ritz eigenvalues (modified on exit)
!> \param ritz_vects              Ritz eigenvectors
!> \param Aop_ritz                approximate action of TDDFPT operator on Ritz vectors
!> \param gs_mos                  molecular orbitals optimised for the ground state
!> \param matrix_s                overlap matrix
!> \param work_fm_ao_mo_occ       work dense matrices with shape [nao x nmo_occ(spin)]
!>                                (modified on exit)
!> \param work_fm_mo_virt_mo_occ  work dense matrices with shape [nmo_virt(spin) x nmo_occ(spin)]
!>                                (modified on exit)
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored to achieve significant performance gain [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_residual_vects(residual_vects, evals, ritz_vects, Aop_ritz, gs_mos, &
                                            matrix_s, work_fm_ao_mo_occ, work_fm_mo_virt_mo_occ)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: residual_vects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: work_fm_ao_mo_occ, work_fm_mo_virt_mo_occ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_residual_vects', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, icol_local, irow_local, irv, &
                                                            ispin, nao, ncols_local, nrows_local, &
                                                            nrvs, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_local, row_indices_local
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: e_occ_plus_lambda, eref, lambda
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: weights_ldata

      CALL timeset(routineN, handle)

      nspins = SIZE(residual_vects, 1)
      nrvs = SIZE(residual_vects, 2)

      CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      IF (nrvs > 0) THEN
         ! *** actually compute residual vectors ***
         DO irv = 1, nrvs
            lambda = evals(irv)

            DO ispin = 1, nspins
               CALL cp_fm_get_info(work_fm_mo_virt_mo_occ(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local, &
                                   row_indices=row_indices_local, col_indices=col_indices_local, local_data=weights_ldata)

               ! work_fm_ao_mo_occ := Ab(ispin, irv) - evals(irv) b(ispin, irv), where 'b' is a Ritz vector
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, ritz_vects(ispin, irv)%matrix, work_fm_ao_mo_occ(ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=-lambda, beta=0.0_dp)
               CALL cp_fm_scale_and_add(1.0_dp, work_fm_ao_mo_occ(ispin)%matrix, 1.0_dp, Aop_ritz(ispin, irv)%matrix)

               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, gs_mos(ispin)%mos_virt, &
                            work_fm_ao_mo_occ(ispin)%matrix, 0.0_dp, work_fm_mo_virt_mo_occ(ispin)%matrix)

               DO icol_local = 1, ncols_local
                  e_occ_plus_lambda = gs_mos(ispin)%evals_occ(col_indices_local(icol_local)) + lambda

                  DO irow_local = 1, nrows_local
                     eref = gs_mos(ispin)%evals_virt(row_indices_local(irow_local)) - e_occ_plus_lambda

                     ! eref = e_virt - e_occ - lambda = e_virt - e_occ - (eref_scale*lambda + (1-eref_scale)*lambda);
                     ! eref_new = e_virt - e_occ - eref_scale*lambda = eref + (1 - eref_scale)*lambda
                     IF (ABS(eref) < threshold) &
                        eref = eref + (1.0_dp - eref_scale)*lambda

                     weights_ldata(irow_local, icol_local) = weights_ldata(irow_local, icol_local)/eref
                  END DO
               END DO

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nmo_virt(ispin), 1.0_dp, gs_mos(ispin)%mos_virt, &
                            work_fm_mo_virt_mo_occ(ispin)%matrix, 0.0_dp, residual_vects(ispin, irv)%matrix)
            END DO
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_residual_vects

! **************************************************************************************************
!> \brief Perform Davidson iterations.
!> \param evects                TDDFPT trial vectors (modified on exit)
!> \param evals                 TDDFPT eigenvalues (modified on exit)
!> \param S_evects              cached matrix product S * evects (modified on exit)
!> \param gs_mos                molecular orbitals optimised for the ground state
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param tddfpt_control        TDDFPT control parameters
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param kernel_env_admm_aux   kernel environment for ADMM correction
!> \param sub_env               parallel (sub)group environment
!> \param logger                CP2K logger
!> \param iter_unit             I/O unit to write basic iteration information
!> \param energy_unit           I/O unit to write detailed energy information
!> \param tddfpt_print_section  TDDFPT print input section (need to write TDDFPT restart files)
!> \param work_matrices         collection of work matrices (modified on exit)
!> \return energy convergence achieved (in Hartree)
!> \par History
!>    * 03.2017 code related to Davidson eigensolver has been moved here from the main subroutine
!>              tddfpt() [Sergey Chulkov]
!> \note Based on the subroutines apply_op() and iterative_solver() originally created by
!>       Thomas Chassaing in 2002.
! **************************************************************************************************
   FUNCTION tddfpt_davidson_solver(evects, evals, S_evects, gs_mos, do_hfx, tddfpt_control, &
                                   qs_env, kernel_env, kernel_env_admm_aux, &
                                   sub_env, logger, iter_unit, energy_unit, &
                                   tddfpt_print_section, work_matrices) RESULT(conv)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: do_hfx
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kernel_env_type), INTENT(in)           :: kernel_env, kernel_env_admm_aux
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, INTENT(in)                                :: iter_unit, energy_unit
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices
      REAL(kind=dp)                                      :: conv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_davidson_solver', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate, iter, &
                                                            max_krylov_vects, nspins, nstates, &
                                                            nstates_conv, nvects_exist, nvects_new
      INTEGER(kind=int_8)                                :: nstates_total
      LOGICAL                                            :: is_nonortho
      REAL(kind=dp)                                      :: t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_last
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Aop_krylov, Aop_ritz, krylov_vects, &
                                                            S_krylov
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)
      nstates = tddfpt_control%nstates
      nstates_total = tddfpt_total_number_of_states(gs_mos)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(evects, 1) == nspins)
         CPASSERT(SIZE(evects, 2) == nstates)
         CPASSERT(SIZE(evals) == nstates)
      END IF

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s)

      ! adjust the number of Krylov vectors
      max_krylov_vects = tddfpt_control%nkvs
      IF (max_krylov_vects < nstates) max_krylov_vects = nstates
      IF (INT(max_krylov_vects, kind=int_8) > nstates_total) max_krylov_vects = INT(nstates_total)

      ALLOCATE (Aop_ritz(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (Aop_ritz(ispin, istate)%matrix)
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_ritz(ispin, istate)%matrix)
         END DO
      END DO

      ALLOCATE (evals_last(max_krylov_vects))
      ALLOCATE (Aop_krylov(nspins, max_krylov_vects), krylov_vects(nspins, max_krylov_vects), S_krylov(nspins, max_krylov_vects))
      DO istate = 1, max_krylov_vects
         DO ispin = 1, nspins
            NULLIFY (Aop_krylov(ispin, istate)%matrix, krylov_vects(ispin, istate)%matrix, S_krylov(ispin, istate)%matrix)
         END DO
      END DO

      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, krylov_vects(ispin, istate)%matrix)
            CALL cp_fm_to_fm(evects(ispin, istate)%matrix, krylov_vects(ispin, istate)%matrix)

            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_krylov(ispin, istate)%matrix)
            CALL cp_fm_to_fm(S_evects(ispin, istate)%matrix, S_krylov(ispin, istate)%matrix)

            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_krylov(ispin, istate)%matrix)
         END DO
      END DO

      nvects_exist = 0
      nvects_new = nstates

      t1 = m_walltime()

      DO
         ! davidson iteration
         CALL cp_iterate(logger%iter_info, iter_nr_out=iter)

         CALL tddfpt_compute_Aop_evects(Aop_evects=Aop_krylov(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                        evects=krylov_vects(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                        S_evects=S_krylov(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                        gs_mos=gs_mos, is_rks_triplets=tddfpt_control%rks_triplets, &
                                        do_hfx=do_hfx, matrix_ks=matrix_ks, &
                                        qs_env=qs_env, kernel_env=kernel_env, &
                                        kernel_env_admm_aux=kernel_env_admm_aux, &
                                        sub_env=sub_env, &
                                        work_matrices=work_matrices)

         CALL tddfpt_compute_ritz_vects(ritz_vects=evects, Aop_ritz=Aop_ritz, &
                                        evals=evals_last(1:nvects_exist + nvects_new), &
                                        krylov_vects=krylov_vects(:, 1:nvects_exist + nvects_new), &
                                        Aop_krylov=Aop_krylov(:, 1:nvects_exist + nvects_new))

         CALL tddfpt_write_restart(evects=evects, evals=evals_last(1:nstates), gs_mos=gs_mos, &
                                   logger=logger, tddfpt_print_section=tddfpt_print_section)

         conv = MAXVAL(ABS(evals_last(1:nstates) - evals(1:nstates)))

         nvects_exist = nvects_exist + nvects_new
         IF (nvects_exist + nvects_new > max_krylov_vects) &
            nvects_new = max_krylov_vects - nvects_exist
         IF (iter >= tddfpt_control%niters) nvects_new = 0

         IF (conv > tddfpt_control%conv .AND. nvects_new > 0) THEN
            ! compute residual vectors for the next iteration
            DO istate = 1, nvects_new
               DO ispin = 1, nspins
                  NULLIFY (Aop_krylov(ispin, nvects_exist + istate)%matrix, &
                           krylov_vects(ispin, nvects_exist + istate)%matrix, &
                           S_krylov(ispin, nvects_exist + istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         krylov_vects(ispin, nvects_exist + istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         S_krylov(ispin, nvects_exist + istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         Aop_krylov(ispin, nvects_exist + istate)%matrix)
               END DO
            END DO

            CALL tddfpt_compute_residual_vects(residual_vects=krylov_vects(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                               evals=evals_last(1:nvects_new), &
                                               ritz_vects=evects(:, 1:nvects_new), Aop_ritz=Aop_ritz(:, 1:nvects_new), &
                                               gs_mos=gs_mos, matrix_s=matrix_s(1)%matrix, &
                                               work_fm_ao_mo_occ=work_matrices%wfm_ao_mo_occ, &
                                               work_fm_mo_virt_mo_occ=work_matrices%wfm_mo_virt_mo_occ)

            CALL tddfpt_orthogonalize_psi1_psi0(krylov_vects(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                                work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)

            CALL tddfpt_orthonormalize_psi1_psi1(krylov_vects(:, 1:nvects_exist + nvects_new), nvects_new, &
                                                 S_krylov(:, 1:nvects_exist + nvects_new), matrix_s(1)%matrix)

            is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(krylov_vects(:, nvects_exist + 1:nvects_exist + nvects_new), &
                                                            work_matrices%S_C0, tddfpt_control%orthogonal_eps, &
                                                            work_matrices%wfm_mo_occ_mo_occ)
         ELSE
            ! convergence or the maximum number of Krylov vectors have been achieved
            nvects_new = 0
            is_nonortho = .FALSE.
         END IF

         t2 = m_walltime()
         IF (energy_unit > 0) THEN
            WRITE (energy_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
            DO istate = 1, nstates
               WRITE (energy_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
                  evals_last(istate)*evolt, (evals_last(istate) - evals(istate))*evolt
            END DO
            WRITE (energy_unit, *)
            CALL m_flush(energy_unit)
         END IF

         IF (iter_unit > 0) THEN
            nstates_conv = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate) - evals(istate)) <= tddfpt_control%conv) &
                  nstates_conv = nstates_conv + 1
            END DO

            WRITE (iter_unit, '(1X,I8,T12,F7.1,T24,ES11.4,T42,I8)') iter, t2 - t1, conv, nstates_conv
            CALL m_flush(iter_unit)
         END IF

         t1 = t2
         evals(1:nstates) = evals_last(1:nstates)

         ! nvects_new == 0 if iter >= tddfpt_control%niters
         IF (nvects_new == 0 .OR. is_nonortho) THEN
            ! restart Davidson iterations
            CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)
            CALL tddfpt_orthonormalize_psi1_psi1(evects, nstates, S_evects, matrix_s(1)%matrix)

            EXIT
         END IF
      END DO

      DO istate = nvects_exist + nvects_new, 1, -1
         DO ispin = nspins, 1, -1
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_krylov(ispin, istate)%matrix)
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_krylov(ispin, istate)%matrix)
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, krylov_vects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_krylov, krylov_vects, S_krylov)
      DEALLOCATE (evals_last)

      DO istate = nstates, 1, -1
         DO ispin = nspins, 1, -1
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_ritz(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_ritz)

      CALL timestop(handle)
   END FUNCTION tddfpt_davidson_solver

! **************************************************************************************************
!> \brief Compute the action of the dipole operator on the ground state wave function.
!> \param dipole_op_mos_occ  2-D array [x,y,z ; spin] of matrices where to put the computed quantity
!>                           (allocated and initialised on exit)
!> \param tddfpt_control     TDDFPT control parameters
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param qs_env             Quickstep environment
!> \par History
!>    * 05.2016 created as 'tddfpt_print_summary' [Sergey Chulkov]
!>    * 06.2018 dipole operator based on the Berry-phase formula [Sergey Chulkov]
!>    * 08.2018 splited of from 'tddfpt_print_summary' and merged with code from 'tddfpt'
!>              [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       The relation between dipole integrals in velocity and length forms are the following:
!>       \f[<\psi_i|\nabla|\psi_a> = <\psi_i|\vec{r}|\hat{H}\psi_a> - <\hat{H}\psi_i|\vec{r}|\psi_a>
!>                                 = (\epsilon_a - \epsilon_i) <\psi_i|\vec{r}|\psi_a> .\f],
!>       due to the commutation identity:
!>       \f[\vec{r}\hat{H} - \hat{H}\vec{r} = [\vec{r},\hat{H}] = [\vec{r},-1/2 \nabla^2] = \nabla\f] .
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control, gs_mos, qs_env)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(inout)                                   :: dipole_op_mos_occ
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_dipole_operator', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cos_sin, icol, ideriv, irow, &
                                                            ispin, jderiv, nao, ncols_local, &
                                                            ndim_periodic, nrows_local, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: eval_occ
      REAL(kind=dp), DIMENSION(3)                        :: kvec, reference_point
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data_ediff, local_data_wfm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: gamma_00, gamma_inv_00
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_mos_virt
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dBerry_mos_occ, gamma_real_imag, opvec
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ediff_inv, rRc_mos_occ, wfm_ao_ao, &
                                                            wfm_mo_virt_mo_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: berry_cossin_xyz, matrix_s, rRc_xyz
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, cell, matrix_s, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, matrix_s=matrix_s, pw_env=pw_env)

      nspins = SIZE(gs_mos)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      ! +++ allocate dipole operator matrices (must be deallocated elsewhere)
      ALLOCATE (dipole_op_mos_occ(nderivs, nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)

         DO ideriv = 1, nderivs
            NULLIFY (dipole_op_mos_occ(ideriv, ispin)%matrix)
            CALL cp_fm_create(dipole_op_mos_occ(ideriv, ispin)%matrix, fm_struct)
         END DO
      END DO

      ! +++ allocate work matrices
      ALLOCATE (S_mos_virt(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_mos_virt(ispin)%matrix)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fm_struct)
         CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, &
                                      gs_mos(ispin)%mos_virt, &
                                      S_mos_virt(ispin)%matrix, &
                                      ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! check that the chosen dipole operator is consistent with the periodic boundary conditions used
      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      ndim_periodic = COUNT(poisson_env%parameters%periodic == 1)

      SELECT CASE (tddfpt_control%dipole_form)
      CASE (tddfpt_dipole_berry)
         IF (ndim_periodic /= 3) THEN
            CALL cp_warn(__LOCATION__, &
                         "Fully periodic Poisson solver (PERIODIC xyz) is needed "// &
                         "for oscillator strengths based on the Berry phase formula")
         END IF

         NULLIFY (berry_cossin_xyz)
         ! index: 1 = Re[exp(-i * G_t * t)],
         !        2 = Im[exp(-i * G_t * t)];
         ! t = x,y,z
         CALL dbcsr_allocate_matrix_set(berry_cossin_xyz, 2)

         DO i_cos_sin = 1, 2
            CALL dbcsr_init_p(berry_cossin_xyz(i_cos_sin)%matrix)
            CALL dbcsr_copy(berry_cossin_xyz(i_cos_sin)%matrix, matrix_s(1)%matrix)
         END DO

         ! +++ allocate berry-phase-related work matrices
         ALLOCATE (gamma_00(nspins), gamma_inv_00(nspins), gamma_real_imag(2, nspins), opvec(2, nspins))
         ALLOCATE (dBerry_mos_occ(nderivs, nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)

            NULLIFY (gamma_00(ispin)%matrix, gamma_inv_00(ispin)%matrix)
            CALL cp_cfm_create(gamma_00(ispin)%matrix, fm_struct)
            CALL cp_cfm_create(gamma_inv_00(ispin)%matrix, fm_struct)

            DO i_cos_sin = 1, 2
               NULLIFY (gamma_real_imag(i_cos_sin, ispin)%matrix)
               CALL cp_fm_create(gamma_real_imag(i_cos_sin, ispin)%matrix, fm_struct)
            END DO
            CALL cp_fm_struct_release(fm_struct)

            ! G_real C_0, G_imag C_0
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)
            DO i_cos_sin = 1, 2
               NULLIFY (opvec(i_cos_sin, ispin)%matrix)
               CALL cp_fm_create(opvec(i_cos_sin, ispin)%matrix, fm_struct)
            END DO

            ! dBerry * C_0
            DO ideriv = 1, nderivs
               NULLIFY (dBerry_mos_occ(ideriv, ispin)%matrix)
               CALL cp_fm_create(dBerry_mos_occ(ideriv, ispin)%matrix, fm_struct)
               CALL cp_fm_set_all(dBerry_mos_occ(ideriv, ispin)%matrix, 0.0_dp)
            END DO
         END DO

         DO ideriv = 1, nderivs
            kvec(:) = twopi*cell%h_inv(ideriv, :)
            DO i_cos_sin = 1, 2
               CALL dbcsr_set(berry_cossin_xyz(i_cos_sin)%matrix, 0.0_dp)
            END DO
            CALL build_berry_moment_matrix(qs_env, berry_cossin_xyz(1)%matrix, berry_cossin_xyz(2)%matrix, kvec)

            DO ispin = 1, nspins
               ! i_cos_sin = 1: cos (real) component; opvec(1) = gamma_real C_0
               ! i_cos_sin = 2: sin (imaginary) component; opvec(2) = gamma_imag C_0
               DO i_cos_sin = 1, 2
                  CALL cp_dbcsr_sm_fm_multiply(berry_cossin_xyz(i_cos_sin)%matrix, &
                                               gs_mos(ispin)%mos_occ, &
                                               opvec(i_cos_sin, ispin)%matrix, &
                                               ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               END DO

               CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, &
                            1.0_dp, gs_mos(ispin)%mos_occ, opvec(1, ispin)%matrix, &
                            0.0_dp, gamma_real_imag(1, ispin)%matrix)

               CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, &
                            -1.0_dp, gs_mos(ispin)%mos_occ, opvec(2, ispin)%matrix, &
                            0.0_dp, gamma_real_imag(2, ispin)%matrix)

               CALL cp_fm_to_cfm(msourcer=gamma_real_imag(1, ispin)%matrix, &
                                 msourcei=gamma_real_imag(2, ispin)%matrix, &
                                 mtarget=gamma_00(ispin)%matrix)

               ! gamma_inv_00 = Q = [C_0^T (gamma_real - i gamma_imag) C_0] ^ {-1}
               CALL cp_cfm_set_all(gamma_inv_00(ispin)%matrix, z_zero, z_one)
               CALL cp_cfm_solve(gamma_00(ispin)%matrix, gamma_inv_00(ispin)%matrix)

               CALL cp_cfm_to_fm(msource=gamma_inv_00(ispin)%matrix, &
                                 mtargetr=gamma_real_imag(1, ispin)%matrix, &
                                 mtargeti=gamma_real_imag(2, ispin)%matrix)

               ! dBerry_mos_occ is identical to dBerry_psi0 from qs_linres_op % polar_operators()
               CALL cp_gemm("N", "N", nao, nmo_occ(ispin), nmo_occ(ispin), &
                            1.0_dp, opvec(1, ispin)%matrix, gamma_real_imag(2, ispin)%matrix, &
                            0.0_dp, dipole_op_mos_occ(1, ispin)%matrix)
               CALL cp_gemm("N", "N", nao, nmo_occ(ispin), nmo_occ(ispin), &
                            -1.0_dp, opvec(2, ispin)%matrix, gamma_real_imag(1, ispin)%matrix, &
                            1.0_dp, dipole_op_mos_occ(1, ispin)%matrix)

               DO jderiv = 1, nderivs
                  CALL cp_fm_scale_and_add(1.0_dp, dBerry_mos_occ(jderiv, ispin)%matrix, &
                                           cell%hmat(jderiv, ideriv), dipole_op_mos_occ(1, ispin)%matrix)
               END DO
            END DO
         END DO

         ! --- release berry-phase-related work matrices
         DO ispin = nspins, 1, -1
            DO i_cos_sin = SIZE(opvec, 1), 1, -1
               CALL cp_fm_release(opvec(i_cos_sin, ispin)%matrix)
            END DO

            DO i_cos_sin = SIZE(gamma_real_imag, 1), 1, -1
               CALL cp_fm_release(gamma_real_imag(i_cos_sin, ispin)%matrix)
            END DO

            CALL cp_cfm_release(gamma_inv_00(ispin)%matrix)
            CALL cp_cfm_release(gamma_00(ispin)%matrix)
         END DO
         DEALLOCATE (gamma_00, gamma_inv_00, gamma_real_imag, opvec)
         CALL dbcsr_deallocate_matrix_set(berry_cossin_xyz)

         NULLIFY (fm_struct, wfm_ao_ao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(wfm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         ! trans_dipole = 2|e|/|G_mu| * Tr Imag(evects^T * (gamma_real - i gamma_imag) * C_0 * gamma_inv_00) +
         !                2|e|/|G_mu| * Tr Imag(C_0^T * (gamma_real - i gamma_imag) * evects * gamma_inv_00) ,
         !
         ! Taking into account the symmetry of the matrices 'gamma_real' and 'gamma_imag' and the fact
         ! that the response wave-function is a real-valued function, the above expression can be simplified as
         ! trans_dipole = 4|e|/|G_mu| * Tr Imag(evects^T * (gamma_real - i gamma_imag) * C_0 * gamma_inv_00)
         !
         ! 1/|G_mu| = |lattice_vector_mu| / (2*pi) .
         DO ispin = 1, nspins
            ! wfm_ao_ao = S * mos_virt * mos_virt^T
            CALL cp_gemm('N', 'T', nao, nao, nmo_virt(ispin), &
                         1.0_dp/twopi, S_mos_virt(ispin)%matrix, gs_mos(ispin)%mos_virt, &
                         0.0_dp, wfm_ao_ao)

            DO ideriv = 1, nderivs
               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, &
                            1.0_dp, wfm_ao_ao, dBerry_mos_occ(ideriv, ispin)%matrix, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO
         END DO

         CALL cp_fm_release(wfm_ao_ao)

         DO ispin = nspins, 1, -1
            DO ideriv = SIZE(dBerry_mos_occ, 1), 1, -1
               CALL cp_fm_release(dBerry_mos_occ(ideriv, ispin)%matrix)
            END DO
         END DO
         DEALLOCATE (dBerry_mos_occ)

      CASE (tddfpt_dipole_length)
         IF (ndim_periodic /= 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Non-periodic Poisson solver (PERIODIC none) is needed "// &
                         "for oscillator strengths based on the length operator")
         END IF

         ! compute components of the dipole operator in the length form
         NULLIFY (rRc_xyz)
         CALL dbcsr_allocate_matrix_set(rRc_xyz, nderivs)

         DO ideriv = 1, nderivs
            CALL dbcsr_init_p(rRc_xyz(ideriv)%matrix)
            CALL dbcsr_copy(rRc_xyz(ideriv)%matrix, matrix_s(1)%matrix)
         END DO

         CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

         CALL rRc_xyz_ao(op=rRc_xyz, qs_env=qs_env, rc=reference_point, order=1, minimum_image=.FALSE., soft=.FALSE.)

         NULLIFY (fm_struct, wfm_ao_ao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(wfm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         DO ispin = 1, nspins
            NULLIFY (rRc_mos_occ)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)
            CALL cp_fm_create(rRc_mos_occ, fm_struct)

            ! wfm_ao_ao = S * mos_virt * mos_virt^T
            CALL cp_gemm('N', 'T', nao, nao, nmo_virt(ispin), &
                         1.0_dp, S_mos_virt(ispin)%matrix, gs_mos(ispin)%mos_virt, &
                         0.0_dp, wfm_ao_ao)

            DO ideriv = 1, nderivs
               CALL cp_dbcsr_sm_fm_multiply(rRc_xyz(ideriv)%matrix, &
                                            gs_mos(ispin)%mos_occ, &
                                            rRc_mos_occ, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, &
                            1.0_dp, wfm_ao_ao, rRc_mos_occ, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO

            CALL cp_fm_release(rRc_mos_occ)
         END DO

         CALL cp_fm_release(wfm_ao_ao)
         CALL dbcsr_deallocate_matrix_set(rRc_xyz)

      CASE (tddfpt_dipole_velocity)
         IF (SIZE(matrix_s) < nderivs + 1) THEN
            CPABORT("Where is the derivative?")
         END IF

         DO ispin = 1, nspins
            NULLIFY (fm_struct, ediff_inv, wfm_mo_virt_mo_occ)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
            CALL cp_fm_create(ediff_inv, fm_struct)
            CALL cp_fm_create(wfm_mo_virt_mo_occ, fm_struct)
            CALL cp_fm_struct_release(fm_struct)

            CALL cp_fm_get_info(ediff_inv, nrow_local=nrows_local, ncol_local=ncols_local, &
                                row_indices=row_indices, col_indices=col_indices, local_data=local_data_ediff)
            CALL cp_fm_get_info(wfm_mo_virt_mo_occ, local_data=local_data_wfm)

!$OMP       PARALLEL DO DEFAULT(NONE), &
!$OMP                PRIVATE(eval_occ, icol, irow), &
!$OMP                SHARED(col_indices, gs_mos, ispin, local_data_ediff, ncols_local, nrows_local, row_indices)
            DO icol = 1, ncols_local
               ! E_occ_i ; imo_occ = col_indices(icol)
               eval_occ = gs_mos(ispin)%evals_occ(col_indices(icol))

               DO irow = 1, nrows_local
                  ! ediff_inv_weights(a, i) = 1.0 / (E_virt_a - E_occ_i)
                  ! imo_virt = row_indices(irow)
                  local_data_ediff(irow, icol) = 1.0_dp/(gs_mos(ispin)%evals_virt(row_indices(irow)) - eval_occ)
               END DO
            END DO
!$OMP       END PARALLEL DO

            DO ideriv = 1, nderivs
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(ideriv + 1)%matrix, &
                                            gs_mos(ispin)%mos_occ, &
                                            dipole_op_mos_occ(ideriv, ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, &
                            1.0_dp, gs_mos(ispin)%mos_virt, dipole_op_mos_occ(ideriv, ispin)%matrix, &
                            0.0_dp, wfm_mo_virt_mo_occ)

               ! in-place element-wise (Schur) product;
               ! avoid allocation of a temporary [nmo_virt x nmo_occ] matrix which is needed for cp_fm_schur_product() subroutine call

!$OMP          PARALLEL DO DEFAULT(NONE), &
!$OMP                   PRIVATE(icol, irow), &
!$OMP                   SHARED(ispin, local_data_ediff, local_data_wfm, ncols_local, nrows_local)
               DO icol = 1, ncols_local
                  DO irow = 1, nrows_local
                     local_data_wfm(irow, icol) = local_data_wfm(irow, icol)*local_data_ediff(irow, icol)
                  END DO
               END DO
!$OMP          END PARALLEL DO

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nmo_virt(ispin), &
                            1.0_dp, S_mos_virt(ispin)%matrix, wfm_mo_virt_mo_occ, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO

            CALL cp_fm_release(wfm_mo_virt_mo_occ)
            CALL cp_fm_release(ediff_inv)
         END DO

      CASE DEFAULT
         CPABORT("Unimplemented form of the dipole operator")
      END SELECT

      ! --- release work matrices
      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_dipole_operator

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and oscillator strengths.
!> \param log_unit           output unit
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param evals              TDDFPT eigenvalues
!> \param mult               multiplicity
!> \param dipole_op_mos_occ  action of the dipole operator on the ground state wave function
!>                           [x,y,z ; spin]
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!>    * 07.2016 spin-unpolarised electron density [Sergey Chulkov]
!>    * 08.2018 compute 'dipole_op_mos_occ' in a separate subroutine [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       Transition dipole moment along direction 'd' is computed as following:
!>       \f[ t_d(spin) = Tr[evects^T dipole\_op\_mos\_occ(d, spin)] .\f]
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_print_summary(log_unit, evects, evals, mult, dipole_op_mos_occ)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: mult
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: dipole_op_mos_occ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_summary', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=1)                                   :: lsd_str
      CHARACTER(len=20)                                  :: mult_str
      INTEGER                                            :: handle, ideriv, ispin, istate, nspins, &
                                                            nstates
      REAL(kind=dp)                                      :: osc_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)

      IF (nspins > 1) THEN
         lsd_str = 'U'
      ELSE
         lsd_str = 'R'
      END IF

      ! *** summary header ***
      IF (log_unit > 0) THEN
         CALL integer_to_string(mult, mult_str)
         WRITE (log_unit, '(/,1X,A1,A,1X,A,/)') lsd_str, "-TDDFPT states of multiplicity", TRIM(mult_str)

         WRITE (log_unit, '(T10,A,T19,A,T37,A,T69,A)') "State", "Excitation", "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T10,A,T19,A,T37,A,T49,A,T61,A,T67,A)') "number", "energy (eV)", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T10,72("-"))')
      END IF

      ! transition dipole moment
      ALLOCATE (trans_dipoles(nstates, nderivs, nspins))
      trans_dipoles(:, :, :) = 0.0_dp

      ! nspins == 1 .AND. mult == 3 : spin-flip transitions are forbidden due to symmetry reasons
      IF (nspins > 1 .OR. mult == 1) THEN
         DO ispin = 1, nspins
            DO ideriv = 1, nderivs
               CALL cp_fm_trace(evects(ispin, :), dipole_op_mos_occ(ideriv, ispin)%matrix, trans_dipoles(:, ideriv, ispin))
            END DO
         END DO

         IF (nspins == 1) THEN
            trans_dipoles(:, :, 1) = SQRT(2.0_dp)*trans_dipoles(:, :, 1)
         ELSE
            trans_dipoles(:, :, 1) = SQRT(trans_dipoles(:, :, 1)**2 + trans_dipoles(:, :, 2)**2)
         END IF
      END IF

      ! *** summary information ***
      DO istate = 1, nstates
         IF (log_unit > 0) THEN
            osc_strength = 2.0_dp/3.0_dp*evals(istate)* &
                           accurate_dot_product(trans_dipoles(istate, :, 1), trans_dipoles(istate, :, 1))

            WRITE (log_unit, '(1X,A,T9,I7,T19,F11.5,T31,3(1X,ES11.4E2),T69,ES12.5E2)') &
               "TDDFPT|", istate, evals(istate)*evolt, trans_dipoles(istate, 1:nderivs, 1), osc_strength
         END IF
      END DO

      DEALLOCATE (trans_dipoles)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary

! **************************************************************************************************
!> \brief Print excitation analysis.
!> \param log_unit           output unit
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param matrix_s           overlap matrix
!> \param min_amplitude      the smallest excitation amplitude to print
!> \par History
!>    * 05.2016 created as 'tddfpt_print_summary' [Sergey Chulkov]
!>    * 08.2018 splited of from 'tddfpt_print_summary' [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_print_excitation_analysis(log_unit, evects, gs_mos, matrix_s, min_amplitude)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_excitation_analysis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=5)                                   :: spin_label
      INTEGER :: handle, icol, iproc, irow, ispin, istate, nao, ncols_local, nrows_local, nspins, &
         nstates, send_handler, send_handler2, state_spin
      INTEGER(kind=int_8)                                :: iexc, imo_occ, imo_virt, ind, nexcs, &
                                                            nexcs_local, nexcs_max_local, &
                                                            nmo_virt_occ, nmo_virt_occ_alpha
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: inds_local, inds_recv, nexcs_recv
      INTEGER(kind=int_8), DIMENSION(1)                  :: nexcs_send
      INTEGER(kind=int_8), DIMENSION(maxspins)           :: nmo_occ8, nmo_virt8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds, recv_handlers, recv_handlers2
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: do_exc_analysis
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_local, weights_neg_abs_recv, &
                                                            weights_recv
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_mos_virt, weights_fm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      do_exc_analysis = min_amplitude < 1.0_dp

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, context=blacs_env, para_env=para_env)
      CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
         nmo_occ8(ispin) = SIZE(gs_mos(ispin)%evals_occ, kind=int_8)
         nmo_virt8(ispin) = SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO

      ! *** excitation analysis ***
      IF (do_exc_analysis) THEN
         CPASSERT(log_unit <= 0 .OR. para_env%ionode)
         nmo_virt_occ_alpha = INT(nmo_virt(1), int_8)*INT(nmo_occ(1), int_8)

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1X,A,/)') "Excitation analysis"

            WRITE (log_unit, '(3X,A,T17,A,T34,A,T49,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(3X,A,T18,A,T34,A,T49,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,57("-"))')

            IF (nspins == 1) THEN
               state_spin = 1
               spin_label = '     '
            END IF
         END IF

         ALLOCATE (S_mos_virt(nspins), weights_fm(nspins))
         DO ispin = 1, nspins
            NULLIFY (S_mos_virt(ispin)%matrix)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fm_struct)
            CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, &
                                         gs_mos(ispin)%mos_virt, &
                                         S_mos_virt(ispin)%matrix, &
                                         ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)

            NULLIFY (fm_struct, weights_fm(ispin)%matrix)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
            CALL cp_fm_create(weights_fm(ispin)%matrix, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         nexcs_max_local = 0
         DO ispin = 1, nspins
            CALL cp_fm_get_info(weights_fm(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local)
            nexcs_max_local = nexcs_max_local + INT(nrows_local, int_8)*INT(ncols_local, int_8)
         END DO

         ALLOCATE (weights_local(nexcs_max_local), inds_local(nexcs_max_local))

         DO istate = 1, nstates
            nexcs_local = 0
            nmo_virt_occ = 0

            ! analyse matrix elements locally and transfer only significant
            ! excitations to the master node for subsequent ordering
            DO ispin = 1, nspins
               ! compute excitation amplitudes
               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin)%matrix, &
                            evects(ispin, istate)%matrix, 0.0_dp, weights_fm(ispin)%matrix)

               CALL cp_fm_get_info(weights_fm(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local, &
                                   row_indices=row_indices, col_indices=col_indices, local_data=local_data)

               ! locate single excitations with significant amplitudes (>= min_amplitude)
               DO icol = 1, ncols_local
                  DO irow = 1, nrows_local
                     IF (ABS(local_data(irow, icol)) >= min_amplitude) THEN
                        ! number of non-negligible excitations
                        nexcs_local = nexcs_local + 1
                        ! excitation amplitude
                        weights_local(nexcs_local) = local_data(irow, icol)
                        ! index of single excitation (ivirt, iocc, ispin) in compressed form
                        inds_local(nexcs_local) = nmo_virt_occ + INT(row_indices(irow), int_8) + &
                                                  INT(col_indices(icol) - 1, int_8)*nmo_virt8(ispin)
                     END IF
                  END DO
               END DO

               nmo_virt_occ = nmo_virt_occ + nmo_virt8(ispin)*nmo_occ8(ispin)
            END DO

            IF (para_env%ionode) THEN
               ! master node
               ALLOCATE (nexcs_recv(para_env%num_pe), recv_handlers(para_env%num_pe), recv_handlers2(para_env%num_pe))

               ! collect number of non-negligible excitations from other nodes
               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) THEN
                     CALL mp_irecv(nexcs_recv(iproc:iproc), iproc - 1, para_env%group, recv_handlers(iproc), 0)
                  ELSE
                     nexcs_recv(iproc) = nexcs_local
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) &
                     CALL mp_wait(recv_handlers(iproc))
               END DO

               ! compute total number of non-negligible excitations
               nexcs = 0
               DO iproc = 1, para_env%num_pe
                  nexcs = nexcs + nexcs_recv(iproc)
               END DO

               ! receive indices and amplitudes of selected excitations
               ALLOCATE (weights_recv(nexcs), weights_neg_abs_recv(nexcs))
               ALLOCATE (inds_recv(nexcs), inds(nexcs))

               nmo_virt_occ = 0
               DO iproc = 1, para_env%num_pe
                  IF (nexcs_recv(iproc) > 0) THEN
                     IF (iproc - 1 /= para_env%mepos) THEN
                        ! excitation amplitudes
                        CALL mp_irecv(weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                      iproc - 1, para_env%group, recv_handlers(iproc), 1)
                        ! compressed indices
                        CALL mp_irecv(inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                      iproc - 1, para_env%group, recv_handlers2(iproc), 2)
                     ELSE
                        ! data on master node
                        weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = weights_local(1:nexcs_recv(iproc))
                        inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = inds_local(1:nexcs_recv(iproc))
                     END IF

                     nmo_virt_occ = nmo_virt_occ + nexcs_recv(iproc)
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos .AND. nexcs_recv(iproc) > 0) THEN
                     CALL mp_wait(recv_handlers(iproc))
                     CALL mp_wait(recv_handlers2(iproc))
                  END IF
               END DO

               DEALLOCATE (nexcs_recv, recv_handlers, recv_handlers2)
            ELSE
               ! working node: send the number of selected excited states to the master node
               nexcs_send(1) = nexcs_local
               CALL mp_isend(nexcs_send, para_env%source, para_env%group, send_handler, 0)
               CALL mp_wait(send_handler)

               IF (nexcs_local > 0) THEN
                  ! send excitation amplitudes
                  CALL mp_isend(weights_local(1:nexcs_local), para_env%source, para_env%group, send_handler, 1)
                  ! send compressed indices
                  CALL mp_isend(inds_local(1:nexcs_local), para_env%source, para_env%group, send_handler2, 2)

                  CALL mp_wait(send_handler)
                  CALL mp_wait(send_handler2)
               END IF
            END IF

            ! sort non-negligible excitations on the master node according to their amplitudes,
            ! uncompress indices and print summary information
            IF (para_env%ionode .AND. log_unit > 0) THEN
               weights_neg_abs_recv(:) = -ABS(weights_recv)
               CALL sort(weights_neg_abs_recv, INT(nexcs), inds)

               WRITE (log_unit, '(1X,I8)') istate

               DO iexc = 1, nexcs
                  ind = inds_recv(inds(iexc)) - 1
                  IF (nspins > 1) THEN
                     IF (ind < nmo_virt_occ_alpha) THEN
                        state_spin = 1
                        spin_label = '(alp)'
                     ELSE
                        state_spin = 2
                        ind = ind - nmo_virt_occ_alpha
                        spin_label = '(bet)'
                     END IF
                  END IF

                  imo_occ = ind/nmo_virt8(state_spin) + 1
                  imo_virt = MOD(ind, nmo_virt8(state_spin)) + 1

                  WRITE (log_unit, '(T14,I8,1X,A5,T30,I8,1X,A5,T50,F9.6)') imo_occ, spin_label, &
                     nmo_occ8(state_spin) + imo_virt, spin_label, weights_recv(inds(iexc))
               END DO
            END IF

            ! deallocate temporary arrays
            IF (para_env%ionode) &
               DEALLOCATE (weights_recv, weights_neg_abs_recv, inds_recv, inds)
         END DO

         DEALLOCATE (weights_local, inds_local)
      END IF

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(weights_fm(ispin)%matrix)
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt, weights_fm)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_excitation_analysis

! **************************************************************************************************
!> \brief Write Ritz vectors to a binary restart file.
!> \param evects               vectors to store
!> \param evals                TDDFPT eigenvalues
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_write_restart(evects, evals, gs_mos, logger, tddfpt_print_section)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_write_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate, nao, nspins, &
                                                            nstates, ounit
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ

      IF (BTEST(cp_print_key_should_output(logger%iter_info, tddfpt_print_section, "RESTART"), cp_p_file)) THEN
         CALL timeset(routineN, handle)

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)

         IF (debug_this_module) THEN
            CPASSERT(SIZE(evals) == nstates)
            CPASSERT(nspins > 0)
            CPASSERT(nstates > 0)
         END IF

         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)

         DO ispin = 1, nspins
            nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         END DO

         ounit = cp_print_key_unit_nr(logger, tddfpt_print_section, "RESTART", &
                                      extension=".tdwfn", file_status="REPLACE", file_action="WRITE", &
                                      do_backup=.TRUE., file_form="UNFORMATTED")

         IF (ounit > 0) THEN
            WRITE (ounit) nstates, nspins, nao
            WRITE (ounit) nmo_occ(1:nspins)
            WRITE (ounit) evals
         END IF

         DO istate = 1, nstates
            DO ispin = 1, nspins
               ! TDDFPT wave function is actually stored as a linear combination of virtual MOs
               ! that replaces the corresponding deoccupied MO. Unfortunately, the phase
               ! of the occupied MOs varies depending on the eigensolver used as well as
               ! how eigenvectors are distributed across computational cores. The phase is important
               ! because TDDFPT wave functions are used to compute a response electron density
               ! \rho^{-} = 1/2 * [C_{0} * evect^T + evect * C_{0}^{-}], where C_{0} is the expansion
               ! coefficients of the reference ground-state wave function. To make the restart file
               ! transferable, TDDFPT wave functions are stored in assumption that all ground state
               ! MOs have a positive phase.
               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)

               CALL cp_fm_write_unformatted(evects(ispin, istate)%matrix, ounit)

               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)
            END DO
         END DO

         CALL cp_print_key_finished_output(ounit, logger, tddfpt_print_section, "RESTART")

         CALL timestop(handle)
      END IF
   END SUBROUTINE tddfpt_write_restart

! **************************************************************************************************
!> \brief Initialise initial guess vectors by reading (un-normalised) Ritz vectors
!>        from a binary restart file.
!> \param evects               vectors to initialise (initialised on exit)
!> \param evals                TDDFPT eigenvalues (initialised on exit)
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_section       TDDFPT input section
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \param fm_pool_ao_mo_occ    pools of dense matrices with shape [nao x nmo_occ(spin)]
!> \param blacs_env_global     BLACS parallel environment involving all the processor
!> \return the number of excited states found in the restart file
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_read_restart(evects, evals, gs_mos, logger, tddfpt_section, tddfpt_print_section, &
                                fm_pool_ao_mo_occ, blacs_env_global) RESULT(nstates_read)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_section, tddfpt_print_section
      TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(in)  :: fm_pool_ao_mo_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      INTEGER                                            :: nstates_read

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_read_restart', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=20)                                  :: read_str, ref_str
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, ispin, istate, iunit, n_rep_val, &
                                                            nao, nao_read, nspins, nspins_read, &
                                                            nstates
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_occ_read
      LOGICAL                                            :: file_exists
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_read
      TYPE(cp_para_env_type), POINTER                    :: para_env_global
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      ! generate restart file name
      CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", c_val=filename)
      ELSE
         print_key => section_vals_get_subs_vals(tddfpt_print_section, "RESTART")
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".tdwfn", my_local=.FALSE.)
      END IF

      CALL get_blacs_info(blacs_env_global, para_env=para_env_global)

      IF (para_env_global%ionode) THEN
         INQUIRE (FILE=filename, exist=file_exists)

         IF (.NOT. file_exists) THEN
            nstates_read = 0
            CALL mp_bcast(nstates_read, para_env_global%source, para_env_global%group)

            CALL cp_warn(__LOCATION__, &
                         "User requested to restart the TDDFPT wave functions from the file '"//TRIM(filename)// &
                         "' which does not exist. Guess wave functions will be constructed using Kohn-Sham orbitals.")
            CALL timestop(handle)
            RETURN
         END IF

         CALL open_file(file_name=filename, file_action="READ", file_form="UNFORMATTED", file_status="OLD", unit_number=iunit)
      END IF

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      IF (para_env_global%ionode) THEN
         READ (iunit) nstates_read, nspins_read, nao_read

         IF (nspins_read /= nspins) THEN
            CALL integer_to_string(nspins, ref_str)
            CALL integer_to_string(nspins_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Restarted TDDFPT wave function contains incompatible number of spin components ("// &
                          TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         IF (nao_read /= nao) THEN
            CALL integer_to_string(nao, ref_str)
            CALL integer_to_string(nao_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Incompatible number of atomic orbitals ("//TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         READ (iunit) nmo_occ_read(1:nspins)

         DO ispin = 1, nspins
            IF (nmo_occ_read(ispin) /= nmo_occ(ispin)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Incompatible number of electrons and/or multiplicity.")
            END IF
         END DO

         IF (nstates_read /= nstates) THEN
            CALL integer_to_string(nstates, ref_str)
            CALL integer_to_string(nstates_read, read_str)
            CALL cp_warn(__LOCATION__, &
                         "TDDFPT restart file contains "//TRIM(read_str)// &
                         " wave function(s) however "//TRIM(ref_str)// &
                         " excited states were requested.")
         END IF
      END IF
      CALL mp_bcast(nstates_read, para_env_global%source, para_env_global%group)

      ! exit if restart file does not exist
      IF (nstates_read <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (para_env_global%ionode) THEN
         ALLOCATE (evals_read(nstates_read))
         READ (iunit) evals_read
         IF (nstates_read <= nstates) THEN
            evals(1:nstates_read) = evals_read(1:nstates_read)
         ELSE
            evals(1:nstates) = evals_read(1:nstates)
         END IF
         DEALLOCATE (evals_read)
      END IF
      CALL mp_bcast(evals, para_env_global%source, para_env_global%group)

      DO istate = 1, nstates_read
         DO ispin = 1, nspins
            IF (istate <= nstates) THEN
               CALL fm_pool_create_fm(fm_pool_ao_mo_occ(ispin)%pool, evects(ispin, istate)%matrix)

               CALL cp_fm_read_unformatted(evects(ispin, istate)%matrix, iunit)

               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)
            END IF
         END DO
      END DO

      IF (para_env_global%ionode) &
         CALL close_file(unit_number=iunit)

      CALL timestop(handle)
   END FUNCTION tddfpt_read_restart
END MODULE qs_tddfpt2_methods
