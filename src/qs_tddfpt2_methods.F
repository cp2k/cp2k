!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_methods
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_deallocate_matrix, &
        cp_dbcsr_deallocate_matrix_set, cp_dbcsr_desymmetrize, cp_dbcsr_distribution_release, &
        cp_dbcsr_get_info, cp_dbcsr_init, cp_dbcsr_init_p, cp_dbcsr_p_type, cp_dbcsr_release, &
        cp_dbcsr_set, cp_dbcsr_type, dbcsr_distribution_obj, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_p_trace,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_schur_product,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_maxabsval, &
        cp_fm_p_type, cp_fm_read_unformatted, cp_fm_release, cp_fm_set_all, cp_fm_to_fm, &
        cp_fm_type, cp_fm_write_unformatted
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
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_2d
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE input_constants,                 ONLY: cholesky_dbcsr,&
                                              cholesky_inverse,&
                                              cholesky_off,&
                                              cholesky_restore,&
                                              xas_dip_len,&
                                              xas_dip_vel
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_split,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_max,&
                                              mp_min,&
                                              mp_sum,&
                                              mp_wait
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types_new,              ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale,&
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
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_operators_ao,                 ONLY: rRc_xyz_ao
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE string_utilities,                ONLY: integer_to_string
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
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

   PUBLIC :: tddfpt
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3

! **************************************************************************************************
!> \brief Ground state molecular orbitals
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_ground_state_mos
      !> occupied MOs stored in a matrix form [nao x nmo_occ]
      TYPE(cp_fm_type), POINTER :: matrix_mos_occ
      !> virtual MOs stored in a matrix form [nao x nmo_virt]
      TYPE(cp_fm_type), POINTER :: matrix_mos_virt
      !> virtual MOs stored as individual vectors [nao x 1]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:) :: vectors_mo_virt
      !> occupied orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_occ
      !> virtual orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals_virt
      !> phase of occupied MOs; +1.0 -- positive, -1.0 -- negative;
      !> it is mainly needed to make the restart file transferable
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: phases_occ
   END TYPE tddfpt_ground_state_mos

! **************************************************************************************************
!> \brief Time-Dependent DFT guess vectors
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_guess_vectors
      !> indices of occupied MOs to deoccupy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imos_occ
      !> indices of virtual MOs to occupy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imos_virt
      !> spin of excited electron (in spin-polarized case)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ispins
      !> transition energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals
   END TYPE tddfpt_guess_vectors

! **************************************************************************************************
!> \brief Allocatable list of full matrices
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_fm_vector
      !> allocatable list of full matrices
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:) :: vect
   END TYPE tddfpt_fm_vector

! **************************************************************************************************
!> \brief Allocatable 2-D matrix
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_allocatable_matrix
      !> allocatable 2-D matrix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :) :: matrix
   END TYPE tddfpt_allocatable_matrix

! **************************************************************************************************
!> \brief Allocatable array of integer numbers
! **************************************************************************************************
   TYPE tddfpt_int_vector
      !> allocatable array of integer numbers
      INTEGER, ALLOCATABLE, DIMENSION(:)          :: vect
   END TYPE tddfpt_int_vector

! **************************************************************************************************
!> \brief Allocatable array of real numbers
! **************************************************************************************************
   TYPE tddfpt_real_vector
      !> allocatable array of real numbers
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)    :: vect
   END TYPE tddfpt_real_vector

! **************************************************************************************************
!> \brief structure to save global multi-grid related parameters
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE mgrid_saved_parameters
      !> create commensurate grids
      LOGICAL                                     :: commensurate_mgrids
      !> create real-space grids
      LOGICAL                                     :: realspace_mgrids
      !> do not perform load balancing
      LOGICAL                                     :: skip_load_balance
      !> cutoff value at the finest grid level
      REAL(KIND=dp)                               :: cutoff
      !> inverse scale factor
      REAL(KIND=dp)                               :: progression_factor
      !> relative cutoff
      REAL(KIND=dp)                               :: relative_cutoff
      !> list of explicitly given cutoff values
      REAL(KIND=dp), DIMENSION(:), POINTER        :: e_cutoff
   END TYPE mgrid_saved_parameters
CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT calculation.
!> \param qs_env  Quickstep environment
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 refactored to be used with Davidson eigensolver [Sergey Chulkov]
!> \note Based on the subroutines apply_op(), iterative_solver(), tddfpt_env_init(), and
!>       tddfpt_env_deallocate().
! **************************************************************************************************
   SUBROUTINE tddfpt(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: max_reortho = 16384

      CHARACTER(len=20)                                  :: navail_str, nstates_str
      INTEGER :: energy_unit, handle, ispin, istate, isubgroup, iter, log_unit, max_krylov_vects, &
         mpi_comm_sub, mult, nao, nao_aux, niters, nspins, nstates, nstates_conv, nstates_total, &
         nsubgroups, nvects_exists, nvects_new
      INTEGER(kind=int_8)                                :: nstates_total_i8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: subgroup
      LOGICAL                                            :: do_hfx, is_nonortho, is_restarted
      REAL(kind=dp)                                      :: C_hf, conv, t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_last, evals_prev
      REAL(kind=dp), DIMENSION(3)                        :: reference_point
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global, blacs_env_sub
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: dipole_op_xyz_sub, matrix_ks_global, &
         matrix_ks_sub, matrix_s_aux_fit_global, matrix_s_global, rho_ij_ao, rRc_xyz
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s_aux_fit_sub, matrix_s_sub
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0_C0T_sub, S_C0_sub
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Aop_evects, Aop_ritz, evects, ritz_vects
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: fm_struct_evects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_aux_fit
      TYPE(cp_fm_type), POINTER                          :: work_fm_aux_fit_global, &
                                                            work_fm_aux_fit_sub, work_fm_global, &
                                                            work_fm_sub
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env_global, para_env_sub
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist_sub
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d_sub
      TYPE(mgrid_saved_parameters)                       :: mgrid_saved
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit_sub, sab_orb_sub
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r, rho_ij_r2, tau_ij_r, tau_ij_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(qs_rho_type), POINTER                         :: rho_ij_struct
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, &
                                                            tddfpt_print_section, tddfpt_section, &
                                                            xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit_sub, task_list_orb_sub
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), ALLOCATABLE, &
         DIMENSION(:)                                    :: gs_mos
      TYPE(tddfpt_guess_vectors)                         :: guess_vectors
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(Iannuzzi2005)

      NULLIFY (blacs_env_global, dft_control, input, matrix_ks_global, matrix_s_global, mos, para_env_global)
      CALL get_qs_env(qs_env, blacs_env=blacs_env_global, dft_control=dft_control, &
                      input=input, matrix_ks=matrix_ks_global, matrix_s=matrix_s_global, &
                      mos=mos, para_env=para_env_global)
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s_aux_fit_global)
      END IF
      tddfpt_control => dft_control%tddfpt2_control
      qs_control => dft_control%qs_control

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
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")

      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=C_hf)
         do_hfx = (C_hf /= 0.0_dp)
      END IF

      ! compute components of the dipole operator in the length form
      IF (tddfpt_control%dipole_form == xas_dip_len) THEN
         NULLIFY (rRc_xyz)
         CALL cp_dbcsr_allocate_matrix_set(rRc_xyz, nderivs)

         DO ispin = 1, nderivs
            CALL cp_dbcsr_init_p(rRc_xyz(ispin)%matrix)
            CALL cp_dbcsr_copy(rRc_xyz(ispin)%matrix, matrix_s_global(1)%matrix)
         END DO

         ! get the reference point
         CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

         CALL rRc_xyz_ao(op=rRc_xyz, qs_env=qs_env, rc=reference_point, order=1, minimum_image=.FALSE., soft=.FALSE.)
      END IF

      ! ++ split mpi communicator

      ! suppress the mp_comm_split() exception when the minimal subgroup size
      ! is larger then the number of available processors
      IF (tddfpt_control%nprocs > para_env_global%num_pe) &
         tddfpt_control%nprocs = para_env_global%num_pe

      ! proc_membership(iproc) stores the index of the mpi subgroup to which the given process belongs
      ALLOCATE (subgroup(0:para_env_global%num_pe-1))

      ! default: assign all available processors to a single subgroup
      IF (tddfpt_control%nprocs == 0) &
         tddfpt_control%nprocs = para_env_global%num_pe

      CALL mp_comm_split(comm=para_env_global%group, sub_comm=mpi_comm_sub, ngroups=nsubgroups, &
                         group_distribution=subgroup, subgroup_min_size=tddfpt_control%nprocs)
      isubgroup = subgroup(para_env_global%mepos)

      NULLIFY (blacs_env_sub, para_env_sub)
      CALL cp_para_env_create(para_env_sub, mpi_comm_sub)
      ! use the default (SQUARE) BLACS grid layout and non-repeatable BLACS collective operations
      ! by omitting optional parameters 'blacs_grid_layout' and 'blacs_repeatable'.
      ! Ideally we should take these parameters from the variables globenv%blacs_grid_layout and
      ! globenv%blacs_repeatable, however the global environment is not available
      ! from the subroutine 'qs_energies_properties'.
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! initialise plain wave environment
      log_unit = cp_logger_get_default_io_unit()

      IF (tddfpt_control%mgrid_is_explicit) &
         CALL init_tddfpt_mgrid(qs_control, tddfpt_control, mgrid_saved)

      NULLIFY (pw_env_sub)
      CALL pw_env_create(pw_env_sub)
      CALL pw_env_rebuild(pw_env_sub, qs_env, para_env_sub)

      CALL tddfpt_build_distribution_2d(distribution_2d=dist_2d_sub, dbcsr_dist=dbcsr_dist_sub, &
                                        blacs_env=blacs_env_sub, qs_env=qs_env, log_unit=log_unit)
      CALL tddfpt_build_tasklist(task_list=task_list_orb_sub, sab=sab_orb_sub, basis_type="ORB", &
                                 distribution_2d=dist_2d_sub, pw_env=pw_env_sub, qs_env=qs_env, &
                                 skip_load_balance=qs_control%skip_load_balance_distributed, reorder_grid_ranks=.TRUE.)

      IF (dft_control%do_admm) THEN
         CALL tddfpt_build_tasklist(task_list=task_list_aux_fit_sub, sab=sab_aux_fit_sub, basis_type="AUX_FIT", &
                                    distribution_2d=dist_2d_sub, pw_env=pw_env_sub, qs_env=qs_env, &
                                    skip_load_balance=qs_control%skip_load_balance_distributed, reorder_grid_ranks=.FALSE.)
      ELSE
         NULLIFY (sab_aux_fit_sub, task_list_aux_fit_sub)
      END IF

      ! recompute ground-state Kohn-Sham orbitals
      ALLOCATE (gs_mos(nspins))
      DO ispin = 1, nspins
         CALL tddfpt_init_ground_state_mos(gs_mos=gs_mos(ispin), mo_set=mos(ispin)%mo_set, &
                                           blacs_env_global=blacs_env_global, blacs_env_sub=blacs_env_sub, &
                                           cholesky_method=cholesky_restore, &
                                           matrix_ks_global=matrix_ks_global(ispin)%matrix, &
                                           matrix_s_global=matrix_s_global(1)%matrix, &
                                           nsubgroups=nsubgroups, subgroup=subgroup)
      END DO

      ! check we have enough molecular orbitals to construct the requested number of excited states
      nstates_total_i8 = 0
      DO ispin = 1, nspins
         nstates_total_i8 = nstates_total_i8+SIZE(gs_mos(ispin)%evals_occ, kind=int_8)* &
                            SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO

      IF (INT(tddfpt_control%nstates, kind=int_8) > nstates_total_i8) THEN
         CALL integer_to_string(tddfpt_control%nstates, nstates_str)
         CALL integer_to_string(INT(nstates_total_i8), navail_str)
         CALL cp_warn(__LOCATION__, TRIM(nstates_str)//" excited states were requested, however only "// &
                      TRIM(navail_str)//" excited states can be constructed using the atomic basis set given.")
         tddfpt_control%nstates = INT(nstates_total_i8)
      END IF

      ! adjust the number of trial vectors
      IF (tddfpt_control%added_states < 0 .OR. &
          INT(tddfpt_control%added_states, kind=int_8)+INT(tddfpt_control%nstates, kind=int_8) > nstates_total_i8) THEN

         tddfpt_control%added_states = INT(nstates_total_i8-INT(tddfpt_control%nstates, kind=int_8))
      END IF
      nstates = tddfpt_control%nstates

      ! reset rks_triplets if UKS is in use
      IF (tddfpt_control%rks_triplets .AND. nspins > 1) THEN
         tddfpt_control%rks_triplets = .FALSE.
         CALL cp_warn(__LOCATION__, "Keyword RKS_TRIPLETS has been ignored for spin-polarized calculations")
      END IF

      ! multiplicity of molecular system
      IF (nspins > 1) THEN
         mult = ABS(SIZE(gs_mos(1)%evals_occ)-SIZE(gs_mos(2)%evals_occ))+1
         IF (mult > 2) &
            CALL cp_warn(__LOCATION__, "There is a convergence issue for multiplicity >= 3")
      ELSE
         IF (tddfpt_control%rks_triplets) THEN
            mult = 3
         ELSE
            mult = 1
         END IF
      END IF

      ! build the list of singly excited states and sort them in ascending order according to excitation energies
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "GUESS_VECTORS", extension=".tddfptLog")
      CALL tddfpt_init_guess_vectors(guess_vectors=guess_vectors, nstates=nstates, &
                                     nstates_active=nstates+tddfpt_control%added_states, &
                                     gs_mos=gs_mos, is_add_degenerate=tddfpt_control%add_degenerate, &
                                     degenerate_eps=tddfpt_control%degenerate_eps, log_unit=log_unit)
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "GUESS_VECTORS")

      nstates_total = SIZE(guess_vectors%evals)

      ! adjust the number of Krylov vectors
      max_krylov_vects = tddfpt_control%nkvs
      IF (max_krylov_vects < nstates) max_krylov_vects = nstates
      IF (max_krylov_vects > nstates_total) max_krylov_vects = nstates_total

      ! matrix structures for eigevectors
      ALLOCATE (fm_struct_evects(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, matrix_struct=fm_struct_evects(ispin)%struct)
      END DO

      ! replicate Kohn-Sham and overlap matrices across subgroups
      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao)

      NULLIFY (fm_struct, work_fm_global, work_fm_sub)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_global)
      CALL cp_fm_create(work_fm_global, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_sub)
      CALL cp_fm_create(work_fm_sub, fm_struct)

      NULLIFY (matrix_ks_sub)
      CALL cp_dbcsr_allocate_matrix_set(matrix_ks_sub, nspins)
      DO ispin = 1, nspins
         CALL tddfpt_dbcsr_create_by_dist(matrix_ks_sub(ispin)%matrix, &
                                          template=matrix_ks_global(ispin)%matrix, &
                                          dbcsr_dist=dbcsr_dist_sub, sab=sab_orb_sub)

         CALL copy_dbcsr_to_fm(matrix_ks_global(ispin)%matrix, work_fm_global)
         CALL tddfpt_fm_replicate_across_subgroups(fm_src=work_fm_global, fm_dest_sub=work_fm_sub, &
                                                   nsubgroups=nsubgroups, subgroup=subgroup)
         CALL copy_fm_to_dbcsr(work_fm_sub, matrix_ks_sub(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      ! overlap matrix
      NULLIFY (matrix_s_sub)
      CALL cp_dbcsr_init_p(matrix_s_sub)

      CALL tddfpt_dbcsr_create_by_dist(matrix_s_sub, &
                                       template=matrix_s_global(1)%matrix, &
                                       dbcsr_dist=dbcsr_dist_sub, &
                                       sab=sab_orb_sub)
      CALL copy_dbcsr_to_fm(matrix_s_global(1)%matrix, work_fm_global)
      CALL tddfpt_fm_replicate_across_subgroups(fm_src=work_fm_global, fm_dest_sub=work_fm_sub, &
                                                nsubgroups=nsubgroups, subgroup=subgroup)
      CALL copy_fm_to_dbcsr(work_fm_sub, matrix_s_sub, keep_sparsity=.TRUE.)

      ! components of the dipole operator
      NULLIFY (dipole_op_xyz_sub)
      CALL cp_dbcsr_allocate_matrix_set(dipole_op_xyz_sub, nderivs)
      DO istate = 1, nderivs
         SELECT CASE (tddfpt_control%dipole_form)
         CASE (xas_dip_len)
            CALL tddfpt_dbcsr_create_by_dist(dipole_op_xyz_sub(istate)%matrix, &
                                             template=rRc_xyz(istate)%matrix, &
                                             dbcsr_dist=dbcsr_dist_sub, &
                                             sab=sab_orb_sub)
            CALL copy_dbcsr_to_fm(rRc_xyz(istate)%matrix, work_fm_global)

         CASE (xas_dip_vel)
            CALL tddfpt_dbcsr_create_by_dist(dipole_op_xyz_sub(istate)%matrix, &
                                             template=matrix_s_global(istate+1)%matrix, &
                                             dbcsr_dist=dbcsr_dist_sub, &
                                             sab=sab_orb_sub)
            CALL copy_dbcsr_to_fm(matrix_s_global(istate+1)%matrix, work_fm_global)

         CASE DEFAULT
            CPABORT("Unimplemented form of the dipole operator")
         END SELECT

         CALL tddfpt_fm_replicate_across_subgroups(fm_src=work_fm_global, fm_dest_sub=work_fm_sub, &
                                                   nsubgroups=nsubgroups, subgroup=subgroup)
         CALL copy_fm_to_dbcsr(work_fm_sub, dipole_op_xyz_sub(istate)%matrix, keep_sparsity=.TRUE.)
      END DO

      CALL cp_fm_release(work_fm_global)
      IF (tddfpt_control%dipole_form == xas_dip_len) &
         CALL cp_dbcsr_deallocate_matrix_set(rRc_xyz)

      NULLIFY (matrix_s_aux_fit_sub)
      IF (dft_control%do_admm) THEN
         CALL cp_dbcsr_get_info(matrix_s_aux_fit_global(1)%matrix, nfullrows_total=nao_aux)
         CALL tddfpt_dbcsr_create_by_dist(matrix_s_aux_fit_sub, &
                                          template=matrix_s_aux_fit_global(1)%matrix, &
                                          dbcsr_dist=dbcsr_dist_sub, &
                                          sab=sab_aux_fit_sub)

         NULLIFY (fm_struct_aux_fit, work_fm_aux_fit_global, work_fm_aux_fit_sub)
         CALL cp_fm_struct_create(fm_struct_aux_fit, nrow_global=nao_aux, ncol_global=nao_aux, context=blacs_env_global)
         CALL cp_fm_create(work_fm_aux_fit_global, fm_struct_aux_fit)
         CALL cp_fm_struct_release(fm_struct_aux_fit)

         CALL cp_fm_struct_create(fm_struct_aux_fit, nrow_global=nao_aux, ncol_global=nao_aux, context=blacs_env_sub)
         CALL cp_fm_create(work_fm_aux_fit_sub, fm_struct_aux_fit)
         CALL cp_fm_struct_release(fm_struct_aux_fit)

         CALL copy_dbcsr_to_fm(matrix_s_aux_fit_global(1)%matrix, work_fm_aux_fit_global)
         CALL tddfpt_fm_replicate_across_subgroups(fm_src=work_fm_aux_fit_global, fm_dest_sub=work_fm_aux_fit_sub, &
                                                   nsubgroups=nsubgroups, subgroup=subgroup)
         CALL copy_fm_to_dbcsr(work_fm_aux_fit_sub, matrix_s_aux_fit_sub, keep_sparsity=.TRUE.)

         CALL cp_fm_release(work_fm_aux_fit_sub)
         CALL cp_fm_release(work_fm_aux_fit_global)
      END IF

      ! *** initialise xc_rho_set and xc_deriv_set (needed for kernel in adiabatic approximation) ***
      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool)

      ! ++ S * C_{0, sigma}
      ALLOCATE (S_C0_sub(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_C0_sub(ispin)%matrix)
         CALL cp_fm_create(S_C0_sub(ispin)%matrix, fm_struct_evects(ispin)%struct)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s_sub, gs_mos(ispin)%matrix_mos_occ, S_C0_sub(ispin)%matrix, &
                                      ncol=SIZE(gs_mos(ispin)%evals_occ), alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! ++ S * P_{0, sigma}
      ALLOCATE (S_C0_C0T_sub(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_C0_C0T_sub(ispin)%matrix)
         CALL cp_fm_create(S_C0_C0T_sub(ispin)%matrix, fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      ! rho_ij_struct
      NULLIFY (rho_ij_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ij_ao, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(rho_ij_ao(ispin)%matrix)
         CALL cp_dbcsr_copy(rho_ij_ao(ispin)%matrix, matrix_s_sub)
      END DO

      NULLIFY (rho_ij_struct)
      CALL qs_rho_create(rho_ij_struct)
      CALL qs_rho_set(rho_ij_struct, rho_ao=rho_ij_ao)
      CALL qs_rho_rebuild(rho_ij_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE., pw_env_external=pw_env_sub)
      CALL qs_rho_get(rho_ij_struct, rho_r=rho_ij_r, tau_r=tau_ij_r)

      DO ispin = 1, nspins
         CALL cp_gemm('N', 'T', nao, nao, SIZE(gs_mos(ispin)%evals_occ), 1.0_dp, &
                      gs_mos(ispin)%matrix_mos_occ, gs_mos(ispin)%matrix_mos_occ, 0.0_dp, work_fm_sub)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s_sub, work_fm_sub, S_C0_C0T_sub(ispin)%matrix, &
                                      ncol=nao, alpha=1.0_dp, beta=0.0_dp)

         ! take into account that all MOs are doubly occupied in spin-restricted case
         IF (nspins == 1 .AND. (.NOT. tddfpt_control%rks_triplets)) &
            CALL cp_fm_scale(2.0_dp, work_fm_sub)

         CALL copy_fm_to_dbcsr(work_fm_sub, rho_ij_ao(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO
      CALL qs_rho_update_rho(rho_ij_struct, qs_env, pw_env_external=pw_env_sub, task_list_external=task_list_orb_sub)

      NULLIFY (xc_rho_set, xc_deriv_set)

      IF (tddfpt_control%rks_triplets) THEN
         ! going to compute triplet states using spin-restricted MOs;
         ! still need spin-beta density in order to compute the TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1)%pw => rho_ij_r(1)%pw
         rho_ij_r2(2)%pw => rho_ij_r(1)%pw
         IF (ASSOCIATED(tau_ij_r)) THEN
            ALLOCATE (tau_ij_r2(2))
            tau_ij_r2(1)%pw => tau_ij_r(1)%pw
            tau_ij_r2(2)%pw => tau_ij_r(1)%pw
         ELSE
            NULLIFY (tau_ij_r2)
         END IF

         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r2, auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r2)

         IF (ASSOCIATED(tau_ij_r2)) &
            DEALLOCATE (tau_ij_r2)

         DEALLOCATE (rho_ij_r2)
      ELSE
         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho_ij_r, auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r)
      END IF

      CALL qs_rho_release(rho_ij_struct)
      CALL cp_fm_release(work_fm_sub)

      ! *** allocate storage for trial vectors and matrix-vector products ***
      ALLOCATE (evects(nspins, max_krylov_vects))
      ALLOCATE (Aop_evects(nspins, max_krylov_vects))
      DO istate = 1, max_krylov_vects
         DO ispin = 1, nspins
            NULLIFY (evects(ispin, istate)%matrix, Aop_evects(ispin, istate)%matrix)
         END DO
      END DO

      ALLOCATE (ritz_vects(nspins, nstates))
      ALLOCATE (Aop_ritz(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (ritz_vects(ispin, istate)%matrix, Aop_ritz(ispin, istate)%matrix)
         END DO
      END DO

      DO istate = nsubgroups-isubgroup, nstates, nsubgroups
         DO ispin = 1, nspins
            ! allocate Ritz vectors;
            CALL cp_fm_create(ritz_vects(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
            CALL cp_fm_create(Aop_ritz(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
         END DO
      END DO

      ! *** initial Krylov vectors ***
      ALLOCATE (evals_last(max_krylov_vects))
      ALLOCATE (evals_prev(nstates))
      evals_prev(:) = 0.0_dp

      DO istate = nsubgroups-isubgroup, nstates, nsubgroups
         ! initial guess : allocate initial guess vectors; other vectors will be allocated
         !                 during the Krylov's space expansion stage
         DO ispin = 1, nspins
            CALL cp_fm_create(evects(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
            CALL cp_fm_set_all(evects(ispin, istate)%matrix, 0.0_dp)
         END DO
      END DO

      ! reuse Ritz vectors from the previous calculation if available
      IF (tddfpt_control%is_restart) THEN
         nstates_conv = tddfpt_read_restart(evects=evects(:, 1:nstates), evals=evals_prev, gs_mos=gs_mos, &
                                            logger=logger, tddfpt_section=tddfpt_section, &
                                            tddfpt_print_section=tddfpt_print_section, blacs_env_global=blacs_env_global)
      ELSE
         nstates_conv = 0
      END IF

      ! restart Davidson iterations at least once if some of the requested
      ! excited states have not been read from the restart file
      is_restarted = nstates_conv >= nstates

      DO istate = nsubgroups-isubgroup, nstates, nsubgroups
         IF (istate > nstates_conv) THEN
            ! guess non-restarted wave functions using appropriate virtual molecular orbitals
            ispin = guess_vectors%ispins(istate)
            CALL cp_fm_to_fm(gs_mos(ispin)%matrix_mos_virt, evects(ispin, istate)%matrix, &
                             ncol=1, source_start=guess_vectors%imos_virt(istate), &
                             target_start=guess_vectors%imos_occ(istate))

            IF (para_env_sub%mepos == 0) &
               evals_prev(istate) = guess_vectors%evals(istate)
         END IF
      END DO

      IF (nstates < nstates_conv) THEN
         CALL mp_sum(evals_prev(nstates_conv+1:nstates), para_env_global%group)
      END IF

      IF (tddfpt_control%is_restart) THEN
         DO istate = 1, nstates
            DO ispin = 1, nspins
               CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=evects(ispin, istate)%matrix, &
                                                         fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                         blacs_env_global=blacs_env_global, &
                                                         nsubgroups=nsubgroups, subgroup=subgroup)
            END DO
         END DO

         DO istate = tddfpt_control%nreortho, 1, -1
            CALL tddfpt_orthogonalize_psi1_psi0(evects=evects(:, 1:nstates), S_C0_C0T=S_C0_C0T_sub)
            CALL tddfpt_orthonormalize_psi1_psi1(evects=evects(:, 1:nstates), nvects_new=nstates, &
                                                 matrix_s=matrix_s_sub)
            is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(evects=evects(:, 1:nstates), S_C0=S_C0_sub, &
                                                            max_norm=tddfpt_control%orthogonal_eps)
            IF (.NOT. is_nonortho) EXIT
         END DO

         DO istate = 1, nstates
            DO ispin = 1, nspins
               IF (.NOT. ASSOCIATED(ritz_vects(ispin, istate)%matrix)) &
                  CALL cp_fm_release(evects(ispin, istate)%matrix)
            END DO
         END DO
      END IF

      nvects_exists = 0
      nvects_new = nstates

      ! *** perform Davidson iterations ***
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
      energy_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "DETAILED_ENERGY", extension=".tddfptLog")

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A,/)') 'TDDFPT WAVEFUNCTION OPTIMIZATION'
         WRITE (log_unit, '(5X,A,T15,A,T24,A,T40,A)') "Step", "Time", "Convergence", "Conv. states"
         WRITE (log_unit, '(1X,50("-"))')
      END IF

      CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

      niters = tddfpt_control%niters
      t1 = m_walltime()

      DO iter = 1, niters
         CALL cp_iterate(logger%iter_info, last=.FALSE., iter_nr=iter)

         CALL tddfpt_compute_Aop_evects(Aop_evects=Aop_evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                        evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                        fm_struct_evects=fm_struct_evects, gs_mos=gs_mos, &
                                        is_rks_triplets=tddfpt_control%rks_triplets, &
                                        matrix_ks_sub=matrix_ks_sub, matrix_s_sub=matrix_s_sub, &
                                        qs_env=qs_env, xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, &
                                        xc_section=xc_section, do_hfx=do_hfx, do_admm=dft_control%do_admm, &
                                        matrix_s_aux_fit_sub=matrix_s_aux_fit_sub, blacs_env_sub=blacs_env_sub, &
                                        pw_env_sub=pw_env_sub, task_list_orb_sub=task_list_orb_sub, &
                                        task_list_aux_fit_sub=task_list_aux_fit_sub, &
                                        nsubgroups=nsubgroups, subgroup=subgroup)

         DO istate = nvects_exists+1, nvects_exists+nvects_new
            DO ispin = 1, nspins
               CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=evects(ispin, istate)%matrix, &
                                                         fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                         blacs_env_global=blacs_env_global, &
                                                         nsubgroups=nsubgroups, subgroup=subgroup)
               CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=Aop_evects(ispin, istate)%matrix, &
                                                         fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                         blacs_env_global=blacs_env_global, &
                                                         nsubgroups=nsubgroups, subgroup=subgroup)
            END DO
         END DO

         CALL tddfpt_compute_ritz_vects(ritz_vects=ritz_vects(:, 1:nstates), Aop_ritz=Aop_ritz(:, 1:nstates), &
                                        evals=evals_last(1:nvects_exists+nvects_new), &
                                        Aop_evects=Aop_evects(:, 1:nvects_exists+nvects_new), &
                                        evects=evects(:, 1:nvects_exists+nvects_new), gs_mos=gs_mos, &
                                        blacs_env_global=blacs_env_global)

         conv = MAXVAL(ABS(evals_last(1:nstates)-evals_prev(1:nstates)))
         IF (conv > tddfpt_control%conv .OR. .NOT. is_restarted) THEN
            ! expand the active space
            nvects_exists = nvects_exists+nvects_new
            IF (nvects_exists < max_krylov_vects .AND. conv > tddfpt_control%conv) THEN
               IF (nvects_exists+nvects_new > max_krylov_vects) &
                  nvects_new = max_krylov_vects-nvects_exists

               CALL tddfpt_compute_residual_vects(residual_vects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                  evals=evals_last(1:nvects_new), Aop_ritz=Aop_ritz(:, 1:nvects_new), &
                                                  ritz_vects=ritz_vects(:, 1:nvects_new), fm_struct_evects=fm_struct_evects, &
                                                  gs_mos=gs_mos, guess_vectors=guess_vectors, matrix_s=matrix_s_sub)

               DO istate = nvects_exists+1, nvects_exists+nvects_new
                  DO ispin = 1, nspins
                     CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=evects(ispin, istate)%matrix, &
                                                               fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                               blacs_env_global=blacs_env_global, &
                                                               nsubgroups=nsubgroups, subgroup=subgroup)
                  END DO
               END DO

               DO istate = tddfpt_control%nreortho, 1, -1
                  CALL tddfpt_orthogonalize_psi1_psi0(evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                      S_C0_C0T=S_C0_C0T_sub)
                  CALL tddfpt_orthonormalize_psi1_psi1(evects=evects(:, 1:nvects_exists+nvects_new), &
                                                       nvects_new=nvects_new, matrix_s=matrix_s_sub)
                  is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(evects=evects(:, nvects_exists+1:nvects_exists+nvects_new), &
                                                                  S_C0=S_C0_sub, max_norm=tddfpt_control%orthogonal_eps)
                  IF (.NOT. is_nonortho) EXIT
               END DO

               DO istate = 1, nvects_new
                  DO ispin = 1, nspins
                     IF (.NOT. ASSOCIATED(ritz_vects(ispin, istate)%matrix)) &
                        CALL cp_fm_release(evects(ispin, nvects_exists+istate)%matrix)
                  END DO
               END DO
            ELSE
               is_nonortho = .FALSE.
            END IF

            ! restart Davidson iterations
            IF (nvects_exists >= max_krylov_vects .OR. is_nonortho .OR. &
                (conv <= tddfpt_control%conv .AND. .NOT. is_restarted)) THEN

               DO istate = SIZE(evects, 2), 1, -1
                  DO ispin = nspins, 1, -1
                     IF (ASSOCIATED(Aop_evects(ispin, istate)%matrix)) &
                        CALL cp_fm_release(Aop_evects(ispin, istate)%matrix)
                     IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
                        CALL cp_fm_release(evects(ispin, istate)%matrix)
                  END DO
               END DO

               DO istate = 1, nstates
                  DO ispin = 1, nspins
                     IF (ASSOCIATED(ritz_vects(ispin, istate)%matrix)) THEN
                        CALL cp_fm_create(evects(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
                        CALL cp_fm_to_fm(ritz_vects(ispin, istate)%matrix, evects(ispin, istate)%matrix)
                     END IF

                     CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=evects(ispin, istate)%matrix, &
                                                               fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                               blacs_env_global=blacs_env_global, &
                                                               nsubgroups=nsubgroups, subgroup=subgroup)
                  END DO
               END DO

               DO istate = 1, max_reortho
                  CALL tddfpt_orthogonalize_psi1_psi0(evects=evects(:, 1:nstates), S_C0_C0T=S_C0_C0T_sub)
                  CALL tddfpt_orthonormalize_psi1_psi1(evects=evects(:, 1:nstates), nvects_new=nstates, &
                                                       matrix_s=matrix_s_sub)
                  is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(evects=evects(:, 1:nstates), &
                                                                  S_C0=S_C0_sub, max_norm=tddfpt_control%orthogonal_eps)
                  IF (.NOT. is_nonortho) EXIT
               END DO

               DO istate = 1, nvects_new
                  DO ispin = 1, nspins
                     IF (.NOT. ASSOCIATED(ritz_vects(ispin, istate)%matrix)) &
                        CALL cp_fm_release(evects(ispin, istate)%matrix)
                  END DO
               END DO

               nvects_exists = 0
               nvects_new = nstates
               IF (conv > tddfpt_control%conv) &
                  is_restarted = .TRUE.
            END IF
         END IF

         t2 = m_walltime()

         IF (energy_unit > 0) THEN
            WRITE (energy_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
            DO istate = 1, nstates
               WRITE (energy_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
                  evals_last(istate)*evolt, (evals_last(istate)-evals_prev(istate))*evolt
            END DO
            WRITE (energy_unit, *)
            CALL m_flush(energy_unit)
         END IF

         IF (log_unit > 0) THEN
            nstates_conv = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate)-evals_prev(istate)) <= tddfpt_control%conv) &
                  nstates_conv = nstates_conv+1
            END DO

            WRITE (log_unit, '(1X,I8,T12,F7.1,T24,ES11.4,T42,I8)') iter, t2-t1, conv, nstates_conv
            CALL m_flush(log_unit)
         END IF

         t1 = t2
         evals_prev(1:nstates) = evals_last(1:nstates)

         IF (conv <= tddfpt_control%conv) THEN
            IF (is_restarted) THEN
               EXIT
            ELSE
               is_restarted = .TRUE.
            END IF
         END IF

         IF (log_unit > 0 .AND. (nvects_exists == 0 .OR. conv <= tddfpt_control%conv)) THEN
            WRITE (log_unit, '(1X,10("-"),1X,A,1X,11("-"))') "Restart Davidson iterations"
            CALL m_flush(log_unit)
         END IF

         CALL tddfpt_write_restart(evects=ritz_vects(:, 1:nstates), evals=evals_prev(1:nstates), &
                                   gs_mos=gs_mos, logger=logger, tddfpt_print_section=tddfpt_print_section, &
                                   blacs_env_global=blacs_env_global)
      END DO

      CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=iter)
      CALL tddfpt_write_restart(evects=ritz_vects(:, 1:nstates), evals=evals_prev(1:nstates), &
                                gs_mos=gs_mos, logger=logger, tddfpt_print_section=tddfpt_print_section, &
                                blacs_env_global=blacs_env_global)

      CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

      IF (log_unit > 0) THEN
         IF (iter <= niters) THEN
            CALL integer_to_string(iter, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run converged in "//TRIM(nstates_str)//" iteration(s) ***"
         ELSE
            CALL integer_to_string(iter-1, nstates_str)
            WRITE (log_unit, '(/,1X,A)') "*** TDDFPT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ***"
         END IF
      END IF

      CALL cp_print_key_finished_output(energy_unit, logger, tddfpt_print_section, "DETAILED_ENERGY")
      CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "ITERATION_INFO")

      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit()

      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL tddfpt_fm_broadcast_to_all_subgroups(fm_sub=ritz_vects(ispin, istate)%matrix, &
                                                      fm_struct_sub=fm_struct_evects(ispin)%struct, &
                                                      blacs_env_global=blacs_env_global, &
                                                      nsubgroups=nsubgroups, subgroup=subgroup)
         END DO
      END DO
      CALL tddfpt_orthogonalize_psi1_psi0(evects=ritz_vects(:, 1:nstates), S_C0_C0T=S_C0_C0T_sub)
      CALL tddfpt_orthonormalize_psi1_psi1(evects=ritz_vects(:, 1:nstates), nvects_new=nstates, &
                                           matrix_s=matrix_s_sub)
      DO istate = 1, nstates
         DO ispin = 1, nspins
            IF (.NOT. ASSOCIATED(Aop_ritz(ispin, istate)%matrix)) &
               CALL cp_fm_release(ritz_vects(ispin, istate)%matrix)
         END DO
      END DO

      CALL tddfpt_print_summary(log_unit=log_unit, evects=ritz_vects(:, 1:nstates), evals=evals_last(1:nstates), &
                                mult=mult, fm_struct_evects=fm_struct_evects, gs_mos=gs_mos, matrix_s=matrix_s_sub, &
                                dipole_op_xyz=dipole_op_xyz_sub, dipole_form=tddfpt_control%dipole_form, &
                                min_amplitude=tddfpt_control%min_excitation_amplitude, &
                                blacs_env_global=blacs_env_global)

      ! -- clean up all useless stuff
      DEALLOCATE (evals_last, evals_prev)

      DO istate = SIZE(ritz_vects, 2), 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(Aop_ritz(ispin, istate)%matrix)) &
               CALL cp_fm_release(Aop_ritz(ispin, istate)%matrix)
            IF (ASSOCIATED(ritz_vects(ispin, istate)%matrix)) &
               CALL cp_fm_release(ritz_vects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (ritz_vects, Aop_ritz)

      DO istate = SIZE(evects, 2), 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(Aop_evects(ispin, istate)%matrix)) &
               CALL cp_fm_release(Aop_evects(ispin, istate)%matrix)
            IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
               CALL cp_fm_release(evects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (evects, Aop_evects)

      CALL xc_dset_release(xc_deriv_set)
      CALL xc_rho_set_release(xc_rho_set)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_C0_C0T_sub(ispin)%matrix)
         CALL cp_fm_release(S_C0_sub(ispin)%matrix)
      END DO
      DEALLOCATE (S_C0_C0T_sub, S_C0_sub)

      CALL cp_dbcsr_deallocate_matrix_set(dipole_op_xyz_sub)
      CALL cp_dbcsr_deallocate_matrix(matrix_s_sub)
      CALL cp_dbcsr_deallocate_matrix_set(matrix_ks_sub)

      DEALLOCATE (fm_struct_evects)

      CALL tddfpt_release_guess_vectors(guess_vectors=guess_vectors)

      DO ispin = nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos=gs_mos(ispin))
      END DO
      DEALLOCATE (gs_mos)

      IF (ASSOCIATED(matrix_s_aux_fit_sub)) &
         CALL cp_dbcsr_deallocate_matrix(matrix_s_aux_fit_sub)

      IF (ASSOCIATED(task_list_aux_fit_sub)) &
         CALL deallocate_task_list(task_list_aux_fit_sub)
      CALL deallocate_task_list(task_list_orb_sub)

      IF (ASSOCIATED(sab_aux_fit_sub)) THEN
         DO istate = SIZE(sab_aux_fit_sub), 1, -1
            CALL deallocate_neighbor_list_set(sab_aux_fit_sub(istate)%neighbor_list_set)
         END DO
         DEALLOCATE (sab_aux_fit_sub)
      END IF

      DO istate = SIZE(sab_orb_sub), 1, -1
         CALL deallocate_neighbor_list_set(sab_orb_sub(istate)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb_sub)

      CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (dbcsr_dist_sub)
      CALL distribution_2d_release(dist_2d_sub)

      CALL pw_env_release(pw_env_sub)

      IF (tddfpt_control%mgrid_is_explicit) &
         CALL restore_qs_mgrid(qs_control, mgrid_saved)

      CALL cp_blacs_env_release(blacs_env_sub)
      CALL cp_para_env_release(para_env_sub)
      ! do not need to call mp_comm_free() as the MPI sub-communicator mpi_comm_sub
      ! has already been released by the subroutine cp_para_env_release()
      DEALLOCATE (subgroup)

      CALL timestop(handle)
   END SUBROUTINE tddfpt

! **************************************************************************************************
!> \brief Replace the global multi-grid related parameters in qs_control by the ones given in the
!>        TDDFPT/MGRID subsection. The original parameters are stored into the 'mgrid_saved'
!>        variable.
!> \param qs_control     Quickstep control parameters (modified on exit)
!> \param tddfpt_control TDDFPT control parameters
!> \param mgrid_saved    structure to hold global MGRID-related parameters (initialised on exit)
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!> \note the code to build the 'e_cutoff' list was taken from the subroutine read_mgrid_section()
! **************************************************************************************************
   SUBROUTINE init_tddfpt_mgrid(qs_control, tddfpt_control, mgrid_saved)
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(mgrid_saved_parameters), INTENT(out)          :: mgrid_saved

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_tddfpt_mgrid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igrid, ngrids

      CALL timeset(routineN, handle)

      ! save global plane-wave grid parameters to the variable 'mgrid_saved'
      mgrid_saved%commensurate_mgrids = qs_control%commensurate_mgrids
      mgrid_saved%realspace_mgrids = qs_control%realspace_mgrids
      mgrid_saved%skip_load_balance = qs_control%skip_load_balance_distributed
      mgrid_saved%cutoff = qs_control%cutoff
      mgrid_saved%progression_factor = qs_control%progression_factor
      mgrid_saved%relative_cutoff = qs_control%relative_cutoff
      mgrid_saved%e_cutoff => qs_control%e_cutoff

      ! set parameters from 'tddfpt_control' as default ones for all newly allocated plane-wave grids
      qs_control%commensurate_mgrids = tddfpt_control%mgrid_commensurate_mgrids
      qs_control%realspace_mgrids = tddfpt_control%mgrid_realspace_mgrids
      qs_control%skip_load_balance_distributed = tddfpt_control%mgrid_skip_load_balance
      qs_control%cutoff = tddfpt_control%mgrid_cutoff
      qs_control%progression_factor = tddfpt_control%mgrid_progression_factor
      qs_control%relative_cutoff = tddfpt_control%mgrid_relative_cutoff

      ALLOCATE (qs_control%e_cutoff(tddfpt_control%mgrid_ngrids))
      ngrids = tddfpt_control%mgrid_ngrids
      IF (ASSOCIATED(tddfpt_control%mgrid_e_cutoff)) THEN
         ! following read_mgrid_section() there is a magic scale factor there (0.5_dp)
         DO igrid = 1, ngrids
            qs_control%e_cutoff(igrid) = tddfpt_control%mgrid_e_cutoff(igrid)*0.5_dp
         END DO
         ! round 'qs_control%cutoff' upward to the nearest sub-grid's cutoff value;
         ! here we take advantage of the fact that the array 'e_cutoff' has been sorted in descending order
         DO igrid = ngrids, 1, -1
            IF (qs_control%cutoff <= qs_control%e_cutoff(igrid)) THEN
               qs_control%cutoff = qs_control%e_cutoff(igrid)
               EXIT
            END IF
         END DO
         ! igrid == 0 if qs_control%cutoff is larger than the largest manually provided cutoff value;
         ! use the largest actual value
         IF (igrid <= 0) &
            qs_control%cutoff = qs_control%e_cutoff(1)
      ELSE
         qs_control%e_cutoff(1) = qs_control%cutoff
         DO igrid = 2, ngrids
            qs_control%e_cutoff(igrid) = qs_control%e_cutoff(igrid-1)/qs_control%progression_factor
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE init_tddfpt_mgrid

! **************************************************************************************************
!> \brief Restore the global multi-grid related parameters stored in the 'mgrid_saved' variable.
!> \param qs_control  Quickstep control parameters (modified on exit)
!> \param mgrid_saved structure that holds global MGRID-related parameters
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE restore_qs_mgrid(qs_control, mgrid_saved)
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(mgrid_saved_parameters), INTENT(in)           :: mgrid_saved

      CHARACTER(LEN=*), PARAMETER :: routineN = 'restore_qs_mgrid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(qs_control%e_cutoff)) &
         DEALLOCATE (qs_control%e_cutoff)

      qs_control%commensurate_mgrids = mgrid_saved%commensurate_mgrids
      qs_control%realspace_mgrids = mgrid_saved%realspace_mgrids
      qs_control%skip_load_balance_distributed = mgrid_saved%skip_load_balance
      qs_control%cutoff = mgrid_saved%cutoff
      qs_control%progression_factor = mgrid_saved%progression_factor
      qs_control%relative_cutoff = mgrid_saved%relative_cutoff
      qs_control%e_cutoff => mgrid_saved%e_cutoff

      CALL timestop(handle)
   END SUBROUTINE restore_qs_mgrid

! **************************************************************************************************
!> \brief Distribute atoms across the two-dimensional grid of processors.
!> \param distribution_2d new two-dimensional distribution of pairs of particles
!>                        (allocated and initialised on exit)
!> \param dbcsr_dist      new DBCSR distribution (allocated and initialised on exit)
!> \param blacs_env       BLACS parallel environment
!> \param qs_env          Quickstep environment
!> \param log_unit        output unit
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_build_distribution_2d(distribution_2d, dbcsr_dist, blacs_env, qs_env, log_unit)
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_build_distribution_2d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, input=input, &
                      molecule_kind_set=molecule_kind_set, molecule_set=molecule_set, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)

      NULLIFY (distribution_2d)
      CALL distribute_molecules_2d(cell=cell, &
                                   atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   qs_kind_set=qs_kind_set, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   distribution_2d=distribution_2d, &
                                   blacs_env=blacs_env, &
                                   force_env_section=input)

      ALLOCATE (dbcsr_dist)
      CALL cp_dbcsr_dist2d_to_dist(distribution_2d, dbcsr_dist, log_unit)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_build_distribution_2d

! **************************************************************************************************
!> \brief Build task and neighbour lists for the given plane wave environment and basis set.
!> \param task_list           new task list (allocated and initialised on exit)
!> \param sab                 new list of neighbours (allocated and initialised on exit)
!> \param basis_type          type of the basis set
!> \param distribution_2d     two-dimensional distribution of pairs of particles
!> \param pw_env              plane wave environment
!> \param qs_env              Quickstep environment
!> \param skip_load_balance   do not perform load balancing
!> \param reorder_grid_ranks  re-optimise grid ranks and re-create the real-space grid descriptor
!>                            as well as grids
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_build_tasklist(task_list, sab, basis_type, distribution_2d, pw_env, qs_env, &
                                    skip_load_balance, reorder_grid_ranks)
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab
      CHARACTER(len=*), INTENT(in)                       :: basis_type
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: skip_load_balance, reorder_grid_ranks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_build_tasklist', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, nkinds
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(kind=dp)                                      :: subcells
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: orb_radius
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, input=input, &
                      ks_env=ks_env, local_particles=local_particles, molecule_set=molecule_set, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)

      nkinds = SIZE(atomic_kind_set)

      ALLOCATE (atom2d(nkinds))
      CALL atom2d_build(atom2d, local_particles, distribution_2d, atomic_kind_set, &
                        molecule_set, molecule_only=.FALSE., particle_set=particle_set)

      ALLOCATE (orb_present(nkinds))
      ALLOCATE (orb_radius(nkinds))
      ALLOCATE (pair_radius(nkinds, nkinds))

      DO ikind = 1, nkinds
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=basis_type)
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         ENDIF
      END DO

      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

      NULLIFY (sab)
      CALL section_vals_val_get(input, "DFT%SUBCELLS", r_val=subcells)
      CALL build_neighbor_lists(sab, particle_set, atom2d, cell, pair_radius, &
                                mic=.FALSE., subcells=subcells, molecular=.FALSE., name="sab_orb")

      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d, orb_present, orb_radius, pair_radius)

      CALL allocate_task_list(task_list)
      CALL generate_qs_task_list(ks_env, task_list, &
                                 reorder_rs_grid_ranks=reorder_grid_ranks, soft_valid=.FALSE., &
                                 basis_type=basis_type, skip_load_balance_distributed=skip_load_balance, &
                                 pw_env_external=pw_env, sab_orb_external=sab)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_build_tasklist

! **************************************************************************************************
!> \brief Create a DBCSR matrix based on a template matrix, distribution object, and the list of
!>        neighbours.
!> \param matrix     matrix to create
!> \param template   template matrix
!> \param dbcsr_dist DBCSR distribution
!> \param sab        list of neighbours
! **************************************************************************************************
   SUBROUTINE tddfpt_dbcsr_create_by_dist(matrix, template, dbcsr_dist, sab)
      TYPE(cp_dbcsr_type), POINTER                       :: matrix, template
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_dbcsr_create_by_dist', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: matrix_type
      CHARACTER(len=default_string_length)               :: matrix_name
      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(template))
      CALL cp_dbcsr_get_info(template, row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                             name=matrix_name, matrix_type=matrix_type)

      IF (ASSOCIATED(matrix)) THEN
         CALL cp_dbcsr_release(matrix)
      ELSE
         ALLOCATE (matrix)
      END IF
      CALL cp_dbcsr_init(matrix)

      CALL cp_dbcsr_create(matrix, matrix_name, dbcsr_dist, matrix_type, row_blk_sizes, col_blk_sizes, nze=0)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix, sab)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_dbcsr_create_by_dist

! **************************************************************************************************
!> \brief Replicate a globally distributed matrix across all sub-groups. As the result,
!>        every sub-group will hold a local subgroup-specific copy of the original globally
!>        distributed matrix.
!>
!>                                 |--------------------|
!>                         fm_src  |  0    1    2    3  |
!>                                 |--------------------|
!>                                    /  MPI  ranks  \
!>                                  |/_              _\|
!>                    |--------------------|    |--------------------|
!> fm_dest_subgroup0  |     0        1     |    |     2        3     |  fm_dest_subgroup1
!>                    |--------------------|    |--------------------|
!>                          subgroup 0                subgroup 1
!>
!> \param fm_src      globally distributed matrix to replicate
!> \param fm_dest_sub subgroup-specific copy of the replicated matrix
!> \param nsubgroups  number of sub-groups (e.g. 2 means there are two sub-group with
!>                    indices 0 and 1)
!> \param subgroup    list that indicates which MPI process belongs to which sub-group
!>                    (e.g. {0, 0, 1, 1} means that the processors with MPI ranks 0 and 1 are
!>                    the members of subgroup 0, while the processors with MPI ranks 2 and 3
!>                    belong to subgroup 1)
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_fm_replicate_across_subgroups(fm_src, fm_dest_sub, nsubgroups, subgroup)
      TYPE(cp_fm_type), POINTER                          :: fm_src, fm_dest_sub
      INTEGER, INTENT(in)                                :: nsubgroups
      INTEGER, DIMENSION(0:), INTENT(in)                 :: subgroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_fm_replicate_across_subgroups', &
         routineP = moduleN//':'//routineN

      INTEGER :: context_global, handle, igroup, igroup_local, ncols_global_dest, &
         ncols_global_src, nrows_global_dest, nrows_global_src
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_type), POINTER                          :: fm_null
      TYPE(cp_para_env_type), POINTER                    :: para_env_global

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_src, nrow_global=nrows_global_src, ncol_global=ncols_global_src, &
                          context=blacs_env_global, para_env=para_env_global)
      CALL cp_fm_get_info(fm_dest_sub, nrow_global=nrows_global_dest, ncol_global=ncols_global_dest)

      CPASSERT(nrows_global_src == nrows_global_dest)
      CPASSERT(ncols_global_src == ncols_global_dest)

      NULLIFY (fm_null)
      context_global = blacs_env_global%group
      igroup_local = subgroup(para_env_global%mepos)

      DO igroup = 0, nsubgroups-1
         IF (igroup == igroup_local) THEN
            CALL cp_fm_copy_general(source=fm_src, destination=fm_dest_sub, global_context=context_global)
         ELSE
            CALL cp_fm_copy_general(source=fm_src, destination=fm_null, global_context=context_global)
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_fm_replicate_across_subgroups

! **************************************************************************************************
!> \brief Broadcast a subgroup-specific dense matrix to all other sub-groups. The number of
!>        processors in each sub-group can be different.
!>
!>    |--------------------|       |--------------------|       |--------------------|
!>    |     0        1     | <---- |     2        3     | ----> |     4        5     |
!>    |--------------------|       |--------------------|       |--------------------|
!>          subgroup 0                   subgroup 1                   subgroup 2
!>
!> \param fm_sub           matrix to replicate. One and only one sub-group should hold the given
!>                         matrix; other sub-groups should pass a nullified pointer here which
!>                         will points to a local copy of the replicated matrix on exit
!> \param fm_struct_sub    dense matrix structure for newly allocated matrices
!> \param blacs_env_global global BLACS parallel environment
!> \param nsubgroups       number of sub-groups
!> \param subgroup         list that indicates which MPI process belongs to which sub-group
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_fm_broadcast_to_all_subgroups(fm_sub, fm_struct_sub, blacs_env_global, &
                                                   nsubgroups, subgroup)
      TYPE(cp_fm_type), POINTER                          :: fm_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      INTEGER, INTENT(in)                                :: nsubgroups
      INTEGER, DIMENSION(0:), INTENT(in)                 :: subgroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_fm_broadcast_to_all_subgroups', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: context_global, handle, igroup, &
                                                            igroup_local, igroup_root
      TYPE(cp_fm_type), POINTER                          :: fm_null
      TYPE(cp_para_env_type), POINTER                    :: para_env_global, para_env_sub

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_get(fm_struct_sub, para_env=para_env_sub)
      CALL get_blacs_info(blacs_env_global, para_env=para_env_global)

      context_global = blacs_env_global%group
      igroup_local = subgroup(para_env_global%mepos)
      igroup_root = 0

      IF (ASSOCIATED(fm_sub)) THEN
         IF (para_env_sub%mepos == para_env_sub%source) &
            igroup_root = igroup_local
      ELSE
         CALL cp_fm_create(fm_sub, fm_struct_sub)
      END IF

      CALL mp_max(igroup_root, para_env_global%group)

      NULLIFY (fm_null)
      DO igroup = 0, nsubgroups-1
         IF (igroup /= igroup_root) THEN
            IF (igroup_local == igroup_root) THEN
               CALL cp_fm_copy_general(source=fm_sub, destination=fm_null, global_context=context_global)
            ELSE
               IF (igroup_local == igroup) THEN
                  CALL cp_fm_copy_general(source=fm_null, destination=fm_sub, global_context=context_global)
               ELSE
                  CALL cp_fm_copy_general(source=fm_null, destination=fm_null, global_context=context_global)
               END IF
            END IF
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE tddfpt_fm_broadcast_to_all_subgroups

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param gs_mos           structure to store occupied and virtual molecular orbitals
!>                         (allocated and initialised on exit)
!> \param mo_set           ground state molecular orbitals for a given spin
!> \param blacs_env_global BLACS parallel environment envolving all processors
!> \param blacs_env_sub    BLACS parallel environment for local subgroup
!> \param cholesky_method  Cholesky method to compute the inverse overlap matrix
!> \param matrix_ks_global Kohn-Sham matrix for a given spin distributed among all processors
!> \param matrix_s_global  overlap matrix distributed among all processors
!> \param nsubgroups       number of MPI subgroups
!> \param subgroup         array of size 'blacs_env_global%para_env%num_pe' that indicates
!>                         which MPI process belongs to which subgroup
!> \par History
!>    * 05.2016 created as tddfpt_lumos() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_init_ground_state_mos(gs_mos, mo_set, blacs_env_global, blacs_env_sub, &
                                           cholesky_method, matrix_ks_global, matrix_s_global, &
                                           nsubgroups, subgroup)
      TYPE(tddfpt_ground_state_mos), INTENT(out)         :: gs_mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global, blacs_env_sub
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_ks_global, matrix_s_global
      INTEGER, INTENT(in)                                :: nsubgroups
      INTEGER, DIMENSION(0:), INTENT(in)                 :: subgroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_ground_state_mos', &
         routineP = moduleN//':'//routineN

      INTEGER :: cholesky_method_inout, handle, icol, icol_global, imo, irow, irow_global, nao, &
         ncol_local, nelectrons, nmo_occ, nmo_scf, nmo_virt, nrow_local, sign_int
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: minrow_neg_array, minrow_pos_array, &
                                                            sum_sign_array
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(kind=dp)                                      :: element, maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block
      TYPE(cp_fm_pool_type), POINTER                     :: ao_ao_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct, ao_vect_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env_global
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      CALL get_blacs_info(blacs_env_global, para_env=para_env_global)
      IF (debug_this_module) THEN
         CPASSERT(SIZE(subgroup) == para_env_global%num_pe)
      END IF

      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=nmo_occ, maxocc=maxocc, &
                      nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_virt = nao-nmo_occ

      IF (nmo_virt <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'Unable to generate virtual molecular orbitals. Please use a larger atomic basis set.')

      ! ++ allocate storage space for gs_mos
      NULLIFY (ao_mo_occ_fm_struct, ao_mo_virt_fm_struct, ao_vect_fm_struct)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env_sub)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt, context=blacs_env_sub)
      CALL cp_fm_struct_create(ao_vect_fm_struct, nrow_global=nao, ncol_global=1, context=blacs_env_sub)

      NULLIFY (gs_mos%matrix_mos_occ, gs_mos%matrix_mos_virt)
      CALL cp_fm_create(gs_mos%matrix_mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(gs_mos%matrix_mos_virt, ao_mo_virt_fm_struct)

      ALLOCATE (gs_mos%vectors_mo_virt(nmo_virt))
      DO imo = 1, nmo_virt
         NULLIFY (gs_mos%vectors_mo_virt(imo)%matrix)
         CALL cp_fm_create(gs_mos%vectors_mo_virt(imo)%matrix, ao_vect_fm_struct)
      END DO

      ALLOCATE (gs_mos%evals_occ(nmo_occ))
      ALLOCATE (gs_mos%evals_virt(nmo_virt))
      ALLOCATE (gs_mos%phases_occ(nmo_occ))

      CALL cp_fm_struct_release(ao_vect_fm_struct)
      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      ! ++ set of molecular orbitals
      NULLIFY (ao_ao_fm_struct, ao_ao_fm_pool)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_global)
      CALL fm_pool_create(ao_ao_fm_pool, ao_ao_fm_struct)

      NULLIFY (mos_extended, mo_coeff_extended, mo_evals_extended, mo_occ_extended)
      CALL allocate_mo_set(mos_extended, nao, nao, nelectrons, &
                           REAL(nelectrons, dp), maxocc, flexible_electron_count=0.0_dp)
      CALL init_mo_set(mos_extended, fm_pool=ao_ao_fm_pool, name="mos-extended")
      CALL fm_pool_release(ao_ao_fm_pool)
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
      mo_occ_extended(nmo_scf+1:) = 0.0_dp

      ! ++ allocate temporary matrices
      NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
      CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
      CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
      CALL cp_fm_create(work_fm, ao_ao_fm_struct)

      CALL cp_fm_struct_release(ao_ao_fm_struct)

      ! some stuff from the subroutine general_eigenproblem()
      CALL copy_dbcsr_to_fm(matrix_s_global, ortho_fm)
      CALL copy_dbcsr_to_fm(matrix_ks_global, matrix_ks_fm)

      IF (cholesky_method == cholesky_dbcsr) THEN
         CPABORT('CHOLESKY DBCSR_INVERSE is not implemented in TDDFT.')
      ELSE IF (cholesky_method == cholesky_off) THEN
         CPABORT('CHOLESKY OFF is not implemented in TDDFT.')
      ELSE
         CALL cp_fm_cholesky_decompose(ortho_fm)
         IF (cholesky_method == cholesky_inverse) THEN
            CALL cp_fm_triangular_invert(ortho_fm)
         END IF

         ! need this temporary variable, as the subroutine eigensolver() is going to update it.
         cholesky_method_inout = cholesky_method
         CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=mos_extended, ortho=ortho_fm, &
                          work=work_fm, cholesky_method=cholesky_method_inout, &
                          do_level_shift=.FALSE., level_shift=0.0_dp, matrix_u_fm=null(), use_jacobi=.FALSE.)
      END IF

      ! -- clean up needless matrices
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(ortho_fm)
      CALL cp_fm_release(matrix_ks_fm)

      ! compute the phase of molecular orbitals;
      ! matrix work_fm holds occupied molecular orbital coefficients distributed among all the processors
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env_global)
      CALL cp_fm_create(work_fm, ao_mo_occ_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      CALL cp_fm_to_fm(mo_coeff_extended, work_fm, ncol=nmo_occ, source_start=1, target_start=1)
      CALL cp_fm_get_info(work_fm, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=my_block)

      ALLOCATE (minrow_neg_array(nmo_occ), minrow_pos_array(nmo_occ), sum_sign_array(nmo_occ))
      minrow_neg_array(:) = nao
      minrow_pos_array(:) = nao
      sum_sign_array(:) = 0
      DO icol = 1, ncol_local
         icol_global = col_indices(icol)

         DO irow = 1, nrow_local
            element = my_block(irow, icol)

            sign_int = 0
            IF (element >= EPSILON(0.0_dp)) THEN
               sign_int = 1
            ELSE IF (element <= -EPSILON(0.0_dp)) THEN
               sign_int = -1
            END IF

            sum_sign_array(icol_global) = sum_sign_array(icol_global)+sign_int

            irow_global = row_indices(irow)
            IF (sign_int > 0) THEN
               IF (minrow_pos_array(icol_global) > irow_global) &
                  minrow_pos_array(icol_global) = irow_global
            ELSE IF (sign_int < 0) THEN
               IF (minrow_neg_array(icol_global) > irow_global) &
                  minrow_neg_array(icol_global) = irow_global
            END IF
         END DO
      END DO

      CALL mp_sum(sum_sign_array, para_env_global%group)
      CALL mp_min(minrow_neg_array, para_env_global%group)
      CALL mp_min(minrow_pos_array, para_env_global%group)

      DO icol = 1, nmo_occ
         IF (sum_sign_array(icol) > 0) THEN
            ! most of the expansion coefficients are positive => MO's phase = +1
            gs_mos%phases_occ(icol) = 1.0_dp
         ELSE IF (sum_sign_array(icol) < 0) THEN
            ! most of the expansion coefficients are negative => MO's phase = -1
            gs_mos%phases_occ(icol) = -1.0_dp
         ELSE
            ! equal number of positive and negative expansion coefficients
            IF (minrow_pos_array(icol) <= minrow_neg_array(icol)) THEN
               ! the first positive expansion coefficient has a lower index then
               ! the first negative expansion coefficient; MO's phase = +1
               gs_mos%phases_occ(icol) = 1.0_dp
            ELSE
               ! MO's phase = -1
               gs_mos%phases_occ(icol) = -1.0_dp
            END IF
         END IF
      END DO

      DEALLOCATE (minrow_neg_array, minrow_pos_array, sum_sign_array)
      CALL cp_fm_release(work_fm)

      ! replicate molecular coefficients across subgroups
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_sub)
      CALL cp_fm_create(work_fm, ao_ao_fm_struct)
      CALL cp_fm_struct_release(ao_ao_fm_struct)

      CALL tddfpt_fm_replicate_across_subgroups(fm_src=mo_coeff_extended, fm_dest_sub=work_fm, &
                                                nsubgroups=nsubgroups, subgroup=subgroup)

      ! return the requested occupied and virtual molecular orbitals and corresponding orbital energies
      CALL cp_fm_to_fm(work_fm, gs_mos%matrix_mos_occ, ncol=nmo_occ, source_start=1, target_start=1)
      gs_mos%evals_occ(1:nmo_occ) = mo_evals_extended(1:nmo_occ)

      CALL cp_fm_to_fm(work_fm, gs_mos%matrix_mos_virt, ncol=nmo_virt, source_start=nmo_occ+1, target_start=1)
      DO imo = 1, nmo_virt
         CALL cp_fm_to_fm(work_fm, gs_mos%vectors_mo_virt(imo)%matrix, ncol=1, source_start=nmo_occ+imo, target_start=1)
      END DO
      gs_mos%evals_virt(1:nmo_virt) = mo_evals_extended(nmo_occ+1:nmo_occ+nmo_virt)

      CALL cp_fm_release(work_fm)

      CALL deallocate_mo_set(mos_extended)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_ground_state_mos

! **************************************************************************************************
!> \brief Release molecular orbitals.
!> \param gs_mos          structure that holds occupied and virtual molecular orbitals
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_ground_state_mos(gs_mos)
      TYPE(tddfpt_ground_state_mos), INTENT(inout)       :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_ground_state_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imo

      CALL timeset(routineN, handle)

      IF (ALLOCATED(gs_mos%phases_occ)) &
         DEALLOCATE (gs_mos%phases_occ)

      IF (ALLOCATED(gs_mos%evals_virt)) &
         DEALLOCATE (gs_mos%evals_virt)

      IF (ALLOCATED(gs_mos%evals_occ)) &
         DEALLOCATE (gs_mos%evals_occ)

      IF (ALLOCATED(gs_mos%vectors_mo_virt)) THEN
         DO imo = SIZE(gs_mos%vectors_mo_virt), 1, -1
            IF (ASSOCIATED(gs_mos%vectors_mo_virt(imo)%matrix)) &
               CALL cp_fm_release(gs_mos%vectors_mo_virt(imo)%matrix)
         END DO
         DEALLOCATE (gs_mos%vectors_mo_virt)
      END IF

      IF (ASSOCIATED(gs_mos%matrix_mos_virt)) &
         CALL cp_fm_release(gs_mos%matrix_mos_virt)

      IF (ASSOCIATED(gs_mos%matrix_mos_occ)) &
         CALL cp_fm_release(gs_mos%matrix_mos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_ground_state_mos

! **************************************************************************************************
!> \brief Generate the list of guess vectors.
!> \param guess_vectors     guess vectors (initialised on exit)
!> \param nstates           number of excited states to converge
!> \param nstates_active    the number of excited states in the active space
!> \param gs_mos            occupied and virtual molecular orbitals optimised for the ground state
!> \param is_add_degenerate indicates that all degenerate excited states should be automatically
!>                          added to the active space
!> \param degenerate_eps    energy threshold which controls when excited states
!>                          are considered to be degenerate
!> \param log_unit          output unit
!> \par History
!>    * 05.2016 created as tddfpt_guess() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype, supports spin-polarised density [Sergey Chulkov]
!> \note Based on the subroutine co_initial_guess() which was originally created by
!>       Thomas Chassaing on 06.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_guess_vectors(guess_vectors, nstates, nstates_active, gs_mos, &
                                        is_add_degenerate, degenerate_eps, log_unit)
      TYPE(tddfpt_guess_vectors), INTENT(out)            :: guess_vectors
      INTEGER, INTENT(in)                                :: nstates, nstates_active
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_add_degenerate
      REAL(kind=dp), INTENT(in)                          :: degenerate_eps
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_guess_vectors', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=5)                                   :: spin_label
      INTEGER :: handle, imo_occ, imo_virt, ind, ispin, istate, jstate, nmo_occ_avail_spin, &
         nmo_occ_selected_spin, nmo_virt_avail_spin, nmo_virt_selected_spin, nspins, &
         nstates_occ_virt_alpha, nstates_selected, nstates_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(2)                              :: nmo_occ_avail, nmo_occ_selected, &
                                                            nmo_virt_selected
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)
      IF (debug_this_module) THEN
         CPASSERT(nspins == 1 .OR. nspins == 2)
      END IF

      DO ispin = 1, nspins
         nmo_occ_avail(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      DO ispin = 1, nspins
         nmo_occ_avail_spin = nmo_occ_avail(ispin)
         nmo_virt_avail_spin = SIZE(gs_mos(ispin)%evals_virt)

         DO nmo_occ_selected_spin = nmo_occ_avail_spin-MIN(nstates_active, nmo_occ_avail_spin), 1, -1
            IF (gs_mos(ispin)%evals_occ(nmo_occ_selected_spin+1)- &
                gs_mos(ispin)%evals_occ(nmo_occ_selected_spin) > degenerate_eps) EXIT
         END DO
         nmo_occ_selected_spin = nmo_occ_avail_spin-nmo_occ_selected_spin

         DO nmo_virt_selected_spin = MIN(nstates_active, nmo_virt_avail_spin), nmo_virt_avail_spin-1
            IF (gs_mos(ispin)%evals_virt(nmo_virt_selected_spin+1)- &
                gs_mos(ispin)%evals_virt(nmo_virt_selected_spin) > degenerate_eps) EXIT
         END DO

         nmo_occ_selected(ispin) = nmo_occ_selected_spin
         nmo_virt_selected(ispin) = nmo_virt_selected_spin
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = DOT_PRODUCT(nmo_occ_selected(1:nspins), nmo_virt_selected(1:nspins))

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO ispin = 1, nspins
         nmo_occ_avail_spin = nmo_occ_avail(ispin)
         nmo_occ_selected_spin = nmo_occ_selected(ispin)
         nmo_virt_selected_spin = nmo_virt_selected(ispin)

         DO imo_occ = 1, nmo_occ_selected_spin
            e_occ = gs_mos(ispin)%evals_occ(nmo_occ_avail_spin-imo_occ+1)

            DO imo_virt = 1, nmo_virt_selected_spin
               istate = istate+1
               e_virt_minus_occ(istate) = gs_mos(ispin)%evals_virt(imo_virt)-e_occ
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

         DO istate = 1, nstates
            ind = inds(istate)-1
            IF (nspins > 1) THEN
               IF (ind < nstates_occ_virt_alpha) THEN
                  ispin = 1
                  spin_label = '(alp)'
               ELSE
                  ispin = 2
                  ind = ind-nstates_occ_virt_alpha
                  spin_label = '(bet)'
               END IF
            END IF

            imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
            imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

            WRITE (log_unit, '(1X,I8,5X,I8,1X,A5,3X,I8,1X,A5,2X,F14.5)') &
               istate, imo_occ, spin_label, imo_virt, spin_label, e_virt_minus_occ(istate)*evolt
         END DO
      END IF

      ! detect degenerate states outside the active space
      nstates_total = nstates_active

      DO istate = nstates_total+1, nstates_selected
         IF (e_virt_minus_occ(istate)-e_virt_minus_occ(istate-1) > degenerate_eps) EXIT
      END DO

      IF (istate > nstates_total+1) THEN
         istate = istate-1

         IF (is_add_degenerate) THEN
            nstates_total = istate
         ELSE
            IF (log_unit > 0) THEN
               DO jstate = nstates_total, 2, -1
                  IF (e_virt_minus_occ(jstate)-e_virt_minus_occ(jstate-1) > degenerate_eps) EXIT
               END DO

               IF (jstate <= nstates) &
                  jstate = nstates+1

               IF (jstate <= nstates_total) THEN
                  WRITE (log_unit, '(1X,19("-"),1X,A,1X,20("-"))') "Other active states"

                  DO jstate = jstate, nstates_total
                     ind = inds(jstate)-1
                     IF (nspins > 1) THEN
                        IF (ind < nstates_occ_virt_alpha) THEN
                           ispin = 1
                           spin_label = '(alp)'
                        ELSE
                           ispin = 2
                           ind = ind-nstates_occ_virt_alpha
                           spin_label = '(bet)'
                        END IF
                     END IF

                     imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
                     imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

                     WRITE (log_unit, '(1X,I8,5X,I8,1X,A5,3X,I8,1X,A5,2X,F14.5)') &
                        istate, imo_occ, spin_label, imo_virt, spin_label, e_virt_minus_occ(jstate)*evolt
                  END DO
               END IF

               WRITE (log_unit, '(1X,21("-"),1X,A,1X,22("-"))') "Inactive states"

               DO jstate = nstates_total+1, istate
                  ind = inds(jstate)-1
                  IF (nspins > 1) THEN
                     IF (ind < nstates_occ_virt_alpha) THEN
                        ispin = 1
                        spin_label = '(alp)'
                     ELSE
                        ispin = 2
                        ind = ind-nstates_occ_virt_alpha
                        spin_label = '(bet)'
                     END IF
                  END IF

                  imo_occ = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
                  imo_virt = nmo_occ_avail(ispin)+MOD(ind, nmo_virt_selected(ispin))+1

                  WRITE (log_unit, '(1X,I8,5X,I8,1X,A5,3X,I8,1X,A5,2X,F14.5)') &
                     istate, imo_occ, spin_label, imo_virt, spin_label, e_virt_minus_occ(jstate)*evolt
               END DO
            END IF

            CALL cp_warn(__LOCATION__, "It is advised to include all proposed"// &
                         " degenerate excited states into TDDFT active space")
         END IF
      END IF

      IF (log_unit > 0) &
         WRITE (log_unit, '(/,1X,A,T46,I8)') 'Number of active states:', nstates_total

      ALLOCATE (guess_vectors%imos_occ(nstates_total))
      ALLOCATE (guess_vectors%imos_virt(nstates_total))
      ALLOCATE (guess_vectors%ispins(nstates_total))
      ALLOCATE (guess_vectors%evals(nstates_total))

      DO istate = 1, nstates_total
         ind = inds(istate)-1
         IF (ind < nstates_occ_virt_alpha) THEN
            ispin = 1
         ELSE
            ispin = 2
            ind = ind-nstates_occ_virt_alpha
         END IF

         guess_vectors%ispins(istate) = ispin
         guess_vectors%imos_occ(istate) = nmo_occ_avail(ispin)-ind/nmo_virt_selected(ispin)
         guess_vectors%imos_virt(istate) = MOD(ind, nmo_virt_selected(ispin))+1
      END DO
      guess_vectors%evals(1:nstates_total) = e_virt_minus_occ(1:nstates_total)

      DEALLOCATE (e_virt_minus_occ)
      DEALLOCATE (inds)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_guess_vectors

! **************************************************************************************************
!> \brief Release initial guess vectors.
!> \param guess_vectors   structure that holds guess vectors
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_guess_vectors(guess_vectors)
      TYPE(tddfpt_guess_vectors), INTENT(inout)          :: guess_vectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_guess_vectors', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(guess_vectors%evals)) &
         DEALLOCATE (guess_vectors%evals)

      IF (ALLOCATED(guess_vectors%ispins)) &
         DEALLOCATE (guess_vectors%ispins)

      IF (ALLOCATED(guess_vectors%imos_virt)) &
         DEALLOCATE (guess_vectors%imos_virt)

      IF (ALLOCATED(guess_vectors%imos_occ)) &
         DEALLOCATE (guess_vectors%imos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_guess_vectors

! **************************************************************************************************
!> \brief Make TDDFPT trial vectors orthogonal to all occupied molecular orbitals.
!> \param evects   trial vectors (modified on exit)
!> \param S_C0_C0T matrix product S * C_0 * C_0^T, where C_0 is the ground state wave function
!>            for each spin expressed in atomic basis set, and S is the corresponding overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note  Based on the subroutine p_preortho() which was created by Thomas Chassaing on 09.2002.
!>        Should be useless when the ground states MOs are computed with extremely high accuracy,
!>        as all virtual orbitals are already orthogonal to the occupied ones by design.
!>        It is vitally important, however, when the norm of residual vectors is relatively small
!>        (e.g. less then SCF_EPS). In this case, new krylov vectors seem to be random and should be
!>        orthogonalised even with respect to the occupied MOs.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthogonalize_psi1_psi0(evects, S_C0_C0T)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0_C0T

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthogonalize_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: C0_C0T_S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0)
         CPASSERT(SIZE(S_C0_C0T) == nspins)
      END IF

      IF (nvects > 0) THEN
         ALLOCATE (nmo_occ(nspins))
         ALLOCATE (C0_C0T_S_C1(nspins))

         DO ispin = 1, nspins
            NULLIFY (fm_struct, C0_C0T_S_C1(ispin)%matrix)
            CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, nrow_global=nao, &
                                ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
            CALL cp_fm_create(C0_C0T_S_C1(ispin)%matrix, fm_struct)
         END DO

         DO ivect = 1, nvects
            DO ispin = 1, nspins
               ! C0 * C0^T * S * C1 == (S * C0 * C0^T)^T * C1
               CALL cp_gemm('T', 'N', nao, nmo_occ(ispin), nao, 1.0_dp, S_C0_C0T(ispin)%matrix, &
                            evects(ispin, ivect)%matrix, 0.0_dp, C0_C0T_S_C1(ispin)%matrix)

               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, ivect)%matrix, -1.0_dp, C0_C0T_S_C1(ispin)%matrix)
            END DO
         END DO

         DO ispin = nspins, 1, -1
            CALL cp_fm_release(C0_C0T_S_C1(ispin)%matrix)
         END DO
         DEALLOCATE (C0_C0T_S_C1)
         DEALLOCATE (nmo_occ)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthogonalize_psi1_psi0

! **************************************************************************************************
!> \brief Check that orthogonalised TDDFPT trial vectors remain orthogonal to
!>        occupied molecular orbitals.
!> \param evects   trial vectors (modified on exit)
!> \param S_C0     matrix product S * C_0, where C_0 is the ground state wave function
!>                 for each spin in atomic basis set, and S is the corresponding overlap matrix
!> \param max_norm the largest possible overlap between the ground state and
!>                 excited state wave functions
!> \retval is_nonortho true if trial vectors are non-orthogonal to occupied molecular orbitals
!> \par History
!>    * 07.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_is_nonorthogonal_psi1_psi0(evects, S_C0, max_norm) RESULT(is_nonortho)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0
      REAL(kind=dp), INTENT(in)                          :: max_norm
      LOGICAL                                            :: is_nonortho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_is_nonorthogonal_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: maxabs_val
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: C0T_S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ALLOCATE (nmo_occ(nspins))
      CALL cp_fm_get_info(matrix=S_C0(1)%matrix, nrow_global=nao, ncol_global=nmo_occ(1), context=blacs_env)

      DO ispin = 2, nspins
         CALL cp_fm_get_info(matrix=S_C0(ispin)%matrix, ncol_global=nmo_occ(ispin))
      END DO

      NULLIFY (fm_struct)
      ALLOCATE (C0T_S_C1(nspins))
      DO ispin = 1, nspins
         NULLIFY (C0T_S_C1(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL cp_fm_create(C0T_S_C1(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      is_nonortho = .FALSE.

      loop: DO ivect = 1, nvects
         DO ispin = 1, nspins
            CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, 1.0_dp, S_C0(ispin)%matrix, &
                         evects(ispin, ivect)%matrix, 0.0_dp, C0T_S_C1(ispin)%matrix)

            CALL cp_fm_maxabsval(C0T_S_C1(ispin)%matrix, maxabs_val)
            is_nonortho = maxabs_val > max_norm
            IF (is_nonortho) EXIT loop
         END DO
      END DO loop

      DO ispin = SIZE(C0T_S_C1), 1, -1
         CALL cp_fm_release(C0T_S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (C0T_S_C1)

      DEALLOCATE (nmo_occ)
      CALL timestop(handle)
   END FUNCTION tddfpt_is_nonorthogonal_psi1_psi0

! **************************************************************************************************
!> \brief Make new TDDFPT trial vectors orthonormal to all previous trial vectors.
!> \param evects      trial vectors (modified on exit)
!> \param nvects_new  number of new trial vectors to orthogonalise
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!> \note Based on the subroutines reorthogonalize() and normalize() which were originally created
!>             by Thomas Chassaing on 03.2003.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthonormalize_psi1_psi1(evects, nvects_new, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      INTEGER, INTENT(in)                                :: nvects_new
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthonormalize_psi1_psi1', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER                           :: threshold = 16*EPSILON(1.0_dp)

      INTEGER                                            :: handle, ispin, ivect, jvect, nspins, &
                                                            nvects_old, nvects_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: norm, norm_spin, weight
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects_total = SIZE(evects, 2)
      nvects_old = nvects_total-nvects_new
      norm_spin = 1.0_dp

      ALLOCATE (nmo_occ(nspins))
      ALLOCATE (S_C1(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, S_C1(ispin)%matrix)
         CALL cp_fm_get_info(evects(ispin, 1)%matrix, ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
         CALL cp_fm_create(S_C1(ispin)%matrix, fm_struct)
      END DO

      ! ensure <psi_{i,alpha} | psi_{j,alpha}> == <psi_{i,beta} | psi_{j,beta}> == 0, for i /= j ,
      ! which leads to <psi_i | psi_j> = 0 for i /= j , where psi_i = psi_{i,alpha} + psi_{i,beta}.
      ! The orthogonality condition <psi_i | psi_j> = delta_{i,j} does not imply that
      ! the spin-components psi_{i,alpha} and psi_{i,beta} are normalised by themselves
      ! (<psi_{i,alpha} | psi_{i,alpha}> in [0..1]; <psi_{i,beta} | psi_{i,beta}> in [0..1]),
      ! only their sum does (<psi_{i,alpha} | psi_{i,alpha}> + <psi_{i,beta} | psi_{i,beta}> == 1)
      DO jvect = nvects_old+1, nvects_total
         DO ivect = 1, jvect-1
            DO ispin = 1, nspins
               ! TO DO: do not recompute S * C1 for old trial wectors
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               IF (nspins > 1) &
                  CALL cp_fm_trace(evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, norm_spin)

               IF (norm_spin > threshold) THEN
                  CALL cp_fm_trace(evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, weight)

                  CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, jvect)%matrix, -weight/norm_spin, evects(ispin, ivect)%matrix)
               END IF
            END DO
         END DO

         norm = 0.0_dp
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, &
                                         ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
            CALL cp_fm_trace(evects(ispin, jvect)%matrix, S_C1(ispin)%matrix, weight)
            norm = norm+weight
         END DO

         norm = 1.0_dp/SQRT(norm)
         DO ispin = 1, nspins
            CALL cp_fm_scale(norm, evects(ispin, ivect)%matrix)
         END DO
      END DO

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (S_C1)
      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************
!> \brief Initialise hat{A} C_{1,i} by energy difference using the following expression:
!>        KS_{sigma} * C_{1,i,sigma} - S * C_{1,i,sigma} * diag(evals_occ) .
!> \param Aop_evects  hat{A} C_{1,i} (initialised on exit)
!> \param evects      trial vectors C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied orbital
!>                    energies [field %evals_occ] are needed)
!> \param matrix_ks   Kohn-Sham matrix
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 initialise all matrix elements in one go [Sergey Chulkov]
!> \note Based on the subroutine p_op_l1() which was originally created by
!>       Thomas Chassaing on 08.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_init_by_energy_diff(Aop_evects, evects, gs_mos, matrix_ks, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_dbcsr_p_type), DIMENSION(:), INTENT(in)    :: matrix_ks
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_by_energy_diff', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ALLOCATE (nmo_occ(nspins))
      ALLOCATE (S_C1(nspins))

      DO ispin = 1, nspins
         NULLIFY (S_C1(ispin)%matrix)
      END DO

      DO ivect = 1, nvects
         DO ispin = 1, nspins
            IF (ASSOCIATED(evects(ispin, ivect)%matrix)) THEN
               IF (.NOT. ASSOCIATED(S_C1(ispin)%matrix)) THEN
                  CALL cp_fm_get_info(evects(ispin, ivect)%matrix, ncol_global=nmo_occ(ispin), matrix_struct=fm_struct)
                  CALL cp_fm_create(S_C1(ispin)%matrix, fm_struct)
               END IF

               CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, evects(ispin, ivect)%matrix, &
                                            Aop_evects(ispin, ivect)%matrix, ncol=nmo_occ(ispin), &
                                            alpha=1.0_dp, beta=0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, ivect)%matrix, S_C1(ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

               CALL cp_fm_column_scale(S_C1(ispin)%matrix, gs_mos(ispin)%evals_occ)
               ! KS * C1 - S * C1 * occupied_orbital_energies
               CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect)%matrix, -1.0_dp, S_C1(ispin)%matrix)
            END IF
         END DO
      END DO

      DO ispin = nspins, 1, -1
         IF (ASSOCIATED(S_C1(ispin)%matrix)) &
            CALL cp_fm_release(S_C1(ispin)%matrix)
      END DO
      DEALLOCATE (S_C1)
      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_init_by_energy_diff

! **************************************************************************************************
!> \brief Update the hat{A} C_{1,i} vectors by adding coulomb terms.
!> \param Aop_evects  hat{A} C_{1,i} (updated on exit)
!> \param evects      trial vectors C_{1,i}
!> \param gs_mos      molecular orbitals optimised for the ground state (only occupied molecular
!>                    orbitals [field %matrix_mos_occ] are needed)
!> \param matrix_s    overlap matrix
!> \param qs_env      Quickstep environment
!> \param blacs_env   BLACS parallel environment
!> \param pw_env      plain wave environment
!> \param task_list   task list
!> \par History
!>    * 05.2016 compute all coulomb terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_coulomb(Aop_evects, evects, gs_mos, matrix_s, qs_env, &
                                   blacs_env, pw_env, task_list)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(task_list_type), POINTER                      :: task_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: alpha, pair_energy
      TYPE(cp_dbcsr_p_type)                              :: J_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_p_type)                                    :: v_gspace, v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      ! *** initialize variables ***
      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (debug_this_module) THEN
         CPASSERT(dft_control%nspins == nspins)
      END IF
      ! a sum J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu can be computed by solving
      ! the Poisson equation for combined density (rho_{ia,alpha} + rho_{ia,beta}) .
      ! Temporary activate spin-restricted mode for acceleration purpose, as we actually
      ! do not need a separate beta-spin component
      dft_control%nspins = 1

      NULLIFY (auxbas_pw_pool, poisson_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      ! alpha-spin electron density
      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, 1)
      CALL cp_dbcsr_init_p(rho_ia_ao(1)%matrix)
      CALL cp_dbcsr_copy(rho_ia_ao(1)%matrix, matrix_s)

      NULLIFY (rho_ia_struct)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE., pw_env_external=pw_env)
      CALL qs_rho_get(rho_ia_struct, rho_g=rho_ia_g)

      NULLIFY (J_ia_munu%matrix)
      CALL cp_dbcsr_init_p(J_ia_munu%matrix)
      CALL cp_dbcsr_copy(J_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (v_gspace%pw, v_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_gspace%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace%pw, use_data=REALDATA3D, in_space=REALSPACE)

      IF (nspins > 1) THEN
         alpha = 1.0_dp
      ELSE
         ! spin-restricted case;
         ! alpha == 2 due to singlet state, for triplet states we should not call this subroutine at all
         alpha = 2.0_dp
      END IF

      ! *** compute coulomb terms ***
      DO ivect = 1, nvects
         IF (ASSOCIATED(evects(1, ivect)%matrix)) THEN
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(1), 0.5_dp, &
                         gs_mos(1)%matrix_mos_occ, evects(1, ivect)%matrix, 0.0_dp, rho_ia_fm)
            CALL cp_gemm('N', 'T', nao, nao, nmo_occ(1), 0.5_dp, &
                         evects(1, ivect)%matrix, gs_mos(1)%matrix_mos_occ, 1.0_dp, rho_ia_fm)

            DO ispin = 2, nspins
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                            gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 1.0_dp, rho_ia_fm)
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                            evects(ispin, ivect)%matrix, gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_fm)
            END DO

            CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(1)%matrix, keep_sparsity=.TRUE.)
            CALL qs_rho_update_rho(rho_ia_struct, qs_env, pw_env_external=pw_env, task_list_external=task_list)

            CALL cp_dbcsr_set(J_ia_munu%matrix, 0.0_dp)

            CALL pw_poisson_solve(poisson_env, rho_ia_g(1)%pw, pair_energy, v_gspace%pw)
            CALL pw_transfer(v_gspace%pw, v_rspace%pw)
            CALL pw_scale(v_rspace%pw, v_rspace%pw%pw_grid%dvol)

            CALL integrate_v_rspace(v_rspace=v_rspace, hmat=J_ia_munu, qs_env=qs_env, &
                                    calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                                    pw_env_external=pw_env, task_list_external=task_list)

            ! (i a || j b) = ( i_alpha a_alpha + i_beta a_beta || j_alpha b_alpha + j_beta b_beta) =
            !                tr (Cj_alpha^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_alpha) +
            !                tr (Cj_beta^T * [J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu] * Cb_beta)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(J_ia_munu%matrix, gs_mos(ispin)%matrix_mos_occ, Aop_evects(ispin, ivect)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=alpha, beta=1.0_dp)
            END DO
         END IF
      END DO

      dft_control%nspins = nspins

      ! *** clean up ***
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_gspace%pw)

      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(J_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_coulomb

! **************************************************************************************************
!> \brief Update action matrix-vector products by adding kernel terms.
!> \param Aop_evects      hat{A} C_{1,i} (updated on exit)
!> \param evects          TDDFPT trial vectors
!> \param gs_mos          molecular orbitals optimised for the ground state (only occupied molecular
!>                        orbitals [field %matrix_mos_occ] are needed)
!> \param is_rks_triplets indicates that the triplet excited states calculation using
!>                        spin-unpolarised molecular orbitals has been requested
!> \param matrix_s        overlap matrix
!> \param qs_env          Quickstep environment
!> \param xc_rho_set      a variable to compute the 2nd derivative of XC-functional
!> \param xc_deriv_set    a variable to compute the 2nd derivative of XC-functional
!> \param xc_section      XC input section
!> \param blacs_env       BLACS parallel environment
!> \param pw_env          plain wave environment
!> \param task_list       task list
!> \par History
!>    * 05.2016 compute all kernel terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_kernel(Aop_evects, evects, gs_mos, is_rks_triplets, matrix_s, qs_env, &
                                  xc_rho_set, xc_deriv_set, xc_section, blacs_env, pw_env, task_list)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(task_list_type), POINTER                      :: task_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, &
                                                            nvects, xc_deriv_method_id, &
                                                            xc_rho_smooth_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      LOGICAL                                            :: lsd
      REAL(kind=dp)                                      :: alpha, tddfpt_fac
      TYPE(cp_dbcsr_p_type)                              :: F_ia_munu
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ia_ao
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: rho_ia_fm
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g2, rho_ia_r, &
                                                            rho_ia_r2, tau_ia_r, tau_ia_r2, v_xc
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_ia_struct
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      lsd = (nspins > 1) .OR. is_rks_triplets

      ! *** initialize variables ***
      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      NULLIFY (rho_ia_ao)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(rho_ia_ao(ispin)%matrix)
         CALL cp_dbcsr_copy(rho_ia_ao(ispin)%matrix, matrix_s)
      END DO

      NULLIFY (rho_ia_struct, rho_ia_r, rho_ia_g)
      CALL qs_rho_create(rho_ia_struct)
      CALL qs_rho_set(rho_ia_struct, rho_ao=rho_ia_ao)

      CALL qs_rho_rebuild(rho_ia_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE., pw_env_external=pw_env)
      CALL qs_rho_get(rho_ia_struct, rho_r=rho_ia_r, rho_g=rho_ia_g, tau_r=tau_ia_r)

      NULLIFY (F_ia_munu%matrix)
      CALL cp_dbcsr_init_p(F_ia_munu%matrix)
      CALL cp_dbcsr_copy(F_ia_munu%matrix, matrix_s)

      NULLIFY (fm_struct, rho_ia_fm)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      CALL cp_fm_create(rho_ia_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

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
         ELSE
            NULLIFY (tau_ia_r2)
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
         ELSE
            NULLIFY (tau_ia_r2)
         END IF
      END IF

      ALLOCATE (v_xc(nspins))
      DO ispin = 1, SIZE(v_xc)
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, in_space=REALSPACE)
      END DO

      NULLIFY (xc_rho1_set)
      CALL xc_rho_set_create(xc_rho1_set, rho_ia_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      xc_deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      xc_rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      ! get arguments needed by the given functional
      needs = xc_functionals_get_needs(functionals=xc_fun_section, lsd=lsd, add_basic_components=.TRUE.)

      alpha = 1.0_dp
      tddfpt_fac = 0.0_dp
      ! tddfpt_fac is taken into account in spin-restricted case only
      IF (nspins == 1) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            tddfpt_fac = -1.0_dp
         ELSE
            !                                                 alpha              tddfpt_fac
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            alpha = 2.0_dp
         END IF
      END IF

      ! *** compute kernel terms ***
      DO ivect = 1, nvects
         IF (ASSOCIATED(evects(1, ivect)%matrix)) THEN
            DO ispin = 1, nspins
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                            gs_mos(ispin)%matrix_mos_occ, evects(ispin, ivect)%matrix, 0.0_dp, rho_ia_fm)
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, &
                            evects(ispin, ivect)%matrix, gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_fm)

               CALL copy_fm_to_dbcsr(rho_ia_fm, rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
            END DO
            CALL qs_rho_update_rho(rho_ia_struct, qs_env, pw_env_external=pw_env, task_list_external=task_list)

            CALL xc_rho_set_update(rho_set=xc_rho1_set, rho_r=rho_ia_r2, rho_g=rho_ia_g2, &
                                   tau=tau_ia_r2, needs=needs, xc_deriv_method_id=xc_deriv_method_id, &
                                   xc_rho_smooth_id=xc_rho_smooth_id, pw_pool=auxbas_pw_pool)

            DO ispin = 1, SIZE(v_xc)
               CALL pw_zero(v_xc(ispin)%pw)
            END DO

            CALL xc_calc_2nd_deriv(v_xc=v_xc, deriv_set=xc_deriv_set, rho_set=xc_rho_set, &
                                   rho1_set=xc_rho1_set, pw_pool=auxbas_pw_pool, &
                                   xc_section=xc_section, gapw=.FALSE., tddfpt_fac=tddfpt_fac)

            DO ispin = 1, nspins
               CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)

               CALL cp_dbcsr_set(F_ia_munu%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=v_xc(ispin), hmat=F_ia_munu, qs_env=qs_env, &
                                       calculate_forces=.FALSE., gapw=.FALSE., &
                                       pw_env_external=pw_env, task_list_external=task_list)

               CALL cp_dbcsr_sm_fm_multiply(F_ia_munu%matrix, gs_mos(ispin)%matrix_mos_occ, Aop_evects(ispin, ivect)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=alpha, beta=1.0_dp)
            END DO
         END IF
      END DO

      ! *** clean up ***
      CALL xc_rho_set_release(xc_rho1_set)

      DO ispin = SIZE(v_xc), 1, -1
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (v_xc)

      IF (ASSOCIATED(tau_ia_r2)) &
         DEALLOCATE (tau_ia_r2)
      DEALLOCATE (rho_ia_r2, rho_ia_g2)

      CALL cp_fm_release(rho_ia_fm)

      CALL cp_dbcsr_deallocate_matrix(F_ia_munu%matrix)
      CALL qs_rho_release(rho_ia_struct)

      DEALLOCATE (nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_kernel

! **************************************************************************************************
!> \brief Update action matrix-vector products by adding exact-exchange terms.
!> \param Aop_evects        hat{A} C_{1,i} (updated on exit)
!> \param evects            trial vectors C_{1,i}
!> \param gs_mos            molecular orbitals optimised for the ground state (only occupied
!>                          molecular orbitals [field %matrix_mos_occ] are needed)
!> \param qs_env            Quickstep environment
!> \param do_admm           perform auxiliary density matrix method calculations;
!>                          all the subsequent arguments are ignored if do_admm == .false.
!> \param blacs_env         BLACS parallel environment
!> \param matrix_s          overlap matrix for the primary basis set
!> \param matrix_s_aux_fit  overlap matrix for the auxiliary basis set
!> \param pw_env            plain wave environment
!> \param task_list_orb     task list for the primary basis set
!> \param task_list_aux_fit task list for the auxiliary basis set
!> \param nsubgroups        number of MPI subgroups
!> \param subgroup          list of sub-group members
!> \par History
!>    * 05.2016 compute all exact-exchange terms in one go [Sergey Chulkov]
!> \note Based on the subroutine kpp1_calc_k_p_p1() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_hfx(Aop_evects, evects, gs_mos, qs_env, do_admm, blacs_env, &
                               matrix_s, matrix_s_aux_fit, pw_env, task_list_orb, &
                               task_list_aux_fit, nsubgroups, subgroup)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: do_admm
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s, matrix_s_aux_fit
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(task_list_type), POINTER                      :: task_list_orb, task_list_aux_fit
      INTEGER, INTENT(in)                                :: nsubgroups
      INTEGER, DIMENSION(0:), INTENT(in)                 :: subgroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_apply_hfx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate, nao, nao_aux, &
                                                            nspins, nstates
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: alpha, exc_aux_fit
      REAL(kind=dp), DIMENSION(:), POINTER               :: tot_rho_aux_fit_r
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: hmat_aux_fit, hmat_global, hmat_orb, &
         matrix_s_aux_fit_global, matrix_s_global, rho_ia_aux_fit, rho_ia_global, rho_ia_orb
      TYPE(cp_dbcsr_type), POINTER                       :: hmat_nosym_global
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: rho_ia_aux_fit_fm, rho_ia_orb_fm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: admm_A_sub, rho_ia_fm_global, &
                                                            work_aux_orb_sub
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_aux_fit_g, rho_aux_fit_r, vxc_rho, &
                                                            vxc_tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit_struct, rho_orb_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)

      IF (nspins > 1) THEN
         alpha = 2.0_dp
      ELSE
         alpha = 4.0_dp
      END IF

      ! *** initialize variables ***
      NULLIFY (admm_env, auxbas_pw_pool, blacs_env_global, ks_env, matrix_s_global, matrix_s_aux_fit_global)
      CALL get_qs_env(qs_env, admm_env=admm_env, blacs_env=blacs_env_global, matrix_s=matrix_s_global)

      IF (do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env, ks_env=ks_env, matrix_s_aux_fit=matrix_s_aux_fit_global)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux)
      END IF

      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao)
      ALLOCATE (nmo_occ(nspins))
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)

      ALLOCATE (rho_ia_orb_fm(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (rho_ia_orb_fm(ispin, istate)%matrix)
            IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
               CALL cp_fm_create(rho_ia_orb_fm(ispin, istate)%matrix, fm_struct)
         END DO
      END DO

      CALL cp_fm_struct_release(fm_struct)

      IF (do_admm) THEN
         NULLIFY (hmat_orb, hmat_aux_fit, rho_ia_orb, rho_ia_aux_fit, rho_orb_struct, rho_aux_fit_struct)

         CALL cp_dbcsr_allocate_matrix_set(hmat_orb, nspins)
         DO ispin = 1, nspins
            CALL cp_dbcsr_init_p(hmat_orb(ispin)%matrix)
            CALL cp_dbcsr_copy(hmat_orb(ispin)%matrix, matrix_s)
         END DO

         CALL cp_dbcsr_allocate_matrix_set(hmat_aux_fit, nspins)
         DO ispin = 1, nspins
            CALL cp_dbcsr_init_p(hmat_aux_fit(ispin)%matrix)
            CALL cp_dbcsr_copy(hmat_aux_fit(ispin)%matrix, matrix_s_aux_fit)
         END DO

         CALL cp_dbcsr_allocate_matrix_set(rho_ia_orb, nspins)
         DO ispin = 1, nspins
            CALL cp_dbcsr_init_p(rho_ia_orb(ispin)%matrix)
            CALL cp_dbcsr_copy(rho_ia_orb(ispin)%matrix, matrix_s)
         END DO

         CALL cp_dbcsr_allocate_matrix_set(rho_ia_aux_fit, nspins)
         DO ispin = 1, nspins
            CALL cp_dbcsr_init_p(rho_ia_aux_fit(ispin)%matrix)
            CALL cp_dbcsr_copy(rho_ia_aux_fit(ispin)%matrix, matrix_s_aux_fit)
         END DO

         CALL qs_rho_create(rho_orb_struct)
         CALL qs_rho_set(rho_orb_struct, rho_ao=rho_ia_orb)

         CALL qs_rho_create(rho_aux_fit_struct)
         CALL qs_rho_set(rho_aux_fit_struct, rho_ao=rho_ia_aux_fit)

         ! grid representation of the electron density is needed to compute ADMM correction
         CALL qs_rho_rebuild(rho_orb_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE., pw_env_external=pw_env)
         CALL qs_rho_rebuild(rho_aux_fit_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE., pw_env_external=pw_env)
         CALL qs_rho_get(rho_aux_fit_struct, rho_g=rho_aux_fit_g, rho_r=rho_aux_fit_r, tot_rho_r=tot_rho_aux_fit_r)

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao_aux, context=blacs_env)

         ALLOCATE (rho_ia_aux_fit_fm(nspins, nstates))
         DO istate = 1, nstates
            DO ispin = 1, nspins
               NULLIFY (rho_ia_aux_fit_fm(ispin, istate)%matrix)
               IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
                  CALL cp_fm_create(rho_ia_aux_fit_fm(ispin, istate)%matrix, fm_struct)
            END DO
         END DO

         CALL cp_fm_struct_release(fm_struct)

         NULLIFY (admm_A_sub, work_aux_orb_sub)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(admm_A_sub, fm_struct)
         CALL cp_fm_create(work_aux_orb_sub, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL tddfpt_fm_replicate_across_subgroups(fm_src=admm_env%A, fm_dest_sub=admm_A_sub, &
                                                   nsubgroups=nsubgroups, subgroup=subgroup)
      END IF

      ! compute the response electron density for all excited state
      DO istate = 1, nstates
         DO ispin = 1, nspins
            IF (ASSOCIATED(evects(ispin, istate)%matrix)) THEN
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, gs_mos(ispin)%matrix_mos_occ, &
                            evects(ispin, istate)%matrix, 0.0_dp, rho_ia_orb_fm(ispin, istate)%matrix)
               CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 0.5_dp, evects(ispin, istate)%matrix, &
                            gs_mos(ispin)%matrix_mos_occ, 1.0_dp, rho_ia_orb_fm(ispin, istate)%matrix)
            END IF
         END DO
      END DO

      IF (do_admm) THEN
         DO istate = 1, nstates
            DO ispin = 1, nspins
               IF (ASSOCIATED(rho_ia_orb_fm(ispin, istate)%matrix)) THEN
                  CALL cp_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, admm_A_sub, &
                               rho_ia_orb_fm(ispin, istate)%matrix, 0.0_dp, work_aux_orb_sub)
                  CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_A_sub, work_aux_orb_sub, &
                               0.0_dp, rho_ia_aux_fit_fm(ispin, istate)%matrix)
               END IF
            END DO
         END DO

         NULLIFY (vxc_rho, vxc_tau)
         DO istate = 1, nstates
            IF (ASSOCIATED(evects(1, istate)%matrix)) THEN

               DO ispin = 1, nspins
                  CALL copy_fm_to_dbcsr(rho_ia_orb_fm(ispin, istate)%matrix, rho_ia_orb(ispin)%matrix, keep_sparsity=.TRUE.)
               END DO
               CALL qs_rho_update_rho(rho_orb_struct, qs_env, pw_env_external=pw_env, task_list_external=task_list_orb)

               DO ispin = 1, nspins
                  CALL copy_fm_to_dbcsr(rho_ia_aux_fit_fm(ispin, istate)%matrix, rho_ia_aux_fit(ispin)%matrix, &
                                        keep_sparsity=.TRUE.)
                  CALL calculate_rho_elec(matrix_p=rho_ia_aux_fit(ispin)%matrix, &
                                          rho=rho_aux_fit_r(ispin), rho_gspace=rho_aux_fit_g(ispin), &
                                          total_rho=tot_rho_aux_fit_r(ispin), ks_env=ks_env, &
                                          soft_valid=.FALSE., basis_type="AUX_FIT", &
                                          pw_env_external=pw_env, task_list_external=task_list_aux_fit)
               END DO

               ! some stuff from qs_ks_build_kohn_sham_matrix
               ! TO DO: add SIC support

               ! E_{x}^{DFT}[\rho]
               CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_orb_struct, xc_section=admm_env%xc_section_aux, &
                                  vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=exc_aux_fit, just_energy=.FALSE., pw_env_external=pw_env)

               IF (ASSOCIATED(vxc_rho)) THEN
                  DO ispin = 1, nspins
                     CALL pw_scale(vxc_rho(ispin)%pw, vxc_rho(ispin)%pw%pw_grid%dvol)

                     CALL cp_dbcsr_set(hmat_orb(ispin)%matrix, 0.0_dp)
                     CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), pmat=rho_ia_orb(ispin), hmat=hmat_orb(ispin), &
                                             qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE., &
                                             pw_env_external=pw_env, task_list_external=task_list_orb)

                     CALL cp_dbcsr_sm_fm_multiply(hmat_orb(ispin)%matrix, gs_mos(ispin)%matrix_mos_occ, &
                                                  Aop_evects(ispin, istate)%matrix, ncol=nmo_occ(ispin), &
                                                  alpha=-1.0_dp, beta=1.0_dp)
                     CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
                  END DO
                  DEALLOCATE (vxc_rho)
               END IF

               IF (ASSOCIATED(vxc_tau)) THEN
                  DO ispin = 1, nspins
                     CALL pw_scale(vxc_tau(ispin)%pw, vxc_tau(ispin)%pw%pw_grid%dvol)

                     CALL cp_dbcsr_set(hmat_orb(ispin)%matrix, 0.0_dp)
                     CALL integrate_v_rspace(v_rspace=vxc_tau(ispin), pmat=rho_ia_orb(ispin), hmat=hmat_orb(ispin), &
                                             qs_env=qs_env, calculate_forces=.FALSE., compute_tau=.TRUE., gapw=.FALSE., &
                                             pw_env_external=pw_env, task_list_external=task_list_orb)

                     CALL cp_dbcsr_sm_fm_multiply(hmat_orb(ispin)%matrix, gs_mos(ispin)%matrix_mos_occ, &
                                                  Aop_evects(ispin, istate)%matrix, ncol=nmo_occ(ispin), &
                                                  alpha=-1.0_dp, beta=1.0_dp)
                     CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_tau(ispin)%pw)
                  END DO
                  DEALLOCATE (vxc_tau)
               END IF

               ! E_{x}^{DFT}[\hat{\rho}]
               CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_aux_fit_struct, xc_section=admm_env%xc_section_aux, &
                                  vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=exc_aux_fit, just_energy=.FALSE., pw_env_external=pw_env)

               IF (ASSOCIATED(vxc_rho)) THEN
                  DO ispin = 1, nspins
                     CALL pw_scale(vxc_rho(ispin)%pw, vxc_rho(ispin)%pw%pw_grid%dvol)

                     CALL cp_dbcsr_set(hmat_aux_fit(ispin)%matrix, 0.0_dp)
                     CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), pmat=rho_ia_aux_fit(ispin), &
                                             hmat=hmat_aux_fit(ispin), qs_env=qs_env, calculate_forces=.FALSE., &
                                             force_adm=.TRUE., ispin=ispin, gapw=.FALSE., basis_type="AUX_FIT", &
                                             pw_env_external=pw_env, task_list_external=task_list_aux_fit)

                     CALL cp_dbcsr_sm_fm_multiply(hmat_aux_fit(ispin)%matrix, admm_A_sub, work_aux_orb_sub, &
                                                  ncol=nao, alpha=1.0_dp, beta=0.0_dp)
                     CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_A_sub, work_aux_orb_sub, &
                                  0.0_dp, rho_ia_orb_fm(ispin, istate)%matrix)
                     CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, 1.0_dp, rho_ia_orb_fm(ispin, istate)%matrix, &
                                  gs_mos(ispin)%matrix_mos_occ, 1.0_dp, Aop_evects(ispin, istate)%matrix)

                     CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
                  END DO
                  DEALLOCATE (vxc_rho)
               END IF

               IF (ASSOCIATED(vxc_tau)) THEN
                  DO ispin = 1, nspins
                     CALL pw_scale(vxc_tau(ispin)%pw, vxc_tau(ispin)%pw%pw_grid%dvol)

                     CALL cp_dbcsr_set(hmat_aux_fit(ispin)%matrix, 0.0_dp)
                     CALL integrate_v_rspace(v_rspace=vxc_tau(ispin), pmat=rho_ia_aux_fit(ispin), &
                                             hmat=hmat_aux_fit(ispin), qs_env=qs_env, calculate_forces=.FALSE., &
                                             force_adm=.TRUE., ispin=ispin, compute_tau=.TRUE., gapw=.FALSE., &
                                             basis_type="AUX_FIT", pw_env_external=pw_env, task_list_external=task_list_aux_fit)

                     CALL cp_dbcsr_sm_fm_multiply(hmat_aux_fit(ispin)%matrix, admm_A_sub, work_aux_orb_sub, &
                                                  ncol=nao, alpha=1.0_dp, beta=0.0_dp)
                     CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_A_sub, work_aux_orb_sub, &
                                  0.0_dp, rho_ia_orb_fm(ispin, istate)%matrix)
                     CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, 1.0_dp, rho_ia_orb_fm(ispin, istate)%matrix, &
                                  gs_mos(ispin)%matrix_mos_occ, 1.0_dp, Aop_evects(ispin, istate)%matrix)

                     CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_tau(ispin)%pw)
                  END DO
                  DEALLOCATE (vxc_tau)
               END IF
            END IF
         END DO

         CALL qs_rho_release(rho_aux_fit_struct)
         CALL qs_rho_release(rho_orb_struct)
         CALL cp_dbcsr_deallocate_matrix_set(hmat_aux_fit)
         CALL cp_dbcsr_deallocate_matrix_set(hmat_orb)
      END IF

      ! E_{x}^{HFX}[\rho] and E_{x}^{HFX}[\hat{\rho}]
      ! create globally distributed sparse matrices to store electron density as well as
      ! exact-exchange contributions to the Hamiltonian
      NULLIFY (hmat_global, rho_ia_global)
      CALL cp_dbcsr_allocate_matrix_set(hmat_global, nspins)
      CALL cp_dbcsr_allocate_matrix_set(rho_ia_global, nspins)
      DO ispin = 1, nspins
         CALL cp_dbcsr_init_p(hmat_global(ispin)%matrix)
         CALL cp_dbcsr_init_p(rho_ia_global(ispin)%matrix)
         IF (do_admm) THEN
            CALL cp_dbcsr_copy(hmat_global(ispin)%matrix, matrix_s_aux_fit_global(1)%matrix)
            CALL cp_dbcsr_copy(rho_ia_global(ispin)%matrix, matrix_s_aux_fit_global(1)%matrix)
         ELSE
            CALL cp_dbcsr_copy(hmat_global(ispin)%matrix, matrix_s_global(1)%matrix)
            CALL cp_dbcsr_copy(rho_ia_global(ispin)%matrix, matrix_s_global(1)%matrix)
         END IF
      END DO

      NULLIFY (hmat_nosym_global)
      CALL cp_dbcsr_init_p(hmat_nosym_global)
      CALL cp_dbcsr_create(hmat_nosym_global, template=hmat_global(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      IF (do_admm) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao_aux, context=blacs_env_global)
      ELSE
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_global)
      END IF
      CALL cp_fm_create(rho_ia_fm_global, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      DO istate = 1, nstates
         DO ispin = 1, nspins
            IF (do_admm) THEN
               CALL cp_fm_copy_general(source=rho_ia_aux_fit_fm(ispin, istate)%matrix, &
                                       destination=rho_ia_fm_global, &
                                       global_context=blacs_env_global%group)
            ELSE
               CALL cp_fm_copy_general(source=rho_ia_orb_fm(ispin, istate)%matrix, &
                                       destination=rho_ia_fm_global, &
                                       global_context=blacs_env_global%group)
            END IF

            CALL copy_fm_to_dbcsr(rho_ia_fm_global, rho_ia_global(ispin)%matrix, keep_sparsity=.TRUE.)
            CALL cp_dbcsr_set(hmat_global(ispin)%matrix, 0.0_dp)
         END DO

         CALL tddft_hfx_matrix(hmat_global, rho_ia_global, qs_env)

         DO ispin = 1, nspins
            CALL cp_dbcsr_desymmetrize(hmat_global(ispin)%matrix, hmat_nosym_global)
            CALL copy_dbcsr_to_fm(hmat_nosym_global, rho_ia_fm_global)

            IF (do_admm) THEN
               CALL cp_fm_copy_general(source=rho_ia_fm_global, &
                                       destination=rho_ia_aux_fit_fm(ispin, istate)%matrix, &
                                       global_context=blacs_env_global%group)
            ELSE
               CALL cp_fm_copy_general(source=rho_ia_fm_global, &
                                       destination=rho_ia_orb_fm(ispin, istate)%matrix, &
                                       global_context=blacs_env_global%group)
            END IF
         END DO
      END DO

      CALL cp_fm_release(rho_ia_fm_global)
      CALL cp_dbcsr_deallocate_matrix_set(rho_ia_global)
      CALL cp_dbcsr_deallocate_matrix(hmat_nosym_global)
      CALL cp_dbcsr_deallocate_matrix_set(hmat_global)

      IF (do_admm) THEN
         DO istate = 1, nstates
            DO ispin = 1, nspins
               IF (ASSOCIATED(evects(ispin, istate)%matrix)) THEN
                  CALL cp_gemm('N', 'N', nao_aux, nao, nao_aux, 1.0_dp, rho_ia_aux_fit_fm(ispin, istate)%matrix, &
                               admm_A_sub, 0.0_dp, work_aux_orb_sub)
                  CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_A_sub, work_aux_orb_sub, &
                               0.0_dp, rho_ia_orb_fm(ispin, istate)%matrix)
               END IF
            END DO
         END DO
      END IF

      DO istate = 1, nstates
         DO ispin = 1, nspins
            IF (ASSOCIATED(evects(ispin, istate)%matrix)) &
               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, alpha, rho_ia_orb_fm(ispin, istate)%matrix, &
                            gs_mos(ispin)%matrix_mos_occ, 1.0_dp, Aop_evects(ispin, istate)%matrix)
         END DO
      END DO

      IF (do_admm) THEN
         CALL cp_fm_release(work_aux_orb_sub)
         CALL cp_fm_release(admm_A_sub)

         DO istate = nstates, 1, -1
            DO ispin = nspins, 1, -1
               IF (ASSOCIATED(rho_ia_aux_fit_fm(ispin, istate)%matrix)) &
                  CALL cp_fm_release(rho_ia_aux_fit_fm(ispin, istate)%matrix)
            END DO
         END DO
         DEALLOCATE (rho_ia_aux_fit_fm)
      END IF

      DO istate = nstates, 1, -1
         DO ispin = nspins, 1, -1
            IF (ASSOCIATED(rho_ia_orb_fm(ispin, istate)%matrix)) &
               CALL cp_fm_release(rho_ia_orb_fm(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (rho_ia_orb_fm)

      DEALLOCATE (nmo_occ)
      CALL timestop(handle)
   END SUBROUTINE tddfpt_apply_hfx

! **************************************************************************************************
!> \brief Compute action matrix-vector products.
!> \param Aop_evects            hat{A} C_{1,i} (allocated and initialised on exit)
!> \param evects                TDDFPT trial vectors
!> \param fm_struct_evects      structures of newly allocated dense matrices
!> \param gs_mos                molecular orbitals optimised for the ground state
!> \param is_rks_triplets       indicates that a triplet excited states calculation using
!>                              spin-unpolarised molecular orbitals has been requested
!> \param matrix_ks_sub         Kohn-Sham matrix
!> \param matrix_s_sub          overlap matrix for the primary basis set
!> \param qs_env                Quickstep environment
!> \param xc_rho_set            a variable to compute the 2nd derivative of XC-functional
!> \param xc_deriv_set          a variable to compute the 2nd derivative of XC-functional
!> \param xc_section            XC input section
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param do_admm               perform auxiliary density matrix method calculations
!> \param matrix_s_aux_fit_sub  overlap matrix for the auxiliary basis set
!>                              (ignored if do_admm == .false.)
!> \param blacs_env_sub         BLACS parallel environment
!> \param pw_env_sub            plain wave environment
!> \param task_list_orb_sub     task list for the primary basis set
!> \param task_list_aux_fit_sub task list for the auxiliary basis set
!>                              (ignored if do_admm == .false.)
!> \param nsubgroups            number of MPI subgroups
!> \param subgroup              list of sub-group members
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_Aop_evects(Aop_evects, evects, fm_struct_evects, gs_mos, &
                                        is_rks_triplets, matrix_ks_sub, matrix_s_sub, qs_env, &
                                        xc_rho_set, xc_deriv_set, xc_section, do_hfx, do_admm, &
                                        matrix_s_aux_fit_sub, blacs_env_sub, pw_env_sub, &
                                        task_list_orb_sub, task_list_aux_fit_sub, nsubgroups, subgroup)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(cp_fm_struct_p_type), DIMENSION(:), &
         INTENT(in)                                      :: fm_struct_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(cp_dbcsr_p_type), DIMENSION(:), INTENT(in)    :: matrix_ks_sub
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s_sub
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(in)                                :: do_hfx, do_admm
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s_aux_fit_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_orb_sub, task_list_aux_fit_sub
      INTEGER, INTENT(in)                                :: nsubgroups
      INTEGER, DIMENSION(0:), INTENT(in)                 :: subgroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_Aop_evects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0)
         CPASSERT(SIZE(Aop_evects, 1) == nspins)
         CPASSERT(SIZE(Aop_evects, 2) == nvects)
         CPASSERT(SIZE(gs_mos) == nspins)
      END IF

      IF (nvects > 0) THEN
         CALL get_qs_env(qs_env, blacs_env=blacs_env_global)

         DO ivect = 1, nvects
            DO ispin = 1, nspins
               IF (ASSOCIATED(evects(ispin, ivect)%matrix) .AND. .NOT. ASSOCIATED(Aop_evects(ispin, ivect)%matrix)) &
                  CALL cp_fm_create(Aop_evects(ispin, ivect)%matrix, fm_struct_evects(ispin)%struct)
            END DO
         END DO

         ! initialise hat{A} C_{1,i} by orbital energy difference
         CALL tddfpt_init_by_energy_diff(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, &
                                         matrix_ks=matrix_ks_sub, matrix_s=matrix_s_sub)

         ! add contributions from coulomb terms
         IF (.NOT. is_rks_triplets) &
            CALL tddfpt_apply_coulomb(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, matrix_s=matrix_s_sub, &
                                      qs_env=qs_env, blacs_env=blacs_env_sub, pw_env=pw_env_sub, task_list=task_list_orb_sub)

         ! add contributions from the adiabatic TDDFT kernel
         CALL tddfpt_apply_kernel(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, is_rks_triplets=is_rks_triplets, &
                                  matrix_s=matrix_s_sub, qs_env=qs_env, xc_rho_set=xc_rho_set, xc_deriv_set=xc_deriv_set, &
                                  xc_section=xc_section, blacs_env=blacs_env_sub, pw_env=pw_env_sub, task_list=task_list_orb_sub)

         ! add contributions from exact-exchange terms
         IF (do_hfx) &
            CALL tddfpt_apply_hfx(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, qs_env=qs_env, do_admm=do_admm, &
                                  blacs_env=blacs_env_sub, matrix_s=matrix_s_sub, matrix_s_aux_fit=matrix_s_aux_fit_sub, &
                                  pw_env=pw_env_sub, task_list_orb=task_list_orb_sub, task_list_aux_fit=task_list_aux_fit_sub, &
                                  nsubgroups=nsubgroups, subgroup=subgroup)
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_Aop_evects

! **************************************************************************************************
!> \brief Solve eigenproblem for the reduced action matrix and find new Ritz eigenvectors and
!>        eigenvalues.
!> \param ritz_vects       Ritz eigenvectors (initialised on exit)
!> \param Aop_ritz         action matrix -- Ritz vector product (initialised on exit)
!> \param evals            Ritz eigenvalues (initialised on exit)
!> \param Aop_evects       action matrix -- trial vector product
!> \param evects           TDDFPT trial vectors
!> \param gs_mos           molecular orbitals optimized for the ground state
!> \param blacs_env_global BLACS parallel environment involving all the processors
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_ritz_vects(ritz_vects, Aop_ritz, evals, Aop_evects, evects, gs_mos, blacs_env_global)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_ritz_vects', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, icol_global, icol_local, iproc, irow_global, irow_local, ispin, &
         mepos_global, mepos_local, ncols_local, nkvs, nprocs_local, nrows_local, nrvs, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ncols_procs, nrows_procs
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_local, row_indices_local
      REAL(kind=dp)                                      :: act, trace
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: evects_Atilde
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: Atilde_fm, evects_Atilde_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(tddfpt_int_vector), ALLOCATABLE, DIMENSION(:) :: col_indices_procs, row_indices_procs

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nkvs = SIZE(evects, 2)
      nrvs = SIZE(ritz_vects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0 .AND. nkvs > 0)
      END IF

      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, context=blacs_env_sub, para_env=para_env_sub)

      mepos_global = blacs_env_global%para_env%mepos
      mepos_local = para_env_sub%mepos
      nprocs_local = para_env_sub%num_pe

      CALL cp_fm_struct_create(fm_struct, nrow_global=nkvs, ncol_global=nkvs, context=blacs_env_global)
      CALL cp_fm_create(Atilde_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(Atilde_fm, nrow_local=nrows_local, ncol_local=ncols_local, &
                          row_indices=row_indices_local, col_indices=col_indices_local, local_data=local_data)

      ALLOCATE (nrows_procs(nprocs_local), ncols_procs(nprocs_local))
      nrows_procs(:) = 0
      ncols_procs(:) = 0
      nrows_procs(mepos_local+1) = nrows_local
      ncols_procs(mepos_local+1) = ncols_local

      CALL mp_max(nrows_procs, para_env_sub%group)
      CALL mp_max(ncols_procs, para_env_sub%group)

      ALLOCATE (row_indices_procs(nprocs_local), col_indices_procs(nprocs_local))
      DO iproc = 1, nprocs_local
         ALLOCATE (row_indices_procs(iproc)%vect(nrows_procs(iproc)))
         ALLOCATE (col_indices_procs(iproc)%vect(ncols_procs(iproc)))

         row_indices_procs(iproc)%vect(:) = 0
         col_indices_procs(iproc)%vect(:) = 0
      END DO

      row_indices_procs(mepos_local+1)%vect(:) = row_indices_local(1:nrows_local)
      col_indices_procs(mepos_local+1)%vect(:) = col_indices_local(1:ncols_local)

      DO iproc = 1, nprocs_local
         CALL mp_max(row_indices_procs(iproc)%vect, para_env_sub%group)
         CALL mp_max(col_indices_procs(iproc)%vect, para_env_sub%group)
      END DO

      ! *** compute upper-diagonal reduced action matrix ***
      DO iproc = 1, nprocs_local
         DO icol_local = 1, ncols_procs(iproc)
            icol_global = col_indices_procs(iproc)%vect(icol_local)

            DO irow_local = 1, nrows_procs(iproc)
               irow_global = row_indices_procs(iproc)%vect(irow_local)

               IF (irow_global <= icol_global) THEN
                  act = 0.0_dp
                  DO ispin = 1, nspins
                     CALL cp_fm_trace(Aop_evects(ispin, irow_global)%matrix, evects(ispin, icol_global)%matrix, trace)
                     act = act+trace
                  END DO

                  IF (mepos_local+1 == iproc) &
                     local_data(irow_local, icol_local) = act
               END IF
            END DO
         END DO
      END DO

      DO iproc = nprocs_local, 1, -1
         DEALLOCATE (row_indices_procs(iproc)%vect, col_indices_procs(iproc)%vect)
      END DO
      DEALLOCATE (nrows_procs, ncols_procs)

      ! *** solve eigenproblem for reduced matrices ***
      CALL choose_eigv_solver(Atilde_fm, evects_Atilde_fm, evals(1:nkvs))

      ALLOCATE (evects_Atilde(nkvs, nrvs))
      CALL cp_fm_get_submatrix(evects_Atilde_fm, evects_Atilde, start_row=1, start_col=1, n_rows=nkvs, n_cols=nrvs)
      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(Atilde_fm)

      ! *** compute Ritz vectors ***
      DO icol_global = 1, nrvs
         IF (ASSOCIATED(ritz_vects(1, icol_global)%matrix)) THEN
            DO ispin = 1, nspins
               CALL cp_fm_set_all(ritz_vects(ispin, icol_global)%matrix, 0.0_dp)
               CALL cp_fm_set_all(Aop_ritz(ispin, icol_global)%matrix, 0.0_dp)
            END DO

            DO irow_global = 1, nkvs
               act = evects_Atilde(irow_global, icol_global)
               DO ispin = 1, nspins
                  CALL cp_fm_scale_and_add(1.0_dp, ritz_vects(ispin, icol_global)%matrix, &
                                           act, evects(ispin, irow_global)%matrix)
                  CALL cp_fm_scale_and_add(1.0_dp, Aop_ritz(ispin, icol_global)%matrix, &
                                           act, Aop_evects(ispin, irow_global)%matrix)
               END DO
            END DO
         END IF
      END DO

      DEALLOCATE (evects_Atilde)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_compute_ritz_vects

! **************************************************************************************************
!> \brief Expand Krylov space by computing residual vectors.
!> \param residual_vects   residual vectors
!> \param evals            Ritz eigenvalues
!> \param Aop_ritz         action matrix -- Ritz vector product
!> \param ritz_vects       Ritz eigenvectors
!> \param fm_struct_evects structures of newly allocated dense matrices
!> \param gs_mos           molecular orbitals optimised for the ground state
!> \param guess_vectors    initial guess vectors
!> \param matrix_s         overlap matrix
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_residual_vects(residual_vects, evals, Aop_ritz, ritz_vects, &
                                            fm_struct_evects, gs_mos, guess_vectors, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: residual_vects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_ritz, ritz_vects
      TYPE(cp_fm_struct_p_type), DIMENSION(:), &
         INTENT(in)                                      :: fm_struct_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(tddfpt_guess_vectors), INTENT(in)             :: guess_vectors
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_residual_vects', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, imo_occ, imo_virt, irv, ispin, &
                                                            istate, nrvs, nspins, nstates_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      REAL(kind=dp)                                      :: eref, minus_lambda, weight
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: Ab_minus_eb, matrix_a, matrix_b
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(tddfpt_fm_vector), ALLOCATABLE, DIMENSION(:)  :: Ab_cols, residual_cols

      CALL timeset(routineN, handle)

      nspins = SIZE(ritz_vects, 1)
      nrvs = SIZE(residual_vects, 2)
      nstates_total = SIZE(guess_vectors%evals)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(ritz_vects, 2) >= nrvs)
      END IF

      IF (nrvs > 0) THEN
         DO istate = 1, nrvs
            DO ispin = 1, nspins
               IF (ASSOCIATED(ritz_vects(ispin, istate)%matrix)) THEN
                  IF (.NOT. ASSOCIATED(residual_vects(ispin, istate)%matrix)) &
                     CALL cp_fm_create(residual_vects(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
               ELSE
                  IF (ASSOCIATED(residual_vects(ispin, istate)%matrix)) &
                     CALL cp_fm_release(residual_vects(ispin, istate)%matrix)
               END IF
            END DO
         END DO

         ALLOCATE (nmo_occ(nspins), Ab_minus_eb(nspins))
         DO ispin = 1, nspins
            NULLIFY (Ab_minus_eb(ispin)%matrix)
            nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
            CALL cp_fm_create(Ab_minus_eb(ispin)%matrix, fm_struct_evects(ispin)%struct)
         END DO

         ALLOCATE (Ab_cols(nspins), residual_cols(nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_struct)
            CALL cp_fm_get_info(gs_mos(ispin)%vectors_mo_virt(1)%matrix, matrix_struct=fm_struct)

            ALLOCATE (Ab_cols(ispin)%vect(nmo_occ(ispin)))
            ALLOCATE (residual_cols(ispin)%vect(nmo_occ(ispin)))
            DO imo_occ = 1, nmo_occ(ispin)
               NULLIFY (Ab_cols(ispin)%vect(imo_occ)%matrix, residual_cols(ispin)%vect(imo_occ)%matrix)
               CALL cp_fm_create(Ab_cols(ispin)%vect(imo_occ)%matrix, fm_struct)
               CALL cp_fm_create(residual_cols(ispin)%vect(imo_occ)%matrix, fm_struct)
            END DO
         END DO

         ALLOCATE (matrix_a(nstates_total), matrix_b(nstates_total), weights(nstates_total))

         ! *** actually compute residual vectors ***
         DO irv = 1, nrvs
            IF (ASSOCIATED(ritz_vects(1, irv)%matrix)) THEN
               minus_lambda = -evals(irv)
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(matrix_s, ritz_vects(ispin, irv)%matrix, Ab_minus_eb(ispin)%matrix, &
                                               ncol=nmo_occ(ispin), alpha=minus_lambda, beta=0.0_dp)
                  CALL cp_fm_scale_and_add(1.0_dp, Ab_minus_eb(ispin)%matrix, 1.0_dp, Aop_ritz(ispin, irv)%matrix)
               END DO

               DO ispin = 1, nspins
                  DO imo_occ = nmo_occ(ispin), 1, -1
                     CALL cp_fm_to_fm(Ab_minus_eb(ispin)%matrix, Ab_cols(ispin)%vect(imo_occ)%matrix, &
                                      ncol=1, source_start=imo_occ, target_start=1)
                     CALL cp_fm_set_all(residual_cols(ispin)%vect(imo_occ)%matrix, 0.0_dp)
                  END DO
               END DO

               DO istate = 1, nstates_total
                  ispin = guess_vectors%ispins(istate)
                  imo_occ = guess_vectors%imos_occ(istate)
                  imo_virt = guess_vectors%imos_virt(istate)

                  matrix_a(istate)%matrix => gs_mos(ispin)%vectors_mo_virt(imo_virt)%matrix
                  matrix_b(istate)%matrix => Ab_cols(ispin)%vect(imo_occ)%matrix
               END DO

               CALL cp_fm_p_trace(matrix_a, matrix_b, weights)

               DO istate = 1, nstates_total
                  ispin = guess_vectors%ispins(istate)
                  imo_occ = guess_vectors%imos_occ(istate)
                  imo_virt = guess_vectors%imos_virt(istate)

                  eref = guess_vectors%evals(istate)+minus_lambda
                  IF (ABS(eref) < threshold) &
                     eref = guess_vectors%evals(istate)+eref_scale*minus_lambda

                  weight = weights(istate)/eref
                  CALL cp_fm_scale_and_add(1.0_dp, residual_cols(ispin)%vect(imo_occ)%matrix, &
                                           weight, gs_mos(ispin)%vectors_mo_virt(imo_virt)%matrix)
               END DO

               DO ispin = 1, nspins
                  DO imo_occ = nmo_occ(ispin), 1, -1
                     CALL cp_fm_to_fm(residual_cols(ispin)%vect(imo_occ)%matrix, residual_vects(ispin, irv)%matrix, &
                                      ncol=1, source_start=1, target_start=imo_occ)
                  END DO
               END DO
            END IF
         END DO

         DEALLOCATE (matrix_a, matrix_b, weights)

         DO ispin = SIZE(Ab_cols), 1, -1
            DO imo_occ = SIZE(Ab_cols(ispin)%vect), 1, -1
               CALL cp_fm_release(residual_cols(ispin)%vect(imo_occ)%matrix)
               CALL cp_fm_release(Ab_cols(ispin)%vect(imo_occ)%matrix)
            END DO
            DEALLOCATE (residual_cols(ispin)%vect, Ab_cols(ispin)%vect)
         END DO
         DEALLOCATE (residual_cols, Ab_cols)

         DO ispin = SIZE(Ab_minus_eb), 1, -1
            CALL cp_fm_release(Ab_minus_eb(ispin)%matrix)
         END DO
         DEALLOCATE (Ab_minus_eb)
         DEALLOCATE (nmo_occ)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_compute_residual_vects

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and analysis.
!> \param log_unit         output unit
!> \param evects           TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                         SIZE(evects,2) -- number of excited states to print)
!> \param evals            TDDFPT eigenvalues
!> \param mult             multiplicity
!> \param fm_struct_evects structures of newly allocated dense matrices
!> \param gs_mos           molecular orbitals optimised for the ground state
!> \param matrix_s         overlap matrix
!> \param dipole_op_xyz    components of the dipole operator in terms of contracted basis functions
!> \param dipole_form      form of the dipole operator (length / vector)
!> \param min_amplitude    the smallest excitation amplitude to print
!> \param blacs_env_global BLACS parallel environment involving all the processors
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!>    * 07.2016 spin-unpolarised electron density [Sergey Chulkov]
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
   SUBROUTINE tddfpt_print_summary(log_unit, evects, evals, mult, fm_struct_evects, gs_mos, &
                                   matrix_s, dipole_op_xyz, dipole_form, min_amplitude, blacs_env_global)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: mult
      TYPE(cp_fm_struct_p_type), DIMENSION(:), &
         INTENT(in)                                      :: fm_struct_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: dipole_op_xyz
      INTEGER, INTENT(in)                                :: dipole_form
      REAL(kind=dp), INTENT(in)                          :: min_amplitude
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_print_summary', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=1)                                   :: lsd_str
      CHARACTER(len=20)                                  :: mult_str
      CHARACTER(len=5)                                   :: spin_label
      INTEGER :: handle, icol, ideriv, iexc, imo_occ, imo_virt, ind, irow, ispin, istate, nao, &
         ncols_local, nmo_virt_occ, nmo_virt_occ_alpha, nmo_virt_occ_spin, nrows_local, nspins, &
         nstates, state_spin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds, nexcitations, nmo_occ, nmo_virt, &
                                                            sendrecv_handlers_amp, &
                                                            sendrecv_handlers_ind, source_nodes
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_exc_analysis
      REAL(kind=dp)                                      :: eval_occ, osc_strength, trans_dipole
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_1d, weights_neg_abs_1d
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: trans_dipoles
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: ediff_inv, ediff_inv_weights, &
                                                            S_mos_virt, scaled_evect, weights_fm
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dS_mos_occ
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env_global, para_env_sub
      TYPE(tddfpt_allocatable_matrix), ALLOCATABLE, &
         DIMENSION(:)                                    :: weights_2d
      TYPE(tddfpt_int_vector), ALLOCATABLE, DIMENSION(:) :: exc_indices
      TYPE(tddfpt_real_vector), ALLOCATABLE, &
         DIMENSION(:)                                    :: exc_amplitudes

      EXTERNAL :: dcopy

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      CALL cp_fm_struct_get(fm_struct_evects(1)%struct, nrow_global=nao, context=blacs_env_sub, para_env=para_env_sub)
      CALL get_blacs_info(blacs_env_global, para_env=para_env_global)
      do_exc_analysis = min_amplitude < 1.0_dp

      IF (nspins > 1) THEN
         lsd_str = 'U'
      ELSE
         lsd_str = 'R'
      END IF

      ALLOCATE (nmo_occ(nspins), nmo_virt(nspins))
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      ALLOCATE (S_mos_virt(nspins), dS_mos_occ(nderivs, nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct, S_mos_virt(ispin)%matrix)
         CALL cp_fm_get_info(gs_mos(ispin)%matrix_mos_virt, matrix_struct=fm_struct)
         CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)

         DO ideriv = 1, nderivs
            NULLIFY (dS_mos_occ(ideriv, ispin)%matrix)
            CALL cp_fm_create(dS_mos_occ(ideriv, ispin)%matrix, fm_struct_evects(ispin)%struct)
         END DO
      END DO

      ALLOCATE (weights_fm(nspins), ediff_inv(nspins), ediff_inv_weights(nspins), scaled_evect(nspins))
      DO ispin = 1, nspins
         NULLIFY (weights_fm(ispin)%matrix, ediff_inv(ispin)%matrix, ediff_inv_weights(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), context=blacs_env_sub)

         CALL cp_fm_create(weights_fm(ispin)%matrix, fm_struct)
         CALL cp_fm_create(ediff_inv(ispin)%matrix, fm_struct)
         CALL cp_fm_create(ediff_inv_weights(ispin)%matrix, fm_struct)

         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_create(scaled_evect(ispin)%matrix, fm_struct_evects(ispin)%struct)
      END DO

      ! compute (E_virt_a - E_occ_i)^{-1}
      DO ispin = 1, nspins
         CALL cp_fm_get_info(ediff_inv(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(eval_occ, icol, irow), &
!$OMP             SHARED(col_indices, gs_mos, ispin, local_data, ncols_local, nrows_local, row_indices)
         DO icol = 1, ncols_local
            ! E_occ_i ; imo_occ = col_indices(icol)
            eval_occ = gs_mos(ispin)%evals_occ(col_indices(icol))

            DO irow = 1, nrows_local
               ! ediff_inv_weights(a, i) = 1.0 / (E_virt_a - E_occ_i)
               ! imo_virt = row_indices(irow)
               local_data(irow, icol) = 1.0_dp/(gs_mos(ispin)%evals_virt(row_indices(irow))-eval_occ)
            END DO
         END DO
!$OMP END PARALLEL DO
      END DO

      DO ispin = 1, nspins
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, gs_mos(ispin)%matrix_mos_virt, S_mos_virt(ispin)%matrix, &
                                      ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)
         DO ideriv = 1, nderivs
            CALL cp_dbcsr_sm_fm_multiply(dipole_op_xyz(ideriv)%matrix, gs_mos(ispin)%matrix_mos_occ, &
                                         dS_mos_occ(ideriv, ispin)%matrix, ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         END DO
      END DO

      ALLOCATE (trans_dipoles(nderivs, nstates))
      trans_dipoles(:, :) = 0.0_dp

      IF (do_exc_analysis) THEN
         nmo_virt_occ_alpha = nmo_virt(1)*nmo_occ(1)
         nmo_virt_occ = DOT_PRODUCT(nmo_virt(:), nmo_occ(:))
         ALLOCATE (weights_1d(nmo_virt_occ), weights_neg_abs_1d(nmo_virt_occ), inds(nmo_virt_occ))
         ALLOCATE (weights_2d(nspins))
         DO ispin = 1, nspins
            ALLOCATE (weights_2d(ispin)%matrix(nmo_virt(ispin), nmo_occ(ispin)))
         END DO

         ALLOCATE (exc_amplitudes(nstates), exc_indices(nstates))
      END IF

      DO istate = 1, nstates
         IF (ASSOCIATED(evects(1, istate)%matrix)) THEN
            nmo_virt_occ = 0
            DO ispin = 1, nspins
               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin)%matrix, &
                            evects(ispin, istate)%matrix, 0.0_dp, weights_fm(ispin)%matrix)

               ! compute electric dipole moments
               IF (dipole_form == xas_dip_vel) THEN
                  CALL cp_fm_schur_product(weights_fm(ispin)%matrix, ediff_inv(ispin)%matrix, ediff_inv_weights(ispin)%matrix)
                  CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nmo_virt(ispin), 1.0_dp, gs_mos(ispin)%matrix_mos_virt, &
                               ediff_inv_weights(ispin)%matrix, 0.0_dp, scaled_evect(ispin)%matrix)
               END IF

               DO ideriv = 1, nderivs
                  SELECT CASE (dipole_form)
                  CASE (xas_dip_len)
                     CALL cp_fm_trace(dS_mos_occ(ideriv, ispin)%matrix, evects(ispin, istate)%matrix, trans_dipole)

                  CASE (xas_dip_vel)
                     CALL cp_fm_trace(dS_mos_occ(ideriv, ispin)%matrix, scaled_evect(ispin)%matrix, trans_dipole)

                  CASE DEFAULT
                     CPABORT("Unimplemented form of the dipole operator")
                  END SELECT

                  IF (para_env_sub%mepos == para_env_sub%source) &
                     trans_dipoles(ideriv, istate) = trans_dipoles(ideriv, istate)+trans_dipole
               END DO

               ! compute excitation amplitudes
               IF (do_exc_analysis) THEN
                  CALL cp_fm_get_submatrix(weights_fm(ispin)%matrix, weights_2d(ispin)%matrix)

                  IF (para_env_sub%mepos == para_env_sub%source) THEN
                     nmo_virt_occ_spin = nmo_virt(ispin)*nmo_occ(ispin)
                     CALL dcopy(nmo_virt_occ_spin, weights_2d(ispin)%matrix, 1, &
                                weights_1d(nmo_virt_occ+1:nmo_virt_occ+nmo_virt_occ_spin), 1)
                     nmo_virt_occ = nmo_virt_occ+nmo_virt_occ_spin
                  END IF
               END IF
            END DO

            IF (do_exc_analysis .AND. para_env_sub%mepos == para_env_sub%source) THEN
               weights_neg_abs_1d(:) = -ABS(weights_1d)
               CALL sort(weights_neg_abs_1d, nmo_virt_occ, inds)

               DO iexc = 1, nmo_virt_occ
                  IF (-weights_neg_abs_1d(iexc) < min_amplitude) EXIT
               END DO
               nrows_local = iexc-1

               IF (nrows_local > 0) THEN
                  ALLOCATE (exc_amplitudes(istate)%vect(nrows_local))
                  ALLOCATE (exc_indices(istate)%vect(nrows_local))

                  exc_indices(istate)%vect(:) = inds(1:nrows_local)
                  DO iexc = 1, nrows_local
                     exc_amplitudes(istate)%vect(iexc) = weights_1d(inds(iexc))
                  END DO
               END IF
            END IF
         END IF
      END DO

      CALL mp_sum(trans_dipoles, para_env_global%group)

      IF (do_exc_analysis) THEN
         ALLOCATE (source_nodes(nstates), nexcitations(nstates))
         ALLOCATE (sendrecv_handlers_amp(nstates), sendrecv_handlers_ind(nstates))

         source_nodes(:) = 0
         nexcitations(:) = 0
         IF (para_env_sub%mepos == para_env_sub%source) THEN
            DO istate = 1, nstates
               IF (ALLOCATED(exc_amplitudes(istate)%vect)) THEN
                  source_nodes(istate) = para_env_global%mepos
                  nexcitations(istate) = SIZE(exc_amplitudes(istate)%vect)
               END IF
            END DO
         END IF

         CALL mp_sum(source_nodes, para_env_global%group)
         CALL mp_sum(nexcitations, para_env_global%group)

         nrows_local = 0
         IF (para_env_global%mepos == para_env_global%source) THEN
            DO istate = 1, nstates
               IF (.NOT. ALLOCATED(exc_amplitudes(istate)%vect)) THEN
                  ALLOCATE (exc_amplitudes(istate)%vect(nexcitations(istate)))
                  ALLOCATE (exc_indices(istate)%vect(nexcitations(istate)))

                  nrows_local = nrows_local+1
                  CALL mp_irecv(exc_amplitudes(istate)%vect, source_nodes(istate), para_env_global%group, &
                                sendrecv_handlers_amp(nrows_local), 2*(istate-1))
                  CALL mp_irecv(exc_indices(istate)%vect, source_nodes(istate), para_env_global%group, &
                                sendrecv_handlers_ind(nrows_local), 2*(istate-1)+1)
               END IF
            END DO
         ELSE
            DO istate = 1, nstates
               IF (ALLOCATED(exc_amplitudes(istate)%vect)) THEN
                  nrows_local = nrows_local+1

                  CALL mp_isend(exc_amplitudes(istate)%vect, para_env_global%source, para_env_global%group, &
                                sendrecv_handlers_amp(nrows_local), 2*(istate-1))
                  CALL mp_isend(exc_indices(istate)%vect, para_env_global%source, para_env_global%group, &
                                sendrecv_handlers_ind(nrows_local), 2*(istate-1)+1)
               END IF
            END DO
         END IF

         DO irow = 1, nrows_local
            CALL mp_wait(sendrecv_handlers_amp(irow))
            CALL mp_wait(sendrecv_handlers_ind(irow))
         END DO

         DEALLOCATE (source_nodes, sendrecv_handlers_amp, sendrecv_handlers_ind)
      END IF

      IF (do_exc_analysis) THEN
         DO ispin = nspins, 1, -1
            DEALLOCATE (weights_2d(ispin)%matrix)
         END DO
         DEALLOCATE (weights_1d, weights_neg_abs_1d, weights_2d, inds)
      END IF

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(scaled_evect(ispin)%matrix)
         CALL cp_fm_release(ediff_inv_weights(ispin)%matrix)
         CALL cp_fm_release(ediff_inv(ispin)%matrix)

         DO ideriv = nderivs, 1, -1
            CALL cp_fm_release(dS_mos_occ(ideriv, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (dS_mos_occ, ediff_inv)

      ! *** summary information ***
      IF (log_unit > 0) THEN
         CALL integer_to_string(mult, mult_str)
         WRITE (log_unit, '(/,1X,A1,A,1X,A,/)') lsd_str, "-TDDFPT states of multiplicity", TRIM(mult_str)

         WRITE (log_unit, '(T11,A,T21,A,T37,A,T69,A)') "State", "Excitation", "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T11,A,T21,A,T39,A,T49,A,T59,A,T67,A)') "number", "energy (eV)", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T10,72("-"))')

         DO istate = 1, nstates
            osc_strength = 2.0_dp/3.0_dp*evals(istate)*DOT_PRODUCT(trans_dipoles(:, istate), trans_dipoles(:, istate))

            WRITE (log_unit, '(1X,A,T9,I8,T21,F11.5,T34,3(1X,F9.4),T69,F10.5)') &
               "TDDFPT|", istate, evals(istate)*evolt, trans_dipoles(1:nderivs, istate), osc_strength
         END DO

         IF (do_exc_analysis) THEN
            WRITE (log_unit, '(/,1X,A,/)') "Excitation analysis"

            WRITE (log_unit, '(3X,A,T17,A,T34,A,T49,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(3X,A,T18,A,T34,A,T49,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,57("-"))')

            IF (nspins == 1) THEN
               state_spin = 1
               spin_label = '     '
            END IF

            DO istate = 1, nstates
               WRITE (log_unit, '(1X,I8)') istate

               IF (nexcitations(istate) > 0) THEN

                  DO iexc = 1, nexcitations(istate)
                     ind = exc_indices(istate)%vect(iexc)-1

                     IF (nspins > 1) THEN
                        IF (ind < nmo_virt_occ_alpha) THEN
                           state_spin = 1
                           spin_label = '(alp)'
                        ELSE
                           state_spin = 2
                           ind = ind-nmo_virt_occ_alpha
                           spin_label = '(bet)'
                        END IF
                     END IF

                     imo_occ = ind/nmo_virt(state_spin)+1
                     imo_virt = MOD(ind, nmo_virt(state_spin))+1

                     WRITE (log_unit, '(T14,I8,1X,A5,T30,I8,1X,A5,T50,F9.6)') imo_occ, spin_label, &
                        nmo_occ(state_spin)+imo_virt, spin_label, exc_amplitudes(istate)%vect(iexc)
                  END DO
               END IF
            END DO
         END IF
      END IF

      IF (do_exc_analysis) THEN
         DO istate = 1, nstates
            IF (ALLOCATED(exc_amplitudes(istate)%vect)) &
               DEALLOCATE (exc_amplitudes(istate)%vect, exc_indices(istate)%vect)
         END DO
         DEALLOCATE (nexcitations, exc_amplitudes, exc_indices)
      END IF

      DEALLOCATE (trans_dipoles)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(weights_fm(ispin)%matrix)
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt, weights_fm)

      DEALLOCATE (nmo_virt, nmo_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary

! **************************************************************************************************
!> \brief Write Ritz vectors to a binary restart file.
!> \param evects               vectors to store
!> \param evals                TDDFPT eigenvalues
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \param blacs_env_global     BLACS parallel environment involving all the processor
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_write_restart(evects, evals, gs_mos, logger, tddfpt_print_section, blacs_env_global)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_write_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate, nao, nspins, &
                                                            nstates, ounit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: evects_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      IF (BTEST(cp_print_key_should_output(logger%iter_info, tddfpt_print_section, "RESTART"), cp_p_file)) THEN
         CALL timeset(routineN, handle)

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)

         IF (debug_this_module) THEN
            CPASSERT(SIZE(evals) == nstates)
            CPASSERT(nspins > 0)
            CPASSERT(nstates > 0)
         END IF

         CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao)

         ALLOCATE (nmo_occ(nspins), evects_global(nspins))
         DO ispin = 1, nspins
!            CALL cp_fm_struct_get(fm_struct_evects(ispin)%struct, nrow_global=nao, ncol_global=nmo_occ(ispin))
            nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)

            NULLIFY (fm_struct, evects_global(ispin)%matrix)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ(ispin), context=blacs_env_global)
            CALL cp_fm_create(evects_global(ispin)%matrix, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         ounit = cp_print_key_unit_nr(logger, tddfpt_print_section, "RESTART", &
                                      extension=".tdwfn", file_status="REPLACE", file_action="WRITE", &
                                      do_backup=.TRUE., file_form="UNFORMATTED")

         IF (ounit > 0) THEN
            WRITE (ounit) nstates, nspins, nao
            WRITE (ounit) nmo_occ
            WRITE (ounit) evals
         END IF

         DO istate = 1, nstates
            DO ispin = 1, nspins
               CALL cp_fm_copy_general(source=evects(ispin, istate)%matrix, &
                                       destination=evects_global(ispin)%matrix, &
                                       global_context=blacs_env_global%group)

               ! TDDFPT wave function is actually stored as a linear combination of virtual MOs
               ! that replaces the corresponding deoccupied MO. Unfortunately, the phase
               ! of the occupied MOs varies depending on the eigensolver used as well as
               ! how eigenvectors are distributed across computational cores. The phase is important
               ! because TDDFPT wave functions are used to compute a response electron density
               ! \rho^{-} = 1/2 * [C_{0} * evect^T + evect * C_{0}^{-}], where C_{0} is the expansion
               ! coefficients of the reference ground-state wave function. To make the restart file
               ! transferable, TDDFPT wave functions are stored in assumption that all ground state
               ! MOs have a positive phase.
               CALL cp_fm_column_scale(evects_global(ispin)%matrix, gs_mos(ispin)%phases_occ)

               CALL cp_fm_write_unformatted(evects_global(ispin)%matrix, ounit)
            END DO
         END DO

         DO ispin = nspins, 1, -1
            CALL cp_fm_release(evects_global(ispin)%matrix)
         END DO
         DEALLOCATE (nmo_occ, evects_global)

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
!> \param blacs_env_global     BLACS parallel environment involving all the processor
!> \retval nstates_read        the number of excited states found in the restart file
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_read_restart(evects, evals, gs_mos, logger, tddfpt_section, &
                                tddfpt_print_section, blacs_env_global) RESULT(nstates_read)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_section, tddfpt_print_section
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      INTEGER                                            :: nstates_read

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_read_restart', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=20)                                  :: read_str, ref_str
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, ispin, istate, iunit, n_rep_val, &
                                                            nao, nao_read, nspins, nspins_read, &
                                                            nstates
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nmo_occ, nmo_occ_read
      LOGICAL                                            :: file_exists
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_read
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: evects_read
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
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
      CALL cp_fm_get_info(gs_mos(1)%matrix_mos_occ, nrow_global=nao)

      ALLOCATE (evects_read(nspins), nmo_occ(nspins), nmo_occ_read(nspins))
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         NULLIFY (fm_struct, evects_read(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ(ispin), context=blacs_env_global)
         CALL cp_fm_create(evects_read(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      IF (para_env_global%ionode) THEN
         READ (iunit) nstates_read, nspins_read, nao_read

         IF (nspins_read /= nspins) THEN
            CALL integer_to_string(nspins, ref_str)
            CALL integer_to_string(nspins_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Restarted TDDFPT wave function contains incompartible number of spin components ("// &
                          TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         IF (nao_read /= nao) THEN
            CALL integer_to_string(nao, ref_str)
            CALL integer_to_string(nao_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Incompartible number of atomic orbitals ("//TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         READ (iunit) nmo_occ_read

         DO ispin = 1, nspins
            IF (nmo_occ_read(ispin) /= nmo_occ(ispin)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Incompartible number of electrons and/or multiplicity.")
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
            CALL cp_fm_read_unformatted(evects_read(ispin)%matrix, iunit)

            IF (istate <= nstates) THEN
               CALL cp_fm_column_scale(evects_read(ispin)%matrix, gs_mos(ispin)%phases_occ)

               CALL cp_fm_copy_general(source=evects_read(ispin)%matrix, &
                                       destination=evects(ispin, istate)%matrix, &
                                       global_context=blacs_env_global%group)
            END IF
         END DO
      END DO

      IF (para_env_global%ionode) &
         CALL close_file(unit_number=iunit)

      DO ispin = nspins, 1, -1
         IF (ASSOCIATED(evects_read(ispin)%matrix)) &
            CALL cp_fm_release(evects_read(ispin)%matrix)
      END DO
      DEALLOCATE (evects_read, nmo_occ, nmo_occ_read)

      CALL timestop(handle)
   END FUNCTION tddfpt_read_restart
END MODULE qs_tddfpt2_methods
