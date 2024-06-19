!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_subgroups
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_get_info,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_2d
   USE hartree_local_methods,           ONLY: init_coulomb_local
   USE hartree_local_types,             ONLY: hartree_local_create,&
                                              hartree_local_release,&
                                              hartree_local_type
   USE input_constants,                 ONLY: tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_release,&
                                              pw_env_retain,&
                                              pw_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_rho0_ggrid,                   ONLY: rho0_s_grid_create
   USE qs_rho0_methods,                 ONLY: init_rho0
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_subgroups'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: tddfpt_subgroup_env_type
   PUBLIC :: tddfpt_sub_env_init, tddfpt_sub_env_release
   PUBLIC :: tddfpt_dbcsr_create_by_dist, tddfpt_fm_replicate_across_subgroups

! **************************************************************************************************
!> \brief Parallel (sub)group environment.
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_subgroup_env_type
      !> indicates that the global MPI communicator has been split into subgroups; if it is .FALSE.
      !> certain components of the structure (blacs_env, para_env, admm_A, and mos_occ)
      !> can still be accessed; in this case they simply point to the corresponding global variables
      LOGICAL                                            :: is_split = .FALSE.
      !> number of parallel groups
      INTEGER                                            :: ngroups = -1
      !> group_distribution(0:ngroups-1) : a process with rank 'i' belongs to the parallel group
      !> with index 'group_distribution(i)'
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: group_distribution
      !> group-specific BLACS parallel environment
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => NULL()
      !> group-specific MPI parallel environment
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      !> occupied MOs stored in a matrix form [nao x nmo_occ(spin)] distributed across processes
      !> in the parallel group
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)      :: mos_occ
      !> group-specific copy of the ADMM A matrix 'admm_type%A'
      TYPE(cp_fm_type), POINTER                          :: admm_A => NULL()
      !
      !> indicates that a set of multi-grids has been allocated; if it is .FALSE. all the components
      !> below point to the corresponding global variables and can be accessed
      LOGICAL                                            :: is_mgrid = .FALSE.
      !> group-specific DBCSR distribution
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist => NULL()
      !> group-specific two-dimensional distribution of pairs of particles
      TYPE(distribution_2d_type), POINTER                :: dist_2d => NULL()
      !> group-specific plane wave environment
      TYPE(pw_env_type), POINTER                         :: pw_env => NULL()
      !> lists of neighbours in auxiliary and primary basis sets
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                           :: sab_aux_fit => NULL(), sab_orb => NULL()
      !> task lists in auxiliary and primary basis sets
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit => NULL(), task_list_orb => NULL()
      !> soft task lists in auxiliary and primary basis sets
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit_soft => NULL(), task_list_orb_soft => NULL()
      !> GAPW local atomic grids
      TYPE(hartree_local_type), POINTER                  :: hartree_local => NULL()
      TYPE(local_rho_type), POINTER                      :: local_rho_set => NULL()
      TYPE(local_rho_type), POINTER                      :: local_rho_set_admm => NULL()
   END TYPE tddfpt_subgroup_env_type

! **************************************************************************************************
!> \brief Structure to save global multi-grid related parameters.
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
! **************************************************************************************************
   TYPE mgrid_saved_parameters
      !> create commensurate grids
      LOGICAL                                     :: commensurate_mgrids = .FALSE.
      !> do not perform load balancing
      LOGICAL                                     :: skip_load_balance = .FALSE.
      !> cutoff value at the finest grid level
      REAL(KIND=dp)                               :: cutoff = 0.0_dp
      !> inverse scale factor
      REAL(KIND=dp)                               :: progression_factor = 0.0_dp
      !> relative cutoff
      REAL(KIND=dp)                               :: relative_cutoff = 0.0_dp
      !> list of explicitly given cutoff values
      REAL(KIND=dp), DIMENSION(:), POINTER        :: e_cutoff => NULL()
   END TYPE mgrid_saved_parameters

CONTAINS

! **************************************************************************************************
!> \brief Split MPI communicator to create a set of parallel (sub)groups.
!> \param sub_env  parallel group environment (initialised on exit)
!> \param qs_env   Quickstep environment
!> \param mos_occ  ground state molecular orbitals in primary atomic basis set
!> \param kernel   Type of kernel (full/sTDA) that will be used
!> \par History
!>    * 01.2017 (sub)group-related code has been moved here from the main subroutine tddfpt()
!>              [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_sub_env_init(sub_env, qs_env, mos_occ, kernel)
      TYPE(tddfpt_subgroup_env_type), INTENT(out)        :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(in)         :: mos_occ
      INTEGER, INTENT(in)                                :: kernel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_sub_env_init'

      INTEGER                                            :: handle, ispin, nao, nao_aux, natom, &
                                                            nmo_occ, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mgrid_saved_parameters)                       :: mgrid_saved
      TYPE(mp_para_env_type), POINTER                    :: para_env_global
      TYPE(pw_env_type), POINTER                         :: pw_env_global
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      nspins = SIZE(mos_occ)

      CALL get_qs_env(qs_env, blacs_env=blacs_env_global, dft_control=dft_control, &
                      para_env=para_env_global, pw_env=pw_env_global)

      tddfpt_control => dft_control%tddfpt2_control
      qs_control => dft_control%qs_control

      ! ++ split mpi communicator if
      !    a) the requested number of processors per group > 0
      !       (means that the split has been requested explicitly), and
      !    b) the number of subgroups is >= 2
      sub_env%is_split = tddfpt_control%nprocs > 0 .AND. tddfpt_control%nprocs*2 <= para_env_global%num_pe

      ALLOCATE (sub_env%mos_occ(nspins))
      NULLIFY (sub_env%admm_A)

      IF (sub_env%is_split) THEN
         ALLOCATE (sub_env%group_distribution(0:para_env_global%num_pe - 1))

         ALLOCATE (sub_env%para_env)
         CALL sub_env%para_env%from_split(comm=para_env_global, ngroups=sub_env%ngroups, &
                                          group_distribution=sub_env%group_distribution, subgroup_min_size=tddfpt_control%nprocs)

         ! ++ create a new parallel environment based on the given sub-communicator)
         NULLIFY (sub_env%blacs_env)

         ! use the default (SQUARE) BLACS grid layout and non-repeatable BLACS collective operations
         ! by omitting optional parameters 'blacs_grid_layout' and 'blacs_repeatable'.
         ! Ideally we should take these parameters from the variables globenv%blacs_grid_layout and
         ! globenv%blacs_repeatable, however the global environment is not available
         ! from the subroutine 'qs_energies_properties'.
         CALL cp_blacs_env_create(sub_env%blacs_env, sub_env%para_env)

         NULLIFY (fm_struct)

         DO ispin = 1, nspins
            CALL cp_fm_get_info(mos_occ(ispin), nrow_global=nao, ncol_global=nmo_occ)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=sub_env%blacs_env)
            CALL cp_fm_create(sub_env%mos_occ(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL tddfpt_fm_replicate_across_subgroups(fm_src=mos_occ(ispin), &
                                                      fm_dest_sub=sub_env%mos_occ(ispin), sub_env=sub_env)
         END DO

         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux, ncol_global=nao)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao, context=sub_env%blacs_env)
            ALLOCATE (sub_env%admm_A)
            CALL cp_fm_create(sub_env%admm_A, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL tddfpt_fm_replicate_across_subgroups(fm_src=admm_env%A, fm_dest_sub=sub_env%admm_A, sub_env=sub_env)
         END IF
      ELSE
         CALL para_env_global%retain()
         sub_env%para_env => para_env_global

         CALL blacs_env_global%retain()
         sub_env%blacs_env => blacs_env_global

         sub_env%mos_occ(:) = mos_occ(:)

         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            sub_env%admm_A => admm_env%A
         END IF
      END IF

      IF (kernel == tddfpt_kernel_full) THEN
         ! ++ allocate a new plane wave environment
         sub_env%is_mgrid = sub_env%is_split .OR. tddfpt_control%mgrid_is_explicit

         NULLIFY (sub_env%dbcsr_dist, sub_env%dist_2d)
         NULLIFY (sub_env%sab_orb, sub_env%sab_aux_fit)
         NULLIFY (sub_env%task_list_orb, sub_env%task_list_aux_fit)
         NULLIFY (sub_env%task_list_orb_soft, sub_env%task_list_aux_fit_soft)

         IF (sub_env%is_mgrid) THEN
            IF (tddfpt_control%mgrid_is_explicit) &
               CALL init_tddfpt_mgrid(qs_control, tddfpt_control, mgrid_saved)

            NULLIFY (sub_env%pw_env)

            CALL pw_env_create(sub_env%pw_env)
            CALL pw_env_rebuild(sub_env%pw_env, qs_env, sub_env%para_env)

            CALL tddfpt_build_distribution_2d(distribution_2d=sub_env%dist_2d, dbcsr_dist=sub_env%dbcsr_dist, &
                                              blacs_env=sub_env%blacs_env, qs_env=qs_env)
            CALL tddfpt_build_tasklist(task_list=sub_env%task_list_orb, sab=sub_env%sab_orb, basis_type="ORB", &
                                       distribution_2d=sub_env%dist_2d, pw_env=sub_env%pw_env, qs_env=qs_env, &
                                       soft_valid=.FALSE., skip_load_balance=qs_control%skip_load_balance_distributed, &
                                       reorder_grid_ranks=.TRUE.)
            IF (qs_control%gapw .OR. qs_control%gapw_xc) THEN
               CALL tddfpt_build_tasklist(task_list=sub_env%task_list_orb_soft, sab=sub_env%sab_orb, basis_type="ORB", &
                                          distribution_2d=sub_env%dist_2d, pw_env=sub_env%pw_env, qs_env=qs_env, &
                                          soft_valid=.TRUE., skip_load_balance=qs_control%skip_load_balance_distributed, &
                                          reorder_grid_ranks=.TRUE.)
            END IF

            IF (dft_control%do_admm) THEN
               CALL tddfpt_build_tasklist(task_list=sub_env%task_list_aux_fit, sab=sub_env%sab_aux_fit, &
                                          basis_type="AUX_FIT", distribution_2d=sub_env%dist_2d, &
                                          pw_env=sub_env%pw_env, qs_env=qs_env, soft_valid=.FALSE., &
                                          skip_load_balance=qs_control%skip_load_balance_distributed, &
                                          reorder_grid_ranks=.FALSE.)
               IF (qs_control%gapw .OR. qs_control%gapw_xc) THEN
                  CALL tddfpt_build_tasklist(task_list=sub_env%task_list_aux_fit_soft, sab=sub_env%sab_aux_fit, &
                                             basis_type="AUX_FIT", distribution_2d=sub_env%dist_2d, &
                                             pw_env=sub_env%pw_env, qs_env=qs_env, soft_valid=.TRUE., &
                                             skip_load_balance=qs_control%skip_load_balance_distributed, &
                                             reorder_grid_ranks=.FALSE.)
               END IF
            END IF

            IF (tddfpt_control%mgrid_is_explicit) &
               CALL restore_qs_mgrid(qs_control, mgrid_saved)
         ELSE
            CALL pw_env_retain(pw_env_global)
            sub_env%pw_env => pw_env_global

            CALL get_qs_env(qs_env, dbcsr_dist=sub_env%dbcsr_dist, &
                            sab_orb=sub_env%sab_orb, task_list=sub_env%task_list_orb)
            IF (dft_control%do_admm) THEN
               CALL get_admm_env(admm_env, sab_aux_fit=sub_env%sab_aux_fit, &
                                 task_list_aux_fit=sub_env%task_list_aux_fit)
               IF (qs_control%gapw .OR. qs_control%gapw_xc) THEN
                  sub_env%task_list_aux_fit_soft => admm_env%admm_gapw_env%task_list
               END IF
            END IF
            IF (qs_control%gapw .OR. qs_control%gapw_xc) THEN
               CALL get_qs_env(qs_env, task_list_soft=sub_env%task_list_orb_soft)
            END IF
         END IF

         ! GAPW initializations
         IF (dft_control%qs_control%gapw) THEN
            CALL get_qs_env(qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            natom=natom, &
                            qs_kind_set=qs_kind_set)

            CALL local_rho_set_create(sub_env%local_rho_set)
            CALL allocate_rho_atom_internals(sub_env%local_rho_set%rho_atom_set, atomic_kind_set, &
                                             qs_kind_set, dft_control, sub_env%para_env)

            CALL init_rho0(sub_env%local_rho_set, qs_env, dft_control%qs_control%gapw_control, &
                           zcore=0.0_dp)
            CALL rho0_s_grid_create(sub_env%pw_env, sub_env%local_rho_set%rho0_mpole)
            CALL hartree_local_create(sub_env%hartree_local)
            CALL init_coulomb_local(sub_env%hartree_local, natom)
         ELSEIF (dft_control%qs_control%gapw_xc) THEN
            CALL get_qs_env(qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            qs_kind_set=qs_kind_set)
            CALL local_rho_set_create(sub_env%local_rho_set)
            CALL allocate_rho_atom_internals(sub_env%local_rho_set%rho_atom_set, atomic_kind_set, &
                                             qs_kind_set, dft_control, sub_env%para_env)
         END IF

         ! ADMM/GAPW
         IF (dft_control%do_admm) THEN
            IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
               CALL local_rho_set_create(sub_env%local_rho_set_admm)
               CALL allocate_rho_atom_internals(sub_env%local_rho_set_admm%rho_atom_set, atomic_kind_set, &
                                                admm_env%admm_gapw_env%admm_kind_set, &
                                                dft_control, sub_env%para_env)
            END IF
         END IF

      ELSE IF (kernel == tddfpt_kernel_stda) THEN
         sub_env%is_mgrid = .FALSE.
         NULLIFY (sub_env%dbcsr_dist, sub_env%dist_2d)
         NULLIFY (sub_env%sab_orb, sub_env%sab_aux_fit)
         NULLIFY (sub_env%task_list_orb, sub_env%task_list_orb_soft, sub_env%task_list_aux_fit)
         NULLIFY (sub_env%pw_env)
         IF (sub_env%is_split) THEN
            CPABORT('Subsys option not available')
         ELSE
            CALL get_qs_env(qs_env, dbcsr_dist=sub_env%dbcsr_dist, sab_orb=sub_env%sab_orb)
         END IF
      ELSE IF (kernel == tddfpt_kernel_none) THEN
         sub_env%is_mgrid = .FALSE.
         NULLIFY (sub_env%dbcsr_dist, sub_env%dist_2d)
         NULLIFY (sub_env%sab_orb, sub_env%sab_aux_fit)
         NULLIFY (sub_env%task_list_orb, sub_env%task_list_orb_soft, sub_env%task_list_aux_fit)
         NULLIFY (sub_env%pw_env)
         IF (sub_env%is_split) THEN
            CPABORT('Subsys option not available')
         ELSE
            CALL get_qs_env(qs_env, dbcsr_dist=sub_env%dbcsr_dist, sab_orb=sub_env%sab_orb)
         END IF
      ELSE
         CPABORT("Unknown kernel type")
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_sub_env_init

! **************************************************************************************************
!> \brief Release parallel group environment
!> \param sub_env  parallel group environment (modified on exit)
!> \par History
!>    * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_sub_env_release(sub_env)
      TYPE(tddfpt_subgroup_env_type), INTENT(inout)      :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_sub_env_release'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      IF (sub_env%is_mgrid) THEN
         IF (ASSOCIATED(sub_env%task_list_aux_fit)) &
            CALL deallocate_task_list(sub_env%task_list_aux_fit)

         IF (ASSOCIATED(sub_env%task_list_orb)) &
            CALL deallocate_task_list(sub_env%task_list_orb)

         IF (ASSOCIATED(sub_env%task_list_orb_soft)) &
            CALL deallocate_task_list(sub_env%task_list_orb_soft)

         CALL release_neighbor_list_sets(sub_env%sab_aux_fit)
         CALL release_neighbor_list_sets(sub_env%sab_orb)

         IF (ASSOCIATED(sub_env%dbcsr_dist)) THEN
            CALL dbcsr_distribution_release(sub_env%dbcsr_dist)
            DEALLOCATE (sub_env%dbcsr_dist)
         END IF

         IF (ASSOCIATED(sub_env%dist_2d)) &
            CALL distribution_2d_release(sub_env%dist_2d)
      END IF

      ! GAPW
      IF (ASSOCIATED(sub_env%local_rho_set)) THEN
         CALL local_rho_set_release(sub_env%local_rho_set)
      END IF
      IF (ASSOCIATED(sub_env%hartree_local)) THEN
         CALL hartree_local_release(sub_env%hartree_local)
      END IF
      IF (ASSOCIATED(sub_env%local_rho_set_admm)) THEN
         CALL local_rho_set_release(sub_env%local_rho_set_admm)
      END IF

      ! if TDDFPT-specific plane-wave environment has not been requested,
      ! the pointers sub_env%dbcsr_dist, sub_env%sab_*, and sub_env%task_list_*
      ! point to the corresponding ground-state variables from qs_env
      ! and should not be deallocated

      CALL pw_env_release(sub_env%pw_env)

      sub_env%is_mgrid = .FALSE.

      IF (sub_env%is_split .AND. ASSOCIATED(sub_env%admm_A)) THEN
         CALL cp_fm_release(sub_env%admm_A)
         DEALLOCATE (sub_env%admm_A)
         NULLIFY (sub_env%admm_A)
      END IF

      IF (sub_env%is_split) THEN
         DO i = SIZE(sub_env%mos_occ), 1, -1
            CALL cp_fm_release(sub_env%mos_occ(i))
         END DO
      END IF
      DEALLOCATE (sub_env%mos_occ)

      CALL cp_blacs_env_release(sub_env%blacs_env)
      CALL mp_para_env_release(sub_env%para_env)

      IF (ALLOCATED(sub_env%group_distribution)) &
         DEALLOCATE (sub_env%group_distribution)

      sub_env%is_split = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE tddfpt_sub_env_release

! **************************************************************************************************
!> \brief Replace the global multi-grid related parameters in qs_control by the ones given in the
!>        TDDFPT/MGRID subsection. The original parameters are stored into the 'mgrid_saved'
!>        variable.
!> \param qs_control     Quickstep control parameters (modified on exit)
!> \param tddfpt_control TDDFPT control parameters
!> \param mgrid_saved    structure to hold global MGRID-related parameters (initialised on exit)
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
!> \note the code to build the 'e_cutoff' list was taken from the subroutine read_mgrid_section()
! **************************************************************************************************
   SUBROUTINE init_tddfpt_mgrid(qs_control, tddfpt_control, mgrid_saved)
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(mgrid_saved_parameters), INTENT(out)          :: mgrid_saved

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_tddfpt_mgrid'

      INTEGER                                            :: handle, igrid, ngrids

      CALL timeset(routineN, handle)

      ! ++ save global plane-wave grid parameters to the variable 'mgrid_saved'
      mgrid_saved%commensurate_mgrids = qs_control%commensurate_mgrids
      mgrid_saved%skip_load_balance = qs_control%skip_load_balance_distributed
      mgrid_saved%cutoff = qs_control%cutoff
      mgrid_saved%progression_factor = qs_control%progression_factor
      mgrid_saved%relative_cutoff = qs_control%relative_cutoff
      mgrid_saved%e_cutoff => qs_control%e_cutoff

      ! ++ set parameters from 'tddfpt_control' as default ones for all newly allocated plane-wave grids
      qs_control%commensurate_mgrids = tddfpt_control%mgrid_commensurate_mgrids
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
         ! ++ round 'qs_control%cutoff' upward to the nearest sub-grid's cutoff value;
         !    here we take advantage of the fact that the array 'e_cutoff' has been sorted in descending order
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
            qs_control%e_cutoff(igrid) = qs_control%e_cutoff(igrid - 1)/qs_control%progression_factor
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'restore_qs_mgrid'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(qs_control%e_cutoff)) &
         DEALLOCATE (qs_control%e_cutoff)

      qs_control%commensurate_mgrids = mgrid_saved%commensurate_mgrids
      qs_control%skip_load_balance_distributed = mgrid_saved%skip_load_balance
      qs_control%cutoff = mgrid_saved%cutoff
      qs_control%progression_factor = mgrid_saved%progression_factor
      qs_control%relative_cutoff = mgrid_saved%relative_cutoff
      qs_control%e_cutoff => mgrid_saved%e_cutoff

      CALL timestop(handle)
   END SUBROUTINE restore_qs_mgrid

! **************************************************************************************************
!> \brief Distribute atoms across the two-dimensional grid of processors.
!> \param distribution_2d  new two-dimensional distribution of pairs of particles
!>                         (allocated and initialised on exit)
!> \param dbcsr_dist       new DBCSR distribution (allocated and initialised on exit)
!> \param blacs_env        BLACS parallel environment
!> \param qs_env           Quickstep environment
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_build_distribution_2d(distribution_2d, dbcsr_dist, blacs_env, qs_env)
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_build_distribution_2d'

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
      CALL cp_dbcsr_dist2d_to_dist(distribution_2d, dbcsr_dist)

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
!> \param soft_valid          generate a task list for soft basis functions (GAPW, GAPW_XC)
!> \param skip_load_balance   do not perform load balancing
!> \param reorder_grid_ranks  re-optimise grid ranks and re-create the real-space grid descriptor
!>                            as well as grids
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_build_tasklist(task_list, sab, basis_type, distribution_2d, pw_env, qs_env, &
                                    soft_valid, skip_load_balance, reorder_grid_ranks)
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab
      CHARACTER(len=*), INTENT(in)                       :: basis_type
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: soft_valid, skip_load_balance, &
                                                            reorder_grid_ranks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_build_tasklist'

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
         END IF
      END DO

      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

      NULLIFY (sab)
      CALL section_vals_val_get(input, "DFT%SUBCELLS", r_val=subcells)
      CALL build_neighbor_lists(sab, particle_set, atom2d, cell, pair_radius, &
                                mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb")

      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d, orb_present, orb_radius, pair_radius)

      CALL allocate_task_list(task_list)
      CALL generate_qs_task_list(ks_env, task_list, &
                                 reorder_rs_grid_ranks=reorder_grid_ranks, soft_valid=soft_valid, &
                                 basis_type=basis_type, skip_load_balance_distributed=skip_load_balance, &
                                 pw_env_external=pw_env, sab_orb_external=sab)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_build_tasklist

! **************************************************************************************************
!> \brief Create a DBCSR matrix based on a template matrix, distribution object, and the list of
!>        neighbours.
!> \param matrix      matrix to create
!> \param template    template matrix
!> \param dbcsr_dist  DBCSR distribution
!> \param sab         list of neighbours
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_dbcsr_create_by_dist(matrix, template, dbcsr_dist, sab)
      TYPE(dbcsr_type), POINTER                          :: matrix, template
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_dbcsr_create_by_dist'

      CHARACTER                                          :: matrix_type
      CHARACTER(len=default_string_length)               :: matrix_name
      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(template))
      CALL dbcsr_get_info(template, row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                          name=matrix_name, matrix_type=matrix_type)

      IF (ASSOCIATED(matrix)) THEN
         CALL dbcsr_release(matrix)
      ELSE
         ALLOCATE (matrix)
      END IF

      CALL dbcsr_create(matrix, matrix_name, dbcsr_dist, matrix_type, row_blk_sizes, col_blk_sizes, nze=0)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix, sab)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_dbcsr_create_by_dist

! **************************************************************************************************
!> \brief Replicate a globally distributed matrix across all sub-groups. At the end
!>        every sub-group will hold a local copy of the original globally distributed matrix.
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
!> \param fm_src       globally distributed matrix to replicate
!> \param fm_dest_sub  subgroup-specific copy of the replicated matrix
!> \param sub_env      subgroup environment
!> \par History
!>   * 09.2016 created [Sergey Chulkov]
!>   * 01.2017 moved from qs_tddfpt2_methods [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_fm_replicate_across_subgroups(fm_src, fm_dest_sub, sub_env)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_src, fm_dest_sub
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_fm_replicate_across_subgroups'

      INTEGER :: handle, igroup, igroup_local, ncols_global_dest, ncols_global_src, ngroups, &
         nrows_global_dest, nrows_global_src
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_type)                                   :: fm_null
      TYPE(mp_para_env_type), POINTER                    :: para_env_global

      IF (sub_env%is_split) THEN
         CALL timeset(routineN, handle)

         CALL cp_fm_get_info(fm_src, nrow_global=nrows_global_src, ncol_global=ncols_global_src, &
                             context=blacs_env_global, para_env=para_env_global)
         CALL cp_fm_get_info(fm_dest_sub, nrow_global=nrows_global_dest, ncol_global=ncols_global_dest)

         IF (debug_this_module) THEN
            CPASSERT(nrows_global_src == nrows_global_dest)
            CPASSERT(ncols_global_src == ncols_global_dest)
         END IF

         igroup_local = sub_env%group_distribution(para_env_global%mepos)
         ngroups = sub_env%ngroups

         DO igroup = 0, ngroups - 1
            IF (igroup == igroup_local) THEN
               CALL cp_fm_copy_general(fm_src, fm_dest_sub, para_env_global)
            ELSE
               CALL cp_fm_copy_general(fm_src, fm_null, para_env_global)
            END IF
         END DO

         CALL timestop(handle)
      END IF
   END SUBROUTINE tddfpt_fm_replicate_across_subgroups
END MODULE qs_tddfpt2_subgroups

