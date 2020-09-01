!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      - rewrote collocate for increased accuracy and speed
!>      - introduced the PGI hack for increased speed with that compiler
!>        (22.02.02)
!>      - Added Multiple Grid feature
!>      - new way to get over the grid (01.03.02)
!>      - removed timing calls since they were getting expensive
!>      - Updated with the new QS data structures (09.04.02,MK)
!>      - introduction of the real space grid type ( prelim. version JVdV 05.02)
!>      - parallel FFT (JGH 22.05.02)
!>      - multigrid arrays independent from density (JGH 30.08.02)
!>      - old density stored in g space (JGH 30.08.02)
!>      - distributed real space code (JGH 17.07.03)
!>      - refactoring and new loop ordering (JGH 23.11.03)
!>      - OpenMP parallelization (JGH 03.12.03)
!>      - Modified to compute tau (Joost 12.03)
!>      - removed the incremental density rebuild (Joost 01.04)
!>      - introduced realspace multigridding (Joost 02.04)
!>      - introduced map_consistent (Joost 02.04)
!>      - Addition of the subroutine calculate_atomic_charge_density (TdK, 08.05)
!>      - rewrite of the collocate/integrate kernels (Joost VandeVondele, 03.07)
!> \author Matthias Krack (03.04.2001)
!>      1) Joost VandeVondele (01.2002)
!>      Thomas D. Kuehne (04.08.2005)
!>      Ole Schuett (2020)
! **************************************************************************************************
MODULE qs_collocate_density
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                              gridlevel_info_type
   USE grid_api,                        ONLY: &
        GRID_FUNC_AB, GRID_FUNC_DABpADB_X, GRID_FUNC_DABpADB_Y, GRID_FUNC_DABpADB_Z, &
        GRID_FUNC_DADB, GRID_FUNC_DX, GRID_FUNC_DXDX, GRID_FUNC_DXDY, GRID_FUNC_DY, &
        GRID_FUNC_DYDY, GRID_FUNC_DYDZ, GRID_FUNC_DZ, GRID_FUNC_DZDX, GRID_FUNC_DZDZ, &
        collocate_pgf_product, grid_collocate_task_list
   USE input_constants,                 ONLY: &
        orb_dx2, orb_dxy, orb_dy2, orb_dyz, orb_dz2, orb_dzx, orb_px, orb_py, orb_pz, orb_s
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lri_environment_types,           ONLY: lri_kind_type
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type,&
                                              pw_pools_create_pws,&
                                              pw_pools_give_back_pws
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              realspace_grid_type,&
                                              rs2pw,&
                                              rs_grid_release,&
                                              rs_grid_retain,&
                                              rs_grid_zero,&
                                              rs_pw_transfer
   USE rs_pw_interface,                 ONLY: density_rs2pw,&
                                              density_rs2pw_basic
   USE task_list_methods,               ONLY: rs_copy_matrix,&
                                              rs_distribute_matrix,&
                                              rs_scatter_matrix
   USE task_list_types,                 ONLY: atom_pair_type,&
                                              task_list_type,&
                                              task_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_collocate_density'
! *** Public subroutines ***

   PUBLIC :: calculate_ppl_grid, &
             calculate_rho_core, &
             calculate_lri_rho_elec, &
             calculate_rho_single_gaussian, &
             calculate_rho_metal, &
             calculate_rho_resp_single, &
             calculate_rho_resp_all, &
             calculate_rho_elec, &
             calculate_drho_elec, &
             calculate_wavefunction, &
             calculate_rho_nlcc

   INTEGER :: debug_count = 0

CONTAINS

! **************************************************************************************************
!> \brief computes the density of the non-linear core correction on the grid
!> \param rho_nlcc ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_rho_nlcc(rho_nlcc, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_nlcc
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_nlcc'

      INTEGER                                            :: atom_a, handle, iatom, iexp_nlcc, ikind, &
                                                            ithread, j, n, natom, nc, nexp_nlcc, &
                                                            ni, npme, nthread, subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, nct_nlcc
      LOGICAL                                            :: nlcc
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_nlcc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_nlcc, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pab, particle_set, atomic_kind_set, &
               qs_kind_set, atom_list, pw_env, rs_rho, auxbas_pw_pool, cores)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is chosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, nlcc_present=nlcc, nexp_nlcc=nexp_nlcc, &
                            alpha_nlcc=alpha_nlcc, nct_nlcc=nct_nlcc, cval_nlcc=cval_nlcc)

         IF (.NOT. nlcc) CYCLE

         DO iexp_nlcc = 1, nexp_nlcc

            alpha = alpha_nlcc(iexp_nlcc)
            nc = nct_nlcc(iexp_nlcc)

            ni = ncoset(2*nc - 2)
            ALLOCATE (pab(ni, 1))
            pab = 0._dp

            nthread = 1
            ithread = 0

            CALL reallocate(cores, 1, natom)
            npme = 0
            cores = 0

            ! prepare core function
            DO j = 1, nc
               SELECT CASE (j)
               CASE (1)
                  pab(1, 1) = cval_nlcc(1, iexp_nlcc)
               CASE (2)
                  n = coset(2, 0, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 2, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 0, 2)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
               CASE (3)
                  n = coset(4, 0, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 4, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 0, 4)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 2, 0)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 0, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 2, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
               CASE (4)
                  n = coset(6, 0, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 6, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 0, 6)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 2, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 0, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 4, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 0, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 4, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 2, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 2, 2)
                  pab(n, 1) = 6._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
               CASE DEFAULT
                  CPABORT("")
               END SELECT
            END DO
            IF (dft_control%nspins == 2) pab = pab*0.5_dp

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
                  ! replicated realspace grid, split the atoms up between procs
                  IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                     npme = npme + 1
                     cores(npme) = iatom
                  ENDIF
               ELSE
                  npme = npme + 1
                  cores(npme) = iatom
               ENDIF
            END DO

            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               subpatch_pattern = 0
               ni = 2*nc - 2
               radius = exp_radius_very_extended(la_min=0, la_max=ni, lb_min=0, lb_max=0, &
                                                 ra=ra, rb=ra, rp=ra, &
                                                 zetp=1/(2*alpha**2), eps=eps_rho_rspace, &
                                                 pab=pab, o1=0, o2=0, &  ! without map_consistent
                                                 prefactor=1.0_dp, cutoff=0.0_dp)

               CALL collocate_pgf_product(ni, 1/(2*alpha**2), 0, 0, 0.0_dp, 0, ra, &
                                          (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, 0, 0, rs_rho, &
                                          cell, cube_info, ga_gb_function=GRID_FUNC_AB, radius=radius, &
                                          use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)

            END DO

            DEALLOCATE (pab)

         END DO

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      CALL rs_pw_transfer(rs_rho, rho_nlcc%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_nlcc

! **************************************************************************************************
!> \brief computes the local pseudopotential (without erf term) on the grid
!> \param vppl ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_ppl_grid(vppl, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: vppl
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ppl_grid'

      INTEGER                                            :: atom_a, handle, iatom, ikind, ithread, &
                                                            j, lppl, n, natom, ni, npme, nthread, &
                                                            subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cexp_ppl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pab, atomic_kind_set, qs_kind_set, particle_set, &
               atom_list, pw_env, rs_rho, auxbas_pw_pool, cores)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is chosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, alpha_ppl=alpha, nexp_ppl=lppl, cexp_ppl=cexp_ppl)

         IF (lppl <= 0) CYCLE

         ni = ncoset(2*lppl - 2)
         ALLOCATE (pab(ni, 1))
         pab = 0._dp

         nthread = 1
         ithread = 0

         CALL reallocate(cores, 1, natom)
         npme = 0
         cores = 0

         ! prepare core function
         DO j = 1, lppl
            SELECT CASE (j)
            CASE (1)
               pab(1, 1) = cexp_ppl(1)
            CASE (2)
               n = coset(2, 0, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 2, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 0, 2)
               pab(n, 1) = cexp_ppl(2)
            CASE (3)
               n = coset(4, 0, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 4, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 0, 4)
               pab(n, 1) = cexp_ppl(3)
               n = coset(2, 2, 0)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(2, 0, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(0, 2, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
            CASE (4)
               n = coset(6, 0, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 6, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 0, 6)
               pab(n, 1) = cexp_ppl(4)
               n = coset(4, 2, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(4, 0, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 4, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 0, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 4, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 2, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 2, 2)
               pab(n, 1) = 6._dp*cexp_ppl(4)
            CASE DEFAULT
               CPABORT("")
            END SELECT
         END DO

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iatom
               ENDIF
            ELSE
               npme = npme + 1
               cores(npme) = iatom
            ENDIF
         END DO

         IF (npme .GT. 0) THEN
            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               subpatch_pattern = 0
               ni = 2*lppl - 2

               radius = exp_radius_very_extended(la_min=0, la_max=ni, &
                                                 lb_min=0, lb_max=0, &
                                                 ra=ra, rb=ra, rp=ra, &
                                                 zetp=alpha, eps=eps_rho_rspace, &
                                                 pab=pab, o1=0, o2=0, &  ! without map_consistent
                                                 prefactor=1.0_dp, cutoff=0.0_dp)

               CALL collocate_pgf_product(ni, alpha, 0, 0, 0.0_dp, 0, ra, &
                                          (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, 0, 0, rs_rho, &
                                          cell, cube_info, radius=radius, ga_gb_function=GRID_FUNC_AB, &
                                          use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)

            END DO
         ENDIF

         DEALLOCATE (pab)

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      CALL rs_pw_transfer(rs_rho, vppl%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL timestop(handle)

   END SUBROUTINE calculate_ppl_grid

! **************************************************************************************************
!> \brief Collocates the fitted lri density on a grid.
!> \param lri_rho_g ...
!> \param lri_rho_r ...
!> \param qs_env ...
!> \param lri_coef ...
!> \param total_rho ...
!> \param basis_type ...
!> \param exact_1c_terms ...
!> \param pmat replicated block diagonal density matrix (optional)
!> \param atomlist list of atoms to be included (optional)
!> \par History
!>      04.2013
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_lri_rho_elec(lri_rho_g, lri_rho_r, qs_env, &
                                     lri_coef, total_rho, basis_type, exact_1c_terms, pmat, atomlist)

      TYPE(pw_p_type), INTENT(INOUT)                     :: lri_rho_g, lri_rho_r
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      CHARACTER(len=*), INTENT(IN)                       :: basis_type
      LOGICAL, INTENT(IN)                                :: exact_1c_terms
      TYPE(dbcsr_type), OPTIONAL                         :: pmat
      INTEGER, DIMENSION(:), OPTIONAL                    :: atomlist

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_lri_rho_elec'

      INTEGER :: atom_a, dir, group_size, handle, iatom, igrid_level, ikind, ipgf, iset, jpgf, &
         jset, m1, maxco, maxsgf_set, my_pos, na1, natom, nb1, ncoa, ncob, nseta, offset, sgfa, &
         sgfb
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: map_it
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: map_it2
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aci
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, pab, sphi_a, work, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: lri_basis_set, orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: mgrid_gspace, mgrid_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(realspace_grid_type), POINTER                 :: rs_grid

      NULLIFY (aci, atomic_kind_set, qs_kind_set, atom_list, cell, &
               cube_info, dft_control, first_sgfa, gridlevel_info, la_max, &
               la_min, lri_basis_set, mgrid_gspace, mgrid_rspace, npgfa, nsgfa, &
               pab, particle_set, pw_env, pw_pools, rs_grid, rs_rho, sphi_a, &
               work, zeta)

      CALL timeset(routineN, handle)

      IF (exact_1c_terms) THEN
         CPASSERT(PRESENT(pmat))
      END IF

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, particle_set=particle_set, &
                      pw_env=pw_env, &
                      dft_control=dft_control)

      cube_info => pw_env%cube_info
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      gridlevel_info => pw_env%gridlevel_info

      ! *** set up the pw multi-grids *** !
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env=pw_env, rs_grids=rs_rho, pw_pools=pw_pools)

      CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                               use_data=REALDATA3D, &
                               in_space=REALSPACE)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                               use_data=COMPLEXDATA1D, &
                               in_space=RECIPROCALSPACE)

      ! *** set up the rs multi-grids *** !
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      !take maxco from the LRI basis set!
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, basis_type=basis_type)

      ALLOCATE (pab(maxco, 1))
      offset = 0
      my_pos = mgrid_rspace(1)%pw%pw_grid%para%my_pos
      group_size = mgrid_rspace(1)%pw%pw_grid%para%group_size

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=lri_basis_set, basis_type=basis_type)

         !Take the lri basis set here!
         CALL get_gto_basis_set(gto_basis_set=lri_basis_set, lmax=la_max, &
                                lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            IF (PRESENT(ATOMLIST)) THEN
               IF (atomlist(atom_a) == 0) CYCLE
            END IF
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            aci => lri_coef(ikind)%acoef(iatom, :)

            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it(m1))
            DO iset = 1, nseta
               ! collocate this set locally?
               map_it = .FALSE.
               DO ipgf = 1, npgfa(iset)
                  igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                  rs_grid => rs_rho(igrid_level)%rs_grid
                  IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                     DO dir = 1, 3
                        ! bounds of local grid (i.e. removing the 'wings'), if periodic
                        tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                        tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                        IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                           lb(dir) = rs_grid%lb_local(dir) + rs_grid%desc%border
                           ub(dir) = rs_grid%ub_local(dir) - rs_grid%desc%border
                        ELSE
                           lb(dir) = rs_grid%lb_local(dir)
                           ub(dir) = rs_grid%ub_local(dir)
                        ENDIF
                        ! distributed grid, only map if it is local to the grid
                        location(dir) = tp(dir) + rs_grid%desc%lb(dir)
                     ENDDO
                     IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                         lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                         lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                        map_it(ipgf) = .TRUE.
                     ENDIF
                  ELSE
                     ! not distributed, just a round-robin distribution over the full set of CPUs
                     IF (MODULO(offset, group_size) == my_pos) map_it(ipgf) = .TRUE.
                  ENDIF
               END DO
               offset = offset + 1

               IF (ANY(map_it(1:npgfa(iset)))) THEN
                  sgfa = first_sgfa(1, iset)
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  m1 = sgfa + nsgfa(iset) - 1
                  ALLOCATE (work(nsgfa(iset), 1))
                  work(1:nsgfa(iset), 1) = aci(sgfa:m1)
                  pab = 0._dp

                  CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), 1.0_dp, lri_basis_set%sphi(1, sgfa), &
                             SIZE(lri_basis_set%sphi, 1), work(1, 1), SIZE(work, 1), 0.0_dp, pab(1, 1), &
                             SIZE(pab, 1))

                  DO ipgf = 1, npgfa(iset)
                     na1 = (ipgf - 1)*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                     rs_grid => rs_rho(igrid_level)%rs_grid
                     IF (map_it(ipgf)) THEN
                        radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                          lb_min=0, lb_max=0, &
                                                          ra=ra, rb=ra, rp=ra, &
                                                          zetp=zeta(ipgf, iset), eps=eps_rho_rspace, &
                                                          prefactor=1.0_dp, cutoff=1.0_dp)

                        CALL collocate_pgf_product(la_max=la_max(iset), &
                                                   zeta=zeta(ipgf, iset), &
                                                   la_min=la_min(iset), &
                                                   lb_max=0, zetb=0.0_dp, lb_min=0, &
                                                   ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                                   scale=1._dp, &
                                                   pab=pab, o1=na1, o2=0, &
                                                   rsgrid=rs_grid, cell=cell, &
                                                   cube_info=cube_info(igrid_level), &
                                                   radius=radius, &
                                                   ga_gb_function=GRID_FUNC_AB)
                     ENDIF
                  END DO
                  DEALLOCATE (work)
               END IF
            END DO
            DEALLOCATE (map_it)
         END DO
      END DO

      DEALLOCATE (pab)

      ! process the one-center terms
      IF (exact_1c_terms) THEN
         ! find maximum numbers
         offset = 0
         CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                              maxco=maxco, &
                              maxsgf_set=maxsgf_set, &
                              basis_type="ORB")
         ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set))

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, lmax=la_max, &
                                   lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                   sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)
            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               CALL dbcsr_get_block_p(matrix=pmat, row=atom_a, col=atom_a, BLOCK=p_block, found=found)
               m1 = MAXVAL(npgfa(1:nseta))
               ALLOCATE (map_it2(m1, m1))
               DO iset = 1, nseta
                  DO jset = 1, nseta
                     ! processor mappint
                     map_it2 = .FALSE.
                     DO ipgf = 1, npgfa(iset)
                        DO jpgf = 1, npgfa(jset)
                           zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                           igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                           rs_grid => rs_rho(igrid_level)%rs_grid
                           IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                              DO dir = 1, 3
                                 ! bounds of local grid (i.e. removing the 'wings'), if periodic
                                 tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                                 tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                                 IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                                    lb(dir) = rs_grid%lb_local(dir) + rs_grid%desc%border
                                    ub(dir) = rs_grid%ub_local(dir) - rs_grid%desc%border
                                 ELSE
                                    lb(dir) = rs_grid%lb_local(dir)
                                    ub(dir) = rs_grid%ub_local(dir)
                                 ENDIF
                                 ! distributed grid, only map if it is local to the grid
                                 location(dir) = tp(dir) + rs_grid%desc%lb(dir)
                              ENDDO
                              IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                                  lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                                  lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                                 map_it2(ipgf, jpgf) = .TRUE.
                              ENDIF
                           ELSE
                              ! not distributed, just a round-robin distribution over the full set of CPUs
                              IF (MODULO(offset, group_size) == my_pos) map_it2(ipgf, jpgf) = .TRUE.
                           ENDIF
                        END DO
                     END DO
                     offset = offset + 1
                     !
                     IF (ANY(map_it2(1:npgfa(iset), 1:npgfa(jset)))) THEN
                        ncoa = npgfa(iset)*ncoset(la_max(iset))
                        sgfa = first_sgfa(1, iset)
                        ncob = npgfa(jset)*ncoset(la_max(jset))
                        sgfb = first_sgfa(1, jset)
                        ! decontract density block
                        CALL dgemm("N", "N", ncoa, nsgfa(jset), nsgfa(iset), &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   p_block(sgfa, sgfb), SIZE(p_block, 1), &
                                   0.0_dp, work(1, 1), maxco)
                        CALL dgemm("N", "T", ncoa, ncob, nsgfa(jset), &
                                   1.0_dp, work(1, 1), maxco, &
                                   sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                   0.0_dp, pab(1, 1), maxco)
                        DO ipgf = 1, npgfa(iset)
                           DO jpgf = 1, npgfa(jset)
                              zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                              igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                              rs_grid => rs_rho(igrid_level)%rs_grid

                              na1 = (ipgf - 1)*ncoset(la_max(iset))
                              nb1 = (jpgf - 1)*ncoset(la_max(jset))

                              IF (map_it2(ipgf, jpgf)) THEN
                                 radius = exp_radius_very_extended(la_min=la_min(iset), &
                                                                   la_max=la_max(iset), &
                                                                   lb_min=la_min(jset), &
                                                                   lb_max=la_max(jset), &
                                                                   ra=ra, rb=ra, rp=ra, &
                                                                   zetp=zetp, eps=eps_rho_rspace, &
                                                                   prefactor=1.0_dp, cutoff=1.0_dp)

                                 CALL collocate_pgf_product( &
                                    la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                    la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                    ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, na1, nb1, &
                                    rs_grid, cell, cube_info(igrid_level), &
                                    radius=radius, ga_gb_function=GRID_FUNC_AB)
                              END IF
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
               DEALLOCATE (map_it2)
               !
            END DO
         END DO
         DEALLOCATE (pab, work)
      END IF

      CALL pw_zero(lri_rho_g%pw)
      CALL pw_zero(lri_rho_r%pw)

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_rspace(igrid_level)%pw)
         CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                             pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid)
      END DO

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_gspace(igrid_level)%pw)
         CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                          mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw, lri_rho_g%pw)
      END DO
      CALL pw_transfer(lri_rho_g%pw, lri_rho_r%pw)
      total_rho = pw_integrate_function(lri_rho_r%pw, isign=-1)

      ! *** give back the multi-grids *** !
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_rho_elec

! **************************************************************************************************
!> \brief computes the density of the core charges on the grid
!> \param rho_core ...
!> \param total_rho ...
!> \param qs_env ...
!> \param only_nopaw ...
! **************************************************************************************************
   SUBROUTINE calculate_rho_core(rho_core, total_rho, qs_env, only_nopaw)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_core
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_nopaw

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_core'

      INTEGER                                            :: atom_a, handle, iatom, ikind, ithread, &
                                                            j, natom, npme, nthread, &
                                                            subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: my_only_nopaw, paw_atom
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pab, atomic_kind_set, qs_kind_set, particle_set, &
               atom_list, pw_env, rs_rho, auxbas_pw_pool, cores)
      ALLOCATE (pab(1, 1))

      my_only_nopaw = .FALSE.
      IF (PRESENT(only_nopaw)) my_only_nopaw = only_nopaw

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is chosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom, &
                          alpha_core_charge=alpha, ccore_charge=pab(1, 1))

         IF (my_only_nopaw .AND. paw_atom) CYCLE
         IF (alpha == 0.0_dp .OR. pab(1, 1) == 0.0_dp) CYCLE

         nthread = 1
         ithread = 0

         CALL reallocate(cores, 1, natom)
         npme = 0
         cores = 0

         DO iatom = 1, natom
            IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iatom
               ENDIF
            ELSE
               npme = npme + 1
               cores(npme) = iatom
            ENDIF
         END DO

         IF (npme .GT. 0) THEN
            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               subpatch_pattern = 0
               radius = exp_radius_very_extended(la_min=0, la_max=0, &
                                                 lb_min=0, lb_max=0, &
                                                 ra=ra, rb=ra, rp=ra, &
                                                 zetp=alpha, eps=eps_rho_rspace, &
                                                 pab=pab, o1=0, o2=0, &  ! without map_consistent
                                                 prefactor=-1.0_dp, cutoff=0.0_dp)

               CALL collocate_pgf_product(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                          (/0.0_dp, 0.0_dp, 0.0_dp/), -1.0_dp, pab, 0, 0, rs_rho, &
                                          cell, cube_info, radius=radius, ga_gb_function=GRID_FUNC_AB, &
                                          use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)

            END DO
         ENDIF

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF
      DEALLOCATE (pab)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      total_rho = pw_integrate_function(rhoc_r%pw, isign=-1)

      CALL pw_transfer(rhoc_r%pw, rho_core%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_core

! **************************************************************************************************
!> \brief collocate a single Gaussian on the grid
!> \param rho_gb charge density generated by a single gaussian
!> \param qs_env qs environment
!> \param iatom_in atom index
!> \par History
!>        12.2011 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_single_gaussian(rho_gb, qs_env, iatom_in)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_gb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom_in

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_single_gaussian'

      INTEGER                                            :: atom_a, handle, iatom, npme, &
                                                            subpatch_pattern
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pab, pw_env, rs_rho, auxbas_pw_pool)

      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      pab(1, 1) = 1.0_dp
      iatom = iatom_in

      npme = 0

      IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
         IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
            npme = npme + 1
         ENDIF
      ELSE
         npme = npme + 1
      ENDIF

      IF (npme .GT. 0) THEN
         atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
         ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r, cell)
         subpatch_pattern = 0
         radius = exp_radius_very_extended(la_min=0, la_max=0, &
                                           lb_min=0, lb_max=0, &
                                           ra=ra, rb=ra, rp=ra, &
                                           zetp=qs_env%qmmm_env_qm%image_charge_pot%eta, &
                                           eps=eps_rho_rspace, &
                                           pab=pab, o1=0, o2=0, &  ! without map_consistent
                                           prefactor=1.0_dp, cutoff=0.0_dp)

         CALL collocate_pgf_product(0, qs_env%qmmm_env_qm%image_charge_pot%eta, &
                                    0, 0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, 0, 0, rs_rho, &
                                    cell, pw_env%cube_info(1), radius=radius, ga_gb_function=GRID_FUNC_AB, &
                                    use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)
      ENDIF

      DEALLOCATE (pab)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_gb%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_single_gaussian

! **************************************************************************************************
!> \brief computes the image charge density on the grid (including coeffcients)
!> \param rho_metal image charge density
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param total_rho_metal total induced image charge density
!> \param qs_env qs environment
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_metal(rho_metal, coeff, total_rho_metal, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_metal
      REAL(KIND=dp), DIMENSION(:), POINTER               :: coeff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: total_rho_metal
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_metal'

      INTEGER                                            :: atom_a, handle, iatom, j, natom, npme, &
                                                            subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: cores
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pab, pw_env, rs_rho, auxbas_pw_pool, cores)

      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      pab(1, 1) = 1.0_dp

      natom = SIZE(qs_env%qmmm_env_qm%image_charge_pot%image_mm_list)

      CALL reallocate(cores, 1, natom)
      npme = 0
      cores = 0

      DO iatom = 1, natom
         IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
               npme = npme + 1
               cores(npme) = iatom
            ENDIF
         ELSE
            npme = npme + 1
            cores(npme) = iatom
         ENDIF
      ENDDO

      IF (npme .GT. 0) THEN
         DO j = 1, npme
            iatom = cores(j)
            atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
            ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r, cell)
            subpatch_pattern = 0
            radius = exp_radius_very_extended(la_min=0, la_max=0, &
                                              lb_min=0, lb_max=0, &
                                              ra=ra, rb=ra, rp=ra, &
                                              zetp=qs_env%qmmm_env_qm%image_charge_pot%eta, &
                                              eps=eps_rho_rspace, &
                                              pab=pab, o1=0, o2=0, &  ! without map_consistent
                                              prefactor=coeff(iatom), cutoff=0.0_dp)

            CALL collocate_pgf_product( &
               0, qs_env%qmmm_env_qm%image_charge_pot%eta, &
               0, 0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), coeff(iatom), pab, 0, 0, rs_rho, &
               cell, pw_env%cube_info(1), radius=radius, ga_gb_function=GRID_FUNC_AB, &
               use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)
         ENDDO
      ENDIF

      DEALLOCATE (pab, cores)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      IF (PRESENT(total_rho_metal)) &
         !minus sign: account for the fact that rho_metal has opposite sign
         total_rho_metal = pw_integrate_function(rhoc_r%pw, isign=-1)

      CALL pw_transfer(rhoc_r%pw, rho_metal%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_metal

! **************************************************************************************************
!> \brief collocate a single Gaussian on the grid for periodic RESP fitting
!> \param rho_gb charge density generated by a single gaussian
!> \param qs_env qs environment
!> \param eta width of single Gaussian
!> \param iatom_in atom index
!> \par History
!>        06.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_resp_single(rho_gb, qs_env, eta, iatom_in)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_gb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: eta
      INTEGER, INTENT(IN)                                :: iatom_in

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_single'

      INTEGER                                            :: handle, iatom, npme, subpatch_pattern
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pab, pw_env, rs_rho, auxbas_pw_pool, &
               particle_set)

      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      pab(1, 1) = 1.0_dp
      iatom = iatom_in

      npme = 0

      IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
         IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
            npme = npme + 1
         ENDIF
      ELSE
         npme = npme + 1
      ENDIF

      IF (npme .GT. 0) THEN
         ra(:) = pbc(particle_set(iatom)%r, cell)
         subpatch_pattern = 0
         radius = exp_radius_very_extended(la_min=0, la_max=0, &
                                           lb_min=0, lb_max=0, &
                                           ra=ra, rb=ra, rp=ra, &
                                           zetp=eta, eps=eps_rho_rspace, &
                                           pab=pab, o1=0, o2=0, &  ! without map_consistent
                                           prefactor=1.0_dp, cutoff=0.0_dp)

         CALL collocate_pgf_product(0, eta, 0, 0, 0.0_dp, 0, ra, &
                                    (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, 0, 0, rs_rho, &
                                    cell, pw_env%cube_info(1), radius=radius, ga_gb_function=GRID_FUNC_AB, &
                                    use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)
      ENDIF

      DEALLOCATE (pab)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_gb%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_resp_single

! **************************************************************************************************
!> \brief computes the RESP charge density on a grid based on the RESP charges
!> \param rho_resp RESP charge density
!> \param coeff RESP charges, take care of normalization factor
!>        (eta/pi)**1.5 later
!> \param natom number of atoms
!> \param eta width of single Gaussian
!> \param qs_env qs environment
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_resp_all(rho_resp, coeff, natom, eta, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_resp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: coeff
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(IN)                          :: eta
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_all'

      INTEGER                                            :: handle, iatom, j, npme, subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: cores
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, cores, dft_control, pab, pw_env, rs_rho, auxbas_pw_pool, &
               particle_set)

      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      pab(1, 1) = 1.0_dp

      CALL reallocate(cores, 1, natom)
      npme = 0
      cores = 0

      DO iatom = 1, natom
         IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
               npme = npme + 1
               cores(npme) = iatom
            ENDIF
         ELSE
            npme = npme + 1
            cores(npme) = iatom
         ENDIF
      ENDDO

      IF (npme .GT. 0) THEN
         DO j = 1, npme
            iatom = cores(j)
            ra(:) = pbc(particle_set(iatom)%r, cell)
            subpatch_pattern = 0
            radius = exp_radius_very_extended(la_min=0, la_max=0, &
                                              lb_min=0, lb_max=0, &
                                              ra=ra, rb=ra, rp=ra, &
                                              zetp=eta, eps=eps_rho_rspace, &
                                              pab=pab, o1=0, o2=0, &  ! without map_consistent
                                              prefactor=coeff(iatom), cutoff=0.0_dp)

            CALL collocate_pgf_product( &
               0, eta, &
               0, 0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), coeff(iatom), pab, 0, 0, rs_rho, &
               cell, pw_env%cube_info(1), radius=radius, ga_gb_function=GRID_FUNC_AB, &
               use_subpatch=.TRUE., subpatch_pattern=subpatch_pattern)
         ENDDO
      ENDIF

      DEALLOCATE (pab, cores)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_resp%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_resp_all

! **************************************************************************************************
!> \brief computes the density corresponding to a given density matrix on the grid
!> \param matrix_p ...
!> \param matrix_p_kp ...
!> \param rho ...
!> \param rho_gspace ...
!> \param total_rho ...
!> \param ks_env ...
!> \param soft_valid ...
!> \param compute_tau ...
!> \param compute_grad ...
!> \param basis_type ...
!> \param der_type ...
!> \param idir ...
!> \param task_list_external ...
!> \param pw_env_external ...
!> \par History
!>      IAB (15-Feb-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      Anything that is not the default ORB basis_type requires an external_task_list 12.2019, (A.Bussy)
!>      Ole Schuett (2020): Migrated to C, see grid_api.F
!> \note
!>      both rho and rho_gspace contain the new rho
!>      (in real and g-space respectively)
! **************************************************************************************************
   SUBROUTINE calculate_rho_elec(matrix_p, matrix_p_kp, rho, rho_gspace, total_rho, &
                                 ks_env, soft_valid, compute_tau, compute_grad, &
                                 basis_type, der_type, idir, task_list_external, pw_env_external)

      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_p_kp
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho, rho_gspace
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: soft_valid, compute_tau, compute_grad
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: der_type, idir
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_elec'

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER                                            :: ga_gb_function, group, handle, ilevel, &
                                                            img, nimages, nlevels
      LOGICAL                                            :: any_distributed, my_compute_grad, &
                                                            my_compute_tau, my_soft_valid
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_images
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (matrix_images, dft_control, pw_env, rs_rho, task_list)

      ! Figure out which function to collocate.
      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
      my_compute_grad = .FALSE.
      IF (PRESENT(compute_grad)) my_compute_grad = compute_grad
      IF (PRESENT(der_type)) THEN
         SELECT CASE (der_type)
         CASE (orb_s)
            ga_gb_function = GRID_FUNC_AB
         CASE (orb_px)
            ga_gb_function = GRID_FUNC_DX
         CASE (orb_py)
            ga_gb_function = GRID_FUNC_DY
         CASE (orb_pz)
            ga_gb_function = GRID_FUNC_DZ
         CASE (orb_dxy)
            ga_gb_function = GRID_FUNC_DXDY
         CASE (orb_dyz)
            ga_gb_function = GRID_FUNC_DYDZ
         CASE (orb_dzx)
            ga_gb_function = GRID_FUNC_DZDX
         CASE (orb_dx2)
            ga_gb_function = GRID_FUNC_DXDX
         CASE (orb_dy2)
            ga_gb_function = GRID_FUNC_DYDY
         CASE (orb_dz2)
            ga_gb_function = GRID_FUNC_DZDZ
         CASE DEFAULT
            CPABORT("Unknown der_type")
         END SELECT
      ELSE IF (my_compute_tau) THEN
         ga_gb_function = GRID_FUNC_DADB
      ELSE IF (my_compute_grad) THEN
         CPASSERT(PRESENT(idir))
         SELECT CASE (idir)
         CASE (1)
            ga_gb_function = GRID_FUNC_DABpADB_X
         CASE (2)
            ga_gb_function = GRID_FUNC_DABpADB_Y
         CASE (3)
            ga_gb_function = GRID_FUNC_DABpADB_Z
         CASE DEFAULT
            CPABORT("invalid idir")
         END SELECT
      ELSE
         ga_gb_function = GRID_FUNC_AB
      ENDIF

      ! Figure out which basis_type to use.
      my_basis_type = "ORB"  ! by default, the full density is calculated
      IF (PRESENT(basis_type)) my_basis_type = basis_type
      CPASSERT(my_basis_type == "ORB" .OR. PRESENT(task_list_external))

      ! Figure out which task_list to use.
      my_soft_valid = .FALSE.
      IF (PRESENT(soft_valid)) my_soft_valid = soft_valid
      IF (my_soft_valid) THEN
         !TODO: Why does soft_valid overrule task_list_external?
         CALL get_ks_env(ks_env, task_list_soft=task_list)
      ELSE IF (PRESENT(task_list_external)) THEN
         task_list => task_list_external
      ELSE
         CALL get_ks_env(ks_env, task_list=task_list)
      ENDIF
      CPASSERT(ASSOCIATED(task_list))

      ! Figure out which pw_env to use.
      IF (PRESENT(pw_env_external)) THEN
         pw_env => pw_env_external
      ELSE
         CALL get_ks_env(ks_env, pw_env=pw_env)
      ENDIF
      CPASSERT(ASSOCIATED(pw_env))

      ! Get grids.
      CALL pw_env_get(pw_env, rs_grids=rs_rho)
      nlevels = SIZE(rs_rho)
      group = rs_rho(1)%rs_grid%desc%group

      ! Zero all grids.
      DO ilevel = 1, nlevels
         CALL rs_grid_retain(rs_rho(ilevel)%rs_grid)
         CALL rs_grid_zero(rs_rho(ilevel)%rs_grid)
      END DO

      ! Check if any of the grids is distributed.
      any_distributed = .FALSE.
      DO ilevel = 1, nlevels
         any_distributed = any_distributed .OR. rs_rho(ilevel)%rs_grid%desc%distributed
      ENDDO

      ! Gather all matrix images in a single array.
      CALL get_ks_env(ks_env, dft_control=dft_control)
      nimages = dft_control%nimages
      ALLOCATE (matrix_images(nimages))
      IF (PRESENT(matrix_p_kp)) THEN
         CPASSERT(.NOT. PRESENT(matrix_p))
         DO img = 1, nimages
            matrix_images(img)%matrix => matrix_p_kp(img)%matrix
         END DO
      ELSE
         CPASSERT(PRESENT(matrix_p) .AND. nimages == 1)
         matrix_images(1)%matrix => matrix_p
      END IF

      ! Distribute matrix blocks.
      IF (any_distributed) THEN
         CALL rs_scatter_matrix(task_list, group, matrix_images)
      ELSE
         CALL rs_copy_matrix(task_list, matrix_images)
      ENDIF
      DEALLOCATE (matrix_images)

      ! Map all tasks onto the grids
      CALL grid_collocate_task_list(task_list=task_list%grid_task_list, &
                                    ga_gb_function=ga_gb_function, &
                                    rs_grids=rs_rho)

      ! Merge realspace multi-grids into single planewave grid.
      CALL density_rs2pw(pw_env, rs_rho, rho, rho_gspace)
      total_rho = pw_integrate_function(rho%pw, isign=-1)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_elec

! **************************************************************************************************
!> \brief computes the gradient of the density corresponding to a given
!>        density matrix on the grid
!> \param matrix_p ...
!> \param matrix_p_kp ...
!> \param drho ...
!> \param drho_gspace ...
!> \param qs_env ...
!> \param soft_valid ...
!> \param basis_type ...
!> \note  this is an alternative to calculate the gradient through FFTs
! **************************************************************************************************
   SUBROUTINE calculate_drho_elec(matrix_p, matrix_p_kp, drho, drho_gspace, qs_env, &
                                  soft_valid, basis_type)

      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_p_kp
      TYPE(pw_p_type), DIMENSION(:), INTENT(INOUT)       :: drho, drho_gspace
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: soft_valid
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_drho_elec'

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER :: bcol, brow, dabqadb_func, handle, i, iatom, iatom_old, idir, igrid_level, ikind, &
         ikind_old, img, img_old, ipgf, iset, iset_old, itask, ithread, jatom, jatom_old, jkind, &
         jkind_old, jpgf, jset, jset_old, maxco, maxsgf_set, na1, na2, natoms, nb1, nb2, ncoa, &
         ncob, nimages, nseta, nsetb, ntasks, nthread, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: atom_pair_changed, distributed_rs_grids, &
                                                            do_kp, found, my_soft, use_subpatch
      REAL(KIND=dp)                                      :: eps_rho_rspace, f, prefactor, radius, &
                                                            scale, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rab_inv, rb, rp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, pab, sphi_a, sphi_b, work, &
                                                            zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: pabt, workt
      TYPE(atom_pair_type), DIMENSION(:), POINTER        :: atom_pair_recv, atom_pair_send
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(task_type), DIMENSION(:), POINTER             :: tasks

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(matrix_p) .OR. PRESENT(matrix_p_kp))
      do_kp = PRESENT(matrix_p_kp)

      NULLIFY (cell, dft_control, orb_basis_set, deltap, qs_kind_set, &
               sab_orb, particle_set, rs_rho, pw_env, rs_descs, la_max, la_min, &
               lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, p_block, sphi_a, &
               sphi_b, zeta, zetb, first_sgfa, first_sgfb, tasks, pabt, workt)

      debug_count = debug_count + 1

      ! by default, the full density is calculated
      my_soft = .FALSE.
      IF (PRESENT(soft_valid)) my_soft = soft_valid

      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      sab_orb=sab_orb, &
                      pw_env=pw_env)

      SELECT CASE (my_basis_type)
      CASE ("ORB")
         CALL get_qs_env(qs_env=qs_env, &
                         task_list=task_list, &
                         task_list_soft=task_list_soft)
      CASE ("AUX_FIT")
         CALL get_qs_env(qs_env=qs_env, &
                         task_list_aux_fit=task_list, &
                         task_list_soft=task_list_soft)
      END SELECT

      ! *** assign from pw_env
      gridlevel_info => pw_env%gridlevel_info
      cube_info => pw_env%cube_info

      !   *** Allocate work storage ***
      nthread = 1
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)
      CALL reallocate(pabt, 1, maxco, 1, maxco, 0, nthread - 1)
      CALL reallocate(workt, 1, maxco, 1, maxsgf_set, 0, nthread - 1)

      ! find maximum numbers
      nimages = dft_control%nimages
      CPASSERT(nimages == 1 .OR. do_kp)

      natoms = SIZE(particle_set)

      ! get the task lists
      IF (my_soft) task_list => task_list_soft
      CPASSERT(ASSOCIATED(task_list))
      tasks => task_list%tasks
      atom_pair_send => task_list%atom_pair_send
      atom_pair_recv => task_list%atom_pair_recv
      ntasks = task_list%ntasks

      ! *** set up the rs multi-grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho)
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         distributed_rs_grids = rs_rho(igrid_level)%rs_grid%desc%distributed
      END DO

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      !   *** Initialize working density matrix ***
      ! distributed rs grids require a matrix that will be changed
      ! whereas this is not the case for replicated grids
      ALLOCATE (deltap(nimages))
      IF (distributed_rs_grids) THEN
         DO img = 1, nimages
         END DO
         ! this matrix has no strict sparsity pattern in parallel
         ! deltap%sparsity_id=-1
         IF (do_kp) THEN
            DO img = 1, nimages
               CALL dbcsr_copy(deltap(img)%matrix, matrix_p_kp(img)%matrix, &
                               name="DeltaP")
            END DO
         ELSE
            CALL dbcsr_copy(deltap(1)%matrix, matrix_p, name="DeltaP")
         END IF
      ELSE
         IF (do_kp) THEN
            DO img = 1, nimages
               deltap(img)%matrix => matrix_p_kp(img)%matrix
            END DO
         ELSE
            deltap(1)%matrix => matrix_p
         END IF
      ENDIF

      ! distribute the matrix
      IF (distributed_rs_grids) THEN
         CALL rs_distribute_matrix(rs_descs=rs_descs, pmats=deltap, &
                                   atom_pair_send=atom_pair_send, atom_pair_recv=atom_pair_recv, &
                                   nimages=nimages, scatter=.TRUE.)
      ENDIF

      ! map all tasks on the grids

      ithread = 0
      pab => pabt(:, :, ithread)
      work => workt(:, :, ithread)

      loop_xyz: DO idir = 1, 3

         DO igrid_level = 1, gridlevel_info%ngrid_levels
            CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
         END DO

         iatom_old = -1; jatom_old = -1; iset_old = -1; jset_old = -1
         ikind_old = -1; jkind_old = -1; img_old = -1
         loop_tasks: DO itask = 1, ntasks

            !decode the atom pair and basis info
            igrid_level = tasks(itask)%grid_level
            img = tasks(itask)%image
            iatom = tasks(itask)%iatom
            jatom = tasks(itask)%jatom
            iset = tasks(itask)%iset
            jset = tasks(itask)%jset
            ipgf = tasks(itask)%ipgf
            jpgf = tasks(itask)%jpgf

            ikind = particle_set(iatom)%atomic_kind%kind_number
            jkind = particle_set(jatom)%atomic_kind%kind_number

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old .OR. img .NE. img_old) THEN

               IF (iatom .NE. iatom_old) ra(:) = pbc(particle_set(iatom)%r, cell)

               IF (iatom <= jatom) THEN
                  brow = iatom
                  bcol = jatom
               ELSE
                  brow = jatom
                  bcol = iatom
               END IF

               IF (ikind .NE. ikind_old) THEN
                  CALL get_qs_kind(qs_kind_set(ikind), softb=my_soft, &
                                   basis_set=orb_basis_set, basis_type=my_basis_type)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfa, &
                                         lmax=la_max, &
                                         lmin=la_min, &
                                         npgf=npgfa, &
                                         nset=nseta, &
                                         nsgf_set=nsgfa, &
                                         sphi=sphi_a, &
                                         zet=zeta)
               ENDIF

               IF (jkind .NE. jkind_old) THEN
                  CALL get_qs_kind(qs_kind_set(jkind), softb=my_soft, &
                                   basis_set=orb_basis_set, basis_type=my_basis_type)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfb, &
                                         lmax=lb_max, &
                                         lmin=lb_min, &
                                         npgf=npgfb, &
                                         nset=nsetb, &
                                         nsgf_set=nsgfb, &
                                         sphi=sphi_b, &
                                         zet=zetb)
               ENDIF

               CALL dbcsr_get_block_p(matrix=deltap(img)%matrix, &
                                      row=brow, col=bcol, BLOCK=p_block, found=found)
               CPASSERT(found)

               iatom_old = iatom
               jatom_old = jatom
               ikind_old = ikind
               jkind_old = jkind
               img_old = img
               atom_pair_changed = .TRUE.

            ELSE

               atom_pair_changed = .FALSE.

            ENDIF

            IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               IF (iatom <= jatom) THEN
                  CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             p_block(sgfa, sgfb), SIZE(p_block, 1), &
                             0.0_dp, work(1, 1), maxco)
                  CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                             1.0_dp, work(1, 1), maxco, &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, pab(1, 1), maxco)
               ELSE
                  CALL dgemm("N", "N", ncob, nsgfa(iset), nsgfb(jset), &
                             1.0_dp, sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             p_block(sgfb, sgfa), SIZE(p_block, 1), &
                             0.0_dp, work(1, 1), maxco)
                  CALL dgemm("N", "T", ncob, ncoa, nsgfa(iset), &
                             1.0_dp, work(1, 1), maxco, &
                             sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             0.0_dp, pab(1, 1), maxco)
               END IF

               iset_old = iset
               jset_old = jset

            ENDIF

            rab(:) = tasks(itask)%rab
            rb(:) = ra(:) + rab(:)
            zetp = zeta(ipgf, iset) + zetb(jpgf, jset)

            f = zetb(jpgf, jset)/zetp
            rp(:) = ra(:) + f*rab(:)
            prefactor = EXP(-zeta(ipgf, iset)*f*DOT_PRODUCT(rab, rab))
            radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                              lb_min=lb_min(jset), lb_max=lb_max(jset), &
                                              ra=ra, rb=rb, rp=rp, &
                                              zetp=zeta(ipgf, iset), eps=eps_rho_rspace, &
                                              prefactor=prefactor, cutoff=1.0_dp)

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))
            nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
            nb2 = jpgf*ncoset(lb_max(jset))

            ! takes the density matrix symmetry in account, i.e. off-diagonal blocks need to be mapped 'twice'
            IF (iatom == jatom .AND. img == 1) THEN
               scale = 1.0_dp
            ELSE
               scale = 2.0_dp
            END IF

            ! check whether we need to use fawzi's generalised collocation scheme
            IF (rs_rho(igrid_level)%rs_grid%desc%distributed) THEN
               !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
               IF (tasks(itask)%dist_type .EQ. 2) THEN
                  use_subpatch = .TRUE.
               ELSE
                  use_subpatch = .FALSE.
               ENDIF
            ELSE
               use_subpatch = .FALSE.
            ENDIF

            SELECT CASE (idir)
            CASE (1)
               dabqadb_func = GRID_FUNC_DABpADB_X
            CASE (2)
               dabqadb_func = GRID_FUNC_DABpADB_Y
            CASE (3)
               dabqadb_func = GRID_FUNC_DABpADB_Z
            CASE DEFAULT
               CPABORT("invalid idir")
            END SELECT

            IF (iatom <= jatom) THEN
               CALL collocate_pgf_product( &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  ra, rab, scale, pab, na1 - 1, nb1 - 1, &
                  rs_rho(igrid_level)%rs_grid, cell, cube_info(igrid_level), &
                  radius=radius, ga_gb_function=dabqadb_func, &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
            ELSE
               rab_inv = -rab
               CALL collocate_pgf_product( &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  rb, rab_inv, scale, pab, nb1 - 1, na1 - 1, &
                  rs_rho(igrid_level)%rs_grid, cell, cube_info(igrid_level), &
                  radius=radius, ga_gb_function=dabqadb_func, &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
            END IF

         END DO loop_tasks

         CALL density_rs2pw_basic(pw_env, rs_rho, drho(idir), drho_gspace(idir))

      END DO loop_xyz

      !   *** Release work storage ***
      IF (ASSOCIATED(rs_rho)) THEN
         DO i = 1, SIZE(rs_rho)
            CALL rs_grid_release(rs_rho(i)%rs_grid)
         END DO
      END IF

      IF (distributed_rs_grids) THEN
         CALL dbcsr_deallocate_matrix_set(deltap)
      ELSE
         DO img = 1, nimages
            NULLIFY (deltap(img)%matrix)
         END DO
         DEALLOCATE (deltap)
      ENDIF

      DEALLOCATE (pabt, workt)

      CALL timestop(handle)

   END SUBROUTINE calculate_drho_elec

! **************************************************************************************************
!> \brief maps a given wavefunction on the grid
!> \param mo_vectors ...
!> \param ivector ...
!> \param rho ...
!> \param rho_gspace ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param dft_control ...
!> \param particle_set ...
!> \param pw_env ...
!> \param basis_type ...
!> \param external_vector ...
!> \par History
!>      08.2002 created [Joost VandeVondele]
!>      03.2006 made independent of qs_env [Joost VandeVondele]
!> \note
!>      modified calculate_rho_elec, should write the wavefunction represented by
!>      it's presumably dominated by the FFT and the rs->pw and back routines
!>
!>    BUGS ???
!>      does it take correctly the periodic images of the basis functions into account
!>      i.e. is only correct if the basis functions are localised enough to be just in 1 cell ?
! **************************************************************************************************
   SUBROUTINE calculate_wavefunction(mo_vectors, ivector, rho, rho_gspace, &
                                     atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                     pw_env, basis_type, external_vector)

      TYPE(cp_fm_type), POINTER                          :: mo_vectors
      INTEGER                                            :: ivector
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho, rho_gspace
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: external_vector

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_wavefunction'

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER :: dir, group, group_size, handle, i, iatom, igrid_level, ikind, ipgf, iset, maxco, &
         maxsgf_set, my_pos, na1, na2, nao, natom, ncoa, nseta, offset, sgfa
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: where_is_the_point
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: local, map_it_here
      REAL(KIND=dp)                                      :: dab, eps_rho_rspace, radius, scale
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvector
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab, sphi_a, work, zeta
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: mgrid_gspace, mgrid_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(REALSPACE_GRID_P_TYPE), DIMENSION(:), POINTER :: rs_rho

      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      CALL timeset(routineN, handle)

      NULLIFY (eigenvector, orb_basis_set, pab, work, la_max, la_min, &
               npgfa, nsgfa, sphi_a, zeta, first_sgfa, rs_rho, pw_pools, &
               mgrid_rspace, mgrid_gspace)

      IF (PRESENT(external_vector)) THEN
         nao = SIZE(external_vector)
         ALLOCATE (eigenvector(nao))
         eigenvector = external_vector
      ELSE
         CALL cp_fm_get_info(matrix=mo_vectors, nrow_global=nao)
         ALLOCATE (eigenvector(nao))
         DO i = 1, nao
            CALL cp_fm_get_element(mo_vectors, i, ivector, eigenvector(i), local)
         ENDDO
      ENDIF

      ! *** set up the pw multi-grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_grids=rs_rho, pw_pools=pw_pools, &
                      cube_info=cube_info, gridlevel_info=gridlevel_info)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                               use_data=COMPLEXDATA1D, &
                               in_space=RECIPROCALSPACE)
      CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                               use_data=REALDATA3D, &
                               in_space=REALSPACE)

      ! *** set up rs multi-grids
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
!   *** Allocate work storage ***
      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)

      ALLOCATE (pab(maxco, 1))
      ALLOCATE (work(maxco, 1))

      offset = 0
      group = mgrid_rspace(1)%pw%pw_grid%para%group
      my_pos = mgrid_rspace(1)%pw%pw_grid%para%my_pos
      group_size = mgrid_rspace(1)%pw%pw_grid%para%group_size
      ALLOCATE (where_is_the_point(0:group_size - 1))

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=my_basis_type)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                first_sgf=first_sgfa, &
                                lmax=la_max, &
                                lmin=la_min, &
                                npgf=npgfa, &
                                nset=nseta, &
                                nsgf_set=nsgfa, &
                                sphi=sphi_a, &
                                zet=zeta)
         ra(:) = pbc(particle_set(iatom)%r, cell)
         dab = 0.0_dp

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO i = 1, nsgfa(iset)
               work(i, 1) = eigenvector(offset + i)
            ENDDO

            CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), &
                       1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                       work(1, 1), SIZE(work, 1), &
                       0.0_dp, pab(1, 1), SIZE(pab, 1))

            DO ipgf = 1, npgfa(iset)

               na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
               na2 = ipgf*ncoset(la_max(iset))

               scale = 1.0_dp
               igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))

               map_it_here = .FALSE.

               IF (.NOT. ALL(rs_rho(igrid_level)%rs_grid%desc%perd == 1)) THEN
                  DO dir = 1, 3
                     ! bounds of local grid (i.e. removing the 'wings'), if periodic
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_rho(igrid_level)%rs_grid%desc%npts(dir))
                     tp(dir) = MODULO(tp(dir), rs_rho(igrid_level)%rs_grid%desc%npts(dir))
                     IF (rs_rho(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local(dir) + rs_rho(igrid_level)%rs_grid%desc%border
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local(dir) - rs_rho(igrid_level)%rs_grid%desc%border
                     ELSE
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local(dir)
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local(dir)
                     ENDIF
                     ! distributed grid, only map if it is local to the grid
                     location(dir) = tp(dir) + rs_rho(igrid_level)%rs_grid%desc%lb(dir)
                  ENDDO
                  IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                      lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                      lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                     map_it_here = .TRUE.
                  ENDIF
               ELSE
                  ! not distributed, just a round-robin distribution over the full set of CPUs
                  IF (MODULO(offset, group_size) == my_pos) map_it_here = .TRUE.
               ENDIF

               IF (map_it_here) THEN
                  radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                    lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                    zetp=zeta(ipgf, iset), eps=eps_rho_rspace, &
                                                    prefactor=1.0_dp, cutoff=1.0_dp)

                  CALL collocate_pgf_product(la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                             0, 0.0_dp, 0, &
                                             ra, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                             scale, pab, na1 - 1, 0, rs_rho(igrid_level)%rs_grid, &
                                             cell, cube_info(igrid_level), &
                                             radius=radius, ga_gb_function=GRID_FUNC_AB)
               ENDIF

            END DO

            offset = offset + nsgfa(iset)

         END DO

      END DO

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid, &
                             mgrid_rspace(igrid_level)%pw, rs2pw)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid)
      ENDDO

      CALL pw_zero(rho_gspace%pw)
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                          mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw, rho_gspace%pw)
      END DO

      CALL pw_transfer(rho_gspace%pw, rho%pw)

      ! Release work storage
      DEALLOCATE (eigenvector)

      DEALLOCATE (pab)

      DEALLOCATE (work)

      ! give back the pw multi-grids
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

      CALL timestop(handle)

   END SUBROUTINE calculate_wavefunction

END MODULE qs_collocate_density
