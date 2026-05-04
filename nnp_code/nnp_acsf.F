!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functionality for atom centered symmetry functions
!>         for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_acsf
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE nnp_cell_list,                   ONLY: nnp_compute_neighbors_cell_list,&
                                              nnp_prepare_cell_list_cache
   USE nnp_environment_types,           ONLY: nnp_acsf_ang_type,&
                                              nnp_acsf_rad_type,&
                                              nnp_cut_cos,&
                                              nnp_cut_tanh,&
                                              nnp_neighbor_type,&
                                              nnp_symfgrp_type,&
                                              nnp_type
   USE nnp_neighbor_interface,          ONLY: nnp_grp_grow_dGdr,&
                                              nnp_neigh_grp_grow,&
                                              nnp_neighbor_interface_prepare,&
                                              nnp_neighbor_interface_reset_neighbor,&
                                              nnp_workspace_grow_caches
   USE periodic_table,                  ONLY: get_ptable_info
!$ USE omp_lib,                          ONLY: omp_get_max_threads
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_acsf'

   ! Persistent radial sym-function spline tables. Built once on the first
   ! nnp_calc_acsf call from the per-(group, sf) eta/rs/cutoff parameters
   ! and reused thereafter. Each radial group s carries its own dense
   ! 2D table workspace(s)%spline_y(sf, i) and spline_dy(sf, i):
   !   y(r)  = EXP(-eta_sf*(r-rs_sf)^2) * fcut(r)
   !   y'(r) = d/dr [ y(r) ]
   ! Within a group all SFs share the cutoff (=> same uniform grid), so
   ! the inner SF loop in nnp_calc_rad computes the Hermite basis params
   ! ONCE and then issues 4 contiguous reads per SF -- the n_symf loop
   ! is fully auto-vectorizable (sf-first layout = unit stride).
   !
   ! Hermite cubic interpolation is O(dx^4). At n=1024 over a ~12 Bohr cutoff
   ! the residual error accumulates through the NN to ~3e-7 relative, which
   ! tripped the regtest tolerance on H2O-64_NNP_MD. Bumping to 8192 reduces
   ! the per-eval error by 4096x (~1e-11) and brings the test back in tolerance
   ! at a memory cost of ~128 KB per (group) table.
   INTEGER, PARAMETER, PRIVATE                           :: nnp_rad_spline_n = 8192

   ! Tolerance on cutoff equality when grouping symmetry functions in
   ! nnp_init_acsf_groups. SFs sharing a cutoff to within this absolute
   ! tolerance are placed in the same group (same spline grid). Tight
   ! enough to keep distinct cutoffs separate at typical n2p2 precision
   ! (~1e-3 Bohr), loose enough to absorb FP-noise on equal cutoffs.
   REAL(KIND=dp), PARAMETER, PRIVATE                     :: cutoff_eq_tol = 1.0e-5_dp

   ! Public subroutines ***
   PUBLIC :: nnp_calc_acsf, &
             nnp_prepare_neighbor_cache, &
             nnp_init_acsf_groups, &
             nnp_sort_acsf, &
             nnp_sort_ele, &
             nnp_write_acsf

CONTAINS

! **************************************************************************************************
!> \brief Calculate atom centered symmetry functions for given atom i
!>
!> n2p2-style storage. The per-atom symmetry-function gradients live in the
!> per-element neighbor workspace as three sparse arrays:
!>   workspace(ind)%self_dGdr        (3, n_input_nodes)             -- atom i itself
!>   workspace(ind)%dGdr_rad   (s)%data(3, n_symf_s, cap_s)         -- radial group s
!>   workspace(ind)%dGdr_ang_jj(s)%data(3, n_symf_s, cap_s)         -- ang group s, j-side
!>   workspace(ind)%dGdr_ang_kk(s)%data(3, n_symf_s, cap_s)         -- ang group s, k-side
!> Each per-neighbor slab is a small contiguous block addressed by the
!> neighbor's slot in workspace(ind)%neighbor%rad(s)%ind / ang1(s)%ind / ang2(s)%ind.
!> No global (3, n_sf, num_atoms) slab and no touched-atom marker is needed:
!> the only atoms that get derivatives are the ones that already live in the
!> per-(ind) neighbor lists for atom i.
!>
!> \param nnp          NNP environment with persistent neighbour caches populated by
!>                     nnp_prepare_neighbor_cache (must be called once before the per-atom loop).
!> \param i            central-atom index in the global atom ordering; selects
!>                     nnp%ele_ind(i) for per-element scratch routing.
!> \param calc_forces  if .TRUE., populate the per-element dGdr workspaces for caller-side
!>                     force assembly via nnp_scatter_dgdr_to_forces.
!> \param stress       optional per-input-node stress accumulator (only valid when calc_forces).
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_calc_acsf(nnp, i, calc_forces, stress)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: i
      LOGICAL, INTENT(IN)                                :: calc_forces
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: stress

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_calc_acsf'

      INTEGER                                            :: handle, handle_nlist, handle_sf, ind, j, &
                                                            k, l, m, off, peak, s, sf, ii, &
                                                            n_ang1_s, n_ang2_s, n_symf_s, &
                                                            n_input_nodes, nthreads_ang, izeta_il
      REAL(KIND=dp)                                      :: r1, r2, r3, r3_sqr, cutoff_s, &
                                                            cutoff_sqr
      REAL(KIND=dp), DIMENSION(3)                        :: rvect1, rvect2, rvect3
      LOGICAL                                            :: do_forces, homo_grp

      ! Inlined angular kernel variables (eliminates nnp_calc_ang call overhead
      ! and intermediate angular_force3tmp array on the serial force path).
      REAL(KIND=dp) :: rsqr1, rsqr2, rsqr3, r2sum_il, f_il, g_il, costheta_il, &
         fcut3_il, dfcut3_il, ftot_il, dfcutdr1_il, dfcutdr2_il, dfcutdr3_il, &
         r1_inv, r2_inv, r3_inv, inv_g2_il, arg_il, tanh_il, &
         tmp_il, angular_il, tmpzeta_il, symtmp_il, sym_il, &
         eta_il, lam_il, zeta_il, prefzeta_il, pref_lam_il, pref_il, &
         dsymdr1_il, dsymdr2_il, dsymdr3_il, tmp1_il, tmp2_il, tmp3_il, &
         dgdx_t1, dgdx_t2
      REAL(KIND=dp), DIMENSION(3) :: dr1dx_il, dr2dx_il, dr3dx_il, &
         dcosbase1_il, dcosbase2_il, dcosbase3_il, f_jj_il, f_kk_il

      CALL timeset(routineN, handle)

      !determine index of atom type
      ind = nnp%ele_ind(i)
      do_forces = calc_forces

      ! Lazy one-shot build of the per-radial-group spline tables. The tables
      ! depend only on eta/rs/cutoff/cut_type which are fixed once
      ! nnp_init_acsf_groups and nnp_sort_acsf have run, so first-touch here
      ! is safe. The per-group %spline_built flag is per-nnp-instance, so
      ! re-init across runs (e.g. regtests in one process) Just Works.
      ! Guard the indexed access in case this element has no radial groups.
      IF (nnp%rad(ind)%n_symfgrp > 0) THEN
         IF (.NOT. nnp%rad(ind)%symfgrp(1)%spline_built) CALL nnp_build_radial_splines(nnp)
      END IF

      ! Persistent, reusable per-element workspace bindings. The per-group
      ! dGdr_* slabs are grown lazily by nnp_grp_grow_dGdr below; the
      ! fc_cache*/dfc_cache* 1D buffers are bound in an inner ASSOCIATE
      ! AFTER nnp_workspace_grow_caches has sized them to the per-element
      ! peak (binding to a possibly-reallocated allocatable here would
      ! invalidate the alias).
      ASSOCIATE (workspace => nnp%neighbor_interface_state%workspace(ind), &
                 neighbor => nnp%neighbor_interface_state%workspace(ind)%neighbor, &
                 radial_symtmp => nnp%neighbor_interface_state%workspace(ind)%radial_sym, &
                 radial_forcetmp => nnp%neighbor_interface_state%workspace(ind)%radial_force, &
                 angular_symtmp => nnp%neighbor_interface_state%workspace(ind)%angular_sym, &
                 angular_force3tmp => nnp%neighbor_interface_state%workspace(ind)%angular_force, &
                 self_dGdr => nnp%neighbor_interface_state%workspace(ind)%self_dGdr)

         n_input_nodes = nnp%neighbor_interface_state%workspace(ind)%n_input_nodes
         IF (do_forces) self_dGdr(:, 1:n_input_nodes) = 0.0_dp

         ! Traverse linked-cell candidates directly, mirroring the one-pass LAMMPS style.
         ! NOTE: the cell-list/neighbor cache must have been prepared once at per-force-eval
         ! scope by nnp_calc_energy_force via nnp_prepare_neighbor_cache -- DO NOT call it here.
         ! Doing so reintroduces O(N) per-atom prep work (=> O(N^2) per step).
         CALL timeset('nnp_acsf_neighbor_fill', handle_nlist)
         neighbor%pbc_copies = nnp%cell_list_cache%exact_pbc_copies
         CALL nnp_neighbor_interface_reset_neighbor(nnp, ind)
         CALL nnp_compute_neighbors_cell_list(nnp, neighbor, i)
         CALL timestop(handle_nlist)

         ! Reset y:
         nnp%rad(ind)%y = 0.0_dp
         nnp%ang(ind)%y = 0.0_dp

         ! Lazily grow the per-element 1D angular cutoff caches to the
         ! peak per-group neighbor count observed for THIS central atom.
         ! Caches are tracked by an independent high-water mark inside
         ! workspace, so the slab settles at the true element-wide peak
         ! after a few MD steps. Bind the (now stable) allocatables in an
         ! inner ASSOCIATE so the angular kernels access them as plain
         ! contiguous arrays.
         ! MAXVAL of a zero-size array is implementation-defined (gfortran returns
         ! -HUGE, ifort 0), so guard each axis independently.
         peak = 0
         IF (SIZE(neighbor%n_ang1) > 0) peak = MAX(peak, MAXVAL(neighbor%n_ang1))
         IF (SIZE(neighbor%n_ang2) > 0) peak = MAX(peak, MAXVAL(neighbor%n_ang2))
         IF (peak > 0) CALL nnp_workspace_grow_caches(workspace, peak)

         ASSOCIATE (fc_cache1 => workspace%fc_cache1, &
                    dfc_cache1 => workspace%dfc_cache1, &
                    fc_cache2 => workspace%fc_cache2, &
                    dfc_cache2 => workspace%dfc_cache2)

         !calc forces
         IF (do_forces) THEN
            !loop over radial sym fnct grps
            CALL timeset('nnp_acsf_radial', handle_sf)
            DO s = 1, nnp%rad(ind)%n_symfgrp
               n_symf_s = nnp%rad(ind)%symfgrp(s)%n_symf
               ! Per-group dense buffer: (3, n_symf_s, cap_s). Grown lazily.
               CALL nnp_grp_grow_dGdr(workspace%dGdr_rad(s), neighbor%n_rad(s))
               ASSOCIATE (rad_buf => workspace%dGdr_rad(s)%data)
                  !loop over associated neighbors
                  DO j = 1, neighbor%n_rad(s)
                     rvect1 = neighbor%rad(s)%dist(1:3, j)
                     r1 = neighbor%rad(s)%dist(4, j)
                     CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, &
                                       radial_symtmp(1:n_symf_s), &
                                       radial_forcetmp(:, 1:n_symf_s))
                     ! Per-group dense write: rad_buf(:, sf, j) holds dG_m/dr_j.
                     DO sf = 1, n_symf_s
                        m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                        self_dGdr(1, m) = self_dGdr(1, m) + radial_forcetmp(1, sf)
                        self_dGdr(2, m) = self_dGdr(2, m) + radial_forcetmp(2, sf)
                        self_dGdr(3, m) = self_dGdr(3, m) + radial_forcetmp(3, sf)
                        rad_buf(1, sf, j) = -radial_forcetmp(1, sf)
                        rad_buf(2, sf, j) = -radial_forcetmp(2, sf)
                        rad_buf(3, sf, j) = -radial_forcetmp(3, sf)
                        IF (PRESENT(stress)) THEN
                           DO l = 1, 3
                              stress(:, l, m) = stress(:, l, m) + rvect1(:)*radial_forcetmp(l, sf)
                           END DO
                        END IF
                        nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + radial_symtmp(sf)
                     END DO
                  END DO
               END ASSOCIATE
            END DO
            CALL timestop(handle_sf)

            !loop over angular sym fnct grps
            CALL timeset('nnp_acsf_angular', handle_sf)
            off = nnp%n_rad(ind)

            ! OpenMP intra-atom parallelism over the angular group index s.
            ! The serial path below is bit-identical to the pre-OMP kernel.
            ! The parallel path runs only when OpenMP gives us >1 thread AND
            ! there is more than one group to spread across them. Angular
            ! groups partition the input-node index m disjointly (see
            ! nnp_sort_acsf), so self_dGdr, stress, and nnp%ang(ind)%y can
            ! be written concurrently without races or reductions.
            nthreads_ang = 1
!$          nthreads_ang = omp_get_max_threads()
            IF (nthreads_ang > 1 .AND. nnp%ang(ind)%n_symfgrp > 1) THEN
               IF (PRESENT(stress)) THEN
                  CALL nnp_acsf_angular_loop_omp(nnp, ind, self_dGdr, off, stress)
               ELSE
                  CALL nnp_acsf_angular_loop_omp(nnp, ind, self_dGdr, off)
               END IF
            ELSE
            DO s = 1, nnp%ang(ind)%n_symfgrp
               cutoff_s = nnp%ang(ind)%symfgrp(s)%cutoff
               cutoff_sqr = cutoff_s*cutoff_s
               n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
               n_ang1_s = neighbor%n_ang1(s)
               homo_grp = (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2))

               ! Grow per-group dense buffers. jj is always indexed in ang1.
               ! kk is indexed in ang1 for homo groups and in ang2 for hetero.
               CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_jj(s), n_ang1_s)
               IF (homo_grp) THEN
                  CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_kk(s), n_ang1_s)
                  n_ang2_s = 0
               ELSE
                  n_ang2_s = neighbor%n_ang2(s)
                  CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_kk(s), n_ang2_s)
               END IF

               ! Per-group reset. Triplets accumulate into the same (sf, j) slot
               ! across multiple k partners, so we MUST zero before the triplet loop.
               IF (n_ang1_s > 0) workspace%dGdr_ang_jj(s)%data(:, 1:n_symf_s, 1:n_ang1_s) = 0.0_dp
               IF (homo_grp) THEN
                  IF (n_ang1_s > 0) workspace%dGdr_ang_kk(s)%data(:, 1:n_symf_s, 1:n_ang1_s) = 0.0_dp
               ELSE
                  IF (n_ang2_s > 0) workspace%dGdr_ang_kk(s)%data(:, 1:n_symf_s, 1:n_ang2_s) = 0.0_dp
               END IF

               ! Precompute cutoff values and derivatives for ang1 neighbors
               CALL nnp_fill_fc_dfc_cache(neighbor%ang1(s)%dist, n_ang1_s, &
                                          nnp%cut_type, cutoff_s, fc_cache1, dfc_cache1)

               ! ---------------------------------------------------------------
               ! Inlined angular kernel: fused compute + scatter.
               ! Eliminates the nnp_calc_ang subroutine call per triple and
               ! the intermediate angular_force3tmp array. The geometry is
               ! computed once per (j,k) triple, then the SF loop scatters
               ! sym values and forces directly into self_dGdr / jj_buf /
               ! kk_buf without touching any temporary slab. This keeps
               ! dcosbase, dr*dx, and all per-triple scalars in registers
               ! across the SF loop iterations.
               ! Math is bit-identical to nnp_calc_ang.
               ! ---------------------------------------------------------------
               IF (homo_grp) THEN
                  ASSOCIATE (jj_buf => workspace%dGdr_ang_jj(s)%data, &
                             kk_buf => workspace%dGdr_ang_kk(s)%data, &
                             grp_il => nnp%ang(ind)%symfgrp(s))
                  DO j = 1, n_ang1_s
                     rvect1 = neighbor%ang1(s)%dist(1:3, j)
                     r1 = neighbor%ang1(s)%dist(4, j)
                     DO k = j + 1, n_ang1_s
                        rvect2 = neighbor%ang1(s)%dist(1:3, k)
                        r2 = neighbor%ang1(s)%dist(4, k)
                        rvect3(1) = rvect2(1) - rvect1(1)
                        rvect3(2) = rvect2(2) - rvect1(2)
                        rvect3(3) = rvect2(3) - rvect1(3)
                        r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                        IF (r3_sqr < cutoff_sqr) THEN
                           r3 = SQRT(r3_sqr)

                           ! -- per-triple geometry (once) --
                           rsqr1 = r1*r1; rsqr2 = r2*r2; rsqr3 = r3*r3
                           r2sum_il = rsqr1 + rsqr2 + rsqr3
                           f_il = rsqr3 - rsqr1 - rsqr2
                           g_il = -2.0_dp*r1*r2
                           costheta_il = f_il/g_il

                           SELECT CASE (nnp%cut_type)
                           CASE (nnp_cut_cos)
                              arg_il = pi*r3/cutoff_s
                              fcut3_il = 0.5_dp*(COS(arg_il) + 1.0_dp)
                              dfcut3_il = -0.5_dp*SIN(arg_il)*(pi/cutoff_s)
                           CASE (nnp_cut_tanh)
                              tanh_il = TANH(1.0_dp - r3/cutoff_s)
                              fcut3_il = tanh_il**3
                              dfcut3_il = (-3.0_dp/cutoff_s)*(tanh_il**2 - tanh_il**4)
                           CASE DEFAULT
                              CPABORT("NNP| Cutoff function unknown")
                           END SELECT

                           ftot_il = fc_cache1(j)*fc_cache1(k)*fcut3_il
                           dfcutdr1_il = dfc_cache1(j)*fc_cache1(k)*fcut3_il
                           dfcutdr2_il = fc_cache1(j)*dfc_cache1(k)*fcut3_il
                           dfcutdr3_il = fc_cache1(j)*fc_cache1(k)*dfcut3_il

                           r1_inv = 1.0_dp/r1; r2_inv = 1.0_dp/r2; r3_inv = 1.0_dp/r3
                           dr1dx_il(:) = rvect1(:)*r1_inv
                           dr2dx_il(:) = rvect2(:)*r2_inv
                           dr3dx_il(:) = rvect3(:)*r3_inv

                           inv_g2_il = 1.0_dp/(g_il*g_il)
                           DO ii = 1, 3
                              dgdx_t1 = 2.0_dp*r2*dr1dx_il(ii)
                              dgdx_t2 = 2.0_dp*r1*dr2dx_il(ii)
                              dcosbase1_il(ii) = -2.0_dp*(rvect1(ii) + rvect2(ii))*g_il &
                                                 - f_il*(-(dgdx_t1 + dgdx_t2))
                              dcosbase2_il(ii) = 2.0_dp*(rvect3(ii) + rvect1(ii))*g_il &
                                                 - f_il*dgdx_t1
                              dcosbase3_il(ii) = 2.0_dp*(rvect2(ii) - rvect3(ii))*g_il &
                                                 - f_il*dgdx_t2
                           END DO

                           ! -- fused SF loop: compute + direct scatter --
                           DO sf = 1, n_symf_s
                              m = off + grp_il%symf(sf)
                              lam_il = grp_il%pack_lam(sf)
                              zeta_il = grp_il%pack_zeta(sf)
                              eta_il = grp_il%pack_eta(sf)
                              prefzeta_il = grp_il%pack_prefzeta(sf)

                              tmp_il = 1.0_dp + lam_il*costheta_il
                              IF (tmp_il <= 0.0_dp) THEN
                                 tmpzeta_il = 0.0_dp
                                 angular_il = 0.0_dp
                              ELSE
                                 IF (grp_il%pack_use_int_zeta(sf)) THEN
                                    izeta_il = grp_il%pack_izeta(sf)
                                    tmpzeta_il = tmp_il**(izeta_il - 1)
                                 ELSE
                                    tmpzeta_il = tmp_il**(zeta_il - 1.0_dp)
                                 END IF
                                 angular_il = tmpzeta_il*tmp_il
                              END IF

                              symtmp_il = EXP(-eta_il*r2sum_il)
                              sym_il = prefzeta_il*angular_il*symtmp_il*ftot_il
                              nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + sym_il

                              pref_lam_il = zeta_il*tmpzeta_il*lam_il*inv_g2_il
                              tmp_il = -2.0_dp*symtmp_il*eta_il
                              dsymdr1_il = tmp_il*r1
                              dsymdr2_il = tmp_il*r2
                              dsymdr3_il = tmp_il*r3

                              pref_il = prefzeta_il*symtmp_il*ftot_il
                              tmp1_il = prefzeta_il*angular_il*(ftot_il*dsymdr1_il + dfcutdr1_il*symtmp_il)
                              tmp2_il = prefzeta_il*angular_il*(ftot_il*dsymdr2_il + dfcutdr2_il*symtmp_il)
                              tmp3_il = prefzeta_il*angular_il*(ftot_il*dsymdr3_il + dfcutdr3_il*symtmp_il)

                              DO ii = 1, 3
                                 f_jj_il(ii) = pref_il*pref_lam_il*dcosbase2_il(ii) &
                                               - tmp1_il*dr1dx_il(ii) + tmp3_il*dr3dx_il(ii)
                                 f_kk_il(ii) = pref_il*pref_lam_il*dcosbase3_il(ii) &
                                               - tmp2_il*dr2dx_il(ii) - tmp3_il*dr3dx_il(ii)
                                 self_dGdr(ii, m) = self_dGdr(ii, m) &
                                    + pref_il*pref_lam_il*dcosbase1_il(ii) &
                                    + tmp1_il*dr1dx_il(ii) + tmp2_il*dr2dx_il(ii)
                                 jj_buf(ii, sf, j) = jj_buf(ii, sf, j) + f_jj_il(ii)
                                 kk_buf(ii, sf, k) = kk_buf(ii, sf, k) + f_kk_il(ii)
                              END DO
                              IF (PRESENT(stress)) THEN
                                 DO l = 1, 3
                                    stress(:, l, m) = stress(:, l, m) &
                                       - rvect1(:)*f_jj_il(l) - rvect2(:)*f_kk_il(l)
                                 END DO
                              END IF
                           END DO

                        END IF
                     END DO
                  END DO
                  END ASSOCIATE
               ELSE
                  ! Precompute cutoff values for ang2 neighbors (different elements)
                  CALL nnp_fill_fc_dfc_cache(neighbor%ang2(s)%dist, n_ang2_s, &
                                             nnp%cut_type, cutoff_s, fc_cache2, dfc_cache2)

                  ASSOCIATE (jj_buf => workspace%dGdr_ang_jj(s)%data, &
                             kk_buf => workspace%dGdr_ang_kk(s)%data, &
                             grp_il => nnp%ang(ind)%symfgrp(s))
                  DO j = 1, n_ang1_s
                     rvect1 = neighbor%ang1(s)%dist(1:3, j)
                     r1 = neighbor%ang1(s)%dist(4, j)
                     DO k = 1, n_ang2_s
                        rvect2 = neighbor%ang2(s)%dist(1:3, k)
                        r2 = neighbor%ang2(s)%dist(4, k)
                        rvect3(1) = rvect2(1) - rvect1(1)
                        rvect3(2) = rvect2(2) - rvect1(2)
                        rvect3(3) = rvect2(3) - rvect1(3)
                        r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                        IF (r3_sqr < cutoff_sqr) THEN
                           r3 = SQRT(r3_sqr)

                           ! -- per-triple geometry (once) --
                           rsqr1 = r1*r1; rsqr2 = r2*r2; rsqr3 = r3*r3
                           r2sum_il = rsqr1 + rsqr2 + rsqr3
                           f_il = rsqr3 - rsqr1 - rsqr2
                           g_il = -2.0_dp*r1*r2
                           costheta_il = f_il/g_il

                           SELECT CASE (nnp%cut_type)
                           CASE (nnp_cut_cos)
                              arg_il = pi*r3/cutoff_s
                              fcut3_il = 0.5_dp*(COS(arg_il) + 1.0_dp)
                              dfcut3_il = -0.5_dp*SIN(arg_il)*(pi/cutoff_s)
                           CASE (nnp_cut_tanh)
                              tanh_il = TANH(1.0_dp - r3/cutoff_s)
                              fcut3_il = tanh_il**3
                              dfcut3_il = (-3.0_dp/cutoff_s)*(tanh_il**2 - tanh_il**4)
                           CASE DEFAULT
                              CPABORT("NNP| Cutoff function unknown")
                           END SELECT

                           ftot_il = fc_cache1(j)*fc_cache2(k)*fcut3_il
                           dfcutdr1_il = dfc_cache1(j)*fc_cache2(k)*fcut3_il
                           dfcutdr2_il = fc_cache1(j)*dfc_cache2(k)*fcut3_il
                           dfcutdr3_il = fc_cache1(j)*fc_cache2(k)*dfcut3_il

                           r1_inv = 1.0_dp/r1; r2_inv = 1.0_dp/r2; r3_inv = 1.0_dp/r3
                           dr1dx_il(:) = rvect1(:)*r1_inv
                           dr2dx_il(:) = rvect2(:)*r2_inv
                           dr3dx_il(:) = rvect3(:)*r3_inv

                           inv_g2_il = 1.0_dp/(g_il*g_il)
                           DO ii = 1, 3
                              dgdx_t1 = 2.0_dp*r2*dr1dx_il(ii)
                              dgdx_t2 = 2.0_dp*r1*dr2dx_il(ii)
                              dcosbase1_il(ii) = -2.0_dp*(rvect1(ii) + rvect2(ii))*g_il &
                                                 - f_il*(-(dgdx_t1 + dgdx_t2))
                              dcosbase2_il(ii) = 2.0_dp*(rvect3(ii) + rvect1(ii))*g_il &
                                                 - f_il*dgdx_t1
                              dcosbase3_il(ii) = 2.0_dp*(rvect2(ii) - rvect3(ii))*g_il &
                                                 - f_il*dgdx_t2
                           END DO

                           ! -- fused SF loop: compute + direct scatter --
                           DO sf = 1, n_symf_s
                              m = off + grp_il%symf(sf)
                              lam_il = grp_il%pack_lam(sf)
                              zeta_il = grp_il%pack_zeta(sf)
                              eta_il = grp_il%pack_eta(sf)
                              prefzeta_il = grp_il%pack_prefzeta(sf)

                              tmp_il = 1.0_dp + lam_il*costheta_il
                              IF (tmp_il <= 0.0_dp) THEN
                                 tmpzeta_il = 0.0_dp
                                 angular_il = 0.0_dp
                              ELSE
                                 IF (grp_il%pack_use_int_zeta(sf)) THEN
                                    izeta_il = grp_il%pack_izeta(sf)
                                    tmpzeta_il = tmp_il**(izeta_il - 1)
                                 ELSE
                                    tmpzeta_il = tmp_il**(zeta_il - 1.0_dp)
                                 END IF
                                 angular_il = tmpzeta_il*tmp_il
                              END IF

                              symtmp_il = EXP(-eta_il*r2sum_il)
                              sym_il = prefzeta_il*angular_il*symtmp_il*ftot_il
                              nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + sym_il

                              pref_lam_il = zeta_il*tmpzeta_il*lam_il*inv_g2_il
                              tmp_il = -2.0_dp*symtmp_il*eta_il
                              dsymdr1_il = tmp_il*r1
                              dsymdr2_il = tmp_il*r2
                              dsymdr3_il = tmp_il*r3

                              pref_il = prefzeta_il*symtmp_il*ftot_il
                              tmp1_il = prefzeta_il*angular_il*(ftot_il*dsymdr1_il + dfcutdr1_il*symtmp_il)
                              tmp2_il = prefzeta_il*angular_il*(ftot_il*dsymdr2_il + dfcutdr2_il*symtmp_il)
                              tmp3_il = prefzeta_il*angular_il*(ftot_il*dsymdr3_il + dfcutdr3_il*symtmp_il)

                              DO ii = 1, 3
                                 f_jj_il(ii) = pref_il*pref_lam_il*dcosbase2_il(ii) &
                                               - tmp1_il*dr1dx_il(ii) + tmp3_il*dr3dx_il(ii)
                                 f_kk_il(ii) = pref_il*pref_lam_il*dcosbase3_il(ii) &
                                               - tmp2_il*dr2dx_il(ii) - tmp3_il*dr3dx_il(ii)
                                 self_dGdr(ii, m) = self_dGdr(ii, m) &
                                    + pref_il*pref_lam_il*dcosbase1_il(ii) &
                                    + tmp1_il*dr1dx_il(ii) + tmp2_il*dr2dx_il(ii)
                                 jj_buf(ii, sf, j) = jj_buf(ii, sf, j) + f_jj_il(ii)
                                 kk_buf(ii, sf, k) = kk_buf(ii, sf, k) + f_kk_il(ii)
                              END DO
                              IF (PRESENT(stress)) THEN
                                 DO l = 1, 3
                                    stress(:, l, m) = stress(:, l, m) &
                                       - rvect1(:)*f_jj_il(l) - rvect2(:)*f_kk_il(l)
                                 END DO
                              END IF
                           END DO

                        END IF
                     END DO
                  END DO
                  END ASSOCIATE
               END IF
            END DO
            END IF
            CALL timestop(handle_sf)
         ELSE
            !loop over radial sym fnct grps
            CALL timeset('nnp_acsf_radial', handle_sf)
            DO s = 1, nnp%rad(ind)%n_symfgrp
               !loop over associated neighbors
               DO j = 1, neighbor%n_rad(s)
                  rvect1 = neighbor%rad(s)%dist(1:3, j)
                  r1 = neighbor%rad(s)%dist(4, j)
                  CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, radial_symtmp(1:nnp%rad(ind)%symfgrp(s)%n_symf))
                  DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
                     m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                     nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + radial_symtmp(sf)
                  END DO
               END DO
            END DO
            CALL timestop(handle_sf)

            !loop over angular sym fnct grps
            CALL timeset('nnp_acsf_angular', handle_sf)
            off = nnp%n_rad(ind)
            DO s = 1, nnp%ang(ind)%n_symfgrp
               cutoff_s = nnp%ang(ind)%symfgrp(s)%cutoff
               cutoff_sqr = cutoff_s*cutoff_s
               n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
               n_ang1_s = neighbor%n_ang1(s)

               ! Precompute cutoff values for ang1 neighbors (no derivatives needed)
               CALL nnp_fill_fc_cache(neighbor%ang1(s)%dist, n_ang1_s, &
                                      nnp%cut_type, cutoff_s, fc_cache1)

               IF (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2)) THEN
                  DO j = 1, n_ang1_s
                     rvect1 = neighbor%ang1(s)%dist(1:3, j)
                     r1 = neighbor%ang1(s)%dist(4, j)
                     DO k = j + 1, n_ang1_s
                        rvect2 = neighbor%ang1(s)%dist(1:3, k)
                        r2 = neighbor%ang1(s)%dist(4, k)
                        rvect3(1) = rvect2(1) - rvect1(1)
                        rvect3(2) = rvect2(2) - rvect1(2)
                        rvect3(3) = rvect2(3) - rvect1(3)
                        r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                        IF (r3_sqr < cutoff_sqr) THEN
                           r3 = SQRT(r3_sqr)
                           CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, &
                                             fc_cache1(j), 0.0_dp, fc_cache1(k), 0.0_dp, &
                                             angular_symtmp(1:n_symf_s))
                           DO sf = 1, n_symf_s
                              m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                              nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + angular_symtmp(sf)
                           END DO
                        END IF
                     END DO
                  END DO
               ELSE
                  ! Precompute cutoff values for ang2 neighbors
                  n_ang2_s = neighbor%n_ang2(s)
                  CALL nnp_fill_fc_cache(neighbor%ang2(s)%dist, n_ang2_s, &
                                         nnp%cut_type, cutoff_s, fc_cache2)

                  DO j = 1, n_ang1_s
                     rvect1 = neighbor%ang1(s)%dist(1:3, j)
                     r1 = neighbor%ang1(s)%dist(4, j)
                     DO k = 1, n_ang2_s
                        rvect2 = neighbor%ang2(s)%dist(1:3, k)
                        r2 = neighbor%ang2(s)%dist(4, k)
                        rvect3(1) = rvect2(1) - rvect1(1)
                        rvect3(2) = rvect2(2) - rvect1(2)
                        rvect3(3) = rvect2(3) - rvect1(3)
                        r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                        IF (r3_sqr < cutoff_sqr) THEN
                           r3 = SQRT(r3_sqr)
                           CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, &
                                             fc_cache1(j), 0.0_dp, fc_cache2(k), 0.0_dp, &
                                             angular_symtmp(1:n_symf_s))
                           DO sf = 1, n_symf_s
                              m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                              nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + angular_symtmp(sf)
                           END DO
                        END IF
                     END DO
                  END DO
               END IF
            END DO
            CALL timestop(handle_sf)
         END IF

         END ASSOCIATE

         ! fc_cache1/2 and dfc_cache1/2 are persistent workspace; nothing to deallocate here.

      END ASSOCIATE

      !check extrapolation
      CALL nnp_check_extrapolation(nnp, ind)

      IF (PRESENT(stress)) THEN
         CALL nnp_scale_acsf(nnp, ind, do_forces, stress)
      ELSE
         CALL nnp_scale_acsf(nnp, ind, do_forces)
      END IF

      CALL timestop(handle)

   END SUBROUTINE nnp_calc_acsf

! **************************************************************************************************
!> \brief Fill the per-neighbor fc and dfc cutoff function caches used by the
!>        force branch of the angular ACSF kernel. The arithmetic here is held
!>        bit-identical to the original inline caches so that extracting this
!>        helper does not perturb per-step forces (the H2O-64 NNP regtest is
!>        Lyapunov-sensitive at ~1e-12 per step over 10 steps).
!> \param dist       (4, :) neighbor array; column 4 holds the scalar distance
!> \param n          number of neighbors to process
!> \param cut_type   cutoff function selector (nnp_cut_cos / nnp_cut_tanh)
!> \param cutoff_s   per-group cutoff radius
!> \param fc_cache   output fc values, sized >= n
!> \param dfc_cache  output dfc values, sized >= n
! **************************************************************************************************
   PURE SUBROUTINE nnp_fill_fc_dfc_cache(dist, n, cut_type, cutoff_s, fc_cache, dfc_cache)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: dist
      INTEGER, INTENT(IN)                                :: n, cut_type
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_s
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: fc_cache, dfc_cache

      INTEGER                                            :: j
      REAL(KIND=dp)                                      :: r_tmp, arg_tmp, tanh_tmp

      DO j = 1, n
         r_tmp = dist(4, j)
         SELECT CASE (cut_type)
         CASE (nnp_cut_cos)
            arg_tmp = pi*r_tmp/cutoff_s
            fc_cache(j) = 0.5_dp*(COS(arg_tmp) + 1.0_dp)
            dfc_cache(j) = -0.5_dp*SIN(arg_tmp)*(pi/cutoff_s)
         CASE (nnp_cut_tanh)
            tanh_tmp = TANH(1.0_dp - r_tmp/cutoff_s)
            fc_cache(j) = tanh_tmp**3
            dfc_cache(j) = (-3.0_dp/cutoff_s)*(tanh_tmp**2 - tanh_tmp**4)
         END SELECT
      END DO

   END SUBROUTINE nnp_fill_fc_dfc_cache

! **************************************************************************************************
!> \brief Fill the per-neighbor fc cutoff function cache used by the sym-only
!>        (no-forces) branch of the angular ACSF kernel. Derivatives are not
!>        needed so the arg_tmp intermediate is folded into the COS call to
!>        preserve bit-identical arithmetic with the original inline cache.
!> \param dist       (4, :) neighbor array; column 4 holds the scalar distance
!> \param n          number of neighbors to process
!> \param cut_type   cutoff function selector (nnp_cut_cos / nnp_cut_tanh)
!> \param cutoff_s   per-group cutoff radius
!> \param fc_cache   output fc values, sized >= n
! **************************************************************************************************
   PURE SUBROUTINE nnp_fill_fc_cache(dist, n, cut_type, cutoff_s, fc_cache)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: dist
      INTEGER, INTENT(IN)                                :: n, cut_type
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_s
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: fc_cache

      INTEGER                                            :: j
      REAL(KIND=dp)                                      :: r_tmp, tanh_tmp

      DO j = 1, n
         r_tmp = dist(4, j)
         SELECT CASE (cut_type)
         CASE (nnp_cut_cos)
            fc_cache(j) = 0.5_dp*(COS(pi*r_tmp/cutoff_s) + 1.0_dp)
         CASE (nnp_cut_tanh)
            tanh_tmp = TANH(1.0_dp - r_tmp/cutoff_s)
            fc_cache(j) = tanh_tmp**3
         END SELECT
      END DO

   END SUBROUTINE nnp_fill_fc_cache

! **************************************************************************************************
!> \brief OpenMP intra-atom parallelization of the forces-branch angular
!>        symmetry-function group loop. Only dispatched when the process
!>        has more than one OpenMP thread AND there is more than one
!>        angular group for the current element -- otherwise the caller
!>        runs the original serial loop which is kept bit-identical to
!>        the pre-OpenMP kernel.
!>
!>        Thread safety notes:
!>
!>        * The outermost loop index is the angular group s. Groups
!>          partition the angular symmetry-function indices disjointly
!>          (see nnp_sort_acsf), so the reductions into self_dGdr(:, m),
!>          stress(:, :, m) and nnp%ang(ind)%y(m-off) -- whose m depends
!>          only on s via symfgrp(s)%symf(sf) -- never collide between
!>          threads. No OpenMP reductions or critical sections are
!>          required for these global outputs.
!>        * The per-group accumulators workspace%dGdr_ang_jj(s)%data and
!>          workspace%dGdr_ang_kk(s)%data are also indexed by s, so each
!>          thread writes to its own group's slab.
!>        * The cutoff caches, angular_sym and angular_force scratch that
!>          live on workspace would otherwise alias across threads. The
!>          helper therefore stages them as automatic arrays sized from
!>          workspace%cache_cap / workspace%max_ang_symf and lists them in
!>          the OMP PRIVATE clause so every thread gets its own stack copy.
!>        * Per-group buffer growth (nnp_grp_grow_dGdr) and the zero-init
!>          of the accumulator slabs are done in a small serial pre-pass
!>          so the parallel region itself never touches ALLOCATABLE state.
!>        * nnp_calc_ang is invoked with INTENT(IN) nnp; the input-only
!>          usage has been verified safe for concurrent calls from
!>          different threads (it writes only into the caller-provided
!>          sym and force output arrays).
!>
!>        The kernel arithmetic is a direct copy of the serial body so
!>        the only behavioural difference between this path and the
!>        serial ELSE branch is the order in which groups are visited;
!>        because groups write to disjoint state that order does not
!>        change any final numerical value at the bit level for n_omp=1
!>        (which also takes the serial ELSE branch unconditionally via
!>        the nthreads_ang>1 guard at the call site).
! **************************************************************************************************
   SUBROUTINE nnp_acsf_angular_loop_omp(nnp, ind, self_dGdr, off, stress)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: ind, off
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: self_dGdr
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: stress

      INTEGER                                            :: s, j, k, sf, m, l
      INTEGER                                            :: n_symf_s, n_ang1_s, n_ang2_s
      REAL(KIND=dp)                                      :: r1, r2, r3, r3_sqr, &
                                                            cutoff_s, cutoff_sqr
      REAL(KIND=dp), DIMENSION(3)                        :: rvect1, rvect2, rvect3
      LOGICAL                                            :: homo_grp
      INTEGER                                            :: cache_cap_loc, max_ang_symf_loc
      ! Per-thread automatic scratch (see PRIVATE clause below). Sizes
      ! are pulled from the pre-sized workspace so they match the peak
      ! usage that the serial kernel would allocate into the workspace
      ! %fc_cache*/angular_* slabs.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fc_c1_loc, dfc_c1_loc, &
                                                            fc_c2_loc, dfc_c2_loc, sym_loc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: force_loc

      ASSOCIATE (workspace => nnp%neighbor_interface_state%workspace(ind), &
                 neighbor => nnp%neighbor_interface_state%workspace(ind)%neighbor)

      cache_cap_loc = MAX(1, workspace%cache_cap)
      max_ang_symf_loc = MAX(1, workspace%max_ang_symf)

      ! Serial pre-pass: grow per-group accumulator slabs and zero them
      ! out. Doing this outside the parallel region guarantees no thread
      ! ever touches ALLOCATABLE components of nnp_dGdr_grp_type.
      DO s = 1, nnp%ang(ind)%n_symfgrp
         n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
         n_ang1_s = neighbor%n_ang1(s)
         homo_grp = (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2))
         CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_jj(s), n_ang1_s)
         IF (homo_grp) THEN
            CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_kk(s), n_ang1_s)
            n_ang2_s = 0
         ELSE
            n_ang2_s = neighbor%n_ang2(s)
            CALL nnp_grp_grow_dGdr(workspace%dGdr_ang_kk(s), n_ang2_s)
         END IF
         IF (n_ang1_s > 0) workspace%dGdr_ang_jj(s)%data(:, 1:n_symf_s, 1:n_ang1_s) = 0.0_dp
         IF (homo_grp) THEN
            IF (n_ang1_s > 0) workspace%dGdr_ang_kk(s)%data(:, 1:n_symf_s, 1:n_ang1_s) = 0.0_dp
         ELSE
            IF (n_ang2_s > 0) workspace%dGdr_ang_kk(s)%data(:, 1:n_symf_s, 1:n_ang2_s) = 0.0_dp
         END IF
      END DO

      ! Each thread allocates its own private scratch inside the
      ! parallel region via the PRIVATE clause. Entering the region an
      ! ALLOCATABLE listed as PRIVATE starts in the unallocated state on
      ! every thread, so the per-thread ALLOCATE below gives each thread
      ! its own independent scratch slab.
      !
      ! DEFAULT(SHARED) is used deliberately rather than DEFAULT(NONE) so
      ! that the OPTIONAL `stress` dummy does not have to appear in any
      ! data-sharing clause (referencing an absent optional in SHARED is
      ! ill-defined across compilers; this way the PRESENT() guard below
      ! is the sole gate on whether stress is ever touched).
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP          PRIVATE(s, j, k, sf, m, l, n_symf_s, n_ang1_s, n_ang2_s, &
!$OMP                  r1, r2, r3, r3_sqr, cutoff_s, cutoff_sqr, homo_grp, &
!$OMP                  rvect1, rvect2, rvect3, &
!$OMP                  fc_c1_loc, dfc_c1_loc, fc_c2_loc, dfc_c2_loc, &
!$OMP                  sym_loc, force_loc)
      ALLOCATE (fc_c1_loc(cache_cap_loc))
      ALLOCATE (dfc_c1_loc(cache_cap_loc))
      ALLOCATE (fc_c2_loc(cache_cap_loc))
      ALLOCATE (dfc_c2_loc(cache_cap_loc))
      ALLOCATE (sym_loc(max_ang_symf_loc))
      ALLOCATE (force_loc(3, 3, max_ang_symf_loc))

!$OMP DO SCHEDULE(dynamic)
      DO s = 1, nnp%ang(ind)%n_symfgrp
         cutoff_s = nnp%ang(ind)%symfgrp(s)%cutoff
         cutoff_sqr = cutoff_s*cutoff_s
         n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
         n_ang1_s = neighbor%n_ang1(s)
         homo_grp = (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2))

         CALL nnp_fill_fc_dfc_cache(neighbor%ang1(s)%dist, n_ang1_s, &
                                    nnp%cut_type, cutoff_s, fc_c1_loc, dfc_c1_loc)

         IF (homo_grp) THEN
            ASSOCIATE (jj_buf => workspace%dGdr_ang_jj(s)%data, &
                       kk_buf => workspace%dGdr_ang_kk(s)%data)
            DO j = 1, n_ang1_s
               rvect1 = neighbor%ang1(s)%dist(1:3, j)
               r1 = neighbor%ang1(s)%dist(4, j)
               DO k = j + 1, n_ang1_s
                  rvect2 = neighbor%ang1(s)%dist(1:3, k)
                  r2 = neighbor%ang1(s)%dist(4, k)
                  rvect3(1) = rvect2(1) - rvect1(1)
                  rvect3(2) = rvect2(2) - rvect1(2)
                  rvect3(3) = rvect2(3) - rvect1(3)
                  r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                  IF (r3_sqr < cutoff_sqr) THEN
                     r3 = SQRT(r3_sqr)
                     CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                       r1, r2, r3, &
                                       fc_c1_loc(j), dfc_c1_loc(j), &
                                       fc_c1_loc(k), dfc_c1_loc(k), &
                                       sym_loc(1:n_symf_s), &
                                       force_loc(:, :, 1:n_symf_s))
                     DO sf = 1, n_symf_s
                        m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                        self_dGdr(1, m) = self_dGdr(1, m) + force_loc(1, 1, sf)
                        self_dGdr(2, m) = self_dGdr(2, m) + force_loc(2, 1, sf)
                        self_dGdr(3, m) = self_dGdr(3, m) + force_loc(3, 1, sf)
                        jj_buf(1, sf, j) = jj_buf(1, sf, j) + force_loc(1, 2, sf)
                        jj_buf(2, sf, j) = jj_buf(2, sf, j) + force_loc(2, 2, sf)
                        jj_buf(3, sf, j) = jj_buf(3, sf, j) + force_loc(3, 2, sf)
                        kk_buf(1, sf, k) = kk_buf(1, sf, k) + force_loc(1, 3, sf)
                        kk_buf(2, sf, k) = kk_buf(2, sf, k) + force_loc(2, 3, sf)
                        kk_buf(3, sf, k) = kk_buf(3, sf, k) + force_loc(3, 3, sf)
                        IF (PRESENT(stress)) THEN
                           DO l = 1, 3
                              stress(:, l, m) = stress(:, l, m) - rvect1(:)*force_loc(l, 2, sf)
                              stress(:, l, m) = stress(:, l, m) - rvect2(:)*force_loc(l, 3, sf)
                           END DO
                        END IF
                        nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + sym_loc(sf)
                     END DO
                  END IF
               END DO
            END DO
            END ASSOCIATE
         ELSE
            n_ang2_s = neighbor%n_ang2(s)
            CALL nnp_fill_fc_dfc_cache(neighbor%ang2(s)%dist, n_ang2_s, &
                                       nnp%cut_type, cutoff_s, fc_c2_loc, dfc_c2_loc)

            ASSOCIATE (jj_buf => workspace%dGdr_ang_jj(s)%data, &
                       kk_buf => workspace%dGdr_ang_kk(s)%data)
            DO j = 1, n_ang1_s
               rvect1 = neighbor%ang1(s)%dist(1:3, j)
               r1 = neighbor%ang1(s)%dist(4, j)
               DO k = 1, n_ang2_s
                  rvect2 = neighbor%ang2(s)%dist(1:3, k)
                  r2 = neighbor%ang2(s)%dist(4, k)
                  rvect3(1) = rvect2(1) - rvect1(1)
                  rvect3(2) = rvect2(2) - rvect1(2)
                  rvect3(3) = rvect2(3) - rvect1(3)
                  r3_sqr = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                  IF (r3_sqr < cutoff_sqr) THEN
                     r3 = SQRT(r3_sqr)
                     CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                       r1, r2, r3, &
                                       fc_c1_loc(j), dfc_c1_loc(j), &
                                       fc_c2_loc(k), dfc_c2_loc(k), &
                                       sym_loc(1:n_symf_s), &
                                       force_loc(:, :, 1:n_symf_s))
                     DO sf = 1, n_symf_s
                        m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                        self_dGdr(1, m) = self_dGdr(1, m) + force_loc(1, 1, sf)
                        self_dGdr(2, m) = self_dGdr(2, m) + force_loc(2, 1, sf)
                        self_dGdr(3, m) = self_dGdr(3, m) + force_loc(3, 1, sf)
                        jj_buf(1, sf, j) = jj_buf(1, sf, j) + force_loc(1, 2, sf)
                        jj_buf(2, sf, j) = jj_buf(2, sf, j) + force_loc(2, 2, sf)
                        jj_buf(3, sf, j) = jj_buf(3, sf, j) + force_loc(3, 2, sf)
                        kk_buf(1, sf, k) = kk_buf(1, sf, k) + force_loc(1, 3, sf)
                        kk_buf(2, sf, k) = kk_buf(2, sf, k) + force_loc(2, 3, sf)
                        kk_buf(3, sf, k) = kk_buf(3, sf, k) + force_loc(3, 3, sf)
                        IF (PRESENT(stress)) THEN
                           DO l = 1, 3
                              stress(:, l, m) = stress(:, l, m) - rvect1(:)*force_loc(l, 2, sf)
                              stress(:, l, m) = stress(:, l, m) - rvect2(:)*force_loc(l, 3, sf)
                           END DO
                        END IF
                        nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + sym_loc(sf)
                     END DO
                  END IF
               END DO
            END DO
            END ASSOCIATE
         END IF
      END DO
!$OMP END DO

      DEALLOCATE (fc_c1_loc, dfc_c1_loc, fc_c2_loc, dfc_c2_loc, sym_loc, force_loc)
!$OMP END PARALLEL

      END ASSOCIATE

   END SUBROUTINE nnp_acsf_angular_loop_omp

! **************************************************************************************************
!> \brief Prepare or update the linked-cell / Verlet cache for the current geometry.
!>        Lazily allocates nnp%cell_list_cache and nnp%neighbor_interface_state on
!>        first use, then delegates to nnp_prepare_cell_list_cache and
!>        nnp_neighbor_interface_prepare. Must be called once per force eval before
!>        the per-atom loop in nnp_calc_energy_force; idempotent and cheap on
!>        re-entry when nothing has changed (Verlet skin not tripped, model
!>        topology stable).
!> \param nnp NNP environment whose persistent caches are to be (re)built.
!> \par History
!>      2026-04 created [dsharma]
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_prepare_neighbor_cache(nnp)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp

      IF (.NOT. ALLOCATED(nnp%cell_list_cache)) ALLOCATE (nnp%cell_list_cache)
      IF (.NOT. ALLOCATED(nnp%neighbor_interface_state)) ALLOCATE (nnp%neighbor_interface_state)
      CALL nnp_prepare_cell_list_cache(nnp)
      CALL nnp_neighbor_interface_prepare(nnp)

   END SUBROUTINE nnp_prepare_neighbor_cache


! **************************************************************************************************
!> \brief Check if the nnp is extrapolating
! **************************************************************************************************
   SUBROUTINE nnp_check_extrapolation(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      REAL(KIND=dp), PARAMETER                           :: threshold = 0.0001_dp

      INTEGER                                            :: j
      LOGICAL                                            :: extrapolate

      extrapolate = nnp%output_expol

      DO j = 1, nnp%n_rad(ind)
         IF (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%rad(ind)%y(j) + nnp%rad(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO
      DO j = 1, nnp%n_ang(ind)
         IF (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%ang(ind)%y(j) + nnp%ang(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO

      nnp%output_expol = extrapolate

   END SUBROUTINE nnp_check_extrapolation

! **************************************************************************************************
!> \brief Scale and center symmetry functions (and gradients)
! **************************************************************************************************
   SUBROUTINE nnp_scale_acsf(nnp, ind, do_forces, stress)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind
      LOGICAL, INTENT(IN)                                :: do_forces
      ! INOUT (not OUT): stress is per-input-node and accumulates across central atoms.
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: stress

      INTEGER                                            :: j, k, m, off, s, sf, &
                                                            n_symf_s, n_ang1_s, n_ang2_s
      REAL(KIND=dp)                                      :: scale
      LOGICAL                                            :: homo_grp

      IF (nnp%center_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j)
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j)
         END DO

         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               nnp%arc(ind)%layer(1)%node(j) = nnp%arc(ind)%layer(1)%node(j)/ &
                  (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))*(nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               nnp%arc(ind)%layer(1)%node(j + off) = nnp%arc(ind)%layer(1)%node(j + off)/ &
                  (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))*(nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
         END IF
      ELSE IF (nnp%scale_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_min(j))/ &
                                            (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
                                            (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_min(j))/ &
                                                  (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
                                                  (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE IF (nnp%scale_sigma_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j))/ &
                                            nnp%rad(ind)%sigma(j)*(nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j))/ &
                                                  nnp%ang(ind)%sigma(j)*(nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = nnp%rad(ind)%y(j)
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = nnp%ang(ind)%y(j)
         END DO
      END IF

      IF (do_forces .AND. (nnp%scale_acsf .OR. nnp%scale_sigma_acsf)) THEN
         ! Scale the per-neighbor dGdr slabs in workspace state. Only the
         ! actually-populated entries in each slab need touching: for each
         ! group s the valid ranges are sf=1..n_symf_s and j=1..n_rad(s)
         ! (radial) or j=1..n_ang1_s, k=1..n_ang2_s (angular). The self
         ! contribution lives in self_dGdr(:, m).
         ASSOCIATE (workspace => nnp%neighbor_interface_state%workspace(ind), &
                    neighbor => nnp%neighbor_interface_state%workspace(ind)%neighbor, &
                    self_dGdr => nnp%neighbor_interface_state%workspace(ind)%self_dGdr)

            ! Radial groups
            DO s = 1, nnp%rad(ind)%n_symfgrp
               n_symf_s = nnp%rad(ind)%symfgrp(s)%n_symf
               ASSOCIATE (rad_buf => workspace%dGdr_rad(s)%data)
               DO sf = 1, n_symf_s
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  IF (nnp%scale_acsf) THEN
                     scale = (nnp%scmax - nnp%scmin)/ &
                             (nnp%rad(ind)%loc_max(m) - nnp%rad(ind)%loc_min(m))
                  ELSE
                     scale = (nnp%scmax - nnp%scmin)/nnp%rad(ind)%sigma(m)
                  END IF
                  self_dGdr(1, m) = self_dGdr(1, m)*scale
                  self_dGdr(2, m) = self_dGdr(2, m)*scale
                  self_dGdr(3, m) = self_dGdr(3, m)*scale
                  DO j = 1, neighbor%n_rad(s)
                     rad_buf(1, sf, j) = rad_buf(1, sf, j)*scale
                     rad_buf(2, sf, j) = rad_buf(2, sf, j)*scale
                     rad_buf(3, sf, j) = rad_buf(3, sf, j)*scale
                  END DO
               END DO
               END ASSOCIATE
            END DO

            ! Angular groups
            off = nnp%n_rad(ind)
            DO s = 1, nnp%ang(ind)%n_symfgrp
               n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
               n_ang1_s = neighbor%n_ang1(s)
               homo_grp = (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2))
               IF (homo_grp) THEN
                  ! kk slab is also indexed in ind_ang1 in the homo case
                  n_ang2_s = n_ang1_s
               ELSE
                  n_ang2_s = neighbor%n_ang2(s)
               END IF
               ASSOCIATE (jj_buf => workspace%dGdr_ang_jj(s)%data, &
                          kk_buf => workspace%dGdr_ang_kk(s)%data)
               DO sf = 1, n_symf_s
                  m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                  IF (nnp%scale_acsf) THEN
                     scale = (nnp%scmax - nnp%scmin)/ &
                             (nnp%ang(ind)%loc_max(m - off) - nnp%ang(ind)%loc_min(m - off))
                  ELSE
                     scale = (nnp%scmax - nnp%scmin)/nnp%ang(ind)%sigma(m - off)
                  END IF
                  self_dGdr(1, m) = self_dGdr(1, m)*scale
                  self_dGdr(2, m) = self_dGdr(2, m)*scale
                  self_dGdr(3, m) = self_dGdr(3, m)*scale
                  DO j = 1, n_ang1_s
                     jj_buf(1, sf, j) = jj_buf(1, sf, j)*scale
                     jj_buf(2, sf, j) = jj_buf(2, sf, j)*scale
                     jj_buf(3, sf, j) = jj_buf(3, sf, j)*scale
                  END DO
                  DO k = 1, n_ang2_s
                     kk_buf(1, sf, k) = kk_buf(1, sf, k)*scale
                     kk_buf(2, sf, k) = kk_buf(2, sf, k)*scale
                     kk_buf(3, sf, k) = kk_buf(3, sf, k)*scale
                  END DO
               END DO
               END ASSOCIATE
            END DO

         END ASSOCIATE
      END IF

      IF (PRESENT(stress)) THEN
         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)/(nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
                                 (nnp%scmax - nnp%scmin)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)/ &
                                       (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
                                       (nnp%scmax - nnp%scmin)
            END DO
         ELSE IF (nnp%scale_sigma_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)/nnp%rad(ind)%sigma(j)*(nnp%scmax - nnp%scmin)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)/nnp%ang(ind)%sigma(j)*(nnp%scmax - nnp%scmin)
            END DO
         END IF
      END IF

   END SUBROUTINE nnp_scale_acsf

! **************************************************************************************************
!> \brief Calculate radial symmetry function and gradient (optional)
! **************************************************************************************************
   SUBROUTINE nnp_calc_rad(nnp, ind, s, rvect, r, sym, force)

      TYPE(nnp_type), INTENT(IN), TARGET                 :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect
      REAL(KIND=dp), INTENT(IN)                          :: r
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: i, n_symf, sf
      REAL(KIND=dp)                                      :: dsymdr_full, r_inv, t, t2, t3, &
                                                            h00, h10, h01, h11, &
                                                            h10_dx, h11_dx, &
                                                            dh00, dh10, dh01, dh11, &
                                                            yi, yi1, dyi, dyi1, &
                                                            drdx_x, drdx_y, drdx_z
      TYPE(nnp_symfgrp_type), POINTER                    :: grp

      grp => nnp%rad(ind)%symfgrp(s)
      n_symf = grp%n_symf

      ! Group-shared Hermite cubic spline. All SFs in this radial group sit on
      ! the same uniform grid (because they share grp%cutoff), so the
      ! interpolation parameters (i, t, h00..h11, dh00..dh11) are computed
      ! ONCE outside the SF loop. The inner loop then issues 4 contiguous
      ! reads per SF from spline_y/spline_dy (sf-first layout) and 2-4
      ! multiply-adds. With n_symf = 8 typical, the inner loop is exactly
      ! one cache line per fetch and auto-vectorizes cleanly under
      ! gfortran -O3 / ifx -O2.
      !
      ! Out-of-range clamping: if r is past spline_x_max (which equals
      ! grp%cutoff), every SF in this group sees the boundary value of
      ! its tabulated y(r), which the build routine pinned to 0 because
      ! fcut(cutoff) = 0. So the clamp branch returns sym = 0 directly.
      IF (r >= grp%spline_x_max) THEN
         DO sf = 1, n_symf
            sym(sf) = 0.0_dp
         END DO
         IF (PRESENT(force)) THEN
            DO sf = 1, n_symf
               force(1, sf) = 0.0_dp
               force(2, sf) = 0.0_dp
               force(3, sf) = 0.0_dp
            END DO
         END IF
         RETURN
      END IF

      i = INT(r*grp%spline_dx_inv) + 1
      IF (i < 1) i = 1
      IF (i > grp%spline_n - 1) i = grp%spline_n - 1

      t = (r - REAL(i - 1, KIND=dp)*grp%spline_dx)*grp%spline_dx_inv
      t2 = t*t
      t3 = t2*t

      h00 = 2.0_dp*t3 - 3.0_dp*t2 + 1.0_dp
      h10 = t3 - 2.0_dp*t2 + t
      h01 = -2.0_dp*t3 + 3.0_dp*t2
      h11 = t3 - t2
      h10_dx = h10*grp%spline_dx
      h11_dx = h11*grp%spline_dx

      IF (PRESENT(force)) THEN
         dh00 = 6.0_dp*(t2 - t)
         dh10 = 3.0_dp*t2 - 4.0_dp*t + 1.0_dp
         dh01 = -dh00
         dh11 = 3.0_dp*t2 - 2.0_dp*t

         r_inv = 1.0_dp/r
         drdx_x = rvect(1)*r_inv
         drdx_y = rvect(2)*r_inv
         drdx_z = rvect(3)*r_inv

         ! Vectorizable inner loop: contiguous reads on sf, no branches,
         ! no function calls, no transcendentals.
         ASSOCIATE (spy_i => grp%spline_y(:, i), spy_i1 => grp%spline_y(:, i + 1), &
                    spdy_i => grp%spline_dy(:, i), spdy_i1 => grp%spline_dy(:, i + 1))
         !$OMP SIMD PRIVATE(yi, yi1, dyi, dyi1, dsymdr_full)
         DO sf = 1, n_symf
            yi = spy_i(sf)
            yi1 = spy_i1(sf)
            dyi = spdy_i(sf)
            dyi1 = spdy_i1(sf)
            sym(sf) = h00*yi + h10_dx*dyi + h01*yi1 + h11_dx*dyi1
            dsymdr_full = (dh00*yi + dh01*yi1)*grp%spline_dx_inv + dh10*dyi + dh11*dyi1
            force(1, sf) = dsymdr_full*drdx_x
            force(2, sf) = dsymdr_full*drdx_y
            force(3, sf) = dsymdr_full*drdx_z
         END DO
         END ASSOCIATE
      ELSE
         ASSOCIATE (spy_i => grp%spline_y(:, i), spy_i1 => grp%spline_y(:, i + 1), &
                    spdy_i => grp%spline_dy(:, i), spdy_i1 => grp%spline_dy(:, i + 1))
         !$OMP SIMD
         DO sf = 1, n_symf
            sym(sf) = h00*spy_i(sf) + h10_dx*spdy_i(sf) + &
                      h01*spy_i1(sf) + h11_dx*spdy_i1(sf)
         END DO
         END ASSOCIATE
      END IF

   END SUBROUTINE nnp_calc_rad

! **************************************************************************************************
!> \brief Build per-radial-group Hermite cubic spline tables.
!>        Each radial group s carries its own dense (n_symf, n_grid) tables
!>        spline_y and spline_dy. Within a group all SFs share the cutoff
!>        => same uniform grid, so the inner SF loop in nnp_calc_rad can
!>        compute the Hermite basis params ONCE and stream the per-SF y/dy
!>        values contiguously (sf-first layout).
!>
!>        Each table stores y(r) = EXP(-eta_sf*(r-rs_sf)^2)*fcut(r) and its
!>        analytic derivative y'(r). Memory cost: 2*n_symf*n_grid doubles
!>        per group, ~128 KB for n_symf=8, n_grid=8192. Built once per
!>        force eval; lives until nnp_env_release.
! **************************************************************************************************
   SUBROUTINE nnp_build_radial_splines(nnp)

      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp

      INTEGER                                            :: ind, k, n_symf, p, s, sf
      REAL(KIND=dp)                                      :: arg, cutoff, dfcutdr, dr, eta, &
                                                            exp_term, fcut, r, rs, tanh_tmp

      DO ind = 1, nnp%n_ele
         DO s = 1, nnp%rad(ind)%n_symfgrp
            ASSOCIATE (grp => nnp%rad(ind)%symfgrp(s))
               n_symf = grp%n_symf
               cutoff = grp%cutoff
               dr = cutoff/REAL(nnp_rad_spline_n - 1, KIND=dp)

               IF (ALLOCATED(grp%spline_y)) DEALLOCATE (grp%spline_y)
               IF (ALLOCATED(grp%spline_dy)) DEALLOCATE (grp%spline_dy)
               ALLOCATE (grp%spline_y(MAX(1, n_symf), nnp_rad_spline_n))
               ALLOCATE (grp%spline_dy(MAX(1, n_symf), nnp_rad_spline_n))
               grp%spline_n = nnp_rad_spline_n
               grp%spline_dx = dr
               grp%spline_dx_inv = 1.0_dp/dr
               grp%spline_x_max = cutoff

               ! Fill the table SF-by-SF. The grid loop is the inner one
               ! during build only; runtime nnp_calc_rad reads sf-first.
               DO sf = 1, n_symf
                  k = grp%symf(sf)
                  eta = nnp%rad(ind)%eta(k)
                  rs = nnp%rad(ind)%rs(k)

                  DO p = 1, nnp_rad_spline_n
                     r = REAL(p - 1, KIND=dp)*dr

                     SELECT CASE (nnp%cut_type)
                     CASE (nnp_cut_cos)
                        arg = pi*r/cutoff
                        fcut = 0.5_dp*(COS(arg) + 1.0_dp)
                        dfcutdr = -0.5_dp*SIN(arg)*(pi/cutoff)
                     CASE (nnp_cut_tanh)
                        tanh_tmp = TANH(1.0_dp - r/cutoff)
                        fcut = tanh_tmp**3
                        dfcutdr = (-3.0_dp/cutoff)*(tanh_tmp**2 - tanh_tmp**4)
                     CASE DEFAULT
                        CPABORT("NNP| Cutoff function unknown")
                     END SELECT

                     exp_term = EXP(-eta*(r - rs)**2)

                     grp%spline_y(sf, p) = exp_term*fcut
                     grp%spline_dy(sf, p) = exp_term*(-2.0_dp*eta*(r - rs))*fcut + &
                                            exp_term*dfcutdr
                  END DO
               END DO

               ! Pin the boundary so the runtime out-of-range branch can
               ! return zeros without re-checking each SF.
               DO sf = 1, n_symf
                  grp%spline_y(sf, nnp_rad_spline_n) = 0.0_dp
                  grp%spline_dy(sf, nnp_rad_spline_n) = 0.0_dp
               END DO

               grp%spline_built = .TRUE.
            END ASSOCIATE
         END DO
      END DO

   END SUBROUTINE nnp_build_radial_splines

! **************************************************************************************************
!> \brief Calculate angular symmetry function and gradient (optional)
!>
!> Vectorized SF-batched form. The original eta-sorted exp-skip loop has
!> been replaced by a sequence of SF passes:
!>   1. Branchless angular base: clamps tmp at the cusp (tmp <= 0), then
!>      computes tmpzeta via integer or real pow per SF. Stays scalar
!>      because the int/real switch and integer pow cannot SIMD.
!>   2. Vectorized EXP: symtmp_arr(sf) = EXP(-eta(sf)*r2sum). Marked
!>      !$OMP SIMD so the compiler maps it to libmvec/SVML vector EXP.
!>      Vector libm hides EXP latency across the SF loop without needing
!>      an eta-dedup pre-pass.
!>   3. SIMD sym scatter: sym(sf) = prefzeta*angular*symtmp*ftot.
!>   4. Force scatter (only when forces requested). The clamp in pass 1
!>      forces tmpzeta_arr(sf) = 0 at the cusp, which makes both angular
!>      AND pref_lam vanish naturally -- no per-iter branch.
!>
!> Caller passes precomputed cutoff values for the j/k legs; the j-k leg
!> is computed inline. Geometry-only derivative bases are hoisted outside
!> the SF loops as before.
! **************************************************************************************************
   SUBROUTINE nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, &
                            fcut_j, dfcut_j, fcut_k, dfcut_k, sym, force)

      TYPE(nnp_type), INTENT(IN), TARGET                 :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect1, rvect2, rvect3
      REAL(KIND=dp), INTENT(IN)                          :: r1, r2, r3
      REAL(KIND=dp), INTENT(IN)                          :: fcut_j, dfcut_j, fcut_k, dfcut_k
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: ii, izeta, n_symf, sf
      REAL(KIND=dp) :: angular, costheta, dfcutdr1, dfcutdr2, dfcutdr3, dsymdr1, dsymdr2, dsymdr3, &
         eta, f, fcut3, dfcut3, fcut_rc, ftot, g, inv_g2, lam, prefzeta, pref_lam, &
         r2sum, rsqr1, rsqr2, rsqr3, symtmp, &
         tmp, tmp1, tmp2, tmp3, tmpzeta, zeta, arg_tmp, tanh_tmp
      REAL(KIND=dp), DIMENSION(3) :: dcosbase1, dcosbase2, dcosbase3, &
         dr1dx, dr2dx, dr3dx, dgdx1, dgdx2, dgdx3
      LOGICAL                                            :: do_forces
      TYPE(nnp_symfgrp_type), POINTER                    :: grp

      ! Per-SF scratch (automatic arrays sized to the current group). For
      ! typical n2p2 ACSF n_symf is 8-32, so this is ~256-1024 bytes on
      ! stack per call. Holds the staged outputs of pass 1 + pass 2 so
      ! pass 3 streams them contiguously.
      REAL(KIND=dp), DIMENSION(nnp%ang(ind)%symfgrp(s)%n_symf) :: &
         tmpzeta_arr, angular_arr, symtmp_arr

      do_forces = PRESENT(force)
      grp => nnp%ang(ind)%symfgrp(s)
      n_symf = grp%n_symf
      fcut_rc = grp%cutoff

      rsqr1 = r1*r1
      rsqr2 = r2*r2
      rsqr3 = r3*r3
      r2sum = rsqr1 + rsqr2 + rsqr3

      f = rsqr3 - rsqr1 - rsqr2
      g = -2.0_dp*r1*r2
      costheta = f/g

      ! Compute fcut3 for r3 (j-k distance -- cannot be precomputed by caller)
      SELECT CASE (nnp%cut_type)
      CASE (nnp_cut_cos)
         arg_tmp = pi*r3/fcut_rc
         fcut3 = 0.5_dp*(COS(arg_tmp) + 1.0_dp)
         IF (do_forces) dfcut3 = -0.5_dp*SIN(arg_tmp)*(pi/fcut_rc)
      CASE (nnp_cut_tanh)
         tanh_tmp = TANH(1.0_dp - r3/fcut_rc)
         fcut3 = tanh_tmp**3
         IF (do_forces) dfcut3 = (-3.0_dp/fcut_rc)*(tanh_tmp**2 - tanh_tmp**4)
      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      ! Use precomputed fcut values for j and k neighbors
      ftot = fcut_j*fcut_k*fcut3

      IF (do_forces) THEN
         ! Combined cutoff derivatives (product rule)
         dfcutdr1 = dfcut_j*fcut_k*fcut3
         dfcutdr2 = fcut_j*dfcut_k*fcut3
         dfcutdr3 = fcut_j*fcut_k*dfcut3

         dr1dx(:) = rvect1(:)/r1
         dr2dx(:) = rvect2(:)/r2
         dr3dx(:) = rvect3(:)/r3

         ! Hoist geometry-only parts of costheta derivatives outside SF loop.
         ! Full angular derivative: dangulardx = (zeta*tmpzeta*lam/g^2) * dcosbase
         ! where dcosbase factors out lam from dfdx and the f*dgdx term.
         ! inv_g2 = 1/g^2 is hoisted: one divide per triplet, one mul per SF.
         inv_g2 = 1.0_dp/(g*g)
         DO ii = 1, 3
            tmp1 = 2.0_dp*r2*dr1dx(ii)
            tmp2 = 2.0_dp*r1*dr2dx(ii)
            dgdx1(ii) = -(tmp1 + tmp2)
            dgdx2(ii) = tmp1
            dgdx3(ii) = tmp2

            dcosbase1(ii) = -2.0_dp*(rvect1(ii) + rvect2(ii))*g - f*dgdx1(ii)
            dcosbase2(ii) = 2.0_dp*(rvect3(ii) + rvect1(ii))*g - f*dgdx2(ii)
            dcosbase3(ii) = 2.0_dp*(rvect2(ii) - rvect3(ii))*g - f*dgdx3(ii)
         END DO
      ELSE
         inv_g2 = 0.0_dp
      END IF

      ! ---- Pass 1: angular base + tmpzeta (scalar tight loop) ----
      ! Branchless cusp handling: when tmp <= 0 we set tmpzeta_arr = 0, which
      ! propagates to both angular_arr and pref_lam in pass 3, giving
      ! sym = 0 and force = 0 without an explicit CYCLE branch.
      DO sf = 1, n_symf
         tmp = 1.0_dp + grp%pack_lam(sf)*costheta
         IF (tmp <= 0.0_dp) THEN
            tmpzeta_arr(sf) = 0.0_dp
            angular_arr(sf) = 0.0_dp
         ELSE
            IF (grp%pack_use_int_zeta(sf)) THEN
               izeta = grp%pack_izeta(sf)
               tmpzeta_arr(sf) = tmp**(izeta - 1)
            ELSE
               tmpzeta_arr(sf) = tmp**(grp%pack_zeta(sf) - 1.0_dp)
            END IF
            angular_arr(sf) = tmpzeta_arr(sf)*tmp
         END IF
      END DO

      ! ---- Pass 2: vectorized EXP via libmvec/SVML ----
      ! Single SIMD loop computing symtmp_arr(sf) = EXP(-eta(sf)*r2sum).
      ! With -fopenmp-simd this maps to vector libm (e.g. _ZGVdN4v_exp on
      ! AVX2, _ZGVeN8v_exp on AVX-512), hiding EXP latency across 4-8 SFs
      ! per iter -- ~10 cycles/element vs ~30 for scalar EXP, which beats
      ! any eta-dedup pre-pass.
      !$OMP SIMD
      DO sf = 1, n_symf
         symtmp_arr(sf) = EXP(-grp%pack_eta(sf)*r2sum)
      END DO

      ! ---- Pass 3: sym scatter (vectorizable) ----
      !$OMP SIMD
      DO sf = 1, n_symf
         sym(sf) = grp%pack_prefzeta(sf)*angular_arr(sf)*symtmp_arr(sf)*ftot
      END DO

      ! ---- Pass 4: force scatter (only when do_forces) ----
      IF (do_forces) THEN
         DO sf = 1, n_symf
            symtmp = symtmp_arr(sf)
            angular = angular_arr(sf)
            tmpzeta = tmpzeta_arr(sf)
            eta = grp%pack_eta(sf)
            lam = grp%pack_lam(sf)
            zeta = grp%pack_zeta(sf)
            prefzeta = grp%pack_prefzeta(sf)

            ! pref_lam carries tmpzeta, so it vanishes at the cusp without
            ! an explicit branch (tmpzeta_arr was clamped to 0 in pass 1).
            pref_lam = zeta*tmpzeta*lam*inv_g2

            tmp = -2.0_dp*symtmp*eta
            dsymdr1 = tmp*r1
            dsymdr2 = tmp*r2
            dsymdr3 = tmp*r3

            tmp = prefzeta*symtmp*ftot
            tmp1 = prefzeta*angular*(ftot*dsymdr1 + dfcutdr1*symtmp)
            tmp2 = prefzeta*angular*(ftot*dsymdr2 + dfcutdr2*symtmp)
            tmp3 = prefzeta*angular*(ftot*dsymdr3 + dfcutdr3*symtmp)
            DO ii = 1, 3
               force(ii, 1, sf) = tmp*pref_lam*dcosbase1(ii) + tmp1*dr1dx(ii) + tmp2*dr2dx(ii)
               force(ii, 2, sf) = tmp*pref_lam*dcosbase2(ii) - tmp1*dr1dx(ii) + tmp3*dr3dx(ii)
               force(ii, 3, sf) = tmp*pref_lam*dcosbase3(ii) - tmp2*dr2dx(ii) - tmp3*dr3dx(ii)
            END DO
         END DO
      END IF

   END SUBROUTINE nnp_calc_ang

! **************************************************************************************************
!> \brief Sort an (ele, nuc_ele) pair of arrays in ascending order of atomic number.
!>        Used to canonicalise element ordering inside the NNP environment so the
!>        same model file produces the same per-element index layout regardless
!>        of input ordering.
!> \param ele      element-symbol array, sorted in place to match nuc_ele.
!> \param nuc_ele  per-element atomic number, sorted in place.
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_ele(ele, nuc_ele)
      CHARACTER(len=2), DIMENSION(:), INTENT(INOUT)      :: ele
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: nuc_ele

      CHARACTER(len=2)                                   :: tmp_ele
      INTEGER                                            :: i, j, loc, minimum, tmp_nuc_ele

      DO i = 1, SIZE(ele)
         CALL get_ptable_info(ele(i), number=nuc_ele(i))
      END DO

      DO i = 1, SIZE(ele) - 1
         minimum = nuc_ele(i)
         loc = i
         DO j = i + 1, SIZE(ele)
            IF (nuc_ele(j) < minimum) THEN
               loc = j
               minimum = nuc_ele(j)
            END IF
         END DO
         tmp_nuc_ele = nuc_ele(i)
         nuc_ele(i) = nuc_ele(loc)
         nuc_ele(loc) = tmp_nuc_ele

         tmp_ele = ele(i)
         ele(i) = ele(loc)
         ele(loc) = tmp_ele
      END DO

   END SUBROUTINE nnp_sort_ele

! **************************************************************************************************
!> \brief Sort radial and angular symmetry functions in canonical order.
!>        Radial SFs are sorted by eta (ascending) then rcut; angular SFs by
!>        eta, lambda, zeta. This is the order downstream code in nnp_init_acsf_groups
!>        relies on for run-length-style group packing.
!> \param nnp NNP environment whose rad/ang SF arrays will be reordered in place.
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_acsf(nnp)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: i, j, k, loc

      DO i = 1, nnp%n_ele
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) > nnp%rad(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) > nnp%rad(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) == nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) > nnp%rad(i)%rs(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) == nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) == nnp%rad(i)%rs(k) .AND. &
                   nnp%rad(i)%nuc_ele(loc) > nnp%rad(i)%nuc_ele(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) > nnp%ang(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) > nnp%ang(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) > nnp%ang(i)%zeta(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) > nnp%ang(i)%lam(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) == nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) > nnp%ang(i)%nuc_ele1(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) == nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) == nnp%ang(i)%nuc_ele1(k) .AND. &
                   nnp%ang(i)%nuc_ele2(loc) > nnp%ang(i)%nuc_ele2(k)) THEN
                  loc = k
               END IF
            END DO
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO
      END DO

   END SUBROUTINE nnp_sort_acsf

! **************************************************************************************************
!> \brief Swap two radial symmetry functions
! **************************************************************************************************
   SUBROUTINE nnp_swaprad(rad, i, j)
      TYPE(nnp_acsf_rad_type), INTENT(INOUT)             :: rad
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = rad%funccut(i)
      rad%funccut(i) = rad%funccut(j)
      rad%funccut(j) = tmpr

      tmpr = rad%eta(i)
      rad%eta(i) = rad%eta(j)
      rad%eta(j) = tmpr

      tmpr = rad%rs(i)
      rad%rs(i) = rad%rs(j)
      rad%rs(j) = tmpr

      tmpc = rad%ele(i)
      rad%ele(i) = rad%ele(j)
      rad%ele(j) = tmpc

      tmpi = rad%nuc_ele(i)
      rad%nuc_ele(i) = rad%nuc_ele(j)
      rad%nuc_ele(j) = tmpi

   END SUBROUTINE nnp_swaprad

! **************************************************************************************************
!> \brief Swap two angular symmetry functions
! **************************************************************************************************
   SUBROUTINE nnp_swapang(ang, i, j)
      TYPE(nnp_acsf_ang_type), INTENT(INOUT)             :: ang
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = ang%funccut(i)
      ang%funccut(i) = ang%funccut(j)
      ang%funccut(j) = tmpr

      tmpr = ang%eta(i)
      ang%eta(i) = ang%eta(j)
      ang%eta(j) = tmpr

      tmpr = ang%zeta(i)
      ang%zeta(i) = ang%zeta(j)
      ang%zeta(j) = tmpr

      tmpr = ang%prefzeta(i)
      ang%prefzeta(i) = ang%prefzeta(j)
      ang%prefzeta(j) = tmpr

      tmpr = ang%lam(i)
      ang%lam(i) = ang%lam(j)
      ang%lam(j) = tmpr

      tmpc = ang%ele1(i)
      ang%ele1(i) = ang%ele1(j)
      ang%ele1(j) = tmpc

      tmpi = ang%nuc_ele1(i)
      ang%nuc_ele1(i) = ang%nuc_ele1(j)
      ang%nuc_ele1(j) = tmpi

      tmpc = ang%ele2(i)
      ang%ele2(i) = ang%ele2(j)
      ang%ele2(j) = tmpc

      tmpi = ang%nuc_ele2(i)
      ang%nuc_ele2(i) = ang%nuc_ele2(j)
      ang%nuc_ele2(j) = tmpi

   END SUBROUTINE nnp_swapang

! **************************************************************************************************
!> \brief Pack symmetry functions into groups that share input parameters.
!>        Builds nnp%rad(i)%symfgrp(:) and nnp%ang(i)%symfgrp(:) so that all
!>        radial / angular SFs with identical (eta, rcut [, lambda, zeta])
!>        live in one group; downstream the descriptor pass evaluates each
!>        group's shared cutoff/exponent once and applies it across the
!>        per-element neighbour slab.
!> \param nnp NNP environment with rad/ang SF metadata already populated by nnp_init_model.
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_acsf_groups(nnp)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: ang, i, j, k, rad, s, sf, m, n_symf, izeta_tmp
      REAL(KIND=dp)                                      :: funccut, eta_tmp, zeta_tmp

      DO i = 1, nnp%n_ele
         nnp%rad(i)%n_symfgrp = 0
         nnp%ang(i)%n_symfgrp = 0
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > cutoff_eq_tol) THEN
                     nnp%rad(i)%n_symfgrp = nnp%rad(i)%n_symfgrp + 1
                     funccut = nnp%rad(i)%funccut(s)
                  END IF
               END IF
            END DO
         END DO
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > cutoff_eq_tol) THEN
                        nnp%ang(i)%n_symfgrp = nnp%ang(i)%n_symfgrp + 1
                        funccut = nnp%ang(i)%funccut(s)
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

      DO i = 1, nnp%n_ele
         ALLOCATE (nnp%rad(i)%symfgrp(nnp%rad(i)%n_symfgrp))
         ALLOCATE (nnp%ang(i)%symfgrp(nnp%ang(i)%n_symfgrp))
         DO j = 1, nnp%rad(i)%n_symfgrp
            nnp%rad(i)%symfgrp(j)%n_symf = 0
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            nnp%ang(i)%symfgrp(j)%n_symf = 0
         END DO
      END DO

      DO i = 1, nnp%n_ele
         rad = 0
         ang = 0
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > cutoff_eq_tol) THEN
                     rad = rad + 1
                     funccut = nnp%rad(i)%funccut(s)
                     nnp%rad(i)%symfgrp(rad)%cutoff = funccut
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele(1))
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele_ind(1))
                     nnp%rad(i)%symfgrp(rad)%ele(1) = nnp%ele(j)
                     nnp%rad(i)%symfgrp(rad)%ele_ind(1) = j
                  END IF
                  nnp%rad(i)%symfgrp(rad)%n_symf = nnp%rad(i)%symfgrp(rad)%n_symf + 1
               END IF
            END DO
         END DO
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > cutoff_eq_tol) THEN
                        ang = ang + 1
                        funccut = nnp%ang(i)%funccut(s)
                        nnp%ang(i)%symfgrp(ang)%cutoff = funccut
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele(2))
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele_ind(2))
                        nnp%ang(i)%symfgrp(ang)%ele(1) = nnp%ele(j)
                        nnp%ang(i)%symfgrp(ang)%ele(2) = nnp%ele(k)
                        nnp%ang(i)%symfgrp(ang)%ele_ind(1) = j
                        nnp%ang(i)%symfgrp(ang)%ele_ind(2) = k
                     END IF
                     nnp%ang(i)%symfgrp(ang)%n_symf = nnp%ang(i)%symfgrp(ang)%n_symf + 1
                  END IF
               END DO
            END DO
         END DO
      END DO

      DO i = 1, nnp%n_ele
         DO j = 1, nnp%rad(i)%n_symfgrp
            ALLOCATE (nnp%rad(i)%symfgrp(j)%symf(nnp%rad(i)%symfgrp(j)%n_symf))
            rad = 0
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%rad(i)%symfgrp(j)%ele(1)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - nnp%rad(i)%symfgrp(j)%cutoff) <= cutoff_eq_tol) THEN
                     rad = rad + 1
                     nnp%rad(i)%symfgrp(j)%symf(rad) = s
                  END IF
               END IF
            END DO
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            ALLOCATE (nnp%ang(i)%symfgrp(j)%symf(nnp%ang(i)%symfgrp(j)%n_symf))
            ang = 0
            DO s = 1, nnp%n_ang(i)
               IF ((nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(1) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(2)) .OR. &
                   (nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(2) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(1))) THEN
                  IF (ABS(nnp%ang(i)%funccut(s) - nnp%ang(i)%symfgrp(j)%cutoff) <= cutoff_eq_tol) THEN
                     ang = ang + 1
                     nnp%ang(i)%symfgrp(j)%symf(ang) = s
                  END IF
               END IF
            END DO
         END DO
      END DO

      ! Populate packed parameter arrays for the angular SF group inner loop
      ! so it reads contiguous memory rather than chasing
      ! ang(i)%{eta,zeta,lam,prefzeta}(symf(sf)) on every iteration.
      ! pack_use_int_zeta and pack_izeta let the inner loop skip NINT and
      ! dispatch a tight integer power when zeta is integral. The radial
      ! group does not need packed parameters since nnp_calc_rad reads
      ! pretabulated splines.
      DO i = 1, nnp%n_ele
         DO j = 1, nnp%ang(i)%n_symfgrp
            n_symf = nnp%ang(i)%symfgrp(j)%n_symf
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_eta(n_symf))
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_zeta(n_symf))
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_lam(n_symf))
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_prefzeta(n_symf))
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_izeta(n_symf))
            ALLOCATE (nnp%ang(i)%symfgrp(j)%pack_use_int_zeta(n_symf))
            DO sf = 1, n_symf
               m = nnp%ang(i)%symfgrp(j)%symf(sf)
               eta_tmp = nnp%ang(i)%eta(m)
               zeta_tmp = nnp%ang(i)%zeta(m)
               nnp%ang(i)%symfgrp(j)%pack_eta(sf) = eta_tmp
               nnp%ang(i)%symfgrp(j)%pack_zeta(sf) = zeta_tmp
               nnp%ang(i)%symfgrp(j)%pack_lam(sf) = nnp%ang(i)%lam(m)
               nnp%ang(i)%symfgrp(j)%pack_prefzeta(sf) = nnp%ang(i)%prefzeta(m)
               izeta_tmp = NINT(zeta_tmp)
               nnp%ang(i)%symfgrp(j)%pack_izeta(sf) = izeta_tmp
               nnp%ang(i)%symfgrp(j)%pack_use_int_zeta(sf) = &
                  (REAL(izeta_tmp, dp) == zeta_tmp)
            END DO
         END DO
      END DO

   END SUBROUTINE nnp_init_acsf_groups

! **************************************************************************************************
!> \brief Print a summary of the active symmetry-function set on the source rank.
!>        Emits one line per element listing the per-element n_rad / n_ang counts
!>        and the per-group breakdown; used at NNP init time and after sort/group
!>        passes for traceability.
!> \param nnp       NNP environment whose SF metadata is to be printed.
!> \param para_env  parallel environment; only the source rank emits output.
!> \param printtag  log-line prefix (typically "NNP" or "HELIUM_NNP").
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_write_acsf(nnp, para_env, printtag)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: printtag

      CHARACTER(len=default_string_length)               :: my_label
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      my_label = TRIM(printtag)//"| "
      IF (para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, '(1X,A,1X,10(I2,1X))') TRIM(my_label)//" Activation functions:", nnp%actfnct(:)
         DO i = 1, nnp%n_ele
            WRITE (unit_nr, *) TRIM(my_label)//" short range atomic symmetry functions element "// &
               nnp%ele(i)//":"
            DO j = 1, nnp%n_rad(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,1X,A2,11X,3(F6.3,1X))') TRIM(my_label), j, nnp%ele(i), 2, &
                  nnp%rad(i)%ele(j), nnp%rad(i)%eta(j), &
                  nnp%rad(i)%rs(j), nnp%rad(i)%funccut(j)
            END DO
            DO j = 1, nnp%n_ang(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,2(1X,A2),1X,4(F6.3,1X))') &
                  TRIM(my_label), j, nnp%ele(i), 3, &
                  nnp%ang(i)%ele1(j), nnp%ang(i)%ele2(j), &
                  nnp%ang(i)%eta(j), nnp%ang(i)%lam(j), &
                  nnp%ang(i)%zeta(j), nnp%ang(i)%funccut(j)
            END DO
         END DO
      END IF

   END SUBROUTINE nnp_write_acsf

END MODULE nnp_acsf
