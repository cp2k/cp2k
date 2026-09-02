!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  GPU offload driver for the NNP force evaluation. Hosts the
!>         ISO_C_BINDING interfaces to the CUDA launchers in src/nnp/gpu/ and
!>         the whole-step device pipeline called from nnp_force. All GPU
!>         entry points compile to inert stubs unless the build enables the
!>         NNP GPU backend (__OFFLOAD without __NO_OFFLOAD_NNP).
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
MODULE nnp_gpu

   USE ISO_C_BINDING,                   ONLY: c_char,&
                                              c_double,&
                                              c_int,&
                                              c_null_ptr
   USE input_constants,                 ONLY: nnp_gpu_prec_mixed
   USE kinds,                           ONLY: dp
   USE nnp_acsf,                        ONLY: nnp_build_radial_splines
   USE nnp_environment_types,           ONLY: nnp_type
   USE offload_api,                     ONLY: offload_activate_chosen_device,&
                                              offload_get_device_count,&
                                              offload_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_gpu'

   ! Initial per-element neighbour-list width for the device walk; grown and
   ! retried on overflow, then reused across MD steps.
   INTEGER, PARAMETER, PRIVATE :: gpu_walk_init_max_n_neigh = 96

   ! nnp_gpu_run_walk, nnp_gpu_build_descriptor and nnp_gpu_predict stay private:
   ! each is only valid after the previous one has run for this step, and
   ! nnp_gpu_energy_force is the entry point that guarantees that order.
   PUBLIC :: nnp_gpu_available
   PUBLIC :: nnp_gpu_init_device
   PUBLIC :: nnp_gpu_model_capacity_ok
   PUBLIC :: nnp_gpu_upload_descriptor_statics
   PUBLIC :: nnp_gpu_upload_network
   PUBLIC :: nnp_gpu_energy_force

   ! Every device entry point below is compiled only for accelerator builds.
   ! The guard tests the backend macros the build system defines globally
   ! (__OFFLOAD_CUDA / __OFFLOAD_HIP) alongside the __OFFLOAD umbrella, so it
   ! reads the same in every translation unit that carries it, including
   ! cp2k_info.F, which does not include base_uses.f90 and therefore never
   ! sees a Fortran-side __OFFLOAD definition.
#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
   INTERFACE
      FUNCTION nnp_gpu_set_precision_c(use_fp32) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_set_precision")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int), VALUE                         :: use_fp32
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_set_precision_c

      ! Counts this environment as a user of the per-process device state and
      ! returns the number now live.
      FUNCTION nnp_gpu_register_env_c() &
         RESULT(n_live) BIND(C, NAME="nnp_gpu_register_env")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int)                                :: n_live

      END FUNCTION nnp_gpu_register_env_c

      ! Undoes a registration the caller decided not to honour, leaving the
      ! state the owning environment is using untouched.
      SUBROUTINE nnp_gpu_unregister_env_c() &
         BIND(C, NAME="nnp_gpu_unregister_env")
      END SUBROUTINE nnp_gpu_unregister_env_c

      ! Compiled-in capacities of the per-(group, element) device buffers and
      ! of the run-constant parameter cache, for the model check at init.
      FUNCTION nnp_gpu_max_se_c() &
         RESULT(cap) BIND(C, NAME="nnp_gpu_max_se_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int)                                :: cap

      END FUNCTION nnp_gpu_max_se_c

      FUNCTION nnp_gpu_spline_cache_max_c() &
         RESULT(cap) BIND(C, NAME="nnp_gpu_spline_cache_max_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int)                                :: cap

      END FUNCTION nnp_gpu_spline_cache_max_c

      FUNCTION nnp_gpu_upload_positions_step_c(num_atoms, h_coord) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_positions_step_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: num_atoms
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_coord
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_positions_step_c

      FUNCTION nnp_gpu_upload_image_table_c(n_images, h_image_table) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_image_table_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_images
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_image_table
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_image_table_c

      FUNCTION nnp_gpu_upload_cell_list_c(n_images, n_cells, h_image_atom, &
                                          h_image_shift, h_head, h_next) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_cell_list_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int), VALUE                         :: n_images, n_cells
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_image_atom, h_image_shift, h_head, &
                                                            h_next
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_cell_list_c

      FUNCTION nnp_gpu_upload_pair_map_c(n_ele, max_pair_grps, max_n_radgrp, &
                                         max_n_anggrp, h_n_rad, h_n_ang1, h_n_ang2, &
                                         h_pair_max_cutoff, h_rad_groups, h_ang1_groups, &
                                         h_ang2_groups, h_rad_cutoff, h_ang_cutoff) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_pair_map_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_ele, max_pair_grps, max_n_radgrp, &
                                                            max_n_anggrp
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_n_rad, h_n_ang1, h_n_ang2
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_pair_max_cutoff
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_rad_groups, h_ang1_groups, &
                                                            h_ang2_groups
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_rad_cutoff, h_ang_cutoff
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_pair_map_c

      FUNCTION nnp_gpu_cell_list_walk_c(num_atoms_global, n_ele, n_images, n_cells, &
                                        n_atoms_in_e_max, max_n_neigh_rad, &
                                        max_n_neigh_ang1, max_n_neigh_ang2, max_n_radgrp, &
                                        max_n_anggrp, h_nbin, h_bin_lower, h_bin_width, &
                                        h_bin_span, h_perd, h_exact_pbc_copies, &
                                        h_coord_scaled, h_ele_ind_global, &
                                        h_atoms_global_per_e, h_n_atoms_in_e) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_cell_list_walk_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE :: num_atoms_global, n_ele, n_images, n_cells, n_atoms_in_e_max, &
         max_n_neigh_rad, max_n_neigh_ang1, max_n_neigh_ang2, max_n_radgrp, max_n_anggrp
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_nbin
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_bin_lower, h_bin_width
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_bin_span, h_perd, h_exact_pbc_copies
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_coord_scaled
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_ele_ind_global, h_atoms_global_per_e, &
                                                            h_n_atoms_in_e
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_cell_list_walk_c

      FUNCTION nnp_gpu_walk_counts_d2h_c(n_atoms_in_e_max, max_n_radgrp, max_n_anggrp, &
                                         n_ele, h_n_rad_pa, h_n_ang1_pa, h_n_ang2_pa) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_walk_counts_d2h_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int), VALUE                         :: n_atoms_in_e_max, max_n_radgrp, &
                                                            max_n_anggrp, n_ele
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(INOUT)   :: h_n_rad_pa, h_n_ang1_pa, h_n_ang2_pa
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_walk_counts_d2h_c

      FUNCTION nnp_gpu_begin_descriptor_step_c(n_atoms, n_input_max) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_begin_descriptor_step")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int), VALUE                         :: n_atoms, n_input_max
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_begin_descriptor_step_c

      FUNCTION nnp_gpu_upload_scaling_tables_c(n_ele, n_input_max, h_loc_av, h_loc_min, &
                                               h_loc_max, h_sigma, h_n_input_nodes) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_scaling_tables_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_ele, n_input_max
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_loc_av, h_loc_min, h_loc_max, h_sigma
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_n_input_nodes
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_scaling_tables_c

      FUNCTION nnp_gpu_upload_se_static_tables_c(n_ele, max_n_radgrp, max_n_anggrp, &
                                                 total_rad_symf, total_ang_symf, &
                                                 h_rad_packed, h_rad_offsets, &
                                                 h_ang_packed, h_ang_offsets) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_se_static_tables_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int), VALUE                         :: n_ele, max_n_radgrp, max_n_anggrp, &
                                                            total_rad_symf, total_ang_symf
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_rad_packed, h_rad_offsets, &
                                                            h_ang_packed, h_ang_offsets
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_se_static_tables_c

      FUNCTION nnp_gpu_scale_G_dev_c(n_atoms, n_input_max, n_ele, scaling_mode, &
                                     scmin, scmax, h_ele_ind_for_atoms) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_scale_G_dev_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_atoms, n_input_max, n_ele, scaling_mode
      REAL(KIND=c_double), VALUE                         :: scmin, scmax
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_ele_ind_for_atoms
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_scale_G_dev_c

      ! Whether the last scaling pass saw a descriptor outside its training
      ! range; valid once the step's device work has been synchronized.
      FUNCTION nnp_gpu_extrapolation_flag_c() &
         RESULT(flag) BIND(C, NAME="nnp_gpu_extrapolation_flag_c")
         USE ISO_C_BINDING, ONLY: c_int
      INTEGER(KIND=c_int)                                :: flag

      END FUNCTION nnp_gpu_extrapolation_flag_c

      FUNCTION nnp_gpu_radial_sf_to_G_dev_c(n_atoms_in_e, max_n_neigh_per_se, n_symf, &
                                            spline_n, spline_dx, spline_dx_inv, &
                                            spline_x_max, h_spline_y, h_spline_dy, &
                                            s_0based, e_0based, max_n_radgrp_walk, &
                                            n_atoms_in_e_max_walk, max_n_neigh_walk, &
                                            istart_1based, n_input_max, se_idx) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_radial_sf_to_G_dev_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_atoms_in_e, max_n_neigh_per_se, &
                                                            n_symf, spline_n
      REAL(KIND=c_double), VALUE                         :: spline_dx, spline_dx_inv, spline_x_max
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_spline_y, h_spline_dy
      INTEGER(KIND=c_int), VALUE :: s_0based, e_0based, max_n_radgrp_walk, n_atoms_in_e_max_walk, &
         max_n_neigh_walk, istart_1based, n_input_max, se_idx
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_radial_sf_to_G_dev_c

      FUNCTION nnp_gpu_angular_sf_to_G_dev_c(homo, cut_type, max_n_ang1_per_se, &
                                             max_n_ang2_per_se, max_n_kk_per_se, n_symf, &
                                             n_atoms_in_e, cutoff_s, cutoff_sqr, &
                                             h_pack_lam, h_pack_zeta, h_pack_eta, &
                                             h_pack_prefzeta, h_pack_izeta, &
                                             h_pack_use_int_zeta, s_0based, e_0based, &
                                             max_n_anggrp_walk, n_atoms_in_e_max_walk, &
                                             max_n_neigh_ang1_walk, max_n_neigh_ang2_walk, &
                                             istart_1based, n_rad_for_e, n_input_max, se_idx) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_angular_sf_to_G_dev_c")
         USE ISO_C_BINDING, ONLY: c_char, c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: homo, cut_type, max_n_ang1_per_se, &
                                                            max_n_ang2_per_se, max_n_kk_per_se, &
                                                            n_symf, n_atoms_in_e
      REAL(KIND=c_double), VALUE                         :: cutoff_s, cutoff_sqr
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_pack_lam, h_pack_zeta, h_pack_eta, &
                                                            h_pack_prefzeta
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_pack_izeta
      CHARACTER(KIND=c_char), DIMENSION(*), INTENT(IN)   :: h_pack_use_int_zeta
      INTEGER(KIND=c_int), VALUE :: s_0based, e_0based, max_n_anggrp_walk, n_atoms_in_e_max_walk, &
         max_n_neigh_ang1_walk, max_n_neigh_ang2_walk, istart_1based, n_rad_for_e, n_input_max, &
         se_idx
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_angular_sf_to_G_dev_c

      FUNCTION nnp_gpu_upload_nn_packed_c(n_committee, n_layer, max_n_nodes, n_ele, &
                                          h_n_nodes, h_act_fnct, h_atom_energies, &
                                          h_W_packed, h_b_packed) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_upload_nn_packed_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_committee, n_layer, max_n_nodes, n_ele
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_n_nodes, h_act_fnct
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_atom_energies, h_W_packed, h_b_packed
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_upload_nn_packed_c

      FUNCTION nnp_gpu_nn_atomwise_with_G_dev_c(n_atoms, n_committee, n_layer, max_n_nodes, &
                                                n_input_max, normnodes, n_ele, h_ele_ind, &
                                                h_atomic_energy, h_dE_dG) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_nn_atomwise_with_G_dev_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int, c_ptr
      INTEGER(KIND=c_int), VALUE                         :: n_atoms, n_committee, n_layer, &
                                                            max_n_nodes, n_input_max, normnodes, &
                                                            n_ele
      INTEGER(KIND=c_int), DIMENSION(*), INTENT(IN)      :: h_ele_ind
      REAL(KIND=c_double), DIMENSION(*), INTENT(INOUT)   :: h_atomic_energy
      TYPE(c_ptr), VALUE                                 :: h_dE_dG
      INTEGER(KIND=c_int)                                :: rc

! Optional descriptor-gradient read-back; C_NULL_PTR leaves the
! gradient on the device for the force scatter to consume.

      END FUNCTION nnp_gpu_nn_atomwise_with_G_dev_c

      FUNCTION nnp_gpu_force_scatter_begin_c(num_atoms_global, n_committee, h_force) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_force_scatter_begin_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: num_atoms_global, n_committee
      REAL(KIND=c_double), DIMENSION(*), INTENT(IN)      :: h_force
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_force_scatter_begin_c

      FUNCTION nnp_gpu_force_scatter_end_c(num_atoms_global, n_committee, h_force) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_force_scatter_end_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: num_atoms_global, n_committee
      REAL(KIND=c_double), DIMENSION(*), INTENT(INOUT)   :: h_force
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_force_scatter_end_c

      ! D2H the per-committee analytical virial accumulator (9 * n_committee
      ! doubles, layout [9*c + 3*a + b]). Call after nnp_gpu_force_scatter_end_c,
      ! whose stream sync has committed the scatter kernels' atomic adds.
      FUNCTION nnp_gpu_force_scatter_virial_d2h_c(n_committee, h_virial) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_force_scatter_virial_d2h_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_committee
      REAL(KIND=c_double), DIMENSION(*), INTENT(INOUT)   :: h_virial
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_force_scatter_virial_d2h_c

      FUNCTION nnp_gpu_force_scatter_radial_c(n_atoms_in_e, n_symf, max_n_neigh_per_se, &
                                              n_committee, n_input_max, num_atoms_global, &
                                              ele_idx, scaling_mode, scmin, scmax, s_0based, &
                                              max_n_radgrp_walk, n_atoms_in_e_max_walk, &
                                              max_n_neigh_walk, istart_1based, se_idx, &
                                              accumulate_virial) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_force_scatter_radial_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE                         :: n_atoms_in_e, n_symf, &
                                                            max_n_neigh_per_se, n_committee, &
                                                            n_input_max, num_atoms_global, &
                                                            ele_idx, scaling_mode
      REAL(KIND=c_double), VALUE                         :: scmin, scmax
      INTEGER(KIND=c_int), VALUE :: s_0based, max_n_radgrp_walk, n_atoms_in_e_max_walk, &
         max_n_neigh_walk, istart_1based, se_idx, accumulate_virial
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_force_scatter_radial_c

      FUNCTION nnp_gpu_force_scatter_angular_c(n_atoms_in_e, n_symf, max_n_ang1_per_se, &
                                               max_n_kk_per_se, n_committee, n_input_max, &
                                               num_atoms_global, ele_idx, n_rad_offset, &
                                               scaling_mode, homo, scmin, scmax, s_0based, &
                                               max_n_anggrp_walk, n_atoms_in_e_max_walk, &
                                               max_n_neigh_ang1_walk, max_n_neigh_ang2_walk, &
                                               istart_1based, se_idx, accumulate_virial) &
         RESULT(rc) BIND(C, NAME="nnp_gpu_force_scatter_angular_c")
         USE ISO_C_BINDING, ONLY: c_double, c_int
      INTEGER(KIND=c_int), VALUE :: n_atoms_in_e, n_symf, max_n_ang1_per_se, max_n_kk_per_se, &
         n_committee, n_input_max, num_atoms_global, ele_idx, n_rad_offset, scaling_mode, homo
      REAL(KIND=c_double), VALUE                         :: scmin, scmax
      INTEGER(KIND=c_int), VALUE :: s_0based, max_n_anggrp_walk, n_atoms_in_e_max_walk, &
         max_n_neigh_ang1_walk, max_n_neigh_ang2_walk, istart_1based, se_idx, accumulate_virial
      INTEGER(KIND=c_int)                                :: rc

      END FUNCTION nnp_gpu_force_scatter_angular_c
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Report whether the NNP GPU backend can run: the build must enable
!>        it and at least one offload device must be visible. Builds without
!>        the backend always report .FALSE., so USE_GPU AUTO resolves to the
!>        CPU path there.
!> \return .TRUE. when the GPU code path is usable on this rank
! **************************************************************************************************
   FUNCTION nnp_gpu_available() RESULT(available)
      LOGICAL                                            :: available

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CALL offload_init()
      available = (offload_get_device_count() >= 1)
#else
      available = .FALSE.
#endif

   END FUNCTION nnp_gpu_available

! **************************************************************************************************
!> \brief Claim the offload device for this NNP environment and push the
!>        requested compute precision to the device side.
!>
!>        The device itself is chosen once per process at startup (see
!>        f77_interface), which is also the id DBCSR is given; this routine only
!>        activates it, exactly as every other CP2K offload consumer does.
!>        Re-selecting here would leave the two out of step for the rest of the
!>        run.
!> \param gpu_precision resolved GPU_PRECISION keyword value
!> \param accepted .FALSE. when another NNP environment already holds the
!>                 device, in which case nothing was initialized and the caller
!>                 must fall back to the CPU path
! **************************************************************************************************
   SUBROUTINE nnp_gpu_init_device(gpu_precision, accepted)
      INTEGER, INTENT(IN)                                :: gpu_precision
      LOGICAL, INTENT(OUT)                               :: accepted

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gpu_init_device'
      INTEGER                                            :: ierr

      CALL timeset(routineN, handle)

      CALL offload_activate_chosen_device()

      ! The device state is one per process, and its once-per-run uploads
      ! (network weights, pair map, scaling and symmetry-function tables) are
      ! keyed by shape rather than by environment. Two environments live at the
      ! same time and sharing an architecture would therefore evaluate the first
      ! one's weights against the second one's descriptor tables, with no
      ! diagnostic. Decline the device for the second one instead and let it run
      ! on the CPU. Sequential environments are fine: releasing one resets the
      ! upload flags.
      IF (nnp_gpu_register_env_c() > 1) THEN
         CALL nnp_gpu_unregister_env_c()
         accepted = .FALSE.
         CALL timestop(handle)
         RETURN
      END IF
      accepted = .TRUE.

      ierr = nnp_gpu_set_precision_c(MERGE(1, 0, gpu_precision == nnp_gpu_prec_mixed))
      IF (ierr /= 0) &
         CPABORT("NNP GPU init failed: cannot push the precision selector to the device.")

      CALL timestop(handle)
#else
      MARK_USED(gpu_precision)
      accepted = .FALSE.
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_init_device

! **************************************************************************************************
!> \brief Test whether the model fits the compiled-in device-table capacities:
!>        NNP_MAX_SE (group, element) slots per symmetry-function kind, and
!>        NNP_SPLINE_CACHE_MAX run-constant parameter tables, of which each
!>        radial group holds two (spline values and derivatives) and each
!>        angular group four (the packed parameter tables). Called before the
!>        device is claimed, so an oversized model is turned away before any
!>        upload rather than in the first force evaluation.
!> \param nnp the NNP environment with its radial and angular groups populated
!> \param reason on a .FALSE. result, a phrase naming the exceeded capacity
!> \return .TRUE. when the model fits the device tables
! **************************************************************************************************
   FUNCTION nnp_gpu_model_capacity_ok(nnp, reason) RESULT(fits)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason
      LOGICAL                                            :: fits

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      INTEGER                                            :: cache_needed, e, mag, mrg, se_cap

      fits = .TRUE.
      reason = ""

      CALL nnp_gpu_max_symfgrp(nnp, mrg, mag)
      se_cap = nnp_gpu_max_se_c()
      IF (mrg*nnp%n_ele > se_cap .OR. mag*nnp%n_ele > se_cap) THEN
         fits = .FALSE.
         WRITE (reason, '(A,I0,A,I0)') "the model needs ", &
            MAX(mrg*nnp%n_ele, mag*nnp%n_ele), &
            " (group, element) device slots per kind and this build holds "// &
            "NNP_MAX_SE = ", se_cap
         RETURN
      END IF

      cache_needed = 0
      DO e = 1, nnp%n_ele
         cache_needed = cache_needed + 2*nnp%rad(e)%n_symfgrp + 4*nnp%ang(e)%n_symfgrp
      END DO
      IF (cache_needed > nnp_gpu_spline_cache_max_c()) THEN
         fits = .FALSE.
         WRITE (reason, '(A,I0,A,I0)') "the model needs ", cache_needed, &
            " device parameter-cache entries and this build holds "// &
            "NNP_SPLINE_CACHE_MAX = ", nnp_gpu_spline_cache_max_c()
      END IF
#else
      MARK_USED(nnp)
      fits = .TRUE.
      reason = ""
#endif

   END FUNCTION nnp_gpu_model_capacity_ok

! **************************************************************************************************
!> \brief Largest per-element radial and angular symmetry-function group
!>        counts, floored to one. Every device-side (group, element) buffer
!>        shape and flattened slot index derives from these two strides, so
!>        every caller must take them from this one definition. The floor keeps
!>        the buffer shapes and kernel launch grids nonzero for a model with no
!>        groups of one kind; the padding slot does no work because the
!>        per-element group counts still drive the iteration.
!> \param nnp the NNP environment with its radial and angular groups populated
!> \param max_rad largest radial group count over the elements, at least one
!> \param max_ang largest angular group count over the elements, at least one
! **************************************************************************************************
   SUBROUTINE nnp_gpu_max_symfgrp(nnp, max_rad, max_ang)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp
      INTEGER, INTENT(OUT)                               :: max_rad, max_ang

      INTEGER                                            :: e

      max_rad = 1
      max_ang = 1
      DO e = 1, nnp%n_ele
         max_rad = MAX(max_rad, nnp%rad(e)%n_symfgrp)
         max_ang = MAX(max_ang, nnp%ang(e)%n_symfgrp)
      END DO

   END SUBROUTINE nnp_gpu_max_symfgrp

! **************************************************************************************************
!> \brief Build the device-resident neighbour lists for one MD step: upload the
!>        positions, image table and cell-list chain, pack and upload the
!>        species-pair routing table, then run the on-device cell-list walk over
!>        this rank's atoms and read back the per-(atom-slot, group, element)
!>        neighbour counts. On a per-list overflow the neighbour-list widths are
!>        doubled and the walk is retried.
!> \param nnp the NNP environment (its cell-list cache and pair map must already
!>            be prepared for this step)
!> \param istart first rank-local atom (1-based, global index)
!> \param mecalc number of rank-local atoms
! **************************************************************************************************
   SUBROUTINE nnp_gpu_run_walk(nnp, istart, mecalc)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: istart, mecalc

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_gpu_run_walk'
      INTEGER :: attempts, e_w, ec, gg, ie, ii, je, mag, mpg, mrg, n_atoms_in_e_max, ne, &
                 pflat, rc_walk, slot_w, ss
      INTEGER(KIND=c_int), DIMENSION(3)                  :: epbc_buf, nbin_buf, perd_buf, span_buf
      REAL(KIND=dp), DIMENSION(3)                        :: lower_buf, width_buf
      INTEGER(KIND=c_int), ALLOCATABLE, DIMENSION(:)     :: atoms_global_per_e_buf, &
                                                            ele_ind_global_buf, n_atoms_in_e_buf, &
                                                            pm_ang1_g, pm_ang2_g, pm_n_ang1, &
                                                            pm_n_ang2, pm_n_rad, pm_rad_g
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pm_ang_cut, pm_max_cut, pm_rad_cut

      CALL timeset(routineN, handle)

      IF (.NOT. ALLOCATED(nnp%gpu_host_state)) ALLOCATE (nnp%gpu_host_state)

      ne = nnp%n_ele
      CALL nnp_gpu_max_symfgrp(nnp, mrg, mag)

      ! Positions every step; the image table only when the translations have
      ! been regenerated since the copy the device is holding. The comparison is
      ! against the cache's rebuild counter rather than a per-step flag: this
      ! routine is reached only on the steps that take the GPU path, while the
      ! cache is rebuilt on any force evaluation, energy-only ones included.
      ec = nnp_gpu_upload_positions_step_c(INT(nnp%num_atoms, c_int), &
                                           nnp%cell_list_cache%coord_primary)
      IF (ec /= 0) CPABORT("NNP GPU: position upload failed.")

      IF (nnp%cell_list_cache%n_images > 0) THEN
         IF (nnp%gpu_host_state%image_table_generation /= &
             nnp%cell_list_cache%rebuild_count) THEN
            ec = nnp_gpu_upload_image_table_c(INT(nnp%cell_list_cache%n_images, c_int), &
                                              nnp%cell_list_cache%image_translation)
            IF (ec /= 0) CPABORT("NNP GPU: image-table upload failed.")
            nnp%gpu_host_state%image_table_generation = nnp%cell_list_cache%rebuild_count
         END IF
         IF (nnp%cell_list_cache%n_cells > 0) THEN
            ec = nnp_gpu_upload_cell_list_c(INT(nnp%cell_list_cache%n_images, c_int), &
                                            INT(nnp%cell_list_cache%n_cells, c_int), &
                                            INT(nnp%cell_list_cache%image_atom, c_int), &
                                            INT(nnp%cell_list_cache%image_shift, c_int), &
                                            INT(nnp%cell_list_cache%head, c_int), &
                                            INT(nnp%cell_list_cache%next, c_int))
            IF (ec /= 0) CPABORT("NNP GPU: cell-list upload failed.")
         END IF
      END IF

      ! Pack and upload the species-pair routing table. The launcher performs
      ! the host-to-device copy once and returns immediately thereafter, so it
      ! is safe to call every step. The walk applies the true symmetry-function
      ! cutoffs (it runs every step), so no Verlet-skin widening is needed.
      mpg = 1
      DO ie = 1, ne
         DO je = 1, ne
            mpg = MAX(mpg, nnp%neighbor_interface_state%pair_map(ie, je)%n_rad, &
                      nnp%neighbor_interface_state%pair_map(ie, je)%n_ang1, &
                      nnp%neighbor_interface_state%pair_map(ie, je)%n_ang2)
         END DO
      END DO
      ALLOCATE (pm_n_rad(ne*ne), pm_n_ang1(ne*ne), pm_n_ang2(ne*ne), pm_max_cut(ne*ne))
      ALLOCATE (pm_rad_g(mpg*ne*ne), pm_ang1_g(mpg*ne*ne), pm_ang2_g(mpg*ne*ne))
      ALLOCATE (pm_rad_cut(mrg*ne), pm_ang_cut(mag*ne))
      pm_rad_g(:) = 0
      pm_ang1_g(:) = 0
      pm_ang2_g(:) = 0
      pm_rad_cut(:) = 0.0_dp
      pm_ang_cut(:) = 0.0_dp
      DO je = 1, ne
         DO ie = 1, ne
            pflat = (ie - 1) + ne*(je - 1)
            ASSOCIATE (pmap => nnp%neighbor_interface_state%pair_map(ie, je))
               pm_n_rad(pflat + 1) = INT(pmap%n_rad, c_int)
               pm_n_ang1(pflat + 1) = INT(pmap%n_ang1, c_int)
               pm_n_ang2(pflat + 1) = INT(pmap%n_ang2, c_int)
               pm_max_cut(pflat + 1) = pmap%max_relevant_cutoff
               DO gg = 1, pmap%n_rad
                  pm_rad_g((gg - 1) + mpg*pflat + 1) = INT(pmap%rad_groups(gg), c_int)
               END DO
               DO gg = 1, pmap%n_ang1
                  pm_ang1_g((gg - 1) + mpg*pflat + 1) = INT(pmap%ang1_groups(gg), c_int)
               END DO
               DO gg = 1, pmap%n_ang2
                  pm_ang2_g((gg - 1) + mpg*pflat + 1) = INT(pmap%ang2_groups(gg), c_int)
               END DO
            END ASSOCIATE
         END DO
      END DO
      DO ie = 1, ne
         DO ss = 1, nnp%rad(ie)%n_symfgrp
            pm_rad_cut((ss - 1) + mrg*(ie - 1) + 1) = nnp%rad(ie)%symfgrp(ss)%cutoff
         END DO
         DO ss = 1, nnp%ang(ie)%n_symfgrp
            pm_ang_cut((ss - 1) + mag*(ie - 1) + 1) = nnp%ang(ie)%symfgrp(ss)%cutoff
         END DO
      END DO
      ec = nnp_gpu_upload_pair_map_c(INT(ne, c_int), INT(mpg, c_int), INT(mrg, c_int), &
                                     INT(mag, c_int), pm_n_rad, pm_n_ang1, pm_n_ang2, &
                                     pm_max_cut, pm_rad_g, pm_ang1_g, pm_ang2_g, &
                                     pm_rad_cut, pm_ang_cut)
      IF (ec /= 0) CPABORT("NNP GPU: pair-map upload failed.")

      ! Count rank-local atoms per element and build the (slot, element) to
      ! global-atom map in the same order the descriptor kernels expect.
      ALLOCATE (n_atoms_in_e_buf(ne))
      n_atoms_in_e_buf(:) = 0
      DO ii = istart, istart + mecalc - 1
         e_w = nnp%ele_ind(ii)
         n_atoms_in_e_buf(e_w) = n_atoms_in_e_buf(e_w) + 1
      END DO
      n_atoms_in_e_max = MAXVAL(n_atoms_in_e_buf)

      ! Publish the per-element counts before any early return: a rank whose
      ! atom window is empty (more ranks than atoms) still has to leave the host
      ! state in a well-defined state, because the caller checks it rather than
      ! recomputing the window.
      nnp%gpu_host_state%n_atoms_in_e_max = n_atoms_in_e_max
      IF (ALLOCATED(nnp%gpu_host_state%n_atoms_in_e)) &
         DEALLOCATE (nnp%gpu_host_state%n_atoms_in_e)
      ALLOCATE (nnp%gpu_host_state%n_atoms_in_e(ne))
      nnp%gpu_host_state%n_atoms_in_e(:) = n_atoms_in_e_buf(:)

      IF (n_atoms_in_e_max <= 0) THEN
         DEALLOCATE (n_atoms_in_e_buf, pm_n_rad, pm_n_ang1, pm_n_ang2, pm_max_cut, &
                     pm_rad_g, pm_ang1_g, pm_ang2_g, pm_rad_cut, pm_ang_cut)
         CALL timestop(handle)
         RETURN
      END IF
      ALLOCATE (atoms_global_per_e_buf(n_atoms_in_e_max*ne))
      atoms_global_per_e_buf(:) = 0
      n_atoms_in_e_buf(:) = 0
      DO ii = istart, istart + mecalc - 1
         e_w = nnp%ele_ind(ii)
         n_atoms_in_e_buf(e_w) = n_atoms_in_e_buf(e_w) + 1
         slot_w = n_atoms_in_e_buf(e_w)
         atoms_global_per_e_buf((slot_w - 1) + n_atoms_in_e_max*(e_w - 1) + 1) = INT(ii, c_int)
      END DO
      ALLOCATE (ele_ind_global_buf(nnp%num_atoms))
      DO ii = 1, nnp%num_atoms
         ele_ind_global_buf(ii) = INT(nnp%ele_ind(ii), c_int)
      END DO

      nbin_buf = INT(nnp%cell_list_cache%nbin, c_int)
      span_buf = INT(nnp%cell_list_cache%bin_span, c_int)
      perd_buf = INT(nnp%cell_list_cache%perd, c_int)
      epbc_buf = INT(nnp%cell_list_cache%exact_pbc_copies, c_int)
      lower_buf = nnp%cell_list_cache%lower
      width_buf = nnp%cell_list_cache%bin_width

      IF (.NOT. nnp%gpu_host_state%walk_initialised) THEN
         nnp%gpu_host_state%max_n_neigh_rad = gpu_walk_init_max_n_neigh
         nnp%gpu_host_state%max_n_neigh_ang1 = gpu_walk_init_max_n_neigh
         nnp%gpu_host_state%max_n_neigh_ang2 = gpu_walk_init_max_n_neigh
         nnp%gpu_host_state%walk_initialised = .TRUE.
      END IF

      attempts = 0
      DO
         attempts = attempts + 1
         rc_walk = nnp_gpu_cell_list_walk_c(INT(nnp%num_atoms, c_int), INT(ne, c_int), &
                                            INT(nnp%cell_list_cache%n_images, c_int), &
                                            INT(nnp%cell_list_cache%n_cells, c_int), &
                                            INT(n_atoms_in_e_max, c_int), &
                                            INT(nnp%gpu_host_state%max_n_neigh_rad, c_int), &
                                            INT(nnp%gpu_host_state%max_n_neigh_ang1, c_int), &
                                            INT(nnp%gpu_host_state%max_n_neigh_ang2, c_int), &
                                            INT(mrg, c_int), INT(mag, c_int), nbin_buf, &
                                            lower_buf, width_buf, span_buf, perd_buf, epbc_buf, &
                                            nnp%cell_list_cache%coord_scaled, ele_ind_global_buf, &
                                            atoms_global_per_e_buf, n_atoms_in_e_buf)
         IF (rc_walk == 0) EXIT
         IF (rc_walk == -2 .AND. attempts < 4) THEN
            nnp%gpu_host_state%max_n_neigh_rad = nnp%gpu_host_state%max_n_neigh_rad*2
            nnp%gpu_host_state%max_n_neigh_ang1 = nnp%gpu_host_state%max_n_neigh_ang1*2
            nnp%gpu_host_state%max_n_neigh_ang2 = nnp%gpu_host_state%max_n_neigh_ang2*2
            CYCLE
         END IF
         IF (rc_walk == -2) THEN
            CPABORT("NNP GPU: a per-atom neighbour list still overflows after three width doublings.")
         END IF
         CPABORT("NNP GPU: cell-list walk failed.")
      END DO

      ! Read back the per-(slot, group, element) neighbour counts, used by the
      ! descriptor kernels to size their per-atom neighbour buffers.
      IF (ALLOCATED(nnp%gpu_host_state%n_rad_pa)) DEALLOCATE (nnp%gpu_host_state%n_rad_pa)
      IF (ALLOCATED(nnp%gpu_host_state%n_ang1_pa)) DEALLOCATE (nnp%gpu_host_state%n_ang1_pa)
      IF (ALLOCATED(nnp%gpu_host_state%n_ang2_pa)) DEALLOCATE (nnp%gpu_host_state%n_ang2_pa)
      ALLOCATE (nnp%gpu_host_state%n_rad_pa(n_atoms_in_e_max, mrg, ne))
      ALLOCATE (nnp%gpu_host_state%n_ang1_pa(n_atoms_in_e_max, mag, ne))
      ALLOCATE (nnp%gpu_host_state%n_ang2_pa(n_atoms_in_e_max, mag, ne))
      nnp%gpu_host_state%n_rad_pa(:, :, :) = 0
      nnp%gpu_host_state%n_ang1_pa(:, :, :) = 0
      nnp%gpu_host_state%n_ang2_pa(:, :, :) = 0
      ec = nnp_gpu_walk_counts_d2h_c(INT(n_atoms_in_e_max, c_int), INT(mrg, c_int), &
                                     INT(mag, c_int), INT(ne, c_int), &
                                     nnp%gpu_host_state%n_rad_pa, nnp%gpu_host_state%n_ang1_pa, &
                                     nnp%gpu_host_state%n_ang2_pa)
      IF (ec /= 0) CPABORT("NNP GPU: walk-count read-back failed.")

      DEALLOCATE (n_atoms_in_e_buf, atoms_global_per_e_buf, ele_ind_global_buf, &
                  pm_n_rad, pm_n_ang1, pm_n_ang2, pm_max_cut, pm_rad_g, pm_ang1_g, pm_ang2_g, &
                  pm_rad_cut, pm_ang_cut)

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(istart)
      MARK_USED(mecalc)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_run_walk

! **************************************************************************************************
!> \brief Upload the model-static descriptor tables to the device once at NNP
!>        initialization: the per-element symmetry-function scaling parameters
!>        (average, minimum, maximum, sigma) and the per-(group, element)
!>        symmetry-function to descriptor-row maps. Both are immutable after the
!>        input parse, so the device copies are valid for the whole environment
!>        lifetime and the per-step descriptor kernels read them without any
!>        further host-to-device traffic.
!> \param nnp the NNP environment, with its radial and angular groups populated
! **************************************************************************************************
   SUBROUTINE nnp_gpu_upload_descriptor_statics(nnp)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gpu_upload_descriptor_statics'
      INTEGER :: e, ec, m, max_n_anggrp, max_n_radgrp, n_ang_se_total, n_input_max, &
                 n_rad_se_total, n_symf_g, s, se_flat, se_off, sf, total_ang_symf, total_rad_symf
      INTEGER(KIND=c_int), ALLOCATABLE, DIMENSION(:) :: ang_offsets, ang_packed, n_in_buf, &
                                                        rad_offsets, rad_packed
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: loc_av_buf, loc_max_buf, loc_min_buf, &
                                                           sigma_buf

      CALL timeset(routineN, handle)

      ASSOCIATE (n_ele => nnp%n_ele)

         ! Per-element scaling tables. The device layout is one slab of
         ! n_input_max doubles per element, with the element's radial groups in
         ! the leading slots and its angular groups next, matching the layout
         ! the scaling kernel indexes.
         n_input_max = 0
         DO e = 1, n_ele
            n_input_max = MAX(n_input_max, nnp%arc(e)%n_nodes(1))
         END DO

         ALLOCATE (loc_av_buf(n_ele*n_input_max))
         ALLOCATE (loc_min_buf(n_ele*n_input_max))
         ALLOCATE (loc_max_buf(n_ele*n_input_max))
         ALLOCATE (sigma_buf(n_ele*n_input_max))
         ALLOCATE (n_in_buf(n_ele))
         loc_av_buf = 0.0_dp
         loc_min_buf = 0.0_dp
         loc_max_buf = 0.0_dp
         sigma_buf = 0.0_dp

         DO e = 1, n_ele
            n_in_buf(e) = INT(nnp%arc(e)%n_nodes(1), c_int)
            DO m = 1, nnp%n_rad(e)
               loc_av_buf((e - 1)*n_input_max + m) = nnp%rad(e)%loc_av(m)
               loc_min_buf((e - 1)*n_input_max + m) = nnp%rad(e)%loc_min(m)
               loc_max_buf((e - 1)*n_input_max + m) = nnp%rad(e)%loc_max(m)
               sigma_buf((e - 1)*n_input_max + m) = nnp%rad(e)%sigma(m)
            END DO
            DO m = 1, nnp%n_ang(e)
               loc_av_buf((e - 1)*n_input_max + nnp%n_rad(e) + m) = nnp%ang(e)%loc_av(m)
               loc_min_buf((e - 1)*n_input_max + nnp%n_rad(e) + m) = nnp%ang(e)%loc_min(m)
               loc_max_buf((e - 1)*n_input_max + nnp%n_rad(e) + m) = nnp%ang(e)%loc_max(m)
               sigma_buf((e - 1)*n_input_max + nnp%n_rad(e) + m) = nnp%ang(e)%sigma(m)
            END DO
         END DO

         ec = nnp_gpu_upload_scaling_tables_c(INT(n_ele, c_int), INT(n_input_max, c_int), &
                                              loc_av_buf, loc_min_buf, loc_max_buf, sigma_buf, &
                                              n_in_buf)
         IF (ec /= 0) CPABORT("NNP GPU: descriptor scaling-table upload failed.")

         DEALLOCATE (loc_av_buf, loc_min_buf, loc_max_buf, sigma_buf, n_in_buf)

         ! Per-(group, element) symmetry-function maps, packed radial then
         ! angular with prefix-sum offsets keyed se_flat = (s-1) + max_grp*(e-1)
         ! (element outermost). Padding slots for elements with fewer groups get
         ! zero-length entries; the final slot carries the sentinel total.
         CALL nnp_gpu_max_symfgrp(nnp, max_n_radgrp, max_n_anggrp)

         n_rad_se_total = max_n_radgrp*n_ele
         n_ang_se_total = max_n_anggrp*n_ele

         total_rad_symf = 0
         DO e = 1, n_ele
            DO s = 1, nnp%rad(e)%n_symfgrp
               total_rad_symf = total_rad_symf + nnp%rad(e)%symfgrp(s)%n_symf
            END DO
         END DO
         total_ang_symf = 0
         DO e = 1, n_ele
            DO s = 1, nnp%ang(e)%n_symfgrp
               total_ang_symf = total_ang_symf + nnp%ang(e)%symfgrp(s)%n_symf
            END DO
         END DO

         ALLOCATE (rad_offsets(MAX(1, n_rad_se_total + 1)))
         ALLOCATE (ang_offsets(MAX(1, n_ang_se_total + 1)))
         ALLOCATE (rad_packed(MAX(1, total_rad_symf)))
         ALLOCATE (ang_packed(MAX(1, total_ang_symf)))
         rad_offsets = 0_c_int
         ang_offsets = 0_c_int
         rad_packed = 0_c_int
         ang_packed = 0_c_int

         se_off = 0
         DO e = 1, n_ele
            DO s = 1, max_n_radgrp
               se_flat = (s - 1) + max_n_radgrp*(e - 1)
               rad_offsets(se_flat + 1) = INT(se_off, c_int)
               IF (s <= nnp%rad(e)%n_symfgrp) THEN
                  n_symf_g = nnp%rad(e)%symfgrp(s)%n_symf
                  DO sf = 1, n_symf_g
                     rad_packed(se_off + sf) = INT(nnp%rad(e)%symfgrp(s)%symf(sf), c_int)
                  END DO
                  se_off = se_off + n_symf_g
               END IF
            END DO
         END DO
         IF (n_rad_se_total > 0) rad_offsets(n_rad_se_total + 1) = INT(se_off, c_int)

         se_off = 0
         DO e = 1, n_ele
            DO s = 1, max_n_anggrp
               se_flat = (s - 1) + max_n_anggrp*(e - 1)
               ang_offsets(se_flat + 1) = INT(se_off, c_int)
               IF (s <= nnp%ang(e)%n_symfgrp) THEN
                  n_symf_g = nnp%ang(e)%symfgrp(s)%n_symf
                  DO sf = 1, n_symf_g
                     ang_packed(se_off + sf) = INT(nnp%ang(e)%symfgrp(s)%symf(sf), c_int)
                  END DO
                  se_off = se_off + n_symf_g
               END IF
            END DO
         END DO
         IF (n_ang_se_total > 0) ang_offsets(n_ang_se_total + 1) = INT(se_off, c_int)

         ec = nnp_gpu_upload_se_static_tables_c(INT(n_ele, c_int), INT(max_n_radgrp, c_int), &
                                                INT(max_n_anggrp, c_int), &
                                                INT(total_rad_symf, c_int), &
                                                INT(total_ang_symf, c_int), &
                                                rad_packed, rad_offsets, ang_packed, ang_offsets)
         IF (ec /= 0) CPABORT("NNP GPU: descriptor symf-map upload failed.")

         DEALLOCATE (rad_packed, rad_offsets, ang_packed, ang_offsets)

      END ASSOCIATE

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_upload_descriptor_statics

! **************************************************************************************************
!> \brief Build the device-resident descriptor for one MD step: zero G_dev,
!>        dispatch the radial then the angular symmetry-function kernels for
!>        every populated (group, element) slot, then scale G_dev into the
!>        network input. The kernels read the neighbour lists, atom maps and
!>        symmetry-function maps from device buffers, so this routine only
!>        supplies per-group sizes, model parameters and index scalars. It must
!>        be called after nnp_gpu_run_walk has populated the host walk-count
!>        state for this step.
!> \param nnp the NNP environment
!> \param istart first rank-local atom (1-based, global index)
!> \param mecalc number of rank-local atoms
!> \param n_input_max descriptor-row stride of G_dev (max over elements of the
!>                    network input-layer width)
! **************************************************************************************************
   SUBROUTINE nnp_gpu_build_descriptor(nnp, istart, mecalc, n_input_max)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: istart, mecalc, n_input_max

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gpu_build_descriptor'
      INTEGER :: cut_int, e, ec, i, i_local, max_n_ang1, max_n_ang2, max_n_anggrp, max_n_kk, &
                 max_n_neigh, max_n_radgrp, n_atoms_e, s, scaling_mode
      LOGICAL :: homo, splines_ready
      REAL(KIND=dp)                                      :: cutoff_s, cutoff_sqr
      INTEGER(KIND=c_int), ALLOCATABLE, DIMENSION(:)     :: ele_ind_buf, izeta_buf
      CHARACTER(KIND=c_char), ALLOCATABLE, DIMENSION(:)  :: use_int_zeta_buf

      CALL timeset(routineN, handle)

      ! Group-count strides used to flatten the (group, element) slot index the
      ! device buffers are keyed by.
      CALL nnp_gpu_max_symfgrp(nnp, max_n_radgrp, max_n_anggrp)

      ! Zero G_dev for this step's rank-local atoms.
      ec = nnp_gpu_begin_descriptor_step_c(INT(mecalc, c_int), INT(n_input_max, c_int))
      IF (ec /= 0) CPABORT("NNP GPU: descriptor-step begin failed.")

      ! The radial kernels sample the Hermite splines built from the model's
      ! radial symmetry-function parameters; build them once if a fresh
      ! environment has not done so yet.
      splines_ready = .TRUE.
      DO e = 1, nnp%n_ele
         DO s = 1, nnp%rad(e)%n_symfgrp
            IF (.NOT. nnp%rad(e)%symfgrp(s)%spline_built) splines_ready = .FALSE.
         END DO
      END DO
      IF (.NOT. splines_ready) CALL nnp_build_radial_splines(nnp)

      ASSOCIATE (hstate => nnp%gpu_host_state)

         ! Radial symmetry functions into G_dev, one dispatch per populated
         ! (group, element) slot.
         DO e = 1, nnp%n_ele
            n_atoms_e = hstate%n_atoms_in_e(e)
            IF (n_atoms_e <= 0) CYCLE
            DO s = 1, nnp%rad(e)%n_symfgrp
               ASSOCIATE (symfgrp => nnp%rad(e)%symfgrp(s))
                  max_n_neigh = MAXVAL(hstate%n_rad_pa(1:n_atoms_e, s, e))
                  IF (max_n_neigh <= 0) CYCLE
                  ec = nnp_gpu_radial_sf_to_G_dev_c( &
                       INT(n_atoms_e, c_int), INT(max_n_neigh, c_int), &
                       INT(symfgrp%n_symf, c_int), INT(symfgrp%spline_n, c_int), &
                       REAL(symfgrp%spline_dx, c_double), &
                       REAL(symfgrp%spline_dx_inv, c_double), &
                       REAL(symfgrp%spline_x_max, c_double), &
                       symfgrp%spline_y, symfgrp%spline_dy, &
                       INT(s - 1, c_int), INT(e - 1, c_int), INT(max_n_radgrp, c_int), &
                       INT(hstate%n_atoms_in_e_max, c_int), &
                       INT(hstate%max_n_neigh_rad, c_int), INT(istart, c_int), &
                       INT(n_input_max, c_int), INT((s - 1) + max_n_radgrp*(e - 1), c_int))
                  IF (ec /= 0) CPABORT("NNP GPU: radial descriptor kernel failed.")
               END ASSOCIATE
            END DO
         END DO

         ! Angular symmetry functions into G_dev. max_n_ang1 / max_n_ang2 are
         ! the per-slot high-water neighbour counts from the walk; max_n_ang2 is
         ! zero for homonuclear groups (their second neighbour list is unused)
         ! and the kernel's neighbour-pair width max_n_kk is the ang1 list width
         ! for homonuclear groups and the ang2 list width otherwise.
         cut_int = nnp%cut_type - 1
         DO e = 1, nnp%n_ele
            n_atoms_e = hstate%n_atoms_in_e(e)
            IF (n_atoms_e <= 0) CYCLE
            DO s = 1, nnp%ang(e)%n_symfgrp
               ASSOCIATE (symfgrp => nnp%ang(e)%symfgrp(s))
                  homo = (symfgrp%ele(1) == symfgrp%ele(2))
                  cutoff_s = symfgrp%cutoff
                  cutoff_sqr = cutoff_s*cutoff_s

                  max_n_ang1 = MAXVAL(hstate%n_ang1_pa(1:n_atoms_e, s, e))
                  IF (max_n_ang1 <= 0) CYCLE
                  IF (homo) THEN
                     max_n_ang2 = 0
                     max_n_kk = MAX(1, max_n_ang1)
                  ELSE
                     max_n_ang2 = MAXVAL(hstate%n_ang2_pa(1:n_atoms_e, s, e))
                     ! No atom of e has an ele(2) neighbour in range: the pair
                     ! list is empty and the group's contribution is exactly
                     ! zero (G_dev is zeroed each step), so skip it, mirroring
                     ! the max_n_ang1 guard above. Happens e.g. for ion-ion
                     ! angular groups in dilute electrolytes whenever no ion
                     ! pair is within the cutoff. The force-scatter loop
                     ! below applies the SAME guard: a skipped group leaves
                     ! its per-(s,e) device buffers unpopulated, and
                     ! launching the scatter on them would read stale data.
                     IF (max_n_ang2 <= 0) CYCLE
                     max_n_kk = MAX(1, max_n_ang2)
                  END IF

                  ALLOCATE (izeta_buf(symfgrp%n_symf))
                  ALLOCATE (use_int_zeta_buf(symfgrp%n_symf))
                  DO i = 1, symfgrp%n_symf
                     izeta_buf(i) = INT(symfgrp%pack_izeta(i), c_int)
                     use_int_zeta_buf(i) = MERGE(ACHAR(1), ACHAR(0), symfgrp%pack_use_int_zeta(i))
                  END DO

                  ec = nnp_gpu_angular_sf_to_G_dev_c( &
                       INT(MERGE(1, 0, homo), c_int), INT(cut_int, c_int), &
                       INT(MAX(1, max_n_ang1), c_int), INT(max_n_ang2, c_int), &
                       INT(max_n_kk, c_int), INT(symfgrp%n_symf, c_int), &
                       INT(n_atoms_e, c_int), REAL(cutoff_s, c_double), &
                       REAL(cutoff_sqr, c_double), &
                       symfgrp%pack_lam, symfgrp%pack_zeta, symfgrp%pack_eta, &
                       symfgrp%pack_prefzeta, izeta_buf, use_int_zeta_buf, &
                       INT(s - 1, c_int), INT(e - 1, c_int), INT(max_n_anggrp, c_int), &
                       INT(hstate%n_atoms_in_e_max, c_int), &
                       INT(hstate%max_n_neigh_ang1, c_int), &
                       INT(hstate%max_n_neigh_ang2, c_int), INT(istart, c_int), &
                       INT(nnp%n_rad(e), c_int), INT(n_input_max, c_int), &
                       INT((s - 1) + max_n_anggrp*(e - 1), c_int))
                  IF (ec /= 0) CPABORT("NNP GPU: angular descriptor kernel failed.")

                  DEALLOCATE (izeta_buf, use_int_zeta_buf)
               END ASSOCIATE
            END DO
         END DO

      END ASSOCIATE

      ! Scale the aggregated descriptor into the network input. The scaling mode
      ! selects which of centering, scaling and sigma normalization the model
      ! requested; the per-atom 0-based element index routes each atom to its
      ! element's scaling parameters.
      IF (nnp%scale_sigma_acsf) THEN
         scaling_mode = 4
      ELSE IF (nnp%center_acsf .AND. nnp%scale_acsf) THEN
         scaling_mode = 3
      ELSE IF (nnp%center_acsf) THEN
         scaling_mode = 2
      ELSE IF (nnp%scale_acsf) THEN
         scaling_mode = 1
      ELSE
         scaling_mode = 0
      END IF

      ALLOCATE (ele_ind_buf(MAX(1, mecalc)))
      DO i = istart, istart + mecalc - 1
         i_local = i - istart + 1
         ele_ind_buf(i_local) = INT(nnp%ele_ind(i) - 1, c_int)
      END DO

      ec = nnp_gpu_scale_G_dev_c(INT(mecalc, c_int), INT(n_input_max, c_int), &
                                 INT(nnp%n_ele, c_int), INT(scaling_mode, c_int), &
                                 REAL(nnp%scmin, c_double), REAL(nnp%scmax, c_double), &
                                 ele_ind_buf)
      IF (ec /= 0) CPABORT("NNP GPU: descriptor scaling failed.")

      DEALLOCATE (ele_ind_buf)

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(istart)
      MARK_USED(mecalc)
      MARK_USED(n_input_max)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_build_descriptor

! **************************************************************************************************
!> \brief Upload the committee network to the device once at NNP initialization:
!>        the per-element layer topology, the shared activation-function codes,
!>        the per-element atomic reference energies and the packed per-committee
!>        weights and biases. These are immutable after the weights are read, so
!>        the device copies serve every MD step and the per-step network launcher
!>        needs no further weight traffic. The weights and biases are packed into
!>        zero-padded max_n_nodes-wide slabs so a single stride addresses every
!>        (element, layer, committee) block regardless of that block's true size.
!> \param nnp the NNP environment, with its architecture and committee weights
!>            already read
! **************************************************************************************************
   SUBROUTINE nnp_gpu_upload_network(nnp)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gpu_upload_network'
      INTEGER :: c, e, ec, i, idx_b, idx_W, j, l, max_n_nodes, n_com, n_ele, n_in, n_layer, n_out
      INTEGER(KIND=c_int), ALLOCATABLE, DIMENSION(:)     :: act_fnct_buf, n_nodes_buf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: b_packed_buf, w_packed_buf

      CALL timeset(routineN, handle)

      n_com = nnp%n_committee
      n_layer = nnp%n_layer
      n_ele = nnp%n_ele

      ! Slab width for the packed weights and biases: the largest layer node
      ! count over all elements and layers. The kernel reads only the leading
      ! n_nodes(L-1) x n_nodes(L) corner of each slab.
      max_n_nodes = 0
      DO e = 1, n_ele
         DO l = 1, n_layer
            IF (nnp%arc(e)%n_nodes(l) > max_n_nodes) max_n_nodes = nnp%arc(e)%n_nodes(l)
         END DO
      END DO

      ! Layer node counts, flattened layer-fastest then element.
      ALLOCATE (n_nodes_buf(n_layer*n_ele))
      DO e = 1, n_ele
         DO l = 1, n_layer
            n_nodes_buf((l - 1) + n_layer*(e - 1) + 1) = INT(nnp%arc(e)%n_nodes(l), c_int)
         END DO
      END DO

      ! Activation-function codes for layers 2..n_layer (shared across elements).
      ALLOCATE (act_fnct_buf(n_layer - 1))
      DO l = 2, n_layer
         act_fnct_buf(l - 1) = INT(nnp%actfnct(l - 1), c_int)
      END DO

      ! Weights and biases packed into zero-padded max_n_nodes-wide slabs, one
      ! slab per (element, layer, committee). Input node i is fastest in the
      ! weight slab; the bias slab is one node column per (element, layer,
      ! committee). Layout matches nnp_gpu_nn_atomwise_with_G_dev_c and the
      ! legacy host packing in nnp_model.F.
      ALLOCATE (w_packed_buf(max_n_nodes*max_n_nodes*n_com*(n_layer - 1)*n_ele))
      ALLOCATE (b_packed_buf(max_n_nodes*n_com*(n_layer - 1)*n_ele))
      w_packed_buf = 0.0_dp
      b_packed_buf = 0.0_dp
      DO e = 1, n_ele
         DO l = 2, n_layer
            n_in = nnp%arc(e)%n_nodes(l - 1)
            n_out = nnp%arc(e)%n_nodes(l)
            DO c = 1, n_com
               DO j = 1, n_out
                  DO i = 1, n_in
                     idx_W = (i - 1) + max_n_nodes*(j - 1) + &
                             max_n_nodes*max_n_nodes*(c - 1) + &
                             max_n_nodes*max_n_nodes*n_com*(l - 2) + &
                             max_n_nodes*max_n_nodes*n_com*(n_layer - 1)*(e - 1) + 1
                     w_packed_buf(idx_W) = nnp%arc(e)%layer(l)%weights(i, j, c)
                  END DO
                  idx_b = (j - 1) + max_n_nodes*(c - 1) + &
                          max_n_nodes*n_com*(l - 2) + &
                          max_n_nodes*n_com*(n_layer - 1)*(e - 1) + 1
                  b_packed_buf(idx_b) = nnp%arc(e)%layer(l)%bweights(j, c)
               END DO
            END DO
         END DO
      END DO

      ec = nnp_gpu_upload_nn_packed_c(INT(n_com, c_int), INT(n_layer, c_int), &
                                      INT(max_n_nodes, c_int), INT(n_ele, c_int), &
                                      n_nodes_buf, act_fnct_buf, nnp%atom_energies, &
                                      w_packed_buf, b_packed_buf)
      IF (ec /= 0) CPABORT("NNP GPU: committee network upload failed.")

      DEALLOCATE (n_nodes_buf, act_fnct_buf, w_packed_buf, b_packed_buf)

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_upload_network

! **************************************************************************************************
!> \brief Evaluate the committee network for one MD step from the scaled
!>        device-resident descriptor. The network topology and weights were
!>        uploaded once by nnp_gpu_upload_network, so this launcher only supplies
!>        the per-atom element index and receives the per-atom per-committee
!>        energies back. The reference atomic energies are added inside the
!>        kernel, so the returned energy is the complete per-atom contribution
!>        and is copied straight into nnp%atomic_energy without any further
!>        offset. The descriptor gradient is left on the device for the force
!>        scatter; its host buffer is only allocated and passed to satisfy the
!>        launcher signature. Must be called after nnp_gpu_build_descriptor.
!> \param nnp the NNP environment; nnp%atomic_energy is filled for this window
!> \param istart first rank-local atom (1-based, global index)
!> \param mecalc number of rank-local atoms
!> \param n_input_max descriptor-row stride of the device descriptor (max over
!>                    elements of the network input-layer width)
! **************************************************************************************************
   SUBROUTINE nnp_gpu_predict(nnp, istart, mecalc, n_input_max)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: istart, mecalc, n_input_max

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_gpu_predict'
      INTEGER :: atom_local, c, e, ec, i, l, max_n_nodes, n_com, n_ele, n_layer
      INTEGER(KIND=c_int), ALLOCATABLE, DIMENSION(:)     :: ele_ind_buf
      REAL(KIND=c_double), ALLOCATABLE, DIMENSION(:)     :: h_atomic_energy

      CALL timeset(routineN, handle)

      n_com = nnp%n_committee
      n_layer = nnp%n_layer
      n_ele = nnp%n_ele

      IF (mecalc <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! Slab width the packed weights were uploaded with; the launcher reads the
      ! device weights with this stride.
      max_n_nodes = 0
      DO e = 1, n_ele
         DO l = 1, n_layer
            IF (nnp%arc(e)%n_nodes(l) > max_n_nodes) max_n_nodes = nnp%arc(e)%n_nodes(l)
         END DO
      END DO

      ! Per-atom 0-based element index over this rank's window [istart, ...].
      ALLOCATE (ele_ind_buf(mecalc))
      DO i = istart, istart + mecalc - 1
         atom_local = i - istart + 1
         ele_ind_buf(atom_local) = INT(nnp%ele_ind(i) - 1, c_int)
      END DO

      ! Network outputs in kernel layout: energy[c + n_com*atom_local]. The
      ! descriptor gradient stays on the device, where the force scatter reads
      ! it, so no host buffer is offered for it.
      ALLOCATE (h_atomic_energy(n_com*mecalc))
      h_atomic_energy = 0.0_c_double

      ec = nnp_gpu_nn_atomwise_with_G_dev_c(INT(mecalc, c_int), INT(n_com, c_int), &
                                            INT(n_layer, c_int), INT(max_n_nodes, c_int), &
                                            INT(n_input_max, c_int), &
                                            MERGE(1_c_int, 0_c_int, nnp%normnodes), &
                                            INT(n_ele, c_int), ele_ind_buf, &
                                            h_atomic_energy, c_null_ptr)
      IF (ec /= 0) CPABORT("NNP GPU: committee network evaluation failed.")

      ! Store the per-atom per-committee energy. The kernel already added the
      ! reference atomic energy atom_energies(ele), so this is the full per-atom
      ! contribution and is copied without any further offset (mirroring
      ! nnp_model.F nnp_predict_gradient_atomwise_gpu, which stores E_out_flat
      ! directly, and nnp_force.F, which copies atomic_energy_batch verbatim).
      DO atom_local = 1, mecalc
         i = istart + atom_local - 1
         DO c = 1, n_com
            nnp%atomic_energy(i, c) = h_atomic_energy((c - 1) + n_com*(atom_local - 1) + 1)
         END DO
      END DO

      DEALLOCATE (ele_ind_buf, h_atomic_energy)

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(istart)
      MARK_USED(mecalc)
      MARK_USED(n_input_max)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_predict

! **************************************************************************************************
!> \brief Run the whole per-step NNP force evaluation on the device: the
!>        cell-list walk, the descriptor build, the committee network and the
!>        force scatter, in that order. The walk builds the device neighbour
!>        lists, the descriptor build fills the scaled device descriptor, the
!>        network fills nnp%atomic_energy and leaves its descriptor gradient on
!>        the device, and the scatter multiplies that gradient by the descriptor
!>        derivatives and reads the total forces back into nnp%myforce. The
!>        scatter accumulator is seeded from nnp%myforce, which the caller has
!>        zeroed, so on return it holds the complete NNP forces for this window.
!> \param nnp the NNP environment; nnp%atomic_energy and nnp%myforce are filled
!> \param istart first rank-local atom (1-based, global index)
!> \param mecalc number of rank-local atoms
!> \param n_input_max descriptor-row stride of the device descriptor (max over
!>                    elements of the network input-layer width)
!> \param calc_stress .TRUE. when the caller wants the analytic virial, which
!>                    the force-scatter kernels then accumulate alongside the
!>                    forces
! **************************************************************************************************
   SUBROUTINE nnp_gpu_energy_force(nnp, istart, mecalc, n_input_max, calc_stress)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: istart, mecalc, n_input_max
      LOGICAL, INTENT(IN)                                :: calc_stress

      INTEGER                                            :: handle

#if (defined(__OFFLOAD) || defined(__OFFLOAD_CUDA) || defined(__OFFLOAD_HIP)) && !defined(__NO_OFFLOAD_NNP)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_gpu_energy_force'
      INTEGER :: accumulate_virial, e, ec, iva, ivb, ivc, max_n_ang1, max_n_ang2, max_n_anggrp, &
                 max_n_kk, max_n_neigh_per_se, max_n_radgrp, n_atoms_e, s, scaling_mode
      LOGICAL                                            :: homo
      REAL(KIND=c_double), ALLOCATABLE                   :: c_virial(:)

      CALL timeset(routineN, handle)

      ! A rank can legitimately own no atoms (more ranks than atoms), in which
      ! case it contributes nothing: the caller has already zeroed nnp%myforce
      ! and nnp%atomic_energy for this window. The device kernels all require a
      ! non-empty window, so return before any of them is launched.
      IF (mecalc <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! When a stress tensor is requested, the from-state scatter kernels also
      ! accumulate the analytical virial on device; the flag keeps the kernel
      ! virial branch off on energy/force-only steps.
      accumulate_virial = MERGE(1, 0, calc_stress)

      ! Device neighbour lists, then the scaled device descriptor, then the
      ! committee energy (which also leaves its descriptor gradient on device).
      CALL nnp_gpu_run_walk(nnp, istart, mecalc)
      CALL nnp_gpu_build_descriptor(nnp, istart, mecalc, n_input_max)
      CALL nnp_gpu_predict(nnp, istart, mecalc, n_input_max)

      ! Force scatter. The scaling mode selects which of centering, scaling and
      ! sigma normalization the model requested; it must match the mode the
      ! descriptor build used, so the kernels undo the same transform when they
      ! propagate the descriptor gradient back to the atoms.
      IF (nnp%scale_sigma_acsf) THEN
         scaling_mode = 4
      ELSE IF (nnp%center_acsf .AND. nnp%scale_acsf) THEN
         scaling_mode = 3
      ELSE IF (nnp%center_acsf) THEN
         scaling_mode = 2
      ELSE IF (nnp%scale_acsf) THEN
         scaling_mode = 1
      ELSE
         scaling_mode = 0
      END IF

      ! Group-count strides used to flatten the (group, element) slot index the
      ! device walk and scatter buffers are keyed by.
      CALL nnp_gpu_max_symfgrp(nnp, max_n_radgrp, max_n_anggrp)

      ! Seed the device force accumulator from the host force array (zeroed by
      ! the caller).
      ec = nnp_gpu_force_scatter_begin_c(INT(nnp%num_atoms, c_int), &
                                         INT(nnp%n_committee, c_int), nnp%myforce)
      IF (ec /= 0) CPABORT("NNP GPU: force-scatter begin failed.")

      ASSOCIATE (hstate => nnp%gpu_host_state)

         ! Radial scatter, one dispatch per populated (group, element) slot.
         DO e = 1, nnp%n_ele
            n_atoms_e = hstate%n_atoms_in_e(e)
            IF (n_atoms_e <= 0) CYCLE
            DO s = 1, nnp%rad(e)%n_symfgrp
               ASSOCIATE (symfgrp => nnp%rad(e)%symfgrp(s))
                  max_n_neigh_per_se = MAXVAL(hstate%n_rad_pa(1:n_atoms_e, s, e))
                  IF (max_n_neigh_per_se <= 0) CYCLE
                  ec = nnp_gpu_force_scatter_radial_c( &
                       INT(n_atoms_e, c_int), INT(symfgrp%n_symf, c_int), &
                       INT(max_n_neigh_per_se, c_int), INT(nnp%n_committee, c_int), &
                       INT(n_input_max, c_int), INT(nnp%num_atoms, c_int), &
                       INT(e - 1, c_int), INT(scaling_mode, c_int), &
                       REAL(nnp%scmin, c_double), REAL(nnp%scmax, c_double), &
                       INT(s - 1, c_int), INT(max_n_radgrp, c_int), &
                       INT(hstate%n_atoms_in_e_max, c_int), &
                       INT(hstate%max_n_neigh_rad, c_int), INT(istart, c_int), &
                       INT((s - 1) + max_n_radgrp*(e - 1), c_int), &
                       INT(accumulate_virial, c_int))
                  IF (ec /= 0) CPABORT("NNP GPU: radial force-scatter kernel failed.")
               END ASSOCIATE
            END DO
         END DO

         ! Angular scatter, one dispatch per populated (group, element) slot. The
         ! per-slot high-water neighbour counts and the neighbour-pair width
         ! max_n_kk are sized exactly as the angular descriptor build: max_n_ang2
         ! is unused for a homonuclear group (its second neighbour list is empty)
         ! and the pair width max_n_kk is the ang1 list width for homonuclear
         ! groups and the ang2 list width otherwise.
         DO e = 1, nnp%n_ele
            n_atoms_e = hstate%n_atoms_in_e(e)
            IF (n_atoms_e <= 0) CYCLE
            DO s = 1, nnp%ang(e)%n_symfgrp
               ASSOCIATE (symfgrp => nnp%ang(e)%symfgrp(s))
                  homo = (symfgrp%ele(1) == symfgrp%ele(2))
                  max_n_ang1 = MAXVAL(hstate%n_ang1_pa(1:n_atoms_e, s, e))
                  IF (max_n_ang1 <= 0) CYCLE
                  IF (homo) THEN
                     max_n_ang2 = 0
                     max_n_kk = MAX(1, max_n_ang1)
                  ELSE
                     max_n_ang2 = MAXVAL(hstate%n_ang2_pa(1:n_atoms_e, s, e))
                     ! Empty pair list: the SF stage skipped this group (same
                     ! guard there), so its per-(s,e) device buffers are
                     ! unpopulated, so launching the scatter would read stale
                     ! data and produce silently wrong forces. Contribution
                     ! is exactly zero; skip.
                     IF (max_n_ang2 <= 0) CYCLE
                     max_n_kk = MAX(1, max_n_ang2)
                  END IF
                  ec = nnp_gpu_force_scatter_angular_c( &
                       INT(n_atoms_e, c_int), INT(symfgrp%n_symf, c_int), &
                       INT(MAX(1, max_n_ang1), c_int), INT(max_n_kk, c_int), &
                       INT(nnp%n_committee, c_int), INT(n_input_max, c_int), &
                       INT(nnp%num_atoms, c_int), INT(e - 1, c_int), &
                       INT(nnp%n_rad(e), c_int), INT(scaling_mode, c_int), &
                       INT(MERGE(1, 0, homo), c_int), REAL(nnp%scmin, c_double), &
                       REAL(nnp%scmax, c_double), INT(s - 1, c_int), &
                       INT(max_n_anggrp, c_int), INT(hstate%n_atoms_in_e_max, c_int), &
                       INT(hstate%max_n_neigh_ang1, c_int), &
                       INT(hstate%max_n_neigh_ang2, c_int), INT(istart, c_int), &
                       INT((s - 1) + max_n_anggrp*(e - 1), c_int), &
                       INT(accumulate_virial, c_int))
                  IF (ec /= 0) CPABORT("NNP GPU: angular force-scatter kernel failed.")
               END ASSOCIATE
            END DO
         END DO

      END ASSOCIATE

      ! Read the accumulated forces back into the host force array.
      ec = nnp_gpu_force_scatter_end_c(INT(nnp%num_atoms, c_int), &
                                       INT(nnp%n_committee, c_int), nnp%myforce)
      IF (ec /= 0) CPABORT("NNP GPU: force-scatter end failed.")

      ! Stress: the scatter kernels accumulated the per-committee analytical
      ! virial into the device buffer when accumulate_virial was set. D2H it and
      ! fold into committee_stress (zeroed by the caller), so the committee
      ! reduction in nnp_calc_energy_force forms virial%pv_virial unchanged. The
      ! device layout virial_dev[9*c + 3*a + b] (a = force component, b = bond
      ! component) maps to committee_stress(b, a, c) (column-major 3x3).
      IF (accumulate_virial /= 0) THEN
         ALLOCATE (c_virial(9*nnp%n_committee))
         ec = nnp_gpu_force_scatter_virial_d2h_c(INT(nnp%n_committee, c_int), c_virial)
         IF (ec /= 0) CPABORT("NNP GPU: force-scatter virial read-back failed.")
         DO ivc = 1, nnp%n_committee
            DO iva = 1, 3
               DO ivb = 1, 3
                  nnp%committee_stress(ivb, iva, ivc) = &
                     nnp%committee_stress(ivb, iva, ivc) + &
                     c_virial(9*(ivc - 1) + 3*(iva - 1) + (ivb - 1) + 1)
               END DO
            END DO
         END DO
         DEALLOCATE (c_virial)
      END IF

      ! Out-of-training-range detection, mirroring nnp_check_extrapolation on
      ! the CPU path: the scaling kernel raised the flag on the device, the
      ! force read-back above synchronized the stream, and the caller resets
      ! output_expol before dispatch and sums it over the ranks.
      nnp%output_expol = (nnp_gpu_extrapolation_flag_c() /= 0)

      CALL timestop(handle)
#else
      MARK_USED(nnp)
      MARK_USED(istart)
      MARK_USED(mecalc)
      MARK_USED(n_input_max)
      MARK_USED(calc_stress)
      MARK_USED(handle)
      CPABORT("NNP GPU support is not compiled into this binary.")
#endif

   END SUBROUTINE nnp_gpu_energy_force

END MODULE nnp_gpu
