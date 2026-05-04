!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Data types for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_environment_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                              atomic_kind_list_release,&
                                              atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_release,&
                                              cell_retain,&
                                              cell_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_create,&
                                              molecule_kind_list_release,&
                                              molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_create,&
                                              molecule_list_release,&
                                              molecule_list_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_list_types,             ONLY: particle_list_create,&
                                              particle_list_release,&
                                              particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_environment_types'

   !> derived data types
   PUBLIC :: nnp_type
   PUBLIC :: nnp_arc_type
   PUBLIC :: nnp_neighbor_type
   PUBLIC :: nnp_neigh_grp_type
   PUBLIC :: nnp_symfgrp_type
   PUBLIC :: nnp_acsf_rad_type
   PUBLIC :: nnp_acsf_ang_type
   PUBLIC :: nnp_cell_list_cache_type
   PUBLIC :: nnp_neighbor_pair_map_type
   PUBLIC :: nnp_dGdr_grp_type
   PUBLIC :: nnp_neighbor_workspace_type
   PUBLIC :: nnp_neighbor_interface_state_type

   ! Public subroutines ***
   PUBLIC :: nnp_env_release, &
             nnp_env_set, &
             nnp_env_get, &
             nnp_cell_list_cache_release, &
             nnp_neighbor_interface_state_release, &
             nnp_workspace_release

   INTEGER, PARAMETER, PUBLIC :: &
      nnp_cut_cos = 1, &
      nnp_cut_tanh = 2

   INTEGER, PARAMETER, PUBLIC :: &
      nnp_actfnct_tanh = 1, &
      nnp_actfnct_gaus = 2, &
      nnp_actfnct_lin = 3, &
      nnp_actfnct_cos = 4, &
      nnp_actfnct_sig = 5, &
      nnp_actfnct_invsig = 6, &
      nnp_actfnct_exp = 7, &
      nnp_actfnct_softplus = 8, &
      nnp_actfnct_quad = 9

! **************************************************************************************************
!> \brief Main data type collecting all relevant data for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_type
      TYPE(nnp_acsf_rad_type), DIMENSION(:), POINTER      :: rad => NULL() ! DIM(n_ele)
      TYPE(nnp_acsf_ang_type), DIMENSION(:), POINTER      :: ang => NULL() ! DIM(n_ele)
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_rad ! # radial symfnct for this element
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_ang ! # angular symfnct for this element
      INTEGER                                             :: n_ele = -1 ! # elements
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele ! elements(n_ele)
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele ! elements(n_ele)
      LOGICAL                                             :: scale_acsf = .FALSE.
      LOGICAL                                             :: scale_sigma_acsf = .FALSE.
      LOGICAL                                             :: center_acsf = .FALSE.
      LOGICAL                                             :: normnodes = .FALSE.
      INTEGER                                             :: n_radgrp = -1
      INTEGER                                             :: n_anggrp = -1
      INTEGER                                             :: cut_type = -1 ! cutofftype
      REAL(KIND=dp)                                       :: eshortmin = -1.0_dp
      REAL(KIND=dp)                                       :: eshortmax = -1.0_dp
      REAL(KIND=dp)                                       :: scmax = -1.0_dp !scale
      REAL(KIND=dp)                                       :: scmin = -1.0_dp !scale
      REAL(KIND=dp)                                       :: max_cut = -1.0_dp !largest cutoff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: atom_energies !DIM(n_ele)
      TYPE(nnp_arc_type), POINTER, DIMENSION(:)           :: arc => NULL() ! DIM(n_ele)
      INTEGER                                             :: n_committee = -1
      INTEGER                                             :: n_hlayer = -1
      INTEGER                                             :: n_layer = -1
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_hnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: actfnct
      INTEGER                                             :: expol = -1 ! extrapolation counter
      LOGICAL                                             :: output_expol = .FALSE. ! output extrapolation
      ! structures for calculation
      INTEGER                                             :: num_atoms = -1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: atomic_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: committee_energy
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: ele_ind, nuc_atoms, sort, sort_inv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coord
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: myforce
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: committee_forces, committee_stress
      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                        :: atoms
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: nnp_forces
      REAL(KIND=dp)                                       :: nnp_potential_energy = -1.0_dp
      TYPE(cp_subsys_type), POINTER                       :: subsys => NULL()
      TYPE(section_vals_type), POINTER                    :: nnp_input => NULL()
      TYPE(section_vals_type), POINTER                    :: force_env_input => NULL()
      TYPE(cell_type), POINTER                            :: cell => NULL()
      TYPE(cell_type), POINTER                            :: cell_ref => NULL()
      LOGICAL                                             :: use_ref_cell = .FALSE.
      ! bias
      LOGICAL                                             :: bias = .FALSE.
      LOGICAL                                             :: bias_align = .FALSE.
      REAL(KIND=dp)                                       :: bias_energy = -1.0_dp
      REAL(KIND=dp)                                       :: bias_kb = -1.0_dp
      REAL(KIND=dp)                                       :: bias_sigma0 = -1.0_dp
      REAL(KIND=dp)                                       :: bias_sigma = -1.0_dp
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: bias_forces
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: bias_e_avrg
      ! Per-nnp persistent neighbour-finder state. Owned by nnp_type so that
      ! multiple &NNP force_evals (committee / MIX) carry independent caches
      ! and so the lifetime tracks nnp_env_release rather than module-global
      ! SAVE state. Operated on by nnp_cell_list and nnp_neighbor_interface.
      ! ALLOCATABLE (rather than embedded) only to satisfy gfortran's rule
      ! that direct-embedded type members must be defined before use; the
      ! component types live further down this module file. Allocated by
      ! nnp_prepare_neighbor_cache on first force eval, freed in
      ! nnp_env_release.
      TYPE(nnp_cell_list_cache_type), ALLOCATABLE         :: cell_list_cache
      TYPE(nnp_neighbor_interface_state_type), ALLOCATABLE :: neighbor_interface_state
   END TYPE nnp_type

! **************************************************************************************************
!> \brief Symmetry functions group type
!> \param n_symf - # of associated sym fncts
!> \param symf   - indices of associated sym fncts       DIM(nsymf)
!> \param ele    - elements indices          rad:DIM(2), ang:DIM(3)
!> \param cutoff - associated cutoff value
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_symfgrp_type
      INTEGER                                             :: n_symf = -1
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: symf
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: ele_ind
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE         :: ele
      REAL(KIND=dp)                                       :: cutoff = -1.0_dp
      ! Packed per-member parameters: contiguous arrays sized n_symf so the
      ! inner SF loop reads sequential memory instead of chasing 4 pointers
      ! into ang(ind)%{eta,zeta,lam,prefzeta} via the symf(sf) index.
      ! pack_use_int_zeta and pack_izeta let us skip NINT inside the inner loop.
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pack_eta
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pack_zeta
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pack_lam
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pack_prefzeta
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: pack_izeta
      LOGICAL, DIMENSION(:), ALLOCATABLE                  :: pack_use_int_zeta
      ! Group-shared Hermite cubic spline tables (radial groups only).
      ! All SFs in a radial group share the same cutoff (= grp%cutoff), so
      ! they all use the same uniform grid: x in [0, spline_x_max] with
      ! spline_n nodes, spacing spline_dx. The per-SF y(r) and y'(r)
      ! tables are packed sf-first into 2D arrays so the inner SF loop in
      ! nnp_calc_rad streams them contiguously after computing the
      ! interpolation params (i, t, h00..h11) ONCE per call.
      LOGICAL                                             :: spline_built = .FALSE.
      INTEGER                                             :: spline_n = 0
      REAL(KIND=dp)                                       :: spline_dx = 1.0_dp
      REAL(KIND=dp)                                       :: spline_dx_inv = 1.0_dp
      REAL(KIND=dp)                                       :: spline_x_max = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: spline_y    ! (n_symf, n_grid)
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: spline_dy   ! (n_symf, n_grid)
   END TYPE nnp_symfgrp_type

! **************************************************************************************************
!> \brief Set of radial symmetry function type
!> \param y       - acsf value                                     - DIM(n_rad)
!> \param funccut - distance cutoff                           bohr - DIM(n_rad)
!> \param eta     - eta parameter of radial sym fncts      bohr^-2 - DIM(n_rad)
!> \param rs      - r shift parameter of radial sym fncts     bohr - DIM(n_rad)
!> \param loc_min - minimum of the sym fnct                          DIM(n_rad)
!> \param loc_max - maximum of the sym fnct                          DIM(n_rad)
!> \param loc_av  - average of the sym fnct                          DIM(n_rad)
!> \param sigma   - SD of the sym fnc                                DIM(n_rad)
!> \param ele     - element associated to the sym fnct               DIM(n_rad)
!> \param nuc_ele - associated atomic number                         DIM(n_rad)
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_acsf_rad_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: y
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: funccut
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: eta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: rs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_av
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: sigma
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele
      INTEGER                                             :: n_symfgrp = -1
      TYPE(nnp_symfgrp_type), DIMENSION(:), ALLOCATABLE   :: symfgrp
   END TYPE nnp_acsf_rad_type

! **************************************************************************************************
!> \brief Set of angular symmetry function type
!> \param y         - acsf value                                  - DIM(n_ang)
!> \param funccut   - distance cutoff                        bohr - DIM(n_ang)
!> \param eta       - eta  param. of angular sym fncts    bohr^-2 - DIM(n_ang)
!> \param zeta      - zeta param. of angular sym fncts              DIM(n_ang)
!> \param lam       - lambda  param. of angular sym fncts           DIM(n_ang)
!> \param loc_min   - minimum of the sym fnct                       DIM(n_ang)
!> \param loc_max   - maximum of the sym fnct                       DIM(n_ang)
!> \param loc_av    - average of the sym fnct                       DIM(n_ang)
!> \param sigma     - SD of the sym fnc                             DIM(n_ang)
!> \param ele1,ele2 - elements associated to the sym fnct           DIM(n_ang)
!> \param nuc_ele2, nuc_ele2 - associated atomic numbers            DIM(n_ang)
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_acsf_ang_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: y
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: funccut
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: eta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: zeta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: prefzeta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: lam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: loc_av
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: sigma
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele1
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele2
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele1
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele2
      INTEGER                                             :: n_symfgrp = -1
      TYPE(nnp_symfgrp_type), DIMENSION(:), ALLOCATABLE   :: symfgrp
   END TYPE nnp_acsf_ang_type

! **************************************************************************************************
!> \brief Per-(symmetry-function-group) dense neighbor container.
!>        cap is the allocated slab size; the live count is in
!>        nnp_neighbor_type%n_rad/n_ang1/n_ang2 for the matching group index.
!>        ind(j) is the absolute atom index of the j-th neighbor in this group;
!>        dist(1:3, j) is the displacement vector and dist(4, j) is its norm.
!>        Sized lazily by nnp_neigh_grp_grow (in nnp_neighbor_interface) so the
!>        slabs settle at the true peak per group rather than the global image-pool
!>        capacity.
! **************************************************************************************************
   TYPE nnp_neigh_grp_type
      INTEGER                                             :: cap = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: ind   ! (cap)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: dist  ! (4, cap)
   END TYPE nnp_neigh_grp_type

! **************************************************************************************************
!> \brief Contains neighbors list of an atom (per-group dense layout).
!> \param n_rad/n_ang1/n_ang2 - live neighbor counts per symfgrp
!> \param rad/ang1/ang2       - per-group dense (ind, dist) containers
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_neighbor_type
      INTEGER, DIMENSION(3)                                :: pbc_copies = -1
      INTEGER, DIMENSION(:), ALLOCATABLE                   :: n_rad
      INTEGER, DIMENSION(:), ALLOCATABLE                   :: n_ang1
      INTEGER, DIMENSION(:), ALLOCATABLE                   :: n_ang2
      TYPE(nnp_neigh_grp_type), ALLOCATABLE, DIMENSION(:)  :: rad
      TYPE(nnp_neigh_grp_type), ALLOCATABLE, DIMENSION(:)  :: ang1
      TYPE(nnp_neigh_grp_type), ALLOCATABLE, DIMENSION(:)  :: ang2
   END TYPE nnp_neighbor_type

! **************************************************************************************************
!> \brief Linked-cell + Verlet-skin cache for the NNP descriptor neighbour walk.
!>        Persisted across MD steps so the image pool, bin geometry, and
!>        reference positions survive between force evaluations and the
!>        head/next chain can be reused as long as no atom has drifted more
!>        than skin/2. Operated on by nnp_cell_list.F.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   TYPE nnp_cell_list_cache_type
      LOGICAL                                             :: initialized = .FALSE.
      INTEGER                                             :: num_atoms = -1
      INTEGER                                             :: n_images = 0
      INTEGER                                             :: n_cells = 1
      INTEGER, DIMENSION(3)                               :: exact_pbc_copies = 0
      INTEGER, DIMENSION(3)                               :: list_pbc_copies = 0
      INTEGER, DIMENSION(3)                               :: image_copies = 0
      INTEGER, DIMENSION(3)                               :: nbin = 1
      INTEGER, DIMENSION(3)                               :: bin_span = 0
      INTEGER, DIMENSION(3)                               :: perd = 0
      LOGICAL                                             :: orthorhombic = .FALSE.
      REAL(KIND=dp)                                       :: exact_cutoff = -1.0_dp
      REAL(KIND=dp)                                       :: list_cutoff = -1.0_dp
      REAL(KIND=dp)                                       :: verlet_skin = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)                         :: lower = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)                         :: upper = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)                         :: bin_width = 1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                      :: h_inv = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coord_primary
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coord_scaled
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: ref_coord_primary
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: image_atom
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: head
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: next
      INTEGER, ALLOCATABLE, DIMENSION(:, :)               :: image_shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: image_translation
   END TYPE nnp_cell_list_cache_type

! **************************************************************************************************
!> \brief Species-pair routing table. For one (central element, neighbour element)
!>        pair, lists which radial / angular symmetry-function groups need that
!>        neighbour and the worst-case relevant cutoff.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   TYPE nnp_neighbor_pair_map_type
      INTEGER                                             :: n_rad = 0
      INTEGER                                             :: n_ang1 = 0
      INTEGER                                             :: n_ang2 = 0
      REAL(KIND=dp)                                       :: max_relevant_cutoff = 0.0_dp
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: rad_groups
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: ang1_groups
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: ang2_groups
   END TYPE nnp_neighbor_pair_map_type

! **************************************************************************************************
!> \brief Per-(element, symmetry-function-group) dG_k/dr buffer. Sized to the
!>        actual n_symf for the group (no max_*_symf padding) and to the
!>        observed maximum number of neighbors that fall in this group (no
!>        image-pool over-allocation). Grown lazily with a 1.5x policy so that
!>        after a handful of MD steps everything settles at the true peak.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   TYPE nnp_dGdr_grp_type
      INTEGER                                             :: cap = 0
      INTEGER                                             :: n_symf = 0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: data    ! (3, n_symf, cap)
   END TYPE nnp_dGdr_grp_type

! **************************************************************************************************
!> \brief Reusable per-element scratch / persistent caches for the ACSF
!>        descriptor and force assembly. See nnp_neighbor_interface for the
!>        routines that populate and grow these buffers.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   TYPE nnp_neighbor_workspace_type
      INTEGER                                             :: max_rad_symf = 0
      INTEGER                                             :: max_ang_symf = 0
      INTEGER                                             :: n_input_nodes = 0
      ! High-water mark of the per-element angular-cache slab. Tracked
      ! independently of the image-pool capacity so the 1D cutoff caches
      ! settle at the true per-element peak rather than the global pool.
      INTEGER                                             :: cache_cap = 0
      TYPE(nnp_neighbor_type)                             :: neighbor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: radial_sym
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: radial_force
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: angular_sym
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: angular_force
      ! Per-element persistent cutoff caches reused across every central
      ! atom of this element. Sized lazily by nnp_workspace_grow_caches
      ! to MAX(n_ang1) / MAX(n_ang2) over the angular groups -- not to
      ! the much larger image-pool capacity.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: fc_cache1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: dfc_cache1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: fc_cache2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: dfc_cache2
      ! n2p2-style per-neighbor dG_k/dr storage. Replaces the global
      ! dsymdxyz(3, n_input_nodes, num_atoms) slab. The per-group buffers
      ! are dense in (3, n_symf_for_group, n_neighbors_in_group) -- no
      ! max_*_symf padding and no image-pool over-allocation -- which
      ! shrinks the working set by ~100x and lets every group walk be a
      ! contiguous L1-friendly stride.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: self_dGdr
      TYPE(nnp_dGdr_grp_type), ALLOCATABLE, DIMENSION(:)  :: dGdr_rad
      TYPE(nnp_dGdr_grp_type), ALLOCATABLE, DIMENSION(:)  :: dGdr_ang_jj
      TYPE(nnp_dGdr_grp_type), ALLOCATABLE, DIMENSION(:)  :: dGdr_ang_kk
   END TYPE nnp_neighbor_workspace_type

! **************************************************************************************************
!> \brief Persistent neighbour-interface state. Owned by nnp_type so that the
!>        per-nnp pair-routing tables and per-element workspaces survive
!>        across MD steps without leaking between independent &NNP
!>        force_evals. Operated on by nnp_neighbor_interface.F.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   TYPE nnp_neighbor_interface_state_type
      LOGICAL                                             :: initialized = .FALSE.
      INTEGER                                             :: n_ele = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_rad
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_ang
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_radgrp
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_anggrp
      TYPE(nnp_neighbor_pair_map_type), ALLOCATABLE, &
         DIMENSION(:, :)                                  :: pair_map
      TYPE(nnp_neighbor_workspace_type), ALLOCATABLE, &
         DIMENSION(:)                                     :: workspace
   END TYPE nnp_neighbor_interface_state_type

! **************************************************************************************************
!> \brief Data type for artificial neural networks
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_arc_type
      TYPE(nnp_arc_layer_type), POINTER, DIMENSION(:)     :: layer => NULL() ! DIM(n_layer)
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_nodes
   END TYPE nnp_arc_type

! **************************************************************************************************
!> \brief Data type for individual layer
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_arc_layer_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: weights ! node weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: bweights ! bias weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: node ! DIM(n_nodes)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: node_grad ! DIM(n_nodes)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: tmp_der ! DIM(n_sym,n_nodes)
   END TYPE nnp_arc_layer_type

CONTAINS

! **************************************************************************************************
!> \brief Release data structure that holds all the information for neural
!>        network potentials
!> \param nnp_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_env_release(nnp_env)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp_env

      INTEGER                                            :: i, j

      IF (ASSOCIATED(nnp_env%rad)) THEN
         DO i = 1, nnp_env%n_ele
            DO j = 1, nnp_env%rad(i)%n_symfgrp
               DEALLOCATE (nnp_env%rad(i)%symfgrp(j)%symf, &
                           nnp_env%rad(i)%symfgrp(j)%ele, &
                           nnp_env%rad(i)%symfgrp(j)%ele_ind)
               IF (ALLOCATED(nnp_env%rad(i)%symfgrp(j)%pack_eta)) &
                  DEALLOCATE (nnp_env%rad(i)%symfgrp(j)%pack_eta)
               IF (ALLOCATED(nnp_env%rad(i)%symfgrp(j)%spline_y)) &
                  DEALLOCATE (nnp_env%rad(i)%symfgrp(j)%spline_y)
               IF (ALLOCATED(nnp_env%rad(i)%symfgrp(j)%spline_dy)) &
                  DEALLOCATE (nnp_env%rad(i)%symfgrp(j)%spline_dy)
               nnp_env%rad(i)%symfgrp(j)%spline_built = .FALSE.
            END DO
            DEALLOCATE (nnp_env%rad(i)%y, &
                        nnp_env%rad(i)%funccut, &
                        nnp_env%rad(i)%eta, &
                        nnp_env%rad(i)%rs, &
                        nnp_env%rad(i)%loc_min, &
                        nnp_env%rad(i)%loc_max, &
                        nnp_env%rad(i)%loc_av, &
                        nnp_env%rad(i)%sigma, &
                        nnp_env%rad(i)%ele, &
                        nnp_env%rad(i)%nuc_ele, &
                        nnp_env%rad(i)%symfgrp)
         END DO
         DEALLOCATE (nnp_env%rad)
      END IF

      IF (ASSOCIATED(nnp_env%ang)) THEN
         DO i = 1, nnp_env%n_ele
            DO j = 1, nnp_env%ang(i)%n_symfgrp
               DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%symf, &
                           nnp_env%ang(i)%symfgrp(j)%ele, &
                           nnp_env%ang(i)%symfgrp(j)%ele_ind)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_eta)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_eta)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_zeta)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_zeta)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_lam)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_lam)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_prefzeta)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_prefzeta)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_izeta)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_izeta)
               IF (ALLOCATED(nnp_env%ang(i)%symfgrp(j)%pack_use_int_zeta)) &
                  DEALLOCATE (nnp_env%ang(i)%symfgrp(j)%pack_use_int_zeta)
            END DO
            DEALLOCATE (nnp_env%ang(i)%y, &
                        nnp_env%ang(i)%funccut, &
                        nnp_env%ang(i)%eta, &
                        nnp_env%ang(i)%zeta, &
                        nnp_env%ang(i)%prefzeta, &
                        nnp_env%ang(i)%lam, &
                        nnp_env%ang(i)%loc_min, &
                        nnp_env%ang(i)%loc_max, &
                        nnp_env%ang(i)%loc_av, &
                        nnp_env%ang(i)%sigma, &
                        nnp_env%ang(i)%ele1, &
                        nnp_env%ang(i)%ele2, &
                        nnp_env%ang(i)%nuc_ele1, &
                        nnp_env%ang(i)%nuc_ele2, &
                        nnp_env%ang(i)%symfgrp)
         END DO
         DEALLOCATE (nnp_env%ang)
      END IF

      IF (ASSOCIATED(nnp_env%arc)) THEN
         DO i = 1, nnp_env%n_ele
            IF (ASSOCIATED(nnp_env%arc(i)%layer)) THEN
               DO j = 1, nnp_env%n_layer
                  IF (ALLOCATED(nnp_env%arc(i)%layer(j)%node)) THEN
                     DEALLOCATE (nnp_env%arc(i)%layer(j)%node)
                  END IF
                  IF (ALLOCATED(nnp_env%arc(i)%layer(j)%node_grad)) THEN
                     DEALLOCATE (nnp_env%arc(i)%layer(j)%node_grad)
                  END IF
                  IF (ALLOCATED(nnp_env%arc(i)%layer(j)%weights)) THEN
                     DEALLOCATE (nnp_env%arc(i)%layer(j)%weights)
                  END IF
                  IF (ALLOCATED(nnp_env%arc(i)%layer(j)%bweights)) THEN
                     DEALLOCATE (nnp_env%arc(i)%layer(j)%bweights)
                  END IF
                  IF (ALLOCATED(nnp_env%arc(i)%layer(j)%tmp_der)) THEN
                     DEALLOCATE (nnp_env%arc(i)%layer(j)%tmp_der)
                  END IF
               END DO
               DEALLOCATE (nnp_env%arc(i)%layer, &
                           nnp_env%arc(i)%n_nodes)
            END IF
         END DO
         DEALLOCATE (nnp_env%arc)
      END IF

      IF (ALLOCATED(nnp_env%ele)) DEALLOCATE (nnp_env%ele)
      IF (ALLOCATED(nnp_env%nuc_ele)) DEALLOCATE (nnp_env%nuc_ele)
      IF (ALLOCATED(nnp_env%n_hnodes)) DEALLOCATE (nnp_env%n_hnodes)
      IF (ALLOCATED(nnp_env%actfnct)) DEALLOCATE (nnp_env%actfnct)
      IF (ALLOCATED(nnp_env%nnp_forces)) DEALLOCATE (nnp_env%nnp_forces)
      IF (ALLOCATED(nnp_env%atomic_energy)) DEALLOCATE (nnp_env%atomic_energy)
      IF (ALLOCATED(nnp_env%committee_energy)) DEALLOCATE (nnp_env%committee_energy)
      IF (ALLOCATED(nnp_env%ele_ind)) DEALLOCATE (nnp_env%ele_ind)
      IF (ALLOCATED(nnp_env%nuc_atoms)) DEALLOCATE (nnp_env%nuc_atoms)
      IF (ALLOCATED(nnp_env%sort)) DEALLOCATE (nnp_env%sort)
      IF (ALLOCATED(nnp_env%sort_inv)) DEALLOCATE (nnp_env%sort_inv)
      IF (ALLOCATED(nnp_env%coord)) DEALLOCATE (nnp_env%coord)
      IF (ALLOCATED(nnp_env%myforce)) DEALLOCATE (nnp_env%myforce)
      IF (ALLOCATED(nnp_env%committee_forces)) DEALLOCATE (nnp_env%committee_forces)
      IF (ALLOCATED(nnp_env%committee_stress)) DEALLOCATE (nnp_env%committee_stress)
      IF (ALLOCATED(nnp_env%atoms)) DEALLOCATE (nnp_env%atoms)

      IF (ALLOCATED(nnp_env%cell_list_cache)) THEN
         CALL nnp_cell_list_cache_release(nnp_env%cell_list_cache)
         DEALLOCATE (nnp_env%cell_list_cache)
      END IF
      IF (ALLOCATED(nnp_env%neighbor_interface_state)) THEN
         CALL nnp_neighbor_interface_state_release(nnp_env%neighbor_interface_state)
         DEALLOCATE (nnp_env%neighbor_interface_state)
      END IF

      IF (ASSOCIATED(nnp_env%subsys)) THEN
         CALL cp_subsys_release(nnp_env%subsys)
      END IF
      IF (ASSOCIATED(nnp_env%cell)) THEN
         CALL cell_release(nnp_env%cell)
      END IF
      IF (ASSOCIATED(nnp_env%cell_ref)) THEN
         CALL cell_release(nnp_env%cell_ref)
      END IF

   END SUBROUTINE nnp_env_release

! **************************************************************************************************
!> \brief Free the allocatable parts of a cell-list cache. Co-located with the
!>        type definition to avoid creating a USE-cycle through nnp_cell_list.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_cache_release(cache)
      TYPE(nnp_cell_list_cache_type), INTENT(INOUT)      :: cache

      IF (ALLOCATED(cache%coord_primary)) DEALLOCATE (cache%coord_primary)
      IF (ALLOCATED(cache%coord_scaled)) DEALLOCATE (cache%coord_scaled)
      IF (ALLOCATED(cache%ref_coord_primary)) DEALLOCATE (cache%ref_coord_primary)
      IF (ALLOCATED(cache%image_atom)) DEALLOCATE (cache%image_atom)
      IF (ALLOCATED(cache%head)) DEALLOCATE (cache%head)
      IF (ALLOCATED(cache%next)) DEALLOCATE (cache%next)
      IF (ALLOCATED(cache%image_shift)) DEALLOCATE (cache%image_shift)
      IF (ALLOCATED(cache%image_translation)) DEALLOCATE (cache%image_translation)
      cache%initialized = .FALSE.
      cache%num_atoms = -1
      cache%n_images = 0
   END SUBROUTINE nnp_cell_list_cache_release

! **************************************************************************************************
!> \brief Free the allocatable parts of a neighbour-interface state. Co-located
!>        with the type definition to avoid creating a USE-cycle through
!>        nnp_neighbor_interface.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_interface_state_release(state)
      TYPE(nnp_neighbor_interface_state_type), INTENT(INOUT) :: state

      INTEGER                                            :: i, j

      IF (ALLOCATED(state%pair_map)) THEN
         DO i = 1, SIZE(state%pair_map, 1)
            DO j = 1, SIZE(state%pair_map, 2)
               IF (ALLOCATED(state%pair_map(i, j)%rad_groups)) DEALLOCATE (state%pair_map(i, j)%rad_groups)
               IF (ALLOCATED(state%pair_map(i, j)%ang1_groups)) DEALLOCATE (state%pair_map(i, j)%ang1_groups)
               IF (ALLOCATED(state%pair_map(i, j)%ang2_groups)) DEALLOCATE (state%pair_map(i, j)%ang2_groups)
            END DO
         END DO
         DEALLOCATE (state%pair_map)
      END IF

      IF (ALLOCATED(state%workspace)) THEN
         DO i = 1, SIZE(state%workspace)
            CALL nnp_workspace_release(state%workspace(i))
         END DO
         DEALLOCATE (state%workspace)
      END IF

      IF (ALLOCATED(state%n_rad)) DEALLOCATE (state%n_rad)
      IF (ALLOCATED(state%n_ang)) DEALLOCATE (state%n_ang)
      IF (ALLOCATED(state%n_radgrp)) DEALLOCATE (state%n_radgrp)
      IF (ALLOCATED(state%n_anggrp)) DEALLOCATE (state%n_anggrp)

      state%initialized = .FALSE.
      state%n_ele = 0
   END SUBROUTINE nnp_neighbor_interface_state_release

! **************************************************************************************************
!> \brief Free the allocatable parts of one per-element workspace.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_workspace_release(workspace)
      TYPE(nnp_neighbor_workspace_type), INTENT(INOUT)   :: workspace

      INTEGER                                            :: s

      IF (ALLOCATED(workspace%neighbor%rad)) THEN
         DO s = 1, SIZE(workspace%neighbor%rad)
            IF (ALLOCATED(workspace%neighbor%rad(s)%ind)) DEALLOCATE (workspace%neighbor%rad(s)%ind)
            IF (ALLOCATED(workspace%neighbor%rad(s)%dist)) DEALLOCATE (workspace%neighbor%rad(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%rad)
      END IF
      IF (ALLOCATED(workspace%neighbor%ang1)) THEN
         DO s = 1, SIZE(workspace%neighbor%ang1)
            IF (ALLOCATED(workspace%neighbor%ang1(s)%ind)) DEALLOCATE (workspace%neighbor%ang1(s)%ind)
            IF (ALLOCATED(workspace%neighbor%ang1(s)%dist)) DEALLOCATE (workspace%neighbor%ang1(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%ang1)
      END IF
      IF (ALLOCATED(workspace%neighbor%ang2)) THEN
         DO s = 1, SIZE(workspace%neighbor%ang2)
            IF (ALLOCATED(workspace%neighbor%ang2(s)%ind)) DEALLOCATE (workspace%neighbor%ang2(s)%ind)
            IF (ALLOCATED(workspace%neighbor%ang2(s)%dist)) DEALLOCATE (workspace%neighbor%ang2(s)%dist)
         END DO
         DEALLOCATE (workspace%neighbor%ang2)
      END IF
      IF (ALLOCATED(workspace%neighbor%n_rad)) DEALLOCATE (workspace%neighbor%n_rad)
      IF (ALLOCATED(workspace%neighbor%n_ang1)) DEALLOCATE (workspace%neighbor%n_ang1)
      IF (ALLOCATED(workspace%neighbor%n_ang2)) DEALLOCATE (workspace%neighbor%n_ang2)
      workspace%neighbor%pbc_copies = -1

      IF (ALLOCATED(workspace%radial_sym)) DEALLOCATE (workspace%radial_sym)
      IF (ALLOCATED(workspace%radial_force)) DEALLOCATE (workspace%radial_force)
      IF (ALLOCATED(workspace%angular_sym)) DEALLOCATE (workspace%angular_sym)
      IF (ALLOCATED(workspace%angular_force)) DEALLOCATE (workspace%angular_force)
      IF (ALLOCATED(workspace%fc_cache1)) DEALLOCATE (workspace%fc_cache1)
      IF (ALLOCATED(workspace%dfc_cache1)) DEALLOCATE (workspace%dfc_cache1)
      IF (ALLOCATED(workspace%fc_cache2)) DEALLOCATE (workspace%fc_cache2)
      IF (ALLOCATED(workspace%dfc_cache2)) DEALLOCATE (workspace%dfc_cache2)
      IF (ALLOCATED(workspace%self_dGdr)) DEALLOCATE (workspace%self_dGdr)

      IF (ALLOCATED(workspace%dGdr_rad)) THEN
         DO s = 1, SIZE(workspace%dGdr_rad)
            IF (ALLOCATED(workspace%dGdr_rad(s)%data)) DEALLOCATE (workspace%dGdr_rad(s)%data)
         END DO
         DEALLOCATE (workspace%dGdr_rad)
      END IF
      IF (ALLOCATED(workspace%dGdr_ang_jj)) THEN
         DO s = 1, SIZE(workspace%dGdr_ang_jj)
            IF (ALLOCATED(workspace%dGdr_ang_jj(s)%data)) DEALLOCATE (workspace%dGdr_ang_jj(s)%data)
         END DO
         DEALLOCATE (workspace%dGdr_ang_jj)
      END IF
      IF (ALLOCATED(workspace%dGdr_ang_kk)) THEN
         DO s = 1, SIZE(workspace%dGdr_ang_kk)
            IF (ALLOCATED(workspace%dGdr_ang_kk(s)%data)) DEALLOCATE (workspace%dGdr_ang_kk(s)%data)
         END DO
         DEALLOCATE (workspace%dGdr_ang_kk)
      END IF

      workspace%max_rad_symf = 0
      workspace%max_ang_symf = 0
      workspace%n_input_nodes = 0
      workspace%cache_cap = 0
   END SUBROUTINE nnp_workspace_release

! **************************************************************************************************
!> \brief Returns various attributes of the nnp environment
!> \param nnp_env ...
!> \param nnp_forces ...
!> \param subsys the particles, molecules,... of this environment
!> \param atomic_kind_set The set of all atomic kinds involved
!> \param particle_set The set of all particles
!> \param local_particles All particles on this particular node
!> \param molecule_kind_set The set of all different molecule kinds involved
!> \param molecule_set The set of all molecules
!> \param local_molecules All molecules on this particular node
!> \param nnp_input ...
!> \param force_env_input Pointer to the force_env input section
!> \param cell The simulation cell
!> \param cell_ref The reference simulation cell
!> \param use_ref_cell Logical which indicates if reference
!>                      simulation cell is used
!> \param nnp_potential_energy ...
!> \param virial Dummy virial pointer
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
!> \note
!>      For possible missing arguments see the attributes of
!>      nnp_type
! **************************************************************************************************
   SUBROUTINE nnp_env_get(nnp_env, nnp_forces, subsys, &
                          atomic_kind_set, particle_set, local_particles, &
                          molecule_kind_set, molecule_set, local_molecules, &
                          nnp_input, force_env_input, cell, cell_ref, &
                          use_ref_cell, nnp_potential_energy, virial)

      TYPE(nnp_type), INTENT(IN)                         :: nnp_env
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: nnp_forces
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles
      TYPE(molecule_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_molecules
      TYPE(section_vals_type), OPTIONAL, POINTER         :: nnp_input, force_env_input
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_ref_cell
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: nnp_potential_energy
      TYPE(virial_type), OPTIONAL, POINTER               :: virial

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      NULLIFY (atomic_kinds, particles, molecules, molecule_kinds)

      IF (PRESENT(nnp_potential_energy)) THEN
         nnp_potential_energy = nnp_env%nnp_potential_energy
      END IF
      IF (PRESENT(nnp_forces)) nnp_forces = nnp_env%nnp_forces

      ! note cell will be overwritten if subsys is associated
      ! helium_env uses nnp without subsys
      IF (PRESENT(cell)) cell => nnp_env%cell

      IF (PRESENT(subsys)) subsys => nnp_env%subsys
      IF (ASSOCIATED(nnp_env%subsys)) THEN
         CALL cp_subsys_get(nnp_env%subsys, &
                            atomic_kinds=atomic_kinds, &
                            particles=particles, &
                            molecule_kinds=molecule_kinds, &
                            molecules=molecules, &
                            local_molecules=local_molecules, &
                            local_particles=local_particles, &
                            virial=virial, &
                            cell=cell)
      END IF
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds%els
      IF (PRESENT(molecule_set)) molecule_set => molecules%els

      IF (PRESENT(nnp_input)) nnp_input => nnp_env%nnp_input
      IF (PRESENT(force_env_input)) force_env_input => nnp_env%force_env_input
      IF (PRESENT(cell_ref)) cell_ref => nnp_env%cell_ref
      IF (PRESENT(use_ref_cell)) use_ref_cell = nnp_env%use_ref_cell

   END SUBROUTINE nnp_env_get

! **************************************************************************************************
!> \brief Sets various attributes of the nnp environment
!> \param nnp_env ...
!> \param nnp_forces ...
!> \param subsys the particles, molecules,... of this environment
!> \param atomic_kind_set The set of all atomic kinds involved
!> \param particle_set The set of all particles
!> \param local_particles All particles on this particular node
!> \param molecule_kind_set The set of all different molecule kinds involved
!> \param molecule_set The set of all molecules
!> \param local_molecules All molecules on this particular node
!> \param nnp_input ...
!> \param force_env_input Pointer to the force_env input section
!> \param cell ...
!> \param cell_ref The reference simulation cell
!> \param use_ref_cell Logical which indicates if reference
!>                      simulation cell is used
!> \param nnp_potential_energy ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
!> \note
!>   For possible missing arguments see the attributes of nnp_type
! **************************************************************************************************
   SUBROUTINE nnp_env_set(nnp_env, nnp_forces, subsys, &
                          atomic_kind_set, particle_set, local_particles, &
                          molecule_kind_set, molecule_set, local_molecules, &
                          nnp_input, force_env_input, cell, cell_ref, &
                          use_ref_cell, nnp_potential_energy)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp_env
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: nnp_forces
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles
      TYPE(molecule_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_molecules
      TYPE(section_vals_type), OPTIONAL, POINTER         :: nnp_input, force_env_input
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_ref_cell
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: nnp_potential_energy

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      IF (PRESENT(nnp_potential_energy)) THEN
         nnp_env%nnp_potential_energy = nnp_potential_energy
      END IF
      IF (PRESENT(nnp_forces)) nnp_env%nnp_forces(:, :) = nnp_forces

      IF (PRESENT(subsys)) THEN
         IF (ASSOCIATED(nnp_env%subsys)) THEN
         IF (.NOT. ASSOCIATED(nnp_env%subsys, subsys)) THEN
            CALL cp_subsys_release(nnp_env%subsys)
         END IF
         END IF
         nnp_env%subsys => subsys
      END IF
      IF (PRESENT(cell)) THEN
         IF (ASSOCIATED(cell)) THEN
            CALL cell_retain(cell)
            CALL cell_release(nnp_env%cell)
            nnp_env%cell => cell
         END IF
         IF (ASSOCIATED(nnp_env%subsys)) THEN
            CALL cp_subsys_set(nnp_env%subsys, cell=cell)
         END IF
      END IF
      IF (PRESENT(atomic_kind_set)) THEN
         CALL atomic_kind_list_create(atomic_kinds, els_ptr=atomic_kind_set)
         CALL cp_subsys_set(nnp_env%subsys, atomic_kinds=atomic_kinds)
         CALL atomic_kind_list_release(atomic_kinds)
      END IF
      IF (PRESENT(particle_set)) THEN
         CALL particle_list_create(particles, els_ptr=particle_set)
         CALL cp_subsys_set(nnp_env%subsys, particles=particles)
         CALL particle_list_release(particles)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
         CALL molecule_kind_list_create(molecule_kinds, els_ptr=molecule_kind_set)
         CALL cp_subsys_set(nnp_env%subsys, molecule_kinds=molecule_kinds)
         CALL molecule_kind_list_release(molecule_kinds)
      END IF
      IF (PRESENT(molecule_set)) THEN
         CALL molecule_list_create(molecules, els_ptr=molecule_set)
         CALL cp_subsys_set(nnp_env%subsys, molecules=molecules)
         CALL molecule_list_release(molecules)
      END IF
      IF (PRESENT(local_particles)) THEN
         CALL cp_subsys_set(nnp_env%subsys, local_particles=local_particles)
      END IF
      IF (PRESENT(local_molecules)) THEN
         CALL cp_subsys_set(nnp_env%subsys, local_molecules=local_molecules)
      END IF

      IF (PRESENT(nnp_input)) nnp_env%nnp_input => nnp_input
      IF (PRESENT(force_env_input)) THEN
         nnp_env%force_env_input => force_env_input
      END IF
      IF (PRESENT(cell_ref)) THEN
         CALL cell_retain(cell_ref)
         CALL cell_release(nnp_env%cell_ref)
         nnp_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(use_ref_cell)) nnp_env%use_ref_cell = use_ref_cell
   END SUBROUTINE nnp_env_set

END MODULE nnp_environment_types
