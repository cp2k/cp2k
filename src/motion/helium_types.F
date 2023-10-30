!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Data types representing superfluid helium
!> \author hforbert
!> \date   2009-01-01
!> \par    History
!>         extracted helium_solvent_type from pint_types.F [lwalewski]
! **************************************************************************************************
MODULE helium_types

   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE input_constants,                 ONLY: helium_sampling_ceperley
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_para_env_type
   USE nnp_environment_types,           ONLY: nnp_type
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE splines_types,                   ONLY: spline_data_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_types'

   !> Energy contributions - symbolic names for indexing energy arrays
   INTEGER, PARAMETER, PUBLIC :: &
      e_id_total = 1, &
      e_id_potential = 2, &
      e_id_kinetic = 3, &
      e_id_interact = 4, &
      e_id_thermo = 5, &
      e_id_virial = 6

   !> Number of energy contributions for static array allocation
   INTEGER, PARAMETER, PUBLIC :: e_num_ids = 10

   !> number of density function identifiers
   INTEGER, PARAMETER, PUBLIC :: rho_num = 5

   !> density function identifier names
   INTEGER, PARAMETER, PUBLIC :: &
      rho_atom_number = 1, &
      rho_projected_area = 2, &
      rho_winding_number = 3, &
      rho_winding_cycle = 4, &
      rho_moment_of_inertia = 5

   !> derived data types
   PUBLIC :: helium_solvent_type
   PUBLIC :: helium_solvent_p_type
   PUBLIC :: int_arr_ptr

   !> functions
   PUBLIC :: helium_destroy_int_arr_ptr

! ***************************************************************************
!> \brief  Vector type useful for averaging
!> \author Lukasz Walewski
!> \date   2014-09-09
! ***************************************************************************
   TYPE helium_vector_type

      !> instantaneous value
      REAL(KIND=dp), DIMENSION(3)            :: inst = 0.0_dp

      !> accumulated value
      REAL(KIND=dp), DIMENSION(3)            :: accu = 0.0_dp

      !> running average
      REAL(KIND=dp), DIMENSION(3)            :: ravr = 0.0_dp

      !> restarted value
      REAL(KIND=dp), DIMENSION(3)            :: rstr = 0.0_dp

   END TYPE helium_vector_type

! ***************************************************************************
!> \brief data structure for solvent helium
!> \author hforbert
! ***************************************************************************
   TYPE helium_solvent_type

      TYPE(section_vals_type), POINTER  :: input => NULL()!< input data structure (the whole tree)
      TYPE(cp_logger_type), POINTER     :: logger => NULL()

      INTEGER       :: num_env = 0!< number of He environments in runtime

      INTEGER :: atoms = 0!< number of atoms
      INTEGER :: beads = 0!< number of beads per atom (needs to be an integer multiple of the solute's number of beads)
      INTEGER :: bead_ratio = 0!< ratio of helium beads to system beads
      REAL(KIND=dp) :: density = 0.0_dp !< helium density for free bulk in box

      ! some useful constants
      !
      REAL(KIND=dp) :: he_mass_au = 0.0_dp! mass of helium 4 in electron masses
      REAL(KIND=dp) :: hb2m = 0.0_dp!< hbar squared over m for 4He in CP2k units
      REAL(KIND=dp) :: tau = 0.0_dp!< 1/(k_B T p) with T - He temperature, p - number of beads
      REAL(KIND=dp) :: wpref = 0.0_dp!< prefactor for calculating superfluid fraction from <(M*W)^2>
      REAL(KIND=dp) :: apref = 0.0_dp!< prefactor for calculating superfluid fraction from <A^2/I_c>

      ! PBC related
      !
      LOGICAL                        :: periodic = .FALSE.!< true if bulk liquid helium in periodic box
      INTEGER                        :: cell_shape = 0!< unit cell shape for PBC calculations
      REAL(KIND=dp)                  :: cell_size = 0.0_dp!< size of the periodic box (helium only)
      REAL(KIND=dp)                   :: cell_size_inv = 0.0_dp!< 1/cell_size (inverse)
      REAL(KIND=dp), DIMENSION(3, 3)  :: cell_m = 0.0_dp!< the unit cell vectors' matrix
      REAL(KIND=dp), DIMENSION(3, 3)  :: cell_m_inv = 0.0_dp!< invrse  of the unit cell vectors' matrix
      REAL(KIND=dp), DIMENSION(3)    :: origin = 0.0_dp!< origin of the cell (first voxel position)
      REAL(KIND=dp)                  :: droplet_radius = 0.0_dp !< radius of the droplet

      REAL(KIND=dp), DIMENSION(3)    :: center = 0.0_dp!< COM of solute (if present) or center of periodic cell (if periodic) or COM of helium

      INTEGER :: sampling_method = helium_sampling_ceperley
      ! worm sampling parameters
      REAL(KIND=dp) :: worm_centroid_drmax = 0.0_dp
      INTEGER       :: worm_nstat = 0
      INTEGER       :: worm_staging_l = 0
      INTEGER       :: worm_repeat_crawl = 0
      INTEGER       :: worm_all_limit = 0
      INTEGER       :: worm_centroid_min = 0, worm_centroid_max = 0
      INTEGER       :: worm_staging_min = 0, worm_staging_max = 0
      INTEGER       :: worm_fcrawl_min = 0, worm_fcrawl_max = 0
      INTEGER       :: worm_bcrawl_min = 0, worm_bcrawl_max = 0
      INTEGER       :: worm_head_min = 0, worm_head_max = 0
      INTEGER       :: worm_tail_min = 0, worm_tail_max = 0
      INTEGER       :: worm_swap_min = 0, worm_swap_max = 0
      INTEGER       :: worm_open_close_min = 0, worm_open_close_max = 0
      INTEGER       :: worm_max_open_cycles = 0
      REAL(KIND=dp) :: worm_open_close_scale = 0.0_dp
      REAL(KIND=dp) :: worm_ln_openclose_scale = 0.0_dp
      LOGICAL       :: worm_allow_open = .FALSE., worm_show_statistics = .FALSE.

      ! worm specific variables
      REAL(KIND=dp), DIMENSION(3) :: worm_xtra_bead = 0.0_dp, worm_xtra_bead_work = 0.0_dp
      INTEGER :: worm_atom_idx = 0, worm_bead_idx = 0
      INTEGER :: worm_atom_idx_work = 0, worm_bead_idx_work = 0
      INTEGER :: iw = 0, it = 0
      LOGICAL :: worm_is_closed = .FALSE.!before isector=1 -> open; isector=0 -> closed

      INTEGER :: iter_norot = 0!< number of iterations to try for a given imaginary time slice rotation (num inner MC loop iters)
      INTEGER :: iter_rot = 0!< number of rotations to try (total number of iterations is iter_norot*iter_rot) (num outer MC loop iters)
      !
      INTEGER       :: maxcycle = 0!< maximum cyclic permutation change to attempt
      INTEGER       :: m_dist_type = 0!< distribution from which the cycle length m is sampled
      INTEGER       :: m_value = 0!< cycle length sampled with different probability than other lengths
      REAL(KIND=dp) :: m_ratio = 0.0_dp!< probability ratio betw m_value and other possible values of m
      !
      INTEGER :: relrot = 0!< relative rotation in imaginary time wrt normal system/starting configuration
      INTEGER :: bisection = 0 !< power of 2 number for bisection algorithm
      INTEGER :: bisctlog2 = 0!< log2(bisection)

      REAL(KIND=dp) :: e_corr = 0.0_dp !< potential correction energy due to finite box
      INTEGER :: pdx = 0!< pair density expansion max exponent

      ! MC step counters
      !
      INTEGER :: num_steps = 0!< number of iterations in the current run
      INTEGER :: first_step = 0!< first step, restarted from MOTION%PINT%ITERATION (default value: 0)
      INTEGER :: last_step = 0
      INTEGER :: current_step = 0 !< first_step + number of steps performed so far

      ! helium variables
      !
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: pos => NULL()!< position of the helium atoms DIM(3,atoms,beads)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: savepos => NULL()!< saved position of the helium atoms DIM(3,atoms,beads)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: work => NULL()!< same dimensions as pos
      !
      INTEGER, DIMENSION(:), POINTER :: permutation => NULL()!< current permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: savepermutation => NULL()!< saved permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: iperm => NULL()!< inverse of the current permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: saveiperm => NULL()!< saved inverse of the current permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: ptable => NULL()!< proposed cyclic permutation, DIM(max_cycle)
      INTEGER(KIND=int_8)              :: accepts = 0_int_8!< number of accepted new configurations
      !
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: tmatrix => NULL()!< ? permutation probability related
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: pmatrix => NULL()!< ? permutation probability related [use might change/new ones added/etc]
      REAL(KIND=dp) :: pweight = 0.0_dp!< ? permutation probability related
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: ipmatrix => NULL()
      INTEGER, DIMENSION(:, :), POINTER :: nmatrix => NULL()

      TYPE(spline_data_type), POINTER :: vij => NULL()!< physical pair potential energy
      TYPE(spline_data_type), POINTER :: u0 => NULL()!< pair density matrix coefficient (action) endpoint approx
      TYPE(spline_data_type), POINTER :: e0 => NULL()!< pair density matrix coefficient (energy) endpoint approx
      !< raw spline data for pair density matrix off diagonal expansion beyond endpoint approx:
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: uoffdiag => NULL()!< (action)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: eoffdiag => NULL()!< (energy)

      ! calculated properties
      !
      REAL(KIND=dp), DIMENSION(e_num_ids)    :: energy_inst = 0.0_dp!< energy contributions (instantaneous)
      REAL(KIND=dp), DIMENSION(e_num_ids)    :: energy_avrg = 0.0_dp!< energy contributions (averaged)
      TYPE(helium_vector_type)               :: wnumber = helium_vector_type()!< winding number
      TYPE(helium_vector_type)               :: wnmber2 = helium_vector_type()!< winding number squared
      TYPE(helium_vector_type)               :: proarea = helium_vector_type()!< projected area
      TYPE(helium_vector_type)               :: prarea2 = helium_vector_type()!< projected area squared
      TYPE(helium_vector_type)               :: mominer = helium_vector_type()!< moment of inertia
      INTEGER                                :: averages_iweight = 0!< weight for restarted averages
      LOGICAL                                :: averages_restarted = .FALSE.!< flag indicating whether the averages have been restarted

      REAL(KIND=dp) :: link_action = 0.0_dp, inter_action = 0.0_dp, pair_action = 0.0_dp

      !
      INTEGER                                :: rdf_nbin = 0!< number of bins for RDF
      INTEGER                                :: rdf_iweight = 0 !< weight for restarted RDF
      INTEGER                                :: rho_iweight = 0!< weight for restarted RHO
      INTEGER                                :: rdf_num = 0!< number of X-He-RDFs
      INTEGER                                :: rdf_num_ctr = 0 !< number of centers for RDF calc
      REAL(KIND=dp)                          :: rdf_delr = 0.0_dp!< delta r for RDF
      REAL(KIND=dp)                          :: rdf_maxr = 0.0_dp!< maximum r for RDF
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rdf_centers => NULL() !< positions of RDF solute  centers
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rdf_inst => NULL()!< RDF (instantaneous/tmp array)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rdf_rstr => NULL()!< RDF (restarted)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rdf_accu => NULL()!< RDF (accumulated for one run)
      LOGICAL :: rdf_present = .FALSE.
      LOGICAL :: rdf_sol_he = .FALSE.
      LOGICAL :: rdf_he_he = .FALSE.
      !
      INTEGER :: rho_nbin = 0
      INTEGER :: rho_num_act = 0!< actual number of density estimators
      INTEGER :: rho_num_min_len_wdg = 0!< number of optional estimators based on winding cycles
      INTEGER :: rho_num_min_len_non = 0!< number of optional estimators based on non-winding cycles
      INTEGER :: rho_num_min_len_all = 0!< number of optional estimators based on all cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_wdg_vals => NULL()!< minimum lengths of winding cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_non_vals => NULL()!< minimum lengths of non-winding cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_all_vals => NULL()!< minimum lengths of all cycles
      REAL(KIND=dp) :: rho_delr = 0.0_dp, rho_maxr = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_inst => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_rstr => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_accu => NULL()
      LOGICAL :: rho_present = .FALSE.
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER ::  rho_incr => NULL()!< increment for density bining

      TYPE(density_properties_type), DIMENSION(:), POINTER  :: rho_property => NULL()

      REAL(KIND=dp), DIMENSION(:, :), POINTER :: num_accepted => NULL()!< average number of accepted permutations of a given length
      !! on a given Levy level, plus one additional level which
      !! counts # of trials, REAL(BISCTLOG2+2, MAX_PERM_CYCLE)
      !! num_accepted(1,l) - # of trials for perm length l
      !! num_accepted(2,l) - # of selected perms of length l
      !! num_accepted(3,l) - # of perms of length l accepted at level 1
      !! average over He environments/processors
      REAL(KIND=dp), DIMENSION(:), POINTER :: plength_avrg => NULL()!< permutation length probability distribution DIM(atoms)
      REAL(KIND=dp), DIMENSION(:), POINTER :: plength_inst => NULL()!< instantaneous permutation length probability DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: atom_plength => NULL()!< length of the permutation cycle the atom belongs to DIM(atoms)

      TYPE(rng_stream_type), POINTER  :: rng_stream_uniform => NULL()!< random number stream with uniform distribution
      TYPE(rng_stream_type), POINTER  :: rng_stream_gaussian => NULL()!< random number stream with gaussian distribution

      ! variables related to solvated molecular system
      !
      LOGICAL :: solute_present = .FALSE.!< switch the interactions with the solute on or off
      INTEGER :: solute_atoms = 0!< number of solute atoms (pint_env%ndim/3)
      INTEGER :: solute_beads = 0!< number of solute beads (pint_env%p)
      INTEGER :: get_helium_forces = 0!< parameter to determine whether the average or last MC force should be taken to MD
      CHARACTER(LEN=2), DIMENSION(:), POINTER :: solute_element => NULL()!< element names of solute atoms (pint_env%ndim/3)
      TYPE(cell_type), POINTER  :: solute_cell => NULL()!< dimensions of the solvated system cell (a,b,c) (should be removed at some point)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: force_avrg => NULL()!< averaged forces exerted by He solvent on the solute DIM(p,ndim)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: force_inst => NULL()!< instantaneous forces exerted by He on the solute (p,ndim)
      CHARACTER(LEN=2), DIMENSION(:), POINTER  :: ename => NULL()
      INTEGER :: enum = 0
      INTEGER :: solute_interaction = 0

      LOGICAL :: interaction_pot_scan = .FALSE.!< whether to perform solute-helium interaction scan

      TYPE(nnp_type), POINTER :: nnp => NULL() !< neural network potential
      REAL(KIND=dp), DIMENSION(:), POINTER :: nnp_sr_cut => NULL() !< hard core cutoff in addition to the nnp

      ! temporary arrays for optimization
      !
      INTEGER, DIMENSION(:), POINTER         :: itmp_atoms_1d => NULL()!< DIM(atoms) - same as permutation
      INTEGER, DIMENSION(:), POINTER         :: itmp_atoms_np_1d => NULL()!< DIM(atoms*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_np_1d => NULL()!< DIM(3*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_p_ndim_1d => NULL()!< DIM(p*ndim)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_p_ndim_np_1d => NULL()!< DIM(p*ndim*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_atoms_beads_1d => NULL()!< DIM(3*atoms*beads)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_atoms_beads_np_1d => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rtmp_p_ndim_2d => NULL()!< DIM(p,ndim)
      LOGICAL, DIMENSION(:, :, :), POINTER     :: ltmp_3_atoms_beads_3d => NULL()!< DIM(3,atoms,beads) - same as pos
      LOGICAL, DIMENSION(:), POINTER         :: ltmp_atoms_1d => NULL()!< DIM(atoms) - for unpacking the permutation

   END TYPE helium_solvent_type

! ***************************************************************************
!> \brief data structure for array of solvent helium environments
!> \author cschran
! ***************************************************************************
   TYPE helium_solvent_p_type
      TYPE(helium_solvent_type), POINTER   :: helium => NULL()
      TYPE(mp_para_env_type), POINTER      :: comm => NULL()
      INTEGER, DIMENSION(:), POINTER       :: env_all => NULL()
   END TYPE helium_solvent_p_type

! ***************************************************************************
!> \brief  Container type for properties of a helium density function
!> \author Lukasz Walewski
!> \date   2014-09-09
! ***************************************************************************
   TYPE density_properties_type

      !> name of this density function
      CHARACTER(len=default_string_length) :: name = ""

      !> flag indicating whether this function should be calculated
      LOGICAL :: is_calculated = .FALSE.

      !> number of components that this function is composed of
      INTEGER :: num_components = 0

      !> suffixes for the filenames storing components of this function
      CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: filename_suffix => NULL()

      !> component names
      CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: component_name => NULL()

      !> indices locating the components of this function in the global density arrays
      INTEGER, DIMENSION(:), POINTER :: component_index => NULL()

   END TYPE density_properties_type

! ***************************************************************************
!> \brief  A pointer to an integer array, data type to be used in arrays of
!>         pointers.
!> \author Lukasz Walewski
!> \date   2013-12-11
! ***************************************************************************
   TYPE int_arr_ptr
      INTEGER, DIMENSION(:), POINTER :: iap => NULL()
   END TYPE int_arr_ptr

! ***************************************************************************
!> \brief  A pointer to a real array, data type to be used in arrays of
!>         pointers.
!> \author Lukasz Walewski
!> \date   2013-12-11
! ***************************************************************************
   TYPE real_arr_ptr
      REAL(KIND=dp), DIMENSION(:), POINTER :: rap => NULL()
   END TYPE real_arr_ptr

CONTAINS

! ***************************************************************************
!> \brief  Deallocate all arrays pointed to by the pointers stored in the
!>         integer pointer array
!> \param int_arr_p ...
!> \date   2013-12-12
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_destroy_int_arr_ptr(int_arr_p)

      TYPE(int_arr_ptr), DIMENSION(:), POINTER           :: int_arr_p

      INTEGER                                            :: ip

! deallocate memory used by each component of the pointer array

      DO ip = 1, SIZE(int_arr_p)
         IF (ASSOCIATED(int_arr_p(ip)%iap)) THEN
            DEALLOCATE (int_arr_p(ip)%iap)
         END IF
      END DO

      ! deallocate the memory used for pointer array
      IF (ASSOCIATED(int_arr_p)) THEN
         DEALLOCATE (int_arr_p)
      END IF

      RETURN
   END SUBROUTINE helium_destroy_int_arr_ptr

! ***************************************************************************
!> \brief  Deallocate all arrays pointed to by the pointers stored in the
!>         real pointer array
!> \param real_arr_p ...
!> \date   2013-12-12
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_destroy_real_arr_ptr(real_arr_p)

      TYPE(real_arr_ptr), DIMENSION(:), POINTER          :: real_arr_p

      INTEGER                                            :: ip

! do not attempt deallocation on null pointer

      IF (.NOT. ASSOCIATED(real_arr_p)) THEN
         RETURN
      END IF

      ! deallocate memory used by each component of the pointer array
      DO ip = 1, SIZE(real_arr_p)
         IF (ASSOCIATED(real_arr_p(ip)%rap)) THEN
            DEALLOCATE (real_arr_p(ip)%rap)
         END IF
      END DO

      ! deallocate the memory used for pointer array itself
      IF (ASSOCIATED(real_arr_p)) THEN
         DEALLOCATE (real_arr_p)
      END IF

      RETURN
   END SUBROUTINE helium_destroy_real_arr_ptr

END MODULE helium_types
