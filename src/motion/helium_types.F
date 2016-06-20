!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE splines_types,                   ONLY: spline_data_p_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_types'

   REAL(KIND=dp), PARAMETER, PUBLIC :: he_mass = 4.00263037059764_dp !< 4He mass in [u]

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

   !> number of atom type identifiers
   INTEGER, PARAMETER, PUBLIC :: hid_num = 4

   !> atom type identifier names (for helium-solute interactions)
   INTEGER, PARAMETER, PUBLIC :: &
      hid_chlorine = 1, &
      hid_oxygen = 2, &
      hid_hydrogen = 3, &
      hid_carbon = 4

   !> derived data types
   PUBLIC :: helium_solvent_type
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
      REAL(KIND=dp), DIMENSION(3)            :: inst

      !> accumulated value
      REAL(KIND=dp), DIMENSION(3)            :: accu

      !> running average
      REAL(KIND=dp), DIMENSION(3)            :: ravr

      !> restarted value
      REAL(KIND=dp), DIMENSION(3)            :: rstr

   END TYPE helium_vector_type

! ***************************************************************************
!> \brief data structure for solvent helium
!> \author hforbert
! ***************************************************************************
   TYPE helium_solvent_type

      INTEGER                           :: id_nr !< identification number of this data structure
      INTEGER                           :: ref_count !< reference count of this data structure
      TYPE(section_vals_type), POINTER  :: input !< input data structure (the whole tree)
      TYPE(cp_logger_type), POINTER     :: logger

      INTEGER       :: num_env !< number of He environments in runtime
      INTEGER       :: num_env_restart !< number of He environments in restart file

      INTEGER :: atoms !< number of atoms
      INTEGER :: beads !< number of beads per atom (needs to be an integer multiple of the solute's number of beads)
      INTEGER :: bead_ratio !< ratio of helium beads to system beads
      REAL(KIND=dp) :: density !< helium density for free bulk in box

      ! some useful constants
      !
      REAL(KIND=dp) :: hb2m !< hbar squared over m for 4He in CP2k units
      REAL(KIND=dp) :: tau !< 1/(k_B T p) with T - He temperature, p - number of beads
      REAL(KIND=dp) :: wpref !< prefactor for calculating superfluid fraction from <(M*W)^2>
      REAL(KIND=dp) :: apref !< prefactor for calculating superfluid fraction from <A^2/I_c>

      ! PBC related
      !
      LOGICAL                        :: periodic !< true if bulk liquid helium in periodic box
      INTEGER                        :: cell_shape !< unit cell shape for PBC calculations
      REAL(KIND=dp)                  :: cell_size !< size of the periodic box (helium only)
      REAL(KIND=dp)                  :: cell_size_inv !< 1/cell_size (inverse)
      REAL(KIND=dp), DIMENSION(3, 3)  :: cell_m !< the unit cell vectors' matrix
      REAL(KIND=dp), DIMENSION(3, 3)  :: cell_m_inv !< invrse  of the unit cell vectors' matrix
      REAL(KIND=dp), DIMENSION(3)    :: origin !< origin of the cell (first voxel position)
      REAL(KIND=dp)                  :: droplet_radius !< radius of the droplet

      REAL(KIND=dp), DIMENSION(3)    :: center !< COM of solute (if present) or center of periodic cell (if periodic) or COM of helium

      ! sampling parameters
      INTEGER :: sampling_method

      INTEGER :: iter_norot !< number of iterations to try for a given imaginary time slice rotation (num inner MC loop iters)
      INTEGER :: iter_rot !< number of rotations to try (total number of iterations is iter_norot*iter_rot) (num outer MC loop iters)
      !
      INTEGER       :: maxcycle !< maximum cyclic permutation change to attempt
      INTEGER       :: m_dist_type !< distribution from which the cycle length m is sampled
      INTEGER       :: m_value !< cycle length sampled with different probability than other lengths
      REAL(KIND=dp) :: m_ratio !< probability ratio betw m_value and other possible values of m
      !
      INTEGER :: relrot !< relative rotation in imaginary time wrt normal system/starting configuration
      INTEGER :: bisection !< power of 2 number for bisection algorithm
      INTEGER :: bisctlog2 !< log2(bisection)

      REAL(KIND=dp) :: e_corr !< potential correction energy due to finite box
      INTEGER :: pdx !< pair density expansion max exponent

      ! MC step counters
      !
      INTEGER :: num_steps !< number of iterations in the current run
      INTEGER :: first_step !< first step, restarted from MOTION%PINT%ITERATION (default value: 0)
      INTEGER :: last_step
      INTEGER :: current_step !< first_step + number of steps performed so far

      ! helium variables
      !
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: pos !< position of the helium atoms DIM(3,atoms,beads)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: work !< same dimensions as pos
      !
      INTEGER, DIMENSION(:), POINTER :: permutation !< current permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: iperm !< inverse of the current permutation state DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: ptable !< proposed cyclic permutation, DIM(max_cycle)
      INTEGER(KIND=int_8)              :: accepts !< number of accepted new configurations
      !
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: tmatrix !< ? permutation probability related
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: pmatrix !< ? permutation probability related [use might change/new ones added/etc]
      REAL(KIND=dp) :: pweight !< ? permutation probability related
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: ipmatrix
      INTEGER, DIMENSION(:, :), POINTER :: nmatrix

      TYPE(spline_data_p_type), DIMENSION(:, :), POINTER :: uij !< pair density matrix coefficients (action)
      TYPE(spline_data_p_type), DIMENSION(:, :), POINTER :: eij !< pair density matrix coefficients (energy)

      ! calculated properties
      !
      REAL(KIND=dp), DIMENSION(e_num_ids)    :: energy_inst !< energy contributions (instantaneous)
      REAL(KIND=dp), DIMENSION(e_num_ids)    :: energy_avrg !< energy contributions (averaged)
      TYPE(helium_vector_type)               :: wnumber !< winding number
      TYPE(helium_vector_type)               :: wnmber2 !< winding number squared
      TYPE(helium_vector_type)               :: proarea !< projected area
      TYPE(helium_vector_type)               :: prarea2 !< projected area squared
      TYPE(helium_vector_type)               :: mominer !< moment of inertia
      INTEGER                                :: averages_iweight !< weight for restarted averages
      LOGICAL                                :: averages_restarted !< flag indicating whether the averages have been restarted

      REAL(KIND=dp) :: link_action, inter_action, pair_action

      !
      INTEGER                                :: rdf_nbin !< number of bins for RDF
      INTEGER                                :: rdf_iweight !< weight for restarted RDF
      INTEGER                                :: rho_iweight !< weight for restarted RHO
      INTEGER                                :: rdf_num !< number of RDF estimators
      INTEGER                                :: rdf_num_ctr !< number of centers for RDF calc
      REAL(KIND=dp)                          :: rdf_delr !< delta r for RDF
      REAL(KIND=dp)                          :: rdf_maxr !< maximum r for RDF
      REAL(KIND=dp), DIMENSION(:), POINTER :: rdf_centers !< positions of RDF centers
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rdf_inst !< RDF (instantaneous/tmp array)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rdf_rstr !< RDF (restarted)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rdf_accu !< RDF (accumulated for one run)
      LOGICAL :: rdf_present
      !
      INTEGER :: rho_nbin
      INTEGER :: rho_num_act !< actual number of density estimators
      INTEGER :: rho_num_min_len_wdg !< number of optional estimators based on winding cycles
      INTEGER :: rho_num_min_len_non !< number of optional estimators based on non-winding cycles
      INTEGER :: rho_num_min_len_all !< number of optional estimators based on all cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_wdg_vals !< minimum lengths of winding cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_non_vals !< minimum lengths of non-winding cycles
      INTEGER, DIMENSION(:), POINTER :: rho_min_len_all_vals !< minimum lengths of all cycles
      REAL(KIND=dp) :: rho_delr, rho_maxr
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_inst
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_rstr
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER :: rho_accu
      LOGICAL :: rho_present
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER ::  rho_incr !< increment for density bining

      TYPE(density_properties_type), DIMENSION(:), POINTER  :: rho_property

      REAL(KIND=dp), DIMENSION(:, :), POINTER :: num_accepted !< average number of accepted permutations of a given length
      !! on a given Levy level, plus one additional level which
      !! counts # of trials, REAL(BISCTLOG2+2, MAX_PERM_CYCLE)
      !! num_accepted(1,l) - # of trials for perm length l
      !! num_accepted(2,l) - # of selected perms of length l
      !! num_accepted(3,l) - # of perms of length l accepted at level 1
      !! average over He environments/processors
      REAL(KIND=dp), DIMENSION(:), POINTER :: plength_avrg !< permutation length probability distribution DIM(atoms)
      REAL(KIND=dp), DIMENSION(:), POINTER :: plength_inst !< instantaneous permutation length probability DIM(atoms)
      INTEGER, DIMENSION(:), POINTER :: atom_plength !< length of the permutation cycle the atom belongs to DIM(atoms)

      TYPE(rng_stream_type), POINTER :: rng_stream_uniform !< random number stream with uniform distribution
      TYPE(rng_stream_type), POINTER :: rng_stream_gaussian !< random number stream with gaussian distribution

      ! variables related to solvated molecular system
      !
      LOGICAL :: solute_present !< switch the interactions with the solute on or off
      INTEGER :: solute_atoms !< number of solute atoms (pint_env%ndim/3)
      INTEGER :: solute_beads !< number of solute beads (pint_env%p)
      INTEGER :: get_helium_forces !< parameter to determine whether the average or last MC force should be taken to MD
      CHARACTER(LEN=2), DIMENSION(:), POINTER :: solute_element !< element names of solute atoms (pint_env%ndim/3)
      TYPE(int_arr_ptr), DIMENSION(:), POINTER :: solute_i !< indices of solute atoms sorted by atomic kinds
      TYPE(cell_type), POINTER  :: solute_cell !< dimensions of the solvated system cell (a,b,c) (should be removed at some point)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: force_avrg !< averaged forces exerted by He solvent on the solute DIM(p,ndim)
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: force_inst !< instantaneous forces exerted by He on the solute (p,ndim)
      CHARACTER(LEN=2), DIMENSION(:), POINTER  :: ename
      INTEGER, DIMENSION(:), POINTER           :: eid
      INTEGER :: enum
      INTEGER :: solute_interaction

      LOGICAL :: interaction_pot_scan !< wheather to perform solute-helium interaction scan

      ! temporary arrays for optimization
      !
      INTEGER, DIMENSION(:), POINTER         :: itmp_atoms_1d !< DIM(atoms) - same as permutation
      INTEGER, DIMENSION(:), POINTER         :: itmp_atoms_np_1d !< DIM(atoms*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_np_1d !< DIM(3*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_p_ndim_1d !< DIM(p*ndim)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_p_ndim_np_1d !< DIM(p*ndim*num_env)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_atoms_beads_1d !< DIM(3*atoms*beads)
      REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_atoms_beads_np_1d
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: rtmp_p_ndim_2d !< DIM(p,ndim)
      LOGICAL, DIMENSION(:, :, :), POINTER     :: ltmp_3_atoms_beads_3d !< DIM(3,atoms,beads) - same as pos
      LOGICAL, DIMENSION(:), POINTER         :: ltmp_atoms_1d !< DIM(atoms) - for unpacking the permutation

   END TYPE helium_solvent_type

! ***************************************************************************
!> \brief  Container type for properties of a helium density function
!> \author Lukasz Walewski
!> \date   2014-09-09
! ***************************************************************************
   TYPE density_properties_type

      !> name of this density function
      CHARACTER(len=default_string_length) :: name

      !> flag indicating whether this function should be calculated
      LOGICAL :: is_calculated

      !> number of components that this function is composed of
      INTEGER :: num_components

      !> suffixes for the filenames storing components of this function
      CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: filename_suffix

      !> component names
      CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: component_name

      !> indices locating the components of this function in the global density arrays
      INTEGER, DIMENSION(:), POINTER :: component_index

   END TYPE density_properties_type

! ***************************************************************************
!> \brief  A pointer to an integer array, data type to be used in arrays of
!>         pointers.
!> \author Lukasz Walewski
!> \date   2013-12-11
! ***************************************************************************
   TYPE int_arr_ptr
      INTEGER, DIMENSION(:), POINTER :: iap
   END TYPE int_arr_ptr

! ***************************************************************************
!> \brief  A pointer to a real array, data type to be used in arrays of
!>         pointers.
!> \author Lukasz Walewski
!> \date   2013-12-11
! ***************************************************************************
   TYPE real_arr_ptr
      REAL(KIND=dp), DIMENSION(:), POINTER :: rap
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

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_destroy_int_arr_ptr', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_destroy_real_arr_ptr', &
         routineP = moduleN//':'//routineN

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
