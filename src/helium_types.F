!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Data types representing superfluid helium
!> \author hforbert
!> \date   2009-01-01
!> \par    History
!>         extracted helium_solvent_type from pint_types.F [lwalewski]
! *****************************************************************************
MODULE helium_types

  USE cell_types,                      ONLY: cell_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE parallel_rng_types,              ONLY: rng_stream_type
  USE splines_types,                   ONLY: spline_data_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  REAL(kind=dp), PARAMETER, PUBLIC :: he_mass = 4.00263037059764_dp ! 4He mass in [u]

  ! Energy contributions - symbolic names for indexing energy arrays
  INTEGER, PARAMETER, PUBLIC :: e_id_total = 1, &
                                e_id_potential = 2, &
                                e_id_kinetic = 3, &
                                e_id_interact = 4, &
                                e_id_thermo = 5, &
                                e_id_virial = 6

  ! Number of energy contributions for static array allocation
  INTEGER, PARAMETER, PUBLIC :: e_num_ids = 10

  PUBLIC :: helium_solvent_type

! *****************************************************************************
!> \brief data of solvent helium
!> \note
!>    ref_count     - reference count of this data structure
!>    id_nr         - identification number of this data structure
!>    input         - input data structure (the whole tree)
!>    num_env       - number of independent He environments
!>    num_env_restart - number of He environments present in the restart file
!>    periodic      - true if bulk liquid helium in periodic box
!>    cell_shape    - unit cell shape for PBC calculations
!>    cell_size     - size of the periodic box for helium only
!>    cell_size_inv - 1 / cell_size (inverse)
!>    cell_m_inv    - invrse  of the unit cell vectors' matrix
!>    tau           - 1/(k_B T p) with T being the He temperature, p - number of beads
!>    density       - helium density for free bulk in box
!>    e_corr        - potential correction energy due to finite box
!>    hb2m          - hbar squared over m for helium
!>    atoms         - number of helium atoms
!>    beads         - number of helium beads per atom (needs to be an integer
!>                      multiple of the solvated molecule bead number)
!>    pdx           - pair density expansion max exponent
!>    maxcycle      - maximum cyclic permutation change to attempt
!>    m_value       - cycle length sampled with different probability than other lengths
!>    m_ratio       - probability ratio betw m_value and other possible values of m
!>    bisection     - power of 2 number for bisection algorithm
!>    bisctlog2     - log2(bisection)
!>    relrot        - relative rotation in imaginary time with
!>                    respect to the normal system/starting configuration
!>    iter_norot    - number of iterations to try for a given
!>                    imaginary time slice rotation
!>    iter_rot      - number of rotations to try
!>                    total number of iterations is iter_norot*iter_rot
!>    iter          - number of iterations tried so far
!>    first_step    - first step, restarted from MOTION%PINT%ITERATION (def val =0)
!>    current_step  - first_step + number of steps performed so far
!>    accepts       - number of accepted new configurations
!>    num_accepted  - average number of accepted permutations of a given length
!>                    on a given Levy level, plus one additional level which
!>                    counts # of trials, REAL(BISCTLOG2+2, MAX_PERM_CYCLE)
!>                    num_accepted(1,l) - # of trials for perm length l
!>                    num_accepted(2,l) - # of selected perms of length l
!>                    num_accepted(3,l) - # of perms of length l accepted at level 1
!>                    average over He environments/processors
!>    wnumber_inst  - winding number (instantaneous) REAL(3)
!>    wnumber_avrg  - winding number (averaged) REAL(3)
!>    sdensity_inst - superfluid density (instantaneous) REAL
!>    sdensity_avrg - superfluid density (averaged) REAL
!>    pos           - position of the helium atoms REAL(3,ATOMS,BEADS)
!>    work          - same dimensions as pos
!>    tmatrix       - ? permutation probability related
!>    pmatrix       - ? permutation probability related
!>                    [use might change/new ones added/etc]
!>    pweight       - ? permutation probability related
!>    ptable        - proposed cyclic permutation
!>                    (dimension max_cycle)
!>    permutation   - current permutation state INTEGER(ATOMS)
!>    iperm         - inverse of the current permutation state INTEGER(ATOMS)
!>    plength_avrg  - permutation length probability distribution REAL(ATOMS)
!>    plength_inst  - instantaneous permutation length probability REAL(ATOMS)
!>    atom_plength  - length of the permutation cycle the atom belongs to INTEGER(ATOMS)
!>    uij           - pair density matrix coefficients (action)
!>    eij           - pair density matrix coefficients (energy)
!>    bead_ratio    - ratio of helium beads to system beads
!>    rng_stream_uniform - random numbers from uniform distribution
!>    rng_stream_gaussian - random numbers from gaussian distribution
!>    rdf_delr      - delta r for RDF
!>    rdf_maxr      - maximum r for RDF
!>    solute_present- switch the interactions with the solute on or off
!>    solute_atoms  - number of solute atoms (=pint_env%ndim/3)
!>    solute_beads  - number of solute beads (=pint_env%p)
!>    solute_element- element names of solute atoms (NDIM/3)
!>    solute_number - number of solute atoms of different atomic kinds
!>    solute_index  - indices of solute atoms sorted by atomic kinds
!>    solute_cell   - dimensions of the solvated system cell (a,b,c)
!>                    (not needed now and should be removed at some point)
!>    force_avrg    - averaged forces exerted by He solvent on the solute (P,NDIM)
!>    force_inst    - instantaneous forces exerted by He on the solute (P, NDIM)
!>    General purpose temporary arrays
!>    (performance measure: no allocation in frequently called subroutines)
!>    rtmp_3_np_1d - real, 1D temp array (3 * NUM_ENV)
!>      rtmp_p_ndim_1d - real, 1D temp array (P * NDIM)
!>    rtmp_p_ndim_np_1d - real, 1D temp array (P * NDIM * NUM_ENV)
!>      rtmp_p_ndim_2d - real, 2D temp array (P, NDIM)
!>      rtmp_3_atoms_beads_1d - real, 1D temp array (3 * ATOMS * BEADS)
!>      ltmp_3_atoms_beads_3d - logical, 3D temp array (3,ATOMS,BEADS) - same as pos
!>      itmp_atoms_1d    - integer, 1D temp array (ATOMS) - same as permutation
!>    itmp_atoms_np_1d - integer, 1D temp array (ATOMS*NUM_ENV)
!>      ltmp_atoms_1d    - logical, 1D temp array (ATOMS) - for unpacking permutation
!> \author hforbert
! *****************************************************************************
  TYPE helium_solvent_type
     INTEGER :: id_nr,ref_count
     TYPE(section_vals_type), POINTER :: input
     TYPE(cp_logger_type), POINTER :: logger
     LOGICAL :: periodic
     INTEGER :: cell_shape, num_env, num_env_restart, m_value
     REAL(kind=dp) :: cell_size, cell_size_inv, m_ratio
     REAL(kind=dp), DIMENSION(3,3) ::  cell_m_inv
     REAL(kind=dp) :: tau, density, e_corr, hb2m, pweight
     INTEGER :: atoms, beads, pdx, maxcycle, relrot, iter_norot, iter_rot
     INTEGER :: bisection, bisctlog2
     INTEGER :: bead_ratio, num_steps, first_step, last_step, current_step
     INTEGER(kind=int_8) :: accepts
     INTEGER, DIMENSION(:), POINTER :: ptable
     INTEGER, DIMENSION(:), POINTER :: permutation
     INTEGER, DIMENSION(:), POINTER :: iperm
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: pos
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: work
     REAL(kind=dp), DIMENSION(:,:), POINTER :: tmatrix
     REAL(kind=dp), DIMENSION(:,:), POINTER :: pmatrix
     REAL(kind=dp), DIMENSION(:,:), POINTER :: ipmatrix
     INTEGER, DIMENSION(:,:), POINTER :: nmatrix
     TYPE (spline_data_p_type), DIMENSION(:,:), POINTER :: uij
     TYPE (spline_data_p_type), DIMENSION(:,:), POINTER :: eij
     REAL(kind=dp), DIMENSION(e_num_ids) :: energy_inst, energy_avrg
     REAL(kind=dp), DIMENSION(3)         :: wnumber_inst, wnumber_avrg, origin
     REAL(kind=dp) :: sdensity_inst, sdensity_avrg
     REAL(kind=dp), DIMENSION(:,:), POINTER :: num_accepted
     REAL(kind=dp), DIMENSION(:), POINTER :: plength_avrg, plength_inst
     INTEGER, DIMENSION(:), POINTER :: atom_plength
     INTEGER :: rdf_nbin
     REAL(kind=dp) :: rdf_delr, rdf_maxr
     REAL(kind=dp), DIMENSION(:), POINTER :: rdf_inst, rdf_avrg
     INTEGER :: rho_nbin, rho_minb, rho_iweight, rho_num
     REAL(kind=dp) :: rho_delr, rho_maxr
     REAL(kind=dp), DIMENSION(:,:,:,:), POINTER :: rho_inst, rho_avrg, rho_rstr
     LOGICAL :: rho_restart, rho_present
     TYPE(rng_stream_type), POINTER :: rng_stream_gaussian, rng_stream_uniform
     ! variables that describe the solvated molecular system
     LOGICAL :: solute_present
     INTEGER :: solute_atoms, solute_beads
     CHARACTER(LEN=2), DIMENSION(:), POINTER :: solute_element
     INTEGER, DIMENSION(:), POINTER :: solute_number
     INTEGER, DIMENSION(:,:), POINTER :: solute_index
     TYPE(cell_type), POINTER  :: solute_cell
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: force_avrg, force_inst
     INTEGER, DIMENSION(:), POINTER         :: itmp_atoms_1d, itmp_atoms_np_1d
     REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_np_1d, rtmp_p_ndim_1d, rtmp_p_ndim_np_1d
     REAL(KIND=dp), DIMENSION(:), POINTER   :: rtmp_3_atoms_beads_1d, rtmp_3_atoms_beads_np_1d
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: rtmp_p_ndim_2d
     LOGICAL, DIMENSION(:,:,:), POINTER :: ltmp_3_atoms_beads_3d
     LOGICAL, DIMENSION(:), POINTER :: ltmp_atoms_1d
  END TYPE helium_solvent_type

END MODULE helium_types
