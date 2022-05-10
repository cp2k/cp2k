!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines CDFT control structures
!> \par   History
!>                 separated from cp_control_types [03.2017]
!> \author Nico Holmberg [03.2017]
! **************************************************************************************************
MODULE qs_cdft_types
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_fm_types,                     ONLY: cp_fm_p_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE hirshfeld_types,                 ONLY: hirshfeld_type,&
                                              release_hirshfeld_type
   USE input_constants,                 ONLY: becke_cutoff_global,&
                                              outer_scf_none,&
                                              radius_single,&
                                              shape_function_gaussian
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE outer_scf_control_types,         ONLY: outer_scf_control_type,&
                                              qs_outer_scf_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_cdft_opt_types,               ONLY: cdft_opt_type_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief some parameters useful for becke_constraints
!> \param aij             pairwise parameters used to adjust the Becke cell boundaries built from atomic radii
!> \param adjust          logical which determines if the Becke potential is adjusted with atomic radii
!> \param cavity          the Gaussian confinement cavity: the constraint is nonzero outside this cavity
!> \param cavity_confine  logical which determines if cavity confinement is active
!> \param cavity_mat      a compacted version of cavity
!> \param cavity_shape    the confinement cavity shape id
!> \param cavity_env      the structure used to build the Gaussian cavity
!> \param confine_bounds  grid point indices outside which the constraint vanishes along Z-axis
!> \param cutoff_type     the cutoff type to use for building the constraint
!> \param cutoffs         element specific cutoffs
!> \param cutoffs_tmp     same as cutoffs but a temporary read during parsing of this type
!> \param eps_cavity      threshold used screen small values of the Gaussian cavity density
!> \param in_memory       logical which determines if the gradients of the Becke potential should be
!> \param print_cavity    logical to print the Gaussian confinement cavity
!> \param radii           permanent copy of radii_tmp
!> \param radii_tmp       temporary list of element specific atomic radii used to adjust the Becke cells
!> \param rcavity         an optional global radius parameter used to define the Gaussian confinement cavity
!> \param rglobal         global cutoff to use for building the constraint
!>                        computed simultaneously with the potential instead of separately
!> \param should_skip     logical which determines is grid points should be skipped if all constraint
!>                        atoms are found to reside beyond the cutoff distance from it
!> \param use_bohr        decides whether to use angstrom or bohr units for the confinement cavity radius
! **************************************************************************************************
   ! Utility vector container for building becke constraint
   TYPE becke_vector_buffer
      LOGICAL                              :: store_vectors
      REAL(kind=dp), ALLOCATABLE, &
         DIMENSION(:)                      :: distances
      REAL(kind=dp), ALLOCATABLE, &
         DIMENSION(:, :)                   :: distance_vecs, &
                                              position_vecs, &
                                              R12
      REAL(kind=dp), ALLOCATABLE, &
         DIMENSION(:, :, :)                :: pair_dist_vecs
   END TYPE becke_vector_buffer

   TYPE becke_constraint_type
      INTEGER                              :: cavity_shape, cutoff_type, &
                                              confine_bounds(2)
      LOGICAL                              :: in_memory, &
                                              adjust, cavity_confine, &
                                              should_skip, print_cavity, &
                                              use_bohr
      REAL(KIND=dp)                        :: rglobal, &
                                              rcavity, eps_cavity
      REAL(KIND=dp), DIMENSION(:), POINTER :: cutoffs, cutoffs_tmp, &
                                              radii_tmp, radii
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                   :: aij
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                :: cavity_mat
      TYPE(becke_vector_buffer)            :: vector_buffer
      TYPE(hirshfeld_type), POINTER        :: cavity_env
      TYPE(pw_p_type)                      :: cavity
   END TYPE becke_constraint_type

! **************************************************************************************************
! \brief control parameters for Hirshfeld constraints
!> \param gaussian_shape  the type of Gaussian to use (shape_function Gaussian)
!> \param radii           list of Gaussian radii for different atomic kinds
!> \param radius          Gaussian radius parameter
!> \param shape_function  the constraint type: atomic density or single Gaussian
!> \param use_bohr        determines whether to use angstrom or bohr units for the radii of Gaussians
!> \param use_atomic_cutoff        Boolean to control atomic_cutoff
!> \param atomic_cutoff        cutoff for Hirshfeld constraint atomic density
!> \param atoms_memory        atoms to store in memory
!> \param eps_cutoff       cutoff for Hirshfeld constraint weight function
!> \param print_density    logical to decide if Hirshfeld densities should be printed
!> \param hirshfeld_env   auxiliary type storing information about the Gaussians
! **************************************************************************************************
   TYPE hirshfeld_constraint_type
      INTEGER                              :: gaussian_shape, shape_function, atoms_memory
      LOGICAL                              :: use_bohr, print_density, use_atomic_cutoff
      REAL(KIND=dp)                        :: radius, eps_cutoff, atomic_cutoff
      REAL(KIND=dp), DIMENSION(:), POINTER :: radii
      TYPE(hirshfeld_type), POINTER        :: hirshfeld_env
   END TYPE hirshfeld_constraint_type

! **************************************************************************************************
!> \brief control parameters for CDFT simulations
!> \param fragment_a_fname      filename of cube file holding the total electron density
!>                              of isolated fragment a
!> \param fragment_b_fname      filename of cube file holding the total electron density
!>                              of isolated fragment b
!> \param fragment_a_spin_fname filename of cube file holding the spin difference density
!>                              of isolated fragment a
!> \param fragment_b_spin_fname filename of cube file holding the spin difference density
!>                              of isolated fragment b
!> \param ref_count             the ref count
!> \param need_pot              logical which determines if the Becke potential needs to be built
!> \param save_pot              logical which determines if the Becke potential should be saved until forces
!>                              have been evaluated
!> \param atomic_charges        flag that determines if atomic CDFT charges should be computed
!> \param total_steps           counter to keep track of the total number of SCF steps
!> \param type                  the type of CDFT constraint to use
!> \param precond_freq          preconditioner can be used if SCF converged in less than precond_freq steps
!> \param nreused               determines how many times the current OT preconditioner has been reused
!> \param max_reuse             the same preconditioner can be used a maximum of max_reuse times
!> \param purge_freq            determines how large nbad_conv can grow before purging the wfn/constraint history
!> \param nbad_conv             a running counter keeping track of the number of CDFT SCF loops when the first
!>                              CDFT SCF iteration required more than 1 outer SCF loop. Reset when convergence is
!>                              smooth
!> \param purge_offset          purging is only allowed when more than purge_offset steps have passed since
!>                              last purge
!> \param istep                 a counter to keep track of how many steps have passed since the last purge
!> \param ienergy               a counter tracking the total number of CDFT energy evaluations
!> \param natoms                the total number of atoms included in constraint/dummy atom groups
!> \param atoms                 list of constraint atoms
!> \param need_pot              logical which determines if the constraint potential needs to be built
!> \param save_pot              logical which determines if the constraint potential should be saved until forces
!>                              have been evaluated
!> \param do_et                 logical which determines if a ET coupling calculation was requested
!> \param reuse_precond         logical which determines if a preconditioner can be reused
!> \param purge_history         logical which determines if the wfn/constraint history can be purged
!> \param should_purge          logical which determines if purging should take place after this CDFT SCF loop
!> \param calculate_metric      logical which determines if the ET coupling reliability metric is computed
!> \param fragment_density      use isolated fragment densities as a reference for the constraint
!> \param fragments_integrated  logical to determine if the fragment densities have been integrated
!> \param flip_fragment         should the spin difference density of the either fragment be flipped
!> \param transfer_pot          logical which determines if constraint should be saved for reuse later
!> \param external_control      logical which determines if the constraint has already been built
!>                              in a mixed_env that holds multiple CDFT states
!> \param first_iteration       a flag to mark the first iteration for printing of additional data
!> \param print_weight          logical which determines if CDFT weight functions should be saved to a file
!> \param in_memory       logical which determines if the gradients of the Becke potential should be
!> \param is_constraint         list of logicals which determines if an atom is included in a constraint group
!> \param strength              Lagrangian multipliers of the constraints
!> \param target                target values of the constraints
!> \param value                 integrated values of the constraints
!> \param charges_fragment      atomic partial charges computed from the isolated fragment densities
!> \param becke_control         control parameters for Becke constraints
!> \param group                 container for atom groups each defining their own constraint
!> \param occupations           occupation numbers in case non-uniform MO occupation (for do_et)
!> \param mo_coeff              save the MO coeffs (for do_et)
!> \param matrix_s              save the overlap matrix (for do_et)
!> \param wmat                  matrix representation of the weight function (for do_et)
!> \param matrix_p              save the density matrix (for calculate_metric)
!> \param hirshfeld_control     control parameters for Hirshfeld constraints
!> \param constraint_control    the outer_scf_control_type for the CDFT constraints
!> \param ot_control            the outer_scf_control_type for OT where data is stashed when outside the OT
!>                              outer loop
!> \param charge                atomic CDFT real space potentials needed to calculate CDFT charges
!> \param fragments             container for isolated fragment densities read from cube files
!> \param constraint            holds information about the CDFT SCF loop
! **************************************************************************************************
   ! To build multiple constraints
   TYPE cdft_group_type
      ! Atoms of this constraint group
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      ! Constraint type: charge constraint, magnetization density constraint, or spin channel specific constraint
      INTEGER                              :: constraint_type
      ! Is the constraint fragment based
      LOGICAL                              :: is_fragment_constraint
      ! Temporary array holding a component of the weight function gradient that only includes
      ! terms defined on constraint atoms
      REAL(kind=dp), ALLOCATABLE, &
         DIMENSION(:, :)                   :: d_sum_const_dR
      ! Coefficients that determine how to sum up the atoms to form the constraint
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff
      ! Result of integration dw/dR * rho_r dr where dw/dR is the weight function gradient
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                   :: integrated
      ! Atomic gradients of the weight function at every grid point
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)             :: gradients
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)             :: gradients_x
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)             :: gradients_y
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)             :: gradients_z
      ! The weight function of this constraint group
      TYPE(pw_p_type)                      :: weight
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: hw_rho_atomic
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: hw_rho_atomic_dr
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: hw_rho_atomic_charge
      TYPE(pw_p_type)             :: hw_rho_total_constraint
      TYPE(pw_p_type)          :: hw_rho_total
   END TYPE cdft_group_type

   TYPE cdft_control_type
      CHARACTER(LEN=default_path_length)   :: fragment_a_fname, &
                                              fragment_b_fname, &
                                              fragment_a_spin_fname, &
                                              fragment_b_spin_fname
      INTEGER                              :: ref_count, total_steps, TYPE, &
                                              precond_freq, nreused, max_reuse, &
                                              purge_freq, nbad_conv, purge_offset, &
                                              istep, ienergy, natoms
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      LOGICAL                              :: need_pot, save_pot, do_et, &
                                              reuse_precond, purge_history, &
                                              should_purge, calculate_metric, &
                                              atomic_charges, fragment_density, &
                                              fragments_integrated, flip_fragment(2), &
                                              transfer_pot, external_control, &
                                              first_iteration, print_weight, in_memory
      LOGICAL, POINTER, DIMENSION(:)       :: is_constraint
      REAL(KIND=dp), DIMENSION(:), POINTER :: strength, TARGET, value
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                   :: charges_fragment
      TYPE(becke_constraint_type), POINTER :: becke_control
      TYPE(cdft_group_type), POINTER, &
         DIMENSION(:)                      :: group
      TYPE(cp_1d_r_p_type), ALLOCATABLE, &
         DIMENSION(:)                      :: occupations
      TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                           :: mo_coeff
      TYPE(dbcsr_p_type)                   :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER                           :: wmat, matrix_p
      TYPE(hirshfeld_constraint_type), &
         POINTER                           :: hirshfeld_control
      TYPE(outer_scf_control_type)         :: constraint_control, ot_control
      TYPE(pw_p_type), POINTER, &
         DIMENSION(:)                      :: charge
      TYPE(pw_p_type), POINTER, &
         DIMENSION(:, :)                   :: fragments
      TYPE(qs_outer_scf_type)              :: constraint
   END TYPE cdft_control_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_types'

   ! Public data types

   PUBLIC :: becke_constraint_type, &
             cdft_control_type, &
             cdft_group_type, &
             hirshfeld_constraint_type

   ! Public subroutines

   PUBLIC :: cdft_control_create, &
             cdft_control_release

CONTAINS

! **************************************************************************************************
!> \brief create the becke_constraint_type
!> \param becke_control the structure to create
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_create(becke_control)
      TYPE(becke_constraint_type), POINTER               :: becke_control

      CPASSERT(.NOT. ASSOCIATED(becke_control))
      ALLOCATE (becke_control)

      becke_control%adjust = .FALSE.
      becke_control%cutoff_type = becke_cutoff_global
      becke_control%cavity_confine = .FALSE.
      becke_control%should_skip = .FALSE.
      becke_control%print_cavity = .FALSE.
      becke_control%in_memory = .FALSE.
      becke_control%use_bohr = .FALSE.
      becke_control%confine_bounds = 0
      becke_control%rcavity = 3.0_dp
      becke_control%rglobal = 6.0_dp
      becke_control%eps_cavity = 1.0e-6_dp
      becke_control%cavity_shape = radius_single
      becke_control%vector_buffer%store_vectors = .TRUE.
      NULLIFY (becke_control%cavity%pw)
      NULLIFY (becke_control%aij)
      NULLIFY (becke_control%cavity_mat)
      NULLIFY (becke_control%cavity_env)
      NULLIFY (becke_control%cutoffs)
      NULLIFY (becke_control%cutoffs_tmp)
      NULLIFY (becke_control%radii)
      NULLIFY (becke_control%radii_tmp)
   END SUBROUTINE becke_control_create

! **************************************************************************************************
!> \brief release the becke_constraint_type
!> \param becke_control the structure to release
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_release(becke_control)
      TYPE(becke_constraint_type), POINTER               :: becke_control

      CPASSERT(ASSOCIATED(becke_control))
      IF (becke_control%vector_buffer%store_vectors) THEN
         IF (ALLOCATED(becke_control%vector_buffer%distances)) &
            DEALLOCATE (becke_control%vector_buffer%distances)
         IF (ALLOCATED(becke_control%vector_buffer%distance_vecs)) &
            DEALLOCATE (becke_control%vector_buffer%distance_vecs)
         IF (ALLOCATED(becke_control%vector_buffer%position_vecs)) &
            DEALLOCATE (becke_control%vector_buffer%position_vecs)
         IF (ALLOCATED(becke_control%vector_buffer%R12)) &
            DEALLOCATE (becke_control%vector_buffer%R12)
         IF (ALLOCATED(becke_control%vector_buffer%pair_dist_vecs)) &
            DEALLOCATE (becke_control%vector_buffer%pair_dist_vecs)
      END IF
      IF (ASSOCIATED(becke_control%cutoffs)) &
         DEALLOCATE (becke_control%cutoffs)
      IF (ASSOCIATED(becke_control%cutoffs_tmp)) &
         DEALLOCATE (becke_control%cutoffs_tmp)
      IF (ASSOCIATED(becke_control%radii_tmp)) &
         DEALLOCATE (becke_control%radii_tmp)
      IF (ASSOCIATED(becke_control%radii)) &
         DEALLOCATE (becke_control%radii)
      IF (ASSOCIATED(becke_control%aij)) &
         DEALLOCATE (becke_control%aij)
      IF (ASSOCIATED(becke_control%cavity_mat)) &
         DEALLOCATE (becke_control%cavity_mat)
      IF (becke_control%cavity_confine) &
         CALL release_hirshfeld_type(becke_control%cavity_env)
      DEALLOCATE (becke_control)

   END SUBROUTINE becke_control_release

! **************************************************************************************************
!> \brief create the cdft_control_type
!> \param cdft_control the structure to create
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_create(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CPASSERT(.NOT. ASSOCIATED(cdft_control))
      ALLOCATE (cdft_control)
      cdft_control%ref_count = 1
      cdft_control%total_steps = 0
      NULLIFY (cdft_control%strength)
      NULLIFY (cdft_control%target)
      NULLIFY (cdft_control%value)
      NULLIFY (cdft_control%atoms)
      NULLIFY (cdft_control%is_constraint)
      NULLIFY (cdft_control%charges_fragment)
      NULLIFY (cdft_control%fragments)
      NULLIFY (cdft_control%group)
      NULLIFY (cdft_control%charge)
      cdft_control%natoms = 0
      cdft_control%type = outer_scf_none
      cdft_control%need_pot = .TRUE.
      cdft_control%save_pot = .FALSE.
      cdft_control%transfer_pot = .FALSE.
      cdft_control%atomic_charges = .FALSE.
      cdft_control%first_iteration = .TRUE.
      cdft_control%fragment_density = .FALSE.
      cdft_control%fragments_integrated = .FALSE.
      cdft_control%flip_fragment = .FALSE.
      cdft_control%external_control = .FALSE.
      cdft_control%do_et = .FALSE.
      cdft_control%reuse_precond = .FALSE.
      cdft_control%nreused = 0
      cdft_control%precond_freq = 0
      cdft_control%max_reuse = 0
      cdft_control%should_purge = .FALSE.
      cdft_control%purge_history = .FALSE.
      cdft_control%calculate_metric = .FALSE.
      cdft_control%in_memory = .FALSE.
      cdft_control%purge_freq = 0
      cdft_control%nbad_conv = 0
      cdft_control%purge_offset = 0
      cdft_control%istep = 0
      cdft_control%ienergy = 0
      NULLIFY (cdft_control%becke_control)
      CALL becke_control_create(cdft_control%becke_control)
      NULLIFY (cdft_control%hirshfeld_control)
      CALL hirshfeld_control_create(cdft_control%hirshfeld_control)
      NULLIFY (cdft_control%wmat)
      NULLIFY (cdft_control%matrix_s%matrix)
      NULLIFY (cdft_control%mo_coeff)
      NULLIFY (cdft_control%matrix_p)
      ! Outer SCF default settings
      cdft_control%ot_control%have_scf = .FALSE.
      cdft_control%ot_control%max_scf = 0
      cdft_control%ot_control%eps_scf = 0.0_dp
      cdft_control%ot_control%step_size = 0.0_dp
      cdft_control%ot_control%type = -1
      cdft_control%ot_control%optimizer = -1
      cdft_control%ot_control%diis_buffer_length = -1
      NULLIFY (cdft_control%ot_control%cdft_opt_control)
      cdft_control%constraint_control%have_scf = .FALSE.
      cdft_control%constraint_control%max_scf = 0
      cdft_control%constraint_control%eps_scf = 0.0_dp
      cdft_control%constraint_control%step_size = 0.0_dp
      cdft_control%constraint_control%type = -1
      cdft_control%constraint_control%optimizer = -1
      cdft_control%constraint_control%diis_buffer_length = -1
      NULLIFY (cdft_control%constraint_control%cdft_opt_control)
      cdft_control%constraint%iter_count = 0
      NULLIFY (cdft_control%constraint%variables)
      NULLIFY (cdft_control%constraint%gradient)
      NULLIFY (cdft_control%constraint%energy)
      NULLIFY (cdft_control%constraint%count)
      NULLIFY (cdft_control%constraint%inv_jacobian)
      cdft_control%constraint%deallocate_jacobian = .TRUE.
   END SUBROUTINE cdft_control_create

! **************************************************************************************************
!> \brief release the cdft_control_type
!> \param cdft_control the structure to release
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_release(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(cdft_control))
      CPASSERT(cdft_control%ref_count > 0)
      cdft_control%ref_count = cdft_control%ref_count - 1
      IF (cdft_control%ref_count == 0) THEN
         ! Constraint settings
         IF (ASSOCIATED(cdft_control%atoms)) &
            DEALLOCATE (cdft_control%atoms)
         IF (ASSOCIATED(cdft_control%strength)) &
            DEALLOCATE (cdft_control%strength)
         IF (ASSOCIATED(cdft_control%target)) &
            DEALLOCATE (cdft_control%target)
         IF (ASSOCIATED(cdft_control%value)) &
            DEALLOCATE (cdft_control%value)
         IF (ASSOCIATED(cdft_control%charges_fragment)) &
            DEALLOCATE (cdft_control%charges_fragment)
         IF (ASSOCIATED(cdft_control%fragments)) &
            DEALLOCATE (cdft_control%fragments)
         IF (ASSOCIATED(cdft_control%is_constraint)) &
            DEALLOCATE (cdft_control%is_constraint)
         IF (ASSOCIATED(cdft_control%charge)) &
            DEALLOCATE (cdft_control%charge)
         ! Constraint atom groups
         IF (ASSOCIATED(cdft_control%group)) THEN
            DO i = 1, SIZE(cdft_control%group)
               IF (ASSOCIATED(cdft_control%group(i)%atoms)) &
                  DEALLOCATE (cdft_control%group(i)%atoms)
               IF (ASSOCIATED(cdft_control%group(i)%coeff)) &
                  DEALLOCATE (cdft_control%group(i)%coeff)
               IF (ALLOCATED(cdft_control%group(i)%d_sum_const_dR)) &
                  DEALLOCATE (cdft_control%group(i)%d_sum_const_dR)
                  IF (ASSOCIATED(cdft_control%group(i)%gradients)) &
                     DEALLOCATE (cdft_control%group(i)%gradients)
                  IF (ASSOCIATED(cdft_control%group(i)%gradients_x)) &
                     DEALLOCATE (cdft_control%group(i)%gradients_x)
                  IF (ASSOCIATED(cdft_control%group(i)%gradients_y)) &
                     DEALLOCATE (cdft_control%group(i)%gradients_y)
                  IF (ASSOCIATED(cdft_control%group(i)%gradients_z)) &
                     DEALLOCATE (cdft_control%group(i)%gradients_z)
               IF (ASSOCIATED(cdft_control%group(i)%integrated)) &
                  DEALLOCATE (cdft_control%group(i)%integrated)
            END DO
            DEALLOCATE (cdft_control%group)
         END IF
         ! Constraint type specific deallocations
         IF (ASSOCIATED(cdft_control%becke_control)) THEN
            CALL becke_control_release(cdft_control%becke_control)
         END IF
         IF (ASSOCIATED(cdft_control%hirshfeld_control)) THEN
            CALL hirshfeld_control_release(cdft_control%hirshfeld_control)
         END IF
         ! Release OUTER_SCF types
         CALL cdft_opt_type_release(cdft_control%ot_control%cdft_opt_control)
         CALL cdft_opt_type_release(cdft_control%constraint_control%cdft_opt_control)
         IF (ASSOCIATED(cdft_control%constraint%variables)) &
            DEALLOCATE (cdft_control%constraint%variables)
         IF (ASSOCIATED(cdft_control%constraint%count)) &
            DEALLOCATE (cdft_control%constraint%count)
         IF (ASSOCIATED(cdft_control%constraint%gradient)) &
            DEALLOCATE (cdft_control%constraint%gradient)
         IF (ASSOCIATED(cdft_control%constraint%energy)) &
            DEALLOCATE (cdft_control%constraint%energy)
         IF (ASSOCIATED(cdft_control%constraint%inv_jacobian)) &
            DEALLOCATE (cdft_control%constraint%inv_jacobian)
         ! Storage for mixed CDFT calculations
         IF (ALLOCATED(cdft_control%occupations)) THEN
            DO i = 1, SIZE(cdft_control%occupations)
               IF (ASSOCIATED(cdft_control%occupations(i)%array)) &
                  DEALLOCATE (cdft_control%occupations(i)%array)
            END DO
            DEALLOCATE (cdft_control%occupations)
         END IF
         ! Release control
         cdft_control%type = outer_scf_none
         DEALLOCATE (cdft_control)
      END IF
   END SUBROUTINE cdft_control_release

! **************************************************************************************************
!> \brief retain the cdft_control_type
!> \param cdft_control the structure to retain
!> \par History
!>      created 12.2015 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_retain(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CPASSERT(ASSOCIATED(cdft_control))
      cdft_control%ref_count = cdft_control%ref_count + 1
   END SUBROUTINE cdft_control_retain

! **************************************************************************************************
!> \brief create the hirshfeld_constraint_type
!> \param hirshfeld_control the structure to create
!> \par History
!>      09.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE hirshfeld_control_create(hirshfeld_control)
      TYPE(hirshfeld_constraint_type), POINTER           :: hirshfeld_control

      CPASSERT(.NOT. ASSOCIATED(hirshfeld_control))
      ALLOCATE (hirshfeld_control)

      hirshfeld_control%use_bohr = .FALSE.
      hirshfeld_control%print_density = .FALSE.
      hirshfeld_control%use_atomic_cutoff = .TRUE.
      hirshfeld_control%radius = 3.0_dp
      hirshfeld_control%eps_cutoff = 1.0e-12_dp
      hirshfeld_control%atomic_cutoff = 1.0e-12_dp
      hirshfeld_control%shape_function = shape_function_gaussian
      hirshfeld_control%atoms_memory = 80
      hirshfeld_control%gaussian_shape = radius_single
      NULLIFY (hirshfeld_control%hirshfeld_env)
      NULLIFY (hirshfeld_control%radii)

   END SUBROUTINE hirshfeld_control_create

! **************************************************************************************************
!> \brief release the hirshfeld_constraint_type
!> \param hirshfeld_control the structure to release
!> \par History
!>      09.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE hirshfeld_control_release(hirshfeld_control)
      TYPE(hirshfeld_constraint_type), POINTER           :: hirshfeld_control

      CPASSERT(ASSOCIATED(hirshfeld_control))

      IF (ASSOCIATED(hirshfeld_control%radii)) &
         DEALLOCATE (hirshfeld_control%radii)
      CALL release_hirshfeld_type(hirshfeld_control%hirshfeld_env)
      DEALLOCATE (hirshfeld_control)

   END SUBROUTINE hirshfeld_control_release

END MODULE qs_cdft_types
