!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines CDFT control structures
!> \par   History
!>                 separated from cp_control_types [03.2017]
!> \author Nico Holmberg [03.2017]
! **************************************************************************************************
MODULE qs_cdft_types
   USE cp_fm_types,                     ONLY: cp_fm_p_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE hirshfeld_types,                 ONLY: hirshfeld_type,&
                                              release_hirshfeld_type
   USE input_constants,                 ONLY: becke_cutoff_global,&
                                              outer_scf_none,&
                                              radius_single
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE pw_types,                        ONLY: pw_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief some parameters useful for becke_constraints
!> \param fragment_a_fname filename of cube file holding the isolated electron density of fragment a
!> \param fragment_b_fname filename of cube file holding the isolated electron density of fragment b
!> \param ref_count the ref count
!> \param rglobal global cutoff to use for building the constraint
!> \param target target values of the constraints
!> \param strength Lagrangian multipliers of the constraints
!> \param becke_order_p integrated values of the constraints
!> \param natoms number of constraint atoms
!> \param cavity_shape the confinement cavity shape id
!> \param cutoff_type the cutoff type to use for building the constraint
!> \param atoms list of constraint atoms
!> \param cavity the Gaussian confinement cavity: the constraint is nonzero outside this cavity
!> \param group container for atom groups each defining their own constraint
!> \param charge atomic Becke real space potentials needed to calculate atomic Becke charges
!> \param fragments container for isolated fragment densities read from cube files
!> \param need_pot logical which determines if the Becke potential needs to be built
!> \param save_pot logical which determines if the Becke potential should be saved until forces
!>                 have been evaluated
!> \param in_memory logical which determines if the gradients of the Becke potential should be
!>                  computed simultaneously with the potential instead of separately
!> \param adjust logical which determines if the Becke potential is adjusted with atomic radii
!> \param atomic_charges logical which determines if atomic Becke charges should be computed
!> \param cavity_confine logical which determines if cavity confinement is active
!> \param should_skip logical which determines is grid points should be skipped if all constraint
!>                    atoms are found to reside beyond the cutoff distance from it
!> \param print_cavity logical to print the Gaussian confinement cavity
!> \param external_control logical which determines if the constraint has already been built
!>                         in a mixed_env that holds multiple CDFT states
!> \param first_iteration a flag to mark the first iteration for printing of additional data
!> \param fragment_density use isolated fragment densities as a reference for the constraint
!> \param fragments_integrated logical to determine if the fragment densities have been integrated
!> \param use_bohr decides whether to use angstrom or bohr units for the confinement cavity radius
!> \param cutoffs element specific cutoffs
!> \param cutoffs_tmp same as cutoffs but a temporary read during parsing of this type
!> \param charges_fragment the values of the integrated isolated fragment densities
!> \param radii_tmp temporary list of element specific atomic radii used to adjust the Becke cells
!> \param radii permanent copy of radii_tmp
!> \param confine_bounds grid point indices outside which the constraint vanishes along Z-axis
!> \param rcavity an optional global radius parameter used to define the Gaussian confinement cavity
!> \param eps_cavity threshold used screen small values of the Gaussian cavity density
!> \param cavity_env the structure used to build the Gaussian cavity
!> \param aij pairwise parameters used to adjust the Becke cell boundaries built from atomic radii
!> \param cavity_mat a compacted version of cavity
! **************************************************************************************************
   ! To build multiple constraints
   TYPE becke_group_type
      ! Atoms of this constraint group
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      ! Constraint type: charge constraint or magnetization density constraint
      LOGICAL                              :: is_charge_constraint
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
      ! The weight function of this constraint group
      TYPE(pw_p_type)                      :: weight
   END TYPE becke_group_type

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
      CHARACTER(LEN=default_path_length)   :: fragment_a_fname, &
                                              fragment_b_fname
      INTEGER                              :: ref_count, &
                                              natoms, &
                                              cavity_shape, cutoff_type, &
                                              confine_bounds(2)
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      LOGICAL                              :: need_pot, save_pot, in_memory, &
                                              adjust, atomic_charges, cavity_confine, &
                                              should_skip, print_cavity, &
                                              external_control, first_iteration, &
                                              fragment_density, fragments_integrated, &
                                              use_bohr, transfer_pot
      LOGICAL, POINTER, DIMENSION(:)       :: is_constraint
      REAL(KIND=dp)                        :: rglobal, &
                                              rcavity, eps_cavity
      REAL(KIND=dp), DIMENSION(:), POINTER :: TARGET, strength, &
                                              becke_order_p, &
                                              cutoffs, cutoffs_tmp, &
                                              charges_fragment, radii_tmp, &
                                              radii
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                   :: aij
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                :: cavity_mat
      TYPE(becke_group_type), POINTER, &
         DIMENSION(:)                      :: group
      TYPE(becke_vector_buffer)            :: vector_buffer
      TYPE(hirshfeld_type), POINTER        :: cavity_env
      TYPE(pw_p_type)                      :: cavity
      TYPE(pw_p_type), POINTER, &
         DIMENSION(:)                      :: fragments, charge
   END TYPE becke_constraint_type

! **************************************************************************************************
! \brief control parameters for hirshfeld constraint
!> \param hirshfeld_env holds the information needed to build the Gaussian Hirshfeld weight
!> \param natoms the number of constraint atoms
!> \param constraint_type the constraint type to use
!> \param atoms list of constraint atoms
!> \param coeff determines how to sum up the constraint atoms to form the value of the constraint
! **************************************************************************************************
   TYPE hirshfeld_constraint_type
      TYPE(hirshfeld_type), POINTER        :: hirshfeld_env
      INTEGER                              :: natoms, constraint_type
      INTEGER, POINTER, DIMENSION(:)       :: atoms
      REAL(KIND=dp), POINTER, DIMENSION(:) :: coeff
   END TYPE hirshfeld_constraint_type

! **************************************************************************************************
! \brief control parameters for CDFT with OT
!        TODO: make cdft_control parent section for Becke restraint type
!> \param ref_count the ref count
!> \param total_steps counter to keep track of the total number of SCF steps
!> \param type type of CDFT constraint
!> \param precond_freq preconditioner can be used if SCF converged in less than precond_freq steps
!> \param nreused determines how many times the current OT preconditioner has been reused
!> \param max_reuse the same preconditioner can be used a maximum of max_reuse times
!> \param purge_freq determines how large nbad_conv can grow before purging the wfn/constraint history
!> \param nbad_conv a running counter keeping track of the number of CDFT SCF loops when the first
!>                  CDFT SCF iteration required more than 1 outer SCF loop. Reset when convergence is
!>                  smooth
!> \param purge_offset purging is only allowed when more than purge_offset steps have passed since
!>                     last purge
!> \param istep a counter to keep track of how many steps have passed since the last purge
!> \param ienergy a counter tracking the total number of CDFT energy evaluations
!> \param constraints holds information about the CDFT SCF loop
!> \param constraint_control the outer_scf_control_type for the constraint
!> \param ot_control the outer_scf_control_type for OT where data is stashed when outside the OT
!>                   outer loop
!> \param hirshfeld_control control parameters for a Hirshfeld constraint
!> \param strength Lagrangian multipliers of the constraints
!> \param target target values of the constraints
!> \param value integrated values of the constraints
!> \param weight the constraint potential
!> \param need_pot logical which determines if the constraint potential needs to be built
!> \param save_pot logical which determines if the constraint potential should be saved until forces
!>                 have been evaluated
!> \param do_et logical which determines if a ET coupling calculation was requested
!> \param reuse_precond logical which determines if a preconditioner can be reused
!> \param purge_history logical which determines if the wfn/constraint history can be purged
!> \param should_purge logical which determines if purging should take place after this CDFT SCF loop
!> \param calculate_metric logical which determines if the ET coupling reliablity metric is computed
!> \param is_charge_constraint list of logicals which determine if the constraint should be applied
!>                             to the total charge density or the magnetization density
!> \param mo_coeff save the MO coeffs (for do_et)
!> \param wmat matrix representation of the weight function (for do_et)
!> \param matrix_s save the overlap matrix (for do_et)
!> \param matrix_p save the density matrix (for calculate_metric)
! **************************************************************************************************
   ! Copied from qs_scf_types to avoid circular dependency
   TYPE qs_outer_scf_type
      INTEGER                              :: iter_count
      LOGICAL                              :: deallocate_jacobian
      REAL(KIND=dp), DIMENSION(:), POINTER :: energy
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: variables
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: gradient
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                           :: inv_jacobian
      INTEGER, DIMENSION(:), POINTER       :: count
   END TYPE qs_outer_scf_type

   ! Copied from scf_control_types to avoid circular dependency
   TYPE outer_scf_control_type
      LOGICAL                              :: have_scf
      LOGICAL                              :: build_jacobian
      LOGICAL                              :: broyden_update
      INTEGER                              :: max_scf, &
                                              ijacobian
      REAL(KIND=dp)                        :: eps_scf, step_size, &
                                              jacobian_step, &
                                              newton_step
      INTEGER                              :: TYPE
      INTEGER                              :: optimizer
      INTEGER                              :: diis_buffer_length
      INTEGER                              :: extrapolation_order
      INTEGER                              :: bisect_trust_count
      INTEGER                              :: jacobian_type
      INTEGER                              :: broyden_type
      INTEGER                              :: jacobian_freq(2)
   END TYPE outer_scf_control_type

   TYPE cdft_control_type
      INTEGER                              :: ref_count, total_steps, TYPE, &
                                              precond_freq, nreused, max_reuse, &
                                              purge_freq, nbad_conv, purge_offset, &
                                              istep, ienergy
      LOGICAL                              :: need_pot, save_pot, do_et, &
                                              reuse_precond, purge_history, &
                                              should_purge, calculate_metric
      LOGICAL, DIMENSION(:), ALLOCATABLE   :: is_charge_constraint
      REAL(KIND=dp), DIMENSION(:), POINTER :: strength, TARGET, value
      TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                           :: mo_coeff
      TYPE(dbcsr_p_type)                   :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER                           :: wmat, matrix_p
      TYPE(hirshfeld_constraint_type), &
         POINTER                           :: hirshfeld_control
      TYPE(outer_scf_control_type)         :: constraint_control, ot_control
      TYPE(pw_p_type)                      :: weight
      TYPE(qs_outer_scf_type)              :: constraint
   END TYPE cdft_control_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_types'

   ! Public data types

   PUBLIC :: becke_constraint_type, &
             cdft_control_type, &
             becke_group_type

   ! Public subroutines

   PUBLIC :: becke_control_create, &
             becke_control_release, &
             cdft_control_create, &
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

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(becke_control))
      ALLOCATE (becke_control)

      becke_control%ref_count = 1
      becke_control%natoms = 0
      becke_control%need_pot = .TRUE.
      becke_control%save_pot = .FALSE.
      becke_control%transfer_pot = .FALSE.
      becke_control%adjust = .FALSE.
      becke_control%atomic_charges = .FALSE.
      becke_control%cutoff_type = becke_cutoff_global
      becke_control%cavity_confine = .FALSE.
      becke_control%should_skip = .FALSE.
      becke_control%print_cavity = .FALSE.
      becke_control%in_memory = .FALSE.
      becke_control%first_iteration = .TRUE.
      becke_control%fragment_density = .FALSE.
      becke_control%fragments_integrated = .FALSE.
      becke_control%use_bohr = .FALSE.
      becke_control%confine_bounds = 0
      becke_control%rcavity = 3.0_dp
      becke_control%rglobal = 6.0_dp
      becke_control%eps_cavity = 1.0e-5_dp
      becke_control%cavity_shape = radius_single
      becke_control%external_control = .FALSE.
      becke_control%vector_buffer%store_vectors = .TRUE.
      NULLIFY (becke_control%cavity%pw)
      NULLIFY (becke_control%strength)
      NULLIFY (becke_control%target)
      NULLIFY (becke_control%becke_order_p)
      NULLIFY (becke_control%atoms)
      NULLIFY (becke_control%group)
      NULLIFY (becke_control%charge)
      NULLIFY (becke_control%aij)
      NULLIFY (becke_control%is_constraint)
      NULLIFY (becke_control%cavity_mat)
      NULLIFY (becke_control%cavity_env)
      NULLIFY (becke_control%cutoffs)
      NULLIFY (becke_control%cutoffs_tmp)
      NULLIFY (becke_control%charges_fragment)
      NULLIFY (becke_control%fragments)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(becke_control))
      CPASSERT(becke_control%ref_count > 0)
      becke_control%ref_count = becke_control%ref_count-1
      IF (becke_control%ref_count == 0) THEN
         IF (ASSOCIATED(becke_control%atoms)) &
            DEALLOCATE (becke_control%atoms)
         IF (ASSOCIATED(becke_control%group)) THEN
            DO i = 1, SIZE(becke_control%group)
               IF (ASSOCIATED(becke_control%group(i)%atoms)) &
                  DEALLOCATE (becke_control%group(i)%atoms)
               IF (ASSOCIATED(becke_control%group(i)%coeff)) &
                  DEALLOCATE (becke_control%group(i)%coeff)
               IF (ALLOCATED(becke_control%group(i)%d_sum_const_dR)) &
                  DEALLOCATE (becke_control%group(i)%d_sum_const_dR)
               IF (ASSOCIATED(becke_control%group(i)%gradients)) &
                  DEALLOCATE (becke_control%group(i)%gradients)
               IF (ASSOCIATED(becke_control%group(i)%integrated)) &
                  DEALLOCATE (becke_control%group(i)%integrated)
            END DO
            DEALLOCATE (becke_control%group)
         END IF
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
         IF (ASSOCIATED(becke_control%charges_fragment)) &
            DEALLOCATE (becke_control%charges_fragment)
         IF (ASSOCIATED(becke_control%fragments)) &
            DEALLOCATE (becke_control%fragments)
         IF (ASSOCIATED(becke_control%charge)) &
            DEALLOCATE (becke_control%charge)
         IF (ASSOCIATED(becke_control%aij)) &
            DEALLOCATE (becke_control%aij)
         IF (ASSOCIATED(becke_control%is_constraint)) &
            DEALLOCATE (becke_control%is_constraint)
         IF (ASSOCIATED(becke_control%cavity_mat)) &
            DEALLOCATE (becke_control%cavity_mat)
         IF (ASSOCIATED(becke_control%target)) &
            DEALLOCATE (becke_control%target)
         IF (ASSOCIATED(becke_control%strength)) &
            DEALLOCATE (becke_control%strength)
         IF (ASSOCIATED(becke_control%becke_order_p)) &
            DEALLOCATE (becke_control%becke_order_p)
         IF (becke_control%cavity_confine) &
            CALL release_hirshfeld_type(becke_control%cavity_env)
         becke_control%ref_count = 0
         becke_control%natoms = 0
         DEALLOCATE (becke_control)
      ENDIF
   END SUBROUTINE becke_control_release

! **************************************************************************************************
!> \brief retain the becke_constraint_type
!> \param becke_control the structure to retain
!> \par History
!>      02.2007 created [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE becke_control_retain(becke_control)
      TYPE(becke_constraint_type), POINTER               :: becke_control

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(becke_control))
      becke_control%ref_count = becke_control%ref_count+1
   END SUBROUTINE becke_control_retain

! **************************************************************************************************
!> \brief create the cdft_control_type
!> \param cdft_control the structure to create
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE cdft_control_create(cdft_control)
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_create', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(cdft_control))
      ALLOCATE (cdft_control)
      cdft_control%ref_count = 1
      cdft_control%total_steps = 0
      NULLIFY (cdft_control%strength)
      NULLIFY (cdft_control%target)
      NULLIFY (cdft_control%value)
      cdft_control%type = outer_scf_none
      cdft_control%need_pot = .TRUE.
      cdft_control%save_pot = .FALSE.
      cdft_control%do_et = .FALSE.
      cdft_control%reuse_precond = .FALSE.
      cdft_control%nreused = 0
      cdft_control%precond_freq = 0
      cdft_control%max_reuse = 0
      cdft_control%should_purge = .FALSE.
      cdft_control%purge_history = .FALSE.
      cdft_control%calculate_metric = .FALSE.
      cdft_control%purge_freq = 0
      cdft_control%nbad_conv = 0
      cdft_control%purge_offset = 0
      cdft_control%istep = 0
      cdft_control%ienergy = 0
      NULLIFY (cdft_control%hirshfeld_control)
      NULLIFY (cdft_control%weight%pw)
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
      cdft_control%ot_control%jacobian_type = -1
      cdft_control%ot_control%jacobian_step = 0.0_dp
      cdft_control%ot_control%newton_step = 1.0_dp
      cdft_control%constraint_control%have_scf = .FALSE.
      cdft_control%constraint_control%max_scf = 0
      cdft_control%constraint_control%eps_scf = 0.0_dp
      cdft_control%constraint_control%step_size = 0.0_dp
      cdft_control%constraint_control%type = -1
      cdft_control%constraint_control%optimizer = -1
      cdft_control%constraint_control%diis_buffer_length = -1
      cdft_control%constraint_control%jacobian_type = -1
      cdft_control%constraint_control%broyden_type = -1
      cdft_control%constraint_control%jacobian_freq = 1
      cdft_control%constraint_control%jacobian_step = 0.0_dp
      cdft_control%constraint_control%newton_step = 1.0_dp
      cdft_control%constraint_control%ijacobian = 0
      cdft_control%constraint_control%build_jacobian = .FALSE.
      cdft_control%constraint_control%broyden_update = .FALSE.
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

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_release', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(cdft_control))
      CPASSERT(cdft_control%ref_count > 0)
      cdft_control%ref_count = cdft_control%ref_count-1
      IF (cdft_control%ref_count == 0) THEN
         IF (ASSOCIATED(cdft_control%strength)) &
            DEALLOCATE (cdft_control%strength)
         IF (ASSOCIATED(cdft_control%target)) &
            DEALLOCATE (cdft_control%target)
         IF (ASSOCIATED(cdft_control%value)) &
            DEALLOCATE (cdft_control%value)
         IF (ALLOCATED(cdft_control%is_charge_constraint)) &
            DEALLOCATE (cdft_control%is_charge_constraint)
         cdft_control%type = outer_scf_none
         NULLIFY (cdft_control%weight%pw)
         IF (ASSOCIATED(cdft_control%hirshfeld_control)) THEN
            CALL release_hirshfeld_type(cdft_control%hirshfeld_control%hirshfeld_env)
            cdft_control%hirshfeld_control%natoms = 0
            IF (ASSOCIATED(cdft_control%hirshfeld_control%atoms)) &
               DEALLOCATE (cdft_control%hirshfeld_control%atoms)
            IF (ASSOCIATED(cdft_control%hirshfeld_control%coeff)) &
               DEALLOCATE (cdft_control%hirshfeld_control%coeff)
            DEALLOCATE (cdft_control%hirshfeld_control)
         END IF
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

      CHARACTER(len=*), PARAMETER :: routineN = 'cdft_control_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(cdft_control))
      cdft_control%ref_count = cdft_control%ref_count+1
   END SUBROUTINE cdft_control_retain

END MODULE qs_cdft_types
