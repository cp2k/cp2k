!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE optimize_basis_types

   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE powell,                          ONLY: opt_state_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis_types'

   PUBLIC :: basis_optimization_type, subset_type, flex_basis_type, &
             derived_basis_info, deallocate_basis_optimization_type

   ! constraint information for a single constraing. boundary is translateed into a fermi function
   ! like setting as for variational limited case
   TYPE exp_constraint_type
      INTEGER                                            :: const_type = -1
      REAL(KIND=dp)                                      :: llim = -1.0_dp, ulim = -1.0_dp
      REAL(KIND=dp)                                      :: init = -1.0_dp, var_fac = -1.0_dp
   END TYPE

   ! Subset of a basis+ additional information on what to optimize.
   ! *_x_ind maps to the index in the optimization vector
   !  opt_* logical whether quantity ahould be optimized
   !  *_const information for exponents used to constrain them
   TYPE subset_type
      INTEGER                                            :: lmin = -1, lmax = -1, nexp = -1
      INTEGER                                            :: n = -1, ncon_tot = -1, nl = -1
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: l
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: coeff
      LOGICAL, DIMENSION(:, :), ALLOCATABLE              :: opt_coeff
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: coeff_x_ind
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: exps
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: opt_exps
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: exp_x_ind
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: exp_has_const
      TYPE(exp_constraint_type), DIMENSION(:), &
         ALLOCATABLE                                      :: exp_const
   END TYPE

   ! Top level information for basis sets+ vector subset with the real information
   TYPE flex_basis_type
      CHARACTER(LEN=default_string_length)               :: basis_name = ""
      INTEGER                                            :: nopt = -1
      INTEGER                                            :: nsets = -1
      TYPE(subset_type), DIMENSION(:), ALLOCATABLE       :: subset
   END TYPE

   ! information for optimization: whether coeff has to be optimized or not
   TYPE use_contr_type
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: in_use
   END TYPE

   ! information about how to generate the derived basis sets
   TYPE derived_basis_info
      CHARACTER(LEN=default_string_length)               :: basis_name = ""
      INTEGER                                            :: reference_set = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: remove_contr
      INTEGER                                            :: nsets = -1, ncontr = -1
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: remove_set
      LOGICAL, DIMENSION(:), ALLOCATABLE                 :: in_use_set
      TYPE(use_contr_type), DIMENSION(:), ALLOCATABLE    :: use_contr
   END TYPE

   ! some usual stuff for basis information and an info type containing the
   ! the translated input on how to genrate the derived basis sets
   ! a flexible basis type for every derived basis
   ! ATTENTION: both vectors go from 0:nbasis_deriv. entry 0 is the one specified
   !            in the template basis file
   TYPE kind_basis_type
      CHARACTER(LEN=default_string_length)               :: basis_name = ""
      CHARACTER(LEN=default_string_length)               :: element = ""
      INTEGER                                            :: nbasis_deriv = -1
      TYPE(derived_basis_info), DIMENSION(:), &
         ALLOCATABLE                                      :: deriv_info
      TYPE(flex_basis_type), DIMENSION(:), ALLOCATABLE   :: flex_basis
   END TYPE

   ! vector of length nparallel_groups containing the id's of the calculations in the group
   TYPE comp_group_type
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: member_list
   END TYPE

! **************************************************************************************************
!> \brief type containing all information needed for basis matching
!> \author Florian Schiffmann
! **************************************************************************************************
   TYPE basis_optimization_type
      TYPE(comp_group_type), DIMENSION(:), ALLOCATABLE  :: comp_group
      INTEGER :: ntraining_sets = -1
      INTEGER :: ncombinations = -1
      LOGICAL :: use_condition_number = .FALSE.
      INTEGER, DIMENSION(:), POINTER   :: group_partition => NULL()
      INTEGER :: n_groups_created = -1
      INTEGER, DIMENSION(:), ALLOCATABLE :: sub_sources
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: combination
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: fval_weight
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: condition_weight
      INTEGER :: nkind = -1
      INTEGER :: write_frequency = -1
      INTEGER :: nbasis_deriv_types = -1
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: x_opt
      TYPE(opt_state_type)  :: powell_param = opt_state_type()
      CHARACTER(LEN=default_path_length), DIMENSION(:), ALLOCATABLE :: training_input
      CHARACTER(LEN=default_path_length), DIMENSION(:), ALLOCATABLE :: training_dir
      CHARACTER(LEN=default_path_length) :: work_basis_file = ""
      CHARACTER(LEN=default_path_length) :: output_basis_file = ""
      CHARACTER(LEN=default_path_length) :: template_basis_file = ""
      TYPE(kind_basis_type), DIMENSION(:), ALLOCATABLE :: kind_basis
      INTEGER :: opt_id = -1
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Deallocate everything which was allocated before.
!>        Note not all arrays are used depending on the type of basis
!>        i.e derived or reference basis set
!> \param opt_bas ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE deallocate_basis_optimization_type(opt_bas)
      TYPE(basis_optimization_type)                      :: opt_bas

      INTEGER                                            :: igroup, ikind

      IF (ASSOCIATED(opt_bas%group_partition)) DEALLOCATE (opt_bas%group_partition)
      IF (ALLOCATED(opt_bas%sub_sources)) DEALLOCATE (opt_bas%sub_sources)
      IF (ALLOCATED(opt_bas%combination)) DEALLOCATE (opt_bas%combination)
      IF (ALLOCATED(opt_bas%x_opt)) DEALLOCATE (opt_bas%x_opt)
      IF (ALLOCATED(opt_bas%training_input)) DEALLOCATE (opt_bas%training_input)
      IF (ALLOCATED(opt_bas%training_dir)) DEALLOCATE (opt_bas%training_dir)
      IF (ALLOCATED(opt_bas%fval_weight)) DEALLOCATE (opt_bas%fval_weight)
      IF (ALLOCATED(opt_bas%condition_weight)) DEALLOCATE (opt_bas%condition_weight)

      IF (ALLOCATED(opt_bas%comp_group)) THEN
         DO igroup = 1, SIZE(opt_bas%comp_group)
            IF (ALLOCATED(opt_bas%comp_group(igroup)%member_list)) DEALLOCATE (opt_bas%comp_group(igroup)%member_list)
         END DO
         DEALLOCATE (opt_bas%comp_group)
      END IF

      IF (ALLOCATED(opt_bas%kind_basis)) THEN
         DO ikind = 1, SIZE(opt_bas%kind_basis)
            CALL deallocate_kind_basis(opt_bas%kind_basis(ikind))
         END DO
         DEALLOCATE (opt_bas%kind_basis)
      END IF

   END SUBROUTINE deallocate_basis_optimization_type

! **************************************************************************************************
!> \brief Some more deallocation of the subtypes of optimize_absis type
!> \param kind ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE deallocate_kind_basis(kind)
      TYPE(kind_basis_type)                              :: kind

      INTEGER                                            :: ibasis, icont, iinfo, iset

      IF (ALLOCATED(kind%deriv_info)) THEN
         DO iinfo = 0, SIZE(kind%deriv_info) - 1
            IF (ALLOCATED(kind%deriv_info(iinfo)%remove_contr)) DEALLOCATE (kind%deriv_info(iinfo)%remove_contr)
            IF (ALLOCATED(kind%deriv_info(iinfo)%remove_set)) DEALLOCATE (kind%deriv_info(iinfo)%remove_set)
            IF (ALLOCATED(kind%deriv_info(iinfo)%in_use_set)) DEALLOCATE (kind%deriv_info(iinfo)%in_use_set)
            IF (ALLOCATED(kind%deriv_info(iinfo)%use_contr)) THEN
               DO icont = 1, SIZE(kind%deriv_info(iinfo)%use_contr)
                  IF (ALLOCATED(kind%deriv_info(iinfo)%use_contr(icont)%in_use)) &
                     DEALLOCATE (kind%deriv_info(iinfo)%use_contr(icont)%in_use)
               END DO
               DEALLOCATE (kind%deriv_info(iinfo)%use_contr)
            END IF
         END DO
         DEALLOCATE (kind%deriv_info)
      END IF

      IF (ALLOCATED(kind%flex_basis)) THEN
         DO ibasis = 0, SIZE(kind%flex_basis) - 1
            IF (ALLOCATED(kind%flex_basis(ibasis)%subset)) THEN
               DO iset = 1, SIZE(kind%flex_basis(ibasis)%subset)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%l)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%l)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%coeff)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%coeff)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%opt_coeff)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%opt_coeff)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%coeff_x_ind)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%coeff_x_ind)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%exps)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%exps)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%opt_exps)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%opt_exps)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%exp_x_ind)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%exp_x_ind)
                  IF (ALLOCATED(kind%flex_basis(ibasis)%subset(iset)%exp_const)) &
                     DEALLOCATE (kind%flex_basis(ibasis)%subset(iset)%exp_const)
               END DO
               DEALLOCATE (kind%flex_basis(ibasis)%subset)
            END IF
         END DO
         DEALLOCATE (kind%flex_basis)
      END IF

   END SUBROUTINE deallocate_kind_basis

END MODULE optimize_basis_types
