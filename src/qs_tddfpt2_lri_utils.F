!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_lri_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE auto_basis,                      ONLY: create_lri_aux_basis_set
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_output,                ONLY: print_basis_set_file
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set,&
                                              sort_gto_basis_set
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lri_environment_init,            ONLY: lri_env_basis,&
                                              lri_env_init
   USE lri_environment_methods,         ONLY: build_lri_matrices
   USE lri_environment_types,           ONLY: lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup,&
                                              write_neighbor_lists
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_lri_utils'

   PUBLIC:: tddfpt2_lri_init, tddfpt2_lri_Amat

CONTAINS

! **************************************************************************************************
!> \brief Initialize LRI environment, basis, neighborlists and matrices
!> \param qs_env ...
!> \param kernel_env ...
!> \param lri_section ...
!> \param tddfpt_print_section ...
! **************************************************************************************************
   SUBROUTINE tddfpt2_lri_init(qs_env, kernel_env, lri_section, tddfpt_print_section)

      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: lri_section, tddfpt_print_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tddfpt2_lri_init'

      INTEGER                                            :: basis_sort, handle, ikind, lmax_sphere, &
                                                            maxl, maxlgto, maxlgto_lri, nkind
      LOGICAL                                            :: explicit, mic, molecule_only, &
                                                            redefine_interaction_radii
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set, p_lri_aux_basis
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: neighbor_list_section, print_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt2_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt2_control => dft_control%tddfpt2_control

      NULLIFY (kernel_env%flex_kernel%lri_env)
! initialize lri_env
      CALL lri_env_init(kernel_env%flex_kernel%lri_env, lri_section)
      NULLIFY (lri_env)
      lri_env => kernel_env%flex_kernel%lri_env
      redefine_interaction_radii = .FALSE.

      !     exact_1c_terms not implemented
      IF (lri_env%exact_1c_terms) THEN
         CPABORT("TDDFT with LRI and exact one-center terms not implemented")
      END IF

!       ! check if LRI_AUX basis is available
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      DO ikind = 1, nkind
         NULLIFY (p_lri_aux_basis)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, basis_set=p_lri_aux_basis, basis_type="P_LRI_AUX")
         IF (.NOT. (ASSOCIATED(p_lri_aux_basis))) THEN
!            ! Generate a default basis
            redefine_interaction_radii = .TRUE.
            CALL cp_warn(__LOCATION__, "Automatic Generation of P_LRI_AUX basis. "// &
                         "This is experimental code.")
            CALL create_lri_aux_basis_set(p_lri_aux_basis, qs_kind, &
                                          dft_control%auto_basis_p_lri_aux, lri_env%exact_1c_terms)
            CALL add_basis_set_to_container(qs_kind%basis_sets, p_lri_aux_basis, "P_LRI_AUX")
         END IF
      END DO !nkind
      ! needs to be done here if p_lri_aux_basis is not specified explicitly
      IF (redefine_interaction_radii) THEN
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=p_lri_aux_basis, basis_type="P_LRI_AUX")
            IF (ASSOCIATED(p_lri_aux_basis)) THEN
               CALL init_orb_basis_set(p_lri_aux_basis) ! standardly done in init_qs_kind
               basis_sort = 0
               CALL sort_gto_basis_set(p_lri_aux_basis, basis_sort)
               CALL init_interaction_radii_orb_basis(p_lri_aux_basis, dft_control%qs_control%eps_pgf_orb)
            END IF
         END DO
      END IF

      print_section => section_vals_get_subs_vals(tddfpt_print_section, "BASIS_SET_FILE")
      CALL section_vals_get(print_section, explicit=explicit)
      IF (explicit) THEN
         CALL print_basis_set_file(qs_env, print_section)
      END IF

      !set maxl as in qs_environment for gs lri
      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto)
      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto_lri, basis_type="P_LRI_AUX")
      !take maxlgto from lri basis if larger (usually)
      maxlgto = MAX(maxlgto, maxlgto_lri)
      lmax_sphere = 2*maxlgto
      maxl = MAX(2*maxlgto, lmax_sphere) + 1

      CALL init_orbital_pointers(maxl)
      CALL init_spherical_harmonics(maxl, 0)

      ! initialize LRI basis sets
      CALL lri_env_basis("P_LRI", qs_env, lri_env, qs_kind_set)

!        check for debugging that automatically generated basis is available
!         DO ikind= 1, nkind
!            qs_kind => qs_kind_set(ikind)
!            CALL get_qs_kind(qs_kind, basis_set=p_lri_aux_basis, basis_type="P_LRI_AUX")
!           CALL get_gto_basis_set(gto_basis_set=p_lri_aux_basis, set_radius=set_radius,&
!                                  pgf_radius=pgf_radius)
!         END DO !nkind

      ! set up LRI neighbor list soo_list => same as in qs_neighbor_lists for ground-state LRI
      NULLIFY (cell, para_env, particle_set)
      CALL get_qs_env(qs_env, para_env=para_env, cell=cell, particle_set=particle_set)

      NULLIFY (distribution_1d, distribution_2d, atomic_kind_set, molecule_set)
      CALL get_qs_env(qs_env, local_particles=distribution_1d, distribution_2d=distribution_2d, &
                      atomic_kind_set=atomic_kind_set, molecule_set=molecule_set)

      ALLOCATE (atom2d(nkind))
      molecule_only = .FALSE. !this still needs to be checked
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, molecule_only, particle_set=particle_set)

      ALLOCATE (orb_present(nkind), orb_radius(nkind), pair_radius(nkind, nkind))
      orb_radius(:) = 0.0_dp
      pair_radius(:, :) = 0.0_dp
      orb_present(:) = .FALSE.
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
         END IF
      END DO ! ikind

      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

      NULLIFY (soo_list)
      soo_list => lri_env%soo_list
      mic = .TRUE. !enforcing minimum image convention
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      CALL build_neighbor_lists(soo_list, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, molecular=molecule_only, subcells=subcells, nlname="soo_list")
!
     !!! make this a TDDFT neighborlist
      neighbor_list_section => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%NEIGHBOR_LISTS")
      CALL write_neighbor_lists(soo_list, particle_set, cell, para_env, neighbor_list_section, &
                                "/SOO_LIST", "soo_list", "ORBITAL ORBITAL (RI)")
      lri_env%soo_list => soo_list

      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (orb_present, orb_radius, pair_radius)

! calculate LRI integrals
      lri_env%ppl_ri = .FALSE. ! make sure that option is not available for ES
      CALL build_lri_matrices(lri_env, qs_env)
      kernel_env%flex_kernel%lri_env => lri_env

!   CALL get_condition_number_of_overlap(lri_env)

      CALL timestop(handle)

   END SUBROUTINE tddfpt2_lri_init
! **************************************************************************************************
!> \brief Calculate contribution to response vector for LRI
!> \param qs_env ...
!> \param sub_env ...
!> \param lri_env ...
!> \param lri_v_int ...
!> \param A_ia_munu_sub ...
! **************************************************************************************************
   SUBROUTINE tddfpt2_lri_Amat(qs_env, sub_env, lri_env, lri_v_int, A_ia_munu_sub)
!
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env
      TYPE(lri_environment_type), INTENT(IN), POINTER    :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: lri_v_int
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: A_ia_munu_sub

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tddfpt2_lri_Amat'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dummymat, matrix_s

!

      CALL timeset(routineN, handle)
!
      NULLIFY (atomic_kind_set)
      NULLIFY (matrix_s)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, matrix_s=matrix_s)
!
      nspins = SIZE(A_ia_munu_sub)
      DO ispin = 1, nspins!
         !no kpoints for excited states => using dummy matrix with no cell index
         NULLIFY (dummymat)
         CALL dbcsr_allocate_matrix_set(dummymat, 1)
         CALL tddfpt_dbcsr_create_by_dist(dummymat(1)%matrix, template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_orb)

         CALL dbcsr_copy(A_ia_munu_sub(ispin)%matrix, dummymat(1)%matrix)

         CALL calculate_lri_ks_matrix(lri_env, lri_v_int, dummymat, atomic_kind_set)

         CALL dbcsr_copy(A_ia_munu_sub(ispin)%matrix, dummymat(1)%matrix)

         CALL dbcsr_deallocate_matrix(dummymat(1)%matrix)
         DEALLOCATE (dummymat)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt2_lri_Amat
! **************************************************************************************************
END MODULE qs_tddfpt2_lri_utils
