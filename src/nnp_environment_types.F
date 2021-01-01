!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Data types for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
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
                                              cp_subsys_retain,&
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
   PUBLIC :: nnp_acsf_rad_type
   PUBLIC :: nnp_acsf_ang_type

   ! Public subroutines ***
   PUBLIC :: nnp_env_retain, &
             nnp_env_release, &
             nnp_env_set, &
             nnp_env_get, &
             nnp_env_create

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
      TYPE(nnp_acsf_rad_type), DIMENSION(:), POINTER      :: rad ! DIM(n_ele)
      TYPE(nnp_acsf_ang_type), DIMENSION(:), POINTER      :: ang ! DIM(n_ele)
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_rad ! # radial symfnct for this element
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_ang ! # angular symfnct for this element
      INTEGER                                             :: n_ele ! # elements
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele ! elements(n_ele)
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele ! elements(n_ele)
      LOGICAL                                             :: scale_acsf
      LOGICAL                                             :: center_acsf
      LOGICAL                                             :: normnodes
      INTEGER                                             :: n_radgrp
      INTEGER                                             :: n_anggrp
      INTEGER                                             :: cut_type ! cutofftype
      REAL(KIND=dp)                                       :: eshortmin
      REAL(KIND=dp)                                       :: eshortmax
      REAL(KIND=dp)                                       :: scmax !scale
      REAL(KIND=dp)                                       :: scmin !scale
      REAL(KIND=dp)                                       :: max_cut !largest cutoff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: atom_energies !DIM(n_ele)
      TYPE(nnp_arc_type), POINTER, DIMENSION(:)           :: arc ! DIM(n_ele)
      INTEGER                                             :: n_committee
      INTEGER                                             :: n_hlayer
      INTEGER                                             :: n_layer
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_hnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: actfnct
      INTEGER                                             :: expol ! extrapolation coutner
      LOGICAL                                             :: output_expol ! output extrapolation
      ! structures for calculation
      INTEGER                                             :: num_atoms
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: atomic_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)            :: committee_energy
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: ele_ind, nuc_atoms, sort, sort_inv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)         :: coord
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: myforce
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)      :: committee_forces
      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                        :: atoms
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: nnp_forces
      REAL(KIND=dp)                                       :: nnp_potential_energy
      TYPE(cp_subsys_type), POINTER                       :: subsys
      TYPE(section_vals_type), POINTER                    :: nnp_input, &
                                                             force_env_input
      TYPE(cell_type), POINTER                            :: cell_ref
      INTEGER                                             :: ref_count
      LOGICAL                                             :: use_ref_cell
      ! bias
      LOGICAL                                             :: bias
      LOGICAL                                             :: bias_align
      REAL(KIND=dp)                                       :: bias_energy
      REAL(KIND=dp)                                       :: bias_kb
      REAL(KIND=dp)                                       :: bias_sigma0
      REAL(KIND=dp)                                       :: bias_sigma
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: bias_forces
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: bias_e_avrg
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
      INTEGER                                             :: n_symf
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: symf
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: ele_ind
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE         :: ele
      REAL(KIND=dp)                                       :: cutoff
   END TYPE

! **************************************************************************************************
!> \brief Set of radial symmetry function type
!> \param y       - acsf value                                     - DIM(n_rad)
!> \param funccut - distance cutoff                           bohr - DIM(n_rad)
!> \param eta     - eta parameter of radial sym fncts      bohr^-2 - DIM(n_rad)
!> \param rs      - r shift parameter of radial sym fncts     bohr - DIM(n_rad)
!> \param loc_min - minimum of the sym fnct                          DIM(n_rad,n_com)
!> \param loc_max - maximum of the sym fnct                          DIM(n_rad,n_com)
!> \param loc_av  - average of the sym fnct                          DIM(n_rad,n_com)
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
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele
      INTEGER                                             :: n_symfgrp
      TYPE(nnp_symfgrp_type), DIMENSION(:), ALLOCATABLE   :: symfgrp
   END TYPE

! **************************************************************************************************
!> \brief Set of angular symmetry function type
!> \param y         - acsf value                                  - DIM(n_ang)
!> \param funccut   - distance cutoff                        bohr - DIM(n_ang)
!> \param eta       - eta  param. of angular sym fncts    bohr^-2 - DIM(n_ang)
!> \param zeta      - zeta param. of angular sym fncts              DIM(n_ang)
!> \param lam       - lambda  param. of angular sym fncts           DIM(n_ang)
!> \param loc_min   - minimum of the sym fnct                       DIM(n_ang,n_com)
!> \param loc_max   - maximum of the sym fnct                       DIM(n_ang,n_com)
!> \param loc_av    - average of the sym fnct                       DIM(n_ang,n_com)
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
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele1
      CHARACTER(len=2), ALLOCATABLE, DIMENSION(:)         :: ele2
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele1
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nuc_ele2
      INTEGER                                             :: n_symfgrp
      TYPE(nnp_symfgrp_type), DIMENSION(:), ALLOCATABLE   :: symfgrp
   END TYPE

! **************************************************************************************************
!> \brief Contains neighbors list of an atom
!> \param dist - distance vectors + norm                      DIM(4,nat)
!> \param n    - number of neighbors
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_neighbor_type
      INTEGER, DIMENSION(3)                               :: pbc_copies
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_rad
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_ang1
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: n_ang2
      INTEGER, DIMENSION(:, :), ALLOCATABLE               :: ind_rad
      INTEGER, DIMENSION(:, :), ALLOCATABLE               :: ind_ang1
      INTEGER, DIMENSION(:, :), ALLOCATABLE               :: ind_ang2
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dist_rad
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dist_ang1
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dist_ang2
   END TYPE

! **************************************************************************************************
!> \brief Data type for artificial neural networks
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
   TYPE nnp_arc_type
      TYPE(nnp_arc_layer_type), POINTER, DIMENSION(:)     :: layer ! DIM(n_layer)
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: n_nodes
   END TYPE

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
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Create data structure that holds all the information for neural
!>        network potentials
!> \param nnp_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_env_create(nnp_env)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env

      NULLIFY (nnp_env)
      ALLOCATE (nnp_env)
      NULLIFY (nnp_env%rad)
      NULLIFY (nnp_env%ang)
      NULLIFY (nnp_env%arc)
      NULLIFY (nnp_env%subsys)
      NULLIFY (nnp_env%nnp_input)
      NULLIFY (nnp_env%force_env_input)
      NULLIFY (nnp_env%cell_ref)

      nnp_env%ref_count = 1
      nnp_env%use_ref_cell = .FALSE.

   END SUBROUTINE nnp_env_create

! **************************************************************************************************
!> \brief Retains a nnp environment (see doc/ReferenceCounting.html)
!> \param nnp_env The nnp environment to retain
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_env_retain(nnp_env)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env

      CPASSERT(ASSOCIATED(nnp_env))
      CPASSERT(nnp_env%ref_count > 0)
      nnp_env%ref_count = nnp_env%ref_count + 1
   END SUBROUTINE nnp_env_retain

! **************************************************************************************************
!> \brief Release data structure that holds all the information for neural
!>        network potentials
!> \param nnp_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_env_release(nnp_env)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env

      INTEGER                                            :: i, j

      IF (ASSOCIATED(nnp_env)) THEN
         CPASSERT(nnp_env%ref_count > 0)
         nnp_env%ref_count = nnp_env%ref_count - 1
         IF (nnp_env%ref_count < 1) THEN
            IF (ASSOCIATED(nnp_env%rad)) THEN
               DO i = 1, nnp_env%n_ele
                  DO j = 1, nnp_env%rad(i)%n_symfgrp
                     DEALLOCATE (nnp_env%rad(i)%symfgrp(j)%symf, &
                                 nnp_env%rad(i)%symfgrp(j)%ele, &
                                 nnp_env%rad(i)%symfgrp(j)%ele_ind)
                  END DO
                  DEALLOCATE (nnp_env%rad(i)%y, &
                              nnp_env%rad(i)%funccut, &
                              nnp_env%rad(i)%eta, &
                              nnp_env%rad(i)%rs, &
                              nnp_env%rad(i)%loc_min, &
                              nnp_env%rad(i)%loc_max, &
                              nnp_env%rad(i)%loc_av, &
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
            IF (ALLOCATED(nnp_env%atoms)) DEALLOCATE (nnp_env%atoms)
            IF (ALLOCATED(nnp_env%nnp_forces)) DEALLOCATE (nnp_env%nnp_forces)

            IF (ASSOCIATED(nnp_env%subsys)) THEN
               CALL cp_subsys_release(nnp_env%subsys)
            END IF
            IF (ASSOCIATED(nnp_env%subsys)) THEN
               CALL cp_subsys_release(nnp_env%subsys)
            END IF
            IF (ASSOCIATED(nnp_env%cell_ref)) THEN
               CALL cell_release(nnp_env%cell_ref)
            END IF
            DEALLOCATE (nnp_env)
         END IF
      END IF

   END SUBROUTINE nnp_env_release

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

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
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

      CPASSERT(ASSOCIATED(nnp_env))
      CPASSERT(nnp_env%ref_count > 0)

      IF (PRESENT(nnp_potential_energy)) THEN
         nnp_potential_energy = nnp_env%nnp_potential_energy
      END IF
      IF (PRESENT(nnp_forces)) nnp_forces = nnp_env%nnp_forces

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
                          nnp_input, force_env_input, cell_ref, &
                          use_ref_cell, nnp_potential_energy)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
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
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_ref
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_ref_cell
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: nnp_potential_energy

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      CPASSERT(ASSOCIATED(nnp_env))
      CPASSERT(nnp_env%ref_count > 0)

      IF (PRESENT(nnp_potential_energy)) THEN
         nnp_env%nnp_potential_energy = nnp_potential_energy
      END IF
      IF (PRESENT(nnp_forces)) nnp_env%nnp_forces(:, :) = nnp_forces

      IF (PRESENT(subsys)) THEN
         CALL cp_subsys_retain(subsys)
         CALL cp_subsys_release(nnp_env%subsys)
         nnp_env%subsys => subsys
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
