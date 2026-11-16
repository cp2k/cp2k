!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for optimizing load balance between processes in HFX calculations
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      11.2019 fixed initial value for potential_id (A. Bussy)
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_pair_list_methods
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE four_center_integral_methods,    ONLY: min_image_gap_sq
   USE hfx_types,                       ONLY: hfx_basis_type,&
                                              hfx_block_range_type,&
                                              hfx_cell_type,&
                                              hfx_screen_coeff_type,&
                                              pair_list_type,&
                                              pair_set_list_type
   USE kinds,                           ONLY: dp
   USE mp2_types,                       ONLY: pair_list_type_mp2
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: build_pair_list, &
             build_pair_list_mp2, &
             build_atomic_pair_list, &
             build_pair_list_pbc_pgf, &
             bra_t, allocate_bra, min_image_gap_sq

   ! Store screened primitive gaussian function pairs between two shells (A and B)
   TYPE :: Bra_t
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)    :: pgf_scr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)     :: pgf_idx, cell_idx
      INTEGER :: cell_cnt = 0, pgf_cnt = 0
   END TYPE Bra_t

CONTAINS

! **************************************************************************************************
!> \brief Allocates all arrays within a Bra_t structure.
!>
!> This subroutine ensures that a given `Bra_t` pointer has allocated memory
!> for storing primitive Gaussian pair indices, screening data, and cell indices.
!> Existing allocations are safely deallocated first.
!>
!> \param[in,out] bra         Pointer to the Bra_t structure to allocate.
!> \param[in]     max_npgf    Maximum number of primitives per shell.
!> \param[in]     max_ncell   Maximum number of neighbor cells.
!>
!> \note The arrays are allocated as follows:
!>       - `pgf_idx(3, max_npgf^2 * max_ncell)`
!>       - `pgf_scr(5, max_npgf^2 * max_ncell)`
!>       - `cell_idx(3, max_ncell)`
!>
! **************************************************************************************************
   SUBROUTINE allocate_Bra(bra, max_npgf, max_ncell)
      TYPE(bra_t), POINTER                               :: bra
      INTEGER                                            :: max_npgf, max_ncell

      IF (ALLOCATED(bra%pgf_idx)) DEALLOCATE (bra%pgf_idx)  ! TODO add if n_prm changed
      IF (ALLOCATED(bra%pgf_scr)) DEALLOCATE (bra%pgf_scr)  ! TODO add if n_prm changed
      IF (ALLOCATED(bra%cell_idx)) DEALLOCATE (bra%cell_idx) ! TODO add if n_prm changed
      ALLOCATE (bra%pgf_idx(3, max_npgf*max_npgf*max_ncell))
      ALLOCATE (bra%pgf_scr(5, max_npgf*max_npgf*max_ncell))
      ALLOCATE (bra%cell_idx(3, max_ncell))
   END SUBROUTINE allocate_Bra

! **************************************************************************************************
!> \brief Builds a screened list of primitives from centers A and B, intersecting with another shell
!>
!> This subroutine populates the `Bra_t` structure with indices and screening data
!> for primitive Gaussian pairs (ipgf,jpgf,n3) that pass the Schwarz screening criterion
!>
!> \param[in]  npgfa              Number of primitives in shell A
!> \param[in]  npgfb              Number of primitives in shell B
!> \param[out] bra                Pointer to Bra_t structure to populate
!> \param[in]  screen1            Screening coeffiecients for the AB pair based on most diffuse pgf
!> \param[in]  screen2            Screening coeffiecients for the CD pair based on most diffuse pgf
!> \param[in]  pgf                Screening coeeficents for primitive gaussians
!> \param[in]  log10_pmax         Log10 of maximum integral prefactor
!> \param[in]  log10_eps_schwarz  Log10 of Schwarz screening threshold
!> \param[in]  ra                 Cartesian coordinates of center A
!> \param[in]  rb                 Cartesian coordinates of center B
!> \param[out] nelements          Total number of valid primitives found
!> \param[in]  neighbor_cells     Array of lattice vectors
!> \param[in]  do_periodic        Logical flag for periodic boundary conditions on B
!>
!> \note Loops over all cells, shifting the B shell position accordingly
!>       Only pairs satisfying the screening condition are stored.
!>
!> \note Each valid cell contributes a block of pairs stored contiguously in
!>       bra%pgf_idx, with cell metadata stored in bra%cell_idx
!> \note screening primitive pairs (ipgf, jpgf) using both coarse and fine screening thresholds
! **************************************************************************************************
   SUBROUTINE build_pair_list_pbc_pgf(npgfa, npgfb, bra, screen1, screen2, pgf, &
                                      log10_pmax, log10_eps_schwarz, ra, rb, nelements, &
                                      neighbor_cells, do_periodic)

      INTEGER, INTENT(IN)                                :: npgfa, npgfb
      TYPE(bra_t), POINTER                               :: bra
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2)
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: pgf
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz, ra(3), &
                                                            rb(3)
      INTEGER, INTENT(OUT)                               :: nelements
      TYPE(hfx_cell_type), DIMENSION(:), POINTER         :: neighbor_cells
      LOGICAL, INTENT(IN)                                :: do_periodic

      INTEGER                                            :: cell_cnt, cell_idx, element_cnt, &
                                                            element_idx, element_off, ipgf, jpgf
      REAL(dp)                                           :: AB(3), im_B(3), pgf_max, rab2

      element_idx = 0
      element_off = 0
      cell_cnt = 0
      DO cell_idx = 1, SIZE(neighbor_cells)
         ! move B to this cell while keeping A fixed
         ! NOTE rb has already been moved by build_pair_list
         !      so that when cell_r is (000) AB is in the pbc cell
         IF (do_periodic) THEN
            im_B = rb + neighbor_cells(cell_idx)%cell_r(:)
         ELSE
            im_B = rb
         END IF
         AB = ra - im_B
         rab2 = AB(1)**2 + AB(2)**2 + AB(3)**2
         ! First screening on the most diffuse AB gaussians along with the most diffuse CD gaussian
         IF (screen1(1)*rab2 + screen1(2) + screen2(2) + log10_pmax < log10_eps_schwarz) CYCLE

         element_cnt = 0
         DO ipgf = 1, npgfa
            DO jpgf = 1, npgfb

               ! Second screening on this actual AB pair and the most diffuse CD from before
               pgf_max = pgf(jpgf, ipgf)%x(1)*rab2 + pgf(jpgf, ipgf)%x(2)
               IF (pgf_max + screen2(2) + log10_pmax < log10_eps_schwarz) CYCLE
               ! This pair passed screening. Add it to the bra.
               element_idx = element_idx + 1
               bra%pgf_idx(:, element_idx) = [ipgf, jpgf, cell_idx]
               bra%pgf_scr(1, element_idx) = pgf_max
               element_cnt = element_cnt + 1
            END DO
         END DO

         ! If this cell produced any pair, add this cell to the bra
         IF (element_cnt == 0) CYCLE
         cell_cnt = cell_cnt + 1
         bra%cell_idx(:, cell_cnt) = [cell_idx, element_cnt, element_off]
         element_off = element_off + element_cnt

      END DO ! cell
      bra%cell_cnt = cell_cnt
      nelements = element_idx
      bra%pgf_cnt = nelements

   END SUBROUTINE build_pair_list_pbc_pgf

! **************************************************************************************************
!> \brief ...
!> \param natom ...
!> \param list ...
!> \param set_list ...
!> \param i_start ...
!> \param i_end ...
!> \param j_start ...
!> \param j_end ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param particle_set ...
!> \param do_periodic ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param coeffs_kind_max0 ...
!> \param log10_eps_schwarz ...
!> \param cell ...
!> \param pmax_blocks ...
!> \param atomic_pair_list ...
! **************************************************************************************************
   SUBROUTINE build_pair_list(natom, list, set_list, i_start, i_end, j_start, j_end, kind_of, basis_parameter, particle_set, &
                              do_periodic, coeffs_set, coeffs_kind, coeffs_kind_max0, log10_eps_schwarz, cell, &
                              pmax_blocks, atomic_pair_list)

      INTEGER, INTENT(IN)                                :: natom
      TYPE(pair_list_type), INTENT(INOUT)                :: list
      TYPE(pair_set_list_type), DIMENSION(:), &
         INTENT(OUT)                                     :: set_list
      INTEGER, INTENT(IN)                                :: i_start, i_end, j_start, j_end, &
                                                            kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: basis_parameter
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :), INTENT(IN), POINTER      :: coeffs_set
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         INTENT(IN)                                      :: coeffs_kind
      REAL(KIND=dp), INTENT(IN)                          :: coeffs_kind_max0, log10_eps_schwarz
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), INTENT(IN)                               :: pmax_blocks
      LOGICAL, DIMENSION(natom, natom), INTENT(IN)       :: atomic_pair_list

      INTEGER                                            :: iatom, ikind, iset, jatom, jkind, jset, &
                                                            n_element, nset_ij, nseta, nsetb
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), DIMENSION(3)                        :: B11, pbc_B, ra, rb, temp

      n_element = 0
      nset_ij = 0

      DO iatom = i_start, i_end
         DO jatom = j_start, j_end
            IF (atomic_pair_list(jatom, iatom) .EQV. .FALSE.) CYCLE

            ikind = kind_of(iatom)
            nseta = basis_parameter(ikind)%nset
            ra = particle_set(iatom)%r(:)

            IF (jatom < iatom) CYCLE
            jkind = kind_of(jatom)
            nsetb = basis_parameter(jkind)%nset
            rb = particle_set(jatom)%r(:)

            IF (do_periodic) THEN
               temp = rb - ra
               pbc_B = pbc(temp, cell)
               B11 = ra + pbc_B
               rab2 = (ra(1) - B11(1))**2 + (ra(2) - B11(2))**2 + (ra(3) - B11(3))**2
            ELSE
               rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2
               B11 = rb ! ra - rb
            END IF
            IF ((coeffs_kind(jkind, ikind)%x(1)*rab2 + &
                 coeffs_kind(jkind, ikind)%x(2)) + coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz) CYCLE

            n_element = n_element + 1
            list%elements(n_element)%pair = [iatom, jatom]
            list%elements(n_element)%kind_pair = [ikind, jkind]
            list%elements(n_element)%r1 = ra
            list%elements(n_element)%r2 = B11
            list%elements(n_element)%dist2 = rab2
            ! build a list of guaranteed overlapping sets
            list%elements(n_element)%set_bounds(1) = nset_ij + 1
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  IF (coeffs_set(jset, iset, jkind, ikind)%x(1)*rab2 + coeffs_set(jset, iset, jkind, ikind)%x(2) + &
                      coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz) CYCLE
                  nset_ij = nset_ij + 1
                  set_list(nset_ij)%pair = [iset, jset]
               END DO
            END DO
            list%elements(n_element)%set_bounds(2) = nset_ij
         END DO
      END DO

      list%n_element = n_element

   END SUBROUTINE build_pair_list

! **************************************************************************************************
!> \brief ...
!> \param natom ...
!> \param atomic_pair_list ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param particle_set ...
!> \param do_periodic ...
!> \param coeffs_kind ...
!> \param coeffs_kind_max0 ...
!> \param log10_eps_schwarz ...
!> \param cell ...
!> \param blocks ...
! **************************************************************************************************
   SUBROUTINE build_atomic_pair_list(natom, atomic_pair_list, kind_of, basis_parameter, particle_set, &
                                     do_periodic, coeffs_kind, coeffs_kind_max0, log10_eps_schwarz, cell, &
                                     blocks)
      INTEGER, INTENT(IN)                                :: natom
      LOGICAL, DIMENSION(natom, natom)                   :: atomic_pair_list
      INTEGER, INTENT(IN)                                :: kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: basis_parameter
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         INTENT(IN)                                      :: coeffs_kind
      REAL(KIND=dp), INTENT(IN)                          :: coeffs_kind_max0, log10_eps_schwarz
      TYPE(cell_type), POINTER                           :: cell
      TYPE(hfx_block_range_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: blocks

      INTEGER                                            :: iatom, iatom_end, iatom_start, iblock, &
                                                            ikind, jatom, jatom_end, jatom_start, &
                                                            jblock, jkind, nseta, nsetb
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), DIMENSION(3)                        :: B11, pbc_B, ra, rb, temp

      atomic_pair_list = .FALSE.

      DO iblock = 1, SIZE(blocks)
         iatom_start = blocks(iblock)%istart
         iatom_end = blocks(iblock)%iend
         DO jblock = 1, SIZE(blocks)
            jatom_start = blocks(jblock)%istart
            jatom_end = blocks(jblock)%iend

            DO iatom = iatom_start, iatom_end
               ikind = kind_of(iatom)
               nseta = basis_parameter(ikind)%nset
               ra = particle_set(iatom)%r(:)
               DO jatom = jatom_start, jatom_end
                  IF (jatom < iatom) CYCLE
                  jkind = kind_of(jatom)
                  nsetb = basis_parameter(jkind)%nset
                  rb = particle_set(jatom)%r(:)

                  IF (do_periodic) THEN
                     temp = rb - ra
                     pbc_B = pbc(temp, cell)
                     B11 = ra + pbc_B
                     rab2 = (ra(1) - B11(1))**2 + (ra(2) - B11(2))**2 + (ra(3) - B11(3))**2
                  ELSE
                     rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2
                     B11 = rb ! ra - rb
                  END IF
                  IF ((coeffs_kind(jkind, ikind)%x(1)*rab2 + &
                       coeffs_kind(jkind, ikind)%x(2)) + coeffs_kind_max0 < log10_eps_schwarz) CYCLE

                  atomic_pair_list(jatom, iatom) = .TRUE.
                  atomic_pair_list(iatom, jatom) = .TRUE.
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE build_atomic_pair_list

! **************************************************************************************************
!> \brief ...
!> \param natom ...
!> \param list ...
!> \param set_list ...
!> \param i_start ...
!> \param i_end ...
!> \param j_start ...
!> \param j_end ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param particle_set ...
!> \param do_periodic ...
!> \param coeffs_set ...
!> \param coeffs_kind ...
!> \param coeffs_kind_max0 ...
!> \param log10_eps_schwarz ...
!> \param cell ...
!> \param pmax_blocks ...
!> \param atomic_pair_list ...
!> \param skip_atom_symmetry ...
! **************************************************************************************************
   SUBROUTINE build_pair_list_mp2(natom, list, set_list, i_start, i_end, j_start, j_end, kind_of, basis_parameter, particle_set, &
                                  do_periodic, coeffs_set, coeffs_kind, coeffs_kind_max0, log10_eps_schwarz, cell, &
                                  pmax_blocks, atomic_pair_list, skip_atom_symmetry)

      INTEGER, INTENT(IN)                                :: natom
      TYPE(pair_list_type_mp2), INTENT(INOUT)            :: list
      TYPE(pair_set_list_type), DIMENSION(:), &
         INTENT(OUT)                                     :: set_list
      INTEGER, INTENT(IN)                                :: i_start, i_end, j_start, j_end, &
                                                            kind_of(*)
      TYPE(hfx_basis_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: basis_parameter
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :), INTENT(IN), POINTER      :: coeffs_set
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         INTENT(IN)                                      :: coeffs_kind
      REAL(KIND=dp), INTENT(IN)                          :: coeffs_kind_max0, log10_eps_schwarz
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), INTENT(IN)                               :: pmax_blocks
      LOGICAL, DIMENSION(natom, natom), INTENT(IN)       :: atomic_pair_list
      LOGICAL, INTENT(IN), OPTIONAL                      :: skip_atom_symmetry

      INTEGER                                            :: iatom, ikind, iset, jatom, jkind, jset, &
                                                            n_element, nset_ij, nseta, nsetb
      LOGICAL                                            :: my_skip_atom_symmetry
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), DIMENSION(3)                        :: B11, pbc_B, ra, rb, temp

      n_element = 0
      nset_ij = 0

      my_skip_atom_symmetry = .FALSE.
      IF (PRESENT(skip_atom_symmetry)) my_skip_atom_symmetry = skip_atom_symmetry

      DO iatom = i_start, i_end
         DO jatom = j_start, j_end
            IF (atomic_pair_list(jatom, iatom) .EQV. .FALSE.) CYCLE

            ikind = kind_of(iatom)
            nseta = basis_parameter(ikind)%nset
            ra = particle_set(iatom)%r(:)

            IF (jatom < iatom .AND. (.NOT. my_skip_atom_symmetry)) CYCLE
            jkind = kind_of(jatom)
            nsetb = basis_parameter(jkind)%nset
            rb = particle_set(jatom)%r(:)

            IF (do_periodic) THEN
               temp = rb - ra
               pbc_B = pbc(temp, cell)
               B11 = ra + pbc_B
               rab2 = (ra(1) - B11(1))**2 + (ra(2) - B11(2))**2 + (ra(3) - B11(3))**2
            ELSE
               rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2
               B11 = rb ! ra - rb
            END IF
            IF ((coeffs_kind(jkind, ikind)%x(1)*rab2 + &
                 coeffs_kind(jkind, ikind)%x(2)) + coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz) CYCLE

            n_element = n_element + 1
            list%elements(n_element)%pair = [iatom, jatom]
            list%elements(n_element)%kind_pair = [ikind, jkind]
            list%elements(n_element)%r1 = ra
            list%elements(n_element)%r2 = B11
            list%elements(n_element)%dist2 = rab2
            ! build a list of guaranteed overlapping sets
            list%elements(n_element)%set_bounds(1) = nset_ij + 1
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  IF (coeffs_set(jset, iset, jkind, ikind)%x(1)*rab2 + coeffs_set(jset, iset, jkind, ikind)%x(2) + &
                      coeffs_kind_max0 + pmax_blocks < log10_eps_schwarz) CYCLE
                  nset_ij = nset_ij + 1
                  set_list(nset_ij)%pair = [iset, jset]
               END DO
            END DO
            list%elements(n_element)%set_bounds(2) = nset_ij
         END DO
      END DO

      list%n_element = n_element

   END SUBROUTINE build_pair_list_mp2

END MODULE hfx_pair_list_methods
