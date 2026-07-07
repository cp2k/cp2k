!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the LibGint-Library.
!> \par History
!>      10.2024 Created
!> \author Marcello Puligheddu
! **************************************************************************************************
MODULE libGint_wrapper

   USE kinds, ONLY: dp
#if(__LIBGINT)
   USE input_constants, ONLY: do_potential_coulomb, do_potential_truncated
   USE libGint, ONLY: libgint_init, libgint_set_Potential_Truncated, libgint_set_hf_fac, libgint_set_max_mem, &
                      libgint_set_P, libgint_set_P_polarized, libgint_set_K, libgint_set_K_polarized, &
                      libgint_get_K, libgint_get_K_polarized, libgint_set_Atom, libgint_set_Atom_L, &
                      libgint_set_cell, libgint_set_neighs, &
                      libgint_add_prm, libgint_add_shell, libgint_add_cell, libgint_add_qrt, &
                      libgint_add_qrtt, libgint_add_set
   USE t_c_g0, ONLY: C0
#endif
   USE hfx_types, ONLY: hfx_type, hfx_memory_type, hfx_potential_type, &
                        hfx_screen_coeff_type, hfx_cell_type, hfx_basis_type

   USE cell_types, ONLY: cell_type
   USE hfx_pair_list_methods, ONLY: build_pair_list_pbc_pgf, bra_t, allocate_bra
   USE particle_types, ONLY: particle_type

   USE iso_C_binding, ONLY: c_ptr

#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

#if(__LIBGINT)
   INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: first_set_of_atom
   TYPE(bra_t), TARGET, SAVE :: bra, ket
   LOGICAL, SAVE :: first_call = .TRUE.
   TYPE(c_ptr), SAVE :: libGint_handle
   !$OMP THREADPRIVATE( first_set_of_atom, first_call )
   !$OMP THREADPRIVATE( bra,ket )
   !$OMP THREADPRIVATE( libGint_handle )
#endif

   PUBLIC :: cp_libGint_init, libGint_update_env, libGint_set_density, libGint_coulomb4, &
             libGint_update_fock_matrix, libGint_get_fock_matrix
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libGint_wrapper'

! Comunicates the current density to the libGint engine.
! Sets the Fock matrix on the device to zero
   INTERFACE libGint_set_density
      MODULE PROCEDURE libGint_set_density_A
      MODULE PROCEDURE libGint_set_density_AB
   END INTERFACE

   INTERFACE libGint_get_fock_matrix
      MODULE PROCEDURE libGint_get_fock_matrix_A
      MODULE PROCEDURE libGint_get_fock_matrix_AB
   END INTERFACE

CONTAINS
! **************************************************************************************************
!> \brief Sets libGint internal enviroment, must be called at least once before libGint can be used
!> \param[in] actual_x_data    pointer to hfx_data
! **************************************************************************************************
   SUBROUTINE cp_libGint_init(actual_x_data)

      TYPE(hfx_type), POINTER, INTENT(in) :: actual_x_data
#if(__LIBGINT)
      ! Init the offload library, creates an handle unique to this omp thread
      CALL libgint_init(libGint_handle)
      ! Comunicate the chosen potential and its parameters to libGint
      IF (actual_x_data%potential_parameter%potential_type == do_potential_truncated) THEN
         ! truncated coulomb needs the C0 coefficients. We do not read or compute them,
         ! they must be already saved in C0 from the t[runcated]_c[oulomb]_g0 module
         CALL libgint_set_Potential_Truncated(libGint_handle, &
                                              actual_x_data%potential_parameter%cutoff_radius, &
                                              C0(:, :))
      ELSE
         CPABORT("The selected interaction potential type is not available with libGint")
      END IF
      first_call = .FALSE.
#else
      MARK_USED(actual_x_data)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE cp_libGint_init
! **************************************************************************************************
!> \brief Initialize and update the libGint computational environment. Must be called at least once
!>        after geo_change and before libGint can be used.
!>
!> This routine sets up data required by the libGint integral engine, including
!> Hartree–Fock scaling factors, memory limits, periodic cell information, and
!> per-atom/basis-set data. It also allocates CPU-side buffers (`bra` and `ket`)
!> used for storing screened Gaussian primitive pairs.
!>
!> \param[in]  fac              Fraction of exact exchange
!> \param[in]  memory_parameter Pointer to memory configuration
!> \param[in]  do_periodic      Logical flag: whether to consider pbc
!> \param[in]  cell             primitive simulation cell
!> \param[in]  actual_x_data    HF exchange data
!> \param[in]  nneighbors       Lattice points
!> \param[in]  max_pgf          Maximum number of primitive Gaussians per shel
!> \param[in]  natom            Number of atoms in the system
!> \param[in]  kind_of          Array mapping atom indices to atomic kinds
!> \param[in]  particle_set     particle set, we extract the positions
!> \param[in]  basis_parameter  gaussian basis set parameters
!>
! **************************************************************************************************
   SUBROUTINE libGint_update_env(fac, memory_parameter, do_periodic, cell, actual_x_data, nneighbors, max_pgf, &
                                 natom, kind_of, particle_set, basis_parameter)

      REAL(dp) :: fac
      TYPE(hfx_memory_type), POINTER :: memory_parameter
      LOGICAL :: do_periodic
      TYPE(cell_type), POINTER :: cell
      TYPE(hfx_type), POINTER :: actual_x_data
      INTEGER :: nneighbors, max_pgf, natom
      INTEGER, ALLOCATABLE, DIMENSION(:) :: kind_of
      TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
      TYPE(hfx_basis_type), DIMENSION(:), POINTER :: basis_parameter
#if(__LIBGINT)
      LOGICAL(1) :: do_pbc
      REAL(dp), DIMENSION(:, :), ALLOCATABLE :: cell_r
      REAL(dp), DIMENSION(:), ALLOCATABLE :: flat_gcc
      INTEGER, DIMENSION(:), POINTER :: la_min, la_max, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER :: nsgfl_a
      INTEGER :: i, l, iset, jset, iatom, ikind, nseta, inla, nla
      REAL(dp) :: ra(3)
      REAL(dp), DIMENSION(:, :), POINTER :: zeta
      REAL(dp), DIMENSION(:, :, :), POINTER :: gcc
      TYPE(bra_t), POINTER                               :: bra_p, ket_p

      ! Set the multiplicative factor fac (the fraction of exact exchange times the spin factor) for libGint
      CALL libGint_set_hf_fac(libGint_handle, fac)
      ! Comunicate max gpu mem per mpi
      CALL libGint_set_max_mem(libGint_handle, memory_parameter%max_memory)

      ! Info about periodic cells and neighbouring cells
      do_pbc = do_periodic
      ALLOCATE (cell_r(3, nneighbors))
      DO i = 1, nneighbors
         cell_r(:, i) = actual_x_data%neighbor_cells(i)%cell_r(:)
      END DO
      CALL libgint_set_cell(libGint_handle, do_pbc, cell%hmat, cell%h_inv)
      CALL libgint_set_neighs(libGint_handle, cell_r, nneighbors)

      ! CPU side temporary arrays for info about the <AB(g) and CD(n)> pairs
      bra_p => bra
      ket_p => ket
      CALL allocate_bra(bra_p, max_pgf, nneighbors)
      CALL allocate_bra(ket_p, max_pgf, nneighbors)

      ! Comunicate atomset and atom info to the engine
      jset = 1
      IF (ALLOCATED(first_set_of_atom)) DEALLOCATE (first_set_of_atom)
      ALLOCATE (first_set_of_atom(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         ra = particle_set(iatom)%r(:)
         nseta = basis_parameter(ikind)%nset
         npgfa => basis_parameter(ikind)%npgf
         la_min => basis_parameter(ikind)%lmin
         la_max => basis_parameter(ikind)%lmax
         zeta => basis_parameter(ikind)%zet
         nsgfa => basis_parameter(ikind)%nsgf
         nsgfl_a => basis_parameter(ikind)%nsgfl
         gcc => basis_parameter(ikind)%gcc
         first_set_of_atom(iatom) = jset
         DO iset = 1, nseta
            CALL libgint_set_Atom(libGint_handle, jset - 1, ra, zeta(:, iset), npgfa(iset))
            inla = 1
            DO l = la_min(iset), la_max(iset)
               nla = nsgfl_a(l, iset)
               IF (ALLOCATED(flat_gcc)) DEALLOCATE (flat_gcc)
               ALLOCATE (flat_gcc(npgfa(iset)*nla))
               flat_gcc(:) = PACK(gcc(1:npgfa(iset), inla:inla + nla - 1, iset), .TRUE.)
               CALL libgint_set_Atom_L(libGint_handle, jset - 1, l, nla, flat_gcc)
               inla = inla + nla
            END DO
            jset = jset + 1
         END DO
      END DO
#else
      MARK_USED(fac)
      MARK_USED(memory_parameter)
      MARK_USED(do_periodic)
      MARK_USED(cell)
      MARK_USED(actual_x_data)
      MARK_USED(nneighbors)
      MARK_USED(max_pgf)
      MARK_USED(natom)
      MARK_USED(kind_of)
      MARK_USED(particle_set)
      MARK_USED(basis_parameter)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_update_env

! **************************************************************************************************
!> \brief communicates the density to libGint, no spin case
!>
!> \param[in]  full_density_alpha       full density matrix array
! **************************************************************************************************
   SUBROUTINE libGint_set_density_A(full_density_alpha)

      REAL(dp), DIMENSION(:, :), POINTER :: full_density_alpha
#if(__LIBGINT)
      CALL libgint_set_P(libGint_handle, full_density_alpha(:, 1))
#else
      MARK_USED(full_density_alpha)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_set_density_A

! **************************************************************************************************
!> \brief communicates the density to libGint, spin case
!>
!> \param[in]  full_density_alpha       full density matrix array on the alpha channel
!> \param[in]  full_density_beta        full density matrix array on the beta  channel
! **************************************************************************************************
   SUBROUTINE libGint_set_density_AB(full_density_alpha, full_density_beta)

      REAL(dp), DIMENSION(:, :), POINTER :: full_density_alpha, full_density_beta
#if(__LIBGINT)
      CALL libgint_set_P_polarized(libGint_handle, full_density_alpha(:, 1), full_density_beta(:, 1))
#else
      MARK_USED(full_density_alpha)
      MARK_USED(full_density_beta)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_set_density_AB

! **************************************************************************************************
!> \brief requests the Fock matrix ( hf_fraction * D @@ I ) from libGint. When this function
!>        returns, full_ks_alpha_from_gpu will contain the fock matrix for this MPI rank.
!>
!> \param[out]  full_ks_alpha_from_gpu     location, already allocated, for the Fock matrix
!> \note assumes full_ks_alpha_from_gpu is already allocated with enough memory
! **************************************************************************************************
   SUBROUTINE libGint_get_fock_matrix_A(full_ks_alpha_from_gpu)

      REAL(dp), DIMENSION(:, :), POINTER :: full_ks_alpha_from_gpu
#if(__LIBGINT)
      CALL libGint_get_K(libGint_handle, full_ks_alpha_from_gpu(:, 1))
#else
      MARK_USED(full_ks_alpha_from_gpu)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_get_fock_matrix_A

! **************************************************************************************************
!> \brief requests the Fock matrices from libGint for both channels. When this function
!>        returns, full_ks_alpha_from_gpu and beta will contain the fock matrix for this MPI rank.
!>
!> \param[out]  full_ks_alpha_from_gpu     location, already allocated, for the Fock matrix alpha
!> \param[out]  full_ks_beta_from_gpu      location, already allocated, for the Fock matrix beta
!> \note assumes full_ks_alpha_from_gpu and full_ks_beta_from_gpu are already allocated
! **************************************************************************************************
   SUBROUTINE libGint_get_fock_matrix_AB(full_ks_alpha_from_gpu, full_ks_beta_from_gpu)

      REAL(dp), DIMENSION(:, :), POINTER :: full_ks_alpha_from_gpu, full_ks_beta_from_gpu
#if(__LIBGINT)
      CALL libgint_get_K_polarized(libGint_handle, full_ks_alpha_from_gpu(:, 1), full_ks_beta_from_gpu(:, 1))
#else
      MARK_USED(full_ks_alpha_from_gpu)
      MARK_USED(full_ks_beta_from_gpu)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_get_fock_matrix_AB

! **************************************************************************************************
!> \brief Assign two-electron integrals of a quartet/shell to libGint
!>
!> First we build a list of ab primitives and cd primitives using build_pair_list_pbc_pgf
!> These lists will contain (before screening) n_cell * npgf * npgf.
!> We loop over both and send ab(g) | cd(h) to libgint, which will loop internally over n
!>
!> \param[in]  iatom                index of atom I / A
!> \param[in]  jatom                index of atom J / B
!> \param[in]  katom                index of atom K / C
!> \param[in]  latom                index of atom L / D
!> \param[in]  iset                 index of the set for atom iatom we are adding
!> \param[in]  jset                 index of the set for atom jatom we are adding
!> \param[in]  kset                 index of the set for atom katom we are adding
!> \param[in]  lset                 index of the set for atom latom we are adding
!> \param[in]  ra                   position of atom iatom
!> \param[in]  rb                   position of atom jatom
!> \param[in]  rc                   position of atom katom
!> \param[in]  rd                   position of atom latom
!> \param[in]  npgfa                number of gaussians in set iset of atom iatom
!> \param[in]  npgfb                number of gaussians in set jset of atom jatom
!> \param[in]  npgfc                number of gaussians in set kset of atom katom
!> \param[in]  npgfd                number of gaussians in set lset of atom latom
!> \param[in]  potential_parameter  information about the potential
!> \param[in]  screen1              screening information for AB pair
!> \param[in]  screen2              screening information for CD pair
!> \param[in]  log10_pmax           density screening factor
!> \param[in]  log10_eps_schwarz    screening tolerance
!> \param[in]  pgf1                 screening information for each AB pair primitive
!> \param[in]  pgf2                 screening information for each CD pair primitive
!> \param[in]  neighbor_cells       array with lattice vectors
!> \param[in]  cell                 information about the simulation box
!> \param[in]  do_periodic          flag for pbc
!> \param[out] screened             whether the whole quartet was screened out
! **************************************************************************************************

   SUBROUTINE libGint_coulomb4(iatom, jatom, katom, latom, iset, jset, kset, lset, &
                               ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
                               potential_parameter, &
                               screen1, screen2, log10_pmax, log10_eps_schwarz, &
                               pgf1, pgf2, &
                               neighbor_cells, cell, do_periodic, screened)

      INTEGER, INTENT(in) :: iatom, jatom, katom, latom, iset, jset, kset, lset
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), rc(3), rd(3)
      INTEGER, INTENT(IN)                                :: npgfa, npgfb, npgfc, npgfd
      TYPE(hfx_potential_type)                           :: potential_parameter
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2)
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: pgf1, pgf2
      TYPE(hfx_cell_type), DIMENSION(:), POINTER         :: neighbor_cells
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: do_periodic
      LOGICAL, INTENT(out) :: screened

#if(__LIBGINT)
      TYPE(bra_t), POINTER                               :: bra_p, ket_p
      LOGICAL :: cell_was_screened
      INTEGER  :: idx_n1, idx_n2, n1, n2, idx_ij, idx_kl, o_ij, n_ij, o_kl, n_kl !, n3
      INTEGER  :: ipgf, jpgf, kpgf, lpgf, iatom_set, jatom_set, katom_set, latom_set
      REAL(dp) :: pgf_max_1, pgf_max_2 ! , R1, R2, rpq2
      INTEGER :: nelements_ij, nelements_kl

      cell = cell
      potential_parameter = potential_parameter
      bra_p => bra
      ket_p => ket

      screened = .TRUE.
      iatom_set = first_set_of_atom(iatom) + iset - 2
      jatom_set = first_set_of_atom(jatom) + jset - 2
      katom_set = first_set_of_atom(katom) + kset - 2
      latom_set = first_set_of_atom(latom) + lset - 2

      CALL build_pair_list_pbc_pgf(npgfa, npgfb, bra_p, screen1, screen2, &
                                   pgf1, log10_pmax, log10_eps_schwarz, ra, rb, &
                                   nelements_ij, neighbor_cells, do_periodic)

      CALL build_pair_list_pbc_pgf(npgfc, npgfd, ket_p, screen2, screen1, &
                                   pgf2, log10_pmax, log10_eps_schwarz, rc, rd, &
                                   nelements_kl, neighbor_cells, do_periodic)

      ! Note: we use 3 numbers n1 n2 and n3 as indices for the lattice traslantion vectors
      ! n1 for the AB pair, n2 for the CD pair and n3 for the PQ pair
      ! so that e.g. B = B0 + ua(n1) means B.x = B0.x + ua(n1).x (and same for y and z)
      ! the ua, saved in neighbor_cells(:)%cell_r(:), are already computed
      ! lattice translation vectors T = i a1 + j a2 + k a3
      ! where a1,a2 and a3 are lattice vectors and i,j and k integers.
      ! So, B.x = B0.x + n1.i * a1.x + n1.j * a2.x + n1.k * a3.x (and same for y and z)
      !

      DO idx_n1 = 1, bra%cell_cnt

         n1 = bra%cell_idx(1, idx_n1)
         n_ij = bra%cell_idx(2, idx_n1)
         o_ij = bra%cell_idx(3, idx_n1)

         DO idx_n2 = 1, ket%cell_cnt

            n2 = ket%cell_idx(1, idx_n2)
            n_kl = ket%cell_idx(2, idx_n2)
            o_kl = ket%cell_idx(3, idx_n2)

            cell_was_screened = .TRUE.
            DO idx_ij = o_ij + 1, o_ij + n_ij

               ipgf = bra%pgf_idx(1, idx_ij)
               jpgf = bra%pgf_idx(2, idx_ij)

               pgf_max_1 = bra%pgf_scr(1, idx_ij)

               DO idx_kl = o_kl + 1, o_kl + n_kl
                  kpgf = ket%pgf_idx(1, idx_kl)
                  lpgf = ket%pgf_idx(2, idx_kl)

                  pgf_max_2 = ket%pgf_scr(1, idx_kl)

                  IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz) CYCLE

                  CALL libGint_add_prm(libGint_handle, ipgf - 1, jpgf - 1, kpgf - 1, lpgf - 1)
                  cell_was_screened = .FALSE.

               END DO ! ket pgf
            END DO ! bra pgf

            IF (.NOT. cell_was_screened) THEN
               CALL libgint_add_shell(libGint_handle, iatom_set, jatom_set, katom_set, latom_set, n1 - 1, n2 - 1)
               screened = .FALSE.
            END IF

         END DO ! ket n2
      END DO ! bra n1
#else
      MARK_USED(iatom)
      MARK_USED(jatom)
      MARK_USED(katom)
      MARK_USED(latom)
      MARK_USED(iset)
      MARK_USED(jset)
      MARK_USED(kset)
      MARK_USED(lset)
      MARK_USED(ra)
      MARK_USED(rb)
      MARK_USED(rc)
      MARK_USED(rd)
      MARK_USED(npgfa)
      MARK_USED(npgfb)
      MARK_USED(npgfc)
      MARK_USED(npgfd)
      MARK_USED(potential_parameter)
      MARK_USED(screen1)
      MARK_USED(screen2)
      MARK_USED(log10_pmax)
      MARK_USED(log10_eps_schwarz)
      MARK_USED(pgf1)
      MARK_USED(pgf2)
      MARK_USED(neighbor_cells)
      MARK_USED(cell)
      MARK_USED(do_periodic)
      MARK_USED(screened)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_coulomb4

! **************************************************************************************************
!> \brief The previous coulomb_4 function assigned an integral between primitive gaussian.
!>        Now we assign the nsgfa(iset)*b*c*d gcc integrals form this set,
!>        along with the Pac Pad Pbc Pbd density to the Kbd Kbc Kad Kac Fock matrix
!>
!> \param[in] symm_fac          simmetry factor from iatom=jatom and the like
!> \param[in] iatom             index of atom A
!> \param[in] jatom             index of atom B
!> \param[in] katom             index of atom C
!> \param[in] latom             index of atom D
!> \param[in] iset              index of set for atom A
!> \param[in] jset              index of set for atom B
!> \param[in] kset              index of set for atom C
!> \param[in] lset              index of set for atom D
!> \param[in] atomic_offset_ac  global offset for the pair of A and C atom in the density ( and Fock) matrix
!> \param[in] atomic_offset_ad  global offset for the pair of A and D atom in the density ( and Fock) matrix
!> \param[in] atomic_offset_bc  global offset for the pair of B and C atom in the density ( and Fock) matrix
!> \param[in] atomic_offset_bd  global offset for the pair of B and D atom in the density ( and Fock) matrix
!> \param[in] offset_ac_set     matrix of sub_offset for sets in atomic_offset_ac
!> \param[in] offset_ad_set     matrix of sub_offset for sets in atomic_offset_ad
!> \param[in] offset_bc_set     matrix of sub_offset for sets in atomic_offset_bc
!> \param[in] offset_bd_set     matrix of sub_offset for sets in atomic_offset_bd
!> \param[in] nsgfa             total number of (spherical, contracted) integrals for iset.
!>                              Used as leading dimension of the AC, AD sublock, depending on transposition
!> \param[in] nsgfb             total number of (spherical, contracted) integrals for jset.
!>                              Used as leading dimension of the BC, BD sublock, depending on transposition
!> \param[in] nsgfc             total number of (spherical, contracted) integrals for kset.
!>                              Used as leading dimension of the AC, BC sublock, depending on transposition
!> \param[in] nsgfd             total number of (spherical, contracted) integrals for lset.
!>                              Used as leading dimension of the AD, BD sublock, depending on transposition
!> \param[in] la_min            minumum angular moment in set iset of atom A
!> \param[in] la_max            maximum angular moment in set iset of atom A
!> \param[in] lb_min            minumum angular moment in set iset of atom B
!> \param[in] lb_max            maximum angular moment in set iset of atom B
!> \param[in] lc_min            minumum angular moment in set iset of atom C
!> \param[in] lc_max            maximum angular moment in set iset of atom C
!> \param[in] ld_min            minumum angular moment in set iset of atom D
!> \param[in] ld_max            maximum angular moment in set iset of atom D
!> \param[in] nsgfl_a           matrix with the number of linear combinations of primitive gaussians
!>                              for each angular moment for each set in atom A. We only read iset
!> \param[in] nsgfl_b           matrix with the number of linear combinations of primitive gaussians
!>                              for each angular moment for each set in atom B. We only read jset
!> \param[in] nsgfl_c           matrix with the number of linear combinations of primitive gaussians
!>                              for each angular moment for each set in atom C. We only read kset
!> \param[in] nsgfl_d           matrix with the number of linear combinations of primitive gaussians
!>                              for each angular moment for each set in atom D. We only read lset
!> \note
!>  The atomic_offset_xy, offset_xy_set matrices provide set offsets which are combined with per-L and
!>    per linear combination offsets to produce the final indices into the dense (but block-sparse)
!>    density and Fock matrices
!>  - The routine assumes Fortran column major order and contiguous storage for the per set
!>    density subblocks as described in the code comments.
!>  - The code computes transposition flags to only use the lower part of P and K
!>
!>
! **************************************************************************************************

   SUBROUTINE libGint_update_fock_matrix( &
      symm_fac, &
      iatom, jatom, katom, latom, &
      iset, jset, kset, lset, &
      atomic_offset_ac, atomic_offset_ad, atomic_offset_bc, atomic_offset_bd, &
      offset_ac_set, offset_ad_set, offset_bc_set, offset_bd_set, &
      nsgfa, nsgfb, nsgfc, nsgfd, &
      la_min, la_max, lb_min, lb_max, &
      lc_min, lc_max, ld_min, ld_max, &
      nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d)

      REAL(dp) :: symm_fac
      INTEGER :: iatom, jatom, katom, latom
      INTEGER :: iset, jset, kset, lset
      INTEGER :: atomic_offset_ac, atomic_offset_ad, atomic_offset_bc, atomic_offset_bd
      INTEGER, DIMENSION(:, :), POINTER :: offset_ac_set, offset_ad_set
      INTEGER, DIMENSION(:, :), POINTER :: offset_bc_set, offset_bd_set
      INTEGER :: nsgfa, nsgfb, nsgfc, nsgfd
      INTEGER :: la_min, la_max, lb_min, lb_max
      INTEGER :: lc_min, lc_max, ld_min, ld_max
      INTEGER, DIMENSION(:, :), POINTER :: nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d

#if(__LIBGINT)
      !! (Hyp)
      ! Let a be a set composed of 2 s and 1 p function.
      ! Let c be a set composed of 1 s and 2 p function.
      ! (1) The density matrix for the ac pair is a 5 x 7 matrix organized as
      !
      !    / -------------------------------------------------------------------------------------------------------------------\
      !    |   a_s1_0@c_s1_0  || a_s1_0@c_p1_0 | a_s1_0@c_p1_1 | a_s1_0@c_p1_2 || a_s1_0@c_p2_0 | a_s1_0@c_p2_1 | a_s1_0@c_p2_2 |
      !    |   a_s2_0@c_s1_0  || a_s2_0@c_p1_0 | a_s2_0@c_p1_1 | a_s2_0@c_p1_2 || a_s2_0@c_p2_0 | a_s2_0@c_p2_1 | a_s2_0@c_p2_2 |
      !    |   a_p1_0@c_s1_0  || a_p1_0@c_p1_0 | a_p1_0@c_p1_1 | a_p1_0@c_p1_2 || a_p1_0@c_p2_0 | a_p1_0@c_p2_1 | a_p1_0@c_p2_2 |
      !    |   a_p1_1@c_s1_0  || a_p1_1@c_p1_0 | a_p1_1@c_p1_1 | a_p1_1@c_p1_2 || a_p1_1@c_p2_0 | a_p1_1@c_p2_1 | a_p1_1@c_p2_2 |
      !    |   a_p1_2@c_s1_0  || a_p1_2@c_p1_0 | a_p1_2@c_p1_1 | a_p1_2@c_p1_2 || a_p1_2@c_p2_0 | a_p1_2@c_p2_1 | a_p1_2@c_p2_2 |
      !    \ -------------------------------------------------------------------------------------------------------------------/
      !
      ! where A_LX_Y means the (Y+1) component of the Xth linear combination of the L angular moment for atom A
      !
      ! (2) This matrix is dense, rectangular and contigous in memory, in fortran column major order.
      ! (3) The big matrix with all pairs is block sparse triangular, only the lower part is valid.

      LOGICAL(1) :: Tac, Tad, Tbc, Tbd
      INTEGER :: offset_ac_L_set, offset_ad_L_set, offset_bc_L_set, offset_bd_L_set
      INTEGER :: s_offset_a, s_offset_b, s_offset_c, s_offset_d
      INTEGER :: s_offset_a_l, s_offset_b_l, s_offset_c_l, s_offset_d_l
      INTEGER :: s_offset_ac, s_offset_ad, s_offset_bc, s_offset_bd
      INTEGER :: ld_ac_set, ld_ad_set, ld_bc_set, ld_bd_set
      INTEGER :: la, lb, lc, ld, nla, nlb, nlc, nld, inla, inlb, inlc, inld
      ! TODO rewrite as update_fock_matrix_gpu(libGint_handle, iatomset,jatomset,katomset,latomset )
      ! AFTER TODO communicate (the pointer to) atomic_offset to libGint AND
      ! AFTER TODO communicate (the pointer to)    set_offset to libGint AND
      ! AFTER TODO check if this idea makes sense in general for other codes
      !
      ! Note: this would not change the need to compute sub offsets and
      ! the 8 loops, it would just transfer them to libGint
      ! Except, if libGint can be sure every set has 1  l, it can collapse the l loops
      ! and/or, if libGint can be sure every l   has 1 nl, it can collapse the n loops
      IF (jatom >= latom) THEN
         offset_bd_L_set = offset_bd_set(jset, lset) + atomic_offset_bd - 2
         ld_bd_set = nsgfb
         Tbd = .FALSE.
      ELSE
         offset_bd_L_set = offset_bd_set(lset, jset) + atomic_offset_bd - 2
         ld_bd_set = nsgfd
         Tbd = .TRUE.
      END IF
      IF (jatom >= katom) THEN
         offset_bc_L_set = offset_bc_set(jset, kset) + atomic_offset_bc - 2
         ld_bc_set = nsgfb
         Tbc = .FALSE.
      ELSE
         offset_bc_L_set = offset_bc_set(kset, jset) + atomic_offset_bc - 2
         ld_bc_set = nsgfc
         Tbc = .TRUE.
      END IF

      IF (iatom >= latom) THEN
         offset_ad_L_set = offset_ad_set(iset, lset) + atomic_offset_ad - 2
         ld_ad_set = nsgfa
         Tad = .FALSE.
      ELSE
         offset_ad_L_set = offset_ad_set(lset, iset) + atomic_offset_ad - 2
         ld_ad_set = nsgfd
         Tad = .TRUE.
      END IF

      IF (iatom >= katom) THEN
         offset_ac_L_set = offset_ac_set(iset, kset) + atomic_offset_ac - 2
         ld_ac_set = nsgfa
         Tac = .FALSE.
      ELSE
         offset_ac_L_set = offset_ac_set(kset, iset) + atomic_offset_ac - 2
         ld_ac_set = nsgfc
         Tac = .TRUE.
      END IF

      s_offset_a_l = 0
      DO la = la_min, la_max
         nla = nsgfl_a(la, iset)
         s_offset_b_l = 0
         DO lb = lb_min, lb_max
            nlb = nsgfl_b(lb, jset)
            s_offset_c_l = 0
            DO lc = lc_min, lc_max
               nlc = nsgfl_c(lc, kset)
               s_offset_d_l = 0
               ld_loop: DO ld = ld_min, ld_max
                  nld = nsgfl_d(ld, lset)
                  CALL libgint_add_qrt(libGint_handle, la, lb, lc, ld, nla, nlb, nlc, nld)
                  DO inla = 1, nla
                     s_offset_a = s_offset_a_l + (inla - 1)*(2*la + 1)
                     DO inlb = 1, nlb
                        s_offset_b = s_offset_b_l + (inlb - 1)*(2*lb + 1)
                        DO inlc = 1, nlc
                           s_offset_c = s_offset_c_l + (inlc - 1)*(2*lc + 1)
                           DO inld = 1, nld
                              s_offset_d = s_offset_d_l + (inld - 1)*(2*ld + 1)
                              IF (.NOT. Tac) THEN
                                 s_offset_ac = offset_ac_L_set + s_offset_c*ld_ac_set + s_offset_a
                              ELSE
                                 s_offset_ac = offset_ac_L_set + s_offset_a*ld_ac_set + s_offset_c
                              END IF

                              IF (.NOT. Tad) THEN
                                 s_offset_ad = offset_ad_L_set + s_offset_d*ld_ad_set + s_offset_a
                              ELSE
                                 s_offset_ad = offset_ad_L_set + s_offset_a*ld_ad_set + s_offset_d
                              END IF

                              IF (.NOT. Tbc) THEN
                                 s_offset_bc = offset_bc_L_set + s_offset_c*ld_bc_set + s_offset_b
                              ELSE
                                 s_offset_bc = offset_bc_L_set + s_offset_b*ld_bc_set + s_offset_c
                              END IF

                              IF (.NOT. Tbd) THEN
                                 s_offset_bd = offset_bd_L_set + s_offset_d*ld_bd_set + s_offset_b
                              ELSE
                                 s_offset_bd = offset_bd_L_set + s_offset_b*ld_bd_set + s_offset_d
                              END IF

                              CALL libgint_add_qrtt(libGint_handle, symm_fac, &
                                                    la, lb, lc, ld, inla - 1, inlb - 1, inlc - 1, inld - 1, &
                                                    ld_ac_set, ld_ad_set, ld_bc_set, ld_bd_set, &
                                                    s_offset_ac, s_offset_ad, s_offset_bc, s_offset_bd, &
                                                    Tac, Tad, Tbc, Tbd)

                           END DO
                        END DO
                     END DO
                  END DO
                  s_offset_d_l = s_offset_d_l + nld*(2*ld + 1)
               END DO ld_loop
               s_offset_c_l = s_offset_c_l + nlc*(2*lc + 1)
            END DO
            s_offset_b_l = s_offset_b_l + nlb*(2*lb + 1)
         END DO
         s_offset_a_l = s_offset_a_l + nla*(2*la + 1)
      END DO

      CALL libgint_add_set(libGint_handle)

#else
      MARK_USED(symm_fac)
      MARK_USED(iatom)
      MARK_USED(jatom)
      MARK_USED(katom)
      MARK_USED(latom)
      MARK_USED(iset)
      MARK_USED(jset)
      MARK_USED(kset)
      MARK_USED(lset)
      MARK_USED(atomic_offset_ac)
      MARK_USED(atomic_offset_ad)
      MARK_USED(atomic_offset_bc)
      MARK_USED(atomic_offset_bd)
      MARK_USED(offset_ac_set)
      MARK_USED(offset_ad_set)
      MARK_USED(offset_bc_set)
      MARK_USED(offset_bd_set)
      MARK_USED(nsgfa)
      MARK_USED(nsgfb)
      MARK_USED(nsgfc)
      MARK_USED(nsgfd)
      MARK_USED(la_min)
      MARK_USED(la_max)
      MARK_USED(lb_min)
      MARK_USED(lb_max)
      MARK_USED(lc_min)
      MARK_USED(lc_max)
      MARK_USED(ld_min)
      MARK_USED(ld_max)
      MARK_USED(nsgfl_a)
      MARK_USED(nsgfl_b)
      MARK_USED(nsgfl_c)
      MARK_USED(nsgfl_d)
      CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   END SUBROUTINE libGint_update_fock_matrix

END MODULE libGint_wrapper
