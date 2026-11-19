!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Band structure unfolding onto a primitive cell for k-point supercell runs
!>        (Quan, Rybin, Scheffler, Carbogno, Phys. Rev. B 113, 085112 (2026)).
!>        Closed-form orbit-Fourier weights
!>        W^k_KN = sum_channels | 1/sqrt(m) sum_{mu in channel} exp(-2 pi i f_k.n_mu) C'_muN |^2
!>        with C' = S(K)^{1/2} C(K) (Loewdin) and displacement-integer labels n_mu, which
!>        live in the L = 0 gauge (no supercell wrap phase needed). The per-band sum rule
!>        sum_q W^{k_q}_KN = C_N^T S(K) C_N = 1 holds strictly for any distortion or defect.
!>        Gauge: atoms are wrapped with the convention of the k-point neighbor
!>        lists. CP2K uses pbc_stable when the SCF ran with k-point symmetry, and
!>        plain pbc otherwise. This matches the rskp_transform anchoring.
! **************************************************************************************************
MODULE kpoint_unfolding

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Quan2026
   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              pbc_stable,&
                                              real_to_scaled
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: lowdin_kp_mo_coeff
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3,&
                                              vector_product
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE reference_manager,               ONLY: cite_reference
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_unfolding'
   INTEGER, PARAMETER, PRIVATE :: msg_length = 256
   ! tolerance for identifying k-points that differ by a reciprocal lattice vector
   REAL(KIND=dp), PARAMETER :: kdifftol = 1.0E-6_dp
   ! bound-inclusion tolerance for unfolded f = M^-T (K + p): coordinates on the
   ! unit-cell edge reach 0 or 1 only within roundoff
   REAL(KIND=dp), PARAMETER :: unit_tol = 1.0E-9_dp
   ! full-atom-selection Loewdin populations equal 1 to the k-point MO orthonormality
   ! accuracy (set by the diagonalization threshold, EPS_DEFAULT-level); violations
   ! indicate an inconsistent S(k) or MO set, well above this guard
   REAL(KIND=dp), PARAMETER :: mo_orthonorm_tol = 1.0E-8_dp
   ! MO block width bounding per-rank channel-reduction scratch in unfolding_kp_weights
   INTEGER, PARAMETER :: unfolding_mo_block = 256

   !> Transformation matrix, derived primitive cell, and the per-AO orbit map
   TYPE :: unfolding_env_type
      INTEGER, DIMENSION(3, 3) :: transmat = 0             ! hmat_SC = hmat_PC . M (column convention)
      INTEGER :: m = 0                                      ! |det(M)|
      REAL(KIND=dp) :: symm_tol = 0.5_dp/angstrom         ! default 0.5 Angstrom in internal Bohr
      REAL(KIND=dp), DIMENSION(3, 3) :: h_pc = 0.0_dp       ! derived primitive lattice (columns)
      INTEGER :: n_sets = 0                                 ! atomic equivalent-site sets
      INTEGER :: n_channels = 0                             ! per-AO equivalence channels
      INTEGER, DIMENSION(:), ALLOCATABLE :: channel_of_ao   ! (nao)
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: n_of_ao      ! (3, nao) primitive translation labels
      INTEGER :: n_incomplete = 0                           ! equivalent-site sets with fewer than m atoms
      INTEGER :: n_groups = 1                               ! atom groups (independent clustering)
      INTEGER, DIMENSION(:), ALLOCATABLE :: channel_group   ! group of each channel (upper bound nao)
      LOGICAL :: full_selection = .TRUE.                    ! no atom screened: Loewdin populations are 1
      LOGICAL :: lowdin_proj = .FALSE.                      ! S(K)^{1/2} prebuilt by the DFT+U projection
   END TYPE unfolding_env_type

   PUBLIC :: unfolding_env_type, unfolding_setup, unfolding_fold_kpath, &
             unfolding_process_batch, unfolding_finish_set, unfolding_env_release

CONTAINS

! **************************************************************************************************
!> \brief Set up band structure unfolding from the UNFOLDING section: primitive cell of
!>        TRANSFORMATION_MATRIX, ATOMS groups, and the supercell-to-primitive orbit map. The
!>        setup is summarized on the output unit and in the header of the band structure file
!> \param uf unfolding environment
!> \param qs_env ...
!> \param u_input UNFOLDING section
!> \param unit_nr default output unit
!> \param bs_data_unit band structure data unit
! **************************************************************************************************
   SUBROUTINE unfolding_setup(uf, qs_env, u_input, unit_nr, bs_data_unit)
      TYPE(unfolding_env_type), INTENT(OUT)              :: uf
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: u_input
      INTEGER, INTENT(IN)                                :: unit_nr, bs_data_unit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'unfolding_setup'

      INTEGER                                            :: g, handle, i, ia, ia_first, ia_last, j, &
                                                            n_atom_groups, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_of_atom
      INTEGER, DIMENSION(3, 3)                           :: transmat
      INTEGER, DIMENSION(:), POINTER                     :: atom_vals, mvals
      LOGICAL                                            :: stable_images
      REAL(KIND=dp)                                      :: symm_tol
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cart_hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL cite_reference(Quan2026)
      NULLIFY (cell, kpoints, particle_set, qs_kind_set)
      CALL get_qs_env(qs_env, cell=cell, kpoints=kpoints, particle_set=particle_set, qs_kind_set=qs_kind_set)

      NULLIFY (mvals)
      CALL section_vals_val_get(u_input, "TRANSFORMATION_MATRIX", i_vals=mvals)
      ! input lists the three matrix rows; RESHAPE fills columns first
      transmat = RESHAPE(mvals, [3, 3], ORDER=[2, 1])
      NULLIFY (mvals)  ! aliases section storage, must not be deallocated
      CALL section_vals_val_get(u_input, "SYMM_TOL", r_val=symm_tol)
      cart_hmat(:, :) = cell%hmat(:, :)
      IF (cell%input_cell_canonicalized) cart_hmat(:, :) = cell%input_hmat(:, :)
      CALL unfolding_env_create(uf, transmat, symm_tol, cart_hmat)

      natom = SIZE(particle_set)
      ALLOCATE (group_of_atom(natom))
      group_of_atom = 0
      CALL section_vals_val_get(u_input, "ATOMS", n_rep_val=n_atom_groups)
      DO g = 1, n_atom_groups
         NULLIFY (atom_vals)
         CALL section_vals_val_get(u_input, "ATOMS", i_rep_val=g, i_vals=atom_vals)
         DO j = 1, SIZE(atom_vals)
            ia = atom_vals(j)
            IF (ia < 1 .OR. ia > natom) THEN
               CPABORT("UNFOLD_BAND| ATOMS index outside 1..natom")
            END IF
            IF (group_of_atom(ia) /= 0) THEN
               CPABORT("UNFOLD_BAND| ATOMS groups must be disjoint")
            END IF
            group_of_atom(ia) = g
         END DO
         NULLIFY (atom_vals)
      END DO
      ! no ATOMS keyword: all atoms form a single implicit group
      IF (n_atom_groups == 0) group_of_atom = 1
      uf%full_selection = ALL(group_of_atom /= 0)
      IF (n_atom_groups > 1 .AND. unit_nr > 0) THEN
         WRITE (unit_nr, FMT="(T2,A,I0,A)") "UNFOLD_BAND| atom groups: ", n_atom_groups, &
            " (weights per group; column 4 = total)"
         DO g = 1, n_atom_groups
            ia_first = 0
            ia_last = 0
            DO j = 1, natom
               IF (group_of_atom(j) == g) THEN
                  IF (ia_first == 0) ia_first = j
                  ia_last = j
               END IF
            END DO
            WRITE (unit_nr, FMT="(T2,A,I0,A,I0,A,I0,A,I0,A)") "UNFOLD_BAND| group ", g, ": ", &
               COUNT(group_of_atom == g), " atoms (first ", ia_first, ", last ", ia_last, ")"
         END DO
      END IF
      IF (n_atom_groups >= 1 .AND. unit_nr > 0) THEN
         WRITE (unit_nr, FMT="(T2,A,I0)") "UNFOLD_BAND| screened atoms (excluded): ", &
            COUNT(group_of_atom == 0)
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, FMT="(T2,A,9(1X,I3))") "UNFOLD_BAND| M (rows) =", TRANSPOSE(uf%transmat)
         WRITE (unit_nr, FMT="(T2,A,I0,A,F8.4)") "UNFOLD_BAND| m = |det M| = ", uf%m, &
            " ; site tolerance [angstrom] = ", uf%symm_tol*angstrom
         DO i = 1, 3
            WRITE (unit_nr, FMT="(T2,A,I0,A,3F12.6)") "UNFOLD_BAND| primitive cell vector ", i, &
               " [Bohr] =", uf%h_pc(1:3, i)
         END DO
      END IF
      ! mirror the neighbor-list anchoring convention (qs_neighbor_lists):
      ! stable images when the SCF ran k-points with symmetry, plain pbc otherwise
      stable_images = .FALSE.
      IF (ASSOCIATED(kpoints)) stable_images = kpoints%nkp > 0 .AND. kpoints%symmetry
      CALL unfolding_env_build_map(uf, cell, particle_set, qs_kind_set, group_of_atom, &
                                   stable_images, unit_nr)

      IF (bs_data_unit > 0 .AND. bs_data_unit /= unit_nr) THEN
         WRITE (UNIT=bs_data_unit, FMT="(A,9(1X,I3))") "# UNFOLD_BAND| M (rows) =", TRANSPOSE(uf%transmat)
         WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,F8.4)") "# UNFOLD_BAND| m = |det M| = ", uf%m, &
            " ; site tolerance [angstrom] = ", uf%symm_tol*angstrom
         DO i = 1, 3
            WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,3F14.8)") "# UNFOLD_BAND| primitive cell vector ", i, &
               " [Bohr] =", uf%h_pc(1:3, i)
         END DO
         WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0)") "# UNFOLD_BAND| equivalent-site sets: ", &
            uf%n_sets, "; incomplete sets: ", uf%n_incomplete
         WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0)") "# UNFOLD_BAND| atom groups: ", uf%n_groups, &
            " ; screened atoms: ", COUNT(group_of_atom == 0)
      END IF
      ! DFT+U inputs (nonzero lowdin_v) build S(K)^{1/2} inside the band structure
      ! diagonalization; rebuilding it for the weights would repeat the eigensolver
      IF (ASSOCIATED(kpoints)) THEN
         IF (ALLOCATED(kpoints%lowdin_v)) uf%lowdin_proj = ANY(kpoints%lowdin_v /= 0.0_dp)
      END IF

      DEALLOCATE (group_of_atom)
      CALL timestop(handle)
   END SUBROUTINE unfolding_setup

! **************************************************************************************************
!> \brief Initialize the unfolding environment from the integer transformation matrix
!> \param uf unfolding environment
!> \param transmat 3x3 integer matrix, hmat_SC = hmat_PC . transmat
!> \param symm_tol site-assignment Cartesian tolerance in Bohr (internal units; the
!>        SYMM_TOL input is given in Angstrom), below half the smallest
!>        primitive-cell lattice plane distance
!> \param cart_hmat supercell cartesian lattice in Bohr (columns are the lattice vectors),
!>        already accounting for input-cell canonicalization
! **************************************************************************************************
   SUBROUTINE unfolding_env_create(uf, transmat, symm_tol, cart_hmat)
      TYPE(unfolding_env_type), INTENT(OUT)              :: uf
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: transmat
      REAL(KIND=dp), INTENT(IN)                          :: symm_tol
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cart_hmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_env_create'

      CHARACTER(LEN=msg_length)                          :: msgstr
      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: plane_d_min
      REAL(KIND=dp), DIMENSION(3)                        :: plane_distance
      REAL(KIND=dp), DIMENSION(3, 3)                     :: m_dp, minv

      CALL timeset(routineN, handle)

      uf%transmat = transmat
      uf%symm_tol = symm_tol
      m_dp = REAL(transmat, KIND=dp)
      uf%m = ABS(NINT(det_3x3(m_dp)))
      IF (uf%m == 0) CPABORT("UNFOLD_BAND| TRANSFORMATION_MATRIX is singular")

      minv = inv_3x3(m_dp)
      DO i = 1, 3
         uf%h_pc(1:3, i) = MATMUL(cart_hmat(1:3, 1:3), minv(1:3, i))
      END DO

      ! plane distance of each cyclic pair of primitive lattice-vector columns:
      ! cell volume divided by the cross-product area
      plane_distance(1) = plane_dist_of(uf%h_pc, 2, 3)
      plane_distance(2) = plane_dist_of(uf%h_pc, 3, 1)
      plane_distance(3) = plane_dist_of(uf%h_pc, 1, 2)

      ! a tolerance below half the smallest primitive plane distance guarantees
      ! that componentwise rounding of ds recovers the nearest ideal site. It
      ! also rules out two ideal sites within tol of one atom
      IF (symm_tol <= 0.0_dp) THEN
         CPABORT("UNFOLD_BAND| SYMM_TOL must be positive")
      END IF
      plane_d_min = MIN(plane_distance(1), plane_distance(2), plane_distance(3))
      IF (symm_tol >= 0.5_dp*plane_d_min) THEN
         WRITE (UNIT=msgstr, FMT="(A,F8.3,A)") &
            "UNFOLD_BAND| SYMM_TOL must be smaller than ", &
            0.5_dp*plane_d_min*angstrom, " angstrom (half the smallest primitive plane distance)"
         CPABORT(TRIM(msgstr))
      END IF

      CALL timestop(handle)
   END SUBROUTINE unfolding_env_create

! **************************************************************************************************
!> \brief Release the unfolding environment
!> \param uf ...
! **************************************************************************************************
   SUBROUTINE unfolding_env_release(uf)
      TYPE(unfolding_env_type), INTENT(INOUT)            :: uf

      IF (ALLOCATED(uf%channel_of_ao)) DEALLOCATE (uf%channel_of_ao)
      IF (ALLOCATED(uf%n_of_ao)) DEALLOCATE (uf%n_of_ao)
      IF (ALLOCATED(uf%channel_group)) DEALLOCATE (uf%channel_group)
      uf%n_sets = 0
      uf%n_channels = 0
      uf%n_incomplete = 0
      uf%m = 0
   END SUBROUTINE unfolding_env_release

! **************************************************************************************************
!> \brief Cluster supercell atoms into primitive-translation orbits (per-AO channels).
!>        Candidates are matched by element (atomic number; GHOST kinds equal their
!>        base element), and inside an orbit the basis set must be identical (nsgf
!>        and the flattened set->shell l sequence) - violations abort. Each atom
!>        joins the same-group, same-kind orbit whose ideal site lies within SYMM_TOL
!>        (Cartesian distance); a vacancy is a missing member (incomplete channel), any
!>        further primitive-cell basis atom or interstitial opens its own orbit. Slots
!>        are translation cosets (n = n' mod M Z^3); a slot claimed twice, or an atom
!>        matching two orbits, aborts - the partition is then atom-order dependent.
!>        Coordinates are wrapped with the convention of the k-point neighbor lists:
!>        pbc_stable when the SCF ran with k-point symmetry, plain pbc otherwise.
!> \param uf unfolding environment
!> \param cell supercell
!> \param particle_set supercell atoms
!> \param qs_kind_set kinds
!> \param group_of_atom (natom) 0 = fully screened, >= 1 = clustering group
!> \param stable_images wrap atoms with pbc_stable (SCF k-point symmetry anchoring)
!> \param unit_nr summary output unit (<= 0 suppresses)
! **************************************************************************************************
   SUBROUTINE unfolding_env_build_map(uf, cell, particle_set, qs_kind_set, group_of_atom, &
                                      stable_images, unit_nr)
      TYPE(unfolding_env_type), INTENT(INOUT)            :: uf
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: group_of_atom
      LOGICAL, INTENT(IN)                                :: stable_images
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_env_build_map'

      CHARACTER(LEN=default_string_length)               :: name_i, name_ref
      CHARACTER(LEN=msg_length)                          :: msgstr
      INTEGER :: handle, iao, iatom, idet, ikind, iset, ishell, j, jatom, kind_i, kind_r, &
         local_ao, lpos, lstore_max, nao, natom, nkind, nmemb, nref, nset, refj, z_i
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: channel_size, first_sgf, kind_lcount, &
                                                            kind_nsgf, kind_z, last_sgf, ref_atom, &
                                                            ref_base, slot_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: kind_lseq, slot_atom
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: slot_n
      INTEGER, DIMENSION(3)                              :: coset_test, dn, n_cand, n_int
      INTEGER, DIMENSION(3, 3)                           :: adj
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      REAL(KIND=dp)                                      :: max_res, res, res_cand
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: s_pbc
      REAL(KIND=dp), DIMENSION(3)                        :: dcart, ds, r_pbc
      REAL(KIND=dp), DIMENSION(3, 3)                     :: m_dp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      ALLOCATE (first_sgf(natom), last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)
      nao = last_sgf(natom)
      ALLOCATE (uf%channel_of_ao(nao), uf%n_of_ao(3, nao))
      uf%channel_of_ao = 0
      ALLOCATE (uf%channel_group(nao))
      uf%channel_group = 0

      ! orbit candidates are matched by element (GHOST kinds carry their base
      ! element); within an orbit the basis must agree, so cache per kind the
      ! atomic number, nsgf, and the flattened set->shell l sequence (m = -l..l
      ! within a shell, hence the sequence pins every spherical function)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (kind_z(nkind), kind_nsgf(nkind), kind_lcount(nkind))
      lstore_max = 0
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zatom=kind_z(ikind))
         NULLIFY (basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
         IF (ASSOCIATED(basis_set)) THEN
            NULLIFY (nshell, lshell)
            CALL get_gto_basis_set(basis_set, nset=nset, nshell=nshell, l=lshell, nsgf=kind_nsgf(ikind))
            kind_lcount(ikind) = SUM(nshell(1:nset))
         ELSE
            kind_nsgf(ikind) = 0
            kind_lcount(ikind) = 0
         END IF
         lstore_max = MAX(lstore_max, kind_lcount(ikind))
      END DO
      ALLOCATE (kind_lseq(lstore_max, nkind))
      kind_lseq = 0
      DO ikind = 1, nkind
         IF (kind_lcount(ikind) == 0) CYCLE
         NULLIFY (basis_set, nshell, lshell)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
         CALL get_gto_basis_set(basis_set, nset=nset, nshell=nshell, l=lshell)
         lpos = 0
         DO iset = 1, nset
            DO ishell = 1, nshell(iset)
               lpos = lpos + 1
               kind_lseq(lpos, ikind) = lshell(ishell, iset)
            END DO
         END DO
         CPASSERT(lpos == kind_lcount(ikind))
      END DO

      ALLOCATE (s_pbc(3, natom))
      DO iatom = 1, natom
         IF (stable_images) THEN
            r_pbc(1:3) = pbc_stable(particle_set(iatom)%r(1:3), cell)
         ELSE
            r_pbc(1:3) = pbc(particle_set(iatom)%r(1:3), cell)
         END IF
         CALL real_to_scaled(s_pbc(1:3, iatom), r_pbc(1:3), cell)
      END DO

      ALLOCATE (ref_atom(natom), ref_base(natom))
      ALLOCATE (slot_n(3, uf%m, natom), slot_atom(uf%m, natom), slot_count(natom))
      slot_count = 0
      max_res = 0.0_dp
      ! integer adjugate and determinant: two translations n, n' belong to the same
      ! coset of the primitive lattice iff adj . (n - n') is divisible by det(M)
      m_dp = REAL(uf%transmat, KIND=dp)
      adj(1, 1:3) = NINT(vector_product(m_dp(1:3, 2), m_dp(1:3, 3)))
      adj(2, 1:3) = NINT(vector_product(m_dp(1:3, 3), m_dp(1:3, 1)))
      adj(3, 1:3) = NINT(vector_product(m_dp(1:3, 1), m_dp(1:3, 2)))
      idet = NINT(det_3x3(m_dp))
      nref = 0
      uf%n_sets = 0
      uf%n_channels = 0
      uf%n_incomplete = 0
      uf%n_groups = MAXVAL(group_of_atom(1:natom))

      DO iatom = 1, natom
         IF (group_of_atom(iatom) == 0) CYCLE
         kind_i = particle_set(iatom)%atomic_kind%kind_number
         z_i = kind_z(kind_i)

         ! the residual norm is frame-free (h_PC may sit in the input frame).
         ! The plane-distance validation rules out two sites within tol. A second
         ! matching orbit or a taken slot gives an atom-order-dependent partition
         ! and aborts
         refj = 0
         DO j = 1, nref
            jatom = ref_atom(j)
            IF (group_of_atom(jatom) /= group_of_atom(iatom)) CYCLE
            IF (kind_z(particle_set(jatom)%atomic_kind%kind_number) /= z_i) CYCLE
            ds(1:3) = MATMUL(REAL(uf%transmat(1:3, 1:3), KIND=dp), &
                             s_pbc(1:3, iatom) - s_pbc(1:3, jatom))
            n_int(1:3) = NINT(ds(1:3))
            dcart(1:3) = MATMUL(uf%h_pc(1:3, 1:3), ds(1:3) - REAL(n_int(1:3), KIND=dp))
            res = NORM2(dcart(1:3))
            IF (res >= uf%symm_tol) CYCLE
            DO nmemb = 1, slot_count(j)
               dn(1:3) = n_int(1:3) - slot_n(1:3, nmemb, j)
               coset_test(1:3) = MATMUL(adj(1:3, 1:3), dn(1:3))
               IF (ALL(MOD(coset_test(1:3), idet) == 0)) THEN
                  WRITE (UNIT=msgstr, FMT="(A,I0,A,I0,A)") "UNFOLD_BAND| assignment conflict: atoms ", &
                     iatom, " and ", slot_atom(nmemb, j), &
                     " map to the same equivalent-site slot. Set a smaller SYMM_TOL."
                  CPABORT(TRIM(msgstr))
               END IF
            END DO
            IF (refj /= 0) THEN
               WRITE (UNIT=msgstr, FMT="(A,I0,A)") "UNFOLD_BAND| ambiguous site for atom ", iatom, &
                  ": it matches several equivalent-site sets. Set a smaller SYMM_TOL or use ATOMS groups."
               CPABORT(TRIM(msgstr))
            END IF
            refj = j
            n_cand(1:3) = n_int(1:3)
            res_cand = res
         END DO

         IF (refj /= 0) THEN
            ! basis-set consistency inside the orbit: same element with a
            ! different basis must be excluded from this orbit (put such atoms
            ! in a separate ATOMS group to unfold them independently)
            kind_r = particle_set(ref_atom(refj))%atomic_kind%kind_number
            IF (kind_nsgf(kind_i) /= kind_nsgf(kind_r) .OR. &
                kind_lcount(kind_i) /= kind_lcount(kind_r) .OR. &
                ANY(kind_lseq(1:kind_lcount(kind_i), kind_i) /= kind_lseq(1:kind_lcount(kind_r), kind_r))) THEN
               CALL get_qs_kind(qs_kind_set(kind_i), name=name_i)
               CALL get_qs_kind(qs_kind_set(kind_r), name=name_ref)
               WRITE (UNIT=msgstr, FMT="(A,I0,A,I0,A,A,A,A,A)") &
                  "UNFOLD_BAND| atoms ", iatom, " and ", ref_atom(refj), &
                  " (kinds ", TRIM(name_i), " / ", TRIM(name_ref), &
                  ") differ in basis set. Put them in separate ATOMS groups."
               CPABORT(TRIM(msgstr))
            END IF
         END IF

         IF (refj == 0) THEN
            nref = nref + 1
            refj = nref
            ref_atom(refj) = iatom
            ref_base(refj) = uf%n_channels + 1
            uf%n_channels = uf%n_channels + (last_sgf(iatom) - first_sgf(iatom) + 1)
            uf%channel_group(ref_base(refj):uf%n_channels) = group_of_atom(iatom)
            n_cand(1:3) = 0
            res_cand = 0.0_dp
         END IF

         slot_count(refj) = slot_count(refj) + 1
         slot_n(1:3, slot_count(refj), refj) = n_cand(1:3)
         slot_atom(slot_count(refj), refj) = iatom
         max_res = MAX(max_res, res_cand)

         DO local_ao = first_sgf(iatom), last_sgf(iatom)
            uf%channel_of_ao(local_ao) = ref_base(refj) + (local_ao - first_sgf(iatom))
            uf%n_of_ao(1:3, local_ao) = n_cand(1:3)
         END DO
      END DO

      ! each coset hosts at most one atom, so channel_size equals the member count
      ALLOCATE (channel_size(uf%n_channels))
      channel_size = 0
      DO iao = 1, nao
         IF (uf%channel_of_ao(iao) > 0) THEN
            channel_size(uf%channel_of_ao(iao)) = channel_size(uf%channel_of_ao(iao)) + 1
         END IF
      END DO
      IF (ANY(channel_size > uf%m)) THEN
         CPABORT("UNFOLD_BAND| equivalent-site channel overflow: check TRANSFORMATION_MATRIX/SYMM_TOL")
      END IF
      uf%n_sets = nref
      uf%n_incomplete = COUNT(slot_count(1:nref) < uf%m)

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T2,A,I0,A,I0,A,I0)") "UNFOLD_BAND| equivalent-site sets: ", &
            uf%n_sets, " (", uf%n_channels, " AO channels); incomplete sets (vacancies): ", &
            uf%n_incomplete
         WRITE (UNIT=unit_nr, FMT="(T2,A,F12.6)") "UNFOLD_BAND| max assignment residual [angstrom] = ", &
            max_res*angstrom
      END IF

      DEALLOCATE (channel_size, first_sgf, last_sgf, kind_lcount, kind_lseq, kind_nsgf, kind_z, ref_atom, &
                  ref_base, slot_n, slot_atom, slot_count, s_pbc)
      CALL timestop(handle)
   END SUBROUTINE unfolding_env_build_map

! **************************************************************************************************
!> \brief Fold primitive-cell k-points onto deduplicated supercell K-points
!> \param uf unfolding environment
!> \param kppath (3, npoints) primitive fractional k-points
!> \param kpfolded (3, nfold) unique supercell K in [-0.5, 0.5)
!> \param path2K (npoints) folded-K index of each path point
!> \param nfold number of unique supercell K-points
! **************************************************************************************************
   SUBROUTINE unfolding_fold_kpath(uf, kppath, kpfolded, path2K, nfold)
      TYPE(unfolding_env_type), INTENT(IN)               :: uf
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kppath
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kpfolded
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: path2K
      INTEGER, INTENT(OUT)                               :: nfold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_fold_kpath'

      INTEGER                                            :: handle, ik, jk, npoints
      REAL(KIND=dp), DIMENSION(3)                        :: dvec, ffold

      CALL timeset(routineN, handle)

      npoints = SIZE(kppath, 2)
      ALLOCATE (kpfolded(3, npoints), path2K(npoints))
      path2K = 0
      nfold = 0
      DO ik = 1, npoints
         ! F = M^T f (paper row convention f.M); wrap into [-0.5, 0.5)
         ffold(1:3) = MATMUL(TRANSPOSE(REAL(uf%transmat(1:3, 1:3), KIND=dp)), kppath(1:3, ik))
         ffold(1:3) = ffold(1:3) - ANINT(ffold(1:3))
         DO jk = 1, nfold
            dvec(1:3) = ffold(1:3) - kpfolded(1:3, jk)
            dvec(1:3) = dvec(1:3) - ANINT(dvec(1:3))
            IF (MAXVAL(ABS(dvec(1:3))) < kdifftol) THEN
               path2K(ik) = jk
               EXIT
            END IF
         END DO
         IF (path2K(ik) == 0) THEN
            nfold = nfold + 1
            kpfolded(1:3, nfold) = ffold(1:3)
            path2K(ik) = nfold
         END IF
      END DO
      CALL timestop(handle)
   END SUBROUTINE unfolding_fold_kpath

! **************************************************************************************************
!> \brief Accumulate unfolding weights, eigenvalues, and occupations for one batch of unique
!>        supercell K-points into the full-set arrays. Releases per-K MO coefficients and
!>        S(K)^{1/2} progressively as each K-point is processed.
!> \param uf unfolding environment
!> \param qs_env ...
!> \param kpoint batch of unique supercell K-points with their MOs
!> \param kpgeneral (3, npoints) primitive-cell path points
!> \param path2K (npoints) global supercell K-point index of each path point
!> \param k_offset global index offset preceding this batch (first batch has k_offset = 0)
!> \param wall (npoints, nmo, 0:n_groups, nspins) accumulated unfolding weights (slice 0 = total)
!> \param eigval_all (nmo, nfold, nspins) accumulated band eigenvalues
!> \param occnum_all (nmo, nfold, nspins) accumulated band occupations
!> \param sumdev running maximum deviation from the subspace sum rule across batches
! **************************************************************************************************
   SUBROUTINE unfolding_process_batch(uf, qs_env, kpoint, kpgeneral, path2K, k_offset, &
                                      wall, eigval_all, occnum_all, sumdev)
      TYPE(unfolding_env_type), INTENT(IN)               :: uf
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoint
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kpgeneral
      INTEGER, DIMENSION(:), INTENT(IN)                  :: path2K
      INTEGER, INTENT(IN)                                :: k_offset
      REAL(KIND=dp), DIMENSION(:, :, 0:, :), &
         INTENT(INOUT)                                   :: wall
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: eigval_all, occnum_all
      REAL(KIND=dp), INTENT(INOUT)                       :: sumdev

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_process_batch'

      INTEGER                                            :: handle, ik, ik_glob, ikpgr, iq, iqc, &
                                                            ispin, jpath, nfold_b, nmo, npoints, &
                                                            nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: is_kp_source
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: psub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fqs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wmat
      REAL(KIND=dp), DIMENSION(3)                        :: dvec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s, scf_env)
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s, scf_env=scf_env)
      npoints = SIZE(path2K)
      nmo = SIZE(wall, 2)
      nspins = SIZE(wall, 4)

      ! S(k)^{1/2} of the batch's local K-points, reused across spins and path points
      IF (.NOT. uf%lowdin_proj) CALL diag_kp_smat(matrix_s, kpoint, scf_env%scf_work1)
      CALL get_kpoint_info(kpoint, nkp=nfold_b, kp_range=kp_range, xkp=xkp)
      is_kp_source = kpoint%para_env_kp%is_source()

      ALLOCATE (wmat(uf%m, nmo, 0:uf%n_groups), psub(nmo))
      DO ik = 1, nfold_b
         IF (ik < kp_range(1) .OR. ik > kp_range(2)) CYCLE
         ik_glob = k_offset + ik
         ikpgr = ik - kp_range(1) + 1
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         CALL unfolding_unfolded_kspace(uf, xkp(1:3, ik), fqs)
         DO ispin = 1, nspins
            CALL unfolding_kp_weights(uf, kp, ispin, fqs, wmat, psub)
            ! generalized subspace sum rule: holds for screened selections too
            sumdev = MAX(sumdev, MAXVAL(ABS(SUM(wmat(:, 1:nmo, 0), 1) - psub(1:nmo))))
            ! with all atoms selected the Loewdin populations are the MO
            ! S-orthonormality row sums and must equal 1
            IF (uf%full_selection .AND. MAXVAL(ABS(psub(1:nmo) - 1.0_dp)) > mo_orthonorm_tol) THEN
               CPABORT("UNFOLD_BAND| the full-selection Loewdin population deviates from 1")
            END IF
            ! Only the source rank of the owning k-point group records each K-point's
            ! entries so a single world para_env%sum works even if group sizes vary across batches.
            IF (is_kp_source) THEN
               CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues, occupation_numbers=occupation_numbers)
               eigval_all(1:nmo, ik_glob, ispin) = eigenvalues(1:nmo)
               occnum_all(1:nmo, ik_glob, ispin) = occupation_numbers(1:nmo)
            END IF
            DO jpath = 1, npoints
               IF (path2K(jpath) /= ik_glob) CYCLE
               iq = 0
               DO iqc = 1, SIZE(fqs, 2)
                  dvec(1:3) = fqs(1:3, iqc) - kpgeneral(1:3, jpath)
                  dvec(1:3) = dvec(1:3) - ANINT(dvec(1:3))
                  IF (MAXVAL(ABS(dvec(1:3))) < kdifftol) THEN
                     iq = iqc
                     EXIT
                  END IF
               END DO
               IF (iq == 0) THEN
                  CPABORT("UNFOLD_BAND| path point not in unfolded set")
               END IF
               IF (is_kp_source) THEN
                  wall(jpath, 1:nmo, 0:uf%n_groups, ispin) = wmat(iq, 1:nmo, 0:uf%n_groups)
               END IF
            END DO
            IF (ASSOCIATED(kp%mos(ispin)%cmo_coeff)) THEN
               CALL cp_cfm_release(kp%mos(ispin)%cmo_coeff)
               DEALLOCATE (kp%mos(ispin)%cmo_coeff)
               NULLIFY (kp%mos(ispin)%cmo_coeff)
            END IF
         END DO
         CALL cp_cfm_release(kp%cshalf)
         DEALLOCATE (fqs)
      END DO
      DEALLOCATE (wmat, psub)

      CALL timestop(handle)
   END SUBROUTINE unfolding_process_batch

! **************************************************************************************************
!> \brief Reduce accumulated unfolding results across MPI ranks and write the unfolded band
!>        structure set and sum-rule check
!> \param para_env global MPI communicator
!> \param kpfolded (3, nfold) unique supercell K-points
!> \param kpgeneral (3, npoints) primitive-cell path points
!> \param path2K (npoints) supercell K-point index of each path point
!> \param eigval_all (nmo, nfold, nspins) band eigenvalues (reduced in place)
!> \param occnum_all (nmo, nfold, nspins) band occupations (reduced in place)
!> \param wall (npoints, nmo, 0:n_groups, nspins) unfolding weights (reduced in place)
!> \param sumdev maximum subspace sum-rule deviation (reduced in place)
!> \param unit_nr default output unit
!> \param bs_data_unit band structure data unit
! **************************************************************************************************
   SUBROUTINE unfolding_finish_set(para_env, kpfolded, kpgeneral, path2K, eigval_all, occnum_all, &
                                   wall, sumdev, unit_nr, bs_data_unit)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kpfolded, kpgeneral
      INTEGER, DIMENSION(:), INTENT(IN)                  :: path2K
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: eigval_all, occnum_all
      REAL(KIND=dp), DIMENSION(:, :, 0:, :), &
         INTENT(INOUT)                                   :: wall
      REAL(KIND=dp), INTENT(INOUT)                       :: sumdev
      INTEGER, INTENT(IN)                                :: unit_nr, bs_data_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_finish_set'

      INTEGER                                            :: g, handle, ik, ikmap, imo, ispin, &
                                                            n_groups, ngrp_out, nmo, npoints, &
                                                            nspins
      REAL(KIND=dp)                                      :: path_weight

      CALL timeset(routineN, handle)

      CALL para_env%sum(wall)
      CALL para_env%sum(eigval_all)
      CALL para_env%sum(occnum_all)
      CALL para_env%max(sumdev)

      IF (bs_data_unit > 0) THEN
         npoints = SIZE(wall, 1)
         nmo = SIZE(wall, 2)
         n_groups = UBOUND(wall, 3)  ! group axis of wall runs 0:n_groups (0 = total)
         nspins = SIZE(wall, 4)
         path_weight = 1.0_dp/REAL(npoints, KIND=dp)
         ! per-group weight columns only for two or more groups: a single group equals the total
         ngrp_out = MERGE(n_groups, 0, n_groups >= 2)

         DO ik = 1, npoints
            ikmap = path2K(ik)
            DO ispin = 1, nspins
               WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0,A,3(1X,F14.8))") &
                  "# UNFOLD_BAND| path point ", ik, " -> SC K #", ikmap, " :", kpfolded(1:3, ikmap)
               WRITE (UNIT=bs_data_unit, FMT="(A,I0,T15,A,I0,A,T24,3(1X,F14.8),3X,F14.8)") &
                  "#  Point ", ik, "  Spin ", ispin, ":", kpgeneral(1:3, ik), path_weight
               WRITE (UNIT=bs_data_unit, FMT="(A)", ADVANCE="NO") &
                  "#   Band    Energy [eV]     Occupation   Unfold_Weight"
               DO g = 1, ngrp_out
                  WRITE (UNIT=bs_data_unit, FMT="(A,I0)", ADVANCE="NO") "   W_grp", g
               END DO
               WRITE (UNIT=bs_data_unit, FMT="(A)") ""
               DO imo = 1, nmo
                  WRITE (UNIT=bs_data_unit, FMT="(T2,I7,*(1X,F14.8))") imo, &
                     eigval_all(imo, ikmap, ispin)*evolt, occnum_all(imo, ikmap, ispin), &
                     (wall(ik, imo, g, ispin), g=0, ngrp_out)
               END DO
            END DO
         END DO
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, FMT="(T2,A,ES9.2)") &
            "UNFOLD_BAND| Set check max|sum_q W - Loewdin population| = ", sumdev
      END IF

      CALL timestop(handle)
   END SUBROUTINE unfolding_finish_set

! **************************************************************************************************
!> \brief The m primitive k-points f (in [0,1)) with M^T f = K (mod 1)
!> \param uf unfolding environment
!> \param K supercell fractional k-point
!> \param fqs (3, m) unfolded primitive k-points
! **************************************************************************************************
   SUBROUTINE unfolding_unfolded_kspace(uf, K, fqs)
      TYPE(unfolding_env_type), INTENT(IN)               :: uf
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: K
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: fqs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_unfolded_kspace'

      INTEGER                                            :: handle, nfound, px, py, pz, r
      INTEGER, DIMENSION(3)                              :: phi, plo
      REAL(KIND=dp), DIMENSION(3)                        :: f, pshift
      REAL(KIND=dp), DIMENSION(3, 3)                     :: minvt

      CALL timeset(routineN, handle)

      ! p = M^T f - K with f in [0,1): component r of M^T f lies between the sums of
      ! the negative and positive entries of column r of M; shifting by the K range
      ! [-0.5, 0.5) keeps these integer bounds exact (no slack needed)
      minvt = TRANSPOSE(inv_3x3(REAL(uf%transmat, KIND=dp)))
      DO r = 1, 3
         plo(r) = SUM(MIN(uf%transmat(1:3, r), 0))
         phi(r) = SUM(MAX(uf%transmat(1:3, r), 0))
      END DO

      ALLOCATE (fqs(3, uf%m))
      nfound = 0
      DO px = plo(1), phi(1)
         DO py = plo(2), phi(2)
            DO pz = plo(3), phi(3)
               pshift(1:3) = REAL([px, py, pz], KIND=dp)
               f(1:3) = MATMUL(minvt(1:3, 1:3), K(1:3) + pshift(1:3))
               IF (ALL(f(1:3) >= -unit_tol) .AND. ALL(f(1:3) < 1.0_dp - unit_tol)) THEN
                  nfound = nfound + 1
                  IF (nfound > uf%m) CPABORT("UNFOLD_BAND| unfolded k set larger than m")
                  fqs(1:3, nfound) = MAX(f(1:3), 0.0_dp)
               END IF
            END DO
         END DO
      END DO
      IF (nfound /= uf%m) THEN
         CPABORT("UNFOLD_BAND| unfolded k set incomplete. Check TRANSFORMATION_MATRIX.")
      END IF
      CALL timestop(handle)
   END SUBROUTINE unfolding_unfolded_kspace

! **************************************************************************************************
!> \brief Unfolding weights for one local k-point environment (complex wavefunctions).
!>        Requires kp%cshalf, i.e. diag_kp_smat must have run for this kpoint object.
!>        Per atom group the weights obey the generalized subspace sum rule
!>        sum_q W^k_q_KN(group) = sum_{mu in group} |C'_muN|^2 (Loewdin population);
!>        slice 0 carries the total over all selected atoms (= 1 for a full selection).
!> \param uf unfolding environment with orbit map
!> \param kp k-point environment (one local K)
!> \param ispin spin index
!> \param fqs (3, nq) unfolded primitive k-points of this K
!> \param weights (nq, nmo, 0:n_groups) spectral weights; slice 0 = total
!> \param psubspace (nmo) Loewdin population of the selected atoms per band
! **************************************************************************************************
   SUBROUTINE unfolding_kp_weights(uf, kp, ispin, fqs, weights, psubspace)
      TYPE(unfolding_env_type), INTENT(IN)               :: uf
      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: fqs
      REAL(KIND=dp), DIMENSION(:, :, 0:), INTENT(OUT)    :: weights
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: psubspace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'unfolding_kp_weights'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: phase_loc
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: tsum
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      INTEGER                                            :: g, handle, iao, ich, icol, imo, iq, &
                                                            irow, j1, j2, jloc, mo_blk, nao, nb, &
                                                            ncol_local, nmo, nq, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: arg
      TYPE(cp_cfm_type)                                  :: cwork
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(mp_para_env_type), POINTER                    :: para_env_kp

      CALL timeset(routineN, handle)

      NULLIFY (mo_set, local_data, row_indices, col_indices, para_env_kp)
      mo_set => kp%mos(ispin)
      CALL get_mo_set(mo_set, nao=nao, nmo=nmo)
      IF (.NOT. ASSOCIATED(mo_set%cmo_coeff)) THEN
         CPABORT("UNFOLD_BAND| complex MO coefficients expected (band structure k-points)")
      END IF
      IF (.NOT. ASSOCIATED(kp%cshalf%matrix_struct)) THEN
         CPABORT("UNFOLD_BAND| kp%cshalf is missing. diag_kp_smat must run before unfolding.")
      END IF

      nq = SIZE(fqs, 2)
      CALL cp_cfm_create(cwork, mo_set%cmo_coeff%matrix_struct)
      CALL lowdin_kp_mo_coeff(kp, ispin, .FALSE., cshalfc=cwork)
      CALL cp_cfm_get_info(cwork, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices, &
                           local_data=local_data, para_env=para_env_kp)

      ! Loewdin population of the selected atoms: the generalized sum-rule reference
      psubspace(1:nmo) = 0.0_dp
      DO icol = 1, ncol_local
         imo = col_indices(icol)
         DO irow = 1, nrow_local
            iao = row_indices(irow)
            IF (uf%channel_of_ao(iao) > 0) THEN
               psubspace(imo) = psubspace(imo) + ABS(local_data(irow, icol))**2
            END IF
         END DO
      END DO
      IF (para_env_kp%num_pe > 1) CALL para_env_kp%sum(psubspace(1:nmo))

      mo_blk = MAX(1, MIN(nmo, unfolding_mo_block))
      ALLOCATE (phase_loc(nrow_local), tsum(uf%n_channels, mo_blk))
      weights = 0.0_dp
      DO iq = 1, nq
         DO irow = 1, nrow_local
            iao = row_indices(irow)
            IF (uf%channel_of_ao(iao) == 0) CYCLE
            arg = -twopi*DOT_PRODUCT(fqs(1:3, iq), REAL(uf%n_of_ao(1:3, iao), KIND=dp))
            phase_loc(irow) = EXP(CMPLX(0.0_dp, arg, KIND=dp))
         END DO
         DO j1 = 1, nmo, mo_blk
            j2 = MIN(nmo, j1 + mo_blk - 1)
            nb = j2 - j1 + 1
            tsum(1:uf%n_channels, 1:nb) = (0.0_dp, 0.0_dp)
            DO icol = 1, ncol_local
               imo = col_indices(icol)
               IF (imo < j1 .OR. imo > j2) CYCLE
               jloc = imo - j1 + 1
               DO irow = 1, nrow_local
                  iao = row_indices(irow)
                  ich = uf%channel_of_ao(iao)
                  IF (ich == 0) CYCLE
                  tsum(ich, jloc) = tsum(ich, jloc) + phase_loc(irow)*local_data(irow, icol)
               END DO
            END DO
            IF (para_env_kp%num_pe > 1) CALL para_env_kp%sum(tsum(1:uf%n_channels, 1:nb))
            DO jloc = 1, nb
               imo = j1 + jloc - 1
               DO ich = 1, uf%n_channels
                  g = uf%channel_group(ich)
                  weights(iq, imo, g) = weights(iq, imo, g) + ABS(tsum(ich, jloc))**2
               END DO
            END DO
         END DO
      END DO
      CALL cp_cfm_release(cwork)
      weights(1:nq, 1:nmo, 0:uf%n_groups) = weights(1:nq, 1:nmo, 0:uf%n_groups)/REAL(uf%m, KIND=dp)
      DO g = 1, uf%n_groups
         weights(1:nq, 1:nmo, 0) = weights(1:nq, 1:nmo, 0) + weights(1:nq, 1:nmo, g)
      END DO

      DEALLOCATE (phase_loc, tsum)
      CALL timestop(handle)
   END SUBROUTINE unfolding_kp_weights

! **************************************************************************************************
!> \brief Plane distance of lattice-vector columns j and k of a cell (columns)
!> \param h ...
!> \param j ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION plane_dist_of(h, j, k) RESULT(d)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h
      INTEGER, INTENT(IN)                                :: j, k
      REAL(KIND=dp)                                      :: d

      d = ABS(det_3x3(h))/NORM2(vector_product(h(1:3, j), h(1:3, k)))
   END FUNCTION plane_dist_of

END MODULE kpoint_unfolding
