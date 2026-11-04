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
   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              pbc_stable,&
                                              real_to_scaled
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: lowdin_kp_mo_coeff
   USE kpoint_types,                    ONLY: kpoint_env_type
   USE mathconstants,                   ONLY: twopi
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_unfolding'
   INTEGER, PARAMETER, PRIVATE :: msg_length = 256
   ! tolerance for identifying k-points that differ by a reciprocal lattice vector
   REAL(KIND=dp), PARAMETER :: kdifftol = 1.0E-6_dp
   REAL(KIND=dp), PARAMETER :: unit_tol = 1.0E-9_dp
   ! full-atom-selection Loewdin populations equal 1 to the k-point MO orthonormality
   ! accuracy (set by the diagonalization threshold, EPS_DEFAULT-level); violations
   ! indicate an inconsistent S(k) or MO set, well above this guard
   REAL(KIND=dp), PARAMETER :: mo_orthonorm_tol = 1.0E-8_dp

   !> Transformation matrix, derived primitive cell, and the per-AO orbit map
   TYPE :: unfolding_env_type
      INTEGER, DIMENSION(3, 3) :: transmat = 0             ! hmat_SC = hmat_PC . M (column convention)
      INTEGER :: m = 0                                      ! |det(M)|
      REAL(KIND=dp) :: symm_tol = 0.5_dp/angstrom         ! default 0.5 Angstrom in internal Bohr
      REAL(KIND=dp), DIMENSION(3, 3) :: h_pc = 0.0_dp       ! derived primitive lattice (columns)
      INTEGER :: n_channels = 0                             ! per-AO equivalence channels
      INTEGER, DIMENSION(:), ALLOCATABLE :: channel_of_ao   ! (nao)
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: n_of_ao      ! (3, nao) primitive translation labels
      INTEGER :: n_incomplete = 0                           ! channels with fewer than m members
      INTEGER :: n_groups = 1                               ! atom groups (independent clustering)
      INTEGER, DIMENSION(:), ALLOCATABLE :: channel_group   ! group of each channel (upper bound nao)
   END TYPE unfolding_env_type

   PUBLIC :: unfolding_env_type, unfolding_env_create, unfolding_env_release, &
             unfolding_env_build_map, unfolding_fold_kpath, unfolding_unfolded_kspace, &
             unfolding_kp_weights, kdifftol, mo_orthonorm_tol

CONTAINS

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
      REAL(KIND=dp), DIMENSION(3, 3)                     :: minv

      CALL timeset(routineN, handle)

      uf%transmat = transmat
      uf%symm_tol = symm_tol
      uf%m = ABS(int_det3(transmat))
      IF (uf%m == 0) CPABORT("UNFOLD_BAND| TRANSFORMATION_MATRIX is singular")

      minv = inverse_int3(transmat)
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
      uf%n_channels = 0
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
!> \param iounit summary output unit (<= 0 suppresses)
! **************************************************************************************************
   SUBROUTINE unfolding_env_build_map(uf, cell, particle_set, qs_kind_set, group_of_atom, &
                                      stable_images, iounit)
      TYPE(unfolding_env_type), INTENT(INOUT)            :: uf
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: group_of_atom
      LOGICAL, INTENT(IN)                                :: stable_images
      INTEGER, INTENT(IN)                                :: iounit

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
      adj = adj_int3(uf%transmat)
      idet = int_det3(uf%transmat)
      nref = 0
      uf%n_channels = 0
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
            res = SQRT(DOT_PRODUCT(dcart(1:3), dcart(1:3)))
            IF (res >= uf%symm_tol) CYCLE
            DO nmemb = 1, slot_count(j)
               dn(1:3) = n_int(1:3) - slot_n(1:3, nmemb, j)
               coset_test(1:3) = MATMUL(adj(1:3, 1:3), dn(1:3))
               IF (ALL(MOD(coset_test(1:3), idet) == 0)) THEN
                  WRITE (UNIT=msgstr, FMT="(A,I0,A,I0,A)") "UNFOLD_BAND| assignment conflict: atoms ", &
                     iatom, " and ", slot_atom(nmemb, j), " map to the same orbit slot. Set a smaller SYMM_TOL."
                  CPABORT(TRIM(msgstr))
               END IF
            END DO
            IF (refj /= 0) THEN
               WRITE (UNIT=msgstr, FMT="(A,I0,A)") "UNFOLD_BAND| ambiguous site for atom ", iatom, &
                  ": it matches several orbits. Set a smaller SYMM_TOL or use ATOMS groups."
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
         CPABORT("UNFOLD_BAND| orbit channel overflow: check TRANSFORMATION_MATRIX/SYMM_TOL")
      END IF
      uf%n_incomplete = COUNT(channel_size < uf%m)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,I0,A,I0)") "UNFOLD_BAND| orbit channels: ", uf%n_channels, &
            "; incomplete channels (vacancies): ", uf%n_incomplete
         WRITE (UNIT=iounit, FMT="(T2,A,F12.6)") "UNFOLD_BAND| max assignment residual [angstrom] = ", &
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
      minvt = TRANSPOSE(inverse_int3(uf%transmat))
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

      COMPLEX(KIND=dp)                                   :: phase
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cbuf, tsum
      INTEGER                                            :: g, handle, iao, ich, iq, nao, nmo, nq
      REAL(KIND=dp)                                      :: arg
      TYPE(cp_cfm_type)                                  :: cwork
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      NULLIFY (mo_set)
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
      ALLOCATE (cbuf(nao, nmo))
      CALL cp_cfm_get_submatrix(cwork, cbuf)
      CALL cp_cfm_release(cwork)

      ALLOCATE (tsum(uf%n_channels, nmo))
      weights = 0.0_dp
      DO iq = 1, nq
         tsum(1:uf%n_channels, 1:nmo) = (0.0_dp, 0.0_dp)
         DO iao = 1, nao
            IF (uf%channel_of_ao(iao) == 0) CYCLE
            arg = twopi*DOT_PRODUCT(fqs(1:3, iq), REAL(uf%n_of_ao(1:3, iao), KIND=dp))
            phase = CMPLX(COS(arg), -SIN(arg), KIND=dp)
            ich = uf%channel_of_ao(iao)
            tsum(ich, 1:nmo) = tsum(ich, 1:nmo) + phase*cbuf(iao, 1:nmo)
         END DO
         DO ich = 1, uf%n_channels
            g = uf%channel_group(ich)
            weights(iq, 1:nmo, g) = weights(iq, 1:nmo, g) + ABS(tsum(ich, 1:nmo))**2
         END DO
      END DO
      weights(1:nq, 1:nmo, 0:uf%n_groups) = weights(1:nq, 1:nmo, 0:uf%n_groups)/REAL(uf%m, KIND=dp)
      DO g = 1, uf%n_groups
         weights(1:nq, 1:nmo, 0) = weights(1:nq, 1:nmo, 0) + weights(1:nq, 1:nmo, g)
      END DO
      ! Loewdin population of the selected atoms: the generalized sum-rule reference
      psubspace(1:nmo) = 0.0_dp
      DO iao = 1, nao
         IF (uf%channel_of_ao(iao) > 0) THEN
            psubspace(1:nmo) = psubspace(1:nmo) + ABS(cbuf(iao, 1:nmo))**2
         END IF
      END DO

      DEALLOCATE (tsum, cbuf)
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

      REAL(KIND=dp)                                      :: det
      REAL(KIND=dp), DIMENSION(3)                        :: cross_jk

      cross_jk(1) = h(2, j)*h(3, k) - h(3, j)*h(2, k)
      cross_jk(2) = h(3, j)*h(1, k) - h(1, j)*h(3, k)
      cross_jk(3) = h(1, j)*h(2, k) - h(2, j)*h(1, k)
      det = h(1, 1)*(h(2, 2)*h(3, 3) - h(3, 2)*h(2, 3)) - h(2, 1)*(h(1, 2)*h(3, 3) - h(3, 2)*h(1, 3)) &
            + h(3, 1)*(h(1, 2)*h(2, 3) - h(2, 2)*h(1, 3))
      d = ABS(det)/SQRT(DOT_PRODUCT(cross_jk, cross_jk))
   END FUNCTION plane_dist_of

! **************************************************************************************************
!> \brief Integer adjugate of a 3x3 matrix (transpose of cofactors)
!> \param a ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION adj_int3(a) RESULT(adj)
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: a
      INTEGER, DIMENSION(3, 3)                           :: adj

      adj(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
      adj(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
      adj(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
      adj(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
      adj(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
      adj(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
      adj(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
      adj(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
      adj(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
   END FUNCTION adj_int3

! **************************************************************************************************
!> \brief Determinant of a 3x3 integer matrix
!> \param a ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION int_det3(a) RESULT(det)
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: a
      INTEGER                                            :: det

      det = a(1, 1)*a(2, 2)*a(3, 3) + a(1, 2)*a(2, 3)*a(3, 1) + a(1, 3)*a(2, 1)*a(3, 2) &
            - a(1, 3)*a(2, 2)*a(3, 1) - a(1, 1)*a(2, 3)*a(3, 2) - a(1, 2)*a(2, 1)*a(3, 3)
   END FUNCTION int_det3

! **************************************************************************************************
!> \brief Real inverse of a 3x3 integer matrix (adjugate / determinant)
!> \param a ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION inverse_int3(a) RESULT(ainv)
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: a
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ainv

      ainv = REAL(adj_int3(a), KIND=dp)/REAL(int_det3(a), KIND=dp)
   END FUNCTION inverse_int3

END MODULE kpoint_unfolding
