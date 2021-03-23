!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to direct methods for electron repulsion integrals for MP2.
! **************************************************************************************************
#:def conditional(n)
$:'' if n else '.NOT.'
#:enddef

MODULE mp2_eri
   USE ai_contraction_sphi, ONLY: ab_contract, &
                                  abc_contract
   USE ai_coulomb, ONLY: coulomb2_new, &
                         coulomb3
   USE atomic_kind_types, ONLY: atomic_kind_type, &
                                get_atomic_kind_set
   USE basis_set_types, ONLY: gto_basis_set_p_type, &
                              gto_basis_set_type
   USE cell_types, ONLY: cell_create, &
                         cell_release, &
                         cell_type, &
                         init_cell, &
                         pbc
   USE cp_eri_mme_interface, ONLY: cp_eri_mme_finalize, &
                                   cp_eri_mme_init_read_input, &
                                   cp_eri_mme_param, &
                                   cp_eri_mme_set_params
   USE cp_para_types, ONLY: cp_para_env_type
   USE dbcsr_api, ONLY: dbcsr_get_block_p, &
                        dbcsr_p_type
   USE eri_mme_integrate, ONLY: eri_mme_2c_integrate, &
                                eri_mme_3c_integrate
   USE eri_mme_test, ONLY: eri_mme_2c_perf_acc_test, &
                           eri_mme_3c_perf_acc_test
   USE eri_mme_types, ONLY: eri_mme_param, &
                            eri_mme_set_potential, eri_mme_coulomb, eri_mme_longrange
   USE input_constants, ONLY: do_eri_gpw, &
                              do_eri_mme, &
                              do_eri_os, &
                              do_potential_coulomb, &
                              do_potential_long
   USE input_section_types, ONLY: section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: dp
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE message_passing, ONLY: mp_sum
   USE orbital_pointers, ONLY: coset, &
                               init_orbital_pointers, &
                               ncoset
   USE particle_types, ONLY: particle_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_integral_utils, ONLY: basis_set_list_setup
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     get_neighbor_list_set_p, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type
   USE util, ONLY: get_limit
   USE cp_eri_mme_interface, ONLY: cp_eri_mme_update_local_counts
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_eri'

   PUBLIC :: &
      mp2_eri_2c_integrate, &
      mp2_eri_3c_integrate, &
      mp2_eri_allocate_forces, &
      mp2_eri_deallocate_forces, &
      mp2_eri_force

   TYPE mp2_eri_force
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: forces
   END TYPE mp2_eri_force

CONTAINS

! **************************************************************************************************
!> \brief high-level integration routine for 2c integrals over CP2K basis sets.
!>        Contiguous column-wise distribution and parallelization over pairs of sets.
!> \param param ...
!> \param para_env mpi environment for local columns
!> \param potential_parameter ...
!> \param qs_env ...
!> \param basis_type_a ...
!> \param basis_type_b ...
!> \param hab columns of ERI matrix
!> \param first_b first column of hab
!> \param last_b last column of hab
!> \param eri_method ...
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \param hdab ...
!> \param hadb ...
!> \param reflection_z_a ...
!> \param reflection_z_b ...
!> \param do_reflection_a ...
!> \param do_reflection_b ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_2c_integrate(param, potential_parameter, para_env, qs_env, basis_type_a, basis_type_b, hab, first_b, &
                                   last_b, eri_method, pab, force_a, force_b, hdab, hadb, &
                                   reflection_z_a, reflection_z_b, do_reflection_a, do_reflection_b)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_a, basis_type_b
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      INTEGER, INTENT(IN)                                :: first_b, last_b
      INTEGER, INTENT(IN), OPTIONAL                      :: eri_method
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pab
      TYPE(mp2_eri_force), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT), OPTIONAL             :: force_a, force_b
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: hdab, hadb
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: reflection_z_a, reflection_z_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_reflection_a, do_reflection_b

      CHARACTER(len=*), PARAMETER :: routineN = 'mp2_eri_2c_integrate'

      INTEGER :: atom_a, atom_b, atom_end, atom_start, first_set, G_count, handle, iatom, ikind, &
                 iset, jatom, jkind, jset, jset_end, jset_start, last_set, my_eri_method, my_setpair, &
                 n_setpair, natom, nkind, nseta, nseta_total, nsetb, nsetb_total, offset_a_end, &
                 offset_a_start, offset_b_end, offset_b_start, R_count, set_end, set_offset_end, &
                 set_offset_start, set_start, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, natom_of_kind
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: eri_offsets
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: map_it_here, my_do_reflection_a, &
                                                            my_do_reflection_b
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sphi_a, sphi_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      my_eri_method = do_eri_mme
      IF (PRESENT(eri_method)) my_eri_method = eri_method

      my_do_reflection_a = .FALSE.
      IF (PRESENT(do_reflection_a) .AND. PRESENT(reflection_z_a)) my_do_reflection_a = do_reflection_a

      my_do_reflection_b = .FALSE.
      IF (PRESENT(do_reflection_b) .AND. PRESENT(reflection_z_b)) my_do_reflection_b = do_reflection_b

      G_count = 0; R_count = 0; 
      ! get mapping between ERIs and atoms, sets, set offsets
      CALL get_eri_offsets(qs_env, basis_type_b, eri_offsets)

      atom_start = eri_offsets(first_b, 1)
      set_start = eri_offsets(first_b, 2)
      set_offset_start = eri_offsets(first_b, 3)

      atom_end = eri_offsets(last_b, 1)
      set_end = eri_offsets(last_b, 2)
      set_offset_end = eri_offsets(last_b, 3)

      ! get QS stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set, natom=natom, nkind=nkind)

      ALLOCATE (atom_of_kind(natom))
      ALLOCATE (kind_of(natom))
      ALLOCATE (natom_of_kind(nkind))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, natom_of_kind=natom_of_kind, atom_of_kind=atom_of_kind)

      IF (PRESENT(force_a)) CALL mp2_eri_allocate_forces(force_a, natom_of_kind)
      IF (PRESENT(force_b)) CALL mp2_eri_allocate_forces(force_b, natom_of_kind)

      ! get total number of local set pairs to integrate
      nseta_total = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)
         nseta_total = nseta_total + basis_set_a%nset
      ENDDO

      nsetb_total = 0
      DO jatom = atom_start, atom_end
         jkind = kind_of(jatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)
         nsetb_total = nsetb_total + basis_set_b%nset
      ENDDO

      n_setpair = nseta_total*nsetb_total

      my_setpair = 0

      offset_a_end = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         atom_a = atom_of_kind(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         npgfa => basis_set_a%npgf

         ra(:) = pbc(particle_set(iatom)%r, cell)

         IF (my_do_reflection_a) THEN
            ra(3) = 2.0_dp*reflection_z_a - ra(3)
         END IF

         DO iset = 1, nseta
            offset_a_start = offset_a_end
            offset_a_end = offset_a_end + nsgfa(iset)
            sgfa = first_sgfa(1, iset)

            offset_b_end = 0
            DO jatom = atom_start, atom_end
               jkind = kind_of(jatom)
               atom_b = atom_of_kind(jatom)
               CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)

               first_sgfb => basis_set_b%first_sgf
               lb_max => basis_set_b%lmax
               lb_min => basis_set_b%lmin
               nsetb = basis_set_b%nset
               nsgfb => basis_set_b%nsgf_set
               sphi_b => basis_set_b%sphi
               zetb => basis_set_b%zet
               npgfb => basis_set_b%npgf

               rb(:) = pbc(particle_set(jatom)%r, cell)

               IF (my_do_reflection_b) THEN
                  rb(3) = 2.0_dp*reflection_z_b - rb(3)
               END IF

               rab(:) = ra(:) - rb(:) ! pbc not needed?
               dab = SQRT(rab(1)**2 + rab(2)**2 + rab(3)**2)

               jset_start = 1; jset_end = nsetb
               IF (jatom == atom_start) jset_start = set_start
               IF (jatom == atom_end) jset_end = set_end

               DO jset = jset_start, jset_end
                  first_set = 1; last_set = nsgfb(jset)
                  IF (jset == jset_start .AND. jatom == atom_start) first_set = set_offset_start
                  IF (jset == jset_end .AND. jatom == atom_end) last_set = set_offset_end

                  offset_b_start = offset_b_end
                  offset_b_end = offset_b_end + last_set + 1 - first_set
                  sgfb = first_sgfb(1, jset)
                  my_setpair = my_setpair + 1
                  map_it_here = MODULO(my_setpair, para_env%num_pe) == para_env%mepos

                  IF (map_it_here) THEN
#!some fypp magic to deal with combinations of optional arguments
#:for doforce_1 in [0, 1]
#:for doforce_2 in [0, 1]
                     IF (${conditional(doforce_1)}$PRESENT(force_a) .AND. &
                         ${conditional(doforce_2)}$PRESENT(force_b)) THEN

                        CALL integrate_set_2c( &
                           param%par, potential_parameter, &
                           la_min(iset), la_max(iset), &
                           lb_min(jset), lb_max(jset), &
                           npgfa(iset), npgfb(jset), &
                           zeta(:, iset), zetb(:, jset), &
                           ra, rb, &
                           hab, nsgfa(iset), last_set - first_set + 1, &
                           offset_a_start, offset_b_start, &
                           0, first_set - 1, &
                           sphi_a, sphi_b, &
                           sgfa, sgfb, nsgfa(iset), nsgfb(jset), &
                           my_eri_method, &
                           pab, &
                           $:                         'force_a=force_a(ikind)%forces(:, atom_a), &'*doforce_1
                           $:                         'force_b=force_b(jkind)%forces(:, atom_b), &'*doforce_2
                           hdab=hdab, hadb=hadb, &
                           G_count=G_count, R_count=R_count, &
                           do_reflection_a=do_reflection_a, do_reflection_b=do_reflection_b)
                     ENDIF
#:endfor
#:endfor
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF (my_eri_method == do_eri_mme) THEN

         CALL cp_eri_mme_update_local_counts(param, para_env, G_count_2c=G_count, R_count_2c=R_count)

      END IF

      CALL mp_sum(hab, para_env%group)
      IF (PRESENT(hdab)) CALL mp_sum(hdab, para_env%group)
      IF (PRESENT(hadb)) CALL mp_sum(hadb, para_env%group)

      CALL timestop(handle)
   END SUBROUTINE mp2_eri_2c_integrate

! **************************************************************************************************
!> \brief Integrate set pair and contract with sphi matrix.
!> \param param ...
!> \param potential_parameter ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfa ...
!> \param npgfb ...
!> \param zeta ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param hab ...
!> \param n_hab_a ...
!> \param n_hab_b ...
!> \param offset_hab_a ...
!> \param offset_hab_b ...
!> \param offset_set_a ...
!> \param offset_set_b ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sgfa ...
!> \param sgfb ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param eri_method ...
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \param hdab ...
!> \param hadb ...
!> \param G_count ...
!> \param R_count ...
!> \param do_reflection_a ...
!> \param do_reflection_b ...
! **************************************************************************************************
   SUBROUTINE integrate_set_2c(param, potential_parameter, la_min, la_max, lb_min, lb_max, npgfa, npgfb, zeta, zetb, &
                               ra, rb, hab, n_hab_a, n_hab_b, offset_hab_a, offset_hab_b, &
                               offset_set_a, offset_set_b, sphi_a, sphi_b, sgfa, sgfb, nsgfa, nsgfb, &
                               eri_method, pab, force_a, force_b, hdab, hadb, G_count, R_count, &
                               do_reflection_a, do_reflection_b)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, npgfa
      REAL(KIND=dp), DIMENSION(npgfa), INTENT(IN)        :: zeta
      INTEGER, INTENT(IN)                                :: npgfb
      REAL(KIND=dp), DIMENSION(npgfb), INTENT(IN)        :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hab
      INTEGER, INTENT(IN)                                :: n_hab_a, n_hab_b, offset_hab_a, &
                                                            offset_hab_b, offset_set_a, &
                                                            offset_set_b
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_a
      INTEGER, INTENT(IN)                                :: sgfa, nsgfa
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_b
      INTEGER, INTENT(IN)                                :: sgfb, nsgfb, eri_method
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: hdab, hadb
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count, R_count
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_reflection_a, do_reflection_b

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_set_2c'

      INTEGER :: ax, ay, az, bx, by, bz, hab_a_end, hab_a_start, hab_b_end, hab_b_start, handle, &
                 i_xyz, ico, icox, icoy, icoz, ipgf, jco, jcox, jcoy, jcoz, jpgf, la, la_max_d, lb, &
                 lb_max_d, na, nb, ncoa, ncob, set_a_end, set_a_start, set_b_end, set_b_start, &
                 sphi_a_start, sphi_b_start
      INTEGER, DIMENSION(3)                              :: la_xyz, lb_xyz
      LOGICAL                                            :: calculate_forces, my_do_reflection_a, &
                                                            my_do_reflection_b, do_force_a, do_force_b
      REAL(KIND=dp)                                      :: rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hab_contr, hab_uncontr, &
                                                            hab_uncontr_d, pab_hh, pab_hs, &
                                                            pab_ss
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: hadb_contr, hadb_uncontr, hdab_contr, &
                                                            hdab_uncontr, v_work

      ! note: tested only for one exponent per pair (npgfa = npgfb = 1)
      CALL timeset(routineN, handle)

      my_do_reflection_a = .FALSE.
      IF (PRESENT(do_reflection_a)) my_do_reflection_a = do_reflection_a

      my_do_reflection_b = .FALSE.
      IF (PRESENT(do_reflection_b)) my_do_reflection_b = do_reflection_b

      do_force_a = PRESENT(force_a) .OR. PRESENT(hdab)
      do_force_b = PRESENT(force_b) .OR. PRESENT(hadb)
      calculate_forces = do_force_a .OR. do_force_b

      IF (PRESENT(force_a) .OR. PRESENT(force_b)) THEN
         CPASSERT(PRESENT(pab))
         CPASSERT(ALL(SHAPE(pab) .EQ. SHAPE(hab)))
      ENDIF

      la_max_d = la_max
      lb_max_d = lb_max

      IF (calculate_forces) THEN
         IF (do_force_a) la_max_d = la_max + 1
         IF (do_force_b) lb_max_d = lb_max + 1
      ENDIF

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)

      rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2

      ALLOCATE (hab_uncontr_d(ncoset(la_max_d), ncoset(lb_max_d))); hab_uncontr_d(:, :) = 0.0_dp
      ALLOCATE (hab_uncontr(ncoa, ncob)); hab_uncontr(:, :) = 0.0_dp
      IF (PRESENT(hdab)) THEN
         ALLOCATE (hdab_uncontr(3, ncoa, ncob)); hdab_uncontr(:, :, :) = 0.0_dp
      ENDIF
      IF (PRESENT(hadb)) THEN
         ALLOCATE (hadb_uncontr(3, ncoa, ncob)); hadb_uncontr(:, :, :) = 0.0_dp
      ENDIF

      hab_a_start = offset_hab_a + 1; hab_a_end = offset_hab_a + n_hab_a
      hab_b_start = offset_hab_b + 1; hab_b_end = offset_hab_b + n_hab_b
      set_a_start = offset_set_a + 1; set_a_end = offset_set_a + n_hab_a
      set_b_start = offset_set_b + 1; set_b_end = offset_set_b + n_hab_b

      IF (eri_method == do_eri_mme) THEN
         CALL eri_mme_set_potential(param, convert_potential_type(potential_parameter%potential_type), potential_parameter%omega)

         IF (calculate_forces .AND. PRESENT(pab)) THEN
            ! uncontracted hermite-gaussian representation of density matrix
            sphi_a_start = sgfa - 1 + set_a_start
            sphi_b_start = sgfb - 1 + set_b_start

            ALLOCATE (pab_ss(n_hab_a, n_hab_b))
            pab_ss(:, :) = pab(hab_a_start:hab_a_end, hab_b_start:hab_b_end)
            ALLOCATE (pab_hs(ncoa, n_hab_b)); ALLOCATE (pab_hh(ncoa, ncob))
            CALL dgemm("N", "N", ncoa, n_hab_b, n_hab_a, 1.0_dp, &
                       sphi_a(:, sphi_a_start), SIZE(sphi_a, 1), pab_ss, n_hab_a, 0.0_dp, pab_hs, ncoa)
            CALL dgemm("N", "T", ncoa, ncob, n_hab_b, 1.0_dp, &
                       pab_hs, ncoa, sphi_b(:, sphi_b_start), SIZE(sphi_b, 1), 0.0_dp, pab_hh, ncoa)
         ENDIF

         DO ipgf = 1, npgfa
            na = (ipgf - 1)*ncoset(la_max)
            DO jpgf = 1, npgfb
               nb = (jpgf - 1)*ncoset(lb_max)
               hab_uncontr_d(:, :) = 0.0_dp
               CALL eri_mme_2c_integrate(param, &
                                         la_min, la_max_d, lb_min, lb_max_d, &
                                         zeta(ipgf), zetb(jpgf), ra - rb, hab_uncontr_d, 0, 0, G_count, R_count)

               hab_uncontr(na + 1:na + ncoset(la_max), nb + 1:nb + ncoset(lb_max)) = &
                  hab_uncontr_d(:ncoset(la_max), :ncoset(lb_max))

               IF (calculate_forces) THEN
                  DO lb = lb_min, lb_max
                  DO bx = 0, lb
                  DO by = 0, lb - bx
                     bz = lb - bx - by
                     jco = coset(bx, by, bz)
                     jcox = coset(bx + 1, by, bz)
                     jcoy = coset(bx, by + 1, bz)
                     jcoz = coset(bx, by, bz + 1)
                     DO la = la_min, la_max
                     DO ax = 0, la
                     DO ay = 0, la - ax
                        az = la - ax - ay
                        la_xyz = [ax, ay, az]
                        lb_xyz = [bx, by, bz]
                        ico = coset(ax, ay, az)
                        icox = coset(ax + 1, ay, az)
                        icoy = coset(ax, ay + 1, az)
                        icoz = coset(ax, ay, az + 1)
                        IF (PRESENT(force_a)) THEN
                           force_a(:) = force_a(:) + 2.0_dp*zeta(ipgf)* &
                                        [pab_hh(na + ico, nb + jco)*hab_uncontr_d(icox, jco), &
                                         pab_hh(na + ico, nb + jco)*hab_uncontr_d(icoy, jco), &
                                         pab_hh(na + ico, nb + jco)*hab_uncontr_d(icoz, jco)]
                        ENDIF
                        IF (PRESENT(force_b)) THEN
                           force_b(:) = force_b(:) + 2.0_dp*zetb(jpgf)* &
                                        [pab_hh(na + ico, nb + jco)*hab_uncontr_d(ico, jcox), &
                                         pab_hh(na + ico, nb + jco)*hab_uncontr_d(ico, jcoy), &
                                         pab_hh(na + ico, nb + jco)*hab_uncontr_d(ico, jcoz)]
                        ENDIF
                        IF (PRESENT(hdab)) THEN
                           hdab_uncontr(1:3, na + ico, nb + jco) = 2.0_dp*zeta(ipgf)* &
                                                                   [hab_uncontr_d(icox, jco), &
                                                                    hab_uncontr_d(icoy, jco), &
                                                                    hab_uncontr_d(icoz, jco)]
                        ENDIF
                        IF (PRESENT(hadb)) THEN
                           hadb_uncontr(1:3, na + ico, nb + jco) = 2.0_dp*zetb(jpgf)* &
                                                                   [hab_uncontr_d(ico, jcox), &
                                                                    hab_uncontr_d(ico, jcoy), &
                                                                    hab_uncontr_d(ico, jcoz)]
                        ENDIF
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF

            ENDDO
         ENDDO

      ELSE IF (eri_method == do_eri_os) THEN

         IF (calculate_forces) CPABORT("NYI")

         ALLOCATE (f_work(0:la_max + lb_max + 2))
         ALLOCATE (v_work(ncoa, ncob, la_max + lb_max + 1))
         v_work(:, :, :) = 0.0_dp
         f_work = 0.0_dp

         CALL coulomb2_new(la_max, npgfa, zeta, la_min, lb_max, npgfb, zetb, lb_min, &
                           rb - ra, rab2, hab_uncontr, v_work, f_work)

         DEALLOCATE (v_work, f_work)

      ELSE IF (eri_method == do_eri_gpw) THEN

         CPABORT("GPW not enabled in the ERI interface.")

      END IF

      ALLOCATE (hab_contr(nsgfa, nsgfb))
      IF (PRESENT(hdab)) THEN
         ALLOCATE (hdab_contr(3, nsgfa, nsgfb))
      ENDIF
      IF (PRESENT(hadb)) THEN
         ALLOCATE (hadb_contr(3, nsgfa, nsgfb))
      ENDIF

      CALL ab_contract(hab_contr, hab_uncontr, sphi_a(:, sgfa:), sphi_b(:, sgfb:), ncoa, ncob, nsgfa, nsgfb)

      IF (calculate_forces) THEN
         DO i_xyz = 1, 3
            IF (PRESENT(hdab)) THEN
               CALL ab_contract(hdab_contr(i_xyz, :, :), hdab_uncontr(i_xyz, :, :), &
                                sphi_a(:, sgfa:), sphi_b(:, sgfb:), ncoa, ncob, nsgfa, nsgfb)
            ENDIF
            IF (PRESENT(hadb)) THEN
               CALL ab_contract(hadb_contr(i_xyz, :, :), hadb_uncontr(i_xyz, :, :), &
                                sphi_a(:, sgfa:), sphi_b(:, sgfb:), ncoa, ncob, nsgfa, nsgfb)
            ENDIF
         ENDDO
      ENDIF

      hab(hab_a_start:hab_a_end, hab_b_start:hab_b_end) = hab_contr(set_a_start:set_a_end, set_b_start:set_b_end)

      IF (calculate_forces) THEN
         IF (PRESENT(hdab)) hdab(:, hab_a_start:hab_a_end, hab_b_start:hab_b_end) = &
            hdab_contr(:, set_a_start:set_a_end, set_b_start:set_b_end)
         IF (PRESENT(hadb)) hadb(:, hab_a_start:hab_a_end, hab_b_start:hab_b_end) = &
            hadb_contr(:, set_a_start:set_a_end, set_b_start:set_b_end)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE integrate_set_2c

! **************************************************************************************************
!> \brief high-level integration routine for 3c integrals (ab|c) over CP2K basis sets.
!>        For each local function of c, (ab|c) is written to a DBCSR matrix mat_ab.
!> \param param ...
!> \param potential_parameter ...
!> \param para_env ...
!> \param qs_env ...
!> \param first_c start index of local range of c
!> \param last_c end index of local range of c
!> \param mat_ab DBCSR matrices for each c
!> \param basis_type_a ...
!> \param basis_type_b ...
!> \param basis_type_c ...
!> \param sab_nl neighbor list for a, b
!> \param eri_method ...
!> \param pabc ...
!> \param force_a ...
!> \param force_b ...
!> \param force_c ...
!> \param mat_dabc ...
!> \param mat_adbc ...
!> \param mat_abdc ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_3c_integrate(param, potential_parameter, para_env, qs_env, &
                                   first_c, last_c, mat_ab, &
                                   basis_type_a, basis_type_b, basis_type_c, &
                                   sab_nl, eri_method, &
                                   pabc, force_a, force_b, force_c, &
                                   mat_dabc, mat_adbc, mat_abdc)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(IN)                                :: first_c, last_c
      TYPE(dbcsr_p_type), DIMENSION(last_c - first_c + 1), &
         INTENT(INOUT)                                   :: mat_ab
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type_a, basis_type_b, basis_type_c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, INTENT(IN), OPTIONAL                      :: eri_method
      TYPE(dbcsr_p_type), DIMENSION(last_c - first_c + 1), &
         INTENT(INOUT), OPTIONAL                         :: pabc
      TYPE(mp2_eri_force), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT), OPTIONAL             :: force_a, force_b, force_c
      TYPE(dbcsr_p_type), &
         DIMENSION(3, last_c - first_c + 1), INTENT(INOUT), &
         OPTIONAL                                        :: mat_dabc, mat_adbc, mat_abdc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_eri_3c_integrate'

      INTEGER :: atom_a, atom_b, atom_c, atom_end, atom_start, first_set, GG_count, GR_count, &
                 handle, i_xyz, iatom, ic, icol, ikind, inode, irow, iset, jatom, jkind, jset, katom, &
                 kkind, kset, kset_end, kset_start, last_jatom, last_set, mepos, my_eri_method, na, natom, &
                 nb, nc, nkind, nseta, nsetb, nsetc, nthread, offset_a_end, offset_a_start, offset_b_end, &
                 offset_b_start, offset_c_end, offset_c_start, RR_count, set_end, set_offset_end, &
                 set_offset_start, set_start, sgfa, sgfb, sgfc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, natom_of_kind
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: eri_offsets
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: calculate_forces, do_symmetric, found, &
                                                            to_be_asserted
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: habc, pabc_block
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: habdc, hadbc, hdabc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: munu_block, pab_block, rpgfb, sphi_a, &
                                                            sphi_b, sphi_c, zeta, zetb, zetc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      calculate_forces = PRESENT(force_a) .OR. PRESENT(force_b) .OR. PRESENT(force_c) .OR. &
                         PRESENT(mat_dabc) .OR. PRESENT(mat_adbc) .OR. PRESENT(mat_abdc)

      my_eri_method = do_eri_mme
      IF (PRESENT(eri_method)) my_eri_method = eri_method

      IF (PRESENT(force_a) .OR. PRESENT(force_b) .OR. PRESENT(force_c)) THEN
         CPASSERT(PRESENT(pabc))
      ENDIF

      GG_count = 0; GR_count = 0; RR_count = 0

      nthread = 1

      ! get mapping between ERIs and atoms, sets, set offsets
      CALL get_eri_offsets(qs_env, basis_type_c, eri_offsets)

      atom_start = eri_offsets(first_c, 1)
      set_start = eri_offsets(first_c, 2)
      set_offset_start = eri_offsets(first_c, 3)

      atom_end = eri_offsets(last_c, 1)
      set_end = eri_offsets(last_c, 2)
      set_offset_end = eri_offsets(last_c, 3)

      ! get QS stuff
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      natom=natom, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      nkind=nkind)

      ALLOCATE (atom_of_kind(natom))
      ALLOCATE (kind_of(natom))
      ALLOCATE (natom_of_kind(nkind))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of, natom_of_kind=natom_of_kind)

      IF (PRESENT(force_a)) CALL mp2_eri_allocate_forces(force_a, natom_of_kind)
      IF (PRESENT(force_b)) CALL mp2_eri_allocate_forces(force_b, natom_of_kind)
      IF (PRESENT(force_c)) CALL mp2_eri_allocate_forces(force_c, natom_of_kind)

      nc = last_c - first_c + 1

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      IF (do_symmetric) THEN
         CPASSERT(basis_type_a == basis_type_b)
      END IF

      ALLOCATE (basis_set_list_a(nkind), basis_set_list_b(nkind))
      CALL basis_set_list_setup(basis_set_list_a, basis_type_a, qs_kind_set)
      CALL basis_set_list_setup(basis_set_list_b, basis_type_b, qs_kind_set)

      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

      mepos = 0

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)

         ! exclude periodic images because method is periodic intrinsically
         IF (inode == 1) last_jatom = 0

         IF (jatom /= last_jatom) THEN
            last_jatom = jatom
         ELSE
            CYCLE
         END IF

         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         na = SUM(nsgfa)

         ra(:) = pbc(particle_set(iatom)%r, cell)

         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         nb = SUM(nsgfb)

         rb(:) = pbc(particle_set(jatom)%r, cell)

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         ALLOCATE (habc(na, nb, nc))
         habc(:, :, :) = 0.0_dp ! needs to be initialized due to screening
         IF (PRESENT(mat_dabc)) THEN
            ALLOCATE (hdabc(3, na, nb, nc))
            hdabc(:, :, :, :) = 0.0_dp
         ENDIF
         IF (PRESENT(mat_adbc)) THEN
            ALLOCATE (hadbc(3, na, nb, nc))
            hadbc(:, :, :, :) = 0.0_dp
         ENDIF
         IF (PRESENT(mat_abdc)) THEN
            ALLOCATE (habdc(3, na, nb, nc))
            habdc(:, :, :, :) = 0.0_dp
         ENDIF

         IF (calculate_forces .AND. PRESENT(pabc)) THEN
            ALLOCATE (pabc_block(na, nb, nc))
            DO ic = 1, nc
               NULLIFY (pab_block)
               CALL dbcsr_get_block_p(matrix=pabc(ic)%matrix, &
                                      row=irow, col=icol, block=pab_block, found=found)
               CPASSERT(found)
               IF (irow .EQ. iatom) THEN
                  to_be_asserted = SIZE(pab_block, 1) .EQ. SIZE(pabc_block, 1) .AND. SIZE(pab_block, 2) .EQ. SIZE(pabc_block, 2)
                  CPASSERT(to_be_asserted)
                  pabc_block(:, :, ic) = pab_block(:, :)
               ELSE
                  to_be_asserted = SIZE(pab_block, 2) .EQ. SIZE(pabc_block, 1) .AND. SIZE(pab_block, 1) .EQ. SIZE(pabc_block, 2)
                  CPASSERT(to_be_asserted)
                  pabc_block(:, :, ic) = TRANSPOSE(pab_block(:, :))
               ENDIF
            ENDDO
         ENDIF

         rab(:) = pbc(rab, cell)
         dab = SQRT(rab(1)**2 + rab(2)**2 + rab(3)**2)

         offset_a_end = 0
         DO iset = 1, nseta
            offset_a_start = offset_a_end
            offset_a_end = offset_a_end + nsgfa(iset)
            sgfa = first_sgfa(1, iset)

            offset_b_end = 0
            DO jset = 1, nsetb
               offset_b_start = offset_b_end
               offset_b_end = offset_b_end + nsgfb(jset)

               sgfb = first_sgfb(1, jset)

               ! Screening
               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               offset_c_end = 0
               DO katom = atom_start, atom_end

                  atom_c = atom_of_kind(katom)

                  kkind = kind_of(katom)
                  CALL get_qs_kind(qs_kind=qs_kind_set(kkind), basis_set=basis_set_c, basis_type=basis_type_c)
                  first_sgfc => basis_set_c%first_sgf
                  lc_max => basis_set_c%lmax
                  lc_min => basis_set_c%lmin
                  nsetc = basis_set_c%nset
                  nsgfc => basis_set_c%nsgf_set
                  sphi_c => basis_set_c%sphi
                  zetc => basis_set_c%zet
                  npgfc => basis_set_c%npgf

                  rc(:) = pbc(particle_set(katom)%r, cell)

                  kset_start = 1; kset_end = nsetc
                  IF (katom == atom_start) kset_start = set_start
                  IF (katom == atom_end) kset_end = set_end

                  DO kset = kset_start, kset_end
                     first_set = 1; last_set = nsgfc(kset)
                     IF (kset == kset_start .AND. katom == atom_start) first_set = set_offset_start
                     IF (kset == kset_end .AND. katom == atom_end) last_set = set_offset_end

                     offset_c_start = offset_c_end
                     offset_c_end = offset_c_end + last_set + 1 - first_set
                     sgfc = first_sgfc(1, kset)

#!some fypp magic to deal with combinations of optional arguments
#:for pabc_present in [0, 1]
#:for doforce_1 in [0, 1]
#:for doforce_2 in [0, 1]
#:for doforce_3 in [0, 1]
#:for dabc in [0, 1]
#:for adbc in [0, 1]
#:for abdc in [0, 1]
                     IF (${conditional(doforce_1)}$PRESENT(force_a) .AND. &
                         ${conditional(doforce_2)}$PRESENT(force_b) .AND. &
                         ${conditional(doforce_3)}$PRESENT(force_c) .AND. &
                         ${conditional(pabc_present)}$PRESENT(pabc) .AND. &
                         ${conditional(dabc)}$PRESENT(mat_dabc) .AND. &
                         ${conditional(adbc)}$PRESENT(mat_adbc) .AND. &
                         ${conditional(abdc)}$PRESENT(mat_abdc)) THEN
                        CALL integrate_set_3c( &
                           param%par, potential_parameter, &
                           la_min(iset), la_max(iset), &
                           lb_min(jset), lb_max(jset), &
                           lc_min(kset), lc_max(kset), &
                           npgfa(iset), npgfb(jset), npgfc(kset), &
                           zeta(:, iset), zetb(:, jset), zetc(:, kset), &
                           ra, rb, rc, &
                           habc, &
                           nsgfa(iset), nsgfb(jset), last_set - first_set + 1, &
                           offset_a_start, offset_b_start, offset_c_start, &
                           0, 0, first_set - 1, &
                           sphi_a, sphi_b, sphi_c, &
                           sgfa, sgfb, sgfc, &
                           nsgfa(iset), nsgfb(jset), nsgfc(kset), &
                           my_eri_method, &
                           $:                         'pabc=pabc_block, &'*pabc_present
                           $:                         'force_a=force_a(ikind)%forces(:, atom_a), &'*doforce_1
                           $:                         'force_b=force_b(jkind)%forces(:, atom_b), &'*doforce_2
                           $:                         'force_c=force_c(kkind)%forces(:, atom_c), &'*doforce_3
                           do_symmetric=do_symmetric, &
                           on_diagonal=iatom .EQ. jatom, &
                           $:                         'hdabc=hdabc, &'*dabc
                           $:                         'hadbc=hadbc, &'*adbc
                           $:                         'habdc=habdc, &'*abdc
                           GG_count=GG_count, GR_count=GR_count, RR_count=RR_count)
                     ENDIF
#:endfor
#:endfor
#:endfor
#:endfor
#:endfor
#:endfor
#:endfor
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF (calculate_forces .AND. PRESENT(pabc)) DEALLOCATE (pabc_block)
         DO ic = 1, nc
            NULLIFY (munu_block)
            CALL dbcsr_get_block_p(matrix=mat_ab(ic)%matrix, &
                                   row=irow, col=icol, block=munu_block, found=found)
            CPASSERT(found)
            munu_block(:, :) = 0.0_dp
            IF (irow .EQ. iatom) THEN
               to_be_asserted = SIZE(munu_block, 1) .EQ. SIZE(habc, 1) .AND. SIZE(munu_block, 2) .EQ. SIZE(habc, 2)
               CPASSERT(to_be_asserted)
               munu_block(:, :) = habc(:, :, ic)
            ELSE
               to_be_asserted = SIZE(munu_block, 2) .EQ. SIZE(habc, 1) .AND. SIZE(munu_block, 1) .EQ. SIZE(habc, 2)
               CPASSERT(to_be_asserted)
               munu_block(:, :) = TRANSPOSE(habc(:, :, ic))
            ENDIF
         ENDDO
         DEALLOCATE (habc)
         IF (calculate_forces) THEN
            DO ic = 1, nc
               DO i_xyz = 1, 3
                  IF (PRESENT(mat_dabc)) THEN
                     NULLIFY (munu_block)
                     CALL dbcsr_get_block_p(matrix=mat_dabc(i_xyz, ic)%matrix, &
                                            row=irow, col=icol, block=munu_block, found=found)
                     CPASSERT(found)
                     munu_block(:, :) = 0.0_dp
                     IF (irow .EQ. iatom) THEN
                        munu_block(:, :) = hdabc(i_xyz, :, :, ic)
                     ELSE
                        munu_block(:, :) = TRANSPOSE(hdabc(i_xyz, :, :, ic))
                     ENDIF
                  ENDIF
                  IF (PRESENT(mat_adbc)) THEN
                     NULLIFY (munu_block)
                     CALL dbcsr_get_block_p(matrix=mat_adbc(i_xyz, ic)%matrix, &
                                            row=irow, col=icol, block=munu_block, found=found)
                     CPASSERT(found)
                     munu_block(:, :) = 0.0_dp
                     IF (irow .EQ. iatom) THEN
                        munu_block(:, :) = hadbc(i_xyz, :, :, ic)
                     ELSE
                        munu_block(:, :) = TRANSPOSE(hadbc(i_xyz, :, :, ic))
                     ENDIF
                  ENDIF
                  IF (PRESENT(mat_abdc)) THEN
                     NULLIFY (munu_block)
                     CALL dbcsr_get_block_p(matrix=mat_abdc(i_xyz, ic)%matrix, &
                                            row=irow, col=icol, block=munu_block, found=found)
                     CPASSERT(found)
                     munu_block(:, :) = 0.0_dp
                     IF (irow .EQ. iatom) THEN
                        munu_block(:, :) = habdc(i_xyz, :, :, ic)
                     ELSE
                        munu_block(:, :) = TRANSPOSE(habdc(i_xyz, :, :, ic))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            IF (PRESENT(mat_dabc)) DEALLOCATE (hdabc)
            IF (PRESENT(mat_adbc)) DEALLOCATE (hadbc)
            IF (PRESENT(mat_abdc)) DEALLOCATE (habdc)
         ENDIF
      ENDDO

      DEALLOCATE (basis_set_list_a, basis_set_list_b)
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL cp_eri_mme_update_local_counts(param, para_env, GG_count_3c=GG_count, GR_count_3c=GR_count, RR_count_3c=RR_count)

      CALL timestop(handle)
   END SUBROUTINE mp2_eri_3c_integrate

! **************************************************************************************************
!> \brief Integrate set triple and contract with sphi matrix
!> \param param ...
!> \param potential_parameter ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param habc ...
!> \param n_habc_a ...
!> \param n_habc_b ...
!> \param n_habc_c ...
!> \param offset_habc_a ...
!> \param offset_habc_b ...
!> \param offset_habc_c ...
!> \param offset_set_a ...
!> \param offset_set_b ...
!> \param offset_set_c ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sgfa ...
!> \param sgfb ...
!> \param sgfc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param eri_method ...
!> \param pabc ...
!> \param force_a ...
!> \param force_b ...
!> \param force_c ...
!> \param do_symmetric ...
!> \param on_diagonal ...
!> \param hdabc ...
!> \param hadbc ...
!> \param habdc ...
!> \param GG_count ...
!> \param GR_count ...
!> \param RR_count ...
!> \note
! **************************************************************************************************
   SUBROUTINE integrate_set_3c(param, potential_parameter, &
                               la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                               npgfa, npgfb, npgfc, &
                               zeta, zetb, zetc, &
                               ra, rb, rc, &
                               habc, &
                               n_habc_a, n_habc_b, n_habc_c, &
                               offset_habc_a, offset_habc_b, offset_habc_c, &
                               offset_set_a, offset_set_b, offset_set_c, &
                               sphi_a, sphi_b, sphi_c, &
                               sgfa, sgfb, sgfc, &
                               nsgfa, nsgfb, nsgfc, &
                               eri_method, &
                               pabc, &
                               force_a, force_b, force_c, &
                               do_symmetric, on_diagonal, &
                               hdabc, hadbc, habdc, &
                               GG_count, GR_count, RR_count)

      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max, npgfa, npgfb, npgfc
      REAL(KIND=dp), DIMENSION(npgfa), INTENT(IN)        :: zeta
      REAL(KIND=dp), DIMENSION(npgfb), INTENT(IN)        :: zetb
      REAL(KIND=dp), DIMENSION(npgfc), INTENT(IN)        :: zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN) :: n_habc_a, n_habc_b, n_habc_c, offset_habc_a, offset_habc_b, &
                             offset_habc_c, offset_set_a, offset_set_b, offset_set_c
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_a, sphi_b, sphi_c
      INTEGER, INTENT(IN)                                :: sgfa, sgfb, sgfc, nsgfa, nsgfb, nsgfc, &
                                                            eri_method
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         OPTIONAL                                        :: pabc
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b, force_c
      LOGICAL, INTENT(IN)                                :: do_symmetric
      LOGICAL, INTENT(IN), OPTIONAL                      :: on_diagonal
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(OUT), OPTIONAL                           :: hdabc, hadbc, habdc
      INTEGER, INTENT(INOUT), OPTIONAL                   :: GG_count, GR_count, RR_count

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_set_3c'

      INTEGER :: ax, ay, az, bx, by, bz, cx, cy, cz, habc_a_end, habc_a_start, habc_b_end, &
                 habc_b_start, habc_c_end, habc_c_start, handle, i_xyz, ico, icoc, icox, icoy, icoz, ipgf, &
                 jco, jcox, jcoy, jcoz, jpgf, kco, kcox, kcoy, kcoz, kpgf, la, la_max_d, lb, &
                 lb_max_d, lc, lc_max_d, na, nb, nc, nc_end, nc_start, ncoa, ncoa_d, ncob, &
                 ncob_d, ncoc, ncoc_d, set_a_end, set_a_start, set_b_end, set_b_start, &
                 set_c_end, set_c_start, sphi_a_start, sphi_b_start, sphi_c_start
      INTEGER, DIMENSION(3)                              :: la_xyz, lb_xyz
      LOGICAL                                            :: calculate_forces, do_force_a, do_force_b, &
                                                            do_force_c
      REAL(KIND=dp)                                      :: rab2, rac2, rbc2, w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f_work, gccc, rpgfa, rpgfb, rpgfc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pab_hh, pab_hs, vabc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: habc_contr, habc_uncontr, &
                                                            habc_uncontr_d, pabc_hhh, &
                                                            pabc_hsh, pabc_hss, pabc_sss
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: habdc_contr, habdc_uncontr, hadbc_contr, &
                                                            hadbc_uncontr, hdabc_contr, &
                                                            hdabc_uncontr, v_work

      CALL timeset(routineN, handle)

      do_force_a = PRESENT(force_a) .OR. PRESENT(hdabc)
      do_force_b = PRESENT(force_b) .OR. PRESENT(hadbc)
      do_force_c = PRESENT(force_c) .OR. PRESENT(habdc)
      calculate_forces = do_force_a .OR. do_force_b .OR. do_force_c

      IF (do_symmetric) THEN
         CPASSERT(PRESENT(on_diagonal))
      END IF

      la_max_d = la_max
      lb_max_d = lb_max
      lc_max_d = lc_max

      IF (calculate_forces) THEN
         IF (do_force_a) la_max_d = la_max + 1
         IF (do_force_b) lb_max_d = lb_max + 1
         IF (do_force_c) lc_max_d = lc_max + 1
      ENDIF

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)
      ncoc = npgfc*ncoset(lc_max)

      ncoa_d = npgfa*ncoset(la_max_d)
      ncob_d = npgfb*ncoset(lb_max_d)
      ncoc_d = npgfc*ncoset(lc_max_d)

      ALLOCATE (habc_uncontr_d(ncoset(la_max_d), ncoset(lb_max_d), ncoset(lc_max_d)))
      habc_uncontr_d(:, :, :) = 0.0_dp
      ALLOCATE (habc_uncontr(ncoa, ncob, ncoc)); habc_uncontr(:, :, :) = 0.0_dp
      IF (PRESENT(hdabc)) THEN
         ALLOCATE (hdabc_uncontr(3, ncoa, ncob, ncoc)); hdabc_uncontr(:, :, :, :) = 0.0_dp
      ENDIF
      IF (PRESENT(hadbc)) THEN
         ALLOCATE (hadbc_uncontr(3, ncoa, ncob, ncoc)); hadbc_uncontr(:, :, :, :) = 0.0_dp
      ENDIF
      IF (PRESENT(habdc)) THEN
         ALLOCATE (habdc_uncontr(3, ncoa, ncob, ncoc)); habdc_uncontr(:, :, :, :) = 0.0_dp
      ENDIF

      habc_a_start = offset_habc_a + 1; habc_a_end = offset_habc_a + n_habc_a
      habc_b_start = offset_habc_b + 1; habc_b_end = offset_habc_b + n_habc_b
      habc_c_start = offset_habc_c + 1; habc_c_end = offset_habc_c + n_habc_c
      set_a_start = offset_set_a + 1; set_a_end = offset_set_a + n_habc_a
      set_b_start = offset_set_b + 1; set_b_end = offset_set_b + n_habc_b
      set_c_start = offset_set_c + 1; set_c_end = offset_set_c + n_habc_c

      IF (eri_method == do_eri_mme) THEN
         CALL eri_mme_set_potential(param, convert_potential_type(potential_parameter%potential_type), potential_parameter%omega)

         IF (calculate_forces .AND. PRESENT(pabc)) THEN
            ! uncontracted hermite-gaussian representation of density matrix
            sphi_a_start = sgfa - 1 + set_a_start
            sphi_b_start = sgfb - 1 + set_b_start
            sphi_c_start = sgfc - 1 + set_c_start

            ALLOCATE (pabc_sss(n_habc_a, n_habc_b, n_habc_c))
            pabc_sss(:, :, :) = pabc(habc_a_start:habc_a_end, habc_b_start:habc_b_end, habc_c_start:habc_c_end)
            ALLOCATE (pabc_hss(ncoa, n_habc_b, n_habc_c))
            ALLOCATE (pabc_hsh(ncoa, n_habc_b, ncoc))
            ALLOCATE (pabc_hhh(ncoa, ncob, ncoc))
            ALLOCATE (pab_hs(ncoa, n_habc_b))
            ALLOCATE (pab_hh(ncoa, ncob))

            CALL dgemm("N", "N", ncoa, n_habc_b*n_habc_c, n_habc_a, 1.0_dp, &
                       sphi_a(:, sphi_a_start), SIZE(sphi_a, 1), pabc_sss, n_habc_a, 0.0_dp, pabc_hss, ncoa)
            CALL dgemm("N", "T", ncoa*n_habc_b, ncoc, n_habc_c, 1.0_dp, &
                       pabc_hss, ncoa*n_habc_b, sphi_c(:, sphi_c_start), SIZE(sphi_c, 1), 0.0_dp, pabc_hsh, ncoa*n_habc_b)

            DO icoc = 1, ncoc
               pab_hs(:, :) = pabc_hsh(:, :, icoc)
               CALL dgemm("N", "T", ncoa, ncob, n_habc_b, 1.0_dp, &
                          pab_hs, ncoa, sphi_b(:, sphi_b_start), SIZE(sphi_b, 1), 0.0_dp, pab_hh, ncoa)
               pabc_hhh(:, :, icoc) = pab_hh(:, :)
            ENDDO
         ENDIF

         DO ipgf = 1, npgfa
            na = (ipgf - 1)*ncoset(la_max)
            DO jpgf = 1, npgfb
               nb = (jpgf - 1)*ncoset(lb_max)
               DO kpgf = 1, npgfc
                  nc = (kpgf - 1)*ncoset(lc_max)
                  habc_uncontr_d(:, :, :) = 0.0_dp
                  CALL eri_mme_3c_integrate(param, &
                                            la_min, la_max_d, lb_min, lb_max_d, lc_min, lc_max_d, &
                                            zeta(ipgf), zetb(jpgf), zetc(kpgf), ra, rb, rc, habc_uncontr_d, 0, 0, 0, &
                                            GG_count, GR_count, RR_count)

                  habc_uncontr(na + 1:na + ncoset(la_max), nb + 1:nb + ncoset(lb_max), nc + 1:nc + ncoset(lc_max)) = &
                     habc_uncontr_d(:ncoset(la_max), :ncoset(lb_max), :ncoset(lc_max))

                  IF (calculate_forces) THEN
                     DO lc = lc_min, lc_max
                     DO cx = 0, lc
                     DO cy = 0, lc - cx
                        cz = lc - cx - cy
                        kco = coset(cx, cy, cz)
                        kcox = coset(cx + 1, cy, cz)
                        kcoy = coset(cx, cy + 1, cz)
                        kcoz = coset(cx, cy, cz + 1)
                        DO lb = lb_min, lb_max
                        DO bx = 0, lb
                        DO by = 0, lb - bx
                           bz = lb - bx - by
                           jco = coset(bx, by, bz)
                           jcox = coset(bx + 1, by, bz)
                           jcoy = coset(bx, by + 1, bz)
                           jcoz = coset(bx, by, bz + 1)
                           DO la = la_min, la_max
                           DO ax = 0, la
                           DO ay = 0, la - ax
                              az = la - ax - ay
                              la_xyz = [ax, ay, az]
                              lb_xyz = [bx, by, bz]
                              ico = coset(ax, ay, az)
                              icox = coset(ax + 1, ay, az)
                              icoy = coset(ax, ay + 1, az)
                              icoz = coset(ax, ay, az + 1)

                              w = 1.0_dp
                              IF (do_symmetric .AND. .NOT. on_diagonal) w = 2.0_dp

                              IF (PRESENT(force_a)) THEN
                                 force_a = force_a + 2.0_dp*w*zeta(ipgf)* &
                                           [pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(icox, jco, kco), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(icoy, jco, kco), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(icoz, jco, kco)]

                              ENDIF
                              IF (PRESENT(force_b)) THEN
                                 force_b = force_b + 2.0_dp*w*zetb(jpgf)* &
                                           [pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jcox, kco), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jcoy, kco), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jcoz, kco)]
                              ENDIF
                              IF (PRESENT(force_c)) THEN
                                 force_c = force_c + 2.0_dp*w*zetc(kpgf)* &
                                           [pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jco, kcox), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jco, kcoy), &
                                            pabc_hhh(na + ico, nb + jco, nc + kco)*habc_uncontr_d(ico, jco, kcoz)]
                              ENDIF

                              IF (PRESENT(hdabc)) THEN
                                 hdabc_uncontr(1:3, na + ico, nb + jco, nc + kco) = 2.0_dp*zeta(ipgf)* &
                                                                                    [habc_uncontr_d(icox, jco, kco), &
                                                                                     habc_uncontr_d(icoy, jco, kco), &
                                                                                     habc_uncontr_d(icoz, jco, kco)]
                              ENDIF
                              IF (PRESENT(hadbc)) THEN
                                 hadbc_uncontr(1:3, na + ico, nb + jco, nc + kco) = 2.0_dp*zetb(jpgf)* &
                                                                                    [habc_uncontr_d(ico, jcox, kco), &
                                                                                     habc_uncontr_d(ico, jcoy, kco), &
                                                                                     habc_uncontr_d(ico, jcoz, kco)]
                              ENDIF
                              IF (PRESENT(habdc)) THEN
                                 habdc_uncontr(1:3, na + ico, nb + jco, nc + kco) = 2.0_dp*zetc(kpgf)* &
                                                                                    [habc_uncontr_d(ico, jco, kcox), &
                                                                                     habc_uncontr_d(ico, jco, kcoy), &
                                                                                     habc_uncontr_d(ico, jco, kcoz)]
                              ENDIF
                           ENDDO
                           ENDDO
                           ENDDO
                        ENDDO
                        ENDDO
                        ENDDO
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

      ELSE IF (eri_method == do_eri_os) THEN

         IF (calculate_forces) CPABORT("NYI")

         ALLOCATE (f_work(0:la_max + lb_max + lc_max + 2))
         f_work(:) = 0.0_dp
         ALLOCATE (v_work(ncoa, ncob, ncoc, la_max + lb_max + lc_max + 1))
         v_work(:, :, :, :) = 0.0_dp
         ! no screening
         ALLOCATE (rpgfa(npgfa))
         ALLOCATE (rpgfb(npgfb))
         ALLOCATE (rpgfc(npgfc))
         rpgfa(:) = 1.0E10_dp
         rpgfb(:) = 1.0E10_dp
         rpgfc(:) = 1.0E10_dp
         ALLOCATE (gccc(ncoc))
         gccc(:) = 0.0_dp
         ALLOCATE (vabc(ncoa, ncob))
         vabc(:, :) = 0.0_dp
         rab2 = (rb(1) - ra(1))**2 + (rb(2) - ra(2))**2 + (rb(3) - ra(3))**2
         rac2 = (rc(1) - ra(1))**2 + (rc(2) - ra(2))**2 + (rc(3) - ra(3))**2
         rbc2 = (rc(1) - rb(1))**2 + (rc(2) - rb(2))**2 + (rc(3) - rb(3))**2

         ! in the RI basis, there is only a single primitive Gaussian
         kpgf = 1

         IF (lc_max == 0) THEN
            nc_start = 1
         ELSE
            nc_start = ncoset(lc_max - 1) + 1
         END IF
         nc_end = ncoset(lc_max)

         CALL coulomb3(la_max, npgfa, zeta(:), rpgfa(:), la_min, &
                       lb_max, npgfb, zetb(:), rpgfb(:), lb_min, &
                       lc_max, zetc(kpgf), rpgfc(kpgf), lc_min, &
                       gccc, rb - ra, rab2, rc - ra, rac2, rbc2, &
                       vabc, habc_uncontr(:, :, nc_start:nc_end), v_work, f_work)

         DEALLOCATE (v_work, f_work, rpgfa, rpgfb, rpgfc, gccc, vabc)

      ELSE IF (eri_method == do_eri_gpw) THEN

         CPABORT("GPW not enabled in the ERI interface.")

      END IF

      ALLOCATE (habc_contr(nsgfa, nsgfb, nsgfc))
      IF (PRESENT(hdabc)) THEN
         ALLOCATE (hdabc_contr(3, nsgfa, nsgfb, nsgfc))
      ENDIF
      IF (PRESENT(hadbc)) THEN
         ALLOCATE (hadbc_contr(3, nsgfa, nsgfb, nsgfc))
      ENDIF
      IF (PRESENT(habdc)) THEN
         ALLOCATE (habdc_contr(3, nsgfa, nsgfb, nsgfc))
      ENDIF

      CALL abc_contract(habc_contr, habc_uncontr, &
                        sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                        ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)

      IF (calculate_forces) THEN
         DO i_xyz = 1, 3
            IF (PRESENT(hdabc)) THEN
               CALL abc_contract(hdabc_contr(i_xyz, :, :, :), hdabc_uncontr(i_xyz, :, :, :), &
                                 sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                 ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)
            ENDIF
            IF (PRESENT(hadbc)) THEN
               CALL abc_contract(hadbc_contr(i_xyz, :, :, :), hadbc_uncontr(i_xyz, :, :, :), &
                                 sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                 ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)
            ENDIF
            IF (PRESENT(habdc)) THEN
               CALL abc_contract(habdc_contr(i_xyz, :, :, :), habdc_uncontr(i_xyz, :, :, :), &
                                 sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                 ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)
            ENDIF
         ENDDO
      ENDIF

      habc(habc_a_start:habc_a_end, habc_b_start:habc_b_end, habc_c_start:habc_c_end) = &
         habc_contr(set_a_start:set_a_end, set_b_start:set_b_end, set_c_start:set_c_end)

      IF (calculate_forces) THEN
         IF (PRESENT(hdabc)) hdabc(:, habc_a_start:habc_a_end, habc_b_start:habc_b_end, habc_c_start:habc_c_end) = &
            hdabc_contr(:, set_a_start:set_a_end, set_b_start:set_b_end, set_c_start:set_c_end)
         IF (PRESENT(hadbc)) hadbc(:, habc_a_start:habc_a_end, habc_b_start:habc_b_end, habc_c_start:habc_c_end) = &
            hadbc_contr(:, set_a_start:set_a_end, set_b_start:set_b_end, set_c_start:set_c_end)
         IF (PRESENT(habdc)) habdc(:, habc_a_start:habc_a_end, habc_b_start:habc_b_end, habc_c_start:habc_c_end) = &
            habdc_contr(:, set_a_start:set_a_end, set_b_start:set_b_end, set_c_start:set_c_end)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE integrate_set_3c

! **************************************************************************************************
!> \brief get pointer to atom, pointer to set and offset in a set for each spherical orbital of a
!>        basis.
!> \param qs_env ...
!> \param basis_type ...
!> \param eri_offsets (:,1) atom numbers
!>                    (:,2) set numbers
!>                    (:,3) set offsets
! **************************************************************************************************
   SUBROUTINE get_eri_offsets(qs_env, basis_type, eri_offsets)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: eri_offsets

      INTEGER                                            :: dimen_basis, iatom, ikind, iset, isgf, &
                                                            natom, nkind, nset, nsgf, offset, &
                                                            set_offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(:), POINTER                     :: nsgf_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, natom=natom, nkind=nkind)

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      dimen_basis = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type=basis_type)
         dimen_basis = dimen_basis + nsgf
      END DO

      ALLOCATE (eri_offsets(dimen_basis, 3))

      offset = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, basis_type=basis_type)
         nset = basis_set%nset
         nsgf_set => basis_set%nsgf_set
         DO iset = 1, nset
            set_offset = 0
            DO isgf = 1, nsgf_set(iset)
               set_offset = set_offset + 1
               eri_offsets(offset + set_offset, :) = [iatom, iset, set_offset]
            ENDDO
            offset = offset + nsgf_set(iset)
         ENDDO
      ENDDO
   END SUBROUTINE get_eri_offsets

! **************************************************************************************************
!> \brief ...
!> \param force ...
!> \param natom_of_kind ...
! **************************************************************************************************
   PURE SUBROUTINE mp2_eri_allocate_forces(force, natom_of_kind)
      TYPE(mp2_eri_force), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: force
      INTEGER, DIMENSION(:), INTENT(IN)                  :: natom_of_kind

      INTEGER                                            :: ikind, n, nkind

      nkind = SIZE(natom_of_kind)

      ALLOCATE (force(nkind))

      DO ikind = 1, nkind
         n = natom_of_kind(ikind)
         ALLOCATE (force(ikind)%forces(3, n))
         force(ikind)%forces(:, :) = 0.0_dp
      ENDDO
   END SUBROUTINE mp2_eri_allocate_forces

! **************************************************************************************************
!> \brief ...
!> \param force ...
! **************************************************************************************************
   PURE SUBROUTINE mp2_eri_deallocate_forces(force)
      TYPE(mp2_eri_force), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: force

      INTEGER                                            :: ikind, nkind

      IF (ALLOCATED(force)) THEN
         nkind = SIZE(force)
         DO ikind = 1, nkind
            IF (ALLOCATED(force(ikind)%forces)) DEALLOCATE (force(ikind)%forces)
         ENDDO

         DEALLOCATE (force)
      ENDIF
   END SUBROUTINE mp2_eri_deallocate_forces

   FUNCTION convert_potential_type(potential_type) RESULT(res)
      INTEGER, INTENT(IN)                                :: potential_type
      INTEGER                                            :: res

      IF (potential_type == do_potential_coulomb) THEN
         res = eri_mme_coulomb
      ELSE IF (potential_type == do_potential_long) THEN
         res = eri_mme_longrange
      ELSE
         CPABORT("MME potential not implemented!")
      END IF

   END FUNCTION convert_potential_type

END MODULE mp2_eri
