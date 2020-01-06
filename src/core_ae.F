!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Calculation of the nuclear attraction contribution to the core Hamiltonian
!>         <a|erfc|b> :we only calculate the non-screened part
!> \par History
!>      - core_ppnl refactored from qs_core_hamiltonian [Joost VandeVondele, 2008-11-01]
!>      - adapted for nuclear attraction [jhu, 2009-02-24]
! **************************************************************************************************
MODULE core_ae
   USE ai_verfc,                        ONLY: verfc
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE external_potential_types,        ONLY: all_potential_type,&
                                              get_potential,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              init_orbital_pointers,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              nl_set_sub_iterator,&
                                              nl_sub_iterate
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ae'

   PUBLIC :: build_core_ae

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_h ...
!> \param matrix_p ...
!> \param force ...
!> \param virial ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param nder ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param sab_orb ...
!> \param sac_ae ...
!> \param nimages ...
!> \param cell_to_index ...
! **************************************************************************************************
   SUBROUTINE build_core_ae(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                            qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ae, nimages, cell_to_index)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL                                            :: use_virial
      INTEGER                                            :: nder
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ae
      INTEGER, INTENT(IN)                                :: nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ae', routineP = moduleN//':'//routineN

      INTEGER :: atom_a, handle, iatom, icol, ikind, img, irow, iset, jatom, jkind, jset, katom, &
         kkind, ldai, ldsab, maxco, maxl, maxnset, maxsgf, mepos, na_plus, natom, nb_plus, ncoa, &
         ncob, nij, nkind, nseta, nsetb, nthread, sgfa, sgfb, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: dokp, found
      REAL(KIND=dp)                                      :: alpha_c, core_charge, core_radius, dab, &
                                                            dac, dbc, f0, rab2, rac2, rbc2, zeta_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: force_thread, habd, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: hab, pab, verf, vnuc
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, rab, rac, rbc
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_loc, pv_thread
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, p_block, rpgfa, rpgfb, sphi_a, &
                                                            sphi_b, zeta, zetb
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ap_iterator
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash, hash1, hash2
!$    INTEGER, PARAMETER                                 :: nlock = 501

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      ENDIF

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      dokp = (nimages > 1)

      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimages
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_p(2, img)%matrix, matrix_p(1, img)%matrix, &
                              alpha_scalar=-2.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
      END IF
      ALLOCATE (force_thread(3, natom))
      force_thread = 0.0_dp

      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF
      pv_thread = 0.0_dp

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, maxlgto=maxl, maxsgf=maxsgf, maxnset=maxnset)
      CALL init_orbital_pointers(maxl + nder + 1)
      ldsab = MAX(maxco, maxsgf)
      ldai = ncoset(maxl + nder + 1)

      nthread = 1
!$    nthread = omp_get_max_threads()

      ! iterator for basis/potential list
      CALL neighbor_list_iterator_create(ap_iterator, sac_ae, search=.TRUE., nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (ap_iterator, basis_set_list, calculate_forces, use_virial, &
!$OMP          matrix_h, matrix_p, atomic_kind_set, qs_kind_set, particle_set, &
!$OMP          sab_orb, sac_ae, nthread, ncoset, nkind, cell_to_index, &
!$OMP          slot, ldsab,  maxnset, ldai, nder, maxl, maxco, dokp, locks, natom) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, rab, basis_set_a, basis_set_b, &
!$OMP          first_sgfa, la_max, la_min, npgfa, nsgfa, sphi_a, &
!$OMP          zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, rpgfb, set_radius_b, sphi_b, &
!$OMP          zetb, zeta_c, alpha_c, core_charge, dab, irow, icol, h_block, found, iset, ncoa, &
!$OMP          sgfa, jset, ncob, sgfb, nsgfb, p_block, work, pab, hab, kkind, nseta, &
!$OMP          rac, dac, rbc, rab2, rac2, rbc2, dbc, na_plus, nb_plus, verf, vnuc, &
!$OMP          set_radius_a,  core_radius, rpgfa, force_a, force_b, mepos, &
!$OMP          habd, f0, katom, cellind, img, nij, ff, &
!$OMP          sgp_potential, all_potential, hash, hash1, hash2) &
!$OMP REDUCTION (+ : pv_thread, force_thread )

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

      mepos = 0
!$    mepos = omp_get_thread_num()

      ALLOCATE (hab(ldsab, ldsab, maxnset*maxnset), work(ldsab, ldsab))
      ALLOCATE (verf(ldai, ldai, 2*maxl + nder + 1), vnuc(ldai, ldai, 2*maxl + nder + 1), ff(0:2*maxl + nder))
      IF (calculate_forces) THEN
         ALLOCATE (pab(maxco, maxco, maxnset*maxnset))
      END IF

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sab_orb(1)%nl_size

         ikind = sab_orb(1)%nlist_task(slot)%ikind
         jkind = sab_orb(1)%nlist_task(slot)%jkind
         iatom = sab_orb(1)%nlist_task(slot)%iatom
         jatom = sab_orb(1)%nlist_task(slot)%jatom
         cellind(:) = sab_orb(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)

         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
!$       hash1 = MOD((iatom - 1)*natom + jatom, nlock) + 1
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
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

         dab = SQRT(SUM(rab*rab))

         IF (dokp) THEN
            img = cell_to_index(cellind(1), cellind(2), cellind(3))
         ELSE
            img = 1
         END IF

         ! *** Use the symmetry of the first derivatives ***
         IF (iatom == jatom) THEN
            f0 = 1.0_dp
         ELSE
            f0 = 2.0_dp
         END IF

         ! *** Create matrix blocks for a new matrix block column ***
         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF
         NULLIFY (h_block)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, img)%matrix, &
                                row=irow, col=icol, BLOCK=h_block, found=found)
         IF (calculate_forces) THEN
            NULLIFY (p_block)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, img)%matrix, &
                                   row=irow, col=icol, BLOCK=p_block, found=found)
            CPASSERT(ASSOCIATED(p_block))
            ! *** Decontract density matrix block ***
            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               DO jset = 1, nsetb
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)
                  nij = jset + (iset - 1)*maxnset
                  ! *** Decontract density matrix block ***
                  IF (iatom <= jatom) THEN
                     work(1:ncoa, 1:nsgfb(jset)) = MATMUL(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1), &
                                                          p_block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1))
                  ELSE
                     work(1:ncoa, 1:nsgfb(jset)) = MATMUL(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1), &
                                                       TRANSPOSE(p_block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1)))
                  END IF
                  pab(1:ncoa, 1:ncob, nij) = MATMUL(work(1:ncoa, 1:nsgfb(jset)), &
                                                    TRANSPOSE(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1)))
               END DO
            END DO
         END IF

         ! loop over all kinds for pseudopotential  atoms
         hab = 0._dp
         DO kkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(kkind), all_potential=all_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(all_potential)) THEN
               CALL get_potential(potential=all_potential, &
                                  alpha_core_charge=alpha_c, zeff=zeta_c, &
                                  ccore_charge=core_charge, core_charge_radius=core_radius)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, &
                                  alpha_core_charge=alpha_c, zeff=zeta_c, &
                                  ccore_charge=core_charge, core_charge_radius=core_radius)
            ELSE
               CYCLE
            END IF

            CALL nl_set_sub_iterator(ap_iterator, ikind, kkind, iatom, mepos=mepos)

            DO WHILE (nl_sub_iterate(ap_iterator, mepos=mepos) == 0)
               CALL get_iterator_info(ap_iterator, jatom=katom, r=rac, mepos=mepos)

               dac = SQRT(SUM(rac*rac))
               rbc(:) = rac(:) - rab(:)
               dbc = SQRT(SUM(rbc*rbc))
               IF ((MAXVAL(set_radius_a(:)) + core_radius < dac) .OR. &
                   (MAXVAL(set_radius_b(:)) + core_radius < dbc)) THEN
                  CYCLE
               END IF

               DO iset = 1, nseta
                  IF (set_radius_a(iset) + core_radius < dac) CYCLE
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1, iset)
                  DO jset = 1, nsetb
                     IF (set_radius_b(jset) + core_radius < dbc) CYCLE
                     ncob = npgfb(jset)*ncoset(lb_max(jset))
                     sgfb = first_sgfb(1, jset)
                     IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                     rab2 = dab*dab
                     rac2 = dac*dac
                     rbc2 = dbc*dbc
                     nij = jset + (iset - 1)*maxnset
                     ! *** Calculate the GTH pseudo potential forces ***
                     IF (calculate_forces) THEN
                        na_plus = npgfa(iset)*ncoset(la_max(iset) + nder)
                        nb_plus = npgfb(jset)*ncoset(lb_max(jset))
                        ALLOCATE (habd(na_plus, nb_plus))
                        habd = 0._dp
                        CALL verfc( &
                           la_max(iset) + nder, npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                           alpha_c, core_radius, zeta_c, core_charge, &
                           rab, rab2, rac, rac2, rbc2, hab(:, :, nij), verf, vnuc, ff(0:), &
                           nder, habd)

                        ! *** The derivatives w.r.t. atomic center c are    ***
                        ! *** calculated using the translational invariance ***
                        ! *** of the first derivatives                      ***
                        CALL verfc_force(habd, pab(:, :, nij), force_a, force_b, nder, &
                                         la_max(iset), la_min(iset), npgfa(iset), zeta(:, iset), &
                                         lb_max(jset), lb_min(jset), npgfb(jset), zetb(:, jset), rab)

                        DEALLOCATE (habd)

                        force_thread(1, iatom) = force_thread(1, iatom) + f0*force_a(1)
                        force_thread(2, iatom) = force_thread(2, iatom) + f0*force_a(2)
                        force_thread(3, iatom) = force_thread(3, iatom) + f0*force_a(3)

                        force_thread(1, jatom) = force_thread(1, jatom) + f0*force_b(1)
                        force_thread(2, jatom) = force_thread(2, jatom) + f0*force_b(2)
                        force_thread(3, jatom) = force_thread(3, jatom) + f0*force_b(3)

                        force_thread(1, katom) = force_thread(1, katom) - f0*force_a(1) - f0*force_b(1)
                        force_thread(2, katom) = force_thread(2, katom) - f0*force_a(2) - f0*force_b(2)
                        force_thread(3, katom) = force_thread(3, katom) - f0*force_a(3) - f0*force_b(3)

                        IF (use_virial) THEN
                           CALL virial_pair_force(pv_thread, f0, force_a, rac)
                           CALL virial_pair_force(pv_thread, f0, force_b, rbc)
                        END IF
                     ELSE
                        CALL verfc( &
                           la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                           alpha_c, core_radius, zeta_c, core_charge, &
                           rab, rab2, rac, rac2, rbc2, hab(:, :, nij), verf, vnuc, ff(0:))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         ! *** Contract nuclear attraction integrals
         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)
               nij = jset + (iset - 1)*maxnset
!$             hash2 = MOD((iset - 1)*nsetb + jset, nlock) + 1
!$             hash = MOD(hash1 + hash2, nlock) + 1
               work(1:ncoa, 1:nsgfb(jset)) = MATMUL(hab(1:ncoa, 1:ncob, nij), &
                                                    sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1))
!$             CALL omp_set_lock(locks(hash))
               IF (iatom <= jatom) THEN
                  h_block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) = &
                     h_block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) + &
                     MATMUL(TRANSPOSE(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1)), work(1:ncoa, 1:nsgfb(jset)))
               ELSE
                  h_block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) = &
                     h_block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) + &
                     MATMUL(TRANSPOSE(work(1:ncoa, 1:nsgfb(jset))), sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1))
               END IF
!$             CALL omp_unset_lock(locks(hash))
            END DO
         END DO

      END DO

      DEALLOCATE (hab, work, verf, vnuc, ff)
      IF (calculate_forces) THEN
         DEALLOCATE (pab)
      END IF

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL neighbor_list_iterator_release(ap_iterator)

      DEALLOCATE (basis_set_list)

      IF (calculate_forces) THEN
         ! *** If LSD, then recover alpha density and beta density     ***
         ! *** from the total density (1) and the spin density (2)     ***
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimages
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=0.5_dp, beta_scalar=0.5_dp)
               CALL dbcsr_add(matrix_p(2, img)%matrix, matrix_p(1, img)%matrix, &
                              alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
      END IF

      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natom), kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, &
                                  kind_of=kind_of)
!$OMP DO
         DO iatom = 1, natom
            atom_a = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) + force_thread(:, iatom)
         END DO
!$OMP END DO
         DEALLOCATE (atom_of_kind, kind_of)
      END IF
      DEALLOCATE (force_thread)

      IF (calculate_forces .AND. use_virial) THEN
         virial%pv_virial = virial%pv_virial + pv_thread
         virial%pv_ppl = virial%pv_virial - pv_loc
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_core_ae

! **************************************************************************************************
!> \brief ...
!> \param habd ...
!> \param pab ...
!> \param fa ...
!> \param fb ...
!> \param nder ...
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param rab ...
! **************************************************************************************************
   SUBROUTINE verfc_force(habd, pab, fa, fb, nder, la_max, la_min, npgfa, zeta, lb_max, lb_min, npgfb, zetb, rab)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: habd, pab
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: fa, fb
      INTEGER, INTENT(IN)                                :: nder, la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'verfc_force', routineP = moduleN//':'//routineN

      INTEGER                                            :: ic_a, ic_b, icam1, icam2, icam3, icap1, &
                                                            icap2, icap3, icax, icbm1, icbm2, &
                                                            icbm3, icbx, icoa, icob, ipgfa, ipgfb, &
                                                            na, nap, nb
      INTEGER, DIMENSION(3)                              :: la, lb
      REAL(KIND=dp)                                      :: zax2, zbx2

      fa = 0.0_dp
      fb = 0.0_dp

      na = ncoset(la_max)
      nap = ncoset(la_max + nder)
      nb = ncoset(lb_max)
      DO ipgfa = 1, npgfa
         zax2 = zeta(ipgfa)*2.0_dp
         DO ipgfb = 1, npgfb
            zbx2 = zetb(ipgfb)*2.0_dp
            DO ic_a = ncoset(la_min - 1) + 1, ncoset(la_max)
               la(1:3) = indco(1:3, ic_a)
               icap1 = coset(la(1) + 1, la(2), la(3))
               icap2 = coset(la(1), la(2) + 1, la(3))
               icap3 = coset(la(1), la(2), la(3) + 1)
               icam1 = coset(la(1) - 1, la(2), la(3))
               icam2 = coset(la(1), la(2) - 1, la(3))
               icam3 = coset(la(1), la(2), la(3) - 1)
               icoa = ic_a + (ipgfa - 1)*na
               icax = (ipgfa - 1)*nap

               DO ic_b = ncoset(lb_min - 1) + 1, ncoset(lb_max)
                  lb(1:3) = indco(1:3, ic_b)
                  icbm1 = coset(lb(1) - 1, lb(2), lb(3))
                  icbm2 = coset(lb(1), lb(2) - 1, lb(3))
                  icbm3 = coset(lb(1), lb(2), lb(3) - 1)
                  icob = ic_b + (ipgfb - 1)*nb
                  icbx = (ipgfb - 1)*nb

                  fa(1) = fa(1) - pab(icoa, icob)*(-zax2*habd(icap1 + icax, icob) + &
                                                   REAL(la(1), KIND=dp)*habd(icam1 + icax, icob))
                  fa(2) = fa(2) - pab(icoa, icob)*(-zax2*habd(icap2 + icax, icob) + &
                                                   REAL(la(2), KIND=dp)*habd(icam2 + icax, icob))
                  fa(3) = fa(3) - pab(icoa, icob)*(-zax2*habd(icap3 + icax, icob) + &
                                                   REAL(la(3), KIND=dp)*habd(icam3 + icax, icob))

                  fb(1) = fb(1) - pab(icoa, icob)*( &
                          -zbx2*(habd(icap1 + icax, icob) - rab(1)*habd(ic_a + icax, icob)) + &
                          REAL(lb(1), KIND=dp)*habd(ic_a + icax, icbm1 + icbx))
                  fb(2) = fb(2) - pab(icoa, icob)*( &
                          -zbx2*(habd(icap2 + icax, icob) - rab(2)*habd(ic_a + icax, icob)) + &
                          REAL(lb(2), KIND=dp)*habd(ic_a + icax, icbm2 + icbx))
                  fb(3) = fb(3) - pab(icoa, icob)*( &
                          -zbx2*(habd(icap3 + icax, icob) - rab(3)*habd(ic_a + icax, icob)) + &
                          REAL(lb(3), KIND=dp)*habd(ic_a + icax, icbm3 + icbx))

               END DO ! ic_b
            END DO ! ic_a
         END DO ! ipgfb
      END DO ! ipgfa

   END SUBROUTINE verfc_force

! **************************************************************************************************

END MODULE core_ae
