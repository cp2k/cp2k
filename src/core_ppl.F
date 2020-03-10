!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Calculation of the local pseudopotential contribution to the core Hamiltonian
!>         <a|V(local)|b> = <a|Sum e^a*rc**2|b>
!> \par History
!>      - core_ppnl refactored from qs_core_hamiltonian [Joost VandeVondele, 2008-11-01]
!>      - adapted for PPL [jhu, 2009-02-23]
!>      - OpenMP added [Iain Bethune, Fiona Reid, 2013-11-13]
!>      - Bug fix: correct orbital pointer range [07.2014,JGH]
!>      - k-point aware [07.2015,JGH]
! **************************************************************************************************
MODULE core_ppl

   USE ai_overlap_ppl,                  ONLY: ppl_integral,&
                                              ppl_integral_ri
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE lri_environment_types,           ONLY: lri_kind_type
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppl'

   PUBLIC :: build_core_ppl, build_core_ppl_ri

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
!> \param sac_ppl ...
!> \param nimages ...
!> \param cell_to_index ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_core_ppl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, basis_type)

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
         POINTER                                         :: sab_orb, sac_ppl
      INTEGER, INTENT(IN)                                :: nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppl', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nexp_max = 30

      INTEGER :: atom_a, handle, i, iatom, icol, ikind, img, irow, iset, jatom, jkind, jset, &
         katom, kkind, ldai, ldsab, maxco, maxder, maxl, maxlgto, maxlppl, maxnset, maxsgf, mepos, &
         n_local, natom, ncoa, ncob, nexp_lpot, nexp_ppl, nkind, nloc, nseta, nsetb, nthread, &
         sgfa, sgfb, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(1:10)                           :: nrloc
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, &
                                                            nct_lpot, npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(nexp_max)                       :: nct_ppl
      LOGICAL                                            :: dokp, ecp_local, found, lpotextended
      REAL(KIND=dp)                                      :: alpha, dab, dac, dbc, f0, ppl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ppl_fwork, ppl_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: hab, pab
      REAL(KIND=dp), DIMENSION(1:10)                     :: aloc, bloc
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, rab, rac, rbc
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_loc, pv_thread
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ap_iterator
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      REAL(KIND=dp), DIMENSION(nexp_max)                 :: alpha_ppl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_lpot, h_block, p_block, rpgfa, &
                                                            rpgfb, sphi_a, sphi_b, zeta, zetb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_local, alpha_lpot, c_local, cexp_ppl, &
                                                            set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(4, nexp_max)              :: cval_ppl
      REAL(KIND=dp), DIMENSION(3, SIZE(particle_set))    :: force_thread
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash, hash1, hash2
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER, PARAMETER                                 :: nlock = 501

      MARK_USED(int_8)

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
      force_thread = 0.0_dp

      maxder = ncoset(nder)

      CALL get_qs_kind_set(qs_kind_set, maxco=maxco, maxlgto=maxlgto, &
                           maxsgf=maxsgf, maxnset=maxnset, maxlppl=maxlppl, &
                           basis_type=basis_type)

      maxl = MAX(maxlgto, maxlppl)
      CALL init_orbital_pointers(2*maxl + 2*nder + 1)

      ldsab = MAX(maxco, ncoset(maxlppl), maxsgf, maxlppl)
      ldai = ncoset(maxl + nder + 1)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF
      pv_thread = 0.0_dp

      nthread = 1
!$    nthread = omp_get_max_threads()

      ! iterator for basis/potential list
      CALL neighbor_list_iterator_create(ap_iterator, sac_ppl, search=.TRUE., nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (ap_iterator, basis_set_list, calculate_forces, use_virial, &
!$OMP          matrix_h, matrix_p, atomic_kind_set, qs_kind_set, particle_set, &
!$OMP          sab_orb, sac_ppl, nthread, ncoset, nkind, cell_to_index, &
!$OMP          ldsab,  maxnset, maxder, &
!$OMP          maxlgto, nder, maxco, dokp, locks, natom) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, rab, basis_set_a, basis_set_b, &
!$OMP          first_sgfa, la_max, la_min, npgfa, nsgfa, sphi_a, &
!$OMP          zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, rpgfb, set_radius_b, sphi_b, &
!$OMP          zetb, dab, irow, icol, h_block, found, iset, ncoa, &
!$OMP          sgfa, jset, ncob, sgfb, nsgfb, p_block, work, pab, hab, kkind, nseta, &
!$OMP          gth_potential, sgp_potential, alpha, cexp_ppl, lpotextended, &
!$OMP          ppl_radius, nexp_lpot, nexp_ppl, alpha_ppl, alpha_lpot, nct_ppl, &
!$OMP          nct_lpot, cval_ppl, cval_lpot, rac, dac, rbc, dbc, &
!$OMP          set_radius_a,  rpgfa, force_a, force_b, ppl_fwork, mepos, &
!$OMP          slot, f0, katom, ppl_work, cellind, img, ecp_local, &
!$OMP          nloc, nrloc, aloc, bloc, n_local, a_local, c_local, ldai, hash, hash1, hash2, iatom8) &
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

      ALLOCATE (hab(ldsab, ldsab, maxnset, maxnset), work(ldsab, ldsab*maxder))
      ldai = ncoset(2*maxlgto + 2*nder)
      ALLOCATE (ppl_work(ldai, ldai, MAX(maxder, 2*maxlgto + 2*nder + 1)))
      IF (calculate_forces) THEN
         ALLOCATE (pab(maxco, maxco, maxnset, maxnset))
         ldai = ncoset(maxlgto)
         ALLOCATE (ppl_fwork(ldai, ldai, maxder))
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

!$       iatom8 = INT(iatom - 1, int_8)*INT(natom, int_8) + INT(jatom, int_8)
!$       hash1 = INT(MOD(iatom8, INT(nlock, int_8)) + 1)

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
         CALL dbcsr_get_block_p(matrix_h(1, img)%matrix, irow, icol, h_block, found)
         CPASSERT(found)
         IF (calculate_forces) THEN
            NULLIFY (p_block)
            CALL dbcsr_get_block_p(matrix_p(1, img)%matrix, irow, icol, p_block, found)
            IF (ASSOCIATED(p_block)) THEN
               DO iset = 1, nseta
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1, iset)
                  DO jset = 1, nsetb
                     ncob = npgfb(jset)*ncoset(lb_max(jset))
                     sgfb = first_sgfb(1, jset)

                     ! *** Decontract density matrix block ***
                     IF (iatom <= jatom) THEN
                        work(1:ncoa, 1:nsgfb(jset)) = MATMUL(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1), &
                                                             p_block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1))
                     ELSE
                        work(1:ncoa, 1:nsgfb(jset)) = MATMUL(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1), &
                                                       TRANSPOSE(p_block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1)))
                     END IF

                     pab(1:ncoa, 1:ncob, iset, jset) = MATMUL(work(1:ncoa, 1:nsgfb(jset)), &
                                                              TRANSPOSE(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1)))
                  END DO
               END DO
            ENDIF
         END IF

         hab = 0._dp

         ! loop over all kinds for pseudopotential atoms
         DO kkind = 1, nkind

            CALL get_qs_kind(qs_kind_set(kkind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               CALL get_potential(potential=gth_potential, &
                                  alpha_ppl=alpha, cexp_ppl=cexp_ppl, &
                                  lpot_present=lpotextended, ppl_radius=ppl_radius)
               nexp_ppl = 1
               alpha_ppl(1) = alpha
               nct_ppl(1) = SIZE(cexp_ppl)
               cval_ppl(1:nct_ppl(1), 1) = cexp_ppl(1:nct_ppl(1))
               IF (lpotextended) THEN
                  CALL get_potential(potential=gth_potential, &
                                     nexp_lpot=nexp_lpot, alpha_lpot=alpha_lpot, nct_lpot=nct_lpot, &
                                     cval_lpot=cval_lpot)
                  CPASSERT(nexp_lpot < nexp_max)
                  nexp_ppl = nexp_lpot + 1
                  alpha_ppl(2:nexp_lpot + 1) = alpha_lpot(1:nexp_lpot)
                  nct_ppl(2:nexp_lpot + 1) = nct_lpot(1:nexp_lpot)
                  DO i = 1, nexp_lpot
                     cval_ppl(1:nct_lpot(i), i + 1) = cval_lpot(1:nct_lpot(i), i)
                  END DO
               END IF
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, ppl_radius=ppl_radius)
               IF (ecp_local) THEN
                  CALL get_potential(potential=sgp_potential, nloc=nloc, nrloc=nrloc, aloc=aloc, bloc=bloc)
                  IF (SUM(ABS(aloc(1:nloc))) < 1.0e-12_dp) CYCLE
                  nexp_ppl = nloc
                  CPASSERT(nexp_ppl <= nexp_max)
                  nct_ppl(1:nloc) = nrloc(1:nloc) - 1
                  alpha_ppl(1:nloc) = bloc(1:nloc)
                  cval_ppl(1, 1:nloc) = aloc(1:nloc)
               ELSE
                  CALL get_potential(potential=sgp_potential, n_local=n_local, a_local=a_local, c_local=c_local)
                  nexp_ppl = n_local
                  CPASSERT(nexp_ppl <= nexp_max)
                  nct_ppl(1:n_local) = 1
                  alpha_ppl(1:n_local) = a_local(1:n_local)
                  cval_ppl(1, 1:n_local) = c_local(1:n_local)
               END IF
            ELSE
               CYCLE
            END IF

            CALL nl_set_sub_iterator(ap_iterator, ikind, kkind, iatom, mepos=mepos)

            DO WHILE (nl_sub_iterate(ap_iterator, mepos=mepos) == 0)

               CALL get_iterator_info(ap_iterator, mepos=mepos, jatom=katom, r=rac)

               dac = SQRT(SUM(rac*rac))
               rbc(:) = rac(:) - rab(:)
               dbc = SQRT(SUM(rbc*rbc))
               IF ((MAXVAL(set_radius_a(:)) + ppl_radius < dac) .OR. &
                   (MAXVAL(set_radius_b(:)) + ppl_radius < dbc)) THEN
                  CYCLE
               END IF

               DO iset = 1, nseta
                  IF (set_radius_a(iset) + ppl_radius < dac) CYCLE
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1, iset)
                  DO jset = 1, nsetb
                     IF (set_radius_b(jset) + ppl_radius < dbc) CYCLE
                     ncob = npgfb(jset)*ncoset(lb_max(jset))
                     sgfb = first_sgfb(1, jset)
                     IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                     ! *** Calculate the GTH pseudo potential forces ***
                     IF (calculate_forces) THEN

                        CALL ppl_integral( &
                           la_max(iset), la_min(iset), npgfa(iset), &
                           rpgfa(:, iset), zeta(:, iset), &
                           lb_max(jset), lb_min(jset), npgfb(jset), &
                           rpgfb(:, jset), zetb(:, jset), &
                           nexp_ppl, alpha_ppl, nct_ppl, cval_ppl, ppl_radius, &
                           rab, dab, rac, dac, rbc, dbc, &
                           hab(:, :, iset, jset), ppl_work, pab(:, :, iset, jset), &
                           force_a, force_b, ppl_fwork)
                        ! *** The derivatives w.r.t. atomic center c are    ***
                        ! *** calculated using the translational invariance ***
                        ! *** of the first derivatives                      ***

                        force_thread(1, iatom) = force_thread(1, iatom) + f0*force_a(1)
                        force_thread(2, iatom) = force_thread(2, iatom) + f0*force_a(2)
                        force_thread(3, iatom) = force_thread(3, iatom) + f0*force_a(3)
                        force_thread(1, katom) = force_thread(1, katom) - f0*force_a(1)
                        force_thread(2, katom) = force_thread(2, katom) - f0*force_a(2)
                        force_thread(3, katom) = force_thread(3, katom) - f0*force_a(3)

                        force_thread(1, jatom) = force_thread(1, jatom) + f0*force_b(1)
                        force_thread(2, jatom) = force_thread(2, jatom) + f0*force_b(2)
                        force_thread(3, jatom) = force_thread(3, jatom) + f0*force_b(3)
                        force_thread(1, katom) = force_thread(1, katom) - f0*force_b(1)
                        force_thread(2, katom) = force_thread(2, katom) - f0*force_b(2)
                        force_thread(3, katom) = force_thread(3, katom) - f0*force_b(3)

                        IF (use_virial) THEN
                           CALL virial_pair_force(pv_thread, f0, force_a, rac)
                           CALL virial_pair_force(pv_thread, f0, force_b, rbc)
                        END IF

                     ELSE
                        CALL ppl_integral( &
                           la_max(iset), la_min(iset), npgfa(iset), &
                           rpgfa(:, iset), zeta(:, iset), &
                           lb_max(jset), lb_min(jset), npgfb(jset), &
                           rpgfb(:, jset), zetb(:, jset), &
                           nexp_ppl, alpha_ppl, nct_ppl, cval_ppl, ppl_radius, &
                           rab, dab, rac, dac, rbc, dbc, hab(:, :, iset, jset), ppl_work)
                     END IF
                  END DO
               END DO
            END DO
         END DO

         ! *** Contract PPL integrals
         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

!$             hash2 = MOD((iset - 1)*nsetb + jset, nlock) + 1
!$             hash = MOD(hash1 + hash2, nlock) + 1

               work(1:ncoa, 1:nsgfb(jset)) = MATMUL(hab(1:ncoa, 1:ncob, iset, jset), &
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

      DEALLOCATE (hab, work, ppl_work)
      IF (calculate_forces) THEN
         DEALLOCATE (pab, ppl_fwork)
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
            force(ikind)%gth_ppl(:, atom_a) = force(ikind)%gth_ppl(:, atom_a) + force_thread(:, iatom)
         END DO
!$OMP END DO
         DEALLOCATE (atom_of_kind, kind_of)
      END IF

      IF (calculate_forces .AND. use_virial) THEN
         virial%pv_virial = virial%pv_virial + pv_thread
         virial%pv_ppl = virial%pv_virial - pv_loc
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_core_ppl

! **************************************************************************************************
!> \brief ...
!> \param lri_ppl_coef ...
!> \param force ...
!> \param virial ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param sac_ppl ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_core_ppl_ri(lri_ppl_coef, force, virial, calculate_forces, use_virial, &
                                qs_kind_set, atomic_kind_set, particle_set, sac_ppl, &
                                basis_type)

      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_ppl_coef
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL                                            :: use_virial
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sac_ppl
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppl_ri', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nexp_max = 30

      INTEGER :: atom_a, handle, i, iatom, ikind, iset, katom, kkind, maxco, maxsgf, n_local, &
         natom, ncoa, nexp_lpot, nexp_ppl, nfun, nkind, nloc, nseta, sgfa, sgfb, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(1:10)                           :: nrloc
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, nct_lpot, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      INTEGER, DIMENSION(nexp_max)                       :: nct_ppl
      LOGICAL                                            :: ecp_local, lpotextended
      REAL(KIND=dp)                                      :: alpha, dac, ppl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: va, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dva, dvas
      REAL(KIND=dp), DIMENSION(1:10)                     :: aloc, bloc
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, rac
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_thread
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      REAL(KIND=dp), DIMENSION(nexp_max)                 :: alpha_ppl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bcon, cval_lpot, rpgfa, sphi_a, zeta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_local, alpha_lpot, c_local, cexp_ppl, &
                                                            set_radius_a
      REAL(KIND=dp), DIMENSION(4, nexp_max)              :: cval_ppl
      REAL(KIND=dp), DIMENSION(3, SIZE(particle_set))    :: force_thread
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      ENDIF

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      force_thread = 0.0_dp
      pv_thread = 0.0_dp
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, &
                               kind_of=kind_of)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, basis_type=basis_type)
         IF (ASSOCIATED(basis_set)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      CALL get_qs_kind_set(qs_kind_set, maxco=maxco, maxsgf=maxsgf, basis_type=basis_type)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (maxco,maxsgf,basis_set_list,calculate_forces,lri_ppl_coef,qs_kind_set,&
!$OMP          locks,natom,use_virial,virial,ncoset,atom_of_kind,sac_ppl) &
!$OMP PRIVATE (ikind,kkind,iatom,katom,atom_a,rac,va,dva,dvas,basis_set,slot,&
!$OMP          first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,rpgfa,set_radius_a,&
!$OMP          sphi_a,zeta,gth_potential,sgp_potential,alpha,cexp_ppl,lpotextended,ppl_radius,&
!$OMP          nexp_ppl,alpha_ppl,nct_ppl,cval_ppl,nloc,n_local,nrloc,a_local,aloc,bloc,c_local,nfun,work,&
!$OMP          hash,dac,force_a,iset,sgfa,sgfb,ncoa,bcon,cval_lpot,nct_lpot,alpha_lpot,nexp_lpot,ecp_local) &
!$OMP REDUCTION (+ : pv_thread, force_thread )

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

      ALLOCATE (va(maxco), work(maxsgf))
      IF (calculate_forces) THEN
         ALLOCATE (dva(maxco, 3), dvas(maxco, 3))
      END IF

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sac_ppl(1)%nl_size

         ikind = sac_ppl(1)%nlist_task(slot)%ikind
         kkind = sac_ppl(1)%nlist_task(slot)%jkind
         iatom = sac_ppl(1)%nlist_task(slot)%iatom
         katom = sac_ppl(1)%nlist_task(slot)%jatom
         rac(1:3) = sac_ppl(1)%nlist_task(slot)%r(1:3)
         atom_a = atom_of_kind(iatom)

         basis_set => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set)) CYCLE

         ! basis ikind
         first_sgfa => basis_set%first_sgf
         la_max => basis_set%lmax
         la_min => basis_set%lmin
         npgfa => basis_set%npgf
         nseta = basis_set%nset
         nsgfa => basis_set%nsgf_set
         nfun = basis_set%nsgf
         rpgfa => basis_set%pgf_radius
         set_radius_a => basis_set%set_radius
         sphi_a => basis_set%sphi
         zeta => basis_set%zet

         CALL get_qs_kind(qs_kind_set(kkind), gth_potential=gth_potential, &
                          sgp_potential=sgp_potential)
         IF (ASSOCIATED(gth_potential)) THEN
            CALL get_potential(potential=gth_potential, &
                               alpha_ppl=alpha, cexp_ppl=cexp_ppl, &
                               lpot_present=lpotextended, ppl_radius=ppl_radius)
            nexp_ppl = 1
            alpha_ppl(1) = alpha
            nct_ppl(1) = SIZE(cexp_ppl)
            cval_ppl(1:nct_ppl(1), 1) = cexp_ppl(1:nct_ppl(1))
            IF (lpotextended) THEN
               CALL get_potential(potential=gth_potential, &
                                  nexp_lpot=nexp_lpot, alpha_lpot=alpha_lpot, nct_lpot=nct_lpot, cval_lpot=cval_lpot)
               CPASSERT(nexp_lpot < nexp_max)
               nexp_ppl = nexp_lpot + 1
               alpha_ppl(2:nexp_lpot + 1) = alpha_lpot(1:nexp_lpot)
               nct_ppl(2:nexp_lpot + 1) = nct_lpot(1:nexp_lpot)
               DO i = 1, nexp_lpot
                  cval_ppl(1:nct_lpot(i), i + 1) = cval_lpot(1:nct_lpot(i), i)
               END DO
            END IF
         ELSE IF (ASSOCIATED(sgp_potential)) THEN
            CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, ppl_radius=ppl_radius)
            IF (ecp_local) THEN
               CALL get_potential(potential=sgp_potential, nloc=nloc, nrloc=nrloc, aloc=aloc, bloc=bloc)
               IF (SUM(ABS(aloc(1:nloc))) < 1.0e-12_dp) CYCLE
               nexp_ppl = nloc
               CPASSERT(nexp_ppl <= nexp_max)
               nct_ppl(1:nloc) = nrloc(1:nloc) - 1
               alpha_ppl(1:nloc) = bloc(1:nloc)
               cval_ppl(1, 1:nloc) = aloc(1:nloc)
            ELSE
               CALL get_potential(potential=sgp_potential, n_local=n_local, a_local=a_local, c_local=c_local)
               nexp_ppl = n_local
               CPASSERT(nexp_ppl <= nexp_max)
               nct_ppl(1:n_local) = 1
               alpha_ppl(1:n_local) = a_local(1:n_local)
               cval_ppl(1, 1:n_local) = c_local(1:n_local)
            END IF
         ELSE
            CYCLE
         END IF

         dac = SQRT(SUM(rac*rac))
         IF ((MAXVAL(set_radius_a(:)) + ppl_radius < dac)) CYCLE
         IF (calculate_forces) force_a = 0.0_dp
         work(1:nfun) = 0.0_dp

         DO iset = 1, nseta
            IF (set_radius_a(iset) + ppl_radius < dac) CYCLE
            ! integrals
            IF (calculate_forces) THEN
               va = 0.0_dp
               dva = 0.0_dp
               CALL ppl_integral_ri( &
                  la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                  nexp_ppl, alpha_ppl, nct_ppl, cval_ppl, ppl_radius, &
                  -rac, dac, va, dva)
            ELSE
               va = 0.0_dp
               CALL ppl_integral_ri( &
                  la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                  nexp_ppl, alpha_ppl, nct_ppl, cval_ppl, ppl_radius, &
                  -rac, dac, va)
            END IF
            ! contraction
            sgfa = first_sgfa(1, iset)
            sgfb = sgfa + nsgfa(iset) - 1
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            bcon => sphi_a(1:ncoa, sgfa:sgfb)
            work(sgfa:sgfb) = MATMUL(TRANSPOSE(bcon), va(1:ncoa))
            IF (calculate_forces) THEN
               dvas(1:nsgfa(iset), 1:3) = MATMUL(TRANSPOSE(bcon), dva(1:ncoa, 1:3))
               force_a(1) = force_a(1) + SUM(lri_ppl_coef(ikind)%acoef(atom_a, sgfa:sgfb)*dvas(1:nsgfa(iset), 1))
               force_a(2) = force_a(2) + SUM(lri_ppl_coef(ikind)%acoef(atom_a, sgfa:sgfb)*dvas(1:nsgfa(iset), 2))
               force_a(3) = force_a(3) + SUM(lri_ppl_coef(ikind)%acoef(atom_a, sgfa:sgfb)*dvas(1:nsgfa(iset), 3))
            END IF
         END DO
!$       hash = MOD(iatom, nlock) + 1
!$       CALL omp_set_lock(locks(hash))
         lri_ppl_coef(ikind)%v_int(atom_a, 1:nfun) = lri_ppl_coef(ikind)%v_int(atom_a, 1:nfun) + work(1:nfun)
!$       CALL omp_unset_lock(locks(hash))
         IF (calculate_forces) THEN
            force_thread(1, iatom) = force_thread(1, iatom) + force_a(1)
            force_thread(2, iatom) = force_thread(2, iatom) + force_a(2)
            force_thread(3, iatom) = force_thread(3, iatom) + force_a(3)
            force_thread(1, katom) = force_thread(1, katom) - force_a(1)
            force_thread(2, katom) = force_thread(2, katom) - force_a(2)
            force_thread(3, katom) = force_thread(3, katom) - force_a(3)
            IF (use_virial) THEN
               CALL virial_pair_force(pv_thread, 1.0_dp, force_a, rac)
            END IF
         END IF
      END DO

      DEALLOCATE (va, work)
      IF (calculate_forces) THEN
         DEALLOCATE (dva, dvas)
      END IF

!$OMP END PARALLEL

      IF (calculate_forces) THEN
         DO iatom = 1, natom
            atom_a = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%gth_ppl(1, atom_a) = force(ikind)%gth_ppl(1, atom_a) + force_thread(1, iatom)
            force(ikind)%gth_ppl(2, atom_a) = force(ikind)%gth_ppl(2, atom_a) + force_thread(2, iatom)
            force(ikind)%gth_ppl(3, atom_a) = force(ikind)%gth_ppl(3, atom_a) + force_thread(3, iatom)
         END DO
      END IF
      DEALLOCATE (atom_of_kind, kind_of)

      IF (calculate_forces .AND. use_virial) THEN
         virial%pv_virial = virial%pv_virial + pv_thread
         virial%pv_ppl = pv_thread
      END IF

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_core_ppl_ri

! **************************************************************************************************

END MODULE core_ppl
