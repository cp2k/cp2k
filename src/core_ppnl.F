!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Calculation of the non-local pseudopotential contribution to the core Hamiltonian
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par History
!>      - refactered from qs_core_hamiltian [Joost VandeVondele, 2008-11-01]
!>      - full rewrite [jhu, 2009-01-23]
! **************************************************************************************************
MODULE core_ppnl
   USE ai_overlap,                      ONLY: overlap
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              clist_type,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppnl'

   PUBLIC :: build_core_ppnl

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
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param nimages ...
!> \param cell_to_index ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_core_ppnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
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
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      INTEGER, INTENT(IN)                                :: nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_core_ppnl'

      INTEGER :: atom_a, first_col, handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, img, &
         irow, iset, j, jatom, jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, &
         ldsab, lppnl, maxco, maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, na, natom, nb, &
         ncoa, ncoc, nkind, nlist, nneighbor, nnl, np, nppnl, nprjc, nseta, nsgfa, prjc, sgfa, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: dogth, dokp, found, ppnl_present
      LOGICAL, DIMENSION(0:9)                            :: is_nonlocal
      REAL(KIND=dp)                                      :: dac, f0, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: radp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sab, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: fa, fb, rab, rac, rbc
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_thread
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(clist_type), POINTER                          :: clist
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint, c_nl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, h_block, h_nl, p_block, rpgfa, &
                                                            sphi_a, vprj_ppnl, zeta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_nl, alpha_ppnl, hprj, set_radius_a
      REAL(KIND=dp), DIMENSION(3, SIZE(particle_set))    :: force_thread
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      MARK_USED(int_8)

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      ENDIF

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

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

         maxder = ncoset(nder)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl, &
                              basis_type=basis_type)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl + nder + 1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl + nder + 1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=basis_type)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, sgp_potential=sgp_potential)
            NULLIFY (gpotential(ikind)%gth_potential)
            NULLIFY (spotential(ikind)%sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
            END IF
         END DO

         !allocate sap int
         DO slot = 1, sap_ppnl(1)%nl_size

            ikind = sap_ppnl(1)%nlist_task(slot)%ikind
            kkind = sap_ppnl(1)%nlist_task(slot)%jkind
            iatom = sap_ppnl(1)%nlist_task(slot)%iatom
            katom = sap_ppnl(1)%nlist_task(slot)%jatom
            nlist = sap_ppnl(1)%nlist_task(slot)%nlist
            ilist = sap_ppnl(1)%nlist_task(slot)%ilist
            nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode

            iac = ikind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(gpotential(kkind)%gth_potential) .AND. &
                .NOT. ASSOCIATED(spotential(kkind)%sgp_potential)) CYCLE
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
         END DO

         !calculate the overlap integrals <a|p>
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, gpotential, spotential, maxder, ncoset, &
!$OMP          sap_ppnl, sap_int, nkind, ldsab, ldai, nder, nco ) &
!$OMP PRIVATE (ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          slot, sphi_a, zeta, cprj, hprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, first_col, vprj_ppnl, i, j, l, dogth, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl, &
!$OMP          na, nb, np, nnl, is_nonlocal, a_nl, h_nl, c_nl, radp)

         ALLOCATE (sab(ldsab, ldsab*maxder), work(ldsab, ldsab*maxder))
         sab = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, ncoset(nder + 1)))
         ai_work = 0.0_dp

!$OMP DO SCHEDULE(GUIDED)
         DO slot = 1, sap_ppnl(1)%nl_size

            ikind = sap_ppnl(1)%nlist_task(slot)%ikind
            kkind = sap_ppnl(1)%nlist_task(slot)%jkind
            iatom = sap_ppnl(1)%nlist_task(slot)%iatom
            katom = sap_ppnl(1)%nlist_task(slot)%jatom
            nlist = sap_ppnl(1)%nlist_task(slot)%nlist
            ilist = sap_ppnl(1)%nlist_task(slot)%ilist
            nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode
            jneighbor = sap_ppnl(1)%nlist_task(slot)%inode
            cell_c(:) = sap_ppnl(1)%nlist_task(slot)%cell(:)
            rac(1:3) = sap_ppnl(1)%nlist_task(slot)%r(1:3)

            iac = ikind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            sphi_a => basis_set(ikind)%gto_basis_set%sphi
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of PP projectors
            IF (ASSOCIATED(gpotential(kkind)%gth_potential)) THEN
               ! GTH potential
               dogth = .TRUE.
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSE IF (ASSOCIATED(spotential(kkind)%sgp_potential)) THEN
               ! SGP potential
               dogth = .FALSE.
               nprjc = spotential(kkind)%sgp_potential%nppnl
               IF (nprjc == 0) CYCLE
               nnl = spotential(kkind)%sgp_potential%n_nonlocal
               lppnl = spotential(kkind)%sgp_potential%lmax
               is_nonlocal = .FALSE.
               is_nonlocal(0:lppnl) = spotential(kkind)%sgp_potential%is_nonlocal(0:lppnl)
               a_nl => spotential(kkind)%sgp_potential%a_nonlocal
               h_nl => spotential(kkind)%sgp_potential%h_nonlocal
               c_nl => spotential(kkind)%sgp_potential%c_nonlocal
               ppnl_radius = spotential(kkind)%sgp_potential%ppnl_radius
               ALLOCATE (radp(nnl))
               radp(:) = ppnl_radius
               cprj => spotential(kkind)%sgp_potential%cprj_ppnl
               hprj => spotential(kkind)%sgp_potential%vprj_ppnl
               nppnl = SIZE(cprj, 2)
            ELSE
               CYCLE
            END IF

            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nppnl, maxder), &
                      clist%achint(nsgfa, nppnl, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)
            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               IF (dogth) THEN
                  ! GTH potential
                  prjc = 1
                  work = 0._dp
                  DO l = 0, lppnl
                     nprjc = nprj_ppnl(l)*nco(l)
                     IF (nprjc == 0) CYCLE
                     rprjc(1) = ppnl_radius
                     IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
                     lc_max = l + 2*(nprj_ppnl(l) - 1)
                     lc_min = l
                     zetc(1) = alpha_ppnl(l)
                     ncoc = ncoset(lc_max)

                     ! *** Calculate the primitive overlap integrals ***
                     CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab, nder, .TRUE., ai_work, ldai)
                     ! *** Transformation step projector functions (cartesian->spherical) ***
                     na = ncoa
                     nb = nprjc
                     np = ncoc
                     DO i = 1, maxder
                        first_col = (i - 1)*ldsab
!                       CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, first_col + 1), SIZE(sab, 1), &
!                                  cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, first_col + prjc), ldsab)
                        work(1:na, first_col + prjc:first_col + prjc + nb - 1) = &
                           MATMUL(sab(1:na, first_col + 1:first_col + np), cprj(1:np, prjc:prjc + nb - 1))
                     END DO
                     prjc = prjc + nprjc
                  END DO
                  na = nsgf_seta(iset)
                  nb = nppnl
                  np = ncoa
                  DO i = 1, maxder
                     first_col = (i - 1)*ldsab + 1
                     ! *** Contraction step (basis functions) ***
!                    CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
!                               work(1, first_col), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                     clist%acint(sgfa:sgfa + na - 1, 1:nb, i) = &
                        MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work(1:np, first_col:first_col + nb - 1))
                     ! *** Multiply with interaction matrix(h) ***
!                    CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
!                               vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
                     clist%achint(sgfa:sgfa + na - 1, 1:nb, i) = &
                        MATMUL(clist%acint(sgfa:sgfa + na - 1, 1:nb, i), vprj_ppnl(1:nb, 1:nb))
                  END DO
               ELSE
                  ! SGP potential
                  ! *** Calculate the primitive overlap integrals ***
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lppnl, 0, nnl, radp, a_nl, rac, dac, sab, nder, .TRUE., ai_work, ldai)
                  na = nsgf_seta(iset)
                  nb = nppnl
                  np = ncoa
                  DO i = 1, maxder
                     first_col = (i - 1)*ldsab + 1
                     ! *** Transformation step projector functions (cartesian->spherical) ***
!                    CALL dgemm("N", "N", ncoa, nppnl, nprjc, 1.0_dp, sab(1, first_col), ldsab, &
!                               cprj(1, 1), SIZE(cprj, 1), 0.0_dp, work(1, 1), ldsab)
                     work(1:np, 1:nb) = MATMUL(sab(1:np, first_col:first_col + nprjc - 1), cprj(1:nprjc, 1:nb))
                     ! *** Contraction step (basis functions) ***
!                    CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
!                               work(1, 1), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                     clist%acint(sgfa:sgfa + na - 1, 1:nb, i) = &
                        MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work(1:np, 1:nb))
                     ! *** Multiply with interaction matrix(h) ***
                     ncoc = sgfa + nsgf_seta(iset) - 1
                     DO j = 1, nppnl
                        clist%achint(sgfa:ncoc, j, i) = clist%acint(sgfa:ncoc, j, i)*hprj(j)
                     END DO
                  END DO
               END IF
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
            IF (.NOT. dogth) DEALLOCATE (radp)
         END DO

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements

         force_thread = 0.0_dp
         pv_thread = 0.0_dp

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (dokp, basis_set, matrix_h, cell_to_index,&
!$OMP          sab_orb, matrix_p, sap_int, nkind, eps_ppnl, force, &
!$OMP          locks, virial, use_virial, calculate_forces) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          slot, iab, atom_a, f0, irow, icol, h_block, &
!$OMP          found,p_block, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, j, fa, fb, rbc, rac, &
!$OMP          kkind, kac, kbc, i, img, hash, iatom8, natom) &
!$OMP REDUCTION (+ : pv_thread, force_thread )

!$OMP SINGLE
!$       ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$       DO lock_num = 1, nlock
!$          call omp_init_lock(locks(lock_num))
!$       END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)
         DO slot = 1, sab_orb(1)%nl_size

            ikind = sab_orb(1)%nlist_task(slot)%ikind
            jkind = sab_orb(1)%nlist_task(slot)%jkind
            iatom = sab_orb(1)%nlist_task(slot)%iatom
            jatom = sab_orb(1)%nlist_task(slot)%jatom
            cell_b(:) = sab_orb(1)%nlist_task(slot)%cell(:)
            rab(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)

            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE

            iab = ikind + nkind*(jkind - 1)

            ! *** Use the symmetry of the first derivatives ***
            IF (iatom == jatom) THEN
               f0 = 1.0_dp
            ELSE
               f0 = 2.0_dp
            END IF

            IF (dokp) THEN
               img = cell_to_index(cell_b(1), cell_b(2), cell_b(3))
            ELSE
               img = 1
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
            IF (calculate_forces) THEN
               NULLIFY (p_block)
               CALL dbcsr_get_block_p(matrix_p(1, img)%matrix, irow, icol, p_block, found)
            END IF

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(h_block)) THEN
               DO kkind = 1, nkind
                  iac = ikind + nkind*(kkind - 1)
                  ibc = jkind + nkind*(kkind - 1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$                         iatom8 = INT(iatom - 1, int_8)*INT(natom, int_8) + INT(jatom, int_8)
!$                         hash = INT(MOD(iatom8, INT(nlock, int_8)) + 1)
!$                         CALL omp_set_lock(locks(hash))
                           IF (iatom <= jatom) THEN
                              h_block(1:na, 1:nb) = h_block(1:na, 1:nb) + &
                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 1)))
                           ELSE
                              h_block(1:nb, 1:na) = h_block(1:nb, 1:na) + &
                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 1)))
                           END IF
!$                         CALL omp_unset_lock(locks(hash))
                           IF (calculate_forces) THEN
                              IF (ASSOCIATED(p_block)) THEN
                                 katom = alist_ac%clist(kac)%catom
                                 DO i = 1, 3
                                    j = i + 1
                                    IF (iatom <= jatom) THEN
                                       fa(i) = SUM(p_block(1:na, 1:nb)* &
                                                   MATMUL(acint(1:na, 1:np, j), TRANSPOSE(bchint(1:nb, 1:np, 1))))
                                       fb(i) = SUM(p_block(1:na, 1:nb)* &
                                                   MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, j))))
                                    ELSE
                                       fa(i) = SUM(p_block(1:nb, 1:na)* &
                                                   MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, j))))
                                       fb(i) = SUM(p_block(1:nb, 1:na)* &
                                                   MATMUL(bcint(1:nb, 1:np, j), TRANSPOSE(achint(1:na, 1:np, 1))))
                                    END IF
                                    force_thread(i, iatom) = force_thread(i, iatom) + f0*fa(i)
                                    force_thread(i, katom) = force_thread(i, katom) - f0*fa(i)
                                    force_thread(i, jatom) = force_thread(i, jatom) + f0*fb(i)
                                    force_thread(i, katom) = force_thread(i, katom) - f0*fb(i)
                                 END DO

                                 IF (use_virial) THEN
                                    rac = alist_ac%clist(kac)%rac
                                    rbc = alist_bc%clist(kbc)%rac
                                    CALL virial_pair_force(pv_thread, f0, fa, rac)
                                    CALL virial_pair_force(pv_thread, f0, fb, rbc)
                                 END IF
                              ENDIF
                           END IF
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO

!$OMP DO
!$       DO lock_num = 1, nlock
!$          call omp_destroy_lock(locks(lock_num))
!$       END DO
!$OMP END DO

!$OMP SINGLE
!$       DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

         CALL release_sap_int(sap_int)

         DEALLOCATE (basis_set, gpotential, spotential)
         IF (calculate_forces) THEN
            ALLOCATE (atom_of_kind(natom), kind_of(natom))
            CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, &
                                     kind_of=kind_of)
!$OMP DO
            DO iatom = 1, natom
               atom_a = atom_of_kind(iatom)
               ikind = kind_of(iatom)
               force(ikind)%gth_ppnl(:, atom_a) = force(ikind)%gth_ppnl(:, atom_a) + force_thread(:, iatom)
            END DO
!$OMP END DO
            DEALLOCATE (atom_of_kind, kind_of)
         END IF

         IF (calculate_forces .AND. use_virial) THEN
            virial%pv_ppnl = virial%pv_ppnl + pv_thread
            virial%pv_virial = virial%pv_virial + pv_thread
         END IF

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

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_core_ppnl

! **************************************************************************************************

END MODULE core_ppnl
