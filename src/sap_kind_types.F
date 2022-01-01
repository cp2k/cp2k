!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief General overlap type integrals containers
!> \par History
!>      - rewrite of PPNL and OCE integrals
! **************************************************************************************************
MODULE sap_kind_types

   USE ai_moments,                      ONLY: moment
   USE ai_overlap,                      ONLY: overlap
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE util,                            ONLY: locate,&
                                              sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sap_kind_types'

   TYPE clist_type
      INTEGER                                    :: catom, nsgf_cnt
      INTEGER, DIMENSION(:), POINTER             :: sgf_list
      INTEGER, DIMENSION(3)                      :: cell
      LOGICAL                                    :: sgf_soft_only
      REAL(KIND=dp)                              :: maxac, maxach
      REAL(KIND=dp), DIMENSION(3)                :: rac
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: acint
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: achint
   END TYPE clist_type

   TYPE alist_type
      INTEGER                                    :: aatom
      INTEGER                                    :: nclist
      TYPE(clist_type), DIMENSION(:), POINTER    :: clist
   END TYPE alist_type

   TYPE sap_int_type
      INTEGER                                    :: a_kind, p_kind
      INTEGER                                    :: nalist
      TYPE(alist_type), DIMENSION(:), POINTER    :: alist
      INTEGER, DIMENSION(:), POINTER             :: asort, aindex
   END TYPE sap_int_type

   PUBLIC :: sap_int_type, clist_type, alist_type, &
             release_sap_int, get_alist, alist_pre_align_blk, &
             alist_post_align_blk, sap_sort, build_sap_ints

CONTAINS

!==========================================================================================================

! **************************************************************************************************
!> \brief ...
!> \param sap_int ...
! **************************************************************************************************
   SUBROUTINE release_sap_int(sap_int)

      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

      INTEGER                                            :: i, j, k
      TYPE(clist_type), POINTER                          :: clist

      CPASSERT(ASSOCIATED(sap_int))

      DO i = 1, SIZE(sap_int)
         IF (ASSOCIATED(sap_int(i)%alist)) THEN
            DO j = 1, SIZE(sap_int(i)%alist)
               IF (ASSOCIATED(sap_int(i)%alist(j)%clist)) THEN
                  DO k = 1, SIZE(sap_int(i)%alist(j)%clist)
                     clist => sap_int(i)%alist(j)%clist(k)
                     IF (ASSOCIATED(clist%acint)) THEN
                        DEALLOCATE (clist%acint)
                     END IF
                     IF (ASSOCIATED(clist%sgf_list)) THEN
                        DEALLOCATE (clist%sgf_list)
                     END IF
                     IF (ASSOCIATED(clist%achint)) THEN
                        DEALLOCATE (clist%achint)
                     END IF
                  END DO
                  DEALLOCATE (sap_int(i)%alist(j)%clist)
               END IF
            END DO
            DEALLOCATE (sap_int(i)%alist)
         END IF
         IF (ASSOCIATED(sap_int(i)%asort)) THEN
            DEALLOCATE (sap_int(i)%asort)
         END IF
         IF (ASSOCIATED(sap_int(i)%aindex)) THEN
            DEALLOCATE (sap_int(i)%aindex)
         END IF
      END DO

      DEALLOCATE (sap_int)

   END SUBROUTINE release_sap_int

! **************************************************************************************************
!> \brief ...
!> \param sap_int ...
!> \param alist ...
!> \param atom ...
! **************************************************************************************************
   SUBROUTINE get_alist(sap_int, alist, atom)

      TYPE(sap_int_type), INTENT(IN)                     :: sap_int
      TYPE(alist_type), INTENT(OUT), POINTER             :: alist
      INTEGER, INTENT(IN)                                :: atom

      INTEGER                                            :: i

      NULLIFY (alist)
      i = locate(sap_int%asort, atom)
      IF (i > 0 .AND. i <= SIZE(sap_int%alist)) THEN
         i = sap_int%aindex(i)
         alist => sap_int%alist(i)
      ELSE IF (i == 0) THEN
         NULLIFY (alist)
      ELSE
         CPABORT("")
      END IF

   END SUBROUTINE get_alist

! **************************************************************************************************
!> \brief ...
!> \param blk_in ...
!> \param ldin ...
!> \param blk_out ...
!> \param ldout ...
!> \param ilist ...
!> \param in ...
!> \param jlist ...
!> \param jn ...
! **************************************************************************************************
   SUBROUTINE alist_pre_align_blk(blk_in, ldin, blk_out, ldout, ilist, in, jlist, jn)
      INTEGER, INTENT(IN)                                :: in, ilist(*), ldout
      REAL(dp), INTENT(INOUT)                            :: blk_out(ldout, *)
      INTEGER, INTENT(IN)                                :: ldin
      REAL(dp), INTENT(IN)                               :: blk_in(ldin, *)
      INTEGER, INTENT(IN)                                :: jlist(*), jn

      INTEGER                                            :: i, i0, i1, i2, i3, inn, inn1, j, j0

      inn = MOD(in, 4)
      inn1 = inn + 1
      DO j = 1, jn
         j0 = jlist(j)
         DO i = 1, inn
            i0 = ilist(i)
            blk_out(i, j) = blk_in(i0, j0)
         END DO
         DO i = inn1, in, 4
            i0 = ilist(i)
            i1 = ilist(i + 1)
            i2 = ilist(i + 2)
            i3 = ilist(i + 3)
            blk_out(i, j) = blk_in(i0, j0)
            blk_out(i + 1, j) = blk_in(i1, j0)
            blk_out(i + 2, j) = blk_in(i2, j0)
            blk_out(i + 3, j) = blk_in(i3, j0)
         END DO
      END DO
   END SUBROUTINE alist_pre_align_blk

! **************************************************************************************************
!> \brief ...
!> \param blk_in ...
!> \param ldin ...
!> \param blk_out ...
!> \param ldout ...
!> \param ilist ...
!> \param in ...
!> \param jlist ...
!> \param jn ...
! **************************************************************************************************
   SUBROUTINE alist_post_align_blk(blk_in, ldin, blk_out, ldout, ilist, in, jlist, jn)
      INTEGER, INTENT(IN)                                :: in, ilist(*), ldout
      REAL(dp), INTENT(INOUT)                            :: blk_out(ldout, *)
      INTEGER, INTENT(IN)                                :: ldin
      REAL(dp), INTENT(IN)                               :: blk_in(ldin, *)
      INTEGER, INTENT(IN)                                :: jlist(*), jn

      INTEGER                                            :: i, i0, i1, i2, i3, inn, inn1, j, j0

      inn = MOD(in, 4)
      inn1 = inn + 1
      DO j = 1, jn
         j0 = jlist(j)
         DO i = 1, inn
            i0 = ilist(i)
            blk_out(i0, j0) = blk_out(i0, j0) + blk_in(i, j)
         END DO
         DO i = inn1, in, 4
            i0 = ilist(i)
            i1 = ilist(i + 1)
            i2 = ilist(i + 2)
            i3 = ilist(i + 3)
            blk_out(i0, j0) = blk_out(i0, j0) + blk_in(i, j)
            blk_out(i1, j0) = blk_out(i1, j0) + blk_in(i + 1, j)
            blk_out(i2, j0) = blk_out(i2, j0) + blk_in(i + 2, j)
            blk_out(i3, j0) = blk_out(i3, j0) + blk_in(i + 3, j)
         END DO
      END DO
   END SUBROUTINE alist_post_align_blk

! **************************************************************************************************
!> \brief ...
!> \param sap_int ...
! **************************************************************************************************
   SUBROUTINE sap_sort(sap_int)
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

      INTEGER                                            :: iac, na

! *** Set up a sorting index

!$OMP PARALLEL DEFAULT(NONE) SHARED(sap_int) PRIVATE(iac,na)
!$OMP DO
      DO iac = 1, SIZE(sap_int)
         IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
         na = SIZE(sap_int(iac)%alist)
         ALLOCATE (sap_int(iac)%asort(na), sap_int(iac)%aindex(na))
         sap_int(iac)%asort(1:na) = sap_int(iac)%alist(1:na)%aatom
         CALL sort(sap_int(iac)%asort, na, sap_int(iac)%aindex)
      END DO
!$OMP END PARALLEL

   END SUBROUTINE sap_sort

!==========================================================================================================

! **************************************************************************************************
!> \brief Calculate overlap and optionally momenta <a|x^n|p> between GTOs and nl. pseudo potential projectors
!>        adapted from core_ppnl.F::build_core_ppnl
!> \param sap_int allocated in parent routine (nkind*nkind)
!> \param sap_ppnl ...
!> \param qs_kind_set ...
!> \param nder Either number of derivatives or order of moments
!> \param moment_mode if present and true, moments are calculated instead of derivatives
!> \param refpoint optionally the reference point for moment calculation
!> \param particle_set needed if refpoint is present
!> \param cell needed if refpoint is present
! **************************************************************************************************
   SUBROUTINE build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder, moment_mode, refpoint, particle_set, cell)
      TYPE(sap_int_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: sap_int
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sap_ppnl
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      INTEGER, INTENT(IN)                                :: nder
      LOGICAL, INTENT(IN), OPTIONAL                      :: moment_mode
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: refpoint
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         OPTIONAL, POINTER                               :: particle_set
      TYPE(cell_type), INTENT(IN), OPTIONAL, POINTER     :: cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_sap_ints'

      INTEGER :: first_col, handle, i, iac, iatom, ikind, ilist, iset, j, jneighbor, katom, kkind, &
         l, lc_max, lc_min, ldai, ldsab, lppnl, maxco, maxder, maxl, maxlgto, maxlppnl, maxppnl, &
         maxsgf, na, nb, ncoa, ncoc, nkind, nlist, nneighbor, nnl, np, nppnl, nprjc, nseta, nsgfa, &
         prjc, sgfa, slot
      INTEGER, DIMENSION(3)                              :: cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: dogth, my_momentmode, my_ref
      LOGICAL, DIMENSION(0:9)                            :: is_nonlocal
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: radp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: sab, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rac, raf, rc, rcf, rf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_nl, alpha_ppnl, hprj, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, h_nl, rpgfa, sphi_a, vprj_ppnl, &
                                                            zeta
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: c_nl
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      nkind = SIZE(qs_kind_set)
      maxder = ncoset(nder)

      ! determine whether moments or derivatives should be calculated
      my_momentmode = .FALSE.
      IF (PRESENT(moment_mode)) THEN
         my_momentmode = moment_mode
      END IF

      my_ref = .FALSE.
      IF (PRESENT(refpoint)) THEN
         CPASSERT((PRESENT(cell) .AND. PRESENT(particle_set))) ! need these as well if refpoint is provided
         rf = refpoint
         my_ref = .TRUE.
      END IF

      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, &
                           maxlgto=maxlgto, &
                           maxsgf=maxsgf, &
                           maxlppnl=maxlppnl, &
                           maxppnl=maxppnl)

      maxl = MAX(maxlgto, maxlppnl)
      CALL init_orbital_pointers(maxl + nder + 1)

      ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
      ldai = ncoset(maxl + nder + 1)

      !set up direct access to basis and potential
      NULLIFY (gpotential, spotential)
      ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
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

      !calculate the overlap integrals <a|pp>
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, gpotential, spotential, maxder, ncoset, my_momentmode, &
!$OMP          sap_ppnl, sap_int, nkind, ldsab, ldai, nder, nco, my_ref, cell, particle_set, rf) &
!$OMP PRIVATE (ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          slot, sphi_a, zeta, cprj, hprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, first_col, vprj_ppnl, i, j, l, dogth, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl, &
!$OMP          na, nb, np, nnl, is_nonlocal, a_nl, h_nl, c_nl, radp, raf, rcf, ra, rc)

      ! allocate work storage
      ALLOCATE (sab(ldsab, ldsab*maxder), work(ldsab, ldsab*maxder))
      sab = 0.0_dp
      ALLOCATE (ai_work(ldai, ldai, ncoset(nder + 1)))
      ai_work = 0.0_dp

!$OMP DO SCHEDULE(GUIDED)
      ! loop over neighbourlist
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

         IF (my_ref) THEN
            ra(:) = pbc(particle_set(iatom)%r(:) - rf, cell) + rf
            rc(:) = pbc(particle_set(katom)%r(:) - rf, cell) + rf
            raf(:) = ra(:) - rf(:)
            rcf(:) = rc(:) - rf(:)
         ELSE
            raf(:) = -rac
            rcf(:) = (/0._dp, 0._dp, 0._dp/)
         END IF

         dac = NORM2(rac)
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
                  IF (my_momentmode) THEN
                     CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab, 0, .FALSE., ai_work, ldai)
                  ELSE
                     CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab, nder, .TRUE., ai_work, ldai)
                  END IF
                  IF (my_momentmode .AND. nder >= 1) THEN
                     CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                 lc_max, 1, zetc, rprjc, nder, raf, rcf, ai_work)
                     ! reduce ai_work to sab
                     na = ncoa
                     np = ncoc
                     DO i = 1, maxder - 1
                        first_col = i*ldsab
                        sab(1:na, first_col + 1:first_col + np) = ai_work(1:na, 1:np, i)
                     END DO
                  END IF

                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  na = ncoa
                  nb = nprjc
                  np = ncoc
                  DO i = 1, maxder
                     first_col = (i - 1)*ldsab
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
                  clist%acint(sgfa:sgfa + na - 1, 1:nb, i) = &
                     MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work(1:np, first_col:first_col + nb - 1))

                  ! *** Multiply with interaction matrix(h) ***
                  clist%achint(sgfa:sgfa + na - 1, 1:nb, i) = &
                     MATMUL(clist%acint(sgfa:sgfa + na - 1, 1:nb, i), vprj_ppnl(1:nb, 1:nb))
               END DO
            ELSE
               ! SGP potential
               ! *** Calculate the primitive overlap integrals ***
               IF (my_momentmode) THEN
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lppnl, 0, nnl, radp, a_nl, rac, dac, sab, 0, .FALSE., ai_work, ldai)
               ELSE
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lppnl, 0, nnl, radp, a_nl, rac, dac, sab, nder, .TRUE., ai_work, ldai)
               END IF
               IF (my_momentmode .AND. nder >= 1) THEN
                  CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lppnl, nnl, a_nl, radp, nder, raf, rcf, ai_work)
                  ! reduce ai_work to sab
                  na = ncoa
                  DO i = 1, maxder - 1
                     first_col = i*ldsab
                     sab(1:na, first_col:first_col + nprjc - 1) = ai_work(1:na, 1:nprjc, i)
                  END DO
               END IF

               na = nsgf_seta(iset)
               nb = nppnl
               np = ncoa
               DO i = 1, maxder
                  first_col = (i - 1)*ldsab + 1
                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  work(1:np, 1:nb) = MATMUL(sab(1:np, first_col:first_col + nprjc - 1), cprj(1:nprjc, 1:nb))

                  ! *** Contraction step (basis functions) ***
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

      DEALLOCATE (basis_set, gpotential, spotential)

      CALL timestop(handle)

   END SUBROUTINE build_sap_ints

END MODULE sap_kind_types
