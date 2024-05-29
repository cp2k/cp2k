!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
MODULE commutator_rpnl
   USE ai_moments,                      ONLY: moment
   USE ai_overlap,                      ONLY: overlap
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
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
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              build_sap_ints,&
                                              clist_type,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'commutator_rpnl'

   PUBLIC :: build_com_rpnl, build_com_mom_nl, build_com_nl_mag, build_com_vnl_giao

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
! **************************************************************************************************
   SUBROUTINE build_com_rpnl(matrix_rv, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_com_rpnl'

      INTEGER :: handle, i, iab, iac, iatom, ibc, icol, ikind, ilist, inode, irow, iset, jatom, &
         jkind, jneighbor, kac, katom, kbc, kkind, l, lc_max, lc_min, ldai, ldsab, lppnl, maxco, &
         maxder, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, mepos, na, nb, ncoa, ncoc, nkind, &
         nlist, nneighbor, nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, gpot, ppnl_present, spot
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work, sab, work
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, x_block, &
                                                            y_block, z_block, zeta
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (ppnl_present) THEN

         nkind = SIZE(qs_kind_set)

         CALL get_qs_kind_set(qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxsgf=maxsgf, &
                              maxlppnl=maxlppnl, &
                              maxppnl=maxppnl)

         maxl = MAX(maxlgto, maxlppnl)
         CALL init_orbital_pointers(maxl + 1)

         ldsab = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
         ldai = ncoset(maxl + 1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to basis and potential
         ALLOCATE (basis_set(nkind), gpotential(nkind), spotential(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               gpotential(ikind)%gth_potential => gth_potential
               NULLIFY (spotential(ikind)%sgp_potential)
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               spotential(ikind)%sgp_potential => sgp_potential
               NULLIFY (gpotential(ikind)%gth_potential)
            ELSE
               NULLIFY (gpotential(ikind)%gth_potential)
               NULLIFY (spotential(ikind)%sgp_potential)
            END IF
         END DO

         maxder = 4
         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sap_ppnl, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, spotential, gpotential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nco ) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, &
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, vprj_ppnl, i, l, gpot, spot, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (sab(ldsab, ldsab, maxder), work(ldsab, ldsab, maxder))
         sab = 0.0_dp
         ALLOCATE (ai_work(ldai, ldai, 1))
         ai_work = 0.0_dp

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            gpot = ASSOCIATED(gpotential(kkind)%gth_potential)
            spot = ASSOCIATED(spotential(kkind)%sgp_potential)
            IF ((.NOT. gpot) .AND. (.NOT. spot)) CYCLE
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
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf

            ! get definition of PP projectors
            IF (gpot) THEN
               alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
               cprj => gpotential(kkind)%gth_potential%cprj
               lppnl = gpotential(kkind)%gth_potential%lppnl
               nppnl = gpotential(kkind)%gth_potential%nppnl
               nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
               ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
               vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
            ELSEIF (spot) THEN
               CPABORT('SGP not implemented')
            ELSE
               CPABORT('PPNL unknown')
            END IF
!$OMP CRITICAL(sap_int_critical)
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
!$OMP END CRITICAL(sap_int_critical)
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
               work = 0._dp
               prjc = 1
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
                  lc_max = l + 2*(nprj_ppnl(l) - 1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)
                  ! Calculate the primitive overlap and dipole moment integrals
                  CALL overlap(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lc_max, lc_min, 1, rprjc, zetc, rac, dac, sab(:, :, 1), 0, .FALSE., ai_work, ldai)
                  CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lc_max, 1, zetc, rprjc, 1, rac, (/0._dp, 0._dp, 0._dp/), sab(:, :, 2:4))
                  ! *** Transformation step projector functions (cartesian->spherical) ***
                  DO i = 1, maxder
                     CALL dgemm("N", "N", ncoa, nprjc, ncoc, 1.0_dp, sab(1, 1, i), ldsab, &
                                cprj(1, prjc), SIZE(cprj, 1), 0.0_dp, work(1, 1, i), ldsab)
                  END DO
                  prjc = prjc + nprjc
               END DO
               DO i = 1, maxder
                  ! Contraction step (basis functions)
                  CALL dgemm("T", "N", nsgf_seta(iset), nppnl, ncoa, 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1, i), ldsab, 0.0_dp, clist%acint(sgfa, 1, i), nsgfa)
                  ! Multiply with interaction matrix(h)
                  CALL dgemm("N", "N", nsgf_seta(iset), nppnl, nppnl, 1.0_dp, clist%acint(sgfa, 1, i), nsgfa, &
                             vprj_ppnl(1, 1), SIZE(vprj_ppnl, 1), 0.0_dp, clist%achint(sgfa, 1, i), nsgfa)
               END DO
            END DO
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         END DO

         DEALLOCATE (sab, ai_work, work)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl ) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, irow, icol, x_block, y_block, z_block, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, rac, kkind, kac, kbc, i)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind + nkind*(jkind - 1)

            ! *** Create matrix blocks for a new matrix block column ***
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
            CALL dbcsr_get_block_p(matrix_rv(1)%matrix, irow, icol, x_block, found)
            CALL dbcsr_get_block_p(matrix_rv(2)%matrix, irow, icol, y_block, found)
            CALL dbcsr_get_block_p(matrix_rv(3)%matrix, irow, icol, z_block, found)

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(x_block) .AND. ASSOCIATED(y_block) .AND. ASSOCIATED(z_block)) THEN
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
!$OMP CRITICAL(h_block_critical)
                           IF (iatom <= jatom) THEN
                              ! Vnl*r
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 2), nb, 1.0_dp, x_block, SIZE(x_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 3), nb, 1.0_dp, y_block, SIZE(y_block, 1))
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 4), nb, 1.0_dp, z_block, SIZE(z_block, 1))
                              ! -r*Vnl
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 2), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, x_block, SIZE(x_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 3), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, y_block, SIZE(y_block, 1))
                              CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 4), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, z_block, SIZE(z_block, 1))
                           ELSE
                              ! Vnl*r
                              CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 2), nb, &
                                         acint(1, 1, 1), na, 1.0_dp, x_block, SIZE(x_block, 1))
                              CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 3), nb, &
                                         acint(1, 1, 1), na, 1.0_dp, y_block, SIZE(y_block, 1))
                              CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 4), nb, &
                                         acint(1, 1, 1), na, 1.0_dp, z_block, SIZE(z_block, 1))
                              ! -r*Vnl
                              CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                                         acint(1, 1, 2), na, 1.0_dp, x_block, SIZE(x_block, 1))
                              CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                                         acint(1, 1, 3), na, 1.0_dp, y_block, SIZE(y_block, 1))
                              CALL dgemm("N", "T", nb, na, np, -1.0_dp, bchint(1, 1, 1), nb, &
                                         acint(1, 1, 4), na, 1.0_dp, z_block, SIZE(z_block, 1))
                           END IF
!$OMP END CRITICAL(h_block_critical)
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END DO
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         CALL release_sap_int(sap_int)

         DEALLOCATE (basis_set, gpotential, spotential)

      END IF !ppnl_present

      CALL timestop(handle)

   END SUBROUTINE build_com_rpnl

! **************************************************************************************************
!> \brief Calculate [r,Vnl] (matrix_rv), r x [r,Vnl] (matrix_rxrv)
!>        or [rr,Vnl] (matrix_rrv) in AO basis.
!>        Reference point is required for the two latter options
!>        Update: Calculate rxVnlxr (matrix_rvr) and rxrxVnl + Vnlxrxr (matrix_rrv_vrr)
!>        in AO basis. Added in the first place for current correction in
!>        the VG formalism (first order wrt vector potential).
!> \param qs_kind_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param cell ...
!> \param matrix_rv ...
!> \param matrix_rxrv ...
!> \param matrix_rrv ...
!> \param matrix_rvr ...
!> \param matrix_rrv_vrr ...
!> \param matrix_r_rxvr ...
!> \param matrix_rxvr_r ...
!> \param matrix_r_doublecom ...
!> \param pseudoatom ...
!> \param ref_point ...
! **************************************************************************************************
   SUBROUTINE build_com_mom_nl(qs_kind_set, sab_all, sap_ppnl, eps_ppnl, particle_set, cell, matrix_rv, matrix_rxrv, &
                    matrix_rrv, matrix_rvr, matrix_rrv_vrr, matrix_r_rxvr, matrix_rxvr_r, matrix_r_doublecom, pseudoatom, ref_point)

      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: matrix_rv, matrix_rxrv, matrix_rrv, &
                                                            matrix_rvr, matrix_rrv_vrr
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: matrix_r_rxvr, matrix_rxvr_r, &
                                                            matrix_r_doublecom
      INTEGER, INTENT(in), OPTIONAL                      :: pseudoatom
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: ref_point

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_com_mom_nl'
      INTEGER, PARAMETER :: i_x = 2, i_xx = 5, i_xy = 6, i_xz = 7, i_y = 3, i_yx = i_xy, i_yy = 8, &
         i_yz = 9, i_z = 4, i_zx = i_xz, i_zy = i_yz, i_zz = 10

      INTEGER                                            :: handle, i, iab, iac, iatom, ibc, icol, &
                                                            ikind, ind, ind2, irow, jatom, jkind, &
                                                            kac, kbc, kkind, na, natom, nb, nkind, &
                                                            np, order, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL :: asso_r_doublecom, asso_r_rxvr, asso_rrv, asso_rrv_vrr, asso_rv, asso_rvr, &
         asso_rxrv, asso_rxvr_r, do_symmetric, found, go, my_r_doublecom, my_r_rxvr, my_ref, &
         my_rrv, my_rrv_vrr, my_rv, my_rvr, my_rxrv, my_rxvr_r, periodic, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rf
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: blocks_rrv, blocks_rrv_vrr, blocks_rv, &
                                                            blocks_rvr, blocks_rxrv
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:, :)   :: blocks_r_doublecom, blocks_r_rxvr, &
                                                            blocks_rxvr_r
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)

      my_r_doublecom = .FALSE.
      my_r_rxvr = .FALSE.
      my_rxvr_r = .FALSE.
      my_rxrv = .FALSE.
      my_rrv = .FALSE.
      my_rv = .FALSE.
      my_rvr = .FALSE.
      my_rrv_vrr = .FALSE.
      IF (PRESENT(matrix_r_doublecom)) my_r_doublecom = .TRUE.
      IF (PRESENT(matrix_r_rxvr)) my_r_rxvr = .TRUE.
      IF (PRESENT(matrix_rxvr_r)) my_rxvr_r = .TRUE.
      IF (PRESENT(matrix_rxrv)) my_rxrv = .TRUE.
      IF (PRESENT(matrix_rrv)) my_rrv = .TRUE.
      IF (PRESENT(matrix_rv)) my_rv = .TRUE.
      IF (PRESENT(matrix_rvr)) my_rvr = .TRUE.
      IF (PRESENT(matrix_rrv_vrr)) my_rrv_vrr = .TRUE.
      IF (.NOT. (my_rv .OR. my_rxrv .OR. my_rrv .OR. my_rvr .OR. my_rrv_vrr .OR. my_r_rxvr .OR. my_rxvr_r .OR. my_r_doublecom)) THEN
         CPABORT('No dbcsr matrix provided for commutator calculation!')
      END IF

      natom = SIZE(particle_set)

      IF (my_rxrv .OR. my_rrv .OR. my_r_rxvr .OR. my_rxvr_r .OR. my_r_doublecom) THEN
         order = 2
         CPASSERT(PRESENT(ref_point)) ! need reference point for r x [r,Vnl] and [rr,Vnl]
      ELSE IF (my_rvr .OR. my_rrv_vrr) THEN
         order = 2
      ELSE
         order = 1
      END IF

      ! When we want the double commutator [[Vnl, r], r], we also want to fix the pseudoatom
      IF (my_r_doublecom) THEN
         CPASSERT(PRESENT(pseudoatom))
      END IF

      periodic = ANY(cell%perd > 0)
      my_ref = .FALSE.
      IF (PRESENT(ref_point)) THEN
         IF (.NOT. periodic) THEN
            rf = ref_point
            my_ref = .TRUE.
         ELSE ! use my_ref = False in periodic case, corresponds to distributed ref point
            IF (order .GT. 1) THEN
               CPWARN("Not clear how to define reference point for order > 1 in periodic cells.")
            END IF
         END IF
      END IF

      nkind = SIZE(qs_kind_set)

      !sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      IF (my_ref) THEN
         ! calculate integrals <a|x^n|p>
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE., refpoint=rf, &
                             particle_set=particle_set, cell=cell)
      ELSE
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE.)
      END IF

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_all, symmetric=do_symmetric)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, matrix_rxrv, matrix_rrv, &
!$OMP          matrix_rvr, matrix_rrv_vrr, matrix_r_doublecom, &
!$OMP          sap_int, natom, nkind, eps_ppnl, locks, sab_all, &
!$OMP          my_rv, my_rxrv, my_rrv, my_rvr, my_rrv_vrr, &
!$OMP          my_r_doublecom, &
!$OMP          matrix_r_rxvr, matrix_rxvr_r, my_r_rxvr, my_rxvr_r, &
!$OMP          pseudoatom, do_symmetric) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, &
!$OMP          blocks_rv, blocks_rxrv, blocks_rrv, blocks_rvr, blocks_rrv_vrr, &
!$OMP          blocks_r_rxvr, blocks_rxvr_r, blocks_r_doublecom, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, &
!$OMP          go, asso_rv, asso_rxrv, asso_rrv, asso_rvr, asso_rrv_vrr, &
!$OMP          asso_r_rxvr, asso_rxvr_r, asso_r_doublecom, hash, &
!$OMP          acint, achint, bcint, bchint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_all(1)%nl_size

         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_all(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

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

         ! allocate blocks
         IF (my_rv) THEN
            ALLOCATE (blocks_rv(3))
         END IF
         IF (my_rxrv) THEN
            ALLOCATE (blocks_rxrv(3))
         END IF
         IF (my_rrv) THEN
            ALLOCATE (blocks_rrv(6))
         END IF
         IF (my_rvr) THEN
            ALLOCATE (blocks_rvr(6))
         END IF
         IF (my_rrv_vrr) THEN
            ALLOCATE (blocks_rrv_vrr(6))
         END IF
         IF (my_r_rxvr) THEN
            ALLOCATE (blocks_r_rxvr(3, 3))
         END IF

         IF (my_rxvr_r) THEN
            ALLOCATE (blocks_rxvr_r(3, 3))
         END IF

         IF (my_r_doublecom) THEN
            ALLOCATE (blocks_r_doublecom(3, 3))
         END IF

         ! get blocks
         IF (my_rv) THEN
            DO ind = 1, 3
               CALL dbcsr_get_block_p(matrix_rv(ind)%matrix, irow, icol, blocks_rv(ind)%block, found)
            END DO
         END IF

         IF (my_rxrv) THEN
            DO ind = 1, 3
               CALL dbcsr_get_block_p(matrix_rxrv(ind)%matrix, irow, icol, blocks_rxrv(ind)%block, found)
               blocks_rxrv(ind)%block(:, :) = 0._dp
            END DO
         END IF

         IF (my_rrv) THEN
            DO ind = 1, 6
               CALL dbcsr_get_block_p(matrix_rrv(ind)%matrix, irow, icol, blocks_rrv(ind)%block, found)
            END DO
         END IF

         IF (my_rvr) THEN
            DO ind = 1, 6
               CALL dbcsr_get_block_p(matrix_rvr(ind)%matrix, irow, icol, blocks_rvr(ind)%block, found)
            END DO
         END IF

         IF (my_rrv_vrr) THEN
            DO ind = 1, 6
               CALL dbcsr_get_block_p(matrix_rrv_vrr(ind)%matrix, irow, icol, blocks_rrv_vrr(ind)%block, found)
            END DO
         END IF

         IF (my_r_rxvr) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  CALL dbcsr_get_block_p(matrix_r_rxvr(ind, ind2)%matrix, irow, icol, &
                                         blocks_r_rxvr(ind, ind2)%block, found)
                  blocks_r_rxvr(ind, ind2)%block(:, :) = 0._dp
               END DO
            END DO
         END IF

         IF (my_rxvr_r) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  CALL dbcsr_get_block_p(matrix_rxvr_r(ind, ind2)%matrix, irow, icol, &
                                         blocks_rxvr_r(ind, ind2)%block, found)
                  blocks_rxvr_r(ind, ind2)%block(:, :) = 0._dp
               END DO
            END DO
         END IF

         IF (my_r_doublecom) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  CALL dbcsr_get_block_p(matrix_r_doublecom(ind, ind2)%matrix, irow, icol, &
                                         blocks_r_doublecom(ind, ind2)%block, found)
                  blocks_r_doublecom(ind, ind2)%block(:, :) = 0._dp
               END DO
            END DO
         END IF

         ! check whether all blocks are associated
         go = .TRUE.
         IF (my_rv) THEN
            asso_rv = (ASSOCIATED(blocks_rv(1)%block) .AND. ASSOCIATED(blocks_rv(2)%block) .AND. &
                       ASSOCIATED(blocks_rv(3)%block))
            go = go .AND. asso_rv
         END IF

         IF (my_rxrv) THEN
            asso_rxrv = (ASSOCIATED(blocks_rxrv(1)%block) .AND. ASSOCIATED(blocks_rxrv(2)%block) .AND. &
                         ASSOCIATED(blocks_rxrv(3)%block))
            go = go .AND. asso_rxrv
         END IF

         IF (my_rrv) THEN
            asso_rrv = (ASSOCIATED(blocks_rrv(1)%block) .AND. ASSOCIATED(blocks_rrv(2)%block) .AND. &
                        ASSOCIATED(blocks_rrv(3)%block) .AND. ASSOCIATED(blocks_rrv(4)%block) .AND. &
                        ASSOCIATED(blocks_rrv(5)%block) .AND. ASSOCIATED(blocks_rrv(6)%block))
            go = go .AND. asso_rrv
         END IF

         IF (my_rvr) THEN
            asso_rvr = (ASSOCIATED(blocks_rvr(1)%block) .AND. ASSOCIATED(blocks_rvr(2)%block) .AND. &
                        ASSOCIATED(blocks_rvr(3)%block) .AND. ASSOCIATED(blocks_rvr(4)%block) .AND. &
                        ASSOCIATED(blocks_rvr(5)%block) .AND. ASSOCIATED(blocks_rvr(6)%block))
            go = go .AND. asso_rvr
         END IF

         IF (my_rrv_vrr) THEN
            asso_rrv_vrr = (ASSOCIATED(blocks_rrv_vrr(1)%block) .AND. ASSOCIATED(blocks_rrv_vrr(2)%block) .AND. &
                            ASSOCIATED(blocks_rrv_vrr(3)%block) .AND. ASSOCIATED(blocks_rrv_vrr(4)%block) .AND. &
                            ASSOCIATED(blocks_rrv_vrr(5)%block) .AND. ASSOCIATED(blocks_rrv_vrr(6)%block))
            go = go .AND. asso_rrv_vrr
         END IF

         IF (my_r_rxvr) THEN
            asso_r_rxvr = .TRUE.
            DO ind = 1, 3
               DO ind2 = 1, 3
                  asso_r_rxvr = asso_r_rxvr .AND. ASSOCIATED(blocks_r_rxvr(ind, ind2)%block)
               END DO
            END DO
            go = go .AND. asso_r_rxvr
         END IF

         IF (my_rxvr_r) THEN
            asso_rxvr_r = .TRUE.
            DO ind = 1, 3
               DO ind2 = 1, 3
                  asso_rxvr_r = asso_rxvr_r .AND. ASSOCIATED(blocks_rxvr_r(ind, ind2)%block)
               END DO
            END DO
            go = go .AND. asso_rxvr_r
         END IF

         IF (my_r_doublecom) THEN
            asso_r_doublecom = .TRUE.
            DO ind = 1, 3
               DO ind2 = 1, 3
                  asso_r_doublecom = asso_r_doublecom .AND. ASSOCIATED(blocks_r_doublecom(ind, ind2)%block)
               END DO
            END DO
            go = go .AND. asso_r_doublecom
         END IF

         ! loop over all kinds for projector atom
         ! < iatom | katom > h < katom | jatom >
         IF (go) THEN
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
                     IF (PRESENT(pseudoatom)) THEN
                        IF (alist_ac%clist(kac)%catom /= pseudoatom) CYCLE
                     END IF

                     IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                        IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                        acint => alist_ac%clist(kac)%acint
                        bcint => alist_bc%clist(kbc)%acint
                        achint => alist_ac%clist(kac)%achint
                        bchint => alist_bc%clist(kbc)%achint
                        na = SIZE(acint, 1)
                        np = SIZE(acint, 2)
                        nb = SIZE(bcint, 1)
!$                      hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                      CALL omp_set_lock(locks(hash))
                        IF (my_rv) THEN
                           ! r*Vnl
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rv(1)%block, SIZE(blocks_rv(1)%block, 1)) ! xV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rv(2)%block, SIZE(blocks_rv(2)%block, 1)) ! yV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 4), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rv(3)%block, SIZE(blocks_rv(3)%block, 1)) ! zV
                           IF (iatom <= jatom) THEN
                              ! with MATMUL
                              blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                                                               MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1))) ! xV
                              blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) + &
                                                               MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 1))) ! yV
                              blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) + &
                                                               MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 1))) ! zV
                           ELSE
                              blocks_rv(1)%block(1:nb, 1:na) = blocks_rv(1)%block(1:nb, 1:na) + &
                                                               MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 1)))
                              blocks_rv(2)%block(1:nb, 1:na) = blocks_rv(2)%block(1:nb, 1:na) + &
                                                               MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 1)))
                              blocks_rv(3)%block(1:nb, 1:na) = blocks_rv(3)%block(1:nb, 1:na) + &
                                                               MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 1)))
                           END IF
                           ! -Vnl r
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 2), nb, 1.0_dp, blocks_rv(1)%block, SIZE(blocks_rv(1)%block, 1)) ! -Vx
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks_rv(2)%block, SIZE(blocks_rv(2)%block, 1)) ! -Vy
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks_rv(3)%block, SIZE(blocks_rv(3)%block, 1)) ! -Vz
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) - &
                                                               MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 2))) ! -Vx
                              blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) - &
                                                               MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 3))) ! -Vy
                              blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) - &
                                                               MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 4))) ! -Vz
                           ELSE
                              blocks_rv(1)%block(1:nb, 1:na) = blocks_rv(1)%block(1:nb, 1:na) - &
                                                               MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 2)))
                              blocks_rv(2)%block(1:nb, 1:na) = blocks_rv(2)%block(1:nb, 1:na) - &
                                                               MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 3)))
                              blocks_rv(3)%block(1:nb, 1:na) = blocks_rv(3)%block(1:nb, 1:na) - &
                                                               MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 4)))
                           END IF

                        END IF

                        IF (my_rxrv) THEN
                           ! x-component (y [z,Vnl] - z [y, Vnl])
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 9), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(1)%block, SIZE(blocks_rxrv(1)%block, 1)) ! yzV
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks_rxrv(1)%block, SIZE(blocks_rxrv(1)%block, 1)) ! -yVz
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 9), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(1)%block, SIZE(blocks_rxrv(1)%block, 1)) ! -zyV
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 4), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks_rxrv(1)%block, SIZE(blocks_rxrv(1)%block, 1)) ! zVy
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rxrv(1)%block(1:na, 1:nb) = blocks_rxrv(1)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 9), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! yzV
                              blocks_rxrv(1)%block(1:na, 1:nb) = blocks_rxrv(1)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 4)))    ! -yVz
                              blocks_rxrv(1)%block(1:na, 1:nb) = blocks_rxrv(1)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 9), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! -zyV
                              blocks_rxrv(1)%block(1:na, 1:nb) = blocks_rxrv(1)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 3)))    ! zVy
                           ELSE
                              blocks_rxrv(1)%block(1:nb, 1:na) = blocks_rxrv(1)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 9), TRANSPOSE(acint(1:na, 1:np, 1)))    ! yzV
                              blocks_rxrv(1)%block(1:nb, 1:na) = blocks_rxrv(1)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 4)))    ! -yVz
                              blocks_rxrv(1)%block(1:nb, 1:na) = blocks_rxrv(1)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 9), TRANSPOSE(acint(1:na, 1:np, 1)))    ! -zyV
                              blocks_rxrv(1)%block(1:nb, 1:na) = blocks_rxrv(1)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 3)))    ! zVy
                           END IF

                           ! y-component (z [x,Vnl] - x [z, Vnl])
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 7), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(2)%block, SIZE(blocks_rxrv(2)%block, 1)) ! zxV
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 4), na, &
                           !            bcint(1, 1, 2), nb, 1.0_dp, blocks_rxrv(2)%block, SIZE(blocks_rxrv(2)%block, 1)) ! -zVx
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 7), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(2)%block, SIZE(blocks_rxrv(2)%block, 1)) ! -xzV
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 4), nb, 1.0_dp, blocks_rxrv(2)%block, SIZE(blocks_rxrv(2)%block, 1)) ! xVz
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rxrv(2)%block(1:na, 1:nb) = blocks_rxrv(2)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 7), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! zxV
                              blocks_rxrv(2)%block(1:na, 1:nb) = blocks_rxrv(2)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 2)))    ! -zVx
                              blocks_rxrv(2)%block(1:na, 1:nb) = blocks_rxrv(2)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 7), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! -xzV
                              blocks_rxrv(2)%block(1:na, 1:nb) = blocks_rxrv(2)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 4)))    ! xVz
                           ELSE
                              blocks_rxrv(2)%block(1:nb, 1:na) = blocks_rxrv(2)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 7), TRANSPOSE(acint(1:na, 1:np, 1)))    ! zxV
                              blocks_rxrv(2)%block(1:nb, 1:na) = blocks_rxrv(2)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 2)))    ! -zVx
                              blocks_rxrv(2)%block(1:nb, 1:na) = blocks_rxrv(2)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 7), TRANSPOSE(acint(1:na, 1:np, 1)))    ! -xzV
                              blocks_rxrv(2)%block(1:nb, 1:na) = blocks_rxrv(2)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 4)))    ! xVz
                           END IF

                           ! z-component (x [y,Vnl] - y [x, Vnl])
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 6), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(3)%block, SIZE(blocks_rxrv(3)%block, 1)) ! xyV
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 2), na, &
                           !            bcint(1, 1, 3), nb, 1.0_dp, blocks_rxrv(3)%block, SIZE(blocks_rxrv(3)%block, 1)) ! -xVy
                           ! CALL dgemm("N", "T", na, nb, np, -1.0_dp, achint(1, 1, 6), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rxrv(3)%block, SIZE(blocks_rxrv(3)%block, 1)) ! -yxV
                           ! CALL dgemm("N", "T", na, nb, np, 1.0_dp, achint(1, 1, 3), na, &
                           !            bcint(1, 1, 2), nb, 1.0_dp, blocks_rxrv(3)%block, SIZE(blocks_rxrv(3)%block, 1)) ! yVx
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rxrv(3)%block(1:na, 1:nb) = blocks_rxrv(3)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 6), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! xyV
                              blocks_rxrv(3)%block(1:na, 1:nb) = blocks_rxrv(3)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 3)))    ! -xVy
                              blocks_rxrv(3)%block(1:na, 1:nb) = blocks_rxrv(3)%block(1:na, 1:nb) - &
                                                                 MATMUL(achint(1:na, 1:np, 6), TRANSPOSE(bcint(1:nb, 1:np, 1)))    ! -yxV
                              blocks_rxrv(3)%block(1:na, 1:nb) = blocks_rxrv(3)%block(1:na, 1:nb) + &
                                                                 MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 2)))    ! zVx
                           ELSE
                              blocks_rxrv(3)%block(1:nb, 1:na) = blocks_rxrv(3)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 6), TRANSPOSE(acint(1:na, 1:np, 1)))    ! xyV
                              blocks_rxrv(3)%block(1:nb, 1:na) = blocks_rxrv(3)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 3)))    ! -xVy
                              blocks_rxrv(3)%block(1:nb, 1:na) = blocks_rxrv(3)%block(1:nb, 1:na) - &
                                                                 MATMUL(bchint(1:nb, 1:np, 6), TRANSPOSE(acint(1:na, 1:np, 1)))    ! -yxV
                              blocks_rxrv(3)%block(1:nb, 1:na) = blocks_rxrv(3)%block(1:nb, 1:na) + &
                                                                 MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 2)))    ! zVx
                           END IF
                        END IF

                        IF (my_rrv) THEN
                           ! r_alpha * r_beta * Vnl
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 5), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(1)%block, SIZE(blocks_rrv(1)%block, 1)) ! xxV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 6), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(2)%block, SIZE(blocks_rrv(2)%block, 1)) ! xyV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 7), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(3)%block, SIZE(blocks_rrv(3)%block, 1)) ! xzV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 8), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(4)%block, SIZE(blocks_rrv(4)%block, 1)) ! yyV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 9), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(5)%block, SIZE(blocks_rrv(5)%block, 1)) ! yzV
                           ! CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 10), na, &
                           !            bcint(1, 1, 1), nb, 1.0_dp, blocks_rrv(6)%block, SIZE(blocks_rrv(6)%block, 1)) ! zzV
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rrv(1)%block(1:na, 1:nb) = blocks_rrv(1)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 5), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xxV
                              blocks_rrv(2)%block(1:na, 1:nb) = blocks_rrv(2)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 6), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xyV
                              blocks_rrv(3)%block(1:na, 1:nb) = blocks_rrv(3)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 7), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xzV
                              blocks_rrv(4)%block(1:na, 1:nb) = blocks_rrv(4)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 8), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! yyV
                              blocks_rrv(5)%block(1:na, 1:nb) = blocks_rrv(5)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 9), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! yzV
                              blocks_rrv(6)%block(1:na, 1:nb) = blocks_rrv(6)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 10), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! zzV
                           ELSE
                              blocks_rrv(1)%block(1:nb, 1:na) = blocks_rrv(1)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 5), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xxV
                              blocks_rrv(2)%block(1:nb, 1:na) = blocks_rrv(2)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 6), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xyV
                              blocks_rrv(3)%block(1:nb, 1:na) = blocks_rrv(3)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 7), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xzV
                              blocks_rrv(4)%block(1:nb, 1:na) = blocks_rrv(4)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 8), TRANSPOSE(acint(1:na, 1:np, 1)))   ! yyV
                              blocks_rrv(5)%block(1:nb, 1:na) = blocks_rrv(5)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 9), TRANSPOSE(acint(1:na, 1:np, 1)))   ! yzV
                              blocks_rrv(6)%block(1:nb, 1:na) = blocks_rrv(6)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 10), TRANSPOSE(acint(1:na, 1:np, 1)))   ! zzV
                           END IF

                           ! - Vnl * r_alpha * r_beta
                           ! with LAPACK
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 5), nb, 1.0_dp, blocks_rrv(1)%block, SIZE(blocks_rrv(1)%block, 1)) ! Vxx
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 6), nb, 1.0_dp, blocks_rrv(2)%block, SIZE(blocks_rrv(2)%block, 1)) ! Vxy
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 7), nb, 1.0_dp, blocks_rrv(3)%block, SIZE(blocks_rrv(3)%block, 1)) ! Vxz
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 8), nb, 1.0_dp, blocks_rrv(4)%block, SIZE(blocks_rrv(4)%block, 1)) ! Vyy
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 9), nb, 1.0_dp, blocks_rrv(5)%block, SIZE(blocks_rrv(5)%block, 1)) ! Vyz
                           ! CALL dgemm("N", "T", na, nb, np, -1._dp, achint(1, 1, 1), na, &
                           !            bcint(1, 1, 10), nb, 1.0_dp, blocks_rrv(6)%block, SIZE(blocks_rrv(6)%block, 1)) ! Vzz
                           ! with MATMUL
                           IF (iatom <= jatom) THEN
                              blocks_rrv(1)%block(1:na, 1:nb) = blocks_rrv(1)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 5)))   ! -Vxx
                              blocks_rrv(2)%block(1:na, 1:nb) = blocks_rrv(2)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 6)))   ! -Vxy
                              blocks_rrv(3)%block(1:na, 1:nb) = blocks_rrv(3)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 7)))   ! -Vxz
                              blocks_rrv(4)%block(1:na, 1:nb) = blocks_rrv(4)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 8)))   ! -Vyy
                              blocks_rrv(5)%block(1:na, 1:nb) = blocks_rrv(5)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 9)))   ! -Vyz
                              blocks_rrv(6)%block(1:na, 1:nb) = blocks_rrv(6)%block(1:na, 1:nb) - &
                                                                MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 10)))   ! -Vzz
                           ELSE
                              blocks_rrv(1)%block(1:nb, 1:na) = blocks_rrv(1)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 5)))   ! -Vxx
                              blocks_rrv(2)%block(1:nb, 1:na) = blocks_rrv(2)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 6)))   ! -Vxy
                              blocks_rrv(3)%block(1:nb, 1:na) = blocks_rrv(3)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 7)))   ! -Vxz
                              blocks_rrv(4)%block(1:nb, 1:na) = blocks_rrv(4)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 8)))   ! -Vyy
                              blocks_rrv(5)%block(1:nb, 1:na) = blocks_rrv(5)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 9)))   ! -Vyz
                              blocks_rrv(6)%block(1:nb, 1:na) = blocks_rrv(6)%block(1:nb, 1:na) - &
                                                                MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 10)))   ! -Vzz
                           END IF
                        END IF

                        IF (my_rvr) THEN
                           ! r_alpha * Vnl * r_beta
                           IF (iatom <= jatom) THEN
                              blocks_rvr(1)%block(1:na, 1:nb) = blocks_rvr(1)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 2)))   ! xVx
                              blocks_rvr(2)%block(1:na, 1:nb) = blocks_rvr(2)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 3)))   ! xVy
                              blocks_rvr(3)%block(1:na, 1:nb) = blocks_rvr(3)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 4)))   ! xVz
                              blocks_rvr(4)%block(1:na, 1:nb) = blocks_rvr(4)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 3)))   ! yVy
                              blocks_rvr(5)%block(1:na, 1:nb) = blocks_rvr(5)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 4)))   ! yVz
                              blocks_rvr(6)%block(1:na, 1:nb) = blocks_rvr(6)%block(1:na, 1:nb) + &
                                                                MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 4)))   ! zVz
                           ELSE
                              blocks_rvr(1)%block(1:nb, 1:na) = blocks_rvr(1)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 2)))   ! xVx
                              blocks_rvr(2)%block(1:nb, 1:na) = blocks_rvr(2)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 3)))   ! xVy
                              blocks_rvr(3)%block(1:nb, 1:na) = blocks_rvr(3)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 4)))   ! xVz
                              blocks_rvr(4)%block(1:nb, 1:na) = blocks_rvr(4)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 3)))   ! yVy
                              blocks_rvr(5)%block(1:nb, 1:na) = blocks_rvr(5)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 4)))   ! yVz
                              blocks_rvr(6)%block(1:nb, 1:na) = blocks_rvr(6)%block(1:nb, 1:na) + &
                                                                MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 4)))   ! zVz
                           END IF
                        END IF

                        IF (my_rrv_vrr) THEN
                           ! r_alpha * r_beta * Vnl
                           IF (iatom <= jatom) THEN
                              blocks_rrv_vrr(1)%block(1:na, 1:nb) = blocks_rrv_vrr(1)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 5), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xxV
                              blocks_rrv_vrr(2)%block(1:na, 1:nb) = blocks_rrv_vrr(2)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 6), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xyV
                              blocks_rrv_vrr(3)%block(1:na, 1:nb) = blocks_rrv_vrr(3)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 7), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! xzV
                              blocks_rrv_vrr(4)%block(1:na, 1:nb) = blocks_rrv_vrr(4)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 8), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! yyV
                              blocks_rrv_vrr(5)%block(1:na, 1:nb) = blocks_rrv_vrr(5)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 9), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! yzV
                              blocks_rrv_vrr(6)%block(1:na, 1:nb) = blocks_rrv_vrr(6)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, 10), TRANSPOSE(bcint(1:nb, 1:np, 1)))   ! zzV
                           ELSE
                              blocks_rrv_vrr(1)%block(1:nb, 1:na) = blocks_rrv_vrr(1)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 5), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xxV
                              blocks_rrv_vrr(2)%block(1:nb, 1:na) = blocks_rrv_vrr(2)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 6), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xyV
                              blocks_rrv_vrr(3)%block(1:nb, 1:na) = blocks_rrv_vrr(3)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 7), TRANSPOSE(acint(1:na, 1:np, 1)))   ! xzV
                              blocks_rrv_vrr(4)%block(1:nb, 1:na) = blocks_rrv_vrr(4)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 8), TRANSPOSE(acint(1:na, 1:np, 1)))   ! yyV
                              blocks_rrv_vrr(5)%block(1:nb, 1:na) = blocks_rrv_vrr(5)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 9), TRANSPOSE(acint(1:na, 1:np, 1)))   ! yzV
                              blocks_rrv_vrr(6)%block(1:nb, 1:na) = blocks_rrv_vrr(6)%block(1:nb, 1:na) + &
                                                             MATMUL(bchint(1:nb, 1:np, 10), TRANSPOSE(acint(1:na, 1:np, 1)))   ! zzV
                           END IF
                           ! + Vnl * r_alpha * r_beta
                           IF (iatom <= jatom) THEN
                              blocks_rrv_vrr(1)%block(1:na, 1:nb) = blocks_rrv_vrr(1)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 5)))   ! +Vxx
                              blocks_rrv_vrr(2)%block(1:na, 1:nb) = blocks_rrv_vrr(2)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 6)))   ! +Vxy
                              blocks_rrv_vrr(3)%block(1:na, 1:nb) = blocks_rrv_vrr(3)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 7)))   ! +Vxz
                              blocks_rrv_vrr(4)%block(1:na, 1:nb) = blocks_rrv_vrr(4)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 8)))   ! +Vyy
                              blocks_rrv_vrr(5)%block(1:na, 1:nb) = blocks_rrv_vrr(5)%block(1:na, 1:nb) + &
                                                                    MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 9)))   ! +Vyz
                              blocks_rrv_vrr(6)%block(1:na, 1:nb) = blocks_rrv_vrr(6)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 10)))   ! +Vzz
                           ELSE
                              blocks_rrv_vrr(1)%block(1:nb, 1:na) = blocks_rrv_vrr(1)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 5)))   ! +Vxx
                              blocks_rrv_vrr(2)%block(1:nb, 1:na) = blocks_rrv_vrr(2)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 6)))   ! +Vxy
                              blocks_rrv_vrr(3)%block(1:nb, 1:na) = blocks_rrv_vrr(3)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 7)))   ! +Vxz
                              blocks_rrv_vrr(4)%block(1:nb, 1:na) = blocks_rrv_vrr(4)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 8)))   ! +Vyy
                              blocks_rrv_vrr(5)%block(1:nb, 1:na) = blocks_rrv_vrr(5)%block(1:nb, 1:na) + &
                                                                    MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 9)))   ! +Vyz
                              blocks_rrv_vrr(6)%block(1:nb, 1:na) = blocks_rrv_vrr(6)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 10)))   ! +Vzz
                           END IF
                        END IF

                        ! The indices are stored in i_1, i_x, ..., i_zzz
                        ! matrix_r_rxvr(alpha, beta)
                        !  = sum_(gamma delta) epsilon_(alpha gamma delta)
                        !      (r_beta * r_gamma * V_nl * r_delta - r_beta * r_gamma * r_delta * V_nl)
                        !  = sum_(gamma delta) epsilon_(alpha gamma delta) r_beta * r_gamma * V_nl * r_delta

                        ! TODO: is this set to zero before?
                        IF (my_r_rxvr) THEN
                           ! beta = 1
                           ! matrix_r_rxvr(x, x) = x * y * V_nl * z  -  x * z * V_nl * y
                           blocks_r_rxvr(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xy), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_rxvr(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xz), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_rxvr(y, x) = x * z * V_nl * x  -  x * x * V_nl * z
                           blocks_r_rxvr(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xz), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_rxvr(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xx), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_rxvr(z, x) = x * x * V_nl * y  -  x * y * V_nl * x
                           blocks_r_rxvr(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xx), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_rxvr(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xy), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                           ! beta = 2
                           ! matrix_r_rxvr(x, y) = y * y * V_nl * z  -  y * z * V_nl * y
                           blocks_r_rxvr(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yy), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_rxvr(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yz), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_rxvr(y, y) = y * z * V_nl * x  -  y * x * V_nl * z
                           blocks_r_rxvr(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yz), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_rxvr(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yx), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_rxvr(z, y) = y * x * V_nl * y  -  y * y * V_nl * x
                           blocks_r_rxvr(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yx), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_rxvr(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yy), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                           ! beta = 3
                           ! matrix_r_rxvr(x, z) = z * y * V_nl * z  -  z * z * V_nl * y
                           blocks_r_rxvr(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zy), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_rxvr(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(1, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zz), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_rxvr(y, z) = z * z * V_nl * x  -  z * x * V_nl * z
                           blocks_r_rxvr(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zz), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_rxvr(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(2, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zx), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_rxvr(z, z) = z * x * V_nl * y  -  z * y * V_nl * x
                           blocks_r_rxvr(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zx), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_rxvr(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_rxvr(3, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zy), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                        END IF ! my_r_rxvr

                        ! The indices are stored in i_1, i_x, ..., i_zzz
                        ! This will put into blocks_rxvr_r
                        ! matrix_rxvr_r(alpha, beta) = sum_(gamma delta) epsilon_(alpha gamma delta)
                        !                              r_gamma * V_nl * r_delta * r_beta
                        IF (my_rxvr_r) THEN
                           ! beta = 1
                           ! matrix_rxvr_r(x, x) = yV zx - zV yx
                           blocks_rxvr_r(1, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_zx)))
                           blocks_rxvr_r(1, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_yx)))

                           ! matrix_rxvr_r(y, x) = zV xx - xV zx
                           blocks_rxvr_r(2, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_xx)))
                           blocks_rxvr_r(2, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_zx)))

                           ! matrix_rxvr_r(z, x) = xV yx - yV xx
                           blocks_rxvr_r(3, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_yx)))
                           blocks_rxvr_r(3, 1)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_xx)))

                           ! beta = 2
                           ! matrix_rxvr_r(x, y) = yV zy - zV yy
                           blocks_rxvr_r(1, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_zy)))
                           blocks_rxvr_r(1, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_yy)))

                           ! matrix_rxvr_r(y, y) = zV xy - xV zy
                           blocks_rxvr_r(2, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_xy)))
                           blocks_rxvr_r(2, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_zy)))

                           ! matrix_rxvr_r(z, y) = xV yy - yV xy
                           blocks_rxvr_r(3, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_yy)))
                           blocks_rxvr_r(3, 2)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_xy)))

                           ! beta = 3
                           ! matrix_rxvr_r(x, z) = yV zz - zV yz
                           blocks_rxvr_r(1, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_zz)))
                           blocks_rxvr_r(1, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(1, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_yz)))

                           ! matrix_rxvr_r(y, z) = zV xz - xV zz
                           blocks_rxvr_r(2, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_xz)))
                           blocks_rxvr_r(2, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(2, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_zz)))

                           ! matrix_rxvr_r(z, z) = xV yz - yV xz
                           blocks_rxvr_r(3, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_yz)))
                           blocks_rxvr_r(3, 3)%block(1:na, 1:nb) = &
                              blocks_rxvr_r(3, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_xz)))

                        END IF ! my_rxvr_r

                        ! matrix_r_doublecom(alpha, beta) = sum_(gamma delta) epsilon_(alpha gamma delta)
                        !                                gamma V^pseudoatom beta delta - gamma beta V^pseudoatom delta

                        IF (my_r_doublecom) THEN
                           ! beta = 1
                           ! matrix_r_doublecom(x, x) = yV xz  -  zV xy  -  yxV z  +  zxV y
                           blocks_r_doublecom(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_xz)))
                           blocks_r_doublecom(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_xy)))
                           blocks_r_doublecom(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yx), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_doublecom(1, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zx), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_doublecom(y, x) = zV xx  -  xV xz  -  zxV x  +  xxV z
                           blocks_r_doublecom(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_xx)))
                           blocks_r_doublecom(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_xz)))
                           blocks_r_doublecom(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zx), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_doublecom(2, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xx), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_doublecom(z, x) = xV xy  -  yV xx  -  xxV y  +  yxV x
                           blocks_r_doublecom(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_xy)))
                           blocks_r_doublecom(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_xx)))
                           blocks_r_doublecom(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 1)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xx), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_doublecom(3, 1)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 1)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yx), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                           ! beta = 2
                           ! matrix_r_doublecom(x, y) = yV yz  -  zV yy  -  yyV z  +  zyV y
                           blocks_r_doublecom(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_yz)))
                           blocks_r_doublecom(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_yy)))
                           blocks_r_doublecom(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yy), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_doublecom(1, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zy), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_doublecom(y, y) = zV yx  -  xV yz  -  zyV x  +  xyV z
                           blocks_r_doublecom(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_yx)))
                           blocks_r_doublecom(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_yz)))
                           blocks_r_doublecom(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zy), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_doublecom(2, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xy), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_doublecom(z, y) = xV yy  -  yV yx  -  xyV y  +  yyV x
                           blocks_r_doublecom(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_yy)))
                           blocks_r_doublecom(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_yx)))
                           blocks_r_doublecom(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 2)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xy), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_doublecom(3, 2)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 2)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yy), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                           ! beta = 3
                           ! matrix_r_doublecom(x, z) = yV zz  -  zV zy  -  yzV z  +  zzV y
                           blocks_r_doublecom(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_zz)))
                           blocks_r_doublecom(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_zy)))
                           blocks_r_doublecom(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_yz), TRANSPOSE(bcint(1:nb, 1:np, i_z)))
                           blocks_r_doublecom(1, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(1, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_zz), TRANSPOSE(bcint(1:nb, 1:np, i_y)))

                           ! matrix_r_doublecom(y, z) = zV zx  -  xV zz  -  zzV x  +  xzV z
                           blocks_r_doublecom(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_z), TRANSPOSE(bcint(1:nb, 1:np, i_zx)))
                           blocks_r_doublecom(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_zz)))
                           blocks_r_doublecom(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_zz), TRANSPOSE(bcint(1:nb, 1:np, i_x)))
                           blocks_r_doublecom(2, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(2, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_xz), TRANSPOSE(bcint(1:nb, 1:np, i_z)))

                           ! matrix_r_doublecom(z, z) = xV zy  -  yV zx  -  xzV y  +  yzV x
                           blocks_r_doublecom(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_x), TRANSPOSE(bcint(1:nb, 1:np, i_zy)))
                           blocks_r_doublecom(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_y), TRANSPOSE(bcint(1:nb, 1:np, i_zx)))
                           blocks_r_doublecom(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 3)%block(1:na, 1:nb) - &
                              MATMUL(achint(1:na, 1:np, i_xz), TRANSPOSE(bcint(1:nb, 1:np, i_y)))
                           blocks_r_doublecom(3, 3)%block(1:na, 1:nb) = &
                              blocks_r_doublecom(3, 3)%block(1:na, 1:nb) + &
                              MATMUL(achint(1:na, 1:np, i_yz), TRANSPOSE(bcint(1:nb, 1:np, i_x)))

                        END IF ! my_r_doublecom
!$                      CALL omp_unset_lock(locks(hash))
                        EXIT ! We have found a match and there can be only one single match
                     END IF
                  END DO
               END DO
            END DO
         END IF
         IF (my_rv) THEN
            DO ind = 1, 3
               NULLIFY (blocks_rv(ind)%block)
            END DO
            DEALLOCATE (blocks_rv)
         END IF
         IF (my_rxrv) THEN
            DO ind = 1, 3
               NULLIFY (blocks_rxrv(ind)%block)
            END DO
            DEALLOCATE (blocks_rxrv)
         END IF
         IF (my_rrv) THEN
            DO ind = 1, 6
               NULLIFY (blocks_rrv(ind)%block)
            END DO
            DEALLOCATE (blocks_rrv)
         END IF
         IF (my_rvr) THEN
            DO ind = 1, 6
               NULLIFY (blocks_rvr(ind)%block)
            END DO
            DEALLOCATE (blocks_rvr)
         END IF
         IF (my_rrv_vrr) THEN
            DO ind = 1, 6
               NULLIFY (blocks_rrv_vrr(ind)%block)
            END DO
            DEALLOCATE (blocks_rrv_vrr)
         END IF
         IF (my_r_rxvr) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  NULLIFY (blocks_r_rxvr(ind, ind2)%block)
               END DO
            END DO
            DEALLOCATE (blocks_r_rxvr)
         END IF
         IF (my_rxvr_r) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  NULLIFY (blocks_rxvr_r(ind, ind2)%block)
               END DO
            END DO
            DEALLOCATE (blocks_rxvr_r)
         END IF
         IF (my_r_doublecom) THEN
            DO ind = 1, 3
               DO ind2 = 1, 3
                  NULLIFY (blocks_r_doublecom(ind, ind2)%block)
               END DO
            END DO
            DEALLOCATE (blocks_r_doublecom)
         END IF
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_com_mom_nl

! **************************************************************************************************
!> \brief calculate \sum_R_ps (R_ps - R_nu) x [V_nl, r] summing over all pseudized atoms R
!> \param qs_kind_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param matrix_mag_nl ...
!> \param refpoint ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE build_com_nl_mag(qs_kind_set, sab_all, sap_ppnl, eps_ppnl, particle_set, matrix_mag_nl, refpoint, cell)

      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_mag_nl
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: refpoint
      TYPE(cell_type), INTENT(IN), OPTIONAL, POINTER     :: cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_com_nl_mag'

      INTEGER                                            :: handle, iab, iac, iatom, ibc, icol, &
                                                            ikind, ind, irow, jatom, jkind, kac, &
                                                            kbc, kkind, na, natom, nb, nkind, np, &
                                                            order, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, go, my_ref, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: r_b, r_ps, rab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: blocks_mag
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)

      my_ref = .FALSE.
      IF (PRESENT(refpoint)) THEN
         my_ref = .TRUE.
         CPASSERT(PRESENT(cell))
      END IF

      natom = SIZE(particle_set)
      nkind = SIZE(qs_kind_set)

      ! allocate integral storage
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO ind = 1, nkind*nkind
         NULLIFY (sap_int(ind)%alist, sap_int(ind)%asort, sap_int(ind)%aindex)
         sap_int(ind)%nalist = 0
      END DO

      ! build integrals over GTO + projector functions, refpoint actually
      order = 1   ! only need first moments (x, y, z)
      ! refpoint actually does not matter in this case, i. e. (order = 1 .and. commutator)
      IF (my_ref) THEN
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE., refpoint=refpoint, &
                             particle_set=particle_set, cell=cell)
      ELSE
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE.)
      END IF

      CALL sap_sort(sap_int)

      ! get access to basis sets
      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_mag_nl, sap_int, natom, nkind, eps_ppnl, locks, sab_all, &
!$OMP          particle_set, my_ref, refpoint) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, blocks_mag, r_ps, r_b, go, hash, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, kkind, kac, kbc)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sab_all(1)%nl_size
         ! get indices
         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_all(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         ! get blocks
         ALLOCATE (blocks_mag(3))
         DO ind = 1, 3
            CALL dbcsr_get_block_p(matrix_mag_nl(ind)%matrix, irow, icol, blocks_mag(ind)%block, found)
         END DO

         go = (ASSOCIATED(blocks_mag(1)%block) .AND. ASSOCIATED(blocks_mag(2)%block) .AND. ASSOCIATED(blocks_mag(3)%block))

         IF (go) THEN
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
                        ! Position of the pseudized atom
                        r_ps = particle_set(alist_ac%clist(kac)%catom)%r
                        r_b = refpoint

!$                      hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                      CALL omp_set_lock(locks(hash))
                        ! assemble integrals
                        IF (iatom <= jatom) THEN
                           blocks_mag(1)%block(1:na, 1:nb) = blocks_mag(1)%block(1:na, 1:nb) + &
                                              (r_ps(2) - r_b(2))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 4))) - &
                                                 MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 1)))) & !   R_y [V_nl, z]
                                            - (r_ps(3) - r_b(3))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 3))) - &
                                                 MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 1))))   ! - R_z [V_nl, y]
                           blocks_mag(2)%block(1:na, 1:nb) = blocks_mag(2)%block(1:na, 1:nb) + &
                                              (r_ps(3) - r_b(3))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 2))) - &
                                                 MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1)))) & !   R_z [V_nl, x]
                                            - (r_ps(1) - r_b(1))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 4))) - &
                                                 MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 1))))   ! - R_x [V_nl, z]
                           blocks_mag(3)%block(1:na, 1:nb) = blocks_mag(3)%block(1:na, 1:nb) + &
                                              (r_ps(1) - r_b(1))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 3))) - &
                                                MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 1)))) &  !   R_x [V_nl, y]
                                            - (r_ps(2) - r_b(2))*(MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 2))) - &
                                                MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1))))    ! - R_y [V_nl, x]
                        ELSE
                           blocks_mag(1)%block(1:nb, 1:na) = blocks_mag(1)%block(1:nb, 1:na) + &
                                              (r_ps(2) - r_b(2))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 4))) - &
                                                 MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 1)))) & !   R_y [V_nl, z]
                                            - (r_ps(3) - r_b(3))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 3))) - &
                                                 MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 1))))   ! - R_z [V_nl, y]
                           blocks_mag(2)%block(1:nb, 1:na) = blocks_mag(2)%block(1:nb, 1:na) + &
                                              (r_ps(3) - r_b(3))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 2))) - &
                                                 MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 1)))) & !   R_z [V_nl, x]
                                            - (r_ps(1) - r_b(1))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 4))) - &
                                                 MATMUL(bchint(1:nb, 1:np, 4), TRANSPOSE(acint(1:na, 1:np, 1))))   ! - R_x [V_nl, z]
                           blocks_mag(3)%block(1:nb, 1:na) = blocks_mag(3)%block(1:nb, 1:na) + &
                                              (r_ps(1) - r_b(1))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 3))) - &
                                                MATMUL(bchint(1:nb, 1:np, 3), TRANSPOSE(acint(1:na, 1:np, 1)))) &  !   R_x [V_nl, y]
                                            - (r_ps(2) - r_b(2))*(MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, 2))) - &
                                                MATMUL(bchint(1:nb, 1:np, 2), TRANSPOSE(acint(1:na, 1:np, 1))))    ! - R_y [V_nl, x]
                        END IF
!$                      CALL omp_unset_lock(locks(hash))
                        EXIT ! We have found a match and there can be only one single match
                     END IF
                  END DO
               END DO
            END DO
         END IF

         DO ind = 1, 3
            NULLIFY (blocks_mag(ind)%block)
         END DO
         DEALLOCATE (blocks_mag)
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      DEALLOCATE (basis_set)
      CALL release_sap_int(sap_int)

      CALL timestop(handle)

   END SUBROUTINE build_com_nl_mag

! **************************************************************************************************
!> \brief Calculate matrix_rv(gamma, delta) = < R^eta_gamma * Vnl * r_delta > for GIAOs
!> \param qs_kind_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param matrix_rv ...
!> \param ref_point ...
!> \param cell ...
!> \param direction_Or If set to true: calculate Vnl * r_delta
!>                     Otherwise       calculate r_delta * Vnl
! **************************************************************************************************
   SUBROUTINE build_com_vnl_giao(qs_kind_set, sab_all, sap_ppnl, eps_ppnl, particle_set, &
                                 matrix_rv, ref_point, cell, direction_Or)

      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL, POINTER                :: matrix_rv
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: ref_point
      TYPE(cell_type), INTENT(IN), OPTIONAL, POINTER     :: cell
      LOGICAL                                            :: direction_Or

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_com_vnl_giao'
      INTEGER, PARAMETER                                 :: i_1 = 1

      INTEGER                                            :: delta, gamma, handle, i, iab, iac, &
                                                            iatom, ibc, icol, ikind, irow, j, &
                                                            jatom, jkind, kac, kbc, kkind, na, &
                                                            natom, nb, nkind, np, order, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, my_ref, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rf
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:, :)   :: blocks_rv
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)

      my_ref = .FALSE.
      IF (PRESENT(ref_point)) THEN
         CPASSERT(PRESENT(cell)) ! need cell as well if refpoint is provided
         rf = ref_point
         my_ref = .TRUE.
      END IF

      nkind = SIZE(qs_kind_set)

      ! sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      order = 1
      IF (my_ref) THEN
         ! calculate integrals <a|x^n|p>
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE., refpoint=rf, &
                             particle_set=particle_set, cell=cell)
      ELSE
         CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, order, moment_mode=.TRUE.)
      END IF

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_all)
      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl, locks, sab_all, &
!$OMP          particle_set, direction_Or) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, blocks_rv, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, &
!$OMP          hash, natom, delta, gamma, achint, bchint, acint, bcint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_all(1)%nl_size

         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_all(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         irow = iatom
         icol = jatom

         ! allocate blocks
         ALLOCATE (blocks_rv(3, 3))

         ! get blocks
         DO i = 1, 3
            DO j = 1, 3
               CALL dbcsr_get_block_p(matrix_rv(i, j)%matrix, irow, icol, &
                                      blocks_rv(i, j)%block, found)
               blocks_rv(i, j)%block(:, :) = 0.0_dp
               CPASSERT(found)
            END DO
         END DO

         ! loop over all kinds for projector atom
         ! < iatom | katom > h < katom | jatom >
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
!$                   hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                   CALL omp_set_lock(locks(hash))

                     !! The atom index is alist_ac%clist(kac)%catom
                     ! The coordinate is particle_set(alist_ac%clist(kac)%catom)%r(:)
                     IF (direction_Or) THEN  ! V * r_delta * (R^eta_gamma - R^nu_gamma)
                        DO delta = 1, 3
                           DO gamma = 1, 3
                              blocks_rv(gamma, delta)%block(1:na, 1:nb) &
                                 = blocks_rv(gamma, delta)%block(1:na, 1:nb) + &
                                   MATMUL(achint(1:na, 1:np, i_1), TRANSPOSE(bcint(1:nb, 1:np, delta + 1))) &
                                   *(particle_set(alist_ac%clist(kac)%catom)%r(gamma) - particle_set(jatom)%r(gamma))
                           END DO
                        END DO
                     ELSE                   ! r_delta * V * (R^eta_gamma - R^nu_gamma)
                        DO delta = 1, 3
                           DO gamma = 1, 3
                              blocks_rv(gamma, delta)%block(1:na, 1:nb) &
                                 = blocks_rv(gamma, delta)%block(1:na, 1:nb) + &
                                   MATMUL(achint(1:na, 1:np, delta + 1), TRANSPOSE(bcint(1:nb, 1:np, i_1))) &
                                   *(particle_set(alist_ac%clist(kac)%catom)%r(gamma) - particle_set(jatom)%r(gamma))
                           END DO
                        END DO
                     END IF

!$                   CALL omp_unset_lock(locks(hash))
                     EXIT ! We have found a match and there can be only one single match
                  END IF
               END DO
            END DO
         END DO
         DO delta = 1, 3
            DO gamma = 1, 3
               NULLIFY (blocks_rv(gamma, delta)%block)
            END DO
         END DO
         DEALLOCATE (blocks_rv)
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_com_vnl_giao

END MODULE commutator_rpnl
