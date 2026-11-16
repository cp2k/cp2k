!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the construction of the coefficients
!>      for the expansion  of the atomic
!>      densities rho1_hard and rho1_soft in terms of primitive spherical gaussians.
!> \par History
!>      05-2004 created
!> \author MI
! **************************************************************************************************
MODULE qs_oce_methods

   USE ai_overlap,                      ONLY: overlap,&
                                              overlap_ab
   USE ao_util,                         ONLY: exp_radius
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset,&
                                              nso
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                              paw_proj_set_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: clist_type,&
                                              sap_int_type,&
                                              sap_sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_methods'

   ! Public subroutines

   PUBLIC :: build_oce_matrices, proj_blk, prj_scatter, prj_gather

CONTAINS

! **************************************************************************************************
!> \brief Overlap and derivative integrals <b|p_a> between orbital basis on B and PAW projectors on A
!> \param oces_t Transposed projector-orbital integrals (nsgf_cnt, nsatbas, maxder)
!> \param atom_ka Atomic kind of projector center A
!> \param atom_kb Atomic kind of orbital basis center B
!> \param rab Distance vector between centers A and B
!> \param nder Derivative order (0 <= nder <= 1)
!> \param sgf_list Indices of active spherical orbital functions on B
!> \param nsgf_cnt Number of active spherical orbital functions on B
!> \param sgf_soft_only True if all active functions on B are soft
!> \param eps_fit Threshold used to classify hard vs. soft functions
!> \param ovs_ws Spherical projector-primitive workspace (>= np_sph, >= maxco, maxder)
!> \param work One-dimensional scratch buffer for overlap_ab and DGEMM contractions
! **************************************************************************************************
   SUBROUTINE build_oce_block(oces_t, atom_ka, atom_kb, rab, nder, sgf_list, nsgf_cnt, &
                              sgf_soft_only, eps_fit, ovs_ws, work)

      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: oces_t
      TYPE(qs_kind_type), POINTER                        :: atom_ka, atom_kb
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      INTEGER, INTENT(IN)                                :: nder
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: sgf_list
      INTEGER, INTENT(OUT)                               :: nsgf_cnt
      LOGICAL, INTENT(OUT)                               :: sgf_soft_only
      REAL(KIND=dp), INTENT(IN)                          :: eps_fit
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: ovs_ws
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:), &
         INTENT(INOUT), TARGET                           :: work

      INTEGER :: ic, ider, igau, ip, ipgf, isgfb, isgfb_cnt, jc, jset, lprj, lsgfb, lshell, m0, &
         maxder, maxlprj, n, na_tot, nb_tot, nco_l, ncob, ncob_max, np_sph, nprj_l, nsatbas, &
         nsetb, nso_l, nso_tot, ofb, ovlp_cnt, rr_off, sgf_hard_only, sph_off
      INTEGER, DIMENSION(:), POINTER                     :: fp_spha, lb_max, lb_min, npgfb, nprjla, &
                                                            nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfb
      LOGICAL                                            :: paw_atom_a, paw_atom_b
      REAL(KIND=dp)                                      :: dab, hard_radius_a, hard_radius_b, radius
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_b, zisominb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: csprj, rpgfb, rzetprja, sab_p, sphi_b, &
                                                            zetb, zetprja
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dab_p
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: sph_p4
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_b
      TYPE(paw_proj_set_type), POINTER                   :: paw_proj_a, paw_proj_b

      CPASSERT(nder >= 0 .AND. nder <= 1)
      nsgf_cnt = 0
      sgf_soft_only = .TRUE.

      NULLIFY (paw_proj_a)
      CALL get_qs_kind(qs_kind=atom_ka, paw_proj_set=paw_proj_a, paw_atom=paw_atom_a, &
                       hard_radius=hard_radius_a)

      NULLIFY (orb_basis_b, paw_proj_b)
      CALL get_qs_kind(qs_kind=atom_kb, basis_set=orb_basis_b, &
                       paw_proj_set=paw_proj_b, paw_atom=paw_atom_b, &
                       hard_radius=hard_radius_b)

      IF (.NOT. paw_atom_a) RETURN

      NULLIFY (first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, sphi_b, set_radius_b, zetb, zisominb)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b, nset=nsetb, &
                             set_radius=set_radius_b, lmax=lb_max, lmin=lb_min, &
                             npgf=npgfb, nsgf_set=nsgfb, pgf_radius=rpgfb, &
                             sphi=sphi_b, zet=zetb, first_sgf=first_sgfb)

      !  Add the block ab
      dab = SQRT(SUM(rab*rab))
      IF (hard_radius_a + MAXVAL(set_radius_b(1:nsetb)) < dab) RETURN

      NULLIFY (nprjla, fp_spha, rzetprja, zetprja, csprj)
      CALL get_paw_proj_set(paw_proj_set=paw_proj_a, csprj=csprj, maxl=maxlprj, &
                            nprj=nprjla, nsgauprj=np_sph, nsatbas=nsatbas, &
                            first_prjs=fp_spha, rzetprj=rzetprja, zetprj=zetprja)
      IF (nsatbas == 0 .OR. np_sph == 0) RETURN
      IF (paw_atom_b) CALL get_paw_proj_set(paw_proj_set=paw_proj_b, zisomin=zisominb)

      maxder = ncoset(nder)
      isgfb_cnt = 1
      sgf_hard_only = 0
      DO jset = 1, nsetb
         !
         ! Set the contribution list
         IF (hard_radius_a + set_radius_b(jset) >= dab) THEN
            isgfb = first_sgfb(1, jset)
            lsgfb = isgfb - 1 + nsgfb(jset)
            DO jc = isgfb, lsgfb
               nsgf_cnt = nsgf_cnt + 1
               sgf_list(nsgf_cnt) = jc
            END DO

            ! check if this function is hard
            radius = exp_radius(lb_max(jset), MAXVAL(zetb(1:npgfb(jset), jset)), eps_fit, 1.0_dp)
            IF (radius <= hard_radius_b) sgf_hard_only = sgf_hard_only + 1

            ! Integral between proj of iatom and primitives of jatom
            ! Calculate the primitives overlap
            ncob_max = ncoset(lb_max(jset))
            ofb = ncoset(lb_min(jset) - 1)
            n = ncob_max - ofb
            ncob = npgfb(jset)*ncob_max
            nb_tot = npgfb(jset)*n
            IF (ofb > 0) ovs_ws(1:np_sph, 1:ncob, 1:maxder) = 0.0_dp

            DO lprj = 0, maxlprj
               nprj_l = nprjla(lprj)
               IF (nprj_l <= 0) CYCLE
               nco_l = nco(lprj)
               nso_l = nso(lprj)
               na_tot = nco_l*nprj_l
               nso_tot = nso_l*nprj_l
               m0 = fp_spha(lprj)
               ovlp_cnt = na_tot*nb_tot*maxder
               sph_off = ovlp_cnt + 1
               rr_off = ovlp_cnt + nso_tot*nb_tot*maxder + 1

               sab_p(1:na_tot, 1:nb_tot) => work(1:na_tot*nb_tot)
               IF (nder == 0) THEN
                  CALL overlap_ab(lprj, lprj, nprj_l, rzetprja(:, lprj), zetprja(:, lprj), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), &
                                  zetb(:, jset), -rab, sab=sab_p, rr_work=work(rr_off:))
               ELSE
                  dab_p(1:na_tot, 1:nb_tot, 1:3) => work(na_tot*nb_tot + 1:ovlp_cnt)
                  CALL overlap_ab(lprj, lprj, nprj_l, rzetprja(:, lprj), zetprja(:, lprj), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), &
                                  zetb(:, jset), -rab, sab=sab_p, dab=dab_p, rr_work=work(rr_off:))
               END IF

               CALL dgemm("N", "N", nso_l, nprj_l*nb_tot*maxder, nco_l, &
                          1.0_dp, orbtramat(lprj)%slm(1, 1), SIZE(orbtramat(lprj)%slm, 1), &
                          work(1), nco_l, &
                          0.0_dp, work(sph_off), nso_l)

               sph_p4(1:nso_tot, 1:n, 1:npgfb(jset), 1:maxder) => &
                  work(sph_off:sph_off + nso_tot*nb_tot*maxder - 1)
               DO ider = 1, maxder
                  DO ip = 1, npgfb(jset)
                     ic = (ip - 1)*ncob_max + ofb
                     ovs_ws(m0:m0 + nso_tot - 1, ic + 1:ic + n, ider) = sph_p4(1:nso_tot, 1:n, ip, ider)
                  END DO
               END DO
            END DO

            IF (paw_atom_b) THEN
               DO ipgf = 1, npgfb(jset)
                  DO lshell = lb_min(jset), lb_max(jset)
                     IF (zetb(ipgf, jset) >= zisominb(lshell)) THEN
                        igau = n*(ipgf - 1) + ncoset(lshell - 1)
                        ovs_ws(1:np_sph, igau + 1:igau + nco(lshell), 1:maxder) = 0.0_dp
                     END IF
                  END DO
               END DO
            END IF

            ! Contraction step (integrals and derivatives)
            DO ider = 1, maxder
               CALL dgemm("N", "N", np_sph, nsgfb(jset), ncob, &
                          1.0_dp, ovs_ws(1, 1, ider), SIZE(ovs_ws, 1), &
                          sphi_b(1, isgfb), SIZE(sphi_b, 1), &
                          0.0_dp, work(1), np_sph)
               CALL dgemm("T", "N", nsgfb(jset), nsatbas, np_sph, &
                          1.0_dp, work(1), np_sph, &
                          csprj(1, 1), SIZE(csprj, 1), &
                          0.0_dp, oces_t(isgfb_cnt, 1, ider), SIZE(oces_t, 1))
            END DO
            isgfb_cnt = isgfb_cnt + nsgfb(jset)
         END IF ! radius
      END DO !jset

      ! Check if the screened functions are all soft
      sgf_soft_only = .FALSE.
      IF (sgf_hard_only == 0) sgf_soft_only = .TRUE.

   END SUBROUTINE build_oce_block

! **************************************************************************************************
!> \brief ...
!> \param oceh ...
!> \param oces ...
!> \param atom_ka ...
!> \param sgf_list ...
!> \param nsgf_cnt ...
! **************************************************************************************************
   SUBROUTINE build_oce_block_local(oceh, oces, atom_ka, sgf_list, nsgf_cnt)

      TYPE(block_p_type), DIMENSION(:), POINTER          :: oceh, oces
      TYPE(qs_kind_type), POINTER                        :: atom_ka
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: sgf_list
      INTEGER, INTENT(OUT)                               :: nsgf_cnt

      INTEGER                                            :: i, iset, isgfa, j, jc, lsgfa, maxlprj, &
                                                            maxso1a, n, nsatbas, nset1a, nseta, &
                                                            nsgfa
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex, nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: paw_atom_a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: prjloc_h, prjloc_s
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_oce_h, local_oce_s
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_a, orb_basis_a
      TYPE(paw_proj_set_type), POINTER                   :: paw_proj_a

      NULLIFY (orb_basis_a, basis_1c_a, paw_proj_a)
      CALL get_qs_kind(qs_kind=atom_ka, basis_set=orb_basis_a, &
                       paw_proj_set=paw_proj_a, paw_atom=paw_atom_a)

      IF (.NOT. paw_atom_a) RETURN

      CALL get_paw_proj_set(paw_proj_set=paw_proj_a, maxl=maxlprj)
      CALL get_qs_kind(qs_kind=atom_ka, basis_set=basis_1c_a, basis_type="GAPW_1C")
      CALL get_gto_basis_set(gto_basis_set=basis_1c_a, nset=nset1a, maxso=maxso1a)
      NULLIFY (n2oindex)
      CALL get_paw_basis_info(basis_1c_a, n2oindex=n2oindex, nsatbas=nsatbas)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, first_sgf=first_sgfa, &
                             nsgf=nsgfa, nsgf_set=nsgf_seta, nset=nseta)

      NULLIFY (local_oce_h, local_oce_s)
      CALL get_paw_proj_set(paw_proj_set=paw_proj_a, &
                            local_oce_sphi_h=local_oce_h, &
                            local_oce_sphi_s=local_oce_s)

      ALLOCATE (prjloc_h(nset1a*maxso1a, nsgfa), prjloc_s(nset1a*maxso1a, nsgfa))
      prjloc_h = 0._dp
      prjloc_s = 0._dp

      nsgf_cnt = 0
      DO iset = 1, nseta
         isgfa = first_sgfa(1, iset)
         lsgfa = isgfa - 1 + nsgf_seta(iset)
         DO jc = isgfa, lsgfa
            nsgf_cnt = nsgf_cnt + 1
            sgf_list(nsgf_cnt) = jc
         END DO
         ! this asumes that the first sets are the same for basis_1c/orb_basis!
         n = maxso1a*(iset - 1)
         prjloc_h(n + 1:n + maxso1a, isgfa:lsgfa) = local_oce_h(1:maxso1a, isgfa:lsgfa)
         prjloc_s(n + 1:n + maxso1a, isgfa:lsgfa) = local_oce_s(1:maxso1a, isgfa:lsgfa)
      END DO

      DO i = 1, nsgfa
         DO j = 1, nsatbas
            jc = n2oindex(j)
            oceh(1)%block(j, i) = prjloc_h(jc, i)
            oces(1)%block(j, i) = prjloc_s(jc, i)
         END DO
      END DO

      DEALLOCATE (prjloc_h, prjloc_s)
      DEALLOCATE (n2oindex)

   END SUBROUTINE build_oce_block_local

! **************************************************************************************************
!> \brief ...
!> \param oceh ...
!> \param oces ...
!> \param atom_ka ...
!> \param sgf_list ...
!> \param nsgf_cnt ...
!> \param eps_fit ...
! **************************************************************************************************
   SUBROUTINE build_oce_block_1c(oceh, oces, atom_ka, sgf_list, nsgf_cnt, eps_fit)

      TYPE(block_p_type), DIMENSION(:), POINTER          :: oceh, oces
      TYPE(qs_kind_type), POINTER                        :: atom_ka
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: sgf_list
      INTEGER, INTENT(OUT)                               :: nsgf_cnt
      REAL(KIND=dp), INTENT(IN)                          :: eps_fit

      INTEGER :: first_col, ic, ig1, igau, ip, ipgf, isgfb, isgfb_cnt, jc, jset, lds, lm, lpoint, &
         lprj, lsgfb, lsgfb_cnt, lshell, m, m1, maxcob, maxlb, maxlprj, maxnprja, maxsoa, msab, n, &
         ncob, np_car, np_sph, nsatbas, nseta, nsetb, nsoatot, ntotsgfb
      INTEGER, DIMENSION(:), POINTER                     :: fp_cara, fp_spha, lb_max, lb_min, npgfb, &
                                                            nprjla, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfb
      LOGICAL                                            :: paw_atom_a
      REAL(KIND=dp)                                      :: radius, rc, rcprja
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ovh, ovs, spa_sb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: s
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_b, zisominb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chprj, csprj, rpgfb, rzetprja, spa_tmp, &
                                                            sphi_b, zetb, zetprja
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_a, orb_basis_b
      TYPE(paw_proj_set_type), POINTER                   :: paw_proj_a

      NULLIFY (orb_basis_b, basis_1c_a, paw_proj_a)
      CALL get_qs_kind(qs_kind=atom_ka, paw_atom=paw_atom_a)

      IF (.NOT. paw_atom_a) RETURN

      CALL get_qs_kind(qs_kind=atom_ka, paw_proj_set=paw_proj_a)
      NULLIFY (nprjla, fp_cara, fp_spha, rzetprja, zetprja)
      CALL get_paw_proj_set(paw_proj_set=paw_proj_a, csprj=csprj, chprj=chprj, maxl=maxlprj, &
                            nprj=nprjla, ncgauprj=np_car, nsgauprj=np_sph, nsatbas=nsatbas, rcprj=rcprja, &
                            first_prj=fp_cara, first_prjs=fp_spha, &
                            rzetprj=rzetprja, zetprj=zetprja)

      CALL get_qs_kind(qs_kind=atom_ka, basis_set=orb_basis_b)
      NULLIFY (first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, sphi_b, set_radius_b, zetb, zisominb)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b, nset=nsetb, nsgf=ntotsgfb, &
                             set_radius=set_radius_b, lmax=lb_max, lmin=lb_min, &
                             npgf=npgfb, nsgf_set=nsgfb, pgf_radius=rpgfb, &
                             sphi=sphi_b, zet=zetb, first_sgf=first_sgfb, &
                             maxco=maxcob, maxl=maxlb)

      CALL get_qs_kind(qs_kind=atom_ka, hard_radius=rc)

      CALL get_qs_kind(qs_kind=atom_ka, basis_set=basis_1c_a, basis_type="GAPW_1C")
      CALL get_gto_basis_set(gto_basis_set=basis_1c_a, nset=nseta, maxso=maxsoa)

      !  Add the block ab
      nsoatot = maxsoa*nseta
      maxnprja = SIZE(zetprja, 1)

      lm = MAX(maxlb, maxlprj)
      lds = ncoset(lm + 1)
      msab = MAX(maxnprja*ncoset(maxlprj), maxcob)

      ALLOCATE (s(lds, lds, 1))
      ALLOCATE (spa_sb(np_car, ntotsgfb))
      ALLOCATE (spa_tmp(msab, msab))
      ALLOCATE (ovs(np_sph, maxcob*nsetb))
      ALLOCATE (ovh(np_sph, maxcob*nsetb))

      m1 = 0
      nsgf_cnt = 0
      isgfb_cnt = 1
      DO jset = 1, nsetb
         !
         ! Set the contribution list
         isgfb = first_sgfb(1, jset)
         lsgfb = isgfb - 1 + nsgfb(jset)
         DO jc = isgfb, lsgfb
            nsgf_cnt = nsgf_cnt + 1
            sgf_list(nsgf_cnt) = jc
         END DO

         ! Integral between proj of iatom and primitives of iatom
         ! Calculate the primitives overlap
         spa_tmp = 0.0_dp
         ovs = 0.0_dp
         ovh = 0.0_dp
         s = 0.0_dp
         ncob = npgfb(jset)*ncoset(lb_max(jset))
         isgfb = first_sgfb(1, jset)
         lsgfb = isgfb - 1 + nsgfb(jset)

         lsgfb_cnt = isgfb_cnt - 1 + nsgfb(jset)

         DO lprj = 0, maxlprj
            CALL overlap(lprj, lprj, nprjla(lprj), &
                         rzetprja(:, lprj), zetprja(:, lprj), &
                         lb_max(jset), lb_min(jset), npgfb(jset), &
                         rpgfb(:, jset), zetb(:, jset), &
                         -[0._dp, 0._dp, 0._dp], 0.0_dp, spa_tmp, &
                         0, .TRUE., s, lds)
            DO ipgf = 1, nprjla(lprj)
               lpoint = ncoset(lprj - 1) + 1 + (ipgf - 1)*ncoset(lprj)
               m = fp_spha(lprj) + (ipgf - 1)*nso(lprj)
               DO ip = 1, npgfb(jset)
                  ic = (ip - 1)*ncoset(lb_max(jset))
                  igau = ic + m1 + ncoset(lb_min(jset) - 1) + 1
                  ig1 = ic + ncoset(lb_min(jset) - 1) + 1
                  n = ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1)
                  ovs(m:m + nso(lprj) - 1, igau:igau + n - 1) = &
                     MATMUL(orbtramat(lprj)%slm(1:nso(lprj), 1:nco(lprj)), &
                            spa_tmp(lpoint:lpoint + nco(lprj) - 1, ig1:ig1 + n - 1))
               END DO
            END DO
         END DO

         ovh(:, :) = ovs(:, :)

         CALL get_paw_proj_set(paw_proj_set=paw_proj_a, zisomin=zisominb)
         DO ipgf = 1, npgfb(jset)
            DO lshell = lb_min(jset), lb_max(jset)
               radius = exp_radius(lshell, zetb(ipgf, jset), eps_fit, 1.0_dp)
               IF (radius < rc) THEN
                  n = ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1)
                  igau = n*(ipgf - 1) + ncoset(lshell - 1)
                  ovs(:, igau + 1 + m1:igau + nco(lshell) + m1) = 0.0_dp
               END IF
            END DO
         END DO

         ! Contraction step (integrals and derivatives)
         first_col = 1 + m1
         spa_sb(1:np_sph, isgfb:isgfb + nsgfb(jset) - 1) = &
            MATMUL(ovs(1:np_sph, first_col:first_col + ncob - 1), &
                   sphi_b(1:ncob, isgfb:isgfb + nsgfb(jset) - 1))

         oces(1)%block(1:nsatbas, isgfb_cnt:isgfb_cnt + nsgfb(jset) - 1) = &
            oces(1)%block(1:nsatbas, isgfb_cnt:isgfb_cnt + nsgfb(jset) - 1) + &
            MATMUL(TRANSPOSE(csprj(1:np_sph, 1:nsatbas)), &
                   spa_sb(1:np_sph, isgfb:isgfb + nsgfb(jset) - 1))

         spa_sb(1:np_sph, isgfb:isgfb + nsgfb(jset) - 1) = &
            MATMUL(ovh(1:np_sph, first_col:first_col + ncob - 1), &
                   sphi_b(1:ncob, isgfb:isgfb + nsgfb(jset) - 1))

         oceh(1)%block(1:nsatbas, isgfb_cnt:isgfb_cnt + nsgfb(jset) - 1) = &
            oceh(1)%block(1:nsatbas, isgfb_cnt:isgfb_cnt + nsgfb(jset) - 1) + &
            MATMUL(TRANSPOSE(chprj(1:np_sph, 1:nsatbas)), &
                   spa_sb(1:np_sph, isgfb:isgfb + nsgfb(jset) - 1))

         isgfb_cnt = isgfb_cnt + nsgfb(jset)
         m1 = m1 + maxcob
      END DO !jset

      DEALLOCATE (s, spa_sb, spa_tmp, ovs)

   END SUBROUTINE build_oce_block_1c

! **************************************************************************************************
!> \brief Set up the sparse matrix for the coefficients of one center expansions
!>      This routine uses the same logic as the nonlocal pseudopotential
!> \param intac TYPE that holds the integrals (a=basis; c=projector)
!> \param calculate_forces ...
!> \param nder Derivative order (0 <= nder <= 1)
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param sap_oce ...
!> \param eps_fit ...
!> \par History
!>      02.2009 created
!> \author jgh
! **************************************************************************************************
   SUBROUTINE build_oce_matrices(intac, calculate_forces, nder, &
                                 qs_kind_set, particle_set, sap_oce, eps_fit)

      TYPE(sap_int_type), DIMENSION(:), POINTER          :: intac
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER                                            :: nder
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sap_oce
      REAL(KIND=dp), INTENT(IN)                          :: eps_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_oce_matrices'

      INTEGER :: atom_a, atom_b, handle, i, iac, ikind, ilist, jkind, jneighbor, lwork, &
         max_nsatbas, maxco, maxder, maxl, maxlgto, maxlprj, maxprj, maxsgf, natom, nkind, nlist, &
         nneighbor, nsgf_cnt, nsgfa, nsobtot, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsatbas_kind, nsgfa_kind, sgf_list
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: local, sgf_soft_only
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: has_1c_basis, has_orb_basis, paw_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oces_t, ovs_ws
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(block_p_type), DIMENSION(:), POINTER          :: oceh, oces
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_paw, orb_basis_set
      TYPE(paw_proj_set_type), POINTER                   :: paw_proj_b
      TYPE(qs_kind_type), POINTER                        :: at_a, at_b

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      END IF

      IF (ASSOCIATED(sap_oce)) THEN

         nkind = SIZE(qs_kind_set)
         natom = SIZE(particle_set)

         CPASSERT(nder >= 0 .AND. nder <= 1)
         maxder = ncoset(nder)

         CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                              maxco=maxco, &
                              maxlgto=maxlgto, &
                              maxlprj=maxlprj, &
                              maxco_proj=maxprj, &
                              maxsgf=maxsgf)

         maxl = MAX(maxlgto, maxlprj)
         CALL init_orbital_pointers(maxl + nder + 1)

         DO i = 1, nkind*nkind
            NULLIFY (intac(i)%alist, intac(i)%asort, intac(i)%aindex)
            intac(i)%nalist = 0
         END DO

         ALLOCATE (has_orb_basis(nkind), paw_kind(nkind), has_1c_basis(nkind), &
                   nsgfa_kind(nkind), nsatbas_kind(nkind))
         nsgfa_kind(:) = 0
         nsatbas_kind(:) = 0
         DO ikind = 1, nkind
            NULLIFY (orb_basis_set, paw_proj_b, orb_basis_paw)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=orb_basis_set, &
                             paw_proj_set=paw_proj_b, paw_atom=paw_kind(ikind))
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=orb_basis_paw, basis_type="GAPW_1C")
            has_orb_basis(ikind) = ASSOCIATED(orb_basis_set)
            has_1c_basis(ikind) = ASSOCIATED(orb_basis_paw)
            IF (has_orb_basis(ikind)) CALL get_gto_basis_set(orb_basis_set, nsgf=nsgfa_kind(ikind))
            IF (paw_kind(ikind) .AND. ASSOCIATED(paw_proj_b)) THEN
               CALL get_paw_proj_set(paw_proj_b, nsatbas=nsatbas_kind(ikind))
            END IF
         END DO
         max_nsatbas = MAXVAL(nsatbas_kind)

         ! Allocate memory list
         DO slot = 1, sap_oce(1)%nl_size
            ikind = sap_oce(1)%nlist_task(slot)%ikind
            jkind = sap_oce(1)%nlist_task(slot)%jkind
            atom_a = sap_oce(1)%nlist_task(slot)%iatom
            nlist = sap_oce(1)%nlist_task(slot)%nlist
            ilist = sap_oce(1)%nlist_task(slot)%ilist
            nneighbor = sap_oce(1)%nlist_task(slot)%nnode

            iac = ikind + nkind*(jkind - 1)

            IF (.NOT. has_orb_basis(ikind)) CYCLE
            IF (.NOT. paw_kind(jkind)) CYCLE
            IF (.NOT. has_1c_basis(jkind)) CYCLE
            IF (.NOT. ASSOCIATED(intac(iac)%alist)) THEN
               intac(iac)%a_kind = ikind
               intac(iac)%p_kind = jkind
               intac(iac)%nalist = nlist
               ALLOCATE (intac(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (intac(iac)%alist(i)%clist)
                  intac(iac)%alist(i)%aatom = 0
                  intac(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(intac(iac)%alist(ilist)%clist)) THEN
               intac(iac)%alist(ilist)%aatom = atom_a
               intac(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (intac(iac)%alist(ilist)%clist(nneighbor))
            END IF
         END DO

         lwork = MAX(1, maxprj*maxsgf, 2*maxprj*maxco*maxder + 3*(maxl + nder + 1)**2)

         !calculate the overlap integrals <a|p>
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (intac, nder, nkind, maxder, sap_oce, qs_kind_set, eps_fit, &
!$OMP         has_orb_basis, paw_kind, has_1c_basis, nsgfa_kind, nsatbas_kind, &
!$OMP         maxco, maxprj, maxsgf, max_nsatbas, lwork) &
!$OMP PRIVATE (oceh, oces, oces_t, ovs_ws, work, sgf_list, slot, ikind, jkind, atom_a, atom_b, &
!$OMP          ilist, jneighbor, rab, cell_b, iac, nsgfa, nsobtot, clist, at_a, at_b, local, i, &
!$OMP          sgf_soft_only, nsgf_cnt)

         ALLOCATE (oces_t(MAX(1, maxsgf), MAX(1, max_nsatbas), maxder), &
                   ovs_ws(MAX(1, maxprj), MAX(1, maxco), maxder), work(lwork), &
                   oceh(maxder), oces(maxder), sgf_list(maxsgf))

!$OMP DO SCHEDULE(GUIDED)
         DO slot = 1, sap_oce(1)%nl_size
            ikind = sap_oce(1)%nlist_task(slot)%ikind
            jkind = sap_oce(1)%nlist_task(slot)%jkind
            atom_a = sap_oce(1)%nlist_task(slot)%iatom
            atom_b = sap_oce(1)%nlist_task(slot)%jatom
            ilist = sap_oce(1)%nlist_task(slot)%ilist
            jneighbor = sap_oce(1)%nlist_task(slot)%inode
            rab(1:3) = sap_oce(1)%nlist_task(slot)%r(1:3)
            cell_b(1:3) = sap_oce(1)%nlist_task(slot)%cell(1:3)

            iac = ikind + nkind*(jkind - 1)

            IF (.NOT. has_orb_basis(ikind)) CYCLE
            IF (.NOT. paw_kind(jkind)) CYCLE
            IF (.NOT. has_1c_basis(jkind)) CYCLE

            nsgfa = nsgfa_kind(ikind)
            nsobtot = nsatbas_kind(jkind)

            clist => intac(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = atom_b
            clist%cell = cell_b
            clist%rac = rab
            clist%nsgf_cnt = 0
            clist%maxac = 0.0_dp
            clist%maxach = 0.0_dp
            NULLIFY (clist%acint, clist%achint, clist%sgf_list)

            at_a => qs_kind_set(jkind)
            at_b => qs_kind_set(ikind)

            local = (atom_a == atom_b .AND. ALL(cell_b == 0))

            IF (local) THEN
               DO i = 1, maxder
                  ALLOCATE (oceh(i)%block(nsobtot, nsgfa), oces(i)%block(nsobtot, nsgfa))
                  oceh(i)%block = 0._dp
                  oces(i)%block = 0._dp
               END DO
               CALL build_oce_block_local(oceh, oces, at_a, sgf_list, nsgf_cnt)
               clist%nsgf_cnt = nsgf_cnt
               clist%sgf_soft_only = .FALSE.
               IF (nsgf_cnt > 0) THEN
                  ALLOCATE (clist%acint(nsgf_cnt, nsobtot, maxder), clist%sgf_list(nsgf_cnt))
                  clist%acint(:, :, :) = 0._dp
                  clist%sgf_list(:) = HUGE(0)
                  CPASSERT(nsgf_cnt == nsgfa)
                  ! *** Special case: A=B
                  ALLOCATE (clist%achint(nsgfa, nsobtot, maxder))
                  clist%achint = 0._dp
                  clist%acint(1:nsgfa, 1:nsobtot, 1) = TRANSPOSE(oces(1)%block(1:nsobtot, 1:nsgfa))
                  clist%achint(1:nsgfa, 1:nsobtot, 1) = TRANSPOSE(oceh(1)%block(1:nsobtot, 1:nsgfa))
                  clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
                  clist%maxach = 0._dp
                  clist%sgf_list(1:nsgf_cnt) = sgf_list(1:nsgf_cnt)
               END IF
               DO i = 1, maxder
                  DEALLOCATE (oceh(i)%block, oces(i)%block)
               END DO
            ELSE
               CALL build_oce_block(oces_t, at_a, at_b, rab, nder, sgf_list, nsgf_cnt, &
                                    sgf_soft_only, eps_fit, ovs_ws, work)
               clist%nsgf_cnt = nsgf_cnt
               clist%sgf_soft_only = sgf_soft_only
               IF (nsgf_cnt > 0) THEN
                  ALLOCATE (clist%acint(nsgf_cnt, nsobtot, maxder), clist%sgf_list(nsgf_cnt))
                  clist%acint(:, :, :) = oces_t(1:nsgf_cnt, 1:nsobtot, 1:maxder)
                  clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
                  clist%maxach = 0._dp
                  clist%sgf_list(1:nsgf_cnt) = sgf_list(1:nsgf_cnt)
               END IF
            END IF

         END DO

         DEALLOCATE (oces_t, ovs_ws, work, sgf_list, oceh, oces)
!$OMP END PARALLEL

         ! Setup sort index
         CALL sap_sort(intac)

         DEALLOCATE (has_orb_basis, paw_kind, has_1c_basis, nsgfa_kind, nsatbas_kind)

      END IF

      CALL timestop(handle)

   END SUBROUTINE build_oce_matrices

! **************************************************************************************************
!> \brief Project a matrix block onto the local atomic functions.
!>
!> \param h_a ...
!> \param s_a ...
!> \param na ...
!> \param h_b ...
!> \param s_b ...
!> \param nb ...
!> \param blk ...
!> \param ldb ...
!> \param proj_h ...
!> \param proj_s ...
!> \param nso ...
!> \param len1 ...
!> \param len2 ...
!> \param fac ...
!> \param distab ...
!> \param work1 ...
!> \param work2 ...
!> \par History
!>      02.2009 created
!>      09.2016 use automatic arrays [M Tucker]
!>      05.2026 move work arrays to callers
! **************************************************************************************************
   SUBROUTINE proj_blk(h_a, s_a, na, h_b, s_b, nb, blk, ldb, proj_h, proj_s, nso, len1, len2, fac, distab, &
                       work1, work2)

      INTEGER, INTENT(IN)                                :: na
      REAL(KIND=dp), INTENT(IN)                          :: s_a(na, *), h_a(na, *)
      INTEGER, INTENT(IN)                                :: nb
      REAL(KIND=dp), INTENT(IN)                          :: s_b(nb, *), h_b(nb, *)
      INTEGER, INTENT(IN)                                :: ldb
      REAL(KIND=dp), INTENT(IN)                          :: blk(ldb, *)
      INTEGER, INTENT(IN)                                :: nso
      REAL(KIND=dp), INTENT(INOUT)                       :: proj_s(nso, *), proj_h(nso, *)
      INTEGER, INTENT(IN)                                :: len1, len2
      REAL(KIND=dp), INTENT(IN)                          :: fac
      LOGICAL, INTENT(IN)                                :: distab
      REAL(KIND=dp), INTENT(INOUT)                       :: work1(len1), work2(len2)

      IF (na == 0 .OR. nb == 0 .OR. nso == 0) RETURN

      ! Handle special cases
      IF (na == 1 .AND. nb == 1) THEN
         ! hard
         CALL dger(nso, nso, fac*blk(1, 1), h_a(1, 1), 1, h_b(1, 1), 1, proj_h(1, 1), nso)
         ! soft
         CALL dger(nso, nso, fac*blk(1, 1), s_a(1, 1), 1, s_b(1, 1), 1, proj_s(1, 1), nso)
      ELSE
         IF (distab) THEN
            ! hard
            CALL dgemm('N', 'N', na, nso, nb, fac, blk(1, 1), ldb, h_b(1, 1), nb, 0.0_dp, work1(1), na)
            CALL dgemm('T', 'N', nso, nso, na, 1.0_dp, h_a(1, 1), na, work1(1), na, 0.0_dp, work2(1), nso)
            CALL daxpy(nso*nso, 1.0_dp, work2(1), 1, proj_h(1, 1), 1)
            ! soft
            CALL daxpy(nso*nso, 1.0_dp, work2(1), 1, proj_s(1, 1), 1)
         ELSE
            ! hard
            CALL dgemm('N', 'N', na, nso, nb, fac, blk(1, 1), ldb, h_b(1, 1), nb, 0.0_dp, work1(1), na)
            CALL dgemm('T', 'N', nso, nso, na, 1.0_dp, h_a(1, 1), na, work1(1), na, 1.0_dp, proj_h(1, 1), nso)
            ! soft
            CALL dgemm('N', 'N', na, nso, nb, fac, blk(1, 1), ldb, s_b(1, 1), nb, 0.0_dp, work1(1), na)
            CALL dgemm('T', 'N', nso, nso, na, 1.0_dp, s_a(1, 1), na, work1(1), na, 1.0_dp, proj_s(1, 1), nso)
         END IF
      END IF

   END SUBROUTINE proj_blk

! **************************************************************************************************
!> \brief ...
!> \param ain matrix in old indexing
!> \param aout matrix in new compressed indexing
!> \param atom ...
! **************************************************************************************************
   SUBROUTINE prj_gather(ain, aout, atom)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ain
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aout
      TYPE(qs_kind_type), INTENT(IN)                     :: atom

      INTEGER                                            :: i, ip, j, jp, nbas
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c

      NULLIFY (basis_1c)
      CALL get_qs_kind(qs_kind=atom, basis_set=basis_1c, basis_type="GAPW_1C")
      NULLIFY (n2oindex)
      CALL get_paw_basis_info(basis_1c, n2oindex=n2oindex, nsatbas=nbas)

      DO i = 1, nbas
         ip = n2oindex(i)
         DO j = 1, nbas
            jp = n2oindex(j)
            aout(j, i) = ain(jp, ip)
         END DO
      END DO

      DEALLOCATE (n2oindex)

   END SUBROUTINE prj_gather

! **************************************************************************************************
!> \brief ...
!> \param ain  matrix in new compressed indexing
!> \param aout matrix in old indexing (addup)
!> \param atom ...
! **************************************************************************************************
   SUBROUTINE prj_scatter(ain, aout, atom)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ain
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aout
      TYPE(qs_kind_type), INTENT(IN)                     :: atom

      INTEGER                                            :: i, ip, j, jp, nbas
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c

      NULLIFY (basis_1c)
      CALL get_qs_kind(qs_kind=atom, basis_set=basis_1c, basis_type="GAPW_1C")
      NULLIFY (n2oindex)
      CALL get_paw_basis_info(basis_1c, n2oindex=n2oindex, nsatbas=nbas)

      DO i = 1, nbas
         ip = n2oindex(i)
         DO j = 1, nbas
            jp = n2oindex(j)
            aout(jp, ip) = aout(jp, ip) + ain(j, i)
         END DO
      END DO

      DEALLOCATE (n2oindex)

   END SUBROUTINE prj_scatter

END MODULE qs_oce_methods
