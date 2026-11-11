!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that build the Kohn-Sham matrix  contributions
!>      coming from local atomic densities
! **************************************************************************************************
MODULE qs_ks_atom

   USE ao_util,                         ONLY: trace_r_AxB
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              neighbor_list_task_type
   USE qs_oce_methods,                  ONLY: prj_gather
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE sap_kind_types,                  ONLY: alist_post_align_blk,&
                                              alist_pre_align_blk,&
                                              alist_type,&
                                              clist_type,&
                                              get_alist
   USE util,                            ONLY: get_limit
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_atom'
   INTEGER, PARAMETER, PRIVATE :: l2_chunk_dp = 16000

   PUBLIC :: update_ks_atom

CONTAINS

! **************************************************************************************************
!> \brief The correction to the KS matrix due to the GAPW local terms to the hartree and
!>      XC contributions is here added. The correspondig forces contribution are also calculated
!>      if required. Each sparse-matrix block A-B is corrected by all the atomic contributions
!>      centered on atoms C for which the triplet A-C-B exists (they are close enough)
!>      To this end special lists are used
!> \param qs_env qs environment, for the lists, the contraction coefficients and the
!>               pre calculated integrals of the potential with the atomic orbitals
!> \param ksmat KS matrix, sparse matrix
!> \param pmat density matrix, sparse matrix, needed only for the forces
!> \param forces switch for the calculation of the forces on atoms
!> \param tddft switch for TDDFT linear response
!> \param rho_atom_external ...
!> \param kind_set_external ...
!> \param oce_external ...
!> \param sab_external ...
!> \param kscale ...
!> \param kintegral ...
!> \param kforce ...
!> \param fscale ...
!> \par History
!>      created [MI]
!>      the loop over the spins is done internally [03-05,MI]
!>      Rewrite using new OCE matrices [08.02.09, jhu]
!>      Add OpenMP [Apr 2016, EPCC]
!>      Allow for external kind_set, rho_atom_set, oce, sab 12.2019 (A. Bussy)
! **************************************************************************************************
   SUBROUTINE update_ks_atom(qs_env, ksmat, pmat, forces, tddft, rho_atom_external, &
                             kind_set_external, oce_external, sab_external, kscale, &
                             kintegral, kforce, fscale)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(*), INTENT(INOUT)    :: ksmat, pmat
      LOGICAL, INTENT(IN)                                :: forces
      LOGICAL, INTENT(IN), OPTIONAL                      :: tddft
      TYPE(rho_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_atom_external
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external
      TYPE(oce_matrix_type), OPTIONAL, POINTER           :: oce_external
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_external
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kscale, kintegral, kforce
      REAL(KIND=dp), DIMENSION(2), INTENT(IN), OPTIONAL  :: fscale

      CHARACTER(len=*), PARAMETER                        :: routineN = 'update_ks_atom'

      INTEGER :: bo(2), dir, handle, i, ia_kind, iac, iat, iatom, ibc, ikind, img, inb, ip, ispin, &
         ja_kind, jatom, jkind, k_total, ka_kind, kac, kac_hi, kac_lo, katom, kbc, kkind, &
         ku_alloc, max_gau, max_nb, max_nclist, max_nsgf, n_cont_a, n_cont_b, n_nb, nat, natom, &
         nimages, nkind, nl_size, nsoctot, nspins, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, ia_ptr, ia_work, k_offset, &
                                                            slot_order
      INTEGER, DIMENSION(3)                              :: cell_b
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, list_a, list_b
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: dista, do_forces, found, my_tddft, &
                                                            paw_atom, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: has_intac, paw_kind
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: u_touched
      REAL(KIND=dp)                                      :: beta, factor1, factor2, ieqj, trace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: C_int_h, C_int_s
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a_matrix, coc, f_katom, k_buf, p_matrix, &
                                                            u_buf
      REAL(KIND=dp), DIMENSION(2)                        :: force_fac
      REAL(KIND=dp), DIMENSION(3)                        :: f_jatom, fb_triple, rac, rbc
      REAL(KIND=dp), DIMENSION(3, 3)                     :: force_tmp, pv_virial_thread
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: C_coeff_hh_a, C_coeff_hh_b, &
                                                            C_coeff_ss_a, C_coeff_ss_b
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(clist_type), POINTER                          :: clist_a, clist_b
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:, :) :: nb_hblk, nb_pblk
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: my_sab
      TYPE(neighbor_list_task_type), POINTER             :: task
      TYPE(oce_matrix_type), POINTER                     :: my_oce
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: my_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_local_h, int_local_s
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: my_rho_atom
      TYPE(virial_type), POINTER                         :: virial

!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num

      CALL timeset(routineN, handle)

      NULLIFY (my_kind_set, atomic_kind_set, force, my_oce, para_env, my_rho_atom, my_sab, dft_control)

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=my_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      force=force, &
                      oce=my_oce, &
                      para_env=para_env, &
                      rho_atom_set=my_rho_atom, &
                      virial=virial, &
                      sab_orb=my_sab, &
                      dft_control=dft_control)

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      factor1 = 1.0_dp
      factor2 = 1.0_dp

      !deal with externals
      my_tddft = .FALSE.
      IF (PRESENT(tddft)) my_tddft = tddft
      IF (my_tddft) THEN
         IF (nspins == 1) factor1 = 2.0_dp
         CPASSERT(nimages == 1)
      END IF
      IF (PRESENT(kscale)) THEN
         factor1 = factor1*kscale
         factor2 = factor2*kscale
      END IF
      IF (PRESENT(kintegral)) factor1 = kintegral
      IF (PRESENT(kforce)) factor2 = kforce
      force_fac = 1.0_dp
      IF (PRESENT(fscale)) force_fac(:) = fscale(:)

      IF (PRESENT(rho_atom_external)) my_rho_atom => rho_atom_external
      IF (PRESENT(kind_set_external)) my_kind_set => kind_set_external
      IF (PRESENT(oce_external)) my_oce => oce_external
      IF (PRESENT(sab_external)) my_sab => sab_external

      ! kpoint images
      NULLIFY (cell_to_index)
      IF (nimages > 1) THEN
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(my_kind_set, maxsgf=max_nsgf, maxgtops=max_gau, basis_type="GAPW_1C")

      IF (forces) THEN
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      ELSE
         use_virial = .FALSE.
      END IF

      pv_virial_thread(:, :) = 0.0_dp ! always initialize to avoid floating point exceptions in OMP REDUCTION

      nkind = SIZE(my_kind_set, 1)
      ! Collect the local potential contributions from all the processors
      ASSOCIATE (mepos => para_env%mepos, num_pe => para_env%num_pe)
      DO ikind = 1, nkind
         NULLIFY (atom_list)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom)
         IF (paw_atom) THEN
            ! gather the atomic block integrals in a more compressed format
            bo = get_limit(nat, num_pe, mepos)
            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               DO ispin = 1, nspins
                  CALL prj_gather(my_rho_atom(iatom)%ga_Vlocal_gb_h(ispin)%r_coef, &
                                  my_rho_atom(iatom)%int_scr_h(ispin)%r_coef, my_kind_set(ikind))
                  CALL prj_gather(my_rho_atom(iatom)%ga_Vlocal_gb_s(ispin)%r_coef, &
                                  my_rho_atom(iatom)%int_scr_s(ispin)%r_coef, my_kind_set(ikind))
               END DO
            END DO
            ! broadcast the CPC arrays to all processors (replicated data)
            DO ip = 0, num_pe - 1
               bo = get_limit(nat, num_pe, ip)
               DO iat = bo(1), bo(2)
                  iatom = atom_list(iat)
                  DO ispin = 1, nspins
                     CALL para_env%bcast(my_rho_atom(iatom)%int_scr_h(ispin)%r_coef, ip)
                     CALL para_env%bcast(my_rho_atom(iatom)%int_scr_s(ispin)%r_coef, ip)
                  END DO
               END DO
            END DO
         END IF
      END DO
      END ASSOCIATE

      ALLOCATE (basis_set_list(nkind))
      ALLOCATE (paw_kind(nkind), has_intac(nkind*nkind))
      paw_kind(:) = .FALSE.
      has_intac(:) = .FALSE.
      DO ikind = 1, nkind
         CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_set_a, paw_atom=paw_kind(ikind))
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      DO ikind = 1, nkind*nkind
         has_intac(ikind) = ASSOCIATED(my_oce%intac(ikind)%alist)
      END DO

      ku_alloc = MAX(1, l2_chunk_dp/nspins, max_nsgf*max_gau) ! 16k-double L2 chunk across spins
      max_nclist = 1
      DO iac = 1, nkind*nkind
         IF (.NOT. has_intac(iac)) CYCLE
         DO i = 1, SIZE(my_oce%intac(iac)%alist)
            max_nclist = MAX(max_nclist, my_oce%intac(iac)%alist(i)%nclist)
         END DO
      END DO

      CPASSERT(ASSOCIATED(my_sab(1)%nlist_task))
      nl_size = my_sab(1)%nl_size
      ALLOCATE (ia_ptr(natom + 1), ia_work(natom), slot_order(MAX(1, nl_size)))
      ia_work(:) = 0
      DO i = 1, 2
         DO slot = 1, nl_size
            task => my_sab(1)%nlist_task(slot)
            IF (.NOT. ASSOCIATED(basis_set_list(task%ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set_list(task%jkind)%gto_basis_set)) CYCLE
            iatom = MIN(task%iatom, task%jatom)
            IF (i == 2) slot_order(ia_work(iatom)) = slot
            ia_work(iatom) = ia_work(iatom) + 1
         END DO
         IF (i == 2) EXIT
         ia_ptr(1) = 1
         max_nb = MAX(1, MAXVAL(ia_work))
         DO iatom = 1, natom
            ia_ptr(iatom + 1) = ia_ptr(iatom) + ia_work(iatom)
            ia_work(iatom) = ia_ptr(iatom)
         END DO
      END DO
      DEALLOCATE (ia_work)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          SHARED(natom, nkind, nspins, nimages, cell_to_index, ksmat, pmat, forces, max_gau, &
!$OMP                 max_nsgf, max_nb, max_nclist, ku_alloc, ia_ptr, slot_order, my_sab, paw_kind, &
!$OMP                 has_intac, my_oce, my_rho_atom, factor1, factor2, force_fac, use_virial, &
!$OMP                 atom_of_kind, force, locks) &
!$OMP          PRIVATE(nb_hblk, nb_pblk, C_int_h, C_int_s, coc, a_matrix, p_matrix, k_buf, k_offset, &
!$OMP                  u_buf, u_touched, f_katom, f_jatom, fb_triple, beta, dir, ieqj, trace, iatom, &
!$OMP                  n_nb, ikind, TASK, inb, jatom, cell_b, img, ispin, found, kkind, iac, &
!$OMP                  alist_ac, kac_lo, kac_hi, k_total, kac, clist_a, n_cont_a, nsoctot, katom, &
!$OMP                  int_local_h, int_local_s, C_coeff_hh_a, C_coeff_ss_a, jkind, ibc, alist_bc, &
!$OMP                  kbc, clist_b, n_cont_b, list_a, list_b, dista, C_coeff_hh_b, C_coeff_ss_b, &
!$OMP                  rac, rbc, do_forces, force_tmp, ia_kind, ja_kind, ka_kind, lock_num) &
!$OMP          REDUCTION(+:pv_virial_thread)

      ALLOCATE (nb_hblk(nspins, max_nb), nb_pblk(nspins, max_nb), k_buf(ku_alloc, nspins), &
                k_offset(max_nclist), C_int_h(MAX(1, max_gau*max_nsgf)), &
                C_int_s(MAX(1, max_gau*max_nsgf)), coc(MAX(1, max_gau), MAX(1, max_gau)), &
                a_matrix(MAX(1, max_nsgf), MAX(1, max_nsgf)), p_matrix(MAX(1, max_nsgf), MAX(1, max_nsgf)))
      IF (forces) ALLOCATE (u_buf(ku_alloc, nspins), u_touched(max_nclist, nspins), f_katom(3, max_nclist))

      IF (forces) THEN
!$OMP SINGLE
!$       ALLOCATE (locks(natom*nkind))
!$OMP END SINGLE

!$OMP DO
!$       do lock_num = 1, natom*nkind
!$          call omp_init_lock(locks(lock_num))
!$       end do
!$OMP END DO
      END IF

!$OMP DO SCHEDULE(DYNAMIC, 1)
      DO iatom = 1, natom
         n_nb = ia_ptr(iatom + 1) - ia_ptr(iatom)
         IF (n_nb == 0) CYCLE
         task => my_sab(1)%nlist_task(slot_order(ia_ptr(iatom)))
         ikind = MERGE(task%jkind, task%ikind, task%iatom > task%jatom)

         DO inb = 1, n_nb
            task => my_sab(1)%nlist_task(slot_order(ia_ptr(iatom) + inb - 1))
            jatom = MAX(task%iatom, task%jatom)
            img = 1
            IF (nimages > 1) img = cell_to_index(task%cell(1), task%cell(2), task%cell(3))
            CPASSERT(img > 0)
            DO ispin = 1, nspins
               CALL dbcsr_get_block_p(ksmat(nspins*(img - 1) + ispin)%matrix, iatom, jatom, nb_hblk(ispin, inb)%array, found)
               CPASSERT(found)
               IF (forces) THEN
                  CALL dbcsr_get_block_p(pmat(nspins*(img - 1) + ispin)%matrix, iatom, jatom, nb_pblk(ispin, inb)%array, found)
                  CPASSERT(found)
               END IF
            END DO
         END DO

         DO kkind = 1, nkind
            IF (.NOT. paw_kind(kkind)) CYCLE
            iac = ikind + nkind*(kkind - 1)
            IF (.NOT. has_intac(iac)) CYCLE
            CALL get_alist(my_oce%intac(iac), alist_ac, iatom)
            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE

            ! Precompute K_a = C_a^h V_h - C_a^s V_s and accumulate U_a = sum_(B/=K) f P_AB C_b^s.
            kac_lo = 1
            DO WHILE (kac_lo <= alist_ac%nclist)
               k_total = 0
               DO kac = kac_lo, alist_ac%nclist
                  clist_a => alist_ac%clist(kac)
                  n_cont_a = clist_a%nsgf_cnt
                  IF (n_cont_a == 0) CYCLE
                  nsoctot = SIZE(clist_a%acint, 2)
                  IF (k_total > 0 .AND. k_total + n_cont_a*nsoctot > ku_alloc) EXIT
                  katom = clist_a%catom
                  k_offset(kac) = k_total + 1
                  k_total = k_total + n_cont_a*nsoctot
                  CALL get_rho_atom(rho_atom=my_rho_atom(katom), int_scr_h=int_local_h, int_scr_s=int_local_s)
                  C_coeff_ss_a => clist_a%acint
                  IF (iatom == katom .AND. ALL(clist_a%cell == 0)) THEN
                     C_coeff_hh_a => clist_a%achint
                     DO ispin = 1, nspins
                        CALL DGEMM('N', 'N', n_cont_a, nsoctot, nsoctot, 1.0_dp, C_coeff_hh_a(:, :, 1), n_cont_a, &
                                   int_local_h(ispin)%r_coef, nsoctot, 0.0_dp, k_buf(k_offset(kac):, ispin), n_cont_a)
                        CALL DGEMM('N', 'N', n_cont_a, nsoctot, nsoctot, -1.0_dp, C_coeff_ss_a(:, :, 1), n_cont_a, &
                                   int_local_s(ispin)%r_coef, nsoctot, 1.0_dp, k_buf(k_offset(kac):, ispin), n_cont_a)
                     END DO
                  ELSE
                     DO ispin = 1, nspins
                        coc(1:nsoctot, 1:nsoctot) = int_local_h(ispin)%r_coef - int_local_s(ispin)%r_coef
                        CALL DGEMM('N', 'N', n_cont_a, nsoctot, nsoctot, 1.0_dp, C_coeff_ss_a(:, :, 1), n_cont_a, &
                                   coc, SIZE(coc, 1), 0.0_dp, k_buf(k_offset(kac):, ispin), n_cont_a)
                     END DO
                  END IF
               END DO
               kac_hi = kac - 1

               IF (k_total > 0) THEN
                  IF (forces) THEN
                     u_touched(kac_lo:kac_hi, 1:nspins) = .FALSE.
                     f_katom(1:3, kac_lo:kac_hi) = 0.0_dp
                  END IF
                  DO inb = 1, n_nb
                     task => my_sab(1)%nlist_task(slot_order(ia_ptr(iatom) + inb - 1))
                     jkind = MERGE(task%ikind, task%jkind, task%iatom > task%jatom)
                     jatom = MAX(task%iatom, task%jatom)
                     cell_b = MERGE(-task%cell, task%cell, task%iatom > task%jatom)
                     ibc = jkind + nkind*(kkind - 1)
                     IF (.NOT. has_intac(ibc)) CYCLE
                     CALL get_alist(my_oce%intac(ibc), alist_bc, jatom)
                     IF (.NOT. ASSOCIATED(alist_bc)) CYCLE

                     f_jatom = 0.0_dp
                     ieqj = MERGE(1.0_dp, 2.0_dp, iatom == jatom)

                     DO kbc = 1, alist_bc%nclist
                        clist_b => alist_bc%clist(kbc)
                        n_cont_b = clist_b%nsgf_cnt
                        IF (n_cont_b == 0) CYCLE
                        katom = clist_b%catom
                        DO kac = kac_lo, kac_hi
                           IF (alist_ac%clist(kac)%catom /= katom) CYCLE
                           IF (ALL(cell_b + clist_b%cell == alist_ac%clist(kac)%cell)) EXIT
                        END DO
                        IF (kac > kac_hi) CYCLE
                        clist_a => alist_ac%clist(kac)
                        n_cont_a = clist_a%nsgf_cnt
                        IF (n_cont_a == 0) CYCLE
                        list_a => clist_a%sgf_list
                        list_b => clist_b%sgf_list

                        IF (jatom /= katom .OR. ANY(clist_b%cell /= 0)) THEN
                           nsoctot = SIZE(clist_a%acint, 2)
                           C_coeff_ss_b => clist_b%acint
                           DO ispin = 1, nspins
                              CALL DGEMM('N', 'T', n_cont_a, n_cont_b, nsoctot, factor1, k_buf(k_offset(kac):, ispin), &
                                         n_cont_a, C_coeff_ss_b(:, :, 1), n_cont_b, 0.0_dp, a_matrix, SIZE(a_matrix, 1))
                              CALL alist_post_align_blk(a_matrix, SIZE(a_matrix, 1), nb_hblk(ispin, inb)%array, &
                                                        SIZE(nb_hblk(ispin, inb)%array, 1), list_a, n_cont_a, list_b, n_cont_b)
                           END DO
                           IF (forces) THEN
                              dista = (iatom /= katom .OR. ANY(clist_a%cell /= 0))
                              fb_triple = 0.0_dp
                              DO ispin = 1, nspins
                                 CALL alist_pre_align_blk(nb_pblk(ispin, inb)%array, SIZE(nb_pblk(ispin, inb)%array, 1), &
                                                          p_matrix, SIZE(p_matrix, 1), list_a, n_cont_a, list_b, n_cont_b)
                                 CALL DGEMM('T', 'N', n_cont_b, nsoctot, n_cont_a, 1.0_dp, p_matrix, SIZE(p_matrix, 1), &
                                            k_buf(k_offset(kac):, ispin), n_cont_a, 0.0_dp, C_int_s, n_cont_b)
                                 IF (dista) THEN
                                    beta = MERGE(1.0_dp, 0.0_dp, u_touched(kac, ispin))
                                    u_touched(kac, ispin) = .TRUE.
                                    CALL DGEMM('N', 'N', n_cont_a, nsoctot, n_cont_b, ieqj*force_fac(ispin), p_matrix, &
                                               SIZE(p_matrix, 1), C_coeff_ss_b(:, :, 1), n_cont_b, beta, &
                                               u_buf(k_offset(kac):, ispin), n_cont_a)
                                 END IF
                                 DO dir = 2, 4
                                    trace = factor2*ieqj*force_fac(ispin)* &
                                            trace_r_AxB(C_int_s, n_cont_b, C_coeff_ss_b(:, :, dir), n_cont_b, n_cont_b, nsoctot)
                                    f_katom(dir - 1, kac) = f_katom(dir - 1, kac) + trace
                                    fb_triple(dir - 1) = fb_triple(dir - 1) - trace
                                 END DO
                              END DO
                              f_jatom = f_jatom + fb_triple
                              IF (use_virial) CALL virial_pair_force(pv_virial_thread, 1.0_dp, fb_triple, clist_b%rac)
                           END IF
                        ELSE
                           C_coeff_ss_a => clist_a%acint
                           IF (iatom == katom .AND. ALL(clist_a%cell == 0)) THEN
                              C_coeff_hh_a => clist_a%achint
                              dista = .FALSE.
                           ELSE
                              C_coeff_hh_a => clist_a%acint
                              dista = .TRUE.
                           END IF
                           C_coeff_hh_b => clist_b%achint
                           C_coeff_ss_b => clist_b%acint
                           nsoctot = SIZE(C_coeff_ss_a, 2)
                           CALL get_rho_atom(rho_atom=my_rho_atom(katom), int_scr_h=int_local_h, int_scr_s=int_local_s)
                           do_forces = forces .AND. dista
                           ! Use the stored block ordering for both KS and force contractions.
                           CALL add_vhxca(nb_hblk(:, inb), nb_pblk(:, inb), C_coeff_hh_a, C_coeff_hh_b, &
                                          C_coeff_ss_a, C_coeff_ss_b, int_local_h, int_local_s, nspins, &
                                          iatom, jatom, nsoctot, factor1, list_a, n_cont_a, list_b, n_cont_b, &
                                          dista, .FALSE., C_int_h, C_int_s, MAX(1, max_nsgf), coc, a_matrix, &
                                          p_matrix, do_forces, force_tmp, force_fac)

                           IF (do_forces) THEN
                              IF (use_virial) THEN
                                 rac = alist_ac%clist(kac)%rac
                                 rbc = alist_bc%clist(kbc)%rac
                              END IF
                              ia_kind = atom_of_kind(iatom)
                              ja_kind = atom_of_kind(jatom)
                              ka_kind = atom_of_kind(katom)
                              force_tmp = factor2*force_tmp
!$                            CALL omp_set_lock(locks((ka_kind - 1)*nkind + kkind))
                              force(kkind)%vhxc_atom(1:3, ka_kind) = force(kkind)%vhxc_atom(1:3, ka_kind) + force_tmp(1:3, 3)
!$                            CALL omp_unset_lock(locks((ka_kind - 1)*nkind + kkind))
!$                            CALL omp_set_lock(locks((ia_kind - 1)*nkind + ikind))
                              force(ikind)%vhxc_atom(1:3, ia_kind) = force(ikind)%vhxc_atom(1:3, ia_kind) + force_tmp(1:3, 1)
!$                            CALL omp_unset_lock(locks((ia_kind - 1)*nkind + ikind))
!$                            CALL omp_set_lock(locks((ja_kind - 1)*nkind + jkind))
                              force(jkind)%vhxc_atom(1:3, ja_kind) = force(jkind)%vhxc_atom(1:3, ja_kind) + force_tmp(1:3, 2)
!$                            CALL omp_unset_lock(locks((ja_kind - 1)*nkind + jkind))
                              IF (use_virial) THEN
                                 CALL virial_pair_force(pv_virial_thread, 1._dp, force_tmp(1:3, 1), rac)
                                 CALL virial_pair_force(pv_virial_thread, 1._dp, force_tmp(1:3, 2), rbc)
                              END IF
                           END IF
                        END IF
                     END DO ! kbc

                     IF (forces .AND. ANY(f_jatom /= 0.0_dp)) THEN
                        ja_kind = atom_of_kind(jatom)
!$                      CALL omp_set_lock(locks((ja_kind - 1)*nkind + jkind))
                        force(jkind)%vhxc_atom(:, ja_kind) = force(jkind)%vhxc_atom(:, ja_kind) + f_jatom
!$                      CALL omp_unset_lock(locks((ja_kind - 1)*nkind + jkind))
                     END IF
                  END DO

                  IF (forces) THEN
                     ia_kind = atom_of_kind(iatom)
                     DO kac = kac_lo, kac_hi
                        clist_a => alist_ac%clist(kac)
                        n_cont_a = clist_a%nsgf_cnt
                        katom = clist_a%catom
                        IF (ANY(u_touched(kac, 1:nspins))) THEN
                           CALL get_rho_atom(rho_atom=my_rho_atom(katom), int_scr_h=int_local_h, int_scr_s=int_local_s)
                           C_coeff_ss_a => clist_a%acint
                           nsoctot = SIZE(C_coeff_ss_a, 2)
                           fb_triple = 0.0_dp
                           DO ispin = 1, nspins
                              IF (.NOT. u_touched(kac, ispin)) CYCLE
                              coc(1:nsoctot, 1:nsoctot) = int_local_h(ispin)%r_coef - int_local_s(ispin)%r_coef
                              CALL DGEMM('N', 'N', n_cont_a, nsoctot, nsoctot, 1.0_dp, u_buf(k_offset(kac):, ispin), &
                                         n_cont_a, coc, SIZE(coc, 1), 0.0_dp, C_int_s, n_cont_a)
                              DO dir = 2, 4
                                 trace = factor2*trace_r_AxB(C_int_s, n_cont_a, C_coeff_ss_a(:, :, dir), &
                                                             n_cont_a, n_cont_a, nsoctot)
                                 f_katom(dir - 1, kac) = f_katom(dir - 1, kac) + trace
                                 fb_triple(dir - 1) = fb_triple(dir - 1) - trace
                              END DO
                           END DO
!$                         CALL omp_set_lock(locks((ia_kind - 1)*nkind + ikind))
                           force(ikind)%vhxc_atom(:, ia_kind) = force(ikind)%vhxc_atom(:, ia_kind) + fb_triple
!$                         CALL omp_unset_lock(locks((ia_kind - 1)*nkind + ikind))
                           IF (use_virial) CALL virial_pair_force(pv_virial_thread, 1.0_dp, fb_triple, clist_a%rac)
                        END IF
                        IF (ANY(f_katom(:, kac) /= 0.0_dp)) THEN
                           ka_kind = atom_of_kind(katom)
!$                         CALL omp_set_lock(locks((ka_kind - 1)*nkind + kkind))
                           force(kkind)%vhxc_atom(:, ka_kind) = force(kkind)%vhxc_atom(:, ka_kind) + f_katom(:, kac)
!$                         CALL omp_unset_lock(locks((ka_kind - 1)*nkind + kkind))
                        END IF
                     END DO
                  END IF
               END IF
               kac_lo = kac_hi + 1
            END DO
         END DO
      END DO
!$OMP END DO

      DEALLOCATE (nb_hblk, nb_pblk, k_buf, k_offset, C_int_h, C_int_s, a_matrix, p_matrix, coc)
      IF (forces) DEALLOCATE (u_buf, u_touched, f_katom)

      IF (forces) THEN
         ! Implicit barrier at end of OMP DO, so locks can be freed
!$OMP DO
!$       DO lock_num = 1, natom*nkind
!$          call omp_destroy_lock(locks(lock_num))
!$       END DO
!$OMP END DO

!$OMP SINGLE
!$       DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT
      END IF

!$OMP END PARALLEL

      IF (use_virial) THEN
         virial%pv_gapw(1:3, 1:3) = virial%pv_gapw(1:3, 1:3) + pv_virial_thread(1:3, 1:3)
         virial%pv_virial(1:3, 1:3) = virial%pv_virial(1:3, 1:3) + pv_virial_thread(1:3, 1:3)
      END IF

      DEALLOCATE (ia_ptr, slot_order)

      DEALLOCATE (basis_set_list, paw_kind, has_intac)

      CALL timestop(handle)

   END SUBROUTINE update_ks_atom

! **************************************************************************************************
!> \brief Add one-center hard-minus-soft KS and coefficient-derivative force contributions.
!>        Real symmetric one-center potentials allow the same unscaled C*V intermediates
!>        to be used for KS assembly and forces within each spin.
!> \param mat_h KS blocks in stored (A <= B) order
!> \param mat_p density blocks in the same order, used only for forces
!> \param C_hh_a hard coefficients and derivatives on A
!> \param C_hh_b hard coefficients and derivatives on B
!> \param C_ss_a soft coefficients and derivatives on A
!> \param C_ss_b soft coefficients and derivatives on B
!> \param int_local_h hard one-center potentials by spin
!> \param int_local_s soft one-center potentials by spin
!> \param nspins number of spins
!> \param ia atom A
!> \param ja atom B
!> \param nsp size of the one-center basis
!> \param factor KS scaling, kept separate from the force scaling
!> \param lista contracted basis indices on A
!> \param nconta number of contracted functions on A
!> \param listb contracted basis indices on B
!> \param ncontb number of contracted functions on B
!> \param dista true if A is offsite relative to the expansion center
!> \param distb true if B is offsite relative to the expansion center
!> \param work_h workspace for unscaled potential-weighted coefficients
!> \param work_s workspace for soft or density-weighted coefficients
!> \param ldwork leading dimension of the coefficient workspaces
!> \param int_diff workspace for the hard-minus-soft potential
!> \param a_matrix workspace for the contracted KS block
!> \param p_matrix workspace for the aligned density block
!> \param do_forces whether nonzero coefficient-derivative contributions are needed
!> \param force forces on A, B, and the expansion center
!> \param force_scaling force scaling by spin
! **************************************************************************************************
   SUBROUTINE add_vhxca(mat_h, mat_p, C_hh_a, C_hh_b, C_ss_a, C_ss_b, &
                        int_local_h, int_local_s, nspins, ia, ja, nsp, factor, &
                        lista, nconta, listb, ncontb, dista, distb, &
                        work_h, work_s, ldwork, int_diff, a_matrix, p_matrix, &
                        do_forces, force, force_scaling)
      TYPE(cp_2d_r_p_type), DIMENSION(:), INTENT(IN)     :: mat_h, mat_p
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: C_hh_a, C_hh_b, C_ss_a, C_ss_b
      TYPE(rho_atom_coeff), DIMENSION(:), INTENT(IN)     :: int_local_h, int_local_s
      INTEGER, INTENT(IN)                                :: nspins, ia, ja, nsp
      REAL(KIND=dp), INTENT(IN)                          :: factor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lista
      INTEGER, INTENT(IN)                                :: nconta
      INTEGER, DIMENSION(:), INTENT(IN)                  :: listb
      INTEGER, INTENT(IN)                                :: ncontb
      LOGICAL, INTENT(IN)                                :: dista, distb
      INTEGER, INTENT(IN)                                :: ldwork
      REAL(KIND=dp), INTENT(OUT)                         :: work_s(ldwork, nsp), &
                                                            work_h(ldwork, nsp), int_diff(nsp, nsp)
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: a_matrix, p_matrix
      LOGICAL, INTENT(IN)                                :: do_forces
      REAL(dp), DIMENSION(3, 3), INTENT(OUT)             :: force
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: force_scaling

      INTEGER                                            :: dir, ispin
      REAL(dp), DIMENSION(:, :), POINTER                 :: h_block, int_hard, int_soft, p_block
      REAL(KIND=dp)                                      :: ieqj, trace

      force(:, :) = 0.0_dp
      ! Factor 2 because only half of the pairs with ia /= ja are considered.
      ieqj = 2.0_dp
      IF (ia == ja) ieqj = 1.0_dp

      DO ispin = 1, nspins
         h_block => mat_h(ispin)%array
         int_hard => int_local_h(ispin)%r_coef
         int_soft => int_local_s(ispin)%r_coef

         IF (dista .AND. distb) int_diff(:, :) = int_hard(1:nsp, 1:nsp) - int_soft(1:nsp, 1:nsp)

         IF (distb .AND. .NOT. dista) THEN
            ! B is offsite: H_ab = (C_a^h V_h - C_a^s V_s) C_b^T.
            CALL DGEMM('N', 'N', nconta, nsp, nsp, 1.0_dp, C_hh_a(:, :, 1), SIZE(C_hh_a, 1), &
                       int_hard, SIZE(int_hard, 1), 0.0_dp, work_h, ldwork)
            CALL DGEMM('N', 'N', nconta, nsp, nsp, -1.0_dp, C_ss_a(:, :, 1), SIZE(C_ss_a, 1), &
                       int_soft, SIZE(int_soft, 1), 1.0_dp, work_h, ldwork)
            CALL DGEMM('N', 'T', nconta, ncontb, nsp, factor, work_h, ldwork, &
                       C_hh_b(:, :, 1), SIZE(C_hh_b, 1), 0.0_dp, a_matrix, SIZE(a_matrix, 1))
         ELSE
            IF (distb) THEN
               ! Both ends are offsite: C_h = C_s, including their derivatives.
               CALL DGEMM('N', 'N', ncontb, nsp, nsp, 1.0_dp, C_hh_b(:, :, 1), SIZE(C_hh_b, 1), &
                          int_diff, nsp, 0.0_dp, work_h, ldwork)
            ELSE
               CALL DGEMM('N', 'N', ncontb, nsp, nsp, 1.0_dp, C_hh_b(:, :, 1), SIZE(C_hh_b, 1), &
                          int_hard, SIZE(int_hard, 1), 0.0_dp, work_h, ldwork)
               IF (dista) THEN
                  CALL DGEMM('N', 'N', ncontb, nsp, nsp, -1.0_dp, C_ss_b(:, :, 1), SIZE(C_ss_b, 1), &
                             int_soft, SIZE(int_soft, 1), 1.0_dp, work_h, ldwork)
               ELSE
                  CALL DGEMM('N', 'N', ncontb, nsp, nsp, 1.0_dp, C_ss_b(:, :, 1), SIZE(C_ss_b, 1), &
                             int_soft, SIZE(int_soft, 1), 0.0_dp, work_s, ldwork)
               END IF
            END IF
            CALL DGEMM('N', 'T', nconta, ncontb, nsp, factor, C_hh_a(:, :, 1), SIZE(C_hh_a, 1), &
                       work_h, ldwork, 0.0_dp, a_matrix, SIZE(a_matrix, 1))
            IF (.NOT. dista) THEN
               CALL DGEMM('N', 'T', nconta, ncontb, nsp, -factor, C_ss_a(:, :, 1), SIZE(C_ss_a, 1), &
                          work_s, ldwork, 1.0_dp, a_matrix, SIZE(a_matrix, 1))
            END IF
         END IF

         CALL alist_post_align_blk(a_matrix, SIZE(a_matrix, 1), h_block, SIZE(h_block, 1), &
                                   lista, nconta, listb, ncontb)

         IF (.NOT. do_forces) CYCLE
         p_block => mat_p(ispin)%array
         CALL alist_pre_align_blk(p_block, SIZE(p_block, 1), p_matrix, SIZE(p_matrix, 1), &
                                  lista, nconta, listb, ncontb)

         ! Onsite coefficient derivatives are zero (see build_oce_matrices).
         ! The offsite derivative is common to hard and soft coefficients.
         IF (dista) THEN
            ! Reuse K_b = C_b^h V_h - C_b^s V_s from KS assembly.
            CALL DGEMM('N', 'N', nconta, nsp, ncontb, 1.0_dp, p_matrix, SIZE(p_matrix, 1), &
                       work_h, ldwork, 0.0_dp, work_s, ldwork)
            ! <(P C_b)^T dC_a, V>_F = <dC_a, P (C_b V)>_F.
            DO dir = 2, 4
               trace = trace_r_AxB(work_s, ldwork, C_hh_a(:, :, dir), SIZE(C_hh_a, 1), nconta, nsp)
               trace = ieqj*trace*force_scaling(ispin)
               force(dir - 1, 3) = force(dir - 1, 3) + trace
               force(dir - 1, 1) = force(dir - 1, 1) - trace
            END DO
         END IF

         IF (distb) THEN
            IF (dista) THEN
               ! A force used K_b; now form K_a = C_a (V_h - V_s).
               CALL DGEMM('N', 'N', nconta, nsp, nsp, 1.0_dp, C_hh_a(:, :, 1), SIZE(C_hh_a, 1), &
                          int_diff, nsp, 0.0_dp, work_h, ldwork)
            END IF
            ! For an onsite A, K_a is already available from KS assembly.
            CALL DGEMM('T', 'N', ncontb, nsp, nconta, 1.0_dp, p_matrix, SIZE(p_matrix, 1), &
                       work_h, ldwork, 0.0_dp, work_s, ldwork)
            DO dir = 2, 4
               trace = trace_r_AxB(work_s, ldwork, C_hh_b(:, :, dir), SIZE(C_hh_b, 1), ncontb, nsp)
               trace = ieqj*trace*force_scaling(ispin)
               force(dir - 1, 3) = force(dir - 1, 3) + trace
               force(dir - 1, 2) = force(dir - 1, 2) - trace
            END DO
         END IF
      END DO
   END SUBROUTINE add_vhxca

END MODULE qs_ks_atom
