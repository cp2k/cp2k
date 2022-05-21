!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief 3-center integrals machinery for the XAS_TDP method
!> \author A. Bussy (03.2020)
! **************************************************************************************************

MODULE xas_tdp_integrals
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              libxsmm_abc_contract
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_2d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_param,&
                                              cp_eri_mme_set_params
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_distribution_get,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type
   USE dbt_api,                         ONLY: &
        dbt_create, dbt_distribution_destroy, dbt_distribution_new, dbt_distribution_type, &
        dbt_finalize, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_put_block, &
        dbt_reserve_blocks, dbt_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                              distribution_2d_type
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate
   USE eri_mme_types,                   ONLY: eri_mme_init,&
                                              eri_mme_release
   USE gamma,                           ONLY: init_md_ftable
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor,&
                                              eri_2center,&
                                              eri_3center,&
                                              libint_potential_type
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_cleanup_3eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE mathlib,                         ONLY: invmat_symm
   USE message_passing,                 ONLY: mp_sum,&
                                              mp_sync
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: &
        get_iterator_info, get_neighbor_list_set_p, neighbor_list_iterate, &
        neighbor_list_iterator_create, neighbor_list_iterator_p_type, &
        neighbor_list_iterator_release, neighbor_list_set_p_type, nl_set_sub_iterator, &
        nl_sub_iterate, release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_o3c_types,                    ONLY: get_o3c_iterator_info,&
                                              init_o3c_container,&
                                              o3c_container_type,&
                                              o3c_iterate,&
                                              o3c_iterator_create,&
                                              o3c_iterator_release,&
                                              o3c_iterator_type,&
                                              release_o3c_container
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init
   USE xas_tdp_types,                   ONLY: xas_tdp_control_type,&
                                              xas_tdp_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_integrals'

   PUBLIC :: create_pqX_tensor, fill_pqX_tensor, compute_ri_3c_coulomb, compute_ri_3c_exchange, &
             build_xas_tdp_3c_nl, build_xas_tdp_ovlp_nl, get_opt_3c_dist2d, &
             compute_ri_coulomb2_int, compute_ri_exchange2_int

CONTAINS

! **************************************************************************************************
!> \brief Prepares a tensor to hold 3-center integrals (pq|X), where p,q are distributed according
!>        to the given 2d dbcsr distribution of the given . The third dimension of the tensor is
!>        iteslf not distributed (i.e. the t_pgrid's third dimension has size 1). The blocks are
!>        reserved according to the neighbor lists
!> \param pq_X the tensor to store the integrals
!> \param ab_nl the 1st and 2nd center neighbor list
!> \param ac_nl the 1st and 3rd center neighbor list
!> \param matrix_dist ...
!> \param blk_size_1 the block size in the first dimension
!> \param blk_size_2 the block size in the second dimension
!> \param blk_size_3 the block size in the third dimension
!> \param only_bc_same_center only keep block if, for the corresponding integral (ab|c), b and c
!>        share the same center, i.e. r_bc = 0.0
! **************************************************************************************************
   SUBROUTINE create_pqX_tensor(pq_X, ab_nl, ac_nl, matrix_dist, blk_size_1, blk_size_2, &
                                blk_size_3, only_bc_same_center)

      TYPE(dbt_type), INTENT(OUT)                        :: pq_X
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_nl, ac_nl
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: matrix_dist
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size_1, blk_size_2, blk_size_3
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_bc_same_center

      CHARACTER(len=*), PARAMETER                        :: routineN = 'create_pqX_tensor'

      INTEGER                                            :: group, handle, i, iatom, ikind, jatom, &
                                                            katom, kkind, nblk, nblk_3, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx1, idx2, idx3
      INTEGER, DIMENSION(3)                              :: pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: mat_pgrid
      LOGICAL                                            :: my_sort_bc, symmetric
      REAL(dp), DIMENSION(3)                             :: rab, rac, rbc
      TYPE(dbt_distribution_type)                        :: t_dist
      TYPE(dbt_pgrid_type)                               :: t_pgrid
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ab_iter, ac_iter

      NULLIFY (ab_iter, ac_iter, col_dist, row_dist, mat_pgrid)

      CALL timeset(routineN, handle)

      my_sort_bc = .FALSE.
      IF (PRESENT(only_bc_same_center)) my_sort_bc = only_bc_same_center

      CALL get_neighbor_list_set_p(ab_nl, symmetric=symmetric)
      CPASSERT(symmetric)

      !get 2D distribution info from matrix_dist
      CALL dbcsr_distribution_get(matrix_dist, pgrid=mat_pgrid, group=group, &
                                  row_dist=row_dist, col_dist=col_dist)

      !create the corresponding tensor proc grid and dist
      pdims(1) = SIZE(mat_pgrid, 1); pdims(2) = SIZE(mat_pgrid, 2); pdims(3) = 1
      CALL dbt_pgrid_create(group, pdims, t_pgrid)

      nblk_3 = SIZE(blk_size_3)
      CALL dbt_distribution_new(t_dist, t_pgrid, nd_dist_1=row_dist, nd_dist_2=col_dist, &
                                nd_dist_3=[(0, i=1, nblk_3)])

      !create the tensor itself.
      CALL dbt_create(pq_X, name="(pq|X)", dist=t_dist, map1_2d=[1, 2], map2_2d=[3], &
                      blk_size_1=blk_size_1, blk_size_2=blk_size_2, blk_size_3=blk_size_3)

      !count the blocks to reserve !note: dbcsr takes care of only keeping unique indices
      CALL neighbor_list_iterator_create(ab_iter, ab_nl)
      CALL neighbor_list_iterator_create(ac_iter, ac_nl, search=.TRUE.)
      nblk = 0
      DO WHILE (neighbor_list_iterate(ab_iter) == 0)
         CALL get_iterator_info(ab_iter, ikind=ikind, iatom=iatom, nkind=nkind, r=rab)

         DO kkind = 1, nkind
            CALL nl_set_sub_iterator(ac_iter, ikind, kkind, iatom)

            DO WHILE (nl_sub_iterate(ac_iter) == 0)

               IF (my_sort_bc) THEN
                  !we check for rbc or rac because of symmetry in ab_nl
                  CALL get_iterator_info(ac_iter, r=rac)
                  rbc(:) = rac(:) - rab(:)
                  IF (.NOT. (ALL(ABS(rbc) .LE. 1.0E-8_dp) .OR. ALL(ABS(rac) .LE. 1.0E-8_dp))) CYCLE

               END IF

               nblk = nblk + 1
            END DO !ac_iter
         END DO !kkind
      END DO !ab_iter
      CALL neighbor_list_iterator_release(ab_iter)
      CALL neighbor_list_iterator_release(ac_iter)

      ALLOCATE (idx1(nblk), idx2(nblk), idx3(nblk))

      !actually reserve the blocks
      CALL neighbor_list_iterator_create(ab_iter, ab_nl)
      CALL neighbor_list_iterator_create(ac_iter, ac_nl, search=.TRUE.)
      nblk = 0
      DO WHILE (neighbor_list_iterate(ab_iter) == 0)
         CALL get_iterator_info(ab_iter, ikind=ikind, iatom=iatom, jatom=jatom, nkind=nkind, r=rab)

         DO kkind = 1, nkind
            CALL nl_set_sub_iterator(ac_iter, ikind, kkind, iatom)

            DO WHILE (nl_sub_iterate(ac_iter) == 0)
               CALL get_iterator_info(ac_iter, jatom=katom, r=rac)

               IF (my_sort_bc) THEN
                  !we check for rbc or rac because of symmetry in ab_nl
                  CALL get_iterator_info(ac_iter, r=rac)
                  rbc(:) = rac(:) - rab(:)
                  IF (.NOT. (ALL(ABS(rbc) .LE. 1.0E-8_dp) .OR. ALL(ABS(rac) .LE. 1.0E-8_dp))) CYCLE

               END IF

               nblk = nblk + 1

               idx1(nblk) = iatom
               idx2(nblk) = jatom
               idx3(nblk) = katom

            END DO !ac_iter
         END DO !kkind
      END DO !ab_iter
      CALL neighbor_list_iterator_release(ab_iter)
      CALL neighbor_list_iterator_release(ac_iter)

!TODO: Parallelize this properly.
!$OMP PARALLEL DEFAULT(NONE) SHARED(pq_X, idx1, idx2, idx3)
      CALL dbt_reserve_blocks(pq_X, idx1, idx2, idx3)
!$OMP END PARALLEL
      CALL dbt_finalize(pq_X)

      !clean-up
      CALL dbt_distribution_destroy(t_dist)
      CALL dbt_pgrid_destroy(t_pgrid)

      CALL timestop(handle)

   END SUBROUTINE create_pqX_tensor

! **************************************************************************************************
!> \brief Fills the given 3 dimensional (pq|X) tensor with integrals
!> \param pq_X the tensor to fill
!> \param ab_nl the neighbor list for the first 2 centers
!> \param ac_nl the neighbor list for the first and third centers
!> \param basis_set_list_a basis sets for first center
!> \param basis_set_list_b basis sets for second center
!> \param basis_set_list_c basis sets for third center
!> \param potential_parameter the operator for the integrals
!> \param qs_env ...
!> \param only_bc_same_center same as in create_pqX_tensor
!> \param eps_screen threshold for possible screening
!> \note The following indices are happily mixed within this routine: First center i,a,p
!>       Second center: j,b,q       Third center: k,c,X
! **************************************************************************************************
   SUBROUTINE fill_pqX_tensor(pq_X, ab_nl, ac_nl, basis_set_list_a, basis_set_list_b, &
                              basis_set_list_c, potential_parameter, qs_env, &
                              only_bc_same_center, eps_screen)

      TYPE(dbt_type)                                     :: pq_X
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_nl, ac_nl
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(libint_potential_type)                        :: potential_parameter
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_bc_same_center
      REAL(dp), INTENT(IN), OPTIONAL                     :: eps_screen

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fill_pqX_tensor'

      INTEGER :: egfa, egfb, egfc, handle, i, iatom, ibasis, ikind, ilist, imax, iset, jatom, &
         jkind, jset, katom, kkind, kset, m_max, max_ncob, max_ncoc, max_nset, max_nsgfa, &
         max_nsgfb, maxli, maxlj, maxlk, mepos, nbasis, ncoa, ncob, ncoc, ni, nj, nk, nseta, &
         nsetb, nsetc, nthread, sgfa, sgfb, sgfc, unit_id
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: do_screen, my_sort_bc
      REAL(dp)                                           :: dij, dik, djk, my_eps_screen, ri(3), &
                                                            rij(3), rik(3), rj(3), rjk(3), rk(3), &
                                                            sabc_ext, screen_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: ccp_buffer, cpp_buffer
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: max_contr, max_contra, max_contrb, &
                                                            max_contrc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: iabc, sabc, work
      REAL(dp), DIMENSION(:), POINTER                    :: set_radius_a, set_radius_b, set_radius_c
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf_a, rpgf_b, rpgf_c, zeta, zetb, zetc
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spb, spc, tspa
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, basis_set_a, basis_set_b, &
                                                            basis_set_c
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      NULLIFY (basis_set, basis_set_list, para_env, la_max, la_min)
      NULLIFY (lb_max, lb_min, lc_max, lc_min, npgfa, npgfb, npgfc, nsgfa, nsgfb, nsgfc)
      NULLIFY (first_sgfa, first_sgfb, first_sgfc, set_radius_a, set_radius_b, set_radius_c)
      NULLIFY (rpgf_a, rpgf_b, rpgf_c, zeta, zetb, zetc)
      NULLIFY (basis_set_a, basis_set_b, basis_set_c, tspa, spb, spc)

      CALL timeset(routineN, handle)

      !Need the max l for each basis for libint (and overall max #of sets for screening)
      nbasis = SIZE(basis_set_list_a)
      max_nsgfa = 0
      max_nset = 0
      maxli = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_a(ibasis)%gto_basis_set, &
                                maxl=imax, nset=iset, nsgf_set=nsgfa)
         maxli = MAX(maxli, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfa = MAX(max_nsgfa, MAXVAL(nsgfa))
      END DO
      max_nsgfb = 0
      max_ncob = 0
      maxlj = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_b(ibasis)%gto_basis_set, &
                                maxl=imax, nset=iset, nsgf_set=nsgfb, npgf=npgfb)
         maxlj = MAX(maxlj, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfb = MAX(max_nsgfb, MAXVAL(nsgfb))
         max_ncob = MAX(max_ncob, MAXVAL(npgfb)*ncoset(maxlj))
      END DO
      maxlk = 0
      max_ncoc = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_set_list_c(ibasis)%gto_basis_set, &
                                maxl=imax, nset=iset, npgf=npgfc)
         maxlk = MAX(maxlk, imax)
         max_nset = MAX(max_nset, iset)
         max_ncoc = MAX(max_ncoc, MAXVAL(npgfc)*ncoset(maxlk))
      END DO
      m_max = maxli + maxlj + maxlk

      !Screening
      do_screen = .FALSE.
      IF (PRESENT(eps_screen)) THEN
         do_screen = .TRUE.
         my_eps_screen = eps_screen
      END IF
      screen_radius = 0.0_dp
      IF (potential_parameter%potential_type == do_potential_truncated .OR. &
          potential_parameter%potential_type == do_potential_short) THEN

         screen_radius = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (potential_parameter%potential_type == do_potential_coulomb) THEN

         screen_radius = 1000000.0_dp
      END IF

      !get maximum contraction values for abc_contract screening
      IF (do_screen) THEN

         !Allocate max_contraction arrays such that we have a specific value for each set/kind
         ALLOCATE (max_contr(max_nset, nbasis), max_contra(max_nset, nbasis), &
                   max_contrb(max_nset, nbasis), max_contrc(max_nset, nbasis))

         !Not the most elegent, but better than copying 3 times the same
         DO ilist = 1, 3

            IF (ilist == 1) basis_set_list => basis_set_list_a
            IF (ilist == 2) basis_set_list => basis_set_list_b
            IF (ilist == 3) basis_set_list => basis_set_list_c

            max_contr = 0.0_dp

            DO ibasis = 1, nbasis
               basis_set => basis_set_list(ibasis)%gto_basis_set

               DO iset = 1, basis_set%nset

                  ncoa = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
                  sgfa = basis_set%first_sgf(1, iset)
                  egfa = sgfa + basis_set%nsgf_set(iset) - 1

                  max_contr(iset, ibasis) = &
                     MAXVAL((/(SUM(ABS(basis_set%sphi(1:ncoa, i))), i=sgfa, egfa)/))

               END DO !iset
            END DO !ibasis

            IF (ilist == 1) max_contra(:, :) = max_contr(:, :)
            IF (ilist == 2) max_contrb(:, :) = max_contr(:, :)
            IF (ilist == 3) max_contrc(:, :) = max_contr(:, :)
         END DO !ilist
         DEALLOCATE (max_contr)
      END IF !do_screen

      !To minimize memory ops in contraction, we need to pre-allocate buffers, pre-tranpose sphi_a
      !and also trim sphi in general to have contiguous arrays
      ALLOCATE (tspa(max_nset, nbasis), spb(max_nset, nbasis), spc(max_nset, nbasis))
      DO ibasis = 1, nbasis
         DO iset = 1, max_nset
            NULLIFY (tspa(iset, ibasis)%array)
            NULLIFY (spb(iset, ibasis)%array)
            NULLIFY (spc(iset, ibasis)%array)
         END DO
      END DO

      DO ilist = 1, 3

         DO ibasis = 1, nbasis
            IF (ilist == 1) basis_set => basis_set_list_a(ibasis)%gto_basis_set
            IF (ilist == 2) basis_set => basis_set_list_b(ibasis)%gto_basis_set
            IF (ilist == 3) basis_set => basis_set_list_c(ibasis)%gto_basis_set

            DO iset = 1, basis_set%nset

               ncoa = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
               sgfa = basis_set%first_sgf(1, iset)
               egfa = sgfa + basis_set%nsgf_set(iset) - 1

               IF (ilist == 1) THEN
                  ALLOCATE (tspa(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoa))
                  tspa(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoa, sgfa:egfa))
               ELSE IF (ilist == 2) THEN
                  ALLOCATE (spb(iset, ibasis)%array(ncoa, basis_set%nsgf_set(iset)))
                  spb(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoa, sgfa:egfa)
               ELSE
                  ALLOCATE (spc(iset, ibasis)%array(ncoa, basis_set%nsgf_set(iset)))
                  spc(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoa, sgfa:egfa)
               END IF

            END DO !iset
         END DO !ibasis
      END DO !ilist

      my_sort_bc = .FALSE.
      IF (PRESENT(only_bc_same_center)) my_sort_bc = only_bc_same_center

      !Init the truncated Coulomb operator
      CALL get_qs_env(qs_env, para_env=para_env)
      IF (potential_parameter%potential_type == do_potential_truncated) THEN

         !open the file only if necessary
         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      !Inint the initial gamma function before the OMP region as it is not thread safe
      CALL init_md_ftable(nmax=m_max)

      !Strategy: we use the o3c iterator because it is OMP parallelized and also takes the
      !          only_bc_same_center argument. Only the dbcsr_put_block is critical

      nthread = 1
!$    nthread = omp_get_max_threads()

      ALLOCATE (o3c)
      CALL init_o3c_container(o3c, 1, basis_set_list_a, basis_set_list_b, basis_set_list_c, &
                              ab_nl, ac_nl, only_bc_same_center=my_sort_bc)
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (pq_X,do_screen,max_nset,basis_set_list_a,max_contra,max_contrb,max_contrc,max_nsgfa,&
!$OMP         basis_set_list_b, basis_set_list_c,ncoset,screen_radius,potential_parameter,max_ncob,&
!$OMP         my_eps_screen,maxli,maxlj,maxlk,my_sort_bc,nthread,o3c,o3c_iterator,tspa,spb,spc,&
!$OMP         max_ncoc,max_nsgfb) &
!$OMP PRIVATE (lib,i,mepos,work,iset,ncoa,sgfa,egfa,nseta,&
!$OMP          iatom,ikind,jatom,jkind,katom,kkind,rij,rik,rjk,basis_set_a,nsetb,&
!$OMP          la_max,la_min,lb_max,lb_min,lc_max,lc_min,npgfa,npgfb,npgfc,nsgfa,nsgfb,nsgfc,ri,rk,&
!$OMP          first_sgfa,first_sgfb,first_sgfc,set_radius_a,set_radius_b,set_radius_c, nsetc,rj,&
!$OMP          rpgf_a,rpgf_b,rpgf_c,zeta,zetb,zetc,basis_set_b,basis_set_c,dij,dik,djk,ni,nj,nk,&
!$OMP          iabc,sabc,jset,kset,ncob,ncoc,sgfb,sgfc,egfb,egfc,sabc_ext,cpp_buffer,ccp_buffer)

      mepos = 0
!$    mepos = omp_get_thread_num()

      !pre-allocate work buffers for LIBXSMM contract in order to avoid memory ops
      ALLOCATE (cpp_buffer(max_nsgfa*max_ncob))
      ALLOCATE (ccp_buffer(max_nsgfa*max_nsgfb*max_ncoc))

      !note: we do not initalize libxsmm here, because we assume that if the flag is there, then it
      !      is done in dbcsr already

      !each thread need its own libint object (internals may change at different rates)
      CALL cp_libint_init_3eri(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, ikind=ikind, jkind=jkind, kkind=kkind, &
                                    iatom=iatom, jatom=jatom, katom=katom, rij=rij, rik=rik)

         !get first center basis info
         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         zeta => basis_set_a%zet
         rpgf_a => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         ni = SUM(nsgfa)
         !second center basis info
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         zetb => basis_set_b%zet
         rpgf_b => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         nj = SUM(nsgfb)
         !third center basis info
         basis_set_c => basis_set_list_c(kkind)%gto_basis_set
         first_sgfc => basis_set_c%first_sgf
         lc_max => basis_set_c%lmax
         lc_min => basis_set_c%lmin
         npgfc => basis_set_c%npgf
         nsetc = basis_set_c%nset
         nsgfc => basis_set_c%nsgf_set
         zetc => basis_set_c%zet
         rpgf_c => basis_set_c%pgf_radius
         set_radius_c => basis_set_c%set_radius
         nk = SUM(nsgfc)

         !position and distances, only relative pos matter for libint
         rjk = rik - rij
         ri = 0.0_dp
         rj = rij ! ri + rij
         rk = rik ! ri + rik

         djk = NORM2(rjk)
         dij = NORM2(rij)
         dik = NORM2(rik)

         !sgf integrals
         ALLOCATE (iabc(ni, nj, nk))
         iabc(:, :, :) = 0.0_dp

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            egfa = sgfa + nsgfa(iset) - 1

            DO jset = 1, nsetb
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)
               egfb = sgfb + nsgfb(jset) - 1

               !screening (overlap)
               IF (set_radius_a(iset) + set_radius_b(jset) < dij) CYCLE

               DO kset = 1, nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1, kset)
                  egfc = sgfc + nsgfc(kset) - 1

                  !screening (potential)
                  IF (set_radius_a(iset) + set_radius_c(kset) + screen_radius < dik) CYCLE
                  IF (set_radius_b(jset) + set_radius_c(kset) + screen_radius < djk) CYCLE

                  !pgf integrals
                  ALLOCATE (sabc(ncoa, ncob, ncoc))
                  sabc(:, :, :) = 0.0_dp

                  IF (do_screen) THEN
                     CALL eri_3center(sabc, la_min(iset), la_max(iset), npgfa(iset), zeta(:, iset), &
                                      rpgf_a(:, iset), ri, lb_min(jset), lb_max(jset), npgfb(jset), &
                                      zetb(:, jset), rpgf_b(:, jset), rj, lc_min(kset), lc_max(kset), &
                                      npgfc(kset), zetc(:, kset), rpgf_c(:, kset), rk, dij, dik, &
                                      djk, lib, potential_parameter, int_abc_ext=sabc_ext)
                     IF (my_eps_screen > sabc_ext*(max_contra(iset, ikind)* &
                                                   max_contrb(jset, jkind)* &
                                                   max_contrc(kset, kkind))) THEN
                        DEALLOCATE (sabc)
                        CYCLE
                     END IF
                  ELSE
                     CALL eri_3center(sabc, la_min(iset), la_max(iset), npgfa(iset), zeta(:, iset), &
                                      rpgf_a(:, iset), ri, lb_min(jset), lb_max(jset), npgfb(jset), &
                                      zetb(:, jset), rpgf_b(:, jset), rj, lc_min(kset), lc_max(kset), &
                                      npgfc(kset), zetc(:, kset), rpgf_c(:, kset), rk, dij, dik, &
                                      djk, lib, potential_parameter)
                  END IF

                  ALLOCATE (work(nsgfa(iset), nsgfb(jset), nsgfc(kset)))

                  CALL libxsmm_abc_contract(work, sabc, tspa(iset, ikind)%array, spb(jset, jkind)%array, &
                                            spc(kset, kkind)%array, ncoa, ncob, ncoc, nsgfa(iset), &
                                            nsgfb(jset), nsgfc(kset), cpp_buffer, ccp_buffer)

                  iabc(sgfa:egfa, sgfb:egfb, sgfc:egfc) = work(:, :, :)
                  DEALLOCATE (sabc, work)

               END DO !kset
            END DO !jset
         END DO !iset

         !Add the integral to the proper tensor block
!$OMP CRITICAL
         CALL dbt_put_block(pq_X, [iatom, jatom, katom], SHAPE(iabc), iabc, summation=.TRUE.)
!$OMP END CRITICAL

         DEALLOCATE (iabc)
      END DO !o3c_iterator

      CALL cp_libint_cleanup_3eri(lib)

!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)
      CALL release_o3c_container(o3c)
      DEALLOCATE (o3c)

      DO iset = 1, max_nset
         DO ibasis = 1, nbasis
            IF (ASSOCIATED(tspa(iset, ibasis)%array)) DEALLOCATE (tspa(iset, ibasis)%array)
            IF (ASSOCIATED(spb(iset, ibasis)%array)) DEALLOCATE (spb(iset, ibasis)%array)
            IF (ASSOCIATED(spc(iset, ibasis)%array)) DEALLOCATE (spc(iset, ibasis)%array)
         END DO
      END DO
      DEALLOCATE (tspa, spb, spc)

      CALL timestop(handle)

   END SUBROUTINE fill_pqX_tensor

! **************************************************************************************************
!> \brief Builds a neighbor lists set for overlaping 2-center S_ab, where b is restricted on a
!>        a given list of atoms. Used for Coulomb RI where (aI|P) = sum_b C_bI (ab|P), where
!>        contraction coeff C_bI is assumed to be non-zero only on excited atoms
!> \param ab_list the neighbor list
!> \param basis_a basis set list for atom a
!> \param basis_b basis set list for atom b
!> \param qs_env ...
!> \param excited_atoms the indices of the excited atoms on which b is centered
!> \param ext_dist2d use an external distribution2d
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_ovlp_nl(ab_list, basis_a, basis_b, qs_env, excited_atoms, ext_dist2d)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_a, basis_b
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: excited_atoms
      TYPE(distribution_2d_type), OPTIONAL, POINTER      :: ext_dist2d

      INTEGER                                            :: ikind, nkind
      LOGICAL                                            :: my_restrictb
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: a_present, b_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a_radius, b_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      my_restrictb = .FALSE.
      IF (PRESENT(excited_atoms)) THEN
         my_restrictb = .TRUE.
      END IF

      ALLOCATE (a_present(nkind), b_present(nkind))
      a_present = .FALSE.
      b_present = .FALSE.
      ALLOCATE (a_radius(nkind), b_radius(nkind))
      a_radius = 0.0_dp
      b_radius = 0.0_dp

!  Set up the radii
      DO ikind = 1, nkind
         IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
            a_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
         END IF

         IF (ASSOCIATED(basis_b(ikind)%gto_basis_set)) THEN
            b_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_b(ikind)%gto_basis_set, kind_radius=b_radius(ikind))
         END IF
      END DO !ikind

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(a_present, b_present, a_radius, b_radius, pair_radius)

!  Set up the nl
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)

      !use an external distribution_2d if required
      IF (PRESENT(ext_dist2d)) distribution_2d => ext_dist2d

      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      IF (my_restrictb) THEN

         CALL build_neighbor_lists(ab_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   atomb_to_keep=excited_atoms, nlname="XAS_TDP_ovlp_nl")

      ELSE

         CALL build_neighbor_lists(ab_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   nlname="XAS_TDP_ovlp_nl")

      END IF
!  Clean-up
      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE build_xas_tdp_ovlp_nl

! **************************************************************************************************
!> \brief Builds a neighbor lists set taylored for 3-center integral within XAS TDP, such that only
!>        excited atoms are taken into account for the list_c
!> \param ac_list the neighbor list ready for 3-center integrals
!> \param basis_a basis set list for atom a
!> \param basis_c basis set list for atom c
!> \param op_type to indicate whther the list should be built with overlap, Coulomb or else in mind
!> \param qs_env ...
!> \param excited_atoms the indices of the excited atoms to consider (if not given, all atoms are taken)
!> \param x_range in case some truncated/screened operator is used, gives its range
!> \param ext_dist2d external distribution_2d to be used
!> \note Based on setup_neighbor_list with added features
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_3c_nl(ac_list, basis_a, basis_c, op_type, qs_env, excited_atoms, &
                                  x_range, ext_dist2d)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ac_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_a, basis_c
      INTEGER, INTENT(IN)                                :: op_type
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: excited_atoms
      REAL(dp), INTENT(IN), OPTIONAL                     :: x_range
      TYPE(distribution_2d_type), OPTIONAL, POINTER      :: ext_dist2d

      INTEGER                                            :: ikind, nkind
      LOGICAL                                            :: sort_atoms
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: a_present, c_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a_radius, c_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      sort_atoms = .FALSE.
      IF ((PRESENT(excited_atoms))) sort_atoms = .TRUE.

      ALLOCATE (a_present(nkind), c_present(nkind))
      a_present = .FALSE.
      c_present = .FALSE.
      ALLOCATE (a_radius(nkind), c_radius(nkind))
      a_radius = 0.0_dp
      c_radius = 0.0_dp

!  Set up the radii, depending on the operator type
      IF (op_type == do_potential_id) THEN

         !overlap => use the kind radius for both a and c
         DO ikind = 1, nkind
            !orbital basis set
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
               a_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
            END IF
            !RI_XAS basis set
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_c(ikind)%gto_basis_set, kind_radius=c_radius(ikind))
            END IF
         END DO !ikind

      ELSE IF (op_type == do_potential_coulomb) THEN

         !Coulomb operator, virtually infinite range => set c_radius to arbitrarily large number
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               c_radius(ikind) = 1000000.0_dp
            END IF
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) a_present(ikind) = .TRUE.
         END DO !ikind

      ELSE IF (op_type == do_potential_truncated .OR. op_type == do_potential_short) THEN

         !Truncated coulomb/short range: set c_radius to x_range + the kind_radii
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
               a_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
            END IF
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_c(ikind)%gto_basis_set, kind_radius=c_radius(ikind))
               c_radius(ikind) = c_radius(ikind) + x_range
            END IF
         END DO !ikind

      ELSE
         CPABORT("Operator not known")
      END IF

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(a_present, c_present, a_radius, c_radius, pair_radius)

!  Actually setup the list
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)

      !use an external distribution_2d if required
      IF (PRESENT(ext_dist2d)) distribution_2d => ext_dist2d

      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      IF (sort_atoms) THEN
         CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   operator_type="ABC", atomb_to_keep=excited_atoms, &
                                   nlname="XAS_TDP_3c_nl")
      ELSE
         CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   operator_type="ABC", nlname="XAS_TDP_3c_nl")
      END IF

!  Clean-up
      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE build_xas_tdp_3c_nl

! **************************************************************************************************
!> \brief Returns an optimized distribution_2d for the given neighbor lists based on an evaluation
!>        of the cost of the corresponding 3-center integrals
!> \param opt_3c_dist2d the optimized distribution_2d
!> \param ab_list ...
!> \param ac_list ...
!> \param basis_set_a ...
!> \param basis_set_b ...
!> \param basis_set_c ...
!> \param qs_env ...
!> \param only_bc_same_center ...
! **************************************************************************************************
   SUBROUTINE get_opt_3c_dist2d(opt_3c_dist2d, ab_list, ac_list, basis_set_a, basis_set_b, &
                                basis_set_c, qs_env, only_bc_same_center)

      TYPE(distribution_2d_type), POINTER                :: opt_3c_dist2d
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_list, ac_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_bc_same_center

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_opt_3c_dist2d'

      INTEGER                                            :: handle, i, iatom, ikind, ip, jatom, &
                                                            jkind, kkind, mypcol, myprow, n, &
                                                            natom, nkind, npcol, nprow, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nparticle_local_col, nparticle_local_row
      INTEGER, DIMENSION(:, :), POINTER                  :: col_dist, row_dist
      LOGICAL                                            :: my_sort_bc
      REAL(dp)                                           :: cost, rab(3), rac(3), rbc(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: col_cost, col_proc_cost, row_cost, &
                                                            row_proc_cost
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: local_particle_col, local_particle_row
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ab_iter, ac_iter
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (para_env, col_dist, row_dist, blacs_env, qs_kind_set, particle_set)
      NULLIFY (local_particle_col, local_particle_row, ab_iter, ac_iter)

      CALL timeset(routineN, handle)

      !Idea: create a,b and a,c nl_iterators in the original dist, then loop over them and compute the
      !      cost of each ab pairs and project on the col/row. Then distribute the atom col/row to
      !      the proc col/row in order to spread out the cost as much as possible

      my_sort_bc = .FALSE.
      IF (PRESENT(only_bc_same_center)) my_sort_bc = only_bc_same_center

      CALL get_qs_env(qs_env, natom=natom, para_env=para_env, blacs_env=blacs_env, &
                      qs_kind_set=qs_kind_set, particle_set=particle_set, nkind=nkind)

      myprow = blacs_env%mepos(1) + 1
      mypcol = blacs_env%mepos(2) + 1
      nprow = blacs_env%num_pe(1)
      npcol = blacs_env%num_pe(2)

      ALLOCATE (col_cost(natom), row_cost(natom))
      col_cost = 0.0_dp; row_cost = 0.0_dp

      ALLOCATE (row_dist(natom, 2), col_dist(natom, 2))
      row_dist = 1; col_dist = 1

      CALL neighbor_list_iterator_create(ab_iter, ab_list)
      CALL neighbor_list_iterator_create(ac_iter, ac_list, search=.TRUE.)
      DO WHILE (neighbor_list_iterate(ab_iter) == 0)
         CALL get_iterator_info(ab_iter, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rab)

         DO kkind = 1, nkind
            CALL nl_set_sub_iterator(ac_iter, ikind, kkind, iatom)

            DO WHILE (nl_sub_iterate(ac_iter) == 0)

               IF (my_sort_bc) THEN
                  !only take a,b,c if b,c (or a,c because of symmetry) share the same center
                  CALL get_iterator_info(ac_iter, r=rac)
                  rbc(:) = rac(:) - rab(:)
                  IF (.NOT. (ALL(ABS(rbc) .LE. 1.0E-8_dp) .OR. ALL(ABS(rac) .LE. 1.0E-8_dp))) CYCLE

               END IF

               !Use the size of integral as measure as contraciton cost seems to dominate
               nsgfa = basis_set_a(ikind)%gto_basis_set%nsgf
               nsgfb = basis_set_b(jkind)%gto_basis_set%nsgf
               nsgfc = basis_set_c(kkind)%gto_basis_set%nsgf

               cost = REAL(nsgfa*nsgfb*nsgfc, dp)

               row_cost(iatom) = row_cost(iatom) + cost
               col_cost(jatom) = col_cost(jatom) + cost

            END DO !ac_iter
         END DO !kkind
      END DO !ab_iter
      CALL neighbor_list_iterator_release(ab_iter)
      CALL neighbor_list_iterator_release(ac_iter)

      CALL mp_sum(row_cost, para_env%group)
      CALL mp_sum(col_cost, para_env%group)

      !Distribute the cost as evenly as possible
      ALLOCATE (col_proc_cost(npcol), row_proc_cost(nprow))
      col_proc_cost = 0.0_dp; row_proc_cost = 0.0_dp
      DO i = 1, natom
         iatom = MAXLOC(row_cost, 1)
         ip = MINLOC(row_proc_cost, 1)
         row_proc_cost(ip) = row_proc_cost(ip) + row_cost(iatom)
         row_dist(iatom, 1) = ip
         row_cost(iatom) = 0.0_dp

         iatom = MAXLOC(col_cost, 1)
         ip = MINLOC(col_proc_cost, 1)
         col_proc_cost(ip) = col_proc_cost(ip) + col_cost(iatom)
         col_dist(iatom, 1) = ip
         col_cost(iatom) = 0.0_dp
      END DO

      !the usual stuff
      ALLOCATE (local_particle_col(nkind), local_particle_row(nkind))
      ALLOCATE (nparticle_local_row(nkind), nparticle_local_col(nkind))
      nparticle_local_row = 0; nparticle_local_col = 0

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (row_dist(iatom, 1) == myprow) nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
         IF (col_dist(iatom, 1) == mypcol) nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
      END DO

      DO ikind = 1, nkind
         n = nparticle_local_row(ikind)
         ALLOCATE (local_particle_row(ikind)%array(n))

         n = nparticle_local_col(ikind)
         ALLOCATE (local_particle_col(ikind)%array(n))
      END DO

      nparticle_local_row = 0; nparticle_local_col = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (row_dist(iatom, 1) == myprow) THEN
            nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
            local_particle_row(ikind)%array(nparticle_local_row(ikind)) = iatom
         END IF
         IF (col_dist(iatom, 1) == mypcol) THEN
            nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
            local_particle_col(ikind)%array(nparticle_local_col(ikind)) = iatom
         END IF
      END DO

      !Finally create the dist_2d
      row_dist(:, 1) = row_dist(:, 1) - 1
      col_dist(:, 1) = col_dist(:, 1) - 1
      CALL distribution_2d_create(opt_3c_dist2d, row_distribution_ptr=row_dist, &
                                  col_distribution_ptr=col_dist, local_rows_ptr=local_particle_row, &
                                  local_cols_ptr=local_particle_col, blacs_env=blacs_env)

      CALL timestop(handle)

   END SUBROUTINE get_opt_3c_dist2d

! **************************************************************************************************
!> \brief Computes the RI exchange 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on excited atoms and kind. The operator used is that of the RI metric
!> \param ex_atoms excited atoms on which the third center is located
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  This routine is called once for each excited atom. Because there are many different a,b
!>        pairs involved, load balance is ok. This allows memory saving
! **************************************************************************************************
   SUBROUTINE compute_ri_3c_exchange(ex_atoms, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: ex_atoms
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_3c_exchange'

      INTEGER                                            :: handle, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_size_orb, blk_size_ri
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: opt_dbcsr_dist
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_list, ac_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_ri, basis_set_orb, ac_list, ab_list, qs_kind_set, para_env, particle_set)

      CALL timeset(routineN, handle)

!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, para_env=para_env, &
                      natom=natom, particle_set=particle_set)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Get the optimized distribution 2d for theses integrals (and store it in xas_tdp_env)
      CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env)
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, &
                               xas_tdp_control%ri_m_potential%potential_type, qs_env, &
                               excited_atoms=ex_atoms, x_range=xas_tdp_control%ri_m_potential%cutoff_radius)

      CALL get_opt_3c_dist2d(xas_tdp_env%opt_dist2d_ex, ab_list, ac_list, basis_set_orb, &
                             basis_set_orb, basis_set_ri, qs_env)
      CALL release_neighbor_list_sets(ab_list)
      CALL release_neighbor_list_sets(ac_list)

!  Build the ab and ac centers neighbor lists based on the optimized distribution
      CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env, &
                                 ext_dist2d=xas_tdp_env%opt_dist2d_ex)
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, &
                               xas_tdp_control%ri_m_potential%potential_type, qs_env, &
                               excited_atoms=ex_atoms, x_range=xas_tdp_control%ri_m_potential%cutoff_radius, &
                               ext_dist2d=xas_tdp_env%opt_dist2d_ex)

!  Allocate, init and compute the integrals.
      ALLOCATE (blk_size_orb(natom), blk_size_ri(natom))
      CALL cp_dbcsr_dist2d_to_dist(xas_tdp_env%opt_dist2d_ex, opt_dbcsr_dist)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_orb, basis=basis_set_orb)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_ri, basis=basis_set_ri)

      ALLOCATE (xas_tdp_env%ri_3c_ex)
      CALL create_pqX_tensor(xas_tdp_env%ri_3c_ex, ab_list, ac_list, opt_dbcsr_dist, blk_size_orb, &
                             blk_size_orb, blk_size_ri)
      CALL fill_pqX_tensor(xas_tdp_env%ri_3c_ex, ab_list, ac_list, basis_set_orb, basis_set_orb, &
                           basis_set_ri, xas_tdp_control%ri_m_potential, qs_env, &
                           eps_screen=xas_tdp_control%eps_screen)

! Clean-up
      CALL release_neighbor_list_sets(ab_list)
      CALL release_neighbor_list_sets(ac_list)
      CALL dbcsr_distribution_release(opt_dbcsr_dist)
      DEALLOCATE (basis_set_ri, basis_set_orb)

      !not strictly necessary but avoid having any load unbalance here being reported in the
      !timings for other routines
      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_3c_exchange

! **************************************************************************************************
!> \brief Computes the RI Coulomb 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on the excited atoms of xas_tdp_env
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note  The ri_3c_coul tensor of xas_tdp_env is defined and allocated here. Only computed once
!>        for the whole system (for optimized load balance). Ok because not too much memory needed
! **************************************************************************************************
   SUBROUTINE compute_ri_3c_coulomb(xas_tdp_env, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_3c_coulomb'

      INTEGER                                            :: handle, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_size_orb, blk_size_ri
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: opt_dbcsr_dist
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(libint_potential_type)                        :: pot
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_list, ac_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_ri, basis_set_orb, ac_list, ab_list, qs_kind_set, para_env, particle_set)

      CALL timeset(routineN, handle)

!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, para_env=para_env, &
                      natom=natom, particle_set=particle_set)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Get the optimized distribution 2d for these integrals (and store it in xas_tdp_env)
      CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env, &
                                 excited_atoms=xas_tdp_env%ex_atom_indices)
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, do_potential_id, &
                               qs_env, excited_atoms=xas_tdp_env%ex_atom_indices)
      CALL get_opt_3c_dist2d(xas_tdp_env%opt_dist2d_coul, ab_list, ac_list, basis_set_orb, &
                             basis_set_orb, basis_set_ri, qs_env, only_bc_same_center=.TRUE.)
      CALL release_neighbor_list_sets(ab_list)
      CALL release_neighbor_list_sets(ac_list)

!  Build a neighbor list for the ab centers. Assume (aI|c) = sum_b c_bI (ab|c), with c_bI only
!  non-zero for b centered on the same atom as c => build overlap nl, but only keeping b if centered
!  on an excited atom
      CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env, &
                                 excited_atoms=xas_tdp_env%ex_atom_indices, &
                                 ext_dist2d=xas_tdp_env%opt_dist2d_coul)

!  Build a neighbor list for the ac centers. Since we later contract as (aI|c) and we assume I is
!  very localized on the same atom as c, we take a,c as neighbors if they overlap
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, do_potential_id, &
                               qs_env, excited_atoms=xas_tdp_env%ex_atom_indices, &
                               ext_dist2d=xas_tdp_env%opt_dist2d_coul)

!  Allocate, init and compute the integrals
      ALLOCATE (blk_size_orb(natom), blk_size_ri(natom))
      CALL cp_dbcsr_dist2d_to_dist(xas_tdp_env%opt_dist2d_coul, opt_dbcsr_dist)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_orb, basis=basis_set_orb)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_ri, basis=basis_set_ri)
      pot%potential_type = do_potential_coulomb

      ALLOCATE (xas_tdp_env%ri_3c_coul)
      CALL create_pqX_tensor(xas_tdp_env%ri_3c_coul, ab_list, ac_list, opt_dbcsr_dist, blk_size_orb, &
                             blk_size_orb, blk_size_ri, only_bc_same_center=.TRUE.)
      CALL fill_pqX_tensor(xas_tdp_env%ri_3c_coul, ab_list, ac_list, basis_set_orb, basis_set_orb, &
                           basis_set_ri, pot, qs_env, only_bc_same_center=.TRUE.)

! Clean-up
      CALL release_neighbor_list_sets(ab_list)
      CALL release_neighbor_list_sets(ac_list)
      CALL dbcsr_distribution_release(opt_dbcsr_dist)
      DEALLOCATE (basis_set_ri, basis_set_orb)

      !not strictly necessary but avoid having any load unbalance here being reported in the
      !timings for other routines
      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_3c_coulomb

! **************************************************************************************************
!> \brief Computes the two-center Exchange integral needed for the RI in kernel calculation. Stores
!>        the integrals in the xas_tdp_env as global (small) arrays. Does that for a given excited
!>        kind. The quantity stored is M^-1 (P|Q) M^-1, where M is the RI metric. If the metric is
!>        the same as the exchange potential, then we end up with the V-approximation (P|Q)^-1
!>        By default (if no metric), the ri_m_potential is a copy of the x_potential
!> \param ex_kind ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Computes all these integrals in non-PBCs as we assume that the range is short enough that
!>       atoms do not exchange with their periodic images
! **************************************************************************************************
   SUBROUTINE compute_ri_exchange2_int(ex_kind, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: ex_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: egfp, egfq, maxl, ncop, ncoq, nset, &
                                                            nsgf, pset, qset, sgfp, sgfq, unit_id
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf_set, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp)                                           :: r(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: metric, pq, work
      REAL(dp), DIMENSION(:, :), POINTER                 :: rpgf, sphi, zet
      TYPE(cp_libint_t)                                  :: lib
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (ri_basis, qs_kind_set, para_env, lmin, lmax, npgf_set, zet, rpgf, first_sgf)
      NULLIFY (sphi, nsgf_set)

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, para_env=para_env)
      IF (ASSOCIATED(xas_tdp_env%ri_inv_ex)) THEN
         DEALLOCATE (xas_tdp_env%ri_inv_ex)
      END IF

!  Get the RI basis of interest and its quantum numbers
      CALL get_qs_kind(qs_kind_set(ex_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nsgf=nsgf, maxl=maxl, npgf=npgf_set, lmin=lmin, &
                             lmax=lmax, zet=zet, pgf_radius=rpgf, first_sgf=first_sgf, &
                             nsgf_set=nsgf_set, sphi=sphi, nset=nset)
      ALLOCATE (metric(nsgf, nsgf))
      metric = 0.0_dp

      r = 0.0_dp

      !init the libint 2-center object
      CALL cp_libint_init_2eri(lib, maxl)
      CALL cp_libint_set_contrdepth(lib, 1)

      !make sure the truncted coulomb is initialized
      IF (xas_tdp_control%ri_m_potential%potential_type == do_potential_truncated) THEN

         IF (2*maxl + 1 > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=xas_tdp_control%ri_m_potential%filename)
            END IF
            CALL init(2*maxl + 1, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

!  Compute the RI metric
      DO pset = 1, nset
         ncop = npgf_set(pset)*ncoset(lmax(pset))
         sgfp = first_sgf(1, pset)
         egfp = sgfp + nsgf_set(pset) - 1

         DO qset = 1, nset
            ncoq = npgf_set(qset)*ncoset(lmax(qset))
            sgfq = first_sgf(1, qset)
            egfq = sgfq + nsgf_set(qset) - 1

            ALLOCATE (work(ncop, ncoq))
            work = 0.0_dp

            CALL eri_2center(work, lmin(pset), lmax(pset), npgf_set(pset), zet(:, pset), rpgf(:, pset), &
                             r, lmin(qset), lmax(qset), npgf_set(qset), zet(:, qset), rpgf(:, qset), &
                             r, 0.0_dp, lib, xas_tdp_control%ri_m_potential)

            CALL ab_contract(metric(sgfp:egfp, sgfq:egfq), work, sphi(:, sgfp:), sphi(:, sgfq:), &
                             ncop, ncoq, nsgf_set(pset), nsgf_set(qset))

            DEALLOCATE (work)
         END DO !qset
      END DO !pset

!  Inverting (to M^-1)
      CALL invmat_symm(metric)

      IF (.NOT. xas_tdp_control%do_ri_metric) THEN

         !If no metric, then x_pot = ri_m_pot and (P|Q)^-1 = M^-1 (V-approximation)
         ALLOCATE (xas_tdp_env%ri_inv_ex(nsgf, nsgf))
         xas_tdp_env%ri_inv_ex(:, :) = metric(:, :)
         CALL cp_libint_cleanup_2eri(lib)
         RETURN

      END IF

      !make sure the truncted coulomb is initialized
      IF (xas_tdp_control%x_potential%potential_type == do_potential_truncated) THEN

         IF (2*maxl + 1 > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=xas_tdp_control%x_potential%filename)
            END IF
            CALL init(2*maxl + 1, unit_id, para_env%mepos, para_env%group)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

!  Compute the proper exchange 2-center
      ALLOCATE (pq(nsgf, nsgf))
      pq = 0.0_dp

      DO pset = 1, nset
         ncop = npgf_set(pset)*ncoset(lmax(pset))
         sgfp = first_sgf(1, pset)
         egfp = sgfp + nsgf_set(pset) - 1

         DO qset = 1, nset
            ncoq = npgf_set(qset)*ncoset(lmax(qset))
            sgfq = first_sgf(1, qset)
            egfq = sgfq + nsgf_set(qset) - 1

            ALLOCATE (work(ncop, ncoq))
            work = 0.0_dp

            CALL eri_2center(work, lmin(pset), lmax(pset), npgf_set(pset), zet(:, pset), rpgf(:, pset), &
                             r, lmin(qset), lmax(qset), npgf_set(qset), zet(:, qset), rpgf(:, qset), &
                             r, 0.0_dp, lib, xas_tdp_control%x_potential)

            CALL ab_contract(pq(sgfp:egfp, sgfq:egfq), work, sphi(:, sgfp:), sphi(:, sgfq:), &
                             ncop, ncoq, nsgf_set(pset), nsgf_set(qset))

            DEALLOCATE (work)
         END DO !qset
      END DO !pset

!  Compute and store M^-1 (P|Q) M^-1
      ALLOCATE (xas_tdp_env%ri_inv_ex(nsgf, nsgf))
      xas_tdp_env%ri_inv_ex = 0.0_dp

      CALL dgemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, metric, nsgf, pq, nsgf, &
                 0.0_dp, xas_tdp_env%ri_inv_ex, nsgf)
      CALL dgemm('N', 'N', nsgf, nsgf, nsgf, 1.0_dp, xas_tdp_env%ri_inv_ex, nsgf, metric, nsgf, &
                 0.0_dp, pq, nsgf)
      xas_tdp_env%ri_inv_ex(:, :) = pq(:, :)

      CALL cp_libint_cleanup_2eri(lib)

   END SUBROUTINE compute_ri_exchange2_int

! **************************************************************************************************
!> \brief Computes the two-center Coulomb integral needed for the RI in kernel calculation. Stores
!>        the integrals (P|Q)^-1 in the xas_tdp_env as global (small) arrays. Does that for a given
!>        excited kind
!> \param ex_kind ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb2_int(ex_kind, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: ex_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: nsgf
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (ri_basis, qs_kind_set)

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      IF (ASSOCIATED(xas_tdp_env%ri_inv_coul)) THEN
         DEALLOCATE (xas_tdp_env%ri_inv_coul)
      END IF

!  Get the RI basis of interest and its quantum numbers
      CALL get_qs_kind(qs_kind_set(ex_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nsgf=nsgf)
      ALLOCATE (xas_tdp_env%ri_inv_coul(nsgf, nsgf))
      xas_tdp_env%ri_inv_coul = 0.0_dp

      IF (.NOT. xas_tdp_control%is_periodic) THEN
         CALL int_operators_r12_ab_os(r12_operator=operator_coulomb, vab=xas_tdp_env%ri_inv_coul, &
                                      rab=(/0.0_dp, 0.0_dp, 0.0_dp/), fba=ri_basis, fbb=ri_basis, &
                                      calculate_forces=.FALSE.)
         CPASSERT(ASSOCIATED(xas_tdp_control))
      ELSE
         CALL periodic_ri_coulomb2(xas_tdp_env%ri_inv_coul, ri_basis, qs_env)
      END IF

!  Inverting
      CALL invmat_symm(xas_tdp_env%ri_inv_coul)

   END SUBROUTINE compute_ri_coulomb2_int

! **************************************************************************************************
!> \brief Computes the two-center inverse coulomb integral in the case of PBCs
!> \param ri_coul2 ...
!> \param ri_basis ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE periodic_ri_coulomb2(ri_coul2, ri_basis, qs_env)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: ri_coul2
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: maxco, ncop, ncoq, nset, op, oq, ppgf, &
                                                            pset, qpgf, qset, sgfp, sgfq
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp)                                           :: r(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: hpq
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi, zet
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_eri_mme_param)                             :: mme_param
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (cell, qs_kind_set, lmin, lmax, nsgf, npgf, zet, sphi, first_sgf)

      ! Use eri_mme for this. Don't want to annoy the user with a full input section just for this
      ! tiny bit => initialize our own eri_mme section with the defaults

      CALL get_qs_env(qs_env, cell=cell, qs_kind_set=qs_kind_set, para_env=para_env)

      CALL eri_mme_init(mme_param%par, n_minimax=20, cutoff=300.0_dp, do_calib_cutoff=.TRUE., &
                        cutoff_min=10.0_dp, cutoff_max=10000.0_dp, cutoff_eps=0.01_dp, &
                        cutoff_delta=0.9_dp, sum_precision=1.0E-12_dp, debug=.FALSE., &
                        debug_delta=1.0E-6_dp, debug_nsum=1000000, unit_nr=0, print_calib=.FALSE., &
                        do_error_est=.FALSE.)
      mme_param%do_calib = .TRUE.

      CALL cp_eri_mme_set_params(mme_param, cell, qs_kind_set, basis_type_1="RI_XAS", para_env=para_env)

      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, lmin=lmin, nset=nset, &
                             nsgf_set=nsgf, sphi=sphi, first_sgf=first_sgf, maxco=maxco)

      r = 0.0_dp
      ALLOCATE (hpq(nset*maxco, nset*maxco))
      hpq = 0.0_dp

      DO pset = 1, nset
         ncop = npgf(pset)*ncoset(lmax(pset))
         sgfp = first_sgf(1, pset)

         DO qset = 1, nset
            ncoq = npgf(qset)*ncoset(lmax(qset))
            sgfq = first_sgf(1, qset)

            DO ppgf = 1, npgf(pset)
               op = (pset - 1)*maxco + (ppgf - 1)*ncoset(lmax(pset))
               DO qpgf = 1, npgf(qset)
                  oq = (qset - 1)*maxco + (qpgf - 1)*ncoset(lmax(qset))

                  CALL eri_mme_2c_integrate(mme_param%par, lmin(pset), lmax(pset), lmin(qset), &
                                            lmax(qset), zet(ppgf, pset), zet(qpgf, qset), r, hpq, &
                                            op, oq)

               END DO !qpgf
            END DO ! ppgf

            !contraction into sgfs
            op = (pset - 1)*maxco + 1
            oq = (qset - 1)*maxco + 1

            CALL ab_contract(ri_coul2(sgfp:sgfp + nsgf(pset) - 1, sgfq:sgfq + nsgf(qset) - 1), &
                             hpq(op:op + ncop - 1, oq:oq + ncoq - 1), sphi(:, sgfp:), sphi(:, sgfq:), &
                             ncop, ncoq, nsgf(pset), nsgf(qset))

         END DO !qset
      END DO !pset

      !celan-up
      CALL eri_mme_release(mme_param%par)

   END SUBROUTINE periodic_ri_coulomb2

END MODULE xas_tdp_integrals
