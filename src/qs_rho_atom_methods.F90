!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_rho_atom_methods

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              gapw_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE gaunt_coefficients,              ONLY: real_gaunt_coefficients
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: create_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: create_harmonics_atom,&
                                              get_maxl_CG,&
                                              get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              neighbor_list_task_type
   USE qs_oce_methods,                  ONLY: proj_blk
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE sap_kind_types,                  ONLY: alist_pre_align_blk,&
                                              alist_type,&
                                              clist_type,&
                                              get_alist
   USE util,                            ONLY: get_limit
   USE whittaker,                       ONLY: whittaker_c0a,&
                                              whittaker_ci

!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_atom_methods'
   INTEGER, PARAMETER, PRIVATE :: l2_chunk_dp = 16000

! *** Public subroutines ***

   PUBLIC :: allocate_rho_atom_internals, &
             calculate_rho_atom, &
             calculate_rho_atom_coeff, &
             init_rho_atom, &
             replicate_rho_atom_radial

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param rho_atom_set ...
!> \param qs_kind ...
!> \param atom_list ...
!> \param natom ...
!> \param nspins ...
!> \param tot_rho1_h ...
!> \param tot_rho1_s ...
!> \param rho1_h_spin ...
!> \param rho1_s_spin ...
!> \param rho1_h_aspin ...
!> \param rho1_s_aspin ...
! **************************************************************************************************
   SUBROUTINE calculate_rho_atom(para_env, rho_atom_set, qs_kind, atom_list, &
                                 natom, nspins, tot_rho1_h, tot_rho1_s, &
                                 rho1_h_spin, rho1_s_spin, rho1_h_aspin, rho1_s_aspin)

      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_list
      INTEGER, INTENT(IN)                                :: natom, nspins
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: tot_rho1_h, tot_rho1_s
      REAL(dp), INTENT(INOUT)                            :: rho1_h_spin, rho1_s_spin, rho1_h_aspin, &
                                                            rho1_s_aspin

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_atom'

      INTEGER :: damax_iso_not0_local, handle, i, i1, i2, iat, iatom, icg, ipgf1, ipgf2, ir, &
         iset1, iset2, iso, iso1, iso1_coeff, iso1_first, iso1_last, iso2, iso2_coeff, iso2_first, &
         iso2_last, j, l, l_iso, l_sub, l_sum, lmax12, lmax_expansion, lmin12, m1s, m2s, &
         max_iso_not0, max_iso_not0_local, max_npgf, max_s_harm, maxl, maxso, mepos, n1s, n2s, na, &
         nr, nset, num_pe, size1, size2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list, dacg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list, dacg_list
      INTEGER, DIMENSION(2)                              :: bo
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, o2nindex
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: done_vgg
      REAL(dp)                                           :: c1, c2, cpc_h, cpc_s, rfun, rho_h, &
                                                            rho_s, root_zet12, zet12
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: erf_zet12, g1, g2, gg0, int1, int2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dgg, gg, gg_lm1, sfun_h, sfun_s
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: g_rad, vgg
      REAL(dp), DIMENSION(:, :), POINTER                 :: coeff_h, coeff_s, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: my_CG_dxyz
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL timeset(routineN, handle)

      !Note: tau is taken care of separately in qs_vxc_atom.F

      NULLIFY (basis_1c)
      NULLIFY (harmonics, grid_atom)
      NULLIFY (lmin, lmax, npgf, zet, my_CG, my_CG_dxyz, coeff_h, coeff_s)

      CALL get_qs_kind(qs_kind, grid_atom=grid_atom, harmonics=harmonics)
      CALL get_qs_kind(qs_kind, basis_set=basis_1c, basis_type="GAPW_1C")

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                             maxl=maxl, npgf=npgf, nset=nset, zet=zet, &
                             maxso=maxso)

      CALL get_paw_basis_info(basis_1c, o2nindex=o2nindex)

      max_iso_not0 = harmonics%max_iso_not0
      max_s_harm = harmonics%max_s_harm

      nr = grid_atom%nr
      max_npgf = MAXVAL(npgf(1:nset))
      lmax_expansion = indso(1, max_iso_not0)
      ! Distribute the atoms of this kind
      num_pe = para_env%num_pe
      mepos = para_env%mepos
      bo = get_limit(natom, num_pe, mepos)

      my_CG => harmonics%my_CG
      my_CG_dxyz => harmonics%my_CG_dxyz

      ALLOCATE (g1(nr), g2(nr), gg0(nr), gg(nr, 0:2*maxl), dgg(nr, 0:2*maxl), gg_lm1(nr, 0:2*maxl))
      ALLOCATE (erf_zet12(nr), vgg(nr, 0:2*maxl, 0:indso(1, max_iso_not0)))
      ALLOCATE (done_vgg(0:2*maxl, 0:indso(1, max_iso_not0)))
      ALLOCATE (int1(nr), int2(nr))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm), &
                dacg_list(2, nsoset(maxl)**2, max_s_harm), dacg_n_list(max_s_harm))
      ALLOCATE (g_rad(nr, max_npgf, nset))

      DO iset1 = 1, nset
         DO ipgf1 = 1, npgf(iset1)
            g_rad(1:nr, ipgf1, iset1) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
         END DO
      END DO

      DO iat = bo(1), bo(2)
         iatom = atom_list(iat)
         DO i = 1, nspins
            IF (.NOT. ASSOCIATED(rho_atom_set(iatom)%rho_rad_h(i)%r_coef)) THEN
               CALL allocate_rho_atom_rad(rho_atom_set, iatom, i, nr, max_iso_not0)
            ELSE
               CALL set2zero_rho_atom_rad(rho_atom_set, iatom, i)
            END IF
         END DO
      END DO

      m1s = 0
      DO iset1 = 1, nset
         m2s = 0
         DO iset2 = 1, nset

            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
            CPASSERT(max_iso_not0_local <= max_iso_not0)
            CALL get_none0_cg_list(my_CG_dxyz, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, dacg_list, dacg_n_list, damax_iso_not0_local)
            n1s = nsoset(lmax(iset1))

            DO ipgf1 = 1, npgf(iset1)
               iso1_first = nsoset(lmin(iset1) - 1) + 1 + n1s*(ipgf1 - 1) + m1s
               iso1_last = nsoset(lmax(iset1)) + n1s*(ipgf1 - 1) + m1s
               size1 = iso1_last - iso1_first + 1
               iso1_first = o2nindex(iso1_first)
               iso1_last = o2nindex(iso1_last)
               i1 = iso1_last - iso1_first + 1
               CPASSERT(size1 == i1)
               i1 = nsoset(lmin(iset1) - 1) + 1

               g1(1:nr) = g_rad(1:nr, ipgf1, iset1)

               n2s = nsoset(lmax(iset2))
               DO ipgf2 = 1, npgf(iset2)
                  iso2_first = nsoset(lmin(iset2) - 1) + 1 + n2s*(ipgf2 - 1) + m2s
                  iso2_last = nsoset(lmax(iset2)) + n2s*(ipgf2 - 1) + m2s
                  size2 = iso2_last - iso2_first + 1
                  iso2_first = o2nindex(iso2_first)
                  iso2_last = o2nindex(iso2_last)
                  i2 = iso2_last - iso2_first + 1
                  CPASSERT(size2 == i2)
                  i2 = nsoset(lmin(iset2) - 1) + 1

                  g2(1:nr) = g_rad(1:nr, ipgf2, iset2)
                  lmin12 = lmin(iset1) + lmin(iset2)
                  lmax12 = lmax(iset1) + lmax(iset2)

                  zet12 = zet(ipgf1, iset1) + zet(ipgf2, iset2)
                  root_zet12 = SQRT(zet(ipgf1, iset1) + zet(ipgf2, iset2))
                  DO ir = 1, nr
                     erf_zet12(ir) = erf(root_zet12*grid_atom%rad(ir))
                  END DO

                  gg = 0.0_dp
                  dgg = 0.0_dp
                  gg_lm1 = 0.0_dp
                  vgg = 0.0_dp
                  done_vgg = .FALSE.
                  ! reduce the number of terms in the expansion local densities
                  IF (lmin12 <= lmax_expansion) THEN
                     IF (lmin12 == 0) THEN
                        gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                        gg_lm1(1:nr, lmin12) = 0.0_dp
                        gg0(1:nr) = gg(1:nr, lmin12)
                     ELSE
                        gg0(1:nr) = g1(1:nr)*g2(1:nr)
                        gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                        gg_lm1(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12 - 1)*g1(1:nr)*g2(1:nr)
                     END IF

                     ! reduce the number of terms in the expansion local densities
                     IF (lmax12 > lmax_expansion) lmax12 = lmax_expansion

                     DO l = lmin12 + 1, lmax12
                        gg(1:nr, l) = grid_atom%rad(1:nr)*gg(1:nr, l - 1)
                        gg_lm1(1:nr, l) = gg(1:nr, l - 1)
                        dgg(1:nr, l - 1) = -2.0_dp*(zet(ipgf1, iset1) + zet(ipgf2, iset2))*gg(1:nr, l)

                     END DO
                     dgg(1:nr, lmax12) = -2.0_dp*(zet(ipgf1, iset1) + &
                                                  zet(ipgf2, iset2))*grid_atom%rad(1:nr)*gg(1:nr, lmax12)

                     c2 = SQRT(pi*pi*pi/(zet12*zet12*zet12))

                     DO iso = 1, max_iso_not0_local
                        l_iso = indso(1, iso)
                        c1 = fourpi/(2._dp*REAL(l_iso, dp) + 1._dp)
                        DO icg = 1, cg_n_list(iso)
                           iso1 = cg_list(1, icg, iso)
                           iso2 = cg_list(2, icg, iso)

                           l = indso(1, iso1) + indso(1, iso2)
                           CPASSERT(l <= lmax_expansion)
                           IF (done_vgg(l, l_iso)) CYCLE
                           L_sum = l + l_iso
                           L_sub = l - l_iso

                           IF (l_sum == 0) THEN
                              vgg(1:nr, l, l_iso) = erf_zet12(1:nr)*grid_atom%oorad2l(1:nr, 1)*c2
                           ELSE
                              CALL whittaker_c0a(int1, grid_atom%rad, gg0, erf_zet12, zet12, l, l_iso, nr)
                              CALL whittaker_ci(int2, grid_atom%rad, gg0, zet12, L_sub, nr)

                              DO ir = 1, nr
                                 int2(ir) = grid_atom%rad2l(ir, l_iso)*int2(ir)
                                 vgg(ir, l, l_iso) = c1*(int1(ir) + int2(ir))
                              END DO
                           END IF
                           done_vgg(l, l_iso) = .TRUE.
                        END DO
                     END DO
                  END IF ! lmax_expansion

                  DO iat = bo(1), bo(2)
                     iatom = atom_list(iat)

                     DO i = 1, nspins
                        coeff_h => rho_atom_set(iatom)%cpc_h(i)%r_coef
                        coeff_s => rho_atom_set(iatom)%cpc_s(i)%r_coef

                        DO iso = 1, max_iso_not0_local
                           l_iso = indso(1, iso)
                           DO icg = 1, cg_n_list(iso)
                              iso1 = cg_list(1, icg, iso)
                              iso2 = cg_list(2, icg, iso)

                              l = indso(1, iso1) + indso(1, iso2)
                              CPASSERT(l <= lmax_expansion)
                              iso1_coeff = iso1_first + iso1 - i1
                              iso2_coeff = iso2_first + iso2 - i2
                              cpc_h = coeff_h(iso1_coeff, iso2_coeff)*my_CG(iso1, iso2, iso)
                              cpc_s = coeff_s(iso1_coeff, iso2_coeff)*my_CG(iso1, iso2, iso)

                              rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr, iso) + &
                                 gg(1:nr, l)*cpc_h

                              rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr, iso) + &
                                 gg(1:nr, l)*cpc_s

                              rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr, iso) + &
                                 dgg(1:nr, l)*cpc_h

                              rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr, iso) + &
                                 dgg(1:nr, l)*cpc_s

                              rho_atom_set(iatom)%vrho_rad_h(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%vrho_rad_h(i)%r_coef(1:nr, iso) + &
                                 vgg(1:nr, l, l_iso)*cpc_h

                              rho_atom_set(iatom)%vrho_rad_s(i)%r_coef(1:nr, iso) = &
                                 rho_atom_set(iatom)%vrho_rad_s(i)%r_coef(1:nr, iso) + &
                                 vgg(1:nr, l, l_iso)*cpc_s

                           END DO ! icg

                        END DO ! iso

                        DO iso = 1, max_iso_not0 !damax_iso_not0_local
                           l_iso = indso(1, iso)
                           DO icg = 1, dacg_n_list(iso)
                              iso1 = dacg_list(1, icg, iso)
                              iso2 = dacg_list(2, icg, iso)
                              l = indso(1, iso1) + indso(1, iso2)
                              CPASSERT(l <= lmax_expansion)
                              iso1_coeff = iso1_first + iso1 - i1
                              iso2_coeff = iso2_first + iso2 - i2
                              cpc_h = coeff_h(iso1_coeff, iso2_coeff)
                              cpc_s = coeff_s(iso1_coeff, iso2_coeff)
                              DO j = 1, 3
                                 rho_atom_set(iatom)%rho_rad_h_d(j, i)%r_coef(1:nr, iso) = &
                                    rho_atom_set(iatom)%rho_rad_h_d(j, i)%r_coef(1:nr, iso) + &
                                    gg_lm1(1:nr, l)*cpc_h*my_CG_dxyz(j, iso1, iso2, iso)

                                 rho_atom_set(iatom)%rho_rad_s_d(j, i)%r_coef(1:nr, iso) = &
                                    rho_atom_set(iatom)%rho_rad_s_d(j, i)%r_coef(1:nr, iso) + &
                                    gg_lm1(1:nr, l)*cpc_s*my_CG_dxyz(j, iso1, iso2, iso)
                              END DO
                           END DO ! icg

                        END DO ! iso

                     END DO ! i
                  END DO ! iat

               END DO ! ipgf2
            END DO ! ipgf1
            m2s = m2s + maxso
         END DO ! iset2
         m1s = m1s + maxso
      END DO ! iset1

      DO iat = bo(1), bo(2)
         iatom = atom_list(iat)

         DO i = 1, nspins

            DO iso = 1, max_iso_not0
               rho_s = 0.0_dp
               rho_h = 0.0_dp
               DO ir = 1, nr
                  rho_h = rho_h + rho_atom_set(iatom)%rho_rad_h(i)%r_coef(ir, iso)*grid_atom%wr(ir)
                  rho_s = rho_s + rho_atom_set(iatom)%rho_rad_s(i)%r_coef(ir, iso)*grid_atom%wr(ir)
               END DO ! ir
               tot_rho1_h(i) = tot_rho1_h(i) + rho_h*harmonics%slm_int(iso)
               tot_rho1_s(i) = tot_rho1_s(i) + rho_s*harmonics%slm_int(iso)
            END DO ! iso

         END DO ! ispin

         IF (nspins == 2) THEN
            na = SIZE(harmonics%slm, 1)
            ALLOCATE (sfun_h(nr, na), sfun_s(nr, na))
            sfun_h = 0.0_dp
            sfun_s = 0.0_dp
            DO iso = 1, max_iso_not0
               DO ir = 1, nr
                  rfun = grid_atom%wr(ir)*(rho_atom_set(iatom)%rho_rad_h(1)%r_coef(ir, iso) - &
                                           rho_atom_set(iatom)%rho_rad_h(2)%r_coef(ir, iso))
                  sfun_h(ir, 1:na) = sfun_h(ir, 1:na) + rfun*harmonics%slm(1:na, iso)*grid_atom%wa(1:na)
                  rfun = grid_atom%wr(ir)*(rho_atom_set(iatom)%rho_rad_s(1)%r_coef(ir, iso) - &
                                           rho_atom_set(iatom)%rho_rad_s(2)%r_coef(ir, iso))
                  sfun_s(ir, 1:na) = sfun_s(ir, 1:na) + rfun*harmonics%slm(1:na, iso)*grid_atom%wa(1:na)
               END DO
            END DO
            rho1_h_spin = rho1_h_spin + SUM(sfun_h(1:nr, 1:na))
            rho1_s_spin = rho1_s_spin + SUM(sfun_s(1:nr, 1:na))
            rho1_h_aspin = rho1_h_aspin + SUM(ABS(sfun_h(1:nr, 1:na)))
            rho1_s_aspin = rho1_s_aspin + SUM(ABS(sfun_s(1:nr, 1:na)))
            DEALLOCATE (sfun_h, sfun_s)
         END IF

      END DO ! iat

      DEALLOCATE (g1, g2, gg0, gg, gg_lm1, dgg, vgg, done_vgg, erf_zet12, int1, int2, g_rad)
      DEALLOCATE (cg_list, cg_n_list, dacg_list, dacg_n_list)
      DEALLOCATE (o2nindex)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_atom

! **************************************************************************************************
!> \brief Replicate the radial hard/soft density data needed to evaluate one-center tails on
!>        rank-local target grids. The compact one-center density matrices are already global;
!>        this routine performs one packed reduction for the derived radial fields of a kind.
!> \param para_env ...
!> \param rho_atom_set ...
!> \param qs_kind ...
!> \param atom_list ...
!> \param natom ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE replicate_rho_atom_radial(para_env, rho_atom_set, qs_kind, atom_list, natom, nspins)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_list
      INTEGER, INTENT(IN)                                :: natom, nspins

      CHARACTER(len=*), PARAMETER :: routineN = 'replicate_rho_atom_radial'

      INTEGER                                            :: block_size, bo(2), cursor, handle, iat, &
                                                            iatom, ispin, j, max_iso_not0, ncoeff, &
                                                            nr
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: buffer
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL timeset(routineN, handle)

      NULLIFY (grid_atom, harmonics)
      CALL get_qs_kind(qs_kind, grid_atom=grid_atom, harmonics=harmonics)
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      nr = grid_atom%nr
      max_iso_not0 = harmonics%max_iso_not0
      ncoeff = nr*max_iso_not0
      block_size = 10*nspins*ncoeff
      ALLOCATE (buffer(block_size*natom))
      buffer = 0.0_dp

      bo = get_limit(natom, para_env%num_pe, para_env%mepos)
      DO iat = bo(1), bo(2)
         iatom = atom_list(iat)
         cursor = (iat - 1)*block_size
         DO ispin = 1, nspins
            buffer(cursor + 1:cursor + ncoeff) = &
               RESHAPE(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef, [ncoeff])
            cursor = cursor + ncoeff
            buffer(cursor + 1:cursor + ncoeff) = &
               RESHAPE(rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef, [ncoeff])
            cursor = cursor + ncoeff
            buffer(cursor + 1:cursor + ncoeff) = &
               RESHAPE(rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef, [ncoeff])
            cursor = cursor + ncoeff
            buffer(cursor + 1:cursor + ncoeff) = &
               RESHAPE(rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef, [ncoeff])
            cursor = cursor + ncoeff
            DO j = 1, 3
               buffer(cursor + 1:cursor + ncoeff) = &
                  RESHAPE(rho_atom_set(iatom)%rho_rad_h_d(j, ispin)%r_coef, [ncoeff])
               cursor = cursor + ncoeff
               buffer(cursor + 1:cursor + ncoeff) = &
                  RESHAPE(rho_atom_set(iatom)%rho_rad_s_d(j, ispin)%r_coef, [ncoeff])
               cursor = cursor + ncoeff
            END DO
         END DO
         CPASSERT(cursor == iat*block_size)
      END DO
      CALL para_env%sum(buffer)

      DO iat = 1, natom
         iatom = atom_list(iat)
         cursor = (iat - 1)*block_size
         DO ispin = 1, nspins
            IF (.NOT. ASSOCIATED(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef)) THEN
               CALL allocate_rho_atom_rad(rho_atom_set, iatom, ispin, nr, max_iso_not0)
            END IF
            rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef = &
               RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
            cursor = cursor + ncoeff
            rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef = &
               RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
            cursor = cursor + ncoeff
            rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = &
               RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
            cursor = cursor + ncoeff
            rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = &
               RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
            cursor = cursor + ncoeff
            DO j = 1, 3
               rho_atom_set(iatom)%rho_rad_h_d(j, ispin)%r_coef = &
                  RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
               cursor = cursor + ncoeff
               rho_atom_set(iatom)%rho_rad_s_d(j, ispin)%r_coef = &
                  RESHAPE(buffer(cursor + 1:cursor + ncoeff), [nr, max_iso_not0])
               cursor = cursor + ncoeff
            END DO
         END DO
         CPASSERT(cursor == iat*block_size)
      END DO
      DEALLOCATE (buffer)

      CALL timestop(handle)

   END SUBROUTINE replicate_rho_atom_radial

! **************************************************************************************************
!> \brief ...
!> \param qs_env        QuickStep environment
!>                      (accessed components: atomic_kind_set, dft_control%nimages,
!>                                            dft_control%nspins, kpoints%cell_to_index)
!> \param rho_ao        density matrix in atomic basis set
!> \param rho_atom_set ...
!> \param qs_kind_set   list of QuickStep kinds
!> \param oce           one-centre expansion coefficients
!> \param sab           neighbour pair list
!> \param para_env      parallel environment
!> \par History
!>      Add OpenMP [Apr 2016, EPCC]
!>      Use automatic arrays [Sep 2016, M Tucker]
!>      Allow for external non-default kind_set, oce and sab [Dec 2019, A Bussy]
!> \note  Consider to declare 'rho_ao' dummy argument as a pointer to the two-dimensional
!>        (1:nspins, 1:nimages) set of matrices.
! **************************************************************************************************
   SUBROUTINE calculate_rho_atom_coeff(qs_env, rho_ao, rho_atom_set, qs_kind_set, oce, sab, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(*)                   :: rho_ao
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_atom_coeff'

      INTEGER :: bo(2), handle, i, iac, iatom, ibc, ikind, img, inb, ispin, jatom, jkind, kac, &
         kac_hi, kac_lo, katom, kbc, kkind, len_PC1, max_gau, max_nb, max_nclist, max_nsgf, mepos, &
         n_cont_a, n_cont_b, n_nb, nat_kind, natom, nimages, nkind, nl_size, nsoctot, nspins, &
         num_pe, slot, w_alloc, w_total
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ia_ptr, ia_work, nsatbas_kind, &
                                                            slot_order, w_offset
      INTEGER, DIMENSION(3)                              :: cell_b
      INTEGER, DIMENSION(:), POINTER                     :: a_list, list_a, list_b
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: dista, distb, found, paw_atom
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: has_intac, paw_kind
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: w_touched
      REAL(KIND=dp)                                      :: beta, eps_cpc, factor, pmax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nb_pmax, proj_work1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: p_matrix, proj_work2, w_buf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: C_coeff_hh_a, C_coeff_hh_b, &
                                                            C_coeff_ss_a, C_coeff_ss_b, r_coef_h, &
                                                            r_coef_s
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(clist_type), POINTER                          :: clist_a, clist_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c, basis_set_a
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_task_type), POINTER             :: task
      TYPE(rho_atom_coeff), ALLOCATABLE, DIMENSION(:, :) :: nb_pblk

!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock, number_of_locks

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set)

      eps_cpc = dft_control%qs_control%gapw_control%eps_cpc

      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(rho_atom_set))
      CPASSERT(ASSOCIATED(oce))
      CPASSERT(ASSOCIATED(sab))

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      NULLIFY (cell_to_index)
      IF (nimages > 1) THEN
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, maxsgf=max_nsgf, maxgtops=max_gau, basis_type='GAPW_1C')

      nkind = SIZE(atomic_kind_set)
      !   Initialize to 0 the CPC coefficients and the local density arrays
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=a_list, natom=nat_kind)
         CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)

         IF (.NOT. paw_atom) CYCLE
         DO i = 1, nat_kind
            iatom = a_list(i)
            DO ispin = 1, nspins
               rho_atom_set(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
               rho_atom_set(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
            END DO ! ispin
         END DO ! i

         num_pe = para_env%num_pe
         mepos = para_env%mepos
         bo = get_limit(nat_kind, num_pe, mepos)
         DO i = bo(1), bo(2)
            iatom = a_list(i)
            DO ispin = 1, nspins
               rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
               rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
            END DO ! ispin
         END DO ! i
      END DO ! ikind

      ALLOCATE (basis_set_list(nkind))
      ALLOCATE (paw_kind(nkind), nsatbas_kind(nkind), has_intac(nkind*nkind))
      paw_kind(:) = .FALSE.
      nsatbas_kind(:) = 0
      has_intac(:) = .FALSE.
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C", &
                          paw_atom=paw_kind(ikind))
         IF (paw_kind(ikind)) CALL get_paw_basis_info(basis_1c, nsatbas=nsatbas_kind(ikind))
      END DO
      DO ikind = 1, nkind*nkind
         has_intac(ikind) = ASSOCIATED(oce%intac(ikind)%alist)
      END DO

      len_PC1 = MAX(1, max_nsgf*max_gau)
      w_alloc = MAX(1, l2_chunk_dp/nspins, len_PC1) ! 16k-double L2 chunk across spins
      max_nclist = 1
      DO iac = 1, nkind*nkind
         IF (.NOT. has_intac(iac)) CYCLE
         DO i = 1, SIZE(oce%intac(iac)%alist)
            max_nclist = MAX(max_nclist, oce%intac(iac)%alist(i)%nclist)
         END DO
      END DO

      CPASSERT(ASSOCIATED(sab(1)%nlist_task))
      nl_size = sab(1)%nl_size
      ALLOCATE (ia_ptr(natom + 1), ia_work(natom), slot_order(MAX(1, nl_size)))
      ia_work(:) = 0
      DO i = 1, 2
         DO slot = 1, nl_size
            task => sab(1)%nlist_task(slot)
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
!$OMP          SHARED(natom, nkind, nspins, nimages, cell_to_index, rho_ao, max_nsgf, max_gau, max_nb, &
!$OMP                 max_nclist, w_alloc, len_PC1, ia_ptr, slot_order, sab, paw_kind, nsatbas_kind, &
!$OMP                 has_intac, oce, eps_cpc, rho_atom_set, locks, number_of_locks) &
!$OMP          PRIVATE(nb_pblk, nb_pmax, p_matrix, proj_work1, proj_work2, w_buf, w_offset, w_touched, &
!$OMP                  iatom, ikind, n_nb, inb, TASK, jatom, cell_b, img, pmax, ispin, found, kkind, &
!$OMP                  nsoctot, iac, alist_ac, kac_lo, kac_hi, w_total, kac, clist_a, n_cont_a, katom, &
!$OMP                  jkind, ibc, alist_bc, factor, kbc, clist_b, n_cont_b, list_a, list_b, dista, &
!$OMP                  distb, C_coeff_hh_a, C_coeff_ss_a, C_coeff_hh_b, C_coeff_ss_b, beta, r_coef_h, r_coef_s)

      ALLOCATE (nb_pblk(nspins, max_nb), nb_pmax(max_nb), w_buf(w_alloc, nspins), w_offset(max_nclist))
      ALLOCATE (p_matrix(MAX(1, max_nsgf), MAX(1, max_nsgf)), proj_work1(len_PC1), &
                proj_work2(MAX(1, max_gau), MAX(1, max_gau)), w_touched(max_nclist, nspins))

!$OMP SINGLE
!$    number_of_locks = nspins*natom
!$    ALLOCATE (locks(number_of_locks))
!$OMP END SINGLE

!$OMP DO
!$    DO lock = 1, number_of_locks
!$       call omp_init_lock(locks(lock))
!$    END DO
!$OMP END DO

      ! Factor sum_(A,B) C_A^T P_AB C_B = sum_A C_A^T (sum_B P_AB C_B) across offsite B /= K.
!$OMP DO SCHEDULE(DYNAMIC, 1)
      DO iatom = 1, natom
         n_nb = ia_ptr(iatom + 1) - ia_ptr(iatom)
         IF (n_nb == 0) CYCLE
         task => sab(1)%nlist_task(slot_order(ia_ptr(iatom)))
         ikind = MERGE(task%jkind, task%ikind, task%iatom > task%jatom)

         DO inb = 1, n_nb
            task => sab(1)%nlist_task(slot_order(ia_ptr(iatom) + inb - 1))
            jatom = MAX(task%iatom, task%jatom)
            img = 1
            IF (nimages > 1) img = cell_to_index(task%cell(1), task%cell(2), task%cell(3))
            CPASSERT(img > 0)
            pmax = 0.0_dp
            DO ispin = 1, nspins
               CALL dbcsr_get_block_p(rho_ao(nspins*(img - 1) + ispin)%matrix, iatom, jatom, &
                                      nb_pblk(ispin, inb)%r_coef, found)
               CPASSERT(found)
               pmax = pmax + MAXVAL(ABS(nb_pblk(ispin, inb)%r_coef))
            END DO
            nb_pmax(inb) = pmax
         END DO

         DO kkind = 1, nkind
            IF (.NOT. paw_kind(kkind)) CYCLE
            nsoctot = nsatbas_kind(kkind)
            iac = ikind + nkind*(kkind - 1)
            IF (.NOT. has_intac(iac)) CYCLE
            CALL get_alist(oce%intac(iac), alist_ac, iatom)
            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE

            kac_lo = 1
            DO WHILE (kac_lo <= alist_ac%nclist)
               w_total = 0
               DO kac = kac_lo, alist_ac%nclist
                  n_cont_a = alist_ac%clist(kac)%nsgf_cnt
                  IF (n_cont_a == 0) CYCLE
                  IF (w_total > 0 .AND. w_total + n_cont_a*nsoctot > w_alloc) EXIT
                  w_offset(kac) = w_total + 1
                  w_total = w_total + n_cont_a*nsoctot
               END DO
               kac_hi = kac - 1

               IF (w_total > 0) THEN
                  w_touched(kac_lo:kac_hi, 1:nspins) = .FALSE.
                  DO inb = 1, n_nb
                     pmax = nb_pmax(inb)
                     IF (pmax <= 0.0_dp) CYCLE
                     task => sab(1)%nlist_task(slot_order(ia_ptr(iatom) + inb - 1))
                     jkind = MERGE(task%ikind, task%jkind, task%iatom > task%jatom)
                     jatom = MAX(task%iatom, task%jatom)
                     cell_b = MERGE(-task%cell, task%cell, task%iatom > task%jatom)
                     ibc = jkind + nkind*(kkind - 1)
                     IF (.NOT. has_intac(ibc)) CYCLE
                     CALL get_alist(oce%intac(ibc), alist_bc, jatom)
                     IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                     factor = MERGE(0.5_dp, 1.0_dp, iatom == jatom)

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
                        IF (n_cont_a == 0 .OR. pmax*clist_b%maxac*clist_a%maxac < eps_cpc) CYCLE

                        list_a => clist_a%sgf_list
                        list_b => clist_b%sgf_list
                        dista = (iatom /= katom .OR. ANY(clist_a%cell /= 0))
                        distb = (jatom /= katom .OR. ANY(clist_b%cell /= 0))
                        C_coeff_ss_b => clist_b%acint(:, :, 1)
                        IF (.NOT. distb) THEN
                           C_coeff_ss_a => clist_a%acint(:, :, 1)
                           C_coeff_hh_a => C_coeff_ss_a
                           IF (.NOT. dista) C_coeff_hh_a => clist_a%achint(:, :, 1)
                           C_coeff_hh_b => clist_b%achint(:, :, 1)
                        END IF

                        DO ispin = 1, nspins
                           CALL alist_pre_align_blk(nb_pblk(ispin, inb)%r_coef, &
                                                    SIZE(nb_pblk(ispin, inb)%r_coef, 1), p_matrix, SIZE(p_matrix, 1), &
                                                    list_a, n_cont_a, list_b, n_cont_b)
                           IF (distb) THEN
                              beta = MERGE(1.0_dp, 0.0_dp, w_touched(kac, ispin))
                              w_touched(kac, ispin) = .TRUE.
                              CALL DGEMM('N', 'N', n_cont_a, nsoctot, n_cont_b, factor, p_matrix, SIZE(p_matrix, 1), &
                                         C_coeff_ss_b, n_cont_b, beta, w_buf(w_offset(kac):, ispin), n_cont_a)
                              CYCLE
                           END IF
                           r_coef_h => rho_atom_set(katom)%cpc_h(ispin)%r_coef
                           r_coef_s => rho_atom_set(katom)%cpc_s(ispin)%r_coef
!$                         CALL omp_set_lock(locks((katom - 1)*nspins + ispin))
                           CALL proj_blk(C_coeff_hh_a, C_coeff_ss_a, n_cont_a, C_coeff_hh_b, C_coeff_ss_b, n_cont_b, &
                                         p_matrix, SIZE(p_matrix, 1), r_coef_h, r_coef_s, nsoctot, len_PC1, &
                                         SIZE(proj_work2), factor, .FALSE., proj_work1, proj_work2)
!$                         CALL omp_unset_lock(locks((katom - 1)*nspins + ispin))
                        END DO
                     END DO
                  END DO

                  DO kac = kac_lo, kac_hi
                     clist_a => alist_ac%clist(kac)
                     n_cont_a = clist_a%nsgf_cnt
                     IF (n_cont_a == 0) CYCLE
                     katom = clist_a%catom
                     dista = (iatom /= katom .OR. ANY(clist_a%cell /= 0))
                     C_coeff_ss_a => clist_a%acint(:, :, 1)
                     DO ispin = 1, nspins
                        IF (.NOT. w_touched(kac, ispin)) CYCLE
                        r_coef_h => rho_atom_set(katom)%cpc_h(ispin)%r_coef
                        r_coef_s => rho_atom_set(katom)%cpc_s(ispin)%r_coef
                        CALL DGEMM('T', 'N', nsoctot, nsoctot, n_cont_a, 1.0_dp, C_coeff_ss_a, n_cont_a, &
                                   w_buf(w_offset(kac):, ispin), n_cont_a, 0.0_dp, proj_work2, SIZE(proj_work2, 1))
!$                      CALL omp_set_lock(locks((katom - 1)*nspins + ispin))
                        r_coef_s = r_coef_s + proj_work2(1:nsoctot, 1:nsoctot)
                        IF (dista) r_coef_h = r_coef_h + proj_work2(1:nsoctot, 1:nsoctot)
                        IF (.NOT. dista) CALL DGEMM('T', 'N', nsoctot, nsoctot, n_cont_a, 1.0_dp, &
                                                    clist_a%achint(:, :, 1), n_cont_a, w_buf(w_offset(kac):, ispin), &
                                                    n_cont_a, 1.0_dp, r_coef_h, SIZE(r_coef_h, 1))
!$                      CALL omp_unset_lock(locks((katom - 1)*nspins + ispin))
                     END DO
                  END DO
               END IF
               kac_lo = kac_hi + 1
            END DO
         END DO
      END DO
!$OMP END DO
      ! Wait for all threads to finish the loop before locks can be freed
!$OMP BARRIER

!$OMP DO
!$    DO lock = 1, number_of_locks
!$       call omp_destroy_lock(locks(lock))
!$    END DO
!$OMP END DO
!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

      DEALLOCATE (nb_pblk, nb_pmax, p_matrix, proj_work1, proj_work2, w_buf, w_offset, w_touched)
!$OMP END PARALLEL

      DEALLOCATE (ia_ptr, slot_order)
      CALL reduce_symmetrize_cpc(rho_atom_set, natom, nspins, max_gau, para_env)

      DEALLOCATE (basis_set_list, paw_kind, nsatbas_kind, has_intac)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_atom_coeff

! **************************************************************************************************
!> \brief ...
!> \param rho_atom_set    the type to initialize
!> \param atomic_kind_set list of atomic kinds
!> \param qs_kind_set     the kind set from which to take quantum numbers and basis info
!> \param dft_control     DFT control type
!> \param para_env        parallel environment
!> \par History:
!>      - Generalised by providing the rho_atom_set and the qs_kind_set 12.2019 (A.Bussy)
! **************************************************************************************************
   SUBROUTINE init_rho_atom(rho_atom_set, atomic_kind_set, qs_kind_set, dft_control, para_env)

      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'init_rho_atom'

      CHARACTER(LEN=4)                                   :: my_ikind, my_la, my_llmax
      INTEGER                                            :: handle, ikind, la, ll, llmax, &
                                                            lmax_sphere, max_s_harm, max_s_set, &
                                                            maxl, maxlgto, maxs, na, nat, natom, &
                                                            nr, nspins, quadrature
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: paw_atom
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      TYPE(gapw_control_type), POINTER                   :: gapw_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_set
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL timeset(routineN, handle)

      NULLIFY (basis_1c_set)
      NULLIFY (my_CG, grid_atom, harmonics, atom_list)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(qs_kind_set))

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)

      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="GAPW_1C")

      nspins = dft_control%nspins
      gapw_control => dft_control%qs_control%gapw_control

      lmax_sphere = gapw_control%lmax_sphere

      llmax = MIN(lmax_sphere, 2*maxlgto)
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)

!   *** calculate the Gaunt coefficients up to the maxl ***
      CALL real_gaunt_coefficients(my_CG, maxlgto, maxlgto, llmax)

!   *** initialize the Lebedev grids ***
      CALL init_lebedev_grids()
      quadrature = gapw_control%quadrature

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          paw_atom=paw_atom, &
                          grid_atom=grid_atom, &
                          harmonics=harmonics, &
                          ngrid_rad=nr, ngrid_ang=na)

!     *** determine the Lebedev grid for this kind ***

         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

         IF (llmax > la) THEN
            IF (para_env%is_source()) THEN
               WRITE (my_la, "(I4)") la
               WRITE (my_ikind, "(I4)") ikind
               WRITE (my_llmax, "(I4)") llmax
               CALL cp_warn(__LOCATION__, &
                            "The Lebedev grid has been built for angular momentum l up "// &
                            "to "//TRIM(my_la)//" for atom kind "//TRIM(my_ikind)//", "// &
                            "but the max l of spherical harmonics "//TRIM(my_llmax)//" is "// &
                            "larger and good integration cannot be guaranteed")
            END IF
         END IF

!     *** calculate the radial grid ***
         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)

!     *** calculate the spherical harmonics on the grid ***

         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_1c_set, basis_type="GAPW_1C")
         CALL get_gto_basis_set(gto_basis_set=basis_1c_set, maxl=maxl)
         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na, llmax, maxs, max_s_harm, ll, grid_atom%wa, &
                                    grid_atom%azi, grid_atom%pol)
         CALL get_maxl_CG(harmonics, basis_1c_set, llmax, max_s_harm)

      END DO

      CALL deallocate_lebedev_grids()
      DEALLOCATE (my_CG)

      CALL allocate_rho_atom_internals(rho_atom_set, atomic_kind_set, qs_kind_set, dft_control, para_env)

      CALL timestop(handle)

   END SUBROUTINE init_rho_atom

! **************************************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param atomic_kind_set list of atomic kinds
!> \param qs_kind_set     the kind set from which to take quantum numbers and basis info
!> \param dft_control     DFT control type
!> \param para_env        parallel environment
! **************************************************************************************************
   SUBROUTINE allocate_rho_atom_internals(rho_atom_set, atomic_kind_set, qs_kind_set, dft_control, para_env)

      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho_atom_internals'

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, ispin, &
                                                            max_iso_not0, maxso, mepos, nat, &
                                                            natom, nsatbas, nset, nsotot, nspins, &
                                                            num_pe
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: paw_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(qs_kind_set))

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)

      nspins = dft_control%nspins

      IF (ASSOCIATED(rho_atom_set)) THEN
         CALL deallocate_rho_atom_set(rho_atom_set)
      END IF
      ALLOCATE (rho_atom_set(natom))

      DO ikind = 1, SIZE(atomic_kind_set)

         NULLIFY (atom_list, harmonics)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          paw_atom=paw_atom, &
                          harmonics=harmonics)

         IF (paw_atom) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
            CALL get_gto_basis_set(gto_basis_set=basis_1c, nset=nset, maxso=maxso)
            nsotot = nset*maxso
            CALL get_paw_basis_info(basis_1c, nsatbas=nsatbas)
         END IF

         max_iso_not0 = harmonics%max_iso_not0
         DO iat = 1, nat
            iatom = atom_list(iat)
            !       *** allocate the radial density for each LM,for each atom ***

            ALLOCATE (rho_atom_set(iatom)%rho_rad_h(nspins))
            ALLOCATE (rho_atom_set(iatom)%rho_rad_s(nspins))
            ALLOCATE (rho_atom_set(iatom)%vrho_rad_h(nspins))
            ALLOCATE (rho_atom_set(iatom)%vrho_rad_s(nspins))

            ALLOCATE (rho_atom_set(iatom)%cpc_h(nspins), &
                      rho_atom_set(iatom)%cpc_s(nspins), &
                      rho_atom_set(iatom)%drho_rad_h(nspins), &
                      rho_atom_set(iatom)%drho_rad_s(nspins), &
                      rho_atom_set(iatom)%rho_rad_h_d(3, nspins), &
                      rho_atom_set(iatom)%rho_rad_s_d(3, nspins))
            ALLOCATE (rho_atom_set(iatom)%int_scr_h(nspins), &
                      rho_atom_set(iatom)%int_scr_s(nspins))

            IF (paw_atom) THEN
               DO ispin = 1, nspins
                  ALLOCATE (rho_atom_set(iatom)%cpc_h(ispin)%r_coef(1:nsatbas, 1:nsatbas), &
                            rho_atom_set(iatom)%cpc_s(ispin)%r_coef(1:nsatbas, 1:nsatbas))
                  ALLOCATE (rho_atom_set(iatom)%int_scr_h(ispin)%r_coef(1:nsatbas, 1:nsatbas), &
                            rho_atom_set(iatom)%int_scr_s(ispin)%r_coef(1:nsatbas, 1:nsatbas))

                  rho_atom_set(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
                  rho_atom_set(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
               END DO
            END IF

         END DO ! iat

         num_pe = para_env%num_pe
         mepos = para_env%mepos
         bo = get_limit(nat, num_pe, mepos)
         DO iat = bo(1), bo(2)
            iatom = atom_list(iat)
            ALLOCATE (rho_atom_set(iatom)%ga_Vlocal_gb_h(nspins), &
                      rho_atom_set(iatom)%ga_Vlocal_gb_s(nspins))
            IF (paw_atom) THEN
               DO ispin = 1, nspins
                  CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef, &
                                  1, nsotot, 1, nsotot)
                  CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef, &
                                  1, nsotot, 1, nsotot)

                  rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
                  rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
               END DO
            END IF

         END DO ! iat

      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_rho_atom_internals

! **************************************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param iatom ...
!> \param ispin ...
!> \param nr ...
!> \param max_iso_not0 ...
! **************************************************************************************************
   SUBROUTINE allocate_rho_atom_rad(rho_atom_set, iatom, ispin, nr, max_iso_not0)

      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      INTEGER, INTENT(IN)                                :: iatom, ispin, nr, max_iso_not0

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho_atom_rad'

      INTEGER                                            :: handle, j

      CALL timeset(routineN, handle)

      ALLOCATE (rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef(1:nr, 1:max_iso_not0), &
                rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef(1:nr, 1:max_iso_not0), &
                rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef(1:nr, 1:max_iso_not0), &
                rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef(1:nr, 1:max_iso_not0))

      rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef = 0.0_dp

      ALLOCATE (rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef(nr, max_iso_not0), &
                rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef(nr, max_iso_not0))
      rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp

      DO j = 1, 3
         ALLOCATE (rho_atom_set(iatom)%rho_rad_h_d(j, ispin)%r_coef(nr, max_iso_not0), &
                   rho_atom_set(iatom)%rho_rad_s_d(j, ispin)%r_coef(nr, max_iso_not0))
         rho_atom_set(iatom)%rho_rad_h_d(j, ispin)%r_coef = 0.0_dp
         rho_atom_set(iatom)%rho_rad_s_d(j, ispin)%r_coef = 0.0_dp
      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_rho_atom_rad

! **************************************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param iatom ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE set2zero_rho_atom_rad(rho_atom_set, iatom, ispin)

      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      INTEGER, INTENT(IN)                                :: iatom, ispin

      INTEGER                                            :: j

      rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef = 0.0_dp

      rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef = 0.0_dp

      rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp

      DO j = 1, 3
         rho_atom_set(iatom)%rho_rad_h_d(j, ispin)%r_coef = 0.0_dp
         rho_atom_set(iatom)%rho_rad_s_d(j, ispin)%r_coef = 0.0_dp
      END DO

   END SUBROUTINE set2zero_rho_atom_rad

! **************************************************************************************************
!> \brief Symmetrize PAW cpc_h/cpc_s blocks and MPI-sum packed lower triangles in L2-sized chunks.
!> \param rho_atom_set atomic PAW density set holding cpc_h and cpc_s blocks
!> \param natom total number of atoms
!> \param nspins number of spin components
!> \param max_gau maximum size of the one-center basis
!> \param para_env MPI parallel environment
! **************************************************************************************************
   SUBROUTINE reduce_symmetrize_cpc(rho_atom_set, natom, nspins, max_gau, para_env)
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      INTEGER, INTENT(IN)                                :: natom, nspins, max_gau
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: b_end, b_lo, iat, ispin, jj, m1, n, off, &
                                                            p
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: r_h, r_s

      ALLOCATE (buf(MAX(l2_chunk_dp, nspins*max_gau*(max_gau + 1))))
      b_lo = 1
      DO WHILE (b_lo <= natom)
         b_end = natom
         DO p = 1, 2
            off = 0
            DO iat = b_lo, b_end
               IF (.NOT. ASSOCIATED(rho_atom_set(iat)%cpc_h(1)%r_coef)) CYCLE
               m1 = SIZE(rho_atom_set(iat)%cpc_h(1)%r_coef, 1)
               IF (p == 1 .AND. off > 0 .AND. off + nspins*m1*(m1 + 1) > SIZE(buf)) EXIT
               DO ispin = 1, nspins
                  r_h => rho_atom_set(iat)%cpc_h(ispin)%r_coef
                  r_s => rho_atom_set(iat)%cpc_s(ispin)%r_coef
                  DO jj = 1, m1
                     n = m1 - jj + 1
                     IF (p == 1) buf(off + 1:off + n) = r_h(jj:m1, jj) + r_h(jj, jj:m1)
                     IF (p == 1) buf(off + n + 1:off + 2*n) = r_s(jj:m1, jj) + r_s(jj, jj:m1)
                     IF (p == 2) r_h(jj:m1, jj) = buf(off + 1:off + n)
                     IF (p == 2) r_h(jj, jj:m1) = buf(off + 1:off + n)
                     IF (p == 2) r_s(jj:m1, jj) = buf(off + n + 1:off + 2*n)
                     IF (p == 2) r_s(jj, jj:m1) = buf(off + n + 1:off + 2*n)
                     off = off + 2*n
                  END DO
               END DO
            END DO
            b_end = iat - 1
            IF (p == 1 .AND. off > 0) CALL para_env%sum(buf(1:off))
         END DO
         b_lo = b_end + 1
      END DO
      DEALLOCATE (buf)
   END SUBROUTINE reduce_symmetrize_cpc

END MODULE qs_rho_atom_methods
