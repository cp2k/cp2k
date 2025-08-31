!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A collection of functions used by CNEO-DFT
!>      (see J. Chem. Theory Comput. 2025, 21, 16, 7865–7877)
!> \par History
!>      08.2025 created [zc62]
!> \author Zehua Chen
! **************************************************************************************************
MODULE qs_cneo_methods
   USE ai_verfc,                        ONLY: verfc
   USE ao_util,                         ONLY: trace_r_AxB
   USE atom_operators,                  ONLY: atom_int_release,&
                                              atom_int_setup
   USE atom_types,                      ONLY: CGTO_BASIS,&
                                              atom_basis_type,&
                                              atom_integrals,&
                                              lmat,&
                                              release_atom_basis
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Chen2025,&
                                              cite_reference
   USE core_ae,                         ONLY: verfc_force
   USE cp_control_types,                ONLY: gapw_control_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fourpi,&
                                              pi
   USE mathlib,                         ONLY: get_pseudo_inverse_svd
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso,&
                                              indso_inv,&
                                              init_orbital_pointers,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: massunit
   USE qs_cneo_types,                   ONLY: allocate_rhoz_cneo_set,&
                                              cneo_potential_type,&
                                              get_cneo_potential,&
                                              rhoz_cneo_type,&
                                              set_cneo_potential
   USE qs_cneo_utils,                   ONLY: atom_solve_cneo,&
                                              cneo_gather,&
                                              create_harmonics_atom_cneo,&
                                              create_my_CG_cneo,&
                                              get_maxl_CG_cneo
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: allocate_harmonics_atom,&
                                              get_none0_cg_list,&
                                              harmonics_atom_type
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
   USE util,                            ONLY: get_limit
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE whittaker,                       ONLY: whittaker_c0a,&
                                              whittaker_ci

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cneo_methods'

   PUBLIC :: allocate_rhoz_cneo_internals, calculate_rhoz_cneo, cneo_core_matrices, &
             init_cneo_potential_internals, Vh_1c_nuc_integrals

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param nuc_basis ...
!> \param nuc_soft_basis ...
!> \param gapw_control ...
!> \param grid_atom ...
! **************************************************************************************************
   SUBROUTINE init_cneo_potential_internals(potential, nuc_basis, nuc_soft_basis, gapw_control, grid_atom)

      TYPE(cneo_potential_type), POINTER                 :: potential
      TYPE(gto_basis_set_type), POINTER                  :: nuc_basis, nuc_soft_basis
      TYPE(gapw_control_type), POINTER                   :: gapw_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CHARACTER(len=*), PARAMETER :: routineN = 'init_cneo_potential_internals'

      INTEGER :: handle, i, icg, ico, ii, ipgf, ipgf1, ipgf2, ir, is1, is2, iset, iset1, iset2, &
         iso, iso1, iso2, iso_pgf, iso_set, j, k, k1, k2, l, l_iso, l_sub, l_sum, ll, llmax, &
         lmax12, lmax_expansion, lmax_sphere, lmin12, m, m1, m2, max_iso_not0, max_iso_not0_local, &
         max_s, max_s_harm, maxl, maxso, n1, n2, nl, nne, npgf2, npgf_sum, npsgf, nr, ns, nset, &
         nsgf, nsotot, nsox
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(0:lmat, 100)                    :: set_index, shell_index
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, n2oindex, npgf, npgf_s, &
                                                            nshell, o2nindex
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, ls
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: done_vgg
      REAL(KIND=dp)                                      :: c1, c2, gcc_tmp, mass, massinv, &
                                                            root_zet12, scal, scal1, zet12
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: erf_zet12, g1, g2, gg0, int1, int2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dist
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kin, my_gcc_h, my_gcc_s, oorad2l, ovlp, &
                                                            rad2l, utrans, zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: distance, gcc_h, gcc_s, gg, my_CG
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: vgg
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(atom_integrals), POINTER                      :: integrals
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CPASSERT(ASSOCIATED(potential))
      CPASSERT(ASSOCIATED(nuc_basis))
      CPASSERT(ASSOCIATED(nuc_soft_basis))

      CALL cite_reference(Chen2025)

      CALL timeset(routineN, handle)

      CALL get_cneo_potential(potential, my_gcc_h=my_gcc_h, my_gcc_s=my_gcc_s, &
                              ovlp=ovlp, kin=kin, utrans=utrans, distance=distance, &
                              harmonics=harmonics, gg=gg, vgg=vgg, n2oindex=n2oindex, &
                              o2nindex=o2nindex, rad2l=rad2l, oorad2l=oorad2l)
      CPASSERT(.NOT. ASSOCIATED(my_gcc_h))
      CPASSERT(.NOT. ASSOCIATED(my_gcc_s))
      CPASSERT(.NOT. ASSOCIATED(ovlp))
      CPASSERT(.NOT. ASSOCIATED(kin))
      CPASSERT(.NOT. ASSOCIATED(utrans))
      CPASSERT(.NOT. ASSOCIATED(distance))
      CPASSERT(.NOT. ASSOCIATED(harmonics))
      CPASSERT(.NOT. ASSOCIATED(gg))
      CPASSERT(.NOT. ASSOCIATED(vgg))
      CPASSERT(.NOT. ASSOCIATED(n2oindex))
      CPASSERT(.NOT. ASSOCIATED(o2nindex))
      CPASSERT(.NOT. ASSOCIATED(rad2l))
      CPASSERT(.NOT. ASSOCIATED(oorad2l))

      ! ovlp, kin and utrans parts are mostly copied from atom_kind_orbitals::calculate_atomic_orbitals
      ! and atom_set_basis::set_kind_basis_atomic
      NULLIFY (basis, integrals, grid)
      ALLOCATE (basis, integrals)
      CALL allocate_grid_atom(grid)
      basis%grid => grid
      NULLIFY (basis%am, basis%cm, basis%as, basis%ns, basis%bf, basis%dbf, basis%ddbf)
      ! fill in the basis data structures
      basis%basis_type = CGTO_BASIS
      basis%eps_eig = 1.e-12_dp

      NULLIFY (nshell, npgf, lmin, lmax, ls, zet, gcc_h, first_sgf)
      CALL get_gto_basis_set(nuc_basis, nset=nset, nshell=nshell, npgf=npgf, lmin=lmin, &
                             lmax=lmax, l=ls, nsgf=nsgf, zet=zet, gcc=gcc_h, first_sgf=first_sgf, &
                             maxl=maxl, maxso=maxso, npgf_sum=npgf_sum)
      NULLIFY (npgf_s, gcc_s)
      CALL get_gto_basis_set(nuc_soft_basis, npgf=npgf_s, gcc=gcc_s)
      ! There is such a limitation because we rely on atomic code to build S, T and U.
      ! Usually l=5 is more than enough, suppoting PB6H basis.
      IF (maxl > lmat) &
         CALL cp_abort(__LOCATION__, "Nuclear basis with angular momentum higher than "// &
                       "atom_types::lmat is not supported yet.")

      set_index = 0
      shell_index = 0
      basis%nprim = 0
      basis%nbas = 0
      DO i = 1, nset
         DO j = lmin(i), MIN(lmax(i), lmat)
            basis%nprim(j) = basis%nprim(j) + npgf(i)
         END DO
         DO j = 1, nshell(i)
            l = ls(j, i)
            IF (l <= lmat) THEN
               basis%nbas(l) = basis%nbas(l) + 1
               k = basis%nbas(l)
               CPASSERT(k <= 100)
               set_index(l, k) = i
               shell_index(l, k) = j
            END IF
         END DO
      END DO

      nl = MAXVAL(basis%nprim)
      ns = MAXVAL(basis%nbas)
      ALLOCATE (basis%am(nl, 0:lmat))
      basis%am = 0._dp
      ALLOCATE (basis%cm(nl, ns, 0:lmat))
      basis%cm = 0._dp
      DO l = 0, lmat
         nl = 0
         ns = 0
         DO i = 1, nset
            IF (l >= lmin(i) .AND. l <= lmax(i)) THEN
               DO ipgf = 1, npgf(i)
                  basis%am(nl + ipgf, l) = zet(ipgf, i)
               END DO
               DO ii = 1, nshell(i)
                  IF (ls(ii, i) == l) THEN
                     ns = ns + 1
                     DO ipgf = 1, npgf(i)
                        basis%cm(nl + ipgf, ns, l) = gcc_h(ipgf, ii, i) ! NOTE: not normalized
                     END DO
                  END IF
               END DO
               nl = nl + npgf(i)
            END IF
         END DO
      END DO

      ! overlap, kinetic and transformation matrices
      CALL atom_int_setup(integrals, basis)

      ! make the integrals full matrix form
      ALLOCATE (ovlp(nsgf, nsgf), kin(nsgf, nsgf), utrans(nsgf, nsgf))
      ovlp = 0.0_dp
      kin = 0.0_dp
      utrans = 0.0_dp
      CALL get_cneo_potential(potential, mass=mass)
      mass = mass*massunit
      massinv = 1._dp/mass
      nne = 0 ! number of linear-independent spherical basis functions
      DO l = 0, lmat
         ll = 2*l
         DO k2 = 1, integrals%nne(l)
            DO m = 0, ll
               nne = nne + 1
               DO k1 = 1, basis%nbas(l)
                  scal1 = SQRT(integrals%ovlp(k1, k1, l))
                  i = first_sgf(shell_index(l, k1), set_index(l, k1))
                  utrans(i + m, nne) = integrals%utrans(k1, k2, l)*scal1
               END DO
            END DO
         END DO
         DO k1 = 1, basis%nbas(l)
            scal1 = 1._dp/SQRT(integrals%ovlp(k1, k1, l))
            i = first_sgf(shell_index(l, k1), set_index(l, k1))
            DO k2 = 1, basis%nbas(l)
               scal = scal1/SQRT(integrals%ovlp(k2, k2, l))
               j = first_sgf(shell_index(l, k2), set_index(l, k2))
               DO m = 0, ll
                  ! normalize the integrals
                  ovlp(i + m, j + m) = integrals%ovlp(k1, k2, l)*scal
                  kin(i + m, j + m) = integrals%kin(k1, k2, l)*scal*massinv
               END DO
            END DO
         END DO
      END DO

      nsotot = maxso*nset
      ALLOCATE (my_gcc_h(nsotot, nsgf), my_gcc_s(nsotot, nsgf))
      my_gcc_h = 0.0_dp
      my_gcc_s = 0.0_dp
      ! create gcc that really 3D-normalize the basis functions
      DO l = 0, MIN(maxl, lmat)
         ns = 0
         m = 0
         ll = 2*l
         k = nsoset(l - 1) + 1
         DO i = 1, nset
            IF (l >= lmin(i) .AND. l <= lmax(i)) THEN
               nsox = nsoset(lmax(i))
               DO ii = 1, nshell(i)
                  IF (ls(ii, i) == l) THEN
                     ns = ns + 1
                     k1 = first_sgf(shell_index(l, ns), set_index(l, ns))
                     scal = 1._dp/SQRT(integrals%ovlp(ns, ns, l))
                     DO ipgf = 1, npgf(i)
                        gcc_tmp = gcc_h(ipgf, ii, i)*scal
                        k2 = (ipgf - 1)*nsox + m
                        DO j = 0, ll
                           my_gcc_h(k + k2 + j, k1 + j) = gcc_tmp
                        END DO
                     END DO
                     DO ipgf = 1, npgf_s(i)
                        gcc_tmp = gcc_s(ipgf, ii, i)*scal
                        k2 = (ipgf - 1)*nsox + m
                        DO j = 0, ll
                           my_gcc_s(k + k2 + j, k1 + j) = gcc_tmp
                        END DO
                     END DO
                  END IF
               END DO
            END IF
            m = m + maxso
         END DO
      END DO

      CALL atom_int_release(integrals)
      CALL set_cneo_potential(potential, nsgf=nsgf, nne=nne, nsotot=nsotot, &
                              my_gcc_h=my_gcc_h, my_gcc_s=my_gcc_s, &
                              ovlp=ovlp, kin=kin, utrans=utrans)
      CALL release_atom_basis(basis)
      DEALLOCATE (basis, integrals)

      ! initialize my_CG
      lmax_sphere = gapw_control%lmax_sphere
      ! make sure llmax is at least 1 such that distance matrices can be generated
      llmax = MAX(1, MIN(lmax_sphere, 2*maxl))
      max_s_harm = nsoset(llmax)
      max_s = nsoset(maxl)
      NULLIFY (my_CG)
      CALL reallocate(my_CG, 1, max_s, 1, max_s, 1, max_s_harm)
      CALL create_my_CG_cneo(my_CG, MAX(llmax, 2*maxl, 1), maxl, llmax)

      ! initialize harmonics
      CALL allocate_harmonics_atom(harmonics)
      CALL create_harmonics_atom_cneo(harmonics, my_CG, llmax, max_s, max_s_harm)
      DEALLOCATE (my_CG)
      CALL get_maxl_CG_cneo(harmonics, nuc_basis, llmax, max_s_harm)

      CALL set_cneo_potential(potential, harmonics=harmonics)

      ! initialize my own rad2l and oorad2l
      ! copied from qs_grid_atom::create_grid_atom
      nr = grid_atom%nr
      NULLIFY (rad2l, oorad2l)
      CALL reallocate(rad2l, 1, nr, 0, llmax + 1)
      CALL reallocate(oorad2l, 1, nr, 0, llmax + 1)
      rad2l(:, 0) = 1._dp
      oorad2l(:, 0) = 1._dp
      DO l = 1, llmax + 1
         rad2l(:, l) = rad2l(:, l - 1)*grid_atom%rad(:)
         oorad2l(:, l) = oorad2l(:, l - 1)/grid_atom%rad(:)
      END DO
      CALL set_cneo_potential(potential, rad2l=rad2l, oorad2l=oorad2l)
      ! still need to bump lmax in grid_atom as qs_rho0_types::calculate_g0 uses it
      IF (SIZE(rad2l, 2) > SIZE(grid_atom%rad2l, 2)) THEN
         CPASSERT(SIZE(rad2l, 1) == SIZE(grid_atom%rad2l, 1))
         DEALLOCATE (grid_atom%rad2l)
         NULLIFY (grid_atom%rad2l)
         CALL reallocate(grid_atom%rad2l, 1, nr, 0, llmax + 1)
         grid_atom%rad2l = rad2l
      END IF
      IF (SIZE(oorad2l, 2) > SIZE(grid_atom%oorad2l, 2)) THEN
         CPASSERT(SIZE(oorad2l, 1) == SIZE(grid_atom%oorad2l, 1))
         DEALLOCATE (grid_atom%oorad2l)
         NULLIFY (grid_atom%oorad2l)
         CALL reallocate(grid_atom%oorad2l, 1, nr, 0, llmax + 1)
         grid_atom%oorad2l = oorad2l
      END IF

      ! distance matrices
      ALLOCATE (distance(nsgf, nsgf, 3), dist(nsotot, nsotot, 3))
      distance = 0.0_dp
      dist = 0.0_dp
      ! initialize gg and vgg
      ! mostly copied from qs_rho_atom_methods::calculate_rho_atom
      max_iso_not0 = harmonics%max_iso_not0
      lmax_expansion = indso(1, max_iso_not0)
      my_CG => harmonics%my_CG
      ALLOCATE (g1(nr), g2(nr), gg0(nr), gg(nr, 0:2*maxl, npgf_sum*(npgf_sum + 1)/2))
      ALLOCATE (erf_zet12(nr), vgg(nr, 0:2*maxl, 0:indso(1, max_iso_not0), npgf_sum*(npgf_sum + 1)/2))
      ALLOCATE (done_vgg(0:2*maxl, 0:indso(1, max_iso_not0)))
      ALLOCATE (int1(nr), int2(nr))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

      j = 0
      m1 = 0
      DO iset1 = 1, nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1, iset1
            n2 = nsoset(lmax(iset2))

            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
            CPASSERT(max_iso_not0_local .LE. max_iso_not0)

            DO ipgf1 = 1, npgf(iset1)
               g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))

               IF (iset2 == iset1) THEN
                  npgf2 = ipgf1
               ELSE
                  npgf2 = npgf(iset2)
               END IF
               DO ipgf2 = 1, npgf2
                  zet12 = zet(ipgf1, iset1) + zet(ipgf2, iset2)

                  ! distance part
                  ! -1 -> y -> 2, 0 -> z -> 3, 1 -> x -> 1
                  DO m = -1, 1
                     k = m + 3
                     IF (m == 1) k = 1
                     iso = indso_inv(1, m)
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2

                        l = indso(1, is1) + indso(1, is2)
                        dist(iso1, iso2, k) = dist(iso1, iso2, k) + my_CG(is1, is2, iso)* &
                                              pi*dfac(l + 2)/ &
                                              ((2.0_dp*zet12)**((l + 3)/2)*SQRT(3.0_dp*zet12))
                        dist(iso2, iso1, k) = dist(iso1, iso2, k) ! symmetric
                     END DO !icg
                  END DO

                  ! gg and vgg part
                  j = j + 1
                  g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                  lmin12 = lmin(iset1) + lmin(iset2)
                  lmax12 = lmax(iset1) + lmax(iset2)

                  root_zet12 = SQRT(zet12)
                  DO ir = 1, nr
                     erf_zet12(ir) = erf(root_zet12*grid_atom%rad(ir))
                  END DO

                  gg(:, :, j) = 0.0_dp
                  vgg(:, :, :, j) = 0.0_dp
                  done_vgg = .FALSE.
                  ! reduce the number of terms in the expansion local densities
                  IF (lmin12 .LE. lmax_expansion) THEN
                     IF (lmin12 == 0) THEN
                        gg(1:nr, lmin12, j) = g1(1:nr)*g2(1:nr)
                        gg0(1:nr) = gg(1:nr, lmin12, j)
                     ELSE
                        gg0(1:nr) = g1(1:nr)*g2(1:nr)
                        gg(1:nr, lmin12, j) = rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                     END IF

                     ! reduce the number of terms in the expansion local densities
                     IF (lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

                     DO l = lmin12 + 1, lmax12
                        gg(1:nr, l, j) = grid_atom%rad(1:nr)*gg(1:nr, l - 1, j)
                     END DO

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
                              vgg(1:nr, l, l_iso, j) = erf_zet12(1:nr)*oorad2l(1:nr, 1)*c2
                           ELSE
                              CALL whittaker_c0a(int1, grid_atom%rad, gg0, erf_zet12, zet12, l, l_iso, nr)
                              CALL whittaker_ci(int2, grid_atom%rad, gg0, zet12, L_sub, nr)

                              DO ir = 1, nr
                                 int2(ir) = rad2l(ir, l_iso)*int2(ir)
                                 vgg(ir, l, l_iso, j) = c1*(int1(ir) + int2(ir))
                              END DO
                           END IF
                           done_vgg(l, l_iso) = .TRUE.
                        END DO
                     END DO
                  END IF ! lmax_expansion

               END DO ! ipgf2
            END DO ! ipgf1
            m2 = m2 + maxso
         END DO ! iset2
         m1 = m1 + maxso
      END DO ! iset1

      DEALLOCATE (g1, g2, gg0, erf_zet12, int1, int2, done_vgg)
      DEALLOCATE (cg_list, cg_n_list)

      ALLOCATE (work(nsotot, nsgf))
      DO k = 1, 3
         CALL dgemm("N", "N", nsotot, nsgf, nsotot, 1.0_dp, dist(:, :, k), nsotot, my_gcc_h, &
                    nsotot, 0.0_dp, work, nsotot)
         CALL dgemm("T", "N", nsgf, nsgf, nsotot, 1.0_dp, my_gcc_h, nsotot, work, &
                    nsotot, 0.0_dp, distance(:, :, k), nsgf)
      END DO
      DEALLOCATE (work, dist)
      CALL set_cneo_potential(potential, distance=distance, gg=gg, vgg=vgg)

      ! Index transformation OLD-NEW
      ! copied from paw_proj_set_types::build_projector
      ALLOCATE (o2nindex(nsotot))
      ALLOCATE (n2oindex(nsotot))
      o2nindex = 0
      n2oindex = 0
      ico = 1
      DO iset = 1, nset
         iso_set = (iset - 1)*maxso + 1
         nsox = nsoset(lmax(iset))
         DO ipgf = 1, npgf(iset)
            iso_pgf = iso_set + (ipgf - 1)*nsox
            iso = iso_pgf + nsoset(lmin(iset) - 1)
            DO l = lmin(iset), lmax(iset)
               DO k = 1, nso(l)
                  n2oindex(ico) = iso
                  o2nindex(iso) = ico
                  iso = iso + 1
                  ico = ico + 1
               END DO
            END DO
         END DO
      END DO
      npsgf = ico - 1
      CALL set_cneo_potential(potential, npsgf=npsgf, n2oindex=n2oindex, o2nindex=o2nindex)

      CALL timestop(handle)

   END SUBROUTINE init_cneo_potential_internals

! **************************************************************************************************
!> \brief ...
!> \param rhoz_cneo_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_rhoz_cneo_internals(rhoz_cneo_set, atomic_kind_set, &
                                           qs_kind_set, qs_env)

      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rhoz_cneo_internals'

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, &
                                                            max_iso_not0, mepos, nat, natom, &
                                                            npsgf, nr, nsgf, nsotot, num_pe
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: paw_atom
      TYPE(cneo_potential_type), POINTER                 :: cneo_potential
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)

      CALL allocate_rhoz_cneo_set(rhoz_cneo_set, natom)

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, para_env=para_env)

      DO ikind = 1, SIZE(atomic_kind_set)

         NULLIFY (cneo_potential)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          ngrid_rad=nr, &
                          paw_atom=paw_atom, &
                          cneo_potential=cneo_potential)

         IF (ASSOCIATED(cneo_potential)) THEN
            CPASSERT(paw_atom)

            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)

            nsgf = cneo_potential%nsgf
            npsgf = cneo_potential%npsgf
            nsotot = cneo_potential%nsotot

            DO iat = 1, nat
               iatom = atom_list(iat)

               ! density matrices, core and soft vmat will be broadcast to all processes
               ALLOCATE (rhoz_cneo_set(iatom)%pmat(1:nsgf, 1:nsgf), &
                         rhoz_cneo_set(iatom)%cpc_h(1:npsgf, 1:npsgf), &
                         rhoz_cneo_set(iatom)%cpc_s(1:npsgf, 1:npsgf), &
                         rhoz_cneo_set(iatom)%core(1:nsgf, 1:nsgf), &
                         rhoz_cneo_set(iatom)%vmat(1:nsgf, 1:nsgf))
               rhoz_cneo_set(iatom)%pmat = 0.0_dp
               rhoz_cneo_set(iatom)%cpc_h = 0.0_dp
               rhoz_cneo_set(iatom)%cpc_s = 0.0_dp
               rhoz_cneo_set(iatom)%core = 0.0_dp
               rhoz_cneo_set(iatom)%vmat = 0.0_dp
            END DO

            max_iso_not0 = cneo_potential%harmonics%max_iso_not0
            num_pe = para_env%num_pe
            mepos = para_env%mepos
            bo = get_limit(nat, num_pe, mepos)
            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)

               ALLOCATE (rhoz_cneo_set(iatom)%fmat(1:nsgf, 1:nsgf), &
                         rhoz_cneo_set(iatom)%wfn(1:nsgf, 1:nsgf))
               rhoz_cneo_set(iatom)%fmat = 0.0_dp
               rhoz_cneo_set(iatom)%wfn = 0.0_dp

               ALLOCATE (rhoz_cneo_set(iatom)%rho_rad_h(1:nr, 1:max_iso_not0), &
                         rhoz_cneo_set(iatom)%rho_rad_s(1:nr, 1:max_iso_not0), &
                         rhoz_cneo_set(iatom)%vrho_rad_h(1:nr, 1:max_iso_not0), &
                         rhoz_cneo_set(iatom)%vrho_rad_s(1:nr, 1:max_iso_not0))
               rhoz_cneo_set(iatom)%rho_rad_h = 0.0_dp
               rhoz_cneo_set(iatom)%rho_rad_s = 0.0_dp
               rhoz_cneo_set(iatom)%vrho_rad_h = 0.0_dp
               rhoz_cneo_set(iatom)%vrho_rad_s = 0.0_dp

               NULLIFY (rhoz_cneo_set(iatom)%ga_Vlocal_gb_h)
               CALL reallocate(rhoz_cneo_set(iatom)%ga_Vlocal_gb_h, 1, nsotot, 1, nsotot)
               rhoz_cneo_set(iatom)%ga_Vlocal_gb_h = 0.0_dp
               NULLIFY (rhoz_cneo_set(iatom)%ga_Vlocal_gb_s)
               CALL reallocate(rhoz_cneo_set(iatom)%ga_Vlocal_gb_s, 1, nsotot, 1, nsotot)
               rhoz_cneo_set(iatom)%ga_Vlocal_gb_s = 0.0_dp
            END DO ! iat
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_rhoz_cneo_internals

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param nder ...
! **************************************************************************************************
   SUBROUTINE cneo_core_matrices(qs_env, calculate_forces, nder)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: nder

      LOGICAL                                            :: use_virial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_cneo
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (rhoz_cneo_set)
      CALL get_qs_env(qs_env=qs_env, rhoz_cneo_set=rhoz_cneo_set)

      IF (ASSOCIATED(rhoz_cneo_set)) THEN
         NULLIFY (force, virial)
         ! force
         IF (calculate_forces) CALL get_qs_env(qs_env=qs_env, force=force)
         ! virial
         CALL get_qs_env(qs_env=qs_env, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         NULLIFY (qs_kind_set, atomic_kind_set, particle_set, distribution_1d, para_env, sab_cneo)
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                         particle_set=particle_set, local_particles=distribution_1d, &
                         para_env=para_env, sab_cneo=sab_cneo)
         CALL build_core_cneo(rhoz_cneo_set, force, virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, distribution_1d, &
                              sab_cneo, para_env)
      END IF

   END SUBROUTINE cneo_core_matrices

! **************************************************************************************************
!> \brief ...
!> \param rhoz_cneo_set ...
!> \param force ...
!> \param virial ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param nder ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param distribution_1d ...
!> \param sab_cneo ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE build_core_cneo(rhoz_cneo_set, force, virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, distribution_1d, &
                              sab_cneo, para_env)
      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: calculate_forces, use_virial
      INTEGER, INTENT(IN)                                :: nder
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_cneo
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_core_cneo'

      INTEGER :: atom_a, handle, iat, iatom, ikind, iset, jatom, jkind, jset, ldai, ldsab, maxco, &
         maxl, maxnset, maxsgf, mepos, na_plus, nat, natom, nb_plus, ncoa, ncob, nij, nkind, nset, &
         nthread, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf, nsgf
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(KIND=dp)                                      :: alpha_c, core_charge, core_radius, dab, &
                                                            f0, rab2, zeta_i, zeta_j
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: habd, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: hab, pab, verf, vnuc
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, force_i, rab
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_thread
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ap_iterator
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(cneo_potential_type), POINTER                 :: cneo_potential, cneo_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: core, pmat, rpgf, sphi, zet
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius
      REAL(KIND=dp), DIMENSION(3, SIZE(particle_set))    :: force_thread

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      END IF

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      force_thread = 0.0_dp
      pv_thread = 0.0_dp

      ! re-initialize core matrices to zero, as later will use para_env%sum to broadcast
      DO ikind = 1, nkind
         NULLIFY (cneo_potential)
         CALL get_qs_kind(qs_kind_set(ikind), cneo_potential=cneo_potential)

         IF (ASSOCIATED(cneo_potential)) THEN
            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
            DO iat = 1, nat
               iatom = atom_list(iat)
               rhoz_cneo_set(iatom)%core = 0.0_dp
            END DO
         END IF
      END DO

      CALL get_qs_kind_set(qs_kind_set, basis_type="NUC", &
                           maxco=maxco, maxlgto=maxl, maxsgf=maxsgf, maxnset=maxnset)
      CALL init_orbital_pointers(maxl + nder + 1)
      ldsab = MAX(maxco, maxsgf)
      ldai = ncoset(maxl + nder + 1)

      nthread = 1
!$    nthread = omp_get_max_threads()

      CALL neighbor_list_iterator_create(ap_iterator, sab_cneo, search=.TRUE., nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (rhoz_cneo_set, ap_iterator, distribution_1d, calculate_forces, use_virial, &
!$OMP          qs_kind_set, nthread, ncoset, nkind, iat, ldsab, maxnset, ldai, nder, maxl, &
!$OMP          maxco, para_env) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, basis_set, first_sgf, lmax, lmin, npgf, nset, &
!$OMP          nsgf, rpgf, sphi, zet, set_radius, zeta_i, zeta_j, alpha_c, core_charge, &
!$OMP          core_radius, rab, rab2, dab, core, pmat, iset, ncoa, sgfa, jset, ncob, sgfb, &
!$OMP          work, pab, hab, na_plus, nb_plus, verf, vnuc, force_a, force_b, force_i, &
!$OMP          mepos, habd, f0, nij, ff, cneo_potential, cneo_tmp) &
!$OMP REDUCTION (+ : pv_thread, force_thread )

      mepos = 0
!$    mepos = omp_get_thread_num()

      ALLOCATE (hab(ldsab, ldsab, maxnset*(maxnset + 1)/2), work(ldsab, ldsab))
      ALLOCATE (verf(ldai, ldai, 2*maxl + nder + 1), vnuc(ldai, ldai, 2*maxl + nder + 1), ff(0:2*maxl + nder))
      IF (calculate_forces) THEN
         ALLOCATE (pab(maxco, maxco, maxnset*(maxnset + 1)/2))
      END IF

      DO ikind = 1, nkind
         NULLIFY (cneo_potential)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, basis_type="NUC", &
                          cneo_potential=cneo_potential, zeff=zeta_i)
         IF (ASSOCIATED(cneo_potential)) THEN
            CPASSERT(ASSOCIATED(basis_set))
            first_sgf => basis_set%first_sgf
            lmax => basis_set%lmax
            lmin => basis_set%lmin
            npgf => basis_set%npgf
            nset = basis_set%nset
            nsgf => basis_set%nsgf_set
            rpgf => basis_set%pgf_radius
            set_radius => basis_set%set_radius
            sphi => basis_set%sphi
            zet => basis_set%zet

!$OMP DO SCHEDULE(GUIDED)
            DO iat = 1, distribution_1d%n_el(ikind)
               iatom = distribution_1d%list(ikind)%array(iat)
               core => rhoz_cneo_set(iatom)%core
               CPASSERT(ASSOCIATED(core))
               core = cneo_potential%kin ! copy kinetic matrix to core
               IF (calculate_forces) THEN
                  CPASSERT(rhoz_cneo_set(iatom)%ready)
                  pmat => rhoz_cneo_set(iatom)%pmat
                  CPASSERT(ASSOCIATED(pmat))
                  ! *** Decontract density matrix ***
                  DO iset = 1, nset
                     ncoa = npgf(iset)*ncoset(lmax(iset))
                     sgfa = first_sgf(1, iset)
                     DO jset = 1, iset
                        ncob = npgf(jset)*ncoset(lmax(jset))
                        sgfb = first_sgf(1, jset)
                        nij = jset + (iset - 1)*iset/2
                        work(1:ncoa, 1:nsgf(jset)) = MATMUL(sphi(1:ncoa, sgfa:sgfa + nsgf(iset) - 1), &
                                                            pmat(sgfa:sgfa + nsgf(iset) - 1, sgfb:sgfb + nsgf(jset) - 1))
                        pab(1:ncoa, 1:ncob, nij) = MATMUL(work(1:ncoa, 1:nsgf(jset)), &
                                                          TRANSPOSE(sphi(1:ncob, sgfb:sgfb + nsgf(jset) - 1)))
                     END DO
                  END DO
               END IF

               hab = 0._dp
               DO jkind = 1, nkind
                  NULLIFY (cneo_tmp)
                  CALL get_qs_kind(qs_kind_set(jkind), cneo_potential=cneo_tmp)
                  IF (.NOT. ASSOCIATED(cneo_tmp)) THEN
                     CALL get_qs_kind(qs_kind_set(jkind), &
                                      alpha_core_charge=alpha_c, zeff=zeta_j, &
                                      ccore_charge=core_charge, core_charge_radius=core_radius)
                     CALL nl_set_sub_iterator(ap_iterator, ikind, jkind, iatom, mepos=mepos)

                     DO WHILE (nl_sub_iterate(ap_iterator, mepos=mepos) == 0)
                        CALL get_iterator_info(ap_iterator, jatom=jatom, r=rab, mepos=mepos)
                        rab2 = SUM(rab*rab)
                        dab = SQRT(rab2)
                        IF (MAXVAL(set_radius(:)) + core_radius < dab) CYCLE
                        DO iset = 1, nset
                           IF (set_radius(iset) + core_radius < dab) CYCLE
                           ncoa = npgf(iset)*ncoset(lmax(iset))
                           sgfa = first_sgf(1, iset)
                           DO jset = 1, iset ! symmetric
                              IF (set_radius(jset) + core_radius < dab) CYCLE
                              ncob = npgf(jset)*ncoset(lmax(jset))
                              sgfb = first_sgf(1, jset)
                              nij = jset + (iset - 1)*iset/2
                              IF (calculate_forces) THEN
                                 IF (jset == iset) THEN
                                    f0 = -zeta_i
                                 ELSE
                                    f0 = -2.0_dp*zeta_i
                                 END IF
                                 na_plus = npgf(iset)*ncoset(lmax(iset) + nder)
                                 nb_plus = npgf(jset)*ncoset(lmax(jset))
                                 ALLOCATE (habd(na_plus, nb_plus))
                                 habd = 0._dp
                                 CALL verfc( &
                                    lmax(iset) + nder, npgf(iset), zet(:, iset), rpgf(:, iset), lmin(iset), &
                                    lmax(jset), npgf(jset), zet(:, jset), rpgf(:, jset), lmin(jset), &
                                    alpha_c, core_radius, zeta_j, core_charge, &
                                    (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, rab, rab2, rab2, &
                                    hab(:, :, nij), verf, vnuc, ff(0:), nder, habd)

                                 ! *** The derivatives w.r.t. atomic center b are    ***
                                 ! *** calculated using the translational invariance ***
                                 ! *** of the first derivatives                      ***
                                 CALL verfc_force(habd, pab(:, :, nij), force_a, force_b, nder, &
                                                  lmax(iset), lmin(iset), npgf(iset), zet(:, iset), &
                                                  lmax(jset), lmin(jset), npgf(jset), zet(:, jset), &
                                                  (/0.0_dp, 0.0_dp, 0.0_dp/))

                                 DEALLOCATE (habd)
                                 force_i = force_a + force_b

                                 force_thread(1, iatom) = force_thread(1, iatom) + f0*force_i(1)
                                 force_thread(2, iatom) = force_thread(2, iatom) + f0*force_i(2)
                                 force_thread(3, iatom) = force_thread(3, iatom) + f0*force_i(3)

                                 force_thread(1, jatom) = force_thread(1, jatom) - f0*force_i(1)
                                 force_thread(2, jatom) = force_thread(2, jatom) - f0*force_i(2)
                                 force_thread(3, jatom) = force_thread(3, jatom) - f0*force_i(3)

                                 IF (use_virial) THEN
                                    CALL virial_pair_force(pv_thread, f0, force_i, rab)
                                 END IF
                              ELSE
                                 CALL verfc( &
                                    lmax(iset), npgf(iset), zet(:, iset), rpgf(:, iset), lmin(iset), &
                                    lmax(jset), npgf(jset), zet(:, jset), rpgf(:, jset), lmin(jset), &
                                    alpha_c, core_radius, zeta_j, core_charge, &
                                    (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, rab, rab2, rab2, &
                                    hab(:, :, nij), verf, vnuc, ff(0:))
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF
               END DO
               ! *** Contract nuclear repulsion integrals
               DO iset = 1, nset
                  ncoa = npgf(iset)*ncoset(lmax(iset))
                  sgfa = first_sgf(1, iset)
                  DO jset = 1, iset
                     ncob = npgf(jset)*ncoset(lmax(jset))
                     sgfb = first_sgf(1, jset)
                     nij = jset + (iset - 1)*iset/2
                     work(1:ncoa, 1:nsgf(jset)) = MATMUL(hab(1:ncoa, 1:ncob, nij), &
                                                         sphi(1:ncob, sgfb:sgfb + nsgf(jset) - 1))
                     core(sgfa:sgfa + nsgf(iset) - 1, sgfb:sgfb + nsgf(jset) - 1) = &
                        core(sgfa:sgfa + nsgf(iset) - 1, sgfb:sgfb + nsgf(jset) - 1) - zeta_i* &
                        MATMUL(TRANSPOSE(sphi(1:ncoa, sgfa:sgfa + nsgf(iset) - 1)), work(1:ncoa, 1:nsgf(jset)))
                     ! symmetrize core matrix
                     IF (iset /= jset) THEN
                        core(sgfb:sgfb + nsgf(jset) - 1, sgfa:sgfa + nsgf(iset) - 1) = &
                           TRANSPOSE(core(sgfa:sgfa + nsgf(iset) - 1, sgfb:sgfb + nsgf(jset) - 1))
                     END IF
                  END DO
               END DO
            END DO
         END IF
      END DO

      DEALLOCATE (hab, work, verf, vnuc, ff)
      IF (calculate_forces) THEN
         DEALLOCATE (pab)
      END IF

!$OMP END PARALLEL

      CALL neighbor_list_iterator_release(ap_iterator)

      IF (calculate_forces) THEN
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, &
                                  kind_of=kind_of)
!$OMP DO
         DO iatom = 1, natom
            atom_a = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%cneo_potential(:, atom_a) = force(ikind)%cneo_potential(:, atom_a) + &
                                                     force_thread(:, iatom)
         END DO
!$OMP END DO
      END IF

      IF (calculate_forces .AND. use_virial) THEN
         virial%pv_ppl = virial%pv_ppl + pv_thread
         virial%pv_virial = virial%pv_virial + pv_thread
      END IF

      ! broadcast core matrices
      DO ikind = 1, nkind
         NULLIFY (cneo_potential)
         CALL get_qs_kind(qs_kind_set(ikind), cneo_potential=cneo_potential)

         IF (ASSOCIATED(cneo_potential)) THEN
            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
            DO iat = 1, nat
               iatom = atom_list(iat)
               CALL para_env%sum(rhoz_cneo_set(iatom)%core)
            END DO
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_core_cneo

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param potential ...
!> \param cg_list ...
!> \param cg_n_list ...
!> \param nset ...
!> \param npgf ...
!> \param lmin ...
!> \param lmax ...
!> \param maxl ...
!> \param maxso ...
! **************************************************************************************************
   SUBROUTINE calculate_rhoz_cneo(rho, potential, cg_list, cg_n_list, nset, npgf, &
                                  lmin, lmax, maxl, maxso)

      TYPE(rhoz_cneo_type), POINTER                      :: rho
      TYPE(cneo_potential_type), POINTER                 :: potential
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: cg_list
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: cg_n_list
      INTEGER, INTENT(IN)                                :: nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, lmin, lmax
      INTEGER, INTENT(IN)                                :: maxl, maxso

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rhoz_cneo'

      INTEGER :: handle, i, i1, i2, icg, ipgf1, ipgf2, iset1, iset2, iso, iso1, iso1_first, &
         iso1_last, iso2, iso2_first, iso2_last, iter, j, l, l1, l2, l_iso, lmax_expansion, m1s, &
         m2s, max_iso_not0, max_iso_not0_local, max_iter, max_s_harm, n1s, n2s, nne, npgf2, npsgf, &
         nsgf, nsotot, size1, size2
      INTEGER, DIMENSION(:), POINTER                     :: n2oindex, o2nindex
      REAL(KIND=dp)                                      :: det, df_norm, factor, g0, g0p, g1, step, &
                                                            zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ener
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: CPCH_sphere, CPCS_sphere, work
      REAL(KIND=dp), DIMENSION(3)                        :: df, f_tmp, r, r_tmp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: jac, jac_inv
      REAL(KIND=dp), DIMENSION(:), POINTER               :: f
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: core, cpc_h, cpc_s, fmat, int_local_h, &
         int_local_s, my_gcc_h, my_gcc_s, pmat, rho_rad_h, rho_rad_s, utrans, vmat, vrho_rad_h, &
         vrho_rad_s, wfn
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: distance, gg, my_CG
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: vgg
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CPASSERT(ASSOCIATED(rho))
      CPASSERT(ASSOCIATED(potential))

      CALL timeset(routineN, handle)

      ! convert ga_Vlocal_gb to compressed form V_Hartree
      ! use fmat to store V_Hartree
      NULLIFY (utrans, my_gcc_h, my_gcc_s, distance, n2oindex, o2nindex)
      CALL get_cneo_potential(potential, zeff=zeff, nsgf=nsgf, nne=nne, npsgf=npsgf, &
                              nsotot=nsotot, my_gcc_h=my_gcc_h, my_gcc_s=my_gcc_s, &
                              utrans=utrans, distance=distance, n2oindex=n2oindex, &
                              o2nindex=o2nindex)
      fmat => rho%fmat
      int_local_h => rho%ga_Vlocal_gb_h
      int_local_s => rho%ga_Vlocal_gb_s
      ALLOCATE (work(nsotot, nsgf))
      CALL dgemm("N", "N", nsotot, nsgf, nsotot, 1.0_dp, int_local_h, nsotot, my_gcc_h, &
                 nsotot, 0.0_dp, work, nsotot)
      CALL dgemm("T", "N", nsgf, nsgf, nsotot, 1.0_dp, my_gcc_h, nsotot, work, &
                 nsotot, 0.0_dp, fmat, nsgf)
      CALL dgemm("N", "N", nsotot, nsgf, nsotot, 1.0_dp, int_local_s, nsotot, my_gcc_s, &
                 nsotot, 0.0_dp, work, nsotot)
      CALL dgemm("T", "N", nsgf, nsgf, nsotot, -1.0_dp, my_gcc_s, nsotot, work, &
                 nsotot, 1.0_dp, fmat, nsgf)
      ! add the soft basis FFT grid part
      vmat => rho%vmat
      fmat = fmat + vmat

      core => rho%core
      wfn => rho%wfn
      pmat => rho%pmat
      f => rho%f
      ALLOCATE (ener(nne))
      ! build the fock matrix: F = T + V_core + V_Hartree
      fmat = fmat + core
      ! conduct the constrained optimization with F + f*x
      ! initial guess of f is taken from the result of last iteration
      CALL atom_solve_cneo(fmat, f, utrans, wfn, ener, pmat, r, distance, nsgf, nne)
      ! test if zero initial guess is better
      IF (SQRT(DOT_PRODUCT(r, r)) > 1.e-12_dp .AND. DOT_PRODUCT(f, f) /= 0.0_dp) THEN
         CALL atom_solve_cneo(fmat, (/0.0_dp, 0.0_dp, 0.0_dp/), utrans, wfn, &
                              ener, pmat, r_tmp, distance, nsgf, nne)
         IF (DOT_PRODUCT(r_tmp, r_tmp) < DOT_PRODUCT(r, r)) THEN
            f = 0.0_dp
            r = r_tmp
         END IF
      END IF
      max_iter = 20
      iter = 0
      ! using Newton's method to solve for f
      DO WHILE (SQRT(DOT_PRODUCT(r, r)) > 1.e-12_dp)
         iter = iter + 1
         ! construct numerical Jacobian with one-side finite difference
         DO i = 1, 3
            f_tmp = f
            f_tmp(i) = f(i) + SIGN(1.e-4_dp, r(i)) ! forward or backward based on the sign of r
            CALL atom_solve_cneo(fmat, f_tmp, utrans, wfn, ener, pmat, r_tmp, distance, nsgf, nne)
            DO j = 1, 3
               jac(j, i) = (r_tmp(j) - r(j))*SIGN(1.e4_dp, r(i))
            END DO
         END DO
         CALL invert_matrix_3x3(jac, jac_inv, det)
         IF (ABS(det) < 1.0E-8_dp) THEN
            CALL cp_warn(__LOCATION__, "Determinant of the CNEO position Jacobian is small! "// &
                         TRIM(cp_to_string(det))//" Trying central difference.")
            ! construct numerical Jacobian with central finite difference
            DO i = 1, 3
               f_tmp = f
               f_tmp(i) = f(i) - SIGN(1.e-4_dp, r(i))
               CALL atom_solve_cneo(fmat, f_tmp, utrans, wfn, ener, pmat, r_tmp, distance, nsgf, nne)
               DO j = 1, 3
                  jac(j, i) = (jac(j, i)*SIGN(1.e-4_dp, r(i)) + r(j) - r_tmp(j)) &
                              /SIGN(2.e-4_dp, r(i))
               END DO
            END DO
            CALL invert_matrix_3x3(jac, jac_inv, det)
            IF (ABS(det) < 1.0E-8_dp) THEN
               CALL cp_warn(__LOCATION__, "Determinant of the CNEO position Jacobian is small! "// &
                            "(Central difference) "//TRIM(cp_to_string(det))//" Using pseudoinverse.")
            END IF
            CALL invert_matrix_3x3(jac, jac_inv, det, try_svd=.TRUE.)
         END IF
         df = -RESHAPE(MATMUL(jac_inv, RESHAPE(r, (/3, 1/))), (/3/))
         df_norm = SQRT(DOT_PRODUCT(df, df))
         f_tmp = f
         r_tmp = r
         g0 = SQRT(DOT_PRODUCT(r_tmp, r_tmp))
         f = f_tmp + df
         CALL atom_solve_cneo(fmat, f, utrans, wfn, ener, pmat, r, distance, nsgf, nne)
         g1 = SQRT(DOT_PRODUCT(r, r))
         step = 1.0_dp
         DO WHILE (g1 >= g0)
            ! line search
            IF (step < 0.0101_dp) THEN
               CPWARN("CNEO nuclear position constraint solver line search failure.")
               EXIT
            END IF
            g0p = -g0/(step*df_norm)
            step = step*MAX(-g0p/(2.0_dp*(g1 - g0 - g0p)), 0.1_dp)
            f = f_tmp + step*df
            CALL atom_solve_cneo(fmat, f, utrans, wfn, ener, pmat, r, distance, nsgf, nne)
            g1 = SQRT(DOT_PRODUCT(r, r))
         END DO
         IF (iter >= max_iter) THEN
            CALL cp_warn(__LOCATION__, "CNEO nuclear position constraint solver failed to "// &
                         "converge in "//TRIM(cp_to_string(max_iter))//" steps. "// &
                         "Nuclear position error (x,y,z): "//TRIM(cp_to_string(r(1)))// &
                         ", "//TRIM(cp_to_string(r(2)))//", "//TRIM(cp_to_string(r(3)))// &
                         ". This does not hurt as long as it is not the final SCF iteration.")
            EXIT
         END IF
      END DO
      DEALLOCATE (ener)
      rho%e_core = trace_r_AxB(core, nsgf, pmat, nsgf, nsgf, nsgf)

      ! decontract the density matrix
      ! first use ga_Vlocal_gb to store the decompressed form
      CALL dgemm("N", "N", nsotot, nsgf, nsgf, 1.0_dp, my_gcc_h, nsotot, pmat, nsgf, &
                 0.0_dp, work, nsotot)
      CALL dgemm("N", "T", nsotot, nsotot, nsgf, 1.0_dp, work, nsotot, my_gcc_h, nsotot, &
                 0.0_dp, int_local_h, nsotot)
      CALL dgemm("N", "N", nsotot, nsgf, nsgf, 1.0_dp, my_gcc_s, nsotot, pmat, nsgf, &
                 0.0_dp, work, nsotot)
      CALL dgemm("N", "T", nsotot, nsotot, nsgf, 1.0_dp, work, nsotot, my_gcc_s, nsotot, &
                 0.0_dp, int_local_s, nsotot)
      DEALLOCATE (work)
      ! compress the density matrix
      cpc_h => rho%cpc_h
      cpc_s => rho%cpc_s
      CALL cneo_gather(int_local_h, cpc_h, npsgf, n2oindex)
      CALL cneo_gather(int_local_s, cpc_s, npsgf, n2oindex)
      ! restore ga_Vlocal_gb to zeros
      int_local_h = 0.0_dp
      int_local_s = 0.0_dp

      ! construct the nuclear density and its Hartree potential
      ! rho_rad_h and vrho_rad_h should contain the -Zeff factor
      ! mostly copied from qs_rho_atom_methods::calculate_rho_atom
      NULLIFY (harmonics, gg, vgg)
      CALL get_cneo_potential(potential, harmonics=harmonics, gg=gg, vgg=vgg)
      rho_rad_h => rho%rho_rad_h
      rho_rad_s => rho%rho_rad_s
      rho_rad_h = 0.0_dp
      rho_rad_s = 0.0_dp
      vrho_rad_h => rho%vrho_rad_h
      vrho_rad_s => rho%vrho_rad_s
      vrho_rad_h = 0.0_dp
      vrho_rad_s = 0.0_dp
      my_CG => harmonics%my_CG
      max_iso_not0 = harmonics%max_iso_not0
      max_s_harm = harmonics%max_s_harm
      lmax_expansion = indso(1, max_iso_not0)

      ALLOCATE (CPCH_sphere(nsoset(maxl), nsoset(maxl)))
      ALLOCATE (CPCS_sphere(nsoset(maxl), nsoset(maxl)))
      j = 0
      m1s = 0
      DO iset1 = 1, nset
         m2s = 0
         n1s = nsoset(lmax(iset1))
         DO iset2 = 1, iset1

            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
            CPASSERT(max_iso_not0_local .LE. max_iso_not0)

            n2s = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               iso1_first = nsoset(lmin(iset1) - 1) + 1 + n1s*(ipgf1 - 1) + m1s
               iso1_last = nsoset(lmax(iset1)) + n1s*(ipgf1 - 1) + m1s
               size1 = iso1_last - iso1_first + 1
               iso1_first = o2nindex(iso1_first)
               iso1_last = o2nindex(iso1_last)
               i1 = iso1_last - iso1_first + 1
               CPASSERT(size1 == i1)
               i1 = nsoset(lmin(iset1) - 1) + 1

               IF (iset2 == iset1) THEN
                  npgf2 = ipgf1
               ELSE
                  npgf2 = npgf(iset2)
               END IF
               DO ipgf2 = 1, npgf2
                  j = j + 1
                  iso2_first = nsoset(lmin(iset2) - 1) + 1 + n2s*(ipgf2 - 1) + m2s
                  iso2_last = nsoset(lmax(iset2)) + n2s*(ipgf2 - 1) + m2s
                  size2 = iso2_last - iso2_first + 1
                  iso2_first = o2nindex(iso2_first)
                  iso2_last = o2nindex(iso2_last)
                  i2 = iso2_last - iso2_first + 1
                  CPASSERT(size2 == i2)
                  i2 = nsoset(lmin(iset2) - 1) + 1

                  IF (iset2 == iset1 .AND. ipgf2 == ipgf1) THEN
                     factor = -zeff
                  ELSE
                     factor = -2.0_dp*zeff
                  END IF

                  CPCH_sphere = 0.0_dp
                  CPCS_sphere = 0.0_dp
                  CPCH_sphere(i1:i1 + size1 - 1, i2:i2 + size2 - 1) = cpc_h(iso1_first:iso1_last, iso2_first:iso2_last)
                  CPCS_sphere(i1:i1 + size1 - 1, i2:i2 + size2 - 1) = cpc_s(iso1_first:iso1_last, iso2_first:iso2_last)
                  DO iso = 1, max_iso_not0_local
                     l_iso = indso(1, iso)
                     DO icg = 1, cg_n_list(iso)
                        iso1 = cg_list(1, icg, iso)
                        iso2 = cg_list(2, icg, iso)

                        l1 = indso(1, iso1)
                        l2 = indso(1, iso2)

                        l = indso(1, iso1) + indso(1, iso2)
                        CPASSERT(l <= lmax_expansion)

                        rho_rad_h(:, iso) = rho_rad_h(:, iso) + gg(:, l, j)* &
                                            CPCH_sphere(iso1, iso2)*my_CG(iso1, iso2, iso)*factor

                        rho_rad_s(:, iso) = rho_rad_s(:, iso) + gg(:, l, j)* &
                                            CPCS_sphere(iso1, iso2)*my_CG(iso1, iso2, iso)*factor

                        vrho_rad_h(:, iso) = vrho_rad_h(:, iso) + vgg(:, l, l_iso, j)* &
                                             CPCH_sphere(iso1, iso2)*my_CG(iso1, iso2, iso)*factor

                        vrho_rad_s(:, iso) = vrho_rad_s(:, iso) + vgg(:, l, l_iso, j)* &
                                             CPCS_sphere(iso1, iso2)*my_CG(iso1, iso2, iso)*factor
                     END DO ! icg
                  END DO ! iso
               END DO ! ipgf2
            END DO ! ipgf1
            m2s = m2s + maxso
         END DO ! iset2
         m1s = m1s + maxso
      END DO ! iset1
      DEALLOCATE (CPCH_sphere, CPCS_sphere)

      CALL timestop(handle)

   END SUBROUTINE calculate_rhoz_cneo

! **************************************************************************************************
!> \brief Mostly copied from hartree_local_methods::Vh_1c_atom_integrals
!> \param rhoz_cneo ...
!> \param zeff ...
!> \param aVh1b_hh ...
!> \param aVh1b_ss ...
!> \param aVh1b_00 ...
!> \param Vh1_h ...
!> \param Vh1_s ...
!> \param max_iso_not0_elec ...
!> \param max_iso_not0_nuc ...
!> \param max_s_harm ...
!> \param llmax ...
!> \param cg_list ...
!> \param cg_n_list ...
!> \param nset ...
!> \param npgf ...
!> \param lmin ...
!> \param lmax ...
!> \param nsotot ...
!> \param maxso ...
!> \param nchan_0 ...
!> \param gsph ...
!> \param g0_h_w ...
!> \param my_CG ...
!> \param Qlm_gg ...
! **************************************************************************************************
   SUBROUTINE Vh_1c_nuc_integrals(rhoz_cneo, zeff, &
                                  aVh1b_hh, aVh1b_ss, aVh1b_00, Vh1_h, Vh1_s, &
                                  max_iso_not0_elec, max_iso_not0_nuc, &
                                  max_s_harm, llmax, cg_list, cg_n_list, &
                                  nset, npgf, lmin, lmax, nsotot, maxso, nchan_0, gsph, &
                                  g0_h_w, my_CG, Qlm_gg)

      TYPE(rhoz_cneo_type), POINTER                      :: rhoz_cneo
      REAL(KIND=dp), INTENT(IN)                          :: zeff
      REAL(KIND=dp), DIMENSION(:, :)                     :: aVh1b_hh, aVh1b_ss, aVh1b_00
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: Vh1_h, Vh1_s
      INTEGER, INTENT(IN)                                :: max_iso_not0_elec, max_iso_not0_nuc, &
                                                            max_s_harm, llmax
      INTEGER, DIMENSION(:, :, :)                        :: cg_list
      INTEGER, DIMENSION(:)                              :: cg_n_list
      INTEGER, INTENT(IN)                                :: nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, lmin, lmax
      INTEGER, INTENT(IN)                                :: nsotot, maxso, nchan_0
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gsph
      REAL(KIND=dp), DIMENSION(:, 0:)                    :: g0_h_w
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: my_CG, Qlm_gg

      INTEGER                                            :: icg, ipgf1, ipgf2, ir, is1, is2, iset1, &
                                                            iset2, iso, iso1, iso2, l_ang, m1, m2, &
                                                            max_iso_not0_local, n1, n2, nr
      REAL(KIND=dp)                                      :: gVg_0, gVg_h, gVg_s

      !       Calculate the integrals of the potential with 2 primitives
      aVh1b_hh = 0.0_dp
      aVh1b_ss = 0.0_dp
      aVh1b_00 = 0.0_dp

      nr = SIZE(gsph, 1)

      m1 = 0
      DO iset1 = 1, nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1, nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, llmax, cg_list, cg_n_list, max_iso_not0_local)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               DO ipgf2 = 1, npgf(iset2)
                  DO iso = 1, MIN(max_iso_not0_elec, max_iso_not0_nuc)
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2
                        gVg_h = 0.0_dp
                        gVg_s = 0.0_dp

                        DO ir = 1, nr
                           gVg_h = gVg_h + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_h(ir, iso)
                           gVg_s = gVg_s + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_s(ir, iso)
                        END DO ! ir

                        aVh1b_hh(iso1, iso2) = aVh1b_hh(iso1, iso2) + gVg_h*my_CG(is1, is2, iso)
                        aVh1b_ss(iso1, iso2) = aVh1b_ss(iso1, iso2) + gVg_s*my_CG(is1, is2, iso)

                     END DO !icg
                  END DO ! iso
                  DO iso = max_iso_not0_elec + 1, max_iso_not0_nuc
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2
                        gVg_s = 0.0_dp

                        DO ir = 1, nr
                           gVg_s = gVg_s + gsph(ir, iso1)*gsph(ir, iso2)*Vh1_s(ir, iso)
                        END DO ! ir

                        aVh1b_ss(iso1, iso2) = aVh1b_ss(iso1, iso2) + gVg_s*my_CG(is1, is2, iso)

                     END DO !icg
                  END DO ! iso
                  DO iso = 1, MIN(nchan_0, max_iso_not0_nuc)
                     l_ang = indso(1, iso)
                     gVg_0 = SUM(Vh1_s(:, iso)*g0_h_w(:, l_ang))
                     DO icg = 1, cg_n_list(iso)
                        is1 = cg_list(1, icg, iso)
                        is2 = cg_list(2, icg, iso)

                        iso1 = is1 + n1*(ipgf1 - 1) + m1
                        iso2 = is2 + n2*(ipgf2 - 1) + m2

                        aVh1b_00(iso1, iso2) = aVh1b_00(iso1, iso2) + gVg_0*Qlm_gg(iso1, iso2, iso)

                     END DO !icg
                  END DO ! iso
               END DO ! ipgf2
            END DO ! ipgf1
            m2 = m2 + maxso
         END DO ! iset2
         m1 = m1 + maxso
      END DO !iset1

      CALL daxpy(nsotot*nsotot, -zeff, aVh1b_hh, 1, rhoz_cneo%ga_Vlocal_gb_h, 1)
      CALL daxpy(nsotot*nsotot, -zeff, aVh1b_ss, 1, rhoz_cneo%ga_Vlocal_gb_s, 1)
      CALL daxpy(nsotot*nsotot, zeff, aVh1b_00, 1, rhoz_cneo%ga_Vlocal_gb_h, 1)
      CALL daxpy(nsotot*nsotot, zeff, aVh1b_00, 1, rhoz_cneo%ga_Vlocal_gb_s, 1)

   END SUBROUTINE Vh_1c_nuc_integrals

! **************************************************************************************************
!> \brief Analytical inversion of a 3x3 matrix
!> \param matrix ...
!> \param inv_matrix ...
!> \param det ...
!> \param try_svd ...
! **************************************************************************************************
   SUBROUTINE invert_matrix_3x3(matrix, inv_matrix, det, try_svd)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: matrix
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: inv_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: det
      LOGICAL, INTENT(IN), OPTIONAL                      :: try_svd

      LOGICAL                                            :: my_try_svd

      my_try_svd = .FALSE.
      IF (PRESENT(try_svd)) my_try_svd = try_svd

      det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) &
            - matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) &
            + matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
      IF (ABS(det) < 1.0E-8_dp) THEN
         IF (my_try_svd) THEN
            ! pseudo inverse using SVD
            CALL get_pseudo_inverse_svd(matrix, inv_matrix, 1.0E-6_dp, det)
         ELSE
            inv_matrix = 0.0_dp
         END IF
      ELSE
         inv_matrix(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
         inv_matrix(1, 2) = matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3)
         inv_matrix(1, 3) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
         inv_matrix(2, 1) = matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3)
         inv_matrix(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
         inv_matrix(2, 3) = matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3)
         inv_matrix(3, 1) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
         inv_matrix(3, 2) = matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2)
         inv_matrix(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
         inv_matrix = inv_matrix/det
      END IF
   END SUBROUTINE invert_matrix_3x3

END MODULE qs_cneo_methods
