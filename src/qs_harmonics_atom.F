!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE qs_harmonics_atom

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: lebedev_grid
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE spherical_harmonics,             ONLY: dy_lm,&
                                              y_lm
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harmonics_atom'

   TYPE harmonics_atom_type
      INTEGER                                :: max_s_harm, llmax, &
                                                max_iso_not0, &
                                                dmax_iso_not0, &
                                                damax_iso_not0, &
                                                ngrid
      REAL(dp), DIMENSION(:, :), POINTER   :: a, slm
      REAL(dp), DIMENSION(:, :, :), POINTER   :: dslm, dslm_dxyz
      REAL(dp), DIMENSION(:, :, :), POINTER   :: my_CG
      REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz
      REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz_asym
      REAL(dp), DIMENSION(:), POINTER        :: slm_int

   END TYPE harmonics_atom_type

   PUBLIC :: allocate_harmonics_atom, &
             create_harmonics_atom, &
             deallocate_harmonics_atom, &
             get_none0_cg_list

   PUBLIC :: harmonics_atom_type, get_maxl_CG

   INTERFACE get_none0_cg_list
      MODULE PROCEDURE get_none0_cg_list3, get_none0_cg_list4
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief   Allocate a spherical harmonics set for the atom grid.
!> \param harmonics ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_harmonics_atom(harmonics)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      IF (ASSOCIATED(harmonics)) CALL deallocate_harmonics_atom(harmonics)

      ALLOCATE (harmonics)

      harmonics%max_s_harm = 0
      harmonics%llmax = 0
      harmonics%max_iso_not0 = 0
      harmonics%dmax_iso_not0 = 0
      harmonics%damax_iso_not0 = 0
      harmonics%ngrid = 0

      NULLIFY (harmonics%slm)
      NULLIFY (harmonics%dslm)
      NULLIFY (harmonics%dslm_dxyz)
      NULLIFY (harmonics%slm_int)
      NULLIFY (harmonics%my_CG)
      NULLIFY (harmonics%my_CG_dxyz)
      NULLIFY (harmonics%my_CG_dxyz_asym)
      NULLIFY (harmonics%a)

   END SUBROUTINE allocate_harmonics_atom

! **************************************************************************************************
!> \brief   Deallocate the spherical harmonics set for the atom grid.
!> \param harmonics ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_harmonics_atom(harmonics)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      IF (ASSOCIATED(harmonics)) THEN

         IF (ASSOCIATED(harmonics%slm)) THEN
            DEALLOCATE (harmonics%slm)
         END IF

         IF (ASSOCIATED(harmonics%dslm)) THEN
            DEALLOCATE (harmonics%dslm)
         END IF

         IF (ASSOCIATED(harmonics%dslm_dxyz)) THEN
            DEALLOCATE (harmonics%dslm_dxyz)
         END IF

         IF (ASSOCIATED(harmonics%slm_int)) THEN
            DEALLOCATE (harmonics%slm_int)
         END IF

         IF (ASSOCIATED(harmonics%my_CG)) THEN
            DEALLOCATE (harmonics%my_CG)
         END IF

         IF (ASSOCIATED(harmonics%my_CG_dxyz)) THEN
            DEALLOCATE (harmonics%my_CG_dxyz)
         END IF

         IF (ASSOCIATED(harmonics%my_CG_dxyz_asym)) THEN
            DEALLOCATE (harmonics%my_CG_dxyz_asym)
         END IF

         IF (ASSOCIATED(harmonics%a)) THEN
            DEALLOCATE (harmonics%a)
         END IF

         DEALLOCATE (harmonics)
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The pointer harmonics is not associated and "// &
                       "cannot be deallocated")
      END IF

   END SUBROUTINE deallocate_harmonics_atom

! **************************************************************************************************
!> \brief ...
!> \param harmonics ...
!> \param my_CG ...
!> \param na ...
!> \param llmax ...
!> \param maxs ...
!> \param max_s_harm ...
!> \param ll ...
!> \param wa ...
!> \param azi ...
!> \param pol ...
!> \note Slight refactoring + OMP parallelized (03.2020 A. Bussy)
! **************************************************************************************************
   SUBROUTINE create_harmonics_atom(harmonics, my_CG, na, llmax, maxs, max_s_harm, ll, wa, azi, pol)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      INTEGER, INTENT(IN)                                :: na, llmax, maxs, max_s_harm, ll
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: wa, azi, pol

      CHARACTER(len=*), PARAMETER :: routineN = 'create_harmonics_atom'

      INTEGER                                            :: handle, i, ia, ic, is, is1, is2, iso, &
                                                            iso1, iso2, l, l1, l2, lx, ly, lz, m, &
                                                            m1, m2, n
      REAL(dp)                                           :: drx, dry, drz, rx, ry, rz
      REAL(dp), DIMENSION(2)                             :: cin, dylm
      REAL(dp), DIMENSION(:), POINTER                    :: slm_int, y
      REAL(dp), DIMENSION(:, :), POINTER                 :: dc, slm
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dslm_dxyz

      CALL timeset(routineN, handle)

      NULLIFY (y, slm, dslm_dxyz, dc)

      CPASSERT(ASSOCIATED(harmonics))

      harmonics%max_s_harm = max_s_harm
      harmonics%llmax = llmax
      harmonics%ngrid = na

      NULLIFY (harmonics%my_CG, harmonics%my_CG_dxyz, harmonics%my_CG_dxyz_asym)
      CALL reallocate(harmonics%my_CG, 1, maxs, 1, maxs, 1, max_s_harm)
      CALL reallocate(harmonics%my_CG_dxyz, 1, 3, 1, maxs, 1, maxs, 1, max_s_harm)
      CALL reallocate(harmonics%my_CG_dxyz_asym, 1, 3, 1, maxs, 1, maxs, 1, max_s_harm)

      DO i = 1, max_s_harm
         DO is1 = 1, maxs
            harmonics%my_CG(1:maxs, is1, i) = my_CG(1:maxs, is1, i)
         END DO
      END DO

      ! allocate and calculate the spherical harmonics LM for this grid
      ! and their derivatives
      NULLIFY (harmonics%slm, harmonics%dslm, harmonics%dslm_dxyz, harmonics%a, harmonics%slm_int)
      CALL reallocate(harmonics%slm, 1, na, 1, max_s_harm)
      CALL reallocate(harmonics%dslm, 1, 2, 1, na, 1, maxs)
      CALL reallocate(harmonics%dslm_dxyz, 1, 3, 1, na, 1, max_s_harm)
      CALL reallocate(harmonics%a, 1, 3, 1, na)
      CALL reallocate(harmonics%slm_int, 1, max_s_harm)

      NULLIFY (slm, dslm_dxyz, slm_int)
      slm => harmonics%slm
      dslm_dxyz => harmonics%dslm_dxyz
      dslm_dxyz = 0.0_dp
      slm_int => harmonics%slm_int
      slm_int = 0.0_dp

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED (slm,dslm_dxyz,slm_int,max_s_harm,ll,lebedev_grid,na,harmonics,wa,indco,orbtramat) &
!$OMP SHARED (nso,nsoset,nco,maxs,indso,ncoset,pol,azi,llmax) &
!$OMP PRIVATE(ia,iso,l,m,i,lx,ly,lz,rx,ry,rz,drx,dry,drz,ic,dc,iso1,iso2,cin,dylm) &
!$OMP PRIVATE(is1,l1,m1,is2,l2,m2,is,n,y)

      ALLOCATE (y(na))
      ALLOCATE (dc(nco(llmax), 3))

!$OMP DO
      DO iso = 1, max_s_harm
         l = indso(1, iso)
         m = indso(2, iso)
         CALL y_lm(lebedev_grid(ll)%r, y, l, m)

         DO ia = 1, na
            slm(ia, iso) = y(ia)
            slm_int(iso) = slm_int(iso) + slm(ia, iso)*wa(ia)
         END DO ! ia
      END DO ! iso
!$OMP END DO

!$OMP DO
      DO ia = 1, na
         harmonics%a(:, ia) = lebedev_grid(ll)%r(:, ia)
      END DO
!$OMP END DO

      !
      ! The derivatives dslm_dxyz and its expansions my_CG_dxyz and my_CG_dxyz_asymm
      ! are NOT the dSlm/dx but the scaled by r**(l-1) derivatives of the monomial
      ! terms x^n1 y^n2 z^n3 transformed by spherical harmonics expansion coefficients
      !

!$OMP DO
      DO ia = 1, na
         DO l = 0, indso(1, max_s_harm)
            DO ic = 1, nco(l)
               lx = indco(1, ic + ncoset(l - 1))
               ly = indco(2, ic + ncoset(l - 1))
               lz = indco(3, ic + ncoset(l - 1))

               IF (lx == 0) THEN
                  rx = 1.0_dp
                  drx = 0.0_dp
               ELSE IF (lx == 1) THEN
                  rx = lebedev_grid(ll)%r(1, ia)
                  drx = 1.0_dp
               ELSE
                  rx = lebedev_grid(ll)%r(1, ia)**lx
                  drx = REAL(lx, dp)*lebedev_grid(ll)%r(1, ia)**(lx - 1)
               END IF
               IF (ly == 0) THEN
                  ry = 1.0_dp
                  dry = 0.0_dp
               ELSE IF (ly == 1) THEN
                  ry = lebedev_grid(ll)%r(2, ia)
                  dry = 1.0_dp
               ELSE
                  ry = lebedev_grid(ll)%r(2, ia)**ly
                  dry = REAL(ly, dp)*lebedev_grid(ll)%r(2, ia)**(ly - 1)
               END IF
               IF (lz == 0) THEN
                  rz = 1.0_dp
                  drz = 0.0_dp
               ELSE IF (lz == 1) THEN
                  rz = lebedev_grid(ll)%r(3, ia)
                  drz = 1.0_dp
               ELSE
                  rz = lebedev_grid(ll)%r(3, ia)**lz
                  drz = REAL(lz, dp)*lebedev_grid(ll)%r(3, ia)**(lz - 1)
               END IF
               dc(ic, 1) = drx*ry*rz
               dc(ic, 2) = rx*dry*rz
               dc(ic, 3) = rx*ry*drz
            END DO
            n = nsoset(l - 1)
            DO is = 1, nso(l)
               iso = is + n
               DO ic = 1, nco(l)
                  dslm_dxyz(:, ia, iso) = dslm_dxyz(:, ia, iso) + &
                                          orbtramat(l)%slm(is, ic)*dc(ic, :)
               END DO
            END DO
         END DO ! l
      END DO !ia
!$OMP END DO

      ! Expansion coefficients of the cartesian derivatives
      ! of the product of two harmonics :
      ! d(Y(l1m1) * Y(l2m2))/dx ; d(Y(l1m1) * Y(l2m2))/dy ; d(Y(l1m1) * Y(l2m2))/dz

!$OMP DO COLLAPSE(3)
      DO iso1 = 1, maxs
         DO iso2 = 1, maxs

            DO iso = 1, max_s_harm
               rx = 0.0_dp
               ry = 0.0_dp
               rz = 0.0_dp

               DO ia = 1, na
                  rx = rx + wa(ia)*slm(ia, iso)* &
                       (dslm_dxyz(1, ia, iso1)*slm(ia, iso2) + slm(ia, iso1)*dslm_dxyz(1, ia, iso2))
                  ry = ry + wa(ia)*slm(ia, iso)* &
                       (dslm_dxyz(2, ia, iso1)*slm(ia, iso2) + slm(ia, iso1)*dslm_dxyz(2, ia, iso2))
                  rz = rz + wa(ia)*slm(ia, iso)* &
                       (dslm_dxyz(3, ia, iso1)*slm(ia, iso2) + slm(ia, iso1)*dslm_dxyz(3, ia, iso2))
               END DO

               harmonics%my_CG_dxyz(1, iso1, iso2, iso) = rx

               harmonics%my_CG_dxyz(2, iso1, iso2, iso) = ry

               harmonics%my_CG_dxyz(3, iso1, iso2, iso) = rz

            END DO
         END DO
      END DO
!$OMP END DO

      ! Expansion coefficients of the cartesian of the combinations
      ! Y(l1m1) * d(Y(l2m2))/dx -  d(Y(l1m1))/dx * Y(l2m2)
      ! Y(l1m1) * d(Y(l2m2))/dy -  d(Y(l1m1))/dy * Y(l2m2)
      ! Y(l1m1) * d(Y(l2m2))/dz -  d(Y(l1m1))/dz * Y(l2m2)

!$OMP DO COLLAPSE(3)
      DO iso1 = 1, maxs
         DO iso2 = 1, maxs
            DO iso = 1, max_s_harm
               drx = 0.0_dp
               dry = 0.0_dp
               drz = 0.0_dp

               DO ia = 1, na
                  drx = drx + wa(ia)*slm(ia, iso)* &
                        (-dslm_dxyz(1, ia, iso1)*slm(ia, iso2) + &
                         slm(ia, iso1)*dslm_dxyz(1, ia, iso2))
                  dry = dry + wa(ia)*slm(ia, iso)* &
                        (-dslm_dxyz(2, ia, iso1)*slm(ia, iso2) + &
                         slm(ia, iso1)*dslm_dxyz(2, ia, iso2))
                  drz = drz + wa(ia)*slm(ia, iso)* &
                        (-dslm_dxyz(3, ia, iso1)*slm(ia, iso2) + &
                         slm(ia, iso1)*dslm_dxyz(3, ia, iso2))
               END DO

               harmonics%my_CG_dxyz_asym(1, iso1, iso2, iso) = drx

               harmonics%my_CG_dxyz_asym(2, iso1, iso2, iso) = dry

               harmonics%my_CG_dxyz_asym(3, iso1, iso2, iso) = drz

            END DO ! iso
         END DO ! iso2
      END DO ! iso1
!$OMP END DO

      ! Calculate the derivatives of the harmonics with respect of the 2 angles
      ! the first angle (polar) is acos(lebedev_grid(ll)%r(3))
      ! the second angle (azimutal) is atan(lebedev_grid(ll)%r(2)/lebedev_grid(ll)%r(1))
!$OMP DO
      DO iso = 1, maxs
         l = indso(1, iso)
         m = indso(2, iso)
         DO ia = 1, na
            cin(1) = pol(ia)
            cin(2) = azi(ia)
            CALL dy_lm(cin, dylm, l, m)
            harmonics%dslm(:, ia, iso) = dylm(:)
         END DO
      END DO
!$OMP END DO

      ! expansion coefficients of product of polar angle derivatives (dslm(1...)) in
      ! spherical harmonics (used for tau functionals)

      DEALLOCATE (y, dc)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE create_harmonics_atom

! **************************************************************************************************
!> \brief ...
!> \param harmonics ...
!> \param orb_basis ...
!> \param llmax ...
!> \param max_s_harm ...
! **************************************************************************************************
   SUBROUTINE get_maxl_CG(harmonics, orb_basis, llmax, max_s_harm)

      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      INTEGER, INTENT(IN)                                :: llmax, max_s_harm

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_maxl_CG'

      INTEGER                                            :: damax_iso_not0, dmax_iso_not0, handle, &
                                                            is1, is2, itmp, max_iso_not0, nset
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(harmonics))

      CALL get_gto_basis_set(gto_basis_set=orb_basis, lmax=lmax, lmin=lmin, nset=nset)

      !   *** Assign indexes for the non null CG coefficients ***
      max_iso_not0 = 0
      dmax_iso_not0 = 0
      damax_iso_not0 = 0
      DO is1 = 1, nset
         DO is2 = 1, nset
            CALL get_none0_cg_list(harmonics%my_CG, &
                                   lmin(is1), lmax(is1), lmin(is2), lmax(is2), &
                                   max_s_harm, llmax, max_iso_not0=itmp)
            max_iso_not0 = MAX(max_iso_not0, itmp)
            CALL get_none0_cg_list(harmonics%my_CG_dxyz, &
                                   lmin(is1), lmax(is1), lmin(is2), lmax(is2), &
                                   max_s_harm, llmax, max_iso_not0=itmp)
            dmax_iso_not0 = MAX(dmax_iso_not0, itmp)
            CALL get_none0_cg_list(harmonics%my_CG_dxyz_asym, &
                                   lmin(is1), lmax(is1), lmin(is2), lmax(is2), &
                                   max_s_harm, llmax, max_iso_not0=itmp)
            damax_iso_not0 = MAX(damax_iso_not0, itmp)
         END DO ! is2
      END DO ! is1
      harmonics%max_iso_not0 = max_iso_not0
      harmonics%dmax_iso_not0 = dmax_iso_not0
      harmonics%damax_iso_not0 = damax_iso_not0

      CALL timestop(handle)

   END SUBROUTINE get_maxl_CG

! **************************************************************************************************
!> \brief ...
!> \param cgc ...
!> \param lmin1 ...
!> \param lmax1 ...
!> \param lmin2 ...
!> \param lmax2 ...
!> \param max_s_harm ...
!> \param llmax ...
!> \param list ...
!> \param n_list ...
!> \param max_iso_not0 ...
! **************************************************************************************************
   SUBROUTINE get_none0_cg_list4(cgc, lmin1, lmax1, lmin2, lmax2, max_s_harm, llmax, &
                                 list, n_list, max_iso_not0)

      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: cgc
      INTEGER, INTENT(IN)                                :: lmin1, lmax1, lmin2, lmax2, max_s_harm, &
                                                            llmax
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), OPTIONAL :: list
      INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL       :: n_list
      INTEGER, INTENT(OUT)                               :: max_iso_not0

      INTEGER                                            :: iso, iso1, iso2, l1, l2, nlist

      CPASSERT(nsoset(lmax1) .LE. SIZE(cgc, 2))
      CPASSERT(nsoset(lmax2) .LE. SIZE(cgc, 3))
      CPASSERT(max_s_harm .LE. SIZE(cgc, 4))
      IF (PRESENT(n_list) .AND. PRESENT(list)) THEN
         CPASSERT(max_s_harm .LE. SIZE(list, 3))
      END IF
      max_iso_not0 = 0
      IF (PRESENT(n_list) .AND. PRESENT(list)) n_list = 0
      DO iso = 1, max_s_harm
         nlist = 0
         DO l1 = lmin1, lmax1
            DO iso1 = nsoset(l1 - 1) + 1, nsoset(l1)
               DO l2 = lmin2, lmax2
                  IF (l1 + l2 > llmax) CYCLE
                  DO iso2 = nsoset(l2 - 1) + 1, nsoset(l2)
                     IF (ABS(cgc(1, iso1, iso2, iso)) + &
                         ABS(cgc(2, iso1, iso2, iso)) + &
                         ABS(cgc(3, iso1, iso2, iso)) > 1.E-8_dp) THEN
                        nlist = nlist + 1
                        IF (PRESENT(n_list) .AND. PRESENT(list)) THEN
                           list(1, nlist, iso) = iso1
                           list(2, nlist, iso) = iso2
                        END IF
                        max_iso_not0 = MAX(max_iso_not0, iso)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (PRESENT(n_list) .AND. PRESENT(list)) n_list(iso) = nlist
      END DO
   END SUBROUTINE get_none0_cg_list4

! **************************************************************************************************
!> \brief ...
!> \param cgc ...
!> \param lmin1 ...
!> \param lmax1 ...
!> \param lmin2 ...
!> \param lmax2 ...
!> \param max_s_harm ...
!> \param llmax ...
!> \param list ...
!> \param n_list ...
!> \param max_iso_not0 ...
! **************************************************************************************************
   SUBROUTINE get_none0_cg_list3(cgc, lmin1, lmax1, lmin2, lmax2, max_s_harm, llmax, &
                                 list, n_list, max_iso_not0)

      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: cgc
      INTEGER, INTENT(IN)                                :: lmin1, lmax1, lmin2, lmax2, max_s_harm, &
                                                            llmax
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), OPTIONAL :: list
      INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL       :: n_list
      INTEGER, INTENT(OUT)                               :: max_iso_not0

      INTEGER                                            :: iso, iso1, iso2, l1, l2, nlist

      CPASSERT(nsoset(lmax1) .LE. SIZE(cgc, 1))
      CPASSERT(nsoset(lmax2) .LE. SIZE(cgc, 2))
      CPASSERT(max_s_harm .LE. SIZE(cgc, 3))
      IF (PRESENT(n_list) .AND. PRESENT(list)) THEN
         CPASSERT(max_s_harm .LE. SIZE(list, 3))
      END IF
      max_iso_not0 = 0
      IF (PRESENT(n_list) .AND. PRESENT(list)) n_list = 0
      DO iso = 1, max_s_harm
         nlist = 0
         DO l1 = lmin1, lmax1
            DO iso1 = nsoset(l1 - 1) + 1, nsoset(l1)
               DO l2 = lmin2, lmax2
                  IF (l1 + l2 > llmax) CYCLE
                  DO iso2 = nsoset(l2 - 1) + 1, nsoset(l2)
                     IF (ABS(cgc(iso1, iso2, iso)) > 1.E-8_dp) THEN
                        nlist = nlist + 1
                        IF (PRESENT(n_list) .AND. PRESENT(list)) THEN
                           list(1, nlist, iso) = iso1
                           list(2, nlist, iso) = iso2
                        END IF
                        max_iso_not0 = MAX(max_iso_not0, iso)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (PRESENT(n_list) .AND. PRESENT(list)) n_list(iso) = nlist
      END DO
   END SUBROUTINE get_none0_cg_list3

END MODULE qs_harmonics_atom
