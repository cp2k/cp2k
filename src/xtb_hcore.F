!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of EHT matrix elements in xTB
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_hcore
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE xtb_parameters,                  ONLY: early3d,&
                                              metal,&
                                              pp_gfn0,&
                                              xtb_set_kab
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_hcore'

   PUBLIC :: gfn0_huckel, gfn1_huckel, gfn0_kpair, gfn1_kpair

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cnumbers ...
!> \param charges ...
!> \param huckel ...
!> \param dhuckel ...
!> \param dqhuckel ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE gfn0_huckel(qs_env, cnumbers, charges, huckel, dhuckel, dqhuckel, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers, charges
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: huckel, dhuckel, dqhuckel
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: i, iatom, ikind, l, natom, nshell
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(25)                             :: lval
      REAL(KIND=dp)                                      :: kqat2
      REAL(KIND=dp), DIMENSION(5)                        :: hena, kcn, kq
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)
      xtb_control => dft_control%qs_control%xtb_control

      CALL get_qs_env(qs_env=qs_env, natom=natom)

      ALLOCATE (huckel(5, natom))
      IF (calculate_forces) THEN
         ALLOCATE (dhuckel(5, natom), dqhuckel(5, natom))
      END IF
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, nshell=nshell, lval=lval, &
                                 kcn=kcn, kq=kq, kqat2=kqat2, hen=hena)
         kcn = kcn/evolt
         kq = kq/evolt
         kqat2 = kqat2/evolt
         huckel(:, iatom) = 0.0_dp
         DO i = 1, nshell
            l = lval(i) + 1
            huckel(i, iatom) = hena(i) - kcn(l)*cnumbers(iatom) &
                               - kq(l)*charges(iatom) - kqat2*charges(iatom)**2
         END DO
         IF (calculate_forces) THEN
            dhuckel(:, iatom) = 0.0_dp
            dqhuckel(:, iatom) = 0.0_dp
            DO i = 1, nshell
               l = lval(i) + 1
               dhuckel(i, iatom) = -kcn(l)
               dqhuckel(i, iatom) = -kq(l) - 2.0_dp*kqat2*charges(iatom)
            END DO
         END IF
      END DO

   END SUBROUTINE gfn0_huckel

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cnumbers ...
!> \param huckel ...
!> \param dhuckel ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE gfn1_huckel(qs_env, cnumbers, huckel, dhuckel, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: huckel, dhuckel
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: i, iatom, ikind, natom, nkind, nshell, za
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(25)                             :: lval
      REAL(KIND=dp)                                      :: kcnd, kcnp, kcns
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: kcnlk
      REAL(KIND=dp), DIMENSION(5)                        :: hena
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)
      xtb_control => dft_control%qs_control%xtb_control

      CALL get_qs_env(qs_env=qs_env, nkind=nkind, natom=natom)

      kcns = xtb_control%kcns
      kcnp = xtb_control%kcnp
      kcnd = xtb_control%kcnd

      ! Calculate Huckel parameters
      ! Eq 12
      ! huckel(nshell,natom)
      ALLOCATE (kcnlk(0:3, nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         IF (metal(za)) THEN
            kcnlk(0:3, ikind) = 0.0_dp
         ELSEIF (early3d(za)) THEN
            kcnlk(0, ikind) = kcns
            kcnlk(1, ikind) = kcnp
            kcnlk(2, ikind) = 0.005_dp
            kcnlk(3, ikind) = 0.0_dp
         ELSE
            kcnlk(0, ikind) = kcns
            kcnlk(1, ikind) = kcnp
            kcnlk(2, ikind) = kcnd
            kcnlk(3, ikind) = 0.0_dp
         END IF
      END DO

      ALLOCATE (huckel(5, natom))
      IF (calculate_forces) THEN
         ALLOCATE (dhuckel(5, natom))
      END IF
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, nshell=nshell, lval=lval, hen=hena)
         huckel(:, iatom) = 0.0_dp
         DO i = 1, nshell
            huckel(i, iatom) = hena(i)*(1._dp + kcnlk(lval(i), ikind)*cnumbers(iatom))
         END DO
         IF (calculate_forces) THEN
            dhuckel(:, iatom) = 0.0_dp
            DO i = 1, nshell
               dhuckel(i, iatom) = hena(i)*kcnlk(lval(i), ikind)
            END DO
         END IF
      END DO

      DEALLOCATE (kcnlk)

   END SUBROUTINE gfn1_huckel

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kijab ...
! **************************************************************************************************
   SUBROUTINE gfn0_kpair(qs_env, kijab)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab

      INTEGER                                            :: i, ikind, j, jkind, la, lb, maxs, na, &
                                                            natorb_a, natorb_b, nb, nkind, za, zb
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      LOGICAL                                            :: defined
      REAL(KIND=dp)                                      :: ben, den, etaa, etab, kab, kd, kden, &
                                                            kdiff, ken, kia, kjb, km, kp, kpen, &
                                                            ks, ksen, ksp, xijab, yijab
      REAL(KIND=dp), DIMENSION(0:3)                      :: ke, kl
      REAL(KIND=dp), DIMENSION(5)                        :: zetaa, zetab
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)
      xtb_control => dft_control%qs_control%xtb_control

      CALL get_qs_env(qs_env=qs_env, nkind=nkind)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxsgf=maxs, basis_type="ORB")

      ks = xtb_control%ks
      kp = xtb_control%kp
      kd = xtb_control%kd
      ksp = xtb_control%ksp
      ksen = xtb_control%ksen
      kpen = xtb_control%kpen
      kden = xtb_control%kden
      ben = xtb_control%ben
      kdiff = xtb_control%k2sh

      kl(0) = ks
      kl(1) = kp
      kl(2) = kd
      kl(3) = 0.0_dp

      ke(0) = ksen
      ke(1) = kpen
      ke(2) = kden
      ke(3) = 0.0_dp

      ! Calculate KAB parameters and electronegativity correction
      ALLOCATE (kijab(maxs, maxs, nkind, nkind))
      kijab = 0.0_dp

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, &
                                 en=etaa, zeta=zetaa)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined, natorb=natorb_b)
            IF (.NOT. defined .OR. natorb_b < 1) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, &
                                    en=etab, zeta=zetab)
            ! Kab
            kab = pp_gfn0(za, zb)
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  kia = kl(la)
                  kjb = kl(lb)
                  km = 0.5_dp*(kia + kjb)*kab
                  IF (za == 1 .AND. na == 2) THEN
                     IF (zb == 1 .AND. nb == 2) THEN
                        km = 0._dp
                     ELSE
                        km = km*kdiff
                     END IF
                  ELSEIF (zb == 1 .AND. nb == 2) THEN
                     km = km*kdiff
                  END IF
                  kijab(i, j, ikind, jkind) = km
               END DO
            END DO
            ! Yab
            DO j = 1, natorb_b
               nb = naob(j)
               kjb = zetab(nb)
               DO i = 1, natorb_a
                  na = naoa(i)
                  kia = zetaa(na)
                  yijab = 2.0_dp*SQRT(kia*kjb)/(kia + kjb)
                  kijab(i, j, ikind, jkind) = kijab(i, j, ikind, jkind)*yijab
               END DO
            END DO
            ! X
            den = etaa - etab
            DO j = 1, natorb_b
               lb = laob(j)
               kjb = ke(lb)
               DO i = 1, natorb_a
                  la = laoa(i)
                  kia = ke(la)
                  ken = 0.5_dp*(kia + kjb)
                  xijab = 1.0_dp + ken*den**2 + ken*ben*den**4
                  kijab(i, j, ikind, jkind) = kijab(i, j, ikind, jkind)*xijab
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE gfn0_kpair

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kijab ...
! **************************************************************************************************
   SUBROUTINE gfn1_kpair(qs_env, kijab)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab

      INTEGER                                            :: i, ikind, j, jkind, la, lb, maxs, na, &
                                                            natorb_a, natorb_b, nb, nkind, za, zb
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      LOGICAL                                            :: defined
      REAL(KIND=dp)                                      :: ena, enb, fen, k2sh, kab, kd, ken, kia, &
                                                            kjb, kp, ks, ksp
      REAL(KIND=dp), DIMENSION(0:3)                      :: kl
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)
      xtb_control => dft_control%qs_control%xtb_control

      CALL get_qs_env(qs_env=qs_env, nkind=nkind)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxsgf=maxs, basis_type="ORB")

      ks = xtb_control%ks
      kp = xtb_control%kp
      kd = xtb_control%kd
      ksp = xtb_control%ksp
      k2sh = xtb_control%k2sh
      ken = xtb_control%ken

      kl(0) = ks
      kl(1) = kp
      kl(2) = kd
      kl(3) = 0.0_dp

      ! Calculate KAB parameters and electronegativity correction
      ! kijab -> K_l_l'[A,B] * X_l_l'[ENa, ENb] * Y[xia, xib]
      ALLOCATE (kijab(maxs, maxs, nkind, nkind))
      kijab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, electronegativity=ena)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined, natorb=natorb_b)
            IF (.NOT. defined .OR. natorb_b < 1) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, electronegativity=enb)
            ! get Fen = (1+ken*deltaEN^2)
            fen = 1.0_dp + ken*(ena - enb)**2
            ! Kab
            kab = xtb_set_kab(za, zb, xtb_control)
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  kia = kl(la)
                  kjb = kl(lb)
                  IF (zb == 1 .AND. nb == 2) kjb = k2sh
                  IF (za == 1 .AND. na == 2) kia = k2sh
                  IF ((zb == 1 .AND. nb == 2) .OR. (za == 1 .AND. na == 2)) THEN
                     kijab(i, j, ikind, jkind) = 0.5_dp*(kia + kjb)
                  ELSE
                     IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                        kijab(i, j, ikind, jkind) = ksp*kab*fen
                     ELSE
                        kijab(i, j, ikind, jkind) = 0.5_dp*(kia + kjb)*kab*fen
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE gfn1_kpair

END MODULE xtb_hcore

