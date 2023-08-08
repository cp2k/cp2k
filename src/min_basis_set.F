!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief generate or use from input minimal basis set
!> \par History
!>      03.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE min_basis_set
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: &
        allocate_sto_basis_set, create_gto_from_sto_basis, deallocate_sto_basis_set, &
        get_gto_basis_set, gto_basis_set_type, init_orb_basis_set, set_sto_basis_set, srules, &
        sto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE periodic_table,                  ONLY: get_ptable_info,&
                                              ptable
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'min_basis_set'

   PUBLIC ::  create_minbas_set

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_minbas_set(qs_env, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bname
      INTEGER                                            :: ikind, lb, mao, ne, ngau, nkind, nsgf, &
                                                            ub, z
      INTEGER, DIMENSION(0:3)                            :: elcon
      INTEGER, DIMENSION(:), POINTER                     :: econf
      REAL(KIND=dp)                                      :: zval
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: refbasis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      ! check for or generate reference basis
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nkind = SIZE(qs_kind_set)
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         NULLIFY (refbasis)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=refbasis, basis_type="MIN")
         IF (.NOT. ASSOCIATED(refbasis)) THEN
            CALL get_qs_kind(qs_kind=qs_kind, mao=mao)
            ! generate a minimal basis set
            elcon = 0
            lb = LBOUND(econf, 1)
            ub = UBOUND(econf, 1)
            ne = ub - lb
            elcon(0:ne) = econf(lb:ub)
            ! STO-3G
            ngau = 3
            CALL create_min_basis(refbasis, z, elcon, mao, ngau)
            CALL init_interaction_radii_orb_basis(refbasis, dft_control%qs_control%eps_pgf_orb)
            CALL add_basis_set_to_container(qs_kind%basis_sets, refbasis, "MIN")
         END IF
         IF (unit_nr > 0) THEN
            CALL get_gto_basis_set(refbasis, name=bname, nsgf=nsgf)
            WRITE (unit_nr, '(T2,A,A,T14,A,I4,T40,A,A24)') "Kind: ", element_symbol, "NBasFun: ", nsgf, &
               "Reference Basis: ", ADJUSTL(TRIM(bname))
         END IF
      END DO

   END SUBROUTINE create_minbas_set

! **************************************************************************************************
!> \brief Creates a minimal basis set based on Slater rules
!> \param min_basis_set ...
!> \param zval ...
!> \param econf ...
!> \param mao ...
!> \param ngau ...
!> \par History
!>      03.2023 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_min_basis(min_basis_set, zval, econf, mao, ngau)
      TYPE(gto_basis_set_type), POINTER                  :: min_basis_set
      INTEGER, INTENT(IN)                                :: zval
      INTEGER, DIMENSION(0:3)                            :: econf
      INTEGER, INTENT(IN)                                :: mao, ngau

      CHARACTER(len=1), DIMENSION(0:3), PARAMETER        :: lnam = (/"S", "P", "D", "F"/)

      CHARACTER(len=6)                                   :: str
      CHARACTER(len=6), DIMENSION(:), POINTER            :: sym
      CHARACTER(LEN=default_string_length)               :: bname
      INTEGER                                            :: i, iss, l, n, nmao, nn, nss
      INTEGER, DIMENSION(0:3)                            :: nae, nco, npe
      INTEGER, DIMENSION(4, 7)                           :: ne
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      CPASSERT(.NOT. ASSOCIATED(min_basis_set))
      NULLIFY (sto_basis_set)

      ! electronic configuration
      ne = 0
      DO l = 1, 4
         nn = 2*(l - 1) + 1
         DO i = l, 7
            ne(l, i) = ptable(zval)%e_conv(l - 1) - 2*nn*(i - l)
            ne(l, i) = MAX(ne(l, i), 0)
            ne(l, i) = MIN(ne(l, i), 2*nn)
         END DO
      END DO
      ! STO definition
      nae = 0
      npe = 0
      DO l = 0, 3
         nn = 2*(2*l + 1)
         nae(l) = ptable(zval)%e_conv(l)/nn
         IF (MOD(ptable(zval)%e_conv(l), nn) /= 0) nae(l) = nae(l) + 1
         npe(l) = econf(l)/nn
         IF (MOD(econf(l), nn) /= 0) npe(l) = npe(l) + 1
      END DO
      CPASSERT(ALL(nae - npe >= 0))
      nco = nae - npe
      ! MAO count
      nmao = 0
      DO l = 0, 3
         nmao = nmao + npe(l)*(2*l + 1)
      END DO

      IF (mao > nmao) THEN
         nmao = mao - nmao
         SELECT CASE (nmao)
         CASE (1)
            npe(0) = npe(0) + 1
         CASE (2)
            npe(0) = npe(0) + 2
         CASE (3)
         CASE (4)
            npe(0) = npe(0) + 1
            npe(1) = npe(1) + 1
         CASE (5)
            IF (npe(2) == 0) THEN
               npe(2) = npe(2) + 1
            ELSE
               npe(0) = npe(0) + 2
               npe(1) = npe(1) + 1
            END IF
         CASE (6)
            npe(0) = npe(0) + 1
            npe(2) = npe(2) + 1
         CASE (7)
            npe(0) = npe(0) + 2
            npe(2) = npe(2) + 1
         CASE (8)
            npe(0) = npe(0) + 3
            npe(2) = npe(2) + 1
         CASE (9)
            npe(0) = npe(0) + 1
            npe(1) = npe(1) + 1
            npe(2) = npe(2) + 1
         CASE DEFAULT
            CPABORT("Default setting of minimal basis failed")
         END SELECT
         CALL cp_warn(__LOCATION__, "Reference basis has been adjusted according to MAO value.")
      END IF

      nss = SUM(npe)
      ALLOCATE (sym(nss), lq(nss), nq(nss), zet(nss))
      iss = 0
      DO l = 0, 3
         DO i = 1, npe(l)
            iss = iss + 1
            lq(iss) = l
            n = nco(l) + l
            nq(iss) = n + i
            str = "      "
            WRITE (str(5:5), FMT='(I1)') nq(iss)
            str(6:6) = lnam(l)
            sym(iss) = str
            zet(iss) = srules(zval, ne, nq(iss), lq(iss))
         END DO
      END DO

      bname = ADJUSTR(ptable(zval)%symbol)//"_MBas"
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, name=bname, nshell=nss, nq=nq, &
                             lq=lq, zet=zet, symbol=sym)
      CALL create_gto_from_sto_basis(sto_basis_set, min_basis_set, ngau)
      min_basis_set%norm_type = 2
      CALL init_orb_basis_set(min_basis_set)
      CALL deallocate_sto_basis_set(sto_basis_set)

      DEALLOCATE (sym, lq, nq, zet)

   END SUBROUTINE create_min_basis

! **************************************************************************************************

END MODULE min_basis_set
