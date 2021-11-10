!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Print basis sets in CP2K format
!> \par History
!> \author JGH (12.2017)
! **************************************************************************************************
MODULE basis_set_output
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp2k_info,                       ONLY: compile_revision,&
                                              cp2k_version,&
                                              r_datx,&
                                              r_host_name,&
                                              r_user_name
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'basis_set_output'
   PUBLIC :: print_basis_set_file

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param base_section ...
! **************************************************************************************************
   SUBROUTINE print_basis_set_file(qs_env, base_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: base_section

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bname, filename
      INTEGER                                            :: ikind, iunit, nkind, ounit
      INTEGER, SAVE                                      :: ncalls = 0
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis, lri_aux_basis, orb_basis, &
                                                            ri_aux_basis, ri_hfx_basis, &
                                                            ri_hxc_basis, ri_xas_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (ncalls > 0) RETURN
      ncalls = ncalls + 1

      logger => cp_get_default_logger()
      ounit = cp_logger_get_default_io_unit(logger)

      CALL section_vals_val_get(base_section, "FILENAME", c_val=filename)

      IF (ounit > 0) THEN
         WRITE (UNIT=ounit, FMT='(/,(T2,A))') REPEAT("-", 79)
         WRITE (UNIT=ounit, FMT='((T2,A,A))') "Print Basis Set File:    ", TRIM(filename)
         WRITE (UNIT=ounit, FMT='((T2,A))') REPEAT("-", 79)
         CALL open_file(filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
         WRITE (UNIT=iunit, FMT="(A8,T11,A)") &
            "# TITLE ", "Basis set file created by "//TRIM(cp2k_version)//" (revision "//TRIM(compile_revision)//")", &
            "# AUTHOR", TRIM(r_user_name)//"@"//TRIM(r_host_name)//" "//r_datx(1:19)

      END IF

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, element_symbol=element_symbol)
         NULLIFY (orb_basis, ri_aux_basis, lri_aux_basis, aux_fit_basis)
         CALL get_qs_kind(qs_kind, basis_set=orb_basis, basis_type="ORB")
         CALL get_qs_kind(qs_kind, basis_set=ri_aux_basis, basis_type="RI_AUX")
         CALL get_qs_kind(qs_kind, basis_set=ri_hxc_basis, basis_type="RI_HXC")
         CALL get_qs_kind(qs_kind, basis_set=ri_hfx_basis, basis_type="RI_HFX")
         CALL get_qs_kind(qs_kind, basis_set=lri_aux_basis, basis_type="LRI_AUX")
         CALL get_qs_kind(qs_kind, basis_set=aux_fit_basis, basis_type="AUX_FIT")
         CALL get_qs_kind(qs_kind, basis_set=ri_xas_basis, basis_type="RI_XAS")
         IF (ounit > 0) THEN
            IF (ASSOCIATED(orb_basis)) THEN
               bname = "local_orbital"
               CALL basis_out(orb_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(ri_aux_basis)) THEN
               bname = "local_ri_aux"
               CALL basis_out(ri_aux_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(ri_hxc_basis)) THEN
               bname = "local_ri_hxc"
               CALL basis_out(ri_hxc_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(lri_aux_basis)) THEN
               bname = "local_lri_aux"
               CALL basis_out(lri_aux_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(aux_fit_basis)) THEN
               bname = "local_aux_fit"
               CALL basis_out(aux_fit_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(ri_xas_basis)) THEN
               bname = "local_ri_xas"
               CALL basis_out(ri_xas_basis, element_symbol, bname, iunit)
            END IF
            IF (ASSOCIATED(ri_hfx_basis)) THEN
               bname = "local_ri_hfx"
               CALL basis_out(ri_hfx_basis, element_symbol, bname, iunit)
            END IF
         END IF
      END DO

      IF (ounit > 0) THEN
         CALL close_file(iunit)
      END IF

   END SUBROUTINE print_basis_set_file

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param element_symbol ...
!> \param bname ...
!> \param iunit ...
! **************************************************************************************************
   SUBROUTINE basis_out(basis, element_symbol, bname, iunit)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol, bname
      INTEGER, INTENT(IN)                                :: iunit

      INTEGER                                            :: ipgf, iset, ishell, ll, nset
      INTEGER, DIMENSION(0:9)                            :: lset
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l, n
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc

      WRITE (iunit, "(A1)") "#"
      WRITE (iunit, "(A2,T5,A)") element_symbol, ADJUSTL(TRIM(bname))

      CALL get_gto_basis_set(basis, nset=nset, npgf=npgf, lmax=lmax, lmin=lmin, &
                             nshell=nshell, n=n, l=l, &
                             gcc=gcc, zet=zet)

      WRITE (iunit, "(I5)") nset
      DO iset = 1, nset
         lset = 0
         DO ishell = 1, nshell(iset)
            ll = l(ishell, iset)
            lset(ll) = lset(ll) + 1
         END DO
         WRITE (iunit, "(I5,2I3,I5,2X,10(I3))") n(1, iset), lmin(iset), lmax(iset), npgf(iset), &
            (lset(ll), ll=lmin(iset), lmax(iset))
         DO ipgf = 1, npgf(iset)
            WRITE (iunit, "(F20.10,50(F15.10))") zet(ipgf, iset), (gcc(ipgf, ishell, iset), ishell=1, nshell(iset))
         END DO
      END DO

   END SUBROUTINE basis_out

! **************************************************************************************************

END MODULE basis_set_output
