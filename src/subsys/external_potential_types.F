!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the atomic potential types.
!> \par History
!>      GT, 22.09.2002: added elp_potential_types
!> \author Matthias Krack (04.07.2000)
! **************************************************************************************************
MODULE external_potential_types
   USE ao_util,                         ONLY: exp_radius
   USE bibliography,                    ONLY: Goedecker1996,&
                                              Hartwigsen1998,&
                                              Krack2000,&
                                              Krack2005,&
                                              cite_reference
   USE cp_linked_list_input,            ONLY: cp_sll_val_next,&
                                              cp_sll_val_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object,&
                                              parser_search_string,&
                                              parser_test_next_token
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_list_get,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE input_val_types,                 ONLY: val_get,&
                                              val_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              pi,&
                                              rootpi
   USE mathlib,                         ONLY: symmetrize_matrix
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: co,&
                                              coset,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset,&
                                              nso
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE periodic_table,                  ONLY: ptable
   USE string_utilities,                ONLY: remove_word,&
                                              uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_types'

   ! Define the all-electron potential type

   ! Literature: M. Krack and M. Parrinello,
   !             Phys. Chem. Chem. Phys. 2, 2105 (2000)

! **************************************************************************************************
   TYPE all_potential_type
      !MK PRIVATE
      CHARACTER(LEN=default_string_length)   :: name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(2)                      :: description
      REAL(KIND=dp)                        :: alpha_core_charge, &
                                              ccore_charge, &
                                              core_charge_radius, &
                                              zeff, zeff_correction
      INTEGER                                :: z
      INTEGER, DIMENSION(:), POINTER         :: elec_conf
   END TYPE all_potential_type

   ! Define the effective charge & inducible dipole potential type (for Fist)
! **************************************************************************************************
   TYPE fist_potential_type
      PRIVATE
      CHARACTER(LEN=default_string_length)     :: name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(1)                        :: description
      REAL(KIND=dp)                          :: apol, cpol, mm_radius, qeff, &
                                                qmmm_corr_radius, qmmm_radius

   END TYPE fist_potential_type

! **************************************************************************************************
! local potential type
! V(r) = SUM_i exp(0.5*(r/rci)**2) * ( C1i + C2i (r/rci)**2 + C3i (r/rci)**4 ...)
! alpha = 0.5/rci**2
   TYPE local_potential_type
      !PRIVATE
      CHARACTER(LEN=default_string_length)       :: name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(4)                          :: description
      INTEGER                                    :: ngau, npol
      REAL(KIND=dp)                            :: radius
      REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: cval
   END TYPE local_potential_type

   ! Define the GTH potential type

   ! Literature: - S. Goedecker, M. Teter and J. Hutter,
   !               Phys. Rev. B 54, 1703 (1996)
   !             - C. Hartwigsen, S. Goedecker and J. Hutter,
   !               Phys. Rev. B 58, 3641 (1998)
   !             - M. Krack,
   !               Theor. Chem. Acc. 114, 145 (2005)

! **************************************************************************************************
   TYPE gth_potential_type
      !PRIVATE
      CHARACTER(LEN=default_string_length)       :: name
      CHARACTER(LEN=default_string_length)       :: aliases = ""
      CHARACTER(LEN=default_string_length), &
         DIMENSION(4)                          :: description
      REAL(KIND=dp)                            :: alpha_core_charge, &
                                                  alpha_ppl, ccore_charge, &
                                                  cerf_ppl, zeff, &
                                                  core_charge_radius, &
                                                  ppl_radius, ppnl_radius, &
                                                  zeff_correction
      INTEGER                                    :: lppnl, lprj_ppnl_max, &
                                                    nexp_ppl, nppnl, &
                                                    nprj_ppnl_max, z
      REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_ppnl, cexp_ppl
      INTEGER, DIMENSION(:), POINTER             :: elec_conf
      ! nonlocal projectors
      INTEGER, DIMENSION(:), POINTER             :: nprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: cprj, cprj_ppnl, vprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: hprj_ppnl
      ! type extensions
      ! NLCC
      LOGICAL                                    :: nlcc
      INTEGER                                    :: nexp_nlcc
      REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_nlcc
      INTEGER, DIMENSION(:), POINTER             :: nct_nlcc
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: cval_nlcc
      ! LSD potential
      LOGICAL                                    :: lsdpot
      INTEGER                                    :: nexp_lsd
      REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lsd
      INTEGER, DIMENSION(:), POINTER             :: nct_lsd
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: cval_lsd
      ! extended local potential
      LOGICAL                                    :: lpotextended
      INTEGER                                    :: nexp_lpot
      REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lpot
      INTEGER, DIMENSION(:), POINTER             :: nct_lpot
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: cval_lpot

   END TYPE gth_potential_type

! **************************************************************************************************
   TYPE sgp_potential_type
      CHARACTER(LEN=default_string_length)       :: name
      CHARACTER(LEN=default_string_length)       :: aliases = ""
      CHARACTER(LEN=default_string_length), &
         DIMENSION(4)                            :: description
      ! CHARGE
      INTEGER                                    :: z
      REAL(KIND=dp)                              :: zeff, &
                                                    zeff_correction
      REAL(KIND=dp)                              :: alpha_core_charge, &
                                                    ccore_charge, &
                                                    core_charge_radius
      REAL(KIND=dp)                              :: ppl_radius, ppnl_radius
      INTEGER, DIMENSION(:), POINTER             :: elec_conf
      ! LOCAL
      LOGICAL                                    :: ecp_local
      INTEGER                                    :: n_local
      REAL(KIND=dp), DIMENSION(:), POINTER       :: a_local => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c_local => Null()
      ! ECP local
      INTEGER                                    :: nloc ! # terms
      INTEGER, DIMENSION(1:10)                   :: nrloc ! r**(n-2)
      REAL(dp), DIMENSION(1:10)                  :: aloc ! coefficient
      REAL(dp), DIMENSION(1:10)                  :: bloc ! exponent
      ! NONLOCAL
      INTEGER                                    :: n_nonlocal
      INTEGER                                    :: nppnl
      INTEGER                                    :: lmax
      LOGICAL, DIMENSION(0:5)                    :: is_nonlocal = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER       :: a_nonlocal => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: h_nonlocal => Null()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: c_nonlocal => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: cprj_ppnl
      REAL(KIND=dp), DIMENSION(:), POINTER       :: vprj_ppnl
      ! NLCC
      LOGICAL                                    :: has_nlcc
      INTEGER                                    :: n_nlcc
      REAL(KIND=dp), DIMENSION(:), POINTER       :: a_nlcc => Null()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c_nlcc => Null()
   END TYPE sgp_potential_type

! **************************************************************************************************
   TYPE all_potential_p_type
      TYPE(all_potential_type), POINTER          :: all_potential
   END TYPE all_potential_p_type

   TYPE gth_potential_p_type
      TYPE(gth_potential_type), POINTER          :: gth_potential
   END TYPE gth_potential_p_type

   TYPE local_potential_p_type
      TYPE(local_potential_type), POINTER        :: local_potential
   END TYPE local_potential_p_type

   TYPE sgp_potential_p_type
      TYPE(sgp_potential_type), POINTER          :: sgp_potential
   END TYPE sgp_potential_p_type
! **************************************************************************************************

   ! Public subroutines
   PUBLIC :: allocate_potential, &
             deallocate_potential, &
             get_potential, &
             init_potential, &
             read_potential, &
             set_potential, &
             set_default_all_potential, &
             write_potential, &
             copy_potential

   ! Public data types

   PUBLIC :: all_potential_type, &
             fist_potential_type, &
             local_potential_type, &
             gth_potential_type, &
             sgp_potential_type
   PUBLIC :: gth_potential_p_type, &
             sgp_potential_p_type

   INTERFACE allocate_potential
      MODULE PROCEDURE allocate_all_potential, &
         allocate_fist_potential, &
         allocate_local_potential, &
         allocate_gth_potential, &
         allocate_sgp_potential
   END INTERFACE

   INTERFACE deallocate_potential
      MODULE PROCEDURE deallocate_all_potential, &
         deallocate_fist_potential, &
         deallocate_local_potential, &
         deallocate_sgp_potential, &
         deallocate_gth_potential
   END INTERFACE

   INTERFACE get_potential
      MODULE PROCEDURE get_all_potential, &
         get_fist_potential, &
         get_local_potential, &
         get_gth_potential, &
         get_sgp_potential
   END INTERFACE

   INTERFACE init_potential
      MODULE PROCEDURE init_all_potential, &
         init_gth_potential, &
         init_sgp_potential
   END INTERFACE

   INTERFACE read_potential
      MODULE PROCEDURE read_all_potential, &
         read_local_potential, &
         read_gth_potential
   END INTERFACE

   INTERFACE set_potential
      MODULE PROCEDURE set_all_potential, &
         set_fist_potential, &
         set_local_potential, &
         set_gth_potential, &
         set_sgp_potential
   END INTERFACE

   INTERFACE write_potential
      MODULE PROCEDURE write_all_potential, &
         write_local_potential, &
         write_gth_potential, &
         write_sgp_potential
   END INTERFACE

   INTERFACE copy_potential
      MODULE PROCEDURE copy_all_potential, &
         copy_gth_potential, &
         copy_sgp_potential
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief   Allocate an atomic all-electron potential data set.
!> \param potential ...
!> \date    25.07.2000,
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_all_potential(potential)
      TYPE(all_potential_type), INTENT(INOUT), POINTER   :: potential

      IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

      ALLOCATE (potential)

      NULLIFY (potential%elec_conf)

      potential%description(1) = "All-electron potential"
      potential%description(2) = "Krack, Parrinello, PCCP 2, 2105 (2000)"

   END SUBROUTINE allocate_all_potential

! **************************************************************************************************
!> \brief   Allocate an effective charge and inducible dipole potential data set.
!> \param potential ...
!> \date    05.03.2010
!> \author  Toon.Verstraelen@gmail.com
! **************************************************************************************************
   SUBROUTINE allocate_fist_potential(potential)
      TYPE(fist_potential_type), INTENT(INOUT), POINTER  :: potential

      IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

      ALLOCATE (potential)

      potential%apol = 0.0_dp
      potential%cpol = 0.0_dp
      potential%mm_radius = 0.0_dp
      potential%qeff = 0.0_dp
      potential%qmmm_radius = 0.0_dp
      potential%qmmm_corr_radius = 0.0_dp

      potential%description(1) = "Effective charge and inducible dipole potential"

   END SUBROUTINE allocate_fist_potential

! **************************************************************************************************
!> \brief   Allocate an atomic local potential data set.
!> \param potential ...
!> \date    24.01.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_local_potential(potential)
      TYPE(local_potential_type), INTENT(INOUT), POINTER :: potential

      IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

      ALLOCATE (potential)

      NULLIFY (potential%alpha)
      NULLIFY (potential%cval)

      potential%description(1) = "Local short-range pseudopotential"
      potential%ngau = 0
      potential%npol = 0
      potential%radius = 0.0_dp

   END SUBROUTINE allocate_local_potential

! **************************************************************************************************
!> \brief   Allocate an atomic GTH potential data set.
!> \param potential ...
!> \date    25.07.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_gth_potential(potential)
      TYPE(gth_potential_type), INTENT(INOUT), POINTER   :: potential

      IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

      ALLOCATE (potential)

      NULLIFY (potential%alpha_ppnl)
      NULLIFY (potential%cexp_ppl)
      NULLIFY (potential%elec_conf)
      NULLIFY (potential%nprj_ppnl)
      NULLIFY (potential%cprj)
      NULLIFY (potential%cprj_ppnl)
      NULLIFY (potential%vprj_ppnl)
      NULLIFY (potential%hprj_ppnl)

      NULLIFY (potential%alpha_lpot)
      NULLIFY (potential%nct_lpot)
      NULLIFY (potential%cval_lpot)
      NULLIFY (potential%alpha_lsd)
      NULLIFY (potential%nct_lsd)
      NULLIFY (potential%cval_lsd)
      NULLIFY (potential%alpha_nlcc)
      NULLIFY (potential%nct_nlcc)
      NULLIFY (potential%cval_nlcc)

      potential%description(1) = "Goedecker-Teter-Hutter pseudopotential"
      potential%description(2) = "Goedecker et al., PRB 54, 1703 (1996)"
      potential%description(3) = "Hartwigsen et al., PRB 58, 3641 (1998)"
      potential%description(4) = "Krack, TCA 114, 145 (2005)"

   END SUBROUTINE allocate_gth_potential

! **************************************************************************************************
!> \brief   Allocate an atomic SGP potential data set.
!> \param potential ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_sgp_potential(potential)
      TYPE(sgp_potential_type), INTENT(INOUT), POINTER   :: potential

      IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

      ALLOCATE (potential)

      NULLIFY (potential%elec_conf)
      NULLIFY (potential%a_local)
      NULLIFY (potential%c_local)
      NULLIFY (potential%a_nonlocal)
      NULLIFY (potential%c_nonlocal)
      NULLIFY (potential%h_nonlocal)
      NULLIFY (potential%a_nlcc)
      NULLIFY (potential%c_nlcc)
      NULLIFY (potential%cprj_ppnl)
      NULLIFY (potential%vprj_ppnl)

      potential%name = ""
      potential%aliases = ""
      potential%description(1) = "Separable Gaussian pseudopotential"
      potential%description(2) = "M. Pelissier, N. Komiha, J.P. Daudey, JCC, 9, 298 (1988)"
      potential%description(3) = "create from"
      potential%description(4) = ""

      potential%z = 0
      potential%zeff = 0.0_dp
      potential%zeff_correction = 0.0_dp
      potential%alpha_core_charge = 0.0_dp
      potential%ccore_charge = 0.0_dp
      potential%core_charge_radius = 0.0_dp
      potential%ppl_radius = 0.0_dp
      potential%ppnl_radius = 0.0_dp
      potential%ecp_local = .FALSE.
      potential%n_local = 0
      potential%nloc = 0
      potential%nrloc = 0
      potential%aloc = 0.0_dp
      potential%bloc = 0.0_dp
      potential%n_nonlocal = 0
      potential%nppnl = 0
      potential%lmax = -1
      potential%is_nonlocal = .FALSE.
      potential%has_nlcc = .FALSE.
      potential%n_nlcc = 0

   END SUBROUTINE allocate_sgp_potential
! **************************************************************************************************
!> \brief   Deallocate an atomic all-electron potential data set.
!> \param potential ...
!> \date    03.11.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_all_potential(potential)
      TYPE(all_potential_type), POINTER                  :: potential

      IF (.NOT. ASSOCIATED(potential)) THEN
         CPABORT("The pointer potential is not associated.")
      END IF

      DEALLOCATE (potential%elec_conf)
      DEALLOCATE (potential)

   END SUBROUTINE deallocate_all_potential

! **************************************************************************************************
!> \brief   Deallocate an effective charge and inducible dipole potential data set.
!> \param potential ...
!> \date    05.03.2010
!> \author  Toon.Verstraelen@gmail.com
! **************************************************************************************************
   SUBROUTINE deallocate_fist_potential(potential)
      TYPE(fist_potential_type), POINTER                 :: potential

      IF (.NOT. ASSOCIATED(potential)) THEN
         CPABORT("The pointer potential is not associated.")
      END IF

      ! Nothing exciting here yet.
      DEALLOCATE (potential)

   END SUBROUTINE deallocate_fist_potential

! **************************************************************************************************
!> \brief   Deallocate an atomic local potential data set.
!> \param potential ...
!> \date    24.01.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_local_potential(potential)
      TYPE(local_potential_type), POINTER                :: potential

      IF (.NOT. ASSOCIATED(potential)) THEN
         CPABORT("The pointer potential is not associated.")
      END IF

      IF (ASSOCIATED(potential%alpha)) THEN
         DEALLOCATE (potential%alpha)
      END IF
      IF (ASSOCIATED(potential%cval)) THEN
         DEALLOCATE (potential%cval)
      END IF

      DEALLOCATE (potential)

   END SUBROUTINE deallocate_local_potential

! **************************************************************************************************
!> \brief   Deallocate an atomic GTH potential data set.
!> \param potential ...
!> \date    03.11.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_gth_potential(potential)
      TYPE(gth_potential_type), POINTER                  :: potential

      IF (.NOT. ASSOCIATED(potential)) THEN
         CPABORT("The pointer potential is not associated.")
      END IF

      DEALLOCATE (potential%elec_conf)
      !     *** Deallocate the parameters of the local part ***

      IF (ASSOCIATED(potential%cexp_ppl)) THEN
         DEALLOCATE (potential%cexp_ppl)
      END IF

      !     *** Deallocate the parameters of the non-local part ***
      IF (ASSOCIATED(potential%alpha_ppnl)) THEN
         DEALLOCATE (potential%alpha_ppnl)
         DEALLOCATE (potential%cprj)
         DEALLOCATE (potential%cprj_ppnl)
         DEALLOCATE (potential%hprj_ppnl)
         DEALLOCATE (potential%nprj_ppnl)
         DEALLOCATE (potential%vprj_ppnl)
      END IF

      IF (ASSOCIATED(potential%alpha_lpot)) THEN
         DEALLOCATE (potential%alpha_lpot)
         DEALLOCATE (potential%nct_lpot)
         DEALLOCATE (potential%cval_lpot)
      END IF

      IF (ASSOCIATED(potential%alpha_lsd)) THEN
         DEALLOCATE (potential%alpha_lsd)
         DEALLOCATE (potential%nct_lsd)
         DEALLOCATE (potential%cval_lsd)
      END IF

      IF (ASSOCIATED(potential%alpha_nlcc)) THEN
         DEALLOCATE (potential%alpha_nlcc)
         DEALLOCATE (potential%nct_nlcc)
         DEALLOCATE (potential%cval_nlcc)
      END IF

      DEALLOCATE (potential)
   END SUBROUTINE deallocate_gth_potential

! **************************************************************************************************
!> \brief   Deallocate an atomic SGP potential data set.
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE deallocate_sgp_potential(potential)
      TYPE(sgp_potential_type), POINTER                  :: potential

      IF (.NOT. ASSOCIATED(potential)) THEN
         CPABORT("The pointer potential is not associated.")
      END IF

      IF (ASSOCIATED(potential%elec_conf)) THEN
         DEALLOCATE (potential%elec_conf)
      END IF
      IF (ASSOCIATED(potential%a_local)) THEN
         DEALLOCATE (potential%a_local)
      END IF
      IF (ASSOCIATED(potential%c_local)) THEN
         DEALLOCATE (potential%c_local)
      END IF

      IF (ASSOCIATED(potential%a_nonlocal)) THEN
         DEALLOCATE (potential%a_nonlocal)
      END IF
      IF (ASSOCIATED(potential%h_nonlocal)) THEN
         DEALLOCATE (potential%h_nonlocal)
      END IF
      IF (ASSOCIATED(potential%c_nonlocal)) THEN
         DEALLOCATE (potential%c_nonlocal)
      END IF
      IF (ASSOCIATED(potential%cprj_ppnl)) THEN
         DEALLOCATE (potential%cprj_ppnl)
      END IF
      IF (ASSOCIATED(potential%vprj_ppnl)) THEN
         DEALLOCATE (potential%vprj_ppnl)
      END IF

      IF (ASSOCIATED(potential%a_nlcc)) THEN
         DEALLOCATE (potential%a_nlcc)
      END IF
      IF (ASSOCIATED(potential%c_nlcc)) THEN
         DEALLOCATE (potential%c_nlcc)
      END IF

      DEALLOCATE (potential)

   END SUBROUTINE deallocate_sgp_potential

! **************************************************************************************************
!> \brief   Get attributes of an all-electron potential data set.
!> \param potential ...
!> \param name ...
!> \param alpha_core_charge ...
!> \param ccore_charge ...
!> \param core_charge_radius ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param elec_conf ...
!> \date    11.01.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_all_potential(potential, name, alpha_core_charge, &
                                ccore_charge, core_charge_radius, z, zeff, &
                                zeff_correction, elec_conf)
      TYPE(all_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: alpha_core_charge, ccore_charge, &
                                                            core_charge_radius
      INTEGER, INTENT(OUT), OPTIONAL                     :: z
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: zeff, zeff_correction
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) name = potential%name
         IF (PRESENT(alpha_core_charge)) &
            alpha_core_charge = potential%alpha_core_charge
         IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
         IF (PRESENT(core_charge_radius)) &
            core_charge_radius = potential%core_charge_radius
         IF (PRESENT(z)) z = potential%z
         IF (PRESENT(zeff)) zeff = potential%zeff
         IF (PRESENT(zeff_correction)) zeff_correction = potential%zeff_correction
         IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_all_potential

! **************************************************************************************************
!> \brief   Get attributes of an effective point charge and inducible dipole
!>          potential.
!> \param potential ...
!> \param name ...
!> \param apol ...
!> \param cpol ...
!> \param mm_radius ...
!> \param qeff ...
!> \param qmmm_corr_radius ...
!> \param qmmm_radius ...
!> \date    05.03-2010
!> \author  Toon.Verstraelen@UGent.be
! **************************************************************************************************
   SUBROUTINE get_fist_potential(potential, name, apol, cpol, mm_radius, qeff, &
                                 qmmm_corr_radius, qmmm_radius)
      TYPE(fist_potential_type), POINTER                 :: potential
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: apol, cpol, mm_radius, qeff, &
                                                            qmmm_corr_radius, qmmm_radius

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) name = potential%name
         IF (PRESENT(apol)) apol = potential%apol
         IF (PRESENT(cpol)) cpol = potential%cpol
         IF (PRESENT(mm_radius)) mm_radius = potential%mm_radius
         IF (PRESENT(qeff)) qeff = potential%qeff
         IF (PRESENT(qmmm_corr_radius)) qmmm_corr_radius = potential%qmmm_corr_radius
         IF (PRESENT(qmmm_radius)) qmmm_radius = potential%qmmm_radius

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_fist_potential

! **************************************************************************************************
!> \brief   Get attributes of an atomic local potential data set.
!> \param potential ...
!> \param name ...
!> \param ngau ...
!> \param npol ...
!> \param alpha ...
!> \param cval ...
!> \param radius ...
!> \date    24.01.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_local_potential(potential, name, ngau, npol, alpha, cval, radius)
      TYPE(local_potential_type), POINTER                :: potential
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name
      INTEGER, INTENT(OUT), OPTIONAL                     :: ngau, npol
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cval
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: radius

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) name = potential%name
         IF (PRESENT(ngau)) ngau = potential%ngau
         IF (PRESENT(npol)) npol = potential%npol
         IF (PRESENT(alpha)) alpha => potential%alpha
         IF (PRESENT(cval)) cval => potential%cval
         IF (PRESENT(radius)) radius = potential%radius

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_local_potential

! **************************************************************************************************
!> \brief   Get attributes of a GTH potential data set.
!> \param potential ...
!> \param name ...
!> \param aliases ...
!> \param alpha_core_charge ...
!> \param alpha_ppl ...
!> \param ccore_charge ...
!> \param cerf_ppl ...
!> \param core_charge_radius ...
!> \param ppl_radius ...
!> \param ppnl_radius ...
!> \param lppnl ...
!> \param lprj_ppnl_max ...
!> \param nexp_ppl ...
!> \param nppnl ...
!> \param nprj_ppnl_max ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param ppl_present ...
!> \param ppnl_present ...
!> \param alpha_ppnl ...
!> \param cexp_ppl ...
!> \param elec_conf ...
!> \param nprj_ppnl ...
!> \param cprj ...
!> \param cprj_ppnl ...
!> \param vprj_ppnl ...
!> \param hprj_ppnl ...
!> \param lpot_present ...
!> \param nexp_lpot ...
!> \param alpha_lpot ...
!> \param nct_lpot ...
!> \param cval_lpot ...
!> \param lsd_present ...
!> \param nexp_lsd ...
!> \param alpha_lsd ...
!> \param nct_lsd ...
!> \param cval_lsd ...
!> \param nlcc_present ...
!> \param nexp_nlcc ...
!> \param alpha_nlcc ...
!> \param nct_nlcc ...
!> \param cval_nlcc ...
!> \date    11.01.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_gth_potential(potential, name, aliases, alpha_core_charge, &
                                alpha_ppl, ccore_charge, cerf_ppl, &
                                core_charge_radius, ppl_radius, ppnl_radius, &
                                lppnl, lprj_ppnl_max, nexp_ppl, nppnl, &
                                nprj_ppnl_max, z, zeff, zeff_correction, &
                                ppl_present, ppnl_present, &
                                alpha_ppnl, cexp_ppl, elec_conf, nprj_ppnl, cprj, &
                                cprj_ppnl, vprj_ppnl, hprj_ppnl, &
                                lpot_present, nexp_lpot, alpha_lpot, nct_lpot, cval_lpot, &
                                lsd_present, nexp_lsd, alpha_lsd, nct_lsd, cval_lsd, &
                                nlcc_present, nexp_nlcc, alpha_nlcc, nct_nlcc, cval_nlcc)

      TYPE(gth_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name, aliases
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: alpha_core_charge, alpha_ppl, &
                                                            ccore_charge, cerf_ppl, &
                                                            core_charge_radius, ppl_radius, &
                                                            ppnl_radius
      INTEGER, INTENT(OUT), OPTIONAL                     :: lppnl, lprj_ppnl_max, nexp_ppl, nppnl, &
                                                            nprj_ppnl_max, z
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: zeff, zeff_correction
      LOGICAL, INTENT(OUT), OPTIONAL                     :: ppl_present, ppnl_present
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha_ppnl, cexp_ppl
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf, nprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cprj, cprj_ppnl, vprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: hprj_ppnl
      LOGICAL, INTENT(OUT), OPTIONAL                     :: lpot_present
      INTEGER, INTENT(OUT), OPTIONAL                     :: nexp_lpot
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha_lpot
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: nct_lpot
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cval_lpot
      LOGICAL, INTENT(OUT), OPTIONAL                     :: lsd_present
      INTEGER, INTENT(OUT), OPTIONAL                     :: nexp_lsd
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha_lsd
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: nct_lsd
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cval_lsd
      LOGICAL, INTENT(OUT), OPTIONAL                     :: nlcc_present
      INTEGER, INTENT(OUT), OPTIONAL                     :: nexp_nlcc
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha_nlcc
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: nct_nlcc
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cval_nlcc

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) name = potential%name
         IF (PRESENT(aliases)) aliases = potential%aliases
         IF (PRESENT(alpha_core_charge)) &
            alpha_core_charge = potential%alpha_core_charge
         IF (PRESENT(alpha_ppl)) alpha_ppl = potential%alpha_ppl
         IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
         IF (PRESENT(cerf_ppl)) cerf_ppl = potential%cerf_ppl
         IF (PRESENT(core_charge_radius)) &
            core_charge_radius = potential%core_charge_radius
         IF (PRESENT(ppl_radius)) ppl_radius = potential%ppl_radius
         IF (PRESENT(ppnl_radius)) ppnl_radius = potential%ppnl_radius
         IF (PRESENT(lppnl)) lppnl = potential%lppnl
         IF (PRESENT(lprj_ppnl_max)) lprj_ppnl_max = potential%lprj_ppnl_max
         IF (PRESENT(nexp_ppl)) nexp_ppl = potential%nexp_ppl
         IF (PRESENT(nppnl)) nppnl = potential%nppnl
         IF (PRESENT(nprj_ppnl_max)) nprj_ppnl_max = potential%nprj_ppnl_max
         IF (PRESENT(z)) z = potential%z
         IF (PRESENT(zeff)) zeff = potential%zeff
         IF (PRESENT(zeff_correction)) zeff_correction = potential%zeff_correction
         IF (PRESENT(ppl_present)) ppl_present = (potential%nexp_ppl > 0)
         IF (PRESENT(ppnl_present)) ppnl_present = (potential%nppnl > 0)
         IF (PRESENT(alpha_ppnl)) alpha_ppnl => potential%alpha_ppnl
         IF (PRESENT(cexp_ppl)) cexp_ppl => potential%cexp_ppl
         IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf
         IF (PRESENT(nprj_ppnl)) nprj_ppnl => potential%nprj_ppnl
         IF (PRESENT(cprj)) cprj => potential%cprj
         IF (PRESENT(cprj_ppnl)) cprj_ppnl => potential%cprj_ppnl
         IF (PRESENT(vprj_ppnl)) vprj_ppnl => potential%vprj_ppnl
         IF (PRESENT(hprj_ppnl)) hprj_ppnl => potential%hprj_ppnl

         IF (PRESENT(lpot_present)) lpot_present = potential%lpotextended
         IF (PRESENT(nexp_lpot)) nexp_lpot = potential%nexp_lpot
         IF (PRESENT(alpha_lpot)) alpha_lpot => potential%alpha_lpot
         IF (PRESENT(nct_lpot)) nct_lpot => potential%nct_lpot
         IF (PRESENT(cval_lpot)) cval_lpot => potential%cval_lpot

         IF (PRESENT(lsd_present)) lsd_present = potential%lsdpot
         IF (PRESENT(nexp_lsd)) nexp_lsd = potential%nexp_lsd
         IF (PRESENT(alpha_lsd)) alpha_lsd => potential%alpha_lsd
         IF (PRESENT(nct_lsd)) nct_lsd => potential%nct_lsd
         IF (PRESENT(cval_lsd)) cval_lsd => potential%cval_lsd

         IF (PRESENT(nlcc_present)) nlcc_present = potential%nlcc
         IF (PRESENT(nexp_nlcc)) nexp_nlcc = potential%nexp_nlcc
         IF (PRESENT(alpha_nlcc)) alpha_nlcc => potential%alpha_nlcc
         IF (PRESENT(nct_nlcc)) nct_nlcc => potential%nct_nlcc
         IF (PRESENT(cval_nlcc)) cval_nlcc => potential%cval_nlcc

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_gth_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param name ...
!> \param description ...
!> \param aliases ...
!> \param elec_conf ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param alpha_core_charge ...
!> \param ccore_charge ...
!> \param core_charge_radius ...
!> \param ppl_radius ...
!> \param ppnl_radius ...
!> \param ppl_present ...
!> \param ppnl_present ...
!> \param ecp_local ...
!> \param n_local ...
!> \param a_local ...
!> \param c_local ...
!> \param nloc ...
!> \param nrloc ...
!> \param aloc ...
!> \param bloc ...
!> \param n_nonlocal ...
!> \param nppnl ...
!> \param lmax ...
!> \param is_nonlocal ...
!> \param a_nonlocal ...
!> \param h_nonlocal ...
!> \param c_nonlocal ...
!> \param cprj_ppnl ...
!> \param vprj_ppnl ...
!> \param has_nlcc ...
!> \param n_nlcc ...
!> \param a_nlcc ...
!> \param c_nlcc ...
! **************************************************************************************************
   SUBROUTINE get_sgp_potential(potential, name, description, aliases, elec_conf, &
                                z, zeff, zeff_correction, alpha_core_charge, &
                                ccore_charge, core_charge_radius, &
                                ppl_radius, ppnl_radius, ppl_present, ppnl_present, &
                                ecp_local, n_local, a_local, c_local, &
                                nloc, nrloc, aloc, bloc, &
                                n_nonlocal, nppnl, lmax, is_nonlocal, a_nonlocal, h_nonlocal, c_nonlocal, &
                                cprj_ppnl, vprj_ppnl, has_nlcc, n_nlcc, a_nlcc, c_nlcc)

      TYPE(sgp_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(4), INTENT(OUT), OPTIONAL             :: description
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: aliases
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf
      INTEGER, INTENT(OUT), OPTIONAL                     :: z
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: zeff, zeff_correction, &
                                                            alpha_core_charge, ccore_charge, &
                                                            core_charge_radius, ppl_radius, &
                                                            ppnl_radius
      LOGICAL, INTENT(OUT), OPTIONAL                     :: ppl_present, ppnl_present, ecp_local
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_local
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_local, c_local
      INTEGER, INTENT(OUT), OPTIONAL                     :: nloc
      INTEGER, DIMENSION(1:10), INTENT(OUT), OPTIONAL    :: nrloc
      REAL(dp), DIMENSION(1:10), INTENT(OUT), OPTIONAL   :: aloc, bloc
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_nonlocal, nppnl, lmax
      LOGICAL, DIMENSION(0:5), OPTIONAL                  :: is_nonlocal
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_nonlocal
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: h_nonlocal
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: c_nonlocal
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cprj_ppnl
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: vprj_ppnl
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_nlcc
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_nlcc
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_nlcc, c_nlcc

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) name = potential%name
         IF (PRESENT(aliases)) aliases = potential%aliases
         IF (PRESENT(description)) description = potential%description

         IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf

         IF (PRESENT(z)) z = potential%z
         IF (PRESENT(zeff)) zeff = potential%zeff
         IF (PRESENT(zeff_correction)) zeff_correction = potential%zeff_correction
         IF (PRESENT(alpha_core_charge)) alpha_core_charge = potential%alpha_core_charge
         IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
         IF (PRESENT(core_charge_radius)) core_charge_radius = potential%core_charge_radius

         IF (PRESENT(ppl_radius)) ppl_radius = potential%ppl_radius
         IF (PRESENT(ppnl_radius)) ppnl_radius = potential%ppnl_radius
         IF (PRESENT(ppl_present)) THEN
            ppl_present = (potential%nloc > 0 .OR. potential%n_local > 0)
         END IF
         IF (PRESENT(ppnl_present)) THEN
            ppnl_present = ANY(potential%is_nonlocal)
         END IF

         IF (PRESENT(ecp_local)) ecp_local = potential%ecp_local
         IF (PRESENT(n_local)) n_local = potential%n_local
         IF (PRESENT(a_local)) a_local => potential%a_local
         IF (PRESENT(c_local)) c_local => potential%c_local

         IF (PRESENT(nloc)) nloc = potential%nloc
         IF (PRESENT(nrloc)) nrloc = potential%nrloc
         IF (PRESENT(aloc)) aloc = potential%aloc
         IF (PRESENT(bloc)) bloc = potential%bloc

         IF (PRESENT(n_nonlocal)) n_nonlocal = potential%n_nonlocal
         IF (PRESENT(nppnl)) nppnl = potential%nppnl
         IF (PRESENT(lmax)) lmax = potential%lmax
         IF (PRESENT(is_nonlocal)) is_nonlocal(:) = potential%is_nonlocal(:)
         IF (PRESENT(a_nonlocal)) a_nonlocal => potential%a_nonlocal
         IF (PRESENT(c_nonlocal)) c_nonlocal => potential%c_nonlocal
         IF (PRESENT(h_nonlocal)) h_nonlocal => potential%h_nonlocal
         IF (PRESENT(cprj_ppnl)) cprj_ppnl => potential%cprj_ppnl
         IF (PRESENT(vprj_ppnl)) vprj_ppnl => potential%vprj_ppnl

         IF (PRESENT(has_nlcc)) has_nlcc = potential%has_nlcc
         IF (PRESENT(n_nlcc)) n_nlcc = potential%n_nlcc
         IF (PRESENT(a_nlcc)) a_nlcc => potential%a_nlcc
         IF (PRESENT(c_nlcc)) c_nlcc => potential%c_nlcc

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE get_sgp_potential

! **************************************************************************************************
!> \brief   Initialise the coefficients of the projectors of the non-local
!>          part of the GTH pseudopotential and the transformation matrices
!>          for Cartesian overlap integrals between the orbital basis
!>          functions and the projector functions.
!> \param potential ...
!> \date    16.10.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_cprj_ppnl(potential)
      TYPE(gth_potential_type), POINTER                  :: potential

      INTEGER                                            :: cpx, cpy, cpz, cx, cy, cz, ico, iprj, &
                                                            iprj_ppnl, l, lp, lprj_ppnl, nprj, px, &
                                                            py, pz
      REAL(KIND=dp)                                      :: alpha_ppnl, cp

      nprj = 0

      DO l = 0, potential%lppnl
         alpha_ppnl = potential%alpha_ppnl(l)
         DO iprj_ppnl = 1, potential%nprj_ppnl(l)
            lp = iprj_ppnl - 1
            lprj_ppnl = l + 2*lp
            cp = SQRT(2.0_dp**(2.0_dp*REAL(lprj_ppnl, dp) + 3.5_dp)* &
                      alpha_ppnl**(REAL(lprj_ppnl, dp) + 1.5_dp)/ &
                      (rootpi*dfac(2*lprj_ppnl + 1)))
            potential%cprj_ppnl(iprj_ppnl, l) = cp
            DO cx = 0, l
               DO cy = 0, l - cx
                  cz = l - cx - cy
                  iprj = nprj + co(cx, cy, cz)
                  DO px = 0, lp
                     DO py = 0, lp - px
                        pz = lp - px - py
                        cpx = cx + 2*px
                        cpy = cy + 2*py
                        cpz = cz + 2*pz
                        ico = coset(cpx, cpy, cpz)
                        potential%cprj(ico, iprj) = cp*fac(lp)/(fac(px)*fac(py)*fac(pz))
                     END DO
                  END DO
               END DO
            END DO
            nprj = nprj + nco(l)
         END DO
      END DO

   END SUBROUTINE init_cprj_ppnl

! **************************************************************************************************
!> \brief   Initialise a GTH potential data set structure.
!> \param potential ...
!> \date    27.10.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_gth_potential(potential)
      TYPE(gth_potential_type), POINTER                  :: potential

      IF (.NOT. ASSOCIATED(potential)) RETURN

      IF (potential%nppnl > 0) THEN

         !     *** Initialise the projector coefficients of the    ***
         !     *** non-local part of the GTH pseudopotential and   ***
         !     *** the transformation matrices "pgf" -> "prj_ppnl" ***

         CALL init_cprj_ppnl(potential)

         !     *** Initialise the h(i,j) projector coefficients of ***
         !     *** the non-local part of the GTH pseudopotential   ***

         CALL init_vprj_ppnl(potential)

      END IF

   END SUBROUTINE init_gth_potential

! **************************************************************************************************
!> \brief   Initialise the h(i,j) projector coefficients of the non-local part
!>          of the GTH pseudopotential.
!> \param potential ...
!> \date    24.10.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_vprj_ppnl(potential)
      TYPE(gth_potential_type), POINTER                  :: potential

      INTEGER                                            :: i, ico, iprj, iprj_ppnl, iso, j, jco, &
                                                            jprj, jprj_ppnl, l, nprj

      nprj = 0

      DO l = 0, potential%lppnl
         DO iprj_ppnl = 1, potential%nprj_ppnl(l)
            iprj = nprj + (iprj_ppnl - 1)*nco(l)
            DO jprj_ppnl = 1, potential%nprj_ppnl(l)
               jprj = nprj + (jprj_ppnl - 1)*nco(l)
               DO ico = 1, nco(l)
                  i = iprj + ico
                  DO jco = 1, nco(l)
                     j = jprj + jco
                     DO iso = 1, nso(l)
                        potential%vprj_ppnl(i, j) = potential%vprj_ppnl(i, j) + &
                                                    orbtramat(l)%slm(iso, ico)* &
                                                    potential%hprj_ppnl(iprj_ppnl, &
                                                                        jprj_ppnl, l)* &
                                                    orbtramat(l)%slm(iso, jco)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         nprj = nprj + potential%nprj_ppnl(l)*nco(l)
      END DO

   END SUBROUTINE init_vprj_ppnl

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param itype ...
!> \param zeff ...
!> \param zeff_correction ...
! **************************************************************************************************
   SUBROUTINE init_all_potential(potential, itype, zeff, zeff_correction)

      TYPE(all_potential_type), POINTER                  :: potential
      CHARACTER(LEN=*), OPTIONAL                         :: itype
      REAL(KIND=dp), OPTIONAL                            :: zeff, zeff_correction

      INTEGER                                            :: dz

      IF (.NOT. ASSOCIATED(potential)) RETURN

      IF (PRESENT(zeff)) potential%zeff = zeff
      IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
      dz = potential%z - INT(potential%zeff - potential%zeff_correction)
      SELECT CASE (dz)
      CASE DEFAULT
      CASE (2)
         potential%elec_conf(0) = potential%elec_conf(0) - 2
      CASE (10)
         potential%elec_conf(0) = potential%elec_conf(0) - 4
         potential%elec_conf(1) = potential%elec_conf(1) - 6
      CASE (18)
         potential%elec_conf(0) = potential%elec_conf(0) - 6
         potential%elec_conf(1) = potential%elec_conf(1) - 12
      CASE (28)
         potential%elec_conf(0) = potential%elec_conf(0) - 6
         potential%elec_conf(1) = potential%elec_conf(1) - 12
         potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (30)
         potential%elec_conf(0) = potential%elec_conf(0) - 8
         potential%elec_conf(1) = potential%elec_conf(1) - 12
         potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (36)
         potential%elec_conf(0) = potential%elec_conf(0) - 8
         potential%elec_conf(1) = potential%elec_conf(1) - 18
         potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (46)
         potential%elec_conf(0) = potential%elec_conf(0) - 8
         potential%elec_conf(1) = potential%elec_conf(1) - 18
         potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (48)
         potential%elec_conf(0) = potential%elec_conf(0) - 10
         potential%elec_conf(1) = potential%elec_conf(1) - 18
         potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (54)
         potential%elec_conf(0) = potential%elec_conf(0) - 10
         potential%elec_conf(1) = potential%elec_conf(1) - 24
         potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (68)
         potential%elec_conf(0) = potential%elec_conf(0) - 10
         potential%elec_conf(1) = potential%elec_conf(1) - 24
         potential%elec_conf(2) = potential%elec_conf(2) - 20
         potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (78)
         potential%elec_conf(0) = potential%elec_conf(0) - 10
         potential%elec_conf(1) = potential%elec_conf(1) - 24
         potential%elec_conf(2) = potential%elec_conf(2) - 30
         potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (80)
         potential%elec_conf(0) = potential%elec_conf(0) - 12
         potential%elec_conf(1) = potential%elec_conf(1) - 24
         potential%elec_conf(2) = potential%elec_conf(2) - 30
         potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (86)
         potential%elec_conf(0) = potential%elec_conf(0) - 12
         potential%elec_conf(1) = potential%elec_conf(1) - 30
         potential%elec_conf(2) = potential%elec_conf(2) - 30
         potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (100)
         potential%elec_conf(0) = potential%elec_conf(0) - 12
         potential%elec_conf(1) = potential%elec_conf(1) - 30
         potential%elec_conf(2) = potential%elec_conf(2) - 30
         potential%elec_conf(3) = potential%elec_conf(3) - 28
      END SELECT

      IF (PRESENT(itype)) THEN
         IF (itype == "BARE") THEN
            potential%description(1) = "Bare Coulomb Potential"
            IF (dz > 0) THEN
               potential%description(2) = "Valence charge only"
            ELSE
               potential%description(2) = "Full atomic charge"
            END IF
         ENDIF
      END IF

   END SUBROUTINE init_all_potential
! **************************************************************************************************
!> \brief   Initialise a SGP potential data set structure.
!> \param potential ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_sgp_potential(potential)
      TYPE(sgp_potential_type), POINTER                  :: potential

      INTEGER                                            :: i1, i2, j1, j2, l, la, lb, n1, n2, nnl, &
                                                            nprj
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ind1, ind2
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, hnl
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cn

      IF (ASSOCIATED(potential)) THEN
         IF (potential%nppnl > 0) THEN
            !
            IF (ASSOCIATED(potential%cprj_ppnl)) THEN
               DEALLOCATE (potential%cprj_ppnl)
            END IF
            nprj = 0
            DO l = 0, potential%lmax
               nnl = potential%n_nonlocal
               IF (potential%is_nonlocal(l)) nprj = nprj + nnl*nso(l)
            END DO
            ALLOCATE (potential%cprj_ppnl(potential%nppnl, nprj))
            cprj => potential%cprj_ppnl
            cprj = 0.0_dp
            cn => potential%c_nonlocal
            nnl = potential%n_nonlocal
            !
            n1 = 0
            DO la = 0, potential%lmax
               n1 = n1 + nnl*nco(la)
            END DO
            ALLOCATE (ind1(n1, 3))
            n1 = 0
            DO i1 = 1, nnl
               DO la = 0, potential%lmax
                  DO j1 = 1, nco(la)
                     n1 = n1 + 1
                     ind1(n1, 1) = la
                     ind1(n1, 2) = j1
                     ind1(n1, 3) = i1
                  END DO
               END DO
            END DO
            !
            n2 = 0
            DO i2 = 1, nnl
               DO lb = 0, potential%lmax
                  IF (.NOT. potential%is_nonlocal(lb)) CYCLE
                  n2 = n2 + nso(lb)
               END DO
            END DO
            ALLOCATE (ind2(n2, 3))
            n2 = 0
            DO i2 = 1, nnl
               DO lb = 0, potential%lmax
                  IF (.NOT. potential%is_nonlocal(lb)) CYCLE
                  DO j2 = 1, nso(lb)
                     n2 = n2 + 1
                     ind2(n2, 1) = lb
                     ind2(n2, 2) = j2
                     ind2(n2, 3) = i2
                  END DO
               END DO
            END DO
            !
            DO n1 = 1, SIZE(ind1, 1)
               la = ind1(n1, 1)
               j1 = ind1(n1, 2)
               i1 = ind1(n1, 3)
               DO n2 = 1, SIZE(ind2, 1)
                  lb = ind2(n2, 1)
                  IF (la /= lb) CYCLE
                  j2 = ind2(n2, 2)
                  i2 = ind2(n2, 3)
                  cprj(n1, n2) = orbtramat(la)%c2s(j2, j1)*cn(i1, i2, la)
               END DO
            END DO
            !
            hnl => potential%h_nonlocal
            IF (ASSOCIATED(potential%vprj_ppnl)) THEN
               DEALLOCATE (potential%vprj_ppnl)
            END IF
            ALLOCATE (potential%vprj_ppnl(nprj))
            potential%vprj_ppnl = 0.0_dp
            DO n2 = 1, SIZE(ind2, 1)
               lb = ind2(n2, 1)
               i2 = ind2(n2, 3)
               potential%vprj_ppnl(n2) = hnl(i2, lb)
            END DO
            !
            DEALLOCATE (ind1, ind2)
         END IF
      END IF

   END SUBROUTINE init_sgp_potential

! **************************************************************************************************
!> \brief   Read an atomic all-electron potential data set.
!> \param element_symbol ...
!> \param potential_name ...
!> \param potential ...
!> \param zeff_correction ...
!> \param para_env ...
!> \param potential_file_name ...
!> \param potential_section ...
!> \param update_input ...
!> \date    14.05.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE read_all_potential(element_symbol, potential_name, potential, zeff_correction, &
                                 para_env, potential_file_name, potential_section, update_input)

      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol, potential_name
      TYPE(all_potential_type), POINTER                  :: potential
      REAL(KIND=dp)                                      :: zeff_correction
      TYPE(cp_para_env_type), POINTER                    :: para_env
      CHARACTER(len=default_path_length), INTENT(IN)     :: potential_file_name
      TYPE(section_vals_type), POINTER                   :: potential_section
      LOGICAL, INTENT(IN)                                :: update_input

      CHARACTER(LEN=240)                                 :: line
      CHARACTER(LEN=242)                                 :: line2
      CHARACTER(len=5*default_string_length)             :: line_att
      CHARACTER(LEN=LEN(element_symbol))                 :: symbol
      CHARACTER(LEN=LEN(element_symbol)+2)               :: symbol2
      CHARACTER(LEN=LEN(potential_name))                 :: apname
      CHARACTER(LEN=LEN(potential_name)+2)               :: apname2
      INTEGER                                            :: irep, l, strlen1, strlen2
      INTEGER, DIMENSION(:), POINTER                     :: elec_conf
      LOGICAL                                            :: found, is_ok, match, read_from_input
      REAL(KIND=dp)                                      :: alpha, r
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(cp_sll_val_type), POINTER                     :: list
      TYPE(val_type), POINTER                            :: val

      line2 = ""
      symbol2 = ""
      apname2 = ""
      NULLIFY (parser)
      CALL cite_reference(Krack2000)

      potential%name = potential_name
      read_from_input = .FALSE.
      CALL section_vals_get(potential_section, explicit=read_from_input)
      IF (.NOT. read_from_input) THEN
         CALL parser_create(parser, potential_file_name, para_env=para_env)
      END IF

      !   *** Search for the requested potential in the potential file   ***
      !   *** until the potential is found or the end of file is reached ***

      apname = potential_name
      symbol = element_symbol
      irep = 0
      search_loop: DO
         IF (read_from_input) THEN
            NULLIFY (list, val)
            found = .TRUE.
            CALL section_vals_list_get(potential_section, "_DEFAULT_KEYWORD_", list=list)
         ELSE
            CALL parser_search_string(parser, TRIM(apname), .TRUE., found, line)
         END IF
         IF (found) THEN
            CALL uppercase(symbol)
            CALL uppercase(apname)

            IF (read_from_input) THEN
               match = .TRUE.
            ELSE
               ! Check both the element symbol and the atomic potential name
               match = .FALSE.
               CALL uppercase(line)
               line2 = " "//line//" "
               symbol2 = " "//TRIM(symbol)//" "
               apname2 = " "//TRIM(apname)//" "
               strlen1 = LEN_TRIM(symbol2) + 1
               strlen2 = LEN_TRIM(apname2) + 1

               IF ((INDEX(line2, symbol2(:strlen1)) > 0) .AND. &
                   (INDEX(line2, apname2(:strlen2)) > 0)) match = .TRUE.
            END IF
            IF (match) THEN
               ! Read the electronic configuration
               NULLIFY (elec_conf)
               l = 0
               CALL reallocate(elec_conf, 0, l)
               IF (read_from_input) THEN
                  is_ok = cp_sll_val_next(list, val)
                  IF (.NOT. is_ok) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  CALL val_get(val, c_val=line_att)
                  READ (line_att, *) elec_conf(l)
                  CALL remove_word(line_att)
                  DO WHILE (LEN_TRIM(line_att) /= 0)
                     l = l + 1
                     CALL reallocate(elec_conf, 0, l)
                     READ (line_att, *) elec_conf(l)
                     CALL remove_word(line_att)
                  END DO
               ELSE
                  CALL parser_get_object(parser, elec_conf(l), newline=.TRUE.)
                  DO WHILE (parser_test_next_token(parser) == "INT")
                     l = l + 1
                     CALL reallocate(elec_conf, 0, l)
                     CALL parser_get_object(parser, elec_conf(l))
                  END DO
                  irep = irep + 1
                  IF (update_input) THEN
                     WRITE (line_att, '(100(1X,I0))') elec_conf
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
               END IF

               CALL reallocate(potential%elec_conf, 0, l)
               potential%elec_conf(:) = elec_conf(:)

               potential%zeff_correction = zeff_correction
               potential%zeff = REAL(SUM(elec_conf), dp) + zeff_correction

               DEALLOCATE (elec_conf)

               ! Read r(loc) to define the exponent of the core charge
               ! distribution and calculate the corresponding coefficient

               IF (read_from_input) THEN
                  is_ok = cp_sll_val_next(list, val)
                  IF (.NOT. is_ok) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  CALL val_get(val, c_val=line_att)
                  READ (line_att, *) r
               ELSE
                  CALL parser_get_object(parser, r, newline=.TRUE.)
                  irep = irep + 1
                  IF (update_input) THEN
                     WRITE (line_att, '(E24.16)') r
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
               END IF
               alpha = 1.0_dp/(2.0_dp*r**2)

               potential%alpha_core_charge = alpha
               potential%ccore_charge = potential%zeff*SQRT((alpha/pi)**3)

               EXIT search_loop
            END IF
         ELSE
            ! Stop program, if the end of file is reached
            CALL cp_abort(__LOCATION__, &
                          "The requested atomic potential <"// &
                          TRIM(potential_name)// &
                          "> for element <"// &
                          TRIM(symbol)// &
                          "> was not found in the potential file <"// &
                          TRIM(potential_file_name)//">")
         END IF
      END DO search_loop

      IF (.NOT. read_from_input) THEN
         ! Dump the potential info the in potential section
         IF (match .AND. update_input) THEN
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential read from the potential filename: "//TRIM(potential_file_name)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
         END IF
         CALL parser_release(parser)
      END IF
   END SUBROUTINE read_all_potential

! **************************************************************************************************
!> \brief   Read an atomic local potential data set.
!> \param element_symbol ...
!> \param potential_name ...
!> \param potential ...
!> \param para_env ...
!> \param potential_file_name ...
!> \param potential_section ...
!> \param update_input ...
!> \date    24.12.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE read_local_potential(element_symbol, potential_name, potential, &
                                   para_env, potential_file_name, potential_section, update_input)

      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol, potential_name
      TYPE(local_potential_type), POINTER                :: potential
      TYPE(cp_para_env_type), POINTER                    :: para_env
      CHARACTER(len=default_path_length), INTENT(IN)     :: potential_file_name
      TYPE(section_vals_type), POINTER                   :: potential_section
      LOGICAL, INTENT(IN)                                :: update_input

      REAL(KIND=dp), PARAMETER                           :: eps_tpot = 1.e-10_dp

      CHARACTER(LEN=240)                                 :: line
      CHARACTER(LEN=242)                                 :: line2
      CHARACTER(len=5*default_string_length)             :: line_att
      CHARACTER(LEN=LEN(element_symbol))                 :: symbol
      CHARACTER(LEN=LEN(element_symbol)+2)               :: symbol2
      CHARACTER(LEN=LEN(potential_name))                 :: apname
      CHARACTER(LEN=LEN(potential_name)+2)               :: apname2
      INTEGER                                            :: igau, ipol, irep, l, ngau, npol, &
                                                            strlen1, strlen2
      LOGICAL                                            :: found, is_ok, match, read_from_input
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(cp_sll_val_type), POINTER                     :: list
      TYPE(val_type), POINTER                            :: val

      line2 = ""
      symbol2 = ""
      apname2 = ""
      NULLIFY (parser, alpha, cval)

      potential%name = potential_name
      read_from_input = .FALSE.
      CALL section_vals_get(potential_section, explicit=read_from_input)
      IF (.NOT. read_from_input) THEN
         CALL parser_create(parser, potential_file_name, para_env=para_env)
      END IF

      !   *** Search for the requested potential in the potential file   ***
      !   *** until the potential is found or the end of file is reached ***

      apname = potential_name
      symbol = element_symbol
      irep = 0
      search_loop: DO
         IF (read_from_input) THEN
            NULLIFY (list, val)
            found = .TRUE.
            CALL section_vals_list_get(potential_section, "_DEFAULT_KEYWORD_", list=list)
         ELSE
            CALL parser_search_string(parser, TRIM(apname), .TRUE., found, line)
         END IF
         IF (found) THEN
            CALL uppercase(symbol)
            CALL uppercase(apname)

            IF (read_from_input) THEN
               match = .TRUE.
            ELSE
               ! Check both the element symbol and the atomic potential name
               match = .FALSE.
               CALL uppercase(line)
               line2 = " "//line//" "
               symbol2 = " "//TRIM(symbol)//" "
               apname2 = " "//TRIM(apname)//" "
               strlen1 = LEN_TRIM(symbol2) + 1
               strlen2 = LEN_TRIM(apname2) + 1

               IF ((INDEX(line2, symbol2(:strlen1)) > 0) .AND. &
                   (INDEX(line2, apname2(:strlen2)) > 0)) match = .TRUE.
            END IF
            IF (match) THEN

               ! Read ngau and npol
               IF (read_from_input) THEN
                  is_ok = cp_sll_val_next(list, val)
                  IF (.NOT. is_ok) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  CALL val_get(val, c_val=line_att)
                  READ (line_att, *) ngau, npol
                  CALL remove_word(line_att)
               ELSE
                  CALL parser_get_object(parser, ngau, newline=.TRUE.)
                  CALL parser_get_object(parser, npol)
                  irep = irep + 1
                  IF (update_input) THEN
                     WRITE (line_att, '(2(1X,I0))') ngau, npol
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
               END IF

               CALL reallocate(alpha, 1, ngau)
               CALL reallocate(cval, 1, ngau, 1, npol)
               DO igau = 1, ngau
                  IF (read_from_input) THEN
                     is_ok = cp_sll_val_next(list, val)
                     IF (.NOT. is_ok) &
                        CALL cp_abort(__LOCATION__, &
                                      "Error reading the Potential from input file!!")
                     CALL val_get(val, c_val=line_att)
                     READ (line_att, *) alpha(igau), (cval(igau, ipol), ipol=1, npol)
                  ELSE
                     CALL parser_get_object(parser, alpha(igau), newline=.TRUE.)
                     DO ipol = 1, npol
                        CALL parser_get_object(parser, cval(igau, ipol), newline=.FALSE.)
                     END DO
                     irep = irep + 1
                     IF (update_input) THEN
                        WRITE (line_att, '(8(E24.16))') alpha(igau), (cval(igau, ipol), ipol=1, npol)
                        CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                  c_val=TRIM(line_att))
                     END IF
                  END IF
               END DO
               alpha = 1.0_dp/(2.0_dp*alpha**2)

               potential%ngau = ngau
               potential%npol = npol

               potential%alpha => alpha
               potential%cval => cval

               potential%radius = 0.0_dp
               DO igau = 1, ngau
                  DO ipol = 1, npol
                     l = 2*(ipol - 1)
                     potential%radius = MAX(potential%radius, &
                                            exp_radius(l, alpha(igau), eps_tpot, cval(igau, ipol), &
                                                       rlow=potential%radius))
                  END DO
               END DO

               EXIT search_loop
            END IF
         ELSE
            ! Stop program, if the end of file is reached
            CALL cp_abort(__LOCATION__, &
                          "The requested local atomic potential <"// &
                          TRIM(potential_name)// &
                          "> for element <"// &
                          TRIM(symbol)// &
                          "> was not found in the potential file <"// &
                          TRIM(potential_file_name)//">")
         END IF
      END DO search_loop

      IF (.NOT. read_from_input) THEN
         ! Dump the potential info in the potential section
         IF (match .AND. update_input) THEN
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential read from the potential filename: "//TRIM(potential_file_name)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
         END IF
         CALL parser_release(parser)
      END IF

   END SUBROUTINE read_local_potential

! **************************************************************************************************
!> \brief   Read an atomic GTH potential data set.
!> \param element_symbol ...
!> \param potential_name ...
!> \param potential ...
!> \param zeff_correction ...
!> \param para_env ...
!> \param potential_file_name ...
!> \param potential_section ...
!> \param update_input ...
!> \date    14.05.2000
!> \par  Literature
!>         - S. Goedecker, M. Teter and J. Hutter,
!>                Phys. Rev. B 54, 1703 (1996)
!>         - C. Hartwigsen, S. Goedecker and J. Hutter,
!>                Phys. Rev. B 58, 3641 (1998)
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE read_gth_potential(element_symbol, potential_name, potential, zeff_correction, &
                                 para_env, potential_file_name, potential_section, update_input)

      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol, potential_name
      TYPE(gth_potential_type), POINTER                  :: potential
      REAL(KIND=dp), INTENT(IN)                          :: zeff_correction
      TYPE(cp_para_env_type), POINTER                    :: para_env
      CHARACTER(len=default_path_length), INTENT(IN)     :: potential_file_name
      TYPE(section_vals_type), POINTER                   :: potential_section
      LOGICAL, INTENT(IN)                                :: update_input

      CHARACTER(LEN=240)                                 :: line
      CHARACTER(LEN=242)                                 :: line2
      CHARACTER(len=5*default_string_length)             :: line_att
      CHARACTER(LEN=LEN(element_symbol))                 :: symbol
      CHARACTER(LEN=LEN(element_symbol)+2)               :: symbol2
      CHARACTER(LEN=LEN(potential_name))                 :: apname
      CHARACTER(LEN=LEN(potential_name)+2)               :: apname2
      INTEGER                                            :: i, ic, ipot, irep, istr, j, l, lppnl, &
                                                            lprj_ppnl_max, maxlppl, n, nppnl, &
                                                            nprj_ppnl, nprj_ppnl_max, strlen1, &
                                                            strlen2
      INTEGER, DIMENSION(:), POINTER                     :: elec_conf
      LOGICAL                                            :: found, is_ok, match, read_from_input
      REAL(KIND=dp)                                      :: alpha, ci, r, rc2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmp_vals
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: hprj_ppnl
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(cp_sll_val_type), POINTER                     :: list
      TYPE(val_type), POINTER                            :: val

      line2 = ""
      symbol2 = ""
      apname2 = ""
      NULLIFY (parser, tmp_vals)
      CALL cite_reference(Goedecker1996)
      CALL cite_reference(Hartwigsen1998)
      CALL cite_reference(Krack2005)

      potential%name = potential_name
      potential%aliases = potential_name
      read_from_input = .FALSE.
      CALL section_vals_get(potential_section, explicit=read_from_input)
      IF (.NOT. read_from_input) THEN
         CALL parser_create(parser, potential_file_name, para_env=para_env)
      END IF

      !initialize extended form
      potential%lpotextended = .FALSE.
      potential%nexp_lpot = 0
      potential%lsdpot = .FALSE.
      potential%nexp_lsd = 0
      potential%nlcc = .FALSE.
      potential%nexp_nlcc = 0

      !   *** Search for the requested potential in the potential file   ***
      !   *** until the potential is found or the end of file is reached ***

      apname = potential_name
      symbol = element_symbol
      irep = 0
      search_loop: DO
         IF (read_from_input) THEN
            NULLIFY (list, val)
            found = .TRUE.
            CALL section_vals_list_get(potential_section, "_DEFAULT_KEYWORD_", list=list)
         ELSE
            CALL parser_search_string(parser, TRIM(apname), .TRUE., found, line)
         END IF
         IF (found) THEN
            CALL uppercase(symbol)
            CALL uppercase(apname)
            IF (read_from_input) THEN
               match = .TRUE.
            ELSE
               !       *** Check both the element symbol and the atomic potential name ***
               match = .FALSE.
               CALL uppercase(line)
               line2 = " "//line//" "
               symbol2 = " "//TRIM(symbol)//" "
               apname2 = " "//TRIM(apname)//" "
               strlen1 = LEN_TRIM(symbol2) + 1
               strlen2 = LEN_TRIM(apname2) + 1

               i = INDEX(line2, symbol2(:strlen1))
               j = INDEX(line2, apname2(:strlen2))
               IF (i > 0 .AND. j > 0) THEN
                  match = .TRUE.
                  i = i + 1 + INDEX(line2(i + 1:), " ")
                  potential%aliases = line2(i:) ! copy all names into aliases field
               ENDIF
            END IF
            IF (match) THEN
               !         *** Read the electronic configuration ***
               NULLIFY (elec_conf)
               l = 0
               CALL reallocate(elec_conf, 0, l)
               IF (read_from_input) THEN
                  is_ok = cp_sll_val_next(list, val)
                  IF (.NOT. is_ok) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  CALL val_get(val, c_val=line_att)
                  READ (line_att, *) elec_conf(l)
                  CALL remove_word(line_att)
                  DO WHILE (LEN_TRIM(line_att) /= 0)
                     l = l + 1
                     CALL reallocate(elec_conf, 0, l)
                     READ (line_att, *) elec_conf(l)
                     CALL remove_word(line_att)
                  END DO
               ELSE
                  CALL parser_get_object(parser, elec_conf(l), newline=.TRUE.)
                  DO WHILE (parser_test_next_token(parser) == "INT")
                     l = l + 1
                     CALL reallocate(elec_conf, 0, l)
                     CALL parser_get_object(parser, elec_conf(l))
                  END DO
                  irep = irep + 1
                  IF (update_input) THEN
                     WRITE (line_att, '(100(1X,I0))') elec_conf
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
               END IF

               CALL reallocate(potential%elec_conf, 0, l)
               potential%elec_conf(:) = elec_conf(:)

               potential%zeff_correction = zeff_correction
               potential%zeff = REAL(SUM(elec_conf), dp) + zeff_correction

               DEALLOCATE (elec_conf)

               !         *** Read r(loc) to define the exponent of the core charge    ***
               !         *** distribution and calculate the corresponding coefficient ***

               IF (read_from_input) THEN
                  is_ok = cp_sll_val_next(list, val)
                  IF (.NOT. is_ok) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  CALL val_get(val, c_val=line_att)
                  READ (line_att, *) r
                  CALL remove_word(line_att)
               ELSE
                  line_att = ""
                  CALL parser_get_object(parser, r, newline=.TRUE.)
                  istr = LEN_TRIM(line_att) + 1
                  WRITE (line_att(istr:), '(E24.16)') r
               END IF
               alpha = 1.0_dp/(2.0_dp*r**2)

               potential%alpha_core_charge = alpha
               potential%ccore_charge = potential%zeff*SQRT((alpha/pi)**3)

               potential%alpha_ppl = alpha
               potential%cerf_ppl = potential%zeff*SQRT((alpha/pi)**3)

               !         *** Read the parameters for the local part ***
               !         *** of the GTH pseudopotential (ppl)       ***

               IF (read_from_input) THEN
                  READ (line_att, *) n
                  CALL remove_word(line_att)
               ELSE
                  CALL parser_get_object(parser, n)
                  istr = LEN_TRIM(line_att) + 1
                  WRITE (line_att(istr:), '(1X,I0)') n
               END IF
               potential%nexp_ppl = n
               CALL reallocate(potential%cexp_ppl, 1, n)

               DO i = 1, n
                  IF (read_from_input) THEN
                     READ (line_att, *) ci
                     CALL remove_word(line_att)
                  ELSE
                     CALL parser_get_object(parser, ci)
                     istr = LEN_TRIM(line_att) + 1
                     WRITE (line_att(istr:), '(E24.16)') ci
                  END IF
                  rc2 = (2.0_dp*potential%alpha_ppl)
                  potential%cexp_ppl(i) = rc2**(i - 1)*ci
               END DO

               IF (.NOT. read_from_input) THEN
                  irep = irep + 1
                  IF (update_input) THEN
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
                  line_att = ""
               ELSE
                  IF (LEN_TRIM(line_att) /= 0) THEN
                     CALL cp_abort(__LOCATION__, &
                                   "Error reading the Potential from input file!!")
                  END IF
               END IF
               maxlppl = 2*(n - 1)

               IF (maxlppl > -1) CALL init_orbital_pointers(maxlppl)

               ! Read extended form of GTH pseudopotential
               ! local potential, NLCC, LSD potential
               IF (read_from_input) THEN
                  DO
                     is_ok = cp_sll_val_next(list, val)
                     CPASSERT(is_ok)
                     CALL val_get(val, c_val=line_att)
                     IF (INDEX(line_att, "LPOT") /= 0) THEN
                        potential%lpotextended = .TRUE.
                        CALL remove_word(line_att)
                        READ (line_att, *) potential%nexp_lpot
                        n = potential%nexp_lpot
                        maxlppl = 2*(n - 1)
                        IF (maxlppl > -1) CALL init_orbital_pointers(maxlppl)
                        NULLIFY (potential%alpha_lpot, potential%nct_lpot, potential%cval_lpot)
                        CALL reallocate(potential%alpha_lpot, 1, n)
                        CALL reallocate(potential%nct_lpot, 1, n)
                        CALL reallocate(potential%cval_lpot, 1, 4, 1, n)
                        DO ipot = 1, potential%nexp_lpot
                           is_ok = cp_sll_val_next(list, val)
                           CPASSERT(is_ok)
                           CALL val_get(val, c_val=line_att)
                           READ (line_att, *) r
                           potential%alpha_lpot(ipot) = 0.5_dp/(r*r)
                           CALL remove_word(line_att)
                           READ (line_att, *) potential%nct_lpot(ipot)
                           CALL remove_word(line_att)
                           DO ic = 1, potential%nct_lpot(ipot)
                              READ (line_att, *) ci
                              rc2 = (2._dp*potential%alpha_lpot(ipot))**(ic - 1)
                              potential%cval_lpot(ic, ipot) = ci*rc2
                              CALL remove_word(line_att)
                           END DO
                        END DO
                     ELSEIF (INDEX(line_att, "NLCC") /= 0) THEN
                        potential%nlcc = .TRUE.
                        CALL remove_word(line_att)
                        READ (line_att, *) potential%nexp_nlcc
                        n = potential%nexp_nlcc
                        NULLIFY (potential%alpha_nlcc, potential%nct_nlcc, potential%cval_nlcc)
                        CALL reallocate(potential%alpha_nlcc, 1, n)
                        CALL reallocate(potential%nct_nlcc, 1, n)
                        CALL reallocate(potential%cval_nlcc, 1, 4, 1, n)
                        DO ipot = 1, potential%nexp_nlcc
                           is_ok = cp_sll_val_next(list, val)
                           CPASSERT(is_ok)
                           CALL val_get(val, c_val=line_att)
                           READ (line_att, *) potential%alpha_nlcc(ipot)
                           CALL remove_word(line_att)
                           READ (line_att, *) potential%nct_nlcc(ipot)
                           CALL remove_word(line_att)
                           DO ic = 1, potential%nct_nlcc(ipot)
                              READ (line_att, *) potential%cval_nlcc(ic, ipot)
                              !make it compatible with bigdft style
                              potential%cval_nlcc(ic, ipot) = potential%cval_nlcc(ic, ipot)/(4.0_dp*pi)
                              CALL remove_word(line_att)
                           END DO
                        END DO
                     ELSEIF (INDEX(line_att, "LSD") /= 0) THEN
                        potential%lsdpot = .TRUE.
                        CALL remove_word(line_att)
                        READ (line_att, *) potential%nexp_lsd
                        n = potential%nexp_lsd
                        NULLIFY (potential%alpha_lsd, potential%nct_lsd, potential%cval_lsd)
                        CALL reallocate(potential%alpha_lsd, 1, n)
                        CALL reallocate(potential%nct_lsd, 1, n)
                        CALL reallocate(potential%cval_lsd, 1, 4, 1, n)
                        DO ipot = 1, potential%nexp_lsd
                           is_ok = cp_sll_val_next(list, val)
                           CPASSERT(is_ok)
                           CALL val_get(val, c_val=line_att)
                           READ (line_att, *) r
                           potential%alpha_lsd(ipot) = 0.5_dp/(r*r)
                           CALL remove_word(line_att)
                           READ (line_att, *) potential%nct_lsd(ipot)
                           CALL remove_word(line_att)
                           DO ic = 1, potential%nct_lsd(ipot)
                              READ (line_att, *) ci
                              rc2 = (2._dp*potential%alpha_lsd(ipot))**(ic - 1)
                              potential%cval_lsd(ic, ipot) = ci*rc2
                              CALL remove_word(line_att)
                           END DO
                        END DO
                     ELSE
                        EXIT
                     END IF
                  END DO
               ELSE
                  DO
                     CALL parser_get_next_line(parser, 1)
                     IF (parser_test_next_token(parser) == "INT") THEN
                        EXIT
                     ELSEIF (parser_test_next_token(parser) == "STR") THEN
                        CALL parser_get_object(parser, line)
                        IF (INDEX(LINE, "LPOT") /= 0) THEN
                           ! local potential
                           potential%lpotextended = .TRUE.
                           CALL parser_get_object(parser, potential%nexp_lpot)
                           n = potential%nexp_lpot
                           NULLIFY (potential%alpha_lpot, potential%nct_lpot, potential%cval_lpot)
                           CALL reallocate(potential%alpha_lpot, 1, n)
                           CALL reallocate(potential%nct_lpot, 1, n)
                           CALL reallocate(potential%cval_lpot, 1, 4, 1, n)
                           ! add to input section
                           irep = irep + 1
                           IF (update_input) THEN
                              WRITE (line_att, '(A,1X,I0)') "LPOT", n
                              CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                        c_val=TRIM(line_att))
                           END IF
                           DO ipot = 1, potential%nexp_lpot
                              CALL parser_get_object(parser, r, newline=.TRUE.)
                              potential%alpha_lpot(ipot) = 0.5_dp/(r*r)
                              CALL parser_get_object(parser, potential%nct_lpot(ipot))
                              CALL reallocate(tmp_vals, 1, potential%nct_lpot(ipot))
                              DO ic = 1, potential%nct_lpot(ipot)
                                 CALL parser_get_object(parser, ci)
                                 tmp_vals(ic) = ci
                                 rc2 = (2._dp*potential%alpha_lpot(ipot))**(ic - 1)
                                 potential%cval_lpot(ic, ipot) = ci*rc2
                              END DO
                              ! add to input section
                              irep = irep + 1
                              IF (update_input) THEN
                                 WRITE (line_att, '(E24.16,1X,I0,100(1X,E24.16))') r, potential%nct_lpot(ipot), &
                                    tmp_vals(1:potential%nct_lpot(ipot))
                                 CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                           c_val=TRIM(line_att))
                              END IF
                           END DO
                        ELSEIF (INDEX(LINE, "NLCC") /= 0) THEN
                           ! NLCC
                           potential%nlcc = .TRUE.
                           CALL parser_get_object(parser, potential%nexp_nlcc)
                           n = potential%nexp_nlcc
                           NULLIFY (potential%alpha_nlcc, potential%nct_nlcc, potential%cval_nlcc)
                           CALL reallocate(potential%alpha_nlcc, 1, n)
                           CALL reallocate(potential%nct_nlcc, 1, n)
                           CALL reallocate(potential%cval_nlcc, 1, 4, 1, n)
                           ! add to input section
                           WRITE (line_att, '(A,1X,I0)') "NLCC", n
                           irep = irep + 1
                           CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                     c_val=TRIM(line_att))
                           DO ipot = 1, potential%nexp_nlcc
                              CALL parser_get_object(parser, potential%alpha_nlcc(ipot), newline=.TRUE.)
                              CALL parser_get_object(parser, potential%nct_nlcc(ipot))
                              CALL reallocate(tmp_vals, 1, potential%nct_nlcc(ipot))
                              DO ic = 1, potential%nct_nlcc(ipot)
                                 CALL parser_get_object(parser, potential%cval_nlcc(ic, ipot))
                                 tmp_vals(ic) = potential%cval_nlcc(ic, ipot)
                                 !make it compatible with bigdft style
                                 potential%cval_nlcc(ic, ipot) = potential%cval_nlcc(ic, ipot)/(4.0_dp*pi)
                              END DO
                              ! add to input section
                              irep = irep + 1
                              IF (update_input) THEN
                                 WRITE (line_att, '(E24.16,1X,I0,100(1X,E24.16))') &
                                    potential%alpha_nlcc(ipot), potential%nct_nlcc(ipot), &
                                    tmp_vals(1:potential%nct_nlcc(ipot))
                                 CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                           c_val=TRIM(line_att))
                              END IF
                           END DO
                        ELSEIF (INDEX(LINE, "LSD") /= 0) THEN
                           ! LSD potential
                           potential%lsdpot = .TRUE.
                           CALL parser_get_object(parser, potential%nexp_lsd)
                           n = potential%nexp_lsd
                           NULLIFY (potential%alpha_lsd, potential%nct_lsd, potential%cval_lsd)
                           CALL reallocate(potential%alpha_lsd, 1, n)
                           CALL reallocate(potential%nct_lsd, 1, n)
                           CALL reallocate(potential%cval_lsd, 1, 4, 1, n)
                           ! add to input section
                           irep = irep + 1
                           IF (update_input) THEN
                              WRITE (line_att, '(A,1X,I0)') "LSD", n
                              CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                        c_val=TRIM(line_att))
                           END IF
                           DO ipot = 1, potential%nexp_lsd
                              CALL parser_get_object(parser, r, newline=.TRUE.)
                              potential%alpha_lsd(ipot) = 0.5_dp/(r*r)
                              CALL parser_get_object(parser, potential%nct_lsd(ipot))
                              CALL reallocate(tmp_vals, 1, potential%nct_lsd(ipot))
                              DO ic = 1, potential%nct_lsd(ipot)
                                 CALL parser_get_object(parser, ci)
                                 tmp_vals(ic) = ci
                                 rc2 = (2._dp*potential%alpha_lsd(ipot))**(ic - 1)
                                 potential%cval_lsd(ic, ipot) = ci*rc2
                              END DO
                              ! add to input section
                              irep = irep + 1
                              IF (update_input) THEN
                                 WRITE (line_att, '(E24.16,1X,I0,100(1X,E24.16))') r, potential%nct_lsd(ipot), &
                                    tmp_vals(1:potential%nct_lsd(ipot))
                                 CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                           c_val=TRIM(line_att))
                              END IF
                           END DO
                        ELSE
                           CALL cp_abort(__LOCATION__, &
                                         "Syntax error for <"// &
                                         TRIM(element_symbol)// &
                                         "> in the atomic potential <"// &
                                         TRIM(potential_name)// &
                                         "> potential file <"// &
                                         TRIM(potential_file_name)//">: "// &
                                         "Expected LPOT/NLCC/LSD keyword, got: <"// &
                                         TRIM(line)//">")
                        END IF
                     ELSE
                        CALL parser_get_object(parser, line)
                        CALL cp_abort(__LOCATION__, &
                                      "Syntax error for <"// &
                                      TRIM(element_symbol)// &
                                      "> in the atomic potential <"// &
                                      TRIM(potential_name)// &
                                      "> potential file <"// &
                                      TRIM(potential_file_name)//">: "// &
                                      "Expected LPOT/NLCC/LSD keyword or INTEGER, got: <"// &
                                      TRIM(line)//">")
                     END IF
                  END DO
               END IF
               !         *** Read the parameters for the non-local  ***
               !         *** part of the GTH pseudopotential (ppnl) ***
               IF (read_from_input) THEN
                  READ (line_att, *) n
                  CALL remove_word(line_att)
               ELSE
                  CALL parser_get_object(parser, n)
                  irep = irep + 1
                  IF (update_input) THEN
                     WRITE (line_att, '(1X,I0)') n
                     CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                               c_val=TRIM(line_att))
                  END IF
               END IF
               potential%lppnl = n - 1
               potential%nppnl = 0

               potential%lprj_ppnl_max = n - 1
               potential%nprj_ppnl_max = 0

               IF (n > 0) THEN

                  lppnl = potential%lppnl
                  nppnl = potential%nppnl

                  CALL init_orbital_pointers(lppnl)

                  NULLIFY (hprj_ppnl)

                  ! Load the parameter for n non-local projectors

                  CALL reallocate(potential%alpha_ppnl, 0, lppnl)
                  CALL reallocate(potential%nprj_ppnl, 0, lppnl)

                  lprj_ppnl_max = -1
                  nprj_ppnl_max = 0

                  DO l = 0, lppnl
                     IF (read_from_input) THEN
                        is_ok = cp_sll_val_next(list, val)
                        IF (.NOT. is_ok) &
                           CALL cp_abort(__LOCATION__, &
                                         "Error reading the Potential from input file!!")
                        CALL val_get(val, c_val=line_att)
                        READ (line_att, *) r
                        CALL remove_word(line_att)
                        READ (line_att, *) nprj_ppnl
                        CALL remove_word(line_att)
                     ELSE
                        line_att = ""
                        CALL parser_get_object(parser, r, newline=.TRUE.)
                        CALL parser_get_object(parser, nprj_ppnl)
                        istr = LEN_TRIM(line_att) + 1
                        WRITE (line_att(istr:), '(E24.16,1X,I0)') r, nprj_ppnl
                     END IF
                     IF (r == 0.0_dp .AND. nprj_ppnl /= 0) THEN
                        CALL cp_abort(__LOCATION__, &
                                      "An error was detected in the atomic potential <"// &
                                      TRIM(potential_name)// &
                                      "> potential file <"// &
                                      TRIM(potential_file_name)//">")
                     END IF
                     potential%alpha_ppnl(l) = 0.0_dp
                     IF (r /= 0.0_dp .AND. n /= 0) potential%alpha_ppnl(l) = 1.0_dp/(2.0_dp*r**2)
                     potential%nprj_ppnl(l) = nprj_ppnl
                     nppnl = nppnl + nprj_ppnl*nco(l)
                     IF (nprj_ppnl > nprj_ppnl_max) THEN
                        nprj_ppnl_max = nprj_ppnl
                        CALL reallocate(hprj_ppnl, 1, nprj_ppnl_max, &
                                        1, nprj_ppnl_max, &
                                        0, lppnl)
                     END IF
                     DO i = 1, nprj_ppnl
                        IF (i == 1) THEN
                           IF (read_from_input) THEN
                              READ (line_att, *) hprj_ppnl(i, i, l)
                              CALL remove_word(line_att)
                           ELSE
                              CALL parser_get_object(parser, hprj_ppnl(i, i, l))
                              istr = LEN_TRIM(line_att) + 1
                              WRITE (line_att(istr:), '(E24.16)') hprj_ppnl(i, i, l)
                           END IF
                        ELSE
                           IF (read_from_input) THEN
                              IF (LEN_TRIM(line_att) /= 0) &
                                 CALL cp_abort(__LOCATION__, &
                                               "Error reading the Potential from input file!!")
                              is_ok = cp_sll_val_next(list, val)
                              IF (.NOT. is_ok) &
                                 CALL cp_abort(__LOCATION__, &
                                               "Error reading the Potential from input file!!")
                              CALL val_get(val, c_val=line_att)
                              READ (line_att, *) hprj_ppnl(i, i, l)
                              CALL remove_word(line_att)
                           ELSE
                              irep = irep + 1
                              CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                        c_val=TRIM(line_att))
                              line_att = ""
                              CALL parser_get_object(parser, hprj_ppnl(i, i, l), newline=.TRUE.)
                              istr = LEN_TRIM(line_att) + 1
                              WRITE (line_att(istr:), '(E24.16)') hprj_ppnl(i, i, l)
                           END IF
                        END IF
                        DO j = i + 1, nprj_ppnl
                           IF (read_from_input) THEN
                              READ (line_att, *) hprj_ppnl(i, j, l)
                              CALL remove_word(line_att)
                           ELSE
                              CALL parser_get_object(parser, hprj_ppnl(i, j, l))
                              istr = LEN_TRIM(line_att) + 1
                              WRITE (line_att(istr:), '(E24.16)') hprj_ppnl(i, j, l)
                           END IF
                        END DO
                     END DO
                     IF (.NOT. read_from_input) THEN
                        irep = irep + 1
                        IF (update_input) THEN
                           CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                                     c_val=TRIM(line_att))
                        END IF
                        line_att = ""
                     ELSE
                        IF (LEN_TRIM(line_att) /= 0) THEN
                           CALL cp_abort(__LOCATION__, &
                                         "Error reading the Potential from input file!!")
                        END IF
                     END IF
                     IF (nprj_ppnl > 1) THEN
                        CALL symmetrize_matrix(hprj_ppnl(:, :, l), "upper_to_lower")
                     END IF
                     lprj_ppnl_max = MAX(lprj_ppnl_max, l + 2*(nprj_ppnl - 1))
                  END DO

                  potential%nppnl = nppnl
                  CALL init_orbital_pointers(lprj_ppnl_max)

                  potential%lprj_ppnl_max = lprj_ppnl_max
                  potential%nprj_ppnl_max = nprj_ppnl_max
                  CALL reallocate(potential%hprj_ppnl, 1, nprj_ppnl_max, &
                                  1, nprj_ppnl_max, &
                                  0, lppnl)
                  potential%hprj_ppnl(:, :, :) = hprj_ppnl(:, :, :)

                  CALL reallocate(potential%cprj, 1, ncoset(lprj_ppnl_max), 1, nppnl)
                  CALL reallocate(potential%cprj_ppnl, 1, nprj_ppnl_max, 0, lppnl)
                  CALL reallocate(potential%vprj_ppnl, 1, nppnl, 1, nppnl)

                  DEALLOCATE (hprj_ppnl)
               END IF
               EXIT search_loop
            END IF
         ELSE
            ! Stop program, if the end of file is reached
            CALL cp_abort(__LOCATION__, &
                          "The requested atomic potential <"// &
                          TRIM(potential_name)// &
                          "> for element <"// &
                          TRIM(symbol)// &
                          "> was not found in the potential file <"// &
                          TRIM(potential_file_name)//">")
         END IF
      END DO search_loop

      IF (.NOT. read_from_input) THEN
         ! Dump the potential info the in potential section
         IF (match .AND. update_input) THEN
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential name: "//apname2(:strlen2)//" for symbol: "//symbol2(:strlen1)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
            irep = irep + 1
            WRITE (line_att, '(A)') "         # Potential read from the potential filename: "//TRIM(potential_file_name)
            CALL section_vals_val_set(potential_section, "_DEFAULT_KEYWORD_", i_rep_val=irep, &
                                      c_val=TRIM(line_att))
         END IF
         CALL parser_release(parser)
      END IF

      IF (ASSOCIATED(tmp_vals)) DEALLOCATE (tmp_vals)

   END SUBROUTINE read_gth_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param z ...
!> \param zeff_correction ...
! **************************************************************************************************
   SUBROUTINE set_default_all_potential(potential, z, zeff_correction)

      TYPE(all_potential_type), POINTER                  :: potential
      INTEGER, INTENT(IN)                                :: z
      REAL(KIND=dp)                                      :: zeff_correction

      CHARACTER(LEN=default_string_length)               :: name
      INTEGER, DIMENSION(:), POINTER                     :: elec_conf
      REAL(KIND=dp)                                      :: alpha, alpha_core_charge, ccore_charge, &
                                                            core_charge_radius, r, zeff

      ALLOCATE (elec_conf(0:3))
      elec_conf(0:3) = ptable(z)%e_conv(0:3)
      zeff = REAL(SUM(elec_conf), dp) + zeff_correction
      name = ptable(z)%name

      r = ptable(z)%covalent_radius*0.5_dp
      r = MAX(r, 0.2_dp)
      r = MIN(r, 1.0_dp)
      alpha = 1.0_dp/(2.0_dp*r**2)

      core_charge_radius = r
      alpha_core_charge = alpha
      ccore_charge = zeff*SQRT((alpha/pi)**3)

      CALL set_all_potential(potential, &
                             name=name, &
                             alpha_core_charge=alpha_core_charge, &
                             ccore_charge=ccore_charge, &
                             core_charge_radius=core_charge_radius, &
                             z=z, &
                             zeff=zeff, &
                             zeff_correction=zeff_correction, &
                             elec_conf=elec_conf)

      DEALLOCATE (elec_conf)

   END SUBROUTINE set_default_all_potential

! **************************************************************************************************
!> \brief   Set the attributes of an all-electron potential data set.
!> \param potential ...
!> \param name ...
!> \param alpha_core_charge ...
!> \param ccore_charge ...
!> \param core_charge_radius ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param elec_conf ...
!> \date    11.01.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE set_all_potential(potential, name, alpha_core_charge, &
                                ccore_charge, core_charge_radius, z, zeff, &
                                zeff_correction, elec_conf)

      TYPE(all_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: name
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha_core_charge, ccore_charge, &
                                                            core_charge_radius
      INTEGER, INTENT(IN), OPTIONAL                      :: z
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: zeff, zeff_correction
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) potential%name = name
         IF (PRESENT(alpha_core_charge)) &
            potential%alpha_core_charge = alpha_core_charge
         IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
         IF (PRESENT(core_charge_radius)) &
            potential%core_charge_radius = core_charge_radius
         IF (PRESENT(z)) potential%z = z
         IF (PRESENT(zeff)) potential%zeff = zeff
         IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
         IF (PRESENT(elec_conf)) THEN
            IF (.NOT. ASSOCIATED(potential%elec_conf)) THEN
               CALL reallocate(potential%elec_conf, 0, SIZE(elec_conf) - 1)
            END IF
            potential%elec_conf(:) = elec_conf(:)
         END IF

      ELSE

         CPABORT("The pointer potential is not associated")

      END IF

   END SUBROUTINE set_all_potential

! **************************************************************************************************
!> \brief   Set the attributes of an atomic local potential data set.
!> \param potential ...
!> \param name ...
!> \param alpha ...
!> \param cval ...
!> \param radius ...
!> \date    24.01.2014
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE set_local_potential(potential, name, alpha, cval, radius)

      TYPE(local_potential_type), POINTER                :: potential
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: name
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cval
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: radius

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) potential%name = name
         IF (PRESENT(alpha)) potential%alpha => alpha
         IF (PRESENT(cval)) potential%cval => cval
         IF (PRESENT(radius)) potential%radius = radius

      ELSE

         CPABORT("The pointer potential is not associated")

      END IF

   END SUBROUTINE set_local_potential

! **************************************************************************************************
!> \brief   Set the attributes of an effective charge and inducible point
!>          dipole potential data set.
!> \param potential ...
!> \param apol ...
!> \param cpol ...
!> \param qeff ...
!> \param mm_radius ...
!> \param qmmm_corr_radius ...
!> \param qmmm_radius ...
!> \date    05.03.2010
!> \author  Toon.Verstraelen@gmail.com
! **************************************************************************************************
   SUBROUTINE set_fist_potential(potential, apol, cpol, qeff, mm_radius, &
                                 qmmm_corr_radius, qmmm_radius)
      TYPE(fist_potential_type), POINTER                 :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: apol, cpol, qeff, mm_radius, &
                                                            qmmm_corr_radius, qmmm_radius

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(apol)) potential%apol = apol
         IF (PRESENT(cpol)) potential%cpol = cpol
         IF (PRESENT(mm_radius)) potential%mm_radius = mm_radius
         IF (PRESENT(qeff)) potential%qeff = qeff
         IF (PRESENT(qmmm_corr_radius)) potential%qmmm_corr_radius = qmmm_corr_radius
         IF (PRESENT(qmmm_radius)) potential%qmmm_radius = qmmm_radius

      ELSE

         CPABORT("The pointer potential is not associated")

      END IF

   END SUBROUTINE set_fist_potential

! **************************************************************************************************
!> \brief   Set the attributes of a GTH potential data set.
!> \param potential ...
!> \param name ...
!> \param alpha_core_charge ...
!> \param alpha_ppl ...
!> \param ccore_charge ...
!> \param cerf_ppl ...
!> \param core_charge_radius ...
!> \param ppl_radius ...
!> \param ppnl_radius ...
!> \param lppnl ...
!> \param lprj_ppnl_max ...
!> \param nexp_ppl ...
!> \param nppnl ...
!> \param nprj_ppnl_max ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param alpha_ppnl ...
!> \param cexp_ppl ...
!> \param elec_conf ...
!> \param nprj_ppnl ...
!> \param cprj ...
!> \param cprj_ppnl ...
!> \param vprj_ppnl ...
!> \param hprj_ppnl ...
!> \date    11.01.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE set_gth_potential(potential, name, alpha_core_charge, alpha_ppl, &
                                ccore_charge, cerf_ppl, core_charge_radius, &
                                ppl_radius, ppnl_radius, lppnl, lprj_ppnl_max, &
                                nexp_ppl, nppnl, nprj_ppnl_max, z, zeff, zeff_correction, &
                                alpha_ppnl, cexp_ppl, elec_conf, nprj_ppnl, cprj, cprj_ppnl, &
                                vprj_ppnl, hprj_ppnl)
      TYPE(gth_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: name
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha_core_charge, alpha_ppl, &
                                                            ccore_charge, cerf_ppl, &
                                                            core_charge_radius, ppl_radius, &
                                                            ppnl_radius
      INTEGER, INTENT(IN), OPTIONAL                      :: lppnl, lprj_ppnl_max, nexp_ppl, nppnl, &
                                                            nprj_ppnl_max, z
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: zeff, zeff_correction
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: alpha_ppnl, cexp_ppl
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf, nprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: cprj, cprj_ppnl, vprj_ppnl
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: hprj_ppnl

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) potential%name = name
         IF (PRESENT(alpha_core_charge)) &
            potential%alpha_core_charge = alpha_core_charge
         IF (PRESENT(alpha_ppl)) potential%alpha_ppl = alpha_ppl
         IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
         IF (PRESENT(cerf_ppl)) potential%cerf_ppl = cerf_ppl
         IF (PRESENT(core_charge_radius)) &
            potential%core_charge_radius = core_charge_radius
         IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius
         IF (PRESENT(ppnl_radius)) potential%ppnl_radius = ppnl_radius
         IF (PRESENT(lppnl)) potential%lppnl = lppnl
         IF (PRESENT(lprj_ppnl_max)) potential%lprj_ppnl_max = lprj_ppnl_max
         IF (PRESENT(nexp_ppl)) potential%nexp_ppl = nexp_ppl
         IF (PRESENT(nppnl)) potential%nppnl = nppnl
         IF (PRESENT(nprj_ppnl_max)) potential%nprj_ppnl_max = nprj_ppnl_max
         IF (PRESENT(z)) potential%z = z
         IF (PRESENT(zeff)) potential%zeff = zeff
         IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
         IF (PRESENT(alpha_ppnl)) potential%alpha_ppnl => alpha_ppnl
         IF (PRESENT(cexp_ppl)) potential%cexp_ppl => cexp_ppl
         IF (PRESENT(elec_conf)) THEN
            IF (ASSOCIATED(potential%elec_conf)) THEN
               DEALLOCATE (potential%elec_conf)
            ENDIF
            ALLOCATE (potential%elec_conf(0:SIZE(elec_conf) - 1))
            potential%elec_conf(:) = elec_conf(:)
         ENDIF
         IF (PRESENT(nprj_ppnl)) potential%nprj_ppnl => nprj_ppnl
         IF (PRESENT(cprj)) potential%cprj => cprj
         IF (PRESENT(cprj_ppnl)) potential%cprj_ppnl => cprj_ppnl
         IF (PRESENT(vprj_ppnl)) potential%vprj_ppnl => vprj_ppnl
         IF (PRESENT(hprj_ppnl)) potential%hprj_ppnl => hprj_ppnl

      ELSE

         CPABORT("The pointer potential is not associated")

      END IF

   END SUBROUTINE set_gth_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param name ...
!> \param description ...
!> \param aliases ...
!> \param elec_conf ...
!> \param z ...
!> \param zeff ...
!> \param zeff_correction ...
!> \param alpha_core_charge ...
!> \param ccore_charge ...
!> \param core_charge_radius ...
!> \param ppl_radius ...
!> \param ppnl_radius ...
!> \param ecp_local ...
!> \param n_local ...
!> \param a_local ...
!> \param c_local ...
!> \param nloc ...
!> \param nrloc ...
!> \param aloc ...
!> \param bloc ...
!> \param n_nonlocal ...
!> \param nppnl ...
!> \param lmax ...
!> \param is_nonlocal ...
!> \param a_nonlocal ...
!> \param h_nonlocal ...
!> \param c_nonlocal ...
!> \param has_nlcc ...
!> \param n_nlcc ...
!> \param a_nlcc ...
!> \param c_nlcc ...
! **************************************************************************************************
   SUBROUTINE set_sgp_potential(potential, name, description, aliases, elec_conf, &
                                z, zeff, zeff_correction, alpha_core_charge, &
                                ccore_charge, core_charge_radius, &
                                ppl_radius, ppnl_radius, &
                                ecp_local, n_local, a_local, c_local, &
                                nloc, nrloc, aloc, bloc, &
                                n_nonlocal, nppnl, lmax, is_nonlocal, a_nonlocal, h_nonlocal, c_nonlocal, &
                                has_nlcc, n_nlcc, a_nlcc, c_nlcc)

      TYPE(sgp_potential_type), POINTER                  :: potential
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: name
      CHARACTER(LEN=default_string_length), &
         DIMENSION(4), INTENT(IN), OPTIONAL              :: description
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: aliases
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: elec_conf
      INTEGER, INTENT(IN), OPTIONAL                      :: z
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: zeff, zeff_correction, &
                                                            alpha_core_charge, ccore_charge, &
                                                            core_charge_radius, ppl_radius, &
                                                            ppnl_radius
      LOGICAL, INTENT(IN), OPTIONAL                      :: ecp_local
      INTEGER, INTENT(IN), OPTIONAL                      :: n_local
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_local, c_local
      INTEGER, INTENT(IN), OPTIONAL                      :: nloc
      INTEGER, DIMENSION(1:10), INTENT(IN), OPTIONAL     :: nrloc
      REAL(dp), DIMENSION(1:10), INTENT(IN), OPTIONAL    :: aloc, bloc
      INTEGER, INTENT(IN), OPTIONAL                      :: n_nonlocal, nppnl, lmax
      LOGICAL, DIMENSION(0:5), INTENT(IN), OPTIONAL      :: is_nonlocal
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_nonlocal
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: h_nonlocal
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: c_nonlocal
      LOGICAL, INTENT(IN), OPTIONAL                      :: has_nlcc
      INTEGER, INTENT(IN), OPTIONAL                      :: n_nlcc
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: a_nlcc, c_nlcc

      IF (ASSOCIATED(potential)) THEN

         IF (PRESENT(name)) potential%name = name
         IF (PRESENT(aliases)) potential%aliases = aliases
         IF (PRESENT(description)) potential%description = description

         IF (PRESENT(elec_conf)) THEN
            IF (ASSOCIATED(potential%elec_conf)) THEN
               DEALLOCATE (potential%elec_conf)
            ENDIF
            ALLOCATE (potential%elec_conf(0:SIZE(elec_conf) - 1))
            potential%elec_conf(:) = elec_conf(:)
         ENDIF

         IF (PRESENT(z)) potential%z = z
         IF (PRESENT(zeff)) potential%zeff = zeff
         IF (PRESENT(zeff_correction)) potential%zeff_correction = zeff_correction
         IF (PRESENT(alpha_core_charge)) potential%alpha_core_charge = alpha_core_charge
         IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
         IF (PRESENT(core_charge_radius)) potential%core_charge_radius = core_charge_radius

         IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius
         IF (PRESENT(ppnl_radius)) potential%ppnl_radius = ppnl_radius

         IF (PRESENT(ecp_local)) potential%ecp_local = ecp_local
         IF (PRESENT(n_local)) potential%n_local = n_local
         IF (PRESENT(a_local)) potential%a_local => a_local
         IF (PRESENT(c_local)) potential%c_local => c_local

         IF (PRESENT(nloc)) potential%nloc = nloc
         IF (PRESENT(nrloc)) potential%nrloc = nrloc
         IF (PRESENT(aloc)) potential%aloc = aloc
         IF (PRESENT(bloc)) potential%bloc = bloc

         IF (PRESENT(n_nonlocal)) potential%n_nonlocal = n_nonlocal
         IF (PRESENT(nppnl)) potential%nppnl = nppnl
         IF (PRESENT(lmax)) potential%lmax = lmax
         IF (PRESENT(is_nonlocal)) potential%is_nonlocal(:) = is_nonlocal(:)
         IF (PRESENT(a_nonlocal)) potential%a_nonlocal => a_nonlocal
         IF (PRESENT(c_nonlocal)) potential%c_nonlocal => c_nonlocal
         IF (PRESENT(h_nonlocal)) potential%h_nonlocal => h_nonlocal

         IF (PRESENT(has_nlcc)) potential%has_nlcc = has_nlcc
         IF (PRESENT(n_nlcc)) potential%n_nlcc = n_nlcc
         IF (PRESENT(a_nlcc)) potential%a_nlcc => a_nlcc
         IF (PRESENT(c_nlcc)) potential%c_nlcc => c_nlcc

      ELSE

         CPABORT("The pointer potential is not associated.")

      END IF

   END SUBROUTINE set_sgp_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_all_potential(potential, output_unit)

      ! Write an atomic all-electron potential data set to the output unit.

      ! - Creation (09.02.2002,MK)

      TYPE(all_potential_type), POINTER                  :: potential
      INTEGER, INTENT(in)                                :: output_unit

      CHARACTER(LEN=20)                                  :: string

      IF (output_unit > 0 .AND. ASSOCIATED(potential)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T6,A,T41,A40,/)") &
            "AE Potential information for", ADJUSTR(TRIM(potential%name))
         WRITE (UNIT=output_unit, FMT="(T8,A,T41,A40)") &
            "Description: ", TRIM(potential%description(1)), &
            "             ", TRIM(potential%description(2))
         WRITE (UNIT=output_unit, FMT="(/,T8,A,T69,F12.6)") &
            "Gaussian exponent of the core charge distribution: ", &
            potential%alpha_core_charge
         WRITE (UNIT=string, FMT="(5I4)") potential%elec_conf
         WRITE (UNIT=output_unit, FMT="(T8,A,T61,A20)") &
            "Electronic configuration (s p d ...):", &
            ADJUSTR(TRIM(string))
      END IF

   END SUBROUTINE write_all_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_local_potential(potential, output_unit)

      ! Write an atomic local potential data set to the output unit.

      ! - Creation (24.01.2014,JGH)

      TYPE(local_potential_type), POINTER                :: potential
      INTEGER, INTENT(in)                                :: output_unit

      INTEGER                                            :: igau, ipol

      IF (output_unit > 0 .AND. ASSOCIATED(potential)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T6,A,T41,A40)") &
            "Local Potential information for", ADJUSTR(TRIM(potential%name))
         WRITE (UNIT=output_unit, FMT="(T8,A,T41,A40)") &
            "Description: ", TRIM(potential%description(1))
         DO igau = 1, potential%ngau
            WRITE (UNIT=output_unit, FMT="(T8,A,F12.6,T50,A,4(T68,I2,F10.4))") &
               "Exponent: ", potential%alpha(igau), &
               "Coefficients: ", (2*ipol - 2, potential%cval(igau, ipol), ipol=1, potential%npol)
         END DO
      END IF

   END SUBROUTINE write_local_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_gth_potential(potential, output_unit)

      ! Write an atomic GTH potential data set to the output unit.
      ! - Creation (09.02.2002,MK)

      TYPE(gth_potential_type), POINTER                  :: potential
      INTEGER, INTENT(in)                                :: output_unit

      CHARACTER(LEN=20)                                  :: string
      INTEGER                                            :: i, j, l
      REAL(KIND=dp)                                      :: r

      IF (output_unit > 0 .AND. ASSOCIATED(potential)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T6,A,T41,A40,/)") &
            "GTH Potential information for", ADJUSTR(TRIM(potential%name))
         WRITE (UNIT=output_unit, FMT="(T8,A,T41,A40)") &
            "Description: ", ADJUSTR(TRIM(potential%description(1))), &
            "             ", ADJUSTR(TRIM(potential%description(2))), &
            "             ", ADJUSTR(TRIM(potential%description(3))), &
            "             ", ADJUSTR(TRIM(potential%description(4)))
         WRITE (UNIT=output_unit, FMT="(/,T8,A,T69,F12.6)") &
            "Gaussian exponent of the core charge distribution: ", &
            potential%alpha_core_charge
         WRITE (UNIT=string, FMT="(5I4)") potential%elec_conf
         WRITE (UNIT=output_unit, FMT="(T8,A,T61,A20)") &
            "Electronic configuration (s p d ...):", &
            ADJUSTR(TRIM(string))

         r = 1.0_dp/SQRT(2.0_dp*potential%alpha_ppl)

         WRITE (UNIT=output_unit, FMT="(/,T8,A,/,/,T27,A,/,T21,5F12.6)") &
            "Parameters of the local part of the GTH pseudopotential:", &
            "rloc        C1          C2          C3          C4", &
            r, (potential%cexp_ppl(i)*r**(2*(i - 1)), i=1, potential%nexp_ppl)

         IF (potential%lppnl > -1) THEN
            WRITE (UNIT=output_unit, FMT="(/,T8,A,/,/,T20,A,/)") &
               "Parameters of the non-local part of the GTH pseudopotential:", &
               "l      r(l)      h(i,j,l)"
            DO l = 0, potential%lppnl
               r = SQRT(0.5_dp/potential%alpha_ppnl(l))
               WRITE (UNIT=output_unit, FMT="(T19,I2,5F12.6)") &
                  l, r, (potential%hprj_ppnl(1, j, l), j=1, potential%nprj_ppnl(l))
               DO i = 2, potential%nprj_ppnl(l)
                  WRITE (UNIT=output_unit, FMT="(T33,4F12.6)") &
                     (potential%hprj_ppnl(i, j, l), j=1, potential%nprj_ppnl(l))
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE write_gth_potential

! **************************************************************************************************
!> \brief ...
!> \param potential ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_sgp_potential(potential, output_unit)

      TYPE(sgp_potential_type), POINTER                  :: potential
      INTEGER, INTENT(in)                                :: output_unit

      CHARACTER(LEN=40)                                  :: string
      INTEGER                                            :: i, l

      IF (output_unit > 0 .AND. ASSOCIATED(potential)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T6,A,T41,A40,/)") &
            "SGP Potential information for", ADJUSTR(TRIM(potential%name))
         WRITE (UNIT=output_unit, FMT="(T8,A,T25,A56)") &
            "Description: ", ADJUSTR(TRIM(potential%description(1))), &
            "             ", ADJUSTR(TRIM(potential%description(2))), &
            "             ", ADJUSTR(TRIM(potential%description(3))), &
            "             ", ADJUSTR(TRIM(potential%description(4)))
         WRITE (UNIT=output_unit, FMT="(/,T8,A,T69,F12.6)") &
            "Gaussian exponent of the core charge distribution: ", &
            potential%alpha_core_charge
         WRITE (UNIT=string, FMT="(10I4)") potential%elec_conf
         WRITE (UNIT=output_unit, FMT="(T8,A,T61,A20)") &
            "Electronic configuration (s p d ...):", &
            ADJUSTR(TRIM(string))
         IF (potential%ecp_local) THEN
            IF (potential%nloc > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T8,'Local pseudopotential')")
               WRITE (UNIT=output_unit, FMT="(T20,'r**(n-2)',T50,'Coefficient',T73,'Exponent')")
               DO i = 1, potential%nloc
                  WRITE (UNIT=output_unit, FMT="(T20,I5,T47,F14.8,T69,F12.6)") &
                     potential%nrloc(i), potential%aloc(i), potential%bloc(i)
               END DO
            END IF
         ELSE
            IF (potential%n_local > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T8,'Local pseudopotential')")
               WRITE (UNIT=output_unit, FMT="(T8,A,10(T21,6F10.4,/))") &
                  'Exponents:', potential%a_local(1:potential%n_local)
               WRITE (UNIT=output_unit, FMT="(T8,A,10(T21,6F10.4,/))") &
                  'Coefficients:', potential%c_local(1:potential%n_local)
            END IF
         END IF
         ! nonlocal PP
         IF (potential%n_nonlocal > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T8,'Nonlocal pseudopotential')")
            WRITE (UNIT=output_unit, FMT="(T8,A,T71,I10)") 'Total number of projectors:', potential%nppnl
            WRITE (UNIT=output_unit, FMT="(T8,A,10(T21,6F10.4,/))") &
               'Exponents:', potential%a_nonlocal(1:potential%n_nonlocal)
            DO l = 0, potential%lmax
               WRITE (UNIT=output_unit, FMT="(T8,'Coupling for l=',I4)") l
               WRITE (UNIT=output_unit, FMT="(10(T21,6F10.4,/))") &
                  potential%h_nonlocal(1:potential%n_nonlocal, l)
            END DO
         END IF
         !
         IF (potential%has_nlcc) THEN
            WRITE (UNIT=output_unit, FMT="(/,T8,'Nonlinear Core Correction')")
            WRITE (UNIT=output_unit, FMT="(T8,A,10(T21,6F10.4,/))") &
               'Exponents:', potential%a_nlcc(1:potential%n_nlcc)
            WRITE (UNIT=output_unit, FMT="(T8,A,10(T21,6F10.4,/))") &
               'Coefficients:', potential%c_nlcc(1:potential%n_nlcc)
         END IF
      END IF

   END SUBROUTINE write_sgp_potential

! **************************************************************************************************
!> \brief Copy an all_potential_type to a new, unallocated variable
!> \param pot_in the input potential to copy
!> \param pot_out the newly copied and allocated potential
!> \par History: created 12.2019 (A. Bussy)
! **************************************************************************************************
   SUBROUTINE copy_all_potential(pot_in, pot_out)

      TYPE(all_potential_type), INTENT(IN), POINTER      :: pot_in
      TYPE(all_potential_type), INTENT(INOUT), POINTER   :: pot_out

      CPASSERT(ASSOCIATED(pot_in))
      CALL allocate_all_potential(pot_out)

      pot_out%name = pot_in%name
      pot_out%alpha_core_charge = pot_in%alpha_core_charge
      pot_out%ccore_charge = pot_in%ccore_charge
      pot_out%core_charge_radius = pot_in%core_charge_radius
      pot_out%zeff = pot_in%zeff
      pot_out%zeff_correction = pot_in%zeff_correction
      pot_out%z = pot_in%z

      IF (ASSOCIATED(pot_in%elec_conf)) THEN
         ALLOCATE (pot_out%elec_conf(LBOUND(pot_in%elec_conf, 1):UBOUND(pot_in%elec_conf, 1)))
         pot_out%elec_conf(:) = pot_in%elec_conf(:)
      END IF

   END SUBROUTINE copy_all_potential

! **************************************************************************************************
!> \brief Copy a gth_potential_type to a new, unallocated variable
!> \param pot_in the input potential to copy
!> \param pot_out the newly copied and allocated potential
!> \par History: created 12.2019 (A. Bussy)
! **************************************************************************************************
   SUBROUTINE copy_gth_potential(pot_in, pot_out)

      TYPE(gth_potential_type), INTENT(IN), POINTER      :: pot_in
      TYPE(gth_potential_type), INTENT(INOUT), POINTER   :: pot_out

      CPASSERT(ASSOCIATED(pot_in))
      CALL allocate_gth_potential(pot_out)

      pot_out%name = pot_in%name
      pot_out%aliases = pot_in%aliases
      pot_out%alpha_core_charge = pot_in%alpha_core_charge
      pot_out%alpha_ppl = pot_in%alpha_ppl
      pot_out%ccore_charge = pot_in%ccore_charge
      pot_out%cerf_ppl = pot_in%cerf_ppl
      pot_out%zeff = pot_in%zeff
      pot_out%core_charge_radius = pot_in%core_charge_radius
      pot_out%ppl_radius = pot_in%ppl_radius
      pot_out%ppnl_radius = pot_in%ppnl_radius
      pot_out%zeff_correction = pot_in%zeff_correction
      pot_out%lppnl = pot_in%lppnl
      pot_out%lprj_ppnl_max = pot_in%lprj_ppnl_max
      pot_out%nexp_ppl = pot_in%nexp_ppl
      pot_out%nppnl = pot_in%nppnl
      pot_out%nprj_ppnl_max = pot_in%nprj_ppnl_max
      pot_out%z = pot_in%z
      pot_out%nlcc = pot_in%nlcc
      pot_out%nexp_nlcc = pot_in%nexp_nlcc
      pot_out%lsdpot = pot_in%lsdpot
      pot_out%nexp_lsd = pot_in%nexp_lsd
      pot_out%lpotextended = pot_in%lpotextended
      pot_out%nexp_lpot = pot_in%nexp_lpot

      IF (ASSOCIATED(pot_in%alpha_ppnl)) THEN
         ALLOCATE (pot_out%alpha_ppnl(LBOUND(pot_in%alpha_ppnl, 1):UBOUND(pot_in%alpha_ppnl, 1)))
         pot_out%alpha_ppnl(:) = pot_in%alpha_ppnl(:)
      END IF
      IF (ASSOCIATED(pot_in%cexp_ppl)) THEN
         ALLOCATE (pot_out%cexp_ppl(LBOUND(pot_in%cexp_ppl, 1):UBOUND(pot_in%cexp_ppl, 1)))
         pot_out%cexp_ppl(:) = pot_in%cexp_ppl(:)
      END IF
      IF (ASSOCIATED(pot_in%elec_conf)) THEN
         ALLOCATE (pot_out%elec_conf(LBOUND(pot_in%elec_conf, 1):UBOUND(pot_in%elec_conf, 1)))
         pot_out%elec_conf(:) = pot_in%elec_conf(:)
      END IF
      IF (ASSOCIATED(pot_in%nprj_ppnl)) THEN
         ALLOCATE (pot_out%nprj_ppnl(LBOUND(pot_in%nprj_ppnl, 1):UBOUND(pot_in%nprj_ppnl, 1)))
         pot_out%nprj_ppnl(:) = pot_in%nprj_ppnl(:)
      END IF
      IF (ASSOCIATED(pot_in%cprj)) THEN
         ALLOCATE (pot_out%cprj(LBOUND(pot_in%cprj, 1):UBOUND(pot_in%cprj, 1), &
                                LBOUND(pot_in%cprj, 2):UBOUND(pot_in%cprj, 2)))
         pot_out%cprj(:, :) = pot_in%cprj(:, :)
      END IF
      IF (ASSOCIATED(pot_in%cprj_ppnl)) THEN
         ALLOCATE (pot_out%cprj_ppnl(LBOUND(pot_in%cprj_ppnl, 1):UBOUND(pot_in%cprj_ppnl, 1), &
                                     LBOUND(pot_in%cprj_ppnl, 2):UBOUND(pot_in%cprj_ppnl, 2)))
         pot_out%cprj_ppnl(:, :) = pot_in%cprj_ppnl(:, :)
      END IF
      IF (ASSOCIATED(pot_in%vprj_ppnl)) THEN
         ALLOCATE (pot_out%vprj_ppnl(LBOUND(pot_in%vprj_ppnl, 1):UBOUND(pot_in%vprj_ppnl, 1), &
                                     LBOUND(pot_in%vprj_ppnl, 2):UBOUND(pot_in%vprj_ppnl, 2)))
         pot_out%vprj_ppnl(:, :) = pot_in%vprj_ppnl(:, :)
      END IF
      IF (ASSOCIATED(pot_in%hprj_ppnl)) THEN
         ALLOCATE (pot_out%hprj_ppnl(LBOUND(pot_in%hprj_ppnl, 1):UBOUND(pot_in%hprj_ppnl, 1), &
                                     LBOUND(pot_in%hprj_ppnl, 2):UBOUND(pot_in%hprj_ppnl, 2), &
                                     LBOUND(pot_in%hprj_ppnl, 3):UBOUND(pot_in%hprj_ppnl, 3)))
         pot_out%hprj_ppnl(:, :, :) = pot_in%hprj_ppnl(:, :, :)
      END IF
      IF (ASSOCIATED(pot_in%alpha_nlcc)) THEN
         ALLOCATE (pot_out%alpha_nlcc(LBOUND(pot_in%alpha_nlcc, 1):UBOUND(pot_in%alpha_nlcc, 1)))
         pot_out%alpha_nlcc(:) = pot_in%alpha_nlcc(:)
      END IF
      IF (ASSOCIATED(pot_in%nct_nlcc)) THEN
         ALLOCATE (pot_out%nct_nlcc(LBOUND(pot_in%nct_nlcc, 1):UBOUND(pot_in%nct_nlcc, 1)))
         pot_out%nct_nlcc(:) = pot_in%nct_nlcc(:)
      END IF
      IF (ASSOCIATED(pot_in%cval_nlcc)) THEN
         ALLOCATE (pot_out%cval_nlcc(LBOUND(pot_in%cval_nlcc, 1):UBOUND(pot_in%cval_nlcc, 1), &
                                     LBOUND(pot_in%cval_nlcc, 2):UBOUND(pot_in%cval_nlcc, 2)))
         pot_out%cval_nlcc(:, :) = pot_in%cval_nlcc(:, :)
      END IF
      IF (ASSOCIATED(pot_in%alpha_lsd)) THEN
         ALLOCATE (pot_out%alpha_lsd(LBOUND(pot_in%alpha_lsd, 1):UBOUND(pot_in%alpha_lsd, 1)))
         pot_out%alpha_lsd(:) = pot_in%alpha_lsd(:)
      END IF
      IF (ASSOCIATED(pot_in%nct_lsd)) THEN
         ALLOCATE (pot_out%nct_lsd(LBOUND(pot_in%nct_lsd, 1):UBOUND(pot_in%nct_lsd, 1)))
         pot_out%nct_lsd(:) = pot_in%nct_lsd(:)
      END IF
      IF (ASSOCIATED(pot_in%cval_lsd)) THEN
         ALLOCATE (pot_out%cval_lsd(LBOUND(pot_in%cval_lsd, 1):UBOUND(pot_in%cval_lsd, 1), &
                                    LBOUND(pot_in%cval_lsd, 2):UBOUND(pot_in%cval_lsd, 2)))
         pot_out%cval_lsd(:, :) = pot_in%cval_lsd(:, :)
      END IF
      IF (ASSOCIATED(pot_in%alpha_lpot)) THEN
         ALLOCATE (pot_out%alpha_lpot(LBOUND(pot_in%alpha_lpot, 1):UBOUND(pot_in%alpha_lpot, 1)))
         pot_out%alpha_lpot(:) = pot_in%alpha_lpot(:)
      END IF
      IF (ASSOCIATED(pot_in%nct_lpot)) THEN
         ALLOCATE (pot_out%nct_lpot(LBOUND(pot_in%nct_lpot, 1):UBOUND(pot_in%nct_lpot, 1)))
         pot_out%nct_lpot(:) = pot_in%nct_lpot(:)
      END IF
      IF (ASSOCIATED(pot_in%cval_lpot)) THEN
         ALLOCATE (pot_out%cval_lpot(LBOUND(pot_in%cval_lpot, 1):UBOUND(pot_in%cval_lpot, 1), &
                                     LBOUND(pot_in%cval_lpot, 2):UBOUND(pot_in%cval_lpot, 2)))
         pot_out%cval_lpot(:, :) = pot_in%cval_lpot(:, :)
      END IF

   END SUBROUTINE copy_gth_potential

! **************************************************************************************************
!> \brief Copy a sgp_potential_type to a new, unallocated variable
!> \param pot_in the input potential to copy
!> \param pot_out the newly copied and allocated potential
!> \par History: created 12.2019 (A. Bussy)
! **************************************************************************************************
   SUBROUTINE copy_sgp_potential(pot_in, pot_out)

      TYPE(sgp_potential_type), INTENT(IN)               :: pot_in
      TYPE(sgp_potential_type), INTENT(INOUT), POINTER   :: pot_out

      CALL allocate_sgp_potential(pot_out)

      pot_out%name = pot_in%name
      pot_out%aliases = pot_in%aliases
      pot_out%z = pot_in%z
      pot_out%zeff = pot_in%zeff
      pot_out%zeff_correction = pot_in%zeff_correction
      pot_out%alpha_core_charge = pot_in%alpha_core_charge
      pot_out%ccore_charge = pot_in%ccore_charge
      pot_out%core_charge_radius = pot_in%core_charge_radius
      pot_out%ppl_radius = pot_in%ppl_radius
      pot_out%ppnl_radius = pot_in%ppnl_radius
      pot_out%ecp_local = pot_in%ecp_local
      pot_out%n_local = pot_in%n_local
      pot_out%nloc = pot_in%nloc
      pot_out%nrloc = pot_in%nrloc
      pot_out%aloc = pot_in%aloc
      pot_out%bloc = pot_in%bloc
      pot_out%n_nonlocal = pot_in%n_nonlocal
      pot_out%nppnl = pot_in%nppnl
      pot_out%lmax = pot_in%lmax
      pot_out%is_nonlocal = pot_in%is_nonlocal
      pot_out%has_nlcc = pot_in%has_nlcc
      pot_out%n_nlcc = pot_in%n_nlcc

      IF (ASSOCIATED(pot_in%elec_conf)) THEN
         ALLOCATE (pot_out%elec_conf(LBOUND(pot_in%elec_conf, 1):UBOUND(pot_in%elec_conf, 1)))
         pot_out%elec_conf(:) = pot_in%elec_conf(:)
      END IF
      IF (ASSOCIATED(pot_in%a_local)) THEN
         ALLOCATE (pot_out%a_local(LBOUND(pot_in%a_local, 1):UBOUND(pot_in%a_local, 1)))
         pot_out%a_local(:) = pot_in%a_local(:)
      END IF
      IF (ASSOCIATED(pot_in%c_local)) THEN
         ALLOCATE (pot_out%c_local(LBOUND(pot_in%c_local, 1):UBOUND(pot_in%c_local, 1)))
         pot_out%c_local(:) = pot_in%c_local(:)
      END IF
      IF (ASSOCIATED(pot_in%a_nonlocal)) THEN
         ALLOCATE (pot_out%a_nonlocal(LBOUND(pot_in%a_nonlocal, 1):UBOUND(pot_in%a_nonlocal, 1)))
         pot_out%a_nonlocal(:) = pot_in%a_nonlocal(:)
      END IF
      IF (ASSOCIATED(pot_in%h_nonlocal)) THEN
         ALLOCATE (pot_out%h_nonlocal(LBOUND(pot_in%h_nonlocal, 1):UBOUND(pot_in%h_nonlocal, 1), &
                                      LBOUND(pot_in%h_nonlocal, 2):UBOUND(pot_in%h_nonlocal, 2)))
         pot_out%h_nonlocal(:, :) = pot_in%h_nonlocal(:, :)
      END IF
      IF (ASSOCIATED(pot_in%c_nonlocal)) THEN
         ALLOCATE (pot_out%c_nonlocal(LBOUND(pot_in%c_nonlocal, 1):UBOUND(pot_in%c_nonlocal, 1), &
                                      LBOUND(pot_in%c_nonlocal, 2):UBOUND(pot_in%c_nonlocal, 2), &
                                      LBOUND(pot_in%c_nonlocal, 3):UBOUND(pot_in%c_nonlocal, 3)))
         pot_out%c_nonlocal(:, :, :) = pot_in%c_nonlocal(:, :, :)
      END IF
      IF (ASSOCIATED(pot_in%cprj_ppnl)) THEN
         ALLOCATE (pot_out%cprj_ppnl(LBOUND(pot_in%cprj_ppnl, 1):UBOUND(pot_in%cprj_ppnl, 1), &
                                     LBOUND(pot_in%cprj_ppnl, 2):UBOUND(pot_in%cprj_ppnl, 2)))
         pot_out%cprj_ppnl(:, :) = pot_in%cprj_ppnl(:, :)
      END IF
      IF (ASSOCIATED(pot_in%vprj_ppnl)) THEN
         ALLOCATE (pot_out%vprj_ppnl(LBOUND(pot_in%vprj_ppnl, 1):UBOUND(pot_in%vprj_ppnl, 1)))
         pot_out%vprj_ppnl(:) = pot_in%vprj_ppnl(:)
      END IF
      IF (ASSOCIATED(pot_in%a_nlcc)) THEN
         ALLOCATE (pot_out%a_nlcc(LBOUND(pot_in%a_nlcc, 1):UBOUND(pot_in%a_nlcc, 1)))
         pot_out%a_nlcc(:) = pot_in%a_nlcc(:)
      END IF
      IF (ASSOCIATED(pot_in%c_nlcc)) THEN
         ALLOCATE (pot_out%c_nlcc(LBOUND(pot_in%c_nlcc, 1):UBOUND(pot_in%c_nlcc, 1)))
         pot_out%c_nlcc(:) = pot_in%c_nlcc(:)
      END IF

   END SUBROUTINE copy_sgp_potential

END MODULE external_potential_types
