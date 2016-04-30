!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the MME method for periodic ERI's to be used in CP2K
!> \note  accuracy may not be sufficient yet because current implementation
!>        can only use up to 15 points in minimax approximation.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE cp_eri_mme_interface

   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_error_control,           ONLY: eri_mme_set_params_custom
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate_low
   USE eri_mme_types,                   ONLY: eri_mme_param
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: ncoset
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_eri_mme_interface'

   PUBLIC :: eri_mme_set_params, eri_mme_2c_integrate

   INTERFACE eri_mme_set_params
      MODULE PROCEDURE eri_mme_set_params_from_basis
      MODULE PROCEDURE eri_mme_set_params_custom
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief High-level integration routine computing 2-center electron repulsion
!>        integrals (P|Q) with P and Q Hermite-Gaussian basis functions.
!> \param param  parameters (cutoff, minimax approx, ...)
!> \param la_min minimum l quantum numbers for P
!> \param la_max maximum l quantum numbers for P
!> \param lb_min minimum l quantum numbers for Q
!> \param lb_max maximum l quantum numbers for Q
!> \param zeta   exponents for P
!> \param zetb   exponents for Q
!> \param ra     Atomic centers of P
!> \param rb     Atomic centers of Q
!> \param hab    Matrix containing integral values
!> \param para_env ...
!> \note  adapted for use in qmmm_image_charge
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate(param, la_min, la_max, lb_min, lb_max, zeta, zetb, &
                                   ra, rb, hab, para_env)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, DIMENSION(:), INTENT(IN)                  :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_2c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: b_first, b_last, handle, i_b, ipgf, &
                                                            jpgf, my_p, n_b, n_b_left, n_bp, n_p, &
                                                            npgfa, npgfb, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: na, nb
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      CALL timeset(routineN, handle)

      param%G_count = 0
      param%R_count = 0
      logger => cp_get_default_logger()

      hab(:, :) = 0.0_dp
      npgfa = SIZE(la_max)
      npgfb = SIZE(lb_max)

      ALLOCATE (na(npgfa), nb(npgfb))

      my_p = para_env%mepos ! mpi rank
      n_p = para_env%num_pe ! total number of processes
      n_b = npgfa*npgfb ! total number of matrix blocks
      n_bp = n_b/n_p ! minimum number of blocks per process
      n_b_left = MOD(n_b, n_p) ! leftovers

      IF (my_p .LT. n_b_left) THEN
         b_first = (n_bp+1)*my_p+1
         b_last = (n_bp+1)*(my_p+1)
      ELSE
         b_first = (n_bp)*my_p+1+n_b_left
         b_last = (n_bp)*(my_p+1)+n_b_left
      ENDIF

      na(1) = 0
      DO ipgf = 1, npgfa-1
         na(ipgf+1) = na(ipgf)+ncoset(la_max(ipgf))-ncoset(la_min(ipgf)-1)
      ENDDO

      nb(1) = 0
      DO jpgf = 1, npgfb-1
         nb(jpgf+1) = nb(jpgf)+ncoset(lb_max(jpgf))-ncoset(lb_min(jpgf)-1)
      ENDDO

      DO i_b = b_first, b_last
         ipgf = (i_b-1)/npgfb+1
         jpgf = MOD(i_b-1, npgfb)+1
         rab(:) = ra(:, ipgf)-rb(:, jpgf)
         CALL eri_mme_2c_integrate_low(param, &
                                       la_min(ipgf), la_max(ipgf), lb_min(jpgf), lb_max(jpgf), &
                                       zeta(ipgf), zetb(jpgf), rab, hab, na(ipgf), nb(jpgf))
      ENDDO
      CALL mp_sum(hab, para_env%group)
      CALL mp_sum(param%G_count, para_env%group)
      CALL mp_sum(param%R_count, para_env%group)

      unit_nr = cp_print_key_unit_nr(logger, param%mme_section, "ERI_MME_INFO", &
                                     extension=".eri_mme")
      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(/T2, A)') "ERI_MME| Percentage of integrals evaluated in"
         WRITE (unit_nr, '(T2, A, T76, F5.1)') "ERI_MME|   G space:", &
            100.0_dp*param%G_count/(param%G_count+param%R_count)
         WRITE (unit_nr, '(T2, A, T76, F5.1/)') "ERI_MME|   R space:", &
            100.0_dp*param%R_count/(param%G_count+param%R_count)
      ENDIF

      CALL cp_print_key_finished_output(unit_nr, logger, param%mme_section, "ERI_MME_INFO")
      CALL timestop(handle)
   END SUBROUTINE eri_mme_2c_integrate

! **************************************************************************************************
!> \brief Set parameters for MME method by deriving basis info from cp2k basis set.
!>        Cutoff can be auto-calibrated to minimize total error.
!>        Switch between G / R space representation can be auto-calibrated to
!>        minimize run time.
!> \param param ...
!> \param cell ...
!> \param qs_kind_set ...
!> \param basis_type ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params_from_basis(param, cell, qs_kind_set, basis_type, para_env)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params_from_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: l_err_cutoff
      REAL(KIND=dp)                                      :: zet_err_minimax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zet_err_cutoff

      CALL timeset(routineN, handle)

      CALL error_est_pgf_params_from_basis(qs_kind_set, basis_type, &
                                           zet_err_minimax, zet_err_cutoff, l_err_cutoff)

      param%l_max = MAXVAL(l_err_cutoff)
      CALL eri_mme_set_params(param, cell, &
                              zet_err_minimax, zet_err_cutoff, l_err_cutoff, &
                              param%l_max, para_env)

      CALL timestop(handle)
   END SUBROUTINE eri_mme_set_params_from_basis

! **************************************************************************************************
!> \brief Get basis parameters for estimating cutoff and minimax error from cp2k basis
!> \param qs_kind_set ...
!> \param basis_type ...
!> \param zet_mm Smallest exponent, used to estimate error due to minimax approx.
!> \param zet_c  contains max. exponent and the largest exponent for max. l,
!>               used to estimate cutoff error
!> \param l_c    contains the largest l for max. exponent and max. l,
!>               used to estimate cutoff error
! **************************************************************************************************
   SUBROUTINE error_est_pgf_params_from_basis(qs_kind_set, basis_type, zet_mm, zet_c, l_c)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      REAL(KIND=dp), INTENT(OUT)                         :: zet_mm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: zet_c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: l_c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'error_est_pgf_params_from_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ikind, j, l_m, l_zet, nkind
      REAL(KIND=dp)                                      :: zet_l, zet_m
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL timeset(routineN, handle)

      l_m = 0
      zet_m = 0.0_dp
      l_zet = -1
      zet_l = -1.0_dp
      zet_mm = -1.0_dp

      nkind = SIZE(qs_kind_set)

      ! 1) get global max l and max zet
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, basis_type=basis_type)
         CPASSERT(ASSOCIATED(basis_set))
         l_m = MAX(l_m, MAXVAL(basis_set%l(:, :)))
         zet_m = MAX(zet_m, MAXVAL(basis_set%zet(:, :)))
         IF (zet_mm .LT. 0.0_dp) THEN
            zet_mm = MINVAL(basis_set%zet(:, :))
         ELSE
            zet_mm = MIN(zet_mm, MINVAL(basis_set%zet(:, :)))
         ENDIF
      ENDDO

      ! 2) get largest zet for max l and largest l for max zet
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, basis_type=basis_type)
         DO i = LBOUND(basis_set%l, 1), UBOUND(basis_set%l, 1)
            DO j = LBOUND(basis_set%l, 2), UBOUND(basis_set%l, 2)
               IF ((basis_set%zet(i, j) .EQ. zet_m) .AND. (basis_set%l(i, j) .GT. l_zet)) THEN
                  l_zet = basis_set%l(i, j)
               ENDIF
               IF ((basis_set%l(i, j) .EQ. l_m) .AND. (basis_set%zet(i, j) .GT. zet_l)) THEN
                  zet_l = basis_set%zet(i, j)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ALLOCATE (l_c(2))
      ALLOCATE (zet_c(2))
      CPASSERT(zet_l .GT. 0.0_dp .AND. l_zet .GT. 0)

      zet_c(:) = [zet_m, zet_l]
      l_c(:) = [l_zet, l_m]

      CALL timestop(handle)
   END SUBROUTINE error_est_pgf_params_from_basis

END MODULE cp_eri_mme_interface
