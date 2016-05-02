!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the MME method for periodic ERI's to be used in CP2K
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE cp_eri_mme_interface

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_error_control,           ONLY: eri_mme_set_params_custom
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate_low,&
                                              eri_mme_2c_integrate_set
   USE eri_mme_types,                   ONLY: eri_mme_param
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_eri_mme_interface'

   PUBLIC :: eri_mme_set_params, &
             eri_mme_2c_integrate, &
             eri_mme_2c_integrate_s

   INTERFACE eri_mme_set_params
      MODULE PROCEDURE eri_mme_set_params_from_basis
      MODULE PROCEDURE eri_mme_set_params_custom
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Reset counters before integration.
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_integrate_prepare(param)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param

      param%G_count = 0
      param%R_count = 0
   END SUBROUTINE eri_mme_integrate_prepare

! **************************************************************************************************
!> \brief print some info after integration is done.
!> \param param ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_integrate_finalize(param, para_env)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, param%mme_section, "ERI_MME_INFO", &
                                     extension=".eri_mme")
      CALL mp_sum(param%G_count, para_env%group)
      CALL mp_sum(param%R_count, para_env%group)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/T2, A)') "ERI_MME| Percentage of integrals evaluated in"
         WRITE (unit_nr, '(T2, A, T76, F5.1)') "ERI_MME|   G space:", &
            100.0_dp*param%G_count/(param%G_count+param%R_count)
         WRITE (unit_nr, '(T2, A, T76, F5.1/)') "ERI_MME|   R space:", &
            100.0_dp*param%R_count/(param%G_count+param%R_count)
      ENDIF
      CALL cp_print_key_finished_output(unit_nr, logger, param%mme_section, "ERI_MME_INFO")
   END SUBROUTINE eri_mme_integrate_finalize

! **************************************************************************************************
!> \brief high-level integration routine for 2c integrals over CP2K basis sets.
!>        Parallelization over pairs of sets.
!> \param param ...
!> \param para_env ...
!> \param qs_env ...
!> \param hab ...
!> \param basis_type_a ...
!> \param basis_type_b ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate(param, para_env, qs_env, hab, basis_type_a, basis_type_b)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_a, basis_type_b

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_2c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iatom, ikind, iset, jatom, jkind, jset, my_setpair, n_setpair, natom, &
         nkind, nseta, nseta_total, nsetb, nsetb_total, offset_a, offset_b, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: map_it_here
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sphi_a, sphi_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL eri_mme_integrate_prepare(param)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      ! get total number of set pairs to integrate
      nseta_total = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)
         nseta_total = nseta_total+basis_set_a%nset
      ENDDO

      nsetb_total = 0
      DO jatom = 1, natom
         jkind = kind_of(jatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)
         nsetb_total = nsetb_total+basis_set_b%nset
      ENDDO

      n_setpair = nseta_total*nsetb_total

      ! integrate
      my_setpair = 0

      offset_a = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type_a)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         npgfa => basis_set_a%npgf

         ra(:) = pbc(particle_set(iatom)%r, cell)

         DO iset = 1, nseta
            offset_a = offset_a+nsgfa(iset)
            sgfa = first_sgfa(1, iset)

            offset_b = 0
            DO jatom = 1, natom
               jkind = kind_of(jatom)
               CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type=basis_type_b)

               first_sgfb => basis_set_b%first_sgf
               lb_max => basis_set_b%lmax
               lb_min => basis_set_b%lmin
               nsetb = basis_set_b%nset
               nsgfb => basis_set_b%nsgf_set
               sphi_b => basis_set_b%sphi
               zetb => basis_set_b%zet
               npgfb => basis_set_b%npgf

               rb(:) = pbc(particle_set(jatom)%r, cell)
               rab(:) = ra(:)-rb(:) ! pbc not needed
               dab = SQRT(rab(1)**2+rab(2)**2+rab(3)**2)

               DO jset = 1, nsetb
                  offset_b = offset_b+nsgfb(jset)
                  sgfb = first_sgfb(1, jset)
                  my_setpair = my_setpair+1
                  map_it_here = MODULO(my_setpair, para_env%num_pe) == para_env%mepos

                  IF (map_it_here) THEN
                     CALL eri_mme_2c_integrate_set(param, la_min(iset), la_max(iset), lb_min(jset), lb_max(jset), &
                                                   npgfa(iset), zeta(:, iset), npgfb(jset), zetb(:, jset), &
                                                   ra, rb, &
                                                   hab, offset_a-nsgfa(iset)+1, offset_b-nsgfb(jset)+1, &
                                                   sphi_a, sgfa, nsgfa(iset), &
                                                   sphi_b, sgfb, nsgfb(jset))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL mp_sum(hab, para_env%group)
      CALL eri_mme_integrate_finalize(param, para_env)

      CALL timestop(handle)
   END SUBROUTINE eri_mme_2c_integrate

! **************************************************************************************************
!> \brief high-level integration routine for 2c integrals over s-type functions.
!>        Parallelization over pairs of functions.
!> \param param ...
!> \param zeta ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param hab ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate_s(param, zeta, zetb, ra, rb, hab, para_env)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_2c_integrate_s', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipgf, ipgf_prod, jpgf, &
                                                            npgf_prod, npgfa, npgfb
      INTEGER, DIMENSION(2)                              :: limits
      REAL(KIND=dp), DIMENSION(3)                        :: rab

      CALL timeset(routineN, handle)
      CALL eri_mme_integrate_prepare(param)

      hab(:, :) = 0.0_dp

      npgfa = SIZE(zeta)
      npgfb = SIZE(zetb)
      npgf_prod = npgfa*npgfb ! total number of integrals

      limits = get_limit(npgf_prod, para_env%num_pe, para_env%mepos)

      DO ipgf_prod = limits(1), limits(2)
         ipgf = (ipgf_prod-1)/npgfb+1
         jpgf = MOD(ipgf_prod-1, npgfb)+1
         rab(:) = ra(:, ipgf)-rb(:, jpgf)
         CALL eri_mme_2c_integrate_low(param, 0, 0, 0, 0, zeta(ipgf), &
                                       zetb(jpgf), rab, hab, ipgf-1, jpgf-1)
      ENDDO

      CALL mp_sum(hab, para_env%group)
      CALL eri_mme_integrate_finalize(param, para_env)
      CALL timestop(handle)

   END SUBROUTINE eri_mme_2c_integrate_s

! **************************************************************************************************
!> \brief Set parameters for MME method by deriving basis info from cp2k basis set.
!>        Cutoff can be auto-calibrated to minimize total error.
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
      CALL eri_mme_set_params_custom(param, cell, &
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

      INTEGER                                            :: handle, i, ikind, ipgf, iset, j, l_m, &
                                                            l_zet, nkind, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf
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
         npgf => basis_set%npgf
         nset = basis_set%nset
         l_m = MAX(l_m, MAXVAL(basis_set%l(:, :)))
         DO iset = 1, nset
            DO ipgf = 1, npgf(iset)
               zet_m = MAX(zet_m, MAXVAL(basis_set%zet(1:ipgf, iset)))
               IF (zet_mm .LT. 0.0_dp) THEN
                  zet_mm = MINVAL(basis_set%zet(1:ipgf, iset))
               ELSE
                  zet_mm = MIN(zet_mm, MINVAL(basis_set%zet(1:ipgf, iset)))
               ENDIF
            ENDDO
         ENDDO
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
      CPASSERT(zet_l .GT. 0.0_dp .AND. l_zet .GE. 0)

      zet_c(:) = [zet_m, zet_l]
      l_c(:) = [l_zet, l_m]

      CALL timestop(handle)
   END SUBROUTINE error_est_pgf_params_from_basis

END MODULE cp_eri_mme_interface
