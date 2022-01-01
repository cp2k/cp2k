!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of the xTB parameter types.
!> \author JGH (10.2018)
! **************************************************************************************************
! To be done:
! 1) Ewald defaults options for GMAX, ALPHA, RCUT
! 2) QM/MM debugging of forces -- done
! 3) Periodic displacement field (debugging)
! 4) Check for RTP and EMD
! 5) Wannier localization
! 6) Charge Mixing methods: Broyden/Pulay (more debugging needed, also add to DFTB)
! **************************************************************************************************
MODULE xtb_types

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_types'

! **************************************************************************************************
   TYPE xtb_atom_type
      ! PRIVATE
      CHARACTER(LEN=default_string_length)   :: typ
      CHARACTER(LEN=default_string_length)   :: aname
      CHARACTER(LEN=2)                       :: symbol
      LOGICAL                                :: defined
      INTEGER                                :: z !atomic number
      REAL(KIND=dp)                          :: zeff !effective core charge
      INTEGER                                :: natorb !number of orbitals
      INTEGER                                :: lmax !max angular momentum
      !
      REAL(KIND=dp)                          :: rcut !cutoff radius for sr-Coulomb
      REAL(KIND=dp)                          :: rcov !covalent radius
      REAL(KIND=dp)                          :: electronegativity !electronegativity
      !
      REAL(KIND=dp)                          :: kx !scaling for halogen term
      !
      REAL(KIND=dp)                          :: eta !Atomic Hubbard parameter
      REAL(KIND=dp)                          :: xgamma !charge derivative of eta
      REAL(KIND=dp)                          :: alpha !exponential scaling parameter for repulsion potential
      REAL(KIND=dp)                          :: zneff !effective core charge for repulsion potential
      ! shell specific parameters
      INTEGER                                :: nshell !number of orbital shells
      INTEGER, DIMENSION(5)                  :: nval ! n-quantum number of shell i
      INTEGER, DIMENSION(5)                  :: lval ! l-quantum number of shell i
      INTEGER, DIMENSION(5)                  :: occupation ! occupation of shell i
      REAL(KIND=dp), DIMENSION(5)            :: kpoly
      REAL(KIND=dp), DIMENSION(5)            :: kappa
      REAL(KIND=dp), DIMENSION(5)            :: hen
      REAL(KIND=dp), DIMENSION(5)            :: zeta
      ! AO to shell pointer
      INTEGER, DIMENSION(25)                 :: nao, lao
      ! Upper limit of Mulliken charge
      REAL(KIND=dp)                          :: chmax
   END TYPE xtb_atom_type

! *** Public data types ***

   PUBLIC :: xtb_atom_type, get_xtb_atom_param, set_xtb_atom_param, write_xtb_atom_param
   PUBLIC :: allocate_xtb_atom_param, deallocate_xtb_atom_param

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param xtb_parameter ...
! **************************************************************************************************
   SUBROUTINE allocate_xtb_atom_param(xtb_parameter)

      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter

      IF (ASSOCIATED(xtb_parameter)) &
         CALL deallocate_xtb_atom_param(xtb_parameter)

      ALLOCATE (xtb_parameter)

      xtb_parameter%defined = .FALSE.
      xtb_parameter%aname = ""
      xtb_parameter%symbol = ""
      xtb_parameter%typ = "NONE"
      xtb_parameter%z = -1
      xtb_parameter%zeff = -1.0_dp
      xtb_parameter%natorb = 0
      xtb_parameter%lmax = -1
      xtb_parameter%rcut = 0.0_dp
      xtb_parameter%rcov = 0.0_dp
      xtb_parameter%electronegativity = 0.0_dp
      xtb_parameter%kx = -100.0_dp
      xtb_parameter%eta = 0.0_dp
      xtb_parameter%xgamma = 0.0_dp
      xtb_parameter%alpha = 0.0_dp
      xtb_parameter%zneff = 0.0_dp
      xtb_parameter%nshell = 0
      xtb_parameter%nval = 0
      xtb_parameter%lval = 0
      xtb_parameter%occupation = 0
      xtb_parameter%kpoly = 0.0_dp
      xtb_parameter%kappa = 0.0_dp
      xtb_parameter%hen = 0.0_dp
      xtb_parameter%zeta = 0.0_dp
      xtb_parameter%nao = 0
      xtb_parameter%lao = 0
      xtb_parameter%chmax = 0.0_dp

   END SUBROUTINE allocate_xtb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param xtb_parameter ...
! **************************************************************************************************
   SUBROUTINE deallocate_xtb_atom_param(xtb_parameter)

      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter

      CPASSERT(ASSOCIATED(xtb_parameter))
      DEALLOCATE (xtb_parameter)

   END SUBROUTINE deallocate_xtb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param xtb_parameter ...
!> \param symbol ...
!> \param aname ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param nao ...
!> \param lao ...
!> \param rcut ...
!> \param rcov ...
!> \param kx ...
!> \param eta ...
!> \param xgamma ...
!> \param alpha ...
!> \param zneff ...
!> \param nshell ...
!> \param nval ...
!> \param lval ...
!> \param kpoly ...
!> \param kappa ...
!> \param hen ...
!> \param zeta ...
!> \param occupation ...
!> \param electronegativity ...
!> \param chmax ...
! **************************************************************************************************
   SUBROUTINE get_xtb_atom_param(xtb_parameter, symbol, aname, typ, defined, z, zeff, natorb, lmax, nao, lao, &
                                 rcut, rcov, kx, eta, xgamma, alpha, zneff, nshell, nval, lval, kpoly, kappa, &
                                 hen, zeta, occupation, electronegativity, chmax)

      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter
      CHARACTER(LEN=2), INTENT(OUT), OPTIONAL            :: symbol
      CHARACTER(LEN=default_string_length), &
         INTENT(OUT), OPTIONAL                           :: aname, typ
      LOGICAL, INTENT(OUT), OPTIONAL                     :: defined
      INTEGER, INTENT(OUT), OPTIONAL                     :: z
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: zeff
      INTEGER, INTENT(OUT), OPTIONAL                     :: natorb, lmax
      INTEGER, DIMENSION(25), INTENT(OUT), OPTIONAL      :: nao, lao
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: rcut, rcov, kx, eta, xgamma, alpha, zneff
      INTEGER, INTENT(OUT), OPTIONAL                     :: nshell
      INTEGER, DIMENSION(5), INTENT(OUT), OPTIONAL       :: nval, lval
      REAL(KIND=dp), DIMENSION(5), INTENT(OUT), OPTIONAL :: kpoly, kappa, hen, zeta
      INTEGER, DIMENSION(5), INTENT(OUT), OPTIONAL       :: occupation
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: electronegativity, chmax

      CPASSERT(ASSOCIATED(xtb_parameter))

      IF (PRESENT(symbol)) symbol = xtb_parameter%symbol
      IF (PRESENT(aname)) aname = xtb_parameter%aname
      IF (PRESENT(typ)) typ = xtb_parameter%typ
      IF (PRESENT(defined)) defined = xtb_parameter%defined
      IF (PRESENT(z)) z = xtb_parameter%z
      IF (PRESENT(zeff)) zeff = xtb_parameter%zeff
      IF (PRESENT(natorb)) natorb = xtb_parameter%natorb
      IF (PRESENT(lmax)) lmax = xtb_parameter%lmax
      IF (PRESENT(nao)) nao = xtb_parameter%nao
      IF (PRESENT(lao)) lao = xtb_parameter%lao
      !
      IF (PRESENT(rcut)) rcut = xtb_parameter%rcut
      IF (PRESENT(rcov)) rcov = xtb_parameter%rcov
      IF (PRESENT(kx)) kx = xtb_parameter%kx
      IF (PRESENT(electronegativity)) electronegativity = xtb_parameter%electronegativity
      IF (PRESENT(eta)) eta = xtb_parameter%eta
      IF (PRESENT(xgamma)) xgamma = xtb_parameter%xgamma
      IF (PRESENT(alpha)) alpha = xtb_parameter%alpha
      IF (PRESENT(zneff)) zneff = xtb_parameter%zneff
      IF (PRESENT(nshell)) nshell = xtb_parameter%nshell
      IF (PRESENT(nval)) nval = xtb_parameter%nval
      IF (PRESENT(lval)) lval = xtb_parameter%lval
      IF (PRESENT(occupation)) occupation = xtb_parameter%occupation
      IF (PRESENT(kpoly)) kpoly = xtb_parameter%kpoly
      IF (PRESENT(kappa)) kappa = xtb_parameter%kappa
      IF (PRESENT(hen)) hen = xtb_parameter%hen
      IF (PRESENT(zeta)) zeta = xtb_parameter%zeta
      IF (PRESENT(chmax)) chmax = xtb_parameter%chmax

   END SUBROUTINE get_xtb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param xtb_parameter ...
!> \param aname ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param lmax ...
!> \param nao ...
!> \param lao ...
!> \param rcut ...
!> \param rcov ...
!> \param kx ...
!> \param eta ...
!> \param xgamma ...
!> \param alpha ...
!> \param zneff ...
!> \param nshell ...
!> \param nval ...
!> \param lval ...
!> \param kpoly ...
!> \param kappa ...
!> \param hen ...
!> \param zeta ...
!> \param electronegativity ...
!> \param occupation ...
!> \param chmax ...
! **************************************************************************************************
   SUBROUTINE set_xtb_atom_param(xtb_parameter, aname, typ, defined, z, zeff, natorb, lmax, nao, lao, &
                                 rcut, rcov, kx, eta, xgamma, alpha, zneff, nshell, nval, lval, kpoly, kappa, &
                                 hen, zeta, electronegativity, occupation, chmax)

      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: aname, typ
      LOGICAL, INTENT(IN), OPTIONAL                      :: defined
      INTEGER, INTENT(IN), OPTIONAL                      :: z
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: zeff
      INTEGER, INTENT(IN), OPTIONAL                      :: natorb, lmax
      INTEGER, DIMENSION(25), INTENT(IN), OPTIONAL       :: nao, lao
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: rcut, rcov, kx, eta, xgamma, alpha, zneff
      INTEGER, INTENT(IN), OPTIONAL                      :: nshell
      INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL        :: nval, lval
      REAL(KIND=dp), DIMENSION(5), INTENT(IN), OPTIONAL  :: kpoly, kappa, hen, zeta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: electronegativity
      INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL        :: occupation
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: chmax

      CPASSERT(ASSOCIATED(xtb_parameter))

      IF (PRESENT(aname)) xtb_parameter%aname = aname
      IF (PRESENT(typ)) xtb_parameter%typ = typ
      IF (PRESENT(defined)) xtb_parameter%defined = defined
      IF (PRESENT(z)) xtb_parameter%z = z
      IF (PRESENT(zeff)) xtb_parameter%zeff = zeff
      IF (PRESENT(natorb)) xtb_parameter%natorb = natorb
      IF (PRESENT(lmax)) xtb_parameter%lmax = lmax
      IF (PRESENT(nao)) xtb_parameter%nao = nao
      IF (PRESENT(lao)) xtb_parameter%lao = lao
      !
      IF (PRESENT(rcut)) xtb_parameter%rcut = rcut
      IF (PRESENT(rcov)) xtb_parameter%rcov = rcov
      IF (PRESENT(kx)) xtb_parameter%kx = kx
      IF (PRESENT(electronegativity)) xtb_parameter%electronegativity = electronegativity
      IF (PRESENT(eta)) xtb_parameter%eta = eta
      IF (PRESENT(xgamma)) xtb_parameter%xgamma = xgamma
      IF (PRESENT(alpha)) xtb_parameter%alpha = alpha
      IF (PRESENT(zneff)) xtb_parameter%zneff = zneff
      IF (PRESENT(nshell)) xtb_parameter%nshell = nshell
      IF (PRESENT(nval)) xtb_parameter%nval = nval
      IF (PRESENT(lval)) xtb_parameter%lval = lval
      IF (PRESENT(occupation)) xtb_parameter%occupation = occupation
      IF (PRESENT(kpoly)) xtb_parameter%kpoly = kpoly
      IF (PRESENT(kappa)) xtb_parameter%kappa = kappa
      IF (PRESENT(hen)) xtb_parameter%hen = hen
      IF (PRESENT(zeta)) xtb_parameter%zeta = zeta
      IF (PRESENT(chmax)) xtb_parameter%chmax = chmax

   END SUBROUTINE set_xtb_atom_param

! **************************************************************************************************
!> \brief ...
!> \param xtb_parameter ...
!> \param subsys_section ...
! **************************************************************************************************
   SUBROUTINE write_xtb_atom_param(xtb_parameter, subsys_section)

      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter
      TYPE(section_vals_type), POINTER                   :: subsys_section

      CHARACTER(LEN=default_string_length)               :: aname, bb
      INTEGER                                            :: i, io_unit, m, natorb, nshell
      INTEGER, DIMENSION(5)                              :: lval, nval, occupation
      LOGICAL                                            :: defined
      REAL(dp)                                           :: zeff
      REAL(KIND=dp)                                      :: alpha, en, eta, xgamma, zneff
      REAL(KIND=dp), DIMENSION(5)                        :: hen, kappa, kpoly, zeta
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (ASSOCIATED(xtb_parameter) .AND. &
          BTEST(cp_print_key_should_output(logger%iter_info, subsys_section, &
                                           "PRINT%KINDS/POTENTIAL"), cp_p_file)) THEN

         io_unit = cp_print_key_unit_nr(logger, subsys_section, "PRINT%KINDS", &
                                        extension=".Log")

         IF (io_unit > 0) THEN
            CALL get_xtb_atom_param(xtb_parameter, aname=aname, defined=defined, zeff=zeff, natorb=natorb)
            CALL get_xtb_atom_param(xtb_parameter, nshell=nshell, lval=lval, nval=nval, occupation=occupation)
            CALL get_xtb_atom_param(xtb_parameter, kpoly=kpoly, kappa=kappa, hen=hen, zeta=zeta)
            CALL get_xtb_atom_param(xtb_parameter, electronegativity=en, xgamma=xgamma, eta=eta, alpha=alpha, zneff=zneff)

            bb = "                                                "
            WRITE (UNIT=io_unit, FMT="(/,A,T67,A14)") " xTB  parameters: ", TRIM(aname)
            IF (defined) THEN
               m = 5 - nshell
               WRITE (UNIT=io_unit, FMT="(T16,A,T71,F10.2)") "Effective core charge:", zeff
               WRITE (UNIT=io_unit, FMT="(T16,A,T71,I10)") "Number of orbitals:", natorb
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5(A4,I1,I2,A1))") "Basis set [nl]", bb(1:8*m), &
                  ("   [", nval(i), lval(i), "]", i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5F8.3)") "Slater Exponent", bb(1:8*m), (zeta(i), i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5I8)") "Ref. occupation", bb(1:8*m), (occupation(i), i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5F8.3)") "Energy levels [au]", bb(1:8*m), (hen(i), i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5F8.3)") "Kpoly", bb(1:8*m), (kpoly(i), i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T71,F10.3)") "Electronegativity", en
               WRITE (UNIT=io_unit, FMT="(T16,A,T71,F10.3)") "Mataga-Nishimoto constant (eta)", eta
               WRITE (UNIT=io_unit, FMT="(T16,A,T41,A,5F8.3)") "Mataga-Nishimoto scaling kappa", bb(1:8*m), (kappa(i), i=1, nshell)
               WRITE (UNIT=io_unit, FMT="(T16,A,T71,F10.3)") "3rd Order constant", xgamma
               WRITE (UNIT=io_unit, FMT="(T16,A,T61,2F10.3)") "Repulsion potential [Z,alpha]", zneff, alpha
            ELSE
               WRITE (UNIT=io_unit, FMT="(T55,A)") "Parameters are not defined"
            END IF
         END IF
         CALL cp_print_key_finished_output(io_unit, logger, subsys_section, "PRINT%KINDS")
      END IF

   END SUBROUTINE write_xtb_atom_param

END MODULE xtb_types

