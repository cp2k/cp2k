!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rotationally invariant parametrization of Fock matrix.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_linpot_rotinv
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE pao_potentials,                  ONLY: calc_gaussian
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_potential_type,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_linpot_rotinv'

   PUBLIC :: linpot_rotinv_count_terms, linpot_rotinv_calc_term

CONTAINS

! **************************************************************************************************
!> \brief Count number of terms for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nterms ...
! **************************************************************************************************
   SUBROUTINE linpot_rotinv_count_terms(qs_env, ikind, nterms)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nterms

      INTEGER                                            :: ipot, iset, ishell, ishell_abs, lmax, &
                                                            lmin, lpot, max_shell, min_shell, &
                                                            npots, nshells, pot_maxl
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: shell_l
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_potentials=pao_potentials)

      nshells = SUM(basis_set%nshell)
      npots = SIZE(pao_potentials)

      IF (npots == 0) CPWARN("Found no PAO_POTENTIAL section")

      ! fill shell_l
      ALLOCATE (shell_l(nshells))
      DO iset = 1, basis_set%nset
      DO ishell = 1, basis_set%nshell(iset)
         ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
         shell_l(ishell_abs) = basis_set%l(ishell, iset)
      ENDDO
      ENDDO

      nterms = 0

      ! terms sensing neighboring atoms
      DO ipot = 1, npots
         pot_maxl = pao_potentials(ipot)%maxl ! maxl is taken from central atom
         IF (pot_maxl < 0) &
            CPABORT("ROTINV parametrization requires non-negative PAO_POTENTIAL%MAXL")
         IF (MOD(pot_maxl, 2) /= 0) &
            CPABORT("ROTINV parametrization requires even-numbered PAO_POTENTIAL%MAXL")
         DO max_shell = 1, nshells
         DO min_shell = 1, max_shell
         DO lpot = 0, pot_maxl, 2
            lmin = shell_l(min_shell)
            lmax = shell_l(max_shell)
            IF (lmin == 0 .AND. lmax == 0) CYCLE ! coverted by central terms
            nterms = nterms+1
         ENDDO
         ENDDO
         ENDDO
      ENDDO

      ! spherical symmetric terms on central atom
      DO max_shell = 1, nshells
      DO min_shell = 1, max_shell
         IF (shell_l(min_shell) /= shell_l(max_shell)) CYCLE ! need quadratic block
         nterms = nterms+1
      ENDDO
      ENDDO

   END SUBROUTINE linpot_rotinv_count_terms

! **************************************************************************************************
!> \brief Fills the given block_V with the requested potential term
!> \param qs_env ...
!> \param iatom ...
!> \param kterm ...
!> \param block_V ...
! **************************************************************************************************
   SUBROUTINE linpot_rotinv_calc_term(qs_env, iatom, kterm, block_V)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom, kterm
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V

      INTEGER :: c, i, ikind, ipot, iset, ishell, ishell_abs, jatom, jkind, jset, jshell, &
         jshell_abs, lmax, lmin, lpot, max_shell, min_shell, N, natoms, npots, nshells, pot_maxl, &
         sgfa1, sgfa2, sgla1, sgla2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: shell_l
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      REAL(KIND=dp)                                      :: pot_beta, pot_weight
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: ipao_potentials, jpao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      cell=cell, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_potentials=ipao_potentials)
      npots = SIZE(ipao_potentials)
      nshells = SUM(basis_set%nshell)

      ! fill shell_l
      ALLOCATE (shell_l(nshells))
      DO iset = 1, basis_set%nset
      DO ishell = 1, basis_set%nshell(iset)
         ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
         shell_l(ishell_abs) = basis_set%l(ishell, iset)
      ENDDO
      ENDDO

      N = basis_set%nsgf ! primary basis-size
      CPASSERT(SIZE(block_V, 1) == N .AND. SIZE(block_V, 2) == N)
      block_V = 0.0_dp

      c = 0

      ! terms sensing neighboring atoms
      outer1: DO ipot = 1, npots
         pot_maxl = ipao_potentials(ipot)%maxl ! taken from central atom
         DO max_shell = 1, nshells
         DO min_shell = 1, max_shell
         DO lpot = 0, pot_maxl, 2
            lmin = shell_l(min_shell)
            lmax = shell_l(max_shell)
            IF (lmin == 0 .AND. lmax == 0) CYCLE ! coverted by central terms
            c = c+1
            IF (c == kterm) EXIT outer1
         ENDDO
         ENDDO
         ENDDO
      ENDDO outer1

      IF (c == kterm) THEN
         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE ! no self-interaction
            CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(qs_kind_set(jkind), pao_potentials=jpao_potentials)
            IF (SIZE(jpao_potentials) /= npots) &
               CPABORT("Not all KINDs have the same number of PAO_POTENTIAL sections")
            pot_beta = jpao_potentials(ipot)%beta ! taken from remote atom
            pot_weight = jpao_potentials(ipot)%weight ! taken from remote atom
            Ra = particle_set(iatom)%r
            Rb = particle_set(jatom)%r
            Rab = pbc(ra, rb, cell)
            CALL calc_gaussian(basis_set, block_V, Rab, lpot=lpot, &
                               beta=pot_beta, weight=pot_weight, &
                               min_shell=min_shell, max_shell=max_shell)
         ENDDO

      ELSE !---------------------------------------------------------------------------------------

         ! spherical symmetric terms on central atom
         outer2: DO max_shell = 1, nshells
            DO min_shell = 1, max_shell
               IF (shell_l(min_shell) /= shell_l(max_shell)) CYCLE ! need quadratic block
               c = c+1
               IF (c == kterm) EXIT outer2
            END DO
         ENDDO outer2
         CPASSERT(c == kterm)

         DO iset = 1, basis_set%nset
         DO jset = 1, basis_set%nset
         DO ishell = 1, basis_set%nshell(iset)
         DO jshell = 1, basis_set%nshell(jset)
            ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
            jshell_abs = SUM(basis_set%nshell(1:jset-1))+jshell
            IF (MIN(ishell_abs, jshell_abs) /= min_shell) CYCLE
            IF (MAX(ishell_abs, jshell_abs) /= max_shell) CYCLE
            sgfa1 = basis_set%first_sgf(ishell, iset)
            sgla1 = basis_set%last_sgf(ishell, iset)
            sgfa2 = basis_set%first_sgf(jshell, jset)
            sgla2 = basis_set%last_sgf(jshell, jset)
            CPASSERT((sgla1-sgfa1) == (sgla2-sgfa2)) ! should be a quadratic block
            ! set diagonal
            DO i = 1, sgla1-sgfa1+1
               block_V(sgfa1-1+i, sgfa2-1+i) = 1.0_dp
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         block_V = block_V/SQRT(SUM(block_V**2)) ! normalize

      ENDIF
   END SUBROUTINE linpot_rotinv_calc_term

END MODULE pao_linpot_rotinv
