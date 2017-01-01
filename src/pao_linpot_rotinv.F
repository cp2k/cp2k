!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rotationally invariant parametrization of Fock matrix.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_linpot_rotinv
   USE ai_overlap,                      ONLY: overlap_aab
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gamma1
   USE mathlib,                         ONLY: multinomial
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
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

   PUBLIC :: linpot_rotinv_count_terms, linpot_rotinv_calc_terms, linpot_rotinv_calc_forces

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

      CHARACTER(len=*), PARAMETER :: routineN = 'linpot_rotinv_count_terms'

      INTEGER                                            :: handle, ipot, iset, ishell, ishell_abs, &
                                                            lmax, lmin, lpot, max_shell, &
                                                            min_shell, npots, nshells, pot_maxl
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: shell_l
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

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

      CALL timestop(handle)

   END SUBROUTINE linpot_rotinv_count_terms

! **************************************************************************************************
!> \brief Calculate all potential terms of the rotinv parametrization
!> \param qs_env ...
!> \param iatom ...
!> \param V_blocks ...
! **************************************************************************************************
   SUBROUTINE linpot_rotinv_calc_terms(qs_env, iatom, V_blocks)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :, :), INTENT(OUT), TARGET  :: V_blocks

      CHARACTER(len=*), PARAMETER :: routineN = 'linpot_rotinv_calc_terms'

      INTEGER :: handle, i, ic, ikind, ipot, iset, ishell, ishell_abs, jatom, jkind, jset, jshell, &
         jshell_abs, kterm, la1_max, la1_min, la2_max, la2_min, lb_max, lb_min, lpot, N, na1, na2, &
         natoms, nb, ncfga1, ncfga2, ncfgb, npgfa1, npgfa2, npgfb, npots, pot_maxl, sgfa1, sgfa2, &
         sgla1, sgla2
      REAL(dp)                                           :: coeff, norm2, pot_beta, pot_weight, &
                                                            rpgfa_max, tab
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      REAL(dp), DIMENSION(:), POINTER                    :: rpgfa1, rpgfa2, rpgfb, zeta1, zeta2, zetb
      REAL(dp), DIMENSION(:, :), POINTER                 :: T1, T2, V12, V21
      REAL(dp), DIMENSION(:, :, :), POINTER              :: block_V_full, saab, saal
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: ipao_potentials, jpao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      cell=cell, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_potentials=ipao_potentials)
      npots = SIZE(ipao_potentials)
      N = basis_set%nsgf ! primary basis-size
      CPASSERT(SIZE(V_blocks, 1) == N .AND. SIZE(V_blocks, 2) == N)
      kterm = 0 ! init counter

      DO ipot = 1, npots
         pot_maxl = ipao_potentials(ipot)%maxl ! taken from central atom

         ! setup description of potential
         lb_min = 0
         lb_max = pot_maxl
         ncfgb = ncoset(lb_max)-ncoset(lb_min-1)
         npgfb = 1 ! number of exponents
         nb = npgfb*ncfgb
         ALLOCATE (rpgfb(npgfb), zetb(npgfb))

         ! build block_V_full
         ALLOCATE (block_V_full(N, N, pot_maxl/2+1))
         block_V_full = 0.0_dp

         DO iset = 1, basis_set%nset
         DO jset = 1, iset

            ! setup iset
            la1_max = basis_set%lmax(iset)
            la1_min = basis_set%lmin(iset)
            npgfa1 = basis_set%npgf(iset)
            ncfga1 = ncoset(la1_max)-ncoset(la1_min-1)
            na1 = npgfa1*ncfga1
            zeta1 => basis_set%zet(:, iset)
            rpgfa1 => basis_set%pgf_radius(:, iset)

            ! setup jset
            la2_max = basis_set%lmax(jset)
            la2_min = basis_set%lmin(jset)
            npgfa2 = basis_set%npgf(jset)
            ncfga2 = ncoset(la2_max)-ncoset(la2_min-1)
            na2 = npgfa2*ncfga2
            zeta2 => basis_set%zet(:, jset)
            rpgfa2 => basis_set%pgf_radius(:, jset)

            ! radius of most diffuse basis-function
            rpgfa_max = MAX(MAXVAL(rpgfa1), MAXVAL(rpgfa2))

            ! allocate space for integrals
            ALLOCATE (saab(na1, na2, nb), saal(na1, na2, pot_maxl/2+1))
            saal = 0.0_dp

            ! loop over neighbors
            DO jatom = 1, natoms
               IF (jatom == iatom) CYCLE ! no self-interaction
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
               CALL get_qs_kind(qs_kind_set(jkind), pao_potentials=jpao_potentials)
               IF (SIZE(jpao_potentials) /= npots) &
                  CPABORT("Not all KINDs have the same number of PAO_POTENTIAL sections")

               ! initialize exponents
               pot_weight = jpao_potentials(ipot)%weight ! taken from remote atom
               pot_beta = jpao_potentials(ipot)%beta ! taken from remote atom
               rpgfb(1) = jpao_potentials(ipot)%beta_radius ! taken from remote atom
               zetb(1) = pot_beta

               ! calculate direction
               Ra = particle_set(iatom)%r
               Rb = particle_set(jatom)%r
               Rab = pbc(ra, rb, cell)

               ! distance screening
               tab = SQRT(SUM(Rab**2))
               IF (rpgfa_max+rpgfb(1) < tab) CYCLE

               ! calculate actual integrals
               saab = 0.0_dp
               CALL overlap_aab(la1_max=la1_max, la1_min=la1_min, npgfa1=npgfa1, rpgfa1=rpgfa1, zeta1=zeta1, &
                                la2_max=la2_max, la2_min=la2_min, npgfa2=npgfa2, rpgfa2=rpgfa2, zeta2=zeta2, &
                                lb_max=lb_max, lb_min=lb_min, npgfb=npgfb, rpgfb=rpgfb, zetb=zetb, &
                                rab=Rab, saab=saab)

               ! sum neighbor contributions according to remote atom's weight and normalization
               DO lpot = 0, pot_maxl, 2
                  norm2 = (2.0_dp*pot_beta)**(-0.5_dp-lpot)*gamma1(lpot)
                  ! sum potential terms: POW(x**2 + y**2 + z**2, lpot/2)
                  DO ic = ncoset(lpot-1)+1, ncoset(lpot)
                     coeff = multinomial(lpot/2, indco(:, ic)/2)
                     saal(:, :, lpot/2+1) = saal(:, :, lpot/2+1)+saab(:, :, ic)*coeff*pot_weight/SQRT(norm2)
                  ENDDO
               ENDDO
            ENDDO ! jatom

            ! find bounds of set-pair and setup transformation matrices
            sgfa1 = basis_set%first_sgf(1, iset)
            sgla1 = sgfa1+basis_set%nsgf_set(iset)-1
            sgfa2 = basis_set%first_sgf(1, jset)
            sgla2 = sgfa2+basis_set%nsgf_set(jset)-1
            T1 => basis_set%scon(1:na1, sgfa1:sgla1)
            T2 => basis_set%scon(1:na2, sgfa2:sgla2)

            ! transform into primary basis
            DO lpot = 0, pot_maxl, 2
               V12 => block_V_full(sgfa1:sgla1, sgfa2:sgla2, lpot/2+1)
               V21 => block_V_full(sgfa2:sgla2, sgfa1:sgla1, lpot/2+1)
               V12 = MATMUL(TRANSPOSE(T1), MATMUL(saal(:, :, lpot/2+1), T2))
               V21 = TRANSPOSE(V12)
            ENDDO
            DEALLOCATE (saab, saal)
         ENDDO ! jset
         ENDDO ! iset
         DEALLOCATE (rpgfb, zetb)

         ! block_V_full is ready -------------------------------------------------------------------
         ! split the full blocks into shell-pair sub-blocks
         DO iset = 1, basis_set%nset
         DO jset = 1, iset
         DO ishell = 1, basis_set%nshell(iset)
         DO jshell = 1, basis_set%nshell(jset)
            IF (basis_set%l(ishell, iset) == 0 .AND. basis_set%l(jshell, jset) == 0) CYCLE ! covered by central terms
            ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
            jshell_abs = SUM(basis_set%nshell(1:jset-1))+jshell
            IF (ishell_abs < jshell_abs) CYCLE

            ! find bounds of shell-pair
            sgfa1 = basis_set%first_sgf(ishell, iset)
            sgla1 = basis_set%last_sgf(ishell, iset)
            sgfa2 = basis_set%first_sgf(jshell, jset)
            sgla2 = basis_set%last_sgf(jshell, jset)

            DO lpot = 0, pot_maxl, 2
               kterm = kterm+1
               V_blocks(:, :, kterm) = 0.0_dp
               V_blocks(sgfa1:sgla1, sgfa2:sgla2, kterm) = block_V_full(sgfa1:sgla1, sgfa2:sgla2, lpot/2+1)
               V_blocks(sgfa2:sgla2, sgfa1:sgla1, kterm) = block_V_full(sgfa2:sgla2, sgfa1:sgla1, lpot/2+1)
            ENDDO ! lpot
         ENDDO ! jshell
         ENDDO ! ishell
         ENDDO ! jset
         ENDDO ! iset
         DEALLOCATE (block_V_full)
      ENDDO ! ipot

      ! terms on central atom ----------------------------------------------------------------------

      DO iset = 1, basis_set%nset
      DO jset = 1, iset
      DO ishell = 1, basis_set%nshell(iset)
      DO jshell = 1, basis_set%nshell(jset)
         IF (basis_set%l(ishell, iset) /= basis_set%l(jshell, jset)) CYCLE ! need quadratic block
         ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
         jshell_abs = SUM(basis_set%nshell(1:jset-1))+jshell
         IF (ishell_abs < jshell_abs) CYCLE
         kterm = kterm+1
         sgfa1 = basis_set%first_sgf(ishell, iset)
         sgla1 = basis_set%last_sgf(ishell, iset)
         sgfa2 = basis_set%first_sgf(jshell, jset)
         sgla2 = basis_set%last_sgf(jshell, jset)
         CPASSERT((sgla1-sgfa1) == (sgla2-sgfa2)) ! should be a quadratic block
         V_blocks(:, :, kterm) = 0.0_dp
         DO i = 1, sgla1-sgfa1+1 ! set diagonal of sub-block
            V_blocks(sgfa1-1+i, sgfa2-1+i, kterm) = 1.0_dp
            V_blocks(sgfa2-1+i, sgfa1-1+i, kterm) = 1.0_dp
         ENDDO
         norm2 = SUM(V_blocks(:, :, kterm)**2)
         V_blocks(:, :, kterm) = V_blocks(:, :, kterm)/SQRT(norm2) ! normalize
      ENDDO ! jshell
      ENDDO ! ishell
      ENDDO ! jset
      ENDDO ! iset

      CPASSERT(SIZE(V_blocks, 3) == kterm) ! ensure we generated all terms

      CALL timestop(handle)
   END SUBROUTINE linpot_rotinv_calc_terms

! **************************************************************************************************
!> \brief Calculate force contribution from rotinv parametrization
!> \param qs_env ...
!> \param iatom ...
!> \param M_blocks ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE linpot_rotinv_calc_forces(qs_env, iatom, M_blocks, forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: M_blocks
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'linpot_rotinv_calc_forces'

      INTEGER :: handle, i, ic, ikind, ipot, iset, ishell, ishell_abs, jatom, jkind, jset, jshell, &
         jshell_abs, kterm, la1_max, la1_min, la2_max, la2_min, lb_max, lb_min, lpot, N, na1, na2, &
         natoms, nb, ncfga1, ncfga2, ncfgb, npgfa1, npgfa2, npgfb, npots, nshells, pot_maxl, &
         sgfa1, sgfa2, sgla1, sgla2
      REAL(dp)                                           :: coeff, f, norm2, pot_beta, pot_weight, &
                                                            rpgfa_max, tab
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      REAL(dp), DIMENSION(:), POINTER                    :: rpgfa1, rpgfa2, rpgfb, zeta1, zeta2, zetb
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_D, T1, T2
      REAL(dp), DIMENSION(:, :, :), POINTER              :: block_M_full, dab
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: daab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: ipao_potentials, jpao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      cell=cell, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_potentials=ipao_potentials)
      npots = SIZE(ipao_potentials)
      nshells = SUM(basis_set%nshell)
      N = basis_set%nsgf ! primary basis-size
      CPASSERT(SIZE(M_blocks, 1) == N .AND. SIZE(M_blocks, 2) == N)
      kterm = 0 ! init counter
      ALLOCATE (block_D(N, N))

      DO ipot = 1, npots
         pot_maxl = ipao_potentials(ipot)%maxl ! taken from central atom

         ! build block_M_full
         ALLOCATE (block_M_full(N, N, pot_maxl/2+1))
         block_M_full = 0.0_dp
         DO iset = 1, basis_set%nset
         DO jset = 1, iset
            DO ishell = 1, basis_set%nshell(iset)
            DO jshell = 1, basis_set%nshell(jset)
               IF (basis_set%l(ishell, iset) == 0 .AND. basis_set%l(jshell, jset) == 0) CYCLE ! covered by central terms
               ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
               jshell_abs = SUM(basis_set%nshell(1:jset-1))+jshell
               IF (ishell_abs < jshell_abs) CYCLE
               ! find bounds of shell-pair
               sgfa1 = basis_set%first_sgf(ishell, iset)
               sgla1 = basis_set%last_sgf(ishell, iset)
               sgfa2 = basis_set%first_sgf(jshell, jset)
               sgla2 = basis_set%last_sgf(jshell, jset)
               DO lpot = 0, pot_maxl, 2
                  kterm = kterm+1
                  block_M_full(sgfa1:sgla1, sgfa2:sgla2, lpot/2+1) = M_blocks(sgfa1:sgla1, sgfa2:sgla2, kterm)
                  block_M_full(sgfa2:sgla2, sgfa1:sgla1, lpot/2+1) = M_blocks(sgfa2:sgla2, sgfa1:sgla1, kterm)
               ENDDO ! lpot
            ENDDO ! jshell
            ENDDO ! ishell
         ENDDO ! jset
         ENDDO ! iset

         ! setup description of potential
         lb_min = 0
         lb_max = pot_maxl
         ncfgb = ncoset(lb_max)-ncoset(lb_min-1)
         npgfb = 1 ! number of exponents
         nb = npgfb*ncfgb
         ALLOCATE (rpgfb(npgfb), zetb(npgfb))

         DO iset = 1, basis_set%nset
         DO jset = 1, iset

            ! setup iset
            la1_max = basis_set%lmax(iset)
            la1_min = basis_set%lmin(iset)
            npgfa1 = basis_set%npgf(iset)
            ncfga1 = ncoset(la1_max)-ncoset(la1_min-1)
            na1 = npgfa1*ncfga1
            zeta1 => basis_set%zet(:, iset)
            rpgfa1 => basis_set%pgf_radius(:, iset)

            ! setup jset
            la2_max = basis_set%lmax(jset)
            la2_min = basis_set%lmin(jset)
            npgfa2 = basis_set%npgf(jset)
            ncfga2 = ncoset(la2_max)-ncoset(la2_min-1)
            na2 = npgfa2*ncfga2
            zeta2 => basis_set%zet(:, jset)
            rpgfa2 => basis_set%pgf_radius(:, jset)

            ! radius of most diffuse basis-function
            rpgfa_max = MAX(MAXVAL(rpgfa1), MAXVAL(rpgfa2))

            ! find bounds of set-pair and setup transformation matrices
            sgfa1 = basis_set%first_sgf(1, iset)
            sgla1 = sgfa1+basis_set%nsgf_set(iset)-1
            sgfa2 = basis_set%first_sgf(1, jset)
            sgla2 = sgfa2+basis_set%nsgf_set(jset)-1
            T1 => basis_set%scon(1:na1, sgfa1:sgla1)
            T2 => basis_set%scon(1:na2, sgfa2:sgla2)

            ! allocate space for integrals
            ALLOCATE (daab(na1, na2, nb, 3), dab(na1, na2, 3))

            ! loop over neighbors
            DO jatom = 1, natoms
               IF (jatom == iatom) CYCLE ! no self-interaction
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
               CALL get_qs_kind(qs_kind_set(jkind), pao_potentials=jpao_potentials)
               IF (SIZE(jpao_potentials) /= npots) &
                  CPABORT("Not all KINDs have the same number of PAO_POTENTIAL sections")

               ! initialize exponents
               pot_weight = jpao_potentials(ipot)%weight ! taken from remote atom
               pot_beta = jpao_potentials(ipot)%beta ! taken from remote atom
               rpgfb(1) = jpao_potentials(ipot)%beta_radius ! taken from remote atom
               zetb(1) = pot_beta

               ! calculate direction
               Ra = particle_set(iatom)%r
               Rb = particle_set(jatom)%r
               Rab = pbc(ra, rb, cell)

               ! distance screening
               tab = SQRT(SUM(Rab**2))
               IF (rpgfa_max+rpgfb(1) < tab) CYCLE

               ! calculate actual integrals
               daab = 0.0_dp
               CALL overlap_aab(la1_max=la1_max, la1_min=la1_min, npgfa1=npgfa1, rpgfa1=rpgfa1, zeta1=zeta1, &
                                la2_max=la2_max, la2_min=la2_min, npgfa2=npgfa2, rpgfa2=rpgfa2, zeta2=zeta2, &
                                lb_max=lb_max, lb_min=lb_min, npgfb=npgfb, rpgfb=rpgfb, zetb=zetb, &
                                rab=Rab, daab=daab)

               ! sum neighbor contributions according to remote atom's weight and normalization
               DO lpot = 0, pot_maxl, 2
                  ! sum potential terms: POW(x**2 + y**2 + z**2, lpot/2)
                  dab = 0.0_dp
                  DO ic = ncoset(lpot-1)+1, ncoset(lpot)
                     norm2 = (2.0_dp*pot_beta)**(-0.5_dp-lpot)*gamma1(lpot)
                     coeff = multinomial(lpot/2, indco(:, ic)/2)
                     dab = dab+coeff*daab(:, :, ic, :)*pot_weight/SQRT(norm2)
                  ENDDO
                  DO i = 1, 3
                     ! transform into primary basis
                     block_D = 0.0_dp
                     block_D(sgfa1:sgla1, sgfa2:sgla2) = MATMUL(TRANSPOSE(T1), MATMUL(dab(:, :, i), T2))
                     block_D(sgfa2:sgla2, sgfa1:sgla1) = TRANSPOSE(block_D(sgfa1:sgla1, sgfa2:sgla2))
                     ! calculate and add forces
                     f = SUM(block_M_full(:, :, lpot/2+1)*block_D)
                     forces(iatom, i) = forces(iatom, i)-f
                     forces(jatom, i) = forces(jatom, i)+f
                  ENDDO
               ENDDO ! lpot
            ENDDO ! jatom
            DEALLOCATE (dab, daab)
         ENDDO ! jset
         ENDDO ! iset
         DEALLOCATE (rpgfb, zetb, block_M_full)
      ENDDO ! ipot
      DEALLOCATE (block_D)

      CALL timestop(handle)
   END SUBROUTINE linpot_rotinv_calc_forces

END MODULE pao_linpot_rotinv
