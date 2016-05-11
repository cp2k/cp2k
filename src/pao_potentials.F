!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Factory routines for potentials used e.g. by pao_param_exp and pao_ml
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_potentials
   USE ai_overlap,                      ONLY: overlap_aab
   USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: multinomial
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_potentials'

   PUBLIC :: pao_build_V_neighbors, pao_build_V_core

CONTAINS

! **************************************************************************************************
!> \brief Construct a potential which mimics the chemical environment created by neighboring atoms
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param block_V ...
! **************************************************************************************************
   SUBROUTINE pao_build_V_neighbors(particle_set, qs_kind_set, cell, iatom, block_V)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: block_V

      INTEGER                                            :: ikind, jatom, jkind, N, natoms
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      REAL(KIND=dp)                                      :: pot_beta, pot_weight
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CPASSERT(.NOT. ALLOCATED(block_V))

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
      N = basis_set%nsgf
      natoms = SIZE(particle_set)

      ALLOCATE (block_V(N, N))
      block_V = 0.0_dp
      DO jatom = 1, natoms
         IF (jatom == iatom) CYCLE
         CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
         CALL get_qs_kind(qs_kind_set(jkind), pao_pot_beta=pot_beta, pao_pot_weight=pot_weight)
         Ra = particle_set(iatom)%r
         Rb = particle_set(jatom)%r
         Rab = pbc(ra, rb, cell)
         CALL calc_gaussian(basis_set, block_V, Rab, lpot=0, beta=pot_beta, weight=pot_weight)
      ENDDO

   END SUBROUTINE pao_build_V_neighbors

! **************************************************************************************************
!> \brief Construct a potential based on the atomic core Hamiltonian
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param iatom ...
!> \param block_V ...
! **************************************************************************************************
   SUBROUTINE pao_build_V_core(particle_set, qs_kind_set, iatom, block_V)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: block_V

      INTEGER                                            :: ikind, N
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_ptr

      CPASSERT(.NOT. ALLOCATED(block_V))

      NULLIFY (block_V_ptr)
      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL calculate_atomic_orbitals(atomic_kind=particle_set(iatom)%atomic_kind, &
                                     qs_kind=qs_kind_set(ikind), &
                                     iunit=-1, &
                                     fmat=block_V_ptr, &
                                     ispin=1)

      ! turn pointer into allocatable
      N = SIZE(block_V_ptr, 1)
      ALLOCATE (block_V(N, N))
      block_V(:, :) = block_V_ptr
      DEALLOCATE (block_V_ptr)
   END SUBROUTINE pao_build_V_core

! **************************************************************************************************
!> \brief Calculates potential term of the form r**lpot * Exp(-beta*r**2)
!> \param basis_set ...
!> \param block_V potential term that is returned
!> \param Rab ...
!> \param lpot polynomial prefactor, r**lpot
!> \param beta exponent of the Gaussian
!> \param weight ...
! **************************************************************************************************
   SUBROUTINE calc_gaussian(basis_set, block_V, Rab, lpot, beta, weight)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V
      REAL(dp), DIMENSION(3)                             :: Rab
      INTEGER, INTENT(IN)                                :: lpot
      REAL(dp), INTENT(IN)                               :: beta, weight

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_gaussian'
      REAL(dp), PARAMETER                                :: inf_dist = 0.1*HUGE(1.0)

      INTEGER :: handle, ic, iset, ishell, jset, jshell, la1_max, la1_min, la2_max, la2_min, &
         lb_max, lb_min, N, na1, na2, nb, ncfga1, ncfga2, ncfgb, npgfa1, npgfa2, npgfb
      REAL(dp)                                           :: coeff
      REAL(dp), DIMENSION(:), POINTER                    :: rpgfa1, rpgfa2, rpgfb, zeta1, zeta2, zetb
      REAL(dp), DIMENSION(:, :), POINTER                 :: new_block, sab
      REAL(dp), DIMENSION(:, :, :), POINTER              :: saab

      CALL timeset(routineN, handle)

      CPASSERT(MOD(lpot, 2) == 0) ! otherwise it's not rotationally invariant
      CALL init_orbital_pointers(lpot)

      N = SIZE(block_V, 1)
      ALLOCATE (new_block(N, N))
      new_block = 0.0_dp

      ! setup description of potential
      lb_min = lpot
      lb_max = lpot
      ncfgb = ncoset(lb_max)-ncoset(lb_min-1)
      npgfb = 1 ! number of exponents
      nb = npgfb*ncfgb

      ! initialize exponents, using R_vdw from neighbour atom
      ALLOCATE (rpgfb(npgfb), zetb(npgfb))
      rpgfb(1) = inf_dist ! range of functions, we dont want screening
      zetb(1) = beta

      ! loop over all set/shell combination and fill block_V
      DO iset = 1, basis_set%nset
      DO jset = 1, basis_set%nset
      DO ishell = 1, basis_set%nshell(iset)
      DO jshell = 1, basis_set%nshell(jset)

         ! setup iset
         la1_max = basis_set%l(ishell, iset)
         la1_min = basis_set%l(ishell, iset)
         npgfa1 = basis_set%npgf(iset)
         ncfga1 = ncoset(la1_max)-ncoset(la1_min-1)
         na1 = npgfa1*ncfga1
         zeta1 => basis_set%zet(:, iset)
         !rpgfa1  => basis_set%pgf_radius(:,iset) !TODO: reenable screening
         ALLOCATE (rpgfa1(npgfa1))
         rpgfa1(:) = inf_dist ! range of functions, we dont want screening

         ! setup jset
         la2_max = basis_set%l(jshell, jset)
         la2_min = basis_set%l(jshell, jset)
         npgfa2 = basis_set%npgf(jset)
         ncfga2 = ncoset(la2_max)-ncoset(la2_min-1)
         na2 = npgfa2*ncfga2
         zeta2 => basis_set%zet(:, jset)
         !rpgfa2  => basis_set%pgf_radius(:,jset) !TODO: reenable screening
         ALLOCATE (rpgfa2(npgfa2))
         rpgfa2(:) = inf_dist ! range of functions, we dont want screening

         ! calculate integrals
         ALLOCATE (saab(na1, na2, nb))
         saab = 0.0_dp
         CALL overlap_aab(la1_max=la1_max, la1_min=la1_min, npgfa1=npgfa1, rpgfa1=rpgfa1, zeta1=zeta1, &
                          la2_max=la2_max, la2_min=la2_min, npgfa2=npgfa2, rpgfa2=rpgfa2, zeta2=zeta2, &
                          lb_max=lb_max, lb_min=lb_min, npgfb=npgfb, rpgfb=rpgfb, zetb=zetb, &
                          rab=Rab, saab=saab)
         DEALLOCATE (rpgfa1, rpgfa2)

         ! sum potential terms: POW(x**2 + y**2 + z**2, lpot/2)
         ALLOCATE (sab(na1, na2))
         sab = 0.0_dp
         DO ic = 1, ncfgb
            coeff = multinomial(lpot/2, indco(:, ncoset(lpot-1)+ic)/2)
            sab = sab+coeff*saab(:, :, ic)
         ENDDO

         CALL my_contract(sab=sab, block_V=new_block, basis_set=basis_set, &
                          iset=iset, ishell=ishell, jset=jset, jshell=jshell)
         DEALLOCATE (sab, saab)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DEALLOCATE (rpgfb, zetb)

      block_V = block_V+weight*new_block
      DEALLOCATE (new_block)

      ! symmetrize
      block_V = 0.5_dp*(block_V+TRANSPOSE(block_V))

      CALL timestop(handle)
   END SUBROUTINE calc_gaussian

! **************************************************************************************************
!> \brief Helper routine, contracts a basis block
!> \param sab ...
!> \param block_V ...
!> \param basis_set ...
!> \param iset ...
!> \param ishell ...
!> \param jset ...
!> \param jshell ...
! **************************************************************************************************
   SUBROUTINE my_contract(sab, block_V, basis_set, iset, ishell, jset, jshell)
      REAL(dp), DIMENSION(:, :), POINTER                 :: sab
      REAL(dp), DIMENSION(:, :), INTENT(OUT), TARGET     :: block_V
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(IN)                                :: iset, ishell, jset, jshell

      INTEGER                                            :: a, b, c, d, ipgf, jpgf, l1, l2, n1, n2, &
                                                            nn1, nn2, sgfa1, sgfa2, sgla1, sgla2
      REAL(dp), DIMENSION(:, :), POINTER                 :: S, T1, T2, V

      ! first and last indices of given shell in block_V.
      ! This matrix is in the contracted spherical basis.
      sgfa1 = basis_set%first_sgf(ishell, iset)
      sgla1 = basis_set%last_sgf(ishell, iset)
      sgfa2 = basis_set%first_sgf(jshell, jset)
      sgla2 = basis_set%last_sgf(jshell, jset)

      ! prepare the result matrix
      V => block_V(sgfa1:sgla1, sgfa2:sgla2)

      ! Calculate strides of sphi matrix.
      ! This matrix is in the uncontraced cartesian basis.
      ! It contains all shells of the set.
      ! It's index runs over all primitive gaussians of the set
      ! and then for each gaussian over all configurations of *the entire set*. (0->lmax)
      nn1 = ncoset(basis_set%lmax(iset))
      nn2 = ncoset(basis_set%lmax(jset))

      ! Calculate strides of sab matrix
      ! This matrix is also in the uncontraced cartensian basis,
      ! however it contains only a single shell.
      ! It's index runs over all primitive gaussians of the set
      ! and then for each gaussian over all configrations of *the given shell*.
      l1 = basis_set%l(ishell, iset)
      l2 = basis_set%l(jshell, jset)
      n1 = ncoset(l1)-ncoset(l1-1)
      n2 = ncoset(l2)-ncoset(l2-1)

      DO ipgf = 1, basis_set%npgf(iset)
      DO jpgf = 1, basis_set%npgf(jset)
         ! prepare first trafo-matrix
         a = (ipgf-1)*nn1+ncoset(l1-1)+1
         T1 => basis_set%sphi(a:a+n1-1, sgfa1:sgla1)

         ! prepare second trafo-matrix
         b = (jpgf-1)*nn2+ncoset(l2-1)+1
         T2 => basis_set%sphi(b:b+n2-1, sgfa2:sgla2)

         ! prepare SAB matrix
         c = (ipgf-1)*n1+1
         d = (jpgf-1)*n2+1
         S => sab(c:c+n1-1, d:d+n2-1)

         ! do the transformation
         V = V+MATMUL(TRANSPOSE(T1), MATMUL(S, T2))
      ENDDO
      ENDDO

   END SUBROUTINE my_contract

END MODULE pao_potentials
