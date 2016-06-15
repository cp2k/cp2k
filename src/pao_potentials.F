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
   USE ao_util,                         ONLY: exp_radius
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gamma1
   USE mathlib,                         ONLY: multinomial
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_potentials'

   PUBLIC :: calc_gaussian

CONTAINS

! **************************************************************************************************
!> \brief Calculates potential term of the form r**lpot * Exp(-beta*r**2)
!> \param basis_set ...
!> \param block_V potential term that is returned
!> \param Rab ...
!> \param lpot polynomial prefactor, r**lpot
!> \param beta exponent of the Gaussian
!> \param weight ...
!> \param min_shell ...
!> \param max_shell ...
!> \param min_l ...
!> \param max_l ...
! **************************************************************************************************
   SUBROUTINE calc_gaussian(basis_set, block_V, Rab, lpot, beta, weight, min_shell, max_shell, min_l, max_l)
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V
      REAL(dp), DIMENSION(3)                             :: Rab
      INTEGER, INTENT(IN)                                :: lpot
      REAL(dp), INTENT(IN)                               :: beta, weight
      INTEGER, INTENT(IN), OPTIONAL                      :: min_shell, max_shell, min_l, max_l

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_gaussian'

      INTEGER :: handle, ic, iset, ishell, ishell_abs, jset, jshell, jshell_abs, la1_max, la1_min, &
         la2_max, la2_min, lb_max, lb_min, N, na1, na2, nb, ncfga1, ncfga2, ncfgb, npgfa1, npgfa2, &
         npgfb
      REAL(dp)                                           :: coeff, norm2
      REAL(dp), DIMENSION(:), POINTER                    :: rpgfa1, rpgfa2, rpgfb, zeta1, zeta2, zetb
      REAL(dp), DIMENSION(:, :), POINTER                 :: new_block, sab
      REAL(dp), DIMENSION(:, :, :), POINTER              :: saab

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(min_shell) .EQV. PRESENT(max_shell))
      CPASSERT(PRESENT(min_l) .EQV. PRESENT(max_l))

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

      ! initialize exponents
      ALLOCATE (rpgfb(npgfb), zetb(npgfb))
      rpgfb(1) = exp_radius(0, beta, 1.0E-12_dp, 1.0_dp) ! TODO get the EPS parameter from somewhere / precompute this elsewhere
      zetb(1) = beta

      ! loop over all set/shell combination and fill block_V
      DO iset = 1, basis_set%nset
      DO jset = 1, basis_set%nset
      DO ishell = 1, basis_set%nshell(iset)
      DO jshell = 1, basis_set%nshell(jset)
         IF (PRESENT(min_shell) .AND. PRESENT(max_shell)) THEN
            ishell_abs = SUM(basis_set%nshell(1:iset-1))+ishell
            jshell_abs = SUM(basis_set%nshell(1:jset-1))+jshell
            IF (MIN(ishell_abs, jshell_abs) /= min_shell) CYCLE
            IF (MAX(ishell_abs, jshell_abs) /= max_shell) CYCLE
         ENDIF
         IF (PRESENT(min_l) .AND. PRESENT(min_l)) THEN
            IF (MIN(basis_set%l(ishell, iset), basis_set%l(jshell, jset)) /= min_l) CYCLE
            IF (MAX(basis_set%l(ishell, iset), basis_set%l(jshell, jset)) /= max_l) CYCLE
         ENDIF

         ! setup iset
         la1_max = basis_set%l(ishell, iset)
         la1_min = basis_set%l(ishell, iset)
         npgfa1 = basis_set%npgf(iset)
         ncfga1 = ncoset(la1_max)-ncoset(la1_min-1)
         na1 = npgfa1*ncfga1
         zeta1 => basis_set%zet(:, iset)
         rpgfa1 => basis_set%pgf_radius(:, iset)

         ! setup jset
         la2_max = basis_set%l(jshell, jset)
         la2_min = basis_set%l(jshell, jset)
         npgfa2 = basis_set%npgf(jset)
         ncfga2 = ncoset(la2_max)-ncoset(la2_min-1)
         na2 = npgfa2*ncfga2
         zeta2 => basis_set%zet(:, jset)
         rpgfa2 => basis_set%pgf_radius(:, jset)

         ! calculate integrals
         ALLOCATE (saab(na1, na2, nb))
         saab = 0.0_dp
         CALL overlap_aab(la1_max=la1_max, la1_min=la1_min, npgfa1=npgfa1, rpgfa1=rpgfa1, zeta1=zeta1, &
                          la2_max=la2_max, la2_min=la2_min, npgfa2=npgfa2, rpgfa2=rpgfa2, zeta2=zeta2, &
                          lb_max=lb_max, lb_min=lb_min, npgfb=npgfb, rpgfb=rpgfb, zetb=zetb, &
                          rab=Rab, saab=saab)

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

      norm2 = (2.0_dp*beta)**(-0.5_dp-lpot)*gamma1(lpot)

      block_V = block_V+weight*new_block/SQRT(norm2)
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
