!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rotationally invariant parametrization of Fock matrix.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_linpot_rotinv
  USE ai_overlap,                      ONLY: overlap_aab
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                             orbtramat
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE spherical_harmonics,             ONLY: y_lm
  USE util,                            ONLY: sort
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_linpot_rotinv'

  PUBLIC :: linpot_rotinv_count_terms, linpot_rotinv_calc_term

CONTAINS

! *****************************************************************************
!> \brief Initialize rotinv parametrization
!> \param qs_env ...
!> \param iatom ...
!> \param nterms ...
! *****************************************************************************
  SUBROUTINE linpot_rotinv_count_terms(qs_env, iatom, nterms)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    INTEGER, INTENT(OUT)                     :: nterms

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_rotinv_count_terms', &
      routineP = moduleN//':'//routineN

    INTEGER :: dim_avail, ikind, iset, ishell, jneighbor, jset, jshell, li, &
      lj, lpot, natoms, pot_maxl, pot_neighbors
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

     CALL get_qs_env(qs_env,&
                     natom=natoms,&
                     qs_kind_set=qs_kind_set,&
                     particle_set=particle_set)

     CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
     CALL get_qs_kind(qs_kind_set(ikind),&
                      pao_potential_maxl=pot_maxl,&
                      pao_potential_neighbors=pot_neighbors,&
                      basis_set=basis_set)

     pot_neighbors = MIN(pot_neighbors, natoms-1)

     nterms = 0
     DO iset=1, basis_set%nset
     DO jset=1, iset
        DO ishell=1, basis_set%nshell(iset)
        DO jshell=1, ishell
           li = basis_set%l(ishell,iset)
           lj = basis_set%l(jshell,jset)

           ! calculate the block's number of dimensions
           IF(li==lj) THEN
              dim_avail = 2*li+1
           ELSE
              dim_avail = (2*li+1) * (2*lj+1)
           ENDIF

           inner:&
           DO lpot=0, pot_maxl
              ! basis function can only detect potentials of lower order
              IF(lpot > MIN(li, lj)) CYCLE
              DO jneighbor=1, pot_neighbors
                 nterms = nterms + 1
                 ! don't generate more terms than available dimensions
                 dim_avail = dim_avail - 1
                 IF(dim_avail==0) EXIT inner
              ENDDO
           ENDDO inner
        ENDDO
        ENDDO
     ENDDO
     ENDDO

  END SUBROUTINE linpot_rotinv_count_terms


! *****************************************************************************
!> \brief Count number of potential terms
!> \param qs_env ...
!> \param iatom ...
!> \param kterm ...
!> \param block_V ...
! *****************************************************************************
  SUBROUTINE linpot_rotinv_calc_term(qs_env, iatom, kterm, block_V)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom, kterm
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: block_V

    CHARACTER(len=*), PARAMETER :: routineN = 'linpot_rotinv_calc_term', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: inf_dist = 0.1*HUGE(1.0) 

    INTEGER :: c, dim_avail, i, idx_i, idx_j, ikind, iset, ishell, jatom, &
      jneighbor, jset, jshell, la1_max, la1_min, la2_max, la2_min, lb_max, &
      lb_min, li, lj, lpot, m, na1, na2, natoms, nb, ncfga1, ncfga2, ncfgb, &
      npgfa1, npgfa2, npgfb, nsgfa1, nsgfa2, pot_maxl, pot_neighbors
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: neighbor_order
    REAL(dp)                                 :: coeff
    REAL(dp), DIMENSION(:), POINTER          :: rpgfa1, rpgfa2, zeta1, zeta2
    REAL(dp), DIMENSION(:, :), POINTER       :: sab
    REAL(KIND=dp)                            :: c2s, pot_beta, pot_exp, &
                                                R_vdw, y
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: neighbor_dist
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rab_dir, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: saab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

     CALL get_qs_env(qs_env,&
                     natom=natoms,&
                     cell=cell,&
                     particle_set=particle_set,&
                     qs_kind_set=qs_kind_set)

     CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
     CALL get_qs_kind(qs_kind_set(ikind),&
                      pao_potential_maxl=pot_maxl,&
                      pao_potential_neighbors=pot_neighbors,&
                      pao_potential_beta=pot_beta,&
                      basis_set=basis_set)

     ! avoid overflows
     pot_neighbors = MIN(pot_neighbors, natoms-1)
     CALL init_orbital_pointers(2*pot_maxl)
     CALL init_spherical_harmonics(2*pot_maxl, output_unit=-1)

     c = 0 ! counter
     outer:&
     DO iset=1, basis_set%nset
     DO jset=1, iset
        DO ishell=1, basis_set%nshell(iset)
        DO jshell=1, ishell
           li = basis_set%l(ishell,iset)
           lj = basis_set%l(jshell,jset)
           IF(li==lj) THEN
              dim_avail = 2*li+1
           ELSE
              dim_avail = (2*li+1) * (2*lj+1)
           ENDIF

           inner:&
           DO lpot=0, pot_maxl
              IF(lpot > MIN(li, lj)) CYCLE
              DO jneighbor=1, pot_neighbors
                 c = c + 1
                 IF(c==kterm) EXIT outer
                 dim_avail = dim_avail - 1
                 IF(dim_avail==0) EXIT inner
              ENDDO
           ENDDO inner
        ENDDO
        ENDDO
     ENDDO
     ENDDO outer


     ! find neighbors
     ! TODO: this is a quadratic algorithm, use a neighbor-list instead
     ALLOCATE(neighbor_dist(natoms), neighbor_order(natoms))
     DO jatom=1, natoms
        ra = particle_set(iatom)%r
        rb = particle_set(jatom)%r
        rab = pbc(ra, rb, cell)
        neighbor_dist(jatom) = SQRT(SUM(rab*rab))
     ENDDO
     CALL sort(neighbor_dist, natoms, neighbor_order)
     CPASSERT(neighbor_order(1) == iatom) !central atom should be closesd to itself

     block_V(:,:) = 0.0_dp
     IF(jneighbor<1) THEN
        ! figure out where to store result
        idx_i = basis_set%first_sgf(ishell, iset)
        idx_j = basis_set%first_sgf(jshell, jset)
        block_V(idx_i, idx_j) = 1.0_dp
     ELSE
        !map jneighbor to jatom
        jatom = neighbor_order(jneighbor + 1)

        ! setup description of potential
        lb_min = lpot
        lb_max = lpot
        ncfgb  = ncoset(lb_max) - ncoset(lb_min-1)
        npgfb  = 1 ! number of exponents
        nb     = npgfb * ncfgb

        ! initialize exponents, using R_vdw from neighbour atom
        ALLOCATE(rpgfb(npgfb), zetb(npgfb))
        rpgfb(1) = inf_dist ! range of functions, we dont want screening
        CALL get_atomic_kind(particle_set(jatom)%atomic_kind, rvdw=R_vdw)
        pot_exp = 0.5_dp / (R_vdw*pot_beta)**2
        !pot_exp = 0.5_dp / pot_beta**2
        zetb(1) =  pot_exp !* 1.6_dp**(kpot-1)

        ! distance vector A-B
        ra = particle_set(iatom)%r
        rb = particle_set(jatom)%r
        rab = pbc(ra, rb, cell)

        ! setup iset
        la1_max = basis_set%l(ishell,iset)
        la1_min = basis_set%l(ishell,iset)
        npgfa1  =  basis_set%npgf(iset)
        ncfga1  =  ncoset(la1_max) - ncoset(la1_min-1)
        na1     =  npgfa1 * ncfga1
        zeta1   => basis_set%zet(:,iset)
        !rpgfa1  => basis_set%pgf_radius(:,iset) !TODO: reenable screening
        ALLOCATE(rpgfa1(npgfa1))
        rpgfa1(:) = inf_dist ! range of functions, we dont want screening


        ! setup jset
        la2_max = basis_set%l(jshell,jset)
        la2_min = basis_set%l(jshell,jset)
        npgfa2  =  basis_set%npgf(jset)
        ncfga2  =  ncoset(la2_max) - ncoset(la2_min-1)
        na2     =  npgfa2 * ncfga2
        zeta2   => basis_set%zet(:,jset)
        !rpgfa2  => basis_set%pgf_radius(:,jset) !TODO: reenable screening
        ALLOCATE(rpgfa2(npgfa2))
        rpgfa2(:) = inf_dist ! range of functions, we dont want screening


        ! allocate storage for integrals in cartesian basis
        ALLOCATE(saab(na1, na2, nb))
        saab = 0.0_dp

        CALL overlap_aab(la1_max, la1_min, npgfa1, rpgfa1, zeta1,&
                         la2_max, la2_min, npgfa2, rpgfa2, zeta2,&
                         lb_max,  lb_min,  npgfb,  rpgfb,  zetb,&
                         rab=rab, saab=saab)

        DEALLOCATE(rpgfb,zetb,rpgfa1,rpgfa2)

        ! number of spherical Gaussian functions in the set
        nsgfa1 = 2*basis_set%l(ishell,iset) + 1
        nsgfa2 = 2*basis_set%l(jshell,jset) + 1

        ! allocate storage for integrals in spherical basis
        ALLOCATE(sab(na1, na2))

        ! map zonal Harmonic in direction of atom_j onto carthesian Harmonics
        rab_dir = rab / SQRT(SUM(rab*rab))
        sab(:,:) = 0.0_dp
        DO i=1, ncfgb
           coeff = 0.0_dp
           DO m=-lpot, +lpot
              CALL Y_lm(rab_dir, y, lpot, m)
              ! with this normalization the y**2 actually sum to 1.
              y = y * SQRT(fourpi / (2*lpot+1))

              ! Empirics showed that one has to use %slm instead of %c2s...
              ! presumably this get's the norms right.
              c2s = orbtramat(lpot)%slm(lpot+m+1, i)

              coeff = coeff + y * c2s
           ENDDO
           sab = sab + coeff * saab(:,:,i)
        ENDDO

        CALL my_contract(sab, block_V, basis_set, iset, ishell, jset, jshell)
        DEALLOCATE(saab, sab)
    ENDIF

    block_V = block_V + TRANSPOSE(block_V)

  END SUBROUTINE linpot_rotinv_calc_term


! *****************************************************************************
!> \brief Helper routine, contracts a basis block
!> \param sab ...
!> \param block_V ...
!> \param basis_set ...
!> \param iset ...
!> \param ishell ...
!> \param jset ...
!> \param jshell ...
! *****************************************************************************
  SUBROUTINE my_contract(sab, block_V, basis_set, iset, ishell, jset, jshell)
    REAL(dp), DIMENSION(:, :), POINTER       :: sab
    REAL(dp), DIMENSION(:, :), INTENT(OUT), &
      TARGET                                 :: block_V
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    INTEGER, INTENT(IN)                      :: iset, ishell, jset, jshell

    INTEGER                                  :: a, b, c, d, ipgf, jpgf, l1, &
                                                l2, n1, n2, nn1, nn2, sgfa1, &
                                                sgfa2, sgla1, sgla2
    REAL(dp), DIMENSION(:, :), POINTER       :: S, T1, T2, V

! first and last indices of given shell in block_V.
! This matrix is in the contracted spherical basis.

     sgfa1 = basis_set%first_sgf(ishell,iset)
     sgla1 = basis_set%last_sgf(ishell,iset)
     sgfa2 = basis_set%first_sgf(jshell,jset)
     sgla2 = basis_set%last_sgf(jshell,jset)

     ! prepare the result matrix
     block_V(:,:) = 0.0_dp
     V => block_V(sgfa1:sgla1, sgfa2:sgla2)

     ! Calculate strides of sphi matrix.
     ! This matrix is in the uncontraced cartesian basis.
     ! It contains all shells of the set.
     ! It's index runs over all primitive gaussians of the set
     ! and then for each gaussian over all configurations of *the entire set*. (0->lmax)
     nn1   =  ncoset(basis_set%lmax(iset))
     nn2   =  ncoset(basis_set%lmax(jset))

     ! Calculate strides of sab matrix
     ! This matrix is also in the uncontraced cartensian basis,
     ! however it contains only a single shell.
     ! It's index runs over all primitive gaussians of the set
     ! and then for each gaussian over all configrations of *the given shell*.
     l1 = basis_set%l(ishell,iset)
     l2 = basis_set%l(jshell,jset)
     n1 =  ncoset(l1) - ncoset(l1-1)
     n2 =  ncoset(l2) - ncoset(l2-1)

     DO ipgf = 1, basis_set%npgf(iset)
     DO jpgf = 1, basis_set%npgf(jset)
        ! prepare first trafo-matrix
        a = (ipgf-1) * nn1 + ncoset(l1 - 1) + 1
        T1 => basis_set%sphi(a:a+n1-1 , sgfa1:sgla1)

        ! prepare second trafo-matrix
        b = (jpgf-1) * nn2 + ncoset(l2 - 1) + 1
        T2 => basis_set%sphi(b:b+n2-1 , sgfa2:sgla2)

        ! prepare SAB matrix
        c = (ipgf-1)*n1 + 1
        d = (jpgf-1)*n2 + 1
        S => sab(c:c+n1-1, d:d+n2-1)

        ! do the transformation
        V = V + MATMUL(TRANSPOSE(T1),MATMUL(S,T2))
     ENDDO
     ENDDO

  END SUBROUTINE my_contract

END MODULE pao_linpot_rotinv
