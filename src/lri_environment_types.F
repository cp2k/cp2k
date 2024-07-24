!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains the types and subroutines for dealing with the lri_env
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) debugged
!> \authors JGH
!>          Dorothea Golze
! **************************************************************************************************
MODULE lri_environment_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: INT_8,&
                                              dp,&
                                              sp
   USE mathlib,                         ONLY: pswitch
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_o3c_types,                    ONLY: o3c_container_type,&
                                              release_o3c_container
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
! Integral container
   TYPE carray
      INTEGER                                                 :: compression = -1
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: cdp => NULL()
      REAL(KIND=sp), DIMENSION(:), POINTER                    :: csp => NULL()
      INTEGER(INT_8), DIMENSION(:), POINTER                   :: cip => NULL()
   END TYPE carray

   TYPE int_container
      INTEGER                                                 :: na = -1, nb = -1, nc = -1
      TYPE(carray), DIMENSION(:), POINTER                     :: ca => NULL()
   END TYPE int_container
! **************************************************************************************************
   TYPE lri_rhoab_type
      ! number of spherical basis functions (a)
      INTEGER                                                 :: nba = -1
      ! number of spherical basis functions (b)
      INTEGER                                                 :: nbb = -1
      ! number of spherical fit basis functions (ai)
      INTEGER                                                 :: nfa = -1
      ! number of spherical fit basis functions (bi)
      INTEGER                                                 :: nfb = -1
      ! expansion coeffs for RI density
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: avec => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: aveca => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: avecb => NULL()
      ! projection coeffs for RI density: SUM_ab (ab,i)*Pab
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvec => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tveca => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvecb => NULL()
      ! integral (ai) * sinv * tvec
      REAL(KIND=dp)                                           :: nst = 0.0_dp
      REAL(KIND=dp)                                           :: nsta = 0.0_dp
      REAL(KIND=dp)                                           :: nstb = 0.0_dp
      ! Lagrange parameter
      REAL(KIND=dp)                                           :: lambda = 0.0_dp
      REAL(KIND=dp)                                           :: lambdaa = 0.0_dp
      REAL(KIND=dp)                                           :: lambdab = 0.0_dp
      ! Charge of pair density
      REAL(KIND=dp)                                           :: charge = 0.0_dp
      REAL(KIND=dp)                                           :: chargea = 0.0_dp
      REAL(KIND=dp)                                           :: chargeb = 0.0_dp
   END TYPE lri_rhoab_type

! **************************************************************************************************

   TYPE lri_int_type
      ! whether to calculate force for pair
      LOGICAL                                                 :: calc_force_pair = .FALSE.
      ! number of spherical basis functions (a)
      INTEGER                                                 :: nba = -1
      ! number of spherical basis functions (b)
      INTEGER                                                 :: nbb = -1
      ! number of spherical fit basis functions (ai)
      INTEGER                                                 :: nfa = -1
      ! number of spherical fit basis functions (bi)
      INTEGER                                                 :: nfb = -1
      ! condition number of overlap matrix
      REAL(KIND=dp)                                           :: cond_num = 0.0_dp
      ! integrals (a,b,ai)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abaint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                :: abascr
      ! integrals (a,b,b)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abbint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                :: abbscr
      ! compressed aba integrals
      TYPE(int_container)                                     :: cabai = int_container()
      ! compressed abb integrals
      TYPE(int_container)                                     :: cabbi = int_container()
      ! integrals (da/dA,b,dai/dA)
      REAL(KIND=dp), DIMENSION(:, :, :, :), ALLOCATABLE       :: dabdaint
      ! integrals (da/dA,b,bi)
      REAL(KIND=dp), DIMENSION(:, :, :, :), ALLOCATABLE       :: dabbint
      ! integrals (a,b)
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE             :: soo
      ! derivative d(a,b)/dA
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: dsoo
      ! integrals (ai,bi)
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE             :: sab
      ! derivative d(ai,bi)/dA
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: dsab
      ! inverse of integrals (ai,bi)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: sinv => NULL()
      ! integral (ai) / (bi), dim(1..nfa,nfa+1..nfa+nfb)
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: n => NULL()
      ! sinv * (ai)
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: sn => NULL()
      ! (ai) * sinv * (ai)
      REAL(KIND=dp)                                           :: nsn = 0.0_dp
      ! distant pair approximation
      LOGICAL                                                 :: lrisr = .FALSE.
      LOGICAL                                                 :: lriff = .FALSE.
      REAL(KIND=dp)                                           :: wsr = 0.0_dp, wff = 0.0_dp, dwsr = 0.0_dp, dwff = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: asinv => NULL(), bsinv => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: na => NULL(), nb => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: sna => NULL(), snb => NULL()
      REAL(KIND=dp)                                           :: nsna = 0.0_dp, nsnb = 0.0_dp
      !
      ! dmax: max deviation for integrals of primitive gtos; for debugging
      ! dmax for overlap integrals (ai,bi); fit bas
      REAL(KIND=dp)                                           :: dmax_ab = 0.0_dp
      ! dmax for overlap integrals (a,b); orb bas
      REAL(KIND=dp)                                           :: dmax_oo = 0.0_dp
      ! dmax for integrals (a,b,ai)
      REAL(KIND=dp)                                           :: dmax_aba = 0.0_dp
      ! dmax for integrals (a,b,bi)
      REAL(KIND=dp)                                           :: dmax_abb = 0.0_dp
   END TYPE lri_int_type

   TYPE lri_int_rho_type
      ! integrals (aa,bb), orb basis
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER           :: soaabb => NULL()
      ! dmax for (aa,bb) integrals; for debugging
      REAL(KIND=dp)                                           :: dmax_aabb = 0.0_dp
   END TYPE lri_int_rho_type

   TYPE lri_node_type
      INTEGER                                                 :: nnode = 0
      TYPE(lri_int_type), DIMENSION(:), POINTER               :: lri_int => NULL()
      TYPE(lri_int_rho_type), DIMENSION(:), POINTER           :: lri_int_rho => NULL()
      TYPE(lri_rhoab_type), DIMENSION(:), POINTER             :: lri_rhoab => NULL()
   END TYPE lri_node_type

   TYPE lri_atom_type
      INTEGER                                                 :: natom = 0
      TYPE(lri_node_type), DIMENSION(:), POINTER              :: lri_node => NULL()
   END TYPE lri_atom_type

   TYPE lri_list_type
      INTEGER                                                 :: nkind = 0
      TYPE(lri_atom_type), DIMENSION(:), POINTER              :: lri_atom => NULL()
   END TYPE lri_list_type

   TYPE lri_list_p_type
      TYPE(lri_list_type), POINTER                             :: lri_list => NULL()
   END TYPE lri_list_p_type

! **************************************************************************************************

   TYPE lri_bas_type
      INTEGER, DIMENSION(:, :, :), POINTER                    :: orb_index => NULL()
      INTEGER, DIMENSION(:, :, :), POINTER                    :: ri_index => NULL()
      ! integral of ri basis fbas
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: int_fbas => NULL()
      ! self overlap ri basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: ri_ovlp => NULL()
      ! inverse of self overlap ri basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: ri_ovlp_inv => NULL()
      ! self overlap orb basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: orb_ovlp => NULL()
      ! self overlap (a,a,fa)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: ovlp3
      ! contraction matrix for SHG integrals ri basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER              :: scon_ri => NULL()
      ! contraction matrix for SHG integrals orb basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER              :: scon_orb => NULL()
      ! contraction matrix for SHG integrals aba/abb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER           :: scon_mix => NULL()
   END TYPE lri_bas_type

! **************************************************************************************************

   TYPE lri_clebsch_gordon_type
      ! Clebsch-Gordon (CG) coefficients
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER                :: cg_coeff => NULL()
      ! list of non-zero CG coefficients
      INTEGER, DIMENSION(:, :, :), POINTER                      :: cg_none0_list => NULL()
      ! number of non-zero CG coefficients
      INTEGER, DIMENSION(:, :), POINTER                         :: ncg_none0 => NULL()
   END TYPE lri_clebsch_gordon_type

! **************************************************************************************************

   TYPE lri_ppl_type
      ! integrals Vppl*fbas (potential*fit basis) dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_int => NULL()
   END TYPE lri_ppl_type

   TYPE lri_ppl_int_type
      TYPE(lri_ppl_type), DIMENSION(:), POINTER                 :: lri_ppl => NULL()
      REAL(KIND=dp)                                             :: ecore_pp_ri = 0.0_dp
   END TYPE lri_ppl_int_type

! **************************************************************************************************

   TYPE ri_fit_type
      INTEGER, DIMENSION(:, :), POINTER                         :: bas_ptr => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: nvec => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: rm1n => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: tvec => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: rm1t => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: avec => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: fout => NULL()
      REAL(KIND=dp)                                             :: ntrm1n = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                               :: ftrm1n = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                               :: echarge = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                               :: lambda = 0.0_dp
   END TYPE ri_fit_type

! **************************************************************************************************
   TYPE wmat_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: mat => NULL()
   END TYPE wmat_type

   TYPE wbas_type
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: vec => NULL()
   END TYPE wbas_type
! **************************************************************************************************
   TYPE stat_type
      REAL(KIND=dp)                                             :: pairs_tt = 0.0_dp
      REAL(KIND=dp)                                             :: pairs_sr = 0.0_dp
      REAL(KIND=dp)                                             :: pairs_ff = 0.0_dp
      REAL(KIND=dp)                                             :: overlap_error = 0.0_dp
      REAL(KIND=dp)                                             :: rho_tt = 0.0_dp
      REAL(KIND=dp)                                             :: rho_sr = 0.0_dp
      REAL(KIND=dp)                                             :: rho_ff = 0.0_dp
      REAL(KIND=dp)                                             :: rho_1c = 0.0_dp
      REAL(KIND=dp)                                             :: coef_mem = 0.0_dp
      REAL(KIND=dp)                                             :: oint_mem = 0.0_dp
      REAL(KIND=dp)                                             :: rhos_mem = 0.0_dp
      REAL(KIND=dp)                                             :: abai_mem = 0.0_dp
      REAL(KIND=dp)                                             :: ppli_mem = 0.0_dp
   END TYPE stat_type
! **************************************************************************************************

   TYPE lri_environment_type
      ! parameter for (pseudo)inverse of overlap
      INTEGER                                                 :: lri_overlap_inv = -1
      ! flag for debugging lri integrals
      LOGICAL                                                 :: debug = .FALSE.
      ! flag for shg (solid haromonic Gaussian) integrals
      LOGICAL                                                 :: use_shg_integrals = .FALSE.
      ! parameter for inversion (autoselect); maximal condition
      ! number up to where inversion is legal
      REAL(KIND=dp)                                           :: cond_max = 0.0_dp
      ! parameter for checking distance between atom pairs
      REAL(KIND=dp)                                           :: delta = 0.0_dp
      ! threshold for aba and abb integrals
      REAL(KIND=dp)                                           :: eps_o3_int = 0.0_dp
      ! orbital basis set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: orb_basis => NULL()
      ! lri (fit) basis set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis => NULL()
      ! orb_basis neighborlist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soo_list => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: saa_list => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soa_list => NULL()
      ! local RI integrals
      TYPE(lri_list_type), POINTER                            :: lri_ints => NULL()
      ! local Vppl integrals
      TYPE(lri_ppl_int_type), POINTER                         :: lri_ppl_ints => NULL()
      ! local integral of rho**2; for optimization
      TYPE(lri_list_type), POINTER                            :: lri_ints_rho => NULL()
      ! properties of orb and aux basis
      TYPE(lri_bas_type), DIMENSION(:), POINTER               :: bas_prop => NULL()
      ! Clebsch-Gordon for solid harmonics
      TYPE(lri_clebsch_gordon_type), POINTER                  :: cg_shg => NULL()
      ! orbital basis overlap
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ob_smat => NULL()
      ! statistics
      LOGICAL                                                 :: statistics = .FALSE.
      TYPE(stat_type)                                         :: stat = stat_type()
      ! exact one-center terms
      LOGICAL                                                 :: exact_1c_terms = .FALSE.
      ! use RI for local pp
      LOGICAL                                                 :: ppl_ri = .FALSE.
      ! store integrals (needed for basis optimization)
      LOGICAL                                                 :: store_integrals = .FALSE.
      ! distant pair approximation
      LOGICAL                                                 :: distant_pair_approximation = .FALSE.
      CHARACTER(len=10)                                       :: distant_pair_method = ""
      REAL(KIND=dp)                                           :: r_in = 0.0_dp
      REAL(KIND=dp)                                           :: r_out = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: aradius => NULL()
      TYPE(wbas_type), DIMENSION(:), POINTER                  :: wbas => NULL()
      TYPE(wmat_type), DIMENSION(:, :), POINTER               :: wmat => NULL()
      ! RI overlap and inverse
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ri_smat => NULL(), &
                                                                 ri_sinv => NULL()
      TYPE(ri_fit_type), POINTER                              :: ri_fit => NULL()
      CHARACTER(len=10)                                       :: ri_sinv_app = ""
      TYPE(o3c_container_type), POINTER                       :: o3c => NULL()
   END TYPE lri_environment_type

! **************************************************************************************************

   TYPE lri_kind_type
      ! expansion coeff for lri density dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: acoef => NULL()
      ! integrals V*fbas (potential*fit basis) dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_int => NULL()
      ! SUM_i integral(V*fbas_i)*davec/dR dim(natom,3)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_dadr => NULL()
      ! integrals V*dfbas/dR
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_dfdr => NULL()
   END TYPE lri_kind_type

   TYPE lri_spin_type
      TYPE(lri_kind_type), DIMENSION(:), POINTER                :: lri_kinds => NULL()
   END TYPE lri_spin_type

! **************************************************************************************************

   TYPE lri_force_type
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: st => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: dssn => NULL(), &
                                                                   dsst => NULL()
      ! derivative dtvec/dR
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: dtvec => NULL()
   END TYPE lri_force_type

! **************************************************************************************************

   TYPE lri_density_type
      INTEGER                                                 :: nspin = 0
      ! pair density expansion (nspin)
      TYPE(lri_list_p_type), DIMENSION(:), POINTER            :: lri_rhos => NULL()
      ! coefficients of RI expansion and gradients (nspin)
      TYPE(lri_spin_type), DIMENSION(:), POINTER              :: lri_coefs => NULL()
   END TYPE lri_density_type

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_types'

   PUBLIC :: lri_environment_type, &
             lri_force_type, lri_list_type, &
             lri_int_type, lri_int_rho_type, lri_density_type, &
             lri_kind_type, lri_rhoab_type
   PUBLIC :: int_container, carray
   PUBLIC :: lri_env_create, lri_env_release, allocate_lri_coefs, &
             allocate_lri_ints, allocate_lri_ints_rho, lri_density_create, &
             allocate_lri_ppl_ints, deallocate_lri_ppl_ints, &
             lri_density_release, allocate_lri_rhos, allocate_lri_force_components, &
             deallocate_lri_ints, deallocate_lri_ints_rho, &
             deallocate_lri_force_components, deallocate_bas_properties

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
! **************************************************************************************************
   SUBROUTINE lri_env_create(lri_env)

      TYPE(lri_environment_type), INTENT(OUT)            :: lri_env

      lri_env%debug = .FALSE.
      lri_env%delta = 1.E-6_dp

      lri_env%store_integrals = .FALSE.

      NULLIFY (lri_env%orb_basis)
      NULLIFY (lri_env%ri_basis)

      NULLIFY (lri_env%soo_list)
      NULLIFY (lri_env%saa_list)
      NULLIFY (lri_env%soa_list)
      NULLIFY (lri_env%lri_ints)
      NULLIFY (lri_env%lri_ppl_ints)
      NULLIFY (lri_env%lri_ints_rho)
      NULLIFY (lri_env%bas_prop)

      NULLIFY (lri_env%ob_smat)
      NULLIFY (lri_env%ri_smat)
      NULLIFY (lri_env%ri_sinv)
      NULLIFY (lri_env%ri_fit)
      NULLIFY (lri_env%o3c)
      NULLIFY (lri_env%aradius)
      NULLIFY (lri_env%wmat)
      NULLIFY (lri_env%wbas)

      NULLIFY (lri_env%cg_shg)
      ALLOCATE (lri_env%cg_shg)
      NULLIFY (lri_env%cg_shg%cg_coeff)
      NULLIFY (lri_env%cg_shg%cg_none0_list)
      NULLIFY (lri_env%cg_shg%ncg_none0)

   END SUBROUTINE lri_env_create

! **************************************************************************************************
!> \brief releases the given lri_env
!> \param lri_env the lri environment to release
! **************************************************************************************************
   SUBROUTINE lri_env_release(lri_env)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env

      INTEGER                                            :: i, ikind, j, nkind

      ! deallocate basis sets
      IF (ASSOCIATED(lri_env%orb_basis)) THEN
         nkind = SIZE(lri_env%orb_basis)
         DO ikind = 1, nkind
            CALL deallocate_gto_basis_set(lri_env%orb_basis(ikind)%gto_basis_set)
         END DO
         DEALLOCATE (lri_env%orb_basis)
      END IF
      IF (ASSOCIATED(lri_env%ri_basis)) THEN
         nkind = SIZE(lri_env%ri_basis)
         DO ikind = 1, nkind
            CALL deallocate_gto_basis_set(lri_env%ri_basis(ikind)%gto_basis_set)
         END DO
         DEALLOCATE (lri_env%ri_basis)
      END IF
      CALL release_neighbor_list_sets(lri_env%soo_list)
      CALL release_neighbor_list_sets(lri_env%saa_list)
      CALL release_neighbor_list_sets(lri_env%soa_list)
      IF (ASSOCIATED(lri_env%lri_ints)) THEN
         CALL deallocate_lri_ints(lri_env%lri_ints)
         DEALLOCATE (lri_env%lri_ints)
         NULLIFY (lri_env%lri_ints)
      END IF
      IF (ASSOCIATED(lri_env%lri_ppl_ints)) THEN
         CALL deallocate_lri_ppl_ints(lri_env%lri_ppl_ints)
         DEALLOCATE (lri_env%lri_ppl_ints)
         NULLIFY (lri_env%lri_ppl_ints)
      END IF
      IF (ASSOCIATED(lri_env%lri_ints_rho)) THEN
         CALL deallocate_lri_ints_rho(lri_env%lri_ints_rho)
         DEALLOCATE (lri_env%lri_ints_rho)
         NULLIFY (lri_env%lri_ints_rho)
      END IF
      CALL deallocate_bas_properties(lri_env)
      IF (ASSOCIATED(lri_env%aradius)) THEN
         DEALLOCATE (lri_env%aradius)
      END IF
      IF (ASSOCIATED(lri_env%wmat)) THEN
         DO i = 1, SIZE(lri_env%wmat, 1)
            DO j = 1, SIZE(lri_env%wmat, 2)
               IF (ASSOCIATED(lri_env%wmat(i, j)%mat)) THEN
                  DEALLOCATE (lri_env%wmat(i, j)%mat)
               END IF
            END DO
         END DO
         DEALLOCATE (lri_env%wmat)
      END IF
      IF (ASSOCIATED(lri_env%wbas)) THEN
         DO i = 1, SIZE(lri_env%wbas, 1)
            IF (ASSOCIATED(lri_env%wbas(i)%vec)) THEN
               DEALLOCATE (lri_env%wbas(i)%vec)
            END IF
         END DO
         DEALLOCATE (lri_env%wbas)
      END IF
      IF (ASSOCIATED(lri_env%cg_shg)) THEN
         IF (ASSOCIATED(lri_env%cg_shg%cg_coeff)) THEN
            DEALLOCATE (lri_env%cg_shg%cg_coeff)
         END IF
         IF (ASSOCIATED(lri_env%cg_shg%cg_none0_list)) THEN
            DEALLOCATE (lri_env%cg_shg%cg_none0_list)
         END IF
         IF (ASSOCIATED(lri_env%cg_shg%ncg_none0)) THEN
            DEALLOCATE (lri_env%cg_shg%ncg_none0)
         END IF
         DEALLOCATE (lri_env%cg_shg)
      END IF
      ! RI
      IF (ASSOCIATED(lri_env%ob_smat)) CALL dbcsr_deallocate_matrix_set(lri_env%ob_smat)
      IF (ASSOCIATED(lri_env%ri_smat)) CALL dbcsr_deallocate_matrix_set(lri_env%ri_smat)
      IF (ASSOCIATED(lri_env%ri_sinv)) CALL dbcsr_deallocate_matrix_set(lri_env%ri_sinv)
      IF (ASSOCIATED(lri_env%ri_fit)) THEN
         IF (ASSOCIATED(lri_env%ri_fit%nvec)) THEN
            DEALLOCATE (lri_env%ri_fit%nvec)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%rm1n)) THEN
            DEALLOCATE (lri_env%ri_fit%rm1n)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%tvec)) THEN
            DEALLOCATE (lri_env%ri_fit%tvec)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%rm1t)) THEN
            DEALLOCATE (lri_env%ri_fit%rm1t)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%avec)) THEN
            DEALLOCATE (lri_env%ri_fit%avec)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%fout)) THEN
            DEALLOCATE (lri_env%ri_fit%fout)
         END IF
         IF (ASSOCIATED(lri_env%ri_fit%bas_ptr)) THEN
            DEALLOCATE (lri_env%ri_fit%bas_ptr)
         END IF
         DEALLOCATE (lri_env%ri_fit)
      END IF
      IF (ASSOCIATED(lri_env%o3c)) THEN
         CALL release_o3c_container(lri_env%o3c)
         DEALLOCATE (lri_env%o3c)
      END IF

   END SUBROUTINE lri_env_release

! **************************************************************************************************
!> \brief creates and initializes an lri_density environment
!> \param lri_density the lri_density environment you want to create
! **************************************************************************************************
   SUBROUTINE lri_density_create(lri_density)

      TYPE(lri_density_type), INTENT(OUT)                :: lri_density

      lri_density%nspin = 0

      NULLIFY (lri_density%lri_rhos)
      NULLIFY (lri_density%lri_coefs)

   END SUBROUTINE lri_density_create

! **************************************************************************************************
!> \brief releases the given lri_density
!> \param lri_density the lri_density to release
! **************************************************************************************************
   SUBROUTINE lri_density_release(lri_density)

      TYPE(lri_density_type), INTENT(INOUT)              :: lri_density

      CALL deallocate_lri_rhos(lri_density%lri_rhos)
      CALL deallocate_lri_coefs(lri_density%lri_coefs)

   END SUBROUTINE lri_density_release

! **************************************************************************************************
!> \brief allocate lri_ints, matrices that store LRI integrals
!> \param lri_env ...
!> \param lri_ints structure storing the LRI integrals
!> \param nkind number of atom kinds
! **************************************************************************************************
   SUBROUTINE allocate_lri_ints(lri_env, lri_ints, nkind)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env
      TYPE(lri_list_type), INTENT(OUT)                   :: lri_ints
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: i, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, nba, nbb, nfa, nfb, &
                                                            nlist, nn, nneighbor
      LOGICAL                                            :: dpa, e1c
      REAL(KIND=dp)                                      :: dab, ra, rab(3), rb
      TYPE(gto_basis_set_type), POINTER                  :: fbasa, fbasb, obasa, obasb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (fbasa, fbasb, nl_iterator, obasa, obasb)

      dpa = lri_env%distant_pair_approximation
      ra = lri_env%r_in
      rb = lri_env%r_out
      lri_env%stat%oint_mem = 0.0_dp

      lri_ints%nkind = nkind
      ALLOCATE (lri_ints%lri_atom(nkind*nkind))

      DO i = 1, nkind*nkind
         NULLIFY (lri_ints%lri_atom(i)%lri_node)
         lri_ints%lri_atom(i)%natom = 0
      END DO

      CALL neighbor_list_iterator_create(nl_iterator, lri_env%soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                iatom=iatom, jatom=jatom, r=rab)

         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set
         fbasa => lri_env%ri_basis(ikind)%gto_basis_set
         fbasb => lri_env%ri_basis(jkind)%gto_basis_set

         IF (.NOT. ASSOCIATED(obasa)) CYCLE
         IF (.NOT. ASSOCIATED(obasb)) CYCLE

         IF (.NOT. ASSOCIATED(lri_ints%lri_atom(iac)%lri_node)) THEN
            lri_ints%lri_atom(iac)%natom = nlist
            ALLOCATE (lri_ints%lri_atom(iac)%lri_node(nlist))
            DO i = 1, nlist
               NULLIFY (lri_ints%lri_atom(iac)%lri_node(i)%lri_int)
               lri_ints%lri_atom(iac)%lri_node(i)%nnode = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int)) THEN
            lri_ints%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
            ALLOCATE (lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(nneighbor))
         END IF

         ASSOCIATE (lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor))

            nba = obasa%nsgf
            nbb = obasb%nsgf
            nfa = fbasa%nsgf
            nfb = fbasb%nsgf
            nn = nfa + nfb

            IF (iatom == jatom .AND. dab < lri_env%delta) THEN
               e1c = lri_env%exact_1c_terms
            ELSE
               e1c = .FALSE.
            END IF

            IF (.NOT. e1c) THEN
               ALLOCATE (lrii%abascr(nfa))
               lrii%abascr = 0._dp
               ALLOCATE (lrii%abbscr(nfb))
               lrii%abbscr = 0._dp
               lri_env%stat%oint_mem = lri_env%stat%oint_mem + nfa + nfb
            END IF

            IF (dpa) THEN
               lrii%wsr = pswitch(dab, ra, rb, 0)
               lrii%wff = 1.0_dp - lrii%wsr
               lrii%dwsr = pswitch(dab, ra, rb, 1)
               lrii%dwff = -lrii%dwsr
               lrii%lrisr = (lrii%wsr > 0.0_dp)
               lrii%lriff = (lrii%wff > 0.0_dp)
               NULLIFY (lrii%asinv, lrii%bsinv)
            ELSE
               lrii%lrisr = .TRUE.
               lrii%lriff = .FALSE.
               lrii%wsr = 1.0_dp
               lrii%wff = 0.0_dp
               lrii%dwsr = 0.0_dp
               lrii%dwff = 0.0_dp
               NULLIFY (lrii%asinv, lrii%bsinv)
            END IF

            ! compressed storage
            NULLIFY (lrii%cabai%ca, lrii%cabbi%ca)

            ! full LRI method term
            IF (lrii%lrisr) THEN
               IF (e1c) THEN
                  NULLIFY (lrii%n, lrii%sn)
                  NULLIFY (lrii%sinv)
               ELSE
                  ALLOCATE (lrii%soo(nba, nbb))
                  lri_env%stat%oint_mem = lri_env%stat%oint_mem + nba*nbb
                  lrii%soo = 0._dp

                  IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                     ALLOCATE (lrii%sinv(nfa, nfa))
                     lri_env%stat%oint_mem = lri_env%stat%oint_mem + nfa*nfa
                  ELSE
                     ALLOCATE (lrii%sinv(nn, nn))
                     lri_env%stat%oint_mem = lri_env%stat%oint_mem + nn*nn
                  END IF
                  lrii%sinv = 0._dp

                  IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                     ALLOCATE (lrii%n(nfa), lrii%sn(nfa))
                     lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nfa
                  ELSE
                     ALLOCATE (lrii%n(nn), lrii%sn(nn))
                     lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nn
                  END IF
                  lrii%n = 0._dp
                  lrii%sn = 0._dp
               END IF
            ELSE
               NULLIFY (lrii%n, lrii%sn)
               NULLIFY (lrii%sinv)
            END IF

            ! far field approximation
            IF (lrii%lriff) THEN
               lrii%asinv => lri_env%bas_prop(ikind)%ri_ovlp_inv
               lrii%bsinv => lri_env%bas_prop(jkind)%ri_ovlp_inv
               ALLOCATE (lrii%na(nfa), lrii%sna(nfa))
               lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nfa
               lrii%na = 0._dp
               lrii%sna = 0._dp
               ALLOCATE (lrii%nb(nfb), lrii%snb(nfb))
               lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nfb
               lrii%nb = 0._dp
               lrii%snb = 0._dp
               IF (.NOT. ALLOCATED(lrii%soo)) THEN
                  ALLOCATE (lrii%soo(nba, nbb))
                  lri_env%stat%oint_mem = lri_env%stat%oint_mem + nba*nbb
                  lrii%soo = 0._dp
               ELSE
                  CPASSERT(SIZE(lrii%soo, 1) == nba .AND. SIZE(lrii%soo, 2) == nbb)
               END IF
            ELSE
               NULLIFY (lrii%na, lrii%sna)
               NULLIFY (lrii%nb, lrii%snb)
            END IF

            lrii%dmax_ab = 0._dp
            lrii%dmax_oo = 0._dp
            lrii%dmax_aba = 0._dp
            lrii%dmax_abb = 0._dp

            lrii%calc_force_pair = .FALSE.

         END ASSOCIATE

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE allocate_lri_ints

! **************************************************************************************************
!> \brief allocate lri_ppl_ints, matrices that store LRI integrals
!> \param lri_env ...
!> \param lri_ppl_ints structure storing the LRI ppl integrals
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_ppl_ints(lri_env, lri_ppl_ints, atomic_kind_set)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env
      TYPE(lri_ppl_int_type), INTENT(OUT)                :: lri_ppl_ints
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ikind, natom, nfa, nkind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gto_basis_set_type), POINTER                  :: fbasa

      lri_env%stat%ppli_mem = 0.0_dp
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (lri_ppl_ints%lri_ppl(nkind))
      DO ikind = 1, nkind
         fbasa => lri_env%ri_basis(ikind)%gto_basis_set
         nfa = fbasa%nsgf
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
         ALLOCATE (lri_ppl_ints%lri_ppl(ikind)%v_int(natom, nfa))
         lri_env%stat%ppli_mem = lri_env%stat%ppli_mem + natom*nfa
      END DO

   END SUBROUTINE allocate_lri_ppl_ints

! **************************************************************************************************
!> \brief allocate lri_ints_rho, storing integral for the exact density
!> \param lri_env ...
!> \param lri_ints_rho structure storing the integrals (aa,bb)
!> \param nkind number of atom kinds
! **************************************************************************************************
   SUBROUTINE allocate_lri_ints_rho(lri_env, lri_ints_rho, nkind)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env
      TYPE(lri_list_type), INTENT(OUT)                   :: lri_ints_rho
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: i, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, nba, nbb, nlist, &
                                                            nneighbor
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      lri_ints_rho%nkind = nkind
      ALLOCATE (lri_ints_rho%lri_atom(nkind*nkind))

      DO i = 1, nkind*nkind
         NULLIFY (lri_ints_rho%lri_atom(i)%lri_node)
         lri_ints_rho%lri_atom(i)%natom = 0
      END DO

      CALL neighbor_list_iterator_create(nl_iterator, lri_env%soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                iatom=iatom, jatom=jatom)

         iac = ikind + nkind*(jkind - 1)

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set

         IF (.NOT. ASSOCIATED(obasa)) CYCLE
         IF (.NOT. ASSOCIATED(obasb)) CYCLE

         IF (.NOT. ASSOCIATED(lri_ints_rho%lri_atom(iac)%lri_node)) THEN
            lri_ints_rho%lri_atom(iac)%natom = nlist
            ALLOCATE (lri_ints_rho%lri_atom(iac)%lri_node(nlist))
            DO i = 1, nlist
               NULLIFY (lri_ints_rho%lri_atom(iac)%lri_node(i)%lri_int_rho)
               lri_ints_rho%lri_atom(iac)%lri_node(i)%nnode = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho)) THEN
            lri_ints_rho%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
            ALLOCATE (lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(nneighbor))
         END IF

         ASSOCIATE (lriir => lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor))

            nba = obasa%nsgf
            nbb = obasb%nsgf

            ALLOCATE (lriir%soaabb(nba, nba, nbb, nbb))
            lriir%soaabb = 0._dp
            lriir%dmax_aabb = 0._dp
         END ASSOCIATE

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE allocate_lri_ints_rho

! **************************************************************************************************
!> \brief creates and initializes lri_rhos
!> \param lri_env ...
!> \param lri_rhos structure storing tvec and avec
!> \param nspin ...
!> \param nkind number of atom kinds
! **************************************************************************************************
   SUBROUTINE allocate_lri_rhos(lri_env, lri_rhos, nspin, nkind)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env
      TYPE(lri_list_p_type), DIMENSION(:), POINTER       :: lri_rhos
      INTEGER, INTENT(IN)                                :: nspin, nkind

      INTEGER                                            :: i, iac, iatom, ikind, ilist, ispin, &
                                                            jatom, jkind, jneighbor, nfa, nfb, &
                                                            nlist, nn, nneighbor
      REAL(KIND=dp)                                      :: dab, rab(3)
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (nl_iterator)

      ALLOCATE (lri_rhos(nspin))

      lri_env%stat%rhos_mem = 0.0_dp

      DO ispin = 1, nspin

         ALLOCATE (lri_rhos(ispin)%lri_list)

         lri_rhos(ispin)%lri_list%nkind = nkind
         ALLOCATE (lri_rhos(ispin)%lri_list%lri_atom(nkind*nkind))

         DO i = 1, nkind*nkind
            NULLIFY (lri_rhos(ispin)%lri_list%lri_atom(i)%lri_node)
            lri_rhos(ispin)%lri_list%lri_atom(i)%natom = 0
         END DO

         ASSOCIATE (lri_rho => lri_rhos(ispin)%lri_list)

            CALL neighbor_list_iterator_create(nl_iterator, lri_env%soo_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, nlist=nlist, ilist=ilist, &
                                      nnode=nneighbor, inode=jneighbor, r=rab)

               iac = ikind + nkind*(jkind - 1)
               dab = SQRT(SUM(rab*rab))

               IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

               IF (.NOT. ASSOCIATED(lri_rho%lri_atom(iac)%lri_node)) THEN
                  lri_rho%lri_atom(iac)%natom = nlist
                  ALLOCATE (lri_rho%lri_atom(iac)%lri_node(nlist))
                  DO i = 1, nlist
                     NULLIFY (lri_rho%lri_atom(iac)%lri_node(i)%lri_rhoab)
                     lri_rho%lri_atom(iac)%lri_node(i)%nnode = 0
                  END DO
               END IF
               IF (.NOT. ASSOCIATED(lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab)) THEN
                  lri_rho%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
                  ALLOCATE (lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(nneighbor))
               END IF

               ASSOCIATE (lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor), &
                          lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor))

                  lrho%nba = lrii%nba
                  lrho%nbb = lrii%nbb
                  lrho%nfa = lrii%nfa
                  lrho%nfb = lrii%nfb

                  nfa = lrho%nfa
                  nfb = lrho%nfb
                  nn = nfa + nfb

                  NULLIFY (lrho%avec, lrho%tvec)
                  IF (lrii%lrisr) THEN
                     IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                        IF (.NOT. lri_env%exact_1c_terms) THEN
                           ALLOCATE (lrho%avec(nfa))
                           ALLOCATE (lrho%tvec(nfa))
                           lri_env%stat%rhos_mem = lri_env%stat%rhos_mem + 2*nfa
                        END IF
                     ELSE
                        ALLOCATE (lrho%avec(nn))
                        ALLOCATE (lrho%tvec(nn))
                        lri_env%stat%rhos_mem = lri_env%stat%rhos_mem + 2*nn
                     END IF
                  END IF
                  NULLIFY (lrho%aveca, lrho%tveca)
                  NULLIFY (lrho%avecb, lrho%tvecb)
                  IF (lrii%lriff) THEN
                     ALLOCATE (lrho%aveca(nfa))
                     ALLOCATE (lrho%avecb(nfb))
                     ALLOCATE (lrho%tveca(nfa))
                     ALLOCATE (lrho%tvecb(nfb))
                     lri_env%stat%rhos_mem = lri_env%stat%rhos_mem + 2*(nfa + nfb)
                  END IF

               END ASSOCIATE

            END DO

         END ASSOCIATE

         CALL neighbor_list_iterator_release(nl_iterator)

      END DO

   END SUBROUTINE allocate_lri_rhos

! **************************************************************************************************
!> \brief creates and initializes lri_coefs
!> \param lri_env ...
!> \param lri_density ...
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_coefs(lri_env, lri_density, atomic_kind_set)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env
      TYPE(lri_density_type), INTENT(INOUT)              :: lri_density
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ikind, ispin, natom, nkind, nsgf, nspin
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gto_basis_set_type), POINTER                  :: fbas
      TYPE(lri_spin_type), DIMENSION(:), POINTER         :: lri_coefs

      NULLIFY (atomic_kind, fbas, lri_coefs)

      nkind = SIZE(atomic_kind_set)
      nspin = lri_density%nspin

      lri_env%stat%coef_mem = 0.0_dp

      ALLOCATE (lri_density%lri_coefs(nspin))
      lri_coefs => lri_density%lri_coefs

      DO ispin = 1, nspin
         ALLOCATE (lri_coefs(ispin)%lri_kinds(nkind))
         DO ikind = 1, nkind
            NULLIFY (lri_coefs(ispin)%lri_kinds(ikind)%acoef)
            NULLIFY (lri_coefs(ispin)%lri_kinds(ikind)%v_int)
            NULLIFY (lri_coefs(ispin)%lri_kinds(ikind)%v_dadr)
            NULLIFY (lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
            fbas => lri_env%ri_basis(ikind)%gto_basis_set
            nsgf = fbas%nsgf
            ALLOCATE (lri_coefs(ispin)%lri_kinds(ikind)%acoef(natom, nsgf))
            lri_coefs(ispin)%lri_kinds(ikind)%acoef = 0._dp
            ALLOCATE (lri_coefs(ispin)%lri_kinds(ikind)%v_int(natom, nsgf))
            lri_coefs(ispin)%lri_kinds(ikind)%v_int = 0._dp
            ALLOCATE (lri_coefs(ispin)%lri_kinds(ikind)%v_dadr(natom, 3))
            lri_coefs(ispin)%lri_kinds(ikind)%v_dadr = 0._dp
            ALLOCATE (lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr(natom, 3))
            lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr = 0._dp
            !
            lri_env%stat%coef_mem = lri_env%stat%coef_mem + 2._dp*natom*(nsgf + 3)
         END DO
      END DO

   END SUBROUTINE allocate_lri_coefs

! **************************************************************************************************
!> \brief creates and initializes lri_force
!> \param lri_force ...
!> \param nfa and nfb number of fit functions on a/b
!> \param nfb ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_force_components(lri_force, nfa, nfb)

      TYPE(lri_force_type), INTENT(OUT)                  :: lri_force
      INTEGER, INTENT(IN)                                :: nfa, nfb

      INTEGER                                            :: nn

      nn = nfa + nfb

      ALLOCATE (lri_force%st(nn))
      lri_force%st = 0._dp
      ALLOCATE (lri_force%dsst(nn, 3))
      lri_force%dsst = 0._dp
      ALLOCATE (lri_force%dssn(nn, 3))
      lri_force%dssn = 0._dp
      ALLOCATE (lri_force%dtvec(nn, 3))
      lri_force%dtvec = 0._dp

   END SUBROUTINE allocate_lri_force_components

! **************************************************************************************************
!> \brief deallocates one-center overlap integrals, integral of ri basis
!>        and scon matrices
!> \param lri_env ...
! **************************************************************************************************
   SUBROUTINE deallocate_bas_properties(lri_env)

      TYPE(lri_environment_type), INTENT(INOUT)          :: lri_env

      INTEGER                                            :: i

      IF (ASSOCIATED(lri_env%bas_prop)) THEN
         DO i = 1, SIZE(lri_env%bas_prop)
            IF (ASSOCIATED(lri_env%bas_prop(i)%int_fbas)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%int_fbas)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%ri_ovlp)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%ri_ovlp)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%ri_ovlp_inv)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%ri_ovlp_inv)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%orb_ovlp)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%orb_ovlp)
            END IF
            IF (ALLOCATED(lri_env%bas_prop(i)%ovlp3)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%ovlp3)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%scon_ri)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%scon_ri)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%scon_orb)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%scon_orb)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%orb_index)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%orb_index)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%ri_index)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%ri_index)
            END IF
            IF (ASSOCIATED(lri_env%bas_prop(i)%scon_mix)) THEN
               DEALLOCATE (lri_env%bas_prop(i)%scon_mix)
            END IF
         END DO
         DEALLOCATE (lri_env%bas_prop)
      END IF

   END SUBROUTINE deallocate_bas_properties

! **************************************************************************************************
!> \brief deallocates the given lri_ints
!> \param lri_ints ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_ints(lri_ints)

      TYPE(lri_list_type), INTENT(INOUT)                 :: lri_ints

      INTEGER                                            :: i, iatom, ijkind, inode, natom, nkind, &
                                                            nnode

      nkind = lri_ints%nkind

      IF (nkind > 0) THEN
         DO ijkind = 1, SIZE(lri_ints%lri_atom)
            natom = lri_ints%lri_atom(ijkind)%natom
            IF (natom > 0) THEN
               DO iatom = 1, natom
                  nnode = lri_ints%lri_atom(ijkind)%lri_node(iatom)%nnode
                  IF (nnode > 0) THEN
                     IF (ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int)) THEN
                        DO inode = 1, nnode
                           ASSOCIATE (lri_int => lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode))
                           IF (ALLOCATED(lri_int%sab)) THEN
                              DEALLOCATE (lri_int%sab)
                           END IF
                           IF (ALLOCATED(lri_int%soo)) THEN
                              DEALLOCATE (lri_int%soo)
                           END IF
                           IF (ALLOCATED(lri_int%abaint)) THEN
                              DEALLOCATE (lri_int%abaint)
                           END IF
                           IF (ALLOCATED(lri_int%abascr)) THEN
                              DEALLOCATE (lri_int%abascr)
                           END IF
                           IF (ALLOCATED(lri_int%abbint)) THEN
                              DEALLOCATE (lri_int%abbint)
                           END IF
                           IF (ALLOCATED(lri_int%abbscr)) THEN
                              DEALLOCATE (lri_int%abbscr)
                           END IF
                           IF (ALLOCATED(lri_int%dsab)) THEN
                              DEALLOCATE (lri_int%dsab)
                           END IF
                           IF (ALLOCATED(lri_int%dsoo)) THEN
                              DEALLOCATE (lri_int%dsoo)
                           END IF
                           IF (ALLOCATED(lri_int%dabbint)) THEN
                              DEALLOCATE (lri_int%dabbint)
                           END IF
                           IF (ALLOCATED(lri_int%dabdaint)) THEN
                              DEALLOCATE (lri_int%dabdaint)
                           END IF
                           IF (ASSOCIATED(lri_int%sinv)) THEN
                              DEALLOCATE (lri_int%sinv)
                           END IF
                           IF (ASSOCIATED(lri_int%n)) THEN
                              DEALLOCATE (lri_int%n)
                           END IF
                           IF (ASSOCIATED(lri_int%sn)) THEN
                              DEALLOCATE (lri_int%sn)
                           END IF
                           IF (ASSOCIATED(lri_int%na)) THEN
                              DEALLOCATE (lri_int%na)
                           END IF
                           IF (ASSOCIATED(lri_int%nb)) THEN
                              DEALLOCATE (lri_int%nb)
                           END IF
                           IF (ASSOCIATED(lri_int%sna)) THEN
                              DEALLOCATE (lri_int%sna)
                           END IF
                           IF (ASSOCIATED(lri_int%snb)) THEN
                              DEALLOCATE (lri_int%snb)
                           END IF
                           ! integral container
                           IF (ASSOCIATED(lri_int%cabai%ca)) THEN
                              DO i = 1, SIZE(lri_int%cabai%ca)
                                 IF (ASSOCIATED(lri_int%cabai%ca(i)%cdp)) DEALLOCATE (lri_int%cabai%ca(i)%cdp)
                                 IF (ASSOCIATED(lri_int%cabai%ca(i)%csp)) DEALLOCATE (lri_int%cabai%ca(i)%csp)
                                 IF (ASSOCIATED(lri_int%cabai%ca(i)%cip)) DEALLOCATE (lri_int%cabai%ca(i)%cip)
                              END DO
                              DEALLOCATE (lri_int%cabai%ca)
                           END IF
                           IF (ASSOCIATED(lri_int%cabbi%ca)) THEN
                              DO i = 1, SIZE(lri_int%cabbi%ca)
                                 IF (ASSOCIATED(lri_int%cabbi%ca(i)%cdp)) DEALLOCATE (lri_int%cabbi%ca(i)%cdp)
                                 IF (ASSOCIATED(lri_int%cabbi%ca(i)%csp)) DEALLOCATE (lri_int%cabbi%ca(i)%csp)
                                 IF (ASSOCIATED(lri_int%cabbi%ca(i)%cip)) DEALLOCATE (lri_int%cabbi%ca(i)%cip)
                              END DO
                              DEALLOCATE (lri_int%cabbi%ca)
                           END IF
                           END ASSOCIATE
                        END DO
                        DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int)
                     END IF
                  END IF
               END DO
               DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node)
            END IF
         END DO
         DEALLOCATE (lri_ints%lri_atom)
      END IF

   END SUBROUTINE deallocate_lri_ints

! **************************************************************************************************
!> \brief deallocates the given lri_ppl_ints
!> \param lri_ppl_ints ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_ppl_ints(lri_ppl_ints)

      TYPE(lri_ppl_int_type), INTENT(INOUT)              :: lri_ppl_ints

      INTEGER                                            :: ikind, nkind

      nkind = SIZE(lri_ppl_ints%lri_ppl)
      DO ikind = 1, nkind
         IF (ASSOCIATED(lri_ppl_ints%lri_ppl(ikind)%v_int)) THEN
            DEALLOCATE (lri_ppl_ints%lri_ppl(ikind)%v_int)
         END IF
      END DO
      DEALLOCATE (lri_ppl_ints%lri_ppl)

   END SUBROUTINE deallocate_lri_ppl_ints

! **************************************************************************************************
!> \brief deallocates the given lri_ints_rho
!> \param lri_ints_rho ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_ints_rho(lri_ints_rho)

      TYPE(lri_list_type), INTENT(INOUT)                 :: lri_ints_rho

      INTEGER                                            :: iatom, ijkind, inode, natom, nkind, nnode

      nkind = lri_ints_rho%nkind

      IF (nkind > 0) THEN
         DO ijkind = 1, SIZE(lri_ints_rho%lri_atom)
            natom = lri_ints_rho%lri_atom(ijkind)%natom
            IF (natom > 0) THEN
               DO iatom = 1, natom
                  nnode = lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%nnode
                  IF (nnode > 0) THEN
                     IF (ASSOCIATED(lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho)) THEN
                        DO inode = 1, nnode
                           IF (ASSOCIATED(lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho(inode)%soaabb)) THEN
                              DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho(inode)%soaabb)
                           END IF
                        END DO
                        DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho)
                     END IF
                  END IF
               END DO
               DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node)
            END IF
         END DO
         DEALLOCATE (lri_ints_rho%lri_atom)
      END IF

   END SUBROUTINE deallocate_lri_ints_rho

! **************************************************************************************************
!> \brief deallocates the given lri_rhos
!> \param lri_rhos ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_rhos(lri_rhos)

      TYPE(lri_list_p_type), DIMENSION(:), POINTER       :: lri_rhos

      INTEGER                                            :: i, iatom, ijkind, inode, natom, nkind, &
                                                            nnode

      IF (ASSOCIATED(lri_rhos)) THEN
         DO i = 1, SIZE(lri_rhos)
            ASSOCIATE (lri_rho => lri_rhos(i)%lri_list)
               nkind = lri_rho%nkind
               IF (nkind > 0) THEN
                  CPASSERT(ASSOCIATED(lri_rho%lri_atom))
                  DO ijkind = 1, SIZE(lri_rho%lri_atom)
                     natom = lri_rho%lri_atom(ijkind)%natom
                     IF (natom > 0) THEN
                        CPASSERT(ASSOCIATED(lri_rho%lri_atom(ijkind)%lri_node))
                        DO iatom = 1, natom
                           nnode = lri_rho%lri_atom(ijkind)%lri_node(iatom)%nnode
                           IF (nnode > 0) THEN
                              IF (ASSOCIATED(lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab)) THEN
                                 DO inode = 1, nnode
                                    ASSOCIATE (lri_rhoab => lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode))
                                       IF (ASSOCIATED(lri_rhoab%avec)) DEALLOCATE (lri_rhoab%avec)
                                       IF (ASSOCIATED(lri_rhoab%tvec)) DEALLOCATE (lri_rhoab%tvec)
                                       IF (ASSOCIATED(lri_rhoab%aveca)) DEALLOCATE (lri_rhoab%aveca)
                                       IF (ASSOCIATED(lri_rhoab%tveca)) DEALLOCATE (lri_rhoab%tveca)
                                       IF (ASSOCIATED(lri_rhoab%avecb)) DEALLOCATE (lri_rhoab%avecb)
                                       IF (ASSOCIATED(lri_rhoab%tvecb)) DEALLOCATE (lri_rhoab%tvecb)
                                    END ASSOCIATE
                                 END DO
                                 DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab)
                              END IF
                           END IF
                        END DO
                        DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node)
                     END IF
                  END DO
                  DEALLOCATE (lri_rho%lri_atom)
               END IF
            END ASSOCIATE
            DEALLOCATE (lri_rhos(i)%lri_list)
         END DO

         DEALLOCATE (lri_rhos)

      END IF
      NULLIFY (lri_rhos)

   END SUBROUTINE deallocate_lri_rhos

! **************************************************************************************************
!> \brief releases the given lri_coefs
!> \param lri_coefs the integral storage environment that is released
! **************************************************************************************************
   SUBROUTINE deallocate_lri_coefs(lri_coefs)
      TYPE(lri_spin_type), DIMENSION(:), POINTER         :: lri_coefs

      INTEGER                                            :: i, j

      IF (ASSOCIATED(lri_coefs)) THEN
         DO i = 1, SIZE(lri_coefs)
            DO j = 1, SIZE(lri_coefs(i)%lri_kinds)
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%acoef)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%acoef)
               END IF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_int)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_int)
               END IF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dadr)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_dadr)
               END IF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dfdr)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_dfdr)
               END IF
            END DO
            DEALLOCATE (lri_coefs(i)%lri_kinds)
         END DO
         DEALLOCATE (lri_coefs)
      END IF

      NULLIFY (lri_coefs)

   END SUBROUTINE deallocate_lri_coefs

! **************************************************************************************************
!> \brief releases the given lri_force_type
!> \param lri_force the integral storage environment that is released
! **************************************************************************************************
   SUBROUTINE deallocate_lri_force_components(lri_force)

      TYPE(lri_force_type), INTENT(INOUT)                :: lri_force

      IF (ASSOCIATED(lri_force%st)) THEN
         DEALLOCATE (lri_force%st)
      END IF
      IF (ASSOCIATED(lri_force%dssn)) THEN
         DEALLOCATE (lri_force%dssn)
      END IF
      IF (ASSOCIATED(lri_force%dsst)) THEN
         DEALLOCATE (lri_force%dsst)
      END IF
      IF (ASSOCIATED(lri_force%dtvec)) THEN
         DEALLOCATE (lri_force%dtvec)
      END IF

   END SUBROUTINE deallocate_lri_force_components

END MODULE lri_environment_types
