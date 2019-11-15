!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kinds,                           ONLY: INT_8,&
                                              dp,&
                                              sp
   USE mathlib,                         ONLY: pswitch
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_o3c_types,                    ONLY: o3c_container_type,&
                                              release_o3c_container
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
! Integral container
   TYPE carray
      INTEGER                                                 :: compression
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: cdp
      REAL(KIND=sp), DIMENSION(:), POINTER                    :: csp
      INTEGER(INT_8), DIMENSION(:), POINTER                   :: cip
   END TYPE carray
   TYPE int_container
      INTEGER                                                 :: na, nb, nc
      TYPE(carray), DIMENSION(:), POINTER                     :: ca
   END TYPE int_container
! **************************************************************************************************
   TYPE lri_rhoab_type
      ! number of spherical basis functions (a)
      INTEGER                                                 :: nba
      ! number of spherical basis functions (b)
      INTEGER                                                 :: nbb
      ! number of spherical fit basis functions (ai)
      INTEGER                                                 :: nfa
      ! number of spherical fit basis functions (bi)
      INTEGER                                                 :: nfb
      ! expansion coeffs for RI density
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: avec
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: aveca
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: avecb
      ! projection coeffs for RI density: SUM_ab (ab,i)*Pab
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvec
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tveca
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvecb
      ! integral (ai) * sinv * tvec
      REAL(KIND=dp)                                           :: nst
      REAL(KIND=dp)                                           :: nsta
      REAL(KIND=dp)                                           :: nstb
      ! Lagrange parameter
      REAL(KIND=dp)                                           :: lambda
      REAL(KIND=dp)                                           :: lambdaa
      REAL(KIND=dp)                                           :: lambdab
      ! Charge of pair density
      REAL(KIND=dp)                                           :: charge
      REAL(KIND=dp)                                           :: chargea
      REAL(KIND=dp)                                           :: chargeb
   END TYPE lri_rhoab_type

! **************************************************************************************************

   TYPE lri_int_type
      ! whether to calculate force for pair
      LOGICAL                                                 :: calc_force_pair
      ! number of spherical basis functions (a)
      INTEGER                                                 :: nba
      ! number of spherical basis functions (b)
      INTEGER                                                 :: nbb
      ! number of spherical fit basis functions (ai)
      INTEGER                                                 :: nfa
      ! number of spherical fit basis functions (bi)
      INTEGER                                                 :: nfb
      ! condition number of overlap matrix
      REAL(KIND=dp)                                           :: cond_num
      ! integrals (a,b,ai)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abaint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                :: abascr
      ! integrals (a,b,b)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: abbint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                :: abbscr
      ! compressed aba integrals
      TYPE(int_container)                                     :: cabai
      ! compressed abb integrals
      TYPE(int_container)                                     :: cabbi
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
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: sinv
      ! integral (ai) / (bi), dim(1..nfa,nfa+1..nfa+nfb)
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: n
      ! sinv * (ai)
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: sn
      ! (ai) * sinv * (ai)
      REAL(KIND=dp)                                           :: nsn
      ! distant pair approximation
      LOGICAL                                                 :: lrisr
      LOGICAL                                                 :: lriff
      REAL(KIND=dp)                                           :: wsr, wff, dwsr, dwff
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: asinv, bsinv
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: na, nb
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: sna, snb
      REAL(KIND=dp)                                           :: nsna, nsnb
      !
      ! dmax: max deviation for integrals of primitive gtos; for debugging
      ! dmax for overlap integrals (ai,bi); fit bas
      REAL(KIND=dp)                                           :: dmax_ab
      ! dmax for overlap integrals (a,b); orb bas
      REAL(KIND=dp)                                           :: dmax_oo
      ! dmax for integrals (a,b,ai)
      REAL(KIND=dp)                                           :: dmax_aba
      ! dmax for integrals (a,b,bi)
      REAL(KIND=dp)                                           :: dmax_abb
   END TYPE lri_int_type

   TYPE lri_int_rho_type
      ! integrals (aa,bb), orb basis
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER              :: soaabb
      ! dmax for (aa,bb) integrals; for debugging
      REAL(KIND=dp)                                           :: dmax_aabb
   END TYPE lri_int_rho_type

   TYPE lri_node_type
      INTEGER                                                 :: nnode
      TYPE(lri_int_type), DIMENSION(:), POINTER               :: lri_int
      TYPE(lri_int_rho_type), DIMENSION(:), POINTER           :: lri_int_rho
      TYPE(lri_rhoab_type), DIMENSION(:), POINTER             :: lri_rhoab
   END TYPE lri_node_type

   TYPE lri_atom_type
      INTEGER                                                 :: natom
      TYPE(lri_node_type), DIMENSION(:), POINTER              :: lri_node
   END TYPE lri_atom_type

   TYPE lri_list_type
      INTEGER                                                 :: nkind
      TYPE(lri_atom_type), DIMENSION(:), POINTER              :: lri_atom
   END TYPE lri_list_type

   TYPE lri_list_p_type
      TYPE(lri_list_type), POINTER                             :: lri_list
   END TYPE lri_list_p_type

! **************************************************************************************************

   TYPE lri_bas_type
      INTEGER, DIMENSION(:, :, :), POINTER                    :: orb_index
      INTEGER, DIMENSION(:, :, :), POINTER                    :: ri_index
      ! integral of ri basis fbas
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: int_fbas
      ! self overlap ri basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: ri_ovlp
      ! inverse of self overlap ri basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: ri_ovlp_inv
      ! self overlap orb basis
      REAL(KIND=dp), DIMENSION(:, :), POINTER                 :: orb_ovlp
      ! self overlap (a,a,fa)
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE          :: ovlp3
      ! contraction matrix for SHG integrals ri basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER              :: scon_ri
      ! contraction matrix for SHG integrals orb basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER              :: scon_orb
      ! contraction matrix for SHG integrals aba/abb
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER           :: scon_mix
   END TYPE lri_bas_type

! **************************************************************************************************

   TYPE lri_clebsch_gordon_type
      ! Clebsch-Gordon (CG) coefficients
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER                :: cg_coeff
      ! list of non-zero CG coefficients
      INTEGER, DIMENSION(:, :, :), POINTER                      :: cg_none0_list
      ! number of non-zero CG coefficients
      INTEGER, DIMENSION(:, :), POINTER                         :: ncg_none0
   END TYPE lri_clebsch_gordon_type

! **************************************************************************************************

   TYPE lri_ppl_type
      ! integrals Vppl*fbas (potential*fit basis) dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_int
   END TYPE lri_ppl_type

   TYPE lri_ppl_int_type
      TYPE(lri_ppl_type), DIMENSION(:), POINTER                 :: lri_ppl
      REAL(KIND=dp)                                             :: ecore_pp_ri
   END TYPE lri_ppl_int_type

! **************************************************************************************************

   TYPE ri_fit_type
      INTEGER, DIMENSION(:, :), POINTER                         :: bas_ptr
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: nvec
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: rm1n
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: tvec
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: rm1t
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: avec
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: fout
      REAL(KIND=dp)                                             :: ntrm1n
      REAL(KIND=dp), DIMENSION(2)                               :: ftrm1n
      REAL(KIND=dp), DIMENSION(2)                               :: echarge
      REAL(KIND=dp), DIMENSION(2)                               :: lambda
   END TYPE ri_fit_type

! **************************************************************************************************
   TYPE wmat_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: mat
   END TYPE wmat_type
   TYPE wbas_type
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: vec
   END TYPE wbas_type
! **************************************************************************************************
   TYPE stat_type
      REAL(KIND=dp)                                             :: pairs_tt
      REAL(KIND=dp)                                             :: pairs_sr
      REAL(KIND=dp)                                             :: pairs_ff
      REAL(KIND=dp)                                             :: overlap_error
      REAL(KIND=dp)                                             :: rho_tt
      REAL(KIND=dp)                                             :: rho_sr
      REAL(KIND=dp)                                             :: rho_ff
      REAL(KIND=dp)                                             :: rho_1c
      REAL(KIND=dp)                                             :: coef_mem
      REAL(KIND=dp)                                             :: oint_mem
      REAL(KIND=dp)                                             :: rhos_mem
      REAL(KIND=dp)                                             :: abai_mem
      REAL(KIND=dp)                                             :: ppli_mem
   END TYPE stat_type
! **************************************************************************************************

   TYPE lri_environment_type
      INTEGER                                                 :: id_nr, ref_count, in_use
      ! parameter for (pseudo)inverse of overlap
      INTEGER                                                 :: lri_overlap_inv
      ! flag for debugging lri integrals
      LOGICAL                                                 :: debug
      ! flag for shg (solid haromonic Gaussian) integrals
      LOGICAL                                                 :: use_shg_integrals
      ! parameter for inversion (autoselect); maximal condition
      ! number up to where inversion is legal
      REAL(KIND=dp)                                           :: cond_max
      ! parameter for checking distance between atom pairs
      REAL(KIND=dp)                                           :: delta
      ! threshold for aba and abb integrals
      REAL(KIND=dp)                                           :: eps_o3_int
      ! orbital basis set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: orb_basis
      ! lri (fit) basis set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis
      ! orb_basis neighborlist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soo_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: saa_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soa_list
      ! local RI integrals
      TYPE(lri_list_type), POINTER                            :: lri_ints
      ! local Vppl integrals
      TYPE(lri_ppl_int_type), POINTER                         :: lri_ppl_ints
      ! local integral of rho**2; for optimization
      TYPE(lri_list_type), POINTER                            :: lri_ints_rho
      ! properties of orb and aux basis
      TYPE(lri_bas_type), DIMENSION(:), POINTER               :: bas_prop
      ! Clebsch-Gordon for solid harmonics
      TYPE(lri_clebsch_gordon_type), POINTER                  :: cg_shg
      ! orbital basis overlap
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ob_smat
      ! statistics
      LOGICAL                                                 :: statistics
      TYPE(stat_type)                                         :: stat
      ! exact one-center terms
      LOGICAL                                                 :: exact_1c_terms
      ! use RI for local pp
      LOGICAL                                                 :: ppl_ri
      ! store integrals (needed for basis optimization)
      LOGICAL                                                 :: store_integrals
      ! distant pair approximation
      LOGICAL                                                 :: distant_pair_approximation
      CHARACTER(len=10)                                       :: distant_pair_method
      REAL(KIND=dp)                                           :: r_in
      REAL(KIND=dp)                                           :: r_out
      REAL(KIND=dp), DIMENSION(:), POINTER                    :: aradius
      TYPE(wbas_type), DIMENSION(:), POINTER                  :: wbas
      TYPE(wmat_type), DIMENSION(:, :), POINTER               :: wmat
      ! RI overlap and inverse
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER               :: ri_smat, &
                                                                 ri_sinv
      TYPE(ri_fit_type), POINTER                              :: ri_fit
      CHARACTER(len=10)                                       :: ri_sinv_app
      TYPE(o3c_container_type), POINTER                       :: o3c
   END TYPE lri_environment_type

! **************************************************************************************************

   TYPE lri_kind_type
      ! expansion coeff for lri density dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: acoef
      ! integrals V*fbas (potential*fit basis) dim(natom,nsgf)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_int
      ! SUM_i integral(V*fbas_i)*davec/dR dim(natom,3)
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_dadr
      ! integrals V*dfbas/dR
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: v_dfdr
   END TYPE lri_kind_type

   TYPE lri_spin_type
      TYPE(lri_kind_type), DIMENSION(:), POINTER                :: lri_kinds
   END TYPE lri_spin_type

! **************************************************************************************************

   TYPE lri_force_type
      REAL(KIND=dp), DIMENSION(:), POINTER                      :: st
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: dssn, &
                                                                   dsst
      ! derivative dtvec/dR
      REAL(KIND=dp), DIMENSION(:, :), POINTER                   :: dtvec
   END TYPE lri_force_type

! **************************************************************************************************

   TYPE lri_density_type
      INTEGER                                                 :: id_nr, ref_count, in_use
      INTEGER                                                 :: nspin
      ! pair density expansion (nspin)
      TYPE(lri_list_p_type), DIMENSION(:), POINTER            :: lri_rhos
      ! coefficients of RI expansion and gradients (nspin)
      TYPE(lri_spin_type), DIMENSION(:), POINTER              :: lri_coefs
   END TYPE lri_density_type

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_types'
   INTEGER, SAVE, PRIVATE :: last_lri_env_id = 0
   INTEGER, SAVE, PRIVATE :: last_lri_density_id = 0

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

      TYPE(lri_environment_type), POINTER                :: lri_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_env_create', routineP = moduleN//':'//routineN

      ALLOCATE (lri_env)

      last_lri_env_id = last_lri_env_id + 1
      lri_env%id_nr = last_lri_env_id
      lri_env%ref_count = 1
      lri_env%in_use = 0

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

      TYPE(lri_environment_type), POINTER                :: lri_env

      CHARACTER(len=*), PARAMETER :: routineN = 'lri_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ikind, j, nkind

      IF (ASSOCIATED(lri_env)) THEN
         lri_env%ref_count = 0

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
         IF (ASSOCIATED(lri_env%soo_list)) THEN
            DO i = 1, SIZE(lri_env%soo_list)
               CALL deallocate_neighbor_list_set(lri_env%soo_list(i)%neighbor_list_set)
            END DO
            DEALLOCATE (lri_env%soo_list)
         END IF
         IF (ASSOCIATED(lri_env%saa_list)) THEN
            DO i = 1, SIZE(lri_env%saa_list)
               CALL deallocate_neighbor_list_set(lri_env%saa_list(i)%neighbor_list_set)
            END DO
            DEALLOCATE (lri_env%saa_list)
         END IF
         IF (ASSOCIATED(lri_env%soa_list)) THEN
            DO i = 1, SIZE(lri_env%soa_list)
               CALL deallocate_neighbor_list_set(lri_env%soa_list(i)%neighbor_list_set)
            END DO
            DEALLOCATE (lri_env%soa_list)
         END IF
         IF (ASSOCIATED(lri_env%lri_ints)) THEN
            CALL deallocate_lri_ints(lri_env%lri_ints)
         END IF
         IF (ASSOCIATED(lri_env%lri_ppl_ints)) THEN
            CALL deallocate_lri_ppl_ints(lri_env%lri_ppl_ints)
         END IF
         IF (ASSOCIATED(lri_env%lri_ints_rho)) THEN
            CALL deallocate_lri_ints_rho(lri_env%lri_ints_rho)
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
         ENDIF
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
         !
         DEALLOCATE (lri_env)
      END IF
      NULLIFY (lri_env)

   END SUBROUTINE lri_env_release

! **************************************************************************************************
!> \brief creates and initializes an lri_density environment
!> \param lri_density the lri_density environment you want to create
! **************************************************************************************************
   SUBROUTINE lri_density_create(lri_density)

      TYPE(lri_density_type), POINTER                    :: lri_density

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_density_create', &
         routineP = moduleN//':'//routineN

      ALLOCATE (lri_density)

      last_lri_density_id = last_lri_density_id + 1
      lri_density%id_nr = last_lri_density_id
      lri_density%ref_count = 1
      lri_density%in_use = 0

      lri_density%nspin = 0

      NULLIFY (lri_density%lri_rhos)
      NULLIFY (lri_density%lri_coefs)

   END SUBROUTINE lri_density_create

! **************************************************************************************************
!> \brief releases the given lri_density
!> \param lri_density the lri_density to release
! **************************************************************************************************
   SUBROUTINE lri_density_release(lri_density)
      TYPE(lri_density_type), POINTER                    :: lri_density

      CHARACTER(len=*), PARAMETER :: routineN = 'lri_density_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(lri_density)) THEN
         lri_density%ref_count = 0

         CALL deallocate_lri_rhos(lri_density%lri_rhos)
         CALL deallocate_lri_coefs(lri_density%lri_coefs)

         DEALLOCATE (lri_density)
      END IF

      NULLIFY (lri_density)

   END SUBROUTINE lri_density_release

! **************************************************************************************************
!> \brief allocate lri_ints, matrices that store LRI integrals
!> \param lri_env ...
!> \param lri_ints structure storing the LRI integrals
!> \param nkind number of atom kinds
! **************************************************************************************************
   SUBROUTINE allocate_lri_ints(lri_env, lri_ints, nkind)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_list_type), POINTER                       :: lri_ints
      INTEGER, INTENT(IN)                                :: nkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_ints', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, nba, nbb, nfa, nfb, &
                                                            nlist, nn, nneighbor
      LOGICAL                                            :: dpa, e1c
      REAL(KIND=dp)                                      :: dab, ra, rab(3), rb
      TYPE(gto_basis_set_type), POINTER                  :: fbasa, fbasb, obasa, obasb
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (fbasa, fbasb, lrii, nl_iterator, obasa, obasb)

      ALLOCATE (lri_ints)

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

         lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf
         nfa = fbasa%nsgf
         nfb = fbasb%nsgf
         nn = nfa + nfb

         IF (iatom == jatom .AND. dab < lri_env%delta) THEN
            e1c = lri_env%exact_1c_terms
         ELSE
            e1c = .FALSE.
         ENDIF

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
               ENDIF
               lrii%sinv = 0._dp

               IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                  ALLOCATE (lrii%n(nfa), lrii%sn(nfa))
                  lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nfa
               ELSE
                  ALLOCATE (lrii%n(nn), lrii%sn(nn))
                  lri_env%stat%oint_mem = lri_env%stat%oint_mem + 2.*nn
               ENDIF
               lrii%n = 0._dp
               lrii%sn = 0._dp
            END IF
         ELSE
            NULLIFY (lrii%n, lrii%sn)
            NULLIFY (lrii%sinv)
         ENDIF

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
         ENDIF

         lrii%dmax_ab = 0._dp
         lrii%dmax_oo = 0._dp
         lrii%dmax_aba = 0._dp
         lrii%dmax_abb = 0._dp

         lrii%calc_force_pair = .FALSE.

      ENDDO

      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE allocate_lri_ints

! **************************************************************************************************
!> \brief allocate lri_ppl_ints, matrices that store LRI integrals
!> \param lri_env ...
!> \param lri_ppl_ints structure storing the LRI ppl integrals
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_ppl_ints(lri_env, lri_ppl_ints, atomic_kind_set)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_ppl_int_type), POINTER                    :: lri_ppl_ints
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_ppl_ints', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, natom, nfa, nkind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gto_basis_set_type), POINTER                  :: fbasa

      lri_env%stat%ppli_mem = 0.0_dp
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (lri_ppl_ints)
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

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_list_type), POINTER                       :: lri_ints_rho
      INTEGER, INTENT(IN)                                :: nkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_ints_rho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iac, iatom, ikind, ilist, jatom, &
                                                            jkind, jneighbor, nba, nbb, nlist, &
                                                            nneighbor
      TYPE(gto_basis_set_type), POINTER                  :: obasa, obasb
      TYPE(lri_int_rho_type), POINTER                    :: lriir
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      ALLOCATE (lri_ints_rho)

      lri_ints_rho%nkind = nkind
      ALLOCATE (lri_ints_rho%lri_atom(nkind*nkind))

      DO i = 1, nkind*nkind
         NULLIFY (lri_ints_rho%lri_atom(i)%lri_node)
         lri_ints_rho%lri_atom(i)%natom = 0
      ENDDO

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

         lriir => lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf

         ALLOCATE (lriir%soaabb(nba, nba, nbb, nbb))
         lriir%soaabb = 0._dp
         lriir%dmax_aabb = 0._dp

      ENDDO

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

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_list_p_type), DIMENSION(:), POINTER       :: lri_rhos
      INTEGER, INTENT(IN)                                :: nspin, nkind

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_rhos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iac, iatom, ikind, ilist, ispin, &
                                                            jatom, jkind, jneighbor, nfa, nfb, &
                                                            nlist, nn, nneighbor
      REAL(KIND=dp)                                      :: dab, rab(3)
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lrho
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      NULLIFY (lri_rho, lrho, lrii, nl_iterator)

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

         lri_rho => lri_rhos(ispin)%lri_list

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

            lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

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
               ENDIF
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

         ENDDO

         CALL neighbor_list_iterator_release(nl_iterator)

      ENDDO

   END SUBROUTINE allocate_lri_rhos

! **************************************************************************************************
!> \brief creates and initializes lri_coefs
!> \param lri_env ...
!> \param lri_density ...
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_coefs(lri_env, lri_density, atomic_kind_set)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_coefs', &
         routineP = moduleN//':'//routineN

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
      ENDDO

   END SUBROUTINE allocate_lri_coefs

! **************************************************************************************************
!> \brief creates and initializes lri_force
!> \param lri_force ...
!> \param nfa and nfb number of fit functions on a/b
!> \param nfb ...
! **************************************************************************************************
   SUBROUTINE allocate_lri_force_components(lri_force, nfa, nfb)

      TYPE(lri_force_type), POINTER                      :: lri_force
      INTEGER, INTENT(IN)                                :: nfa, nfb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_force_components', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nn

      nn = nfa + nfb

      CPASSERT(.NOT. ASSOCIATED(lri_force))

      ALLOCATE (lri_force)

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

      TYPE(lri_environment_type), POINTER                :: lri_env

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_bas_properties', &
         routineP = moduleN//':'//routineN

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

      TYPE(lri_list_type), POINTER                       :: lri_ints

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, ijkind, inode, natom, nkind, &
                                                            nnode
      TYPE(lri_int_type), POINTER                        :: lri_int

      CPASSERT(ASSOCIATED(lri_ints))
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
                           lri_int => lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)
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
      DEALLOCATE (lri_ints)

   END SUBROUTINE deallocate_lri_ints

! **************************************************************************************************
!> \brief deallocates the given lri_ppl_ints
!> \param lri_ppl_ints ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_ppl_ints(lri_ppl_ints)

      TYPE(lri_ppl_int_type), POINTER                    :: lri_ppl_ints

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ppl_ints', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nkind

      CPASSERT(ASSOCIATED(lri_ppl_ints))
      IF (ASSOCIATED(lri_ppl_ints%lri_ppl)) THEN
         nkind = SIZE(lri_ppl_ints%lri_ppl)
         DO ikind = 1, nkind
            IF (ASSOCIATED(lri_ppl_ints%lri_ppl(ikind)%v_int)) THEN
               DEALLOCATE (lri_ppl_ints%lri_ppl(ikind)%v_int)
            END IF
         END DO
         DEALLOCATE (lri_ppl_ints%lri_ppl)
      END IF
      DEALLOCATE (lri_ppl_ints)

   END SUBROUTINE deallocate_lri_ppl_ints

! **************************************************************************************************
!> \brief deallocates the given lri_ints_rho
!> \param lri_ints_rho ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_ints_rho(lri_ints_rho)

      TYPE(lri_list_type), POINTER                       :: lri_ints_rho

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints_rho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ijkind, inode, natom, nkind, nnode

      CPASSERT(ASSOCIATED(lri_ints_rho))
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
      DEALLOCATE (lri_ints_rho)

   END SUBROUTINE deallocate_lri_ints_rho

! **************************************************************************************************
!> \brief deallocates the given lri_rhos
!> \param lri_rhos ...
! **************************************************************************************************
   SUBROUTINE deallocate_lri_rhos(lri_rhos)

      TYPE(lri_list_p_type), DIMENSION(:), POINTER       :: lri_rhos

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_rhos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, ijkind, inode, natom, nkind, &
                                                            nnode
      TYPE(lri_list_type), POINTER                       :: lri_rho
      TYPE(lri_rhoab_type), POINTER                      :: lri_rhoab

      NULLIFY (lri_rho)

      IF (ASSOCIATED(lri_rhos)) THEN

         DO i = 1, SIZE(lri_rhos)

            lri_rho => lri_rhos(i)%lri_list
            CPASSERT(ASSOCIATED(lri_rho))
            nkind = lri_rho%nkind

            IF (nkind > 0) THEN
               DO ijkind = 1, SIZE(lri_rho%lri_atom)
                  natom = lri_rho%lri_atom(ijkind)%natom
                  IF (natom > 0) THEN
                     DO iatom = 1, natom
                        nnode = lri_rho%lri_atom(ijkind)%lri_node(iatom)%nnode
                        IF (nnode > 0) THEN
                           IF (ASSOCIATED(lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab)) THEN
                              DO inode = 1, nnode
                                 lri_rhoab => lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode)
                                 IF (ASSOCIATED(lri_rhoab%avec)) DEALLOCATE (lri_rhoab%avec)
                                 IF (ASSOCIATED(lri_rhoab%tvec)) DEALLOCATE (lri_rhoab%tvec)
                                 IF (ASSOCIATED(lri_rhoab%aveca)) DEALLOCATE (lri_rhoab%aveca)
                                 IF (ASSOCIATED(lri_rhoab%tveca)) DEALLOCATE (lri_rhoab%tveca)
                                 IF (ASSOCIATED(lri_rhoab%avecb)) DEALLOCATE (lri_rhoab%avecb)
                                 IF (ASSOCIATED(lri_rhoab%tvecb)) DEALLOCATE (lri_rhoab%tvecb)
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
            DEALLOCATE (lri_rho)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_coefs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j

      IF (ASSOCIATED(lri_coefs)) THEN
         DO i = 1, SIZE(lri_coefs)
            DO j = 1, SIZE(lri_coefs(i)%lri_kinds)
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%acoef)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%acoef)
               ENDIF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_int)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_int)
               ENDIF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dadr)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_dadr)
               ENDIF
               IF (ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dfdr)) THEN
                  DEALLOCATE (lri_coefs(i)%lri_kinds(j)%v_dfdr)
               ENDIF
            ENDDO
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

      TYPE(lri_force_type), POINTER                      :: lri_force

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_force_components', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(lri_force)) THEN

         IF (ASSOCIATED(lri_force%st)) THEN
            DEALLOCATE (lri_force%st)
         ENDIF
         IF (ASSOCIATED(lri_force%dssn)) THEN
            DEALLOCATE (lri_force%dssn)
         ENDIF
         IF (ASSOCIATED(lri_force%dsst)) THEN
            DEALLOCATE (lri_force%dsst)
         ENDIF
         IF (ASSOCIATED(lri_force%dtvec)) THEN
            DEALLOCATE (lri_force%dtvec)
         ENDIF

         DEALLOCATE (lri_force)

         NULLIFY (lri_force)
      ENDIF

   END SUBROUTINE deallocate_lri_force_components

END MODULE lri_environment_types

