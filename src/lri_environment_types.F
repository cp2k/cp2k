!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the types and subroutines for dealing with the lri_env
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) debugged
!> \authors JGH 
!>          Dorothea Golze  
! *****************************************************************************
MODULE lri_environment_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac,&
                                             pi
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
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
     ! projection coeffs for RI density: SUM_ab (ab,i)*Pab
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: tvec
     ! Lagrange parameter
     REAL(KIND=dp)                                           :: lambda
     ! Charge of pair density
     REAL(KIND=dp)                                           :: charge
  END TYPE lri_rhoab_type

  TYPE lri_int_type
     ! number of spherical basis functions (a)
     INTEGER                                                 :: nba
     ! number of spherical basis functions (b)
     INTEGER                                                 :: nbb
     ! number of spherical fit basis functions (ai)
     INTEGER                                                 :: nfa
     ! number of spherical fit basis functions (bi)
     INTEGER                                                 :: nfb
     ! integrals (a,b,ai)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abaint
     ! integrals (a,b,bi)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abbint
     ! integrals (a,b)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: soo
     ! integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sab
     ! derivative of fit coeff with respect to pmatrix
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: dacoef
     ! inverse of integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sinv
     ! sinv * (ai)
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: sn
     ! (ai) * sinv * (ai)
     REAL(KIND=dp)                                           :: nsn
  END TYPE lri_int_type     

  TYPE lri_node_type
     INTEGER                                                 :: nnode
     TYPE(lri_int_type), DIMENSION(:), POINTER               :: lri_int
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
     TYPE(lri_list_type),POINTER                             :: lri_list
  END TYPE lri_list_p_type
  
! *****************************************************************************

  TYPE lri_bas_norm_type
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: norm_radial
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: norm_s_func
  END TYPE lri_bas_norm_type

  TYPE lri_bas_int_type
     ! integral of ri basis fbas
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: int_fbas
  END TYPE lri_bas_int_type

  TYPE lri_bas_overlap_type
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: ovlp
  END TYPE lri_bas_overlap_type

! *****************************************************************************

  TYPE lri_environment_type
     INTEGER :: id_nr, ref_count, in_use
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: orb_basis
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis
     ! orb_basis neighborlist, LRI integrals
     TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soo_list
     ! local RI integrals
     TYPE(lri_list_type), POINTER                            :: lri_ints
     ! norm of ri basis
     TYPE(lri_bas_norm_type), DIMENSION(:), POINTER          :: bas_norm
     ! integral of a single contracted gaussian
     TYPE(lri_bas_int_type), DIMENSION(:), POINTER           :: bas_int
     ! self overlap of ri basis
     TYPE(lri_bas_overlap_type), DIMENSION(:), POINTER       :: bas_ovlp
  END TYPE lri_environment_type

! *****************************************************************************

  TYPE lri_kind_type
    ! expansion coeff for lri density dim(natom,nsgf)
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: acoef
    ! integrals V*fbas (potentia*fit basis) dim(natom,nsgf)
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: v_int
  END TYPE lri_kind_type

  TYPE lri_spin_type
     TYPE(lri_kind_type), DIMENSION(:), POINTER              :: lri_kinds
  END TYPE lri_spin_type

! *****************************************************************************

  TYPE lri_density_type
     INTEGER :: id_nr, ref_count, in_use
     INTEGER                                                 :: nspin
     ! pair density expansion (nspin)
     TYPE(lri_list_p_type), DIMENSION(:),POINTER             :: lri_rhos
     ! coefficients of RI expansion and gradients (nspin)
     TYPE(lri_spin_type), DIMENSION(:),POINTER               :: lri_coefs
  END TYPE lri_density_type

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_types'
  INTEGER, SAVE, PRIVATE :: last_lri_env_id=0
  INTEGER, SAVE, PRIVATE :: last_lri_density_id=0

  PUBLIC :: lri_environment_type, lri_bas_norm_type, lri_bas_int_type,&
            lri_bas_overlap_type, lri_list_type, lri_atom_type, lri_node_type,&
            lri_int_type, lri_density_type, lri_list_p_type, lri_spin_type,&
            lri_kind_type, lri_rhoab_type 
  PUBLIC :: lri_env_create, lri_env_release, allocate_lri_coefs,&
            lri_density_create, lri_density_release, allocate_lri_list_type,&
            lri_basis_init, deallocate_lri_ints      

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_env_create (lri_env, error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( lri_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_lri_env_id=last_lri_env_id+1
    lri_env%id_nr=last_lri_env_id
    lri_env%ref_count=1
    lri_env%in_use=0
  
    NULLIFY(lri_env%orb_basis)
    NULLIFY(lri_env%ri_basis)

    NULLIFY(lri_env%soo_list)
    NULLIFY(lri_env%lri_ints)
    NULLIFY(lri_env%bas_norm)
    NULLIFY(lri_env%bas_int)
    NULLIFY(lri_env%bas_ovlp)

  END SUBROUTINE lri_env_create

! *****************************************************************************
!> \brief releases the given lri_env
!> \param lri_env the lri environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_env_release(lri_env, error)
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ikind, nkind, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_env)) THEN
       lri_env%ref_count=0

       ! deallocate basis sets
       IF(ASSOCIATED(lri_env%orb_basis)) THEN
          nkind = SIZE(lri_env%orb_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%orb_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%orb_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF(ASSOCIATED(lri_env%ri_basis)) THEN
          nkind = SIZE(lri_env%ri_basis)
          DO ikind=1,nkind
             CALL deallocate_gto_basis_set(lri_env%ri_basis(ikind)%gto_basis_set,error)
          END DO
          DEALLOCATE(lri_env%ri_basis,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%soo_list)) THEN
          DO i=1,SIZE(lri_env%soo_list)
             CALL deallocate_neighbor_list_set(lri_env%soo_list(i)%neighbor_list_set)
          END DO
          DEALLOCATE(lri_env%soo_list,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%lri_ints)) THEN
          CALL deallocate_lri_ints(lri_env%lri_ints,error)
       END IF
       IF (ASSOCIATED(lri_env%bas_norm)) THEN
          DO i=1,SIZE(lri_env%bas_norm)
             DEALLOCATE(lri_env%bas_norm(i)%norm_radial,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             DEALLOCATE(lri_env%bas_norm(i)%norm_s_func,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO
          DEALLOCATE(lri_env%bas_norm,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%bas_int)) THEN
          DO i=1,SIZE(lri_env%bas_int)
             DEALLOCATE(lri_env%bas_int(i)%int_fbas,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO
          DEALLOCATE(lri_env%bas_int,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       IF (ASSOCIATED(lri_env%bas_ovlp)) THEN
          DO i=1,SIZE(lri_env%bas_ovlp)
             DEALLOCATE(lri_env%bas_ovlp(i)%ovlp,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO
          DEALLOCATE(lri_env%bas_ovlp,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       DEALLOCATE(lri_env,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_env)

  END SUBROUTINE lri_env_release

! *****************************************************************************
!> \brief creates and initializes an lri_density environment
!> \param lri_density the lri_density environment you want to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_density_create(lri_density, error)

    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_density_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( lri_density, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    last_lri_density_id=last_lri_density_id+1
    lri_density%id_nr=last_lri_density_id
    lri_density%ref_count=1
    lri_density%in_use=0

    lri_density%nspin=0
  
    NULLIFY(lri_density%lri_rhos)
    NULLIFY(lri_density%lri_coefs)

  END SUBROUTINE lri_density_create

! *****************************************************************************
!> \brief releases the given lri_env
!> \param lri_env the lri environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_density_release(lri_density, error)
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_density_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_density)) THEN
       lri_density%ref_count=0

       IF (ASSOCIATED(lri_density%lri_rhos)) THEN
          DO i=1,SIZE(lri_density%lri_rhos)
             CALL deallocate_lri_rho(lri_density%lri_rhos(i)%lri_list,error)
          END DO
          DEALLOCATE(lri_density%lri_rhos,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF

       CALL deallocate_lri_coefs(lri_density%lri_coefs, error)

       DEALLOCATE(lri_density,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_density)

  END SUBROUTINE lri_density_release

! *****************************************************************************
!> \brief initializes the lri basis: calculates the norm, self-overlap
!>        and integral of the ri basis
! *****************************************************************************
  SUBROUTINE lri_basis_init(lri_env,atomic_kind_set,error)
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_basis_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ri_int_fbas, ri_norm_r, &
                                                ri_norm_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ri_ovlp
    TYPE(gto_basis_set_type), POINTER        :: orb_basis, ri_basis

    failure=.FALSE.

    IF (ASSOCIATED(lri_env)) THEN
       IF (ASSOCIATED(lri_env%orb_basis)) THEN
          CPPreconditionNoFail(ASSOCIATED(lri_env%ri_basis),cp_warning_level,routineP,error)
          nkind = SIZE(lri_env%orb_basis)
          ALLOCATE (lri_env%bas_norm(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%bas_int(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (lri_env%bas_ovlp(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ikind = 1, nkind
             NULLIFY(orb_basis, ri_basis)
             orb_basis => lri_env%orb_basis(ikind)%gto_basis_set
             IF (ASSOCIATED(orb_basis)) THEN
                ri_basis => lri_env%ri_basis(ikind)%gto_basis_set
                CPPreconditionNoFail(ASSOCIATED(ri_basis),cp_warning_level,routineP,error)
                NULLIFY(ri_norm_r)
                CALL basis_norm_radial(ri_basis,ri_norm_r,error)
                lri_env%bas_norm(ikind)%norm_radial => ri_norm_r
                NULLIFY(ri_norm_s)
                CALL basis_norm_s_func(ri_basis,ri_norm_s,error)
                lri_env%bas_norm(ikind)%norm_s_func => ri_norm_s
                NULLIFY(ri_int_fbas)
                CALL basis_int(ri_basis,ri_int_fbas,ri_norm_s,error)
                lri_env%bas_int(ikind)%int_fbas => ri_int_fbas
                NULLIFY(ri_ovlp)
                CALL basis_ovlp(ri_basis,ri_ovlp,ri_norm_r,error)
                lri_env%bas_ovlp(ikind)%ovlp => ri_ovlp
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE lri_basis_init

!*****************************************************************************
!> \brief normalization for a contracted Gaussian s-function,
!>        spherical = cartesian Gaussian for s-functions
!*****************************************************************************
  SUBROUTINE basis_norm_s_func(basis,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:), POINTER          :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_s_func', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, &
                                                jpgf, l, nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aai, aaj, cci, ccj, expa, ppl

    NULLIFY(norm)

    nbas = basis%nsgf
    ALLOCATE (norm(nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    norm = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          l = basis%l(ishell,iset)
          IF (l/=0) CYCLE
          expa = 0.5_dp*REAL(2*l + 3,dp)
          ppl = pi**(3._dp/2._dp)
          DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
             DO ipgf=1,basis%npgf(iset)
                cci = basis%gcc(ipgf,ishell,iset)
                aai = basis%zet(ipgf,iset)
                DO jpgf=1,basis%npgf(iset)
                   ccj = basis%gcc(jpgf,ishell,iset)
                   aaj = basis%zet(jpgf,iset)
                   norm(isgf) =norm(isgf) + cci*ccj*ppl/(aai+aaj)**expa
                END DO
             END DO
             norm(isgf)=1.0_dp/SQRT(norm(isgf)) 
          END DO
       END DO
    END DO

  END SUBROUTINE basis_norm_s_func

!*****************************************************************************
!> \brief normalization for radial part of contracted spherical Gaussian 
!>        functions
!*****************************************************************************
  SUBROUTINE basis_norm_radial(basis,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:), POINTER          :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_radial', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, &
                                                jpgf, l, nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aai, aaj, cci, ccj, expa, ppl

    NULLIFY(norm)

    nbas = basis%nsgf
    ALLOCATE (norm(nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    norm = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          l = basis%l(ishell,iset)
          expa = 0.5_dp*REAL(2*l + 3,dp)
          ppl = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3,dp)/fac(l+1)
          DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
             DO ipgf=1,basis%npgf(iset)
                cci = basis%gcc(ipgf,ishell,iset)
                aai = basis%zet(ipgf,iset)
                DO jpgf=1,basis%npgf(iset)
                   ccj = basis%gcc(jpgf,ishell,iset)
                   aaj = basis%zet(jpgf,iset)
                   norm(isgf) =norm(isgf) + cci*ccj*ppl/(aai+aaj)**expa
                END DO
             END DO
             norm(isgf)=1.0_dp/SQRT(norm(isgf)) 
          END DO
       END DO
    END DO

  END SUBROUTINE basis_norm_radial

!*****************************************************************************
!> \brief integral over a single (contracted) lri auxiliary basis function,
!>        integral is zero for all but s-functions               
!*****************************************************************************
  SUBROUTINE basis_int(basis,int_aux,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:), POINTER          :: int_aux, norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, l, &
                                                nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, cc, pp

    nbas = basis%nsgf
    ALLOCATE (int_aux(nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    int_aux = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          l = basis%l(ishell,iset)
          IF (l/=0) CYCLE
          DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
             DO ipgf=1,basis%npgf(iset)
                cc = basis%gcc(ipgf,ishell,iset)
                aa = basis%zet(ipgf,iset)
                pp = (pi/aa)**(3._dp/2._dp)
                int_aux(isgf) = int_aux(isgf) + norm(isgf)*cc*pp
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE basis_int

!*****************************************************************************
!> \brief self-overlap of lri basis for contracted spherical Gaussians.
!>        Overlap of radial part. Norm contains only normalization of radial
!>        part. Norm and overlap of spherical harmonics not explicitly 
!>        calculated since this cancels for the self-overlap anyway. 
!*****************************************************************************
  SUBROUTINE basis_ovlp(basis,ovlp,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:, :), POINTER       :: ovlp
    REAL(dp), DIMENSION(:), POINTER          :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_ovlp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, isgf, ishell, &
                                                jpgf, jset, jsgf, jshell, l, &
                                                li, lj, nbas, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aai, aaj, cci, ccj, expa, &
                                                m_i, m_j, norm_i, norm_j, oo, &
                                                ppl

    nbas = basis%nsgf
    ALLOCATE (ovlp(nbas,nbas),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ovlp = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          li = basis%l(ishell,iset)
          DO jset=1,basis%nset
             DO jshell=1,basis%nshell(jset)
                lj = basis%l(jshell,jset)
                IF (li == lj) THEN
                   l = li
                   expa = 0.5_dp*REAL(2*l + 3,dp)
                   ppl = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3,dp)/fac(l+1)
                   DO isgf=basis%first_sgf(ishell,iset),basis%last_sgf(ishell,iset)
                      m_i=basis%m(isgf)
                      DO jsgf=basis%first_sgf(jshell,jset),basis%last_sgf(jshell,jset)
                         m_j=basis%m(jsgf)
                         IF(m_i == m_j) THEN
                           DO ipgf=1,basis%npgf(iset)
                              cci = basis%gcc(ipgf,ishell,iset)
                              aai = basis%zet(ipgf,iset)
                              norm_i= norm(isgf) 
                              DO jpgf=1,basis%npgf(jset)
                                 ccj = basis%gcc(jpgf,jshell,jset)
                                 aaj = basis%zet(jpgf,jset)
                                 oo = 1._dp/(aai+aaj)**expa 
                                 norm_j= norm(jsgf) 
                                 ovlp(isgf,jsgf) = ovlp(isgf,jsgf) + norm_i*norm_j*ppl*cci*ccj*oo
                              END DO
                           END DO
                         ENDIF
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
    END DO
    
  END SUBROUTINE basis_ovlp

! *****************************************************************************
!> \brief creates and initializes lri_list
! *****************************************************************************
  SUBROUTINE allocate_lri_list_type(lri_list, nkind, error)

    TYPE(lri_list_type), POINTER             :: lri_list
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_list_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure =.FALSE.

    ALLOCATE(lri_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    lri_list%nkind = nkind
    ALLOCATE(lri_list%lri_atom(nkind*nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    DO i=1,nkind*nkind
      NULLIFY(lri_list%lri_atom(i)%lri_node)
      lri_list%lri_atom(i)%natom = 0
    END DO

  END SUBROUTINE allocate_lri_list_type

! *****************************************************************************
!> \brief creates and initializes lri_coefs
! *****************************************************************************
  SUBROUTINE allocate_lri_coefs(lri_env, lri_density, atomic_kind_set, error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_coefs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, ispin, natom, nkind, &
                                                nsgf, nspin, stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: fbas
    TYPE(lri_spin_type), DIMENSION(:), &
      POINTER                                :: lri_coefs

    failure =.FALSE.
    NULLIFY(atomic_kind, fbas, lri_coefs)
    nkind = SIZE(atomic_kind_set)
    nspin = lri_density%nspin  
 
    ALLOCATE(lri_density%lri_coefs(nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    lri_coefs => lri_density%lri_coefs   

    DO ispin =1, nspin
     ALLOCATE(lri_density%lri_coefs(ispin)%lri_kinds(nkind),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ikind=1,nkind
       NULLIFY(lri_coefs(ispin)%lri_kinds(ikind)%acoef)
       NULLIFY(lri_coefs(ispin)%lri_kinds(ikind)%v_int)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
       fbas => lri_env%ri_basis(ikind)%gto_basis_set
       nsgf = fbas%nsgf
       ALLOCATE(lri_coefs(ispin)%lri_kinds(ikind)%acoef(natom,nsgf),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lri_coefs(ispin)%lri_kinds(ikind)%acoef = 0._dp
       ALLOCATE(lri_coefs(ispin)%lri_kinds(ikind)%v_int(natom,nsgf),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lri_coefs(ispin)%lri_kinds(ikind)%v_int = 0._dp
     END DO
    ENDDO

  END SUBROUTINE allocate_lri_coefs

! *****************************************************************************
!> \brief releases the given lri_coefs
!> \param lri_coefs the integral storage environment that is released 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_coefs(lri_coefs, error)
    TYPE(lri_spin_type), DIMENSION(:), &
      POINTER                                :: lri_coefs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_coefs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_coefs)) THEN
       DO i=1,SIZE(lri_coefs)
        DO j=1,SIZE(lri_coefs(i)%lri_kinds)
           IF(ASSOCIATED(lri_coefs(i)%lri_kinds(j)%acoef)) THEN
            DEALLOCATE(lri_coefs(i)%lri_kinds(j)%acoef,stat=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           ENDIF
           IF(ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_int)) THEN
            DEALLOCATE(lri_coefs(i)%lri_kinds(j)%v_int,stat=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           ENDIF
        ENDDO
        DEALLOCATE(lri_coefs(i)%lri_kinds,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END DO
       DEALLOCATE(lri_coefs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_coefs)
 
  END SUBROUTINE deallocate_lri_coefs

! *****************************************************************************
!> \brief deallocates the given lri_ints
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_ints (lri_ints,error)

    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ijkind, inode, istat, &
                                                natom, nkind, nnode
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(lri_ints),cp_failure_level,routineP,error,failure)
    nkind = lri_ints%nkind

    IF(nkind > 0) THEN
       DO ijkind = 1,SIZE(lri_ints%lri_atom)
          natom = lri_ints%lri_atom(ijkind)%natom
          IF(natom > 0) THEN
             DO iatom=1,natom
                nnode = lri_ints%lri_atom(ijkind)%lri_node(iatom)%nnode
                IF(nnode > 0) THEN
                   IF(ASSOCIATED(lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int)) THEN
                      DO inode = 1,nnode
                         DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%abaint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%abbint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%soo,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sab,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dacoef,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sinv,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sn,&
                                     STAT=istat)
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      END DO
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int,STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
             DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF
       END DO
       DEALLOCATE (lri_ints%lri_atom,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (lri_ints,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_lri_ints

! *****************************************************************************
!> \brief deallocates the given lri_ints
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_rho(lri_rho,error)

    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ijkind, inode, istat, &
                                                natom, nkind, nnode
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(lri_rho),cp_failure_level,routineP,error,failure)
    nkind = lri_rho%nkind

    IF(nkind > 0) THEN
       DO ijkind = 1,SIZE(lri_rho%lri_atom)
          natom = lri_rho%lri_atom(ijkind)%natom
          IF(natom > 0) THEN
             DO iatom=1,natom
                nnode = lri_rho%lri_atom(ijkind)%lri_node(iatom)%nnode
                IF(nnode > 0) THEN
                   IF(ASSOCIATED(lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab)) THEN
                     DO inode = 1,nnode
                        DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode)%avec,&
                                    lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab(inode)%tvec,&
                                    STAT=istat)
                        CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                     END DO
                     DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab,STAT=istat)
                     CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
             DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF
       END DO
       DEALLOCATE (lri_rho%lri_atom,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (lri_rho,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_lri_rho

END MODULE lri_environment_types

