!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

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
     ! integral (ai) * sinv * tvec
     REAL(KIND=dp)                                           :: nst
     ! Lagrange parameter
     REAL(KIND=dp)                                           :: lambda
     ! Charge of pair density
     REAL(KIND=dp)                                           :: charge
  END TYPE lri_rhoab_type

! *****************************************************************************

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
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abaint
     ! integrals (a,b,bi)
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: abbint
     ! integrals (da/dA,b,dai/dA)
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER              :: dabdaint
     ! integrals (da/dA,b,bi)
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER              :: dabbint
     ! integrals (a,b)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: soo 
     ! derivative d(a,b)/dA
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: dsoo
     ! integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sab
     ! derivative d(ai,bi)/dA
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: dsab
     ! derivative of fit coeff dacoef/dpmatrix
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: dacoef
     ! inverse of integrals (ai,bi)
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: sinv
     ! integral (ai) / (bi), dim(1..nfa,nfa+1..nfa+nfb)
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: n
     ! sinv * (ai)
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: sn
     ! (ai) * sinv * (ai)
     REAL(KIND=dp)                                           :: nsn
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
     REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER              :: soaabb
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
     TYPE(lri_list_type),POINTER                             :: lri_list
  END TYPE lri_list_p_type
  
! *****************************************************************************

  TYPE lri_bas_type
     ! integral of ri basis fbas
     REAL(KIND=dp), DIMENSION(:), POINTER                    :: int_fbas
     ! self overlap ri basis
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: ri_ovlp
     ! self overlap orb basis
     REAL(KIND=dp), DIMENSION(:,:), POINTER                  :: orb_ovlp
     ! contraction matrix for SHG integrals ri basis
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: scon_ri
     ! contraction matrix for SHG integrals orb basis
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: scon_orb
  END TYPE lri_bas_type

! *****************************************************************************

  TYPE lri_clebsch_gordon_type
     ! Clebsch-Gordon (CG) coefficients
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER                :: cg_coeff
     ! list of non-zero CG coefficients
     INTEGER, DIMENSION(:,:,:), POINTER                      :: cg_none0_list
     ! number of non-zero CG coefficients
     INTEGER, DIMENSION(:,:), POINTER                        :: ncg_none0
  END TYPE lri_clebsch_gordon_type

! *****************************************************************************

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
     ! orbital basis set
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: orb_basis
     ! lri (fit) basis set 
     TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER       :: ri_basis
     ! orb_basis neighborlist, LRI integrals
     TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER   :: soo_list
     ! local RI integrals
     TYPE(lri_list_type), POINTER                            :: lri_ints
     ! local integral of rho**2; for optimization 
     TYPE(lri_list_type), POINTER                            :: lri_ints_rho
     ! properties of orb and aux basis
     TYPE(lri_bas_type), DIMENSION(:), POINTER               :: bas_prop
     ! Clebsch-Gordon for solid harmonics
     TYPE(lri_clebsch_gordon_type), POINTER                  :: cg_shg
  END TYPE lri_environment_type

! *****************************************************************************

  TYPE lri_kind_type
    ! expansion coeff for lri density dim(natom,nsgf)
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: acoef
    ! integrals V*fbas (potential*fit basis) dim(natom,nsgf)
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: v_int
    ! SUM_i integral(V*fbas_i)*davec/dR dim(natom,3) 
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: v_dadr
    ! integrals V*dfbas/dR
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: v_dfdr
  END TYPE lri_kind_type

  TYPE lri_spin_type
     TYPE(lri_kind_type), DIMENSION(:), POINTER              :: lri_kinds
  END TYPE lri_spin_type

! *****************************************************************************

  TYPE lri_force_type
    REAL(KIND=dp), DIMENSION(:), POINTER                     :: st
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: dssn,&
                                                                sdssn,&
                                                                dsst,&
                                                                sdsst,&
                                                                sdt
    ! derivative dtvec/dR 
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: dtvec
    ! derivative davec/dR 
    REAL(KIND=dp), DIMENSION(:,:), POINTER                   :: davec
    ! derivative overlap matrix dS/dR 
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER                 :: ds
  END TYPE lri_force_type

! *****************************************************************************

  TYPE lri_density_type
     INTEGER :: id_nr, ref_count, in_use
     INTEGER                                                 :: nspin
     ! pair density expansion (nspin)
     TYPE(lri_list_p_type), DIMENSION(:),POINTER             :: lri_rhos
     ! coefficients of RI expansion and gradients (nspin)
     TYPE(lri_spin_type), DIMENSION(:),POINTER               :: lri_coefs
     TYPE(lri_force_type), POINTER                           :: lri_force
  END TYPE lri_density_type

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_types'
  INTEGER, SAVE, PRIVATE :: last_lri_env_id=0
  INTEGER, SAVE, PRIVATE :: last_lri_density_id=0

  PUBLIC :: lri_environment_type,&
            lri_force_type, lri_list_type,&
            lri_int_type, lri_int_rho_type, lri_density_type,&
            lri_kind_type, lri_rhoab_type
  PUBLIC :: lri_env_create, lri_env_release, allocate_lri_coefs,&
            allocate_lri_ints, allocate_lri_ints_rho, lri_density_create,&
            lri_density_release,allocate_lri_rhos, allocate_lri_force_components,&
            deallocate_lri_ints, deallocate_lri_ints_rho,&
            deallocate_lri_force_components, deallocate_bas_properties

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

    INTEGER                                  :: stat

    ALLOCATE ( lri_env, stat=stat )
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    last_lri_env_id=last_lri_env_id+1
    lri_env%id_nr=last_lri_env_id
    lri_env%ref_count=1
    lri_env%in_use=0

    lri_env%debug = .FALSE.
    lri_env%delta = 1.E-6_dp

    NULLIFY(lri_env%orb_basis)
    NULLIFY(lri_env%ri_basis)

    NULLIFY(lri_env%soo_list)
    NULLIFY(lri_env%lri_ints)
    NULLIFY(lri_env%lri_ints_rho)
    NULLIFY(lri_env%bas_prop)
    NULLIFY(lri_env%cg_shg)

    ALLOCATE ( lri_env%cg_shg, stat=stat )
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

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
       IF (ASSOCIATED(lri_env%lri_ints_rho)) THEN
          CALL deallocate_lri_ints_rho(lri_env%lri_ints_rho,error)
       END IF
       CALL deallocate_bas_properties(lri_env,error)
       IF (ASSOCIATED(lri_env%cg_shg)) THEN
          DEALLOCATE(lri_env%cg_shg%cg_coeff,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%cg_shg%cg_none0_list,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%cg_shg%ncg_none0,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%cg_shg,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
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

    INTEGER                                  :: stat

    ALLOCATE ( lri_density, stat=stat )
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    last_lri_density_id=last_lri_density_id+1
    lri_density%id_nr=last_lri_density_id
    lri_density%ref_count=1
    lri_density%in_use=0

    lri_density%nspin=0
  
    NULLIFY(lri_density%lri_rhos)
    NULLIFY(lri_density%lri_coefs)
    NULLIFY(lri_density%lri_force)

  END SUBROUTINE lri_density_create

! *****************************************************************************
!> \brief releases the given lri_density
!> \param lri_density the lri_density to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_density_release(lri_density, error)
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_density_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_density)) THEN
       lri_density%ref_count=0

       CALL deallocate_lri_rhos(lri_density%lri_rhos,error)
       CALL deallocate_lri_coefs(lri_density%lri_coefs, error)
       CALL deallocate_lri_force_components(lri_density%lri_force, error)

       DEALLOCATE(lri_density,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(lri_density)

  END SUBROUTINE lri_density_release

! *****************************************************************************
!> \brief allocate lri_ints, matrices that store LRI integrals
!> \param lri_env ...
!> \param lri_ints structure storing the LRI integrals
!> \param nkind number of atom kinds 
!> \param calculate_forces ...
!> \param virial ...
!> \param error variable to control error logging, stopping,...
! *****************************************************************************
  SUBROUTINE allocate_lri_ints(lri_env,lri_ints,nkind,calculate_forces,&
                               virial,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_list_type), POINTER             :: lri_ints
    INTEGER, INTENT(IN)                      :: nkind
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(virial_type), POINTER               :: virial
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_ints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iac, iatom, ikind, ilist, &
                                                jatom, jkind, jneighbor, nba, &
                                                nbb, nfa, nfb, nlist, nn, &
                                                nneighbor, stat
    LOGICAL                                  :: failure, use_virial
    REAL(KIND=dp)                            :: dab, rab(3)
    TYPE(gto_basis_set_type), POINTER        :: fbasa, fbasb, obasa, obasb
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    failure =.FALSE.
    NULLIFY(fbasa, fbasb, lrii, nl_iterator, obasa, obasb)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    ALLOCATE(lri_ints,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    lri_ints%nkind = nkind
    ALLOCATE(lri_ints%lri_atom(nkind*nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    DO i=1,nkind*nkind
      NULLIFY(lri_ints%lri_atom(i)%lri_node)
      lri_ints%lri_atom(i)%natom = 0
    END DO

    CALL neighbor_list_iterator_create(nl_iterator,lri_env%soo_list)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)

       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
            iatom=iatom,jatom=jatom,r=rab)
      
       iac = ikind + nkind*(jkind - 1)
       dab = SQRT(SUM(rab*rab))
 
       obasa => lri_env%orb_basis(ikind)%gto_basis_set
       obasb => lri_env%orb_basis(jkind)%gto_basis_set
       fbasa => lri_env%ri_basis(ikind)%gto_basis_set
       fbasb => lri_env%ri_basis(jkind)%gto_basis_set

       IF (.NOT.ASSOCIATED(obasa)) CYCLE
       IF (.NOT.ASSOCIATED(obasb)) CYCLE

       IF(.NOT.ASSOCIATED(lri_ints%lri_atom(iac)%lri_node)) THEN
          lri_ints%lri_atom(iac)%natom = nlist
          ALLOCATE(lri_ints%lri_atom(iac)%lri_node(nlist),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,nlist
             NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%lri_int)
             lri_ints%lri_atom(iac)%lri_node(i)%nnode = 0
          END DO
       END IF
       IF(.NOT.ASSOCIATED(lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int)) THEN
          lri_ints%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
          ALLOCATE(lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(nneighbor),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

       nba = obasa%nsgf
       nbb = obasb%nsgf
       nfa = fbasa%nsgf
       nfb = fbasb%nsgf
       nn  = nfa + nfb
       
       ALLOCATE(lrii%abaint(nba,nbb,nfa),&
                lrii%abbint(nba,nbb,nfb),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lrii%abaint = 0._dp
       lrii%abbint = 0._dp

       ALLOCATE(lrii%dabdaint(nba,nbb,nfa,3),&
                lrii%dabbint(nba,nbb,nfb,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lrii%dabdaint = 0._dp
       lrii%dabbint = 0._dp

       ALLOCATE(lrii%sab(nfa,nfb),lrii%dsab(nfa,nfb,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lrii%sab  = 0._dp
       lrii%dsab = 0._dp

       IF(iatom == jatom.AND.dab < lri_env%delta) THEN
        ALLOCATE(lrii%sinv(nfa,nfa),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
        ALLOCATE(lrii%sinv(nn,nn),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       lrii%sinv = 0._dp

       IF(iatom == jatom.AND.dab < lri_env%delta) THEN
        ALLOCATE(lrii%n(nfa),lrii%sn(nfa),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
        ALLOCATE(lrii%n(nn),lrii%sn(nn),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       lrii%n    = 0._dp
       lrii%sn   = 0._dp

       ALLOCATE(lrii%soo(nba,nbb),lrii%dsoo(nba,nbb,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lrii%soo  = 0._dp
       lrii%dsoo = 0._dp

       IF(iatom == jatom.AND.dab < lri_env%delta) THEN
        ALLOCATE(lrii%dacoef(nba,nbb,nfa),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
        ALLOCATE(lrii%dacoef(nba,nbb,nn),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
       lrii%dacoef = 0._dp
 
       lrii%dmax_ab   = 0._dp
       lrii%dmax_oo   = 0._dp
       lrii%dmax_aba  = 0._dp
       lrii%dmax_abb  = 0._dp

       lrii%calc_force_pair = .FALSE.

       ! forces: not every pair is giving contributions
       IF(calculate_forces) THEN
         !no forces for self-pair aa
         IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           lrii%calc_force_pair = .FALSE.
         ELSE
           !forces for periodic self-pair aa' required for virial 
           IF(iatom ==  jatom.AND..NOT.use_virial) THEN
             lrii%calc_force_pair = .FALSE.       
           ELSE
             lrii%calc_force_pair = .TRUE.
           ENDIF 
         ENDIF
       ENDIF

    ENDDO

    CALL neighbor_list_iterator_release(nl_iterator)

  END SUBROUTINE allocate_lri_ints

! *****************************************************************************
!> \brief allocate lri_ints_rho, storing integral for the exact density
!> \param lri_env ...
!> \param lri_ints_rho structure storing the integrals (aa,bb)
!> \param nkind number of atom kinds 
!> \param error variable to control error logging, stopping,...
! *****************************************************************************
  SUBROUTINE allocate_lri_ints_rho(lri_env,lri_ints_rho,nkind,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_list_type), POINTER             :: lri_ints_rho
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_ints_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iac, iatom, ikind, ilist, &
                                                jatom, jkind, jneighbor, nba, &
                                                nbb, nlist, nneighbor, stat
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: obasa, obasb
    TYPE(lri_int_rho_type), POINTER          :: lriir
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    ALLOCATE(lri_ints_rho,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    lri_ints_rho%nkind = nkind
    ALLOCATE(lri_ints_rho%lri_atom(nkind*nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    DO i=1,nkind*nkind
      NULLIFY(lri_ints_rho%lri_atom(i)%lri_node)
      lri_ints_rho%lri_atom(i)%natom = 0
    ENDDO

    CALL neighbor_list_iterator_create(nl_iterator,lri_env%soo_list)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)

       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
            nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
            iatom=iatom,jatom=jatom)
      
       iac = ikind + nkind*(jkind - 1)

       obasa => lri_env%orb_basis(ikind)%gto_basis_set
       obasb => lri_env%orb_basis(jkind)%gto_basis_set

       IF (.NOT.ASSOCIATED(obasa)) CYCLE
       IF (.NOT.ASSOCIATED(obasb)) CYCLE

       IF(.NOT.ASSOCIATED(lri_ints_rho%lri_atom(iac)%lri_node)) THEN
          lri_ints_rho%lri_atom(iac)%natom = nlist
          ALLOCATE(lri_ints_rho%lri_atom(iac)%lri_node(nlist),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,nlist
             NULLIFY(lri_ints_rho%lri_atom(iac)%lri_node(i)%lri_int_rho)
             lri_ints_rho%lri_atom(iac)%lri_node(i)%nnode = 0
          END DO
       END IF
       IF(.NOT.ASSOCIATED(lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho)) THEN
          lri_ints_rho%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
          ALLOCATE(lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(nneighbor),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       lriir => lri_ints_rho%lri_atom(iac)%lri_node(ilist)%lri_int_rho(jneighbor)

       nba = obasa%nsgf
       nbb = obasb%nsgf

       ALLOCATE(lriir%soaabb(nba,nba,nbb,nbb),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lriir%soaabb = 0._dp
       lriir%dmax_aabb = 0._dp

    ENDDO

    CALL neighbor_list_iterator_release(nl_iterator)
  
  END SUBROUTINE allocate_lri_ints_rho

! *****************************************************************************
!> \brief creates and initializes lri_rhos
!> \param lri_env ...
!> \param lri_rhos structure storing tvec and avec
!> \param nspin ...
!> \param nkind number of atom kinds 
!> \param error variable to control error logging, stopping,...
! *****************************************************************************
  SUBROUTINE allocate_lri_rhos(lri_env,lri_rhos,nspin,nkind,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_list_p_type), DIMENSION(:), &
      POINTER                                :: lri_rhos
    INTEGER, INTENT(IN)                      :: nspin, nkind
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_lri_rhos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iac, iatom, ikind, ilist, &
                                                ispin, jatom, jkind, &
                                                jneighbor, nfa, nfb, nlist, &
                                                nn, nneighbor, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, rab(3)
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    failure =.FALSE.
    NULLIFY(lri_rho, lrho, lrii, nl_iterator)

    ALLOCATE(lri_rhos(nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1, nspin

       ALLOCATE(lri_rhos(ispin)%lri_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       lri_rhos(ispin)%lri_list%nkind = nkind
       ALLOCATE(lri_rhos(ispin)%lri_list%lri_atom(nkind*nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       DO i=1,nkind*nkind
         NULLIFY(lri_rhos(ispin)%lri_list%lri_atom(i)%lri_node)
         lri_rhos(ispin)%lri_list%lri_atom(i)%natom = 0
       END DO

       lri_rho => lri_rhos(ispin)%lri_list

       CALL neighbor_list_iterator_create(nl_iterator,lri_env%soo_list)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
               iatom=iatom,jatom=jatom,nlist=nlist,ilist=ilist,&
               nnode=nneighbor,inode=jneighbor,r=rab)
  
          iac = ikind + nkind*(jkind - 1)
          dab = SQRT(SUM(rab*rab))

          IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

          IF(.NOT.ASSOCIATED(lri_rho%lri_atom(iac)%lri_node)) THEN
             lri_rho%lri_atom(iac)%natom = nlist
             ALLOCATE(lri_rho%lri_atom(iac)%lri_node(nlist),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO i=1,nlist
                NULLIFY(lri_rho%lri_atom(iac)%lri_node(i)%lri_rhoab)
                lri_rho%lri_atom(iac)%lri_node(i)%nnode = 0
             END DO
          END IF
          IF(.NOT.ASSOCIATED(lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab)) THEN
             lri_rho%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
             ALLOCATE(lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(nneighbor),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF

          lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
          lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

          lrho%nba = lrii%nba
          lrho%nbb = lrii%nbb
          lrho%nfa = lrii%nfa
          lrho%nfb = lrii%nfb
 
          nfa = lrho%nfa
          nfb = lrho%nfb
          nn  = nfa + nfb

          IF(iatom == jatom.AND.dab < lri_env%delta) THEN
           ALLOCATE(lrho%avec(nfa), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(lrho%tvec(nfa), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
           ALLOCATE(lrho%avec(nn), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(lrho%tvec(nn), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          lrho%avec = 0._dp
          lrho%tvec = 0._dp
          
       ENDDO

       CALL neighbor_list_iterator_release(nl_iterator)

    ENDDO

  END SUBROUTINE allocate_lri_rhos

! *****************************************************************************
!> \brief creates and initializes lri_coefs
!> \param lri_env ...
!> \param lri_density ...
!> \param atomic_kind_set ...
!> \param error variable to control error logging, stopping,...
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
       NULLIFY(lri_coefs(ispin)%lri_kinds(ikind)%v_dadr)
       NULLIFY(lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr)
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
       ALLOCATE(lri_coefs(ispin)%lri_kinds(ikind)%v_dadr(natom,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lri_coefs(ispin)%lri_kinds(ikind)%v_dadr = 0._dp
       ALLOCATE(lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr(natom,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       lri_coefs(ispin)%lri_kinds(ikind)%v_dfdr = 0._dp
     END DO
    ENDDO

  END SUBROUTINE allocate_lri_coefs

! *****************************************************************************
!> \brief creates and initializes lri_force
!> \param lri_force ...
!> \param nfa and nfb number of fit functions on a/b
!> \param nfb ...
!> \param error variable to control error logging, stopping,...
! *****************************************************************************
  SUBROUTINE allocate_lri_force_components(lri_force,nfa,nfb,error)

    TYPE(lri_force_type), POINTER            :: lri_force
    INTEGER, INTENT(IN)                      :: nfa, nfb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'allocate_lri_force_components', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nn, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    nn = nfa + nfb

    IF(.NOT.ASSOCIATED(lri_force)) THEN
      ALLOCATE(lri_force,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE(lri_force%ds(nn,nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%ds    = 0._dp
      ALLOCATE(lri_force%st(nn),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%st    = 0._dp
      ALLOCATE(lri_force%dsst(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%dsst  = 0._dp
      ALLOCATE(lri_force%sdsst(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%sdsst = 0._dp
      ALLOCATE(lri_force%dssn(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%dssn  = 0._dp
      ALLOCATE(lri_force%sdssn(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%sdssn = 0._dp
      ALLOCATE(lri_force%sdt(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%sdt   = 0._dp
      ALLOCATE(lri_force%davec(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%davec = 0._dp
      ALLOCATE(lri_force%dtvec(nn,3),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_force%dtvec = 0._dp
    ENDIF

  END SUBROUTINE allocate_lri_force_components
 
! *****************************************************************************
!> \brief deallocates one-center overlap integrals, integral of ri basis
!>        and scon matrices
!> \param lri_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_bas_properties (lri_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_bas_properties', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat

    IF (ASSOCIATED(lri_env%bas_prop)) THEN
       DO i=1,SIZE(lri_env%bas_prop)
          DEALLOCATE(lri_env%bas_prop(i)%int_fbas,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%bas_prop(i)%ri_ovlp,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%bas_prop(i)%orb_ovlp,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%bas_prop(i)%scon_ri,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(lri_env%bas_prop(i)%scon_orb,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END DO
       DEALLOCATE(lri_env%bas_prop,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

  END SUBROUTINE deallocate_bas_properties

! *****************************************************************************
!> \brief deallocates the given lri_ints
!> \param lri_ints ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_ints (lri_ints,error)

    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ijkind, inode, natom, &
                                                nkind, nnode, stat
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
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dabdaint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dabbint,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%soo,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dsoo,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sab,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dsab,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%dacoef,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sinv,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%n,&
                                     lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int(inode)%sn,&
                                     STAT=stat)
                         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                      END DO
                      DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node(iatom)%lri_int,STAT=stat)
                      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
             DEALLOCATE (lri_ints%lri_atom(ijkind)%lri_node,STAT=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF
       END DO
       DEALLOCATE (lri_ints%lri_atom,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (lri_ints,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_lri_ints

! *****************************************************************************
!> \brief deallocates the given lri_ints_rho
!> \param lri_ints_rho ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_ints_rho (lri_ints_rho,error)

    TYPE(lri_list_type), POINTER             :: lri_ints_rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_ints_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ijkind, inode, natom, &
                                                nkind, nnode, stat
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(lri_ints_rho),cp_failure_level,routineP,error,failure)
    nkind = lri_ints_rho%nkind

    IF(nkind > 0) THEN
       DO ijkind = 1,SIZE(lri_ints_rho%lri_atom)
          natom = lri_ints_rho%lri_atom(ijkind)%natom
          IF(natom > 0) THEN
             DO iatom=1,natom
                nnode = lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%nnode
                IF(nnode > 0) THEN
                   IF(ASSOCIATED(lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho)) THEN
                      DO inode = 1,nnode
                         DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho(inode)%soaabb,&
                                     STAT=stat)
                         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                      END DO
                      DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node(iatom)%lri_int_rho,STAT=stat)
                      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
             DEALLOCATE (lri_ints_rho%lri_atom(ijkind)%lri_node,STAT=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
           END IF
       END DO
       DEALLOCATE (lri_ints_rho%lri_atom,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE (lri_ints_rho,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_lri_ints_rho

! *****************************************************************************
!> \brief deallocates the given lri_rhos
!> \param lri_rhos ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_rhos(lri_rhos,error)

    TYPE(lri_list_p_type), DIMENSION(:), &
      POINTER                                :: lri_rhos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_lri_rhos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ijkind, inode, &
                                                natom, nkind, nnode, stat
    LOGICAL                                  :: failure
    TYPE(lri_list_type), POINTER             :: lri_rho

    failure = .FALSE.
    NULLIFY(lri_rho)

    IF(ASSOCIATED(lri_rhos)) THEN

      DO i=1,SIZE(lri_rhos)

         lri_rho => lri_rhos(i)%lri_list
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
                                      STAT=stat)
                          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                       END DO
                       DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node(iatom)%lri_rhoab,STAT=stat)
                       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                     END IF
                   END IF
                END DO
                DEALLOCATE (lri_rho%lri_atom(ijkind)%lri_node,STAT=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
               END IF
           END DO
           DEALLOCATE (lri_rho%lri_atom,STAT=stat)
           CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
         END IF
         DEALLOCATE (lri_rho,STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO

      DEALLOCATE(lri_rhos,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    END IF

  END SUBROUTINE deallocate_lri_rhos

! *****************************************************************************
!> \brief releases the given lri_coefs
!> \param lri_coefs the integral storage environment that is released 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_coefs(lri_coefs,error)
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
           IF(ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dadr)) THEN
            DEALLOCATE(lri_coefs(i)%lri_kinds(j)%v_dadr,stat=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           ENDIF
           IF(ASSOCIATED(lri_coefs(i)%lri_kinds(j)%v_dfdr)) THEN
            DEALLOCATE(lri_coefs(i)%lri_kinds(j)%v_dfdr,stat=stat)
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
!> \brief releases the given lri_force_type
!> \param lri_force the integral storage environment that is released 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE deallocate_lri_force_components(lri_force,error)

    TYPE(lri_force_type), POINTER            :: lri_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'deallocate_lri_force_components', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(lri_force)) THEN

       IF(ASSOCIATED(lri_force%st)) THEN
        DEALLOCATE(lri_force%st,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%dssn)) THEN
        DEALLOCATE(lri_force%dssn,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%sdssn)) THEN
        DEALLOCATE(lri_force%sdssn,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%dsst)) THEN
        DEALLOCATE(lri_force%dsst,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%sdsst)) THEN
        DEALLOCATE(lri_force%sdsst,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%sdt)) THEN
        DEALLOCATE(lri_force%sdt,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%dtvec)) THEN
        DEALLOCATE(lri_force%dtvec,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%davec)) THEN
        DEALLOCATE(lri_force%davec,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF(ASSOCIATED(lri_force%ds)) THEN
        DEALLOCATE(lri_force%ds,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       
       DEALLOCATE(lri_force, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       NULLIFY(lri_force)
    ENDIF

  END SUBROUTINE deallocate_lri_force_components

END MODULE lri_environment_types

