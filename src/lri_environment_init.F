!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initializes the environment for lri
!>        lri : local resolution of the identity
!> \par History
!>      created [06.2015]
!> \author Dorothea Golze  
! *****************************************************************************
MODULE lri_environment_init
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                             gto_basis_set_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: deallocate_bas_properties,&
                                             lri_env_create,&
                                             lri_environment_type
  USE mathconstants,                   ONLY: fac,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_deallocate,&
                                             clebsch_gordon_init
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_init'

  PUBLIC :: lri_env_init, lri_basis_init 


! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief initializes the lri env
!> \param qs_env ...
!> \param force_env_section ...
!> \param qs_kind_set ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE lri_env_init(qs_env,force_env_section,qs_kind_set,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, lmax_ikind_orb, &
                                                lmax_ikind_ri, maxl_orb, &
                                                maxl_ri, nkind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gto_basis_set_type), POINTER        :: lri_basis_set, orb_basis_set
    TYPE(lri_environment_type), POINTER      :: lri_env

    NULLIFY(atomic_kind_set,lri_basis_set,lri_env,orb_basis_set)
    CALL lri_env_create(lri_env,error)

    ! initialize the basic basis sets (orb and ri)
    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
    nkind = SIZE(atomic_kind_set)
    ALLOCATE(lri_env%orb_basis(nkind),lri_env%ri_basis(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    maxl_orb = 0
    maxl_ri = 0
    DO ikind=1,nkind
       NULLIFY(orb_basis_set,lri_basis_set)
       CALL get_qs_kind(qs_kind_set(ikind),&
                        orb_basis_set=orb_basis_set,&
                        lri_basis_set=lri_basis_set)
       NULLIFY(lri_env%orb_basis(ikind)%gto_basis_set)
       NULLIFY(lri_env%ri_basis(ikind)%gto_basis_set)
       IF (ASSOCIATED(orb_basis_set)) THEN
          CALL copy_gto_basis_set(orb_basis_set,lri_env%orb_basis(ikind)%gto_basis_set,error)
          CALL copy_gto_basis_set(lri_basis_set,lri_env%ri_basis(ikind)%gto_basis_set,error)
       END IF
       lmax_ikind_orb = MAXVAL(lri_env%orb_basis(ikind)%gto_basis_set%lmax)
       lmax_ikind_ri  = MAXVAL(lri_env%ri_basis(ikind)%gto_basis_set%lmax)
       maxl_orb = MAX(maxl_orb,lmax_ikind_orb)
       maxl_ri  = MAX(maxl_ri,lmax_ikind_ri)
    END DO
    CALL lri_basis_init(lri_env,error)

    ! CG coefficients needed for lri integrals
    IF(ASSOCIATED(lri_env%cg_shg)) THEN
      CALL get_clebsch_gordon_coefficients(lri_env%cg_shg%cg_coeff,&
                                           lri_env%cg_shg%cg_none0_list,&
                                           lri_env%cg_shg%ncg_none0,&
                                           maxl_orb,maxl_ri,error)
    ENDIF
    ! init keywords

    ! check for debug (OS scheme)
    CALL section_vals_val_get(force_env_section,"DFT%QS%LRIGPW%DEBUG_LRI_INTEGRALS",&
                              l_val=lri_env%debug,error=error)
    ! integrals based on solid harmonic Gaussians
    CALL section_vals_val_get(force_env_section,"DFT%QS%LRIGPW%SHG_LRI_INTEGRALS",&
                              l_val=lri_env%use_shg_integrals,error=error)
    ! how to calculate inverse/pseuodinverse of overlap
    CALL section_vals_val_get(force_env_section,"DFT%QS%LRIGPW%LRI_OVERLAP_MATRIX",&
                              i_val=lri_env%lri_overlap_inv,error=error)
    CALL section_vals_val_get(force_env_section,"DFT%QS%LRIGPW%MAX_CONDITION_NUM",&
                              r_val=lri_env%cond_max,error=error)
    !
    CALL set_qs_env (qs_env,lri_env=lri_env,error=error)

  END SUBROUTINE lri_env_init

! *****************************************************************************
!> \brief initializes the lri basis: calculates the norm, self-overlap
!>        and integral of the ri basis
!> \param lri_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_basis_init(lri_env,error)
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lri_basis_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: orb_norm_r, ri_int_fbas, &
                                                ri_norm_r, ri_norm_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: orb_ovlp, ri_ovlp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scon_orb, scon_ri
    TYPE(gto_basis_set_type), POINTER        :: orb_basis, ri_basis

    failure=.FALSE.

    IF (ASSOCIATED(lri_env)) THEN
       IF (ASSOCIATED(lri_env%orb_basis)) THEN
          CPPreconditionNoFail(ASSOCIATED(lri_env%ri_basis),cp_warning_level,routineP,error)
          nkind = SIZE(lri_env%orb_basis)
          CALL deallocate_bas_properties(lri_env,error)
          ALLOCATE (lri_env%bas_prop(nkind),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ikind = 1, nkind
             NULLIFY(orb_basis, ri_basis)
             orb_basis => lri_env%orb_basis(ikind)%gto_basis_set
             IF (ASSOCIATED(orb_basis)) THEN
                ri_basis => lri_env%ri_basis(ikind)%gto_basis_set
                CPPreconditionNoFail(ASSOCIATED(ri_basis),cp_warning_level,routineP,error)
                NULLIFY(ri_norm_r)
                CALL basis_norm_radial(ri_basis,ri_norm_r,error)
                NULLIFY(orb_norm_r)
                CALL basis_norm_radial(orb_basis,orb_norm_r,error)
                NULLIFY(ri_norm_s)
                CALL basis_norm_s_func(ri_basis,ri_norm_s,error)
                NULLIFY(ri_int_fbas)
                CALL basis_int(ri_basis,ri_int_fbas,ri_norm_s,error)
                lri_env%bas_prop(ikind)%int_fbas => ri_int_fbas
                NULLIFY(ri_ovlp)
                CALL basis_ovlp(ri_basis,ri_ovlp,ri_norm_r,error)
                lri_env%bas_prop(ikind)%ri_ovlp => ri_ovlp
                NULLIFY(orb_ovlp)
                CALL basis_ovlp(orb_basis,orb_ovlp,orb_norm_r,error)
                lri_env%bas_prop(ikind)%orb_ovlp => orb_ovlp
                NULLIFY(scon_ri)
                CALL contraction_matrix_shg(ri_basis,scon_ri,error)
                lri_env%bas_prop(ikind)%scon_ri => scon_ri
                NULLIFY(scon_orb)
                CALL contraction_matrix_shg(orb_basis,scon_orb,error)
                lri_env%bas_prop(ikind)%scon_orb => scon_orb
                DEALLOCATE (orb_norm_r,ri_norm_r,ri_norm_s,STAT=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE lri_basis_init

!*****************************************************************************
!> \brief normalization for a contracted Gaussian s-function,
!>        spherical = cartesian Gaussian for s-functions
!> \param basis ...
!> \param norm ...
!> \param error ...
! *****************************************************************************
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
!> \param basis ...
!> \param norm ...
!> \param error ...
! *****************************************************************************
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
!> \brief normalization solid harmonic Gaussians (SHG)
!> \param basis ...
!> \param norm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE basis_norm_shg(basis,norm,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(dp), DIMENSION(:, :), POINTER       :: norm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_shg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, ishell, jpgf, l
    REAL(KIND=dp)                            :: aai, aaj, cci, ccj, expa, ppl

    norm = 0._dp

    DO iset=1,basis%nset
       DO ishell=1,basis%nshell(iset)
          l = basis%l(ishell,iset)
          expa = 0.5_dp*REAL(2*l + 3,dp)
          ppl = fac(2*l+2)*pi**(1.5_dp)/fac(l+1)
          ppl = ppl/(2._dp**REAL(2*l+1,dp))
          ppl = ppl/REAL(2*l+1,dp)
          DO ipgf=1,basis%npgf(iset)
             cci = basis%gcc(ipgf,ishell,iset)
             aai = basis%zet(ipgf,iset)
             DO jpgf=1,basis%npgf(iset)
                ccj = basis%gcc(jpgf,ishell,iset)
                aaj = basis%zet(jpgf,iset)
                norm(iset,ishell) =norm(iset,ishell) + cci*ccj*ppl/(aai+aaj)**expa
             END DO
          END DO
          norm(iset,ishell)=1.0_dp/SQRT(norm(iset,ishell)) 
       END DO
    END DO

  END SUBROUTINE basis_norm_shg

!*****************************************************************************
!> \brief integral over a single (contracted) lri auxiliary basis function,
!>        integral is zero for all but s-functions               
!> \param basis ...
!> \param int_aux ...
!> \param norm ...
!> \param error ...
! *****************************************************************************
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
!> \param basis ...
!> \param ovlp ...
!> \param norm ...
!> \param error ...
! *****************************************************************************
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

!*****************************************************************************
!> \brief contraction matrix for SHG integrals
!> \param basis ...
!> \param scon_shg contraction matrix
!> \param error ...
! *****************************************************************************
  SUBROUTINE contraction_matrix_shg(basis,scon_shg,error)

    TYPE(gto_basis_set_type), POINTER        :: basis
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scon_shg
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contraction_matrix_shg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, iset, ishell, l, &
                                                maxpgf, maxshell, nset, stat
    INTEGER, DIMENSION(:), POINTER           :: npgf, nshell
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: gcc, norm_w
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: norm, zet

    nset        =   basis%nset
    npgf        =>  basis%npgf
    nshell      =>  basis%nshell
    zet         =>  basis%zet

    maxpgf   = SIZE(basis%gcc,1)
    maxshell = SIZE(basis%gcc,2)
    ALLOCATE(norm(basis%nset,maxshell),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(scon_shg(maxpgf,maxshell,nset),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL basis_norm_shg(basis,norm,error)

    DO iset=1,nset
       DO ishell=1,nshell(iset)
          l = basis%l(ishell,iset)
          norm_w = norm(iset,ishell)
          DO ipgf = 1,npgf(iset)
             gcc  = basis%gcc(ipgf,ishell,iset)
             scon_shg(ipgf,ishell,iset) = norm(iset,ishell) * gcc 
          END DO
       END DO
    END DO
   
   DEALLOCATE(norm,STAT=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE contraction_matrix_shg

!*****************************************************************************
!> \brief calculate the Clebsch-Gordon (CG) coefficients for expansion of the
!>        product of two spherical harmonic Gaussians
!> \param my_cg matrix storing CG coefficients
!> \param cg_none0_list list of none-zero CG coefficients
!> \param ncg_none0 number of none-zero CG coefficients
!> \param maxl1 maximal l quantum number of 1st spherical function
!> \param maxl2 maximal l quantum number of 2nd spherical function
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_clebsch_gordon_coefficients(my_cg,cg_none0_list,ncg_none0,&
                                             maxl1,maxl2,error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: my_cg
    INTEGER, DIMENSION(:, :, :), POINTER     :: cg_none0_list
    INTEGER, DIMENSION(:, :), POINTER        :: ncg_none0
    INTEGER, INTENT(IN)                      :: maxl1, maxl2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'get_clebsch_gordon_coefficients', &
      routineP = moduleN//':'//routineN

    INTEGER :: il, ilist, iso, iso1, iso2, l1, l1l2, l2, lc1, lc2, lp, m1, &
      m2, maxl, mm, mp, nlist, nlist_max, nsfunc, nsfunc1, nsfunc2, stat
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rga

    failure = .FALSE.

    nlist_max = 6
    nsfunc1 = nsoset(maxl1)
    nsfunc2 = nsoset(maxl2)
    maxl = maxl1+maxl2
    nsfunc  = nsoset(maxl)

    CALL clebsch_gordon_init(maxl)
   
    ALLOCATE(my_cg(nsfunc1,nsfunc2,nsfunc),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    my_cg =0.0_dp
    ALLOCATE(ncg_none0(nsfunc1,nsfunc2),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ncg_none0 = 0
    ALLOCATE(cg_none0_list(nsfunc1,nsfunc2,nlist_max),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    cg_none0_list = 0

    ALLOCATE(rga(maxl,2),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    rga = 0.0_dp
    DO lc1 = 0,maxl1
       DO iso1=nsoset(lc1-1)+1,nsoset(lc1)
          l1 = indso(1,iso1)
          m1 = indso(2,iso1)
          DO lc2 = 0,maxl2
             DO iso2 = nsoset(lc2-1)+1,nsoset(lc2)
                l2 = indso(1,iso2)
                m2 = indso(2,iso2)
                CALL clebsch_gordon (l1,m1,l2,m2,rga)
                l1l2 = l1+l2
                mp = m1 + m2
                mm = m1 - m2
                IF ( m1*m2 < 0  .OR. (m1*m2==0 .AND. (m1<0 .OR. m2<0))) THEN
                  mp = -ABS(mp)
                  mm = -ABS(mm)
                ELSE
                  mp = ABS(mp)
                  mm = ABS(mm)
                END IF
                DO  lp=MOD(l1+l2,2),l1l2,2
                    il = lp/2 + 1
                    IF( ABS(mp) <= lp) THEN
                      IF(mp >= 0) THEN
                        iso = nsoset(lp-1) + lp + 1 + mp
                      ELSE
                        iso = nsoset(lp-1) + lp + 1 - ABS(mp)
                      END IF
                      my_cg(iso1,iso2,iso) =  rga(il,1)
                    ENDIF
                    IF(mp /= mm .AND. ABS(mm) <= lp) THEN
                       IF(mm >= 0) THEN
                         iso = nsoset(lp-1) + lp + 1 + mm
                       ELSE
                          iso = nsoset(lp-1) + lp + 1 - ABS(mm)
                       END IF
                       my_cg(iso1,iso2,iso) =  rga(il,2)
                    ENDIF
                END DO
                nlist = 0
                DO ilist=1,nsfunc
                  IF(ABS(my_cg(iso1,iso2,ilist)) > 1.E-8_dp) THEN
                   nlist = nlist + 1
                   IF(nlist > nlist_max) THEN
                     CALL reallocate(cg_none0_list,1,nsfunc1,1,nsfunc2,1,nlist)
                     nlist_max = nlist 
                   ENDIF
                   cg_none0_list(iso1,iso2,nlist) = ilist
                  ENDIF
                ENDDO
                ncg_none0(iso1,iso2) = nlist
             ENDDO  ! iso2
          ENDDO  ! lc2
       ENDDO    ! iso1
    ENDDO      ! lc1

    DEALLOCATE(rga,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL clebsch_gordon_deallocate()

  END SUBROUTINE get_clebsch_gordon_coefficients

END MODULE lri_environment_init
