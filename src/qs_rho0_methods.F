!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho0_types
!!
!!   NAME
!!     qs_multipole_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_rho0_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,pi,rootpi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_rho0_types,                   ONLY: allocate_multipoles,&
                                             allocate_rho0_atom,&
                                             allocate_rho0_atom_rad,&
                                             allocate_rho0_mpole,&
                                             initialize_mpole_rho,&
                                             mpole_gau_overlap,&
                                             mpole_rho_atom,&
                                             rho0_mpole_type, &
                                             rho0_atom_type,&
                                             write_rho0_info
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: exp_radius,&
                                             gaussint_sph,&
                                             trace_r_AxB
  USE termination,                     ONLY: stop_memory


  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_rho0_methods"

! *** Public subroutines ***

  PUBLIC :: calculate_rho0_atom,init_rho0

!******************************************************************************

 CONTAINS

! *****************************************************************************

  SUBROUTINE  calculate_mpole_gau(mp_gau,orb_basis,harmonics,nchannels,nsotot)

    TYPE(mpole_gau_overlap)           :: mp_gau
    TYPE(gto_basis_set_type), POINTER :: orb_basis
    TYPE(harmonics_atom_type),POINTER :: harmonics
    INTEGER, INTENT(IN)               :: nchannels, nsotot

    CHARACTER(LEN=*), PARAMETER :: routine_name = "initialize_mpole_rho"

    INTEGER   :: icg, ig1, ig2, ipgf1, ipgf2, iset1, iset2, iso, iso1, &
                 iso2, l, l1, l2, m1, m2, maxso, n1, n2, nset
    INTEGER, DIMENSION(:), POINTER    :: lmax,lmin,npgf
    REAL(dp)                          :: zet1, zet2
    REAL(dp), DIMENSION(:,:), POINTER :: zet
    REAL(dp), DIMENSION(:,:,:), &
                           POINTER    :: my_CG

!   ---------------------------------------------------------------------------

    NULLIFY(lmax,lmin,npgf,my_CG,zet)

    CALL reallocate(mp_gau%Qlm_gg,1,nsotot,1,nsotot,1,nchannels)
    
    CALL get_gto_basis_set(gto_basis_set = orb_basis, &
                               lmax=lmax,lmin=lmin,maxso=maxso,&
                               npgf=npgf,nset=nset,zet=zet)
    my_CG => harmonics%my_CG

    m1 = 0
    m2 = 0
    DO iset1 = 1,nset
      n1 = nsoset(lmax(iset1))
      DO ipgf1 = 1,npgf(iset1)
        zet1 = zet(ipgf1,iset1)
        DO iset2 = 1,nset
          n2 = nsoset(lmax(iset2))
          DO ipgf2 = 1,npgf(iset2)
            zet2 = zet(ipgf2,iset2)
          
            DO iso = 1,nchannels
              l = indso(1,iso)
              DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                 iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                 iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                 l1 = indso(1,iso1)
                 l2 = indso(1,iso2)
                 ig1 = iso1 + n1*(ipgf1-1) + m1
                 ig2 = iso2 + n2*(ipgf2-1) + m2
                 mp_gau%Qlm_gg(ig1,ig2,iso) = fourpi/(2._dp*l+1._dp)*&
                        my_CG(iso1,iso2,iso)*gaussint_sph(zet1+zet2,l+l1+l2)
              END DO  ! icg
            END DO  ! iso

          END DO  ! ipgf2 
          m2 = m2 + maxso
        END DO  ! iset2
      END DO ! ipgf1
      m1 = m1 + maxso
    END DO  ! iset1

  END SUBROUTINE calculate_mpole_gau

! *****************************************************************************

  SUBROUTINE calculate_rho0_atom(qs_env,rho0_mp,a_list,g_atom,&
                                 paw_atom,natom,ikind)

    TYPE(qs_environment_type), INTENT(inout) :: qs_env
    INTEGER, DIMENSION(:), INTENT(IN)        :: a_list
    TYPE(grid_atom_type), INTENT(IN)         :: g_atom
    TYPE(rho0_mpole_type), INTENT(INOUT)     :: rho0_mp
    INTEGER, INTENT(IN)                      :: ikind, natom
    LOGICAL, INTENT(IN)                      :: paw_atom

    CHARACTER(len=*), PARAMETER :: routine_name = 'calculate_rho0_atom', &
      routineP = module_name//':'//routine_name

    REAL(dp)                                 :: z_h, z_s
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: gh_tmp,gs_tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)    :: g0_h, g0_s
    TYPE(mpole_rho_atom),  POINTER           :: mpole_rho
    TYPE(mpole_gau_overlap), POINTER         :: mpole_gau
    TYPE(rho0_atom_type), DIMENSION(:),&
                           POINTER           :: rho0_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(rho_atom_type), DIMENSION(:), &
                           POINTER           :: rho_atom_set

    INTEGER :: iat, iatom, iso, istat, l, lmax0, nchannels, nr
!   ---------------------------------------------------------------------------

    NULLIFY(mpole_gau)
    NULLIFY(mpole_rho)
    NULLIFY(rho0_atom_set)
    NULLIFY(rho_atom_set)

    CALL get_qs_env(qs_env=qs_env,rho_atom_set=rho_atom_set,&
                    rho0_atom_set=rho0_atom_set)

    mpole_gau => rho0_mp%mp_gau(ikind)
    lmax0 = rho0_mp%lmax0_kind(ikind)
    z_h = rho0_mp%zet0_h
    z_s = rho0_mp%zet0_s
    nchannels = nsoset(lmax0) 
    nr = g_atom%nr

! Set density coefficient to zero befor the calculation
    DO iat = 1,natom
      iatom = a_list(iat)
      rho0_atom_set(iatom)%rho0_rad_h%r_coef = 0.0_dp
      rho0_atom_set(iatom)%rho0_rad_s%r_coef = 0.0_dp
      rho0_mp%mp_rho(iatom)%Qlm_tot = 0.0_dp
      rho0_mp%mp_rho(iatom)%Qlm_tot(1) = rho0_mp%mp_rho(iatom)%Qlm_z
    ENDDO

    ALLOCATE(g0_h(nr,0:lmax0),g0_s(nr,0:lmax0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "g0_h,g0_s")
    ALLOCATE(gh_tmp(nr),gs_tmp(nr),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gh_tmp,gs_tmp")
    gh_tmp = 0.0_dp
    gs_tmp = 0.0_dp



    gh_tmp(1:nr) = EXP(-z_h*g_atom%rad2(1:nr))
    g0_h(1:nr,0) = gh_tmp(1:nr)*rho0_mp%norm_g0l_h(0)
    gs_tmp(1:nr) = EXP(-z_s*g_atom%rad2(1:nr))
    g0_s(1:nr,0) = gs_tmp(1:nr)*rho0_mp%norm_g0l_s(0)

    DO l = 1,lmax0
      gh_tmp(1:nr) = gh_tmp(1:nr)*g_atom%rad(1:nr)
      g0_h(1:nr,l) = gh_tmp(1:nr)*rho0_mp%norm_g0l_h(l)
      gs_tmp(1:nr) = gs_tmp(1:nr)*g_atom%rad(1:nr)
      g0_s(1:nr,l) = gs_tmp(1:nr)*rho0_mp%norm_g0l_s(l)
    END DO  ! l

    DO iso = 1, nsoset(lmax0)
      l = indso(1,iso)
      DO iat = 1,natom
         iatom = a_list(iat)
         mpole_rho => rho0_mp%mp_rho(iatom)
         rho_atom => rho_atom_set(iatom)
         IF(paw_atom) CALL calculate_Qlm_atom(mpole_gau,mpole_rho,rho_atom,iso)
         rho0_atom_set(iatom)%rho0_rad_h%r_coef(1:nr,iso) = &
              g0_h(1:nr,l) * mpole_rho%Qlm_tot(iso)  
      END DO  ! iat

    END DO  ! iso
    
    DEALLOCATE(g0_h,g0_s,gh_tmp,gs_tmp,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,&
                         "g0_h,g0_s,gh_tmp,gs_tmp")
    ENDIF

  END SUBROUTINE calculate_rho0_atom


! *****************************************************************************

  SUBROUTINE calculate_Qlm_atom(mpole_gau,mpole_rho,rho_atom,iso)

    TYPE(mpole_rho_atom),  POINTER           :: mpole_rho
    TYPE(mpole_gau_overlap), POINTER         :: mpole_gau
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      ::  iso

    TYPE(rho_atom_coeff), DIMENSION(:),&
                                    POINTER  :: cpc_h, cpc_s
    INTEGER                                  :: ispin, nsotot, nspins

!   ---------------------------------------------------------------------------

    CALL get_rho_atom(rho_atom=rho_atom,cpc_h=cpc_h,cpc_s=cpc_s)

    nspins = SIZE(cpc_h)
    nsotot = SIZE(mpole_gau%Qlm_gg,1)

    mpole_rho%Qlm_h(iso) = 0.0_dp
    mpole_rho%Qlm_s(iso) = 0.0_dp

    DO ispin = 1,nspins
      mpole_rho%Qlm_h(iso) =  mpole_rho%Qlm_h(iso) + &
            trace_r_AxB(mpole_gau%Qlm_gg(:,:,iso),nsotot,&
                        cpc_h(ispin)%r_coef,nsotot,&
                        nsotot,nsotot)
      mpole_rho%Qlm_s(iso) =  mpole_rho%Qlm_s(iso) + &
            trace_r_AxB(mpole_gau%Qlm_gg(:,:,iso),nsotot,&
                        cpc_s(ispin)%r_coef,nsotot,&
                        nsotot,nsotot)
    END DO  ! ispin
    mpole_rho%Qlm_tot(iso) = mpole_rho%Qlm_tot(iso) + &
                             mpole_rho%Qlm_h(iso)-mpole_rho%Qlm_s(iso)

  END SUBROUTINE calculate_Qlm_atom

! *****************************************************************************

  SUBROUTINE init_rho0(qs_env,gapw_control,globenv)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(global_environment_type)            :: globenv    

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "init_rho0"

    TYPE(atomic_kind_type), DIMENSION(:), &
                            POINTER          :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rho0_atom_type), DIMENSION(:), &
                           POINTER           :: rho0_atom_set
    REAL(dp)                                 :: ecut, eps_fit, gcut, &
                                                radius, rc_orb, rc_min
    INTEGER  ::  iat, iatom, ikind, l, laddg, lmaxg, maxl, maxso, nat, &
                 natom, nchannels, nkind, nr, nset, nsotot,  zeff
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set)
    NULLIFY(atom_kind)
    NULLIFY(dft_control)
    NULLIFY(harmonics)
    NULLIFY(orb_basis)
    NULLIFY(rho0_mpole)
    NULLIFY(rho0_atom_set)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                    dft_control=dft_control)

    ecut = dft_control%qs_control%cutoff 
    nkind = size(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom)

!   Initialize the multipole and the compensation charge type 
    CALL allocate_rho0_mpole(rho0_mpole)
    CALL allocate_rho0_atom(rho0_atom_set,natom)

!   Allocate the multipole set
    CALL allocate_multipoles(rho0_mpole%mp_rho,natom,rho0_mpole%mp_gau,nkind)

!   For each kind, determine the max l for the compensation charge density
    lmaxg = gapw_control%lmax_rho0
    laddg = gapw_control%ladd_rho0
    eps_fit = gapw_control%eps_fit

    CALL reallocate(rho0_mpole%lmax0_kind,1,nkind)

    rho0_mpole%lmax_0 = 0
    rc_min = 100.0_dp
    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                          orb_basis_set=orb_basis, &
                          atom_list=atom_list,natom=nat, &
                          harmonics=harmonics,ngrid_rad=nr,&
                          paw_atom=paw_atom,zeff=zeff)

      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             kind_radius=rc_orb,maxl=maxl,&
                             maxso=maxso,nset=nset)

      IF(paw_atom) THEN
        rho0_mpole%lmax0_kind(ikind) = min(2*maxl, maxl+laddg,lmaxg)
      ELSE
        rho0_mpole%lmax0_kind(ikind) = 0
      END IF

      rho0_mpole%lmax_0 = max(rho0_mpole%lmax_0,rho0_mpole%lmax0_kind(ikind))
      rc_min = min(rc_min,rc_orb)

      nchannels = nsoset(rho0_mpole%lmax0_kind(ikind))
      nsotot = maxso * nset

      DO iat = 1,nat
         iatom = atom_list(iat)

!        Allocate the multipole for rho1_h rho1_s and rho_z
         CALL initialize_mpole_rho(rho0_mpole%mp_rho(iatom),nchannels,zeff)

!        Allocate the radial part of rho0_h and rho0_s 
!        This is calculated on the radial grid centered at the atomic position
         CALL allocate_rho0_atom_rad(rho0_atom_set(iatom),nr,nchannels)
      END DO

      IF(paw_atom) THEN
         CALL calculate_mpole_gau(rho0_mpole%mp_gau(ikind),&
                                  orb_basis,harmonics,nchannels,nsotot)
      END IF
        
    END DO   ! ikind

!   Calculate the hard and soft esponents for the spherical gaussians
!   used for the construction of the compensation charge densities (g0_lm_h,g0_lm_s)
    rho0_mpole%zet0_h =  0.1_dp
    DO
      radius=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_h,1.0_dp,eps_fit)
      IF (radius<=rc_min) EXIT
      rho0_mpole%zet0_h=rho0_mpole%zet0_h+0.1_dp 
    END DO

    rho0_mpole%zet0_s = 0.25_dp/rho0_mpole%zet0_h

    DO
      gcut=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_s,1.0_dp,eps_fit)
      If (gcut<=SQRT(ecut)) EXIT
      rho0_mpole%zet0_s=rho0_mpole%zet0_s*1.1_dp 
    END DO
    rho0_mpole%zet0_s=0.25_dp/rho0_mpole%zet0_s

! check whether this is necessary to build the proper grid
! expmax=max(expmax,rho0_mpole%zet0_h)
!  expmin=min(expmin,rho0_mpole%zet0_s)

!   Allocate and calculate the normalization factors for g0_lm_h and g0_lm_s
    CALL reallocate(rho0_mpole%norm_g0l_h,0,rho0_mpole%lmax_0)
    CALL reallocate(rho0_mpole%norm_g0l_s,0,rho0_mpole%lmax_0)
    DO l = 0,rho0_mpole%lmax_0
      rho0_mpole%norm_g0l_h(l) = (2._dp*l+1._dp)/&
                                 ( fourpi*gaussint_sph(rho0_mpole%zet0_h,2*l)) 
      rho0_mpole%norm_g0l_s(l) = (2._dp*l+1._dp)/&
                                 ( fourpi*gaussint_sph(rho0_mpole%zet0_s,2*l))
    END DO

    CALL set_qs_env(qs_env=qs_env, rho0_mpole = rho0_mpole, &
                                   rho0_atom_set=rho0_atom_set)

    IF(globenv%print%rho0_information) THEN
      CALL write_rho0_info(rho0_mpole,globenv)
    END IF

  END SUBROUTINE init_rho0
! *****************************************************************************

! *****************************************************************************
END MODULE qs_rho0_methods
!******************************************************************************
