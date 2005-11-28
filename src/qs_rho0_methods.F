!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho0_methods
!!
!!   NAME
!!     qs_rho0_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_rho0_methods

  USE ai_coulomb,                      ONLY: coulomb3
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             set_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE external_potential_types,        ONLY: all_potential_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,fourpi,pi,rootpi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco, indso,&
                                             nco, ncoset,&
                                             nso, nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_local_rho_types,              ONLY: local_rho_type
  USE qs_rho0_ggrid,                   ONLY: rho0_s_grid_create
  USE qs_rho0_types,                   ONLY: allocate_multipoles,&
                                             allocate_rho0_atom,&
                                             allocate_rho0_atom_rad,&
                                             allocate_rho0_mpole,&
                                             calculate_g0,&
                                             get_rho0_mpole, &
                                             initialize_mpole_rho,&
                                             mpole_gau_overlap,&
                                             mpole_rho_atom,&
                                             rho0_mpole_type, &
                                             rho0_atom_type,&
                                             write_rho0_info
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_local_rho_types,              ONLY: allocate_rhoz,&
                                             calculate_rhoz,&
                                             rhoz_type
  USE qs_util,                         ONLY: exp_radius,&
                                             gaussint_sph,&
                                             trace_r_AxB
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

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
    DO iset1 = 1,nset
      n1 = nsoset(lmax(iset1))
      DO ipgf1 = 1,npgf(iset1)
        zet1 = zet(ipgf1,iset1)
        m2 = 0
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

  SUBROUTINE calculate_rho0_atom(rho_atom_set, rho0_atom_set, &
                                 rho0_mp,a_list,g_atom,&
                                 paw_atom,natom,ikind,harmonics,&
                                 rho0_h_tot,rho0_s_tot)

    TYPE(rho_atom_type), DIMENSION(:), &
                           POINTER           :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:),&
                           POINTER           :: rho0_atom_set
    INTEGER, DIMENSION(:), INTENT(IN)        :: a_list
    TYPE(grid_atom_type), INTENT(IN)         :: g_atom
    TYPE(rho0_mpole_type), POINTER           :: rho0_mp
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: ikind, natom
    LOGICAL, INTENT(IN)                      :: paw_atom
    REAL(dp), INTENT(INOUT)                  :: rho0_h_tot,rho0_s_tot

    CHARACTER(len=*), PARAMETER :: routine_name = 'calculate_rho0_atom', &
      routineP = module_name//':'//routine_name

    REAL(dp), DIMENSION(:,:), POINTER        :: g0_h, g0_s,vg0_h
    TYPE(mpole_rho_atom),  POINTER           :: mpole_rho
    TYPE(mpole_gau_overlap), POINTER         :: mpole_gau
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(:), POINTER          :: norm_g0l_h, norm_g0l_s

    INTEGER                                  :: handle, iat, iatom, ic, ico, ir, is, iso,&
                                                istat, l, lshell, lx, ly, lz,&
                                                lmax0, nr,i,j
    REAL(dp) :: sum1, sum2
!   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho0_atom","I"," ",handle)

    NULLIFY(mpole_gau)
    NULLIFY(mpole_rho)
    NULLIFY(g0_h, g0_s, vg0_h)
    NULLIFY(norm_g0l_h, norm_g0l_s)

    CALL get_rho0_mpole(rho0_mpole=rho0_mp,ikind=ikind, &
                        l0_ikind=lmax0,mp_gau_ikind=mpole_gau,&
                        g0_h=g0_h,g0_s=g0_s, vg0_h=vg0_h,&
                        norm_g0l_h=norm_g0l_h,&
                        norm_g0l_s=norm_g0l_s)

    nr = g_atom%nr

! Set density coefficient to zero befor the calculation
    DO iat = 1,natom
      iatom = a_list(iat)
      rho0_atom_set(iatom)%rho0_rad_h%r_coef = 0.0_dp
      rho0_atom_set(iatom)%rho0_rad_s%r_coef = 0.0_dp
      rho0_mp%mp_rho(iatom)%Qlm_tot = 0.0_dp
      rho0_mp%mp_rho(iatom)%Qlm_tot(1) = rho0_mp%mp_rho(iatom)%Qlm_z
    ENDDO


    DO iso = 1, nsoset(lmax0)
      l = indso(1,iso)
      DO iat = 1,natom
         iatom = a_list(iat)
         mpole_rho => rho0_mp%mp_rho(iatom)
         rho_atom => rho_atom_set(iatom)
         IF(paw_atom) CALL calculate_Qlm_atom(mpole_gau,mpole_rho,rho_atom,iso)
         rho0_atom_set(iatom)%rho0_rad_h%r_coef(1:nr,iso) = &
              g0_h(1:nr,l) * mpole_rho%Qlm_tot(iso)  
         rho0_atom_set(iatom)%rho0_rad_s%r_coef(1:nr,iso) = &
              g0_s(1:nr,l) * mpole_rho%Qlm_tot(iso)
         rho0_atom_set(iatom)%vrho0_rad_h%r_coef(1:nr,iso) = &
              vg0_h(1:nr,l) * mpole_rho%Qlm_tot(iso)  

         sum1 = 0.0_dp
         sum2 = 0.0_dp
         DO ir = 1,nr
           sum1 = sum1 + g_atom%wr(ir)*&
                  rho0_atom_set(iatom)%rho0_rad_h%r_coef(ir,iso)
           sum2 = sum2 + g_atom%wr(ir)*&
                  rho0_atom_set(iatom)%rho0_rad_s%r_coef(ir,iso)
         ENDDO
         rho0_h_tot = rho0_h_tot + sum1* harmonics%slm_int(iso)
         rho0_s_tot = rho0_s_tot + sum2* harmonics%slm_int(iso)
      END DO  ! iat
    END DO  ! iso

!   Transform the coefficinets from spherical to cartesian

    DO iat = 1,natom
       iatom = a_list(iat)
       mpole_rho => rho0_mp%mp_rho(iatom)

       DO lshell = 0,lmax0
         DO ic = 1,nco(lshell)
           ico = ic + ncoset(lshell-1)
           mpole_rho%Qlm_car_h(ico) = 0.0_dp
           mpole_rho%Qlm_car_s(ico) = 0.0_dp
           lx = indco(1,ico)
           ly = indco(2,ico)
           lz = indco(3,ico)

           DO is = 1,nso(lshell)
             iso = is + nsoset(lshell-1)

             mpole_rho%Qlm_car_s(ico) = mpole_rho%Qlm_car_s(ico) + &
                  orbtramat(lshell)%c2s(is,ic)*mpole_rho%Qlm_tot(iso)*&
                  norm_g0l_s(lshell)&
                 /SQRT(dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)*fourpi/dfac(2*lshell+1))

             mpole_rho%Qlm_car_h(ico) = mpole_rho%Qlm_car_h(ico) + &
                  orbtramat(lshell)%c2s(is,ic)*mpole_rho%Qlm_tot(iso)*&
                  norm_g0l_h(lshell)&
                 /SQRT(dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)*fourpi/dfac(2*lshell+1))

           END DO 
         END DO
       END DO  ! lshell
    END DO  ! iat

   CALL timestop(0.0_dp,handle)

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

    NULLIFY(cpc_h,cpc_s)
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

  SUBROUTINE init_rho0(qs_env, gapw_control, cell, &
                       tddft, tddft_local_rho_set, error)

    ! ARGUMENTS
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(local_rho_type), POINTER, OPTIONAL  :: tddft_local_rho_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error   

    ! PARAMETERS
    CHARACTER(LEN=*), PARAMETER              :: routine_name = "init_rho0"

    ! LOCALS
    TYPE(atomic_kind_type), DIMENSION(:), &
                            POINTER          :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rho0_atom_type), DIMENSION(:), &
                           POINTER           :: rho0_atom_set
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(rhoz_type), DIMENSION(:), POINTER   :: rhoz_set
    REAL(dp)                                 :: alpha_core,&
                                                ecut, eps_fit, eps_Vrho0,&
                                                gcut, radius,&
                                                rc_orb, rc_min, &
                                                max_rpgf0_s, &
                                                total_rho_core_rspace, &
                                                zetsoft_max
    INTEGER  ::  iat, iatom, ikind, istat, l, l_rho1_max, laddg, lmaxg, max_nr, maxl, maxnset, maxso, nat, &
                 natom, nchan_c, nchan_s, nkind, nr, nset, nsotot
    REAL(KIND=dp)                            :: zeff
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom, my_tddft

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set)
    NULLIFY(atom_kind)
    NULLIFY(dft_control)
    NULLIFY(harmonics)
    NULLIFY(orb_basis)
    NULLIFY(rho0_mpole)
    NULLIFY(rho0_atom_set)
    NULLIFY(rhoz_set)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                    dft_control=dft_control)

    ecut = dft_control%qs_control%cutoff 
    nkind = SIZE(atomic_kind_set,1)
    eps_Vrho0 = gapw_control%eps_Vrho0

!   Initialize rhoz total to zero
!   in gapw rhoz is calculated on local the lebedev grids
    total_rho_core_rspace = 0.0_dp

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,max_ngrid_rad=max_nr,&
                             zetsoft_max=zetsoft_max)

!   Initialize the multipole and the compensation charge type 
    CALL allocate_rho0_mpole(rho0_mpole)
    CALL allocate_rho0_atom(rho0_atom_set,natom)

!   Allocate the multipole set
    CALL allocate_multipoles(rho0_mpole%mp_rho,natom,rho0_mpole%mp_gau,nkind)

!   Allocate the core density on the radial grid for each kind: rhoz_set
    CALL allocate_rhoz(rhoz_set,nkind)

!   For each kind, determine the max l for the compensation charge density
    lmaxg = gapw_control%lmax_rho0
    laddg = gapw_control%ladd_rho0
    eps_fit = gapw_control%eps_fit

    CALL reallocate(rho0_mpole%lmax0_kind,1,nkind)

    rho0_mpole%lmax_0 = 0
    rc_min = 100.0_dp
    maxnset = 0
    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                          orb_basis_set=orb_basis, &
                          atom_list=atom_list,natom=nat, &
                          harmonics=harmonics,&
                          grid_atom=grid_atom,&
                          ngrid_rad=nr,&
                          paw_atom=paw_atom,zeff=zeff,&
                          alpha_core_charge=alpha_core,&
                          hard0_radius=rc_orb)

      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             maxl=maxl,&
                             maxso=maxso,nset=nset)

      maxnset = MAX(maxnset,nset)

      l_rho1_max = indso(1,harmonics%max_iso_not0)
      IF(paw_atom) THEN
        rho0_mpole%lmax0_kind(ikind) = MIN(2*maxl,l_rho1_max,maxl+laddg,lmaxg)
      ELSE
        rho0_mpole%lmax0_kind(ikind) = 0
      END IF

      CALL set_atomic_kind(atomic_kind=atom_kind,&
                           lmax_rho0 = rho0_mpole%lmax0_kind(ikind))

      IF(gapw_control%lrho1_eq_lrho0) harmonics%max_iso_not0 = &
           nsoset(rho0_mpole%lmax0_kind(ikind))

      rho0_mpole%lmax_0 = max(rho0_mpole%lmax_0,rho0_mpole%lmax0_kind(ikind))
      rc_min = min(rc_min,rc_orb)

      nchan_s = nsoset(rho0_mpole%lmax0_kind(ikind))
      nchan_c = ncoset(rho0_mpole%lmax0_kind(ikind))
      nsotot = maxso * nset

      DO iat = 1,nat
         iatom = atom_list(iat)

!        Allocate the multipole for rho1_h rho1_s and rho_z
         CALL initialize_mpole_rho(rho0_mpole%mp_rho(iatom),nchan_s,nchan_c,zeff,my_tddft)

!        Allocate the radial part of rho0_h and rho0_s 
!        This is calculated on the radial grid centered at the atomic position
         CALL allocate_rho0_atom_rad(rho0_atom_set(iatom),nr,nchan_s)
      END DO
!     
      IF(paw_atom) THEN
!        Calculate multipoles given by the product of 2 primitives Qlm_gg
         CALL calculate_mpole_gau(rho0_mpole%mp_gau(ikind),&
                                  orb_basis,harmonics,nchan_s,nsotot)
      END IF

!     Calculate the core density rhoz 
!                  exp(-alpha_c**2 r**2)Z(alpha_c**2/pi)**(3/2)
!     on the logarithmic radial grid
!     WARNING: alpha_core_charge = alpha_c**2 
      CALL calculate_rhoz(rhoz_set(ikind),grid_atom,alpha_core,zeff,&
            nat,total_rho_core_rspace,harmonics)

    END DO   ! ikind
    total_rho_core_rspace = - total_rho_core_rspace

    IF(gapw_control%alpha0_hard_from_input) THEN
!   The Exponent for the compensation charge rho0_hard is read from input

      rho0_mpole%zet0_h = gapw_control%alpha0_hard

    ELSE
!   Calculate the exponent for the compensation charge rho0_hard 

      rho0_mpole%zet0_h =  0.1_dp
      DO
!        radius=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_h,eps_fit,1.0_dp)
        radius=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_h,eps_Vrho0,1.0_dp)
        IF (radius<=rc_min) EXIT
        rho0_mpole%zet0_h=rho0_mpole%zet0_h+0.1_dp 
      END DO

    END IF

    IF(gapw_control%alpha0_soft_from_input) THEN
!   The Exponent for the compensation charge rho0_soft is read from input

      rho0_mpole%zet0_s = gapw_control%alpha0_soft

    ELSE
!   Calculate the exponent for the compensation charge rho0_soft

      rho0_mpole%zet0_s = 0.25_dp/rho0_mpole%zet0_h
      DO
!        gcut=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_s,eps_fit,1.0_dp)
        gcut=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_s,eps_Vrho0,1.0_dp)
        IF (gcut<=SQRT(2.0_dp*ecut)) EXIT
        rho0_mpole%zet0_s=rho0_mpole%zet0_s*1.1_dp 
      END DO
      rho0_mpole%zet0_s=0.25_dp/rho0_mpole%zet0_s

!     rho0_mpole%zet0_s can be as big as the largest exponent of the soft density
!     that is collocated on the grid
      IF(rho0_mpole%zet0_s < zetsoft_max) rho0_mpole%zet0_s = zetsoft_max

      IF(rho0_mpole%zet0_s > 0.9_dp*rho0_mpole%zet0_h) &
           rho0_mpole%zet0_s=0.9_dp*rho0_mpole%zet0_h
    END IF

!   Allocate and calculate the normalization factors for g0_lm_h and g0_lm_s
    CALL reallocate(rho0_mpole%norm_g0l_h,0,rho0_mpole%lmax_0)
    CALL reallocate(rho0_mpole%norm_g0l_s,0,rho0_mpole%lmax_0)
    DO l = 0,rho0_mpole%lmax_0
      rho0_mpole%norm_g0l_h(l) = (2._dp*l+1._dp)/&
                                 ( fourpi*gaussint_sph(rho0_mpole%zet0_h,2*l)) 
      rho0_mpole%norm_g0l_s(l) = (2._dp*l+1._dp)/&
                                 ( fourpi*gaussint_sph(rho0_mpole%zet0_s,2*l))
    END DO


!   Allocate and Initialize the g0 gaussians used to build the compensation density
!   an calculate the interaction radii 
    max_rpgf0_s = 0.0_dp
    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           grid_atom=grid_atom)

      CALL calculate_g0(rho0_mpole,grid_atom,ikind)

      CALL interaction_radii_g0(rho0_mpole,grid_atom,ikind,eps_Vrho0,max_rpgf0_s,.TRUE.)

    END DO
    rho0_mpole%max_rpgf0_s = max_rpgf0_s
    ALLOCATE(rho0_mpole%orb_radius_3c(nkind),STAT=istat)
    ALLOCATE(rho0_mpole%set_radius_3c(maxnset,nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "orb_radius_3c,set_radius_3c")


!   Calculate the interaction radii for the 3 centers lists
    CALL radii_3c_orbxrho0(qs_env,rho0_mpole)

    IF (.NOT.my_tddft) THEN
       CALL set_qs_env(qs_env=qs_env, rho0_mpole = rho0_mpole, &
                                      rho0_atom_set=rho0_atom_set,&
                                      rhoz_set = rhoz_set,&
                                      rhoz_tot=total_rho_core_rspace)
    ELSE
       tddft_local_rho_set%rho0_mpole => rho0_mpole
       tddft_local_rho_set%rho0_atom_set => rho0_atom_set
       tddft_local_rho_set%rhoz_set => rhoz_set
       tddft_local_rho_set%rhoz_tot = total_rho_core_rspace
       CALL rho0_s_grid_create(qs_env, rho0_mpole, .TRUE., error)
    END IF


    CALL write_rho0_info(rho0_mpole,cell,error)

    END SUBROUTINE init_rho0

! *****************************************************************************

  SUBROUTINE interaction_radii_g0(rho0_mpole,grid_atom,ik,&
                                  eps_Vrho0,max_rpgf0_s,calculate_V0hmV0s)

    TYPE(rho0_mpole_type), POINTER      :: rho0_mpole
    TYPE(grid_atom_type), POINTER       :: grid_atom
    INTEGER , INTENT(IN)                :: ik
    REAL(dp), INTENT(IN)                :: eps_Vrho0
    REAL(dp), INTENT(INOUT)             :: max_rpgf0_s
    LOGICAL , INTENT(IN)                :: calculate_V0hmV0s

    CHARACTER(LEN=*), PARAMETER :: routine_name = "interaction_radii_g0"

    INTEGER                             :: ir, l, lmax, nr
    REAL(dp)                            :: z0_h, z0_s
    REAL(dp)                            :: r_h, r_s, r_hms
    REAL(dp)                            :: I_down_h, I_down_s
    REAL(dp)                            :: I_up_h, I_up_s
    REAL(dp)                            :: V_h, V_s
    REAL(dp), DIMENSION(:), POINTER     :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:,:), POINTER   :: g0_h, g0_s

!   ---------------------------------------------------------------------------


    CALL get_rho0_mpole(rho0_mpole, ikind =ik, l0_ikind=lmax,&
                        zet0_h=z0_h, zet0_s=z0_s, norm_g0l_h=ng0_h,&
                        norm_g0l_s=ng0_s, g0_h=g0_h, g0_s=g0_s)
    nr = grid_atom%nr
    
    r_h = 0.0_dp
    r_s = 0.0_dp
    DO l = 0,lmax
      r_h = max(r_h, exp_radius(l,z0_h,eps_Vrho0,ng0_h(l)))
      r_s = max(r_s, exp_radius(l,z0_s,eps_Vrho0,ng0_s(l)))
    END DO

    rho0_mpole%mp_gau(ik)%rpgf0_h = r_h
    rho0_mpole%mp_gau(ik)%rpgf0_s = r_s
    max_rpgf0_s = MAX(max_rpgf0_s,r_s)

    IF(calculate_V0hmV0s) THEN

      r_hms = 0.0_dp
      DO l = 0,lmax
        
        I_up_h = g0_h(nr,l)*grid_atom%rad2l(nr,l)*grid_atom%wr(nr)  
        I_up_s = g0_s(nr,l)*grid_atom%rad2l(nr,l)*grid_atom%wr(nr)  
        I_down_h = 0.0_dp
        I_down_s = 0.0_dp
        DO ir = nr-1,1,-1
          I_down_h = I_down_h + g0_h(ir,l)*&
                       grid_atom%oorad2l(ir,l+1)*grid_atom%wr(ir)
          I_down_s = I_down_s + g0_s(ir,l)*&
                       grid_atom%oorad2l(ir,l+1)*grid_atom%wr(ir)
        END DO
 
        V_h = grid_atom%oorad2l(nr,l+1)*I_up_h + &
              grid_atom%rad2l(nr,l)*I_down_h
        V_s = grid_atom%oorad2l(nr,l+1)*I_up_s + &
              grid_atom%rad2l(nr,l)*I_down_s

        DO ir = nr-1,1,-1

          V_h = grid_atom%oorad2l(ir,l+1)*I_up_h + &
              grid_atom%rad2l(ir,l)*I_down_h

          V_s = grid_atom%oorad2l(ir,l+1)*I_up_s + &
              grid_atom%rad2l(ir,l)*I_down_s

          IF(abs(V_h-V_s)<eps_Vrho0 .AND.grid_atom%rad(ir)>0.7_dp ) THEN
            r_hms = MAX(r_hms,grid_atom%rad(ir))
            exit
          END IF

          I_up_h = I_up_h + g0_h(ir,l)*&
                    grid_atom%rad2l(ir,l)*grid_atom%wr(ir)
          I_down_h = I_down_h - g0_h(ir,l)* &
                      grid_atom%oorad2l(ir,l+1)*grid_atom%wr(ir)
          I_up_s = I_up_s + g0_s(ir,l)*&
                    grid_atom%rad2l(ir,l)*grid_atom%wr(ir)
          I_down_s = I_down_s - g0_s(ir,l)* &
                      grid_atom%oorad2l(ir,l+1)*grid_atom%wr(ir)

        END DO  ! ir

      END DO  ! l
      rho0_mpole%mp_gau(ik)%r_V0hmV0s = r_hms

    ELSE
      rho0_mpole%mp_gau(ik)%r_V0hmV0s = r_s
    END IF
 
  END SUBROUTINE  interaction_radii_g0

!******************************************************************************

  SUBROUTINE radii_3c_orbxrho0(qs_env,rho0_mpole)

    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(rho0_mpole_type),        POINTER      :: rho0_mpole

    TYPE(atomic_kind_type), DIMENSION(:),&
                                  POINTER      :: atomic_kind_set
    TYPE(all_potential_type),     POINTER      :: all_potential
    TYPE(atomic_kind_type),       POINTER      :: atom_kind,atom_kindj
    TYPE(gto_basis_set_type),     POINTER      :: orb_basis

    CHARACTER(LEN=*), PARAMETER :: routine_name = "radii_3c_orbxrho0"
    REAL(dp), PARAMETER         :: eps_radii_3c = 1.E-4_dp ! XXXX Marcella should this be input ?

    INTEGER :: i, ic, ico, ico0, icoa, ikind, ipgf, iset, is, iso0, istat, &
               j, jkind, l0, la, lmax_0, lx, ly, lz, maxlgto, n, nco0, nso0, ncoa, nkind, nseta
    INTEGER, DIMENSION(:),        POINTER      :: la_max, la_min, npgfa 

    REAL(dp) :: intmax0, r(3,0:60), r2(0:60), rab(3), rab2, rac(3), rac2, ri, &
                rmax, rp0, z0_h, z0_s
    REAL(dp), DIMENSION(:),       POINTER      :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:),       POINTER      :: orb_rad
    REAL(dp), DIMENSION(:,:),     POINTER      :: set_rad
    REAL(dp), DIMENSION(:,:),     POINTER      :: rpgfa, sphi_a, zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: Qlm
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: r_cut_set 
    REAL(dp), DIMENSION(:), ALLOCATABLE       :: ff
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: vabc
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: int_abc, int_tmp 
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE  :: vv

    NULLIFY(atomic_kind_set,orb_rad,set_rad)

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set)
    nkind = SIZE(atomic_kind_set,1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)

    NULLIFY(ng0_h,ng0_s)
    CALL get_rho0_mpole(rho0_mpole, lmax_0=lmax_0,&
                        zet0_h=z0_h, zet0_s=z0_s, norm_g0l_h=ng0_h,&
                        norm_g0l_s=ng0_s ,max_rpgf0_s=rp0)

    ALLOCATE(vv(ncoset(maxlgto),ncoset(maxlgto),ncoset(lmax_0),&
         maxlgto+maxlgto+lmax_0+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vv")
    ALLOCATE(ff(0:maxlgto+maxlgto+lmax_0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "ff")

    nco0 = ncoset(lmax_0)
    nso0 = nsoset(lmax_0)
    ALLOCATE(Qlm(nco0),STAT=istat)
    Qlm= 0.0_dp

    orb_rad => rho0_mpole%orb_radius_3c
    set_rad => rho0_mpole%set_radius_3c

    rmax = 20.0_dp
    do i=0,60
      ri=REAL(i,dp)*rmax/60._dp+0.01_dp
      r2(i)=ri**2
      r(:,i)=(/ri,0._dp,0._dp/)
    enddo
    rab(1:3) = 0.0_dp
    rab2 = 0.0_dp

    DO ikind = 1,nkind
      NULLIFY(atom_kind,orb_basis,la_max,la_min,npgfa,rpgfa,zeta)
      atom_kind => atomic_kind_set(ikind)
!      CALL get_atomic_kind(atomic_kind=atom_kind,soft_basis_set=orb_basis)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             pgf_radius=rpgfa,&
                             sphi=sphi_a,&
                             zet=zeta)
      orb_rad(ikind) = 0.0_dp
      set_rad(:,ikind) = 0.0_dp
      ALLOCATE(r_cut_set(nseta),STAT=istat)
      r_cut_set = 0.0_dp

      DO iset = 1,nseta

        ncoa = npgfa(iset)*ncoset(la_max(iset))
        ALLOCATE(int_tmp(ncoa,ncoa,nco0),int_abc(ncoa,2,nco0),STAT=istat)
        int_abc = 0.0_dp
        
        ALLOCATE(vabc(ncoa,ncoa),STAT=istat)
        vabc=0.0_dp

        DO i = 0,60 
          rac(1:3) = r(1:3,i)
          rac2 = r2(i)

          int_tmp = 0.0_dp
          CALL coulomb3(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset),&
                        la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset),&
                        lmax_0,z0_h,rp0,0,Qlm,&
                        rab,rab2,rac,rac2,rac,rac2,&
                        vabc,int_tmp,vv,ff(0:))
          DO ico0 = 1,nco0
            l0 = indco(1,ico0)+indco(2,ico0)+indco(3,ico0)
            DO icoa = 1,ncoa
              int_abc(icoa,1,ico0)=ng0_h(l0)*int_tmp(icoa,icoa,ico0)
            END DO
          END DO  ! ico0

          int_tmp = 0.0_dp
          CALL coulomb3(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset),&
                        la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset),&
                        lmax_0,z0_s,rp0,0,Qlm,&
                        rab,rab2,rac,rac2,rac,rac2,&
                        vabc,int_tmp,vv,ff(0:))
          intmax0 = 0.0_dp
          l0 = 0

          DO ipgf = 1,npgfa(iset)
             n=(ipgf-1)*ncoset(la_max(iset))
             DO ico = 1,ncoset(la_max(iset))
               la = indco(1,ico)+indco(2,ico)+indco(3,ico)
               icoa = ico + n          
               int_abc(icoa,2,1)= ng0_s(l0)*int_tmp(icoa,icoa,1)

               intmax0 = MAX(intmax0,ABS(int_abc(icoa,1,1)-int_abc(icoa,2,1)))
             END DO
          END DO

          DO ico0 = 2,nco0
            l0 = indco(1,ico0)+indco(2,ico0)+indco(3,ico0)
            DO icoa = 1,ncoa
              int_abc(icoa,2,ico0)= ng0_s(l0)*int_tmp(icoa,icoa,ico0)
            END DO
          END DO  ! ico0

          IF( intmax0 < eps_radii_3c ) THEN
            r_cut_set(iset) = r(1,i)
            EXIT
          END IF

  
        END DO  ! i

        orb_rad(ikind) = MAX(orb_rad(ikind) ,r_cut_set(iset) ) 
        set_rad(iset,ikind) = r_cut_set(iset)*r_cut_set(iset)  

        DEALLOCATE(vabc,STAT=istat)
        DEALLOCATE(int_tmp,int_abc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
              "vabc,int_tmp,int_abc")
 
      END DO  ! iset
      
      DEALLOCATE(r_cut_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
            "r_cut_set")
!  write(*,'(A,I4,5f12.5)') 'prima ',ikind , orb_rad(ikind),set_rad(1:nseta,ikind)

    DO jkind = 1,nkind
      NULLIFY(atom_kindj,all_potential)
      atom_kindj => atomic_kind_set(jkind)
      CALL get_atomic_kind(atomic_kind=atom_kindj,all_potential=all_potential)
      IF(ASSOCIATED(all_potential)) THEN
        orb_rad(ikind) = MAX(orb_rad(ikind),all_potential%orb_radius_3c(ikind))
      END IF
    END DO
!  write(*,'(A,I4,f12.5)') 'dopo ',ikind , orb_rad(ikind)

    END DO  ! ikind

    DEALLOCATE(Qlm,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "Qlm")
    DEALLOCATE(vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vv,ff")

  END SUBROUTINE radii_3c_orbxrho0


!******************************************************************************

END MODULE qs_rho0_methods

!******************************************************************************
