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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,fourpi,pi,rootpi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco, indso,&
                                             nco, ncoset,&
                                             nso, nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE qs_charges_types,                ONLY: qs_charges_type
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
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rho0_atom_type), DIMENSION(:), &
                           POINTER           :: rho0_atom_set
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(rhoz_type), DIMENSION(:), POINTER   :: rhoz_set
    REAL(dp)                                 :: alpha_core,&
                                                ecut, eps_fit, eps_Vrho0,&
                                                gcut, radius, rhoz_tot,&
                                                rc_orb, rc_min
    INTEGER  ::  iat, iatom, ikind, l, laddg, lmaxg, max_nr, maxl, maxso, nat, &
                 natom, nchan_c, nchan_s, nkind, nr, nset, nsotot,  zeff
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
    NULLIFY(qs_charges)
    NULLIFY(rhoz_set)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                    dft_control=dft_control,qs_charges=qs_charges)

    ecut = dft_control%qs_control%cutoff 
    nkind = SIZE(atomic_kind_set,1)
    eps_Vrho0 = dft_control%qs_control%gapw_control%eps_Vrho0

!   Initialize rhoz total to zero
!   in gapw rhoz is calculated on local the lebedev grids
    qs_charges%total_rho_core_rspace = 0.0_dp

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,max_ngrid_rad=max_nr)

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

      IF(paw_atom) THEN
        rho0_mpole%lmax0_kind(ikind) = min(2*maxl, maxl+laddg,lmaxg)
      ELSE
        rho0_mpole%lmax0_kind(ikind) = 0
      END IF

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
!        Calculate multipoles given by the prodct of 2 primitives Qlm_gg
         CALL calculate_mpole_gau(rho0_mpole%mp_gau(ikind),&
                                  orb_basis,harmonics,nchan_s,nsotot)
      END IF

!     Calculate the core density rhoz 
!                  exp(-alpha_c**2 r**2)Z(alpha_c**2/pi)**(3/2)
!     on the logarithmic radial grid
!     WARNING: alpha_core_charge = alpha_c**2 
      CALL calculate_rhoz(rhoz_set(ikind),grid_atom,alpha_core,zeff,&
            nat,qs_charges%total_rho_core_rspace,harmonics)

    END DO   ! ikind
    qs_charges%total_rho_core_rspace = - qs_charges%total_rho_core_rspace

!   Calculate the hard and soft exponents for the spherical gaussians
!   used for the construction of the compensation charge densities (g0_lm_h,g0_lm_s)

    rho0_mpole%zet0_h =  0.1_dp
    DO
      radius=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_h,eps_fit,1.0_dp)
      IF (radius<=rc_min) EXIT
      rho0_mpole%zet0_h=rho0_mpole%zet0_h+0.1_dp 
    END DO

!dbg
    rho0_mpole%zet0_h = 10._dp
!dbg

    rho0_mpole%zet0_s = 0.25_dp/rho0_mpole%zet0_h
    DO
      gcut=exp_radius(rho0_mpole%lmax_0,rho0_mpole%zet0_s,eps_fit,1.0_dp)
      IF (gcut<=SQRT(2.0_dp*ecut)) EXIT
      rho0_mpole%zet0_s=rho0_mpole%zet0_s*1.1_dp 
    END DO
    rho0_mpole%zet0_s=0.25_dp/rho0_mpole%zet0_s
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
    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           grid_atom=grid_atom)

      CALL calculate_g0(rho0_mpole,grid_atom,ikind)

      CALL interaction_radii_g0(rho0_mpole,grid_atom,ikind,eps_Vrho0,.TRUE.)

    END DO

    rhoz_tot = qs_charges%total_rho_core_rspace

    IF (.NOT.my_tddft) THEN
       CALL set_qs_env(qs_env=qs_env, rho0_mpole = rho0_mpole, &
                                      rho0_atom_set=rho0_atom_set,&
                                      rhoz_set = rhoz_set,&
                                      rhoz_tot=rhoz_tot)
!      Allocate and Initialie rho0 soft on the global grid
       CALL rho0_s_grid_create(qs_env, qs_env%local_rho_set, error)
    ELSE
       tddft_local_rho_set%rho0_mpole => rho0_mpole
       tddft_local_rho_set%rho0_atom_set => rho0_atom_set
       tddft_local_rho_set%rhoz_set => rhoz_set
       tddft_local_rho_set%rhoz_tot = rhoz_tot
       write (*,*) rhoz_tot
       CALL rho0_s_grid_create(qs_env, tddft_local_rho_set, error)
    END IF

    CALL write_rho0_info(rho0_mpole,cell,error)

    END SUBROUTINE init_rho0

! *****************************************************************************

  SUBROUTINE interaction_radii_g0(rho0_mpole,grid_atom,ik,&
                                  eps_Vrho0,calculate_V0hmV0s)

    TYPE(rho0_mpole_type), POINTER      :: rho0_mpole
    TYPE(grid_atom_type), POINTER       :: grid_atom
    INTEGER , INTENT(IN)                :: ik
    REAL(dp), INTENT(IN)                :: eps_Vrho0
    LOGICAL , INTENT(IN)                :: calculate_V0hmV0s

    CHARACTER(LEN=*), PARAMETER :: routine_name = "interaction_radii_g0"

    INTEGER                             :: ir, l, lmax, nr
    REAL(dp)                            :: z_h, z_s
    REAL(dp)                            :: r_h, r_s, r_hms
    REAL(dp)                            :: I_down_h, I_down_s
    REAL(dp)                            :: I_up_h, I_up_s
    REAL(dp)                            :: V_h, V_s
    REAL(dp), DIMENSION(:), POINTER     :: n_h, n_s
    REAL(dp), DIMENSION(:,:), POINTER   :: g0_h, g0_s

!   ---------------------------------------------------------------------------


    CALL get_rho0_mpole(rho0_mpole, ikind =ik, l0_ikind=lmax,&
                        zet0_h=z_h, zet0_s=z_s, norm_g0l_h=n_h,&
                        norm_g0l_s=n_s, g0_h=g0_h, g0_s=g0_s)
    nr = grid_atom%nr
    
    r_h = 0.0_dp
    r_s = 0.0_dp
    DO l = 0,lmax
      r_h = max(r_h, exp_radius(l,z_h,eps_Vrho0,n_h(l)))
      r_s = max(r_s, exp_radius(l,z_s,eps_Vrho0,n_s(l)))
    END DO

    rho0_mpole%mp_gau(ik)%rpgf0_h = r_h
    rho0_mpole%mp_gau(ik)%rpgf0_s = r_s

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

END MODULE qs_rho0_methods

!******************************************************************************
