!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/hartree_local_methods
!!
!!   NAME
!!     hartree_local_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE  hartree_local_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE global_types,                    ONLY: global_environment_type
  USE hartree_local_types,             ONLY: allocate_ecoul_1center,&
                                             ecoul_1center_type,&
                                             get_ecoul_1c,&
                                             set_ecoul_1c
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size, int_size
  USE mathconstants,                   ONLY: fourpi,pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_oce_types,                    ONLY: vtriple_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole, &
                                             rho0_mpole_type, &
                                             rho0_atom_type
  USE qs_local_rho_types,              ONLY: rhoz_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "hartree_local_methods"

! *** Define multipole type ***

! Public Subroutine
 
  PUBLIC :: init_coulomb_local, Vh_1c_gg_integrals

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************
 
  SUBROUTINE init_coulomb_local(qs_env,gapw_control,tddft,p_env)  
 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gapw_control_type), POINTER         :: gapw_control
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL   :: p_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "init_rho0"

    TYPE(atomic_kind_type), DIMENSION(:), &
                                   POINTER   :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(ecoul_1center_type), DIMENSION(:),&
                                POINTER      :: ecoul_1c
    TYPE(vtriple_type), DIMENSION(:,:),&
                                POINTER      :: vtriple
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom, my_tddft

    INTEGER                                  :: iat, iatom, ikind, istat,&
                                                max_iso_not0, maxso, nat,&
                                                natom, nkind, nr, nset, nsotot 

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set)
    NULLIFY(atom_kind)
    NULLIFY(orb_basis)
    NULLIFY(atom_list)
    NULLIFY(ecoul_1c)
    NULLIFY(vtriple)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
    nkind = SIZE(atomic_kind_set,1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom)

!   Allocate and Initialize 1-center Potentials and Integrals
    CALL allocate_ecoul_1center(ecoul_1c,natom)

    DO ikind = 1,nkind

      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                           max_iso_not0=max_iso_not0, ngrid_rad=nr, natom=nat,&
                           orb_basis_set=orb_basis, paw_atom=paw_atom)
 
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             maxso=maxso,nset=nset)
      nsotot = maxso*nset
 
      DO iat = 1,nat
        iatom = atom_list(iat)

        CALL reallocate(ecoul_1c(iatom)%Vh1_h%r_coef,1,nr,1,max_iso_not0)
        CALL reallocate(ecoul_1c(iatom)%Vh1_s%r_coef,1,nr,1,max_iso_not0)
        IF(.NOT. paw_atom) CYCLE
        CALL reallocate(ecoul_1c(iatom)%aVh1b_h%r_coef,1,nsotot,1,nsotot)
        CALL reallocate(ecoul_1c(iatom)%aVh1b_s%r_coef,1,nsotot,1,nsotot)
        CALL reallocate(ecoul_1c(iatom)%aVh1b_0%r_coef,1,nsotot,1,nsotot)
        CALL reallocate(ecoul_1c(iatom)%aVhsb_0s%r_coef,1,nsotot,1,nsotot)
        CALL reallocate(ecoul_1c(iatom)%aVh2c0b%r_coef,1,nsotot,1,nsotot)
        CALL reallocate(ecoul_1c(iatom)%aVh3c0b%r_coef,1,nsotot,1,nsotot)

      END DO
    END DO

    ALLOCATE (vtriple(nkind,natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,"vtriple",nkind*natom)
    DO iat = 1,natom
      DO ikind = 1,nkind
        NULLIFY(vtriple(ikind,iat)%neighbor)
      ENDDO
    ENDDO

    IF (.NOT.my_tddft) THEN
       CALL set_qs_env(qs_env=qs_env,ecoul_1c=ecoul_1c,VH_3c_list=vtriple)
    ELSE
       p_env%hartree_local%ecoul_1c => ecoul_1c
       p_env%hartree_local%vtriple => vtriple
    END IF

  END SUBROUTINE init_coulomb_local

!******************************************************************************

  SUBROUTINE calculate_Vh_1center(vrad_h,vrad_s,rrad_h,rrad_s,rrad_0,rrad_z,grid_atom)

    REAL(dp), DIMENSION(:,:), INTENT(OUT)      :: vrad_h, vrad_s
    REAL(dp), DIMENSION(:,:), INTENT(IN)       :: rrad_0
    REAL(dp), DIMENSION(:), INTENT(IN)         :: rrad_z
    TYPE (rho_atom_coeff), DIMENSION(:),&
                                INTENT(IN)     :: rrad_h, rrad_s
    TYPE (grid_atom_type), POINTER             :: grid_atom

    REAL(dp)                                   :: I1_up, I1_down,&
                                                  I2_up, I2_down, &
                                                  prefactor, sum
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)     :: rho_1, rho_2
    REAL(dp), DIMENSION(:), POINTER             :: wr
    REAL(dp), DIMENSION(:,:), POINTER          :: r2l, oor2l
    INTEGER                                    :: ir, iso, ispin, l_ang, &
                                                  max_s_harm, nchannels, &
                                                  nr, nspins
!   ---------------------------------------------------------------------------

    nr = grid_atom%nr
    max_s_harm = SIZE(vrad_h,2)
    nspins = SIZE(rrad_h,1)
    nchannels = SIZE(rrad_0,2)

    r2l => grid_atom%rad2l
    oor2l => grid_atom%oorad2l
    wr => grid_atom%wr

    ALLOCATE(rho_1(nr,max_s_harm),rho_2(nr,max_s_harm))
    rho_1 = 0.0_dp
    rho_2 = 0.0_dp

!   Case lm = 0
!    rho_1(:,1) = rrad_z(:)
    rho_2(:,1) = rrad_0(:,1)

    DO iso = 2,nchannels
      rho_2(:,iso) = rrad_0(:,iso)      
    END DO

    DO iso = 1,max_s_harm
       DO ispin = 1,nspins
         rho_1(:,iso) = rho_1(:,iso) + rrad_h(ispin)%r_coef(:,iso)
         rho_2(:,iso) = rho_2(:,iso) + rrad_s(ispin)%r_coef(:,iso)
       END DO

       l_ang = indso(1,iso)
       prefactor = fourpi/(2._dp*l_ang+1._dp)

       rho_1(:,iso) = rho_1(:,iso)*wr(:)
       rho_2(:,iso) = rho_2(:,iso)*wr(:)
 
       I1_up = 0.0_dp
       I1_down = 0.0_dp
       I2_up = 0.0_dp
       I2_down = 0.0_dp

       I1_up = r2l(nr,l_ang)*rho_1(nr,iso)
       I2_up = r2l(nr,l_ang)*rho_2(nr,iso)

       DO ir = nr-1,1,-1
         I1_down = I1_down + oor2l(ir,l_ang+1)*rho_1(ir,iso)
         I2_down = I2_down + oor2l(ir,l_ang+1)*rho_2(ir,iso)
       END DO

       vrad_h(nr,iso) = vrad_h(nr,iso) + prefactor*&
                       (oor2l(nr,l_ang+1)*I1_up + r2l(nr,l_ang)*I1_down)
       vrad_s(nr,iso) = vrad_s(nr,iso) + prefactor*&
                       (oor2l(nr,l_ang+1)*I2_up + r2l(nr,l_ang)*I2_down)

       DO ir = nr-1,1,-1
         I1_up = I1_up + r2l(ir,l_ang)*rho_1(ir,iso)
         I1_down = I1_down -oor2l(ir,l_ang+1)*rho_1(ir,iso)
         I2_up = I2_up + r2l(ir,l_ang)*rho_2(ir,iso)
         I2_down = I2_down -oor2l(ir,l_ang+1)*rho_2(ir,iso)

         vrad_h(ir,iso) = vrad_h(ir,iso) + prefactor*&
                       (oor2l(ir,l_ang+1)*I1_up + r2l(ir,l_ang)*I1_down)
         vrad_s(ir,iso) = vrad_s(ir,iso) + prefactor*&
                       (oor2l(ir,l_ang+1)*I2_up + r2l(ir,l_ang)*I2_down)

       END DO

    END DO

    DEALLOCATE(rho_1,rho_2)

  END SUBROUTINE calculate_Vh_1center

!******************************************************************************

  SUBROUTINE Vh_1c_gg_integrals(qs_env,energy_hartree_1c,tddft,p_env)

    TYPE(qs_environment_type), POINTER         :: qs_env
    REAL(kind=dp), INTENT(out)                 :: energy_hartree_1c
    LOGICAL, INTENT(IN), OPTIONAL              :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL     :: p_env

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "Vh_1c_gg_integrals",&
                routineP = module_name//':'//routine_name

    TYPE(atomic_kind_type), DIMENSION(:), &
                            POINTER            :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER            :: atom_kind
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE(gto_basis_set_type), POINTER          :: orb_basis
    TYPE(grid_atom_type), POINTER              :: grid_atom
    TYPE(harmonics_atom_type), &
                           POINTER             :: harmonics
    TYPE(dft_control_type), POINTER            :: dft_control
    TYPE(ecoul_1center_type), DIMENSION(:), & 
                                 POINTER       :: ecoul_1c
    TYPE(rho_atom_type), DIMENSION(:), &
                                     POINTER   :: rho_atom_set
    TYPE(rho_atom_type),       POINTER         :: rho_atom
    TYPE(rho0_mpole_type), POINTER             :: rho0_mpole
    TYPE(rho0_atom_type), DIMENSION(:), &
                                  POINTER      :: rho0_atom_set
    TYPE(rhoz_type), DIMENSION(:), POINTER     :: rhoz_set
    TYPE(rho_atom_coeff),DIMENSION(:),&
                                 POINTER       :: r_h,r_s
    TYPE(rho_atom_coeff),DIMENSION(:),&
                                 POINTER       :: vr_h,vr_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                                      POINTER  :: cpc_h, cpc_s
    REAL(dp)                                   :: ec_tmp,es_tmp,&
                                                  e_1_h,e_1_s,&
                                                  gVg_h, gVg_s,&
                                                  gVg_0,&
                                                  ecoul_1_0,&
                                                  ecoul_1_z,&
                                                  ecoul_1_h,&
                                                  ecoul_1_s
    REAL(dp), DIMENSION(:),      POINTER       :: rrad_z, vrrad_z
    REAL(dp), DIMENSION(:),  ALLOCATABLE       :: gexp, sqrtwr
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: g0_h_w
    REAL(dp), DIMENSION(:,:), POINTER          :: gsph
    REAL(dp), DIMENSION(:,:),    POINTER       :: Vh1_h, Vh1_s
    REAL(dp), DIMENSION(:,:),    POINTER       :: rrad_0, vrrad_0
    REAL(dp), DIMENSION(:,:),    POINTER       :: aVh1b_h, aVh1b_s,&
                                                  aVh1b_0
    REAL(dp), DIMENSION(:,:), POINTER          :: zet
    REAL(dp), DIMENSION(:,:), POINTER          :: g0_h
    REAL(dp), DIMENSION(:,:,:),  POINTER       :: Qlm_gg
    REAL(dp), DIMENSION(:,:,:), POINTER        :: my_CG

    INTEGER                                    :: bo(2), handle, iat, iatom, icg, ikind, &
                                                  ipgf1, ipgf2, ir, iset1, iset2, &
                                                  is1, is2, iso, iso1, iso2, ispin, istat, &
                                                  l_ang, lmax_0, lmax0, m1, m2, max_nr, &
                                                  max_iso, max_iso_not0, maxso, mepos, n1,&
                                                  n2, nat, nchan_0, nkind, nr, nset, &
                                                  nsotot, nspins, num_pe
    INTEGER, DIMENSION(:),        POINTER      :: atom_list, lmin, lmax, npgf
    LOGICAL                                    :: my_tddft, paw_atom

!   -----------------------------------------------------------------------------

    CALL timeset("Vh_1c_gg_integrals","I","",handle)

    NULLIFY(dft_control,para_env)
    NULLIFY(atomic_kind_set,rho_atom_set,rho_atom,rho0_atom_set)
    NULLIFY(rho0_mpole,rhoz_set,ecoul_1c)
    NULLIFY(atom_kind,atom_list,grid_atom,harmonics)
    NULLIFY(orb_basis,lmin,lmax,npgf,zet)
    NULLIFY(gsph)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    atomic_kind_set=atomic_kind_set)


    IF (my_tddft) THEN
       rho_atom_set => p_env%local_rho_set%rho_atom_set
       rho0_atom_set=> p_env%local_rho_set%rho0_atom_set
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       ecoul_1c     => p_env%hartree_local%ecoul_1c
    ELSE
       CALL get_qs_env(qs_env=qs_env, &
                       rho_atom_set= rho_atom_set,&
                       rho0_atom_set=rho0_atom_set, &
                       rho0_mpole=rho0_mpole,&
                       rhoz_set=rhoz_set,&
                       ecoul_1c=ecoul_1c)
    END IF


    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             max_ngrid_rad=max_nr,maxg_iso_not0=max_iso)
    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,lmax_0=lmax_0)

    ALLOCATE(gexp(max_nr),STAT=istat)
    ALLOCATE(sqrtwr(max_nr),g0_h_w(max_nr,0:lmax_0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "gexp,sqrtwr,g0_h_w",(2*max_nr+max_nr*(1+lmax_0))*dp_size)    

!   Put to 0 the local hartree energy contribution from 1 center integrals
    energy_hartree_1c = 0.0_dp

!   Here starts the loop over all the atoms
    DO ikind = 1,nkind

      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                           orb_basis_set=orb_basis,&
                           natom=nat,grid_atom=grid_atom,&
                           harmonics=harmonics,ngrid_rad=nr,&
                           max_iso_not0=max_iso_not0,paw_atom=paw_atom)

      CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                             maxso=maxso,npgf=npgf,&
                             nset=nset,zet=zet)

      nsotot = maxso*nset
      CALL reallocate(gsph,1,max_nr,1,nsotot)

      NULLIFY(Qlm_gg,g0_h)
      CALL get_rho0_mpole(rho0_mpole=rho0_mpole, ikind=ikind, &
                           l0_ikind=lmax0, &
                           Qlm_gg=Qlm_gg, g0_h=g0_h)


      nchan_0 = nsoset(lmax0)

      IF(nchan_0 > max_iso_not0) CALL stop_program(routine_name,module_name,__LINE__,&
                        "channels for rho0 > # max of spherical harmonics")

      NULLIFY(rrad_z,my_CG)
      my_CG  => harmonics%my_CG

!     set to zero temporary arrays
      sqrtwr=0.0_dp
      g0_h_w=0.0_dp
      gexp=0.0_dp
      gsph=0.0_dp

      sqrtwr(1:nr) = sqrt(grid_atom%wr(1:nr))
      DO l_ang = 0,lmax0
        g0_h_w(1:nr,l_ang) = g0_h(1:nr,l_ang)*grid_atom%wr(1:nr)
      END DO 

      m1 = 0
      DO iset1 = 1,nset
        n1 = nsoset(lmax(iset1))
        DO ipgf1  = 1,npgf(iset1)
          gexp(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))*sqrtwr(1:nr)
          DO is1 = nsoset(lmin(iset1)-1)+1,nsoset(lmax(iset1))
            iso = is1 + (ipgf1-1)*n1 + m1
            l_ang = indso(1,is1)
            gsph(1:nr,iso) = grid_atom%rad2l(1:nr,l_ang)*gexp(1:nr) 
          END DO  ! is1
        END DO  ! ipgf1
        m1 = m1 + maxso
      END DO  ! iset1

!     Initialize to zero the local ks integrals
      DO iat = 1,nat
        iatom = atom_list(iat)
        ecoul_1c(iatom)%Vh1_h%r_coef = 0.0_dp
        ecoul_1c(iatom)%Vh1_s%r_coef = 0.0_dp
        ecoul_1c(iatom)%aVh1b_h%r_coef = 0.0_dp
        ecoul_1c(iatom)%aVh1b_s%r_coef = 0.0_dp
        ecoul_1c(iatom)%aVh1b_0%r_coef = 0.0_dp
      END DO

!     Distribute the atoms of this kind
      num_pe = para_env%num_pe
      mepos  = para_env%mepos
      bo = get_limit( nat, num_pe, mepos )

      DO iat = bo(1), bo(2) !1,nat
        iatom = atom_list(iat)

        NULLIFY(rho_atom,r_h,r_s,vr_h,vr_s,cpc_h,cpc_s)
        rho_atom => rho_atom_set(iatom)
        CALL get_rho_atom(rho_atom=rho_atom,&
                   rho_rad_h=r_h,rho_rad_s=r_s,&
               vrho_rad_h=vr_h,vrho_rad_s=vr_s,&
                        cpc_h=cpc_h,cpc_s=cpc_s)

        NULLIFY(Vh1_h,Vh1_s,rrad_0)
        CALL  get_ecoul_1c(ecoul_1c,iatom=iatom,Vh1_h=Vh1_h,Vh1_s=Vh1_s)

        rrad_0 => rho0_atom_set(iatom)%rho0_rad_h%r_coef
        vrrad_0 => rho0_atom_set(iatom)%vrho0_rad_h%r_coef

        IF (.NOT.my_tddft) THEN
             rrad_z  => rhoz_set(ikind)%r_coef
             vrrad_z => rhoz_set(ikind)%vr_coef
             Vh1_h(:,1) = vrrad_z(:) 
        END IF

        DO iso = 1,nchan_0
          Vh1_s(:,iso) = vrrad_0(:,iso)
        END DO

        DO ispin = 1,nspins
          DO iso =1,max_iso_not0
            Vh1_h(:,iso) = Vh1_h(:,iso) + vr_h(ispin)%r_coef(:,iso)
            Vh1_s(:,iso) = Vh1_s(:,iso) + vr_s(ispin)%r_coef(:,iso)
          END DO
        END DO

!        CALL calculate_Vh_1center(Vh1_h,Vh1_s,r_h,r_s,rrad_0,rrad_z,grid_atom)

!       Calculate the contributions to Ecoul coming from Vh1_h*rhoz 

        ecoul_1_z = 0.0_dp
        IF(.NOT. my_tddft) THEN
           DO ir = 1,nr
              ecoul_1_z = ecoul_1_z + Vh1_h(ir,1)*rrad_z(ir)*grid_atom%wr(ir)
           END DO
           ecoul_1_z = 0.5_dp*ecoul_1_z
        END IF

!       Calculate the contributions to Ecoul coming from  Vh1_s*rho0
        ecoul_1_0 = 0.0_dp
        DO iso = 1,nchan_0
          ec_tmp = 0.0_dp
          DO ir = 1,nr
            ec_tmp = ec_tmp + Vh1_s(ir,iso)*rrad_0(ir,iso)*&
                              grid_atom%wr(ir)
          END DO  ! ir
          ecoul_1_0 = ecoul_1_0 + ec_tmp
        END DO  ! iso
        ecoul_1_0 = 0.5_dp*ecoul_1_0

        CALL set_ecoul_1c(ecoul_1c,iatom,ecoul_1_z=ecoul_1_z,ecoul_1_0=ecoul_1_0)

        energy_hartree_1c = energy_hartree_1c +  ecoul_1_z -  ecoul_1_0

        IF(.NOT. paw_atom) CYCLE

!       Calculate the integrals of the potential with 2 primitives
        NULLIFY(aVh1b_h,aVh1b_s,aVh1b_0)
        CALL  get_ecoul_1c(ecoul_1c,iatom=iatom,&
                           aVh1b_h=aVh1b_h,aVh1b_s=aVh1b_s,aVh1b_0=aVh1b_0)

        m1 = 0
        DO iset1 = 1,nset
          n1 = nsoset(lmax(iset1))
          DO ipgf1  = 1,npgf(iset1)
            m2 = 0
            DO iset2 = 1,nset !iset1,nset 
              n2 = nsoset(lmax(iset2))
              DO ipgf2  = 1,npgf(iset2)
!               with contributions to  V1_s*rho0
                DO iso = 1,nchan_0
                  l_ang = indso(1,iso)
                  DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                    is1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                    is2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                    iso1 = is1 + n1*(ipgf1-1) + m1
                    iso2 = is2 + n2*(ipgf2-1) + m2 
                    gVg_h = 0.0_dp
                    gVg_s = 0.0_dp
                    gVg_0 = 0.0_dp

                    DO ir = 1,nr
                       gVg_h = gVg_h + gsph(ir,iso1)*gsph(ir,iso2)* &
                                       Vh1_h(ir,iso)
                       gVg_s = gVg_s + gsph(ir,iso1)*gsph(ir,iso2)* &
                                       Vh1_s(ir,iso)
                       gVg_0 = gVg_0 + Vh1_s(ir,iso)*g0_h_w(ir,l_ang)
                    END DO  ! ir

                    aVh1b_h(iso1,iso2) = aVh1b_h(iso1,iso2) + &
                                         gVg_h*my_CG(is1,is2,iso)
                    aVh1b_s(iso1,iso2) = aVh1b_s(iso1,iso2) + &
                                         gVg_s*my_CG(is1,is2,iso)
                    aVh1b_0(iso1,iso2) = aVh1b_0(iso1,iso2) + &
                                         Qlm_gg(iso1,iso2,iso)*gVg_0

                  END DO  !icg
                END DO  ! iso
!               without contributions to  V1_s*rho0 
                DO iso = nchan_0+1,max_iso_not0
                  DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                    is1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                    is2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                    iso1 = is1 + n1*(ipgf1-1) + m1
                    iso2 = is2 + n2*(ipgf2-1) + m2 
                    gVg_h = 0.0_dp
                    gVg_s = 0.0_dp

                    DO ir = 1,nr
                       gVg_h = gVg_h + gsph(ir,iso1)*gsph(ir,iso2)* &
                                       Vh1_h(ir,iso)
                       gVg_s = gVg_s + gsph(ir,iso1)*gsph(ir,iso2)* &
                                       Vh1_s(ir,iso)
                    END DO  ! ir

                    aVh1b_h(iso1,iso2) = aVh1b_h(iso1,iso2) + &
                                         gVg_h*my_CG(is1,is2,iso)
                    aVh1b_s(iso1,iso2) = aVh1b_s(iso1,iso2) + &
                                         gVg_s*my_CG(is1,is2,iso)

                  END DO  !icg
                END DO  ! iso 
              END DO  ! ipgf2
              m2 = m2 + maxso
            END DO  ! iset2
          END DO  ! ipgf1
          m1 = m1 + maxso
        END DO  !iset1

!       Symmetry
!        DO iset1 = 1,nset
!          DO is1 = 1,maxso
!           iso1 = is1 + (iset1 - 1)*maxso
!            DO iset2 = iset1+1,nset
!              DO is2 = 1,maxso
!                iso2 = is2 + (iset2 - 1)*maxso
!                aVh1b_h(iso2,iso1) = aVh1b_h(iso1,iso2)
!                aVh1b_s(iso2,iso1) = aVh1b_s(iso1,iso2)
!                aVh1b_0(iso2,iso1) = aVh1b_0(iso1,iso2)
!              END DO
!            END DO
!          END DO
!        END DO  

!       Calculate the contributions to Ecoul coming from Vh1_h*rho1_h and Vh1_s*rho1_s

        ecoul_1_h = 0.0_dp
        ecoul_1_s = 0.0_dp
        DO ispin = 1,nspins
           ecoul_1_h  =  ecoul_1_h + &
                      trace_r_AxB(aVh1b_h,nsotot,&
                          cpc_h(ispin)%r_coef,nsotot,&
                          nsotot,nsotot) 

           ecoul_1_s  =  ecoul_1_s + &
                      trace_r_AxB(aVh1b_s,nsotot,&
                          cpc_s(ispin)%r_coef,nsotot,&
                          nsotot,nsotot) 
        END DO  ! ispin

        ecoul_1_h = 0.5_dp*ecoul_1_h
        ecoul_1_s = 0.5_dp*ecoul_1_s


        CALL set_ecoul_1c(ecoul_1c=ecoul_1c,iatom=iatom,&
                          ecoul_1_h=ecoul_1_h,ecoul_1_s=ecoul_1_s)

        energy_hartree_1c = energy_hartree_1c + ecoul_1_h - ecoul_1_s 

      END DO  ! iat

    END DO  ! ikind 

    nat = SIZE(ecoul_1c,1)
    DO iatom = 1,nat
      CALL mp_sum(ecoul_1c(iatom)%aVh1b_h%r_coef,para_env%group)
      CALL mp_sum(ecoul_1c(iatom)%aVh1b_s%r_coef,para_env%group)
      CALL mp_sum(ecoul_1c(iatom)%aVh1b_0%r_coef,para_env%group)
    END DO
    CALL mp_sum(energy_hartree_1c,para_env%group)

    DEALLOCATE(gexp,gsph,g0_h_w,sqrtwr,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,&
            "gexp,gsph,g0_h_w,sqrtwr")
    ENDIF

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE Vh_1c_gg_integrals

!******************************************************************************

END MODULE hartree_local_methods

!******************************************************************************
