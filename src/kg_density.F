!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_density [1.0] *
!!
!!   NAME
!!     kg_density
!!
!!   FUNCTION
!!     Calculate the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_density

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   dbl,&
                   wp_size => dp_size

  USE atomic_kind_types,      ONLY: atomic_kind_type,&
                                    get_atomic_kind,&
                                    get_atomic_kind_set
  USE basis_set_types,        ONLY: get_gto_basis_set,&
                                    gto_basis_set_type
  USE coefficient_types,      ONLY: coeff_allocate, coeff_deallocate, &
                                    coeff_add,coeff_sumup,&
                                    coeff_copy,&
                                    coeff_transform_space,&
                                    coeff_type,&
                                    coeff_zero
  USE realspace_grid_types,   ONLY: realspace_grid_type, rs_grid_allocate,  &
                                    rs_grid_deallocate, &
                                    rs_grid_setup, rs_pw_transfer, rs_grid_zero
  USE cube_utils,             ONLY: cube_info_type,&
                                    return_cube
  USE dft_types,              ONLY: dft_control_type
  USE gaussian_gridlevels,    ONLY: gaussian_gridlevel,&
                                    gridlevel_info_type
  USE global_types,           ONLY: global_environment_type
  USE l_utils,                ONLY: l_info_type,&
                                    return_l_info
  USE mathconstants,          ONLY: pi,twopi
  USE message_passing,        ONLY: mp_sum
  USE orbital_pointers,       ONLY: coset,indco,ncoset,nso
  USE particle_types,         ONLY: particle_type
  USE pw_types,               ONLY: COMPLEXDATA1D,COMPLEXDATA3D,&
                                    REALDATA3D,REALSPACE,&
                                    RECIPROCALSPACE
  USE kg_environment_types,   ONLY: get_kg_env,&
                                    kg_environment_type
  USE kg_force_types,              ONLY: kg_force_type
  USE qs_interactions,        ONLY: exp_radius_very_extended
  USE simulation_cell,        ONLY: cell_type,&
                                    pbc
  USE termination,            ONLY: stop_memory,&
                                    stop_program
  USE timings,                ONLY: timeset,&
                                    timestop

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: calculate_density, calculate_epc_density, &
            calculate_epc_rspace_forces, calculate_p_density, &
            calculate_v_rspace_forces, calculate_vp_rspace_forces

! *** Public functions ***

  PUBLIC :: calculate_total_density

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_density(rho,rho_gspace,&
                                    mgrid_rspace,mgrid_gspace,&
                                    gridlevel_info,cube_info,l_info,&
                                    total_rho, kg_env,globenv)

!   REAL(dbl), DIMENSION(:), POINTER             :: eigenvector   
    REAL(dbl), INTENT(OUT)                     :: total_rho
    TYPE(coeff_type), INTENT(INOUT)             :: rho,rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(kg_environment_type), INTENT(IN)       :: kg_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_density (MODULE kg_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER            :: part

    REAL(dbl) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ipgf,iset,ishell,istat,j,jatom,jco,jkind,jpgf,jset,&
                jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,ncob,&
                nkind,nseta,nsetb,sgfa,sgfb

    REAL(dbl), DIMENSION(3) :: ra,rab,rb

    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pab,work

    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(dbl), DIMENSION(:,:), POINTER :: p_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    INTEGER, DIMENSION(3)             :: nsmax
    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER  :: rs_rho
    INTEGER :: offset,natom

!   ---------------------------------------------------------------------------
    CALL timeset("calculate_density","I","",handle)

    ! only some processors have non-zero elements. The sum of eigenvector on 
    ! all cpus is the real eigenvector. combined with the rs->pw this will give the
    ! real density

    ALLOCATE(rs_rho(gridlevel_info%ngrid_levels),STAT=istat)
    IF (istat.ne.0) then
       CALL stop_program("calculate_rho_elec","rs_rho")
    ENDIF
    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=1000

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=part ) 

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*wp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*wp_size)

    offset=0

    DO iatom=1,natom

      CALL get_atomic_kind(atomic_kind=part(iatom)%atomic_kind,&
                               orb_basis_set=orb_basis_set)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      ra(:) = pbc(part(iatom)%r,cell)
      rb(:) = 0.0_dbl
      rab(:) = 0.0_dbl
      rab2  = 0.0_dbl
      dab   = 0.0_dbl

      DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

! In KG simulation, density is frozen. The PDM effective charges are
! already included in sphi and are the ones read from the POTENTIAL file in
! the initialization. Thus all prefactors (defined in work(:,:)) are set to 1.0  

         do i=1,nsgfa(iset)
            work(i,1)=1.0_dbl 
         enddo

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dbl,pab(1,1),SIZE(pab,1))

         DO ipgf=1,npgfa(iset)


            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))


            scale = -1.0_dbl
            zetp = zeta(ipgf,iset)
            igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

            CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_dbl,0,&
                        ra,rab,rab2,scale,pab(na1:na2,1:1),&
                        rs_rho(igrid_level),cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_pw_transfer(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw,"FORWARD")
       call rs_grid_deallocate(rs_rho(igrid_level))
    ENDDO

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),mgrid_gspace(igrid_level))
      CALL coeff_add(rho_gspace,mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    deallocate(rs_rho)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    total_rho = calculate_total_density(rho)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_density
! *****************************************************************************

  FUNCTION calculate_total_density(rho) RESULT(total_rho)

    TYPE(coeff_type), INTENT(IN), TARGET :: rho

    REAL(dbl) :: total_rho

!   *** Local variables ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION calculate_total_density (MODULE kg_density)"

!   *** Local variables ***

    INTEGER, DIMENSION(:), POINTER :: lb_grid

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho%pw%cc3d)) THEN
      lb_grid => rho%pw%pw_grid%bounds(1,:)
      total_rho = rho%pw%pw_grid%vol*rho%pw%cc3d(lb_grid(1),&
                                                  lb_grid(2),&
                                                  lb_grid(3))
    ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
      total_rho = rho%pw%pw_grid%dvol*SUM(rho%pw%cr3d)
    ELSE IF (ASSOCIATED(rho%pw%cc)) THEN
      IF ( rho%pw%pw_grid%have_g0 ) THEN
        total_rho = rho%pw%pw_grid%vol*rho%pw%cc(1)
      ELSE
        total_rho = 0._dbl
      END IF
    ELSE
      CALL stop_program(routine,"No density coefficients available")
    END IF
    IF (rho%pw%pw_grid%para%mode.ne.0) THEN
       CALL mp_sum(total_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_density

! *****************************************************************************

  SUBROUTINE calculate_epc_density (rho_core, cube_info, l_info, total_rho, &
                                    kg_env, globenv)
! *****************************************************************************

  IMPLICIT NONE

!   *** Arguments ***

    TYPE(coeff_type), INTENT(OUT)             :: rho_core
    TYPE(cube_info_type), INTENT(IN)          :: cube_info
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(l_info_type), INTENT(IN)             :: l_info
    TYPE(kg_environment_type), INTENT(IN)     :: kg_env
    REAL(dbl), INTENT(OUT)                     :: total_rho

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    TYPE(dft_control_type), POINTER :: dft_control

    REAL(dbl) :: alpha,eps_rho_rspace,scale
    INTEGER  :: atom_a,handle,iatom,ikind,natom

    REAL(dbl), DIMENSION(3)   :: ra
    REAL(dbl), DIMENSION(1,1) :: pab

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: part
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    TYPE(realspace_grid_type)                     :: rs_rho
    INTEGER, DIMENSION(3)                         :: nsmax
    TYPE(coeff_type)                              :: rhoc_r

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_epc_density","I","",handle)

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=part)

    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=100
    call rs_grid_setup(rs_rho,rho_core%pw%pw_grid,nsmax)
    call rs_grid_allocate(rs_rho)
    call rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace
    scale = -1.0_dbl

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      DO iatom=1,natom

        atom_a = atom_list(iatom)
        ra(:) = pbc(part(atom_a)%r,cell)
        IF(rs_rho%group_head) &
        CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dbl,0,ra,&
                                          (/0.0_dbl,0.0_dbl,0.0_dbl/),&
                                          0.0_dbl,scale,pab,rs_rho,&
                                          cube_info,l_info,eps_rho_rspace)

      END DO

    END DO

    CALL coeff_allocate(rhoc_r,rho_core%pw%pw_grid,REALDATA3D)
    rhoc_r%pw%in_space = REALSPACE

    call rs_pw_transfer(rs_rho,rhoc_r%pw,"FORWARD")
    call rs_grid_deallocate(rs_rho)

    total_rho = calculate_total_density(rhoc_r)

    CALL coeff_transform_space(rhoc_r,rho_core)

    CALL coeff_deallocate(rhoc_r)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_epc_density

! *****************************************************************************

  SUBROUTINE calculate_p_density(eigenvector,rho,rho_gspace,&
                                    mgrid_rspace,mgrid_gspace,&
                                    gridlevel_info,cube_info,l_info,&
                                    total_rho, kg_env,globenv)

    REAL(dbl), DIMENSION(:), POINTER             :: eigenvector   
    REAL(dbl), INTENT(OUT)                     :: total_rho
    TYPE(coeff_type), INTENT(INOUT)             :: rho,rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(kg_environment_type), INTENT(IN)       :: kg_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_p_density (MODULE kg_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER            :: part

    REAL(dbl) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ipgf,iset,ishell,istat,j,jatom,jco,jkind,jpgf,jset,&
                jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,ncob,&
                nkind,nseta,nsetb,sgfa,sgfb

    REAL(dbl), DIMENSION(3) :: ra,rab,rb

    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: pab,work

    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(dbl), DIMENSION(:,:), POINTER :: p_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    INTEGER, DIMENSION(3)             :: nsmax
    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER  :: rs_rho
    INTEGER :: offset,natom

!   ---------------------------------------------------------------------------
    CALL timeset("calculate_p_density","I","",handle)

    ! only some processors have non-zero elements. The sum of eigenvector on 
    ! all cpus is the real eigenvector. combined with the rs->pw this will give the
    ! real density

    ALLOCATE(rs_rho(gridlevel_info%ngrid_levels),STAT=istat)
    IF (istat.ne.0) then
       CALL stop_program("calculate_rho_elec","rs_rho")
    ENDIF
    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=1000

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=part ) 

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*wp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*wp_size)

    offset=0

    DO iatom=1,natom

      CALL get_atomic_kind(atomic_kind=part(iatom)%atomic_kind,&
                               aux_basis_set=aux_basis_set)

      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      ra(:) = pbc(part(iatom)%r,cell)
      rb(:) = 0.0_dbl
      rab(:) = 0.0_dbl
      rab2  = 0.0_dbl
      dab   = 0.0_dbl

      DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

! polarization density is not frozen. The coefficients of each contracted
! gaussian polarization functions (i.e. the dynamical variables associated to
! polarization), stored in eigenvector(:), are here multiplied for the contraction
! coefficients in sphi and stored in  work(:,:). Then a decontraction is performed
! in order to operate on each pgf.

         do i=1,nsgfa(iset)
            work(i,1)=eigenvector(offset+i)
         enddo

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dbl,pab(1,1),SIZE(pab,1))

         DO ipgf=1,npgfa(iset)


            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))


            scale = 1.0_dbl
            zetp = zeta(ipgf,iset)
            igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

            CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_dbl,0,&
                        ra,rab,rab2,scale,pab(na1:na2,1:1),&
                        rs_rho(igrid_level),cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_pw_transfer(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw,"FORWARD")
       call rs_grid_deallocate(rs_rho(igrid_level))
    ENDDO

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),mgrid_gspace(igrid_level))
      CALL coeff_add(rho_gspace,mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    deallocate(rs_rho)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    total_rho = calculate_total_density(rho)

    CALL timestop(0.0_dbl,handle)

   END SUBROUTINE calculate_p_density

! *****************************************************************************

  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,&
                                          rsgrid,cube_info,l_info,eps_rho_rspace)

    TYPE(realspace_grid_type), INTENT(INOUT), TARGET :: rsgrid
    TYPE(cube_info_type), INTENT(IN)     :: cube_info
    TYPE(l_info_type), INTENT(IN)        :: l_info
    REAL(dbl), INTENT(IN)                 :: eps_rho_rspace,rab2,scale,zeta,zetb
    INTEGER, INTENT(IN)                  :: la_max,la_min,lb_max,lb_min
    REAL(dbl), DIMENSION(3), INTENT(IN)   :: ra,rab

    REAL(dbl), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab

!   *** Local variables ***

    REAL(dbl) :: f,prefactor,radius,zetp
    INTEGER  :: i,ico,ig,j,jco,jg,k,kg

    REAL(dbl), DIMENSION(3) :: dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(dbl), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER :: lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER :: coef_max,lxa,lxyz,l,gridbounds(2,3)
    INTEGER, pointer,dimension(:) :: ly_max,lz_max
    REAL(dbl) :: a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,&
                rpg,s,pg

    REAL(dbl), pointer :: dpx(:,:)
    REAL(dbl), pointer :: pzyx(:), alpha(:,:)
    integer,  pointer :: ipzyx(:,:,:,:,:,:)
    REAL(dbl), pointer :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    INTEGER, pointer, dimension(:)  :: sphere_bounds
    REAL(dbl)  :: xa,xb,xap,xbp,polypart,cutoff
    INTEGER  :: start,length,offset,cmax
    integer  handle

!   ---------------------------------------------------------------------------


    coef_max=la_max+lb_max+1

    zetp      = zeta + zetb
    f         = zetb/zetp
    prefactor = scale*EXP(-zeta*f*rab2)
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)  
    rb(:)     = ra(:)+rab(:)
    cutoff    = 0.0_dbl

    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp,&
                                    zetp,eps_rho_rspace,prefactor,cutoff)

    IF (radius .eq. 0.0_dbl ) THEN
      RETURN
    END IF

!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    !dr(:) = pw%pw_grid%dr(:)
    !ng(:) = pw%pw_grid%npts(:)
    !grid => pw%cr3d(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min,la_max,lb_min,lb_max,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** position of the gaussian product

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),dbl)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    DO i=1,3
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         if (start+length.ge.ub_cube(i)) EXIT
         start=start+length+1
        END DO
    ENDDO

!   *** initialise the p terms and loop logic
    lxyz=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             pzyx(lxyz)=prefactor*pab(ico,jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

 
!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,dbl)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max
       zbp=1.0_dbl
       DO lzb=0,lb_max
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO 
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,dbl)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max
       ybp=1.0_dbl
       DO lyb=0,lb_max
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

 
!   *** make the alpha matrix ***
    alpha(:,:)=0.0_dbl
    lx=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       lx=lx+1
       binomial_k_lxa=1.0_dbl
       a=1.0_dbl
       DO k=0,lxa
        binomial_l_lxb=1.0_dbl
        b=1.0_dbl
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*real(lxb-l,dbl)/real(l+1,dbl)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*real(lxa-k,dbl)/real(k+1,dbl)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO

    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,dbl)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

!   *** do the loop over the grid
    CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax), &
                        polz(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
                        lxyz_max,coef_max, &
                        cmax,ly_max(1),lz_max(1),gridbounds(1,1),&
                        map(-cmax,1),sphere_bounds(1))

  END SUBROUTINE collocate_pgf_product_rspace

! *****************************************************************************

  SUBROUTINE calculate_epc_rspace_forces(v_rspace,cube_info,l_info,kg_env,force)
   
    IMPLICIT NONE

!   *** Arguments ***

    TYPE(coeff_type), INTENT(INOUT)          :: v_rspace
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(INOUT) :: kg_env
    TYPE(kg_force_type), DIMENSION(:), POINTER    :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_epc_rspace_forces (MODULE kg_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    TYPE(dft_control_type), POINTER :: dft_control

    TYPE(realspace_grid_type) :: rs_v
    REAL(dbl)                  :: alpha_core_charge,ccore_charge,eps_rho_rspace
    INTEGER                   :: atom_a,handle,i,iatom,ikind,istat,&
                                 natom_of_kind

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list

    REAL(dbl), DIMENSION(3)     :: force_a,force_b,ra
    INTEGER, DIMENSION(3)      :: nsmax
    REAL(dbl), DIMENSION(1,1)   :: hab,pab

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_epc_rspace_forces","I","",handle)

    nsmax(:) = 100
    CALL rs_grid_setup(rs_v,v_rspace%pw%pw_grid,nsmax)
    CALL rs_grid_allocate(rs_v)
    CALL rs_pw_transfer(rs_v,v_rspace%pw,"BACKWARD")

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=particle_set) 

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

      pab(1,1) = ccore_charge

      DO iatom=1,natom_of_kind

        atom_a = atom_list(iatom)
        ra(:)  = pbc(particle_set(atom_a)%r,cell)

        hab(1,1) = 0.0_dbl
        force_a(:) = 0.0_dbl
        force_b(:) = 0.0_dbl

        CALL integrate_pgf_product_rspace(0,alpha_core_charge,0,&
                                          0,0.0_dbl,0,&
                                          ra,(/0.0_dbl,0.0_dbl,0.0_dbl/),0.0_dbl,&
                                          rs_v,cube_info,l_info,&
                                          hab(1,1),pab(1,1),1.0_dbl,&
                                          eps_rho_rspace,&
                                          .TRUE.,force_a,force_b)

        force(ikind)%f_hartree_core(:,iatom) =&
          force(ikind)%f_hartree_core(:,iatom) + force_a(:)

      END DO

    END DO

    CALL rs_grid_deallocate(rs_v)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_epc_rspace_forces

! *****************************************************************************

  SUBROUTINE calculate_v_rspace_forces(v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,&
                                kg_env,force,globenv,force_type)
    IMPLICIT NONE

!   ***  Arguments   ***

    TYPE(coeff_type)                            :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    INTEGER, INTENT (IN)                        :: auxbas_grid
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(kg_force_type),DIMENSION (:), POINTER  :: force 
    TYPE(kg_environment_type), INTENT(IN)       :: kg_env
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL      :: force_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_v_rspace_forces (MODULE kg_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    REAL(dbl) :: dab,eps_gvg_rspace,&
                rab2,scale,vratio,zetp
    INTEGER  :: atom_a,atom_b,handle,iatom,ico,igrid_level,ikind,&
                ipgf,iset,ishell,istat,i,&
                maxco,maxsgf_set,na1,na2,natom,&
                natom_of_kind,nb1,nb2,ncoa,nkind,&
                nseta,sgfa

    REAL(dbl), DIMENSION(3) :: force_a,force_b,ra,rab,rb

    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: hab,pab,work

    INTEGER, DIMENSION(:), POINTER                :: atom_list

    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,npgfa,nsgfa
    REAL(dbl), DIMENSION(:,:), POINTER :: sphi_a,rpgfa,zeta
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa

    INTEGER, DIMENSION(3)             :: nsmax
    LOGICAL:: calculate_forces

    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER :: rs_v

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_v_rspace_forces","I","",handle)

    vratio = 1.0_dbl
    calculate_forces =.TRUE.

!   *** Get the potential on the subgrids in real space, via fft ***

    CALL coeff_transform_space(v_rspace,mgrid_gspace(auxbas_grid))

    DO igrid_level=1,gridlevel_info%ngrid_levels

      IF ( igrid_level /= auxbas_grid )  &
           CALL coeff_copy(mgrid_gspace(auxbas_grid),mgrid_gspace(igrid_level))
      CALL coeff_transform_space(mgrid_gspace(igrid_level),mgrid_rspace(igrid_level))

!     *** Multiply by the grid volume element ratio ***
      IF ( igrid_level /= auxbas_grid ) THEN
         scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                 mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
         mgrid_rspace(igrid_level)%pw%cr3d = scale*mgrid_rspace(igrid_level)%pw%cr3d
      END IF

    END DO

    ALLOCATE(rs_v(gridlevel_info%ngrid_levels),STAT=istat)
    IF (istat.ne.0) then
        CALL stop_memory("calculate_v_rspace_forces","rs_v")
    ENDIF

    nsmax(:)=100

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_setup(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, &
            nsmax)
       CALL rs_grid_allocate(rs_v(igrid_level))
       CALL rs_pw_transfer(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw,&
                           "BACKWARD")
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=particle_set )

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_gvg_rspace = dft_control%kg_control%eps_gvg_rspace

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*wp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*wp_size)

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*wp_size)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      DO iatom = 1, natom_of_kind

        atom_a = atom_list(iatom)
        ra(:)  = pbc(particle_set(atom_a)%r,cell)

        force_a(:) = 0.0_dbl
        force_b(:) = 0.0_dbl
        rb(:) = 0.0_dbl
        rab(:) = 0.0_dbl
        rab2  = 0.0_dbl
        dab   = 0.0_dbl

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          do i=1,nsgfa(iset)
            work(i,1)=1.0_dbl
          enddo

          CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dbl,pab(1,1),SIZE(pab,1))

           DO ipgf=1,npgfa(iset)

             hab(:,:) = 0.0_dbl

             na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
             na2 = ipgf*ncoset(la_max(iset))

             zetp = zeta(ipgf,iset) 

             igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

             CALL integrate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0, 0.0_dbl,0,&
                        ra,rab,rab2,rs_v(igrid_level),&
                        cube_info(igrid_level),l_info,&
                        hab(na1:na2,1:1),pab(na1:na2,1:1),&
                        vratio,eps_gvg_rspace,calculate_forces,&
                        force_a,force_b)

           END DO

        END DO

!   *** Update forces ***

        IF (PRESENT (force_type)) THEN
          IF (force_type=='hartree') THEN
            force(ikind)%f_hartree(:,iatom) =  force_a(:)
          ELSE IF (force_type=='xc') THEN
            force(ikind)%f_xc(:,iatom) = force_a(:)
          END IF
        ELSE
          force(ikind)%f_rho(:,iatom) =&
               force(ikind)%f_rho(:,iatom) + 1.0_dbl*force_a(:)
        END IF
      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_deallocate(rs_v(igrid_level))
    ENDDO

!   *** Release work storage ***

    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab")

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_v_rspace_forces

! *****************************************************************************

  SUBROUTINE calculate_vp_rspace_forces(eigenvector,eigenforce, &
                                v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,kg_env,&
                                force,globenv,force_type)
    IMPLICIT NONE

!   ***  Arguments   ***

    REAL(dbl), DIMENSION(:), POINTER            :: eigenvector
    REAL(dbl), DIMENSION(:), POINTER            :: eigenforce
    TYPE(coeff_type)                            :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    INTEGER, INTENT (IN)                        :: auxbas_grid
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(kg_force_type),DIMENSION (:), POINTER  :: force 
    TYPE(kg_environment_type), INTENT(IN)       :: kg_env
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL      :: force_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_vp_rspace_forces (MODULE kg_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    REAL(dbl) :: dab,eps_gvg_rspace,&
                rab2,scale,vratio,zetp
    INTEGER  :: atom_a,atom_b,handle,iatom,ico,igrid_level,ikind,&
                ipgf,iset,ishell,istat,i,&
                maxco,maxsgf_set,na1,na2,natom,&
                natom_of_kind,nb1,nb2,ncoa,nkind,&
                nseta,sgfa
    INTEGER :: offset

    REAL(dbl), DIMENSION(3) :: force_a,force_b,ra,rab,rb

    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: hab,pab,work

    INTEGER, DIMENSION(:), POINTER                :: atom_list

    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,npgfa,nsgfa
    REAL(dbl), DIMENSION(:,:), POINTER :: sphi_a,rpgfa,zeta
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa

    INTEGER, DIMENSION(3)             :: nsmax
    LOGICAL:: calculate_forces

    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER :: rs_v

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_v_rspace_forces","I","",handle)

    vratio = 1.0_dbl
    calculate_forces =.TRUE.

!   *** Get the potential on the subgrids in real space, via fft ***

    CALL coeff_transform_space(v_rspace,mgrid_gspace(auxbas_grid))

    DO igrid_level=1,gridlevel_info%ngrid_levels

      IF ( igrid_level /= auxbas_grid )  &
           CALL coeff_copy(mgrid_gspace(auxbas_grid),mgrid_gspace(igrid_level))
      CALL coeff_transform_space(mgrid_gspace(igrid_level),mgrid_rspace(igrid_level))

!     *** Multiply by the grid volume element ratio ***
      IF ( igrid_level /= auxbas_grid ) THEN
         scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                 mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
         mgrid_rspace(igrid_level)%pw%cr3d = scale*mgrid_rspace(igrid_level)%pw%cr3d
      END IF

    END DO

    ALLOCATE(rs_v(gridlevel_info%ngrid_levels),STAT=istat)
    IF (istat.ne.0) then
        CALL stop_memory("calculate_v_rspace_forces","rs_v")
    ENDIF

    nsmax(:)=100

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_setup(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, &
            nsmax)
       CALL rs_grid_allocate(rs_v(igrid_level))
       CALL rs_pw_transfer(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw,&
                           "BACKWARD")
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=cell,&
                    dft_control=dft_control,&
                    part=particle_set )

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_gvg_rspace = dft_control%kg_control%eps_gvg_rspace

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*wp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*wp_size)

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*wp_size)
    
    offset = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
                           aux_basis_set=aux_basis_set)

      IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      DO iatom = 1, natom_of_kind

        atom_a = atom_list(iatom)
        ra(:)  = pbc(particle_set(atom_a)%r,cell)

        force_a(:) = 0.0_dbl
        force_b(:) = 0.0_dbl
        rb(:) = 0.0_dbl
        rab(:) = 0.0_dbl
        rab2  = 0.0_dbl
        dab   = 0.0_dbl

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          do i=1,nsgfa(iset)
             work(i,1)=eigenvector(offset+i)
          enddo

          CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dbl,pab(1,1),SIZE(pab,1))

           hab(:,:) = 0.0_dbl

           DO ipgf=1,npgfa(iset)


             na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
             na2 = ipgf*ncoset(la_max(iset))

             zetp = zeta(ipgf,iset) 

             igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

             CALL integrate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0, 0.0_dbl,0,&
                        ra,rab,rab2,rs_v(igrid_level),&
                        cube_info(igrid_level),l_info,&
                        hab(na1:na2,1:1),pab(na1:na2,1:1),&
                        vratio,eps_gvg_rspace,calculate_forces,&
                        force_a,force_b)

           END DO
! in work are stored the eigenforces (i.e. forces on coefficients)
           CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))
 
           DO i=1,nsgfa(iset)
            eigenforce(i+offset) = work(i,1)
           ENDDO

           offset=offset+nsgfa(iset)

        END DO

!   *** Update forces ***

        IF (PRESENT (force_type)) THEN
          IF (force_type=='hartree') THEN
            force(ikind)%f_hartree(:,iatom) = force(ikind)%f_hartree(:,iatom) +  force_a(:)
          ELSE IF (force_type=='xc') THEN
            force(ikind)%f_xc(:,iatom) =  force(ikind)%f_xc(:,iatom) + force_a(:)
          END IF
        ELSE
          force(ikind)%f_rho(:,iatom) =&
               force(ikind)%f_rho(:,iatom) + force_a(:)
        END IF
      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_deallocate(rs_v(igrid_level))
    ENDDO

!   *** Release work storage ***

    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab")

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_vp_rspace_forces

! *****************************************************************************

    SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                            lb_max,zetb,lb_min,&
                                            ra,rab,rab2,rsgrid,&
                                            cube_info,l_info,hab,pab,&
                                            vratio,eps_gvg_rspace,&
                                            calculate_forces,force_a,force_b)

    TYPE(realspace_grid_type), TARGET, INTENT(IN) :: rsgrid
    REAL(dbl), INTENT(IN)                  :: eps_gvg_rspace,rab2,zeta,zetb,vratio
    INTEGER, INTENT(IN)                   :: la_max,la_min,lb_max,lb_min
    REAL(dbl), DIMENSION(3), INTENT(IN)    :: ra,rab
    TYPE(cube_info_type),INTENT(IN)       :: cube_info
    TYPE(l_info_type),INTENT(IN)          :: l_info
    REAL(dbl), DIMENSION(3), INTENT(INOUT) :: force_a,force_b
    REAL(dbl), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(INOUT) :: hab
    REAL(dbl), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab
    LOGICAL, INTENT(IN)                   :: calculate_forces

!   *** Local variables ***

    REAL(dbl) :: f,ftza,ftzb,pij,prefactor,radius,yz,z,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,handle,i,ico,ig,j,jco,jg,k,kg,la,lb,&
                lb_cube_min,ub_cube_max,gridbounds(2,3)

    REAL(dbl), DIMENSION(3) :: dhab,dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(dbl), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER coef_max,lxa,lxyz,l
    REAL(dbl) a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,rpg,s,pg
    REAL(dbl)  :: xa,xb,xap,xbp,radius2,maxpab
    REAL(dbl) :: dx2,dy2,dz2,dxi,dyi,dzi,cutoff
    INTEGER  :: kgmax,kgmin,jgmax,jgmin,igmax,igmin,j2,jg2,lar(3),lbr(3), &
                start,length,offset
    INTEGER  :: la_max_l,la_min_l,lb_min_l

    INTEGER, pointer, dimension(:) :: ly_max,lz_max
    REAL(dbl), pointer :: pzyx(:), alpha(:,:)
    integer, pointer :: ipzyx(:,:,:,:,:,:)
    REAL(dbl), pointer :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    REAL(dbl), pointer :: dpx(:,:)
    INTEGER, pointer, dimension(:)  :: sphere_bounds
    INTEGER :: cmax

!   ---------------------------------------------------------------------------

    IF (calculate_forces) THEN
      la_max_l=la_max+1  ! needed for the derivative of the gaussian, unimportant which one
      la_min_l=MAX(la_min-1,0) ! just in case the la_min,lb_min is not zero
      lb_min_l=MAX(lb_min-1,0)
    ELSE
      la_max_l=la_max
      la_min_l=la_min
      lb_min_l=lb_min
    END IF

    coef_max=la_max_l+lb_max+1
    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
!   *** position of the gaussian product
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)  ! this is the gaussian center in real coordinates
    rb(:) = ra(:) + rab(:)


! the cutoff could be vratio instead of 1.0_dbl if not rebuild
! in this case the minimization goes fine up to the point where the energy is
! noisy
    cutoff=1.0_dbl
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp,&
                                    zetp,eps_gvg_rspace,prefactor*vratio,cutoff)

    IF (radius == 0.0_dbl) RETURN
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    CALL return_l_info(l_info,la_min_l,la_max_l,lb_min_l,lb_max,lx_max,lxy_max, &
                         lxyz_max,ly_max,lz_max,&
                         map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** properties of the grid ***

    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),dbl)*dr(:)
    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    DO i=1,3
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
          map(ig,i) = ig+offset
         END DO
         if (start+length.ge.ub_cube(i)) EXIT
         start=start+length+1
        END DO
    ENDDO

!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,dbl)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max_l
       zbp=1.0_dbl
       DO lzb=0,lb_max
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max_l
      DO lxb=0,lb_max
       DO lya=0,la_max_l-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min_l-lxa-lya,0),la_max_l-lxa-lya
          DO lzb=max(lb_min_l-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,dbl)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max_l
       ybp=1.0_dbl
       DO lyb=0,lb_max
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max_l
      DO lxb=0,lb_max
       DO lya=0,la_max_l-lxa
       DO lyb=0,lb_max-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO


!   *** make the alpha matrix ***

    alpha(:,:)=0.0_dbl
    lx=0
    DO lxa=0,la_max_l
    DO lxb=0,lb_max
       lx=lx+1
       binomial_k_lxa=1.0_dbl
       a=1.0_dbl
       DO k=0,lxa
        binomial_l_lxb=1.0_dbl
        b=1.0_dbl
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*real(lxb-l,dbl)/real(l+1,dbl)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*real(lxa-k,dbl)/real(k+1,dbl)
        a=a*(rp(1)-ra(1))
       ENDDO
    ENDDO
    ENDDO


    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,dbl)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

!   *** do the loop over the grid

    call integrate_core(polx(1,-cmax),poly(1,-cmax),polz(1,-cmax),grid(1,1,1),alpha(1,1), &
                        lx_max,lxy_max,lxyz_max,coef_max,cmax, &
                        ly_max(1),lz_max(1),gridbounds(1,1),map(-cmax,1),pzyx(1), &
                        la_min_l,la_max_l,lb_min_l,lb_max,sphere_bounds(1))

    do i=1,lxyz_max
       pzyx(i)=pzyx(i)*prefactor
    enddo

!   *** pzyx contains all the information needed to find the elements of hab
!   *** and optionally of derivatives of these elements

    ftza = 2.0_dbl*zeta
    ftzb = 2.0_dbl*zetb

    DO la=la_min,la_max
      DO ax=0,la
        DO ay=0,la-ax
          az = la - ax - ay
          ico=coset(ax,ay,az)
          DO lb=lb_min,lb_max
            DO bx=0,lb
              DO by=0,lb-bx
                bz = lb - bx - by
                jco=coset(bx,by,bz)
                hab(ico,jco) = hab(ico,jco) + pzyx(ipzyx(ax,ay,az,bx,by,bz))
                IF (calculate_forces) THEN
                  dhab(1) = ftza*pzyx(ipzyx(ax+1,ay,az,bx,by,bz)) -&
                            REAL(ax,dbl)*pzyx(ipzyx(MAX(0,ax-1),ay,az,bx,by,bz))
                  dhab(2) = ftza*pzyx(ipzyx(ax,ay+1,az,bx,by,bz)) -&
                            REAL(ay,dbl)*pzyx(ipzyx(ax,MAX(0,ay-1),az,bx,by,bz))
                  dhab(3) = ftza*pzyx(ipzyx(ax,ay,az+1,bx,by,bz)) -&
                            REAL(az,dbl)*pzyx(ipzyx(ax,ay,MAX(0,az-1),bx,by,bz))
                  force_a(:) = force_a(:) + pab(ico,jco)*dhab(:)
                  dhab(1) = ftzb*(pzyx(ipzyx(ax+1,ay,az,bx,by,bz)) -&
                                  rab(1)*pzyx(ipzyx(ax,ay,az,bx,by,bz))) -&
                            REAL(bx,dbl)*pzyx(ipzyx(ax,ay,az,MAX(0,bx-1),by,bz))
                  dhab(2) = ftzb*(pzyx(ipzyx(ax,ay+1,az,bx,by,bz)) -&
                                  rab(2)*pzyx(ipzyx(ax,ay,az,bx,by,bz))) -&
                            REAL(by,dbl)*pzyx(ipzyx(ax,ay,az,bx,MAX(0,by-1),bz))
                  dhab(3) = ftzb*(pzyx(ipzyx(ax,ay,az+1,bx,by,bz)) -&
                                  rab(3)*pzyx(ipzyx(ax,ay,az,bx,by,bz))) -&
                            REAL(bz,dbl)*pzyx(ipzyx(ax,ay,az,bx,by,MAX(0,bz-1)))
                  force_b(:) = force_b(:) + pab(ico,jco)*dhab(:)
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE integrate_pgf_product_rspace

! *****************************************************************************

END MODULE kg_density
