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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cube_utils,                      ONLY: cube_info_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE kinds,                           ONLY: dbl,&
                                             dp_size,&
                                             int_size
  USE l_utils,                         ONLY: l_info_type
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_allocate,&
                                             rs_grid_deallocate,&
                                             rs_grid_setup,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: calculate_density, calculate_epc_density, &
            calculate_epc_rspace_forces, calculate_p_density, &
            calculate_v_rspace_forces, calculate_vp_rspace_forces,&
            calculate_density_single, calculate_epc_density_single, &
            calculate_vp_rspace_forces_s

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
    TYPE(coeff_type), INTENT(INOUT)          :: rho, rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(dbl), INTENT(OUT)                   :: total_rho
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_density (MODULE kg_density)"

    INTEGER :: handle, i, iatom, igrid_level, ikind, iparticle_local, ipgf, &
      iset, istat, maxco, maxsgf_set, na1, na2, ncoa, nkind, nparticle_local, &
      nseta, nsmax, sgfa
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    REAL(dbl)                                :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(dbl), DIMENSION(3)                  :: ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho

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
    nsmax=-1

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    local_particles=local_particles ) 

    nkind=size(atomic_kind_set)
    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*dp_size)

    DO ikind=1,nkind

      atomic_kind=> atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local=1,nparticle_local
         iatom = local_particles%list(ikind)%array(iparticle_local)
         ra(:) = pbc(particle_set(iatom)%r,cell)
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
                           ra,rab,rab2,scale,pab,na1-1,0,&
                           rs_rho(igrid_level),cube_info(igrid_level),&
                           l_info,eps_rho_rspace)
            END DO

         END DO                     

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_pw_transfer(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw,"FORWARD")
       call rs_grid_deallocate(rs_rho(igrid_level))
    ENDDO

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),mgrid_gspace(igrid_level))
      CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
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

  SUBROUTINE calculate_density_single(rho,rho_gspace,&
                                    mgrid_rspace,mgrid_gspace,&
                                    gridlevel_info,cube_info,l_info,&
                                    kg_env,akind, globenv)

!   REAL(dbl), DIMENSION(:), POINTER             :: eigenvector   
    TYPE(coeff_type), INTENT(INOUT)          :: rho, rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    INTEGER, INTENT(IN)                      :: akind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_density_single (MODULE kg_density)"

    INTEGER :: handle, i, iatom, igrid_level, ikind, ipgf, iset, istat, &
      maxco, maxsgf_set, na1, na2, natom, ncoa, nseta, nsmax, offset, sgfa
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    REAL(dbl)                                :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(dbl), DIMENSION(3)                  :: ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho

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
    nsmax=-1

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set ) 

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*dp_size)

    offset=0

    DO ikind=1,SIZE(atomic_kind_set)

      if (ikind/=akind) CYCLE

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               orb_basis_set=orb_basis_set)

      iatom=atom_list(1)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      ra(:) = pbc(particle_set(iatom)%r,cell)
      rb(:) = 0.0_dbl
      rab(:) = 0.0_dbl
      rab2  = 0.0_dbl
      dab   = 0.0_dbl

      DO iset=1,1

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
                        ra,rab,rab2,scale,pab,na1-1,0,&
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
      CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    deallocate(rs_rho)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_density_single


! *****************************************************************************

  FUNCTION calculate_total_density(rho) RESULT(total_rho)

    TYPE(coeff_type), INTENT(IN), TARGET     :: rho
    REAL(dbl)                                :: total_rho

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "FUNCTION calculate_total_density (MODULE kg_density)"

    INTEGER, DIMENSION(:), POINTER           :: lb_grid

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


!   *** Arguments ***

    TYPE(coeff_type), INTENT(OUT)            :: rho_core
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(dbl), INTENT(OUT)                   :: total_rho
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: handle, iatom, ikind, &
                                                iparticle_local, &
                                                nparticle_local, nsmax
    REAL(dbl)                                :: alpha, eps_rho_rspace, scale
    REAL(dbl), DIMENSION(3)                  :: ra
    REAL(dbl), DIMENSION(:, :), POINTER      :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rhoc_r
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(realspace_grid_type)                :: rs_rho

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_epc_density","I","",handle)

    ALLOCATE(pab(1,1))

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    local_particles=local_particles,&
                    particle_set=particle_set)

    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax=-1   ! should be temporary fix
    call rs_grid_setup(rs_rho,rho_core%pw%pw_grid,nsmax)
    call rs_grid_allocate(rs_rho)
    call rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace
    scale = -1.0_dbl

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local=1,nparticle_local
        iatom = local_particles%list(ikind)%array(iparticle_local)
        ra(:) = pbc(particle_set(iatom)%r,cell)
        CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dbl,0,ra,&
                                          (/0.0_dbl,0.0_dbl,0.0_dbl/),&
                                          0.0_dbl,scale,pab,0,0,rs_rho,&
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

    DEALLOCATE(pab)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_epc_density

! *****************************************************************************

  SUBROUTINE calculate_epc_density_single (rho_core, cube_info, l_info, akind, &
                                    kg_env, globenv)
! *****************************************************************************


!   *** Arguments ***

    TYPE(coeff_type), INTENT(OUT)            :: rho_core
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    INTEGER, INTENT(IN)                      :: akind
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: atom_a, handle, ikind, nsmax
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    REAL(dbl)                                :: alpha, eps_rho_rspace, scale
    REAL(dbl), DIMENSION(3)                  :: ra
    REAL(dbl), DIMENSION(:, :), POINTER      :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rhoc_r
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(realspace_grid_type)                :: rs_rho

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_epc_density","I","",handle)

    ALLOCATE(pab(1,1))

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set)

    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax=-1   ! should be temporary fix
    call rs_grid_setup(rs_rho,rho_core%pw%pw_grid,nsmax)
    call rs_grid_allocate(rs_rho)
    call rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace
    scale = -1.0_dbl

    DO ikind=1,SIZE(atomic_kind_set)

      if (ikind /= akind) CYCLE
      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      atom_a = atom_list(1)
      ra(:) = pbc(particle_set(atom_a)%r,cell)
      IF(rs_rho%group_head) &
      CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dbl,0,ra,&
                                          (/0.0_dbl,0.0_dbl,0.0_dbl/),&
                                          0.0_dbl,scale,pab,0,0,rs_rho,&
                                          cube_info,l_info,eps_rho_rspace)


    END DO

    CALL coeff_allocate(rhoc_r,rho_core%pw%pw_grid,REALDATA3D)
    rhoc_r%pw%in_space = REALSPACE

    call rs_pw_transfer(rs_rho,rhoc_r%pw,"FORWARD")
    call rs_grid_deallocate(rs_rho)

    CALL coeff_transform_space(rhoc_r,rho_core)

    CALL coeff_deallocate(rhoc_r)

    DEALLOCATE(pab)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_epc_density_single

! *****************************************************************************

  SUBROUTINE calculate_p_density(eigenvector,rho,rho_gspace,&
                                    mgrid_rspace,mgrid_gspace,&
                                    gridlevel_info,cube_info,l_info,&
                                    total_rho, kg_env,globenv)

    REAL(dbl), DIMENSION(:), POINTER         :: eigenvector
    TYPE(coeff_type), INTENT(INOUT)          :: rho, rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(dbl), INTENT(OUT)                   :: total_rho
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_p_density (MODULE kg_density)"

    INTEGER :: handle, i, iatom, igrid_level, ipgf, iset, istat, maxco, &
      maxsgf_set, na1, na2, natom, ncoa, nseta, nsmax, offset, sgfa
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    REAL(dbl)                                :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(dbl), DIMENSION(3)                  :: ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho

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
    nsmax=-1

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set ) 

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*dp_size)

    offset=0

    DO iatom=1,natom

      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               aux_basis_set=aux_basis_set)

      IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      ra(:) = pbc(particle_set(iatom)%r,cell)
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
                        ra,rab,rab2,scale,pab,na1-1,0,&
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
      CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
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

  SUBROUTINE calculate_epc_rspace_forces(v_rspace,cube_info,l_info,kg_env,force)
   

!   *** Arguments ***

    TYPE(coeff_type), INTENT(INOUT)          :: v_rspace
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(INOUT) :: kg_env
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_epc_rspace_forces (MODULE kg_density)"

    INTEGER                                  :: atom_a, handle, iatom, &
                                                ikind, iparticle_local, &
                                                istat, natom, nkind, &
                                                nparticle_local, nsmax
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    REAL(dbl)                                :: alpha_core_charge, &
                                                ccore_charge, eps_rho_rspace
    REAL(dbl), DIMENSION(3)                  :: force_a, force_b, ra
    REAL(dbl), DIMENSION(:, :), POINTER      :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(realspace_grid_type), POINTER       :: rs_pot
    TYPE(realspace_grid_type), target        :: rs_v

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_epc_rspace_forces","I","",handle)

    ALLOCATE(hab(1,1),pab(1,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab,pab",1*int_size)

    nsmax = -1
    CALL rs_grid_setup(rs_v,v_rspace%pw%pw_grid,nsmax)
    CALL rs_grid_allocate(rs_v)
    CALL rs_pw_transfer(rs_v,v_rspace%pw,"BACKWARD")

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    local_particles=local_particles,&
                    particle_set=particle_set) 

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_rho_rspace = dft_control%kg_control%eps_rho_rspace

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge)

      pab(1,1) = ccore_charge
      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local=1,nparticle_local

        iatom=local_particles%list(ikind)%array(iparticle_local)
        ra(:)  = pbc(particle_set(iatom)%r,cell)

        hab(1,1) = 0.0_dbl
        force_a(:) = 0.0_dbl
        force_b(:) = 0.0_dbl

        rs_pot=> rs_v
        CALL integrate_pgf_product_rspace(0,alpha_core_charge,0,&
                                          0,0.0_dbl,0,&
                                          ra,(/0.0_dbl,0.0_dbl,0.0_dbl/),0.0_dbl,&
                                          rs_pot,cube_info,l_info,&
                                          hab,pab,0,0,1.0_dbl,&
                                          eps_rho_rspace,&
                                          .TRUE.,force_a,force_b)

        atom_a=atom_of_kind(iatom)
        force(ikind)%f_hartree_core(:,atom_a) =&
          force(ikind)%f_hartree_core(:,atom_a) + force_a(:)

      END DO

    END DO

    CALL rs_grid_deallocate(rs_v)

    DEALLOCATE(pab,hab,atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab,pab,atom_of_kind")

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_epc_rspace_forces

! *****************************************************************************

  SUBROUTINE calculate_v_rspace_forces(v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,&
                                kg_env,force,globenv,force_type)

!   ***  Arguments   ***

    TYPE(coeff_type)                         :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    INTEGER, INTENT(IN)                      :: auxbas_grid
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: force_type

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_v_rspace_forces (MODULE kg_density)"

    INTEGER :: atom_a, handle, i, iatom, igrid_level, ikind, iparticle_local, &
      ipgf, iset, istat, maxco, maxsgf_set, na1, na2, natom, ncoa, nkind, &
      nparticle_local, nseta, nsmax, sgfa
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: calculate_forces
    REAL(dbl)                                :: dab, eps_gvg_rspace, rab2, &
                                                scale, vratio, zetp
    REAL(dbl), DIMENSION(3)                  :: force_a, force_b, ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: hab, pab, rpgfa, sphi_a, &
                                                work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(REALSPACE_GRID_TYPE), POINTER       :: rs_pot

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

    nsmax=-1   ! should be temporary fix

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_setup(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, &
            nsmax)
       CALL rs_grid_allocate(rs_v(igrid_level))
       CALL rs_pw_transfer(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw,&
                           "BACKWARD")
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    local_particles=local_particles,&
                    particle_set=particle_set )

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_gvg_rspace = dft_control%kg_control%eps_gvg_rspace

!   *** Allocate work storage ***

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*dp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*dp_size)

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*dp_size)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
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

      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local = 1, nparticle_local

        iatom = local_particles%list(ikind)%array(iparticle_local)
        ra(:)  = pbc(particle_set(iatom)%r,cell)

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

             rs_pot=>rs_v(igrid_level)
             CALL integrate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0, 0.0_dbl,0,&
                        ra,rab,rab2,rs_pot,&
                        cube_info(igrid_level),l_info,&
                        hab,pab=pab,o1=na1-1,o2=0,&
                        vratio=vratio,eps_gvg_rspace=eps_gvg_rspace,&
                        calculate_forces=calculate_forces,&
                        force_a=force_a,force_b=force_b)

           END DO

        END DO

!   *** Update forces ***

        atom_a=atom_of_kind(iatom)
        IF (PRESENT (force_type)) THEN
          IF (force_type=='hartree') THEN
            force(ikind)%f_hartree(:,atom_a) =  force_a(:)
          ELSE IF (force_type=='xc') THEN
            force(ikind)%f_xc(:,atom_a) = force_a(:)
          END IF
        ELSE
          force(ikind)%f_rho(:,atom_a) =&
               force(ikind)%f_rho(:,atom_a) + 1.0_dbl*force_a(:)
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

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_v_rspace_forces

! *****************************************************************************

  SUBROUTINE calculate_vp_rspace_forces(eigenvector,eigenforce, &
                                v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,kg_env,&
                                force,globenv,force_type)

!   ***  Arguments   ***

    REAL(dbl), DIMENSION(:), POINTER         :: eigenvector, eigenforce
    TYPE(coeff_type)                         :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    INTEGER, INTENT(IN)                      :: auxbas_grid
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: force_type

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_vp_rspace_forces (MODULE kg_density)"

    INTEGER :: atom_a, handle, i, iatom, igrid_level, ikind, ipgf, iset, &
      istat, maxco, maxsgf_set, na1, na2, natom, natom_of_kind, ncoa, nkind, &
      nseta, nsmax, offset, sgfa
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: calculate_forces
    REAL(dbl)                                :: dab, eps_gvg_rspace, rab2, &
                                                scale, vratio, zetp
    REAL(dbl), DIMENSION(3)                  :: force_a, force_b, ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: hab, pab, rpgfa, sphi_a, &
                                                work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(REALSPACE_GRID_TYPE), POINTER       :: rs_pot

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

    nsmax=-1   ! should be temporary fix

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_setup(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, &
            nsmax)
       CALL rs_grid_allocate(rs_v(igrid_level))
       CALL rs_pw_transfer(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw,&
                           "BACKWARD")
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set )

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_gvg_rspace = dft_control%kg_control%eps_gvg_rspace

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*dp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*dp_size)

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*dp_size)
    
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

             rs_pot=>rs_v(igrid_level)
             CALL integrate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0, 0.0_dbl,0,&
                        ra,rab,rab2,rs_pot,&             
                        cube_info(igrid_level),l_info,&
                        hab,pab,na1-1,0,&
                        vratio,eps_gvg_rspace,calculate_forces,&
                        force_a,force_b)

           END DO
! in work are stored the eigenforces (i.e. forces on coefficients)
           CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))
 
           DO i=1,nsgfa(iset)
            eigenforce(i+offset) = eigenforce(i+offset)+ work(i,1)
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

  SUBROUTINE calculate_vp_rspace_forces_s(eigenforce, &
                                v_rspace,mgrid_rspace,mgrid_gspace,&
                                auxbas_grid,gridlevel_info,&
                                cube_info,l_info,kg_env,&
                                akind,globenv)

!   ***  Arguments   ***

    REAL(dbl), DIMENSION(:, :), &
      intent(inout)                          :: eigenforce
    TYPE(coeff_type)                         :: v_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_rspace, mgrid_gspace
    INTEGER, INTENT(IN)                      :: auxbas_grid
    TYPE(gridlevel_info_type), INTENT(IN)    :: gridlevel_info
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    INTEGER, INTENT(IN)                      :: akind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_vp_rspace_forces (MODULE kg_density)"

    INTEGER :: atom_a, handle, igrid_level, ikind, ipgf, iset, isgf, &
      istat, maxco, maxsgf_set, na1, na2, natom, ncoa, nkind, nseta, nsmax, &
      offset, sgfa
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: calculate_forces
    REAL(dbl)                                :: dab, eps_gvg_rspace, rab2, &
                                                scale, vratio, zetp
    REAL(dbl), DIMENSION(3)                  :: force_a, force_b, ra, rab, rb
    REAL(dbl), DIMENSION(:, :), POINTER      :: hab, pab, rpgfa, sphi_a, &
                                                work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(REALSPACE_GRID_TYPE), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(REALSPACE_GRID_TYPE), POINTER       :: rs_pot

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_v_rspace_forces","I","",handle)

    vratio = 1.0_dbl
    calculate_forces =.FALSE.

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

    nsmax=-1   ! should be temporary fix

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_setup(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, &
            nsmax)
       CALL rs_grid_allocate(rs_v(igrid_level))
       CALL rs_pw_transfer(rs_v(igrid_level),mgrid_rspace(igrid_level)%pw,&
                           "BACKWARD")
    ENDDO

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set )

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    eps_gvg_rspace = dft_control%kg_control%eps_gvg_rspace

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*dp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*dp_size)

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*dp_size)
    
    offset = 0

    DO ikind=1,nkind

      if (ikind /= akind) cycle
      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
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

      atom_a = atom_list(1)
      ra(:)  = pbc(particle_set(atom_a)%r,cell)

      force_a(:) = 0.0_dbl
      force_b(:) = 0.0_dbl
      rb(:) = 0.0_dbl
      rab(:) = 0.0_dbl
      rab2  = 0.0_dbl
      dab   = 0.0_dbl

      DO iset=1,1

        ncoa = npgfa(iset)*ncoset(la_max(iset))
        sgfa = first_sgfa(1,iset)

         hab(:,:) = 0.0_dbl
         pab(:,:) = 0.0_dbl

         DO ipgf=1,npgfa(iset)


           na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
           na2 = ipgf*ncoset(la_max(iset))

           zetp = zeta(ipgf,iset) 

           igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

           rs_pot=>rs_v(igrid_level)
           CALL integrate_pgf_product_rspace(&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      0, 0.0_dbl,0,&
                      ra,rab,rab2,rs_pot,&             
                      cube_info(igrid_level),l_info,&
                      hab,pab,na1-1,0,&
                      vratio,eps_gvg_rspace,calculate_forces,&
                      force_a, force_b)  

         END DO
! in work are stored the eigenforces (i.e. forces on coefficients)
         CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                        1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                        hab(1,1),SIZE(hab,1),&
                        0.0_dbl,work(1,1),SIZE(work,1))
 
         DO isgf=1,nsgfa(iset)
          eigenforce(akind,isgf) = eigenforce(akind,isgf)+ work(isgf,1)
         ENDDO

         offset=offset+nsgfa(iset)

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

  END SUBROUTINE calculate_vp_rspace_forces_s

! *****************************************************************************

END MODULE kg_density
! *****************************************************************************
