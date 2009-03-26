!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the SPC plane wave density by collocating the primitive Gaussian
!>      functions (pgf). Also computes the integral of the potential with a density
!>      in the primitive gaussian basis
! *****************************************************************************
MODULE scp_density_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sync
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_fft_wrap,&
                                             pw_integrate_function,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace,&
                                             potential_pw2rs
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs2pw,&
                                             rs_grid_create,&
                                             rs_grid_release,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE scp_coeff_types,                 ONLY: aux_coeff_distributed,&
                                             aux_coeff_set_type,&
                                             aux_coeff_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_force_types,                 ONLY: scp_force_type
  USE scp_rho_types,                   ONLY: scp_rho_get,&
                                             scp_rho_type
  USE scp_rspw_types,                  ONLY: scp_rspw_get,&
                                             scp_rspw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_density_methods'
 
! *** Public subroutines ***

  PUBLIC :: collocate_scp_density, integrate_rhoscp_vrspace, update_rhoscp

CONTAINS

! *****************************************************************************
  SUBROUTINE update_rhoscp ( qs_env,  error )
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_rhoscp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: t1
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      POINTER                                :: rs_descs( : )
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(scp_rspw_type), POINTER             :: rspw

    CALL timeset ( routineN,handle )
    NULLIFY(atomic_kind_set,gridlevel_info,scp_env)
    NULLIFY(pw_pools,cube_info)
    NULLIFY(cube_info, rs_descs, pw_env, rspw, local_particles)
    NULLIFY(particle_set, cell, dft_control)

    stat = 0
    t1 = m_walltime()

    CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error=error )
    CALL get_scp_env ( scp_env = scp_env, rspw = rspw, error = error )
    CALL scp_rspw_get ( scp_rspw=rspw, qs_pw_env = pw_env, error = error ) 
    para_env=>qs_env%para_env
    ionode=para_env%ionode

    CALL get_qs_env ( qs_env = qs_env, &
                    atomic_kind_set=atomic_kind_set,    &
                    local_particles=local_particles, &
                    particle_set = particle_set, cell = cell,  &
                    dft_control = dft_control, error = error )

    CALL pw_env_get ( pw_env = pw_env, gridlevel_info=gridlevel_info, &
                      rs_descs = rs_descs, pw_pools=pw_pools, cube_info = cube_info, &
                      error = error )

! Compute SCP density
    CALL collocate_scp_density ( scp_env, atomic_kind_set, particle_set,        &
                                 local_particles, cell, rs_descs, pw_pools,     &
                                 cube_info, gridlevel_info,                     &
                                 dft_control % qs_control % eps_rho_rspace,     &
                                 error=error )
  t1 = m_walltime()
  CALL timestop(handle)
  END SUBROUTINE update_rhoscp 

! *****************************************************************************
  SUBROUTINE collocate_scp_density( scp_env, atomic_kind_set, particle_set, &
                                    local_particles, cell, rs_descs, pw_pools, &
                                    cube_info, gridlevel_info, eps_rho_rspace, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    REAL(dp), INTENT(IN), OPTIONAL           :: eps_rho_rspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocate_scp_density', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, igrid_level, ikind, iparticle_local, ipgf, &
      iset, istat, maxco, na1, na2, natom, ncoa, nkind, nparticle_local, &
      nseta, offset, sgfa, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    REAL(KIND=dp)                            :: dab, my_eps_rho_rspace, rab2, &
                                                scale, total_rho, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_p_type), POINTER                 :: rhop_g, rhop_r
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(scp_rho_type), POINTER              :: rho_scp

    CALL timeset(routineN,handle)

    NULLIFY(aux_basis_set,aux_coeff_set,&
           local_coeffs,rs_rho, mgrid_rspace, &
           mgrid_gspace,rhop_r, rhop_g)
    NULLIFY(sphi_a,zeta,pab,work)
    NULLIFY(la_max,la_min,npgfa,nsgfa,first_sgfa)

    IF ( PRESENT ( eps_rho_rspace ) ) THEN
      my_eps_rho_rspace = eps_rho_rspace
    ELSE 
      CALL get_scp_env(scp_env=scp_env, eps_rho_rspace = my_eps_rho_rspace, error= error )
    ENDIF 

    CALL get_scp_env( scp_env=scp_env, rho_scp = rho_scp, &
                      aux_coeff_set=aux_coeff_set, error=error)
    CALL scp_rho_get ( rho_scp = rho_scp, rhop_r = rhop_r, &
                       rhop_g = rhop_g, error=error )
    
    ! *** set up the pw multi-grids *** !

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)
    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)

    ! *** set up the rs multi-grids
    ALLOCATE (rs_rho(gridlevel_info%ngrid_levels),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"rs_rho")
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_create(rs_rho(igrid_level)%rs_grid, rs_descs(igrid_level)%rs_desc, error=error)
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    nkind=SIZE(atomic_kind_set)
    IF(aux_coeff_set%distribution_method/=aux_coeff_distributed) THEN
      CALL stop_program ( 'collocate_scp_density','replicated coefs not yet implemented')
    END IF

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom )
                             
    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"pab",maxco*1*dp_size)
    pab = 0.0_dp

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",maxco*1*dp_size)
    work = 0.0_dp

    DO ikind=1,nkind

      atomic_kind=> atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           aux_basis_set=aux_basis_set)

      IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE

      local_coeffs=>aux_coeff_set%coeffs_of_kind(ikind)%coeffs
      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
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
         rb(:) = 0.0_dp
         rab(:) = 0.0_dp
         rab2  = 0.0_dp
         dab   = 0.0_dp
         offset = 0

         DO iset=1,nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1,iset)

! polarization density is not frozen. The coefficients of each contracted
! gaussian polarization functions , stored in eigenvector(:), are here 
! multiplied for the contraction coefficients in sphi and stored in work(:,:). 
! Then a decontraction is performed in order to operate on each pgf.

            DO i=1,nsgfa(iset)
              work(i,1)=local_coeffs%c(iparticle_local,offset+i)
            ENDDO

            CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dp,pab(1,1),SIZE(pab,1))

            DO ipgf=1,npgfa(iset)

               na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
               na2 = ipgf*ncoset(la_max(iset))

               scale = 1.0_dp  ! SCP density is electronic
               zetp = zeta(ipgf,iset)
               igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

               CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_dp,0,&
                        ra,rab,rab2,scale,pab,na1-1,0,&
                        rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                        eps_rho_rspace,ga_gb_function=401,error=error)
            END DO

            offset=offset+nsgfa(iset)

         END DO
      END DO
    END DO

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

    IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rhop_r%pw,rs2pw,error=error)
       CALL rs_grid_release(rs_rho(1)%rs_grid,error=error)
       DEALLOCATE (rs_rho,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"rs_rho")
       CALL pw_fft_wrap(rhop_r%pw,rhop_g%pw)
       IF (rhop_r%pw%pw_grid%spherical) THEN ! rhop_g = rhop_r
          CALL pw_fft_wrap(rhop_g%pw,rhop_r%pw)
       ENDIF
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw,error=error)
          CALL rs_grid_release(rs_rho(igrid_level)%rs_grid,error=error)
       ENDDO
       DEALLOCATE (rs_rho,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"rs_rho")

       ! we want both rho an rho_gspace, the latter for Hartree and co-workers.
       CALL pw_zero(rhop_g%pw)
       DO igrid_level=1,gridlevel_info%ngrid_levels
         CALL pw_fft_wrap(mgrid_rspace(igrid_level)%pw,&
              mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw,rhop_g%pw)
       END DO
       CALL pw_fft_wrap(rhop_g%pw,rhop_r%pw)
    END IF

    total_rho = pw_integrate_function(rhop_r%pw)

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,&
                                   error=error)
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,&
                                   error=error)

    CALL timestop(handle)

   END SUBROUTINE collocate_scp_density

! *****************************************************************************
  SUBROUTINE integrate_rhoscp_vrspace ( scp_env, v_rspace, pw_env, atomic_kind_set, &
                                        particle_set, local_particles, cell,  &
                                        eps_gvg_rspace, calculate_forces, just_energy, debug, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(pw_p_type)                          :: v_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN), OPTIONAL           :: eps_gvg_rspace
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces, &
                                                just_energy, debug
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_rhoscp_vrspace', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, auxbas_grid, handle, i, iatom, igrid_level, ikind, &
      iparticle_local, ipgf, iset, istat, maxco, maxsgf_set, na1, na2, natom, &
      natom_of_kind, ncoa, nkind, nparticle_local, nseta, offset, sgfa, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: energy_only, &
                                                my_calculate_forces, my_debug
    REAL(KIND=dp)                            :: dab, my_eps_gvg_rspace, rab2, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, rb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab, rpgfa, sphi_a, &
                                                work, zeta
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(scp_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(section_vals_type), POINTER         :: dft_section, interp_section

    CALL timeset(routineN,handle)
    NULLIFY(aux_basis_set,aux_coeff_set)
    NULLIFY(local_coeffs)
    NULLIFY(rs_v, rs_descs, pw_pools)
    NULLIFY(sphi_a,rpgfa,zeta,hab,pab,work)
    NULLIFY(la_max,la_min,npgfa,nsgfa,first_sgfa)
    NULLIFY(para_env, interp_section )

    my_calculate_forces = .FALSE.
    IF ( PRESENT ( calculate_forces ) ) my_calculate_forces = calculate_forces
    my_debug= .FALSE.
    IF ( PRESENT ( debug ) ) my_debug = debug
    energy_only = .FALSE.
    IF ( PRESENT ( just_energy ) ) energy_only = just_energy
    CALL get_scp_env(scp_env=scp_env,&
                     aux_coeff_set=aux_coeff_set,&
                     para_env=para_env,  &
                     input = dft_section, &
                     error=error)

    IF ( my_calculate_forces ) CALL get_scp_env ( scp_env, force = force, error = error )

    CALL pw_env_get ( pw_env = pw_env, auxbas_grid=auxbas_grid, &
                      rs_descs = rs_descs, pw_pools=pw_pools, &
                      cube_info = cube_info, gridlevel_info=gridlevel_info, &
                      error = error ) 

    IF ( PRESENT ( eps_gvg_rspace ) ) THEN
      my_eps_gvg_rspace = eps_gvg_rspace 
    ELSE
      CALL get_scp_env(scp_env=scp_env,&
                       eps_gvg_rspace = my_eps_gvg_rspace, &
                       error=error)
    ENDIF

    CALL mp_sync ( para_env%group )
    ! *** set up of the potential on the multigrids using the pw_env
    interp_section => section_vals_get_subs_vals(dft_section,"MGRID%INTERPOLATOR",&
         error=error)
    ALLOCATE (rs_v(SIZE(rs_descs)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"rs_v")
    DO i=1, SIZE(rs_v)
      CALL rs_grid_create(rs_v(i)%rs_grid, rs_descs(i)%rs_desc, error=error) 
    END DO
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

!   *** having the potential on the rs_multigrids, just integrate ...

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"atom_of_kind",natom*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (hab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"hab",maxco*dp_size)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"pab",maxco*dp_size)
    pab=0.0_dp

    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",maxco*dp_size)
    work=0.0_dp

    offset = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
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

      local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local = 1, nparticle_local

        iatom = local_particles%list(ikind)%array(iparticle_local)
        ra(:)  = pbc(particle_set(iatom)%r,cell)

        force_a(:) = 0.0_dp
        force_b(:) = 0.0_dp

        rb(:) = 0.0_dp
        rab(:) = 0.0_dp
        rab2  = 0.0_dp
        dab   = 0.0_dp
        offset = 0

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          DO i=1,nsgfa(iset)
             work(i,1)=local_coeffs%c(iparticle_local,offset+i)
          ENDDO

          CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dp,pab(1,1),SIZE(pab,1))

           hab(:,:) = 0.0_dp

           DO ipgf=1,npgfa(iset)

             na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
             na2 = ipgf*ncoset(la_max(iset))

             zetp = zeta(ipgf,iset)

             igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

             CALL integrate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0, 0.0_dp,0,&
                        ra,rab,rab2,rs_v(igrid_level)%rs_grid,&
                        cell,cube_info(igrid_level),&
                        hab,pab,na1-1,0,&
                        eps_gvg_rspace=eps_gvg_rspace,&
                        calculate_forces=calculate_forces,  &
                        force_a = force_a, force_b = force_b, &
                        error=error)

           END DO
           ! in work are stored the eigenforces (i.e. forces on coefficients)
           IF ( .NOT. energy_only ) THEN
             CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                            1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                            hab(1,1),SIZE(hab,1),&
                            0.0_dp,work(1,1),SIZE(work,1))

             DO i=1,nsgfa(iset)
              local_coeffs%fc(iparticle_local,offset+i) = &
                    local_coeffs%fc(iparticle_local,offset+i) -  work(i,1)
              IF ( my_debug ) WRITE ( *, * ) 'F_ANALYTICAL', -work ( i, 1 )
             ENDDO
           END IF ! energy only

           offset=offset+nsgfa(iset)

        END DO
        !   *** Update forces ***

        IF(my_calculate_forces) THEN
          atom_a=atom_of_kind(iatom)
          force(ikind)%f_scp (:,atom_a) = force(ikind)%f_scp(:,atom_a) - force_a(:)
        END IF

      END DO

    END DO

    IF (ASSOCIATED(rs_v)) THEN
      DO i=1, SIZE(rs_v)
        CALL rs_grid_release(rs_v(i)%rs_grid, error=error)
      END DO
      DEALLOCATE (rs_v,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"rs_v") 
    END IF

    !   *** Release work storage ***

    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"hab")

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work")

    CALL timestop(handle)

  END SUBROUTINE integrate_rhoscp_vrspace

END MODULE scp_density_methods
