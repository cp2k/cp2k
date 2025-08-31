!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CNEO soft nuclear densities on the global grid
!>      (see J. Chem. Theory Comput. 2025, 21, 16, 7865–7877)
!> \par History
!>      08.2025 created [zc62]
!> \author Zehua Chen
! **************************************************************************************************
MODULE qs_cneo_ggrid
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                              gridlevel_info_type
   USE grid_api,                        ONLY: GRID_FUNC_AB,&
                                              collocate_pgf_product
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type,&
                                              pw_pools_create_pws,&
                                              pw_pools_give_back_pws
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_cneo_types,                   ONLY: cneo_potential_type,&
                                              get_cneo_potential,&
                                              rhoz_cneo_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              rho0_mpole_type
   USE realspace_grid_types,            ONLY: map_gaussian_here,&
                                              realspace_grid_type,&
                                              rs_grid_zero,&
                                              transfer_rs2pw
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cneo_ggrid'

   PUBLIC :: put_rhoz_cneo_s_on_grid, rhoz_cneo_s_grid_create, integrate_vhgg_rspace

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0 ...
!> \param rhoz_cneo_set ...
!> \param tot_rs_int ...
! **************************************************************************************************
   SUBROUTINE put_rhoz_cneo_s_on_grid(qs_env, rho0, rhoz_cneo_set, tot_rs_int)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho0_mpole_type), POINTER                     :: rho0
      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      REAL(KIND=dp), INTENT(OUT)                         :: tot_rs_int

      CHARACTER(LEN=*), PARAMETER :: routineN = 'put_rhoz_cneo_s_on_grid'

      INTEGER :: atom_a, group_size, handle, iatom, igrid_level, ikind, ipgf, iset, jpgf, jset, &
         m1, maxco, maxsgf_set, my_pos, na1, natom, nb1, ncoa, ncob, npgf2, nseta, offset, sgfa, &
         sgfb
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: paw_atom
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: map_it2
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius, scale, zeff, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, pab, sphi_a, work, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cneo_potential_type), POINTER                 :: cneo_potential
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: nuc_soft_basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_gspace
      TYPE(pw_c1d_gs_type), POINTER                      :: rhoz_cneo_s_gs
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_r3d_rs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_rspace
      TYPE(pw_r3d_rs_type), POINTER                      :: rhoz_cneo_s_rs
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), DIMENSION(:), POINTER   :: rs_rho
      TYPE(realspace_grid_type), POINTER                 :: rs_grid

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, qs_kind_set, atom_list, cell, dft_control, &
               first_sgfa, gridlevel_info, la_max, la_min, nuc_soft_basis, &
               npgfa, nsgfa, p_block, pab, particle_set, pw_env, pw_pools, &
               rs_grid, rs_rho, sphi_a, work, zeta)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, particle_set=particle_set, &
                      pw_env=pw_env, &
                      dft_control=dft_control)

      ! find maximum numbers
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type="NUC_SOFT")
      IF (maxco > 0) THEN
         ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set))

         eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

         ! *** set up the pw multi-grids *** !
         CPASSERT(ASSOCIATED(pw_env))
         CALL pw_env_get(pw_env=pw_env, rs_grids=rs_rho, pw_pools=pw_pools, &
                         gridlevel_info=gridlevel_info)

         CALL pw_pools_create_pws(pw_pools, mgrid_rspace)

         CALL pw_pools_create_pws(pw_pools, mgrid_gspace)

         ! *** set up the rs multi-grids *** !
         DO igrid_level = 1, gridlevel_info%ngrid_levels
            CALL rs_grid_zero(rs_rho(igrid_level))
         END DO

         offset = 0
         my_pos = mgrid_rspace(1)%pw_grid%para%group%mepos
         group_size = mgrid_rspace(1)%pw_grid%para%group%num_pe

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
            IF (.NOT. paw_atom) CYCLE

            NULLIFY (cneo_potential)
            CALL get_qs_kind(qs_kind_set(ikind), cneo_potential=cneo_potential)
            IF (.NOT. ASSOCIATED(cneo_potential)) CYCLE

            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
            NULLIFY (nuc_soft_basis)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=nuc_soft_basis, basis_type="NUC_SOFT")
            CALL get_gto_basis_set(gto_basis_set=nuc_soft_basis, lmax=la_max, &
                                   lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                   sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)
            CALL get_cneo_potential(cneo_potential, zeff=zeff)
            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it2(m1, m1))

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               IF (rhoz_cneo_set(atom_a)%ready) THEN
                  ra(:) = pbc(particle_set(atom_a)%r, cell)
                  p_block => rhoz_cneo_set(atom_a)%pmat
                  DO iset = 1, nseta
                     DO jset = 1, iset
                        ! processor mapping
                        map_it2 = .FALSE.
                        DO ipgf = 1, npgfa(iset)
                           IF (jset == iset) THEN
                              npgf2 = ipgf
                           ELSE
                              npgf2 = npgfa(jset)
                           END IF
                           DO jpgf = 1, npgf2
                              zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                              igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                              rs_grid => rs_rho(igrid_level)
                              map_it2(ipgf, jpgf) = map_gaussian_here(rs_grid, cell%h_inv, ra, offset, group_size, my_pos)
                           END DO
                        END DO
                        ! skip empty sets (not uncommon for soft nuclear basis)
                        IF (npgfa(iset) > 0 .AND. npgfa(jset) > 0) THEN
                           offset = offset + 1
                        END IF
                        !
                        IF (ANY(map_it2(1:npgfa(iset), 1:npgfa(jset)))) THEN
                           ncoa = npgfa(iset)*ncoset(la_max(iset))
                           sgfa = first_sgfa(1, iset)
                           ncob = npgfa(jset)*ncoset(la_max(jset))
                           sgfb = first_sgfa(1, jset)
                           ! decontract density block
                           CALL dgemm("N", "N", ncoa, nsgfa(jset), nsgfa(iset), &
                                      1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                      p_block(sgfa, sgfb), SIZE(p_block, 1), &
                                      0.0_dp, work(1, 1), maxco)
                           CALL dgemm("N", "T", ncoa, ncob, nsgfa(jset), &
                                      1.0_dp, work(1, 1), maxco, &
                                      sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                      0.0_dp, pab(1, 1), maxco)
                           DO ipgf = 1, npgfa(iset)
                              IF (jset == iset) THEN
                                 npgf2 = ipgf
                              ELSE
                                 npgf2 = npgfa(jset)
                              END IF
                              DO jpgf = 1, npgf2
                                 IF (map_it2(ipgf, jpgf)) THEN
                                    zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                                    igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                                    rs_grid => rs_rho(igrid_level)

                                    na1 = (ipgf - 1)*ncoset(la_max(iset))
                                    nb1 = (jpgf - 1)*ncoset(la_max(jset))

                                    radius = exp_radius_very_extended(la_min=la_min(iset), &
                                                                      la_max=la_max(iset), &
                                                                      lb_min=la_min(jset), &
                                                                      lb_max=la_max(jset), &
                                                                      ra=ra, rb=ra, rp=ra, &
                                                                      zetp=zetp, eps=eps_rho_rspace, &
                                                                      prefactor=zeff, cutoff=1.0_dp)

                                    IF (jset == iset .AND. jpgf == ipgf) THEN
                                       scale = -zeff ! nuclear charge density is positive
                                    ELSE
                                       scale = -2.0_dp*zeff ! symmetric density matrix
                                    END IF
                                    CALL collocate_pgf_product( &
                                       la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                       la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                       ra, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                       scale, pab, na1, nb1, rs_grid, &
                                       radius=radius, ga_gb_function=GRID_FUNC_AB)
                                 END IF
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END IF
            END DO
            DEALLOCATE (map_it2)
         END DO
         DEALLOCATE (pab, work)

         NULLIFY (rhoz_cneo_s_gs, rhoz_cneo_s_rs)
         CALL get_rho0_mpole(rho0_mpole=rho0, &
                             rhoz_cneo_s_gs=rhoz_cneo_s_gs, &
                             rhoz_cneo_s_rs=rhoz_cneo_s_rs)

         CALL pw_zero(rhoz_cneo_s_gs)
         CALL pw_zero(rhoz_cneo_s_rs)

         DO igrid_level = 1, gridlevel_info%ngrid_levels
            CALL pw_zero(mgrid_rspace(igrid_level))
            CALL transfer_rs2pw(rs=rs_rho(igrid_level), &
                                pw=mgrid_rspace(igrid_level))
         END DO

         DO igrid_level = 1, gridlevel_info%ngrid_levels
            CALL pw_zero(mgrid_gspace(igrid_level))
            CALL pw_transfer(mgrid_rspace(igrid_level), &
                             mgrid_gspace(igrid_level))
            CALL pw_axpy(mgrid_gspace(igrid_level), rhoz_cneo_s_gs)
         END DO
         CALL pw_transfer(rhoz_cneo_s_gs, rhoz_cneo_s_rs)
         tot_rs_int = pw_integrate_function(rhoz_cneo_s_rs, isign=-1)

         ! *** give back the multi-grids *** !
         CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
         CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)
      ELSE
         tot_rs_int = 0.0_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE put_rhoz_cneo_s_on_grid

! **************************************************************************************************
!> \brief ...
!> \param pw_env ...
!> \param rho0_mpole ...
! **************************************************************************************************
   SUBROUTINE rhoz_cneo_s_grid_create(pw_env, rho0_mpole)

      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole

      CHARACTER(len=*), PARAMETER :: routineN = 'rhoz_cneo_s_grid_create'

      INTEGER                                            :: handle
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CPASSERT(ASSOCIATED(auxbas_pw_pool))

      ! reallocate rho0 on the global grid in real and reciprocal space
      CPASSERT(ASSOCIATED(rho0_mpole))

      ! soft nuclear charge density in real space
      IF (ASSOCIATED(rho0_mpole%rhoz_cneo_s_rs)) THEN
         CALL rho0_mpole%rhoz_cneo_s_rs%release()
      ELSE
         ALLOCATE (rho0_mpole%rhoz_cneo_s_rs)
      END IF
      CALL auxbas_pw_pool%create_pw(rho0_mpole%rhoz_cneo_s_rs)

      ! soft nuclear charge density in reciprocal space
      IF (ASSOCIATED(rho0_mpole%rhoz_cneo_s_gs)) THEN
         CALL rho0_mpole%rhoz_cneo_s_gs%release()
      ELSE
         ALLOCATE (rho0_mpole%rhoz_cneo_s_gs)
      END IF
      CALL auxbas_pw_pool%create_pw(rho0_mpole%rhoz_cneo_s_gs)

      CALL timestop(handle)

   END SUBROUTINE rhoz_cneo_s_grid_create

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param v_rspace ...
!> \param para_env ...
!> \param calculate_forces ...
!> \param rhoz_cneo_set ...
!> \param kforce ...
! **************************************************************************************************
   SUBROUTINE integrate_vhgg_rspace(qs_env, v_rspace, para_env, calculate_forces, rhoz_cneo_set, &
                                    kforce)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: v_rspace
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(rhoz_cneo_type), DIMENSION(:), POINTER        :: rhoz_cneo_set
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: kforce

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_vhgg_rspace'

      INTEGER :: atom_a, group_size, handle, iatom, igrid_level, ikind, ipgf, iset, jpgf, jset, &
         m1, maxco, maxsgf_set, my_pos, na1, natom, nb1, ncoa, ncob, npgf2, nseta, offset, sgfa, &
         sgfb
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: paw_atom, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: map_it2
      REAL(KIND=dp)                                      :: eps_rho_rspace, f0, fscale, radius, &
                                                            zeff, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, p_block, pab, sphi_a, vmat, work, &
                                                            zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cneo_potential_type), POINTER                 :: cneo_potential
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: nuc_soft_basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), DIMENSION(:), POINTER   :: rs_rho
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, qs_kind_set, atom_list, cell, dft_control, &
               first_sgfa, force, gridlevel_info, hab, la_max, la_min, &
               nuc_soft_basis, npgfa, nsgfa, p_block, pab, particle_set, &
               pw_env, rs_grid, rs_rho, sphi_a, virial, vmat, work, zeta)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      force=force, pw_env=pw_env, &
                      particle_set=particle_set, &
                      virial=virial)

      ! find maximum numbers
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type="NUC_SOFT")
      IF (maxco > 0) THEN
         fscale = 1.0_dp
         IF (PRESENT(kforce)) THEN
            fscale = kforce
         END IF
         ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set), hab(maxco, maxco))
         pab = 0.0_dp

         eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         CPASSERT(ASSOCIATED(pw_env))
         CALL pw_env_get(pw_env=pw_env, rs_grids=rs_rho, gridlevel_info=gridlevel_info)

         ! transform the potential on the rs_multigrids
         CALL potential_pw2rs(rs_rho, v_rspace, pw_env)

         offset = 0
         my_pos = rs_rho(1)%desc%my_pos
         group_size = rs_rho(1)%desc%group_size

         DO ikind = 1, SIZE(atomic_kind_set, 1)
            CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
            IF (.NOT. paw_atom) CYCLE

            NULLIFY (cneo_potential)
            CALL get_qs_kind(qs_kind_set(ikind), cneo_potential=cneo_potential)
            IF (.NOT. ASSOCIATED(cneo_potential)) CYCLE

            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            NULLIFY (nuc_soft_basis)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=nuc_soft_basis, basis_type="NUC_SOFT")
            CALL get_gto_basis_set(gto_basis_set=nuc_soft_basis, lmax=la_max, &
                                   lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                   sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)
            CALL get_cneo_potential(cneo_potential, zeff=zeff)
            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it2(m1, m1))

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               IF (rhoz_cneo_set(atom_a)%ready) THEN
                  p_block => rhoz_cneo_set(atom_a)%pmat
               ELSE
                  NULLIFY (p_block)
               END IF
               vmat => rhoz_cneo_set(atom_a)%vmat
               vmat = 0.0_dp
               DO iset = 1, nseta
                  DO jset = 1, iset
                     ! processor mapping
                     map_it2 = .FALSE.
                     DO ipgf = 1, npgfa(iset)
                        IF (jset == iset) THEN
                           npgf2 = ipgf
                        ELSE
                           npgf2 = npgfa(jset)
                        END IF
                        DO jpgf = 1, npgf2
                           zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                           igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                           rs_grid => rs_rho(igrid_level)
                           map_it2(ipgf, jpgf) = map_gaussian_here(rs_grid, cell%h_inv, ra, offset, group_size, my_pos)
                        END DO
                     END DO
                     ! skip empty sets (not uncommon for soft nuclear basis)
                     IF (npgfa(iset) > 0 .AND. npgfa(jset) > 0) THEN
                        offset = offset + 1
                     END IF
                     !
                     IF (ANY(map_it2(1:npgfa(iset), 1:npgfa(jset)))) THEN
                        hab = 0.0_dp
                        IF (calculate_forces) THEN
                           force_a = 0.0_dp
                           force_b = 0.0_dp
                        END IF
                        IF (use_virial) THEN
                           my_virial_a = 0.0_dp
                           my_virial_b = 0.0_dp
                        END IF
                        ncoa = npgfa(iset)*ncoset(la_max(iset))
                        sgfa = first_sgfa(1, iset)
                        ncob = npgfa(jset)*ncoset(la_max(jset))
                        sgfb = first_sgfa(1, jset)
                        IF (calculate_forces .AND. ASSOCIATED(p_block)) THEN
                           ! decontract density block
                           CALL dgemm("N", "N", ncoa, nsgfa(jset), nsgfa(iset), &
                                      1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                      p_block(sgfa, sgfb), SIZE(p_block, 1), &
                                      0.0_dp, work(1, 1), maxco)
                           CALL dgemm("N", "T", ncoa, ncob, nsgfa(jset), &
                                      1.0_dp, work(1, 1), maxco, &
                                      sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                      0.0_dp, pab(1, 1), maxco)
                        END IF
                        DO ipgf = 1, npgfa(iset)
                           IF (jset == iset) THEN
                              npgf2 = ipgf
                           ELSE
                              npgf2 = npgfa(jset)
                           END IF
                           DO jpgf = 1, npgf2
                              IF (map_it2(ipgf, jpgf)) THEN
                                 zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                                 igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                                 rs_grid => rs_rho(igrid_level)

                                 na1 = (ipgf - 1)*ncoset(la_max(iset))
                                 nb1 = (jpgf - 1)*ncoset(la_max(jset))

                                 radius = exp_radius_very_extended(la_min=la_min(iset), &
                                                                   la_max=la_max(iset), &
                                                                   lb_min=la_min(jset), &
                                                                   lb_max=la_max(jset), &
                                                                   ra=ra, rb=ra, rp=ra, &
                                                                   zetp=zetp, eps=eps_rho_rspace, &
                                                                   prefactor=zeff, cutoff=1.0_dp)

                                 CALL integrate_pgf_product( &
                                    la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                    la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                    ra, (/0.0_dp, 0.0_dp, 0.0_dp/), rs_grid, &
                                    hab, pab=pab, o1=na1, o2=nb1, &
                                    radius=radius, &
                                    calculate_forces=calculate_forces, force_a=force_a, force_b=force_b, &
                                    use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b)
                                 f0 = 0.0_dp
                                 IF (calculate_forces .OR. use_virial) THEN
                                    IF (jset == iset .AND. jpgf == ipgf) THEN
                                       f0 = -fscale*zeff
                                    ELSE
                                       f0 = -2.0_dp*fscale*zeff
                                    END IF
                                 END IF
                                 IF (calculate_forces) THEN
                                    force(ikind)%rho_cneo_nuc(1:3, iatom) = &
                                       force(ikind)%rho_cneo_nuc(1:3, iatom) + f0*(force_a + force_b)
                                 END IF
                                 IF (use_virial) THEN
                                    virial%pv_virial = virial%pv_virial + f0*(my_virial_a + my_virial_b)
                                 END IF
                                 ! symmetrize
                                 IF (iset == jset .AND. ipgf /= jpgf) THEN
                                    hab(nb1 + 1:nb1 + ncoset(la_max(jset)), &
                                        na1 + 1:na1 + ncoset(la_max(iset))) = &
                                       TRANSPOSE(hab(na1 + 1:na1 + ncoset(la_max(iset)), &
                                                     nb1 + 1:nb1 + ncoset(la_max(jset))))
                                 END IF
                              END IF
                           END DO
                        END DO
                        ! contract the soft basis V_Hartree integral
                        work(1:ncoa, 1:nsgfa(jset)) = MATMUL(hab(1:ncoa, 1:ncob), &
                                                             sphi_a(1:ncob, sgfb:sgfb + nsgfa(jset) - 1))
                        vmat(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfa(jset) - 1) = &
                           vmat(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfa(jset) - 1) - zeff* &
                           MATMUL(TRANSPOSE(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1)), &
                                  work(1:ncoa, 1:nsgfa(jset)))
                        ! symmetrize
                        IF (iset /= jset) THEN
                           vmat(sgfb:sgfb + nsgfa(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) = &
                              TRANSPOSE(vmat(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfa(jset) - 1))
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            DEALLOCATE (map_it2)
         END DO
         DEALLOCATE (pab, work, hab)

         DO ikind = 1, SIZE(atomic_kind_set, 1)
            CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
            IF (.NOT. paw_atom) CYCLE

            NULLIFY (cneo_potential)
            CALL get_qs_kind(qs_kind_set(ikind), cneo_potential=cneo_potential)
            IF (.NOT. ASSOCIATED(cneo_potential)) CYCLE

            NULLIFY (atom_list)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               vmat => rhoz_cneo_set(atom_a)%vmat
               CALL para_env%sum(vmat)
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE integrate_vhgg_rspace

END MODULE qs_cneo_ggrid
