!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE qs_rho0_ggrid
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
   USE grid_api,                        ONLY: GRID_FUNC_AB,&
                                              collocate_pgf_product
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_release,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_local_rho_types,              ONLY: get_local_rho,&
                                              local_rho_type
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              rho0_mpole_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_desc_type,&
                                              realspace_grid_type,&
                                              rs_grid_create,&
                                              rs_grid_release,&
                                              rs_grid_zero,&
                                              transfer_pw2rs,&
                                              transfer_rs2pw
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho0_ggrid'

   ! Public subroutines

   PUBLIC :: put_rho0_on_grid, rho0_s_grid_create, integrate_vhg0_rspace

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0 ...
!> \param tot_rs_int ...
! **************************************************************************************************
   SUBROUTINE put_rho0_on_grid(qs_env, rho0, tot_rs_int)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho0_mpole_type), POINTER                     :: rho0
      REAL(KIND=dp), INTENT(OUT)                         :: tot_rs_int

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'put_rho0_on_grid'

      INTEGER                                            :: auxbas_grid, handle, iat, iatom, igrid, &
                                                            ikind, ithread, j, l0_ikind, lmax0, &
                                                            nat, nch_ik, nch_max, npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: paw_atom
      REAL(KIND=dp)                                      :: eps_rho_rspace, rpgf0, zet0
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: Qlm_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_type)                                      :: coeff_gspace, coeff_rspace, rho0_r_tmp
      TYPE(pw_type), POINTER                             :: rho0_s_gs, rho0_s_rs
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: descs
      TYPE(realspace_grid_desc_type), POINTER            :: desc
      TYPE(realspace_grid_type), DIMENSION(:), POINTER   :: grids
      TYPE(realspace_grid_type), POINTER                 :: rs_grid

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, qs_kind_set, cores, pab, Qlm_c)

      NULLIFY (dft_control, pw_env, particle_set, para_env, cell, rho0_s_gs, rho0_s_rs)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      pw_env=pw_env, cell=cell)
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      NULLIFY (descs, pw_pools)
      CALL pw_env_get(pw_env=pw_env, rs_descs=descs, rs_grids=grids, pw_pools=pw_pools)
      auxbas_grid = pw_env%auxbas_grid

      NULLIFY (rho0_s_gs, rho0_s_rs)
      CALL get_rho0_mpole(rho0_mpole=rho0, lmax_0=lmax0, &
                          zet0_h=zet0, igrid_zet0_s=igrid, &
                          rho0_s_gs=rho0_s_gs, &
                          rho0_s_rs=rho0_s_rs)

      ! *** set up the rs grid at level igrid
      NULLIFY (rs_grid, desc, pw_pool)
      desc => descs(igrid)%rs_desc
      rs_grid => grids(igrid)
      pw_pool => pw_pools(igrid)%pool

      CPASSERT(ASSOCIATED(desc))
      CPASSERT(ASSOCIATED(pw_pool))

      IF (igrid /= auxbas_grid) THEN
         CALL pw_pool_create_pw(pw_pool, coeff_rspace, use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(pw_pool, coeff_gspace, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
      END IF
      CALL rs_grid_zero(rs_grid)

      nch_max = ncoset(lmax0)

      ALLOCATE (pab(nch_max, 1))

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)

         IF (.NOT. paw_atom .AND. dft_control%qs_control%gapw_control%nopaw_as_gpw) CYCLE

         CALL get_rho0_mpole(rho0_mpole=rho0, ikind=ikind, l0_ikind=l0_ikind, &
                             rpgf0_s=rpgf0)

         nch_ik = ncoset(l0_ikind)
         pab = 0.0_dp

         CALL reallocate(cores, 1, nat)
         npme = 0
         cores = 0

         DO iat = 1, nat
            iatom = atom_list(iat)
            ra(:) = pbc(particle_set(iatom)%r, cell)
            IF (rs_grid%desc%parallel .AND. .NOT. rs_grid%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(nat, rs_grid%desc%group_size) == rs_grid%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iat
               END IF
            ELSE
               npme = npme + 1
               cores(npme) = iat
            END IF

         END DO

         ithread = 0
         DO j = 1, npme

            iat = cores(j)
            iatom = atom_list(iat)

            CALL get_rho0_mpole(rho0_mpole=rho0, iat=iatom, Qlm_car=Qlm_c)

            pab(1:nch_ik, 1) = Qlm_c(1:nch_ik)

            ra(:) = pbc(particle_set(iatom)%r, cell)

            CALL collocate_pgf_product( &
               l0_ikind, zet0, 0, 0, 0.0_dp, 0, &
               ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 1.0_dp, pab, 0, 0, &
               rs_grid, ga_gb_function=GRID_FUNC_AB, radius=rpgf0, &
               use_subpatch=.TRUE., subpatch_pattern=0)

         END DO ! j

      END DO ! ikind

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      DEALLOCATE (pab)

      IF (igrid /= auxbas_grid) THEN
         CALL transfer_rs2pw(rs_grid, coeff_rspace)
         CALL pw_zero(rho0_s_gs)
         CALL pw_transfer(coeff_rspace, coeff_gspace)
         CALL pw_axpy(coeff_gspace, rho0_s_gs)

         tot_rs_int = pw_integrate_function(coeff_rspace, isign=-1)

         CALL pw_pool_give_back_pw(pw_pool, coeff_rspace)
         CALL pw_pool_give_back_pw(pw_pool, coeff_gspace)
      ELSE

         CALL pw_pool_create_pw(pw_pool, rho0_r_tmp, &
                                use_data=REALDATA3D, in_space=REALSPACE)

         CALL transfer_rs2pw(rs_grid, rho0_r_tmp)

         tot_rs_int = pw_integrate_function(rho0_r_tmp, isign=-1)

         CALL pw_transfer(rho0_r_tmp, rho0_s_rs)
         CALL pw_pool_give_back_pw(pw_pool, rho0_r_tmp)

         CALL pw_zero(rho0_s_gs)
         CALL pw_transfer(rho0_s_rs, rho0_s_gs)
      END IF
      CALL timestop(handle)

   END SUBROUTINE put_rho0_on_grid

! **************************************************************************************************
!> \brief ...
!> \param pw_env ...
!> \param rho0_mpole ...
! **************************************************************************************************
   SUBROUTINE rho0_s_grid_create(pw_env, rho0_mpole)

      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole

      CHARACTER(len=*), PARAMETER :: routineN = 'rho0_s_grid_create'

      INTEGER                                            :: handle
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CPASSERT(ASSOCIATED(auxbas_pw_pool))

      ! reallocate rho0 on the global grid in real and reciprocal space
      CPASSERT(ASSOCIATED(rho0_mpole))

      ! rho0 density in real space
      IF (ASSOCIATED(rho0_mpole%rho0_s_rs)) THEN
         CALL pw_release(rho0_mpole%rho0_s_rs)
      ELSE
         ALLOCATE (rho0_mpole%rho0_s_rs)
      END IF
      CALL pw_pool_create_pw(auxbas_pw_pool, rho0_mpole%rho0_s_rs, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      ! rho0 density in reciprocal space
      IF (ASSOCIATED(rho0_mpole%rho0_s_gs)) THEN
         CALL pw_release(rho0_mpole%rho0_s_gs)
      ELSE
         ALLOCATE (rho0_mpole%rho0_s_gs)
      END IF
      CALL pw_pool_create_pw(auxbas_pw_pool, rho0_mpole%rho0_s_gs, &
                             use_data=COMPLEXDATA1D)
      rho0_mpole%rho0_s_gs%in_space = RECIPROCALSPACE

      ! Find the grid level suitable for rho0_soft
      rho0_mpole%igrid_zet0_s = gaussian_gridlevel(pw_env%gridlevel_info, 2.0_dp*rho0_mpole%zet0_h)

      CALL timestop(handle)

   END SUBROUTINE rho0_s_grid_create

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param v_rspace ...
!> \param para_env ...
!> \param calculate_forces ...
!> \param local_rho_set ...
!> \param local_rho_set_2nd ...
!> \param atener ...
! **************************************************************************************************
   SUBROUTINE integrate_vhg0_rspace(qs_env, v_rspace, para_env, calculate_forces, local_rho_set, &
                                    local_rho_set_2nd, atener)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_type), INTENT(IN)                          :: v_rspace
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(local_rho_type), OPTIONAL, POINTER            :: local_rho_set, local_rho_set_2nd
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atener

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_vhg0_rspace'

      INTEGER :: auxbas_grid, bo(2), handle, i, iat, iatom, ic, icg, ico, ig1, ig2, igrid, ii, &
         ikind, ipgf1, ipgf2, is, iset1, iset2, iso, iso1, iso2, ispin, j, l0_ikind, llmax, lmax0, &
         lshell, lx, ly, lz, m1, m2, max_iso_not0_local, max_s_harm, maxl, maxso, mepos, n1, n2, &
         nat, nch_ik, nch_max, ncurr, nset, nsotot, nspins, num_pe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf
      LOGICAL                                            :: grid_distributed, paw_atom, use_virial
      REAL(KIND=dp)                                      :: eps_rho_rspace, force_tmp(3), ra(3), &
                                                            rpgf0, scale, zet0
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: hab_sph, norm_l, Qlm
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, hdab_sph, intloc, pab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: a_hdab_sph, hdab, Qlm_gg
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: a_hdab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_set
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: pw_aux, pw_pool
      TYPE(pw_type)                                      :: coeff_gaux, coeff_gspace, coeff_raux, &
                                                            coeff_rspace
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_type)                          :: rs_v
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_local_h, int_local_s
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! In case of the external density computed forces probably also
      ! need to be stored outside qs_env. We can then remove the
      ! attribute 'OPTIONAL' from the argument 'local_rho_set'.

      ! CPASSERT(.NOT. (calculate_forces .AND. PRESENT(local_rho_set)))
!      IF (calculate_forces .AND. PRESENT(local_rho_set)) THEN
!         CPWARN("Forces and External Density!")
!      END IF

      NULLIFY (atomic_kind_set, qs_kind_set, dft_control, particle_set)
      NULLIFY (cell, force, pw_env, rho0_mpole, rho_atom_set)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      force=force, pw_env=pw_env, &
                      rho0_mpole=rho0_mpole, &
                      rho_atom_set=rho_atom_set, &
                      particle_set=particle_set, &
                      virial=virial)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      nspins = dft_control%nspins

      ! The aim of the following code was to return immediately if the subroutine
      ! was called for triplet excited states in spin-restricted case. This check
      ! is also performed before invocation of this subroutine. It should be save
      ! to remove the optional argument 'do_triplet' from the subroutine interface.
      !my_tddft = PRESENT(local_rho_set)
      !IF (my_tddft) THEN
      !   IF (PRESENT(do_triplet)) THEN
      !      IF (nspins == 1 .AND. do_triplet) RETURN
      !   ELSE
      !      IF (nspins == 1 .AND. dft_control%tddfpt_control%res_etype /= tddfpt_singlet) RETURN
      !   END IF
      !END IF

      IF (PRESENT(local_rho_set)) &
         CALL get_local_rho(local_rho_set, rho0_mpole=rho0_mpole, rho_atom_set=rho_atom_set)
      ! Q from rho0_mpole of local_rho_set
      ! for TDDFT forces we need mixed potential / integral space
      ! potential stored on local_rho_set_2nd
      IF (PRESENT(local_rho_set_2nd)) THEN
         CALL get_local_rho(local_rho_set_2nd, rho_atom_set=rho_atom_set)
      END IF
      CALL get_rho0_mpole(rho0_mpole=rho0_mpole, lmax_0=lmax0, &
                          zet0_h=zet0, igrid_zet0_s=igrid, &
                          norm_g0l_h=norm_l)

      ! Setup of the potential on the multigrids
      NULLIFY (rs_descs, pw_pools)
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_descs=rs_descs, pw_pools=pw_pools)

      ! Assign from pw_env
      auxbas_grid = pw_env%auxbas_grid

      ! Get the potential on the right grid
      NULLIFY (rs_desc, pw_pool, pw_aux)
      rs_desc => rs_descs(igrid)%rs_desc
      pw_pool => pw_pools(igrid)%pool

      CALL pw_pool_create_pw(pw_pool, coeff_gspace, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL pw_pool_create_pw(pw_pool, coeff_rspace, use_data=REALDATA3D, &
                             in_space=REALSPACE)

      IF (igrid /= auxbas_grid) THEN
         pw_aux => pw_pools(auxbas_grid)%pool
         CALL pw_pool_create_pw(pw_aux, coeff_gaux, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_transfer(v_rspace, coeff_gaux)
         CALL pw_copy(coeff_gaux, coeff_gspace)
         CALL pw_transfer(coeff_gspace, coeff_rspace)
         CALL pw_pool_give_back_pw(pw_aux, coeff_gaux)
         CALL pw_pool_create_pw(pw_aux, coeff_raux, use_data=REALDATA3D, &
                                in_space=REALSPACE)
         scale = coeff_rspace%pw_grid%dvol/coeff_raux%pw_grid%dvol
         coeff_rspace%cr3d = scale*coeff_rspace%cr3d
         CALL pw_pool_give_back_pw(pw_aux, coeff_raux)
      ELSE

         IF (coeff_gspace%pw_grid%spherical) THEN
            CALL pw_transfer(v_rspace, coeff_gspace)
            CALL pw_transfer(coeff_gspace, coeff_rspace)
         ELSE
            CALL pw_copy(v_rspace, coeff_rspace)
         END IF
      END IF
      CALL pw_pool_give_back_pw(pw_pool, coeff_gspace)

      ! Setup the rs grid at level igrid
      CALL rs_grid_create(rs_v, rs_desc)
      CALL rs_grid_zero(rs_v)
      CALL transfer_pw2rs(rs_v, coeff_rspace)

      CALL pw_pool_give_back_pw(pw_pool, coeff_rspace)

      ! Now the potential is on the right grid => integration

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      ! Allocate work storage

      NULLIFY (hab, hab_sph, hdab, hdab_sph, pab, a_hdab, a_hdab_sph)
      nch_max = ncoset(lmax0)
      CALL reallocate(hab, 1, nch_max, 1, 1)
      CALL reallocate(hab_sph, 1, nch_max)
      CALL reallocate(hdab, 1, 3, 1, nch_max, 1, 1)
      CALL reallocate(hdab_sph, 1, 3, 1, nch_max)
      CALL reallocate(a_hdab, 1, 3, 1, 3, 1, nch_max, 1, 1)
      CALL reallocate(a_hdab_sph, 1, 3, 1, 3, 1, nch_max)
      CALL reallocate(pab, 1, nch_max, 1, 1)

      ncurr = -1

      grid_distributed = rs_v%desc%distributed

      DO ikind = 1, SIZE(atomic_kind_set, 1)
         NULLIFY (basis_1c_set, atom_list, harmonics)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=basis_1c_set, basis_type="GAPW_1C", &
                          paw_atom=paw_atom, &
                          harmonics=harmonics)

         IF (.NOT. paw_atom) CYCLE

         NULLIFY (Qlm_gg, lmax, npgf)
         CALL get_rho0_mpole(rho0_mpole=rho0_mpole, ikind=ikind, &
                             l0_ikind=l0_ikind, Qlm_gg=Qlm_gg, &  ! Qs different
                             rpgf0_s=rpgf0)

         CALL get_gto_basis_set(gto_basis_set=basis_1c_set, &
                                lmax=lmax, lmin=lmin, &
                                maxso=maxso, maxl=maxl, &
                                nset=nset, npgf=npgf)

         nsotot = maxso*nset
         ALLOCATE (intloc(nsotot, nsotot))

         ! Initialize the local KS integrals

         nch_ik = ncoset(l0_ikind)
         pab = 1.0_dp
         max_s_harm = harmonics%max_s_harm
         llmax = harmonics%llmax

         ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

         num_pe = para_env%num_pe
         mepos = para_env%mepos
         DO j = 0, num_pe - 1
            bo = get_limit(nat, num_pe, j)
            IF (.NOT. grid_distributed .AND. j /= mepos) CYCLE

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               ra(:) = pbc(particle_set(iatom)%r, cell)

               NULLIFY (Qlm)
               CALL get_rho0_mpole(rho0_mpole=rho0_mpole, iat=iatom, Qlm_tot=Qlm)

               hab = 0.0_dp
               hdab = 0.0_dp
               intloc = 0._dp
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
                  a_hdab = 0.0_dp
               END IF

               CALL integrate_pgf_product( &
                  l0_ikind, zet0, 0, 0, 0.0_dp, 0, &
                  ra, (/0.0_dp, 0.0_dp, 0.0_dp/), rs_v, &
                  hab, pab, o1=0, o2=0, &
                  radius=rpgf0, &
                  calculate_forces=calculate_forces, &
                  use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b, &
                  hdab=hdab, a_hdab=a_hdab, use_subpatch=.TRUE., subpatch_pattern=0)

               ! Convert from cartesian to spherical
               DO lshell = 0, l0_ikind
                  DO is = 1, nso(lshell)
                     iso = is + nsoset(lshell - 1)
                     hab_sph(iso) = 0.0_dp
                     hdab_sph(1:3, iso) = 0.0_dp
                     a_hdab_sph(1:3, 1:3, iso) = 0.0_dp
                     DO ic = 1, nco(lshell)
                        ico = ic + ncoset(lshell - 1)
                        lx = indco(1, ico)
                        ly = indco(2, ico)
                        lz = indco(3, ico)
                        hab_sph(iso) = hab_sph(iso) + &
                                       norm_l(lshell)* &
                                       orbtramat(lshell)%slm(is, ic)* &
                                       hab(ico, 1)
                        IF (calculate_forces) THEN
                           hdab_sph(1:3, iso) = hdab_sph(1:3, iso) + &
                                                norm_l(lshell)* &
                                                orbtramat(lshell)%slm(is, ic)* &
                                                hdab(1:3, ico, 1)
                        END IF
                        IF (use_virial) THEN
                           DO ii = 1, 3
                           DO i = 1, 3
                              a_hdab_sph(i, ii, iso) = a_hdab_sph(i, ii, iso) + &
                                                       norm_l(lshell)* &
                                                       orbtramat(lshell)%slm(is, ic)* &
                                                       a_hdab(i, ii, ico, 1)
                           END DO
                           END DO
                        END IF

                     END DO ! ic
                  END DO ! is
               END DO ! lshell

               m1 = 0
               DO iset1 = 1, nset

                  m2 = 0
                  DO iset2 = 1, nset
                     CALL get_none0_cg_list(harmonics%my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                            max_s_harm, llmax, cg_list, cg_n_list, max_iso_not0_local)
                     n1 = nsoset(lmax(iset1))
                     DO ipgf1 = 1, npgf(iset1)
                        n2 = nsoset(lmax(iset2))
                        DO ipgf2 = 1, npgf(iset2)

                           DO iso = 1, MIN(nsoset(l0_ikind), max_iso_not0_local)
                              DO icg = 1, cg_n_list(iso)
                                 iso1 = cg_list(1, icg, iso)
                                 iso2 = cg_list(2, icg, iso)

                                 ig1 = iso1 + n1*(ipgf1 - 1) + m1
                                 ig2 = iso2 + n2*(ipgf2 - 1) + m2

                                 intloc(ig1, ig2) = intloc(ig1, ig2) + Qlm_gg(ig1, ig2, iso)*hab_sph(iso) ! potential times Q

                              END DO ! icg
                           END DO ! iso

                        END DO ! ipgf2
                     END DO ! ipgf1
                     m2 = m2 + maxso
                  END DO ! iset2
                  m1 = m1 + maxso
               END DO ! iset1

               IF (grid_distributed) THEN
                  ! Sum result over all processors
                  CALL para_env%sum(intloc)
               END IF

               IF (j == mepos) THEN
                  rho_atom => rho_atom_set(iatom)
                  CALL get_rho_atom(rho_atom=rho_atom, ga_Vlocal_gb_h=int_local_h, ga_Vlocal_gb_s=int_local_s)
                  DO ispin = 1, nspins
                     int_local_h(ispin)%r_coef = int_local_h(ispin)%r_coef + intloc
                     int_local_s(ispin)%r_coef = int_local_s(ispin)%r_coef + intloc
                  END DO
               END IF

               IF (PRESENT(atener)) THEN
                  DO iso = 1, nsoset(l0_ikind)
                     atener(iatom) = atener(iatom) + 0.5_dp*Qlm(iso)*hab_sph(iso)
                  END DO
               END IF

               IF (calculate_forces) THEN
                  force_tmp(1:3) = 0.0_dp
                  DO iso = 1, nsoset(l0_ikind)
                     force_tmp(1) = force_tmp(1) + Qlm(iso)*hdab_sph(1, iso) ! Q here is from local_rho_set
                     force_tmp(2) = force_tmp(2) + Qlm(iso)*hdab_sph(2, iso)
                     force_tmp(3) = force_tmp(3) + Qlm(iso)*hdab_sph(3, iso)
                  END DO
                  force(ikind)%g0s_Vh_elec(1:3, iat) = force(ikind)%g0s_Vh_elec(1:3, iat) + force_tmp(1:3)
               END IF
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  DO iso = 1, nsoset(l0_ikind)
                     DO ii = 1, 3
                     DO i = 1, 3
                        ! Q from local_rho_set
                        virial%pv_gapw(i, ii) = virial%pv_gapw(i, ii) + Qlm(iso)*a_hdab_sph(i, ii, iso)
                        virial%pv_virial(i, ii) = virial%pv_virial(i, ii) + Qlm(iso)*a_hdab_sph(i, ii, iso)
                     END DO
                     END DO
                  END DO
               END IF

            END DO
         END DO

         DEALLOCATE (intloc)
         DEALLOCATE (cg_list, cg_n_list)

      END DO ! ikind

      CALL rs_grid_release(rs_v)

      DEALLOCATE (hab, hdab, hab_sph, hdab_sph, pab, a_hdab, a_hdab_sph)

      CALL timestop(handle)

   END SUBROUTINE integrate_vhg0_rspace

END MODULE qs_rho0_ggrid
