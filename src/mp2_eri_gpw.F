!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate 2c- and 3c-integrals for RI with GPW
!> \par History
!>      07.2019 separated from mp2_integrals.F [Frederick Stein]
! **************************************************************************************************
MODULE mp2_eri_gpw
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_long
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_gauss_damp,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              rs_grid_release,&
                                              rs_grid_retain
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_eri_gpw'

   PUBLIC :: mp2_eri_2c_integrate_gpw, mp2_eri_3c_integrate_gpw, calc_potential_gpw, cleanup_gpw, prepare_gpw

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mo_coeff ...
!> \param psi_L ...
!> \param rho_g ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param dft_control ...
!> \param particle_set ...
!> \param pw_env_sub ...
!> \param external_vector ...
!> \param poisson_env ...
!> \param rho_r ...
!> \param pot_g ...
!> \param ri_metric ...
!> \param omega_metric ...
!> \param mat_munu ...
!> \param qs_env ...
!> \param task_list_sub ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_3c_integrate_gpw(mo_coeff, psi_L, rho_g, atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                       pw_env_sub, external_vector, poisson_env, rho_r, pot_g, ri_metric, omega_metric, mat_munu, &
                                       qs_env, task_list_sub)

      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(pw_p_type), INTENT(INOUT)                     :: psi_L, rho_g
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: external_vector
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r, pot_g
      INTEGER, INTENT(IN)                                :: ri_metric
      REAL(KIND=dp), INTENT(IN)                          :: omega_metric
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_eri_3c_integrate_gpw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! pseudo psi_L
      CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                  basis_type="RI_AUX", &
                                  external_vector=external_vector)

      rho_r%pw%cr3d = psi_L%pw%cr3d

      CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, ri_metric, omega_metric)

      ! and finally (K|mu nu)
      CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
      CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, &
                              calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                              pw_env_external=pw_env_sub, task_list_external=task_list_sub)

      CALL timestop(handle)

   END SUBROUTINE mp2_eri_3c_integrate_gpw

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env_sub ...
!> \param dimen_RI ...
!> \param mo_coeff ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param ri_metric ...
!> \param natom ...
!> \param omega ...
!> \param sab_orb_sub ...
!> \param L_local_col ...
!> \param kind_of ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_2c_integrate_gpw(qs_env, para_env_sub, dimen_RI, mo_coeff, my_group_L_start, my_group_L_end, ri_metric, &
                                       natom, omega, sab_orb_sub, L_local_col, kind_of)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(IN)                                :: my_group_L_start, my_group_L_end, &
                                                            ri_metric, natom
      REAL(KIND=dp), INTENT(IN)                          :: omega
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: L_local_col
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_eri_2c_integrate_gpw'

      INTEGER :: dir, handle, handle2, i, i_counter, iatom, igrid_level, ikind, ipgf, iset, lb(3), &
         LLL, location(3), na1, na2, ncoa, nseta, offset, sgfa, tp(3), ub(3)
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: map_it_here
      REAL(KIND=dp)                                      :: cutoff_old, radius, relative_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old, wf_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: I_ab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: I_tmp2, rpgfa, sphi_a, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: pot_g, psi_L, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                       auxbas_pw_pool, poisson_env, task_list_sub, rho_r, rho_g, pot_g, psi_L, sab_orb_sub)

      ALLOCATE (wf_vector(dimen_RI))

      CALL get_qs_env(qs_env, cell=cell, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      L_local_col = 0.0_dp

      i_counter = 0
      DO LLL = my_group_L_start, my_group_L_end
         i_counter = i_counter + 1

         wf_vector = 0.0_dp
         wf_vector(LLL) = 1.0_dp

         ! pseudo psi_L
         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                     basis_type="RI_AUX", &
                                     external_vector=wf_vector)

         CALL timeset(routineN//"_pot_lm", handle2)
         rho_r%pw%cr3d = psi_L%pw%cr3d

         CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, ri_metric, omega)

         NULLIFY (rs_v)
         NULLIFY (rs_descs)
         CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
         DO i = 1, SIZE(rs_v)
            CALL rs_grid_retain(rs_v(i)%rs_grid)
         END DO
         CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

         CALL timestop(handle2)

         ! integrate the little bastards
         offset = 0
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")

            first_sgfa => basis_set_a%first_sgf
            la_max => basis_set_a%lmax
            la_min => basis_set_a%lmin
            npgfa => basis_set_a%npgf
            nseta = basis_set_a%nset
            nsgfa => basis_set_a%nsgf_set
            rpgfa => basis_set_a%pgf_radius
            set_radius_a => basis_set_a%set_radius
            sphi_a => basis_set_a%sphi
            zeta => basis_set_a%zet

            ra(:) = pbc(particle_set(iatom)%r, cell)

            DO iset = 1, nseta
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)

               ALLOCATE (I_tmp2(ncoa, 1))
               I_tmp2 = 0.0_dp
               ALLOCATE (I_ab(nsgfa(iset), 1))
               I_ab = 0.0_dp

               igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(zeta(:, iset)))
               map_it_here = .FALSE.
               IF (.NOT. ALL(rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
                  DO dir = 1, 3
                     ! bounds of local grid (i.e. removing the 'wings'), if periodic
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                     tp(dir) = MODULO(tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir))
                     IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                        lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir) + rs_v(igrid_level)%rs_grid%desc%border
                        ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir) - rs_v(igrid_level)%rs_grid%desc%border
                     ELSE
                        lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir)
                        ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir)
                     ENDIF
                     ! distributed grid, only map if it is local to the grid
                     location(dir) = tp(dir) + rs_v(igrid_level)%rs_grid%desc%lb(dir)
                  ENDDO
                  IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                      lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                      lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                     map_it_here = .TRUE.
                  ENDIF
               ELSE
                  ! not distributed, just a round-robin distribution over the full set of CPUs
                  IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
               ENDIF

               offset = offset + nsgfa(iset)

               IF (map_it_here) THEN
                  DO ipgf = 1, npgfa(iset)
                     sgfa = first_sgfa(1, iset)
                     na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                     na2 = ipgf*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, zeta(ipgf, iset))

                     radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                       lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                       zetp=zeta(ipgf, iset), &
                                                       eps=dft_control%qs_control%eps_gvg_rspace, &
                                                       prefactor=1.0_dp, cutoff=1.0_dp)

                     CALL integrate_pgf_product( &
                        la_max=la_max(iset), zeta=zeta(ipgf, iset)/2.0_dp, la_min=la_min(iset), &
                        lb_max=0, zetb=zeta(ipgf, iset)/2.0_dp, lb_min=0, &
                        ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                        rsgrid=rs_v(igrid_level)%rs_grid, &
                        cell=cell, &
                        cube_info=pw_env_sub%cube_info(igrid_level), &
                        hab=I_tmp2, &
                        o1=na1 - 1, &
                        o2=0, &
                        radius=radius, &
                        calculate_forces=.FALSE.)
                  END DO

                  CALL dgemm("T", "N", nsgfa(iset), 1, ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             I_tmp2(1, 1), SIZE(I_tmp2, 1), &
                             1.0_dp, I_ab(1, 1), SIZE(I_ab, 1))

                  L_local_col(offset - nsgfa(iset) + 1:offset, i_counter) = I_ab(1:nsgfa(iset), 1)
               END IF

               DEALLOCATE (I_tmp2)
               DEALLOCATE (I_ab)

            END DO
         END DO

         DO i = 1, SIZE(rs_v)
            CALL rs_grid_release(rs_v(i)%rs_grid)
         END DO

      END DO

      DEALLOCATE (wf_vector)

      CALL mp_sum(L_local_col, para_env_sub%group)

      CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, pw_env_sub, &
                       task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dft_control ...
!> \param e_cutoff_old ...
!> \param cutoff_old ...
!> \param relative_cutoff_old ...
!> \param para_env_sub ...
!> \param pw_env_sub ...
!> \param auxbas_pw_pool ...
!> \param poisson_env ...
!> \param task_list_sub ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param psi_L ...
!> \param sab_orb_sub ...
! **************************************************************************************************
   SUBROUTINE prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                          auxbas_pw_pool, poisson_env, task_list_sub, rho_r, rho_g, pot_g, psi_L, sab_orb_sub)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: e_cutoff_old
      REAL(KIND=dp), INTENT(OUT)                         :: cutoff_old, relative_cutoff_old
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(pw_p_type), INTENT(OUT)                       :: rho_r, rho_g, pot_g, psi_L
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_orb_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_gpw'

      INTEGER                                            :: handle, i_multigrid, n_multigrid
      LOGICAL                                            :: skip_load_balance_distributed
      REAL(KIND=dp)                                      :: progression_factor
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, ks_env=ks_env)

      ! hack hack hack XXXXXXXXXXXXXXX rebuilds the pw_en with the new cutoffs
      progression_factor = dft_control%qs_control%progression_factor
      n_multigrid = SIZE(dft_control%qs_control%e_cutoff)
      ALLOCATE (e_cutoff_old(n_multigrid))
      e_cutoff_old(:) = dft_control%qs_control%e_cutoff
      cutoff_old = dft_control%qs_control%cutoff

      dft_control%qs_control%cutoff = qs_env%mp2_env%mp2_gpw%cutoff*0.5_dp
      dft_control%qs_control%e_cutoff(1) = dft_control%qs_control%cutoff
      DO i_multigrid = 2, n_multigrid
         dft_control%qs_control%e_cutoff(i_multigrid) = dft_control%qs_control%e_cutoff(i_multigrid - 1) &
                                                        /progression_factor
      END DO

      relative_cutoff_old = dft_control%qs_control%relative_cutoff
      dft_control%qs_control%relative_cutoff = qs_env%mp2_env%mp2_gpw%relative_cutoff*0.5_dp

      ! a pw_env
      NULLIFY (pw_env_sub)
      CALL pw_env_create(pw_env_sub)
      CALL pw_env_rebuild(pw_env_sub, qs_env, para_env_sub)

      CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      ! hack hack hack XXXXXXXXXXXXXXX

      ! now we need a task list, hard code skip_load_balance_distributed
      NULLIFY (task_list_sub)
      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      CALL allocate_task_list(task_list_sub)
      CALL generate_qs_task_list(ks_env, task_list_sub, &
                                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                                 skip_load_balance_distributed=skip_load_balance_distributed, &
                                 pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub)

      ! get some of the grids ready
      NULLIFY (rho_r%pw, rho_g%pw, pot_g%pw, psi_L%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, psi_L%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! run the FFT once, to set up buffers and to take into account the memory
      rho_r%pw%cr3d = 0.0D0
      CALL pw_transfer(rho_r%pw, rho_g%pw)

      CALL timestop(handle)

   END SUBROUTINE prepare_gpw

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param e_cutoff_old ...
!> \param cutoff_old ...
!> \param relative_cutoff_old ...
!> \param pw_env_sub ...
!> \param task_list_sub ...
!> \param auxbas_pw_pool ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param psi_L ...
! **************************************************************************************************
   SUBROUTINE cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, pw_env_sub, &
                          task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: e_cutoff_old
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_old, relative_cutoff_old
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r, rho_g, pot_g, psi_L

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cleanup_gpw'

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      ! and now free the whole lot
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, psi_L%pw)

      CALL deallocate_task_list(task_list_sub)
      CALL pw_env_release(pw_env_sub)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      ! restore the initial value of the cutoff
      dft_control%qs_control%e_cutoff = e_cutoff_old
      dft_control%qs_control%cutoff = cutoff_old
      dft_control%qs_control%relative_cutoff = relative_cutoff_old

      CALL timestop(handle)

   END SUBROUTINE cleanup_gpw

! **************************************************************************************************
!> \brief ...
!> \param rho_r ...
!> \param rho_g ...
!> \param poisson_env ...
!> \param pot_g ...
!> \param potential_type ...
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_type, omega)
      TYPE(pw_p_type), INTENT(IN)                        :: rho_r, rho_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_p_type), INTENT(IN)                        :: pot_g
      INTEGER, INTENT(IN), OPTIONAL                      :: potential_type
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_potential_gpw'

      INTEGER                                            :: handle, my_potential_type

      CALL timeset(routineN, handle)

      my_potential_type = do_potential_coulomb
      IF (PRESENT(potential_type)) THEN
         my_potential_type = potential_type
         IF (my_potential_type /= do_potential_coulomb .AND. .NOT. PRESENT(omega)) THEN
            CPABORT("Need omega when longrange potential is requested.")
         END IF
      END IF

      ! in case we do Coulomb metric for RI, we need the Coulomb operator, but for RI with the
      ! overlap metric, we do not need the Coulomb operator
      IF (my_potential_type == do_potential_coulomb .OR. my_potential_type == do_potential_long) THEN
         CALL pw_transfer(rho_r%pw, rho_g%pw)
         CALL pw_poisson_solve(poisson_env, rho_g%pw, vhartree=pot_g%pw)
         IF (my_potential_type == do_potential_long) CALL pw_gauss_damp(pot_g%pw, omega)
         CALL pw_transfer(pot_g%pw, rho_r%pw)
      END IF

      CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)
      CALL timestop(handle)

   END SUBROUTINE calc_potential_gpw

END MODULE mp2_eri_gpw
