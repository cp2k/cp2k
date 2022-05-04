!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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
                                              do_potential_id,&
                                              do_potential_long,&
                                              do_potential_mix_cl,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_methods,                      ONLY: &
        pw_compl_gauss_damp, pw_copy, pw_derive, pw_gauss_damp, pw_gauss_damp_mix, pw_integral_ab, &
        pw_scale, pw_transfer, pw_truncated, pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                              calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
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

   PUBLIC :: mp2_eri_2c_integrate_gpw, mp2_eri_3c_integrate_gpw, calc_potential_gpw, cleanup_gpw, prepare_gpw, &
             integrate_potential_forces_2c, integrate_potential_forces_3c_1c, integrate_potential_forces_3c_2c, &
             virial_gpw_potential

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
!> \param potential_parameter ...
!> \param mat_munu ...
!> \param qs_env ...
!> \param task_list_sub ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_3c_integrate_gpw(mo_coeff, psi_L, rho_g, atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                       pw_env_sub, external_vector, poisson_env, rho_r, pot_g, &
                                       potential_parameter, mat_munu, qs_env, task_list_sub)

      TYPE(cp_fm_type), INTENT(IN), POINTER              :: mo_coeff
      TYPE(pw_p_type), INTENT(INOUT)                     :: psi_L, rho_g
      TYPE(atomic_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(pw_env_type), INTENT(IN), POINTER             :: pw_env_sub
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: external_vector
      TYPE(pw_poisson_type), INTENT(IN), POINTER         :: poisson_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r, pot_g
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(task_list_type), INTENT(IN), POINTER          :: task_list_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_eri_3c_integrate_gpw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! pseudo psi_L
      CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                  basis_type="RI_AUX", &
                                  external_vector=external_vector)

      CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_parameter)

      ! and finally (K|mu nu)
      CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
      CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, &
                              calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                              pw_env_external=pw_env_sub, task_list_external=task_list_sub)

      CALL timestop(handle)

   END SUBROUTINE mp2_eri_3c_integrate_gpw

! **************************************************************************************************
!> \brief Integrates the potential of an RI function
!> \param qs_env ...
!> \param para_env_sub ...
!> \param dimen_RI ...
!> \param mo_coeff ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param natom ...
!> \param potential_parameter ...
!> \param sab_orb_sub ...
!> \param L_local_col ...
!> \param kind_of ...
! **************************************************************************************************
   SUBROUTINE mp2_eri_2c_integrate_gpw(qs_env, para_env_sub, dimen_RI, mo_coeff, my_group_L_start, my_group_L_end, &
                                       natom, potential_parameter, sab_orb_sub, L_local_col, kind_of)

      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: mo_coeff
      INTEGER, INTENT(IN)                                :: my_group_L_start, my_group_L_end, natom
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_orb_sub
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: L_local_col
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_eri_2c_integrate_gpw'

      INTEGER :: dir, handle, handle2, i, i_counter, iatom, igrid_level, ikind, ipgf, iset, lb(3), &
         LLL, location(3), max_nseta, na1, na2, ncoa, nseta, offset, sgfa, tp(3), ub(3)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: offset_2d
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: map_it_here, use_subpatch
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

         CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_parameter)

         NULLIFY (rs_v)
         NULLIFY (rs_descs)
         CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
         DO i = 1, SIZE(rs_v)
            CALL rs_grid_retain(rs_v(i)%rs_grid)
         END DO
         CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

         CALL timestop(handle2)

         offset = 0
         ! Prepare offset ahead of OMP parallel loop
         CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxnset=max_nseta, basis_type="RI_AUX")
         ALLOCATE (offset_2d(max_nseta, natom))

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")
            nseta = basis_set_a%nset
            nsgfa => basis_set_a%nsgf_set
            DO iset = 1, nseta
               offset = offset + nsgfa(iset)
               offset_2d(iset, iatom) = offset
            END DO
         END DO

         ! integrate the little bastards
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(natom, particle_set, cell, pw_env_sub, rs_v, offset_2d, &
         !$OMP        qs_kind_set, ncoset, para_env_sub, dft_control, i_counter, &
         !$OMP        kind_of, l_local_col) &
         !$OMP PRIVATE(iatom, ikind, basis_set_a, first_sgfa, la_max, la_min, npgfa, &
         !$OMP         nseta, nsgfa, rpgfa, set_radius_a, sphi_a, zeta, &
         !$OMP         ra, iset, ncoa, I_tmp2, I_ab, igrid_level, &
         !$OMP         map_it_here, dir, tp, lb, ub, location, ipgf, &
         !$OMP         sgfa, na1, na2, radius, offset, use_subpatch)
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

               offset = offset_2d(iset, iatom)

               igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(zeta(:, iset)))
               map_it_here = .FALSE.
               use_subpatch = .FALSE.
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
                     END IF
                     ! distributed grid, only map if it is local to the grid
                     location(dir) = tp(dir) + rs_v(igrid_level)%rs_grid%desc%lb(dir)
                  END DO
                  IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                      lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                      lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                     map_it_here = .TRUE.
                  END IF
                  use_subpatch = .TRUE.
               ELSE
                  ! not distributed, just a round-robin distribution over the full set of CPUs
                  IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
               END IF

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
                        hab=I_tmp2, &
                        o1=na1 - 1, &
                        o2=0, &
                        radius=radius, &
                        calculate_forces=.FALSE., &
                        use_subpatch=use_subpatch, &
                        subpatch_pattern=0)
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
         !$OMP END PARALLEL DO
         DEALLOCATE (offset_2d)

         DO i = 1, SIZE(rs_v)
            CALL rs_grid_release(rs_v(i)%rs_grid)
         END DO

      END DO

      DEALLOCATE (wf_vector)

      CALL mp_sum(L_local_col, para_env_sub%group)

      CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                       task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Integrates the potential of a RI function obtaining the forces and stress tensor
!> \param rho_r ...
!> \param wf_vector ...
!> \param matrix ...
!> \param rho_g ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param pot_g ...
!> \param potential_parameter ...
!> \param use_virial ...
!> \param rho_g_copy ...
!> \param dvg ...
!> \param kind_of ...
!> \param atom_of_kind ...
!> \param G_PQ_local ...
!> \param force ...
!> \param h_stress ...
!> \param para_env_sub ...
!> \param dft_control ...
!> \param psi_L ...
! **************************************************************************************************
   SUBROUTINE integrate_potential_forces_2c(rho_r, wf_vector, matrix, rho_g, atomic_kind_set, &
                                            qs_kind_set, particle_set, cell, pw_env_sub, poisson_env, pot_g, &
                                            potential_parameter, use_virial, rho_g_copy, dvg, &
                                            kind_of, atom_of_kind, G_PQ_local, force, h_stress, para_env_sub, &
                                            dft_control, psi_L)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wf_vector
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: matrix
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_g
      TYPE(atomic_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(pw_env_type), INTENT(IN), POINTER             :: pw_env_sub
      TYPE(pw_poisson_type), INTENT(IN), POINTER         :: poisson_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: pot_g
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_g_copy, dvg(3)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of, atom_of_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: G_PQ_local
      TYPE(qs_force_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: force
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: h_stress
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      TYPE(pw_p_type), INTENT(INOUT)                     :: psi_L

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_potential_forces_2c'

      INTEGER                                            :: handle, handle2

      CALL timeset(routineN, handle)

      ! calculate potential associated to the single aux function
      CALL timeset(routineN//"_wf_pot", handle2)
      ! pseudo psi_L
      CALL pw_zero(rho_r%pw)
      CALL pw_zero(rho_g%pw)
      CALL calculate_wavefunction(matrix, 1, rho_r, rho_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, &
                                  pw_env_sub, basis_type='RI_AUX', &
                                  external_vector=wf_vector)
      IF (use_virial) THEN
         CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_parameter, dvg)
      ELSE
         CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_parameter)
      END IF
      CALL timestop(handle2)

      IF (use_virial) THEN
         ! make a copy of the density in G space
         CALL pw_copy(rho_g%pw, rho_g_copy%pw)

         ! add the volume contribution to the virial due to
         ! the (P|Q) integrals, first we put the full gamme_PQ
         ! pseudo wave-function into grid in order to calculate the
         ! hartree potential derivatives
         CALL pw_zero(psi_L%pw)
         CALL pw_zero(rho_g%pw)
         CALL calculate_wavefunction(matrix, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                     basis_type="RI_AUX", &
                                     external_vector=-2.0_dp*G_PQ_local)
         ! transfer to reciprocal space and calculate potential
         CALL calc_potential_gpw(psi_L, rho_g, poisson_env, pot_g, potential_parameter, no_transfer=.TRUE.)
         ! update virial with volume term (first calculate hartree like energy (diagonal part of the virial))
         CALL virial_gpw_potential(rho_g_copy, pot_g, rho_g, dvg, h_stress, potential_parameter, para_env_sub)
      END IF

      ! integrate the potential of the single gaussian and update
      ! 2-center forces with Gamma_PQ
      CALL integrate_potential(pw_env_sub, rho_r, kind_of, atom_of_kind, particle_set, qs_kind_set, &
                               G_PQ_local, cell, force, use_virial, h_stress, para_env_sub, dft_control)

      CALL timestop(handle)
   END SUBROUTINE integrate_potential_forces_2c

   ! If virial, then extract it from virial%pv_virial
! **************************************************************************************************
!> \brief Takes the precomputed potential of an RI wave-function and determines matrix element and
!>        gradients with product of Gaussians
!> \param mat_munu ...
!> \param rho_r ...
!> \param matrix_P_munu ...
!> \param qs_env ...
!> \param pw_env_sub ...
!> \param task_list_sub ...
! **************************************************************************************************
   SUBROUTINE integrate_potential_forces_3c_1c(mat_munu, rho_r, matrix_P_munu, qs_env, pw_env_sub, &
                                               task_list_sub)

      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r
      TYPE(dbcsr_p_type), INTENT(IN)                     :: matrix_P_munu
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(pw_env_type), INTENT(IN), POINTER             :: pw_env_sub
      TYPE(task_list_type), INTENT(INOUT), POINTER       :: task_list_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_potential_forces_3c_1c'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! integrate the potential of the single gaussian and update
      ! 3-center forces
      CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
      CALL integrate_v_rspace(rho_r, hmat=mat_munu, pmat=matrix_P_munu, &
                              qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.FALSE., gapw=.FALSE., &
                              pw_env_external=pw_env_sub, &
                              task_list_external=task_list_sub)

      CALL timestop(handle)
   END SUBROUTINE integrate_potential_forces_3c_1c

! **************************************************************************************************
!> \brief Integrates potential of two Gaussians to a potential
!> \param matrix_P_munu ...
!> \param rho_r ...
!> \param rho_g ...
!> \param task_list_sub ...
!> \param pw_env_sub ...
!> \param potential_parameter ...
!> \param ks_env ...
!> \param poisson_env ...
!> \param pot_g ...
!> \param use_virial ...
!> \param rho_g_copy ...
!> \param dvg ...
!> \param h_stress ...
!> \param para_env_sub ...
!> \param kind_of ...
!> \param atom_of_kind ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param LLL ...
!> \param force ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE integrate_potential_forces_3c_2c(matrix_P_munu, rho_r, rho_g, task_list_sub, pw_env_sub, &
                                               potential_parameter, &
                                               ks_env, poisson_env, pot_g, use_virial, rho_g_copy, dvg, &
                                               h_stress, para_env_sub, kind_of, atom_of_kind, &
                                               qs_kind_set, particle_set, cell, LLL, force, dft_control)

      TYPE(dbcsr_p_type), INTENT(IN)                     :: matrix_P_munu
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_r, rho_g
      TYPE(task_list_type), INTENT(IN), POINTER          :: task_list_sub
      TYPE(pw_env_type), INTENT(IN), POINTER             :: pw_env_sub
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(qs_ks_env_type), INTENT(IN), POINTER          :: ks_env
      TYPE(pw_poisson_type), INTENT(IN), POINTER         :: poisson_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: pot_g
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_g_copy, dvg(3)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: h_stress
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_sub
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of, atom_of_kind
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      INTEGER, INTENT(IN)                                :: LLL
      TYPE(qs_force_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: force
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_potential_forces_3c_2c'

      INTEGER :: atom_a, dir, handle, handle2, i, iatom, igrid_level, ikind, iorb, ipgf, iset, &
         lb(3), location(3), na1, na2, ncoa, nseta, offset, sgfa, tp(3), ub(3)
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: map_it_here, skip_shell, use_subpatch
      REAL(KIND=dp)                                      :: radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: I_ab
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: I_tmp2, pab, rpgfa, sphi_a, zeta
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v

      CALL timeset(routineN, handle)

      ! put the gamma density on grid
      CALL timeset(routineN//"_Gpot", handle2)
      CALL pw_zero(rho_r%pw)
      CALL pw_zero(rho_g%pw)
      CALL calculate_rho_elec(matrix_p=matrix_P_munu%matrix, &
                              rho=rho_r, &
                              rho_gspace=rho_g, &
                              task_list_external=task_list_sub, &
                              pw_env_external=pw_env_sub, &
                              ks_env=ks_env)
      ! calculate associated hartree potential
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_parameter)
      CALL timestop(handle2)

      IF (use_virial) CALL virial_gpw_potential(rho_g_copy, pot_g, rho_g, dvg, h_stress, potential_parameter, para_env_sub)

      ! integrate potential with auxiliary basis function derivatives
      NULLIFY (rs_v)
      NULLIFY (rs_descs)
      CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO
      CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

      offset = 0
      DO iatom = 1, SIZE(kind_of)
         ikind = kind_of(iatom)
         atom_a = atom_of_kind(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, &
                          basis_type="RI_AUX")

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

         force_a(:) = 0.0_dp
         force_b(:) = 0.0_dp
         IF (use_virial) THEN
            my_virial_a = 0.0_dp
            my_virial_b = 0.0_dp
         END IF

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            ALLOCATE (I_tmp2(ncoa, 1))
            I_tmp2 = 0.0_dp
            ALLOCATE (I_ab(nsgfa(iset), 1))
            I_ab = 0.0_dp
            ALLOCATE (pab(ncoa, 1))
            pab = 0.0_dp

            skip_shell = .TRUE.
            DO iorb = 1, nsgfa(iset)
               IF (iorb + offset == LLL) THEN
                  I_ab(iorb, 1) = 1.0_dp
                  skip_shell = .FALSE.
               END IF
            END DO

            IF (skip_shell) THEN
               offset = offset + nsgfa(iset)
               DEALLOCATE (I_tmp2)
               DEALLOCATE (I_ab)
               DEALLOCATE (pab)
               CYCLE
            END IF

            CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), &
                       1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                       I_ab(1, 1), SIZE(I_ab, 1), &
                       0.0_dp, pab(1, 1), SIZE(pab, 1))
            DEALLOCATE (I_ab)

            igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(zeta(:, iset)))
            map_it_here = .FALSE.
            use_subpatch = .FALSE.

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
                  END IF
                  ! distributed grid, only map if it is local to the grid
                  location(dir) = tp(dir) + rs_v(igrid_level)%rs_grid%desc%lb(dir)
               END DO
               IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                   lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                   lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                  map_it_here = .TRUE.
               END IF
               use_subpatch = .TRUE.
            ELSE
               ! not distributed, just a round-robin distribution over the full set of CPUs
               IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
            END IF

            offset = offset + nsgfa(iset)

            IF (map_it_here) THEN
               DO ipgf = 1, npgfa(iset)
                  na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                  na2 = ipgf*ncoset(la_max(iset))
                  igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, zeta(ipgf, iset))

                  radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                    lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                    zetp=zeta(ipgf, iset), &
                                                    eps=dft_control%qs_control%eps_gvg_rspace, &
                                                    prefactor=1.0_dp, cutoff=1.0_dp)

                  CALL integrate_pgf_product(la_max=la_max(iset), zeta=zeta(ipgf, iset)/2.0_dp, la_min=la_min(iset), &
                                             lb_max=0, zetb=zeta(ipgf, iset)/2.0_dp, lb_min=0, &
                                             ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                             rsgrid=rs_v(igrid_level)%rs_grid, &
                                             hab=I_tmp2, &
                                             pab=pab, &
                                             o1=na1 - 1, &
                                             o2=0, &
                                             radius=radius, &
                                             calculate_forces=.TRUE., &
                                             force_a=force_a, force_b=force_b, &
                                             use_virial=use_virial, &
                                             my_virial_a=my_virial_a, &
                                             my_virial_b=my_virial_b, &
                                             use_subpatch=use_subpatch, &
                                             subpatch_pattern=0)
               END DO
            END IF

            DEALLOCATE (I_tmp2)
            DEALLOCATE (pab)

         END DO

         force(ikind)%rho_elec(:, atom_a) = &
            force(ikind)%rho_elec(:, atom_a) + force_a(:) + force_b(:)
         IF (use_virial) THEN
            h_stress = h_stress + my_virial_a + my_virial_b
         END IF
      END DO

      DO i = 1, SIZE(rs_v)
         CALL rs_grid_release(rs_v(i)%rs_grid)
      END DO

      CALL timestop(handle)

   END SUBROUTINE integrate_potential_forces_3c_2c

! **************************************************************************************************
!> \brief Calculates stress tensor contribution from the operator
!> \param rho_g_copy ...
!> \param pot_g ...
!> \param rho_g ...
!> \param dvg ...
!> \param h_stress ...
!> \param potential_parameter ...
!> \param para_env_sub ...
! **************************************************************************************************
   SUBROUTINE virial_gpw_potential(rho_g_copy, pot_g, rho_g, dvg, h_stress, potential_parameter, para_env_sub)

      TYPE(pw_p_type)                                    :: rho_g_copy, pot_g, rho_g, dvg(3)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: h_stress
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_gpw_potential'

      INTEGER                                            :: alpha, beta, handle
      INTEGER, DIMENSION(3)                              :: comp
      REAL(KIND=dp)                                      :: e_hartree

      ! add the volume contribution to the virial coming from
      ! the 3-center integrals (mu nu|P)
      CALL timeset(routineN, handle)
      e_hartree = 0.0_dp
      e_hartree = pw_integral_ab(rho_g_copy%pw, pot_g%pw)

      DO alpha = 1, 3
         comp = 0
         comp(alpha) = 1
         CALL pw_copy(pot_g%pw, rho_g%pw)
         CALL pw_derive(rho_g%pw, comp)
         CALL factor_virial_gpw(rho_g%pw, potential_parameter)
         h_stress(alpha, alpha) = h_stress(alpha, alpha) - e_hartree/REAL(para_env_sub%num_pe, dp)
         DO beta = alpha, 3
            h_stress(alpha, beta) = h_stress(alpha, beta) &
                                    - 2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw)/fourpi/REAL(para_env_sub%num_pe, dp)
            h_stress(beta, alpha) = h_stress(alpha, beta)
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE virial_gpw_potential

! **************************************************************************************************
!> \brief Multiplies a potential in g space with the factor ln(V(g)/Vc(g))' with Vc(g) being the
!>        Fourier-transformed of the Coulomg potential
!> \param pw ...
!> \param potential_parameter parameters of potential V(g)
! **************************************************************************************************
   SUBROUTINE factor_virial_gpw(pw, potential_parameter)
      TYPE(pw_type), INTENT(IN)                          :: pw
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: omega_2, rchalf, scale_coul, scale_long, &
                                                            tmp

      IF (.NOT. (pw%in_space == RECIPROCALSPACE .AND. pw%in_use == COMPLEXDATA1D)) &
         CPABORT("pw in wrong space or wrong data type")
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_coulomb)
         RETURN
      CASE (do_potential_long)
         omega_2 = 1.0_dp/(2.0_dp*potential_parameter%omega)**2
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(pw,potential_parameter,omega_2)
         pw%cc(:) = pw%cc(:)*(1.0_dp + omega_2*pw%pw_grid%gsq(:))
!$OMP END PARALLEL WORKSHARE
      CASE (do_potential_short)
         omega_2 = 1.0_dp/(2.0_dp*potential_parameter%omega)**2
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2)
         DO i = 1, SIZE(pw%cc)
            tmp = omega_2*pw%pw_grid%gsq(i)
            ! For too small arguments, use the Taylor polynomial to prevent division by zero
            IF (ABS(tmp) >= 0.003_dp) THEN
               pw%cc(i) = pw%cc(i)*(1.0_dp - tmp*EXP(-tmp)/(1.0_dp - EXP(-tmp)))
            ELSE
               pw%cc(i) = pw%cc(i)*(0.5_dp*tmp - tmp**2/12.0_dp)
            END IF
         END DO
!$OMP END PARALLEL DO
      CASE (do_potential_mix_cl)
         omega_2 = 1.0_dp/(2.0_dp*potential_parameter%omega)**2
         scale_coul = potential_parameter%scale_coulomb
         scale_long = potential_parameter%scale_longrange
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,omega_2,scale_long,scale_coul)
         DO i = 1, SIZE(pw%cc)
            tmp = omega_2*pw%pw_grid%gsq(i)
            pw%cc(i) = pw%cc(i)*(1.0_dp + scale_long*tmp*EXP(-tmp)/(scale_coul + scale_long*EXP(-tmp)))
         END DO
!$OMP END PARALLEL DO
      CASE (do_potential_truncated)
         rchalf = 0.5_dp*potential_parameter%cutoff_radius
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,tmp) &
!$OMP             SHARED(pw,rchalf)
         DO i = 1, SIZE(pw%cc)
            tmp = rchalf*SQRT(pw%pw_grid%gsq(i))
            ! For too small arguments, use the Taylor polynomial to prevent division by zero
            IF (ABS(tmp) >= 0.0001_dp) THEN
               pw%cc(i) = pw%cc(i)*(1.0_dp - tmp/TAN(tmp))
            ELSE
               pw%cc(i) = pw%cc(i)*tmp**2*(1.0_dp + tmp**2/15.0_dp)/3.0_dp
            END IF
         END DO
!$OMP END PARALLEL DO
      CASE DEFAULT
         CPABORT("Unknown potential type")
      END SELECT

   END SUBROUTINE factor_virial_gpw

! **************************************************************************************************
!> \brief Integrate potential of RI function with RI function
!> \param pw_env_sub ...
!> \param pot_r ...
!> \param kind_of ...
!> \param atom_of_kind ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param G_PQ_local ...
!> \param cell ...
!> \param force ...
!> \param use_virial ...
!> \param h_stress ...
!> \param para_env_sub ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE integrate_potential(pw_env_sub, pot_r, kind_of, atom_of_kind, particle_set, qs_kind_set, &
                                  G_PQ_local, cell, force, use_virial, h_stress, para_env_sub, dft_control)

      TYPE(pw_env_type), INTENT(IN), POINTER             :: pw_env_sub
      TYPE(pw_p_type), INTENT(IN)                        :: pot_r
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of, atom_of_kind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: G_PQ_local
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      LOGICAL, INTENT(IN)                                :: use_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: h_stress
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_potential'

      INTEGER                                            :: atom_a, dir, handle, i, iatom, &
                                                            igrid_level, ikind, ipgf, iset, lb(3), &
                                                            location(3), na1, na2, ncoa, nseta, &
                                                            offset, sgfa, tp(3), ub(3)
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: map_it_here, use_subpatch
      REAL(KIND=dp)                                      :: radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: I_ab
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: I_tmp2, pab, rpgfa, sphi_a, zeta
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v

      CALL timeset(routineN, handle)

      NULLIFY (rs_v)
      NULLIFY (rs_descs)
      CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO
      CALL potential_pw2rs(rs_v, pot_r, pw_env_sub)

      offset = 0
      DO iatom = 1, SIZE(kind_of)
         ikind = kind_of(iatom)
         atom_a = atom_of_kind(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, &
                          basis_type="RI_AUX")

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

         force_a(:) = 0.0_dp
         force_b(:) = 0.0_dp
         IF (use_virial) THEN
            my_virial_a = 0.0_dp
            my_virial_b = 0.0_dp
         END IF

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            ALLOCATE (I_tmp2(ncoa, 1))
            I_tmp2 = 0.0_dp
            ALLOCATE (I_ab(nsgfa(iset), 1))
            I_ab = 0.0_dp
            ALLOCATE (pab(ncoa, 1))
            pab = 0.0_dp

            I_ab(1:nsgfa(iset), 1) = -4.0_dp*G_PQ_local(offset + 1:offset + nsgfa(iset))

            CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), &
                       1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                       I_ab(1, 1), SIZE(I_ab, 1), &
                       0.0_dp, pab(1, 1), SIZE(pab, 1))

            I_ab = 0.0_dp

            igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(zeta(:, iset)))
            map_it_here = .FALSE.
            use_subpatch = .FALSE.
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
                  END IF
                  ! distributed grid, only map if it is local to the grid
                  location(dir) = tp(dir) + rs_v(igrid_level)%rs_grid%desc%lb(dir)
               END DO
               IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                   lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                   lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                  map_it_here = .TRUE.
               END IF
               use_subpatch = .TRUE.
            ELSE
               ! not distributed, just a round-robin distribution over the full set of CPUs
               IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
            END IF

            offset = offset + nsgfa(iset)

            IF (map_it_here) THEN
               DO ipgf = 1, npgfa(iset)
                  na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                  na2 = ipgf*ncoset(la_max(iset))
                  igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, zeta(ipgf, iset))

                  radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                    lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                    zetp=zeta(ipgf, iset), &
                                                    eps=dft_control%qs_control%eps_gvg_rspace, &
                                                    prefactor=1.0_dp, cutoff=1.0_dp)

                  CALL integrate_pgf_product(la_max=la_max(iset), zeta=zeta(ipgf, iset)/2.0_dp, la_min=la_min(iset), &
                                             lb_max=0, zetb=zeta(ipgf, iset)/2.0_dp, lb_min=0, &
                                             ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                             rsgrid=rs_v(igrid_level)%rs_grid, &
                                             hab=I_tmp2, &
                                             pab=pab, &
                                             o1=na1 - 1, &
                                             o2=0, &
                                             radius=radius, &
                                             calculate_forces=.TRUE., &
                                             force_a=force_a, force_b=force_b, &
                                             use_virial=use_virial, &
                                             my_virial_a=my_virial_a, &
                                             my_virial_b=my_virial_b, &
                                             use_subpatch=use_subpatch, &
                                             subpatch_pattern=0)

               END DO

            END IF

            DEALLOCATE (I_tmp2)
            DEALLOCATE (I_ab)
            DEALLOCATE (pab)

         END DO

         force(ikind)%rho_elec(:, atom_a) = &
            force(ikind)%rho_elec(:, atom_a) + force_a(:) + force_b
         IF (use_virial) THEN
            h_stress = h_stress + my_virial_a + my_virial_b
         END IF
      END DO

      DO i = 1, SIZE(rs_v)
         CALL rs_grid_release(rs_v(i)%rs_grid)
      END DO
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Prepares GPW calculation for RI-MP2/RI-RPA
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
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: e_cutoff_old
      REAL(KIND=dp), INTENT(OUT)                         :: cutoff_old, relative_cutoff_old
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: auxbas_pw_pool
      TYPE(pw_poisson_type), INTENT(IN), POINTER         :: poisson_env
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
!> \brief Cleanup GPW integration for RI-MP2/RI-RPA
!> \param qs_env ...
!> \param e_cutoff_old ...
!> \param cutoff_old ...
!> \param relative_cutoff_old ...
!> \param para_env_sub ...
!> \param pw_env_sub ...
!> \param task_list_sub ...
!> \param auxbas_pw_pool ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param psi_L ...
! **************************************************************************************************
   SUBROUTINE cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                          task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: e_cutoff_old
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_old, relative_cutoff_old
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: auxbas_pw_pool
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
      CALL pw_env_release(pw_env_sub, para_env=para_env_sub)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      ! restore the initial value of the cutoff
      dft_control%qs_control%e_cutoff = e_cutoff_old
      dft_control%qs_control%cutoff = cutoff_old
      dft_control%qs_control%relative_cutoff = relative_cutoff_old

      CALL timestop(handle)

   END SUBROUTINE cleanup_gpw

! **************************************************************************************************
!> \brief Calculates potential from a given density in g-space
!> \param pot_r on output: potential in r space
!> \param rho_g on input: rho in g space
!> \param poisson_env ...
!> \param pot_g on output: potential in g space
!> \param potential_parameter Potential parameters, if not provided, assume Coulomb potential
!> \param dvg in output: first derivatives of the corresponding Coulomb potential
!> \param no_transfer whether NOT to transform potential from g-space to r-space (default: do it)
! **************************************************************************************************
   SUBROUTINE calc_potential_gpw(pot_r, rho_g, poisson_env, pot_g, potential_parameter, dvg, no_transfer)
      TYPE(pw_p_type), INTENT(IN)                        :: pot_r, rho_g
      TYPE(pw_poisson_type), INTENT(IN), POINTER         :: poisson_env
      TYPE(pw_p_type), INTENT(IN)                        :: pot_g
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: potential_parameter
      TYPE(pw_p_type), DIMENSION(3), INTENT(IN), &
         OPTIONAL                                        :: dvg
      LOGICAL, INTENT(IN), OPTIONAL                      :: no_transfer

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_potential_gpw'

      INTEGER                                            :: handle, my_potential_type
      LOGICAL                                            :: my_transfer

      CALL timeset(routineN, handle)

      my_potential_type = do_potential_coulomb
      IF (PRESENT(potential_parameter)) THEN
         my_potential_type = potential_parameter%potential_type
      END IF

      my_transfer = .TRUE.
      IF (PRESENT(no_transfer)) my_transfer = .NOT. no_transfer

      ! in case we do Coulomb metric for RI, we need the Coulomb operator, but for RI with the
      ! overlap metric, we do not need the Coulomb operator
      IF (my_potential_type /= do_potential_id) THEN
         CALL pw_poisson_solve(poisson_env, rho_g%pw, vhartree=pot_g%pw, dvhartree=dvg)
         IF (my_potential_type == do_potential_long) CALL pw_gauss_damp(pot_g%pw, potential_parameter%omega)
         IF (my_potential_type == do_potential_short) CALL pw_compl_gauss_damp(pot_g%pw, potential_parameter%omega)
         IF (my_potential_type == do_potential_truncated) CALL pw_truncated(pot_g%pw, potential_parameter%cutoff_radius)
         IF (my_potential_type == do_potential_mix_cl) CALL pw_gauss_damp_mix(pot_g%pw, potential_parameter%omega, &
                                                                              potential_parameter%scale_coulomb, &
                                                                              potential_parameter%scale_longrange)
         IF (my_transfer) CALL pw_transfer(pot_g%pw, pot_r%pw)
      ELSE
         ! If we use an overlap metric, make sure to use the correct potential=density on output
         IF (my_transfer) CALL pw_transfer(rho_g%pw, pot_r%pw)
      END IF

      IF (my_transfer) CALL pw_scale(pot_r%pw, pot_r%pw%pw_grid%dvol)
      CALL timestop(handle)

   END SUBROUTINE calc_potential_gpw

END MODULE mp2_eri_gpw
