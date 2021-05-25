!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate gradients of RI-GPW-MP2 energy using pw
!> \par History
!>      10.2013 created [Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_ri_grad
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_param
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_indxg2l, cp_fm_indxg2p, cp_fm_indxl2g, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_all, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_split
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_copy_into_existing, dbcsr_create, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_transposed, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: &
        mp_allgather, mp_alltoall, mp_irecv, mp_isend, mp_request_null, mp_sendrecv, mp_sum, &
        mp_wait, mp_waitall
   USE mp2_eri,                         ONLY: mp2_eri_2c_integrate,&
                                              mp2_eri_3c_integrate,&
                                              mp2_eri_deallocate_forces,&
                                              mp2_eri_force
   USE mp2_eri_gpw,                     ONLY: calc_potential_gpw,&
                                              cleanup_gpw,&
                                              prepare_gpw
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              integ_mat_buffer_type_2D,&
                                              mp2_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                              calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
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
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_grad'

   PUBLIC :: calc_ri_mp2_nonsep

CONTAINS

! **************************************************************************************************
!> \brief Calculate the non-separable part of the gradients and update the
!>        Lagrangian
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param Eigenval ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param sab_orb_sub ...
!> \param mat_munu ...
!> \param blacs_env_sub ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE calc_ri_mp2_nonsep(qs_env, mp2_env, para_env, para_env_sub, cell, particle_set, &
                                 atomic_kind_set, qs_kind_set, mo_coeff, nmo, homo, dimen_RI, Eigenval, &
                                 my_group_L_start, my_group_L_end, my_group_L_size, sab_orb_sub, mat_munu, &
                                 blacs_env_sub)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_coeff
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: my_group_L_start, my_group_L_end, &
                                                            my_group_L_size
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ri_mp2_nonsep'

      INTEGER :: alpha, atom_a, beta, dimen, dir, eri_method, handle, handle2, handle3, i, iatom, &
         igrid_level, ikind, iorb, ipgf, iset, ispin, itmp(2), L_counter, lb(3), LLL, location(3), &
         my_P_end, my_P_size, my_P_start, na1, na2, natom, ncoa, nseta, nspins, offset, &
         potential_type, sgfa, tp(3), ub(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, virtual
      INTEGER, DIMENSION(3)                              :: comp
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: alpha_beta, map_it_here, skip_shell, &
                                                            use_virial
      REAL(KIND=dp)                                      :: cutoff_old, e_hartree, eps_filter, &
                                                            factor, omega, pair_energy, radius, &
                                                            relative_cutoff_old, total_rho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old, wf_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G_PQ_local, I_ab
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_stress, my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: I_tmp2, pab, rpgfa, sphi_a, zeta
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: L1_mu_i, L2_nu_a
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(dbcsr_p_type)                                 :: matrix_P_munu, matrix_P_munu_nosym
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: Lag_mu_i_1, Lag_nu_a_2, matrix_P_inu, &
                                                            mo_coeff_o, mo_coeff_v
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:, :)   :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_munu_local, matrix_P_munu_local
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(mp2_eri_force), ALLOCATABLE, DIMENSION(:)     :: force_2c, force_3c_aux, force_3c_orb_mu, &
                                                            force_3c_orb_nu
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: dvg(3), pot_g, psi_L, rho_g, rho_r, &
                                                            temp_pw_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      eri_method = mp2_env%eri_method
      eri_param => mp2_env%eri_mme_param

      ! Find out whether we have a closed or open shell
      nspins = SIZE(homo)
      alpha_beta = (nspins == 2)

      dimen = nmo
      ALLOCATE (virtual(nspins))
      virtual(:) = dimen - homo(:)
      potential_type = mp2_env%potential_parameter%potential_type
      omega = mp2_env%potential_parameter%omega
      eps_filter = mp2_env%mp2_gpw%eps_filter
      ALLOCATE (mo_coeff_o(nspins), mo_coeff_v(nspins), G_P_ia(nspins, my_group_L_size))
      DO ispin = 1, nspins
         mo_coeff_o(ispin)%matrix => mp2_env%ri_grad%mo_coeff_o(ispin)%matrix
         mo_coeff_v(ispin)%matrix => mp2_env%ri_grad%mo_coeff_v(ispin)%matrix
         DO LLL = 1, my_group_L_size
            G_P_ia(ispin, LLL)%matrix => mp2_env%ri_grad%G_P_ia(LLL, ispin)%matrix
         END DO
      ENDDO
      DEALLOCATE (mp2_env%ri_grad%G_P_ia)

      itmp = get_limit(dimen_RI, para_env_sub%num_pe, para_env_sub%mepos)
      my_P_start = itmp(1)
      my_P_end = itmp(2)
      my_P_size = itmp(2) - itmp(1) + 1
      ALLOCATE (G_PQ_local(dimen_RI, my_group_L_size))

      G_PQ_local = 0.0_dp
      IF (.NOT. alpha_beta) THEN
         G_PQ_local(my_P_start:my_P_end, 1:my_group_L_size) = mp2_env%ri_grad%Gamma_PQ(1)%array
      ELSE
         G_PQ_local(my_P_start:my_P_end, 1:my_group_L_size) = &
            0.50_dp*(mp2_env%ri_grad%Gamma_PQ(1)%array + mp2_env%ri_grad%Gamma_PQ(2)%array)
      ENDIF
      DEALLOCATE (mp2_env%ri_grad%Gamma_PQ(1)%array)
      IF (alpha_beta) THEN
         DEALLOCATE (mp2_env%ri_grad%Gamma_PQ(2)%array)
      ENDIF
      CALL mp_sum(G_PQ_local, para_env_sub%group)

      ! deallocate here PQ_half, maybe useful in the future
      ! This is really bad style
      DEALLOCATE (mp2_env%ri_grad%PQ_half)

      ! create matrix holding the back transformation (G_P_inu)
      ALLOCATE (matrix_P_inu(nspins))
      DO ispin = 1, nspins
         ALLOCATE (matrix_P_inu(ispin)%matrix)
         CALL dbcsr_create(matrix_P_inu(ispin)%matrix, template=mo_coeff_o(ispin)%matrix)
      ENDDO

      ! non symmetric matrix
      ALLOCATE (matrix_P_munu_nosym%matrix)
      CALL dbcsr_create(matrix_P_munu_nosym%matrix, template=mat_munu%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! create Lagrangian matrices in mixed AO/MO formalism
      ALLOCATE (Lag_mu_i_1(nspins))
      DO ispin = 1, nspins
         ALLOCATE (Lag_mu_i_1(ispin)%matrix)
         CALL dbcsr_create(Lag_mu_i_1(ispin)%matrix, template=mo_coeff_o(ispin)%matrix)
         CALL dbcsr_set(Lag_mu_i_1(ispin)%matrix, 0.0_dp)
      ENDDO

      ALLOCATE (Lag_nu_a_2(nspins))
      DO ispin = 1, nspins
         ALLOCATE (Lag_nu_a_2(ispin)%matrix)
         CALL dbcsr_create(Lag_nu_a_2(ispin)%matrix, template=mo_coeff_v(ispin)%matrix)
         CALL dbcsr_set(Lag_nu_a_2(ispin)%matrix, 0.0_dp)
      ENDDO

      ! get forces
      NULLIFY (force, virial)
      CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)

      ! prepare integral derivatives with mme method
      IF (eri_method .EQ. do_eri_mme) THEN
         ALLOCATE (matrix_P_munu_local(my_group_L_size))
         ALLOCATE (mat_munu_local(my_group_L_size))
         L_counter = 0
         DO LLL = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1
            ALLOCATE (matrix_P_munu_local(L_counter)%matrix)
            ALLOCATE (mat_munu_local(L_counter)%matrix)
            CALL dbcsr_create(matrix_P_munu_local(L_counter)%matrix, template=mat_munu%matrix, &
                              matrix_type=dbcsr_type_symmetric)
            CALL dbcsr_create(mat_munu_local(L_counter)%matrix, template=mat_munu%matrix, &
                              matrix_type=dbcsr_type_symmetric)
            CALL dbcsr_copy(mat_munu_local(L_counter)%matrix, mat_munu%matrix)
            CALL dbcsr_set(mat_munu_local(L_counter)%matrix, 0.0_dp)

            CALL G_P_transform_MO_to_AO(matrix_P_munu_local(L_counter), matrix_P_munu_nosym, mat_munu, &
                                        G_P_ia(:, L_counter), matrix_P_inu, &
                                        mo_coeff_v, mo_coeff_o, eps_filter)
         ENDDO

         ALLOCATE (I_tmp2(dimen_RI, my_group_L_size))
         I_tmp2(:, :) = 0.0_dp
         CALL mp2_eri_2c_integrate(eri_param, potential_type, omega, para_env_sub, qs_env, &
                                   basis_type_a="RI_AUX", basis_type_b="RI_AUX", &
                                   hab=I_tmp2, first_b=my_group_L_start, last_b=my_group_L_end, &
                                   eri_method=eri_method, pab=G_PQ_local, force_a=force_2c)

         DEALLOCATE (I_tmp2)
         CALL mp2_eri_3c_integrate(eri_param, potential_type, omega, para_env_sub, qs_env, &
                                   first_c=my_group_L_start, last_c=my_group_L_end, mat_ab=mat_munu_local, &
                                   basis_type_a="ORB", basis_type_b="ORB", basis_type_c="RI_AUX", &
                                   sab_nl=sab_orb_sub, eri_method=eri_method, &
                                   pabc=matrix_P_munu_local, &
                                   force_a=force_3c_orb_mu, force_b=force_3c_orb_nu, force_c=force_3c_aux)

         L_counter = 0
         DO LLL = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1
            ! we recompute matrix_P_inu
            DO ispin = 1, nspins
               CALL dbcsr_multiply("N", "T", 1.0_dp, mo_coeff_v(ispin)%matrix, G_P_ia(ispin, L_counter)%matrix, &
                                   0.0_dp, matrix_P_inu(ispin)%matrix, filter_eps=eps_filter)
            ENDDO

            ! The matrices of G_P_ia are deallocated here
            CALL update_lagrangian(mat_munu_local(L_counter), matrix_P_inu, Lag_mu_i_1, &
                                   G_P_ia(:, L_counter), mo_coeff_o, Lag_nu_a_2, &
                                   eps_filter)
         ENDDO

         DO ikind = 1, SIZE(force)
            force(ikind)%mp2_non_sep(:, :) = -4.0_dp*force_2c(ikind)%forces(:, :) + &
                                             force_3c_orb_mu(ikind)%forces(:, :) + &
                                             force_3c_orb_nu(ikind)%forces(:, :) + &
                                             force_3c_aux(ikind)%forces(:, :)
         END DO

         CALL mp2_eri_deallocate_forces(force_2c)
         CALL mp2_eri_deallocate_forces(force_3c_aux)
         CALL mp2_eri_deallocate_forces(force_3c_orb_mu)
         CALL mp2_eri_deallocate_forces(force_3c_orb_nu)
         CALL dbcsr_deallocate_matrix_set(matrix_P_munu_local)
         CALL dbcsr_deallocate_matrix_set(mat_munu_local)

      ELSEIF (eri_method == do_eri_gpw) THEN
         ALLOCATE (matrix_P_munu%matrix)
         CALL dbcsr_create(matrix_P_munu%matrix, template=mat_munu%matrix, &
                           matrix_type=dbcsr_type_symmetric)

         CALL get_qs_env(qs_env, ks_env=ks_env)

         natom = SIZE(particle_set)

         ALLOCATE (kind_of(natom))
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

         ! Supporting stuff for GPW
         CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                          auxbas_pw_pool, poisson_env, task_list_sub, rho_r, rho_g, pot_g, psi_L, sab_orb_sub)

         ! wave function vector and supporting stuff
         ALLOCATE (wf_vector(dimen_RI))

         ! check if we want to calculate the virial
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         ! in case virial is required we need auxiliary pw
         ! for calculate the MP2-volume contribution to the virial
         ! (hartree potential derivatives)
         IF (use_virial) THEN
            NULLIFY (temp_pw_g%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, temp_pw_g%pw, &
                                   use_data=COMPLEXDATA1D, &
                                   in_space=RECIPROCALSPACE)
            DO i = 1, 3
               NULLIFY (dvg(i)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, dvg(i)%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE)
            END DO
         END IF

         ! start main loop over auxiliary basis functions
         CALL timeset(routineN//"_loop", handle2)

         L_counter = 0
         DO LLL = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1

            CALL G_P_transform_MO_to_AO(matrix_P_munu, matrix_P_munu_nosym, mat_munu, &
                                        G_P_ia(:, L_counter), matrix_P_inu, &
                                        mo_coeff_v, mo_coeff_o, eps_filter)

            ! calculate potential associted to the single aux function
            CALL timeset(routineN//"_wf_pot", handle3)
            wf_vector = 0.0_dp
            wf_vector(LLL) = 1.0_dp
            ! pseudo psi_L
            CALL pw_zero(rho_r%pw)
            DO ispin = 1, nspins
               CALL calculate_wavefunction(mo_coeff(ispin)%matrix, 1, psi_L, rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, &
                                           pw_env_sub, basis_type='RI_AUX', &
                                           external_vector=wf_vector)
               rho_r%pw%cr3d = rho_r%pw%cr3d+psi_L%pw%cr3d
            ENDDO
            IF (alpha_beta) rho_r%pw%cr3d = rho_r%pw%cr3d*0.5_dp

            CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, potential_type, omega)
            CALL timestop(handle3)

            IF (use_virial) THEN
               ! make a copy of the density in G space
               ! calculate the potential derivatives in G space
               CALL timeset(routineN//"_Virial", handle3)
               CALL pw_copy(rho_g%pw, temp_pw_g%pw)
               DO i = 1, 3
                  comp = 0
                  comp(i) = 1
                  CALL pw_copy(pot_g%pw, dvg(i)%pw)
                  CALL pw_derive(dvg(i)%pw, comp)
               END DO
               CALL timestop(handle3)
            END IF

            ! integrate the potential of the single gaussian and update
            ! 2-center forces with Gamma_PQ
            CALL timeset(routineN//"_int_PQ", handle3)
            NULLIFY (rs_v)
            NULLIFY (rs_descs)
            CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
            DO i = 1, SIZE(rs_v)
               CALL rs_grid_retain(rs_v(i)%rs_grid)
            END DO
            CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

            offset = 0
            DO iatom = 1, natom
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

                  I_ab(1:nsgfa(iset), 1) = -4.0_dp*G_PQ_local(offset + 1:offset + nsgfa(iset), L_counter)

                  CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             I_ab(1, 1), SIZE(I_ab, 1), &
                             0.0_dp, pab(1, 1), SIZE(pab, 1))

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
                                                   cell=cell, &
                                                   cube_info=pw_env_sub%cube_info(igrid_level), &
                                                   hab=I_tmp2, &
                                                   pab=pab, &
                                                   o1=na1 - 1, &
                                                   o2=0, &
                                                   radius=radius, &
                                                   calculate_forces=.TRUE., &
                                                   force_a=force_a, force_b=force_b, &
                                                   use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b)

                     END DO

                  END IF

                  DEALLOCATE (I_tmp2)
                  DEALLOCATE (I_ab)
                  DEALLOCATE (pab)

               END DO

               force(ikind)%rho_elec(:, atom_a) = &
                  force(ikind)%rho_elec(:, atom_a) + force_a(:) + force_b
               IF (use_virial) THEN
                  virial%pv_mp2 = virial%pv_mp2 + my_virial_a + my_virial_b
                  virial%pv_virial = virial%pv_virial + my_virial_a + my_virial_b
               END IF
            END DO

            DO i = 1, SIZE(rs_v)
               CALL rs_grid_release(rs_v(i)%rs_grid)
            END DO
            CALL timestop(handle3)
            ! here we are done with the 2 centers

            ! integrate the potential of the single gaussian and update
            ! 3-center forces
            CALL timeset(routineN//"_int", handle3)
            CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
            CALL integrate_v_rspace(rho_r, hmat=mat_munu, pmat=matrix_P_munu, &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.FALSE., gapw=.FALSE., &
                                    pw_env_external=pw_env_sub, &
                                    task_list_external=task_list_sub)
            CALL timestop(handle3)

            ! The matrices of G_P_ia are deallocated here
            CALL update_lagrangian(mat_munu, matrix_P_inu, Lag_mu_i_1, &
                                   G_P_ia(:, L_counter), mo_coeff_o, Lag_nu_a_2, &
                                   eps_filter)

            IF (use_virial) THEN
               ! add the volume contribution to the virial due to
               ! the (P|Q) integrals, first we put the full gamme_PQ
               ! pseudo wave-function into grid in order to calculate the
               ! hartree potential derivatives
               CALL timeset(routineN//"_Virial", handle3)
               wf_vector = 0.0_dp
               wf_vector(:) = -2.0_dp*G_PQ_local(:, L_counter)
               CALL calculate_wavefunction(mo_coeff(1)%matrix, 1, psi_L, rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                           basis_type="RI_AUX", &
                                           external_vector=wf_vector)
               ! transfer to reciprocal space and calculate potential
               rho_r%pw%cr3d = psi_L%pw%cr3d
               CALL pw_transfer(rho_r%pw, rho_g%pw)
               CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
               ! update virial with volume term (first calculate hartree like energy (diagonal part of the virial))
               e_hartree = 0.0_dp
               h_stress = 0.0_dp
               e_hartree = pw_integral_ab(temp_pw_g%pw, pot_g%pw)
               DO alpha = 1, 3
                  comp = 0
                  comp(alpha) = 1
                  CALL pw_copy(pot_g%pw, rho_g%pw)
                  CALL pw_derive(rho_g%pw, comp)
                  h_stress(alpha, alpha) = -e_hartree
                  DO beta = alpha, 3
                     h_stress(alpha, beta) = h_stress(alpha, beta) &
                                             - 2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw)/fourpi
                     h_stress(beta, alpha) = h_stress(alpha, beta)
                  END DO
               END DO
               virial%pv_mp2 = virial%pv_mp2 + h_stress/REAL(para_env_sub%num_pe, dp)
               virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env_sub%num_pe, dp)
               CALL timestop(handle3)
            END IF

            ! put the gamma density on grid
            CALL timeset(routineN//"_Gpot", handle3)
            CALL calculate_rho_elec(matrix_p=matrix_P_munu%matrix, &
                                    rho=rho_r, &
                                    rho_gspace=rho_g, &
                                    total_rho=total_rho, &
                                    task_list_external=task_list_sub, &
                                    pw_env_external=pw_env_sub, &
                                    ks_env=ks_env)
            ! calculate associated hartree potential
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
            CALL pw_transfer(pot_g%pw, rho_r%pw)
            CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)
            CALL timestop(handle3)

            IF (use_virial) THEN
               ! add the volume contribution to the virial coming from
               ! the 3-center integrals (mu nu|P)
               CALL timeset(routineN//"_Virial", handle3)
               e_hartree = 0.0_dp
               h_stress = 0.0_dp
               e_hartree = pw_integral_ab(temp_pw_g%pw, pot_g%pw)
               DO alpha = 1, 3
                  comp = 0
                  comp(alpha) = 1
                  CALL pw_copy(pot_g%pw, rho_g%pw)
                  CALL pw_derive(rho_g%pw, comp)
                  h_stress(alpha, alpha) = -e_hartree
                  DO beta = alpha, 3
                     h_stress(alpha, beta) = h_stress(alpha, beta) &
                                             - 2.0_dp*pw_integral_ab(rho_g%pw, dvg(beta)%pw)/fourpi
                     h_stress(beta, alpha) = h_stress(alpha, beta)
                  END DO
               END DO
               virial%pv_mp2 = virial%pv_mp2 + h_stress/REAL(para_env_sub%num_pe, dp)
               virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env_sub%num_pe, dp)
               CALL timestop(handle3)
            END IF

            ! integrate potential with auxiliary basis function derivatives
            NULLIFY (rs_v)
            NULLIFY (rs_descs)
            CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
            DO i = 1, SIZE(rs_v)
               CALL rs_grid_retain(rs_v(i)%rs_grid)
            END DO
            CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

            offset = 0
            DO iatom = 1, natom
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
                                                   cell=cell, &
                                                   cube_info=pw_env_sub%cube_info(igrid_level), &
                                                   hab=I_tmp2, &
                                                   pab=pab, &
                                                   o1=na1 - 1, &
                                                   o2=0, &
                                                   radius=radius, &
                                                   calculate_forces=.TRUE., &
                                                   force_a=force_a, force_b=force_b, &
                                                   use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b)
                     END DO
                  END IF

                  DEALLOCATE (I_tmp2)
                  DEALLOCATE (I_ab)
                  DEALLOCATE (pab)

               END DO

               force(ikind)%rho_elec(:, atom_a) = &
                  force(ikind)%rho_elec(:, atom_a) + force_a(:) + force_b(:)
               IF (use_virial) THEN
                  virial%pv_mp2 = virial%pv_mp2 + my_virial_a + my_virial_b
                  virial%pv_virial = virial%pv_virial + my_virial_a + my_virial_b
               END IF
            END DO

            DO i = 1, SIZE(rs_v)
               CALL rs_grid_release(rs_v(i)%rs_grid)
            END DO

         END DO

         CALL timestop(handle2)

         DEALLOCATE (kind_of)
         DEALLOCATE (atom_of_kind)
         DEALLOCATE (wf_vector)

         IF (use_virial) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, temp_pw_g%pw)
            DO i = 1, 3
               CALL pw_pool_give_back_pw(auxbas_pw_pool, dvg(i)%pw)
            END DO
         END IF

         CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, pw_env_sub, &
                          task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)

         CALL dbcsr_release(matrix_P_munu%matrix)
         DEALLOCATE (matrix_P_munu%matrix)

      ENDIF

      DEALLOCATE (G_PQ_local)

      CALL dbcsr_release(matrix_P_munu_nosym%matrix)
      DEALLOCATE (matrix_P_munu_nosym%matrix)

      DO ispin = 1, nspins
         CALL dbcsr_release(matrix_P_inu(ispin)%matrix)
         DEALLOCATE (matrix_P_inu(ispin)%matrix)
      ENDDO
      DEALLOCATE (matrix_P_inu, G_P_ia)

      ! move the forces in the correct place
      IF (eri_method .EQ. do_eri_gpw) THEN
         DO ikind = 1, SIZE(force)
            force(ikind)%mp2_non_sep(:, :) = force(ikind)%rho_elec(:, :)
            force(ikind)%rho_elec(:, :) = 0.0_dp
         END DO
      ENDIF

      ! Now we move from the local matrices to the global ones
      ! defined over all MPI tasks
      ! Start with moving from the DBCSR to FM for the lagrangians

      ALLOCATE (L1_mu_i(nspins), L2_nu_a(nspins))
      DO ispin = 1, nspins
         ! Now we move from the local matrices to the global ones
         ! defined over all MPI tasks
         ! Start with moving from the DBCSR to FM for the lagrangians
         NULLIFY (L1_mu_i(ispin)%matrix, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env_sub, context=blacs_env_sub, &
                                  nrow_global=dimen, ncol_global=homo(ispin))
         CALL cp_fm_create(L1_mu_i(ispin)%matrix, fm_struct_tmp, name="Lag_mu_i")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(L1_mu_i(ispin)%matrix, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=Lag_mu_i_1(ispin)%matrix, fm=L1_mu_i(ispin)%matrix)

         ! release Lag_mu_i_1
         CALL dbcsr_release(Lag_mu_i_1(ispin)%matrix)
         DEALLOCATE (Lag_mu_i_1(ispin)%matrix)

         NULLIFY (L2_nu_a(ispin)%matrix, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env_sub, context=blacs_env_sub, &
                                  nrow_global=dimen, ncol_global=virtual(ispin))
         CALL cp_fm_create(L2_nu_a(ispin)%matrix, fm_struct_tmp, name="Lag_nu_a")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(L2_nu_a(ispin)%matrix, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=Lag_nu_a_2(ispin)%matrix, fm=L2_nu_a(ispin)%matrix)

         ! release Lag_nu_a_2
         CALL dbcsr_release(Lag_nu_a_2(ispin)%matrix)
         DEALLOCATE (Lag_nu_a_2(ispin)%matrix)
      ENDDO
      DEALLOCATE (Lag_mu_i_1, Lag_nu_a_2)

      ! Set the factor to multiply P_ij (depends on the open or closed shell)
      factor = 1.0_dp
      IF (alpha_beta) factor = 0.50_dp

      DO ispin = 1, nspins
         CALL create_W_P(qs_env, mp2_env, mo_coeff(ispin)%matrix, homo(ispin), virtual(ispin), dimen, para_env, &
                         para_env_sub, Eigenval(:, ispin), L1_mu_i(ispin)%matrix, L2_nu_a(ispin)%matrix, &
                         factor, ispin /= 1)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE calc_ri_mp2_nonsep

! **************************************************************************************************
!> \brief ...
!> \param G_P_munu ...
!> \param G_P_munu_nosym ...
!> \param mat_munu ...
!> \param G_P_ia ...
!> \param G_P_inu ...
!> \param mo_coeff_v ...
!> \param mo_coeff_o ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE G_P_transform_MO_to_AO(G_P_munu, G_P_munu_nosym, mat_munu, G_P_ia, G_P_inu, &
                                     mo_coeff_v, mo_coeff_o, eps_filter)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: G_P_munu, G_P_munu_nosym, mat_munu
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: G_P_inu
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff_v, mo_coeff_o
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'G_P_transform_MO_to_AO'

      INTEGER                                            :: handle

      CALL G_P_transform_alpha_beta(G_P_ia, G_P_inu, G_P_munu_nosym, mo_coeff_v, mo_coeff_o, eps_filter)

      ! symmetrize
      CALL timeset(routineN//"_symmetrize", handle)
      CALL dbcsr_set(G_P_munu%matrix, 0.0_dp)
      CALL dbcsr_transposed(G_P_munu%matrix, G_P_munu_nosym%matrix)
      CALL dbcsr_add(G_P_munu%matrix, G_P_munu_nosym%matrix, &
                     alpha_scalar=2.0_dp, beta_scalar=2.0_dp)
      ! this is a trick to avoid that integrate_v_rspace starts to cry
      CALL dbcsr_copy_into_existing(mat_munu%matrix, G_P_munu%matrix)
      CALL dbcsr_copy(G_P_munu%matrix, mat_munu%matrix)

      CALL timestop(handle)

   END SUBROUTINE G_P_transform_MO_to_AO

! **************************************************************************************************
!> \brief ...
!> \param G_P_ia ...
!> \param G_P_inu ...
!> \param G_P_munu ...
!> \param mo_coeff_v ...
!> \param mo_coeff_o ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE G_P_transform_alpha_beta(G_P_ia, G_P_inu, G_P_munu, mo_coeff_v, mo_coeff_o, eps_filter)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: G_P_inu
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: G_P_munu
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff_v, mo_coeff_o
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'G_P_transform_alpha_beta'

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: factor

      CALL timeset(routineN, handle)

      factor = 1.0_dp/REAL(SIZE(G_P_ia), dp)

      CALL dbcsr_set(G_P_munu%matrix, 0.0_dp)

      DO ispin = 1, SIZE(G_P_ia)
         ! first back-transformation a->nu
         CALL dbcsr_multiply("N", "T", 1.0_dp, mo_coeff_v(ispin)%matrix, G_P_ia(ispin)%matrix, &
                             0.0_dp, G_P_inu(ispin)%matrix, filter_eps=eps_filter)

         ! second back-transformation i->mu
         CALL dbcsr_multiply("N", "T", factor, G_P_inu(ispin)%matrix, mo_coeff_o(ispin)%matrix, &
                             1.0_dp, G_P_munu%matrix, filter_eps=eps_filter)
      END DO

      CALL timestop(handle)

   END SUBROUTINE G_P_transform_alpha_beta

! **************************************************************************************************
!> \brief ...
!> \param mat_munu ...
!> \param matrix_P_inu ...
!> \param Lag_mu_i_1 ...
!> \param G_P_ia ...
!> \param mo_coeff_o ...
!> \param Lag_nu_a_2 ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE update_lagrangian(mat_munu, matrix_P_inu, Lag_mu_i_1, &
                                G_P_ia, mo_coeff_o, Lag_nu_a_2, &
                                eps_filter)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: matrix_P_inu, Lag_mu_i_1, G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: mo_coeff_o
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: Lag_nu_a_2
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_lagrangian'

      INTEGER                                            :: handle, ispin

      ! update lagrangian
      CALL timeset(routineN, handle)

      DO ispin = 1, SIZE(G_P_ia)
         ! first contract mat_munu with the half back transformed Gamma_i_nu
         ! in order to update Lag_mu_i_1
         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, matrix_P_inu(ispin)%matrix, &
                             1.0_dp, Lag_mu_i_1(ispin)%matrix, filter_eps=eps_filter)

         ! transform first index of mat_munu and store the result into matrix_P_inu
         CALL dbcsr_set(matrix_P_inu(ispin)%matrix, 0.0_dp)
         CALL dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o(ispin)%matrix, &
                             0.0_dp, matrix_P_inu(ispin)%matrix, filter_eps=eps_filter)

         ! contract the transformend matrix_P_inu with the untransformend Gamma_i_a
         ! in order to update Lag_nu_a_2
         CALL dbcsr_multiply("N", "N", -1.0_dp, matrix_P_inu(ispin)%matrix, G_P_ia(ispin)%matrix, &
                             1.0_dp, Lag_nu_a_2(ispin)%matrix, filter_eps=eps_filter)

         ! release the actual gamma_P_ia
         CALL dbcsr_release(G_P_ia(ispin)%matrix)
         DEALLOCATE (G_P_ia(ispin)%matrix)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE update_lagrangian

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param mo_coeff ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param Eigenval ...
!> \param L1_mu_i ...
!> \param L2_nu_a ...
!> \param factor ...
!> \param alpha_beta ...
! **************************************************************************************************
   SUBROUTINE create_W_P(qs_env, mp2_env, mo_coeff, homo, virtual, dimen, para_env, para_env_sub, &
                         Eigenval, L1_mu_i, L2_nu_a, factor, alpha_beta)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      TYPE(cp_fm_type), POINTER                          :: L1_mu_i, L2_nu_a
      REAL(KIND=dp), INTENT(IN)                          :: factor
      LOGICAL, INTENT(IN)                                :: alpha_beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_W_P'

      INTEGER :: color_exchange, dummy_proc, handle, handle2, handle3, i, i_global, i_local, iiB, &
         iii, iproc, ispin, itmp(2), j_global, j_local, jjB, max_col_size, max_row_size, &
         my_B_virtual_end, my_B_virtual_start, my_spin, mypcol, myprow, ncol_block, ncol_block_1i, &
         ncol_block_2a, ncol_local, ncol_local_1i, ncol_local_2a, npcol, npcol_1i, npcol_2a, &
         nprow, nprow_1i, nprow_2a, nrow_block, nrow_block_1i, nrow_block_2a, nrow_local, &
         nrow_local_1i, nrow_local_2a, number_of_rec, number_of_send, proc_receive, &
         proc_receive_static, proc_send, proc_send_ex, proc_send_static, proc_send_sub, proc_shift
      INTEGER :: rec_col_size, rec_counter, rec_row_size, send_col_size, send_counter, send_pcol, &
         send_prow, send_row_size, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iii_vet, map_rec_size, map_send_size, pos_info, &
         pos_info_ex, proc_2_send_pos, proc_map_ex, req_send, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: grid_2_mepos, mepos_2_grid, my_col_indeces_info_1i, &
         my_col_indeces_info_2a, my_row_indeces_info_1i, my_row_indeces_info_2a, sizes, sizes_1i, &
         sizes_2a
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: col_indeces_info_1i, &
                                                            col_indeces_info_2a, &
                                                            row_indeces_info_1i, &
                                                            row_indeces_info_2a
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_1i, &
                                                            col_indices_2a, row_indices, &
                                                            row_indices_1i, row_indices_2a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ab_rec, ab_send, mat_rec, mat_send
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: fm_P_ij, L_mu_q
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(integ_mat_buffer_type_2D), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_cyclic

      CALL timeset(routineN, handle)

      my_spin = 1
      IF (alpha_beta) my_spin = 2

      ! create the globally distributed mixed lagrangian
      NULLIFY (blacs_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env)

      NULLIFY (L_mu_q, fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(L_mu_q, fm_struct_tmp, name="Lag_mu_q")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(L_mu_q, 0.0_dp)

      ! create all information array
      ALLOCATE (pos_info(0:para_env%num_pe - 1))
      CALL mp_allgather(para_env_sub%mepos, pos_info, para_env%group)

      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
      sub_proc_map = 0
      DO i = 0, para_env_sub%num_pe - 1
         sub_proc_map(i) = i
         sub_proc_map(-i - 1) = para_env_sub%num_pe - i - 1
         sub_proc_map(para_env_sub%num_pe + i) = i
      END DO

      ! get matrix information for the global
      CALL cp_fm_get_info(matrix=L_mu_q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          nrow_block=nrow_block, &
                          ncol_block=ncol_block)
      myprow = L_mu_q%matrix_struct%context%mepos(1)
      mypcol = L_mu_q%matrix_struct%context%mepos(2)
      nprow = L_mu_q%matrix_struct%context%num_pe(1)
      npcol = L_mu_q%matrix_struct%context%num_pe(2)

      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      grid_2_mepos(myprow, mypcol) = para_env%mepos
      CALL mp_sum(grid_2_mepos, para_env%group)

      ! get matrix information for L1_mu_i
      CALL cp_fm_get_info(matrix=L1_mu_i, &
                          nrow_local=nrow_local_1i, &
                          ncol_local=ncol_local_1i, &
                          row_indices=row_indices_1i, &
                          col_indices=col_indices_1i, &
                          nrow_block=nrow_block_1i, &
                          ncol_block=ncol_block_1i)
      nprow_1i = L1_mu_i%matrix_struct%context%num_pe(1)
      npcol_1i = L1_mu_i%matrix_struct%context%num_pe(2)

      ALLOCATE (sizes_1i(2, 0:para_env_sub%num_pe - 1))
      CALL mp_allgather([nrow_local_1i, ncol_local_1i], sizes_1i, para_env_sub%group)

      ! get matrix information for L2_nu_a
      CALL cp_fm_get_info(matrix=L2_nu_a, &
                          nrow_local=nrow_local_2a, &
                          ncol_local=ncol_local_2a, &
                          row_indices=row_indices_2a, &
                          col_indices=col_indices_2a, &
                          nrow_block=nrow_block_2a, &
                          ncol_block=ncol_block_2a)
      nprow_2a = L2_nu_a%matrix_struct%context%num_pe(1)
      npcol_2a = L2_nu_a%matrix_struct%context%num_pe(2)

      ALLOCATE (sizes_2a(2, 0:para_env_sub%num_pe - 1))
      CALL mp_allgather([nrow_local_2a, ncol_local_2a], sizes_2a, para_env_sub%group)

      ! Here we perform a ring communication scheme taking into account
      ! for the sub-group distribution of the source matrices.
      ! as a first step we need to redistribute the data within
      ! the subgroup.
      ! In order to do so we have to allocate the structure
      ! that will hold the local data involved in communication, this
      ! structure will be the same for processes in different subgroups
      ! sharing the same position in the subgroup.
      ! -1) create the exchange para_env
      color_exchange = para_env_sub%mepos
      CALL cp_para_env_split(para_env_exchange, para_env, color_exchange)
      ! crate the proc maps exchange and info
      ALLOCATE (proc_map_ex(-para_env_exchange%num_pe:2*para_env_exchange%num_pe - 1))
      DO i = 0, para_env_exchange%num_pe - 1
         proc_map_ex(i) = i
         proc_map_ex(-i - 1) = para_env_exchange%num_pe - i - 1
         proc_map_ex(para_env_exchange%num_pe + i) = i
      END DO
      ALLOCATE (pos_info_ex(0:para_env%num_pe - 1))
      CALL mp_allgather(para_env_exchange%mepos, pos_info_ex, para_env%group)
      ALLOCATE (sizes(2, 0:para_env_exchange%num_pe - 1))
      CALL mp_allgather([nrow_local, ncol_local], sizes, para_env_exchange%group)

      ! 0) store some info about indeces of the fm matrices (subgroup)
      CALL timeset(routineN//"_inx", handle2)
      ! matrix L1_mu_i
      max_row_size = MAXVAL(sizes_1i(1, :))
      max_col_size = MAXVAL(sizes_1i(2, :))
      ALLOCATE (row_indeces_info_1i(2, max_row_size, 0:para_env_sub%num_pe - 1))
      ALLOCATE (col_indeces_info_1i(2, max_col_size, 0:para_env_sub%num_pe - 1))
      ALLOCATE (my_row_indeces_info_1i(2, max_row_size))
      ALLOCATE (my_col_indeces_info_1i(2, max_col_size))
      row_indeces_info_1i = 0
      col_indeces_info_1i = 0
      dummy_proc = 0
      ! row
      DO iiB = 1, nrow_local_1i
         i_global = row_indices_1i(iiB)
         send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                   L_mu_q%matrix_struct%first_p_pos(1), nprow)
         i_local = cp_fm_indxg2l(i_global, nrow_block, dummy_proc, &
                                 L_mu_q%matrix_struct%first_p_pos(1), nprow)
         my_row_indeces_info_1i(1, iiB) = send_prow
         my_row_indeces_info_1i(2, iiB) = i_local
      END DO
      ! col
      DO jjB = 1, ncol_local_1i
         j_global = col_indices_1i(jjB)
         send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                   L_mu_q%matrix_struct%first_p_pos(2), npcol)
         j_local = cp_fm_indxg2l(j_global, ncol_block, dummy_proc, &
                                 L_mu_q%matrix_struct%first_p_pos(2), npcol)
         my_col_indeces_info_1i(1, jjB) = send_pcol
         my_col_indeces_info_1i(2, jjB) = j_local
      END DO
      CALL mp_allgather(my_row_indeces_info_1i, row_indeces_info_1i, para_env_sub%group)
      CALL mp_allgather(my_col_indeces_info_1i, col_indeces_info_1i, para_env_sub%group)
      DEALLOCATE (my_row_indeces_info_1i, my_col_indeces_info_1i)

      ! matrix L2_nu_a
      max_row_size = MAXVAL(sizes_2a(1, :))
      max_col_size = MAXVAL(sizes_2a(2, :))
      ALLOCATE (row_indeces_info_2a(2, max_row_size, 0:para_env_sub%num_pe - 1))
      ALLOCATE (col_indeces_info_2a(2, max_col_size, 0:para_env_sub%num_pe - 1))
      ALLOCATE (my_row_indeces_info_2a(2, max_row_size))
      ALLOCATE (my_col_indeces_info_2a(2, max_col_size))
      row_indeces_info_2a = 0
      col_indeces_info_2a = 0
      ! row
      DO iiB = 1, nrow_local_2a
         i_global = row_indices_2a(iiB)
         send_prow = cp_fm_indxg2p(i_global, nrow_block, dummy_proc, &
                                   L_mu_q%matrix_struct%first_p_pos(1), nprow)
         i_local = cp_fm_indxg2l(i_global, nrow_block, dummy_proc, &
                                 L_mu_q%matrix_struct%first_p_pos(1), nprow)
         my_row_indeces_info_2a(1, iiB) = send_prow
         my_row_indeces_info_2a(2, iiB) = i_local
      END DO
      ! col
      DO jjB = 1, ncol_local_2a
         j_global = col_indices_2a(jjB) + homo
         send_pcol = cp_fm_indxg2p(j_global, ncol_block, dummy_proc, &
                                   L_mu_q%matrix_struct%first_p_pos(2), npcol)
         j_local = cp_fm_indxg2l(j_global, ncol_block, dummy_proc, &
                                 L_mu_q%matrix_struct%first_p_pos(2), npcol)
         my_col_indeces_info_2a(1, jjB) = send_pcol
         my_col_indeces_info_2a(2, jjB) = j_local
      END DO
      CALL mp_allgather(my_row_indeces_info_2a, row_indeces_info_2a, para_env_sub%group)
      CALL mp_allgather(my_col_indeces_info_2a, col_indeces_info_2a, para_env_sub%group)
      DEALLOCATE (my_row_indeces_info_2a, my_col_indeces_info_2a)
      CALL timestop(handle2)

      ! 1) define the map for sending data in the subgroup starting with L1_mu_i
      CALL timeset(routineN//"_subinfo", handle2)
      ALLOCATE (map_send_size(0:para_env_sub%num_pe - 1))
      map_send_size = 0
      DO jjB = 1, ncol_local_1i
         send_pcol = col_indeces_info_1i(1, jjB, para_env_sub%mepos)
         DO iiB = 1, nrow_local_1i
            send_prow = row_indeces_info_1i(1, iiB, para_env_sub%mepos)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            proc_send_sub = pos_info(proc_send)
            map_send_size(proc_send_sub) = map_send_size(proc_send_sub) + 1
         END DO
      END DO
      ! and the same for L2_nu_a
      DO jjB = 1, ncol_local_2a
         send_pcol = col_indeces_info_2a(1, jjB, para_env_sub%mepos)
         DO iiB = 1, nrow_local_2a
            send_prow = row_indeces_info_2a(1, iiB, para_env_sub%mepos)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            proc_send_sub = pos_info(proc_send)
            map_send_size(proc_send_sub) = map_send_size(proc_send_sub) + 1
         END DO
      END DO
      ! and exchange data in order to create map_rec_size
      ALLOCATE (map_rec_size(0:para_env_sub%num_pe - 1))
      map_rec_size = 0
      CALL mp_alltoall(map_send_size, map_rec_size, 1, para_env_sub%group)
      CALL timestop(handle2)

      ! 2) reorder data in sending buffer
      CALL timeset(routineN//"_sub_Bsend", handle2)
      ! count the number of messages (include myself)
      number_of_send = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      send_counter = 0
      ALLOCATE (proc_2_send_pos(0:para_env_sub%num_pe - 1))
      proc_2_send_pos = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%msg = 0.0_dp
            buffer_send(send_counter)%proc = proc_send
            proc_2_send_pos(proc_send) = send_counter
         END IF
      END DO
      ! loop over the locally held data and fill the buffer_send
      ! for doing that we need an array that keep track if the
      ! sequential increase of the index for each message
      ALLOCATE (iii_vet(number_of_send))
      iii_vet = 0
      DO jjB = 1, ncol_local_1i
         send_pcol = col_indeces_info_1i(1, jjB, para_env_sub%mepos)
         DO iiB = 1, nrow_local_1i
            send_prow = row_indeces_info_1i(1, iiB, para_env_sub%mepos)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            proc_send_sub = pos_info(proc_send)
            send_counter = proc_2_send_pos(proc_send_sub)
            iii_vet(send_counter) = iii_vet(send_counter) + 1
            iii = iii_vet(send_counter)
            buffer_send(send_counter)%msg(iii) = L1_mu_i%local_data(iiB, jjB)
         END DO
      END DO
      ! release the local data of L1_mu_i
      DEALLOCATE (L1_mu_i%local_data)
      ! and the same for L2_nu_a
      DO jjB = 1, ncol_local_2a
         send_pcol = col_indeces_info_2a(1, jjB, para_env_sub%mepos)
         DO iiB = 1, nrow_local_2a
            send_prow = row_indeces_info_2a(1, iiB, para_env_sub%mepos)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            proc_send_sub = pos_info(proc_send)
            send_counter = proc_2_send_pos(proc_send_sub)
            iii_vet(send_counter) = iii_vet(send_counter) + 1
            iii = iii_vet(send_counter)
            buffer_send(send_counter)%msg(iii) = L2_nu_a%local_data(iiB, jjB)
         END DO
      END DO
      DEALLOCATE (L2_nu_a%local_data)
      DEALLOCATE (proc_2_send_pos)
      DEALLOCATE (iii_vet)
      CALL timestop(handle2)

      ! 3) create the buffer for receive, post the message with irecv
      !    and send the messages with mp_isend
      CALL timeset(routineN//"_sub_isendrecv", handle2)
      ! count the number of messages to be received
      number_of_rec = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO
      ALLOCATE (buffer_rec(number_of_rec))
      rec_counter = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%msg = 0.0_dp
            buffer_rec(rec_counter)%proc = proc_receive
            ! post the message to be received (not need to send to myself)
            IF (proc_receive /= para_env_sub%mepos) THEN
               CALL mp_irecv(buffer_rec(rec_counter)%msg, proc_receive, para_env_sub%group, &
                             buffer_rec(rec_counter)%msg_req)
            END IF
         END IF
      END DO
      ! send messages
      ALLOCATE (req_send(number_of_send))
      req_send = mp_request_null
      send_counter = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            IF (proc_send == para_env_sub%mepos) THEN
               buffer_rec(send_counter)%msg(:) = buffer_send(send_counter)%msg
            ELSE
               CALL mp_isend(buffer_send(send_counter)%msg, proc_send, para_env_sub%group, &
                             buffer_send(send_counter)%msg_req)
               req_send(send_counter) = buffer_send(send_counter)%msg_req
            END IF
         END IF
      END DO
      DEALLOCATE (map_send_size)
      CALL timestop(handle2)

      ! 4) (if memory is a problem we should move this part after point 5)
      !    Here we create the new buffer for cyclic(ring) communication and
      !    we fill it with the data received from the other member of the
      !    subgroup
      CALL timeset(routineN//"_Bcyclic", handle2)
      ! first allocata new structure
      ALLOCATE (buffer_cyclic(0:para_env_exchange%num_pe - 1))
      DO iproc = 0, para_env_exchange%num_pe - 1
         rec_row_size = sizes(1, iproc)
         rec_col_size = sizes(2, iproc)
         ALLOCATE (buffer_cyclic(iproc)%msg(rec_row_size, rec_col_size))
         buffer_cyclic(iproc)%msg = 0.0_dp
      END DO
      ! now collect data from other member of the subgroup and fill
      ! buffer_cyclic
      rec_counter = 0
      DO proc_shift = 0, para_env_sub%num_pe - 1
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1

            ! wait for the message
            IF (proc_receive /= para_env_sub%mepos) CALL mp_wait(buffer_rec(rec_counter)%msg_req)

            CALL timeset(routineN//"_fill", handle3)
            iii = 0
            DO jjB = 1, sizes_1i(2, proc_receive)
               send_pcol = col_indeces_info_1i(1, jjB, proc_receive)
               j_local = col_indeces_info_1i(2, jjB, proc_receive)
               DO iiB = 1, sizes_1i(1, proc_receive)
                  send_prow = row_indeces_info_1i(1, iiB, proc_receive)
                  proc_send = grid_2_mepos(send_prow, send_pcol)
                  proc_send_sub = pos_info(proc_send)
                  IF (proc_send_sub /= para_env_sub%mepos) CYCLE
                  iii = iii + 1
                  i_local = row_indeces_info_1i(2, iiB, proc_receive)
                  proc_send_ex = pos_info_ex(proc_send)
                  buffer_cyclic(proc_send_ex)%msg(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
               END DO
            END DO
            ! and the same for L2_nu_a
            DO jjB = 1, sizes_2a(2, proc_receive)
               send_pcol = col_indeces_info_2a(1, jjB, proc_receive)
               j_local = col_indeces_info_2a(2, jjB, proc_receive)
               DO iiB = 1, sizes_2a(1, proc_receive)
                  send_prow = row_indeces_info_2a(1, iiB, proc_receive)
                  proc_send = grid_2_mepos(send_prow, send_pcol)
                  proc_send_sub = pos_info(proc_send)
                  IF (proc_send_sub /= para_env_sub%mepos) CYCLE
                  iii = iii + 1
                  i_local = row_indeces_info_2a(2, iiB, proc_receive)
                  proc_send_ex = pos_info_ex(proc_send)
                  buffer_cyclic(proc_send_ex)%msg(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
               END DO
            END DO
            CALL timestop(handle3)

            ! deallocate the received message
            DEALLOCATE (buffer_rec(rec_counter)%msg)
         END IF
      END DO
      DEALLOCATE (row_indeces_info_1i)
      DEALLOCATE (col_indeces_info_1i)
      DEALLOCATE (row_indeces_info_2a)
      DEALLOCATE (col_indeces_info_2a)
      DEALLOCATE (buffer_rec)
      DEALLOCATE (map_rec_size)
      CALL timestop(handle2)

      ! 5)  Wait for all messeges to be sent in the subgroup
      CALL timeset(routineN//"_sub_waitall", handle2)
      CALL mp_waitall(req_send(:))
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE (buffer_send)
      DEALLOCATE (req_send)
      CALL timestop(handle2)

      ! 6) Start with ring communication
      CALL timeset(routineN//"_ring", handle2)
      proc_send_static = proc_map_ex(para_env_exchange%mepos + 1)
      proc_receive_static = proc_map_ex(para_env_exchange%mepos - 1)
      max_row_size = MAXVAL(sizes(1, :))
      max_col_size = MAXVAL(sizes(2, :))
      ALLOCATE (mat_send(max_row_size, max_col_size))
      ALLOCATE (mat_rec(max_row_size, max_col_size))
      mat_send = 0.0_dp
      mat_send(1:nrow_local, 1:ncol_local) = buffer_cyclic(para_env_exchange%mepos)%msg(:, :)
      DEALLOCATE (buffer_cyclic(para_env_exchange%mepos)%msg)
      DO proc_shift = 1, para_env_exchange%num_pe - 1
         proc_send = proc_map_ex(para_env_exchange%mepos + proc_shift)
         proc_receive = proc_map_ex(para_env_exchange%mepos - proc_shift)

         rec_row_size = sizes(1, proc_receive)
         rec_col_size = sizes(2, proc_receive)

         mat_rec = 0.0_dp
         CALL mp_sendrecv(mat_send, proc_send_static, &
                          mat_rec, proc_receive_static, &
                          para_env_exchange%group)

         mat_send = 0.0_dp
         mat_send(1:rec_row_size, 1:rec_col_size) = mat_rec(1:rec_row_size, 1:rec_col_size) + &
                                                    buffer_cyclic(proc_receive)%msg(:, :)

         DEALLOCATE (buffer_cyclic(proc_receive)%msg)
      END DO
      ! and finally
      CALL mp_sendrecv(mat_send, proc_send_static, &
                       mat_rec, proc_receive_static, &
                       para_env_exchange%group)
      L_mu_q%local_data(1:nrow_local, 1:ncol_local) = mat_rec(1:nrow_local, 1:ncol_local)
      DEALLOCATE (buffer_cyclic)
      DEALLOCATE (mat_send)
      DEALLOCATE (mat_rec)
      CALL timestop(handle2)

      DEALLOCATE (proc_map_ex)
      ! release para_env_exchange
      CALL cp_para_env_release(para_env_exchange)

      CALL cp_fm_release(L1_mu_i)
      CALL cp_fm_release(L2_nu_a)
      DEALLOCATE (pos_info_ex)
      DEALLOCATE (grid_2_mepos)
      DEALLOCATE (sizes)
      DEALLOCATE (sizes_1i)
      DEALLOCATE (sizes_2a)

      ! update the P_ij block of P_MP2 with the
      ! non-singular ij pairs
      CALL timeset(routineN//"_Pij", handle2)
      NULLIFY (fm_P_ij, fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(fm_P_ij, fm_struct_tmp, name="fm_P_ij")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(fm_P_ij, 0.0_dp)

      CALL cp_gemm('T', 'N', homo, homo, dimen, 1.0_dp, &
                   mo_coeff, L_mu_q, 0.0_dp, fm_P_ij, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)
      ! don't know if it is better to transpose (for communication reasons)
      ! or just recompute the transposed matrix
      CALL cp_gemm('T', 'N', homo, homo, dimen, -2.0_dp, &
                   L_mu_q, mo_coeff, 2.0_dp, fm_P_ij, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)
      ! we have it, update P_ij local
      CALL cp_fm_get_info(matrix=fm_P_ij, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            ! diagonal elements alreasy updated
            IF (j_global == i_global) CYCLE
            ! check if the given element is above the threshold
            IF (ABS(Eigenval(j_global) - Eigenval(i_global)) < mp2_env%ri_mp2%eps_canonical) CYCLE
            mp2_env%ri_grad%P_ij(my_spin)%array(i_global, j_global) = &
               factor*fm_P_ij%local_data(iiB, jjB)/(Eigenval(j_global) - Eigenval(i_global))
         END DO
      END DO
      ! release fm_P_ij
      CALL cp_fm_release(fm_P_ij)
      ! mp_sum it (we can avoid mp_sum, but for now let's keep it easy)
      CALL mp_sum(mp2_env%ri_grad%P_ij(my_spin)%array, para_env%group)
      CALL timestop(handle2)

      ! Now create and fill the P matrix (MO)
      ! FOR NOW WE ASSUME P_ab AND P_ij ARE REPLICATED OVER EACH MPI RANK
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%P_mo)) THEN
         ALLOCATE (mp2_env%ri_grad%P_mo(SIZE(mp2_env%ri_grad%mo_coeff_o)))
         DO ispin = 1, SIZE(mp2_env%ri_grad%mo_coeff_o)
            NULLIFY (mp2_env%ri_grad%P_mo(ispin)%matrix)
         END DO
      END IF

      CALL timeset(routineN//"_PMO", handle2)
      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(mp2_env%ri_grad%P_mo(my_spin)%matrix, fm_struct_tmp, name="P_MP2_MO")
      CALL cp_fm_set_all(mp2_env%ri_grad%P_mo(my_spin)%matrix, 0.0_dp)

      ! start with the (easy) occ-occ block and locally held P_ab elements
      itmp = get_limit(virtual, para_env_sub%num_pe, para_env_sub%mepos)
      my_B_virtual_start = itmp(1)
      my_B_virtual_end = itmp(2)

      CALL cp_fm_get_info(mp2_env%ri_grad%P_mo(my_spin)%matrix, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          nrow_block=nrow_block, &
                          ncol_block=ncol_block)
      myprow = mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%context%mepos(1)
      mypcol = mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%context%mepos(2)
      nprow = mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%context%num_pe(1)
      npcol = mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%context%num_pe(2)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (i_global <= homo .AND. j_global <= homo) THEN
               mp2_env%ri_grad%P_mo(my_spin)%matrix%local_data(iiB, jjB) = mp2_env%ri_grad%P_ij(my_spin)%array(i_global, j_global)
            END IF
            IF ((my_B_virtual_start <= i_global - homo .AND. i_global - homo <= my_B_virtual_end) .AND. (j_global > homo)) THEN
               mp2_env%ri_grad%P_mo(my_spin)%matrix%local_data(iiB, jjB) = &
                  mp2_env%ri_grad%P_ab(my_spin)%array(i_global - homo - my_B_virtual_start + 1, j_global - homo)
            END IF
         END DO
      END DO
      ! deallocate the local P_ij
      DEALLOCATE (mp2_env%ri_grad%P_ij(my_spin)%array)

      ! send around the sub_group the local data and check if we
      ! have to update our block with external elements
      ALLOCATE (mepos_2_grid(2, 0:para_env_sub%num_pe - 1))
      CALL mp_allgather([myprow, mypcol], mepos_2_grid, para_env_sub%group)

      ALLOCATE (sizes(2, 0:para_env_sub%num_pe - 1))
      CALL mp_allgather([nrow_local, ncol_local], sizes, para_env_sub%group)

      ALLOCATE (ab_rec(nrow_local, ncol_local))
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

         send_prow = mepos_2_grid(1, proc_send)
         send_pcol = mepos_2_grid(2, proc_send)

         send_row_size = sizes(1, proc_send)
         send_col_size = sizes(2, proc_send)

         ALLOCATE (ab_send(send_row_size, send_col_size))
         ab_send = 0.0_dp

         ! first loop over row since in this way we can cycle
         DO iiB = 1, send_row_size
            i_global = cp_fm_indxl2g(iiB, nrow_block, send_prow, &
                                     mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%first_p_pos(1), nprow)
            IF (i_global <= homo) CYCLE
            i_global = i_global - homo
            IF (.NOT. (my_B_virtual_start <= i_global .AND. i_global <= my_B_virtual_end)) CYCLE
            DO jjB = 1, send_col_size
               j_global = cp_fm_indxl2g(jjB, ncol_block, send_pcol, &
                                        mp2_env%ri_grad%P_mo(my_spin)%matrix%matrix_struct%first_p_pos(2), npcol)
               IF (j_global <= homo) CYCLE
               j_global = j_global - homo
               ab_send(iiB, jjB) = mp2_env%ri_grad%P_ab(my_spin)%array(i_global - my_B_virtual_start + 1, j_global)
            END DO
         END DO

         ab_rec = 0.0_dp
         CALL mp_sendrecv(ab_send, proc_send, &
                          ab_rec, proc_receive, &
                          para_env_sub%group)
         mp2_env%ri_grad%P_mo(my_spin)%matrix%local_data(1:nrow_local, 1:ncol_local) = &
            mp2_env%ri_grad%P_mo(my_spin)%matrix%local_data(1:nrow_local, 1:ncol_local) + &
            ab_rec(1:nrow_local, 1:ncol_local)

         DEALLOCATE (ab_send)
      END DO
      DEALLOCATE (ab_rec)
      DEALLOCATE (mepos_2_grid)
      DEALLOCATE (sizes)

      ! deallocate the local P_ab
      DEALLOCATE (mp2_env%ri_grad%P_ab(my_spin)%array)
      CALL timestop(handle2)

      ! create also W_MP2_MO
      CALL timeset(routineN//"_WMO", handle2)
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%W_mo)) THEN
         ALLOCATE (mp2_env%ri_grad%W_mo(SIZE(mp2_env%ri_grad%mo_coeff_o)))
         DO ispin = 1, SIZE(mp2_env%ri_grad%mo_coeff_o)
            NULLIFY (mp2_env%ri_grad%W_mo(ispin)%matrix)
         END DO
      END IF

      CALL cp_fm_create(mp2_env%ri_grad%W_mo(my_spin)%matrix, fm_struct_tmp, name="W_MP2_MO")
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! all block
      CALL cp_gemm('T', 'N', dimen, dimen, dimen, 2.0_dp*factor, &
                   L_mu_q, mo_coeff, 0.0_dp, mp2_env%ri_grad%W_mo(my_spin)%matrix, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)

      ! occ-occ block
      CALL cp_gemm('T', 'N', homo, homo, dimen, -2.0_dp*factor, &
                   L_mu_q, mo_coeff, 0.0_dp, mp2_env%ri_grad%W_mo(my_spin)%matrix, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)

      ! occ-virt block
      CALL cp_gemm('T', 'N', homo, virtual, dimen, 2.0_dp*factor, &
                   mo_coeff, L_mu_q, 0.0_dp, mp2_env%ri_grad%W_mo(my_spin)%matrix, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=homo + 1, &
                   b_first_row=1, &
                   c_first_col=homo + 1, &
                   c_first_row=1)
      CALL timestop(handle2)

      ! Calculate occ-virt block of the lagrangian in MO
      CALL timeset(routineN//"_Ljb", handle2)
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%L_jb)) THEN
         ALLOCATE (mp2_env%ri_grad%L_jb(SIZE(mp2_env%ri_grad%mo_coeff_o)))
         DO ispin = 1, SIZE(mp2_env%ri_grad%mo_coeff_o)
            NULLIFY (mp2_env%ri_grad%L_jb(ispin)%matrix)
         END DO
      END IF

      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=homo, ncol_global=virtual)
      CALL cp_fm_create(mp2_env%ri_grad%L_jb(my_spin)%matrix, fm_struct_tmp, name="fm_L_jb")
      CALL cp_fm_struct_release(fm_struct_tmp)

      ! first Virtual
      CALL cp_gemm('T', 'N', homo, virtual, dimen, 2.0_dp*factor, &
                   L_mu_q, mo_coeff, 0.0_dp, mp2_env%ri_grad%L_jb(my_spin)%matrix, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=homo + 1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)
      ! then occupied
      CALL cp_gemm('T', 'N', homo, virtual, dimen, 2.0_dp*factor, &
                   mo_coeff, L_mu_q, 1.0_dp, mp2_env%ri_grad%L_jb(my_spin)%matrix, &
                   a_first_col=1, &
                   a_first_row=1, &
                   b_first_col=homo + 1, &
                   b_first_row=1, &
                   c_first_col=1, &
                   c_first_row=1)

      ! finally release L_mu_q
      CALL cp_fm_release(L_mu_q)
      CALL timestop(handle2)

      ! here we should be done next CPHF

      DEALLOCATE (sub_proc_map)
      DEALLOCATE (pos_info)

      CALL timestop(handle)

   END SUBROUTINE create_W_P

END MODULE mp2_ri_grad
