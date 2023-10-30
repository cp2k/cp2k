!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to perform the RTP in the velocity gauge
! **************************************************************************************************

MODULE rt_propagation_velocity_gauge
   USE ai_moments,                      ONLY: cossin
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Mattiat2022,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_get_block_p,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE efield_utils,                    ONLY: make_field
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: one,&
                                              zero
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              clist_type,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_velocity_gauge'

   PUBLIC :: velocity_gauge_ks_matrix, update_vector_potential, velocity_gauge_nl_force

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param subtract_nl_term ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_ks_matrix(qs_env, subtract_nl_term)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: subtract_nl_term

      CHARACTER(len=*), PARAMETER :: routineN = 'velocity_gauge_ks_matrix'

      INTEGER                                            :: handle, idir, image, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: calculate_forces, my_subtract_nl_term, &
                                                            ppnl_present, use_virial
      REAL(KIND=dp)                                      :: eps_ppnl, factor
      REAL(KIND=dp), DIMENSION(3)                        :: vec_pot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: momentum, nl_term
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_h_im, matrix_nl, &
                                                            matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL cite_reference(Mattiat2022)

      my_subtract_nl_term = .FALSE.
      IF (PRESENT(subtract_nl_term)) my_subtract_nl_term = subtract_nl_term

      NULLIFY (dft_control, matrix_s, sab_orb, matrix_h, cell, input, matrix_h_im, kpoints, cell_to_index, &
               sap_ppnl, particle_set, qs_kind_set, atomic_kind_set, virial, force, matrix_p, rho, matrix_nl)

      CALL get_qs_env(qs_env, &
                      rho=rho, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      sap_ppnl=sap_ppnl, &
                      matrix_s_kp=matrix_s, &
                      matrix_h_kp=matrix_h, &
                      cell=cell, &
                      input=input, &
                      matrix_h_im_kp=matrix_h_im)

      nimages = dft_control%nimages
      ppnl_present = ASSOCIATED(sap_ppnl)

      IF (nimages > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      IF (my_subtract_nl_term) THEN
         IF (ppnl_present) THEN
            CALL get_qs_env(qs_env, &
                            qs_kind_set=qs_kind_set, &
                            particle_set=particle_set, &
                            atomic_kind_set=atomic_kind_set, &
                            virial=virial, &
                            rho=rho, &
                            force=force)

            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            calculate_forces = .FALSE.
            use_virial = .FALSE.
            nder = 1
            eps_ppnl = dft_control%qs_control%eps_ppnl

            CALL dbcsr_allocate_matrix_set(matrix_nl, 1, nimages)
            DO image = 1, nimages
               ALLOCATE (matrix_nl(1, image)%matrix)
               CALL dbcsr_create(matrix_nl(1, image)%matrix, template=matrix_s(1, 1)%matrix)
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_nl(1, image)%matrix, sab_orb)
               CALL dbcsr_set(matrix_nl(1, image)%matrix, zero)
            END DO

            CALL build_core_ppnl(matrix_nl, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                                 qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                                 nimages, cell_to_index, "ORB")

            DO image = 1, nimages
               CALL dbcsr_add(matrix_h(1, image)%matrix, matrix_nl(1, image)%matrix, one, -one)
            END DO

            CALL dbcsr_deallocate_matrix_set(matrix_nl)
         END IF
      END IF

      !get vector potential
      vec_pot = dft_control%rtp_control%vec_pot

      ! allocate and build matrices for linear momentum term
      NULLIFY (momentum)
      CALL dbcsr_allocate_matrix_set(momentum, 3)
      DO idir = 1, 3
         CALL dbcsr_init_p(momentum(idir)%matrix)
         CALL dbcsr_create(momentum(idir)%matrix, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(momentum(idir)%matrix, sab_orb)
         CALL dbcsr_set(momentum(idir)%matrix, zero)
      END DO
      CALL build_lin_mom_matrix(qs_env, momentum)

      ! set imaginary part of KS matrix to zero
      DO image = 1, nimages
         CALL dbcsr_set(matrix_h_im(1, image)%matrix, zero)
      END DO

      ! add linear term in vector potential to imaginary part of KS-matrix
      DO image = 1, nimages
         DO idir = 1, 3
            CALL dbcsr_add(matrix_h_im(1, image)%matrix, momentum(idir)%matrix, one, -vec_pot(idir))
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(momentum)

      ! add quadratic term to real part of KS matrix
      factor = 0._dp
      DO idir = 1, 3
         factor = factor + vec_pot(idir)**2
      END DO

      DO image = 1, nimages
         CALL dbcsr_add(matrix_h(1, image)%matrix, matrix_s(1, image)%matrix, one, 0.5*factor)
      END DO

      ! add Non local term
      IF (ppnl_present) THEN
         IF (dft_control%rtp_control%nl_gauge_transform) THEN
            NULLIFY (nl_term)
            CALL dbcsr_allocate_matrix_set(nl_term, 2)

            CALL dbcsr_init_p(nl_term(1)%matrix)
            CALL dbcsr_create(nl_term(1)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_symmetric, name="nl gauge term real part")
            CALL cp_dbcsr_alloc_block_from_nbl(nl_term(1)%matrix, sab_orb)
            CALL dbcsr_set(nl_term(1)%matrix, zero)

            CALL dbcsr_init_p(nl_term(2)%matrix)
            CALL dbcsr_create(nl_term(2)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric, name="nl gauge term imaginary part")
            CALL cp_dbcsr_alloc_block_from_nbl(nl_term(2)%matrix, sab_orb)
            CALL dbcsr_set(nl_term(2)%matrix, zero)

            CALL velocity_gauge_nl_term(qs_env, nl_term, vec_pot)

            DO image = 1, nimages
               CALL dbcsr_add(matrix_h(1, image)%matrix, nl_term(1)%matrix, one, one)
               CALL dbcsr_add(matrix_h_im(1, image)%matrix, nl_term(2)%matrix, one, one)
            END DO
            CALL dbcsr_deallocate_matrix_set(nl_term)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_ks_matrix

! **************************************************************************************************
!> \brief Update the vector potential in the case where a time-dependant
!>        electric field is apply.
!> \param qs_env ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE update_vector_potential(qs_env, dft_control)
      TYPE(qs_environment_type), INTENT(INOUT), POINTER  :: qs_env
      TYPE(dft_control_type), INTENT(INOUT), POINTER     :: dft_control

      REAL(kind=dp)                                      :: field(3)

      CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
      dft_control%rtp_control%vec_pot = dft_control%rtp_control%vec_pot - field*qs_env%rtp%dt
      ! Update the vec_pot_initial value for RTP restart:
      dft_control%efield_fields(1)%efield%vec_pot_initial = dft_control%rtp_control%vec_pot

   END SUBROUTINE update_vector_potential

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nl_term ...
!> \param vec_pot ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_nl_term(qs_env, nl_term, vec_pot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: nl_term
      REAL(KIND=dp), DIMENSION(3), INTENT(in)            :: vec_pot

      CHARACTER(len=*), PARAMETER :: routiuneN = "velocity_gauge_nl_term"

      INTEGER                                            :: handle, i, iac, iatom, ibc, icol, ikind, &
                                                            irow, jatom, jkind, kac, kbc, kkind, &
                                                            maxl, maxlgto, maxlppnl, na, natom, &
                                                            nb, nkind, np, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found
      REAL(dp)                                           :: eps_ppnl
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: imag_block, real_block
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: achint_cos, achint_sin, acint_cos, &
                                                            acint_sin, bchint_cos, bchint_sin, &
                                                            bcint_cos, bcint_sin
      TYPE(alist_type), POINTER                          :: alist_cos_ac, alist_cos_bc, &
                                                            alist_sin_ac, alist_sin_bc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int_cos, sap_int_sin

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      MARK_USED(int_8)

      CALL timeset(routiuneN, handle)

      NULLIFY (sap_ppnl, sab_orb)
      CALL get_qs_env(qs_env, &
                      sap_ppnl=sap_ppnl, &
                      sab_orb=sab_orb)

      IF (ASSOCIATED(sap_ppnl)) THEN
         NULLIFY (qs_kind_set, particle_set, cell, dft_control)
         CALL get_qs_env(qs_env, &
                         dft_control=dft_control, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         cell=cell, &
                         atomic_kind_set=atomic_kind_set)

         nkind = SIZE(atomic_kind_set)
         natom = SIZE(particle_set)
         eps_ppnl = dft_control%qs_control%eps_ppnl

         CALL get_qs_kind_set(qs_kind_set, &
                              maxlgto=maxlgto, &
                              maxlppnl=maxlppnl)

         maxl = MAX(maxlppnl, maxlgto)
         CALL init_orbital_pointers(maxl + 1)

         ! initalize sab_int types to store the integrals
         NULLIFY (sap_int_cos, sap_int_sin)
         ALLOCATE (sap_int_cos(nkind*nkind), sap_int_sin(nkind*nkind))
         DO i = 1, SIZE(sap_int_cos)
            NULLIFY (sap_int_cos(i)%alist, sap_int_cos(i)%asort, sap_int_cos(i)%aindex)
            sap_int_cos(i)%nalist = 0
            NULLIFY (sap_int_sin(i)%alist, sap_int_sin(i)%asort, sap_int_sin(i)%aindex)
            sap_int_sin(i)%nalist = 0
         END DO

         ! get basis set
         ALLOCATE (basis_set(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
         END DO

         ! calculate exponential integrals
         CALL build_sap_exp_ints(sap_int_cos, sap_int_sin, sap_ppnl, qs_kind_set, particle_set, &
                                 cell, kvec=vec_pot, basis_set=basis_set, nkind=nkind, &
                                 derivative=.FALSE.)

         CALL sap_sort(sap_int_cos)
         CALL sap_sort(sap_int_sin)

         ! assemble the integrals for the gauge term
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (basis_set, nl_term, sab_orb, sap_int_cos, sap_int_sin, eps_ppnl, locks, nkind, natom) &
!$OMP PRIVATE (real_block, imag_block, acint_cos, achint_cos, bcint_cos, bchint_cos, acint_sin,&
!$OMP          achint_sin, bcint_sin, bchint_sin, slot, ikind, jkind, iatom, jatom, cell_b, rab, irow, icol,&
!$OMP          found, kkind, iac, ibc, alist_cos_ac, alist_cos_bc, alist_sin_ac, alist_sin_bc, kac, kbc, &
!$OMP          na, np, nb, iatom8, hash)

!$OMP SINGLE
!$       ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$       DO lock_num = 1, nlock
!$          call omp_init_lock(locks(lock_num))
!$       END DO
!$OMP END DO

         NULLIFY (real_block, imag_block)
         NULLIFY (acint_cos, bcint_cos, achint_cos, bchint_cos)
         NULLIFY (acint_sin, bcint_sin, achint_sin, bchint_sin)

         ! loop over atom pairs
!$OMP DO SCHEDULE(GUIDED)
         DO slot = 1, sab_orb(1)%nl_size
            ikind = sab_orb(1)%nlist_task(slot)%ikind
            jkind = sab_orb(1)%nlist_task(slot)%jkind
            iatom = sab_orb(1)%nlist_task(slot)%iatom
            jatom = sab_orb(1)%nlist_task(slot)%jatom
            cell_b(:) = sab_orb(1)%nlist_task(slot)%cell
            rab(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)

            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE

            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF

            CALL dbcsr_get_block_p(nl_term(1)%matrix, irow, icol, real_block, found)
            CALL dbcsr_get_block_p(nl_term(2)%matrix, irow, icol, imag_block, found)

            IF (ASSOCIATED(real_block) .AND. ASSOCIATED(imag_block)) THEN
               ! loop over the <gto_a|ppln_c>h_ij<ppnl_c|gto_b> pairs
               DO kkind = 1, nkind
                  iac = ikind + nkind*(kkind - 1)
                  ibc = jkind + nkind*(kkind - 1)
                  IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_cos(ibc)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_sin(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int_cos(iac), alist_cos_ac, iatom)
                  CALL get_alist(sap_int_cos(ibc), alist_cos_bc, jatom)
                  CALL get_alist(sap_int_sin(iac), alist_sin_ac, iatom)
                  CALL get_alist(sap_int_sin(ibc), alist_sin_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_cos_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_cos_bc)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_sin_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_sin_bc)) CYCLE

                  ! only use cos for indexing, as cos and sin integrals are constructed by the same routine
                  ! in the same way
                  DO kac = 1, alist_cos_ac%nclist
                     DO kbc = 1, alist_cos_bc%nclist
                        ! the next two ifs should be the same for sine integrals
                        IF (alist_cos_ac%clist(kac)%catom /= alist_cos_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b + alist_cos_bc%clist(kbc)%cell - alist_cos_ac%clist(kac)%cell == 0)) THEN
                           ! screening
                           IF (alist_cos_ac%clist(kac)%maxac*alist_cos_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_cos_ac%clist(kac)%maxac*alist_sin_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_sin_ac%clist(kac)%maxac*alist_cos_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_sin_ac%clist(kac)%maxac*alist_sin_bc%clist(kbc)%maxach < eps_ppnl) CYCLE

                           acint_cos => alist_cos_ac%clist(kac)%acint
                           bcint_cos => alist_cos_bc%clist(kbc)%acint
                           achint_cos => alist_cos_ac%clist(kac)%achint
                           bchint_cos => alist_cos_bc%clist(kbc)%achint
                           acint_sin => alist_sin_ac%clist(kac)%acint
                           bcint_sin => alist_sin_bc%clist(kbc)%acint
                           achint_sin => alist_sin_ac%clist(kac)%achint
                           bchint_sin => alist_sin_bc%clist(kbc)%achint

                           na = SIZE(acint_cos, 1)
                           np = SIZE(acint_cos, 2)
                           nb = SIZE(bcint_cos, 1)
!$                         iatom8 = INT(iatom - 1, int_8)*INT(natom, int_8) + INT(jatom, int_8)
!$                         hash = INT(MOD(iatom8, INT(nlock, int_8)) + 1)
!$                         CALL omp_set_lock(locks(hash))
                           IF (iatom <= jatom) THEN
                              ! cos*cos + sin*sin
                              real_block(1:na, 1:nb) = real_block(1:na, 1:nb) + &
                                 MATMUL(achint_cos(1:na, 1:np, 1), TRANSPOSE(bcint_cos(1:nb, 1:np, 1))) + &
                                 MATMUL(achint_sin(1:na, 1:np, 1), TRANSPOSE(bcint_sin(1:nb, 1:np, 1)))
                              ! sin * cos - cos * sin
                              imag_block(1:na, 1:nb) = imag_block(1:na, 1:nb) - &
                                                       MATMUL(achint_sin(1:na, 1:np, 1), TRANSPOSE(bcint_cos(1:nb, 1:np, 1))) + &
                                                       MATMUL(achint_cos(1:na, 1:np, 1), TRANSPOSE(bcint_sin(1:nb, 1:np, 1)))
                           ELSE
                              ! cos*cos + sin*sin
                              real_block(1:nb, 1:na) = real_block(1:nb, 1:na) + &
                                 MATMUL(bchint_cos(1:nb, 1:np, 1), TRANSPOSE(acint_cos(1:na, 1:np, 1))) + &
                                 MATMUL(bchint_sin(1:nb, 1:np, 1), TRANSPOSE(acint_sin(1:na, 1:np, 1)))
                              ! sin * cos - cos * sin
                              imag_block(1:nb, 1:na) = imag_block(1:nb, 1:na) - &
                                                       MATMUL(bchint_sin(1:nb, 1:np, 1), TRANSPOSE(acint_cos(1:na, 1:np, 1))) + &
                                                       MATMUL(bchint_cos(1:nb, 1:np, 1), TRANSPOSE(acint_sin(1:na, 1:np, 1)))

                           END IF
!$                         CALL omp_unset_lock(locks(hash))
                           EXIT
                        END IF
                     END DO
                  END DO
               END DO
            END IF

         END DO

!$OMP DO
!$       DO lock_num = 1, nlock
!$          call omp_destroy_lock(locks(lock_num))
!$       END DO
!$OMP END DO

!$OMP SINGLE
!$       DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL
         CALL release_sap_int(sap_int_cos)
         CALL release_sap_int(sap_int_sin)

         DEALLOCATE (basis_set)
      END IF

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_nl_term

! **************************************************************************************************
!> \brief calculate <a|sin/cos|p> integrals and store in sap_int_type
!>        adapted from build_sap_ints
!>        Do this on each MPI task as the integrals need to be available globally.
!>        Might be faster than communicating as the integrals are obtained analytically.
!>        If asked, compute <da/dRa|sin/cos|p>
!> \param sap_int_cos ...
!> \param sap_int_sin ...
!> \param sap_ppnl ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param kvec ...
!> \param basis_set ...
!> \param nkind ...
!> \param derivative ...
! **************************************************************************************************
   SUBROUTINE build_sap_exp_ints(sap_int_cos, sap_int_sin, sap_ppnl, qs_kind_set, particle_set, cell, &
                                 kvec, basis_set, nkind, derivative)
      TYPE(sap_int_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: sap_int_cos, sap_int_sin
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sap_ppnl
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(in)            :: kvec
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: basis_set
      INTEGER, INTENT(IN)                                :: nkind
      LOGICAL, INTENT(IN)                                :: derivative

      CHARACTER(len=*), PARAMETER :: routiuneN = "build_sap_exp_ints"

      INTEGER :: handle, i, iac, iatom, idir, ikind, ilist, iset, jneighbor, katom, kkind, l, &
         lc_max, lc_min, ldai, ldints, lppnl, maxco, maxl, maxlgto, maxlppnl, maxppnl, maxsgf, na, &
         nb, ncoa, ncoc, nlist, nneighbor, np, nppnl, nprjc, nseta, nsgfa, prjc, sgfa, slot
      INTEGER, DIMENSION(3)                              :: cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nprj_ppnl, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: dogth
      REAL(KIND=dp)                                      :: dac, ppnl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ai_work_cos, ai_work_sin, work_cos, &
                                                            work_sin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ai_work_dcos, ai_work_dsin, work_dcos, &
                                                            work_dsin
      REAL(KIND=dp), DIMENSION(1)                        :: rprjc, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rac, raf, rc, rcf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_ppnl, set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cprj, rpgfa, sphi_a, vprj_ppnl, zeta
      TYPE(clist_type), POINTER                          :: clist, clist_sin
      TYPE(gth_potential_p_type), DIMENSION(:), POINTER  :: gpotential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_p_type), DIMENSION(:), POINTER  :: spotential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routiuneN, handle)

      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, &
                           maxlppnl=maxlppnl, &
                           maxppnl=maxppnl, &
                           maxsgf=maxsgf, &
                           maxlgto=maxlgto)

      ! maximum dimensions for allocations
      maxl = MAX(maxlppnl, maxlgto)
      ldints = MAX(maxco, ncoset(maxlppnl), maxsgf, maxppnl)
      ldai = ncoset(maxl + 1)

      !set up direct access to basis and potential
      NULLIFY (gpotential, spotential)
      ALLOCATE (gpotential(nkind), spotential(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential, sgp_potential=sgp_potential)
         NULLIFY (gpotential(ikind)%gth_potential)
         NULLIFY (spotential(ikind)%sgp_potential)
         IF (ASSOCIATED(gth_potential)) THEN
            gpotential(ikind)%gth_potential => gth_potential
         ELSE IF (ASSOCIATED(sgp_potential)) THEN
            spotential(ikind)%sgp_potential => sgp_potential
         END IF
      END DO

      !allocate sap int
      NULLIFY (clist)
      DO slot = 1, sap_ppnl(1)%nl_size

         ikind = sap_ppnl(1)%nlist_task(slot)%ikind
         kkind = sap_ppnl(1)%nlist_task(slot)%jkind
         iatom = sap_ppnl(1)%nlist_task(slot)%iatom
         katom = sap_ppnl(1)%nlist_task(slot)%jatom
         nlist = sap_ppnl(1)%nlist_task(slot)%nlist
         ilist = sap_ppnl(1)%nlist_task(slot)%ilist
         nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode

         iac = ikind + nkind*(kkind - 1)
         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(gpotential(kkind)%gth_potential) .AND. &
             .NOT. ASSOCIATED(spotential(kkind)%sgp_potential)) CYCLE
         IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist)) THEN
            sap_int_cos(iac)%a_kind = ikind
            sap_int_cos(iac)%p_kind = kkind
            sap_int_cos(iac)%nalist = nlist
            ALLOCATE (sap_int_cos(iac)%alist(nlist))
            DO i = 1, nlist
               NULLIFY (sap_int_cos(iac)%alist(i)%clist)
               sap_int_cos(iac)%alist(i)%aatom = 0
               sap_int_cos(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist(ilist)%clist)) THEN
            sap_int_cos(iac)%alist(ilist)%aatom = iatom
            sap_int_cos(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE (sap_int_cos(iac)%alist(ilist)%clist(nneighbor))
            DO i = 1, nneighbor
               clist => sap_int_cos(iac)%alist(ilist)%clist(i)
               clist%catom = 0
               NULLIFY (clist%acint)
               NULLIFY (clist%achint)
               NULLIFY (clist%sgf_list)
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist)) THEN
            sap_int_sin(iac)%a_kind = ikind
            sap_int_sin(iac)%p_kind = kkind
            sap_int_sin(iac)%nalist = nlist
            ALLOCATE (sap_int_sin(iac)%alist(nlist))
            DO i = 1, nlist
               NULLIFY (sap_int_sin(iac)%alist(i)%clist)
               sap_int_sin(iac)%alist(i)%aatom = 0
               sap_int_sin(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist(ilist)%clist)) THEN
            sap_int_sin(iac)%alist(ilist)%aatom = iatom
            sap_int_sin(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE (sap_int_sin(iac)%alist(ilist)%clist(nneighbor))
            DO i = 1, nneighbor
               clist => sap_int_sin(iac)%alist(ilist)%clist(i)
               clist%catom = 0
               NULLIFY (clist%acint)
               NULLIFY (clist%achint)
               NULLIFY (clist%sgf_list)
            END DO
         END IF
      END DO

      ! actual calculation of the integrals <a|cos|p> and <a|sin|p>
      ! allocate temporary storage using maximum dimensions

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (basis_set, gpotential, ncoset, sap_ppnl, sap_int_cos, sap_int_sin, nkind, &
!$OMP         ldints, maxco, nco, cell, particle_set, kvec, derivative) &
!$OMP PRIVATE (slot, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, dac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta,&
!$OMP          rpgfa, set_radius_a, sphi_a, zeta, alpha_ppnl, cprj, lppnl, nppnl, nprj_ppnl,&
!$OMP          ppnl_radius, vprj_ppnl, clist, clist_sin, ra, rc, ncoa, sgfa, prjc, work_cos, work_sin,&
!$OMP          nprjc, rprjc, lc_max, lc_min, zetc, ncoc, ai_work_sin, ai_work_cos, na, nb, np, dogth, &
!$OMP          raf, rcf,  work_dcos, work_dsin, ai_work_dcos, ai_work_dsin, idir)

      ALLOCATE (work_cos(ldints, ldints), work_sin(ldints, ldints))
      ALLOCATE (ai_work_cos(maxco, maxco), ai_work_sin(maxco, maxco))
      IF (derivative) THEN
         ALLOCATE (work_dcos(ldints, ldints, 3), work_dsin(ldints, ldints, 3))
         ALLOCATE (ai_work_dcos(maxco, maxco, 3), ai_work_dsin(maxco, maxco, 3))
      END IF
      work_cos = 0.0_dp
      work_sin = 0.0_dp
      ai_work_cos = 0.0_dp
      ai_work_sin = 0.0_dp
      IF (derivative) THEN
         ai_work_dcos = 0.0_dp
         ai_work_dsin = 0.0_dp
      END IF
      dogth = .FALSE.

      NULLIFY (first_sgfa, la_max, la_min, npgfa, nsgf_seta, rpgfa, set_radius_a, sphi_a, zeta)
      NULLIFY (alpha_ppnl, cprj, nprj_ppnl, vprj_ppnl)
      NULLIFY (clist, clist_sin)

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sap_ppnl(1)%nl_size
         ikind = sap_ppnl(1)%nlist_task(slot)%ikind
         kkind = sap_ppnl(1)%nlist_task(slot)%jkind
         iatom = sap_ppnl(1)%nlist_task(slot)%iatom
         katom = sap_ppnl(1)%nlist_task(slot)%jatom
         nlist = sap_ppnl(1)%nlist_task(slot)%nlist
         ilist = sap_ppnl(1)%nlist_task(slot)%ilist
         nneighbor = sap_ppnl(1)%nlist_task(slot)%nnode
         jneighbor = sap_ppnl(1)%nlist_task(slot)%inode
         cell_c(:) = sap_ppnl(1)%nlist_task(slot)%cell(:)
         rac(1:3) = sap_ppnl(1)%nlist_task(slot)%r(1:3)
         dac = NORM2(rac)

         iac = ikind + nkind*(kkind - 1)
         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         ! get definition of gto basis set
         first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
         la_max => basis_set(ikind)%gto_basis_set%lmax
         la_min => basis_set(ikind)%gto_basis_set%lmin
         npgfa => basis_set(ikind)%gto_basis_set%npgf
         nseta = basis_set(ikind)%gto_basis_set%nset
         nsgfa = basis_set(ikind)%gto_basis_set%nsgf
         nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
         rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
         set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
         sphi_a => basis_set(ikind)%gto_basis_set%sphi
         zeta => basis_set(ikind)%gto_basis_set%zet

         IF (ASSOCIATED(gpotential(kkind)%gth_potential)) THEN
            ! GTH potential
            dogth = .TRUE.
            alpha_ppnl => gpotential(kkind)%gth_potential%alpha_ppnl
            cprj => gpotential(kkind)%gth_potential%cprj
            lppnl = gpotential(kkind)%gth_potential%lppnl
            nppnl = gpotential(kkind)%gth_potential%nppnl
            nprj_ppnl => gpotential(kkind)%gth_potential%nprj_ppnl
            ppnl_radius = gpotential(kkind)%gth_potential%ppnl_radius
            vprj_ppnl => gpotential(kkind)%gth_potential%vprj_ppnl
         ELSE
            CYCLE
         END IF

         clist => sap_int_cos(iac)%alist(ilist)%clist(jneighbor)
         clist_sin => sap_int_sin(iac)%alist(ilist)%clist(jneighbor)

         clist%catom = katom
         clist%cell = cell_c
         clist%rac = rac
         clist_sin%catom = katom
         clist_sin%cell = cell_c
         clist_sin%rac = rac

         IF (.NOT. derivative) THEN
            ALLOCATE (clist%acint(nsgfa, nppnl, 1), clist%achint(nsgfa, nppnl, 1))
         ELSE
            ALLOCATE (clist%acint(nsgfa, nppnl, 4), clist%achint(nsgfa, nppnl, 4))
         END IF
         clist%acint = 0.0_dp
         clist%achint = 0.0_dp
         clist%nsgf_cnt = 0

         IF (.NOT. derivative) THEN
            ALLOCATE (clist_sin%acint(nsgfa, nppnl, 1), clist_sin%achint(nsgfa, nppnl, 1))
         ELSE
            ALLOCATE (clist_sin%acint(nsgfa, nppnl, 4), clist_sin%achint(nsgfa, nppnl, 4))
         END IF
         clist_sin%acint = 0.0_dp
         clist_sin%achint = 0.0_dp
         clist_sin%nsgf_cnt = 0

         ! reference point at zero
         ra(:) = pbc(particle_set(iatom)%r(:), cell)
         rc(:) = ra + rac

         ! reference point at pseudized atom
         raf(:) = ra - rc
         rcf(:) = 0._dp

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            IF (dogth) THEN
               prjc = 1
               work_cos = 0.0_dp
               work_sin = 0.0_dp
               DO l = 0, lppnl
                  nprjc = nprj_ppnl(l)*nco(l)
                  IF (nprjc == 0) CYCLE
                  rprjc(1) = ppnl_radius
                  IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
                  lc_max = l + 2*(nprj_ppnl(l) - 1)
                  lc_min = l
                  zetc(1) = alpha_ppnl(l)
                  ncoc = ncoset(lc_max)

                  IF (.NOT. derivative) THEN
                     CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                 lc_max, 1, zetc, rprjc, lc_min, raf, rcf, kvec, ai_work_cos, ai_work_sin)
                  ELSE
                     CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                 lc_max, 1, zetc, rprjc, lc_min, raf, rcf, kvec, ai_work_cos, ai_work_sin, &
                                 dcosab=ai_work_dcos, dsinab=ai_work_dsin)
                  END IF
                  ! projector functions: Cartesian -> spherical
                  na = ncoa
                  nb = nprjc
                  np = ncoc
                  work_cos(1:na, prjc:prjc + nb - 1) = &
                     MATMUL(ai_work_cos(1:na, 1:np), cprj(1:np, prjc:prjc + nb - 1))
                  work_sin(1:na, prjc:prjc + nb - 1) = &
                     MATMUL(ai_work_sin(1:na, 1:np), cprj(1:np, prjc:prjc + nb - 1))

                  IF (derivative) THEN
                     DO idir = 1, 3
                        work_dcos(1:na, prjc:prjc + nb - 1, idir) = &
                           MATMUL(ai_work_dcos(1:na, 1:np, idir), cprj(1:np, prjc:prjc + nb - 1))
                        work_dsin(1:na, prjc:prjc + nb - 1, idir) = &
                           MATMUL(ai_work_dsin(1:na, 1:np, idir), cprj(1:np, prjc:prjc + nb - 1))
                     END DO
                  END IF

                  prjc = prjc + nprjc
               END DO

               ! contract gto basis set into acint
               na = nsgf_seta(iset)
               nb = nppnl
               np = ncoa
               clist%acint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_cos(1:np, 1:nb))
               clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_sin(1:np, 1:nb))
               IF (derivative) THEN
                  DO idir = 1, 3
                     clist%acint(sgfa:sgfa + na - 1, 1:nb, 1 + idir) = &
                        MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_dcos(1:np, 1:nb, idir))
                     clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1 + idir) = &
                        MATMUL(TRANSPOSE(sphi_a(1:np, sgfa:sgfa + na - 1)), work_dsin(1:np, 1:nb, idir))
                  END DO
               END IF

               ! multiply with interaction matrix h_ij of the nl pp
               clist%achint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(clist%acint(sgfa:sgfa + na - 1, 1:nb, 1), vprj_ppnl(1:nb, 1:nb))
               clist_sin%achint(sgfa:sgfa + na - 1, 1:nb, 1) = &
                  MATMUL(clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1), vprj_ppnl(1:nb, 1:nb))
               IF (derivative) THEN
                  DO idir = 1, 3
                     clist%achint(sgfa:sgfa + na - 1, 1:nb, 1 + idir) = &
                        MATMUL(clist%acint(sgfa:sgfa + na - 1, 1:nb, 1 + idir), vprj_ppnl(1:nb, 1:nb))
                     clist_sin%achint(sgfa:sgfa + na - 1, 1:nb, 1 + idir) = &
                        MATMUL(clist_sin%acint(sgfa:sgfa + na - 1, 1:nb, 1 + idir), vprj_ppnl(1:nb, 1:nb))
                  END DO
               END IF
            END IF

         END DO
         clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
         clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         clist_sin%maxac = MAXVAL(ABS(clist_sin%acint(:, :, 1)))
         clist_sin%maxach = MAXVAL(ABS(clist_sin%achint(:, :, 1)))
      END DO

      DEALLOCATE (work_cos, work_sin, ai_work_cos, ai_work_sin)
      IF (derivative) DEALLOCATE (work_dcos, work_dsin, ai_work_dcos, ai_work_dsin)

!$OMP END PARALLEL

      DEALLOCATE (gpotential, spotential)

      CALL timestop(handle)

   END SUBROUTINE build_sap_exp_ints

! **************************************************************************************************
!> \brief Calculate the force associated to non-local pseudo potential in the velocity gauge
!> \param qs_env ...
!> \param particle_set ...
!> \date    09.2023
!> \author  Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE velocity_gauge_nl_force(qs_env, particle_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CHARACTER(len=*), PARAMETER :: routiuneN = "velocity_gauge_nl_force"

      INTEGER :: handle, i, iac, iatom, ibc, icol, idir, ikind, irow, jatom, jkind, kac, katom, &
         kbc, kkind, maxl, maxlgto, maxlppnl, na, natom, nb, nkind, np, slot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found_imag, found_real
      REAL(dp)                                           :: eps_ppnl, f0, sign_imag
      REAL(KIND=dp), DIMENSION(3)                        :: fa, fb, rab, vec_pot
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_im
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(alist_type), POINTER                          :: alist_cos_ac, alist_cos_bc, &
                                                            alist_sin_ac, alist_sin_bc
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: achint_cos, achint_sin, acint_cos, &
                                                            acint_sin, bchint_cos, bchint_sin, &
                                                            bcint_cos, bcint_sin
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix_p_imag, matrix_p_real
      REAL(KIND=dp), DIMENSION(3, SIZE(particle_set))    :: force_thread
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int_cos, sap_int_sin

      CALL timeset(routiuneN, handle)

      NULLIFY (sap_ppnl)

      CALL get_qs_env(qs_env, &
                      sap_ppnl=sap_ppnl)

      IF (ASSOCIATED(sap_ppnl)) THEN
         NULLIFY (qs_kind_set, cell, dft_control, force, sab_orb, atomic_kind_set, &
                  sap_int_cos, sap_int_sin)
         ! Load and initialized the required quantities

         CALL get_qs_env(qs_env, &
                         sab_orb=sab_orb, &
                         force=force, &
                         dft_control=dft_control, &
                         qs_kind_set=qs_kind_set, &
                         cell=cell, &
                         atomic_kind_set=atomic_kind_set, &
                         rho=rho)

         nkind = SIZE(atomic_kind_set)
         natom = SIZE(particle_set)
         eps_ppnl = dft_control%qs_control%eps_ppnl

         CALL get_qs_kind_set(qs_kind_set, &
                              maxlgto=maxlgto, &
                              maxlppnl=maxlppnl)

         maxl = MAX(maxlppnl, maxlgto)
         CALL init_orbital_pointers(maxl + 1)

         ! initalize sab_int types to store the integrals
         ALLOCATE (sap_int_cos(nkind*nkind), sap_int_sin(nkind*nkind))
         DO i = 1, SIZE(sap_int_cos)
            NULLIFY (sap_int_cos(i)%alist, sap_int_cos(i)%asort, sap_int_cos(i)%aindex)
            sap_int_cos(i)%nalist = 0
            NULLIFY (sap_int_sin(i)%alist, sap_int_sin(i)%asort, sap_int_sin(i)%aindex)
            sap_int_sin(i)%nalist = 0
         END DO

         ! get basis set
         ALLOCATE (basis_set(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
         END DO

         !get vector potential
         vec_pot = dft_control%rtp_control%vec_pot

         force_thread = 0.0_dp

         CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao, rho_ao_im=rho_ao_im)
         ! To avoid FOR loop over spin, sum the 2 spin into the first one directly. Undone later on
         IF (SIZE(rho_ao) == 2) THEN
            CALL dbcsr_add(rho_ao(1)%matrix, rho_ao(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            CALL dbcsr_add(rho_ao_im(1)%matrix, rho_ao_im(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END IF

         ! Compute cosap = <a|cos kr|p>, sindap = <a|sin kr|p>, cosdap = <da/dRA|cos kr|p>, and sindap = <da/dRA|sin kr|p>
         CALL build_sap_exp_ints(sap_int_cos, sap_int_sin, sap_ppnl, qs_kind_set, particle_set, &
                                 cell, kvec=vec_pot, basis_set=basis_set, nkind=nkind, derivative=.TRUE.)
         CALL sap_sort(sap_int_cos)
         CALL sap_sort(sap_int_sin)

         ! Compute the force, on nuclei A it is given by: Re(P_ab) Re(dV_ab/dRA) - Im(P_ab) Im(dV_ab/dRA)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (basis_set, sab_orb, sap_int_cos, sap_int_sin, eps_ppnl, nkind, natom,&
!$OMP         rho_ao, rho_ao_im) &
!$OMP PRIVATE (matrix_p_real, matrix_p_imag, acint_cos, achint_cos, bcint_cos, bchint_cos, acint_sin,&
!$OMP          achint_sin, bcint_sin, bchint_sin, slot, ikind, jkind, iatom, jatom,&
!$OMP          cell_b, rab, irow, icol, fa, fb, f0, found_real, found_imag, sign_imag, &
!$OMP          kkind, iac, ibc, alist_cos_ac, alist_cos_bc, alist_sin_ac, alist_sin_bc, kac, kbc,&
!$OMP          na, np, nb,  katom) &
!$OMP REDUCTION (+ : force_thread )

         NULLIFY (acint_cos, bcint_cos, achint_cos, bchint_cos)
         NULLIFY (acint_sin, bcint_sin, achint_sin, bchint_sin)

         ! loop over atom pairs
!$OMP DO SCHEDULE(GUIDED)
         DO slot = 1, sab_orb(1)%nl_size
            ikind = sab_orb(1)%nlist_task(slot)%ikind
            jkind = sab_orb(1)%nlist_task(slot)%jkind
            iatom = sab_orb(1)%nlist_task(slot)%iatom
            jatom = sab_orb(1)%nlist_task(slot)%jatom
            cell_b(:) = sab_orb(1)%nlist_task(slot)%cell
            rab(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)

            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE

            ! Use the symmetry of the first derivatives
            IF (iatom == jatom) THEN
               f0 = 1.0_dp
            ELSE
               f0 = 2.0_dp
            END IF

            fa = 0.0_dp
            fb = 0.0_dp

            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
               sign_imag = +1.0_dp
            ELSE
               irow = jatom
               icol = iatom
               sign_imag = -1.0_dp
            END IF
            NULLIFY (matrix_p_real, matrix_p_imag)
            CALL dbcsr_get_block_p(rho_ao(1)%matrix, irow, icol, matrix_p_real, found_real)
            CALL dbcsr_get_block_p(rho_ao_im(1)%matrix, irow, icol, matrix_p_imag, found_imag)

            IF (found_real .OR. found_imag) THEN
               ! loop over the <gto_a|ppln_c>h_ij<ppnl_c|gto_b> pairs
               DO kkind = 1, nkind
                  iac = ikind + nkind*(kkind - 1)
                  ibc = jkind + nkind*(kkind - 1)
                  IF (.NOT. ASSOCIATED(sap_int_cos(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_cos(ibc)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_sin(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int_sin(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int_cos(iac), alist_cos_ac, iatom)
                  CALL get_alist(sap_int_cos(ibc), alist_cos_bc, jatom)
                  CALL get_alist(sap_int_sin(iac), alist_sin_ac, iatom)
                  CALL get_alist(sap_int_sin(ibc), alist_sin_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_cos_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_cos_bc)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_sin_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_sin_bc)) CYCLE

                  ! only use cos for indexing, as cos and sin integrals are constructed by the same routine
                  ! in the same way
                  DO kac = 1, alist_cos_ac%nclist
                     DO kbc = 1, alist_cos_bc%nclist
                        ! the next two ifs should be the same for sine integrals
                        IF (alist_cos_ac%clist(kac)%catom /= alist_cos_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b + alist_cos_bc%clist(kbc)%cell - alist_cos_ac%clist(kac)%cell == 0)) THEN
                           ! screening
                           IF (alist_cos_ac%clist(kac)%maxac*alist_cos_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_cos_ac%clist(kac)%maxac*alist_sin_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_sin_ac%clist(kac)%maxac*alist_cos_bc%clist(kbc)%maxach < eps_ppnl &
                               .AND. alist_sin_ac%clist(kac)%maxac*alist_sin_bc%clist(kbc)%maxach < eps_ppnl) CYCLE

                           acint_cos => alist_cos_ac%clist(kac)%acint
                           bcint_cos => alist_cos_bc%clist(kbc)%acint
                           achint_cos => alist_cos_ac%clist(kac)%achint
                           bchint_cos => alist_cos_bc%clist(kbc)%achint
                           acint_sin => alist_sin_ac%clist(kac)%acint
                           bcint_sin => alist_sin_bc%clist(kbc)%acint
                           achint_sin => alist_sin_ac%clist(kac)%achint
                           bchint_sin => alist_sin_bc%clist(kbc)%achint

                           na = SIZE(acint_cos, 1)
                           np = SIZE(acint_cos, 2)
                           nb = SIZE(bcint_cos, 1)
                           ! Re(dV_ab/dRA) = <da/dRA|cos kr|p><p|cos kr|b> + <db/dRA|cos kr|p><p|cos kr|a> + <da/dRA|sin kr|p><p|sin kr|b> + <db/dRA|sin kr|p><p|sin|a>
                           ! Im(dV_ab/dRA) = <da/dRA|sin kr|p><p|cos kr|b> - <db/dRA|sin kr|p><p|cos kr|a> - <da/dRA|cos kr|p><p|sin kr|b> + <db/dRA|cos kr|p><p|sin|a>
                           katom = alist_cos_ac%clist(kac)%catom
                           DO idir = 1, 3
                              IF (iatom <= jatom) THEN
                                 ! For fa:
                                 IF (found_real) &
                                    fa(idir) = SUM(matrix_p_real(1:na, 1:nb)* &
                                                   (+MATMUL(acint_cos(1:na, 1:np, 1 + idir), TRANSPOSE(bchint_cos(1:nb, 1:np, 1))) &
                                                   + MATMUL(acint_sin(1:na, 1:np, 1 + idir), TRANSPOSE(bchint_sin(1:nb, 1:np, 1)))))
                                 IF (found_imag) &
                                    fa(idir) = fa(idir) - sign_imag*SUM(matrix_p_imag(1:na, 1:nb)* &
                                                   (+MATMUL(acint_sin(1:na, 1:np, 1 + idir), TRANSPOSE(bchint_cos(1:nb, 1:np, 1))) &
                                                   - MATMUL(acint_cos(1:na, 1:np, 1 + idir), TRANSPOSE(bchint_sin(1:nb, 1:np, 1)))))
                                 ! For fb:
                                 IF (found_real) &
                                    fb(idir) = SUM(matrix_p_real(1:na, 1:nb)* &
                                                   (+MATMUL(achint_cos(1:na, 1:np, 1), TRANSPOSE(bcint_cos(1:nb, 1:np, 1 + idir))) &
                                                   + MATMUL(achint_sin(1:na, 1:np, 1), TRANSPOSE(bcint_sin(1:nb, 1:np, 1 + idir)))))
                                 IF (found_imag) &
                                    fb(idir) = fb(idir) - sign_imag*SUM(matrix_p_imag(1:na, 1:nb)* &
                                                   (-MATMUL(achint_cos(1:na, 1:np, 1), TRANSPOSE(bcint_sin(1:nb, 1:np, 1 + idir))) &
                                                   + MATMUL(achint_sin(1:na, 1:np, 1), TRANSPOSE(bcint_cos(1:nb, 1:np, 1 + idir)))))
                              ELSE
                                 ! For fa:
                                 IF (found_real) &
                                    fa(idir) = SUM(matrix_p_real(1:nb, 1:na)* &
                                                   (+MATMUL(bchint_cos(1:nb, 1:np, 1), TRANSPOSE(acint_cos(1:na, 1:np, 1 + idir))) &
                                                   + MATMUL(bchint_sin(1:nb, 1:np, 1), TRANSPOSE(acint_sin(1:na, 1:np, 1 + idir)))))
                                 IF (found_imag) &
                                    fa(idir) = fa(idir) - sign_imag*SUM(matrix_p_imag(1:nb, 1:na)* &
                                                   (+MATMUL(bchint_sin(1:nb, 1:np, 1), TRANSPOSE(acint_cos(1:na, 1:np, 1 + idir))) &
                                                   - MATMUL(bchint_cos(1:nb, 1:np, 1), TRANSPOSE(acint_sin(1:na, 1:np, 1 + idir)))))
                                 ! For fb
                                 IF (found_real) &
                                    fb(idir) = SUM(matrix_p_real(1:nb, 1:na)* &
                                                   (+MATMUL(bcint_cos(1:nb, 1:np, 1 + idir), TRANSPOSE(achint_cos(1:na, 1:np, 1))) &
                                                   + MATMUL(bcint_sin(1:nb, 1:np, 1 + idir), TRANSPOSE(achint_sin(1:na, 1:np, 1)))))
                                 IF (found_imag) &
                                    fb(idir) = fb(idir) - sign_imag*SUM(matrix_p_imag(1:nb, 1:na)* &
                                                   (-MATMUL(bcint_cos(1:nb, 1:np, 1 + idir), TRANSPOSE(achint_sin(1:na, 1:np, 1))) &
                                                   + MATMUL(bcint_sin(1:nb, 1:np, 1 + idir), TRANSPOSE(achint_cos(1:na, 1:np, 1)))))
                              END IF
                              force_thread(idir, iatom) = force_thread(idir, iatom) + f0*fa(idir)
                              force_thread(idir, katom) = force_thread(idir, katom) - f0*fa(idir)
                              force_thread(idir, jatom) = force_thread(idir, jatom) + f0*fb(idir)
                              force_thread(idir, katom) = force_thread(idir, katom) - f0*fb(idir)
                           END DO
                           EXIT
                        END IF
                     END DO
                  END DO
               END DO
            END IF

         END DO

!$OMP END PARALLEL

         ! Update the force
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
!$OMP DO
         DO iatom = 1, natom
            i = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%gth_ppnl(:, i) = force(ikind)%gth_ppnl(:, i) + force_thread(:, iatom)
         END DO
!$OMP END DO

         ! Clean up
         IF (SIZE(rho_ao) == 2) THEN
            CALL dbcsr_add(rho_ao(1)%matrix, rho_ao(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_add(rho_ao_im(1)%matrix, rho_ao_im(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         END IF
         CALL release_sap_int(sap_int_cos)
         CALL release_sap_int(sap_int_sin)

         DEALLOCATE (basis_set, atom_of_kind, kind_of)

      END IF

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_nl_force

END MODULE rt_propagation_velocity_gauge
