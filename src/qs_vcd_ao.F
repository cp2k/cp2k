!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_vcd_ao
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_kinetic,                      ONLY: kinetic
   USE ai_overlap_ppl,                  ONLY: ppl_integral
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_desymmetrize, dbcsr_distribution_get, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_init_p, &
        dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_work_create
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: coset,&
                                              init_orbital_pointers,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_integrate_potential,          ONLY: integrate_pgf_product
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: vcd_env_type
   USE qs_moments,                      ONLY: build_dsdv_moments
   USE qs_neighbor_list_types,          ONLY: &
        get_iterator_info, get_neighbor_list_set_p, neighbor_list_iterate, &
        neighbor_list_iterator_create, neighbor_list_iterator_p_type, &
        neighbor_list_iterator_release, neighbor_list_set_p_type, nl_set_sub_iterator, &
        nl_sub_iterate
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE realspace_grid_types,            ONLY: realspace_grid_type
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE sap_kind_types,                  ONLY: alist_type,&
                                              build_sap_ints,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE task_list_types,                 ONLY: atom_pair_type,&
                                              task_list_type,&
                                              task_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vcd_ao'
   INTEGER, PARAMETER :: bi_1 = 1, bi_x = 2, bi_y = 3, bi_z = 4, bi_xx = 5, &
                         bi_xy = 6, bi_xz = 7, bi_yy = 8, bi_yz = 9, bi_zz = 10

! *** Public subroutines ***

   PUBLIC :: build_dSdV_matrix, build_com_rpnl_r, &
             hr_mult_by_delta_3d, hr_mult_by_delta_1d, build_dcom_rpnl, build_drpnl_matrix, &
             build_tr_matrix, &
             build_rcore_matrix, &
             build_rpnl_matrix, &
             build_matrix_r_vhxc, &
             build_matrix_hr_rh

CONTAINS

! **************************************************************************************************
!> \brief Build the matrix Hr*delta_nu^\lambda - rH*delta_mu^\lambda
!> \param vcd_env ...
!> \param qs_env ...
!> \param rc ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE build_matrix_hr_rh(vcd_env, qs_env, rc)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(3)                             :: rc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_matrix_hr_rh'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: eps_ppnl
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      sab_all=sab_all, &
                      sap_ppnl=sap_ppnl, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell)

      eps_ppnl = dft_control%qs_control%eps_ppnl
      DO i = 1, 3
         CALL dbcsr_set(vcd_env%matrix_hr(ispin, i)%matrix, 0._dp)
         CALL dbcsr_set(vcd_env%matrix_rh(ispin, i)%matrix, 0._dp)
      END DO

      ASSOCIATE (my_matrix_hr_1d => vcd_env%matrix_hr(ispin, 1:3))
         CALL build_rpnl_matrix(my_matrix_hr_1d, qs_kind_set, particle_set, sab_all, sap_ppnl, &
                                dft_control%qs_control%eps_ppnl, cell, rc, &
                                direction_Or=.TRUE.)
         CALL build_tr_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, &
                              direction_Or=.TRUE., rc=rc)
         CALL build_rcore_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, rc)
         CALL build_matrix_r_vhxc(vcd_env%matrix_hr, qs_env, rc)
      END ASSOCIATE

      ASSOCIATE (my_matrix_hr_1d => vcd_env%matrix_rh(ispin, 1:3))
         CALL build_rpnl_matrix(my_matrix_hr_1d, qs_kind_set, particle_set, sab_all, sap_ppnl, &
                                dft_control%qs_control%eps_ppnl, cell, rc, &
                                direction_Or=.FALSE.)
         CALL build_tr_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, &
                              direction_Or=.FALSE., rc=rc)
         CALL build_rcore_matrix(my_matrix_hr_1d, qs_env, qs_kind_set, "ORB", sab_all, rc)
         CALL build_matrix_r_vhxc(vcd_env%matrix_rh, qs_env, rc)
      END ASSOCIATE

      CALL timestop(handle)
   END SUBROUTINE build_matrix_hr_rh

! **************************************************************************************************
!> \brief Product of r with V_nl. Adapted from build_com_rpnl.
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param cell ...
!> \param ref_point ...
!> \param direction_Or ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE build_rpnl_matrix(matrix_rv, qs_kind_set, particle_set, sab_all, sap_ppnl, eps_ppnl, &
                                cell, ref_point, direction_Or)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(cell_type), INTENT(IN), OPTIONAL, POINTER     :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: ref_point
      LOGICAL                                            :: direction_Or

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_rpnl_matrix'

      INTEGER                                            :: handle, i, iab, iac, iatom, ibc, icol, &
                                                            ikind, irow, jatom, jkind, kac, kbc, &
                                                            kkind, na, natom, nb, nkind, np, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, ppnl_present
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), DIMENSION(3)                   :: blocks_rv
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)
      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      MARK_USED(ref_point)
      ! "nder" in moment_mode is "order"
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=1, moment_mode=.TRUE., &
                          particle_set=particle_set, cell=cell, refpoint=ref_point)

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl, locks, sab_all, direction_Or) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, &
!$OMP          iab, irow, icol, blocks_rv, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, &
!$OMP          hash, natom, acint, bcint, achint, bchint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_all(1)%nl_size

         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         ! *** Create matrix blocks for a new matrix block column ***
         irow = iatom
         icol = jatom
         DO i = 1, 3
            CALL dbcsr_get_block_p(matrix_rv(i)%matrix, irow, icol, blocks_rv(i)%block, found)
            CPASSERT(found)
         END DO

         ! loop over all kinds for projector atom
         DO kkind = 1, nkind
            iac = ikind + nkind*(kkind - 1)
            ibc = jkind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
            IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
            CALL get_alist(sap_int(iac), alist_ac, iatom)
            CALL get_alist(sap_int(ibc), alist_bc, jatom)

            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
            IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
            DO kac = 1, alist_ac%nclist
               DO kbc = 1, alist_bc%nclist
                  IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                  IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                     IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                     acint => alist_ac%clist(kac)%acint
                     bcint => alist_bc%clist(kbc)%acint
                     achint => alist_ac%clist(kac)%achint
                     bchint => alist_bc%clist(kbc)%achint
                     na = SIZE(acint, 1)
                     np = SIZE(acint, 2)
                     nb = SIZE(bcint, 1)
!$                   hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                   CALL omp_set_lock(locks(hash))
                     IF (direction_Or) THEN
                        ! Vnl*r
                        blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 2)))
                        blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 3)))
                        blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, 4)))
                     ELSE
                        ! r*Vnl
                        blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1)))
                        blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 3), TRANSPOSE(bcint(1:nb, 1:np, 1)))
                        blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) + &
                                                         MATMUL(achint(1:na, 1:np, 4), TRANSPOSE(bcint(1:nb, 1:np, 1)))
                     END IF
!$                   CALL omp_unset_lock(locks(hash))
                     EXIT ! We have found a match and there can be only one single match
                  END IF
               END DO
            END DO
         END DO
         DO i = 1, 3
            NULLIFY (blocks_rv(i)%block)
         END DO
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_rpnl_matrix

! **************************************************************************************************
!> \brief   Calculation of the product Tr or rT over Cartesian Gaussian functions.
!> \param   matrix_tr ...
!> \param   qs_env ...
!> \param   qs_kind_set ...
!> \param   basis_type basis set to be used
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   direction_Or ...
!> \param   rc ...
!> \date    11.10.2010
!> \par     History
!>          Ported from qs_overlap, replaces code in build_core_hamiltonian
!>          Refactoring [07.2014] JGH
!>          Simplify options and use new kinetic energy integral routine
!>          Adapted from qs_kinetic [07.2016]
!>          Adapted from build_com_tr_matrix [2021] by ED
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE build_tr_matrix(matrix_tr, qs_env, qs_kind_set, basis_type, sab_nl, direction_Or, rc)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix_tr
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      LOGICAL                                            :: direction_Or
      REAL(KIND=dp), DIMENSION(3), OPTIONAL              :: rc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_tr_matrix'

      INTEGER                                            :: handle, iatom, icol, ikind, ir, irow, &
                                                            iset, jatom, jkind, jset, ldsab, ltab, &
                                                            natom, ncoa, ncob, nkind, nseta, &
                                                            nsetb, sgfa, sgfb, slot
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found, trans
      REAL(KIND=dp)                                      :: tab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: qab, tkab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: kab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kx_block, ky_block, kz_block, rpgfa, &
                                                            rpgfb, scon_a, scon_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: qs_cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash, hash1, hash2
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER, PARAMETER                                 :: nlock = 501

      MARK_USED(int_8)

      CALL timeset(routineN, handle)

      nkind = SIZE(qs_kind_set)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      ! prepare basis set
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      ! *** Allocate work storage ***
      ldsab = get_memory_usage(qs_kind_set, basis_type)

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=qs_cell, &
                      natom=natom)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (ldsab,do_symmetric, sab_nl, rc,&
!$OMP         ncoset,matrix_tr,basis_set_list,qs_cell,natom,locks)  &
!$OMP PRIVATE (kx_block,ky_block,kz_block,kab,qab,tab,ikind,jkind,iatom,jatom,rab,rac,rbc,cell, &
!$OMP          basis_set_a, basis_set_b, nseta, ncoa, ncob, ltab, nsetb, tkab, &
!$OMP          irow, icol, found, trans, sgfa, sgfb, iset, jset, &
!$OMP          hash, hash1, hash2, iatom8, slot) &
!$OMP PRIVATE (first_sgfa, la_max, la_min, npgfa, nsgfa, rpgfa, set_radius_a, zeta, scon_a) &
!$OMP PRIVATE (first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, set_radius_b, zetb, scon_b) &
!$OMP SHARED(particle_set, direction_or) &
!$OMP PRIVATE(ra, rb)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

      ALLOCATE (kab(ldsab, ldsab, 3), qab(ldsab, ldsab))

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sab_nl(1)%nl_size

         ikind = sab_nl(1)%nlist_task(slot)%ikind
         jkind = sab_nl(1)%nlist_task(slot)%jkind
         iatom = sab_nl(1)%nlist_task(slot)%iatom
         jatom = sab_nl(1)%nlist_task(slot)%jatom
         cell(:) = sab_nl(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_nl(1)%nlist_task(slot)%r(1:3)

         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE

!$       iatom8 = INT(iatom - 1, int_8)*INT(natom, int_8) + INT(jatom, int_8)
!$       hash1 = INT(MOD(iatom8, INT(nlock, int_8)) + 1)

         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         nseta = basis_set_a%nset
         nsetb = basis_set_b%nset

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF
         NULLIFY (kx_block)
         CALL dbcsr_get_block_p(matrix=matrix_tr(1)%matrix, &
                                row=irow, col=icol, BLOCK=kx_block, found=found)
         CPASSERT(found)
         NULLIFY (ky_block)
         CALL dbcsr_get_block_p(matrix=matrix_tr(2)%matrix, &
                                row=irow, col=icol, BLOCK=ky_block, found=found)
         CPASSERT(found)
         NULLIFY (kz_block)
         CALL dbcsr_get_block_p(matrix=matrix_tr(3)%matrix, &
                                row=irow, col=icol, BLOCK=kz_block, found=found)
         CPASSERT(found)

         ! The kinetic integrals depend only on rab (also for the screening)
         tab = NORM2(rab)

         ! With and without MIC, rab(:) is the vector giving us the coordinates of rb.
         ra = pbc(particle_set(iatom)%r(:), qs_cell)
         rb(:) = ra(:) + rab(:)
         rac = pbc(rc, ra, qs_cell)
         rbc = rac + rab

         trans = do_symmetric .AND. (iatom > jatom)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < tab) CYCLE

!$             hash2 = MOD((iset - 1)*nsetb + jset, nlock) + 1
!$             hash = MOD(hash1 + hash2, nlock) + 1

               ncob = npgfb(jset)*(ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1))
               sgfb = first_sgfb(1, jset)

               ! calculate integrals
               ltab = MAX(npgfa(iset)*ncoset(la_max(iset) + 1), npgfb(jset)*ncoset(lb_max(jset) + 1))
               ALLOCATE (tkab(ltab, ltab))
               CALL kinetic(la_max(iset) + 1, la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                            lb_max(jset) + 1, lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                            rab, tkab)
               ! commutator
               CALL ab_opr(la_max(iset), npgfa(iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), rpgfb(:, jset), lb_min(jset), &
                           tab, tkab, kab, rac, rbc, direction_Or)

               DEALLOCATE (tkab)
               ! Contraction step
               DO ir = 1, 3
                  CALL contraction(kab(:, :, ir), qab, ca=scon_a(:, sgfa:), na=ncoa, ma=nsgfa(iset), &
                                   cb=scon_b(:, sgfb:), nb=ncob, mb=nsgfb(jset), trans=trans)

!$                CALL omp_set_lock(locks(hash))
                  SELECT CASE (ir)
                  CASE (1)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), kx_block, sgfa, sgfb, trans=trans)
                  CASE (2)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), ky_block, sgfa, sgfb, trans=trans)
                  CASE (3)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), kz_block, sgfa, sgfb, trans=trans)
                  END SELECT
!$                CALL omp_unset_lock(locks(hash))
               END DO

            END DO
         END DO
      END DO !iterator
      DEALLOCATE (kab, qab)
!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      ! Release work storage
      DEALLOCATE (basis_set_list)
      CALL timestop(handle)

   END SUBROUTINE build_tr_matrix

! **************************************************************************************************
!> \brief Commutator of the of the local part of the pseudopotential with r
!> \param matrix_rcore ...
!> \param qs_env ...
!> \param qs_kind_set ...
!> \param basis_type ...
!> \param sab_nl ...
!> \param rf ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE build_rcore_matrix(matrix_rcore, qs_env, qs_kind_set, basis_type, sab_nl, rf)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix_rcore
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      CHARACTER(LEN=*)                                   :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      REAL(KIND=dp), DIMENSION(3), OPTIONAL              :: rf

      CHARACTER(len=*), PARAMETER :: routineN = 'build_rcore_matrix'
      INTEGER, PARAMETER                                 :: nexp_max = 30

      INTEGER :: atom_a, atom_b, handle, i, iatom, icol, idir, ikind, inode, irow, iset, jatom, &
         jkind, jset, katom, kkind, ldai, ldsab, maxco, maxder, maxl, maxlgto, maxlppl, maxnset, &
         maxsgf, mepos, n_local, ncoa, ncob, nder, nexp_lpot, nexp_ppl, nimages, nkind, nloc, &
         nseta, nsetb, nthread, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(1:10)                           :: nrloc
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, &
                                                            nct_lpot, npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(nexp_max)                       :: nct_ppl
      LOGICAL                                            :: do_symmetric, dokp, ecp_local, found, &
                                                            lpotextended
      REAL(KIND=dp)                                      :: alpha, dab, dac, dbc, ppl_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hab, qab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ppl_work, rhab, work
      REAL(KIND=dp), DIMENSION(1:10)                     :: aloc, bloc
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, raf, rb, rbc, rbf
      REAL(KIND=dp), DIMENSION(4, nexp_max)              :: cval_ppl
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_local, alpha_lpot, c_local, cexp_ppl, &
                                                            set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_lpot, hx_block, hy_block, hz_block, &
                                                            rpgfa, rpgfb, scon_a, scon_b, sphi_a, &
                                                            sphi_b, zeta, zetb
      REAL(KIND=dp), DIMENSION(nexp_max)                 :: alpha_ppl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ap_iterator, nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sac_ppl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      sac_ppl=sac_ppl, &
                      cell=cell)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      nkind = SIZE(qs_kind_set)

      ! prepare basis set
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      nder = 0
      nimages = 1

      alpha_ppl = 0
      nct_ppl = 0
      cval_ppl = 0

      nkind = SIZE(atomic_kind_set)

      dokp = (nimages > 1)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      maxder = ncoset(nder)

      CALL get_qs_kind_set(qs_kind_set, maxco=maxco, maxlgto=maxlgto, &
                           maxsgf=maxsgf, maxnset=maxnset, maxlppl=maxlppl, &
                           basis_type=basis_type)

      maxl = MAX(maxlgto, maxlppl)
      CALL init_orbital_pointers(2*maxl + 2*nder + 2)

      !tz: maxco in maxco*ncoset(maxlgto+1) is an overkill,
      !    properly there should be maxpgf*ncoset(maxlgto+1), but maxpgf is difficult to get
      ldsab = MAX(maxco, ncoset(maxlppl), maxsgf, maxlppl, maxco*ncoset(maxlgto + 1))
      ldai = ncoset(2*maxlgto + 2)

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      nthread = 1
!$    nthread = omp_get_max_threads()

      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

      ! iterator for basis/potential list
      CALL neighbor_list_iterator_create(ap_iterator, sac_ppl, search=.TRUE., nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, ap_iterator, basis_set_list, &
!$OMP          atomic_kind_set, qs_kind_set, particle_set, &
!$OMP          sab_nl, sac_ppl, nthread, ncoset, nkind, &
!$OMP          atom_of_kind, ldsab,  maxnset, maxder, &
!$OMP          maxlgto, nder, maxco, dokp, cell) &
!$OMP SHARED (matrix_rcore, rf) &
!$OMP PRIVATE (ikind, jkind, inode, iatom, jatom, rab, basis_set_a, basis_set_b, atom_a) &
!$OMP PRIVATE (atom_b) &
!$OMP PRIVATE (nsetb) &
!$OMP PRIVATE (dab, irow, icol, hx_block, hy_block, hz_block, found, iset, ncoa) &
!$OMP PRIVATE (sgfa, jset, ncob, sgfb, work, hab, rhab, kkind, nseta) &
!$OMP PRIVATE (gth_potential, sgp_potential, alpha, cexp_ppl, lpotextended) &
!$OMP PRIVATE (ppl_radius, nexp_lpot, nexp_ppl, alpha_ppl, alpha_lpot, nct_ppl) &
!$OMP PRIVATE (nct_lpot, cval_ppl, cval_lpot, rac, dac, rbc, dbc) &
!$OMP PRIVATE (mepos) &
!$OMP PRIVATE (katom, ppl_work, cellind, ecp_local) &
!$OMP PRIVATE (first_sgfa, la_max, la_min, npgfa, nsgfa, rpgfa, set_radius_a, sphi_a, zeta, scon_a) &
!$OMP PRIVATE (first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, set_radius_b, sphi_b, zetb, scon_b) &
!$OMP PRIVATE (nloc, nrloc, aloc, bloc, n_local, a_local, c_local, ldai) &
!$OMP PRIVATE (ra, rb, qab, raf, rbf)

      mepos = 0
!$    mepos = omp_get_thread_num()

      ALLOCATE (hab(ldsab, ldsab), rhab(ldsab, ldsab, 3), work(ldsab, ldsab*(nder + 1), 3))
      ALLOCATE (qab(ldsab, ldsab))

      ldai = ncoset(2*maxlgto + 2)
      ALLOCATE (ppl_work(ldai, ldai, MAX(maxder, 2*maxlgto + 2 + 1)))

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cellind)

         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE

         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         scon_a => basis_set_a%scon
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         scon_b => basis_set_b%scon

         nseta = basis_set_a%nset
         nsetb = basis_set_b%nset

         ! *** Create matrix blocks for a new matrix block column ***
         irow = iatom
         icol = jatom

         NULLIFY (hx_block)
         NULLIFY (hy_block)
         NULLIFY (hz_block)
         CALL dbcsr_get_block_p(matrix=matrix_rcore(1)%matrix, &
                                row=irow, col=icol, BLOCK=hx_block, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_rcore(2)%matrix, &
                                row=irow, col=icol, BLOCK=hy_block, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_rcore(3)%matrix, &
                                row=irow, col=icol, BLOCK=hz_block, found=found)
         CPASSERT(found)

         dab = NORM2(rab)
         ra = pbc(particle_set(iatom)%r(:), cell)
         rb(:) = ra(:) + rab(:)

         raf = pbc(rf, ra, cell)
         rbf = raf + rab

         ! loop over all kinds for pseudopotential atoms
         DO kkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(kkind), gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            IF (ASSOCIATED(gth_potential)) THEN
               CALL get_potential(potential=gth_potential, &
                                  alpha_ppl=alpha, cexp_ppl=cexp_ppl, &
                                  lpot_present=lpotextended, ppl_radius=ppl_radius)
               nexp_ppl = 1
               alpha_ppl(1) = alpha
               nct_ppl(1) = SIZE(cexp_ppl)
               cval_ppl(1:nct_ppl(1), 1) = cexp_ppl(1:nct_ppl(1))
               IF (lpotextended) THEN
                  CALL get_potential(potential=gth_potential, &
                                     nexp_lpot=nexp_lpot, alpha_lpot=alpha_lpot, nct_lpot=nct_lpot, cval_lpot=cval_lpot)
                  CPASSERT(nexp_lpot < nexp_max)
                  nexp_ppl = nexp_lpot + 1
                  alpha_ppl(2:nexp_lpot + 1) = alpha_lpot(1:nexp_lpot)
                  nct_ppl(2:nexp_lpot + 1) = nct_lpot(1:nexp_lpot)
                  DO i = 1, nexp_lpot
                     cval_ppl(1:nct_lpot(i), i + 1) = cval_lpot(1:nct_lpot(i), i)
                  END DO
               END IF
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, ppl_radius=ppl_radius)
               IF (ecp_local) THEN
                  CALL get_potential(potential=sgp_potential, nloc=nloc, nrloc=nrloc, aloc=aloc, bloc=bloc)
                  IF (SUM(ABS(aloc(1:nloc))) < 1.0e-12_dp) CYCLE
                  nexp_ppl = nloc
                  CPASSERT(nexp_ppl <= nexp_max)
                  nct_ppl(1:nloc) = nrloc(1:nloc) - 1
                  alpha_ppl(1:nloc) = bloc(1:nloc)
                  cval_ppl(1, 1:nloc) = aloc(1:nloc)
               ELSE
                  CALL get_potential(potential=sgp_potential, n_local=n_local, a_local=a_local, c_local=c_local)
                  nexp_ppl = n_local
                  CPASSERT(nexp_ppl <= nexp_max)
                  nct_ppl(1:n_local) = 1
                  alpha_ppl(1:n_local) = a_local(1:n_local)
                  cval_ppl(1, 1:n_local) = c_local(1:n_local)
               END IF
            ELSE
               CYCLE
            END IF

            CALL nl_set_sub_iterator(ap_iterator, ikind, kkind, iatom, mepos=mepos)

            DO WHILE (nl_sub_iterate(ap_iterator, mepos=mepos) == 0)

               CALL get_iterator_info(ap_iterator, mepos=mepos, jatom=katom, r=rac)
               dac = SQRT(SUM(rac*rac))
               rbc(:) = rac(:) - rab(:)
               dbc = SQRT(SUM(rbc*rbc))
               IF ((MAXVAL(set_radius_a(:)) + ppl_radius < dac) .OR. &
                   (MAXVAL(set_radius_b(:)) + ppl_radius < dbc)) THEN
                  CYCLE
               END IF
               DO iset = 1, nseta
                  IF (set_radius_a(iset) + ppl_radius < dac) CYCLE
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  !                  ncoa = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
                  sgfa = first_sgfa(1, iset)
                  DO jset = 1, nsetb
                     IF (set_radius_b(jset) + ppl_radius < dbc) CYCLE
                     ncob = npgfb(jset)*ncoset(lb_max(jset))
                     !                     ncob = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
                     sgfb = first_sgfb(1, jset)
                     IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                     ! *** Calculate the GTH pseudo potential forces ***
                     hab = 0
                     rhab = 0
                     ppl_work = 0
                     work = 0

                     CALL ppl_integral( &
                        la_max(iset) + 1, la_min(iset), npgfa(iset), &
                        rpgfa(:, iset), zeta(:, iset), &
                        lb_max(jset) + 1, lb_min(jset), npgfb(jset), &
                        rpgfb(:, jset), zetb(:, jset), &
                        nexp_ppl, alpha_ppl, nct_ppl, cval_ppl, ppl_radius, &
                        rab, dab, rac, dac, rbc, dbc, hab(:, :), ppl_work)

                     ! product with r
                     CALL ab_opr(la_max(iset), npgfa(iset), rpgfa(:, iset), 0, &
                                 lb_max(jset), npgfb(jset), rpgfb(:, jset), 0, &
                                 dab, hab(:, :), rhab(:, :, :), raf, rbf, &
                                 direction_Or=.FALSE.)

                     DO idir = 1, 3
                        CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                                   1.0_dp, rhab(1, 1, idir), SIZE(rhab, 1), &
                                   sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                                   0.0_dp, work(1, 1, idir), SIZE(work, 1))
!$OMP CRITICAL(h_block_critical)
                        SELECT CASE (idir)
                        CASE (1)
                           CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                      1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                      work(1, 1, idir), SIZE(work, 1), &
                                      1.0_dp, hx_block(sgfa, sgfb), SIZE(hx_block, 1))
                        CASE (2)
                           CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                      1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                      work(1, 1, idir), SIZE(work, 1), &
                                      1.0_dp, hy_block(sgfa, sgfb), SIZE(hy_block, 1))
                        CASE (3)
                           CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                      1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                      work(1, 1, idir), SIZE(work, 1), &
                                      1.0_dp, hz_block(sgfa, sgfb), SIZE(hz_block, 1))
                        END SELECT
!$OMP END CRITICAL(h_block_critical)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO ! iterator

      DEALLOCATE (hab, rhab, work, ppl_work)

!$OMP END PARALLEL

      CALL neighbor_list_iterator_release(ap_iterator)
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (atom_of_kind, basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_rcore_matrix

! **************************************************************************************************
!> \brief Commutator of the Hartree+XC potentials with r
!> \param matrix_rv ...
!> \param qs_env ...
!> \param rc ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE build_matrix_r_vhxc(matrix_rv, qs_env, rc)
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), POINTER                          :: matrix_rv
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3)                        :: rc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_matrix_r_vhxc'
      INTEGER, PARAMETER                                 :: nspins = 1

      INTEGER                                            :: handle, idir, ispin
      REAL(kind=dp)                                      :: edisp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_rvxc, matrix_rvxc_desymm
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), DIMENSION(:), POINTER               :: v_rspace, v_tau_rspace
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, &
                      ks_env=ks_env, &
                      pw_env=pw_env, &
                      input=input, &
                      v_hartree_rspace=v_hartree_rspace, &
                      energy=energy)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! matrix_rv(nspin, 3) was allocated as
      ! CALL dbcsr_init_p(vcd_env%matrix_hr(ispin, i)%matrix)
      ! CALL dbcsr_copy(vcd_env%matrix_hr(ispin, i)%matrix, vcd_env%matrix_nosym_temp(1)%matrix)

      NULLIFY (matrix_rvxc, matrix_rvxc_desymm)
      CALL dbcsr_allocate_matrix_set(matrix_rvxc, nspins, 3)
      CALL dbcsr_allocate_matrix_set(matrix_rvxc_desymm, nspins, 3)

      DO ispin = 1, nspins
         DO idir = 1, 3
            CALL dbcsr_init_p(matrix_rvxc(ispin, idir)%matrix)
            CALL dbcsr_init_p(matrix_rvxc_desymm(ispin, idir)%matrix)

            CALL dbcsr_copy(matrix_rvxc_desymm(ispin, idir)%matrix, matrix_rv(1, 1)%matrix)
            CALL dbcsr_set(matrix_rvxc_desymm(ispin, idir)%matrix, 0._dp)

            CALL dbcsr_copy(matrix_rvxc(ispin, idir)%matrix, matrix_ks(ispin)%matrix)
            CALL dbcsr_set(matrix_rvxc(ispin, idir)%matrix, 0.0_dp)
         END DO
      END DO

      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL get_qs_env(qs_env=qs_env, rho=rho_struct)

      NULLIFY (v_rspace)
      NULLIFY (v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=energy%exc, &
                         edisp=edisp, dispersion_env=qs_env%dispersion_env, &
                         just_energy=.FALSE.)

      IF (.NOT. ASSOCIATED(v_rspace)) THEN
         ALLOCATE (v_rspace(nspins))
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%create_pw(v_rspace(ispin), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(v_rspace(ispin))
         END DO
      END IF

      DO ispin = 1, nspins
         CALL pw_axpy(v_hartree_rspace, v_rspace(ispin), 1.0_dp/v_hartree_rspace%pw_grid%dvol)
         CALL integrate_rv_rspace(v_rspace=v_rspace(ispin), hmat=matrix_rvxc(ispin, :), qs_env=qs_env, &
                                  rc=rc)

         DO idir = 1, 3
            CALL dbcsr_scale(matrix_rvxc(ispin, idir)%matrix, v_rspace(ispin)%pw_grid%dvol)
            CALL dbcsr_desymmetrize(matrix_rvxc(ispin, idir)%matrix, matrix_rvxc_desymm(ispin, idir)%matrix)
            CALL dbcsr_add(matrix_rv(ispin, idir)%matrix, matrix_rvxc_desymm(ispin, idir)%matrix, 1.0_dp, 1.0_dp)
         END DO
      END DO

      ! return pw grids
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(v_rspace(ispin))
      END DO
      DEALLOCATE (v_rspace)

      CALL dbcsr_deallocate_matrix_set(matrix_rvxc)
      CALL dbcsr_deallocate_matrix_set(matrix_rvxc_desymm)

      CALL timestop(handle)

   END SUBROUTINE build_matrix_r_vhxc

! **************************************************************************************************
!> \brief Calculates the integrals < mu | r * V | nu >
!>        There is no direction_Or argument, because the potentials commute with r
!>        This routine uses integrate_pgf_product directly. It could probably be rewritten to use
!>          the new task_list interface.
!> \param v_rspace ...
!> \param hmat ...
!> \param qs_env ...
!> \param rc ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE integrate_rv_rspace(v_rspace, hmat, qs_env, rc)

      TYPE(pw_type)                                      :: v_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: hmat
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3)                        :: rc

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_rv_rspace'

      CHARACTER(len=default_string_length)               :: my_basis_type
      INTEGER :: bcol, brow, handle, iatom, idir, igrid_level, ikind, ikind_old, ilevel, img, &
         ipair, ipgf, ipgf_new, iset, iset_new, iset_old, itask, ithread, jatom, jkind, jkind_old, &
         jpgf, jpgf_new, jset, jset_new, jset_old, ldsab, maxco, maxlgto, maxpgf, maxset, &
         maxsgf_set, na1, na2, natom, nb1, nb2, ncoa, ncoa_full, ncob, ncob_full, nkind, nseta, &
         nsetb, nthread, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL :: atom_pair_changed, atom_pair_done, distributed_grids, found, has_threads, &
         my_compute_tau, my_gapw, new_set_pair_coming
      REAL(KIND=dp)                                      :: dab, eps_rho_rspace, f, prefactor, &
                                                            radius, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rab_inv, rac, rb, rbc, rp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, rpgfa, rpgfb, sphi_a, sphi_b, work, &
                                                            zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: habt, rhab, workt
      TYPE(atom_pair_type), DIMENSION(:), POINTER        :: atom_pair_recv, atom_pair_send
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: h_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(mp_comm_type)                                 :: group
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), DIMENSION(:), POINTER   :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(task_type), DIMENSION(:), POINTER             :: tasks
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      my_compute_tau = .FALSE.
      my_gapw = .FALSE.
      my_basis_type = "ORB"

      ! get the task lists
      CALL get_qs_env(qs_env=qs_env, &
                      task_list=task_list, &
                      task_list_soft=task_list_soft)
      CPASSERT(ASSOCIATED(task_list))

      ! the information on the grids is provided through pw_env
      ! pw_env has to be the parent env for the potential grid (input)
      ! there is an option to provide an external grid
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      ! get all the general information on the system we are working on
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      virial=virial, &
                      natom=natom)

      rab = 0._dp
      rac = 0._dp
      rbc = 0._dp

      ! short cuts to task list variables
      tasks => task_list%tasks
      atom_pair_send => task_list%atom_pair_send
      atom_pair_recv => task_list%atom_pair_recv

      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_grids=rs_v)

      ! get mpi group from rs_v
      group = rs_v(1)%desc%group

      ! assign from pw_env
      gridlevel_info => pw_env%gridlevel_info

      ! transform the potential on the rs_multigrids
      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      nkind = SIZE(qs_kind_set)

      ! needs to be consistent with rho_rspace
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxlgto=maxlgto, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)

      distributed_grids = .FALSE.
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF (rs_v(igrid_level)%desc%distributed) THEN
            distributed_grids = .TRUE.
         END IF
      END DO

      nthread = 1
!$    nthread = omp_get_max_threads()

      ! get maximum numbers
      maxset = 0
      maxpgf = 0
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=orb_basis_set, basis_type=my_basis_type)

         IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE

         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                npgf=npgfa, nset=nseta)

         maxset = MAX(nseta, maxset)
         maxpgf = MAX(MAXVAL(npgfa), maxpgf)
      END DO

      ldsab = MAX(maxco, maxsgf_set, maxpgf*ncoset(maxlgto + 1))

      !   *** Allocate work storage ***
      NULLIFY (habt, workt)
      CALL reallocate(habt, 1, ldsab, 1, ldsab, 0, nthread)
      CALL reallocate(workt, 1, ldsab, 1, maxsgf_set, 0, nthread)
      ALLOCATE (rhab(ldsab, ldsab, 3))

      ALLOCATE (h_block(3))

      ithread = 0
!$    ithread = omp_get_thread_num()
      work => workt(:, :, ithread)
      hab => habt(:, :, ithread)
      hab(:, :) = 0._dp

      iset_old = -1; jset_old = -1
      ikind_old = -1; jkind_old = -1

      ! Here we loop over gridlevels first, finalising the matrix after each grid level is
      ! completed.  On each grid level, we loop over atom pairs, which will only access
      ! a single block of each matrix, so with OpenMP, each matrix block is only touched
      ! by a single thread for each grid level
      loop_gridlevels: DO igrid_level = 1, gridlevel_info%ngrid_levels
         DO idir = 1, 3
            CALL dbcsr_work_create(hmat(idir)%matrix, work_mutable=.TRUE., n=nthread)
            CALL dbcsr_get_info(hmat(idir)%matrix, distribution=dist)
            CALL dbcsr_distribution_get(dist, has_threads=has_threads)
!$          IF (.NOT. has_threads) &
!$             CPABORT("No thread distribution defined.")
         END DO

         loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
         loop_tasks: DO itask = task_list%taskstart(ipair, igrid_level), task_list%taskstop(ipair, igrid_level)
            ilevel = tasks(itask)%grid_level
            img = tasks(itask)%image
            iatom = tasks(itask)%iatom
            jatom = tasks(itask)%jatom
            iset = tasks(itask)%iset
            jset = tasks(itask)%jset
            ipgf = tasks(itask)%ipgf
            jpgf = tasks(itask)%jpgf
            CPASSERT(img == 1)

            ! At the start of a block of tasks, get atom data (and kind data, if needed)
            IF (itask .EQ. task_list%taskstart(ipair, igrid_level)) THEN

               ikind = particle_set(iatom)%atomic_kind%kind_number
               jkind = particle_set(jatom)%atomic_kind%kind_number

               IF (iatom <= jatom) THEN
                  brow = iatom
                  bcol = jatom
               ELSE
                  brow = jatom
                  bcol = iatom
               END IF

               IF (ikind .NE. ikind_old) THEN
                  CALL get_qs_kind(qs_kind_set(ikind), &
                                   basis_set=orb_basis_set, basis_type=my_basis_type)

                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfa, &
                                         lmax=la_max, &
                                         lmin=la_min, &
                                         npgf=npgfa, &
                                         nset=nseta, &
                                         nsgf_set=nsgfa, &
                                         pgf_radius=rpgfa, &
                                         set_radius=set_radius_a, &
                                         sphi=sphi_a, &
                                         zet=zeta)
               END IF

               IF (jkind .NE. jkind_old) THEN
                  CALL get_qs_kind(qs_kind_set(jkind), &
                                   basis_set=orb_basis_set, basis_type=my_basis_type)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfb, &
                                         lmax=lb_max, &
                                         lmin=lb_min, &
                                         npgf=npgfb, &
                                         nset=nsetb, &
                                         nsgf_set=nsgfb, &
                                         pgf_radius=rpgfb, &
                                         set_radius=set_radius_b, &
                                         sphi=sphi_b, &
                                         zet=zetb)

               END IF

               DO idir = 1, 3
                  NULLIFY (h_block(idir)%block)
                  CALL dbcsr_get_block_p(hmat(idir)%matrix, brow, bcol, h_block(idir)%block, found)
                  CPASSERT(found)
               END DO

               ikind_old = ikind
               jkind_old = jkind

               atom_pair_changed = .TRUE.

            ELSE

               atom_pair_changed = .FALSE.

            END IF

            IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN
               ! We reuse the hab(:, :) array to put the new integrals in.

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               ncoa_full = npgfa(iset)*ncoset(la_max(iset) + 1)
               sgfa = first_sgfa(1, iset)
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               ncob_full = npgfb(jset)*ncoset(lb_max(jset) + 1)
               sgfb = first_sgfb(1, jset)

               IF (iatom <= jatom) THEN
                  hab(1:ncoa_full, 1:ncob_full) = 0._dp
               ELSE
                  hab(1:ncob_full, 1:ncoa_full) = 0._dp
               END IF

               iset_old = iset
               jset_old = jset

            END IF

            rab = tasks(itask)%rab
            dab = norm2(rab)
            ! With and without MIC, rab(:) is the vector giving us the coordinates of rb.
            ra = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra(:) + rab(:)
            rac = pbc(rc, ra, cell)
            rbc = rac + rab

            zetp = zeta(ipgf, iset) + zetb(jpgf, jset)
            f = zetb(jpgf, jset)/zetp
            rp(:) = ra(:) + f*rab(:)

            prefactor = EXP(-zeta(ipgf, iset)*f*DOT_PRODUCT(rab, rab))
            radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                              lb_min=lb_min(jset), lb_max=lb_max(jset), &
                                              ra=ra, rb=rb, rp=rp, &
                                              zetp=zetp, eps=eps_rho_rspace, &
                                              prefactor=prefactor, cutoff=1.0_dp)

            na1 = (ipgf - 1)*ncoset(la_max(iset) + 1) + 1
            na2 = ipgf*ncoset(la_max(iset) + 1)
            nb1 = (jpgf - 1)*ncoset(lb_max(jset) + 1) + 1
            nb2 = jpgf*ncoset(lb_max(jset) + 1)

            IF (iatom <= jatom) THEN
               CALL integrate_pgf_product( &
                  la_max(iset) + 1, zeta(ipgf, iset), la_min(iset), &
                  lb_max(jset) + 1, zetb(jpgf, jset), lb_min(jset), &
                  ra, rab, rs_v(igrid_level), &
                  hab, o1=na1 - 1, o2=nb1 - 1, &
                  radius=radius, &
                  calculate_forces=.FALSE.)
            ELSE
               rab_inv = -rab
               CALL integrate_pgf_product( &
                  lb_max(jset) + 1, zetb(jpgf, jset), lb_min(jset), &
                  la_max(iset) + 1, zeta(ipgf, iset), la_min(iset), &
                  rb, rab_inv, rs_v(igrid_level), &
                  hab, o1=nb1 - 1, o2=na1 - 1, &
                  radius=radius, &
                  calculate_forces=.FALSE.)
            END IF

            new_set_pair_coming = .FALSE.
            atom_pair_done = .FALSE.
            IF (itask < task_list%taskstop(ipair, igrid_level)) THEN
               ilevel = tasks(itask + 1)%grid_level
               img = tasks(itask + 1)%image
               iatom = tasks(itask + 1)%iatom
               jatom = tasks(itask + 1)%jatom
               iset_new = tasks(itask + 1)%iset
               jset_new = tasks(itask + 1)%jset
               ipgf_new = tasks(itask + 1)%ipgf
               jpgf_new = tasks(itask + 1)%jpgf
               IF (iset_new .NE. iset .OR. jset_new .NE. jset) THEN
                  new_set_pair_coming = .TRUE.
               END IF
            ELSE
               ! do not forget the last block
               new_set_pair_coming = .TRUE.
               atom_pair_done = .TRUE.
            END IF

            IF (new_set_pair_coming) THEN
               ! Increase lx, ly, lz by one to account for the | r * b >
               IF (iatom <= jatom) THEN
                  ! direction_Or = .false. so that we use rac
                  CALL ab_opr(la_max(iset), npgfa(iset), rpgfa(:, iset), 0, &
                              lb_max(jset), npgfb(jset), rpgfb(:, jset), 0, &
                              dab, hab(:, :), rhab(:, :, :), rac, rbc, direction_Or=.FALSE.)

               ELSE
                  ! direction_Or = .true. so that we use rac
                  CALL ab_opr(lb_max(jset), npgfb(jset), rpgfb(:, jset), 0, &
                              la_max(iset), npgfa(iset), rpgfa(:, iset), 0, &
                              dab, hab(:, :), rhab(:, :, :), rbc, rac, direction_Or=.TRUE.)
               END IF

               ! contract the block into h if we're done with the current set pair
               DO idir = 1, 3
                  IF (iatom <= jatom) THEN
                     work = 0._dp
                     work(1:ncoa, 1:nsgfb(jset)) = MATMUL(rhab(1:ncoa, 1:ncob, idir), sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1))
                     h_block(idir)%block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) = &
                        h_block(idir)%block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) + &
                        MATMUL(TRANSPOSE(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1)), work(1:ncoa, 1:nsgfb(jset)))
                  ELSE
                     work(1:ncob, 1:nsgfa(iset)) = MATMUL(rhab(1:ncob, 1:ncoa, idir), sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1))
                     h_block(idir)%block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) = &
                        h_block(idir)%block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) + &
                        MATMUL(TRANSPOSE(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1)), work(1:ncob, 1:nsgfa(iset)))
                  END IF
               END DO
            END IF

         END DO loop_tasks
         END DO loop_pairs

         DO idir = 1, 3
            CALL dbcsr_finalize(hmat(idir)%matrix)
         END DO

      END DO loop_gridlevels

      !   *** Release work storage ***
      DEALLOCATE (habt, rhab, workt, h_block)

      CALL timestop(handle)

   END SUBROUTINE integrate_rv_rspace

! **************************************************************************************************
!> \brief Builds the overlap derivative wrt nuclear velocities
!>         dS/dV = < mu | r | nu > * (nu - mu)
!> \param qs_env ...
!> \param matrix_dsdv ...
!> \param deltaR ...
!> \param rcc ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE build_dSdV_matrix(qs_env, matrix_dsdv, deltaR, rcc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix_dsdv
      REAL(KIND=dp), DIMENSION(:, :)                     :: deltaR
      REAL(KIND=dp), DIMENSION(3)                        :: rcc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_dSdV_matrix'

      INTEGER                                            :: handle, i
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, my_matrix_dsdv, &
                                                            my_matrix_dsdv2, my_matrix_mom
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (my_matrix_mom, my_matrix_dsdv, my_matrix_dsdv2, sab_all, qs_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      sab_all=sab_all, &
                      qs_kind_set=qs_kind_set, &
                      matrix_ks=matrix_ks)

      ! my_matrix_mom is needed because build_local_moment only works correctly with symmetric matrices
      CALL dbcsr_allocate_matrix_set(my_matrix_dsdv, 3)
      CALL dbcsr_allocate_matrix_set(my_matrix_dsdv2, 3)
      CALL dbcsr_allocate_matrix_set(my_matrix_mom, 3)

      DO i = 1, 3
         ALLOCATE (my_matrix_dsdv(i)%matrix)
         ALLOCATE (my_matrix_dsdv2(i)%matrix)
         ALLOCATE (my_matrix_mom(i)%matrix)

         CALL dbcsr_copy(my_matrix_dsdv(i)%matrix, matrix_dsdv(i)%matrix)
         CALL dbcsr_copy(my_matrix_dsdv2(i)%matrix, matrix_dsdv(i)%matrix)
         CALL dbcsr_copy(my_matrix_mom(i)%matrix, matrix_ks(1)%matrix)

         CALL dbcsr_set(my_matrix_dsdv2(i)%matrix, 0.0_dp)
         CALL dbcsr_set(my_matrix_dsdv(i)%matrix, 0.0_dp)
         CALL dbcsr_set(my_matrix_mom(i)%matrix, 0.0_dp)
         CALL dbcsr_set(matrix_dsdv(i)%matrix, 0.0_dp)
      END DO

      CALL build_dsdv_moments(qs_env, my_matrix_mom, 1, rcc)

      DO i = 1, 3
         CALL dbcsr_desymmetrize(my_matrix_mom(i)%matrix, my_matrix_dsdv(i)%matrix)
         CALL dbcsr_copy(my_matrix_dsdv2(i)%matrix, my_matrix_dsdv(i)%matrix)
      END DO

      ! delta_nu^A <mu|r|nu>
      CALL hr_mult_by_delta_3d(my_matrix_dsdv, qs_kind_set, "ORB", sab_all, &
                               deltaR, direction_Or=.TRUE.)
      ! -delta_mu^A <mu|r|nu>
      CALL hr_mult_by_delta_3d(my_matrix_dsdv2, qs_kind_set, "ORB", sab_all, &
                               deltaR, direction_Or=.FALSE.)
      DO i = 1, 3
         CALL dbcsr_copy(matrix_dsdv(i)%matrix, my_matrix_dsdv(i)%matrix)
         CALL dbcsr_add(matrix_dsdv(i)%matrix, my_matrix_dsdv2(i)%matrix, 1.0_dp, -1.0_dp)
      END DO

      CALL dbcsr_deallocate_matrix_set(my_matrix_dsdv)
      CALL dbcsr_deallocate_matrix_set(my_matrix_dsdv2)
      CALL dbcsr_deallocate_matrix_set(my_matrix_mom)

      CALL timestop(handle)

   END SUBROUTINE build_dSdV_matrix

! **************************************************************************************************
!> \brief Builds the  [Vnl, r] * r from either side
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param cell ...
!> \param direction_Or ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE build_com_rpnl_r(matrix_rv, qs_kind_set, sab_all, sap_ppnl, eps_ppnl, &
                               particle_set, cell, direction_Or)

      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: direction_Or

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_com_rpnl_r'

      INTEGER                                            :: handle, i, iab, iac, iatom, ibc, icol, &
                                                            ikind, irow, j, jatom, jkind, kac, &
                                                            kbc, kkind, na, natom, nb, nkind, np, &
                                                            slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), DIMENSION(3, 3)                :: blocks_rvr
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)
      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      ! We put zero as a reference point, because actually in the integrals we need two different ones:
      !  < a | (r - R^\lambda_\beta) * [V, r_\alpha - R^\eta_\alpha] | b >
      !  The first reference point can be added in a seperate step as the term will be
      !        - R^\lambda_\beta * < a | [V, r_\alpha] | b >
      !      = + R^\lambda_\beta * < a | [r_\alpha, V] | b >
      !  The second reference point is not important, because it disappears in the commutator
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=2, moment_mode=.TRUE., &
                          particle_set=particle_set, cell=cell, refpoint=[0._dp, 0._dp, 0._dp])

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, direction_Or, &
!$OMP          sap_int, nkind, eps_ppnl, locks, sab_all) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, blocks_rvr, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, j, &
!$OMP          hash, natom, acint, bcint, achint, bchint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_all(1)%nl_size

         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_all(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         ! *** Create matrix blocks for a new matrix block column ***
         irow = iatom
         icol = jatom
         DO i = 1, 3
            DO j = 1, 3
               ! t_alpha = MOD(i - 1, 3) + 1
               ! t_beta = FLOOR(REAL(i - 1, dp)/3._dp) + 1

               CALL dbcsr_get_block_p(matrix_rv(i, j)%matrix, irow, icol, blocks_rvr(i, j)%block, found)
               CPASSERT(found)
            END DO
         END DO

         ! loop over all kinds for projector atom
         DO kkind = 1, nkind
            iac = ikind + nkind*(kkind - 1)
            ibc = jkind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
            IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
            CALL get_alist(sap_int(iac), alist_ac, iatom)
            CALL get_alist(sap_int(ibc), alist_bc, jatom)
            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
            IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
            DO kac = 1, alist_ac%nclist
               DO kbc = 1, alist_bc%nclist
                  IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE

                  ! Some documention, considering the C1 O1 O2 molecule
                  ! The integrals are <a|p>
                  ! sap_int(1:2, 1:2) -> [(a=C, p=C), (a=C, p=O), (a=O, p=C), (a=O, p=0)]
                  ! the (a=O, p=C) entry has an alist
                  !  alist has two elements: O1, O2
                  !     alist(O1) -> clist(C1)%acint are the integrals
                  !     alist(O2) -> clist(C1)%acint are the integrals

                  IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                     IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                     acint => alist_ac%clist(kac)%acint
                     bcint => alist_bc%clist(kbc)%acint
                     achint => alist_ac%clist(kac)%achint
                     bchint => alist_bc%clist(kbc)%achint
                     na = SIZE(acint, 1)
                     np = SIZE(acint, 2)
                     nb = SIZE(bcint, 1)
!$                   hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                   CALL omp_set_lock(locks(hash))

                     ! Template:
                     ! blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                     !                                  MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1))) ! xV
                     ! 1, x, y, z, xx, xy, xz, yy, yz, zz
                     IF (direction_Or) THEN
                        ! (Vnl*r_alpha)*r_beta

                        ! This part is symmetric because [r_beta, r_alpha] = 0
                        ! matrix_rvr(x, gamma) += <mu|p>h_ij * <p|x*gamma|nu>
                        blocks_rvr(1, 1)%block(1:na, 1:nb) = blocks_rvr(1, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xx)))
                        blocks_rvr(1, 2)%block(1:na, 1:nb) = blocks_rvr(1, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xy)))
                        blocks_rvr(1, 3)%block(1:na, 1:nb) = blocks_rvr(1, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xz)))

                        ! matrix_rvr(y, gamma) += <mu|p>h_ij * <p|y*gamma|nu>
                        blocks_rvr(2, 1)%block(1:na, 1:nb) = blocks_rvr(2, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xy)))
                        blocks_rvr(2, 2)%block(1:na, 1:nb) = blocks_rvr(2, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yy)))
                        blocks_rvr(2, 3)%block(1:na, 1:nb) = blocks_rvr(2, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yz)))

                        ! matrix_rvr(z, gamma) += <mu|p>h_ij * <p|z*gamma|nu>
                        blocks_rvr(3, 1)%block(1:na, 1:nb) = blocks_rvr(3, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xz)))
                        blocks_rvr(3, 2)%block(1:na, 1:nb) = blocks_rvr(3, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yz)))
                        blocks_rvr(3, 3)%block(1:na, 1:nb) = blocks_rvr(3, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_zz)))

                        ! -(r_alpha*Vnl)*r_beta
                        ! matrix_rvr(x, gamma) += <mu|gamma|p>h_ij * <p|x|nu>
                        blocks_rvr(1, 1)%block(1:na, 1:nb) = blocks_rvr(1, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                        blocks_rvr(1, 2)%block(1:na, 1:nb) = blocks_rvr(1, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                        blocks_rvr(1, 3)%block(1:na, 1:nb) = blocks_rvr(1, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        ! matrix_rvr(y, gamma) += <mu|gamma|p>h_ij * <p|y|nu>
                        blocks_rvr(2, 1)%block(1:na, 1:nb) = blocks_rvr(2, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                        blocks_rvr(2, 2)%block(1:na, 1:nb) = blocks_rvr(2, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                        blocks_rvr(2, 3)%block(1:na, 1:nb) = blocks_rvr(2, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        ! matrix_rvr(z, gamma) += <mu|gamma|p>h_ij * <p|z|nu>
                        blocks_rvr(3, 1)%block(1:na, 1:nb) = blocks_rvr(3, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                        blocks_rvr(3, 2)%block(1:na, 1:nb) = blocks_rvr(3, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                        blocks_rvr(3, 3)%block(1:na, 1:nb) = blocks_rvr(3, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        ! This means that we have
                        ! blocks_rvr(1:3, 1:3) = blocks_rvr(alpha, beta)
                     ELSE
                        ! r_beta*(Vnl*r_alpha)
                        blocks_rvr(1, 1)%block(1:na, 1:nb) = blocks_rvr(1, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                        blocks_rvr(1, 2)%block(1:na, 1:nb) = blocks_rvr(1, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                        blocks_rvr(1, 3)%block(1:na, 1:nb) = blocks_rvr(1, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rvr(2, 1)%block(1:na, 1:nb) = blocks_rvr(2, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                        blocks_rvr(2, 2)%block(1:na, 1:nb) = blocks_rvr(2, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                        blocks_rvr(2, 3)%block(1:na, 1:nb) = blocks_rvr(2, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rvr(3, 1)%block(1:na, 1:nb) = blocks_rvr(3, 1)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))
                        blocks_rvr(3, 2)%block(1:na, 1:nb) = blocks_rvr(3, 2)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))
                        blocks_rvr(3, 3)%block(1:na, 1:nb) = blocks_rvr(3, 3)%block(1:na, 1:nb) + &
                                                             MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        ! -r_beta*(r_alpha*Vnl)
                        blocks_rvr(1, 1)%block(1:na, 1:nb) = blocks_rvr(1, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_xx), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(1, 2)%block(1:na, 1:nb) = blocks_rvr(1, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_xy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(1, 3)%block(1:na, 1:nb) = blocks_rvr(1, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_xz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rvr(2, 1)%block(1:na, 1:nb) = blocks_rvr(2, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_xy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(2, 2)%block(1:na, 1:nb) = blocks_rvr(2, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_yy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(2, 3)%block(1:na, 1:nb) = blocks_rvr(2, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_yz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rvr(3, 1)%block(1:na, 1:nb) = blocks_rvr(3, 1)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_xz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(3, 2)%block(1:na, 1:nb) = blocks_rvr(3, 2)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_yz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                        blocks_rvr(3, 3)%block(1:na, 1:nb) = blocks_rvr(3, 3)%block(1:na, 1:nb) - &
                                                             MATMUL(achint(1:na, 1:np, bi_zz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                     END IF

!$                   CALL omp_unset_lock(locks(hash))
                     EXIT ! We have found a match and there can be only one single match
                  END IF
               END DO
            END DO
         END DO
         DO i = 1, 3
            NULLIFY (blocks_rvr(i, 1)%block)
            NULLIFY (blocks_rvr(i, 2)%block)
            NULLIFY (blocks_rvr(i, 3)%block)
         END DO
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_com_rpnl_r

! **************************************************************************************************
!> \brief Calculate the double commutator [[Vnl, r], r]
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param pseudoatom Only consider pseudopotentials on atom lambda
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE build_dcom_rpnl(matrix_rv, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, pseudoatom)

      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      INTEGER, INTENT(IN)                                :: pseudoatom

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_dcom_rpnl'

      INTEGER                                            :: handle, i, iab, iac, iatom, ibc, icol, &
                                                            ikind, irow, j, jatom, jkind, kac, &
                                                            kbc, kkind, na, natom, nb, nkind, np, &
                                                            slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), DIMENSION(3, 3)                :: blocks_rv
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)
      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      ! "nder" in moment_mode is "order"
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=2, moment_mode=.TRUE., &
                          particle_set=particle_set, pseudoatom=pseudoatom)

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl, locks, sab_orb) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, blocks_rv, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, j, &
!$OMP          hash, natom, acint, bcint, achint, bchint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_orb(1)%nl_size

         ikind = sab_orb(1)%nlist_task(slot)%ikind
         jkind = sab_orb(1)%nlist_task(slot)%jkind
         iatom = sab_orb(1)%nlist_task(slot)%iatom
         jatom = sab_orb(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_orb(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         ! *** Create matrix blocks for a new matrix block column ***
         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         DO i = 1, 3
            DO j = 1, 3
               CALL dbcsr_get_block_p(matrix_rv(i, j)%matrix, irow, icol, blocks_rv(i, j)%block, found)
               blocks_rv(i, j)%block = 0._dp
               CPASSERT(found)
            END DO
         END DO

         ! loop over all kinds for projector atom
         DO kkind = 1, nkind
            iac = ikind + nkind*(kkind - 1)
            ibc = jkind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
            IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
            CALL get_alist(sap_int(iac), alist_ac, iatom)
            CALL get_alist(sap_int(ibc), alist_bc, jatom)
            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
            IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
            DO kac = 1, alist_ac%nclist
               DO kbc = 1, alist_bc%nclist
                  IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                  IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                     IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                     acint => alist_ac%clist(kac)%acint
                     bcint => alist_bc%clist(kbc)%acint
                     achint => alist_ac%clist(kac)%achint
                     bchint => alist_bc%clist(kbc)%achint
                     na = SIZE(acint, 1)
                     np = SIZE(acint, 2)
                     nb = SIZE(bcint, 1)
!$                   hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                   CALL omp_set_lock(locks(hash))
                     ! Template:
                     ! blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                     !                                  MATMUL(achint(1:na, 1:np, 2), TRANSPOSE(bcint(1:nb, 1:np, 1))) ! xV
                     IF (iatom <= jatom) THEN
                        ! r_alpha*Vnl*r_beta
                        blocks_rv(1, 1)%block(1:na, 1:nb) = blocks_rv(1, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(1, 2)%block(1:na, 1:nb) = blocks_rv(1, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(1, 3)%block(1:na, 1:nb) = blocks_rv(1, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        blocks_rv(2, 1)%block(1:na, 1:nb) = blocks_rv(2, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(2, 2)%block(1:na, 1:nb) = blocks_rv(2, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(2, 3)%block(1:na, 1:nb) = blocks_rv(2, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        blocks_rv(3, 1)%block(1:na, 1:nb) = blocks_rv(3, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(3, 2)%block(1:na, 1:nb) = blocks_rv(3, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(3, 3)%block(1:na, 1:nb) = blocks_rv(3, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        ! -r_alpha*r_beta*Vnl
                        blocks_rv(1, 1)%block(1:na, 1:nb) = blocks_rv(1, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_xx), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(1, 2)%block(1:na, 1:nb) = blocks_rv(1, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_xy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(1, 3)%block(1:na, 1:nb) = blocks_rv(1, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_xz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(2, 1)%block(1:na, 1:nb) = blocks_rv(2, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_xy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(2, 2)%block(1:na, 1:nb) = blocks_rv(2, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_yy), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(2, 3)%block(1:na, 1:nb) = blocks_rv(2, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_yz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(3, 1)%block(1:na, 1:nb) = blocks_rv(3, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_xz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(3, 2)%block(1:na, 1:nb) = blocks_rv(3, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_yz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        blocks_rv(3, 3)%block(1:na, 1:nb) = blocks_rv(3, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_zz), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

                        ! -Vnl*r_beta*r_alpha
                        blocks_rv(1, 1)%block(1:na, 1:nb) = blocks_rv(1, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xx)))

                        blocks_rv(1, 2)%block(1:na, 1:nb) = blocks_rv(1, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xy)))

                        blocks_rv(1, 3)%block(1:na, 1:nb) = blocks_rv(1, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xz)))

                        blocks_rv(2, 1)%block(1:na, 1:nb) = blocks_rv(2, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xy)))

                        blocks_rv(2, 2)%block(1:na, 1:nb) = blocks_rv(2, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yy)))

                        blocks_rv(2, 3)%block(1:na, 1:nb) = blocks_rv(2, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yz)))

                        blocks_rv(3, 1)%block(1:na, 1:nb) = blocks_rv(3, 1)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_xz)))

                        blocks_rv(3, 2)%block(1:na, 1:nb) = blocks_rv(3, 2)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_yz)))

                        blocks_rv(3, 3)%block(1:na, 1:nb) = blocks_rv(3, 3)%block(1:na, 1:nb) - &
                                                            MATMUL(achint(1:na, 1:np, bi_1), TRANSPOSE(bcint(1:nb, 1:np, bi_zz)))

                        ! +r_beta*Vnl*r_alpha
                        blocks_rv(1, 1)%block(1:na, 1:nb) = blocks_rv(1, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(1, 2)%block(1:na, 1:nb) = blocks_rv(1, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(1, 3)%block(1:na, 1:nb) = blocks_rv(1, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))

                        blocks_rv(2, 1)%block(1:na, 1:nb) = blocks_rv(2, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(2, 2)%block(1:na, 1:nb) = blocks_rv(2, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(2, 3)%block(1:na, 1:nb) = blocks_rv(2, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))

                        blocks_rv(3, 1)%block(1:na, 1:nb) = blocks_rv(3, 1)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        blocks_rv(3, 2)%block(1:na, 1:nb) = blocks_rv(3, 2)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                        blocks_rv(3, 3)%block(1:na, 1:nb) = blocks_rv(3, 3)%block(1:na, 1:nb) + &
                                                            MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))
                     ELSE
                        ! r_alpha*Vnl*r_beta
                        blocks_rv(1, 1)%block(1:nb, 1:na) = blocks_rv(1, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(1, 2)%block(1:nb, 1:na) = blocks_rv(1, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(1, 3)%block(1:nb, 1:na) = blocks_rv(1, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                        blocks_rv(2, 1)%block(1:nb, 1:na) = blocks_rv(2, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(2, 2)%block(1:nb, 1:na) = blocks_rv(2, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(2, 3)%block(1:nb, 1:na) = blocks_rv(2, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                        blocks_rv(3, 1)%block(1:nb, 1:na) = blocks_rv(3, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(3, 2)%block(1:nb, 1:na) = blocks_rv(3, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(3, 3)%block(1:nb, 1:na) = blocks_rv(3, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                        ! -r_alpha*r_beta*Vnl
                        blocks_rv(1, 1)%block(1:nb, 1:na) = blocks_rv(1, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_xx), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(1, 2)%block(1:nb, 1:na) = blocks_rv(1, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_xy), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(1, 3)%block(1:nb, 1:na) = blocks_rv(1, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_xz), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(2, 1)%block(1:nb, 1:na) = blocks_rv(2, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_xy), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(2, 2)%block(1:nb, 1:na) = blocks_rv(2, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_yy), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(2, 3)%block(1:nb, 1:na) = blocks_rv(2, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_yz), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(3, 1)%block(1:nb, 1:na) = blocks_rv(3, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_xz), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(3, 2)%block(1:nb, 1:na) = blocks_rv(3, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_yz), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        blocks_rv(3, 3)%block(1:nb, 1:na) = blocks_rv(3, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_zz), TRANSPOSE(acint(1:na, 1:np, bi_1)))

                        ! -Vnl*r_beta*r_alpha
                        blocks_rv(1, 1)%block(1:nb, 1:na) = blocks_rv(1, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_xx)))

                        blocks_rv(1, 2)%block(1:nb, 1:na) = blocks_rv(1, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_xy)))

                        blocks_rv(1, 3)%block(1:nb, 1:na) = blocks_rv(1, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_xz)))

                        blocks_rv(2, 1)%block(1:nb, 1:na) = blocks_rv(2, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_xy)))

                        blocks_rv(2, 2)%block(1:nb, 1:na) = blocks_rv(2, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_yy)))

                        blocks_rv(2, 3)%block(1:nb, 1:na) = blocks_rv(2, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_yz)))

                        blocks_rv(3, 1)%block(1:nb, 1:na) = blocks_rv(3, 1)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_xz)))

                        blocks_rv(3, 2)%block(1:nb, 1:na) = blocks_rv(3, 2)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_yz)))

                        blocks_rv(3, 3)%block(1:nb, 1:na) = blocks_rv(3, 3)%block(1:nb, 1:na) - &
                                                            MATMUL(bchint(1:nb, 1:np, bi_1), TRANSPOSE(acint(1:na, 1:np, bi_zz)))

                        ! +r_beta*Vnl*r_alpha
                        blocks_rv(1, 1)%block(1:nb, 1:na) = blocks_rv(1, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(1, 2)%block(1:nb, 1:na) = blocks_rv(1, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(1, 3)%block(1:nb, 1:na) = blocks_rv(1, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_x)))

                        blocks_rv(2, 1)%block(1:nb, 1:na) = blocks_rv(2, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(2, 2)%block(1:nb, 1:na) = blocks_rv(2, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(2, 3)%block(1:nb, 1:na) = blocks_rv(2, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_y)))

                        blocks_rv(3, 1)%block(1:nb, 1:na) = blocks_rv(3, 1)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_x), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                        blocks_rv(3, 2)%block(1:nb, 1:na) = blocks_rv(3, 2)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_y), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                        blocks_rv(3, 3)%block(1:nb, 1:na) = blocks_rv(3, 3)%block(1:nb, 1:na) + &
                                                            MATMUL(bchint(1:nb, 1:np, bi_z), TRANSPOSE(acint(1:na, 1:np, bi_z)))

                     END IF
!$                   CALL omp_unset_lock(locks(hash))
                     EXIT ! We have found a match and there can be only one single match
                  END IF
               END DO
            END DO
         END DO
         DO i = 1, 3
            NULLIFY (blocks_rv(i, 1)%block)
            NULLIFY (blocks_rv(i, 2)%block)
            NULLIFY (blocks_rv(i, 3)%block)
         END DO
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_dcom_rpnl

! **************************************************************************************************
!> \brief dV_nl/dV. Adapted from build_com_rpnl.
!> \param matrix_rv ...
!> \param qs_kind_set ...
!> \param sab_all ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param particle_set ...
!> \param pseudoatom ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE build_drpnl_matrix(matrix_rv, qs_kind_set, sab_all, sap_ppnl, eps_ppnl, particle_set, pseudoatom)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix_rv
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      INTEGER                                            :: pseudoatom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_drpnl_matrix'

      INTEGER                                            :: handle, i, iab, iac, iatom, ibc, icol, &
                                                            ikind, irow, jatom, jkind, kac, kbc, &
                                                            kkind, na, natom, nb, nkind, np, slot
      INTEGER, DIMENSION(3)                              :: cell_b
      LOGICAL                                            :: found, ppnl_present
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(block_p_type), DIMENSION(3)                   :: blocks_rv
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

!$    INTEGER(kind=omp_lock_kind), &
!$       ALLOCATABLE, DIMENSION(:) :: locks
!$    INTEGER                                            :: lock_num, hash
!$    INTEGER, PARAMETER                                 :: nlock = 501

      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (.NOT. ppnl_present) RETURN

      CALL timeset(routineN, handle)
      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! sap_int needs to be shared as multiple threads need to access this
      NULLIFY (sap_int)
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      ! "nder" in moment_mode is "order"
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=1, moment_mode=.TRUE., pseudoatom=pseudoatom)

      ! *** Set up a sorting index
      CALL sap_sort(sap_int)

      ALLOCATE (basis_set(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
            NULLIFY (basis_set(ikind)%gto_basis_set)
         END IF
      END DO

      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the commutator matrix elements

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (basis_set, matrix_rv, &
!$OMP          sap_int, nkind, eps_ppnl, locks, sab_all) &
!$OMP PRIVATE (ikind, jkind, iatom, jatom, cell_b, rab, &
!$OMP          iab, irow, icol, blocks_rv, &
!$OMP          found, iac, ibc, alist_ac, alist_bc, &
!$OMP          na, np, nb, kkind, kac, kbc, i, &
!$OMP          hash, natom, acint, bcint, achint, bchint)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED)

      DO slot = 1, sab_all(1)%nl_size

         ikind = sab_all(1)%nlist_task(slot)%ikind
         jkind = sab_all(1)%nlist_task(slot)%jkind
         iatom = sab_all(1)%nlist_task(slot)%iatom
         jatom = sab_all(1)%nlist_task(slot)%jatom
         cell_b(:) = sab_all(1)%nlist_task(slot)%cell(:)
         rab(1:3) = sab_all(1)%nlist_task(slot)%r(1:3)

         IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)

         ! *** Create matrix blocks for a new matrix block column ***
         irow = iatom
         icol = jatom
         DO i = 1, 3
            CALL dbcsr_get_block_p(matrix_rv(i)%matrix, irow, icol, blocks_rv(i)%block, found)
            CPASSERT(found)
         END DO

         ! loop over all kinds for projector atom
         DO kkind = 1, nkind
            iac = ikind + nkind*(kkind - 1)
            ibc = jkind + nkind*(kkind - 1)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
            IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
            CALL get_alist(sap_int(iac), alist_ac, iatom)
            CALL get_alist(sap_int(ibc), alist_bc, jatom)

            IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
            IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
            DO kac = 1, alist_ac%nclist
               DO kbc = 1, alist_bc%nclist
                  IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                  IF (ALL(cell_b + alist_bc%clist(kbc)%cell - alist_ac%clist(kac)%cell == 0)) THEN
                     IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                     acint => alist_ac%clist(kac)%acint
                     bcint => alist_bc%clist(kbc)%acint
                     achint => alist_ac%clist(kac)%achint
                     bchint => alist_bc%clist(kbc)%achint
                     na = SIZE(acint, 1)
                     np = SIZE(acint, 2)
                     nb = SIZE(bcint, 1)
!$                   hash = MOD((iatom - 1)*natom + jatom, nlock) + 1
!$                   CALL omp_set_lock(locks(hash))
                     ! Vnl*r
                     blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) + &
                                                      MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, bi_x)))
                     blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) + &
                                                      MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, bi_y)))
                     blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) + &
                                                      MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, bi_z)))

                     ! r*Vnl
                     blocks_rv(1)%block(1:na, 1:nb) = blocks_rv(1)%block(1:na, 1:nb) - &
                                                      MATMUL(achint(1:na, 1:np, bi_x), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                     blocks_rv(2)%block(1:na, 1:nb) = blocks_rv(2)%block(1:na, 1:nb) - &
                                                      MATMUL(achint(1:na, 1:np, bi_y), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))
                     blocks_rv(3)%block(1:na, 1:nb) = blocks_rv(3)%block(1:na, 1:nb) - &
                                                      MATMUL(achint(1:na, 1:np, bi_z), TRANSPOSE(bcint(1:nb, 1:np, bi_1)))

!$                   CALL omp_unset_lock(locks(hash))
                     EXIT ! We have found a match and there can be only one single match
                  END IF
               END DO
            END DO
         END DO
         DO i = 1, 3
            NULLIFY (blocks_rv(i)%block)
         END DO
      END DO

!$OMP DO
!$    DO lock_num = 1, nlock
!$       call omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO

!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT

!$OMP END PARALLEL

      CALL release_sap_int(sap_int)

      DEALLOCATE (basis_set)

      CALL timestop(handle)

   END SUBROUTINE build_drpnl_matrix

! **************************************************************************************************
!> \brief   Adapted from comab_opr. Calculate the product O*r or r*O from the integrals [a|O|b].
!>          We assume that on input all integrals [a+1|O|b+1] are available.
!> \param la_max ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param dab ...
!> \param ab ...
!> \param comabr ...
!>
!> \param ra ...
!> \param rb ...
!> \param direction_Or ...
!> \par Literature
!>          S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - ncoset    : Number of orbitals in a Cartesian orbital set.
!>      - npgf{a,b} : Degree of contraction of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rab2      : Square of the distance between the atomic centers a and b.
!>      - rac       : Distance vector between the atomic centers a and c.
!>      - rac2      : Square of the distance between the atomic centers a and c.
!>      - rbc       : Distance vector between the atomic centers b and c.
!>      - rbc2      : Square of the distance between the atomic centers b and c.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
!>
!> \author  Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE ab_opr(la_max, npgfa, rpgfa, la_min, lb_max, npgfb, rpgfb, lb_min, &
                     dab, ab, comabr, ra, rb, direction_Or)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: comabr
      REAL(KIND=dp), DIMENSION(1:3), INTENT(IN)          :: ra, rb
      LOGICAL                                            :: direction_Or

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, coap, &
                                                            coapx, coapy, coapz, cob, cobp, cobpx, &
                                                            cobpy, cobpz, ipgf, jpgf, la, lb, na, &
                                                            nap, nb, nbp, ofa, ofb

      comabr = 0.0_dp

      ofa = ncoset(la_min - 1)
      ofb = ncoset(lb_min - 1)

      na = 0
      nap = 0
      DO ipgf = 1, npgfa
         nb = 0
         nbp = 0
         DO jpgf = 1, npgfb
            IF (rpgfa(ipgf) + rpgfb(jpgf) > dab) THEN
               DO la = la_min, la_max
                  DO ax = 0, la
                     DO ay = 0, la - ax
                        az = la - ax - ay
                        coa = na + coset(ax, ay, az) - ofa
                        coap = nap + coset(ax, ay, az) - ofa
                        coapx = nap + coset(ax + 1, ay, az) - ofa
                        coapy = nap + coset(ax, ay + 1, az) - ofa
                        coapz = nap + coset(ax, ay, az + 1) - ofa
                        DO lb = lb_min, lb_max
                           DO bx = 0, lb
                              DO by = 0, lb - bx
                                 bz = lb - bx - by
                                 cob = nb + coset(bx, by, bz) - ofb
                                 cobp = nbp + coset(bx, by, bz) - ofb
                                 cobpx = nbp + coset(bx + 1, by, bz) - ofb
                                 cobpy = nbp + coset(bx, by + 1, bz) - ofb
                                 cobpz = nbp + coset(bx, by, bz + 1) - ofb
                                 IF (direction_Or) THEN
                                    ! [a|O * x|b] = [a|O|b(x+1)] + [a|O|b] * X_b
                                    !             = [a|O * (x - X_b)|b] + [a|O|b] * X_b
                                    ! So the second term makes sure that we actually calculate
                                    !  <O*r> and not <O*(r-R)>
                                    comabr(coa, cob, 1) = ab(coap, cobpx) + ab(coap, cobp)*rb(1)
                                    comabr(coa, cob, 2) = ab(coap, cobpy) + ab(coap, cobp)*rb(2)
                                    comabr(coa, cob, 3) = ab(coap, cobpz) + ab(coap, cobp)*rb(3)
                                 ELSE
                                    comabr(coa, cob, 1) = ab(coapx, cobp) + ab(coap, cobp)*ra(1)
                                    comabr(coa, cob, 2) = ab(coapy, cobp) + ab(coap, cobp)*ra(2)
                                    comabr(coa, cob, 3) = ab(coapz, cobp) + ab(coap, cobp)*ra(3)
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            nb = nb + ncoset(lb_max) - ofb
            nbp = nbp + ncoset(lb_max + 1) - ofb
         END DO
         na = na + ncoset(la_max) - ofa
         nap = nap + ncoset(la_max + 1) - ofa
      END DO

   END SUBROUTINE ab_opr

! **************************************************************************************************
!> \brief Apply the operator \delta_\mu^\lambda to zero out all elements of the matrix
!>         which don't fulfill the condition.
!> \param matrix ...
!> \param qs_kind_set ...
!> \param basis_type ...
!> \param sab_nl ...
!> \param lambda ...
!> \param direction_Or ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE hr_mult_by_delta_1d(matrix, qs_kind_set, basis_type, sab_nl, lambda, direction_Or)

      TYPE(dbcsr_type)                                   :: matrix
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, INTENT(IN)                                :: lambda
      LOGICAL                                            :: direction_Or

      CHARACTER(len=*), PARAMETER :: routineN = 'hr_mult_by_delta_1d'

      INTEGER                                            :: handle, iatom, icol, ikind, irow, jatom, &
                                                            jkind, ldsab, mepos, nkind, nseta, &
                                                            nsetb, nthread
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: k_block, rpgfa, rpgfb, scon_a, scon_b, &
                                                            sphi_a, sphi_b, zeta, zetb
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      nkind = SIZE(qs_kind_set)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      ! prepare basis set
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      ! *** Allocate work storage ***
      ldsab = get_memory_usage(qs_kind_set, basis_type)

      nthread = 1
!$    nthread = omp_get_max_threads()
      ! Iterate of neighbor list
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,ldsab,nl_iterator, do_symmetric) &
!$OMP SHARED (ncoset,matrix,basis_set_list) &
!$OMP SHARED (direction_or, lambda) &
!$OMP PRIVATE (k_block,mepos,ikind,jkind,iatom,jatom,rab,cell) &
!$OMP PRIVATE (basis_set_a,basis_set_b) &
!$OMP PRIVATE (first_sgfa, la_max, la_min, npgfa, nsgfa, nseta, rpgfa, set_radius_a) &
!$OMP PRIVATE (zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, rpgfb, set_radius_b, nsgfb) &
!$OMP PRIVATE (zetb, scon_a, scon_b, irow, icol, found, sphi_a, sphi_b)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cell)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         scon_a => basis_set_a%scon
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         scon_b => basis_set_b%scon

         nseta = basis_set_a%nset
         nsetb = basis_set_b%nset

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         NULLIFY (k_block)
         CALL dbcsr_get_block_p(matrix, irow, icol, k_block, found)
         CPASSERT(found)

         IF (direction_Or) THEN
            IF (jatom /= lambda) k_block(:, :) = 0._dp
         ELSE IF (.NOT. direction_Or) THEN
            IF (iatom /= lambda) k_block(:, :) = 0._dp
         END IF
      END DO
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE hr_mult_by_delta_1d

! **************************************************************************************************
!> \brief Apply the operator \delta_\mu^\lambda to zero out all elements of the matrix
!>         which don't fulfill the condition.
!>        Operates on matrix_hr(1:3) instead of a single matrix
!> \param matrix_hr ...
!> \param qs_kind_set ...
!> \param basis_type ...
!> \param sab_nl ...
!> \param deltaR ...
!> \param direction_Or ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE hr_mult_by_delta_3d(matrix_hr, qs_kind_set, basis_type, sab_nl, deltaR, direction_Or)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hr
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      REAL(KIND=dp), DIMENSION(:, :)                     :: deltaR
      LOGICAL                                            :: direction_Or

      CHARACTER(len=*), PARAMETER :: routineN = 'hr_mult_by_delta_3d'

      INTEGER                                            :: handle, iatom, icol, ikind, ir, irow, &
                                                            jatom, jkind, ldsab, mepos, nkind, &
                                                            nseta, nsetb, nthread
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kx_block, ky_block, kz_block, rpgfa, &
                                                            rpgfb, scon_a, scon_b, sphi_a, sphi_b, &
                                                            zeta, zetb
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      nkind = SIZE(qs_kind_set)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      ! prepare basis set
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      ! *** Allocate work storage ***
      ldsab = get_memory_usage(qs_kind_set, basis_type)

      nthread = 1
!$    nthread = omp_get_max_threads()
      ! Iterate of neighbor list
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,ldsab,nl_iterator, do_symmetric) &
!$OMP SHARED (ncoset,matrix_hr,basis_set_list) &
!$OMP SHARED (direction_or, deltar) &
!$OMP PRIVATE (kx_block,ky_block,kz_block,mepos,ikind,jkind,iatom,jatom,rab,cell) &
!$OMP PRIVATE (basis_set_a,basis_set_b) &
!$OMP PRIVATE (nseta, nsetb) &
!$OMP PRIVATE (first_sgfa, la_max, la_min, npgfa, nsgfa, rpgfa, set_radius_a, sphi_a, zeta, scon_a) &
!$OMP PRIVATE (first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, set_radius_b, sphi_b, zetb, scon_b) &
!$OMP PRIVATE (irow, icol, found)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cell)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         scon_a => basis_set_a%scon
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         scon_b => basis_set_b%scon

         nseta = basis_set_a%nset
         nsetb = basis_set_b%nset

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         NULLIFY (kx_block, ky_block, kz_block)
         CALL dbcsr_get_block_p(matrix_hr(1)%matrix, irow, icol, kx_block, found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix_hr(2)%matrix, irow, icol, ky_block, found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix_hr(3)%matrix, irow, icol, kz_block, found)
         CPASSERT(found)

         IF (direction_Or) THEN
            DO ir = 1, 3
!$OMP CRITICAL(blockadd)
               SELECT CASE (ir)
               CASE (1)
                  kx_block(:, :) = kx_block(:, :)*deltaR(ir, jatom)
               CASE (2)
                  ky_block(:, :) = ky_block(:, :)*deltaR(ir, jatom)
               CASE (3)
                  kz_block(:, :) = kz_block(:, :)*deltaR(ir, jatom)
               END SELECT
!$OMP END CRITICAL(blockadd)
            END DO
         ELSE
            DO ir = 1, 3
!$OMP CRITICAL(blockadd)
               SELECT CASE (ir)
               CASE (1)
                  kx_block(:, :) = kx_block(:, :)*deltaR(ir, iatom)
               CASE (2)
                  ky_block(:, :) = ky_block(:, :)*deltaR(ir, iatom)
               CASE (3)
                  kz_block(:, :) = kz_block(:, :)*deltaR(ir, iatom)
               END SELECT
!$OMP END CRITICAL(blockadd)
            END DO
         END IF
      END DO
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE hr_mult_by_delta_3d

END MODULE qs_vcd_ao

