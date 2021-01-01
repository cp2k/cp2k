!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the moment integrals <a|r^m|b> and <a|r x d/dr|b>
!> \par History
!>      added angular moments (JGH 11.2012)
!> \author JGH (20.07.2006)
! **************************************************************************************************
MODULE qs_moments
   USE ai_angmom,                       ONLY: angmom
   USE ai_moments,                      ONLY: contract_cossin,&
                                              cossin,&
                                              diff_momop,&
                                              moment
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_double,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, &
        dbcsr_distribution_type, dbcsr_dot, dbcsr_get_block_p, dbcsr_init_p, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_sum
   USE moments_utils,                   ONLY: get_reference_point
   USE orbital_pointers,                ONLY: current_maxl,&
                                              indco,&
                                              ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: bohr,&
                                              debye
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_moments'

   ! Public subroutines
   PUBLIC :: build_berry_moment_matrix, build_local_moment_matrix
   PUBLIC :: build_berry_kpoint_matrix
   PUBLIC :: qs_moment_berry_phase, qs_moment_locop

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param moments ...
!> \param nmoments ...
!> \param ref_point ...
!> \param ref_points ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_local_moment_matrix(qs_env, moments, nmoments, ref_point, ref_points, basis_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      INTEGER, INTENT(IN)                                :: nmoments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ref_points
      CHARACTER(len=*), OPTIONAL                         :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_moment_matrix'

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         maxco, maxsgf, natom, ncoa, ncob, nkind, nm, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, npgfa, npgfb, &
                                                            nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: mab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: mint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (nmoments < 1) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, cell, particle_set, sab_orb)

      nm = (6 + 11*nmoments + 6*nmoments**2 + nmoments**3)/6 - 1
      CPASSERT(SIZE(moments) >= nm)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! Allocate work storage
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, maxsgf=maxsgf, &
                           basis_type=basis_type)

      ALLOCATE (mab(maxco, maxco, nm))
      mab(:, :, :) = 0.0_dp

      ALLOCATE (work(maxco, maxsgf))
      work(:, :) = 0.0_dp

      ALLOCATE (mint(nm))
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
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
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         IF (inode == 1) last_jatom = 0

         ! this guarantees minimum image convention
         ! anything else would not make sense
         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         DO i = 1, nm
            NULLIFY (mint(i)%block)
            CALL dbcsr_get_block_p(matrix=moments(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mint(i)%block, found=found)
            mint(i)%block = 0._dp
         END DO

         ! fold atomic position back into unit cell
         IF (PRESENT(ref_points)) THEN
            rc(:) = 0.5_dp*(ref_points(:, iatom) + ref_points(:, jatom))
         ELSE IF (PRESENT(ref_point)) THEN
            rc(:) = ref_point(:)
         ELSE
            rc(:) = 0._dp
         END IF
         ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
         ! by folding around the center, such screwing can be avoided for a proper choice of center.
         ra(:) = pbc(particle_set(iatom)%r(:) - rc, cell) + rc
         rb(:) = pbc(particle_set(jatom)%r(:) - rc, cell) + rc
         ! we dont use PBC at this point
         rab(:) = ra(:) - rb(:)
         rac(:) = ra(:) - rc(:)
         rbc(:) = rb(:) - rc(:)
         dab = NORM2(rab)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), &
                           rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), &
                           rpgfb(:, jset), nmoments, rac, rbc, mab)

               ! Contraction step
               DO i = 1, nm

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, mab(1, 1, i), SIZE(mab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  IF (iatom <= jatom) THEN

                     CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1.0_dp, mint(i)%block(sgfa, sgfb), &
                                SIZE(mint(i)%block, 1))

                  ELSE

                     CALL dgemm("T", "N", nsgfb(jset), nsgfa(iset), ncoa, &
                                1.0_dp, work(1, 1), SIZE(work, 1), &
                                sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                1.0_dp, mint(i)%block(sgfb, sgfa), &
                                SIZE(mint(i)%block, 1))

                  END IF

               END DO

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (mab, basis_set_list)
      DEALLOCATE (work)
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO
      DEALLOCATE (mint)

      CALL timestop(handle)

   END SUBROUTINE build_local_moment_matrix

! **************************************************************************************************
!> \brief Calculate derivatives of multipole moments assuming the moments are already calculated
!>        Note that the multipole moments are symmetric while their derivatives are anti-symmetric
!> \param qs_env ...
!> \param moments contains the multipole moments, already calculated
!> \param moments_der will contain the derivatives of the multipole moments
!> \param nmoments order of the multipole moments
!> \param nders ...
!> \param ref_point ...
!> \note
!>        Adapted from rRc_xyz_der_ao in qs_operators_ao
! **************************************************************************************************
   SUBROUTINE build_local_moments_der_matrix(qs_env, moments, moments_der, nmoments, nders, ref_point)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: moments_der
      INTEGER, INTENT(IN)                                :: nmoments, nders
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_moments_der_matrix'

      INTEGER :: handle, i, iatom, icol, idir, ii, ikind, inode, ipgf, irow, iset, j, jatom, &
         jkind, jpgf, jset, last_jatom, M_dim, maxco, maxsgf, na, nb, ncoa, ncob, nda, ndb, nkind, &
         nm, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: difmab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: mab, mab_tmp
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: mom_block
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:, :)   :: mom_block_der
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      IF (nmoments < 1) RETURN

      nm = (6 + 11*nmoments + 6*nmoments**2 + nmoments**3)/6 - 1
      CPASSERT(SIZE(moments) >= nm)
      CPASSERT(nmoments >= nders)
      NULLIFY (qs_kind_set, particle_set, sab_all, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      sab_all=sab_all)

      nkind = SIZE(qs_kind_set)

      ! Work storage
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, maxsgf=maxsgf)

      NULLIFY (mab)
      ALLOCATE (mab(maxco, maxco, nm))
      mab(:, :, :) = 0.0_dp

      M_dim = ncoset(nders) - 1
      ! CPASSERT(M_dim <= SIZE(moments, 1))
      ALLOCATE (difmab(maxco, maxco, M_dim, 3))
      difmab(:, :, :, :) = 0.0_dp

      ALLOCATE (work(maxco, maxsgf))
      work(:, :) = 0.0_dp

      ALLOCATE (mom_block(nm))
      ALLOCATE (mom_block_der(M_dim, 3))
      DO i = 1, nm
         NULLIFY (mom_block(i)%block)
      ENDDO
      DO i = 1, M_dim
         DO idir = 1, 3
            NULLIFY (mom_block_der(i, idir)%block)
         ENDDO
      END DO

      NULLIFY (basis_set_a, basis_set_b, basis_set_list)
      NULLIFY (qs_kind)
      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      ! Calculate derivatives looping over neighbour list
      NULLIFY (first_sgfa, la_max, la_min, npgfa, nsgfa, set_radius_a, sphi_a, zeta)
      NULLIFY (first_sgfb, lb_max, lb_min, npgfb, nsgfb, set_radius_b, sphi_b, zetb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_all)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
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
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         ! reference point
         IF (PRESENT(ref_point)) THEN
            rc(:) = ref_point(:)
         ELSE
            rc(:) = 0._dp
         ENDIF
         ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
         ! by folding around the center, such screwing can be avoided for a proper choice of center.
         ra(:) = pbc(particle_set(iatom)%r(:) - rc, cell) + rc
         rb(:) = pbc(particle_set(jatom)%r(:) - rc, cell) + rc
         ! we dont use PBC at this point
         rab(:) = ra(:) - rb(:)
         rac(:) = ra(:) - rc(:)
         rbc(:) = rb(:) - rc(:)
         dab = NORM2(rab)

         ! get blocks
         IF (inode == 1) last_jatom = 0

         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         irow = iatom
         icol = jatom

         DO i = 1, nm
            NULLIFY (mom_block(i)%block)
            ! get block from pre calculated overlap matrix
            CALL dbcsr_get_block_p(matrix=moments(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mom_block(i)%block, found=found)
            CPASSERT(found .AND. ASSOCIATED(mom_block(i)%block))
            mom_block(i)%block = 0._dp
         ENDDO
         DO i = 1, M_dim
            DO idir = 1, 3
               NULLIFY (mom_block_der(i, idir)%block)
               CALL dbcsr_get_block_p(matrix=moments_der(i, idir)%matrix, &
                                      row=irow, col=icol, &
                                      block=mom_block_der(i, idir)%block, &
                                      found=found)
               CPASSERT(found .AND. ASSOCIATED(mom_block_der(i, idir)%block))
               mom_block_der(i, idir)%block = 0._dp
            ENDDO
         END DO

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ALLOCATE (mab_tmp(npgfa(iset)*ncoset(la_max(iset) + 1), &
                                 npgfb(jset)*ncoset(lb_max(jset) + 1), ncoset(nmoments) - 1))

               ! Calculate the primitive integrals (need l+1 for derivatives)
               CALL moment(la_max(iset) + 1, npgfa(iset), zeta(:, iset), &
                           rpgfa(:, iset), la_min(iset), &
                           lb_max(jset) + 1, npgfb(jset), zetb(:, jset), &
                           rpgfb(:, jset), nmoments, rac, rbc, mab_tmp)

               CALL diff_momop(la_max(iset), npgfa(iset), zeta(:, iset), &
                               rpgfa(:, iset), la_min(iset), &
                               lb_max(jset), npgfb(jset), zetb(:, jset), &
                               rpgfb(:, jset), lb_min(jset), &
                               nders, rac, rbc, difmab, mab_ext=mab_tmp)
               ! copy subset of mab_tmp (l+1) to mab (l)
               mab = 0.0_dp
               DO ii = 1, nm
                  na = 0
                  nda = 0
                  DO ipgf = 1, npgfa(iset)
                     nb = 0
                     ndb = 0
                     DO jpgf = 1, npgfb(jset)
                        DO j = 1, ncoset(lb_max(jset))
                           DO i = 1, ncoset(la_max(iset))
                              mab(i + na, j + nb, ii) = mab_tmp(i + nda, j + ndb, ii)
                           END DO ! i
                        END DO ! j
                        nb = nb + ncoset(lb_max(jset))
                        ndb = ndb + ncoset(lb_max(jset) + 1)
                     END DO ! jpgf
                     na = na + ncoset(la_max(iset))
                     nda = nda + ncoset(la_max(iset) + 1)
                  END DO ! ipgf

               ENDDO
               ! Contraction step
               DO i = 1, nm

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, mab(1, 1, i), SIZE(mab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, mom_block(i)%block(sgfa, sgfb), &
                             SIZE(mom_block(i)%block, 1))
               ENDDO
               DO i = 1, M_dim
                  DO idir = 1, 3
                     CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                                1.0_dp, difmab(1, 1, i, idir), SIZE(difmab, 1), &
                                sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                                0._dp, work(1, 1), SIZE(work, 1))

                     CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1._dp, mom_block_der(i, idir)%block(sgfa, sgfb), &
                                SIZE(mom_block_der(i, idir)%block, 1))
                  ENDDO
               END DO
               DEALLOCATE (mab_tmp)
            END DO
         END DO
      ENDDO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! deallocations
      DEALLOCATE (basis_set_list)
      DEALLOCATE (mab)
      DEALLOCATE (difmab)
      DEALLOCATE (work)
      DO i = 1, nm
         NULLIFY (mom_block(i)%block)
         DO idir = 1, 3
            NULLIFY (mom_block_der(i, idir)%block)
         ENDDO
      END DO
      DEALLOCATE (mom_block)
      DEALLOCATE (mom_block_der)

      CALL timestop(handle)

   END SUBROUTINE build_local_moments_der_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magmom ...
!> \param nmoments ...
!> \param ref_point ...
!> \param ref_points ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_local_magmom_matrix(qs_env, magmom, nmoments, ref_point, ref_points, basis_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom
      INTEGER, INTENT(IN)                                :: nmoments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ref_points
      CHARACTER(len=*), OPTIONAL                         :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_magmom_matrix'

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         maxco, maxsgf, natom, ncoa, ncob, nkind, nm, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, npgfa, npgfb, &
                                                            nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: mab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: mint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (nmoments < 1) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, cell, particle_set, sab_all)
      NULLIFY (matrix_s)

      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)

      ! magnetic dipoles/angular moments only
      nm = 3
      CALL dbcsr_allocate_matrix_set(magmom, nm)

      DO i = 1, nm
         ALLOCATE (magmom(i)%matrix)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, magmom(i)%matrix)
         CALL dbcsr_set(magmom(i)%matrix, 0.0_dp)
      END DO

      NULLIFY (qs_kind_set, particle_set, sab_all, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_all=sab_all)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! Allocate work storage
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, maxsgf=maxsgf)

      ALLOCATE (mab(maxco, maxco, nm))
      mab(:, :, :) = 0.0_dp

      ALLOCATE (work(maxco, maxsgf))
      work(:, :) = 0.0_dp

      ALLOCATE (mint(nm))
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_all)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
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
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         IF (inode == 1) last_jatom = 0

         ! this guarantees minimum image convention
         ! anything else would not make sense
         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         irow = iatom
         icol = jatom

         DO i = 1, nm
            NULLIFY (mint(i)%block)
            CALL dbcsr_get_block_p(matrix=magmom(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mint(i)%block, found=found)
            mint(i)%block = 0._dp
            CPASSERT(ASSOCIATED(mint(i)%block))
         END DO

         ! fold atomic position back into unit cell
         IF (PRESENT(ref_points)) THEN
            rc(:) = 0.5_dp*(ref_points(:, iatom) + ref_points(:, jatom))
         ELSE IF (PRESENT(ref_point)) THEN
            rc(:) = ref_point(:)
         ELSE
            rc(:) = 0._dp
         END IF
         ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
         ! by folding around the center, such screwing can be avoided for a proper choice of center.
         ra(:) = pbc(particle_set(iatom)%r(:) - rc, cell) + rc
         rb(:) = pbc(particle_set(jatom)%r(:) - rc, cell) + rc
         ! we dont use PBC at this point
         rab(:) = ra(:) - rb(:)
         rac(:) = ra(:) - rc(:)
         rbc(:) = rb(:) - rc(:)
         dab = NORM2(rab)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL angmom(la_max(iset), npgfa(iset), zeta(:, iset), &
                           rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), &
                           rpgfb(:, jset), rac, rbc, mab)

               ! Contraction step
               DO i = 1, nm

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, mab(1, 1, i), SIZE(mab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, mint(i)%block(sgfa, sgfb), &
                             SIZE(mint(i)%block, 1))

               END DO

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (mab, basis_set_list)
      DEALLOCATE (work)
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO
      DEALLOCATE (mint)

      CALL timestop(handle)

   END SUBROUTINE build_local_magmom_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec, sab_orb_external, basis_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix'

      INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, ldsa, &
         ldsb, ldwork, natom, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, cosab, sblock, sinab, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      IF (PRESENT(sab_orb_external)) sab_orb => sab_orb_external

      CALL dbcsr_set(sinmat, 0.0_dp)
      CALL dbcsr_set(cosmat, 0.0_dp)

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork
      ALLOCATE (cosab(ldab, ldab))
      ALLOCATE (sinab(ldab, ldab))
      ALLOCATE (work(ldwork, ldwork))

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
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
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         IF (ASSOCIATED(cblock) .AND. .NOT. ASSOCIATED(sblock) .OR. &
             .NOT. ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN
            CPABORT("")
         ENDIF

         IF (ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN

            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra + rab
            dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

            DO iset = 1, nseta

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)

               DO jset = 1, nsetb

                  IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)

                  ! Calculate the primitive integrals
                  CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                              ra, rb, kvec, cosab, sinab)
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosab, sinab, ldab, work, ldwork)

               END DO
            END DO

         ENDIF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosab)
      DEALLOCATE (sinab)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_berry_moment_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_berry_kpoint_matrix(qs_env, cosmat, sinmat, kvec, basis_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      CHARACTER(len=*), OPTIONAL                         :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_kpoint_matrix'

      INTEGER :: handle, i, iatom, ic, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, &
         ldsa, ldsb, ldwork, natom, ncoa, ncob, nimg, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: icell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb, row_blk_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_cell_mapping
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, cosab, sblock, sinab, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control)
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         use_cell_mapping = .TRUE.
      ELSE
         use_cell_mapping = .FALSE.
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)
      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type=basis_type)
         IF (ASSOCIATED(basis_set)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      ALLOCATE (row_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, &
                            basis=basis_set_list)
      CALL get_ks_env(ks_env, dbcsr_dist=dbcsr_dist)
      ! (re)allocate matrix sets
      CALL dbcsr_allocate_matrix_set(sinmat, 1, nimg)
      CALL dbcsr_allocate_matrix_set(cosmat, 1, nimg)
      DO i = 1, nimg
         ! sin
         ALLOCATE (sinmat(1, i)%matrix)
         CALL dbcsr_create(matrix=sinmat(1, i)%matrix, &
                           name="SINMAT", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(sinmat(1, i)%matrix, sab_orb)
         CALL dbcsr_set(sinmat(1, i)%matrix, 0.0_dp)
         ! cos
         ALLOCATE (cosmat(1, i)%matrix)
         CALL dbcsr_create(matrix=cosmat(1, i)%matrix, &
                           name="COSMAT", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(cosmat(1, i)%matrix, sab_orb)
         CALL dbcsr_set(cosmat(1, i)%matrix, 0.0_dp)
      END DO

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork)
      ldab = ldwork
      ALLOCATE (cosab(ldab, ldab))
      ALLOCATE (sinab(ldab, ldab))
      ALLOCATE (work(ldwork, ldwork))

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab, cell=icell)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
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
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         IF (use_cell_mapping) THEN
            ic = cell_to_index(icell(1), icell(2), icell(3))
            CPASSERT(ic > 0)
         ELSE
            ic = 1
         END IF

         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         CPASSERT(found)

         ra(:) = pbc(particle_set(iatom)%r(:), cell)
         rb(:) = ra + rab
         dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                           ra, rb, kvec, cosab, sinab)
               CALL contract_cossin(cblock, sblock, &
                                    iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                    jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                    cosab, sinab, ldab, work, ldwork)

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosab)
      DEALLOCATE (sinab)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)
      DEALLOCATE (row_blk_sizes)

      CALL timestop(handle)

   END SUBROUTINE build_berry_kpoint_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
! **************************************************************************************************
   SUBROUTINE qs_moment_berry_phase(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase'

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      COMPLEX(dp)                                        :: xphase(3), zdet, zdeta, zi(3), &
                                                            zij(3, 3), zijk(3, 3, 3), &
                                                            zijkl(3, 3, 3, 3), zphase(3), zz
      INTEGER                                            :: handle, i, ia, idim, ikind, ispin, ix, &
                                                            iy, iz, j, k, l, nao, nm, nmo, nmom, &
                                                            nmotot, tmp_dim
      LOGICAL                                            :: floating, ghost, uniform
      REAL(dp)                                           :: charge, ci(3), cij(3, 3), dd, occ, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: kvec, qq, rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: eigrmat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, opvec
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp

      CPASSERT(ASSOCIATED(qs_env))

      IF (ASSOCIATED(qs_env%ls_scf_env)) THEN
         IF (unit_number > 0) WRITE (unit_number, *) "Periodic moment calculation not implemented in linear scaling code"
         RETURN
      ENDIF

      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 2)

      nm = (6 + 11*nmom + 6*nmom**2 + nmom**3)/6 - 1
      ! rmom(:,1)=electronic
      ! rmom(:,2)=nuclear
      ! rmom(:,1)=total
      ALLOCATE (rmom(nm + 1, 3))
      ALLOCATE (rlab(nm + 1))
      rmom = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

      NULLIFY (dft_control, rho, cell, particle_set, results, para_env, &
               local_particles, matrix_s, mos, rho_ao)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      cell=cell, &
                      results=results, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      local_particles=local_particles, &
                      matrix_s=matrix_s, &
                      mos=mos)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))
      nmotot = 0
      DO ispin = 1, dft_control%nspins
         NULLIFY (tmp_fm_struct, mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         nmotot = nmotot + nmo
         CALL cp_fm_create(opvec(ispin)%matrix, mo_coeff%matrix_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, &
                                  ncol_global=nmo, para_env=para_env, context=mo_coeff%matrix_struct%context)
         DO i = 1, SIZE(op_fm_set, 1)
            NULLIFY (op_fm_set(i, ispin)%matrix)
            CALL cp_fm_create(op_fm_set(i, ispin)%matrix, tmp_fm_struct)
         END DO
         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      ! occupation
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=occ, uniform_occupation=uniform)
         IF (.NOT. uniform) THEN
            CPWARN("Berry phase moments for non uniform MOs' occupation numbers not implemented")
         END IF
      END DO

      ! reference point
      CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      rcc = pbc(rcc, cell)

      ! label
      DO l = 1, nm
         ix = indco(1, l + 1)
         iy = indco(2, l + 1)
         iz = indco(3, l + 1)
         CALL set_label(rlab(l + 1), ix, iy, iz)
      END DO

      ! nuclear contribution
      DO ia = 1, SIZE(particle_set)
         atomic_kind => particle_set(ia)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
         IF (.NOT. ghost .AND. .NOT. floating) THEN
            rmom(1, 2) = rmom(1, 2) - charge
         ENDIF
      END DO
      ria = twopi*MATMUL(cell%h_inv, rcc)
      zphase = CMPLX(COS(ria), SIN(ria), dp)**rmom(1, 2)

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            zi(:) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria, cell)
                  DO i = 1, 3
                     kvec(:) = twopi*cell%h_inv(i, :)
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zi(i) = zi(i)*zdeta
                  END DO
               ENDIF
            END DO
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
            qq = AIMAG(LOG(zi))
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)
         CASE (2)
            ! Quadrupole
            CPABORT("Berry phase moments bigger than 1 not implemented")
            zij(:, :) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge)
               ria = particle_set(ia)%r
               ria = pbc(ria, cell)
               DO i = 1, 3
                  DO j = i, 3
                     kvec(:) = twopi*(cell%h_inv(i, :) + cell%h_inv(j, :))
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zij(i, j) = zij(i, j)*zdeta
                     zij(j, i) = zij(i, j)
                  END DO
               END DO
            END DO
            DO i = 1, 3
               DO j = 1, 3
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
            DO k = 4, 9
               ix = indco(1, k + 1)
               iy = indco(2, k + 1)
               iz = indco(3, k + 1)
               IF (ix == 0) THEN
                  rmom(k + 1, 2) = cij(iy, iz)
               ELSE IF (iy == 0) THEN
                  rmom(k + 1, 2) = cij(ix, iz)
               ELSE IF (iz == 0) THEN
                  rmom(k + 1, 2) = cij(ix, iy)
               END IF
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      ! electronic contribution

      ria = twopi*REAL(nmotot, dp)*occ*MATMUL(cell%h_inv, rcc)
      xphase = CMPLX(COS(ria), SIN(ria), dp)

      ! charge
      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_dot(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1) + trace
      END DO

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            DO i = 1, 3
               kvec(:) = twopi*cell%h_inv(i, :)
               CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
               IF (qs_env%run_rtp) THEN
                  CALL get_qs_env(qs_env, rtp=rtp)
                  CALL get_rtp(rtp, mos_new=mos_new)
                  CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               ELSE
                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               END IF
               zdet = CMPLX(1._dp, 0._dp, dp)
               DO ispin = 1, dft_control%nspins
                  CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                  DO idim = 1, tmp_dim
                     eigrmat(ispin)%matrix%local_data(:, idim) = &
                        CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                              -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                  END DO
                  CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                  zdet = zdet*zdeta
                  IF (dft_control%nspins == 1) THEN
                     zdet = zdet*zdeta
                  ENDIF
               END DO
               zi(i) = zdet
            END DO
            zi = zi*xphase
         CASE (2)
            ! Quadrupole
            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :) + cell%h_inv(j, :))
                  CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
                  IF (qs_env%run_rtp) THEN
                     CALL get_qs_env(qs_env, rtp=rtp)
                     CALL get_rtp(rtp, mos_new=mos_new)
                     CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
                  ELSE
                     CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  END IF
                  zdet = CMPLX(1._dp, 0._dp, dp)
                  DO ispin = 1, dft_control%nspins
                     CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                     DO idim = 1, tmp_dim
                        eigrmat(ispin)%matrix%local_data(:, idim) = &
                           CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                                 -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                     END DO
                     CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                     zdet = zdet*zdeta
                     IF (dft_control%nspins == 1) THEN
                        zdet = zdet*zdeta
                     ENDIF
                  END DO
                  zij(i, j) = zdet*xphase(i)*xphase(j)
                  zij(j, i) = zdet*xphase(i)*xphase(j)
               END DO
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole (apply periodic (2 Pi) boundary conditions)
            ci = AIMAG(LOG(zi))
            DO i = 1, 3
               IF (qq(i) + ci(i) > pi) ci(i) = ci(i) - twopi
               IF (qq(i) + ci(i) < -pi) ci(i) = ci(i) + twopi
            END DO
            rmom(2:4, 1) = MATMUL(cell%hmat, ci)/twopi
         CASE (2)
            ! Quadrupole
            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = 1, 3
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
            DO k = 4, 9
               ix = indco(1, k + 1)
               iy = indco(2, k + 1)
               iz = indco(3, k + 1)
               IF (ix == 0) THEN
                  rmom(k + 1, 1) = cij(iy, iz)
               ELSE IF (iy == 0) THEN
                  rmom(k + 1, 1) = cij(ix, iz)
               ELSE IF (iz == 0) THEN
                  rmom(k + 1, 1) = cij(ix, iy)
               END IF
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      rmom(:, 3) = rmom(:, 1) + rmom(:, 2)
      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))
      IF (magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE., mmom=mmom)
      ELSE
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE.)
      END IF

      DEALLOCATE (rmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (op_fm_set)
      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)

      CALL timestop(handle)

   END SUBROUTINE qs_moment_berry_phase

! **************************************************************************************************
!> \brief ...
!> \param cosmat ...
!> \param sinmat ...
!> \param mos ...
!> \param op_fm_set ...
!> \param opvec ...
! **************************************************************************************************
   SUBROUTINE op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)

      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: opvec

      INTEGER                                            :: i, nao, nmo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff

      DO i = 1, SIZE(op_fm_set, 2) ! spin
         CALL get_mo_set(mo_set=mos(i)%mo_set, nao=nao, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_dbcsr_sm_fm_multiply(cosmat, mo_coeff, opvec(i)%matrix, ncol=nmo)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, opvec(i)%matrix, 0.0_dp, &
                      op_fm_set(1, i)%matrix)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, mo_coeff, opvec(i)%matrix, ncol=nmo)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, opvec(i)%matrix, 0.0_dp, &
                      op_fm_set(2, i)%matrix)
      ENDDO

   END SUBROUTINE op_orbbas

! **************************************************************************************************
!> \brief ...
!> \param cosmat ...
!> \param sinmat ...
!> \param mos ...
!> \param op_fm_set ...
!> \param mos_new ...
! **************************************************************************************************
   SUBROUTINE op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)

      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new

      INTEGER                                            :: i, icol, lcol, nao, newdim, nmo
      LOGICAL                                            :: double_col, double_row
      TYPE(cp_fm_struct_type), POINTER                   :: newstruct, newstruct1
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, work, work1, work2

      DO i = 1, SIZE(op_fm_set, 2) ! spin
         CALL get_mo_set(mo_set=mos(i)%mo_set, nao=nao, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_fm_get_info(mos_new(2*i)%matrix, ncol_local=lcol, ncol_global=nmo)
         double_col = .TRUE.
         double_row = .FALSE.
         CALL cp_fm_struct_double(newstruct, &
                                  mos_new(2*i)%matrix%matrix_struct, &
                                  mos_new(2*i)%matrix%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work, matrix_struct=newstruct)
         CALL cp_fm_create(work1, matrix_struct=newstruct)
         CALL cp_fm_create(work2, matrix_struct=newstruct)
         CALL cp_fm_get_info(work, ncol_global=newdim)

         CALL cp_fm_set_all(work, 0.0_dp, 0.0_dp)
         DO icol = 1, lcol
            work%local_data(:, icol) = mos_new(2*i - 1)%matrix%local_data(:, icol)
            work%local_data(:, icol + lcol) = mos_new(2*i)%matrix%local_data(:, icol)
         END DO

         CALL cp_dbcsr_sm_fm_multiply(cosmat, work, work1, ncol=newdim)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, work, work2, ncol=newdim)

         DO icol = 1, lcol
            work%local_data(:, icol) = work1%local_data(:, icol) - work2%local_data(:, icol + lcol)
            work%local_data(:, icol + lcol) = work1%local_data(:, icol + lcol) + work2%local_data(:, icol)
         END DO

         CALL cp_fm_release(work1)
         CALL cp_fm_release(work2)

         CALL cp_fm_struct_double(newstruct1, &
                                  op_fm_set(1, i)%matrix%matrix_struct, &
                                  op_fm_set(1, i)%matrix%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work1, matrix_struct=newstruct1)

         CALL cp_gemm("T", "N", nmo, newdim, nao, 1.0_dp, mos_new(2*i - 1)%matrix, &
                      work, 0.0_dp, work1)

         DO icol = 1, lcol
            op_fm_set(1, i)%matrix%local_data(:, icol) = work1%local_data(:, icol)
            op_fm_set(2, i)%matrix%local_data(:, icol) = work1%local_data(:, icol + lcol)
         END DO

         CALL cp_gemm("T", "N", nmo, newdim, nao, 1.0_dp, mos_new(2*i)%matrix, &
                      work, 0.0_dp, work1)

         DO icol = 1, lcol
            op_fm_set(1, i)%matrix%local_data(:, icol) = &
               op_fm_set(1, i)%matrix%local_data(:, icol) + work1%local_data(:, icol + lcol)
            op_fm_set(2, i)%matrix%local_data(:, icol) = &
               op_fm_set(2, i)%matrix%local_data(:, icol) - work1%local_data(:, icol)
         END DO

         CALL cp_fm_release(work)
         CALL cp_fm_release(work1)
         CALL cp_fm_struct_release(newstruct)
         CALL cp_fm_struct_release(newstruct1)

      ENDDO

   END SUBROUTINE op_orbbas_rtp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
!> \param vel_reprs ...
!> \param com_nl ...
! **************************************************************************************************
   SUBROUTINE qs_moment_locop(qs_env, magnetic, nmoments, reference, ref_point, unit_number, vel_reprs, com_nl)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), INTENT(IN), POINTER        :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number
      LOGICAL, INTENT(IN), OPTIONAL                      :: vel_reprs, com_nl

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_moment_locop'

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: akind, handle, i, ia, iatom, idir, &
                                                            ikind, ispin, ix, iy, iz, l, nm, nmm, &
                                                            nmom, order
      LOGICAL                                            :: my_com_nl, my_velreprs
      REAL(dp)                                           :: charge, dd, strace, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom, nlcom_rrv, nlcom_rv, nlcom_rxrv, &
                                                            qupole_der, rmom_vel
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom, matrix_s, moments, momentum, &
                                                            rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: moments_der
      TYPE(dbcsr_type), POINTER                          :: tmp_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CPASSERT(ASSOCIATED(qs_env))

      CALL timeset(routineN, handle)

      my_velreprs = .FALSE.
      IF (PRESENT(vel_reprs)) my_velreprs = vel_reprs
      IF (PRESENT(com_nl)) my_com_nl = com_nl

      ! reference point
      CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)

      ! only allow for moments up to maxl set by basis
      nmom = MIN(nmoments, current_maxl)
      ! electronic contribution
      NULLIFY (dft_control, rho, cell, particle_set, qs_kind_set, results, para_env, matrix_s, rho_ao, sab_orb)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      cell=cell, &
                      results=results, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      matrix_s=matrix_s, &
                      sab_all=sab_orb)

      IF (my_com_nl) THEN
         IF ((nmom >= 1) .AND. my_velreprs) THEN
            ALLOCATE (nlcom_rv(3))
            nlcom_rv(:) = 0._dp
         ENDIF
         IF ((nmom >= 2) .AND. my_velreprs) THEN
            ALLOCATE (nlcom_rrv(6))
            nlcom_rrv(:) = 0._dp
         ENDIF
         IF (magnetic) THEN
            ALLOCATE (nlcom_rxrv(3))
            nlcom_rxrv = 0._dp
         ENDIF
         ! Calculate non local correction terms
         CALL calculate_commutator_nl_terms(qs_env, nlcom_rv, nlcom_rxrv, nlcom_rrv, rcc)
      ENDIF

      NULLIFY (moments)
      nm = (6 + 11*nmom + 6*nmom**2 + nmom**3)/6 - 1
      CALL dbcsr_allocate_matrix_set(moments, nm)
      DO i = 1, nm
         ALLOCATE (moments(i)%matrix)
         IF (my_velreprs .AND. (nmom >= 2)) THEN
            CALL dbcsr_create(moments(i)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(moments(i)%matrix, sab_orb)
         ELSE
            CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
         ENDIF
         CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
      END DO

      ! calculate derivatives if quadrupole in vel. reprs. is requested
      IF (my_velreprs .AND. (nmom >= 2)) THEN
         NULLIFY (moments_der)
         CALL dbcsr_allocate_matrix_set(moments_der, 3, 3)
         DO i = 1, 3 ! x, y, z
            DO idir = 1, 3 ! d/dx, d/dy, d/dz
               CALL dbcsr_init_p(moments_der(i, idir)%matrix)
               CALL dbcsr_create(moments_der(i, idir)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_alloc_block_from_nbl(moments_der(i, idir)%matrix, sab_orb)
               CALL dbcsr_set(moments_der(i, idir)%matrix, 0.0_dp)
            ENDDO
         END DO
         CALL build_local_moments_der_matrix(qs_env, moments, moments_der, 2, 1, rcc)
      ELSE
         CALL build_local_moment_matrix(qs_env, moments, nmom, ref_point=rcc)
      ENDIF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ALLOCATE (rmom(nm + 1, 3))
      ALLOCATE (rlab(nm + 1))
      rmom = 0.0_dp
      rlab = ""

      IF ((my_velreprs .AND. (nmoments >= 1)) .OR. magnetic) THEN
         ! Allocate matrix to store the matrix product to be traced (dbcsr_dot only works for products of
         ! symmetric matrices)
         NULLIFY (tmp_ao)
         ALLOCATE (tmp_ao)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, tmp_ao)
      ENDIF

      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_dot(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1) + trace
      END DO

      DO i = 1, SIZE(moments)
         strace = 0._dp
         DO ispin = 1, dft_control%nspins
            IF (my_velreprs .AND. nmoments >= 2) THEN
               CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, moments(i)%matrix, &
                                   0.0_dp, tmp_ao)
               CALL dbcsr_trace(tmp_ao, trace)
            ELSE
               CALL dbcsr_dot(rho_ao(ispin)%matrix, moments(i)%matrix, trace)
            ENDIF
            strace = strace + trace
         END DO
         rmom(i + 1, 1) = strace
      END DO

      CALL dbcsr_deallocate_matrix_set(moments)

      ! nuclear contribution
      CALL get_qs_env(qs_env=qs_env, &
                      local_particles=local_particles)
      DO ikind = 1, SIZE(local_particles%n_el)
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            ! fold atomic positions back into unit cell
            ria = pbc(particle_set(iatom)%r - rcc, cell) + rcc
            ria = ria - rcc
            atomic_kind => particle_set(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind, kind_number=akind)
            CALL get_qs_kind(qs_kind_set(akind), core_charge=charge)
            rmom(1, 2) = rmom(1, 2) - charge
            DO l = 1, nm
               ix = indco(1, l + 1)
               iy = indco(2, l + 1)
               iz = indco(3, l + 1)
               dd = 1._dp
               IF (ix > 0) dd = dd*ria(1)**ix
               IF (iy > 0) dd = dd*ria(2)**iy
               IF (iz > 0) dd = dd*ria(3)**iz
               rmom(l + 1, 2) = rmom(l + 1, 2) - charge*dd
               CALL set_label(rlab(l + 1), ix, iy, iz)
            END DO
         END DO
      END DO
      CALL mp_sum(rmom(:, 2), para_env%group)
      rmom(:, :) = -rmom(:, :)
      rmom(:, 3) = rmom(:, 1) + rmom(:, 2)

      ! magnetic moments
      IF (magnetic) THEN
         NULLIFY (magmom)
         CALL build_local_magmom_matrix(qs_env, magmom, nmom, ref_point=rcc)
         nmm = SIZE(magmom)
         ALLOCATE (mmom(nmm))

         IF (qs_env%run_rtp) THEN
            ! get imaginary part of the density in rho_ao (the real part is not needed since the trace of the product
            ! of a symmetric (REAL(rho_ao)) and an anti-symmetric (L_AO) matrix is zero)
            ! There may be other cases, where the imaginary part of the density is relevant
            NULLIFY (rho_ao)
            CALL qs_rho_get(rho, rho_ao_im=rho_ao)
         ENDIF
         ! if the density is purely real this is an expensive way to calculate zero
         DO i = 1, SIZE(magmom)
            strace = 0._dp
            DO ispin = 1, dft_control%nspins
               CALL dbcsr_set(tmp_ao, 0.0_dp)
               CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, magmom(i)%matrix, &
                                   0.0_dp, tmp_ao)
               CALL dbcsr_trace(tmp_ao, trace)
               strace = strace + trace
            END DO
            mmom(i) = strace
         END DO

         CALL dbcsr_deallocate_matrix_set(magmom)
      END IF

      ! velocity representations
      IF (my_velreprs) THEN
         ALLOCATE (rmom_vel(nm))
         rmom_vel = 0.0_dp

         DO order = 1, nmom
            SELECT CASE (order)

            CASE (1) ! expectation value of momentum
               NULLIFY (sab_orb)
               CALL get_qs_env(qs_env, sab_orb=sab_orb)
               NULLIFY (momentum)
               CALL dbcsr_allocate_matrix_set(momentum, 3)
               DO i = 1, 3
                  CALL dbcsr_init_p(momentum(i)%matrix)
                  CALL dbcsr_create(momentum(i)%matrix, template=matrix_s(1)%matrix, &
                                    matrix_type=dbcsr_type_antisymmetric)
                  CALL cp_dbcsr_alloc_block_from_nbl(momentum(i)%matrix, sab_orb)
                  CALL dbcsr_set(momentum(i)%matrix, 0.0_dp)
               ENDDO
               CALL build_lin_mom_matrix(qs_env, momentum)

               ! imaginary part of the density for RTP, real part gives 0 since momentum is antisymmetric
               IF (qs_env%run_rtp) THEN
                  NULLIFY (rho_ao)
                  CALL qs_rho_get(rho, rho_ao_im=rho_ao)
                  DO idir = 1, SIZE(momentum)
                     strace = 0._dp
                     DO ispin = 1, dft_control%nspins
                        CALL dbcsr_set(tmp_ao, 0.0_dp)
                        CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, momentum(idir)%matrix, &
                                            0.0_dp, tmp_ao)
                        CALL dbcsr_trace(tmp_ao, trace)
                        strace = strace + trace
                     END DO
                     rmom_vel(idir) = rmom_vel(idir) + strace
                  ENDDO
               ENDIF

               CALL dbcsr_deallocate_matrix_set(momentum)

            CASE (2) ! expectation value of quadrupole moment in vel. reprs.
               ALLOCATE (qupole_der(9)) ! will contain the expectation values of r_\alpha * d/d r_\beta
               qupole_der = 0._dp

               NULLIFY (rho_ao)
               CALL qs_rho_get(rho, rho_ao=rho_ao)

               ! Calculate expectation value over real part
               trace = 0._dp
               DO i = 1, 3
                  DO idir = 1, 3
                     strace = 0._dp
                     DO ispin = 1, dft_control%nspins
                        CALL dbcsr_set(tmp_ao, 0._dp)
                        CALL dbcsr_multiply("T", "N", 1._dp, rho_ao(ispin)%matrix, moments_der(i, idir)%matrix, 0._dp, tmp_ao)
                        CALL dbcsr_trace(tmp_ao, trace)
                        strace = strace + trace
                     ENDDO
                     qupole_der((i - 1)*3 + idir) = qupole_der((i - 1)*3 + idir) + strace
                  ENDDO
               ENDDO

               IF (qs_env%run_rtp) THEN
                  NULLIFY (rho_ao)
                  CALL qs_rho_get(rho, rho_ao_im=rho_ao)

                  ! Calculate expectation value over imaginary part
                  trace = 0._dp
                  DO i = 1, 3
                     DO idir = 1, 3
                        strace = 0._dp
                        DO ispin = 1, dft_control%nspins
                           CALL dbcsr_set(tmp_ao, 0._dp)
                           CALL dbcsr_multiply("T", "N", 1._dp, rho_ao(ispin)%matrix, moments_der(i, idir)%matrix, 0._dp, tmp_ao)
                           CALL dbcsr_trace(tmp_ao, trace)
                           strace = strace + trace
                        ENDDO
                        qupole_der((i - 1)*3 + idir) = qupole_der((i - 1)*3 + idir) + strace
                     ENDDO
                  ENDDO
               ENDIF

               ! calculate vel. reprs. of quadrupole moment from derivatives
               rmom_vel(4) = -2*qupole_der(1) - rmom(1, 1)
               rmom_vel(5) = -qupole_der(2) - qupole_der(4)
               rmom_vel(6) = -qupole_der(3) - qupole_der(7)
               rmom_vel(7) = -2*qupole_der(5) - rmom(1, 1)
               rmom_vel(8) = -qupole_der(6) - qupole_der(8)
               rmom_vel(9) = -2*qupole_der(9) - rmom(1, 1)

               DEALLOCATE (qupole_der)
            CASE DEFAULT
            END SELECT
         ENDDO
      ENDIF

      IF ((my_velreprs .AND. (nmoments >= 1)) .OR. magnetic) THEN
         CALL dbcsr_deallocate_matrix(tmp_ao)
      ENDIF
      IF (my_velreprs .AND. (nmoments >= 2)) THEN
         CALL dbcsr_deallocate_matrix_set(moments_der)
      ENDIF

      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))

      IF (magnetic .AND. my_velreprs) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., mmom=mmom, rmom_vel=rmom_vel)
      ELSE IF (magnetic .AND. .NOT. my_velreprs) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., mmom=mmom)
      ELSE IF (my_velreprs .AND. .NOT. magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., rmom_vel=rmom_vel)
      ELSE
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE.)
      END IF

      IF (my_com_nl) THEN
         IF (magnetic) THEN
            mmom(:) = mmom(:) + nlcom_rxrv(:)
         ENDIF
         IF (my_velreprs .AND. (nmom >= 1)) THEN
            rmom_vel(1:3) = rmom_vel(1:3) + nlcom_rv
         ENDIF
         IF (my_velreprs .AND. (nmom >= 2)) THEN
            rmom_vel(4:9) = rmom_vel(4:9) + nlcom_rrv
         ENDIF
         IF (magnetic .AND. .NOT. my_velreprs) THEN
            CALL print_moments_nl(unit_number, nmom, rlab, mmom=mmom)
         ELSE IF (my_velreprs .AND. .NOT. magnetic) THEN
            CALL print_moments_nl(unit_number, nmom, rlab, rmom_vel=rmom_vel)
         ELSE IF (my_velreprs .AND. magnetic) THEN
            CALL print_moments_nl(unit_number, nmom, rlab, mmom=mmom, rmom_vel=rmom_vel)
         ENDIF

      ENDIF

      IF (my_com_nl) THEN
         IF (nmom >= 1 .AND. my_velreprs) DEALLOCATE (nlcom_rv)
         IF (nmom >= 2 .AND. my_velreprs) DEALLOCATE (nlcom_rrv)
         IF (magnetic) DEALLOCATE (nlcom_rxrv)
      ENDIF

      DEALLOCATE (rmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF
      IF (my_velreprs) THEN
         DEALLOCATE (rmom_vel)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE qs_moment_locop

! **************************************************************************************************
!> \brief ...
!> \param label ...
!> \param ix ...
!> \param iy ...
!> \param iz ...
! **************************************************************************************************
   SUBROUTINE set_label(label, ix, iy, iz)
      CHARACTER(LEN=*), INTENT(OUT)                      :: label
      INTEGER, INTENT(IN)                                :: ix, iy, iz

      INTEGER                                            :: i

      label = ""
      DO i = 1, ix
         WRITE (label(i:), "(A1)") "X"
      END DO
      DO i = ix + 1, ix + iy
         WRITE (label(i:), "(A1)") "Y"
      END DO
      DO i = ix + iy + 1, ix + iy + iz
         WRITE (label(i:), "(A1)") "Z"
      END DO

   END SUBROUTINE set_label

! **************************************************************************************************
!> \brief ...
!> \param unit_number ...
!> \param nmom ...
!> \param rmom ...
!> \param rlab ...
!> \param rcc ...
!> \param cell ...
!> \param periodic ...
!> \param mmom ...
!> \param rmom_vel ...
! **************************************************************************************************
   SUBROUTINE print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic, mmom, rmom_vel)
      INTEGER, INTENT(IN)                                :: unit_number, nmom
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rmom
      CHARACTER(LEN=8), DIMENSION(:)                     :: rlab
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rcc
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: periodic
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: mmom, rmom_vel

      INTEGER                                            :: i, i0, i1, j, l
      REAL(dp)                                           :: dd

      IF (unit_number > 0) THEN
         DO l = 0, nmom
            SELECT CASE (l)
            CASE (0)
               WRITE (unit_number, "(T3,A,T33,3F16.8)") "Reference Point [Bohr]", rcc
               WRITE (unit_number, "(T3,A)") "Charges"
               WRITE (unit_number, "(T5,A,T18,F14.8,T36,A,T42,F14.8,T60,A,T67,F14.8)") &
                  "Electronic=", rmom(1, 1), "Core=", rmom(1, 2), "Total=", rmom(1, 3)
            CASE (1)
               IF (periodic) THEN
                  WRITE (unit_number, "(T3,A)") "Dipole vectors are based on the periodic (Berry phase) operator."
                  WRITE (unit_number, "(T3,A)") "They are defined modulo integer multiples of the cell matrix [Debye]."
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[X] [", cell%hmat(1, :)*debye, "] [i]"
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[Y]=[", cell%hmat(2, :)*debye, "]*[j]"
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[Z] [", cell%hmat(3, :)*debye, "] [k]"
               ELSE
                  WRITE (unit_number, "(T3,A)") "Dipoles are based on the traditional operator."
               ENDIF
               dd = SQRT(SUM(rmom(2:4, 3)**2))*debye
               WRITE (unit_number, "(T3,A)") "Dipole moment [Debye]"
               WRITE (unit_number, "(T5,3(A,A,F14.8,1X),T60,A,T67,F14.8)") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye, i=2, 4), "Total=", dd
            CASE (2)
               WRITE (unit_number, "(T3,A)") "Quadrupole moment [Debye*Angstrom]"
               WRITE (unit_number, "(T17,3(A,A,F14.8,9X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr, i=5, 7)
               WRITE (unit_number, "(T17,3(A,A,F14.8,9X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr, i=8, 10)
            CASE (3)
               WRITE (unit_number, "(T3,A)") "Octapole moment [Debye*Angstrom**2]"
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=11, 14)
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=15, 18)
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=19, 20)
            CASE (4)
               WRITE (unit_number, "(T3,A)") "Hexadecapole moment [Debye*Angstrom**3]"
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=21, 24)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=25, 28)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=29, 32)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=32, 35)
            CASE DEFAULT
               WRITE (unit_number, "(T3,A,A,I2)") "Higher moment [Debye*Angstrom**(L-1)]", &
                  "  L=", l
               i0 = (6 + 11*(l - 1) + 6*(l - 1)**2 + (l - 1)**3)/6
               i1 = (6 + 11*l + 6*l**2 + l**3)/6 - 1
               dd = debye/(bohr)**(l - 1)
               DO i = i0, i1, 3
                  WRITE (unit_number, "(T18,3(A,A,F14.8,4X))") &
                     (TRIM(rlab(j + 1)), "=", rmom(j + 1, 3)*dd, j=i, MIN(i1, i + 2))
               END DO
            END SELECT
         END DO
         IF (PRESENT(mmom)) THEN
            IF (nmom >= 1) THEN
               dd = SQRT(SUM(mmom(1:3)**2))
               WRITE (unit_number, "(T3,A)") "Orbital angular momentum [a. u.]"
               WRITE (unit_number, "(T5,3(A,A,F14.8,1X),T60,A,T67,F14.8)") &
                  (TRIM(rlab(i + 1)), "=", mmom(i), i=1, 3), "Total=", dd
            END IF
         END IF
         IF (PRESENT(rmom_vel)) THEN
            DO l = 1, nmom
               SELECT CASE (l)
               CASE (1)
                  dd = SQRT(SUM(rmom_vel(1:3)**2))
                  WRITE (unit_number, "(T3,A)") "Expectation value of momentum operator [a. u.]"
                  WRITE (unit_number, "(T5,3(A,A,D14.8,1X),T60,A,T67,F14.8)") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=1, 3), "Total=", dd
               CASE (2)
                  WRITE (unit_number, "(T3,A)") "Expectation value of quadrupole operator in vel. repr. [a. u.]"
                  WRITE (unit_number, "(T17,3(A,A,D14.8,9X))") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=4, 6)
                  WRITE (unit_number, "(T17,3(A,A,D14.8,9X))") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=7, 9)
               CASE DEFAULT
               END SELECT
            ENDDO
         ENDIF
      END IF

   END SUBROUTINE print_moments

! **************************************************************************************************
!> \brief ...
!> \param unit_number ...
!> \param nmom ...
!> \param rlab ...
!> \param mmom ...
!> \param rmom_vel ...
! **************************************************************************************************
   SUBROUTINE print_moments_nl(unit_number, nmom, rlab, mmom, rmom_vel)
      INTEGER, INTENT(IN)                                :: unit_number, nmom
      CHARACTER(LEN=8), DIMENSION(:)                     :: rlab
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: mmom, rmom_vel

      INTEGER                                            :: i, l
      REAL(dp)                                           :: dd

      IF (unit_number > 0) THEN
         IF (PRESENT(mmom)) THEN
            IF (nmom >= 1) THEN
               dd = SQRT(SUM(mmom(1:3)**2))
               WRITE (unit_number, "(T3,A)") "Orbital angular momentum incl. com_nl [a. u.]"
               WRITE (unit_number, "(T5,3(A,A,D14.8,1X),T60,A,T67,F14.8)") &
                  (TRIM(rlab(i + 1)), "=", mmom(i), i=1, 3), "Total=", dd
            END IF
         END IF
         IF (PRESENT(rmom_vel)) THEN
            DO l = 1, nmom
               SELECT CASE (l)
               CASE (1)
                  dd = SQRT(SUM(rmom_vel(1:3)**2))
                  WRITE (unit_number, "(T3,A)") "Expectation value of momentum operator incl. com_nl [a. u.]"
                  WRITE (unit_number, "(T5,3(A,A,D14.8,1X),T60,A,T67,F14.8)") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=1, 3), "Total=", dd
               CASE (2)
                  WRITE (unit_number, "(T3,A)") "Expectation value of quadrupole operator in vel. repr. incl. com_nl [a. u.]"
                  WRITE (unit_number, "(T17,3(A,A,D14.8,9X))") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=4, 6)
                  WRITE (unit_number, "(T17,3(A,A,D14.8,9X))") &
                     (TRIM(rlab(i + 1)), "=", rmom_vel(i), i=7, 9)
               CASE DEFAULT
               END SELECT
            ENDDO
         ENDIF
      END IF

   END SUBROUTINE print_moments_nl

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nlcom_rv ...
!> \param nlcom_rxrv ...
!> \param nlcom_rrv ...
!> \param ref_point ...
! **************************************************************************************************
   SUBROUTINE calculate_commutator_nl_terms(qs_env, nlcom_rv, nlcom_rxrv, nlcom_rrv, ref_point)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), ALLOCATABLE, DIMENSION(:), OPTIONAL      :: nlcom_rv, nlcom_rxrv, nlcom_rrv
      REAL(dp), DIMENSION(3)                             :: ref_point

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_commutator_nl_terms'

      INTEGER                                            :: handle, ind, ispin
      LOGICAL                                            :: calc_rrv, calc_rv, calc_rxrv
      REAL(dp)                                           :: eps_ppnl, strace, trace
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_rrv, matrix_rv, matrix_rxrv, &
                                                            matrix_s, rho_ao
      TYPE(dbcsr_type), POINTER                          :: tmp_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      calc_rv = .FALSE.
      calc_rxrv = .FALSE.
      calc_rrv = .FALSE.

      IF (ALLOCATED(nlcom_rv)) calc_rv = .TRUE.
      IF (ALLOCATED(nlcom_rrv)) calc_rrv = .TRUE.
      IF (ALLOCATED(nlcom_rxrv)) calc_rxrv = .TRUE.

      IF (.NOT. (calc_rv .OR. calc_rrv .OR. calc_rxrv)) RETURN

      ! calculate expectation values
      ! real part
      ! This evaluates to zero, because all commutator matrices are anti-symmetric while
      ! the real part of the density matrix rho_ao is symmetric. Tr[A*S] = 0
      IF (calc_rv) THEN
         nlcom_rv(:) = 0._dp
      ENDIF

      IF (calc_rrv) THEN
         nlcom_rrv(:) = 0._dp
      ENDIF

      IF (calc_rxrv) THEN
         nlcom_rxrv(:) = 0._dp
      ENDIF

      ! imagninary part of the density matrix
      IF (qs_env%run_rtp) THEN
         NULLIFY (cell, matrix_s, particle_set, qs_kind_set, rho, sab_orb, sap_ppnl)
         CALL get_qs_env(qs_env, &
                         cell=cell, &
                         dft_control=dft_control, &
                         matrix_s=matrix_s, &
                         particle_set=particle_set, &
                         qs_kind_set=qs_kind_set, &
                         rho=rho, &
                         sab_orb=sab_orb, &
                         sap_ppnl=sap_ppnl)

         eps_ppnl = dft_control%qs_control%eps_ppnl
         ! Calculate commutators, only needed of there is an imaginary part of the density matrix
         ! Allocate storage
         NULLIFY (matrix_rv, matrix_rrv, matrix_rxrv)
         IF (calc_rv) THEN
            CALL dbcsr_allocate_matrix_set(matrix_rv, 3)
            DO ind = 1, 3
               CALL dbcsr_init_p(matrix_rv(ind)%matrix)
               CALL dbcsr_create(matrix_rv(ind)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_rv(ind)%matrix, sab_orb)
               CALL dbcsr_set(matrix_rv(ind)%matrix, 0._dp)
            ENDDO
         ENDIF

         IF (calc_rrv) THEN
            CALL dbcsr_allocate_matrix_set(matrix_rrv, 6)
            DO ind = 1, 6
               CALL dbcsr_init_p(matrix_rrv(ind)%matrix)
               CALL dbcsr_create(matrix_rrv(ind)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_rrv(ind)%matrix, sab_orb)
               CALL dbcsr_set(matrix_rrv(ind)%matrix, 0._dp)
            ENDDO
         ENDIF

         IF (calc_rxrv) THEN
            CALL dbcsr_allocate_matrix_set(matrix_rxrv, 3)
            DO ind = 1, 3
               CALL dbcsr_init_p(matrix_rxrv(ind)%matrix)
               CALL dbcsr_create(matrix_rxrv(ind)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_rxrv(ind)%matrix, sab_orb)
               CALL dbcsr_set(matrix_rxrv(ind)%matrix, 0._dp)
            ENDDO
         ENDIF

         ! calculate commutators in AO basis
         IF (calc_rv .AND. calc_rrv .AND. .NOT. calc_rxrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rv=matrix_rv, &
                                  matrix_rrv=matrix_rrv, ref_point=ref_point, cell=cell)
         ELSE IF (calc_rv .AND. calc_rxrv .AND. .NOT. calc_rrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rv=matrix_rv, &
                                  matrix_rxrv=matrix_rxrv, ref_point=ref_point, cell=cell)
         ELSE IF (calc_rrv .AND. calc_rxrv .AND. .NOT. calc_rv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rxrv=matrix_rxrv, &
                                  matrix_rrv=matrix_rrv, ref_point=ref_point, cell=cell)
         ELSE IF (calc_rv .AND. calc_rxrv .AND. calc_rrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rv=matrix_rv, &
                                  matrix_rxrv=matrix_rxrv, matrix_rrv=matrix_rrv, ref_point=ref_point, cell=cell)
         ELSE IF (calc_rv .AND. .NOT. calc_rrv .AND. .NOT. calc_rxrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rv=matrix_rv, &
                                  ref_point=ref_point, cell=cell)
         ELSE IF (calc_rrv .AND. .NOT. calc_rv .AND. .NOT. calc_rxrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rrv=matrix_rrv, &
                                  ref_point=ref_point, cell=cell)
         ELSE IF (calc_rxrv .AND. .NOT. calc_rv .AND. .NOT. calc_rrv) THEN
            CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rxrv=matrix_rxrv, &
                                  ref_point=ref_point, cell=cell)
         ENDIF

         NULLIFY (rho_ao)
         CALL qs_rho_get(rho, rho_ao_im=rho_ao)

         NULLIFY (tmp_ao)
         ALLOCATE (tmp_ao)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, tmp_ao)

         IF (calc_rv) THEN
            trace = 0._dp
            DO ind = 1, SIZE(matrix_rv)
               strace = 0._dp
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_set(tmp_ao, 0.0_dp)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, matrix_rv(ind)%matrix, &
                                      0.0_dp, tmp_ao)
                  CALL dbcsr_trace(tmp_ao, trace)
                  strace = strace + trace
               END DO
               nlcom_rv(ind) = nlcom_rv(ind) + strace
            ENDDO
         ENDIF

         IF (calc_rrv) THEN
            trace = 0._dp
            DO ind = 1, SIZE(matrix_rrv)
               strace = 0._dp
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_set(tmp_ao, 0.0_dp)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, matrix_rrv(ind)%matrix, &
                                      0.0_dp, tmp_ao)
                  CALL dbcsr_trace(tmp_ao, trace)
                  strace = strace + trace
               END DO
               nlcom_rrv(ind) = nlcom_rrv(ind) + strace
            ENDDO
         ENDIF

         IF (calc_rxrv) THEN
            trace = 0._dp
            DO ind = 1, SIZE(matrix_rxrv)
               strace = 0._dp
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_set(tmp_ao, 0.0_dp)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, matrix_rxrv(ind)%matrix, &
                                      0.0_dp, tmp_ao)
                  CALL dbcsr_trace(tmp_ao, trace)
                  strace = strace + trace
               END DO
               nlcom_rxrv(ind) = nlcom_rxrv(ind) + strace
            ENDDO
         ENDIF
         CALL dbcsr_deallocate_matrix(tmp_ao)
         IF (calc_rv) CALL dbcsr_deallocate_matrix_set(matrix_rv)
         IF (calc_rrv) CALL dbcsr_deallocate_matrix_set(matrix_rrv)
         IF (calc_rxrv) CALL dbcsr_deallocate_matrix_set(matrix_rxrv)
      ENDIF ! qs_env%run_rtp

      CALL timestop(handle)
   END SUBROUTINE calculate_commutator_nl_terms

END MODULE qs_moments
