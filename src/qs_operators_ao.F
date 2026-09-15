!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      created 07.2005
!> \author MI (07.2005)
! **************************************************************************************************
MODULE qs_operators_ao
   USE ai_angmom,                       ONLY: angmom
   USE ai_moments,                      ONLY: contract_cossin,&
                                              cossin
   USE ai_overlap,                      ONLY: overlap_ab
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_get_matrix_type,&
                                              dbcsr_has_symmetry,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_operators_ao'

! *** Public subroutines ***

   PUBLIC :: build_exp_ikr_matrix
   PUBLIC :: build_lin_mom_matrix, build_ang_mom_matrix

CONTAINS

! **************************************************************************************************
!> \brief Build real and imaginary AO matrices for exp(i*k*r).
!> \param qs_env ...
!> \param op_sm_set op_sm_set(1, i) is the cosine matrix and op_sm_set(2, i) is the sine matrix
!>        for kvec(:, i)
!> \param kvec Cartesian wave vectors, one column for each output matrix pair
!> \param sab_orb_external optional neighbor list used instead of the default orbital list
!> \param basis_type optional basis-set name
!> \param force_periodic temporarily use all three periodic directions
!> \param cell_external optional cell used for coordinate folding and periodicity
!> \author CP2K developers
!> \note The cosine and sine matrices must be preallocated with matching symmetric
!>       DBCSR topology. The matrices are overwritten from zero and only canonical
!>       atom-pair blocks are written.
! **************************************************************************************************
   SUBROUTINE build_exp_ikr_matrix(qs_env, op_sm_set, kvec, sab_orb_external, basis_type, force_periodic, cell_external)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: op_sm_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(LEN=*), OPTIONAL                         :: basis_type
      LOGICAL, OPTIONAL                                  :: force_periodic
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_external

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_exp_ikr_matrix'

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         ldsa, ldsb, ldwork, ncoa, ncob, nkind, nkvec, nseta, nsetb, reim, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: perd0
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found, my_force_periodic, new_atom_b
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cosab, rpgfa, rpgfb, sinab, sphi_a, &
                                                            sphi_b, work, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:, :)   :: op_cossin
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

      CPASSERT(SIZE(kvec, 1) == 3)
      nkvec = SIZE(kvec, 2)
      CPASSERT(nkvec > 0)
      CPASSERT(SIZE(op_sm_set, 1) == 2)
      CPASSERT(SIZE(op_sm_set, 2) == nkvec)
      DO i = 1, nkvec
         DO reim = 1, 2
            CPASSERT(ASSOCIATED(op_sm_set(reim, i)%matrix))
            CPASSERT(dbcsr_has_symmetry(op_sm_set(reim, i)%matrix))
            CALL dbcsr_set(op_sm_set(reim, i)%matrix, 0.0_dp)
         END DO
      END DO

      NULLIFY (qs_kind, qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, sab_orb=sab_orb)
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(particle_set))
      IF (PRESENT(cell_external)) THEN
         CPASSERT(ASSOCIATED(cell_external))
         cell => cell_external
      END IF

      IF (PRESENT(sab_orb_external)) THEN
         CPASSERT(ASSOCIATED(sab_orb_external))
         sab_orb => sab_orb_external
      END IF
      CPASSERT(ASSOCIATED(sab_orb))

      my_force_periodic = .FALSE.
      IF (PRESENT(force_periodic)) my_force_periodic = force_periodic
      IF (my_force_periodic) THEN
         perd0(:) = cell%perd(:)
         cell%perd(:) = 1
      END IF

      nkind = SIZE(qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ALLOCATE (cosab(ldwork, ldwork), sinab(ldwork, ldwork), work(ldwork, ldwork))

      ALLOCATE (op_cossin(2, nkvec))

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

      last_jatom = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE

         ra(:) = pbc(particle_set(iatom)%r(:), cell)
         rb(:) = ra(:) + rab(:)

         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet

         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         nseta = basis_set_a%nset
         nsetb = basis_set_b%nset
         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)
         IF (inode == 1) last_jatom = 0

         IF (jatom /= last_jatom) THEN
            new_atom_b = .TRUE.
            last_jatom = jatom
         ELSE
            new_atom_b = .FALSE.
         END IF

         IF (new_atom_b) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF

            DO i = 1, nkvec
               DO reim = 1, 2
                  NULLIFY (op_cossin(reim, i)%block)
                  CALL dbcsr_get_block_p(matrix=op_sm_set(reim, i)%matrix, &
                                         row=irow, col=icol, block=op_cossin(reim, i)%block, found=found)
               END DO
               IF (ASSOCIATED(op_cossin(1, i)%block) .NEQV. ASSOCIATED(op_cossin(2, i)%block)) THEN
                  CPABORT("cosine and sine blocks should have the same topology")
               END IF
            END DO
         END IF

         dab = NORM2(rab)
         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               DO i = 1, nkvec
                  IF (.NOT. ASSOCIATED(op_cossin(1, i)%block)) CYCLE
                  CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), &
                              la_min(iset), lb_max(jset), npgfb(jset), zetb(:, jset), &
                              rpgfb(:, jset), lb_min(jset), ra, rb, kvec(:, i), cosab, sinab)
                  CALL contract_cossin(op_cossin(1, i)%block, op_cossin(2, i)%block, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosab, sinab, ldwork, work, ldwork)
               END DO
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (op_cossin, cosab, sinab, work, basis_set_list)

      IF (my_force_periodic) cell%perd(:) = perd0(:)

      CALL timestop(handle)
   END SUBROUTINE build_exp_ikr_matrix

! **************************************************************************************************
!> \brief   Calculation of the linear momentum matrix <mu|∂|nu> over
!>          Cartesian Gaussian functions.
!> \param qs_env ...
!> \param matrix ...
!> \param minimum_image take into account only the first neighbors in the lists
!> \date    27.02.2009
!> \author  VW
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE build_lin_mom_matrix(qs_env, matrix, minimum_image)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: minimum_image

      CHARACTER(len=*), PARAMETER :: routineN = 'build_lin_mom_matrix'

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         ldai, maxco, maxlgto, maxsgf, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found, my_minimum_image, &
                                                            new_atom_b
      REAL(KIND=dp)                                      :: alpha, dab, Lxo2, Lyo2, Lzo2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: rr_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         TARGET                                          :: intab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: integral
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      NULLIFY (cell, sab_nl, qs_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell)

      nkind = SIZE(qs_kind_set)

      my_minimum_image = .FALSE.
      IF (PRESENT(minimum_image)) THEN
         my_minimum_image = minimum_image
         Lxo2 = NORM2(cell%hmat(:, 1))/2.0_dp
         Lyo2 = NORM2(cell%hmat(:, 2))/2.0_dp
         Lzo2 = NORM2(cell%hmat(:, 3))/2.0_dp
      END IF

      ! Take into account the symmetry of the input matrix
      do_symmetric = dbcsr_has_symmetry(matrix(1)%matrix)
      IF (do_symmetric) THEN
         CALL get_qs_env(qs_env=qs_env, sab_orb=sab_nl)
      ELSE
         CALL get_qs_env(qs_env=qs_env, sab_all=sab_nl)
      END IF
!   *** Allocate work storage ***

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxlgto=maxlgto, &
                           maxsgf=maxsgf)

      ldai = ncoset(maxlgto + 1)
      CALL init_orbital_pointers(ldai)

      ALLOCATE (rr_work(ldai*ldai*3), intab(maxco, maxco, 3), work(maxco, maxsgf), integral(3))
      rr_work(:) = 0.0_dp
      intab(:, :, :) = 0.0_dp
      work(:, :) = 0.0_dp

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
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
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

         IF (my_minimum_image) THEN
            IF (ANY(ABS(rab(:)) > [Lxo2, Lyo2, Lzo2])) CYCLE
         END IF

         IF (jatom /= last_jatom) THEN
            new_atom_b = .TRUE.
            last_jatom = jatom
         ELSE
            new_atom_b = .FALSE.
         END IF

         IF (new_atom_b) THEN
            alpha = 1.0_dp
            IF (do_symmetric) THEN
               IF (iatom <= jatom) THEN
                  irow = iatom
                  icol = jatom
               ELSE
                  irow = jatom
                  icol = iatom
                  IF (dbcsr_get_matrix_type(matrix(1)%matrix) == dbcsr_type_antisymmetric) THEN
                     alpha = -1.0_dp
                  END IF
               END IF
            ELSE
               irow = iatom
               icol = jatom
            END IF

            DO i = 1, 3
               NULLIFY (integral(i)%block)
               CALL dbcsr_get_block_p(matrix=matrix(i)%matrix, &
                                      row=irow, col=icol, BLOCK=integral(i)%block, found=found)
               CPASSERT(ASSOCIATED(INTEGRAL(i)%block))
            END DO
         END IF

         dab = NORM2(rab)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! *** Calculate the primitive fermi contact integrals ***

               ! Keep the full primitive layout used by the contraction below.
               CALL overlap_ab(la_max(iset), 0, npgfa(iset), &
                               rpgfa(:, iset), zeta(:, iset), &
                               lb_max(jset), 0, npgfb(jset), &
                               rpgfb(:, jset), zetb(:, jset), &
                               rab, dab=intab, rr_work=rr_work)

               ! *** Contraction step ***

               DO i = 1, 3

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, intab(1, 1, i), SIZE(intab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  IF (do_symmetric) THEN
                     IF (iatom <= jatom) THEN

                        CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   work(1, 1), SIZE(work, 1), &
                                   1.0_dp, integral(i)%block(sgfa, sgfb), &
                                   SIZE(integral(i)%block, 1))

                     ELSE

                        CALL dgemm("T", "N", nsgfb(jset), nsgfa(iset), ncoa, &
                                   alpha, work(1, 1), SIZE(work, 1), &
                                   sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   1.0_dp, integral(i)%block(sgfb, sgfa), &
                                   SIZE(integral(i)%block, 1))

                     END IF
                  ELSE
                     CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1.0_dp, integral(i)%block(sgfa, sgfb), &
                                SIZE(integral(i)%block, 1))
                  END IF

               END DO

            END DO

         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! *** Release work storage ***

      DEALLOCATE (intab, rr_work, work, integral, basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_lin_mom_matrix

! **************************************************************************************************
!> \brief   Calculation of the angular momentum matrix over
!>          Cartesian Gaussian functions.
!> \param qs_env ...
!> \param matrix ...
!> \param rc ...
!> \date    27.02.2009
!> \author  VW
!> \version 1.0
! **************************************************************************************************

   SUBROUTINE build_ang_mom_matrix(qs_env, matrix, rc)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: rc

      CHARACTER(len=*), PARAMETER :: routineN = 'build_ang_mom_matrix'

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         maxco, maxsgf, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found, new_atom_b
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: intab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rbc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: integral
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

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      sab_all=sab_all, &
                      cell=cell)

      nkind = SIZE(qs_kind_set)

!   *** Allocate work storage ***

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf=maxsgf)

      ALLOCATE (intab(maxco, maxco, 3), work(maxco, maxsgf), integral(3))
      intab(:, :, :) = 0.0_dp
      work(:, :) = 0.0_dp

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
      CALL neighbor_list_iterator_create(nl_iterator, sab_all)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ra = pbc(particle_set(iatom)%r, cell)
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

         IF (inode == 1) last_jatom = 0

         IF (jatom /= last_jatom) THEN
            new_atom_b = .TRUE.
            last_jatom = jatom
         ELSE
            new_atom_b = .FALSE.
         END IF

         IF (new_atom_b) THEN
            !IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            !ELSE
            !   irow = jatom
            !   icol = iatom
            !END IF

            DO i = 1, 3
               NULLIFY (INTEGRAL(i)%block)
               CALL dbcsr_get_block_p(matrix=matrix(i)%matrix, &
                                      row=irow, col=icol, BLOCK=INTEGRAL(i)%block, found=found)
               CPASSERT(found .AND. ASSOCIATED(INTEGRAL(i)%block))
            END DO
         END IF

         dab = NORM2(rab)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

               !IF(PRESENT(wancen)) THEN
               !   rc = wancen
               rac = pbc(rc, ra, cell)
               rbc = rac + rab
               !ELSE
               !   rc(1:3) = rb(1:3)
               !   rac(1:3) = -rab(1:3)
               !   rbc(1:3) = 0.0_dp
               !ENDIF

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! *** Calculate the primitive angular momentum integrals ***

               CALL angmom(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), rac, rbc, intab)

               ! *** Contraction step ***

               DO i = 1, 3

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, intab(1, 1, i), SIZE(intab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  !IF (iatom <= jatom) THEN

                  CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, integral(i)%block(sgfa, sgfb), &
                             SIZE(integral(i)%block, 1))

                  !ELSE
                  !
                  !   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                  !              -1.0_dp,work(1,1),SIZE(work,1),&
                  !              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  !              1.0_dp,integral(i)%block(sgfb,sgfa),&
                  !              SIZE(integral(i)%block,1))
                  !
                  !ENDIF

               END DO

            END DO

         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! *** Release work storage ***

      DEALLOCATE (intab, work, integral, basis_set_list)

!   *** Print the spin orbit matrix, if requested ***

      !IF (BTEST(cp_print_key_should_output(logger%iter_info,&
      !     qs_env%input,"DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM"),cp_p_file)) THEN
      !   iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM",&
      !        extension=".Log")
      !   CALL cp_dbcsr_write_sparse_matrix(matrix(1)%matrix,4,6,qs_env,para_env,output_unit=iw)
      !   CALL cp_dbcsr_write_sparse_matrix(matrix(2)%matrix,4,6,qs_env,para_env,output_unit=iw)
      !   CALL cp_dbcsr_write_sparse_matrix(matrix(3)%matrix,4,6,qs_env,para_env,output_unit=iw)
      !   CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
      !        "DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM")
      !END IF

      CALL timestop(handle)

   END SUBROUTINE build_ang_mom_matrix

! **************************************************************************************************
END MODULE qs_operators_ao
