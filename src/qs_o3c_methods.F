!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods used with 3-center overlap type integrals containers
!> \par History
!>      - none
!>      - 11.2018 fixed OMP race condition in contract3_o3c routine (A.Bussy)
! **************************************************************************************************
MODULE qs_o3c_methods
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE ai_overlap3,                     ONLY: overlap3
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: ncoset
   USE qs_o3c_types,                    ONLY: &
        get_o3c_container, get_o3c_iterator_info, get_o3c_vec, o3c_container_type, o3c_iterate, &
        o3c_iterator_create, o3c_iterator_release, o3c_iterator_type, o3c_vec_type, &
        set_o3c_container

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_o3c_methods'

   PUBLIC :: calculate_o3c_integrals, contract12_o3c, contract3_o3c

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param o3c ...
!> \param calculate_forces ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE calculate_o3c_integrals(o3c, calculate_forces, matrix_p)
      TYPE(o3c_container_type), POINTER                  :: o3c
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_p

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_o3c_integrals'

      INTEGER :: egfa, egfb, egfc, handle, i, iatom, icol, ikind, irow, iset, ispin, j, jatom, &
         jkind, jset, katom, kkind, kset, mepos, ncoa, ncob, ncoc, ni, nj, nk, nseta, nsetb, &
         nsetc, nspin, nthread, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lc_max, &
                                                            lc_min, npgfa, npgfb, npgfc, nsgfa, &
                                                            nsgfb, nsgfc
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfc
      LOGICAL                                            :: do_force, found, trans
      REAL(KIND=dp)                                      :: dij, dik, djk, fpre
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sabc, sabc_contr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: iabdc, iadbc, idabc, sabdc, sdabc
      REAL(KIND=dp), DIMENSION(3)                        :: rij, rik, rjk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fi, fj, fk, pblock, rpgfa, rpgfb, rpgfc, &
                                                            sphi_a, sphi_b, sphi_c, tvec, zeta, &
                                                            zetb, zetc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_c
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      do_force = .FALSE.
      IF (PRESENT(calculate_forces)) do_force = calculate_forces
      CALL get_o3c_container(o3c, nspin=nspin)

      ! basis sets
      CALL get_o3c_container(o3c, basis_set_list_a=basis_set_list_a, &
                             basis_set_list_b=basis_set_list_b, basis_set_list_c=basis_set_list_c)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,ncoset,nspin,basis_set_list_a,basis_set_list_b,&
!$OMP         basis_set_list_c,do_force,matrix_p)&
!$OMP PRIVATE (mepos,ikind,jkind,kkind,basis_set_a,basis_set_b,basis_set_c,rij,rik,rjk,&
!$OMP          first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,rpgfa,set_radius_a,sphi_a,zeta,&
!$OMP          first_sgfb,lb_max,lb_min,npgfb,nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,&
!$OMP          first_sgfc,lc_max,lc_min,npgfc,nsetc,nsgfc,rpgfc,set_radius_c,sphi_c,zetc,&
!$OMP          iset,jset,kset,dij,dik,djk,ni,nj,nk,iabc,idabc,iadbc,iabdc,tvec,fi,fj,fk,ncoa,&
!$OMP          ncob,ncoc,sabc,sabc_contr,sdabc,sabdc,sgfa,sgfb,sgfc,egfa,egfb,egfc,i,j,&
!$OMP          pblock,pmat,ispin,iatom,jatom,katom,irow,icol,found,trans,fpre)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, &
                                    ikind=ikind, jkind=jkind, kkind=kkind, rij=rij, rik=rik, &
                                    integral=iabc, tvec=tvec, force_i=fi, force_j=fj, force_k=fk)
         CPASSERT(.NOT. ASSOCIATED(iabc))
         CPASSERT(.NOT. ASSOCIATED(tvec))
         CPASSERT(.NOT. ASSOCIATED(fi))
         CPASSERT(.NOT. ASSOCIATED(fj))
         CPASSERT(.NOT. ASSOCIATED(fk))
         ! basis
         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         basis_set_c => basis_set_list_c(kkind)%gto_basis_set
         ! center A
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
         ! center B
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
         ! center C (RI)
         first_sgfc => basis_set_c%first_sgf
         lc_max => basis_set_c%lmax
         lc_min => basis_set_c%lmin
         npgfc => basis_set_c%npgf
         nsetc = basis_set_c%nset
         nsgfc => basis_set_c%nsgf_set
         rpgfc => basis_set_c%pgf_radius
         set_radius_c => basis_set_c%set_radius
         sphi_c => basis_set_c%sphi
         zetc => basis_set_c%zet

         ni = SUM(nsgfa)
         nj = SUM(nsgfb)
         nk = SUM(nsgfc)

         ALLOCATE (iabc(ni, nj, nk))
         iabc(:, :, :) = 0.0_dp
         IF (do_force) THEN
            ALLOCATE (fi(nk, 3), fj(nk, 3), fk(nk, 3))
            fi(:, :) = 0.0_dp
            fj(:, :) = 0.0_dp
            fk(:, :) = 0.0_dp
            ALLOCATE (idabc(ni, nj, nk, 3))
            idabc(:, :, :, :) = 0.0_dp
            ALLOCATE (iadbc(ni, nj, nk, 3))
            iadbc(:, :, :, :) = 0.0_dp
            ALLOCATE (iabdc(ni, nj, nk, 3))
            iabdc(:, :, :, :) = 0.0_dp
         ELSE
            NULLIFY (fi, fj, fk)
         END IF
         ALLOCATE (tvec(nk, nspin))
         tvec(:, :) = 0.0_dp

         rjk(1:3) = rik(1:3) - rij(1:3)
         dij = NORM2(rij)
         dik = NORM2(rik)
         djk = NORM2(rjk)

         DO iset = 1, nseta
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dij) CYCLE
               DO kset = 1, nsetc
                  IF (set_radius_a(iset) + set_radius_c(kset) < dik) CYCLE
                  IF (set_radius_b(jset) + set_radius_c(kset) < djk) CYCLE

                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))

                  sgfa = first_sgfa(1, iset)
                  sgfb = first_sgfb(1, jset)
                  sgfc = first_sgfc(1, kset)

                  egfa = sgfa + nsgfa(iset) - 1
                  egfb = sgfb + nsgfb(jset) - 1
                  egfc = sgfc + nsgfc(kset) - 1

                  IF (ncoa*ncob*ncoc > 0) THEN
                     ALLOCATE (sabc(ncoa, ncob, ncoc))
                     sabc(:, :, :) = 0.0_dp
                     IF (do_force) THEN
                        ALLOCATE (sdabc(ncoa, ncob, ncoc, 3))
                        sdabc(:, :, :, :) = 0.0_dp
                        ALLOCATE (sabdc(ncoa, ncob, ncoc, 3))
                        sabdc(:, :, :, :) = 0.0_dp
                        CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                      lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                      lc_max(kset), npgfc(kset), zetc(:, kset), rpgfc(:, kset), lc_min(kset), &
                                      rij, dij, rik, dik, rjk, djk, sabc, sdabc, sabdc)
                     ELSE
                        CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                      lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                      lc_max(kset), npgfc(kset), zetc(:, kset), rpgfc(:, kset), lc_min(kset), &
                                      rij, dij, rik, dik, rjk, djk, sabc)
                     END IF
                     ALLOCATE (sabc_contr(nsgfa(iset), nsgfb(jset), nsgfc(kset)))

                     CALL abc_contract(sabc_contr, sabc, &
                                       sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                       ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                     iabc(sgfa:egfa, sgfb:egfb, sgfc:egfc) = &
                        sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                     IF (do_force) THEN
                        DO i = 1, 3
                           CALL abc_contract(sabc_contr, sdabc(:, :, :, i), &
                                             sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                             ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                           idabc(sgfa:egfa, sgfb:egfb, sgfc:egfc, i) = &
                              sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                           CALL abc_contract(sabc_contr, sabdc(:, :, :, i), &
                                             sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                                             ncoa, ncob, ncoc, nsgfa(iset), nsgfb(jset), nsgfc(kset))
                           iabdc(sgfa:egfa, sgfb:egfb, sgfc:egfc, i) = &
                              sabc_contr(1:nsgfa(iset), 1:nsgfb(jset), 1:nsgfc(kset))
                        END DO
                     END IF

                     DEALLOCATE (sabc_contr)
                     DEALLOCATE (sabc)
                  END IF
                  IF (do_force) THEN
                     DEALLOCATE (sdabc, sabdc)
                  END IF
               END DO
            END DO
         END DO
         IF (do_force) THEN
            ! translational invariance
            iadbc(:, :, :, :) = -idabc(:, :, :, :) - iabdc(:, :, :, :)
            !
            ! get the atom indices
            CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, &
                                       iatom=iatom, jatom=jatom, katom=katom)
            !
            ! contract over i and j to get forces
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
               trans = .FALSE.
            ELSE
               irow = jatom
               icol = iatom
               trans = .TRUE.
            END IF
            IF (iatom == jatom) THEN
               fpre = 1.0_dp
            ELSE
               fpre = 2.0_dp
            END IF
            ALLOCATE (pmat(ni, nj))
            pmat(:, :) = 0.0_dp
            DO ispin = 1, nspin
               CALL dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix, &
                                      row=irow, col=icol, BLOCK=pblock, found=found)
               IF (found) THEN
                  IF (trans) THEN
                     pmat(:, :) = pmat(:, :) + TRANSPOSE(pblock(:, :))
                  ELSE
                     pmat(:, :) = pmat(:, :) + pblock(:, :)
                  END IF
               END IF
            END DO
            DO i = 1, 3
               DO j = 1, nk
                  fi(j, i) = fpre*SUM(pmat(:, :)*idabc(:, :, j, i))
                  fj(j, i) = fpre*SUM(pmat(:, :)*iadbc(:, :, j, i))
                  fk(j, i) = fpre*SUM(pmat(:, :)*iabdc(:, :, j, i))
               END DO
            END DO
            DEALLOCATE (pmat)
            !
            DEALLOCATE (idabc, iadbc, iabdc)
         END IF
         !
         CALL set_o3c_container(o3c_iterator, mepos=mepos, &
                                integral=iabc, tvec=tvec, force_i=fi, force_j=fj, force_k=fk)

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE calculate_o3c_integrals

! **************************************************************************************************
!> \brief Contraction of 3-tensor over indices 1 and 2 (assuming symmetry)
!>        t(k) = sum_ij (ijk)*p(ij)
!> \param o3c ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE contract12_o3c(o3c, matrix_p)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract12_o3c'

      INTEGER                                            :: handle, iatom, icol, ik, irow, ispin, &
                                                            jatom, mepos, nk, nspin, nthread
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(KIND=dp)                                      :: fpre
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, tvec
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      nspin = SIZE(matrix_p, 1)
      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,matrix_p,nspin)&
!$OMP PRIVATE (mepos,ispin,iatom,jatom,ik,nk,irow,icol,iabc,tvec,found,pblock,trans,fpre)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, &
                                    integral=iabc, tvec=tvec)
         nk = SIZE(tvec, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF
         IF (iatom == jatom) THEN
            fpre = 1.0_dp
         ELSE
            fpre = 2.0_dp
         END IF

         DO ispin = 1, nspin
            CALL dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix, &
                                   row=irow, col=icol, BLOCK=pblock, found=found)
            IF (found) THEN
               IF (trans) THEN
                  DO ik = 1, nk
                     tvec(ik, ispin) = fpre*SUM(TRANSPOSE(pblock(:, :))*iabc(:, :, ik))
                  END DO
               ELSE
                  DO ik = 1, nk
                     tvec(ik, ispin) = fpre*SUM(pblock(:, :)*iabc(:, :, ik))
                  END DO
               END IF
            END IF
         END DO

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract12_o3c

! **************************************************************************************************
!> \brief Contraction of 3-tensor over index 3
!>        h(ij) = h(ij) + sum_k (ijk)*v(k)
!> \param o3c ...
!> \param vec ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE contract3_o3c(o3c, vec, matrix)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: vec
      TYPE(dbcsr_type)                                   :: matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract3_o3c'

      INTEGER                                            :: handle, iatom, icol, ik, irow, jatom, &
                                                            katom, mepos, nk, nthread, s1, s2
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,vec,matrix)&
!$OMP PRIVATE (mepos,iabc,iatom,jatom,katom,irow,icol,trans,pblock,v,found,ik,nk,work,s1,s2)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, katom=katom, &
                                    integral=iabc)

         CALL get_o3c_vec(vec, katom, v)
         nk = SIZE(v)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF

         CALL dbcsr_get_block_p(matrix=matrix, row=irow, col=icol, BLOCK=pblock, found=found)

         IF (found) THEN
            s1 = SIZE(pblock, 1); s2 = SIZE(pblock, 2)
            ALLOCATE (work(s1, s2))
            work(:, :) = 0.0_dp

            IF (trans) THEN
               DO ik = 1, nk
                  CALL daxpy(s1*s2, v(ik), TRANSPOSE(iabc(:, :, ik)), 1, work(:, :), 1)
               END DO
            ELSE
               DO ik = 1, nk
                  CALL daxpy(s1*s2, v(ik), iabc(:, :, ik), 1, work(:, :), 1)
               END DO
            END IF

            ! Multiple threads with same irow, icol but different katom (same even in PBCs) can try
            ! to access the dbcsr block at the same time. Prevent that by CRITICAL section but keep
            ! computations before hand in order to retain speed

!$OMP CRITICAL
            CALL dbcsr_get_block_p(matrix=matrix, row=irow, col=icol, BLOCK=pblock, found=found)
            CALL daxpy(s1*s2, 1.0_dp, work(:, :), 1, pblock(:, :), 1)
!$OMP END CRITICAL

            DEALLOCATE (work)
         END IF

      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract3_o3c

END MODULE qs_o3c_methods
