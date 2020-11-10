!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that build the Kohn-Sham matrix for the LRIGPW
!>      and xc parts
!> \par History
!>      09.2013 created [Dorothea Golze]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE lri_ks_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_add_block_node,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE lri_compression,                 ONLY: lri_decomp_i
   USE lri_environment_types,           ONLY: lri_environment_type,&
                                              lri_int_type,&
                                              lri_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_o3c_methods,                  ONLY: contract3_o3c
   USE qs_o3c_types,                    ONLY: get_o3c_vec,&
                                              o3c_vec_create,&
                                              o3c_vec_release,&
                                              o3c_vec_type
   USE ri_environment_methods,          ONLY: ri_metric_solver

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_ks_methods'

   PUBLIC :: calculate_lri_ks_matrix, calculate_ri_ks_matrix

CONTAINS

!*****************************************************************************
!> \brief update of LRIGPW KS matrix
!> \param lri_env ...
!> \param lri_v_int integrals of potential * ri basis set
!> \param h_matrix KS matrix, on entry containing the core hamiltonian
!> \param atomic_kind_set ...
!> \param cell_to_index ...
!> \note including this in lri_environment_methods?
! **************************************************************************************************
   SUBROUTINE calculate_lri_ks_matrix(lri_env, lri_v_int, h_matrix, atomic_kind_set, cell_to_index)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: h_matrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER, DIMENSION(:, :, :), OPTIONAL, POINTER     :: cell_to_index

      CHARACTER(*), PARAMETER :: routineN = 'calculate_lri_ks_matrix'

      INTEGER :: atom_a, atom_b, col, handle, i, iac, iatom, ic, ikind, ilist, jatom, jkind, &
         jneighbor, mepos, natom, nba, nbb, nfa, nfb, nkind, nlist, nm, nn, nthread, row
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      LOGICAL                                            :: found, trans, use_cell_mapping
      REAL(KIND=dp)                                      :: dab, fw, isn, isna, isnb, rab(3), &
                                                            threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vi, via, vib
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hf_work, hs_work, int3, wab, wbb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block
      TYPE(dbcsr_type), POINTER                          :: hmat
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CALL timeset(routineN, handle)
      NULLIFY (atom_of_kind, h_block, lrii, nl_iterator, soo_list)

      threshold = lri_env%eps_o3_int

      use_cell_mapping = (SIZE(h_matrix, 1) > 1)
      IF (use_cell_mapping) THEN
         CPASSERT(PRESENT(cell_to_index))
      END IF

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nkind = lri_env%lri_ints%nkind
         nthread = 1
!$       nthread = omp_get_max_threads()

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)
         CALL neighbor_list_iterator_create(nl_iterator, soo_list, nthread=nthread)
!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (nthread,nl_iterator,nkind,atom_of_kind,threshold,lri_env,lri_v_int,&
!$OMP         h_matrix,use_cell_mapping,cell_to_index)&
!$OMP PRIVATE (mepos,ikind,jkind,iatom,jatom,nlist,ilist,jneighbor,rab,iac,dab,lrii,&
!$OMP          nfa,nfb,nba,nbb,nn,hs_work,hf_work,h_block,row,col,trans,found,wab,wbb,&
!$OMP          atom_a,atom_b,isn,nm,vi,isna,isnb,via,vib,fw,int3,cell,ic,hmat)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, inode=jneighbor, &
                                   r=rab, cell=cell)

            iac = ikind + nkind*(jkind - 1)
            dab = SQRT(SUM(rab*rab))

            IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE
            IF (lri_env%exact_1c_terms) THEN
               IF (iatom == jatom .AND. dab < lri_env%delta) CYCLE
            END IF

            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

            nfa = lrii%nfa
            nfb = lrii%nfb
            nba = lrii%nba
            nbb = lrii%nbb
            nn = nfa + nfb

            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)

            IF (use_cell_mapping) THEN
               ic = cell_to_index(cell(1), cell(2), cell(3))
               CPASSERT(ic > 0)
            ELSE
               ic = 1
            END IF
            hmat => h_matrix(ic)%matrix

            ALLOCATE (int3(nba, nbb))
            IF (lrii%lrisr) THEN
               ALLOCATE (hs_work(nba, nbb))
               IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                  nm = nfa
                  ALLOCATE (vi(nfa))
                  vi(1:nfa) = lri_v_int(ikind)%v_int(atom_a, 1:nfa)
               ELSE
                  nm = nn
                  ALLOCATE (vi(nn))
                  vi(1:nfa) = lri_v_int(ikind)%v_int(atom_a, 1:nfa)
                  vi(nfa + 1:nn) = lri_v_int(jkind)%v_int(atom_b, 1:nfb)
               END IF
               isn = SUM(lrii%sn(1:nm)*vi(1:nm))/lrii%nsn
               vi(1:nm) = MATMUL(lrii%sinv(1:nm, 1:nm), vi(1:nm)) - isn*lrii%sn(1:nm)
               hs_work(1:nba, 1:nbb) = isn*lrii%soo(1:nba, 1:nbb)
               IF (iatom == jatom .AND. dab < lri_env%delta) THEN
                  DO i = 1, nfa
                     CALL lri_decomp_i(int3, lrii%cabai, i)
                     hs_work(1:nba, 1:nbb) = hs_work(1:nba, 1:nbb) + vi(i)*int3(1:nba, 1:nbb)
                  END DO
               ELSE
                  DO i = 1, nfa
                     CALL lri_decomp_i(int3, lrii%cabai, i)
                     hs_work(1:nba, 1:nbb) = hs_work(1:nba, 1:nbb) + vi(i)*int3(1:nba, 1:nbb)
                  END DO
                  DO i = 1, nfb
                     CALL lri_decomp_i(int3, lrii%cabbi, i)
                     hs_work(1:nba, 1:nbb) = hs_work(1:nba, 1:nbb) + vi(nfa + i)*int3(1:nba, 1:nbb)
                  END DO
               END IF
               DEALLOCATE (vi)
            END IF

            IF (lrii%lriff) THEN
               ALLOCATE (hf_work(nba, nbb), wab(nba, nbb), wbb(nba, nbb))
               wab(1:nba, 1:nbb) = lri_env%wmat(ikind, jkind)%mat(1:nba, 1:nbb)
               wbb(1:nba, 1:nbb) = 1.0_dp - lri_env%wmat(ikind, jkind)%mat(1:nba, 1:nbb)
               !
               ALLOCATE (via(nfa), vib(nfb))
               via(1:nfa) = lri_v_int(ikind)%v_int(atom_a, 1:nfa)
               vib(1:nfb) = lri_v_int(jkind)%v_int(atom_b, 1:nfb)
               !
               isna = SUM(lrii%sna(1:nfa)*via(1:nfa))/lrii%nsna
               isnb = SUM(lrii%snb(1:nfb)*vib(1:nfb))/lrii%nsnb
               via(1:nfa) = MATMUL(lrii%asinv(1:nfa, 1:nfa), via(1:nfa)) - isna*lrii%sna(1:nfa)
               vib(1:nfb) = MATMUL(lrii%bsinv(1:nfb, 1:nfb), vib(1:nfb)) - isnb*lrii%snb(1:nfb)
               !
               hf_work(1:nba, 1:nbb) = (isna*wab(1:nba, 1:nbb) + isnb*wbb(1:nba, 1:nbb))*lrii%soo(1:nba, 1:nbb)
               !
               DO i = 1, nfa
                  IF (lrii%abascr(i) > threshold) THEN
                     CALL lri_decomp_i(int3, lrii%cabai, i)
                     hf_work(1:nba, 1:nbb) = hf_work(1:nba, 1:nbb) + &
                                             via(i)*int3(1:nba, 1:nbb)*wab(1:nba, 1:nbb)
                  END IF
               END DO
               DO i = 1, nfb
                  IF (lrii%abbscr(i) > threshold) THEN
                     CALL lri_decomp_i(int3, lrii%cabbi, i)
                     hf_work(1:nba, 1:nbb) = hf_work(1:nba, 1:nbb) + &
                                             vib(i)*int3(1:nba, 1:nbb)*wbb(1:nba, 1:nbb)
                  END IF
               END DO
               !
               DEALLOCATE (via, vib, wab, wbb)
            END IF
            DEALLOCATE (int3)

            ! add h_work to core hamiltonian
            IF (iatom <= jatom) THEN
               row = iatom
               col = jatom
               trans = .FALSE.
            ELSE
               row = jatom
               col = iatom
               trans = .TRUE.
            END IF
!$OMP CRITICAL(addhamiltonian)
            NULLIFY (h_block)
            CALL dbcsr_get_block_p(hmat, row, col, h_block, found)
            IF (.NOT. ASSOCIATED(h_block)) THEN
               CALL dbcsr_add_block_node(hmat, row, col, h_block)
            END IF
            IF (lrii%lrisr) THEN
               fw = lrii%wsr
               IF (trans) THEN
                  h_block(1:nbb, 1:nba) = h_block(1:nbb, 1:nba) + fw*TRANSPOSE(hs_work(1:nba, 1:nbb))
               ELSE
                  h_block(1:nba, 1:nbb) = h_block(1:nba, 1:nbb) + fw*hs_work(1:nba, 1:nbb)
               END IF
            END IF
            IF (lrii%lriff) THEN
               fw = lrii%wff
               IF (trans) THEN
                  h_block(1:nbb, 1:nba) = h_block(1:nbb, 1:nba) + fw*TRANSPOSE(hf_work(1:nba, 1:nbb))
               ELSE
                  h_block(1:nba, 1:nbb) = h_block(1:nba, 1:nbb) + fw*hf_work(1:nba, 1:nbb)
               END IF
            END IF
!$OMP END CRITICAL(addhamiltonian)

            IF (lrii%lrisr) DEALLOCATE (hs_work)
            IF (lrii%lriff) DEALLOCATE (hf_work)
         END DO
!$OMP END PARALLEL

         DO ic = 1, SIZE(h_matrix, 1)
            CALL dbcsr_finalize(h_matrix(ic)%matrix)
         END DO

         CALL neighbor_list_iterator_release(nl_iterator)
         DEALLOCATE (atom_of_kind)

      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_ks_matrix

!*****************************************************************************
!> \brief update of RIGPW KS matrix
!> \param lri_env ...
!> \param lri_v_int integrals of potential * ri basis set
!> \param h_matrix KS matrix, on entry containing the core hamiltonian
!> \param s_matrix overlap matrix
!> \param atomic_kind_set ...
!> \param ispin ...
!> \note including this in lri_environment_methods?
! **************************************************************************************************
   SUBROUTINE calculate_ri_ks_matrix(lri_env, lri_v_int, h_matrix, s_matrix, &
                                     atomic_kind_set, ispin)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(dbcsr_type), POINTER                          :: h_matrix, s_matrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(*), PARAMETER :: routineN = 'calculate_ri_ks_matrix'

      INTEGER                                            :: atom_a, handle, i1, i2, iatom, ikind, n, &
                                                            natom, nbas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, nsize
      INTEGER, DIMENSION(:, :), POINTER                  :: bas_ptr
      REAL(KIND=dp)                                      :: fscal, ftrm1n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fout, fvec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: o3c_vec

      CALL timeset(routineN, handle)

      bas_ptr => lri_env%ri_fit%bas_ptr
      natom = SIZE(bas_ptr, 2)
      nbas = bas_ptr(2, natom)
      ALLOCATE (fvec(nbas), fout(nbas))
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         atom_a = atom_of_kind(iatom)
         i1 = bas_ptr(1, iatom)
         i2 = bas_ptr(2, iatom)
         n = i2 - i1 + 1
         fvec(i1:i2) = lri_v_int(ikind)%v_int(atom_a, 1:n)
      END DO
      DEALLOCATE (atom_of_kind, kind_of)
      ! f(T) * R^(-1)*n
      ftrm1n = SUM(fvec(:)*lri_env%ri_fit%rm1n(:))
      lri_env%ri_fit%ftrm1n(ispin) = ftrm1n
      fscal = ftrm1n/lri_env%ri_fit%ntrm1n
      ! renormalize fvec -> fvec - fscal * n
      fvec(:) = fvec(:) - fscal*lri_env%ri_fit%nvec(:)
      ! solve Rx=f'
      CALL ri_metric_solver(mat=lri_env%ri_smat(1)%matrix, &
                            vecr=fvec(:), &
                            vecx=fout(:), &
                            matp=lri_env%ri_sinv(1)%matrix, &
                            solver=lri_env%ri_sinv_app, &
                            ptr=bas_ptr)
      lri_env%ri_fit%fout(:, ispin) = fout(:)

      ! add overlap matrix contribution
      CALL dbcsr_add(h_matrix, s_matrix, 1.0_dp, fscal)

      ! create a o3c_vec from fout
      ALLOCATE (nsize(natom), o3c_vec(natom))
      DO iatom = 1, natom
         i1 = bas_ptr(1, iatom)
         i2 = bas_ptr(2, iatom)
         n = i2 - i1 + 1
         nsize(iatom) = n
      END DO
      CALL o3c_vec_create(o3c_vec, nsize)
      DEALLOCATE (nsize)
      DO iatom = 1, natom
         i1 = bas_ptr(1, iatom)
         i2 = bas_ptr(2, iatom)
         n = i2 - i1 + 1
         CALL get_o3c_vec(o3c_vec, iatom, v)
         v(1:n) = fout(i1:i2)
      END DO
      ! add <T.f'>
      CALL contract3_o3c(lri_env%o3c, o3c_vec, h_matrix)
      !
      CALL o3c_vec_release(o3c_vec)
      DEALLOCATE (o3c_vec, fvec, fout)

      CALL timestop(handle)

   END SUBROUTINE calculate_ri_ks_matrix

END MODULE lri_ks_methods
