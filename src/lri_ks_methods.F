!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
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
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE lri_environment_methods,         ONLY: ri_metric_solver
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
!> \note including this in lri_environment_methods?
! **************************************************************************************************
   SUBROUTINE calculate_lri_ks_matrix(lri_env, lri_v_int, h_matrix, &
                                      atomic_kind_set)

      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(dbcsr_type), POINTER                          :: h_matrix
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CHARACTER(*), PARAMETER :: routineN = 'calculate_lri_ks_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, col, handle, i, iac, iatom, ikind, ilist, jatom, jkind, &
         jneighbor, natom, nba, nbb, nfa, nfb, nkind, nlist, nn, nneighbor, row
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      LOGICAL                                            :: found, trans
      REAL(KIND=dp)                                      :: dab, rab(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block
      TYPE(lri_int_type), POINTER                        :: lrii
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list

      CALL timeset(routineN, handle)
      NULLIFY (atom_of_kind, h_block, lrii, nl_iterator, soo_list)

      IF (ASSOCIATED(lri_env%soo_list)) THEN
         soo_list => lri_env%soo_list

         nkind = lri_env%lri_ints%nkind

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)
         CALL neighbor_list_iterator_create(nl_iterator, soo_list)

         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nneighbor, inode=jneighbor, &
                                   r=rab)

            iac = ikind+nkind*(jkind-1)
            dab = SQRT(SUM(rab*rab))

            IF (.NOT. ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

            nfa = lrii%nfa
            nfb = lrii%nfb
            nba = lrii%nba
            nbb = lrii%nbb
            nn = nfa+nfb

            ALLOCATE (h_work(nba, nbb))
            h_work(:, :) = 0.0_dp

            IF (iatom <= jatom) THEN
               row = iatom
               col = jatom
               trans = .FALSE.
            ELSE
               row = jatom
               col = iatom
               trans = .TRUE.
            END IF

            NULLIFY (h_block)
            CALL dbcsr_get_block_p(h_matrix, row, col, h_block, found)
            IF (.NOT. ASSOCIATED(h_block)) THEN
               CALL dbcsr_add_block_node(h_matrix, row, col, h_block)
            END IF

            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)

            ! sum: Hab(a)+Hab(b)
            ! Hab(a)
            DO i = 1, nfa
               h_work(1:nba, 1:nbb) = h_work(1:nba, 1:nbb)+ &
                                      lri_v_int(ikind)%v_int(atom_a, i)*lrii%dacoef(1:nba, 1:nbb, i)
            ENDDO
            ! Hab(b)
            IF (dab > lri_env%delta) THEN
               DO i = 1, nfb
                  h_work(1:nba, 1:nbb) = h_work(1:nba, 1:nbb)+ &
                                         lri_v_int(jkind)%v_int(atom_b, i)*lrii%dacoef(1:nba, 1:nbb, nfa+i)
               ENDDO
            ENDIF

            ! add h_work to core hamiltonian
            IF (trans) THEN
               h_work(1:nba, 1:nbb) = h_work(1:nba, 1:nbb)+TRANSPOSE(h_block(1:nbb, 1:nba))
               h_block(1:nbb, 1:nba) = TRANSPOSE(h_work(1:nba, 1:nbb))
            ELSE
               h_block(1:nba, 1:nbb) = h_block(1:nba, 1:nbb)+h_work(1:nba, 1:nbb)
            ENDIF

            DEALLOCATE (h_work)
         END DO

         CALL dbcsr_finalize(h_matrix)

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

      CHARACTER(*), PARAMETER :: routineN = 'calculate_ri_ks_matrix', &
         routineP = moduleN//':'//routineN

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
         n = i2-i1+1
         fvec(i1:i2) = lri_v_int(ikind)%v_int(atom_a, 1:n)
      END DO
      DEALLOCATE (atom_of_kind, kind_of)
      ! f(T) * R^(-1)*n
      ftrm1n = SUM(fvec(:)*lri_env%ri_fit%rm1n(:))
      lri_env%ri_fit%ftrm1n(ispin) = ftrm1n
      fscal = ftrm1n/lri_env%ri_fit%ntrm1n
      ! renormalize fvec -> fvec - fscal * n
      fvec(:) = fvec(:)-fscal*lri_env%ri_fit%nvec(:)
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
         n = i2-i1+1
         nsize(iatom) = n
      END DO
      CALL o3c_vec_create(o3c_vec, nsize)
      DEALLOCATE (nsize)
      DO iatom = 1, natom
         i1 = bas_ptr(1, iatom)
         i2 = bas_ptr(2, iatom)
         n = i2-i1+1
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
