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
   USE dbcsr_api,                       ONLY: dbcsr_add_block_node,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_environment_type,&
                                              lri_int_type,&
                                              lri_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_ks_methods'

   PUBLIC :: calculate_lri_ks_matrix

CONTAINS

!*****************************************************************************
!> \brief update of LRIPGW KS matrix
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

END MODULE lri_ks_methods
