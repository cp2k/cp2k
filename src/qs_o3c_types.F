!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief 3-center overlap type integrals containers
!> \par History
!>      - Added options to only keep (abc) triplet if b and c share the same center (2019 A.Bussy)
! **************************************************************************************************
MODULE qs_o3c_types

   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE kinds,                           ONLY: dp
   USE qs_neighbor_list_types,          ONLY: &
        get_iterator_info, get_neighbor_list_set_p, neighbor_list_iterate, &
        neighbor_list_iterator_create, neighbor_list_iterator_p_type, &
        neighbor_list_iterator_release, neighbor_list_set_p_type, nl_set_sub_iterator, &
        nl_sub_iterate
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_o3c_types'

! **************************************************************************************************
! O3C Integrals
! **************************************************************************************************

   TYPE o3c_int_type
      PRIVATE
      INTEGER                                    :: katom, kkind
      INTEGER                                    :: ni, nj, nk
      REAL(KIND=dp), DIMENSION(3)                :: rik
      INTEGER, DIMENSION(3)                      :: cellk
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: integral
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: tvec
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: force_i
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: force_j
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: force_k
   END TYPE o3c_int_type

   TYPE o3c_pair_type
      PRIVATE
      INTEGER                                    :: iatom, ikind
      INTEGER                                    :: jatom, jkind
      REAL(KIND=dp), DIMENSION(3)                :: rij
      INTEGER, DIMENSION(3)                      :: cellj
      INTEGER                                    :: nklist
      TYPE(o3c_int_type), DIMENSION(:), POINTER  :: ijk
   END TYPE o3c_pair_type

   TYPE o3c_container_type
      PRIVATE
      LOGICAL                                    :: ijsymmetric
      INTEGER                                    :: nijpairs
      INTEGER                                    :: nspin
      TYPE(o3c_pair_type), DIMENSION(:), POINTER :: ijpair
      ! basis sets and neighbor lists are pointing to other resources
      ! we don't keep track if the data is available and correct
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         POINTER                                 :: basis_set_list_a, basis_set_list_b, &
                                                    basis_set_list_c
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                   :: sab_nl, sac_nl
   END TYPE o3c_container_type

! **************************************************************************************************
! O3C Iterator
! **************************************************************************************************

   TYPE o3c_iterator_type
      PRIVATE
      TYPE(o3c_container_type), POINTER     :: o3c
      INTEGER                               :: ijp_last, k_last
      INTEGER, DIMENSION(:), POINTER        :: ijp_thread, k_thread
   END TYPE o3c_iterator_type

! **************************************************************************************************
! O3C vector
! **************************************************************************************************

   TYPE o3c_vec_type
      PRIVATE
      INTEGER                               :: n
      REAL(KIND=dp), DIMENSION(:), POINTER  :: v
   END TYPE o3c_vec_type

! **************************************************************************************************

   PUBLIC :: o3c_container_type
   PUBLIC :: release_o3c_container, init_o3c_container, get_o3c_container, set_o3c_container
   PUBLIC :: o3c_iterator_type
   PUBLIC :: o3c_iterator_create, o3c_iterator_release, get_o3c_iterator_info, o3c_iterate
   PUBLIC :: o3c_vec_type
   PUBLIC :: o3c_vec_create, o3c_vec_release, get_o3c_vec

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param o3c ...
!> \param nspin ...
!> \param basis_set_list_a ...
!> \param basis_set_list_b ...
!> \param basis_set_list_c ...
!> \param sab_nl ...
!> \param sac_nl ...
!> \param only_bc_same_center only consider a,b,c atoms if b and c share the same center
!> \par History: only_bc_same_cetner added by A.Bussy for XAS_TDP (04.2019)
! **************************************************************************************************
   SUBROUTINE init_o3c_container(o3c, nspin, basis_set_list_a, basis_set_list_b, basis_set_list_c, &
                                 sab_nl, sac_nl, only_bc_same_center)
      TYPE(o3c_container_type)                           :: o3c
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl, sac_nl
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_bc_same_center

      INTEGER                                            :: kkind, nij, nk, nkind
      LOGICAL                                            :: my_sort_bc, symmetric
      REAL(dp)                                           :: rik(3), rjk(3)
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: ac_iterator, nl_iterator
      TYPE(o3c_int_type), POINTER                        :: ijk
      TYPE(o3c_pair_type), POINTER                       :: ijpair

      CALL get_neighbor_list_set_p(sab_nl, symmetric=symmetric)
      o3c%ijsymmetric = symmetric
      CPASSERT(symmetric)

      o3c%nspin = nspin

      o3c%basis_set_list_a => basis_set_list_a
      o3c%basis_set_list_b => basis_set_list_b
      o3c%basis_set_list_c => basis_set_list_c

      o3c%sab_nl => sab_nl
      o3c%sac_nl => sac_nl

      nkind = SIZE(basis_set_list_a)

      my_sort_bc = .FALSE.
      IF (PRESENT(only_bc_same_center)) my_sort_bc = only_bc_same_center

      ! determine the number of ij pairs
      nij = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         nij = nij + 1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      o3c%nijpairs = nij
      NULLIFY (o3c%ijpair)
      ALLOCATE (o3c%ijpair(nij))

      ! for each pair set up the ijk lists
      nij = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      CALL neighbor_list_iterator_create(ac_iterator, sac_nl, search=.TRUE.)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         nij = nij + 1
         ijpair => o3c%ijpair(nij)
         CALL get_iterator_info(nl_iterator, ikind=ijpair%ikind, jkind=ijpair%jkind, &
                                iatom=ijpair%iatom, jatom=ijpair%jatom, &
                                r=ijpair%rij, cell=ijpair%cellj)
         NULLIFY (ijpair%ijk)
         nk = 0
         DO kkind = 1, nkind
            CALL nl_set_sub_iterator(ac_iterator, ijpair%ikind, kkind, ijpair%iatom)
            DO WHILE (nl_sub_iterate(ac_iterator) == 0)
               IF (my_sort_bc) THEN
                  !we only take ijk if rjk = 0 OR rik = 0 (because of symmetry)
                  CALL get_iterator_info(ac_iterator, r=rik)
                  rjk(:) = rik(:) - ijpair%rij(:)
                  IF (.NOT. (ALL(ABS(rjk) .LE. 1.0E-4_dp) .OR. ALL(ABS(rik) .LE. 1.0E-4_dp))) CYCLE
               END IF
               nk = nk + 1
            END DO
         END DO
         ! ijk lists
         ijpair%nklist = nk
         ALLOCATE (ijpair%ijk(nk))
         ! fill the ijk lists
         nk = 0
         DO kkind = 1, nkind
            CALL nl_set_sub_iterator(ac_iterator, ijpair%ikind, kkind, ijpair%iatom)
            DO WHILE (nl_sub_iterate(ac_iterator) == 0)
               IF (my_sort_bc) THEN
                  !we only take ijk if rjk = 0 OR rik = 0 (because of symmetry)
                  CALL get_iterator_info(ac_iterator, r=rik)
                  rjk(:) = rik(:) - ijpair%rij(:)
                  IF (.NOT. (ALL(ABS(rjk) .LE. 1.0E-4_dp) .OR. ALL(ABS(rik) .LE. 1.0E-4_dp))) CYCLE
               END IF

               nk = nk + 1
               ijk => ijpair%ijk(nk)
               CALL get_iterator_info(ac_iterator, jatom=ijk%katom, r=ijk%rik, cell=ijk%cellk)
               ijk%kkind = kkind
               ijk%ni = 0
               ijk%nj = 0
               ijk%nk = 0
               NULLIFY (ijk%integral)
               NULLIFY (ijk%tvec)
               NULLIFY (ijk%force_i)
               NULLIFY (ijk%force_j)
               NULLIFY (ijk%force_k)
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(ac_iterator)
      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE init_o3c_container
! **************************************************************************************************
!> \brief ...
!> \param o3c_container ...
! **************************************************************************************************
   SUBROUTINE release_o3c_container(o3c_container)

      TYPE(o3c_container_type)                           :: o3c_container

      o3c_container%ijsymmetric = .FALSE.
      o3c_container%nijpairs = 0

      NULLIFY (o3c_container%basis_set_list_a)
      NULLIFY (o3c_container%basis_set_list_b)
      NULLIFY (o3c_container%basis_set_list_c)

      NULLIFY (o3c_container%sab_nl)
      NULLIFY (o3c_container%sac_nl)

      IF (ASSOCIATED(o3c_container%ijpair)) THEN
         CALL release_ijpair(o3c_container%ijpair)
         DEALLOCATE (o3c_container%ijpair)
      END IF

   END SUBROUTINE release_o3c_container

! **************************************************************************************************
!> \brief ...
!> \param ijpair ...
! **************************************************************************************************
   SUBROUTINE release_ijpair(ijpair)

      TYPE(o3c_pair_type), DIMENSION(:)                  :: ijpair

      INTEGER                                            :: i

      DO i = 1, SIZE(ijpair)
         ijpair(i)%iatom = 0
         ijpair(i)%ikind = 0
         ijpair(i)%jatom = 0
         ijpair(i)%jkind = 0
         ijpair(i)%nklist = 0
         ijpair(i)%rij = 0.0_dp
         ijpair(i)%cellj = 0
         IF (ASSOCIATED(ijpair(i)%ijk)) THEN
            CALL release_ijk(ijpair(i)%ijk)
            DEALLOCATE (ijpair(i)%ijk)
         END IF
      END DO

   END SUBROUTINE release_ijpair

! **************************************************************************************************
!> \brief ...
!> \param ijk ...
! **************************************************************************************************
   SUBROUTINE release_ijk(ijk)

      TYPE(o3c_int_type), DIMENSION(:)                   :: ijk

      INTEGER                                            :: i

      DO i = 1, SIZE(ijk)
         ijk(i)%katom = 0
         ijk(i)%kkind = 0
         ijk(i)%ni = 0
         ijk(i)%nj = 0
         ijk(i)%nk = 0
         ijk(i)%rik = 0.0_dp
         ijk(i)%cellk = 0
         IF (ASSOCIATED(ijk(i)%integral)) THEN
            DEALLOCATE (ijk(i)%integral)
         END IF
         IF (ASSOCIATED(ijk(i)%tvec)) THEN
            DEALLOCATE (ijk(i)%tvec)
         END IF
         IF (ASSOCIATED(ijk(i)%force_i)) THEN
            DEALLOCATE (ijk(i)%force_i)
         END IF
         IF (ASSOCIATED(ijk(i)%force_j)) THEN
            DEALLOCATE (ijk(i)%force_j)
         END IF
         IF (ASSOCIATED(ijk(i)%force_k)) THEN
            DEALLOCATE (ijk(i)%force_k)
         END IF
      END DO

   END SUBROUTINE release_ijk

! **************************************************************************************************
!> \brief ...
!> \param o3c ...
!> \param ijsymmetric ...
!> \param nspin ...
!> \param nijpairs ...
!> \param ijpair ...
!> \param basis_set_list_a ...
!> \param basis_set_list_b ...
!> \param basis_set_list_c ...
!> \param sab_nl ...
!> \param sac_nl ...
! **************************************************************************************************

   SUBROUTINE get_o3c_container(o3c, ijsymmetric, nspin, nijpairs, ijpair, &
                                basis_set_list_a, basis_set_list_b, basis_set_list_c, &
                                sab_nl, sac_nl)
      TYPE(o3c_container_type)                           :: o3c
      LOGICAL, OPTIONAL                                  :: ijsymmetric
      INTEGER, OPTIONAL                                  :: nspin, nijpairs
      TYPE(o3c_pair_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ijpair
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: basis_set_list_a, basis_set_list_b, &
                                                            basis_set_list_c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_nl, sac_nl

      IF (PRESENT(ijsymmetric)) ijsymmetric = o3c%ijsymmetric
      IF (PRESENT(nspin)) nspin = o3c%nspin
      IF (PRESENT(nijpairs)) nijpairs = o3c%nijpairs
      IF (PRESENT(ijpair)) ijpair => o3c%ijpair
      IF (PRESENT(basis_set_list_a)) basis_set_list_a => o3c%basis_set_list_a
      IF (PRESENT(basis_set_list_b)) basis_set_list_b => o3c%basis_set_list_b
      IF (PRESENT(basis_set_list_c)) basis_set_list_c => o3c%basis_set_list_c
      IF (PRESENT(sab_nl)) sab_nl => o3c%sab_nl
      IF (PRESENT(sac_nl)) sac_nl => o3c%sac_nl

   END SUBROUTINE get_o3c_container

! **************************************************************************************************
! O3C Iterator
! **************************************************************************************************
!> \brief ...
!> \param o3c ...
!> \param o3c_iterator ...
!> \param nthread ...
! **************************************************************************************************
   SUBROUTINE o3c_iterator_create(o3c, o3c_iterator, nthread)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      INTEGER, OPTIONAL                                  :: nthread

      INTEGER                                            :: n

      IF (PRESENT(nthread)) THEN
         n = nthread
      ELSE
         n = 1
      END IF

      o3c_iterator%o3c => o3c
      o3c_iterator%ijp_last = 0
      o3c_iterator%k_last = 0
      ALLOCATE (o3c_iterator%ijp_thread(0:n - 1))
      ALLOCATE (o3c_iterator%k_thread(0:n - 1))
      o3c_iterator%ijp_thread = 0
      o3c_iterator%k_thread = 0

   END SUBROUTINE o3c_iterator_create

! **************************************************************************************************
!> \brief ...
!> \param o3c_iterator ...
! **************************************************************************************************
   SUBROUTINE o3c_iterator_release(o3c_iterator)
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      NULLIFY (o3c_iterator%o3c)
      o3c_iterator%ijp_last = 0
      o3c_iterator%k_last = 0
      DEALLOCATE (o3c_iterator%ijp_thread)
      DEALLOCATE (o3c_iterator%k_thread)

   END SUBROUTINE o3c_iterator_release

! **************************************************************************************************
!> \brief ...
!> \param o3c_iterator ...
!> \param mepos ...
!> \param iatom ...
!> \param jatom ...
!> \param katom ...
!> \param ikind ...
!> \param jkind ...
!> \param kkind ...
!> \param rij ...
!> \param rik ...
!> \param cellj ...
!> \param cellk ...
!> \param integral ...
!> \param tvec ...
!> \param force_i ...
!> \param force_j ...
!> \param force_k ...
! **************************************************************************************************
   SUBROUTINE get_o3c_iterator_info(o3c_iterator, mepos, &
                                    iatom, jatom, katom, ikind, jkind, kkind, &
                                    rij, rik, cellj, cellk, &
                                    integral, tvec, force_i, force_j, force_k)
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      INTEGER, OPTIONAL                                  :: mepos, iatom, jatom, katom, ikind, &
                                                            jkind, kkind
      REAL(KIND=dp), DIMENSION(3), OPTIONAL              :: rij, rik
      INTEGER, DIMENSION(3), OPTIONAL                    :: cellj, cellk
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: integral
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: tvec, force_i, force_j, force_k

      INTEGER                                            :: ij, k, me
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_int_type), POINTER                        :: ijk
      TYPE(o3c_pair_type), POINTER                       :: ijp

      IF (PRESENT(mepos)) THEN
         me = mepos
      ELSE
         me = 0
      END IF

      ij = o3c_iterator%ijp_thread(me)
      k = o3c_iterator%k_thread(me)

      o3c => o3c_iterator%o3c
      ijp => o3c%ijpair(ij)
      ijk => ijp%ijk(k)

      IF (PRESENT(iatom)) iatom = ijp%iatom
      IF (PRESENT(jatom)) jatom = ijp%jatom
      IF (PRESENT(ikind)) ikind = ijp%ikind
      IF (PRESENT(jkind)) jkind = ijp%jkind
      IF (PRESENT(katom)) katom = ijk%katom
      IF (PRESENT(kkind)) kkind = ijk%kkind

      IF (PRESENT(rij)) rij(1:3) = ijp%rij(1:3)
      IF (PRESENT(rik)) rik(1:3) = ijk%rik(1:3)

      IF (PRESENT(cellj)) cellj(1:3) = ijp%cellj(1:3)
      IF (PRESENT(cellk)) cellk(1:3) = ijk%cellk(1:3)

      IF (PRESENT(integral)) integral => ijk%integral
      IF (PRESENT(tvec)) tvec => ijk%tvec
      IF (PRESENT(force_i)) force_i => ijk%force_i
      IF (PRESENT(force_j)) force_j => ijk%force_j
      IF (PRESENT(force_k)) force_k => ijk%force_k

   END SUBROUTINE get_o3c_iterator_info

! **************************************************************************************************
!> \brief ...
!> \param o3c_iterator ...
!> \param mepos ...
!> \return ...
! **************************************************************************************************
   FUNCTION o3c_iterate(o3c_iterator, mepos) RESULT(istat)
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      INTEGER, OPTIONAL                                  :: mepos
      INTEGER                                            :: istat

      INTEGER                                            :: ij, ijpair, klist, me
      TYPE(o3c_container_type), POINTER                  :: o3c

      IF (PRESENT(mepos)) THEN
         me = mepos
      ELSE
         me = 0
      END IF

      !If the neighbors lists are restricted (XAS_TDP), might have nijpairs = 0 on some procs
      IF (o3c_iterator%o3c%nijpairs == 0) THEN
         istat = 1
         RETURN
      END IF

!$OMP CRITICAL(o3c_iterate_critical)
      o3c => o3c_iterator%o3c
      ! we iterate from the last position
      ijpair = o3c_iterator%ijp_last
      klist = o3c_iterator%k_last

      IF (ijpair == 0 .AND. klist == 0) THEN
         ! first step
         istat = 1
         DO ij = 1, o3c%nijpairs
            IF (o3c%ijpair(ij)%nklist > 0) THEN
               o3c_iterator%ijp_thread(me) = ij
               o3c_iterator%k_thread(me) = 1
               istat = 0
               EXIT
            END IF
         END DO
      ELSE IF (ijpair == o3c%nijpairs .AND. klist == o3c%ijpair(ijpair)%nklist) THEN
         ! last step reached
         istat = 1
      ELSE IF (klist == o3c%ijpair(ijpair)%nklist) THEN
         ! last step in this ij list
         istat = 1
         DO ij = ijpair + 1, o3c%nijpairs
            IF (o3c%ijpair(ij)%nklist > 0) THEN
               o3c_iterator%ijp_thread(me) = ij
               o3c_iterator%k_thread(me) = 1
               istat = 0
               EXIT
            END IF
         END DO
      ELSE
         ! increase klist
         o3c_iterator%ijp_thread(me) = ijpair
         o3c_iterator%k_thread(me) = klist + 1
         istat = 0
      END IF

      IF (istat == 0) THEN
         ! set last to this thread
         o3c_iterator%ijp_last = o3c_iterator%ijp_thread(me)
         o3c_iterator%k_last = o3c_iterator%k_thread(me)
      ELSE
         ! set last to final position
         o3c_iterator%ijp_last = o3c%nijpairs
         o3c_iterator%k_last = o3c%ijpair(o3c%nijpairs)%nklist
      END IF
!$OMP END CRITICAL(o3c_iterate_critical)

   END FUNCTION o3c_iterate

! **************************************************************************************************
!> \brief ...
!> \param o3c_iterator ...
!> \param mepos ...
!> \param integral ...
!> \param tvec ...
!> \param force_i ...
!> \param force_j ...
!> \param force_k ...
! **************************************************************************************************
   SUBROUTINE set_o3c_container(o3c_iterator, mepos, integral, tvec, force_i, force_j, force_k)
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      INTEGER, OPTIONAL                                  :: mepos
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: integral
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: tvec, force_i, force_j, force_k

      INTEGER                                            :: ij, k, me
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_int_type), POINTER                        :: ijk
      TYPE(o3c_pair_type), POINTER                       :: ijp

      IF (PRESENT(mepos)) THEN
         me = mepos
      ELSE
         me = 0
      END IF

      ij = o3c_iterator%ijp_thread(me)
      k = o3c_iterator%k_thread(me)

      o3c => o3c_iterator%o3c
      ijp => o3c%ijpair(ij)
      ijk => ijp%ijk(k)

      IF (PRESENT(integral)) ijk%integral => integral
      IF (PRESENT(tvec)) ijk%tvec => tvec
      IF (PRESENT(force_i)) ijk%force_i => force_i
      IF (PRESENT(force_j)) ijk%force_j => force_j
      IF (PRESENT(force_k)) ijk%force_k => force_k

   END SUBROUTINE set_o3c_container

! **************************************************************************************************
!> \brief ...
!> \param o3c_vec ...
!> \param nsize ...
! **************************************************************************************************
   SUBROUTINE o3c_vec_create(o3c_vec, nsize)
      TYPE(o3c_vec_type), DIMENSION(:)                   :: o3c_vec
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nsize

      INTEGER                                            :: i, m, n

      m = SIZE(o3c_vec)
      CPASSERT(SIZE(nsize) == m)

      DO i = 1, m
         n = nsize(i)
         ALLOCATE (o3c_vec(i)%v(n))
         o3c_vec(i)%v = 0.0_dp
         o3c_vec(i)%n = n
      END DO

   END SUBROUTINE o3c_vec_create

! **************************************************************************************************
!> \brief ...
!> \param o3c_vec ...
! **************************************************************************************************
   SUBROUTINE o3c_vec_release(o3c_vec)
      TYPE(o3c_vec_type), DIMENSION(:)                   :: o3c_vec

      INTEGER                                            :: i

      DO i = 1, SIZE(o3c_vec)
         IF (ASSOCIATED(o3c_vec(i)%v)) THEN
            DEALLOCATE (o3c_vec(i)%v)
         END IF
      END DO

   END SUBROUTINE o3c_vec_release

! **************************************************************************************************
!> \brief ...
!> \param o3c_vec ...
!> \param i ...
!> \param vec ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE get_o3c_vec(o3c_vec, i, vec, n)
      TYPE(o3c_vec_type), DIMENSION(:)                   :: o3c_vec
      INTEGER, INTENT(IN)                                :: i
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: vec
      INTEGER, OPTIONAL                                  :: n

      CPASSERT(i > 0 .AND. i <= SIZE(o3c_vec))

      IF (PRESENT(vec)) vec => o3c_vec(i)%v
      IF (PRESENT(n)) n = o3c_vec(i)%n

   END SUBROUTINE get_o3c_vec

! **************************************************************************************************

END MODULE qs_o3c_types
