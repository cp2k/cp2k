!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      HAF (16-Apr-2025) : Import into CP2K
!> \author HAF and yury-lysogorskiy and ralf-drautz
! **************************************************************************************************

MODULE ace_nlist

   USE ace_wrapper,                     ONLY: ace_model_compute
   USE cell_types,                      ONLY: cell_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: ace_data_type,&
                                              fist_nonbond_env_get,&
                                              fist_nonbond_env_type,&
                                              pos_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: angstrom
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC ace_interface

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ace_nlist'

CONTAINS

!
!-------------------------------------------------------------------------------------

! **************************************************************************************************
!> \brief ...
!> \param ace_natom ...
!> \param ace_atype ...
!> \param pot_ace ...
!> \param ace_force ...
!> \param ace_virial ...
!> \param fist_nonbond_env ...
!> \param cell ...
!> \param ace_data ...
! **************************************************************************************************
   SUBROUTINE ace_interface(ace_natom, ace_atype, pot_ace, ace_force, ace_virial, &
                            fist_nonbond_env, cell, ace_data)

      INTEGER, INTENT(IN)                                :: ace_natom, ace_atype(1:ace_natom)
      REAL(kind=dp), INTENT(OUT)                         :: pot_ace, ace_force(1:3, 1:ace_natom), &
                                                            ace_virial(1:6)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(ace_data_type), POINTER                       :: ace_data

#if defined(__ACE)
      INTEGER                                            :: atom_a, atom_b, counter, ilist, k, m, n, &
                                                            natom, nghost, num_update
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ghostidx, listidx
      REAL(KIND=8), ALLOCATABLE                          :: forceunroll(:)
      REAL(kind=dp)                                      :: cell_v(3), dv(1:3), energy(1:ace_natom)
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update_pbc

      natom = ace_natom

      CALL fist_nonbond_env_get(fist_nonbond_env, nonbonded=nonbonded, &
                                r_last_update_pbc=r_last_update_pbc, &
                                num_update=num_update, counter=counter)

      IF ((counter == 1) .OR. (ace_data%refupdate /= num_update)) THEN
         ! fist neighborlist are new:
         ace_data%refupdate = num_update

         IF (.NOT. ALLOCATED(ace_data%neiat)) THEN
            ALLOCATE (ace_data%neiat(0:natom))
         ELSE
            CPASSERT(SIZE(ace_data%neiat) > natom)
         END IF

         !if more than one mpi task, some neighorlists might be empty
         !do we need to check for lone atoms?
         ALLOCATE (ghostidx(natom), listidx(natom))
         nghost = 0
         ace_data%neiat(:) = 0
         ace_data%nei = 0
         DO n = 1, SIZE(nonbonded%neighbor_kind_pairs)
            neighbor_kind_pair => nonbonded%neighbor_kind_pairs(n)
            IF (neighbor_kind_pair%npairs > 0) THEN
               IF ((neighbor_kind_pair%cell_vector(1) == 0) .AND. &
                   (neighbor_kind_pair%cell_vector(2) == 0) .AND. &
                   (neighbor_kind_pair%cell_vector(3) == 0)) THEN
                  DO ilist = 1, natom
                     ghostidx(ilist) = ilist
                  END DO
               ELSE
                  ghostidx = 0
               END IF
               DO ilist = 1, neighbor_kind_pair%npairs
                  atom_a = ace_data%inverse_index_map(neighbor_kind_pair%list(1, ilist))
                  atom_b = ace_data%inverse_index_map(neighbor_kind_pair%list(2, ilist))
                  IF ((atom_a == 0) .OR. (atom_b == 0)) CYCLE
                  ace_data%neiat(atom_a) = ace_data%neiat(atom_a) + 1
                  IF (ghostidx(atom_b) == 0) THEN
                     ! new ghost
                     nghost = nghost + 1
                     ghostidx(atom_b) = nghost + natom
                  END IF
               END DO
            END IF
         END DO
         ! build running sum
         DO n = 1, natom
            ace_data%neiat(n) = ace_data%neiat(n) + ace_data%neiat(n - 1)
         END DO
         ace_data%nei = ace_data%neiat(natom)
         ace_data%nghost = nghost

         IF (ALLOCATED(ace_data%nlist)) THEN
            IF (SIZE(ace_data%nlist) < ace_data%nei) THEN
               DEALLOCATE (ace_data%nlist)
               ALLOCATE (ace_data%nlist(1:ace_data%nei))
            END IF
         ELSE
            ALLOCATE (ace_data%nlist(1:ace_data%nei))
         END IF

         IF (ALLOCATED(ace_data%attype)) THEN
            IF (SIZE(ace_data%attype) < natom + nghost) THEN
               DEALLOCATE (ace_data%atpos)
               DEALLOCATE (ace_data%attype)
               DEALLOCATE (ace_data%origin)
               DEALLOCATE (ace_data%shift)
               ALLOCATE (ace_data%atpos(1:3, 1:natom + nghost))
               ALLOCATE (ace_data%attype(1:natom + nghost))
               ALLOCATE (ace_data%origin(1:natom + nghost))
               ALLOCATE (ace_data%shift(1:3, 1:natom + nghost))
            END IF
         ELSE
            ALLOCATE (ace_data%atpos(1:3, 1:natom + nghost))
            ALLOCATE (ace_data%attype(1:natom + nghost))
            ALLOCATE (ace_data%origin(1:natom + nghost))
            ALLOCATE (ace_data%shift(1:3, 1:natom + nghost))
         END IF
         ace_data%attype(1:natom) = ace_atype(:)

         DO n = 1, natom
            ace_data%atpos(:, n) = r_last_update_pbc(ace_data%use_indices(n))%r*angstrom
            ace_data%origin(n) = n
         END DO
         ace_data%shift(:, :) = 0

         k = natom
         listidx(1:natom) = ace_data%neiat(0:natom - 1)
         DO n = 1, SIZE(nonbonded%neighbor_kind_pairs)
            neighbor_kind_pair => nonbonded%neighbor_kind_pairs(n)
            IF (neighbor_kind_pair%npairs > 0) THEN
               IF ((neighbor_kind_pair%cell_vector(1) == 0) .AND. &
                   (neighbor_kind_pair%cell_vector(2) == 0) .AND. &
                   (neighbor_kind_pair%cell_vector(3) == 0)) THEN
                  DO m = 1, natom
                     ghostidx(m) = m
                  END DO
               ELSE
                  ghostidx = 0
               END IF
               dv = REAL(neighbor_kind_pair%cell_vector, KIND=dp)
               ! dimensions it's not periodic with should be zero in cell_vector
               ! so should always be valid:
               cell_v = MATMUL(cell%hmat, dv)*angstrom
               DO ilist = 1, neighbor_kind_pair%npairs
                  atom_a = ace_data%inverse_index_map(neighbor_kind_pair%list(1, ilist))
                  atom_b = ace_data%inverse_index_map(neighbor_kind_pair%list(2, ilist))
                  IF ((atom_a == 0) .OR. (atom_b == 0)) CYCLE
                  IF (ghostidx(atom_b) == 0) THEN
                     k = k + 1
                     ace_data%atpos(:, k) = ace_data%atpos(:, atom_b) + cell_v
                     ace_data%shift(:, k) = neighbor_kind_pair%cell_vector
                     ace_data%origin(k) = atom_b
                     ace_data%attype(k) = ace_atype(atom_b)
                     ghostidx(atom_b) = k
                  END IF
                  listidx(atom_a) = listidx(atom_a) + 1
                  ace_data%nlist(listidx(atom_a)) = ghostidx(atom_b)
               END DO
            END IF
         END DO

         DEALLOCATE (ghostidx)
         DEALLOCATE (listidx)

!         transforming to c call
!     -> atom index will change from 1...n to 0...n-1 -> subtract 1 from neighborlist
         ace_data%nlist(1:ace_data%nei) = ace_data%nlist(1:ace_data%nei) - 1
         ace_data%origin(1:natom + nghost) = ace_data%origin(1:natom + nghost) - 1
!-----------------------------------------------

      ELSE ! no changes in neighbor list just update positions:
         nghost = ace_data%nghost
         DO n = 1, natom
            ace_data%atpos(:, n) = r_last_update_pbc(ace_data%use_indices(n))%r*angstrom
         END DO
         DO n = natom + 1, natom + nghost
            dv = REAL(ace_data%shift(:, n), KIND=dp)*angstrom
            !origin+1 since we already changed to C notation for origin:
            ace_data%atpos(:, n) = ace_data%atpos(:, ace_data%origin(n) + 1) + MATMUL(cell%hmat, dv)
         END DO
      END IF

! -> force array
      ALLOCATE (forceunroll(1:3*natom))
      forceunroll = 0.0
      pot_ace = 0.0
      ace_virial = 0.0

      CALL ace_model_compute( &
         natomc=natom, &
         nghostc=nghost, &
         neic=ace_data%nei, &
         neiatc=ace_data%neiat, &
         originc=ace_data%origin, &
         nlistc=ace_data%nlist, &
         attypec=ace_data%attype, &
         atposc=RESHAPE(ace_data%atpos, (/3*(natom + nghost)/)), &
         forcec=forceunroll, &
         virialc=ace_virial, &
         energyc=energy, &
         model=ace_data%model)

      ace_force = RESHAPE(forceunroll, (/3, natom/))
      pot_ace = SUM(energy(1:natom))

      DEALLOCATE (forceunroll)

#else
      MARK_USED(ace_natom)
      MARK_USED(ace_atype)
      MARK_USED(pot_ace)
      MARK_USED(ace_force)
      MARK_USED(ace_virial)
      MARK_USED(fist_nonbond_env)
      MARK_USED(cell)
      MARK_USED(ace_data)
      CPABORT("CP2K was compiled without ACE library.")
#endif

   END SUBROUTINE ace_interface

!----------------------------------------------------------------------------------

END MODULE ace_nlist
