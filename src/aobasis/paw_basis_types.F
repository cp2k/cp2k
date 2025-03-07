!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
!> \author JGH (9.2022)
! **************************************************************************************************
MODULE paw_basis_types

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE orbital_pointers,                ONLY: nso,&
                                              nsoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_basis_types'

   PUBLIC :: get_paw_basis_info

CONTAINS

! **************************************************************************************************
!> \brief Return some info on the PAW basis derived from a GTO basis set
!> \param basis_1c The parent GTO basis set
!> \param o2nindex ...
!> \param n2oindex ...
!> \param nsatbas ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_paw_basis_info(basis_1c, o2nindex, n2oindex, nsatbas)

      TYPE(gto_basis_set_type), INTENT(IN)               :: basis_1c
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: o2nindex, n2oindex
      INTEGER, INTENT(OUT), OPTIONAL                     :: nsatbas

      INTEGER                                            :: ico, ipgf, iset, iso, iso_pgf, iso_set, &
                                                            k, lx, maxso, nset, nsox
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      LOGICAL                                            :: n2o, nsa, o2n

      o2n = PRESENT(o2nindex)
      n2o = PRESENT(n2oindex)
      nsa = PRESENT(nsatbas)

      IF (o2n .OR. n2o .OR. nsa) THEN
         CALL get_gto_basis_set(gto_basis_set=basis_1c, &
                                nset=nset, lmax=lmax, lmin=lmin, npgf=npgf, maxso=maxso)

         ! Index transformation OLD-NEW
         IF (o2n) THEN
            ALLOCATE (o2nindex(maxso*nset))
            o2nindex = 0
         END IF
         IF (n2o) THEN
            ALLOCATE (n2oindex(maxso*nset))
            n2oindex = 0
         END IF

         ico = 1
         DO iset = 1, nset
            iso_set = (iset - 1)*maxso + 1
            nsox = nsoset(lmax(iset))
            DO ipgf = 1, npgf(iset)
               iso_pgf = iso_set + (ipgf - 1)*nsox
               iso = iso_pgf + nsoset(lmin(iset) - 1)
               DO lx = lmin(iset), lmax(iset)
                  DO k = 1, nso(lx)
                     IF (n2o) n2oindex(ico) = iso
                     IF (o2n) o2nindex(iso) = ico
                     iso = iso + 1
                     ico = ico + 1
                  END DO
               END DO
            END DO
         END DO
         IF (nsa) nsatbas = ico - 1
      END IF

   END SUBROUTINE get_paw_basis_info

END MODULE paw_basis_types
