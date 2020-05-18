!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
!> \author JGH (11.2017)
! **************************************************************************************************
MODULE aux_basis_set

   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lapack,                          ONLY: lapack_spotrf
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset,&
                                              nso
   USE orbital_symbols,                 ONLY: cgf_symbol,&
                                              sgf_symbol
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'aux_basis_set'

! *** Public subroutines ***

   PUBLIC :: create_aux_basis

CONTAINS

! **************************************************************************************************
!> \brief create a basis in GTO form
!> \param aux_basis ...
!> \param bsname ...
!> \param nsets ...
!> \param lmin ...
!> \param lmax ...
!> \param nl ...
!> \param npgf ...
!> \param zet ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE create_aux_basis(aux_basis, bsname, nsets, lmin, lmax, nl, npgf, zet)

      TYPE(gto_basis_set_type), POINTER                  :: aux_basis
      CHARACTER(LEN=default_string_length)               :: bsname
      INTEGER, INTENT(IN)                                :: nsets
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lmin, lmax
      INTEGER, DIMENSION(0:, :), INTENT(IN)              :: nl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npgf
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: zet

      CHARACTER(len=*), PARAMETER :: routineN = 'create_aux_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ico, info, iset, ishell, j, l, &
                                                            lshell, m, maxco, maxpgf, maxshell, &
                                                            ncgf, ns, nsgf, nx
      REAL(KIND=dp)                                      :: za, zb, zetab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: so

      CPASSERT(.NOT. ASSOCIATED(aux_basis))
      ALLOCATE (aux_basis)
      !
      aux_basis%name = bsname
      aux_basis%aliases = bsname
      aux_basis%nset = nsets
      !
      ALLOCATE (aux_basis%npgf(nsets), aux_basis%nshell(nsets), &
                aux_basis%lmax(nsets), aux_basis%lmin(nsets))
      aux_basis%lmax(1:nsets) = lmax(1:nsets)
      aux_basis%lmin(1:nsets) = lmin(1:nsets)
      aux_basis%npgf(1:nsets) = npgf(1:nsets)
      DO iset = 1, nsets
         aux_basis%nshell(iset) = 0
         DO l = lmin(iset), lmax(iset)
            aux_basis%nshell(iset) = aux_basis%nshell(iset) + nl(l, iset)
         END DO
      END DO
      maxpgf = MAXVAL(npgf(1:nsets))
      maxshell = MAXVAL(aux_basis%nshell(1:nsets))
      ALLOCATE (aux_basis%zet(maxpgf, nsets))
      aux_basis%zet(1:maxpgf, 1:nsets) = zet(1:maxpgf, 1:nsets)

      ALLOCATE (aux_basis%n(maxshell, nsets))
      ALLOCATE (aux_basis%l(maxshell, nsets))
      ALLOCATE (aux_basis%gcc(maxpgf, maxshell, nsets))

      DO iset = 1, nsets
         ns = 0
         DO l = lmin(iset), lmax(iset)
            DO i = 1, nl(l, iset)
               ns = ns + 1
               aux_basis%l(ns, iset) = l
               aux_basis%n(ns, iset) = l + i
            END DO
         END DO
      END DO

      ! contraction
      aux_basis%gcc = 0.0_dp
      DO iset = 1, nsets
         ns = 0
         DO l = lmin(iset), lmax(iset)
            nx = aux_basis%npgf(iset)
            ALLOCATE (so(nx, nx))
            CPASSERT(nx >= nl(l, iset))
            DO i = 1, nx
               za = (2.0_dp*zet(i, iset))**(0.25_dp*(2*l + 3))
               DO j = i, nx
                  zb = (2.0_dp*zet(j, iset))**(0.25_dp*(2*l + 3))
                  zetab = zet(i, iset) + zet(j, iset)
                  so(i, j) = za*zb/zetab**(l + 1.5_dp)
                  so(j, i) = so(i, j)
               END DO
            END DO
            info = 0
            CALL lapack_spotrf("U", nx, so, nx, info)
            CPASSERT(info == 0)
            CALL dtrtri("U", "N", nx, so, nx, info)
            CPASSERT(info == 0)
            DO i = ns + 1, ns + nl(l, iset)
               DO j = 1, i - ns
                  aux_basis%gcc(j, i, iset) = so(j, i - ns)
               END DO
            END DO
            IF (nl(l, iset) < nx) THEN
               i = ns + nl(l, iset)
               DO j = nl(l, iset) + 1, nx
                  aux_basis%gcc(j, i, iset) = 1.0_dp
               END DO
            END IF
            ns = ns + nl(l, iset)
            DEALLOCATE (so)
         END DO
      END DO

      ! Initialise the depending aux_basis structures
      ALLOCATE (aux_basis%first_cgf(maxshell, nsets))
      ALLOCATE (aux_basis%first_sgf(maxshell, nsets))
      ALLOCATE (aux_basis%last_cgf(maxshell, nsets))
      ALLOCATE (aux_basis%last_sgf(maxshell, nsets))
      ALLOCATE (aux_basis%ncgf_set(nsets))
      ALLOCATE (aux_basis%nsgf_set(nsets))

      maxco = 0
      ncgf = 0
      nsgf = 0
      DO iset = 1, nsets
         aux_basis%ncgf_set(iset) = 0
         aux_basis%nsgf_set(iset) = 0
         DO ishell = 1, aux_basis%nshell(iset)
            lshell = aux_basis%l(ishell, iset)
            aux_basis%first_cgf(ishell, iset) = ncgf + 1
            ncgf = ncgf + nco(lshell)
            aux_basis%last_cgf(ishell, iset) = ncgf
            aux_basis%ncgf_set(iset) = &
               aux_basis%ncgf_set(iset) + nco(lshell)
            aux_basis%first_sgf(ishell, iset) = nsgf + 1
            nsgf = nsgf + nso(lshell)
            aux_basis%last_sgf(ishell, iset) = nsgf
            aux_basis%nsgf_set(iset) = &
               aux_basis%nsgf_set(iset) + nso(lshell)
         END DO
         maxco = MAX(maxco, npgf(iset)*ncoset(lmax(iset)))
      END DO
      aux_basis%ncgf = ncgf
      aux_basis%nsgf = nsgf

      ALLOCATE (aux_basis%lx(ncgf))
      ALLOCATE (aux_basis%ly(ncgf))
      ALLOCATE (aux_basis%lz(ncgf))
      ALLOCATE (aux_basis%m(nsgf))
      ALLOCATE (aux_basis%cgf_symbol(ncgf))
      ALLOCATE (aux_basis%sgf_symbol(nsgf))

      ncgf = 0
      nsgf = 0

      DO iset = 1, nsets
         DO ishell = 1, aux_basis%nshell(iset)
            lshell = aux_basis%l(ishell, iset)
            DO ico = ncoset(lshell - 1) + 1, ncoset(lshell)
               ncgf = ncgf + 1
               aux_basis%lx(ncgf) = indco(1, ico)
               aux_basis%ly(ncgf) = indco(2, ico)
               aux_basis%lz(ncgf) = indco(3, ico)
               aux_basis%cgf_symbol(ncgf) = &
                  cgf_symbol(aux_basis%n(ishell, iset), (/aux_basis%lx(ncgf), &
                                                          aux_basis%ly(ncgf), &
                                                          aux_basis%lz(ncgf)/))
            END DO
            DO m = -lshell, lshell
               nsgf = nsgf + 1
               aux_basis%m(nsgf) = m
               aux_basis%sgf_symbol(nsgf) = &
                  sgf_symbol(aux_basis%n(ishell, iset), lshell, m)
            END DO
         END DO
      ENDDO

      ! orbital radii (initialize later)
      aux_basis%kind_radius = 0.0_dp
      aux_basis%short_kind_radius = 0.0_dp
      ALLOCATE (aux_basis%set_radius(nsets))
      ALLOCATE (aux_basis%pgf_radius(maxpgf, nsets))
      aux_basis%set_radius = 0.0_dp
      aux_basis%pgf_radius = 0.0_dp

      ! basis transformation matrices
      ALLOCATE (aux_basis%cphi(maxco, ncgf))
      ALLOCATE (aux_basis%sphi(maxco, nsgf))
      ALLOCATE (aux_basis%scon(maxco, nsgf))
      ALLOCATE (aux_basis%norm_cgf(ncgf))
      aux_basis%norm_type = 2
!     CALL init_orb_basis_set(aux_basis)

   END SUBROUTINE create_aux_basis

END MODULE aux_basis_set
