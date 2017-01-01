!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialization for solid harmonic Gaussian (SHG) integral scheme. Scheme for calculation
!>        of contracted, spherical Gaussian integrals using the solid harmonics. Initialization of
!>        the contraction matrices
!> \par History
!>      created [08.2016]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE generic_shg_integrals_init
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              ifac,&
                                              pi
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                              clebsch_gordon_deallocate,&
                                              clebsch_gordon_init
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'generic_shg_integrals_init'

   PUBLIC :: contraction_matrix_shg, contraction_matrix_shg_mix, contraction_matrix_shg_rx2m, &
             get_clebsch_gordon_coefficients

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief contraction matrix for SHG integrals
!> \param basis ...
!> \param scon_shg contraction matrix
! **************************************************************************************************
   SUBROUTINE contraction_matrix_shg(basis, scon_shg)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_shg

      CHARACTER(len=*), PARAMETER :: routineN = 'contraction_matrix_shg', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, ishell, l, maxpgf, maxshell, &
                                                            nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nshell
      REAL(KIND=dp)                                      :: aif, gcc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: norm
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet

      nset = basis%nset
      npgf => basis%npgf
      nshell => basis%nshell
      zet => basis%zet

      maxpgf = SIZE(basis%gcc, 1)
      maxshell = SIZE(basis%gcc, 2)
      ALLOCATE (norm(basis%nset, maxshell))
      ALLOCATE (scon_shg(maxpgf, maxshell, nset))
      scon_shg = 0.0_dp

      CALL basis_norm_shg(basis, norm)

      DO iset = 1, nset
         DO ishell = 1, nshell(iset)
            l = basis%l(ishell, iset)
            DO ipgf = 1, npgf(iset)
               aif = 1.0_dp/((2._dp*zet(ipgf, iset))**l)
               gcc = basis%gcc(ipgf, ishell, iset)
               scon_shg(ipgf, ishell, iset) = norm(iset, ishell)*gcc*aif
            END DO
         END DO
      END DO

      DEALLOCATE (norm)

   END SUBROUTINE contraction_matrix_shg

!***************************************************************************************************
!> \brief normalization solid harmonic Gaussians (SHG)
!> \param basis ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_norm_shg(basis, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_shg', routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, ishell, jpgf, l
      REAL(KIND=dp)                                      :: aai, aaj, cci, ccj, expa, ppl

      norm = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            expa = 0.5_dp*REAL(2*l+3, dp)
            ppl = fac(2*l+2)*pi**(1.5_dp)/fac(l+1)
            ppl = ppl/(2._dp**REAL(2*l+1, dp))
            ppl = ppl/REAL(2*l+1, dp)
            DO ipgf = 1, basis%npgf(iset)
               cci = basis%gcc(ipgf, ishell, iset)
               aai = basis%zet(ipgf, iset)
               DO jpgf = 1, basis%npgf(iset)
                  ccj = basis%gcc(jpgf, ishell, iset)
                  aaj = basis%zet(jpgf, iset)
                  norm(iset, ishell) = norm(iset, ishell)+cci*ccj*ppl/(aai+aaj)**expa
               END DO
            END DO
            norm(iset, ishell) = 1.0_dp/SQRT(norm(iset, ishell))
         END DO
      END DO

   END SUBROUTINE basis_norm_shg

! **************************************************************************************************
!> \brief mixed contraction matrix for SHG integrals [aba] and [abb] for orbital and ri basis
!>        at the same atom
!> \param orb_basis orbital basis
!> \param ri_basis   ...
!> \param orb_index index for orbital basis
!> \param ri_index index for ri basis
!> \param scon_mix mixed contraction matrix
! **************************************************************************************************
   SUBROUTINE contraction_matrix_shg_mix(orb_basis, ri_basis, orb_index, ri_index, scon_mix)

      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      INTEGER, DIMENSION(:, :, :), POINTER               :: orb_index, ri_index
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: scon_mix

      CHARACTER(len=*), PARAMETER :: routineN = 'contraction_matrix_shg_mix', &
         routineP = moduleN//':'//routineN

      INTEGER :: forb, fri, iil, il, ipgf, iset, ishell, jpgf, jset, jshell, l, l1, l2, lmax_orb, &
         lmax_ri, maxpgf_orb, maxpgf_ri, maxshell_orb, maxshell_ri, nf_orb, nf_ri, nl, nl_max, &
         nset_orb, nset_ri
      INTEGER, DIMENSION(:), POINTER                     :: npgf_orb, npgf_ri, nshell_orb, nshell_ri
      REAL(KIND=dp)                                      :: cjf, const, const1, const2, gcc_orb, &
                                                            gcc_ri, prefac, scon1, scon2, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: shg_fac
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: norm_orb, norm_ri, zet_orb, zet_ri

      nset_orb = orb_basis%nset
      npgf_orb => orb_basis%npgf
      nshell_orb => orb_basis%nshell
      zet_orb => orb_basis%zet
      nset_ri = ri_basis%nset
      npgf_ri => ri_basis%npgf
      nshell_ri => ri_basis%nshell
      zet_ri => ri_basis%zet

      maxpgf_orb = SIZE(orb_basis%gcc, 1)
      maxshell_orb = SIZE(orb_basis%gcc, 2)
      ALLOCATE (norm_orb(nset_orb, maxshell_orb))
      maxpgf_ri = SIZE(ri_basis%gcc, 1)
      maxshell_ri = SIZE(ri_basis%gcc, 2)
      ALLOCATE (norm_ri(nset_ri, maxshell_ri))

      CALL basis_norm_shg(orb_basis, norm_orb)
      CALL basis_norm_shg(ri_basis, norm_ri)

      ALLOCATE (orb_index(maxpgf_orb, maxshell_orb, nset_orb))
      ALLOCATE (ri_index(maxpgf_ri, maxshell_ri, nset_ri))

      !** index orbital basis set
      nf_orb = 0
      DO iset = 1, nset_orb
         DO ishell = 1, nshell_orb(iset)
            DO ipgf = 1, npgf_orb(iset)
               nf_orb = nf_orb+1
               orb_index(ipgf, ishell, iset) = nf_orb
            END DO
         END DO
      END DO

      !** index ri basis set
      nf_ri = 0
      DO iset = 1, nset_ri
         DO ishell = 1, nshell_ri(iset)
            DO ipgf = 1, npgf_ri(iset)
               nf_ri = nf_ri+1
               ri_index(ipgf, ishell, iset) = nf_ri
            END DO
         END DO
      END DO

      lmax_orb = MAXVAL(orb_basis%lmax)
      lmax_ri = MAXVAL(ri_basis%lmax)
      nl_max = INT((lmax_orb+lmax_ri)/2)+1
      ALLOCATE (scon_mix(nl_max, nf_ri, nf_orb, nl_max))
      scon_mix = 0.0_dp

      ALLOCATE (shg_fac(0:nl_max-1))
      shg_fac(0) = 1.0_dp

      DO iset = 1, nset_orb
         DO ishell = 1, nshell_orb(iset)
            l1 = orb_basis%l(ishell, iset)
            const1 = SQRT(1.0_dp/REAL(2*l1+1, dp))
            DO jset = 1, nset_ri
               DO jshell = 1, nshell_ri(jset)
                  l2 = ri_basis%l(jshell, jset)
                  const2 = SQRT(1.0_dp/REAL(2*l2+1, dp))
                  nl = INT((l1+l2)/2)
                  IF (l1 == 0 .OR. l2 == 0) nl = 0
                  DO il = 0, nl
                     l = l1+l2-2*il
                     const = const1*const2*2.0_dp*SQRT(pi*REAL(2*l+1, dp))
                     DO iil = 1, il
                        shg_fac(iil) = fac(l+iil-1)*ifac(l)*REAL(l, dp) &
                                       *fac(il)/fac(il-iil)/fac(iil)
                     ENDDO
                     DO ipgf = 1, npgf_orb(iset)
                        forb = orb_index(ipgf, ishell, iset)
                        gcc_orb = orb_basis%gcc(ipgf, ishell, iset)
                        scon1 = norm_orb(iset, ishell)*gcc_orb
                        DO jpgf = 1, npgf_ri(jset)
                           fri = ri_index(jpgf, jshell, jset)
                           gcc_ri = ri_basis%gcc(jpgf, jshell, jset)
                           scon2 = norm_ri(jset, jshell)*gcc_ri
                           zet = zet_orb(ipgf, iset)+zet_ri(jpgf, jset)
                           cjf = 1.0_dp/((2._dp*zet)**l)
                           prefac = const*cjf*scon1*scon2
                           DO iil = 0, il
                              scon_mix(iil+1, fri, forb, il+1) = prefac*shg_fac(iil)/zet**iil
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (norm_orb, norm_ri, shg_fac)

   END SUBROUTINE contraction_matrix_shg_mix

! **************************************************************************************************
!> \brief ...
!> \param basis ...
!> \param m ...
!> \param scon_shg ...
!> \param scon_rx2m ...
! **************************************************************************************************
   SUBROUTINE contraction_matrix_shg_rx2m(basis, m, scon_shg, scon_rx2m)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: scon_shg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: scon_rx2m

      CHARACTER(len=*), PARAMETER :: routineN = 'contraction_matrix_shg_rx2m', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, ishell, j, l, maxpgf, &
                                                            maxshell, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nshell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: shg_fac
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet

      npgf => basis%npgf
      nshell => basis%nshell
      zet => basis%zet
      nset = basis%nset

      maxpgf = SIZE(basis%gcc, 1)
      maxshell = SIZE(basis%gcc, 2)
      ALLOCATE (scon_rx2m(maxpgf, m+1, maxshell, nset))
      scon_rx2m = 0.0_dp
      ALLOCATE (shg_fac(0:m))
      shg_fac(0) = 1.0_dp

      DO iset = 1, nset
         DO ishell = 1, nshell(iset)
            l = basis%l(ishell, iset)
            DO j = 1, m
               shg_fac(j) = fac(l+j-1)*ifac(l)*REAL(l, dp) &
                            *fac(m)/fac(m-j)/fac(j)
            ENDDO
            DO ipgf = 1, npgf(iset)
               DO j = 0, m
                  scon_rx2m(ipgf, j+1, ishell, iset) = scon_shg(ipgf, ishell, iset) &
                                                       *shg_fac(j)/zet(ipgf, iset)**j
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (shg_fac)

   END SUBROUTINE contraction_matrix_shg_rx2m

! **************************************************************************************************
!> \brief calculate the Clebsch-Gordon (CG) coefficients for expansion of the
!>        product of two spherical harmonic Gaussians
!> \param my_cg matrix storing CG coefficients
!> \param cg_none0_list list of none-zero CG coefficients
!> \param ncg_none0 number of none-zero CG coefficients
!> \param maxl1 maximal l quantum number of 1st spherical function
!> \param maxl2 maximal l quantum number of 2nd spherical function
! **************************************************************************************************
   SUBROUTINE get_clebsch_gordon_coefficients(my_cg, cg_none0_list, ncg_none0, &
                                              maxl1, maxl2)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: my_cg
      INTEGER, DIMENSION(:, :, :), POINTER               :: cg_none0_list
      INTEGER, DIMENSION(:, :), POINTER                  :: ncg_none0
      INTEGER, INTENT(IN)                                :: maxl1, maxl2

      CHARACTER(len=*), PARAMETER :: routineN = 'get_clebsch_gordon_coefficients', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: il, ilist, iso, iso1, iso2, l1, l1l2, &
                                                            l2, lc1, lc2, lp, m1, m2, maxl, mm, &
                                                            mp, nlist, nlist_max, nsfunc, nsfunc1, &
                                                            nsfunc2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rga

      nlist_max = 6
      nsfunc1 = nsoset(maxl1)
      nsfunc2 = nsoset(maxl2)
      maxl = maxl1+maxl2
      nsfunc = nsoset(maxl)

      CALL clebsch_gordon_init(maxl)

      ALLOCATE (my_cg(nsfunc1, nsfunc2, nsfunc))
      my_cg = 0.0_dp
      ALLOCATE (ncg_none0(nsfunc1, nsfunc2))
      ncg_none0 = 0
      ALLOCATE (cg_none0_list(nsfunc1, nsfunc2, nlist_max))
      cg_none0_list = 0

      ALLOCATE (rga(maxl, 2))
      rga = 0.0_dp
      DO lc1 = 0, maxl1
         DO iso1 = nsoset(lc1-1)+1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxl2
               DO iso2 = nsoset(lc2-1)+1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  l1l2 = l1+l2
                  mp = m1+m2
                  mm = m1-m2
                  IF (m1*m2 < 0 .OR. (m1*m2 == 0 .AND. (m1 < 0 .OR. m2 < 0))) THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1+l2, 2), l1l2, 2
                     il = lp/2+1
                     IF (ABS(mp) <= lp) THEN
                        IF (mp >= 0) THEN
                           iso = nsoset(lp-1)+lp+1+mp
                        ELSE
                           iso = nsoset(lp-1)+lp+1-ABS(mp)
                        END IF
                        my_cg(iso1, iso2, iso) = rga(il, 1)
                     ENDIF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                        IF (mm >= 0) THEN
                           iso = nsoset(lp-1)+lp+1+mm
                        ELSE
                           iso = nsoset(lp-1)+lp+1-ABS(mm)
                        END IF
                        my_cg(iso1, iso2, iso) = rga(il, 2)
                     ENDIF
                  END DO
                  nlist = 0
                  DO ilist = 1, nsfunc
                     IF (ABS(my_cg(iso1, iso2, ilist)) > 1.E-8_dp) THEN
                        nlist = nlist+1
                        IF (nlist > nlist_max) THEN
                           CALL reallocate(cg_none0_list, 1, nsfunc1, 1, nsfunc2, 1, nlist)
                           nlist_max = nlist
                        ENDIF
                        cg_none0_list(iso1, iso2, nlist) = ilist
                     ENDIF
                  ENDDO
                  ncg_none0(iso1, iso2) = nlist
               ENDDO ! iso2
            ENDDO ! lc2
         ENDDO ! iso1
      ENDDO ! lc1

      DEALLOCATE (rga)
      CALL clebsch_gordon_deallocate()

   END SUBROUTINE get_clebsch_gordon_coefficients

END MODULE generic_shg_integrals_init
