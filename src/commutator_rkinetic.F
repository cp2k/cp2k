!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of commutator of kinetic energy and position operator
!> \par History
!>      JGH: from qs_kinetic
!> \author Juerg Hutter
! **************************************************************************************************
MODULE commutator_rkinetic
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_kinetic,                      ONLY: kinetic
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                              cp_dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'commutator_rkinetic'

! *** Public subroutines ***

   PUBLIC :: build_com_tr_matrix

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of commutator [T,r] over Cartesian Gaussian functions.
!> \param   matrix_tr ...
!> \param   qs_kind_set ...
!> \param   basis_type basis set to be used
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \date    11.10.2010
!> \par     History
!>          Ported from qs_overlap, replaces code in build_core_hamiltonian
!>          Refactoring [07.2014] JGH
!>          Simplify options and use new kinetic energy integral routine
!>          Adapted from qs_kinetic [07.2016]
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE build_com_tr_matrix(matrix_tr, qs_kind_set, basis_type, sab_nl)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_tr
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CHARACTER(len=*), PARAMETER :: routineN = 'build_com_tr_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, icol, ikind, ir, irow, &
                                                            iset, jatom, jkind, jset, ldsab, ltab, &
                                                            mepos, ncoa, ncob, nkind, nseta, &
                                                            nsetb, nthread, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found, trans
      REAL(KIND=dp)                                      :: rab2, tab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: qab, tkab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: kab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kx_block, ky_block, kz_block, rpgfa, &
                                                            rpgfb, scon_a, scon_b, zeta, zetb
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      nkind = SIZE(qs_kind_set)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      ! prepare basis set
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      ! *** Allocate work storage ***
      ldsab = get_memory_usage(qs_kind_set, basis_type)

      nthread = 1
!$    nthread = omp_get_max_threads()
      ! Iterate of neighbor list
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,ldsab,nl_iterator, do_symmetric,&
!$OMP         ncoset,matrix_tr,basis_set_list)  &
!$OMP PRIVATE (kx_block,ky_block,kz_block,mepos,kab,qab,tab,ikind,jkind,iatom,jatom,rab,cell,&
!$OMP          basis_set_a,basis_set_b,&
!$OMP          first_sgfa, la_max, la_min, npgfa, nsgfa, nseta, rpgfa, set_radius_a, ncoa, ncob, &
!$OMP          zeta, first_sgfb, lb_max, lb_min, ltab, npgfb, nsetb, rpgfb, set_radius_b, nsgfb, tkab, &
!$OMP          zetb, scon_a, scon_b, irow, icol, found, trans, rab2, sgfa, sgfb, iset, jset)

      mepos = 0
!$    mepos = omp_get_thread_num()

      ALLOCATE (kab(ldsab, ldsab, 3), qab(ldsab, ldsab))

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rab, cell=cell)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
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
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF
         NULLIFY (kx_block)
         CALL cp_dbcsr_get_block_p(matrix=matrix_tr(1)%matrix, &
                                   row=irow, col=icol, BLOCK=kx_block, found=found)
         CPASSERT(found)
         NULLIFY (ky_block)
         CALL cp_dbcsr_get_block_p(matrix=matrix_tr(2)%matrix, &
                                   row=irow, col=icol, BLOCK=ky_block, found=found)
         CPASSERT(found)
         NULLIFY (kz_block)
         CALL cp_dbcsr_get_block_p(matrix=matrix_tr(3)%matrix, &
                                   row=irow, col=icol, BLOCK=kz_block, found=found)
         CPASSERT(found)

         rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
         tab = SQRT(rab2)
         trans = do_symmetric .AND. (iatom > jatom)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset)+set_radius_b(jset) < tab) CYCLE

               ncob = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
               sgfb = first_sgfb(1, jset)

               ! calclulate integrals
               ltab = MAX(npgfa(iset)*ncoset(la_max(iset)+1), npgfb(jset)*ncoset(lb_max(jset)+1))
               ALLOCATE (tkab(ltab, ltab))
               CALL kinetic(la_max(iset)+1, la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                            lb_max(jset)+1, lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                            rab, tkab)
               ! commutator
               CALL comab_opr(la_max(iset), npgfa(iset), rpgfa(:, iset), la_min(iset), &
                              lb_max(jset), npgfb(jset), rpgfb(:, jset), lb_min(jset), &
                              tab, tkab, kab)
               DEALLOCATE (tkab)
               ! Contraction step
               DO ir = 1, 3
                  CALL contraction(kab(:, :, ir), qab, ca=scon_a(:, sgfa:), na=ncoa, ma=nsgfa(iset), &
                                   cb=scon_b(:, sgfb:), nb=ncob, mb=nsgfb(jset), trans=trans)
!$OMP CRITICAL(blockadd)
                  SELECT CASE (ir)
                  CASE (1)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), kx_block, sgfa, sgfb, trans=trans)
                  CASE (2)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), ky_block, sgfa, sgfb, trans=trans)
                  CASE (3)
                     CALL block_add("IN", qab, nsgfa(iset), nsgfb(jset), kz_block, sgfa, sgfb, trans=trans)
                  END SELECT
!$OMP END CRITICAL(blockadd)
               END DO

            END DO
         END DO

      END DO
      DEALLOCATE (kab, qab)
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_com_tr_matrix

! **************************************************************************************************
!> \brief   Calculate the commutator [O,r] from the integrals [a|O|b].
!>          We assume that on input all integrals [a+1|O|b+1] are available.
!>          [a|[O,ri]|b] = [a|O|b+1i] - [a+1i|O|b]
!> \param la_max ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param dab ...
!> \param ab ...
!> \param comabr ...
!>
!> \date    25.07.2016
!> \par Literature
!>          S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - ncoset    : Number of orbitals in a Cartesian orbital set.
!>      - npgf{a,b} : Degree of contraction of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rab2      : Square of the distance between the atomic centers a and b.
!>      - rac       : Distance vector between the atomic centers a and c.
!>      - rac2      : Square of the distance between the atomic centers a and c.
!>      - rbc       : Distance vector between the atomic centers b and c.
!>      - rbc2      : Square of the distance between the atomic centers b and c.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
!>
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE comab_opr(la_max, npgfa, rpgfa, la_min, lb_max, npgfb, rpgfb, lb_min, &
                        dab, ab, comabr)
      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: comabr

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, coap, &
                                                            coapx, coapy, coapz, cob, cobp, cobpx, &
                                                            cobpy, cobpz, ipgf, jpgf, la, lb, na, &
                                                            nap, nb, nbp, ofa, ofb

      comabr = 0.0_dp

      ofa = ncoset(la_min-1)
      ofb = ncoset(lb_min-1)

      na = 0
      nap = 0
      DO ipgf = 1, npgfa
         nb = 0
         nbp = 0
         DO jpgf = 1, npgfb
            IF (rpgfa(ipgf)+rpgfb(jpgf) > dab) THEN
               DO la = la_min, la_max
                  DO ax = 0, la
                     DO ay = 0, la-ax
                        az = la-ax-ay
                        coa = na+coset(ax, ay, az)-ofa
                        coap = nap+coset(ax, ay, az)-ofa
                        coapx = nap+coset(ax+1, ay, az)-ofa
                        coapy = nap+coset(ax, ay+1, az)-ofa
                        coapz = nap+coset(ax, ay, az+1)-ofa
                        DO lb = lb_min, lb_max
                           DO bx = 0, lb
                              DO by = 0, lb-bx
                                 bz = lb-bx-by
                                 cob = nb+coset(bx, by, bz)-ofb
                                 cobp = nbp+coset(bx, by, bz)-ofb
                                 cobpx = nbp+coset(bx+1, by, bz)-ofb
                                 cobpy = nbp+coset(bx, by+1, bz)-ofb
                                 cobpz = nbp+coset(bx, by, bz+1)-ofb
                                 ! [a|[O,ri]|b] = [a|O|b+1i] - [a+1i|O|b]
                                 comabr(coa, cob, 1) = ab(coap, cobpx)-ab(coapx, cobp)
                                 comabr(coa, cob, 2) = ab(coap, cobpy)-ab(coapy, cobp)
                                 comabr(coa, cob, 3) = ab(coap, cobpz)-ab(coapz, cobp)
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            nb = nb+ncoset(lb_max)-ofb
            nbp = nbp+ncoset(lb_max+1)-ofb
         END DO
         na = na+ncoset(la_max)-ofa
         nap = nap+ncoset(la_max+1)-ofa
      END DO

   END SUBROUTINE comab_opr

END MODULE commutator_rkinetic

