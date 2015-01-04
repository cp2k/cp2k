!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Calculation of the nuclear attraction contribution to the core Hamiltonian
!>         <a|erfc|b> :we only calculate the non-screened part
!> \par History
!>      - core_ppnl refactored from qs_core_hamiltonian [Joost VandeVondele, 2008-11-01]
!>      - adapted for nuclear attraction [jhu, 2009-02-24]
! *****************************************************************************
MODULE core_ae
  USE ai_verfc,                        ONLY: verfc
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type,&
                                             nl_set_sub_iterator,&
                                             nl_sub_iterate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ae'

  PUBLIC :: build_core_ae

CONTAINS

!==========================================================================================================

! *****************************************************************************
!> \brief ...
!> \param matrix_h ...
!> \param matrix_p ...
!> \param force ...
!> \param virial ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param nder ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param sab_orb ...
!> \param sac_ae ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_core_ae(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ae, error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL                                  :: use_virial
    INTEGER                                  :: nder
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_ae
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ae', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, handle, iatom, icol, ikind, inode, &
      irow, iset, jatom, jkind, jset, katom, kkind, last_jatom, ldai, ldsab, &
      maxco, maxder, maxl, maxnset, maxsgf, na_plus, natom, nb_plus, ncoa, &
      ncob, nij, nkind, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b
    REAL(KIND=dp)                            :: alpha_c, core_charge, &
                                                core_radius, dab, dac, dbc, &
                                                f0, rab2, rac2, rbc2, zeta_c
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: habd, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hab, pab, verf, vnuc
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, rac, &
                                                rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, rpgfa, &
                                                rpgfb, sphi_a, sphi_b, zeta, &
                                                zetb
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: ap_iterator, nl_iterator

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//"_forces",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    IF (calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-2.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    END IF

    maxder = ncoset(nder)

    CALL get_qs_kind_set(qs_kind_set,&
            maxco=maxco,maxlgto=maxl,maxsgf=maxsgf,maxnset=maxnset)

    CALL init_orbital_pointers(maxl+nder+1)

    ldsab = MAX(maxco,maxsgf)
    ldai = ncoset(maxl+nder+1)
    ALLOCATE(hab(ldsab,ldsab,maxnset*maxnset),work(ldsab,ldsab),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (verf(ldai,ldai,2*maxl+nder+1),vnuc(ldai,ldai,2*maxl+nder+1),ff(0:2*maxl+nder),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (calculate_forces) THEN
       ALLOCATE(pab(maxco,maxco,maxnset*maxnset),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! iterator for basis/potential list
    CALL neighbor_list_iterator_create(ap_iterator,sac_ae,search=.TRUE.)

    ALLOCATE (basis_set_list(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      CALL get_qs_kind(qs_kind_set(ikind),orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       atom_a = atom_of_kind(iatom)
       atom_b = atom_of_kind(jatom)
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       IF(inode==1) last_jatom = 0
       dab = SQRT(SUM(rab*rab))

       IF (jatom /= last_jatom) THEN
          new_atom_b = .TRUE.
          last_jatom = jatom
       ELSE
          new_atom_b = .FALSE.
       END IF

       ! *** Use the symmetry of the first derivatives ***
       IF (iatom == jatom) THEN
          f0 = 1.0_dp
       ELSE
          f0 = 2.0_dp
       END IF

       ! *** Create matrix blocks for a new matrix block column ***
       IF (new_atom_b) THEN
          IF (iatom <= jatom) THEN
             irow = iatom
             icol = jatom
          ELSE
             irow = jatom
             icol = iatom
          END IF
          NULLIFY(h_block)
          CALL cp_dbcsr_get_block_p(matrix=matrix_h(1)%matrix,&
               row=irow,col=icol,BLOCK=h_block,found=found)
          IF (calculate_forces) THEN
             NULLIFY(p_block)
             CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
                  row=irow,col=icol,BLOCK=p_block,found=found)
             CPPostcondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)
             ! *** Decontract density matrix block ***
             DO iset=1,nseta
                ncoa = npgfa(iset)*ncoset(la_max(iset))
                sgfa = first_sgfa(1,iset)
                DO jset=1,nsetb
                   ncob = npgfb(jset)*ncoset(lb_max(jset))
                   sgfb = first_sgfb(1,jset)
                   nij=jset+(iset-1)*maxnset
                   IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              p_block(sgfa,sgfb),SIZE(p_block,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                   ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              p_block(sgfb,sgfa),SIZE(p_block,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                   END IF
                   CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                           1.0_dp,work(1,1),SIZE(work,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,pab(1,1,nij),SIZE(pab,1))
                END DO
             END DO
          END IF

       END IF

       ! loop over all kinds for pseudopotential  atoms
       hab = 0._dp
       DO kkind=1,nkind
          CALL get_qs_kind(qs_kind_set(kkind),all_potential=all_potential)
          IF (.NOT.ASSOCIATED(all_potential)) CYCLE
          CALL get_potential(potential=all_potential,&
                             alpha_core_charge=alpha_c,zeff=zeta_c,&
                             ccore_charge=core_charge,core_charge_radius=core_radius)

          CALL nl_set_sub_iterator(ap_iterator,ikind,kkind,iatom)
          DO WHILE (nl_sub_iterate(ap_iterator)==0)
             CALL get_iterator_info(ap_iterator,jatom=katom,r=rac)
             dac = SQRT(SUM(rac*rac))
             rbc(:) = rac(:) - rab(:)
             dbc = SQRT(SUM(rbc*rbc))
             IF ( (MAXVAL(set_radius_a(:)) + core_radius < dac)  .OR. &
                  (MAXVAL(set_radius_b(:)) + core_radius < dbc) ) THEN
               CYCLE
             END IF

             DO iset=1,nseta
                IF (set_radius_a(iset) + core_radius < dac) CYCLE
                ncoa = npgfa(iset)*ncoset(la_max(iset))
                sgfa = first_sgfa(1,iset)
                DO jset=1,nsetb
                   IF (set_radius_b(jset) + core_radius < dbc) CYCLE
                   ncob = npgfb(jset)*ncoset(lb_max(jset))
                   sgfb = first_sgfb(1,jset)
                   IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                   rab2 = dab*dab
                   rac2 = dac*dac
                   rbc2 = dbc*dbc
                   nij=jset+(iset-1)*maxnset
                   ! *** Calculate the GTH pseudo potential forces ***
                   IF (calculate_forces) THEN
                      na_plus = npgfa(iset)*ncoset(la_max(iset)+nder)
                      nb_plus = npgfb(jset)*ncoset(lb_max(jset))
                      ALLOCATE(habd(na_plus,nb_plus),STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      habd = 0._dp
                      CALL verfc(&
                              la_max(iset)+nder,npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              alpha_c,core_radius,zeta_c,core_charge,&
                              rab,rab2,rac,rac2,rbc2,hab(:,:,nij),verf,vnuc,ff(0:),&
                              nder,habd)

                      ! *** The derivatives w.r.t. atomic center c are    ***
                      ! *** calculated using the translational invariance ***
                      ! *** of the first derivatives                      ***
                      CALL verfc_force(habd,pab(:,:,nij),force_a,force_b,nder,&
                                      la_max(iset),la_min(iset),npgfa(iset),zeta(:,iset),&
                                      lb_max(jset),lb_min(jset),npgfb(jset),zetb(:,jset),rab,error)

                      DEALLOCATE(habd,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                      atom_c = atom_of_kind(katom)
                      force(ikind)%all_potential(1,atom_a) =force(ikind)%all_potential(1,atom_a) + f0*force_a(1)
                      force(ikind)%all_potential(2,atom_a) =force(ikind)%all_potential(2,atom_a) + f0*force_a(2)
                      force(ikind)%all_potential(3,atom_a) =force(ikind)%all_potential(3,atom_a) + f0*force_a(3)

                      force(jkind)%all_potential(1,atom_b) =force(jkind)%all_potential(1,atom_b) + f0*force_b(1)
                      force(jkind)%all_potential(2,atom_b) =force(jkind)%all_potential(2,atom_b) + f0*force_b(2)
                      force(jkind)%all_potential(3,atom_b) =force(jkind)%all_potential(3,atom_b) + f0*force_b(3)

                      force(kkind)%all_potential(1,atom_c) =force(kkind)%all_potential(1,atom_c) &
                                                            - f0*force_a(1) - f0*force_b(1)
                      force(kkind)%all_potential(2,atom_c) =force(kkind)%all_potential(2,atom_c) &
                                                            - f0*force_a(2) - f0*force_b(2)
                      force(kkind)%all_potential(3,atom_c) =force(kkind)%all_potential(3,atom_c) &
                                                            - f0*force_a(3) - f0*force_b(3)

                      IF (use_virial) THEN
                         CALL virial_pair_force ( virial%pv_virial, f0, force_a, rac, error)
                         CALL virial_pair_force ( virial%pv_virial, f0, force_b, rbc, error)
                      END IF
                   ELSE
                      CALL verfc(&
                              la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              alpha_c,core_radius,zeta_c,core_charge,&
                              rab,rab2,rac,rac2,rbc2,hab(:,:,nij),verf,vnuc,ff(0:))
                   END IF
                END DO
             END DO
!            sac_ae_neighbor_node => next(sac_ae_neighbor_node)
          END DO
       END DO
       ! *** Contract nuclear attraction integrals
       DO iset=1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          DO jset=1,nsetb
             ncob = npgfb(jset)*ncoset(lb_max(jset))
             sgfb = first_sgfb(1,jset)
             nij=jset+(iset-1)*maxnset
             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                  1.0_dp,hab(1,1,nij),SIZE(hab,1),&
                  sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  0.0_dp,work(1,1),SIZE(work,1))
             IF (iatom <= jatom) THEN
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     work(1,1),SIZE(work,1),&
                     1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
             ELSE
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                     1.0_dp,work(1,1),SIZE(work,1),&
                     sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
             END IF
          END DO
       END DO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    CALL neighbor_list_iterator_release(ap_iterator)

    DEALLOCATE (atom_of_kind,basis_set_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(hab,work,verf,vnuc,ff,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (calculate_forces) THEN
       DEALLOCATE(pab,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (calculate_forces) THEN
       ! *** If LSD, then recover alpha density and beta density     ***
       ! *** from the total density (1) and the spin density (2)     ***
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 0.5_dp, beta_scalar=0.5_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-1.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_core_ae

!==========================================================================================================
! *****************************************************************************
!> \brief ...
!> \param habd ...
!> \param pab ...
!> \param fa ...
!> \param fb ...
!> \param nder ...
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param rab ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE verfc_force(habd,pab,fa,fb,nder,la_max,la_min,npgfa,zeta,lb_max,lb_min,npgfb,zetb,rab,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: habd, pab
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: fa, fb
    INTEGER, INTENT(IN)                      :: nder, la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'verfc_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: ic_a, ic_b, icam1, icam2, icam3, icap1, icap2, icap3, icax, &
      icbm1, icbm2, icbm3, icbx, icoa, icob, ipgfa, ipgfb, na, nap, nb
    INTEGER, DIMENSION(3)                    :: la, lb
    REAL(KIND=dp)                            :: zax2, zbx2

    fa = 0.0_dp
    fb = 0.0_dp

    na = ncoset(la_max)
    nap = ncoset(la_max+nder)
    nb = ncoset(lb_max)
    DO ipgfa = 1,npgfa
       zax2 = zeta(ipgfa)*2.0_dp
       DO ipgfb = 1,npgfb
          zbx2 = zetb(ipgfb)*2.0_dp
          DO ic_a = ncoset(la_min-1)+1,ncoset(la_max)
             la(1:3) = indco(1:3,ic_a)
             icap1 = coset(la(1)+1,la(2),la(3))
             icap2 = coset(la(1),la(2)+1,la(3))
             icap3 = coset(la(1),la(2),la(3)+1)
             icam1 = coset(la(1)-1,la(2),la(3))
             icam2 = coset(la(1),la(2)-1,la(3))
             icam3 = coset(la(1),la(2),la(3)-1)
             icoa =  ic_a + (ipgfa-1)*na
             icax =  (ipgfa-1)*nap

             DO ic_b = ncoset(lb_min-1)+1,ncoset(lb_max)
                lb(1:3) = indco(1:3,ic_b)
                icbm1 = coset(lb(1)-1,lb(2),lb(3))
                icbm2 = coset(lb(1),lb(2)-1,lb(3))
                icbm3 = coset(lb(1),lb(2),lb(3)-1)
                icob = ic_b + (ipgfb-1)*nb
                icbx = (ipgfb-1)*nb

                fa(1) = fa(1) - pab(icoa,icob)*(-zax2*habd(icap1+icax,icob)+&
                                                REAL(la(1),KIND=dp)*habd(icam1+icax,icob))
                fa(2) = fa(2) - pab(icoa,icob)*(-zax2*habd(icap2+icax,icob)+&
                                                REAL(la(2),KIND=dp)*habd(icam2+icax,icob))
                fa(3) = fa(3) - pab(icoa,icob)*(-zax2*habd(icap3+icax,icob)+&
                                                REAL(la(3),KIND=dp)*habd(icam3+icax,icob))

                fb(1) = fb(1) - pab(icoa,icob)*(-zbx2*(habd(icap1+icax,icob)-&
                                rab(1)*habd(ic_a+icax,icob))+REAL(lb(1),KIND=dp)*habd(ic_a+icax,icbm1+icbx))
                fb(2) = fb(2) - pab(icoa,icob)*(-zbx2*(habd(icap2+icax,icob)-&
                                rab(2)*habd(ic_a+icax,icob))+REAL(lb(2),KIND=dp)*habd(ic_a+icax,icbm2+icbx))
                fb(3) = fb(3) - pab(icoa,icob)*(-zbx2*(habd(icap3+icax,icob)-&
                                rab(3)*habd(ic_a+icax,icob))+REAL(lb(3),KIND=dp)*habd(ic_a+icax,icbm3+icbx))

             END DO  ! ic_b
          END DO  ! ic_a
       END DO  ! ipgfb
    END DO  ! ipgfa

  END SUBROUTINE verfc_force

!==========================================================================================================

END MODULE core_ae
