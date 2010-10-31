!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Calculation of the local pseudopotential contribution to the core Hamiltonian 
!>         <a|V(local)|b> = <a|Sum e^a*rc**2|b>
!> \par History
!>      - core_ppnl refactored from qs_core_hamiltonian [Joost VandeVondele, 2008-11-01]
!>      - adapted for PPL [jhu, 2009-02-23]
! *****************************************************************************
MODULE core_ppl

  USE ai_overlap_ppl,                  ONLY: overlap_ppl
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_force_types,                  ONLY: qs_force_type
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
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppl'

  PUBLIC :: build_core_ppl

CONTAINS

!==========================================================================================================

  SUBROUTINE build_core_ppl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    atomic_kind_set, particle_set, sab_orb, sac_ppl, error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL                                  :: use_virial
    INTEGER                                  :: nder
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_ppl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppl', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, handle, iatom, icol, ikind, inode, &
      irow, iset, jatom, jkind, jset, katom, kkind, last_jatom, ldai, ldsab, &
      maxco, maxder, maxl, maxlgto, maxlppl, maxnset, maxsgf, natom, ncoa, &
      ncob, nkind, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b
    REAL(KIND=dp)                            :: alpha_ppl, dab, dac, dbc, f0, &
                                                ppl_radius
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ai_work, hab, pab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, rac, &
                                                rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cexp_ppl, set_radius_a, &
                                                set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, rpgfa, &
                                                rpgfb, sphi_a, sphi_b, zeta, &
                                                zetb
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
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

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            maxco=maxco,maxlgto=maxlgto,maxlppl=maxlppl,maxsgf=maxsgf,maxnset=maxnset)

    maxl = MAX(maxlgto,maxlppl)
    CALL init_orbital_pointers(maxl+nder+1)

    ldsab = MAX(maxco,ncoset(maxlppl),maxsgf,maxlppl)
    ldai = ncoset(maxl+nder+1)
    ALLOCATE(hab(ldsab,ldsab,maxnset,maxnset),work(ldsab,ldsab*maxder),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ai_work(ldai,ldai,MAX(1,ncoset(maxlppl)),ncoset(nder+1)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (calculate_forces) THEN
       ALLOCATE(pab(maxco,maxco,maxnset,maxnset),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! iterator for basis/potential list
    CALL neighbor_list_iterator_create(ap_iterator,sac_ppl,search=.TRUE.)

    ALLOCATE (basis_set_list(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
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
          CALL cp_dbcsr_get_block_p(matrix_h(1)%matrix,irow,icol,h_block,found)
          !CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)
          IF(ASSOCIATED(h_block)) THEN
          IF (calculate_forces) THEN
             NULLIFY(p_block)
             CALL cp_dbcsr_get_block_p(matrix_p(1)%matrix,irow,icol,p_block,found)
             IF(ASSOCIATED(p_block)) THEN
                DO iset=1,nseta
                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   DO jset=1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      ! *** Decontract density matrix block ***
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
                           0.0_dp,pab(1,1,iset,jset),SIZE(pab,1))
                   END DO
                END DO
             ENDIF
          END IF
          END IF
       END IF

       hab = 0._dp

       ! loop over all kinds for pseudopotential  atoms
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,gth_potential=gth_potential)
          IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alpha_ppl,cexp_ppl=cexp_ppl,ppl_radius=ppl_radius)
          CALL nl_set_sub_iterator(ap_iterator,ikind,kkind,iatom)
          DO WHILE (nl_sub_iterate(ap_iterator)==0)
             CALL get_iterator_info(ap_iterator,jatom=katom,r=rac)

             dac = SQRT(SUM(rac*rac))
             rbc(:) = rac(:) - rab(:)
             dbc = SQRT(SUM(rbc*rbc))
             IF ( (MAXVAL(set_radius_a(:)) + ppl_radius < dac)  .OR. &
                  (MAXVAL(set_radius_b(:)) + ppl_radius < dbc) ) THEN
               CYCLE
             END IF

             DO iset=1,nseta
                IF (set_radius_a(iset) + ppl_radius < dac) CYCLE
                ncoa = npgfa(iset)*ncoset(la_max(iset))
                sgfa = first_sgfa(1,iset)
                DO jset=1,nsetb
                   IF (set_radius_b(jset) + ppl_radius < dbc) CYCLE
                   ncob = npgfb(jset)*ncoset(lb_max(jset))
                   sgfb = first_sgfb(1,jset)
                   IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                   ! *** Calculate the GTH pseudo potential forces ***
                   IF (calculate_forces) THEN

                      CALL overlap_ppl(&
                              la_max(iset),la_min(iset),npgfa(iset),&
                              rpgfa(:,iset),zeta(:,iset),&
                              lb_max(jset),lb_min(jset),npgfb(jset),&
                              rpgfb(:,jset),zetb(:,jset),&
                              cexp_ppl,alpha_ppl,ppl_radius,&
                              rab,dab,rac,dac,rbc,dbc,&
                              hab(:,:,iset,jset),nder,nder,.FALSE.,ai_work,&
                              pab(:,:,iset,jset),force_a,force_b)

                      ! *** The derivatives w.r.t. atomic center c are    ***
                      ! *** calculated using the translational invariance ***
                      ! *** of the first derivatives                      ***
                      atom_c = atom_of_kind(katom)
                      force(ikind)%gth_ppl(1,atom_a) =force(ikind)%gth_ppl(1,atom_a) + f0*force_a(1)
                      force(ikind)%gth_ppl(2,atom_a) =force(ikind)%gth_ppl(2,atom_a) + f0*force_a(2)
                      force(ikind)%gth_ppl(3,atom_a) =force(ikind)%gth_ppl(3,atom_a) + f0*force_a(3)
                      force(kkind)%gth_ppl(1,atom_c) =force(kkind)%gth_ppl(1,atom_c) - f0*force_a(1)
                      force(kkind)%gth_ppl(2,atom_c) =force(kkind)%gth_ppl(2,atom_c) - f0*force_a(2)
                      force(kkind)%gth_ppl(3,atom_c) =force(kkind)%gth_ppl(3,atom_c) - f0*force_a(3)

                      force(jkind)%gth_ppl(1,atom_b) =force(jkind)%gth_ppl(1,atom_b) + f0*force_b(1)
                      force(jkind)%gth_ppl(2,atom_b) =force(jkind)%gth_ppl(2,atom_b) + f0*force_b(2)
                      force(jkind)%gth_ppl(3,atom_b) =force(jkind)%gth_ppl(3,atom_b) + f0*force_b(3)
                      force(kkind)%gth_ppl(1,atom_c) =force(kkind)%gth_ppl(1,atom_c) - f0*force_b(1)
                      force(kkind)%gth_ppl(2,atom_c) =force(kkind)%gth_ppl(2,atom_c) - f0*force_b(2)
                      force(kkind)%gth_ppl(3,atom_c) =force(kkind)%gth_ppl(3,atom_c) - f0*force_b(3)

                      IF (use_virial) THEN
                         CALL virial_pair_force ( virial%pv_virial, f0, force_a, rac, error)
                         CALL virial_pair_force ( virial%pv_virial, f0, force_b, rbc, error)
                      END IF
                   ELSE
                      CALL overlap_ppl(&
                           la_max(iset),la_min(iset),npgfa(iset),&
                           rpgfa(:,iset),zeta(:,iset),&
                           lb_max(jset),lb_min(jset),npgfb(jset),&
                           rpgfb(:,jset),zetb(:,jset),&
                           cexp_ppl,alpha_ppl,ppl_radius,&
                           rab,dab,rac,dac,rbc,dbc,&
                           hab(:,:,iset,jset),0,0,.FALSE.,ai_work)
                   END IF
                END DO
             END DO
          END DO
       END DO

       ! *** Contract PPL integrals
       DO iset=1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          DO jset=1,nsetb
             ncob = npgfb(jset)*ncoset(lb_max(jset))
             sgfb = first_sgfb(1,jset)
             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                  1.0_dp,hab(1,1,iset,jset),SIZE(hab,1),&
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
    DEALLOCATE(hab,work,ai_work,STAT=stat)
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
    
  END SUBROUTINE build_core_ppl

!==========================================================================================================
  
END MODULE core_ppl
