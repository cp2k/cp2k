!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE dbcsr_types,                     ONLY: dbcsr_1d_array_obj
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       find_neighbor_list, first_list, first_node, get_neighbor_list, &
       get_neighbor_list_set, get_neighbor_node, neighbor_list_set_p_type, &
       neighbor_list_type, neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             get_block_node,&
                                             real_matrix_p_type
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

    !TYPE(dbcsr_1d_array_obj)                   :: matrix_h
    TYPE(real_matrix_p_type), DIMENSION(:), &
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

    INTEGER :: atom_a, atom_b, atom_c, handle, iab, iac, iatom, icol, ikind, &
      ilist, inode, irow, iset, jatom, jkind, jset, katom, kkind, kneighbor, &
      last_jatom, ldai, ldsab, maxco, maxder, maxl, maxlgto, maxlppl, &
      maxnset, maxsgf, natom, ncoa, ncob, nkind, nlist, nneighbor, nnode, &
      nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, new_atom_b
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
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sac_ppl_neighbor_node

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//" (forces)",handle)
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
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,matrix_p(2)%matrix,1.0_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-2.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
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

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                              first_sgf=first_sgfa,&
                              lmax=la_max,&
                              lmin=la_min,&
                              npgf=npgfa,&
                              nset=nseta,&
                              nsgf_set=nsgfa,&
                              pgf_radius=rpgfa,&
                              set_radius=set_radius_a,&
                              sphi=sphi_a,&
                              zet=zeta)

       DO jkind=1,nkind
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 first_sgf=first_sgfb,&
                                 lmax=lb_max,&
                                 lmin=lb_min,&
                                 npgf=npgfb,&
                                 nset=nsetb,&
                                 nsgf_set=nsgfb,&
                                 pgf_radius=rpgfb,&
                                 set_radius=set_radius_b,&
                                 sphi=sphi_b,&
                                 zet=zetb)

          iab = ikind + nkind*(jkind - 1)
          IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
          CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)

          NULLIFY ( sab_orb_neighbor_list )
          DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,atom=iatom,nnode=nnode)
             atom_a = atom_of_kind(iatom)

             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

             DO inode=1,nnode
                CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)
                dab = SQRT(SUM(rab*rab))
                atom_b = atom_of_kind(jatom)

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
                   !CALL dbcsr_get_block_node(matrix=matrix_h%s(1),&
                   CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                             block_row=irow,&
                                             block_col=icol,&
                                             BLOCK=h_block)

                   CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)
                   IF (calculate_forces) THEN
                      CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=p_block)
                      CPPostcondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)
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

                   iac= ikind + nkind*(kkind - 1)
                   IF (.NOT.ASSOCIATED(sac_ppl(iac)%neighbor_list_set)) CYCLE
                   sac_ppl_neighbor_list => find_neighbor_list(sac_ppl(iac)%neighbor_list_set,atom=iatom)

                   CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,nnode=nneighbor)

                   sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)

                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,neighbor=katom,r=rac)

                      dac = SQRT(SUM(rac*rac))
                      rbc(:) = rac(:) - rab(:)
                      dbc = SQRT(SUM(rbc*rbc))
                      IF ( (MAXVAL(set_radius_a(:)) + ppl_radius < dac)  .OR. &
                           (MAXVAL(set_radius_b(:)) + ppl_radius < dbc) ) THEN
                        sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
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
                      sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
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

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO

          END DO
       END DO
    END DO

    DEALLOCATE (atom_of_kind,STAT=stat)
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
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 0.5_dp,matrix_p(2)%matrix,0.5_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-1.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
       END IF
    END IF

    CALL timestop(handle)
    
  END SUBROUTINE build_core_ppl

!==========================================================================================================
  
END MODULE core_ppl
