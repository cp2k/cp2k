!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Calculation of the non-local pseudopotential contribution to the core Hamiltonian 
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par History
!>      - refactered from qs_core_hamiltian [Joost VandeVondele, 2008-11-01]
! *****************************************************************************
MODULE core_ppnl

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       find_neighbor_list, first_list, first_node, get_neighbor_list, &
       get_neighbor_list_set, get_neighbor_node, neighbor_list_set_p_type, &
       neighbor_list_type, neighbor_node_type, next, qlist_type
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             get_block_node,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppnl'

  TYPE vppnl_type
     TYPE(qlist_type), DIMENSION(:,:,:), POINTER :: neighbor
      REAL(KIND = dp)                            :: ppnl_radius
      REAL(KIND = dp), DIMENSION(:), POINTER     :: alpha_ppnl
      REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cprj,vprj_ppnl
      INTEGER                                    :: lppnl,nppnl
      INTEGER, DIMENSION(:), POINTER             :: nprj_ppnl
   END TYPE vppnl_type

  PUBLIC :: build_core_ppnl

CONTAINS

  SUBROUTINE build_core_ppnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    atomic_kind_set, particle_set, sab_orb, sac_ppnl, sbc_ppnl, error)


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
      DIMENSION(:), POINTER                  :: sab_orb, sac_ppnl, sbc_ppnl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppnl', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, first_col, handle, i, iab, iac, iatom, &
      ibc, icell, ico, icol, idx, ikind, ilist, inode, irow, iset, j, jatom, &
      jcell, jco, jkind, jset, katom, kcell, kkind, kneighbor, knode, l, &
      last_jatom, lc_max, lc_min, ldai, ldsab, lppnl, maxco, maxder, maxl, &
      maxlgto, maxlppnl, maxppnl, maxsgf, n, natom, ncoa, ncoa_sum, ncob, &
      ncoc, ncol, nkind, nlist, nneighbor, nnode, nnode_c, nppnl, nprjc, &
      nrow, nseta, nsetb, nsgf, prjc, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cell_c, cell_c_max, cell_c_min
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, new_atom_b, &
                                                ppnl_present
    REAL(KIND=dp)                            :: dab, dac, dbc, f, f0, rab2, &
                                                rac2, rbc2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: hab, pab, sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(1)              :: rprjc, zetc
    REAL(KIND=dp), DIMENSION(3)              :: ffa, ffb, rab, rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block, rpgfa, &
                                                rpgfb, sphi_a, sphi_b, zeta, &
                                                zetb
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list, &
      sab_orb_neighbor_list_local, sac_ppnl_neighbor_list, &
      sbc_ppnl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sac_ppnl_neighbor_node, &
                                                sbc_ppnl_neighbor_node
    TYPE(qlist_type), POINTER                :: sac_ppnl_neighbor
    TYPE(vppnl_type), ALLOCATABLE, &
      DIMENSION(:)                           :: vppnl

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//" (forces)",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF


    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    
    IF (calculate_forces) THEN
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,matrix_p(2)%matrix,1.0_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-2.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
       END IF
    END IF
    
    maxder = ncoset(nder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)
    
   ! *** Allocate work storage ***
   CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                            maxco=maxco,&
                            maxlgto=maxlgto,&
                            maxlppnl=maxlppnl,&
                            maxppnl=maxppnl,&
                            maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppnl)
    CALL init_orbital_pointers(maxl+nder+1)

    ldsab = MAX(maxco,ncoset(maxlppnl),maxsgf,maxppnl)
    ldai = ncoset(maxl+nder+1)
    
    ALLOCATE(hab(ldsab,ldsab*maxder),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(ldsab,ldsab*maxder),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(sab(ldsab,ldsab*maxder),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    hab=0.0_dp ; work=0.0_dp ; sab=0.0_dp
    IF (calculate_forces) THEN
       ALLOCATE(pab(maxco,maxco),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       pab=0.0_dp
    END IF
    
    ! *** Load GTH pseudo potential data (non-local part -> PPNL) ***
    ppnl_present = ASSOCIATED(sac_ppnl)
    IF (ppnl_present) THEN
       ALLOCATE (vppnl(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               gth_potential=gth_potential)

             NULLIFY (vppnl(kkind)%alpha_ppnl)
             NULLIFY (vppnl(kkind)%cprj)
             NULLIFY (vppnl(kkind)%neighbor)
             NULLIFY (vppnl(kkind)%nprj_ppnl)
             NULLIFY (vppnl(kkind)%vprj_ppnl)
             IF (ASSOCIATED(gth_potential)) THEN
                CALL get_potential(potential=gth_potential,&
                                   alpha_ppnl=vppnl(kkind)%alpha_ppnl,&
                                   cprj=vppnl(kkind)%cprj,&
                                   lppnl=vppnl(kkind)%lppnl,&
                                   nppnl=vppnl(kkind)%nppnl,&
                                   nprj_ppnl=vppnl(kkind)%nprj_ppnl,&
                                   ppnl_radius=vppnl(kkind)%ppnl_radius,&
                                   vprj_ppnl=vppnl(kkind)%vprj_ppnl)
             END IF
       END DO
    END IF
   
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
       
       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                              first_sgf=first_sgfa,&
                              lmax=la_max,&
                              lmin=la_min,&
                              nco_sum=ncoa_sum,&
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
          ALLOCATE (ai_work(ldai,ldai,ncoset(nder+1)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,atom=iatom,nnode=nnode)

             IF (calculate_forces) atom_a = atom_of_kind(iatom)

             IF (ppnl_present) THEN
                DO kkind=1,nkind
                   iac = ikind + nkind*(kkind - 1)
                   IF (.NOT.ASSOCIATED(sac_ppnl(iac)%neighbor_list_set)) CYCLE
                   lppnl = vppnl(kkind)%lppnl
                   nppnl = vppnl(kkind)%nppnl
                   rprjc(1) = vppnl(kkind)%ppnl_radius
                   sac_ppnl_neighbor_list => find_neighbor_list(sac_ppnl(iac)%neighbor_list_set,atom=iatom)
                   
                   CALL get_neighbor_list(neighbor_list=sac_ppnl_neighbor_list,nnode=nneighbor)
                   
                   ! *** Find the proper cell index ranges ***
                   cell_c_max(:) = 0
                   cell_c_min(:) = 0
                   
                   sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)
                   
                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,cell=cell_c)
                      DO i=1,3
                         cell_c_max(i) = MAX(cell_c_max(i),cell_c(i))
                         cell_c_min(i) = MIN(cell_c_min(i),cell_c(i))
                      END DO
                      sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
                   END DO
                   
                   ! *** Allocate and initialize the sac_ppnl neighbor lists ***
                   ALLOCATE (vppnl(kkind)%neighbor(cell_c_min(1):cell_c_max(1),&
                                                   cell_c_min(2):cell_c_max(2),&
                                                   cell_c_min(3):cell_c_max(3)),&
                                                   STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   vppnl(kkind)%neighbor(:,:,:)%n = 0
                   
                   ! *** Find proper array size for each cell ***
                   sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)
                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,cell=cell_c)
                      vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n =&
                           vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n + 1
                      sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
                   END DO
                   
                   ! *** Allocate the sac_ppnl neighbor list ***
                   DO kcell=cell_c_min(3),cell_c_max(3)
                      DO jcell=cell_c_min(2),cell_c_max(2)
                         DO icell=cell_c_min(1),cell_c_max(1)
                            sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,jcell,kcell)
                            NULLIFY (sac_ppnl_neighbor%sac)
                            NULLIFY (sac_ppnl_neighbor%r2)
                            NULLIFY (sac_ppnl_neighbor%r)
                            NULLIFY (sac_ppnl_neighbor%index_list)
                            NULLIFY (sac_ppnl_neighbor%list)
                            n = sac_ppnl_neighbor%n
                            IF (n > 0) THEN
                               CALL reallocate(sac_ppnl_neighbor%r2,1,n)
                               CALL reallocate(sac_ppnl_neighbor%r,1,3,1,n)
                               CALL reallocate(sac_ppnl_neighbor%index_list,1,n)
                               CALL reallocate(sac_ppnl_neighbor%list,1,n)
                               ALLOCATE (sac_ppnl_neighbor%sac(maxder),STAT=stat)
                               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                               DO i=1,maxder
                                  ncol = n*nppnl
                                  ALLOCATE (sac_ppnl_neighbor%sac(i)%block(ncoa_sum,ncol),STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                               END DO
                            END IF
                         END DO
                      END DO
                   END DO
                   
                   ! *** Fill sac_ppnl neighbor lists ***
                   vppnl(kkind)%neighbor(:,:,:)%n = 0
                   sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)
                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,&
                           neighbor=katom,cell=cell_c,r=rac)
                      sac_ppnl_neighbor => vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))
                      sac_ppnl_neighbor%n = sac_ppnl_neighbor%n + 1
                      sac_ppnl_neighbor%list(sac_ppnl_neighbor%n) = katom
                      rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                      sac_ppnl_neighbor%r2(sac_ppnl_neighbor%n) = rac2
                      sac_ppnl_neighbor%r(:,sac_ppnl_neighbor%n) = rac(:)
                      sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
                   END DO
                   
                   ! *** Build sac overlap integral blocks ***
                   DO kcell=cell_c_min(3),cell_c_max(3)
                      DO jcell=cell_c_min(2),cell_c_max(2)
                         DO icell=cell_c_min(1),cell_c_max(1)
                            
                            sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,jcell,kcell)
                            ! *** Sort sac_ppnl neighbor lists ***
                            IF (sac_ppnl_neighbor%n > 0) THEN
                               CALL sort(sac_ppnl_neighbor%list,sac_ppnl_neighbor%n,sac_ppnl_neighbor%index_list)
                            ELSE
                               CYCLE
                            END IF
                            ncol = 1
                            DO kneighbor=1,sac_ppnl_neighbor%n
                               idx = sac_ppnl_neighbor%index_list(kneighbor)
                               rac(:) = sac_ppnl_neighbor%r(:,idx)
                               rac2 = sac_ppnl_neighbor%r2(idx)
                               dac = SQRT(rac2)
                               nrow = 1
                               DO iset=1,nseta
                                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                                  prjc = 1
                                  DO l=0,lppnl
                                     nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)
                                     IF (nprjc == 0) CYCLE
                                     IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
                                     lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                                     lc_min = l
                                     zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                                     ncoc = ncoset(lc_max)
                                     
                                     ! *** Calculate the primitive overlap integrals ***
                                     CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                          rpgfa(:,iset),zeta(:,iset),&
                                          lc_max,lc_min,1,rprjc,zetc,&
                                          rac,dac,sab,nder,.TRUE.,ai_work,ldai)
                                     
                                     ! *** Contraction step (projector functions) ***
                                     DO i=1,maxder
                                        first_col = (i - 1)*SIZE(work,1)
                                        CALL dgemm("N","N",ncoa,nprjc,ncoc,&
                                             1.0_dp,sab(1,first_col+1),SIZE(sab,1),&
                                             vppnl(kkind)%cprj(1,prjc),&
                                             SIZE(vppnl(kkind)%cprj,1),0.0_dp,&
                                             work(1,first_col+prjc),SIZE(work,1))
                                     END DO
                                     prjc = prjc + nprjc
                                  END DO
                                  DO i=1,maxder
                                     first_col = (i - 1)*SIZE(work,1) + 1
                                     CALL dgemm("N","N",ncoa,nppnl,nppnl,&
                                          1.0_dp,work(1,first_col),SIZE(work,1),&
                                          vppnl(kkind)%vprj_ppnl(1,1),&
                                          SIZE(vppnl(kkind)%vprj_ppnl,1),0.0_dp,&
                                          sac_ppnl_neighbor%sac(i)%block(nrow,ncol),&
                                          SIZE(sac_ppnl_neighbor%sac(i)%block,1))
                                  END DO
                                  nrow = nrow + ncoa
                               END DO
                               ncol = ncol + nppnl
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
             END IF
             
             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)
             
             DO inode=1,nnode
                CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)
                
                IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                END IF
                
                IF (calculate_forces) atom_b = atom_of_kind(jatom)
                
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
                   CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                        block_row=irow,&
                                        block_col=icol,&
                                        BLOCK=h_block)
                   IF (calculate_forces) THEN
                      CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=p_block)
                   END IF
                END IF
                
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)
                nrow = 1
                DO iset=1,nseta
                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   DO jset=1,nsetb
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      DO jco=1,ncob
                         DO ico=1,ncoa
                            hab(ico,jco) = 0.0_dp
                         END DO
                      END DO
                      IF (calculate_forces) THEN
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
                              0.0_dp,pab(1,1),SIZE(pab,1))
                      END IF
                      ! *** Loop over the GTH pseudo potential atoms (non-local part) ***
                      IF (ppnl_present.AND.new_atom_b) THEN
                         DO kkind=1,nkind
                            IF (.NOT.ASSOCIATED(vppnl(kkind)%neighbor)) CYCLE
                            ibc = jkind + nkind*(kkind - 1)
                            IF (.NOT.ASSOCIATED(sbc_ppnl(ibc)%neighbor_list_set)) CYCLE
                            
                            lppnl = vppnl(kkind)%lppnl
                            nppnl = vppnl(kkind)%nppnl
                            rprjc(1) = vppnl(kkind)%ppnl_radius
                            
                            sbc_ppnl_neighbor_list => find_neighbor_list(sbc_ppnl(ibc)%neighbor_list_set,atom=jatom)
                            
                            CALL get_neighbor_list(neighbor_list=sbc_ppnl_neighbor_list,nnode=nnode_c)
                            
                            sbc_ppnl_neighbor_node => first_node(sbc_ppnl_neighbor_list)
                            
                            knodeloop: DO knode=1,nnode_c
                               CALL get_neighbor_node(sbc_ppnl_neighbor_node,neighbor=katom,cell=cell_c,r=rbc)
                               
                               ! The < and > operators don't work correctly on arrays in
                               ! the openmp version
                               DO i=1, 3
                                  IF ( (cell_c(i) < LBOUND(vppnl(kkind)%neighbor,i)) .OR. &
                                       (cell_c(i) > UBOUND(vppnl(kkind)%neighbor,i))) THEN
                                     sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                                     CYCLE knodeloop
                                  END IF
                               END DO
                               
                               sac_ppnl_neighbor => vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))
                               IF (sac_ppnl_neighbor%n == 0) THEN
                                  sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                                  CYCLE
                               END IF
                               
                               ! *** Locate operator atom in the sac_ppnl neighbor list ***
                               kneighbor = locate(sac_ppnl_neighbor%list,katom)
                               
                               IF (kneighbor == 0) THEN
                                  sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                                  CYCLE
                               END IF
                               
                               idx = sac_ppnl_neighbor%index_list(kneighbor)
                               rac(:) = sac_ppnl_neighbor%r(:,idx)
                               
                               rac2 = sac_ppnl_neighbor%r2(idx)
                               dac = SQRT(rac2)
                               IF (set_radius_a(iset) + rprjc(1) < dac) THEN
                                  sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                                  CYCLE
                               END IF
                               
                               rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                               dbc = SQRT(rbc2)
                               IF (set_radius_b(jset) + rprjc(1) < dbc) THEN
                                  sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                                  CYCLE
                               END IF
                               
                               prjc = 1
                               DO l=0,lppnl
                                  nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)
                                  IF (nprjc == 0) CYCLE
                                  lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                                  lc_min = l
                                  zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                                  ncoc = ncoset(lc_max)
                                  
                                  ! *** Calculate the primitive overlap integrals ***
                                  CALL overlap(lb_max(jset),lb_min(jset),npgfb(jset),&
                                       rpgfb(:,jset),zetb(:,jset),&
                                       lc_max,lc_min,1,rprjc,zetc,&
                                       rbc,dbc,sab,nder,.TRUE.,ai_work,ldai)
                                  
                                  ! *** Contraction step (projector functions) ***
                                  DO i=1,maxder
                                     first_col = (i - 1)*SIZE(work,1)
                                     CALL dgemm("N","N",ncob,nprjc,ncoc,1.0_dp,&
                                          sab(1,first_col+1),SIZE(sab,1),&
                                          vppnl(kkind)%cprj(1,prjc),&
                                          SIZE(vppnl(kkind)%cprj,1),&
                                          0.0_dp,work(1,first_col+prjc),&
                                          SIZE(work,1))
                                  END DO
                                  prjc = prjc + nprjc
                               END DO
                               ncol = (kneighbor - 1)*nppnl + 1
                               
                               CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                    sac_ppnl_neighbor%sac(1)%block(nrow,ncol),&
                                    SIZE(sac_ppnl_neighbor%sac(1)%block,1),&
                                    work(1,1),SIZE(work,1),&
                                    1.0_dp,hab(1,1),SIZE(hab,1))
                               
                               IF (calculate_forces) THEN
                                  atom_c = atom_of_kind(katom)
                                  DO i=1,3
                                     first_col = i*SIZE(work,1) + 1
                                     CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                          sac_ppnl_neighbor%sac(i+1)%block(nrow,ncol),&
                                          SIZE(sac_ppnl_neighbor%sac(i+1)%block,1),&
                                          work(1,1),SIZE(work,1),&
                                          0.0_dp,hab(1,first_col),SIZE(hab,1))
                                     f = 0.0_dp
                                     DO jco=1,ncob
                                        j = first_col + jco - 1
                                        DO ico=1,ncoa
                                           f = f + pab(ico,jco)*hab(ico,j)
                                        END DO
                                     END DO
                                     force(ikind)%gth_ppnl(i,atom_a) =force(ikind)%gth_ppnl(i,atom_a) + f0*f
                                     force(kkind)%gth_ppnl(i,atom_c) =force(kkind)%gth_ppnl(i,atom_c) - f0*f
                                     ffa(i) = f
                                     CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                          sac_ppnl_neighbor%sac(1)%block(nrow,ncol),&
                                          SIZE(sac_ppnl_neighbor%sac(1)%block,1),&
                                          work(1,first_col),SIZE(work,1),&
                                          0.0_dp,hab(1,first_col),SIZE(hab,1))
                                    
                                     f = 0.0_dp
                                     DO jco=1,ncob
                                        j = first_col + jco - 1
                                        DO ico=1,ncoa
                                           f = f + pab(ico,jco)*hab(ico,j)
                                        END DO
                                     END DO
                                     force(jkind)%gth_ppnl(i,atom_b) =force(jkind)%gth_ppnl(i,atom_b) + f0*f
                                     force(kkind)%gth_ppnl(i,atom_c) =force(kkind)%gth_ppnl(i,atom_c) - f0*f
                                     ffb(i) = f
                                  END DO
                                  IF (use_virial) THEN
                                     CALL virial_pair_force ( virial%pv_virial, f0, ffa, rac, error)
                                     CALL virial_pair_force ( virial%pv_virial, f0, ffb, rbc, error)
                                  END IF
                               END IF
                               sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                            END DO knodeloop
                         END DO
                      END IF
                      
                      ! *** Contraction step (core Hamiltonian matrix) ***
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
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
                   nrow = nrow + ncoa
                END DO
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             END DO
             IF (ppnl_present) THEN
                DO kkind=1,nkind
                   IF (ASSOCIATED(vppnl(kkind)%neighbor)) THEN
                      DO kcell=LBOUND(vppnl(kkind)%neighbor,3),&
                           UBOUND(vppnl(kkind)%neighbor,3)
                         DO jcell=LBOUND(vppnl(kkind)%neighbor,2),&
                              UBOUND(vppnl(kkind)%neighbor,2)
                            DO icell=LBOUND(vppnl(kkind)%neighbor,1),&
                                 UBOUND(vppnl(kkind)%neighbor,1)
                               sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,jcell,kcell)
                               IF (sac_ppnl_neighbor%n > 0) THEN
                                  DEALLOCATE (sac_ppnl_neighbor%r2,STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                                  DEALLOCATE (sac_ppnl_neighbor%r,STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                                  DEALLOCATE (sac_ppnl_neighbor%list,STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                                  DEALLOCATE (sac_ppnl_neighbor%index_list,STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                                  DO i=1,maxder
                                     DEALLOCATE (sac_ppnl_neighbor%sac(i)%block,STAT=stat)
                                     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                                  END DO
                                  DEALLOCATE (sac_ppnl_neighbor%sac,STAT=stat)
                                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                               END IF
                            END DO
                         END DO
                      END DO
                      DEALLOCATE (vppnl(kkind)%neighbor,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END IF
          END DO
          IF (ALLOCATED(ai_work)) THEN
             DEALLOCATE (ai_work,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO
    END DO

    IF (ppnl_present) THEN
       DEALLOCATE (vppnl,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    DEALLOCATE (hab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (pab,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! *** If LSD, then recover alpha density and beta density     ***
       ! *** from the total density (1) and the spin density (2)     ***
       ! *** The W matrix is neglected, since it will be destroyed   ***
       ! *** in the calling force routine after leaving this routine ***
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 0.5_dp,matrix_p(2)%matrix,0.5_dp,error=error)
          CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-1.0_dp,matrix_p(1)%matrix,1.0_dp,error=error)
       END IF
    END IF

    CALL timestop(handle)
    
  END SUBROUTINE build_core_ppnl
  
END MODULE core_ppnl
