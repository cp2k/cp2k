!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Calculation of the non-local pseudopotential contribution to the core Hamiltonian
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par History
!>      - refactered from qs_core_hamiltian [Joost VandeVondele, 2008-11-01]
!>      - full rewrite [jhu, 2009-01-23]
! *****************************************************************************
MODULE core_ppnl

  USE ai_overlap,                      ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE external_potential_types,        ONLY: gth_potential_p_type,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE sap_kind_types,                  ONLY: alist_type,&
                                             clist_type,&
                                             get_alist,&
                                             release_sap_int,&
                                             sap_int_type,&
                                             sap_sort
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppnl'

  PUBLIC :: build_core_ppnl

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
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param sab_orb ...
!> \param sap_ppnl ...
!> \param eps_ppnl ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_core_ppnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
                    atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, error)

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
      DIMENSION(:), POINTER                  :: sab_orb, sap_ppnl
    REAL(KIND=dp), INTENT(IN)                :: eps_ppnl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppnl', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, atom_c, first_col, handle, i, iab, iac, iatom, &
      ibc, icol, ikind, ilist, inode, irow, iset, j, jatom, jkind, jneighbor, &
      kac, katom, kbc, kkind, l, last_iatom, last_jatom, lc_max, lc_min, &
      ldai, ldsab, lppnl, maxco, maxder, maxl, maxlgto, maxlppnl, maxppnl, &
      maxsgf, mepos, na, natom, nb, ncoa, ncoc, nkind, nlist, nneighbor, &
      nnode, np, nppnl, nprjc, nseta, nsgfa, nthread, prjc, sgfa, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cell_b, cell_c
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, &
                                                nprj_ppnl, nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, found, &
                                                new_atom_pair, ppnl_present
    REAL(KIND=dp)                            :: dac, f0, ppnl_radius
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(1)              :: rprjc, zetc
    REAL(KIND=dp), DIMENSION(3)              :: fa, fb, rab, rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_ppnl, set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cprj, h_block, p_block, &
                                                rpgfa, sphi_a, vprj_ppnl, zeta
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: achint, acint, bchint, bcint
    TYPE(alist_type), POINTER                :: alist_ac, alist_bc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(clist_type), POINTER                :: clist
    TYPE(gth_potential_p_type), &
      DIMENSION(:), POINTER                  :: potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: sap_int

!$  INTEGER :: omp_get_thread_num, omp_get_max_threads

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//"_forces",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    ppnl_present = ASSOCIATED(sap_ppnl)

    IF ( ppnl_present ) THEN

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
                               maxco=maxco,&
                               maxlgto=maxlgto,&
                               maxlppnl=maxlppnl,&
                               maxppnl=maxppnl,&
                               maxsgf=maxsgf)

      maxl = MAX(maxlgto,maxlppnl)
      CALL init_orbital_pointers(maxl+nder+1)

      ldsab = MAX(maxco,ncoset(maxlppnl),maxsgf,maxppnl)
      ldai = ncoset(maxl+nder+1)

      !sap_int needs to be shared as multiple threads need to access this
      ALLOCATE(sap_int(nkind*nkind),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO i=1,nkind*nkind
        NULLIFY(sap_int(i)%alist,sap_int(i)%asort,sap_int(i)%aindex)
        sap_int(i)%nalist=0
      END DO

      !set up direct access to basis and potential
      ALLOCATE(basis_set(nkind),potential(nkind),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
           basis_set(ikind)%gto_basis_set => orb_basis_set
         ELSE
           NULLIFY(basis_set(ikind)%gto_basis_set)
         END IF
         CALL get_atomic_kind(atomic_kind=atomic_kind,gth_potential=gth_potential)
         IF (ASSOCIATED(gth_potential)) THEN
           potential(ikind)%gth_potential => gth_potential
         ELSE
           NULLIFY(potential(ikind)%gth_potential)
         END IF
      END DO

      nthread=1
!$    nthread=omp_get_max_threads()

      !calculate the overlap integrals <a|p>
      CALL neighbor_list_iterator_create(nl_iterator,sap_ppnl,nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &                                                                                      
!$OMP SHARED  (nl_iterator, basis_set, potential, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab,  ldai, nder, nco ) &
!$OMP PRIVATE (mepos, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          sphi_a, zeta, cprj, lppnl, nppnl, nprj_ppnl, & 
!$OMP          clist, iset, ncoa, sgfa, prjc, work, sab, ai_work, nprjc,  ppnl_radius, &
!$OMP          ncoc, rpgfa, first_col, vprj_ppnl, failure, i, l, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, zetc, alpha_ppnl, error, stat)
      mepos=0
!$    mepos=omp_get_thread_num()

      ALLOCATE(sab(ldsab,ldsab*maxder),work(ldsab,ldsab*maxder),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      sab=0.0_dp
      ALLOCATE (ai_work(ldai,ldai,ncoset(nder+1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ai_work=0.0_dp

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos)==0)
         CALL get_iterator_info(nl_iterator,mepos=mepos,ikind=ikind,jkind=kkind,iatom=iatom,&
              jatom=katom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,cell=cell_c,r=rac)
         iac = ikind + nkind*(kkind - 1)
         IF (.NOT.ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT.ASSOCIATED(potential(kkind)%gth_potential)) CYCLE
         ! get definition of basis set
         first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
         la_max => basis_set(ikind)%gto_basis_set%lmax
         la_min => basis_set(ikind)%gto_basis_set%lmin
         npgfa => basis_set(ikind)%gto_basis_set%npgf
         nseta = basis_set(ikind)%gto_basis_set%nset
         nsgfa = basis_set(ikind)%gto_basis_set%nsgf
         nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
         rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
         set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
         sphi_a => basis_set(ikind)%gto_basis_set%sphi
         zeta => basis_set(ikind)%gto_basis_set%zet
         ! get definition of PP projectors
         alpha_ppnl => potential(kkind)%gth_potential%alpha_ppnl
         cprj => potential(kkind)%gth_potential%cprj
         lppnl = potential(kkind)%gth_potential%lppnl
         nppnl = potential(kkind)%gth_potential%nppnl
         nprj_ppnl => potential(kkind)%gth_potential%nprj_ppnl
         ppnl_radius = potential(kkind)%gth_potential%ppnl_radius
         vprj_ppnl => potential(kkind)%gth_potential%vprj_ppnl
!$OMP CRITICAL(sap_int_critical)
         IF (.NOT.ASSOCIATED(sap_int(iac)%alist)) THEN
            sap_int(iac)%a_kind = ikind
            sap_int(iac)%p_kind = kkind
            sap_int(iac)%nalist = nlist
            ALLOCATE(sap_int(iac)%alist(nlist),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            DO i=1,nlist
               NULLIFY(sap_int(iac)%alist(i)%clist)
               sap_int(iac)%alist(i)%aatom = 0
               sap_int(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT.ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
            sap_int(iac)%alist(ilist)%aatom  = iatom
            sap_int(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE(sap_int(iac)%alist(ilist)%clist(nneighbor),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            DO i=1,nneighbor
               sap_int(iac)%alist(ilist)%clist(i)%catom = 0
            END DO
         END IF
!$OMP END CRITICAL(sap_int_critical)
         dac = SQRT(SUM(rac*rac))
         clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
         clist%catom = katom
         clist%cell  = cell_c
         clist%rac   = rac
         ALLOCATE(clist%acint(nsgfa,nppnl,maxder),&
                  clist%achint(nsgfa,nppnl,maxder),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         clist%acint=0._dp
         clist%achint=0._dp
         clist%nsgf_cnt = 0
         NULLIFY(clist%sgf_list)
         DO iset=1,nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1,iset)
            prjc = 1
           ! XXX fix, use correct bounds 
            work = 0._dp
            DO l=0,lppnl
               nprjc = nprj_ppnl(l)*nco(l)
               IF (nprjc == 0) CYCLE
               rprjc(1) = ppnl_radius
               IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE
               lc_max = l + 2*(nprj_ppnl(l) - 1)
               lc_min = l
               zetc(1) = alpha_ppnl(l)
               ncoc = ncoset(lc_max)

               ! *** Calculate the primitive overlap integrals ***
               CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                            lc_max,lc_min,1,rprjc,zetc,rac,dac,sab,nder,.TRUE.,ai_work,ldai)
               ! *** Transformation step projector functions (cartesian->spherical) ***
               DO i=1,maxder
                  first_col = (i - 1)*SIZE(work,1)
                  CALL dgemm("N","N",ncoa,nprjc,ncoc,1.0_dp,sab(1,first_col+1),SIZE(sab,1),&
                       cprj(1,prjc),SIZE(cprj,1),0.0_dp,work(1,first_col+prjc),SIZE(work,1))
               END DO
               prjc = prjc + nprjc
            END DO
            DO i=1,maxder
               first_col = (i - 1)*SIZE(work,1) + 1
               ! *** Contraction step (basis functions) ***
               CALL dgemm("T","N",nsgf_seta(iset),nppnl,ncoa,1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,first_col),SIZE(work,1),0.0_dp,clist%acint(sgfa,1,i),nsgfa)
               ! *** Multiply with interaction matrix(h) ***
               CALL dgemm("N","N",nsgf_seta(iset),nppnl,nppnl,1.0_dp,clist%acint(sgfa,1,i),nsgfa,&
                    vprj_ppnl(1,1),SIZE(vprj_ppnl,1),0.0_dp,clist%achint(sgfa,1,i),nsgfa)
            END DO
         END DO
         clist%maxac=MAXVAL(ABS(clist%acint(:,:,1)))
         clist%maxach=MAXVAL(ABS(clist%achint(:,:,1)))
      END DO

      DEALLOCATE(sab,ai_work,work,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      ! *** Set up a sorting index
      CALL sap_sort(sap_int,error)
      ! *** All integrals needed have been calculated and stored in sap_int
      ! *** We now calculate the Hamiltonian matrix elements
      CALL neighbor_list_iterator_create(nl_iterator,sab_orb,nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &                                                                                      
!$OMP SHARED  (nl_iterator, basis_set, atom_of_kind, matrix_h, &
!$OMP          matrix_p, sap_int, nkind, eps_ppnl, force, virial, use_virial, calculate_forces) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, atom_a, atom_b, last_iatom, last_jatom, f0, irow, icol, h_block, &
!$OMP          found, error, p_block, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, atom_c, j, fa, fb, rbc, rac, new_atom_pair, &
!$OMP          kkind, kac, kbc, i)

      mepos=0
!$    mepos=omp_get_thread_num()
      last_iatom = 0
      last_jatom = 0

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos)==0)
         CALL get_iterator_info(nl_iterator,mepos=mepos,ikind=ikind,jkind=jkind,iatom=iatom,&
              jatom=jatom,nlist=nlist,ilist=ilist,nnode=nnode,inode=inode,cell=cell_b,r=rab)
         IF (.NOT.ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
         IF (.NOT.ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
         iab = ikind + nkind*(jkind - 1)
         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         IF (iatom == last_iatom .AND. jatom == last_jatom) THEN
            new_atom_pair = .FALSE.
         ELSE
            new_atom_pair = .TRUE.
            last_jatom = jatom
            last_iatom = iatom
         END IF

         ! *** Use the symmetry of the first derivatives ***
         IF (iatom == jatom) THEN
            f0 = 1.0_dp
         ELSE
            f0 = 2.0_dp
         END IF

         ! *** Create matrix blocks for a new matrix block column ***
         IF (new_atom_pair) THEN
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
            IF (calculate_forces) THEN
               NULLIFY(p_block)
               CALL cp_dbcsr_get_block_p(matrix_p(1)%matrix,irow,icol,p_block,found)
            END IF
         END IF

         ! loop over all kinds for projector atom
         IF(ASSOCIATED(h_block)) THEN
            DO kkind=1,nkind
               iac= ikind + nkind*(kkind - 1)
               ibc= jkind + nkind*(kkind - 1)
               IF (.NOT.ASSOCIATED(sap_int(iac)%alist)) CYCLE
               IF (.NOT.ASSOCIATED(sap_int(ibc)%alist)) CYCLE
               CALL get_alist(sap_int(iac), alist_ac, iatom, error)
               CALL get_alist(sap_int(ibc), alist_bc, jatom, error)
               IF (.NOT.ASSOCIATED(alist_ac)) CYCLE
               IF (.NOT.ASSOCIATED(alist_bc)) CYCLE
               DO kac=1,alist_ac%nclist
                  DO kbc=1,alist_bc%nclist
                     IF ( alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom ) CYCLE
                     IF ( ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0) ) THEN
                        IF ( alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl ) CYCLE
                        acint => alist_ac%clist(kac)%acint
                        bcint => alist_bc%clist(kbc)%acint
                        achint => alist_ac%clist(kac)%achint
                        bchint => alist_bc%clist(kbc)%achint
                        na = SIZE(acint,1)
                        np = SIZE(acint,2)
                        nb = SIZE(bcint,1)
!$OMP CRITICAL(h_block_critical)
                        IF (iatom <= jatom) THEN
                           CALL dgemm("N","T",na,nb,np,1._dp,achint(1,1,1),na,&
                                bcint(1,1,1),nb,1.0_dp,h_block,SIZE(h_block,1))
                        ELSE
                           CALL dgemm("N","T",nb,na,np,1.0_dp,bchint(1,1,1),nb,&
                                acint(1,1,1),na,1.0_dp,h_block,SIZE(h_block,1))
                        END IF
!$OMP END CRITICAL(h_block_critical)
                        IF (calculate_forces) THEN
                           IF(ASSOCIATED(p_block)) THEN
                              katom = alist_ac%clist(kac)%catom
                              atom_c = atom_of_kind(katom)
                              DO i=1,3
                                 j=i+1
                                 IF (iatom <= jatom) THEN
                                    fa(i)=SUM(p_block(1:na,1:nb)*&
                                         MATMUL(acint(1:na,1:np,j),TRANSPOSE(bchint(1:nb,1:np,1))))
                                    fb(i)=SUM(p_block(1:na,1:nb)*&
                                         MATMUL(achint(1:na,1:np,1),TRANSPOSE(bcint(1:nb,1:np,j))))
                                 ELSE
                                    fa(i)=SUM(p_block(1:nb,1:na)*&
                                         MATMUL(bchint(1:nb,1:np,1),TRANSPOSE(acint(1:na,1:np,j))))
                                    fb(i)=SUM(p_block(1:nb,1:na)*&
                                         MATMUL(bcint(1:nb,1:np,j),TRANSPOSE(achint(1:na,1:np,1))))
                                 END IF
!$OMP CRITICAL(force_critical)
                                 force(ikind)%gth_ppnl(i,atom_a) =force(ikind)%gth_ppnl(i,atom_a) + f0*fa(i)
                                 force(kkind)%gth_ppnl(i,atom_c) =force(kkind)%gth_ppnl(i,atom_c) - f0*fa(i)
                                 force(jkind)%gth_ppnl(i,atom_b) =force(jkind)%gth_ppnl(i,atom_b) + f0*fb(i)
                                 force(kkind)%gth_ppnl(i,atom_c) =force(kkind)%gth_ppnl(i,atom_c) - f0*fb(i)
!$OMP END CRITICAL(force_critical)
                              END DO  

                              IF (use_virial) THEN
                                 rac = alist_ac%clist(kac)%rac
                                 rbc = alist_bc%clist(kbc)%rac
!$OMP CRITICAL(virial_critical)
                                 CALL virial_pair_force ( virial%pv_virial, f0, fa, rac, error)
                                 CALL virial_pair_force ( virial%pv_virial, f0, fb, rbc, error)
!$OMP END CRITICAL(virial_critical)
                              END IF
                           ENDIF
                        END IF
                        EXIT ! We have found a match and there can be only one single match
                     END IF
                  END DO  
               END DO  
            END DO 
         ENDIF
      END DO  
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL release_sap_int(sap_int, error)

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(basis_set,potential,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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

    END IF !ppnl_present

    CALL timestop(handle)

  END SUBROUTINE build_core_ppnl

!==========================================================================================================

END MODULE core_ppnl
