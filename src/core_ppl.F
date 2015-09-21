!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Calculation of the local pseudopotential contribution to the core Hamiltonian 
!>         <a|V(local)|b> = <a|Sum e^a*rc**2|b>
!> \par History
!>      - core_ppnl refactored from qs_core_hamiltonian [Joost VandeVondele, 2008-11-01]
!>      - adapted for PPL [jhu, 2009-02-23]
!>      - OpenMP added [Iain Bethune, Fiona Reid, 2013-11-13]
!>      - Bug fix: correct orbital pointer range [07.2014,JGH]
!>      - k-point aware [07.2015,JGH]
! *****************************************************************************
MODULE core_ppl

  USE ai_overlap_ppl,                  ONLY: ppl_integral
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
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

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_ppl'

  PUBLIC :: build_core_ppl

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
!> \param sac_ppl ...
!> \param nimages ...
!> \param cell_to_index ...
! *****************************************************************************
  SUBROUTINE build_core_ppl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder,&
             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, nimages, cell_to_index)

    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
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
      DIMENSION(:), POINTER                  :: sab_orb, sac_ppl
    INTEGER, INTENT(IN)                      :: nimages
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_ppl', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nexp_max = 10

    INTEGER :: atom_a, atom_b, atom_c, handle, i, iatom, icol, ikind, img, &
      inode, irow, iset, jatom, jkind, jset, katom, kkind, ldai, ldsab, &
      maxco, maxder, maxl, maxlgto, maxlppl, maxnset, maxsgf, mepos, natom, &
      ncoa, ncob, nexp_lpot, nexp_ppl, nkind, nseta, nsetb, nthread, sgfa, &
      sgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cellind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, nct_lpot, npgfa, &
                                                npgfb, nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    INTEGER, DIMENSION(nexp_max)             :: nct_ppl
    LOGICAL                                  :: dokp, failure, found, &
                                                lpotextended
    REAL(KIND=dp)                            :: alpha, dab, dac, dbc, f0, &
                                                ppl_radius
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ppl_fwork, ppl_work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: hab, pab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, rab, rac, &
                                                rbc
    REAL(KIND=dp), DIMENSION(4, nexp_max)    :: cval_ppl
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lpot, cexp_ppl, &
                                                set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cval_lpot, h_block, p_block, &
                                                rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    REAL(KIND=dp), DIMENSION(nexp_max)       :: alpha_ppl
    TYPE(gth_potential_type), POINTER        :: gth_potential
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

    dokp = (nimages > 1)

    ALLOCATE (atom_of_kind(natom))
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    IF (calculate_forces) THEN
       IF (SIZE(matrix_p,1) == 2) THEN
          DO img=1,nimages
             CALL cp_dbcsr_add(matrix_p(1,img)%matrix, matrix_p(2,img)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp)
             CALL cp_dbcsr_add(matrix_p(2,img)%matrix, matrix_p(1,img)%matrix, &
                         alpha_scalar=-2.0_dp, beta_scalar=1.0_dp)
          END DO
       END IF
    END IF

    maxder = ncoset(nder)

    CALL get_qs_kind_set(qs_kind_set, maxco=maxco,maxlgto=maxlgto,&
                         maxsgf=maxsgf,maxnset=maxnset, maxlppl=maxlppl)

    maxl = MAX(maxlgto,maxlppl)
    CALL init_orbital_pointers(2*maxl+2*nder+1)

    ldsab = MAX(maxco,ncoset(maxlppl),maxsgf,maxlppl)
    ldai = ncoset(maxl+nder+1)

    ALLOCATE (basis_set_list(nkind))
    DO ikind=1,nkind
      CALL get_qs_kind(qs_kind_set(ikind),basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO

    nthread=1
!$  nthread=omp_get_max_threads()

    CALL neighbor_list_iterator_create(nl_iterator,sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &                                                                                      
!$OMP SHARED  (nl_iterator, basis_set_list, calculate_forces, force, use_virial, &                             
!$OMP          matrix_h, matrix_p, virial, atomic_kind_set, qs_kind_set, particle_set, &
!$OMP          sab_orb, sac_ppl, nthread, ncoset, nkind, cell_to_index, &
!$OMP          atom_of_kind, ldsab,  maxnset, maxder, &                                                       
!$OMP          maxlgto, nder, maxco, dokp) &                                                                        
!$OMP PRIVATE (ikind, jkind, inode, iatom, jatom, rab, basis_set_a, basis_set_b, atom_a, &                   
!$OMP          atom_b, first_sgfa, la_max, la_min, npgfa, nsgfa, sphi_a, &                                    
!$OMP          zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, rpgfb, set_radius_b, sphi_b, &                 
!$OMP          zetb, dab, irow, icol, h_block, found, iset, ncoa, &                  
!$OMP          sgfa, jset, ncob, sgfb, nsgfb, p_block, work, pab, hab, kkind, nseta, &                       
!$OMP          gth_potential, alpha, cexp_ppl, lpotextended, &
!$OMP          ppl_radius, nexp_lpot, nexp_ppl, alpha_ppl, alpha_lpot, nct_ppl, &                            
!$OMP          nct_lpot, cval_ppl, cval_lpot, ap_iterator, rac, dac, rbc, dbc, &                             
!$OMP          set_radius_a,  rpgfa, force_a, force_b, ppl_fwork, mepos, &                                   
!$OMP          f0, failure, katom, ppl_work, atom_c, cellind, img, &                                                        
!$OMP          ldai)                 

        mepos=0
!$      mepos=omp_get_thread_num()

    ! iterator for basis/potential list
    CALL neighbor_list_iterator_create(ap_iterator,sac_ppl,search=.TRUE.)

    ALLOCATE(hab(ldsab,ldsab,maxnset,maxnset),work(ldsab,ldsab*maxder))
      ldai = ncoset(2*maxlgto+2*nder)
      ALLOCATE (ppl_work(ldai,ldai,MAX(maxder,2*maxlgto+2*nder+1)))
      IF (calculate_forces) THEN
         ALLOCATE(pab(maxco,maxco,maxnset,maxnset))
         ldai = ncoset(maxlgto)
         ALLOCATE (ppl_fwork(ldai,ldai,maxder))
      END IF
      
    DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos)==0)

       CALL get_iterator_info(nl_iterator,mepos=mepos,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab,cell=cellind)

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

       dab = SQRT(SUM(rab*rab))

       IF(dokp) THEN
          img = cell_to_index(cellind(1),cellind(2),cellind(3))
       ELSE
          img = 1
       END IF

       ! *** Use the symmetry of the first derivatives ***
       IF (iatom == jatom) THEN
          f0 = 1.0_dp
       ELSE
          f0 = 2.0_dp
       END IF

       ! *** Create matrix blocks for a new matrix block column ***
       IF (iatom <= jatom) THEN
          irow = iatom
          icol = jatom
       ELSE
          irow = jatom
          icol = iatom
       END IF
       NULLIFY(h_block)
       CALL cp_dbcsr_get_block_p(matrix_h(1,img)%matrix,irow,icol,h_block,found)
       IF(ASSOCIATED(h_block)) THEN
          IF (calculate_forces) THEN
             NULLIFY(p_block)
             CALL cp_dbcsr_get_block_p(matrix_p(1,img)%matrix,irow,icol,p_block,found)
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

       hab = 0._dp

       ! loop over all kinds for pseudopotential atoms
       DO kkind=1,nkind

          CALL get_qs_kind(qs_kind_set(kkind),gth_potential=gth_potential)
          IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
          CALL get_potential(potential=gth_potential,&
               alpha_ppl=alpha,cexp_ppl=cexp_ppl,&
               lpot_present=lpotextended,ppl_radius=ppl_radius)
               nexp_ppl=1
               alpha_ppl(1)=alpha
               nct_ppl(1)=SIZE(cexp_ppl)
               cval_ppl(1:nct_ppl(1),1)=cexp_ppl(1:nct_ppl(1))
          IF (lpotextended) THEN
            CALL get_potential(potential=gth_potential,&
                 nexp_lpot=nexp_lpot,alpha_lpot=alpha_lpot,nct_lpot=nct_lpot,cval_lpot=cval_lpot)
            CPASSERT(nexp_lpot<nexp_max)
            nexp_ppl=nexp_lpot+1
            alpha_ppl(2:nexp_lpot+1) = alpha_lpot(1:nexp_lpot)
            nct_ppl(2:nexp_lpot+1) = nct_lpot(1:nexp_lpot)
            DO i=1,nexp_lpot
              cval_ppl(1:nct_lpot(i),i+1) = cval_lpot(1:nct_lpot(i),i)
            END DO
          ENDIF

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

                      CALL ppl_integral(&
                              la_max(iset),la_min(iset),npgfa(iset),&
                              rpgfa(:,iset),zeta(:,iset),&
                              lb_max(jset),lb_min(jset),npgfb(jset),&
                              rpgfb(:,jset),zetb(:,jset),&
                              nexp_ppl,alpha_ppl,nct_ppl,cval_ppl,ppl_radius,&
                              rab,dab,rac,dac,rbc,dbc,&
                              hab(:,:,iset,jset),ppl_work,pab(:,:,iset,jset),&
                              force_a,force_b,ppl_fwork)
                      ! *** The derivatives w.r.t. atomic center c are    ***
                      ! *** calculated using the translational invariance ***
                      ! *** of the first derivatives                      ***
                      atom_c = atom_of_kind(katom)

!$OMP CRITICAL(force_critical)
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
                         CALL virial_pair_force ( virial%pv_virial, f0, force_a, rac)
                         CALL virial_pair_force ( virial%pv_virial, f0, force_b, rbc)
                      END IF
!$OMP END CRITICAL(force_critical)

                   ELSE
                      CALL ppl_integral(&
                           la_max(iset),la_min(iset),npgfa(iset),&
                           rpgfa(:,iset),zeta(:,iset),&
                           lb_max(jset),lb_min(jset),npgfb(jset),&
                           rpgfb(:,jset),zetb(:,jset),&
                           nexp_ppl,alpha_ppl,nct_ppl,cval_ppl,ppl_radius,&
                           rab,dab,rac,dac,rbc,dbc,hab(:,:,iset,jset),ppl_work)
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

!$OMP CRITICAL(h_block_critical)
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
!$OMP END CRITICAL(h_block_critical)

          END DO
       END DO  
    END DO  

    CALL neighbor_list_iterator_release(ap_iterator)
        
    DEALLOCATE(hab,work,ppl_work)
        
    IF (calculate_forces) THEN
       DEALLOCATE(pab,ppl_fwork)
    END IF

!$OMP END PARALLEL

    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (atom_of_kind,basis_set_list)

    IF (calculate_forces) THEN
       ! *** If LSD, then recover alpha density and beta density     ***
       ! *** from the total density (1) and the spin density (2)     ***
       IF (SIZE(matrix_p,1) == 2) THEN
          DO img=1,nimages
             CALL cp_dbcsr_add(matrix_p(1,img)%matrix, matrix_p(2,img)%matrix, &
                         alpha_scalar= 0.5_dp, beta_scalar=0.5_dp)
             CALL cp_dbcsr_add(matrix_p(2,img)%matrix, matrix_p(1,img)%matrix, &
                         alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
          END DO
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_core_ppl

!==========================================================================================================

END MODULE core_ppl
