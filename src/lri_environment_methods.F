!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates integral matrices for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created JGH [08.2012]
!>      Dorothea Golze [02.2014] (1) extended, re-structured, cleaned
!>                               (2) heavily debugged
!> \authors JGH 
!>          Dorothea Golze  
! *****************************************************************************
MODULE lri_environment_methods

  USE ai_overlap3,                     ONLY: overlap3
  USE ai_overlap3_debug,               ONLY: init_os_overlap3,&
                                             os_overlap3
  USE ai_overlap_debug,                ONLY: init_os_overlap2,&
                                             os_overlap2
  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: &
       allocate_lri_coefs, allocate_lri_ints, allocate_lri_rhos, &
       deallocate_lri_ints, lri_density_create, lri_density_release, &
       lri_density_type, lri_environment_type, lri_int_type, lri_kind_type, &
       lri_list_type, lri_rhoab_type
  USE mathlib,                         ONLY: invmat
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_collocate_density,            ONLY: calculate_lri_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_methods'

  PUBLIC :: build_lri_matrices, calculate_lri_densities

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
! *****************************************************************************
  SUBROUTINE build_lri_matrices(lri_env,qs_env,calculate_forces,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_lri_matrices', &
      routineP = moduleN//':'//routineN

! calculate the integrals needed to do the local (2-center) expansion
! of the (pair) densities

    CALL calculate_lri_integrals(lri_env,qs_env,calculate_forces,error)

  END SUBROUTINE build_lri_matrices

! *****************************************************************************
!> \brief calculates integrals needed for the LRI density fitting,
!>        integrals are calculated once, before the SCF starts
! *****************************************************************************
  SUBROUTINE calculate_lri_integrals(lri_env,qs_env,calculate_forces,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iac, iatom, ikind, ilist, jatom, jkind, jneighbor, &
      nba, nbb, nfa, nfb, nkind, nlist, nn, nneighbor, stat
    LOGICAL                                  :: calc_abb, failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(3)              :: rab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: fbasa, fbasb, obasa, obasb
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_list_type), POINTER             :: lri_ints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(cell, fbasa, fbasb, lrii, lri_ints, nl_iterator, obasa, obasb,&
            particle_set,soo_list)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      CALL get_qs_env(qs_env=qs_env,nkind=nkind,&
                      cell=cell,particle_set=particle_set,error=error)

      IF ( ASSOCIATED(lri_env%lri_ints) ) THEN
        CALL deallocate_lri_ints (lri_env%lri_ints,error)
      END IF

      ! allocate matrices storing the LRI integrals
      CALL allocate_lri_ints(lri_env,lri_env%lri_ints,nkind,error)
      lri_ints => lri_env%lri_ints

      CALL neighbor_list_iterator_create(nl_iterator,soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator)==0)

         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
              iatom=iatom,jatom=jatom,r=rab)
        
         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set
         fbasa => lri_env%ri_basis(ikind)%gto_basis_set
         fbasb => lri_env%ri_basis(jkind)%gto_basis_set

         IF (.NOT.ASSOCIATED(obasa)) CYCLE
         IF (.NOT.ASSOCIATED(obasb)) CYCLE

         lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf
         nfa = fbasa%nsgf
         nfb = fbasb%nsgf

         lrii%nba = nba
         lrii%nbb = nbb
         lrii%nfa = nfa
         lrii%nfb = nfb

         ! calculate integrals (a,b,fa) and (a,b,fb)
         IF(iatom == jatom) THEN
          CALL lri_int_abc(lrii%abaint,lrii%abbint,rab=rab,oba=obasa,obb=obasb,&
                           fba=fbasa,fbb=fbasb,calculate_forces=.FALSE.,&
                           calc_abb=.FALSE.,error=error)
          lrii%dabdaint = 0.0_dp
          lrii%dabbint  = 0.0_dp
         ELSE
          CALL lri_int_abc(lrii%abaint,lrii%abbint,lrii%dabdaint,&
                           lrii%dabbint,rab,obasa,obasb,fbasa,fbasb,&
                           calculate_forces,calc_abb=.TRUE.,error=error)
         ENDIF

         ! calculate integrals (fa,fb); for iatom=jatom this is the self-overlap
         IF(iatom == jatom) THEN
           lrii%sab(1:nfa,1:nfa)=lri_env%bas_ovlp(ikind)%ri_ovlp(1:nfa,1:nfa)
           lrii%dsab = 0._dp
         ELSE
           CALL lri_int_ab(lrii%sab,lrii%dsab,rab,fbasa,fbasb,calculate_forces,error)
         ENDIF

         ! construct and invert S matrix 
         lrii%sinv(1:nfa,1:nfa) = lri_env%bas_ovlp(ikind)%ri_ovlp(1:nfa,1:nfa)
         IF(iatom /= jatom) THEN
          nn = nfa+nfb
          lrii%sinv(1:nfa,nfa+1:nn) = lrii%sab(1:nfa,1:nfb)
          lrii%sinv(nfa+1:nn,1:nfa) = TRANSPOSE(lrii%sab(1:nfa,1:nfb))
          lrii%sinv(nfa+1:nn,nfa+1:nn) = lri_env%bas_ovlp(jkind)%ri_ovlp(1:nfb,1:nfb)
         ENDIF
         CALL invmat(lrii%sinv,stat,error)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ! calculate Sinv*n and n*Sinv*n
         lrii%n(1:nfa) = lri_env%bas_int(ikind)%int_fbas(1:nfa)
         IF(iatom == jatom) THEN
          lrii%sn(1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),lrii%n(1:nfa))
          lrii%nsn = SUM(lrii%sn(1:nfa)*lrii%n(1:nfa))
         ELSE
          lrii%n(nfa+1:nn) = lri_env%bas_int(jkind)%int_fbas(1:nfb)
          lrii%sn(1:nn) = MATMUL(lrii%sinv(1:nn,1:nn),lrii%n(1:nn))
          lrii%nsn = SUM(lrii%sn(1:nn)*lrii%n(1:nn))
         ENDIF

         ! calculate integrals (a,b), overlap of primary basis
         IF(iatom == jatom) THEN
          lrii%soo(1:nba,1:nba) =  lri_env%bas_ovlp(ikind)%orb_ovlp(1:nba,1:nba)
          lrii%dsoo = 0._dp
         ELSE
          CALL lri_int_ab(lrii%soo,lrii%dsoo,rab,obasa,obasb,calculate_forces,&
                          error=error)
         ENDIF
         
         ! calculate derivative of fit coefficients, needed for update of KS matrix
         CALL lri_calculate_derivative_acoef(lrii,iatom,jatom,nba,nbb,nfa,nfb,error)
 
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_integrals

! *****************************************************************************
!> \brief performs the fitting of the density and distributes the fitted 
!>        density on the grid
!> \param lri_env the lri environment
!>        lri_density the environment for the fitting
!>        pmatrix density matrix
!>        lri_rho_struct where the fitted density is stored 
! *****************************************************************************
  SUBROUTINE calculate_lri_densities(lri_env,lri_density,qs_env,pmatrix,&
                                     lri_rho_struct,atomic_kind_set,para_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(qs_rho_type), POINTER               :: lri_rho_struct
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_densities', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iac, iatom, ikind, ilist, ispin, &
      jatom, jkind, jneighbor, nat, natom, nba, nbb, nfa, nfb, nkind, nlist, &
      nn, nneighbor, nspin, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: aci, acj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pbij
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_type), POINTER             :: pmat
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_coef
    TYPE(lri_list_type), POINTER             :: lri_rho
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(aci, acj, atomic_kind, atom_of_kind, lri_coef, lrii, lri_rho,&
            nl_iterator, pbij, pmat, soo_list)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nspin = SIZE(pmatrix)
      nkind = lri_env%lri_ints%nkind

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=nat)
      ALLOCATE(atom_of_kind(nat),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL lri_density_release(lri_density,error)
      CALL lri_density_create(lri_density,error)
      lri_density%nspin = nspin

      ! allocate structure lri_rhos and vectors tvec and avec
      CALL allocate_lri_rhos(lri_env,lri_density%lri_rhos,nspin,nkind,error)

      DO ispin = 1, nspin
         pmat => pmatrix(ispin)%matrix
         lri_rho => lri_density%lri_rhos(ispin)%lri_list
 
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
                 jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)
  
            iac = ikind + nkind*(jkind - 1)

            IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

            ! get the density matrix Pab
            NULLIFY (pbij)
            IF (iatom <= jatom) THEN
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=iatom,col=jatom,block=pbij,found=found)
              trans = .FALSE.
            ELSE
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=jatom,col=iatom,block=pbij,found=found)
              trans = .TRUE.
            END IF
            CPPostcondition(found,cp_failure_level,routineP,error,failure)

            lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
  
            nba = lrii%nba
            nbb = lrii%nbb
            nfa = lrii%nfa
            nfb = lrii%nfb

            nn = nfa + nfb

            ! compute tvec = SUM_ab Pab *(a,b,x) and charge contraint 
            IF (trans) THEN
               lrho%charge = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%soo(1:nba,1:nbb)) 
               DO i=1,nfa
                  lrho%tvec(i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abaint(1:nba,1:nbb,i))
               END DO
               IF(iatom /= jatom) THEN
                DO i=1,nfb
                   lrho%tvec(nfa+i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abbint(1:nba,1:nbb,i))
                END DO
               ENDIF
            ELSE
               lrho%charge = SUM(pbij(1:nba,1:nbb)*lrii%soo(1:nba,1:nbb))
               DO i=1,nfa
                  lrho%tvec(i) = SUM(pbij(1:nba,1:nbb)*lrii%abaint(1:nba,1:nbb,i))
               END DO
               IF(iatom /= jatom) THEN
                DO i=1,nfb
                   lrho%tvec(nfa+i) = SUM(pbij(1:nba,1:nbb)*lrii%abbint(1:nba,1:nbb,i))
                END DO
               ENDIF
            END IF

            IF(iatom == jatom) THEN
             lrho%nst = SUM(lrho%tvec(1:nfa) * lrii%sn(1:nfa))
            ELSE
             lrho%nst = SUM(lrho%tvec(1:nn) * lrii%sn(1:nn))
            ENDIF
            lrho%lambda = (lrho%charge - lrho%nst)/lrii%nsn

            ! solve the linear system of equations
            ALLOCATE(m(nn),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            m = 0._dp
            IF(iatom == jatom) THEN
              m(1:nfa) = lrho%tvec(1:nfa) + lrho%lambda * lrii%n(1:nfa)
              lrho%avec(1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),m(1:nfa))
            ELSE
              m(1:nn) = lrho%tvec(1:nn) + lrho%lambda * lrii%n(1:nn)
              lrho%avec(1:nn)  = MATMUL(lrii%sinv(1:nn,1:nn),m(1:nn))
            ENDIF
            DEALLOCATE(m,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

      END DO

      ! allocate the arrays to hold RI expansion coefficients lri_coefs
      CALL allocate_lri_coefs(lri_env,lri_density,atomic_kind_set,error)

      DO ispin = 1, nspin

         lri_coef => lri_density%lri_coefs(ispin)%lri_kinds  
         lri_rho  => lri_density%lri_rhos(ispin)%lri_list

         ! sum up expansion coefficients
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
                 iatom=iatom,jatom=jatom,ilist=ilist,inode=jneighbor)
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            aci => lri_coef(ikind)%acoef(atom_a,:)
            acj => lri_coef(jkind)%acoef(atom_b,:)
            iac = ikind + nkind*(jkind - 1)
            lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            nfa = lrho%nfa
            nfb = lrho%nfb
            IF(iatom == jatom) THEN
             aci(1:nfa) = aci(1:nfa) + lrho%avec(1:nfa)
            ELSE
             aci(1:nfa) = aci(1:nfa) + 2.0_dp*lrho%avec(1:nfa)
             acj(1:nfb) = acj(1:nfb) + 2.0_dp*lrho%avec(nfa+1:nfa+nfb)
            ENDIF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         ! replicate the acoef infomation 
         DO ikind=1,nkind
           atomic_kind => atomic_kind_set(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
           DO iatom=1,natom
              aci => lri_coef(ikind)%acoef(iatom,:)
              CALL mp_sum(aci,para_env%group)
           END DO
         END DO

      END DO

      !distribute fitted density on the grid 
      DO ispin=1,nspin
       CALL calculate_lri_rho_elec(lri_rho_struct%rho_g(ispin),&
                                   lri_rho_struct%rho_r(ispin), qs_env, lri_env,&
                                   lri_density%lri_coefs(ispin)%lri_kinds,&
                                   lri_rho_struct%tot_rho_r(ispin), error)
      ENDDO

      CALL set_qs_env(qs_env, lri_density=lri_density, error=error)

      DEALLOCATE(atom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF
    CALL timestop(handle)

  END SUBROUTINE calculate_lri_densities

! *****************************************************************************
!> \brief calculate integrals (a,b,fa) and (a,b,fb)
! *****************************************************************************
  SUBROUTINE lri_int_abc(abaint,abbint,dabdaint,dabbint,rab,oba,obb,&
                         fba,fbb,calculate_forces,calc_abb,error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abaint, abbint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: dabdaint, dabbint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces, calc_abb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_abc', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iset, jset, kaset, kbset, m1, m2, m3, ncoa, ncob, &
      ncoc, nseta, nsetb, nsetca, nsetcb, sgfa, sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, &
      lca_max, lca_min, lcb_max, lcb_min, npgfa, npgfb, npgfca, npgfcb, &
      nsgfa, nsgfb, nsgfca, nsgfcb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfca, first_sgfcb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: saba, sabb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sabda, sdaba, sdabb, sdabda
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfca, rpgfcb, &
                                                sphi_a, sphi_b, sphi_ca, &
                                                sphi_cb, zeta, zetb, zetca, &
                                                zetcb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, lca_max, lca_min, lcb_max, lcb_min,&
            npgfa, npgfb, npgfca, npgfcb, nsgfa, nsgfb, nsgfca, nsgfcb)
    NULLIFY(first_sgfa, first_sgfb, first_sgfca, first_sgfcb, set_radius_a,&
            set_radius_b, rpgfa, rpgfb, rpgfca, rpgfcb, sphi_a, sphi_b, sphi_ca,&
            sphi_cb, zeta, zetb, zetca, zetcb) 

    ! basis ikind
    first_sgfa   =>  oba%first_sgf
    la_max       =>  oba%lmax
    la_min       =>  oba%lmin
    npgfa        =>  oba%npgf
    nseta        =   oba%nset
    nsgfa        =>  oba%nsgf_set
    rpgfa        =>  oba%pgf_radius
    set_radius_a =>  oba%set_radius
    sphi_a       =>  oba%sphi
    zeta         =>  oba%zet
    ! basis jkind
    first_sgfb   =>  obb%first_sgf
    lb_max       =>  obb%lmax
    lb_min       =>  obb%lmin
    npgfb        =>  obb%npgf
    nsetb        =   obb%nset
    nsgfb        =>  obb%nsgf_set
    rpgfb        =>  obb%pgf_radius
    set_radius_b =>  obb%set_radius
    sphi_b       =>  obb%sphi
    zetb         =>  obb%zet

    ! basis RI A
    first_sgfca  =>  fba%first_sgf
    lca_max      =>  fba%lmax
    lca_min      =>  fba%lmin
    npgfca       =>  fba%npgf
    nsetca       =   fba%nset
    nsgfca       =>  fba%nsgf_set
    rpgfca       =>  fba%pgf_radius
    sphi_ca      =>  fba%sphi
    zetca        =>  fba%zet
    ! basis RI B
    first_sgfcb  =>  fbb%first_sgf
    lcb_max      =>  fbb%lmax
    lcb_min      =>  fbb%lmin
    npgfcb       =>  fbb%npgf
    nsetcb       =   fbb%nset
    nsgfcb       =>  fbb%nsgf_set
    rpgfcb       =>  fbb%pgf_radius
    sphi_cb      =>  fbb%sphi
    zetcb        =>  fbb%zet

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1

          ! calculate integrals abaint and derivative [d(a,b,a)/dA] dabdaint if requested
          rac = 0._dp
          dac = 0._dp
          rbc = -rab
          dbc = dab
          DO kaset=1,nsetca
             ncoc = npgfca(kaset)*ncoset(lca_max(kaset))
             sgfc = first_sgfca(1,kaset)
             m3=sgfc+nsgfca(kaset)-1
             IF(ncoa*ncob*ncoc > 0) THEN
                ALLOCATE(saba(ncoa,ncob,ncoc),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                saba(:,:,:)    = 0._dp
                ! integrals
                IF(calculate_forces) THEN
                 ALLOCATE(sdaba(ncoa,ncob,ncoc,3),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 ALLOCATE(sabda(ncoa,ncob,ncoc,3),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 ALLOCATE(sdabda(ncoa,ncob,ncoc,3),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 sdaba(:,:,:,:) = 0._dp
                 sabda(:,:,:,:) = 0._dp
                 sdabda(:,:,:,:) = 0._dp
                 CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                               lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                               lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                               rab,dab,rac,dac,rbc,dbc,saba,sdaba,sabda,error=error)
                 !d(a,b,a)/dA = (da/dA,b,a) + (a,b,da/dA)
                 sdabda = sdaba + sabda

                 DO i=1,3
                  CALL abc_contract(dabdaint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabda(:,:,:,i),& 
                       sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                       ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset),error)
                 ENDDO

                 DEALLOCATE(sdaba,sabda,sdabda,STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ELSE
                 CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                               lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                               lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                               rab,dab,rac,dac,rbc,dbc,saba,error=error)
                ENDIF
                CALL abc_contract(abaint(sgfa:m1,sgfb:m2,sgfc:m3),saba,& 
                     sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset),error)
                DEALLOCATE(saba,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO

          ! calculate integrals abbint and derivative [d(a,b,b)/dA] dabbint if requested
          IF(calc_abb) THEN
            rac = rab
            dac = dab
            rbc = 0._dp
            dbc = 0._dp
            DO kbset=1,nsetcb
               ncoc = npgfcb(kbset)*ncoset(lcb_max(kbset))
               sgfc = first_sgfcb(1,kbset)
               m3=sgfc+nsgfcb(kbset)-1
               IF(ncoa*ncob*ncoc > 0) THEN
                  ALLOCATE(sabb(ncoa,ncob,ncoc),STAT=stat)
                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                  sabb(:,:,:) = 0._dp
                  IF(calculate_forces) THEN
                   ALLOCATE(sdabb(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   sdabb(:,:,:,:) = 0._dp
                  ! integrals
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb,sdabc=sdabb,error=error)
                   DO i=1,3
                    CALL abc_contract(dabbint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabb(:,:,:,i),&
                         sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                         ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset),error)
                   ENDDO
                   DEALLOCATE(sdabb,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                  ELSE
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb,error=error)
                  ENDIF
                  CALL abc_contract(abbint(sgfa:m1,sgfb:m2,sgfc:m3),sabb,&
                       sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                       ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset),error)
                  DEALLOCATE(sabb,STAT=stat)
                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               END IF
            END DO
          ENDIF

       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE lri_int_abc

! *****************************************************************************
!> \brief calculate overlap integrals (a,b)
! *****************************************************************************
  SUBROUTINE lri_int_ab(sab,dsab,rab,fba,fbb,calculate_forces,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dsab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: fba, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_ab', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iset, jset, lds, m1, m2, maxco, maxcoa, maxcob, &
      maxl, maxla, maxlb, ncoa, ncob, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: devab, swork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb,&
            first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb,&
            sphi_a, sphi_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  fba%first_sgf
    la_max       =>  fba%lmax
    la_min       =>  fba%lmin
    npgfa        =>  fba%npgf
    nseta        =   fba%nset
    nsgfa        =>  fba%nsgf_set
    rpgfa        =>  fba%pgf_radius
    set_radius_a =>  fba%set_radius
    sphi_a       =>  fba%sphi
    zeta         =>  fba%zet
    ! basis jkind
    first_sgfb   =>  fbb%first_sgf
    lb_max       =>  fbb%lmax
    lb_min       =>  fbb%lmin
    npgfb        =>  fbb%npgf
    nsetb        =   fbb%nset
    nsgfb        =>  fbb%nsgf_set
    rpgfb        =>  fbb%pgf_radius
    set_radius_b =>  fbb%set_radius
    sphi_b       =>  fbb%sphi
    zetb         =>  fbb%zet

    CALL get_gto_basis_set(fba,maxco=maxcoa,maxl=maxla)
    CALL get_gto_basis_set(fbb,maxco=maxcob,maxl=maxlb)
    maxco = MAX(maxcoa,maxcob)
    IF(calculate_forces) THEN
     maxl = MAX(maxla+1,maxlb)
    ELSE
     maxl = MAX(maxla,maxlb)
    ENDIF
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(calculate_forces) THEN
     !derivatives will be stored in devab(:,:,2:4)
     ALLOCATE(swork(lds,lds,4),devab(maxco,maxco,4),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     devab = 0._dp
    ELSE
     ALLOCATE(swork(lds,lds,1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    sint  = 0._dp
    swork = 0._dp

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1
          sint=0._dp
          swork = 0._dp

          ! calculate integrals
          IF(calculate_forces) THEN
             devab = 0._dp 
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds,sdab=devab)
             DO i=1,3
              !NOTE: devab(:,:,2:4) contains all derivatives for lmin=0 to lmax=lmax 
              !      correct after contraction (multiply with zero for elements l < lmin) 
              CALL ab_contract(dsab(sgfa:m1,sgfb:m2,i),devab(:,:,i+1),sphi_a(:,sgfa:),&
                   sphi_b(:,sgfb:),ncoa,ncob,nsgfa(iset),nsgfb(jset),error)
             ENDDO

          ELSE
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds)
          ENDIF

          CALL ab_contract(sab(sgfa:m1,sgfb:m2),sint,sphi_a(:,sgfa:),sphi_b(:,sgfb:),&
               ncoa,ncob,nsgfa(iset),nsgfb(jset),error)
       END DO
    END DO

    IF(calculate_forces) THEN
     DEALLOCATE(devab,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_ab

! *****************************************************************************
!> \brief contract three-center overlap integrals (a,b,c) and transfer
!>        to spherical Gaussians
! *****************************************************************************
  SUBROUTINE abc_contract(abcint,sabc,sphi_a,sphi_b,sphi_c,ncoa,ncob,ncoc,&
                          nsgfa,nsgfb,nsgfc,error)

    REAL(KIND=dp), DIMENSION(:, :, :)        :: abcint, sabc
    REAL(KIND=dp), DIMENSION(:, :)           :: sphi_a, sphi_b, sphi_c
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, nsgfa, &
                                                nsgfb, nsgfc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'abc_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isgfc, m1, m2, m3, msphia, &
                                                msphib, msphic, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: m_work, temp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cpc, cpp

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)
    msphic = SIZE(sphi_c,1)

    m1 = SIZE(sabc,1)
    m2 = SIZE(sabc,2)
    m3 = SIZE(sabc,3)

    ALLOCATE(cpp(nsgfa,m2,m3),cpc(nsgfa,m2,nsgfc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL dgemm("T","N",nsgfa,m2*m3,ncoa,1._dp,sphi_a,msphia,sabc,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa*m2,nsgfc,ncoc,1._dp,cpp,nsgfa*m2,sphi_c,msphic,0.0_dp,&
               cpc,nsgfa*m2)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(m_work(nsgfa,ncob),temp(nsgfa,nsgfb),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    m_work(:,:)=0.0_dp
    temp(:,:)=0.0_dp

    DO isgfc=1,nsgfc
     m_work(1:nsgfa,1:ncob)=cpc(:,:,isgfc)
     CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,m_work,nsgfa,sphi_b,msphib,&
                0.0_dp,temp,nsgfa)
     abcint(:,:,isgfc)=temp(:,:)
    END DO

    DEALLOCATE(cpc,m_work,temp,stat=STAT)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE abc_contract

! *****************************************************************************
!> \brief contract overlap integrals (a,b) and transfer to spherical Gaussians
! *****************************************************************************
  SUBROUTINE ab_contract(abint,sab,sphi_a,sphi_b,ncoa,ncob,nsgfa,nsgfb,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: abint
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab, sphi_a, sphi_b
    INTEGER, INTENT(IN)                      :: ncoa, ncob, nsgfa, nsgfb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ab_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m1, m2, msphia, msphib, nn, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cpp

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)

    m1 = SIZE(sab,1)
    m2 = SIZE(sab,2)

    nn = SIZE(abint,1)

    ALLOCATE(cpp(nsgfa,m2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL dgemm("T","N",nsgfa,m2,ncoa,1._dp,sphi_a,msphia,sab,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,cpp,nsgfa,sphi_b,msphib,0.0_dp,&
               abint,nn)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE ab_contract

! *****************************************************************************
!> \brief calculate derivate of fit coefficients acoef with respect to 
!>        density matrix pmatrix 
!>        R = (a,b)/nsn  - SUM_i (a,b,ai)*sn(i)/nsn 
!>        Q = SUM_i sinv*(a,b,ai)
!>        derviate_aci = R + Q
!> \param nba, nbb number of primary basis functions on a/b 
!> \param nfa, nfb number of ri basis functions on a/b 
! *****************************************************************************
  SUBROUTINE lri_calculate_derivative_acoef(lrii,iatom,jatom,nba,nbb,nfa,nfb,error)

    TYPE(lri_int_type), POINTER              :: lrii
    INTEGER, INTENT(IN)                      :: iatom, jatom, nba, nbb, nfa, &
                                                nfb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'lri_calculate_derivative_acoef', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nn, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: abxint_work, Q, R

    failure = .FALSE.
    CALL timeset(routineN,handle)

    nn = nfa + nfb
 
    ALLOCATE(R(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(Q(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(abxint_work(nba,nbb,nn),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    R(:,:,:)=0._dp
    Q(:,:,:)=0._dp
    abxint_work(:,:,:) =0._dp

    abxint_work(1:nba,1:nbb,1:nfa)=lrii%abaint(1:nba,1:nbb,1:nfa)    !abaint
    abxint_work(1:nba,1:nbb,nfa+1:nn)=lrii%abbint(1:nba,1:nbb,1:nfb) !abbint
    IF(iatom == jatom) THEN
     DO i=1,nba
      DO j=1,nbb
       R(i,j,1:nfa) = MATMUL(lrii%sinv(1:nfa,1:nfa),abxint_work(i,j,1:nfa))
       Q(i,j,1:nfa) = lrii%soo(i,j)/lrii%nsn*lrii%sn(1:nfa) - &
                      DOT_PRODUCT(lrii%sn(1:nfa),abxint_work(i,j,1:nfa))/lrii%nsn*lrii%sn(1:nfa)
       lrii%dacoef(1:nba,1:nbb,1:nfa)= R(1:nba,1:nbb,1:nfa)+Q(1:nba,1:nbb,1:nfa)
      ENDDO
     ENDDO
    ELSE
     DO i=1,nba
      DO j=1,nbb
       R(i,j,1:nn) =  MATMUL(lrii%sinv(1:nn,1:nn),abxint_work(i,j,1:nn))
       Q(i,j,1:nn) =  lrii%soo(i,j)/lrii%nsn*lrii%sn(1:nn) - &
                      DOT_PRODUCT(lrii%sn(1:nn),abxint_work(i,j,1:nn))/lrii%nsn*lrii%sn(1:nn)
       lrii%dacoef(1:nba,1:nbb,1:nn)= R(1:nba,1:nbb,1:nn)+Q(1:nba,1:nbb,1:nn)
      ENDDO
     ENDDO
    ENDIF

    DEALLOCATE(abxint_work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(R,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(Q,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
  
  END SUBROUTINE lri_calculate_derivative_acoef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *****************************************************************************
!  DEBUG ROUTINES ==> RECURSIVE INTEGRAL EVALUATION
!  TO BE DELETED LATER
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE overlap3_test(error)
   
    !INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap3_test', &
      routineP = moduleN//':'//routineN

    INTEGER :: ia1, iax, iay, iaz, ib1, ibx, iby, ibz, ic1, icx, icy, icz, &
      la_max, la_min, lb_max, lb_min, lc_max, lc_min, ma, mb, mc, stat
    INTEGER, DIMENSION(3)                    :: na, nb, nc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc, dmax, res1, &
                                                xa, xb, xc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: sabc
    REAL(KIND=dp), DIMENSION(1)              :: rpgfa, rpgfb, rpgfc, xa_work, &
                                                xb_work, xc_work
    REAL(KIND=dp), DIMENSION(3)              :: A, B, C, rab, rac, rbc

    xa = 0.783300000000   ! exponents
    xb = 1.239648746700
    xc = 0.548370000000
 
    A = (/0.329309000000,0.28408240000,0.28408240000/)!* bohr  !positions
    B = (/0.983983000000,0.00453720000,0.00432740000/)!* bohr
    C = (/0.032380000000,1.23470000000,0.11137400000/)!* bohr
   
    la_min = 0
    lb_min = 0
    lc_min = 0

    la_max = 0
    lb_max = 0
    lc_max = 1
  
    !---------------------------------------
    rab(:)= B(:) - A(:)
    dab = SQRT(DOT_PRODUCT(rab,rab))
    rac(:)= C(:) - A(:)
    dac = SQRT(DOT_PRODUCT(rac,rac))
    rbc(:)= C(:) - B(:)
    dbc = SQRT(DOT_PRODUCT(rbc,rbc))
    ALLOCATE(sabc(ncoset(la_max),ncoset(lb_max),ncoset(lc_max)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    xa_work(1)=xa
    xb_work(1)=xb
    xc_work(1)=xc
    rpgfa =20._dp
    rpgfb =20._dp
    rpgfc =20._dp
    sabc=0._dp

    CALL overlap3(la_max_set=la_max,npgfa=1,zeta=xa_work,rpgfa=rpgfa,la_min_set=la_min,&
                  lb_max_set=lb_max,npgfb=1,zetb=xb_work,rpgfb=rpgfb,lb_min_set=lb_min,&
                  lc_max_set=lc_max,npgfc=1,zetc=xc_work,rpgfc=rpgfc,lc_min_set=lc_min,&
                  rab=rab,dab=dab,rac=rac,dac=dac,rbc=rbc,dbc=dbc,sabc=sabc,error=error)
   ! write(*,*) "sabc", sabc 
    !sabc(:,:,:)=0._dp 
   !CALL overlap3_old(la_max=la_max,npgfa=1,zeta=xa_work,rpgfa=rpgfa,la_min=la_min,&
   !              lb_max=lb_max,npgfb=1,zetb=xb_work,rpgfb=rpgfb,lb_min=lb_min,&
   !              lc_max=lc_max,npgfc=1,zetc=xc,rpgfc=20._dp,lc_min=lc_min,&
   !              rab=rab,dab=dab,rac=rac,dac=dac,rbc=rbc,dbc=dbc,sabc=sabc)
  
  !---------------------------------------
 
    CALL init_os_overlap3(xa,xb,xc,A,B,C)

    dmax=0._dp
    DO ma=la_min,la_max
      DO mc=lc_min,lc_max
        DO mb=lb_min,lb_max
          DO iax=0,ma
            DO iay=0,ma-iax
               iaz=ma-iax-iay
               na(1)=iax; na(2)=iay; na(3)=iaz
               ia1=coset(iax,iay,iaz)
               DO icx=0,mc
                 DO icy=0,mc-icx
                    icz=mc-icx-icy
                    nc(1)=icx; nc(2)=icy; nc(3)=icz
                    ic1=coset(icx,icy,icz)
                    DO ibx=0,mb
                      DO iby=0,mb-ibx
                        ibz=mb-ibx-iby
                        nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                        ib1=coset(ibx,iby,ibz)
                        res1=os_overlap3(na,nc,nb)
                        !write(*,*) "la, lc, lb,na,nc, nb, res1", ma, mc, mb, na, nc, nb, res1
                       !write(*,*) "sabc ia1, ib1, ic1", ia1, ib1, ic1, sabc(ia1,ib1,ic1)
                        dmax=MAX(dmax,ABS(res1-sabc(ia1,ib1,ic1)))
                      END DO
                    END DO
                 END DO
               END DO
            END DO
          END DO
        END DO
      END DO
    END DO
   
    WRITE(*,*) "dmax", dmax

   DEALLOCATE(sabc,STAT=stat) 
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
  END SUBROUTINE overlap3_test
! ******************************************************************************
  SUBROUTINE overlap3_test_2(la_max,npgfa,zeta,la_min,&
                           lb_max,npgfb,zetb,lb_min,&
                           lc_max,npgfc,zetc,lc_min,&
                           ra,rb,rc,sabc,error)
   
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lb_min, lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: sabc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap3_test_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coa, cob, coc, ia1, iax, iay, &
                                                iaz, ib1, ibx, iby, ibz, ic1, &
                                                icx, icy, icz, ipgf, jpgf, &
                                                kpgf, ma, mb, mc
    INTEGER, DIMENSION(3)                    :: na, nb, nc
    REAL(KIND=dp)                            :: dmax, res1, res2, xa, xb, xc
    REAL(KIND=dp), DIMENSION(3)              :: A, B, C

    dmax=0._dp
    coa=0
    DO ipgf = 1, npgfa
     cob=0
     DO jpgf = 1, npgfb
      coc=0
      DO kpgf = 1, npgfc
 
       xa = zeta(ipgf)   ! exponents
       xb = zetb(jpgf)
       xc = zetc(kpgf)
 
       A = Ra  !positions
       B = Rb
       C = Rc
       
       CALL init_os_overlap3(xa,xb,xc,A,B,C)

       DO ma=la_min,la_max
         DO mc=lc_min,lc_max
           DO mb=lb_min,lb_max
             DO iax=0,ma
               DO iay=0,ma-iax
                  iaz=ma-iax-iay
                  na(1)=iax; na(2)=iay; na(3)=iaz
                  ia1=coset(iax,iay,iaz)    
                  DO icx=0,mc
                    DO icy=0,mc-icx
                       icz=mc-icx-icy
                       nc(1)=icx; nc(2)=icy; nc(3)=icz
                       ic1=coset(icx,icy,icz)
                       DO ibx=0,mb
                         DO iby=0,mb-ibx
                           ibz=mb-ibx-iby
                           nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                           ib1=coset(ibx,iby,ibz)
                           res1=os_overlap3(na,nc,nb)
                           res2=sabc(coa+ia1,cob+ib1,coc+ic1)
                          ! write(*,*) "la, na", ma, na
                          ! write(*,*) "lb, nb", mb, nb
                          ! write(*,*) "lc, nc", mc, nc
                          ! write(*,*) "sabc coa+ia1, cob+ib1, coc+ic1 ", &
                          !             coa+ia1, cob+ib1, coc+ic1, res2
                           dmax=MAX(dmax,ABS(res1-res2))
                           IF(dmax > 1.E-10) WRITE(*,*) "dmax in loop", dmax
                         END DO
                       END DO
                    END DO
                  END DO
               END DO
             END DO
           END DO
         END DO
       END DO

       coc = coc + ncoset(lc_max)
      ENDDO 
      cob = cob + ncoset(lb_max)
     ENDDO
     coa = coa + ncoset(la_max)
    ENDDO
    WRITE(*,*) "dmax", dmax
    
  END SUBROUTINE overlap3_test_2
! *****************************************************************************
! Two-center overlap testing routines
! *****************************************************************************
  SUBROUTINE overlap2_test(error)
   
    !INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap2_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia1, iax, iay, iaz, ib1, ibx, &
                                                iby, ibz, la_max, la_min, &
                                                lb_max, lb_min, lds, ma, &
                                                maxl, mb, stat
    INTEGER, DIMENSION(3)                    :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dmax, res1, xa, xb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: swork
    REAL(KIND=dp), DIMENSION(1)              :: rpgfa, rpgfb, xa_work, xb_work
    REAL(KIND=dp), DIMENSION(3)              :: A, B, rab

    xa = 0.783300000000   ! exponents
    xb = 1.239648746700
 
    A = (/0.329309000000,0.28408240000,0.28408240000/)!* bohr  !positions
    B = (/0.983983000000,0.00453720000,0.00432740000/)!* bohr
   
    la_min = 0
    lb_min = 0

    la_max = 3
    lb_max = 4
  

    !---------------------------------------
    ALLOCATE(sab(ncoset(la_max),ncoset(lb_max)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    maxl=MAX(la_max, lb_max)
    lds = ncoset(maxl)
    ALLOCATE(swork(lds,lds,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sab=0._dp
    rab(:)= B(:) - A(:)
    dab = SQRT(DOT_PRODUCT(rab,rab))
    xa_work(1)=xa
    xb_work(1)=xb
    rpgfa =20._dp
    rpgfb =20._dp
    CALL overlap(la_max_set=la_max,la_min_set=la_min,npgfa=1,rpgfa=rpgfa,zeta=xa_work,&
               lb_max_set=lb_max,lb_min_set=lb_min,npgfb=1,rpgfb=rpgfb,zetb=xb_work,&
               rab=rab,dab=dab,sab=sab,da_max_set=0,return_derivatives=.FALSE.,s=swork,lds=lds)
    !---------------------------------------

    CALL init_os_overlap2(xa,xb,A,B)

    dmax=0._dp
    DO ma=la_min,la_max
      DO mb=lb_min,lb_max
        DO iax=0,ma
          DO iay=0,ma-iax
             iaz=ma-iax-iay
             na(1)=iax; na(2)=iay; na(3)=iaz
             ia1=coset(iax,iay,iaz)   
              DO ibx=0,mb
                DO iby=0,mb-ibx
                  ibz=mb-ibx-iby
                  nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                  ib1=coset(ibx,iby,ibz)         
                  res1=os_overlap2(na,nb)
                 ! write(*,*) "la, lb,na, nb, res1", ma, mb, na, nb, res1
                 ! write(*,*) "sab ia1, ib1", ia1, ib1, sab(ia1,ib1)
                  dmax=MAX(dmax,ABS(res1-sab(ia1,ib1)))
                END DO
              END DO
          END DO
        END DO
      END DO
    END DO
   
    WRITE(*,*) "dmax", dmax

   DEALLOCATE(sab,swork,STAT=stat) 
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE overlap2_test
! ******************************************************************************
  SUBROUTINE overlap2_test_2(la_max,npgfa,zeta,la_min,&
                           lb_max,npgfb,zetb,lb_min,&
                           ra,rb,sab,error)
   
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap2_test_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coa, cob, ia1, iax, iay, iaz, &
                                                ib1, ibx, iby, ibz, ipgf, &
                                                jpgf, ma, mb
    INTEGER, DIMENSION(3)                    :: na, nb
    REAL(KIND=dp)                            :: dmax, res1, res2, xa, xb
    REAL(KIND=dp), DIMENSION(3)              :: A, B

    dmax=0._dp
    coa=0
    DO ipgf = 1, npgfa
     cob=0
     DO jpgf = 1, npgfb
 
       xa = zeta(ipgf)   ! exponents
       xb = zetb(jpgf)
 
       A = ra  !positions
       B = rb
       
       CALL init_os_overlap2(xa,xb,A,B)

       DO ma=la_min,la_max
         DO mb=lb_min,lb_max
           DO iax=0,ma
             DO iay=0,ma-iax
                iaz=ma-iax-iay
                na(1)=iax; na(2)=iay; na(3)=iaz
                ia1=coset(iax,iay,iaz)    
                 DO ibx=0,mb
                   DO iby=0,mb-ibx
                     ibz=mb-ibx-iby
                     nb(1)=ibx; nb(2)=iby; nb(3)=ibz
                     ib1=coset(ibx,iby,ibz)
                     res1=os_overlap2(na,nb)
                    ! write(*,*) "la, lb,na, nb, res1", ma, mb, na, nb, res1
                    ! write(*,*) "sab ia1, ib1", ia1, ib1, sab(coa+ia1,cob+ib1)
                     res2=sab(coa+ia1,cob+ib1)
                     dmax=MAX(dmax,ABS(res1-res2))
                     IF(dmax > 1.E-13) WRITE(*,*) "dmax in loop", dmax
                   END DO
                 END DO
             END DO
           END DO
         END DO
       END DO
       cob = cob + ncoset(lb_max)
     ENDDO
     coa = coa + ncoset(la_max)
    ENDDO
    WRITE(*,*) "dmax", dmax
    
  END SUBROUTINE overlap2_test_2

END MODULE lri_environment_methods
