!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates integral matrices for lripaw method
!>        lri : local resolution of the identity
!> \author JGH August 2012
! *****************************************************************************
MODULE lri_environment_methods

  USE ai_overlap3,                     ONLY: overlap3
  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: deallocate_lri_ints,&
                                             lri_density_release,&
                                             lri_density_type,&
                                             lri_environment_type,&
                                             lri_int_list,&
                                             lri_int_type,&
                                             lri_rhoab_type
  USE mathlib,                         ONLY: invmat
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: ncoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix_simple
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_environment_methods'

  PUBLIC :: build_lri_matrices, calculate_lri_densities, calculate_lri_potential

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an lri_env
!> \param lri_env the lri_environment you want to create
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_lri_matrices(lri_env,qs_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_lri_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nder
    REAL(KIND=dp)                            :: eps_filter

    nder = 0
    eps_filter = qs_env%dft_control%qs_control%eps_filter_matrix

    ! calculate the integrals needed to do the local (2-center) expansion
    ! of the (pair) densities
    CALL calculate_lri_integrals(lri_env,qs_env,error)
    ! build the overlap matrix between the RI basis and the PAW projectors
    ! this is needed for the PAW expansion of the fitted (RI) density
    CALL build_overlap_matrix_simple(qs_env,lri_env%ri_basis_soft,lri_env%proj_basis,&
       nder,lri_env%rip_matrix,"RI-P Overlap",lri_env%sip_list,eps_filter,error)

  END SUBROUTINE build_lri_matrices

! *****************************************************************************

  SUBROUTINE calculate_lri_integrals(lri_env,qs_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iac, ikind, ilist, &
                                                jkind, jneighbor, nba, nbb, &
                                                nfa, nfb, nkind, nlist, nn, &
                                                nneighbor, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, ffi
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: n
    REAL(KIND=dp), DIMENSION(3)              :: rab
    TYPE(gto_basis_set_type), POINTER        :: fbasa, fbasb, obasa, obasb
    TYPE(lri_int_list), POINTER              :: lri_ints
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      CALL get_qs_env(qs_env=qs_env,nkind=nkind,error=error)

      IF ( ASSOCIATED(lri_env%lri_ints) ) THEN
        CALL deallocate_lri_ints (lri_env%lri_ints,error)
      END IF
      NULLIFY(lri_ints)
      ALLOCATE(lri_ints,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_ints%nkind = nkind
      ALLOCATE(lri_ints%lri_atom(nkind*nkind),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i=1,nkind*nkind
        NULLIFY(lri_ints%lri_atom(i)%lri_node)
        lri_ints%lri_atom(i)%natom = 0
      END DO

      !calculate the overlap integrals <a|b|ai> and <a|b|bi>
      CALL neighbor_list_iterator_create(nl_iterator,soo_list)
      DO WHILE (neighbor_list_iterate(nl_iterator)==0)
         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
              nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,r=rab)

         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         obasa => lri_env%orb_basis(ikind)%gto_basis_set
         obasb => lri_env%orb_basis(jkind)%gto_basis_set
         fbasa => lri_env%ri_basis(ikind)%gto_basis_set
         fbasb => lri_env%ri_basis(jkind)%gto_basis_set

         IF (.NOT.ASSOCIATED(obasa)) CYCLE
         IF (.NOT.ASSOCIATED(obasb)) CYCLE

         IF(.NOT.ASSOCIATED(lri_ints%lri_atom(iac)%lri_node)) THEN
            lri_ints%lri_atom(iac)%natom = nlist
            ALLOCATE(lri_ints%lri_atom(iac)%lri_node(nlist),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            DO i=1,nlist
               NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%lri_int)
               NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%lri_rhoab)
               NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%lri_paw)
               NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%acoef)
               NULLIFY(lri_ints%lri_atom(iac)%lri_node(i)%fcoef)
               lri_ints%lri_atom(iac)%lri_node(i)%nnode = 0
            END DO
         END IF
         IF(.NOT.ASSOCIATED(lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int)) THEN
            lri_ints%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
            ALLOCATE(lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(nneighbor),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END IF

         lrii => lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)

         nba = obasa%nsgf
         nbb = obasb%nsgf
         nfa = fbasa%nsgf
         nfb = fbasb%nsgf
         ALLOCATE(lrii%abaint(nba,nbb,nfa),&
                  lrii%abbint(nba,nbb,nfb),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         lrii%abaint = 0._dp
         lrii%abbint = 0._dp

         lrii%nba = nba
         lrii%nbb = nbb
         lrii%nfa = nfa
         lrii%nfb = nfb

         ! calculate integrals (a,b,fa) and (a,b,fb)
         CALL lri_int_ab(lrii%abaint,lrii%abbint,rab,obasa,obasb,fbasa,fbasb,error)

         nn = nfa+nfb
         ALLOCATE(lrii%sab(nfa,nfb),lrii%sinv(nn,nn),lrii%sn(nn),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         lrii%sab = 0._dp

         ! calculate integrals (fa,fb)
         CALL lri_int_oo(lrii%sab,rab,fbasa,fbasb,error)
         ! invert S matrix 
         lrii%sinv(1:nfa,1:nfa) = lri_env%bas_ovlp(ikind)%ovlp(1:nfa,1:nfa)
         lrii%sinv(1:nfa,nfa+1:nn) = lrii%sab(1:nfa,1:nfb)
         lrii%sinv(nfa+1:nn,1:nfa) = TRANSPOSE(lrii%sab(1:nfa,1:nfb))
         lrii%sinv(nfa+1:nn,nfa+1:nn) = lri_env%bas_ovlp(jkind)%ovlp(1:nfb,1:nfb)
         CALL invmat(lrii%sinv,stat,error)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ! calculate Sinv*n and n*Sinv*n
         ALLOCATE(n(nn),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         n(1:nfa)    = lri_env%bas_norm(ikind)%norm(1:nfa)
         n(nfa+1:nn) = lri_env%bas_norm(jkind)%norm(1:nfb)
         lrii%sn(1:nn) = MATMUL(lrii%sinv(1:nn,1:nn),n(1:nn))
         lrii%nsn = SUM(lrii%sn(1:nn)*n(1:nn))
         DEALLOCATE(n,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ! calculate qab = -SUM_i (a,b,ai)*sn(i)/nsn
         ALLOCATE(lrii%qab(nba,nbb),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         lrii%qab(nba,nbb) = 0._dp
         DO i=1,nfa
            ffi = lrii%sn(i)/lrii%nsn
            lrii%qab(1:nba,1:nbb) = lrii%qab(1:nba,1:nbb) - ffi * lrii%abaint(1:nba,1:nbb,i)
         END DO
         DO i=1,nfb
            ffi = lrii%sn(i+nfa)/lrii%nsn
            lrii%qab(1:nba,1:nbb) = lrii%qab(1:nba,1:nbb) - ffi * lrii%abbint(1:nba,1:nbb,i)
         END DO

         ! calculate integrals (a,b), overlap of primery basis
         ALLOCATE(lrii%soo(nba,nbb),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         lrii%soo = 0._dp
         CALL lri_int_oo(lrii%soo,rab,obasa,obasb,error)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      lri_env%lri_ints => lri_ints

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_integrals

! *****************************************************************************

  SUBROUTINE calculate_lri_densities(lri_env,lri_density,pmatrix,&
                                     atomic_kind_set,para_env,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_densities', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, blk, handle, i, iac, iatom, ikind, ilist, &
      ispin, jatom, jkind, jneighbor, nat, natom, nba, nbb, nfa, nfb, nkind, &
      nlist, nn, nneighbor, nspin, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp)                            :: nst
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: n
    REAL(KIND=dp), DIMENSION(:), POINTER     :: aci, acj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pbij
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: pmat
    TYPE(gto_basis_set_type), POINTER        :: fbas, ribas
    TYPE(lri_int_list), POINTER              :: lri_coef, lri_paw, lri_rhos
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nspin = SIZE(pmatrix)
      nkind = lri_env%lri_ints%nkind

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=nat)
      ALLOCATE(atom_of_kind(nat),kind_of(nat),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of,atom_of_kind=atom_of_kind)

      CALL lri_density_release(lri_density,error)

      ALLOCATE(lri_density,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      lri_density%nspin = nspin
      ALLOCATE(lri_density%lri_rhos(nspin),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO ispin = 1, nspin
         pmat => pmatrix(ispin)%matrix

         ALLOCATE(lri_rhos,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         lri_rhos%nkind = nkind
         ALLOCATE(lri_rhos%lri_atom(nkind*nkind),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
         DO i=1,nkind*nkind
           NULLIFY(lri_rhos%lri_atom(i)%lri_node)
           lri_rhos%lri_atom(i)%natom = 0
         END DO
   
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,&
                 nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)
   
            iac = ikind + nkind*(jkind - 1)

            IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

            IF(.NOT.ASSOCIATED(lri_rhos%lri_atom(iac)%lri_node)) THEN
               lri_rhos%lri_atom(iac)%natom = nlist
               ALLOCATE(lri_rhos%lri_atom(iac)%lri_node(nlist),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               DO i=1,nlist
                  NULLIFY(lri_rhos%lri_atom(iac)%lri_node(i)%lri_int)
                  NULLIFY(lri_rhos%lri_atom(iac)%lri_node(i)%lri_rhoab)
                  NULLIFY(lri_rhos%lri_atom(iac)%lri_node(i)%lri_paw)
                  NULLIFY(lri_rhos%lri_atom(iac)%lri_node(i)%acoef)
                  NULLIFY(lri_rhos%lri_atom(iac)%lri_node(i)%fcoef)
                  lri_rhos%lri_atom(iac)%lri_node(i)%nnode = 0
               END DO
            END IF
            IF(.NOT.ASSOCIATED(lri_rhos%lri_atom(iac)%lri_node(ilist)%lri_rhoab)) THEN
               lri_rhos%lri_atom(iac)%lri_node(ilist)%nnode = nneighbor
               ALLOCATE(lri_rhos%lri_atom(iac)%lri_node(ilist)%lri_rhoab(nneighbor),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END IF

            NULLIFY (pbij)
            IF (iatom <= jatom) THEN
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=iatom,col=jatom,block=pbij,found=found)
              trans = .FALSE.
            ELSE
              CALL cp_dbcsr_get_block_p(matrix=pmat,row=jatom,col=iatom,block=pbij,found=found)
              trans = .TRUE.
            END IF
            CPPostcondition(found,cp_failure_level,routineP,error,failure)

            lrho => lri_rhos%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
   
            lrho%nba = lrii%nba
            lrho%nbb = lrii%nbb
            lrho%nfa = lrii%nfa
            lrho%nfb = lrii%nfb
            nba = lrii%nba
            nbb = lrii%nbb
            nfa = lrii%nfa
            nfb = lrii%nfb

            nn = lrho%nfa + lrho%nfb
            ALLOCATE(lrho%avec(nn),&
                     lrho%tvec(nn),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

            IF (trans) THEN
               lrho%charge = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%soo(1:nba,1:nbb))
               DO i=1,nfa
                  lrho%tvec(i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abaint(1:nba,1:nbb,i))
               END DO
               DO i=1,nfb
                  lrho%tvec(nfa+i) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%abbint(1:nba,1:nbb,i))
               END DO
            ELSE
               lrho%charge = SUM(pbij(1:nba,1:nbb)*lrii%soo(1:nba,1:nbb))
               DO i=1,nfa
                  lrho%tvec(i) = SUM(pbij(1:nba,1:nbb)*lrii%abaint(1:nba,1:nbb,i))
               END DO
               DO i=1,nfb
                  lrho%tvec(nfa+i) = SUM(pbij(1:nba,1:nbb)*lrii%abbint(1:nba,1:nbb,i))
               END DO
            END IF
            nst = SUM(lrho%tvec(1:nn) * lrii%sn(1:nn))
            lrho%lambda = (lrho%charge - nst)/lrii%nsn

            ALLOCATE(n(nn),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            n(1:nfa)    = lri_env%bas_norm(ikind)%norm(1:nfa)
            n(nfa+1:nn) = lri_env%bas_norm(jkind)%norm(1:nfb)
            n(1:nn) = lrho%tvec(1:nn) + lrho%lambda * n(1:nn)
            lrho%avec(1:nn) = MATMUL(lrii%sinv(1:nn,1:nn),n(1:nn))
            DEALLOCATE(n,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         lri_density%lri_rhos(ispin)%int_list => lri_rhos

      END DO

      ! allocate the arrays to hold RI expansion coefficients and forces
      ALLOCATE(lri_density%lri_coef(nspin),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO ispin = 1, nspin

         ALLOCATE(lri_coef,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         lri_coef%nkind = nkind
         ALLOCATE(lri_coef%lri_atom(nkind),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
         DO ikind=1,nkind
           NULLIFY(lri_coef%lri_atom(ikind)%lri_node)
           atomic_kind => atomic_kind_set(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
           lri_coef%lri_atom(ikind)%natom = natom
           ALLOCATE(lri_coef%lri_atom(ikind)%lri_node(natom),STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           fbas => lri_env%ri_basis(ikind)%gto_basis_set
           nfa = fbas%nsgf
           DO iatom=1,natom
              lri_coef%lri_atom(ikind)%lri_node(iatom)%nnode = 0
              NULLIFY(lri_coef%lri_atom(ikind)%lri_node(iatom)%lri_int)
              NULLIFY(lri_coef%lri_atom(ikind)%lri_node(iatom)%lri_rhoab)
              NULLIFY(lri_coef%lri_atom(ikind)%lri_node(iatom)%lri_paw)
              ALLOCATE(lri_coef%lri_atom(ikind)%lri_node(iatom)%acoef(nfa),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_coef%lri_atom(ikind)%lri_node(iatom)%acoef = 0._dp
              ALLOCATE(lri_coef%lri_atom(ikind)%lri_node(iatom)%fcoef(nfa),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_coef%lri_atom(ikind)%lri_node(iatom)%fcoef = 0._dp
           END DO
         END DO

         ! sum up expansion coefficients
         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,&
                 ilist=ilist,inode=jneighbor)

            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            aci => lri_coef%lri_atom(ikind)%lri_node(atom_a)%acoef
            acj => lri_coef%lri_atom(jkind)%lri_node(atom_b)%acoef
            iac = ikind + nkind*(jkind - 1)
            lrho => lri_rhos%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            nfa = lrho%nfa
            nfb = lrho%nfb
            aci(1:nfa) = aci(1:nfa) + lrho%avec(1:nfa)
            acj(1:nfb) = acj(1:nfb) + lrho%avec(nfa+1:nfa+nfb)

         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         ! replicate the acoef infomation 
         DO ikind=1,nkind
           atomic_kind => atomic_kind_set(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
           DO iatom=1,natom
              aci => lri_coef%lri_atom(ikind)%lri_node(iatom)%acoef
              CALL mp_sum(aci,para_env%group)
           END DO
         END DO

         lri_density%lri_coef(ispin)%int_list => lri_coef

      END DO

      ! calculate PAW expansion of RI density
      ALLOCATE(lri_density%lri_paw(nspin),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO ispin = 1, nspin

         ALLOCATE(lri_paw,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         lri_paw%nkind = nkind
         ALLOCATE(lri_paw%lri_atom(nkind),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
         DO ikind=1,nkind
           NULLIFY(lri_paw%lri_atom(ikind)%lri_node)
           atomic_kind => atomic_kind_set(ikind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
           lri_paw%lri_atom(ikind)%natom = natom
           ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(natom),STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           fbas => lri_env%proj_basis(ikind)%gto_basis_set
           nfa = fbas%nsgf
           ribas => lri_env%ri_basis(ikind)%gto_basis_set
           nfb = ribas%nsgf
           DO iatom=1,natom
              lri_paw%lri_atom(ikind)%lri_node(iatom)%nnode = 0
              NULLIFY(lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_int)
              NULLIFY(lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_rhoab)
              ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(iatom)%acoef(nfa),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_paw%lri_atom(ikind)%lri_node(iatom)%acoef = 0._dp
              ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(iatom)%fcoef(nfa),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_paw%lri_atom(ikind)%lri_node(iatom)%fcoef = 0._dp
              ! PAW force expansions
              ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%nb = nfa
              ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fhh(nfa),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fhs(nfa),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fxh(nfa),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fxs(nfa),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fhh = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fhs = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fxh = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%fxs = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%na = nfb
              ALLOCATE(lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhh(nfb),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhs(nfb),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxh(nfb),&
                       lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxs(nfb),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhh = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhs = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxh = 0._dp
              lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxs = 0._dp
           END DO
         END DO

         ! perform projection of the RI coefficients
         pmat => lri_env%rip_matrix(1)%matrix
         CALL cp_dbcsr_iterator_start(iter, pmat)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, pbij, blk)
            IF (iatom == jatom) CYCLE
            ! RI expansion vector
            ikind  = kind_of(iatom)
            atom_a = atom_of_kind(iatom)
            lri_coef => lri_density%lri_coef(ispin)%int_list
            aci => lri_coef%lri_atom(ikind)%lri_node(atom_a)%acoef
            ! PAW expansion vector
            jkind  = kind_of(jatom)
            atom_b = atom_of_kind(jatom)
            acj => lri_paw%lri_atom(jkind)%lri_node(atom_b)%acoef
            ! PAW = PAW + PROJECTOR*RI
            acj = acj + MATMUL(aci,pbij)
         ENDDO
         CALL cp_dbcsr_iterator_stop(iter)

         lri_density%lri_paw(ispin)%int_list => lri_paw

      END DO

      DEALLOCATE(atom_of_kind,kind_of,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_densities

! *****************************************************************************

  SUBROUTINE calculate_lri_potential(lri_env,lri_density,vmatrix,&
                                     atomic_kind_set,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: vmatrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_potential', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iac, iatom, ikind, ilist, ispin, &
      jatom, jkind, jneighbor, nat, nba, nbb, nfa, nfb, nkind, nlist, nn, &
      nneighbor, nspin, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp)                            :: nsf
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: fcij, fs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fci, fcj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vij
    TYPE(cp_dbcsr_type), POINTER             :: vmat
    TYPE(lri_int_list), POINTER              :: lri_coef, lri_rhos
    TYPE(lri_int_type), POINTER              :: lrii
    TYPE(lri_rhoab_type), POINTER            :: lrho
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: soo_list

    failure = .FALSE.
    CALL timeset(routineN,handle)

    ! get the contribution of one-center terms from PAW expansion
    CALL paw_potential(lri_env,lri_density,atomic_kind_set,error)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nspin = SIZE(vmatrix)
      nkind = lri_env%lri_ints%nkind

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=nat)
      ALLOCATE(atom_of_kind(nat),kind_of(nat),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of,atom_of_kind=atom_of_kind)

      DO ispin = 1, nspin
         vmat => vmatrix(ispin)%matrix
         lri_rhos => lri_density%lri_rhos(ispin)%int_list

         CALL neighbor_list_iterator_create(nl_iterator,soo_list)
         DO WHILE (neighbor_list_iterate(nl_iterator)==0)
            CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,&
                 nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor)

            NULLIFY (vij)
            IF (iatom <= jatom) THEN
              CALL cp_dbcsr_get_block_p(matrix=vmat,row=iatom,col=jatom,block=vij,found=found)
              trans = .FALSE.
            ELSE
              CALL cp_dbcsr_get_block_p(matrix=vmat,row=jatom,col=iatom,block=vij,found=found)
              trans = .TRUE.
            END IF
            CPPostcondition(found,cp_failure_level,routineP,error,failure)

            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            fci => lri_coef%lri_atom(ikind)%lri_node(atom_a)%fcoef
            fcj => lri_coef%lri_atom(jkind)%lri_node(atom_b)%fcoef
            lrho => lri_rhos%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)
            lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
            nba = lrii%nba 
            nbb = lrii%nbb 
            nfa = lrii%nfa 
            nfb = lrii%nfb 
            nn = nfa+nfb
            ALLOCATE(fcij(nn),fs(nn),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            fcij(1:nfa) = fci(1:nfa)
            fcij(nfa+1:nn) = fcj(1:nfb)
             
            nsf = SUM(lrii%sn*fcij)
            fs = MATMUL(lrii%sinv,fcij)

            IF(trans) THEN
               vij(1:nbb,1:nba) = vij(1:nbb,1:nba) + nsf*TRANSPOSE(lrii%qab(1:nba,1:nbb))
               DO i=1,nfa
                  vij(1:nbb,1:nba) = vij(1:nbb,1:nba) + fs(i)*TRANSPOSE(lrii%abaint(1:nba,1:nbb,i))
               END DO
               DO i=1,nfb
                  vij(1:nbb,1:nba) = vij(1:nbb,1:nba) + fs(nfa+i)*TRANSPOSE(lrii%abbint(1:nba,1:nbb,i))
               END DO
            ELSE
               vij(1:nba,1:nbb) = vij(1:nba,1:nbb) + nsf*lrii%qab(1:nba,1:nbb)
               DO i=1,nfa
                  vij(1:nba,1:nbb) = vij(1:nba,1:nbb) + fs(i)*lrii%abaint(1:nba,1:nbb,i)
               END DO
               DO i=1,nfb
                  vij(1:nba,1:nbb) = vij(1:nba,1:nbb) + fs(nfa+i)*lrii%abbint(1:nba,1:nbb,i)
               END DO
            END IF

            DEALLOCATE(fcij,fs,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

      END DO

      DEALLOCATE(atom_of_kind,kind_of,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_potential

! *****************************************************************************

  SUBROUTINE paw_potential(lri_env,lri_density,atomic_kind_set,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'paw_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, blk, handle, &
                                                iatom, ikind, ispin, jatom, &
                                                jkind, nat, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: fci, fhh, fhs, fxh, fxs
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vij
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: vmat
    TYPE(lri_int_list), POINTER              :: lri_coef, lri_paw

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=nat)
    ALLOCATE(atom_of_kind(nat),kind_of(nat),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of,atom_of_kind=atom_of_kind)

    DO ispin = 1, lri_density%nspin
       lri_coef => lri_density%lri_coef(ispin)%int_list
       lri_paw => lri_density%lri_paw(ispin)%int_list
       ! on-site expansion coefficients
       DO ikind = 1,lri_paw%nkind
          DO iatom = 1,lri_paw%lri_atom(ikind)%natom
             fhh => lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhh
             fhs => lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flhs
             fxh => lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxh
             fxs => lri_paw%lri_atom(ikind)%lri_node(iatom)%lri_paw%flxs
             fci => lri_coef%lri_atom(ikind)%lri_node(iatom)%fcoef
             fci = fci + fhh - fhs + fxh - fxs
          END DO
       END DO
       ! perform back-projection of the RI coefficients
       vmat => lri_env%rip_matrix(1)%matrix
       CALL cp_dbcsr_iterator_start(iter, vmat)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, vij, blk)
          IF (iatom == jatom) CYCLE
          ! RI force expansion vector
          ikind  = kind_of(iatom)
          atom_a = atom_of_kind(iatom)
          fci => lri_coef%lri_atom(ikind)%lri_node(atom_a)%fcoef
          ! PAW force expansion vectors
          jkind  = kind_of(jatom)
          atom_b = atom_of_kind(jatom)
          fhh => lri_paw%lri_atom(jkind)%lri_node(atom_b)%lri_paw%fhh
          fhs => lri_paw%lri_atom(jkind)%lri_node(atom_b)%lri_paw%fhs
          fxh => lri_paw%lri_atom(jkind)%lri_node(atom_b)%lri_paw%fxh
          fxs => lri_paw%lri_atom(jkind)%lri_node(atom_b)%lri_paw%fxs
          ! PAW = PAW + PROJECTOR*RI
          fci = fci + MATMUL(vij,fhh)
          fci = fci - MATMUL(vij,fhs)
          fci = fci + MATMUL(vij,fxh)
          fci = fci - MATMUL(vij,fxs)
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)

    END DO

    DEALLOCATE(atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE paw_potential

! *****************************************************************************

  SUBROUTINE lri_int_ab(abaint,abbint,rab,oba,obb,fba,fbb,error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abaint, abbint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba, fbb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iset, jset, kaset, kbset, m1, &
                                                m2, m3, ncoa, ncob, ncoc, &
                                                nseta, nsetb, nsetca, nsetcb, &
                                                sgfa, sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, &
      lca_max, lca_min, lcb_max, lcb_min, npgfa, npgfb, npgfca, npgfcb, &
      nsgfa, nsgfb, nsgfca, nsgfcb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfca, first_sgfcb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: saba, sabb
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfca, rpgfcb, &
                                                sphi_a, sphi_b, sphi_ca, &
                                                sphi_cb, zeta, zetb, zetca, &
                                                zetcb

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

          ! calculate integrals
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
                ! integrals
                CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                              rab,dab,rac,dac,rbc,dbc,saba)
                CALL ricontract(abaint(sgfa:m1,sgfb:m2,sgfc:m3),&
                     saba,sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset),error)
                DEALLOCATE(saba,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO

          ! calculate integrals
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
                ! integrals
                CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                              rab,dab,rac,dac,rbc,dbc,sabb)
                CALL ricontract(abbint(sgfa:m1,sgfb:m2,sgfc:m3),&
                     sabb,sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset),error)
                DEALLOCATE(sabb,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO

       END DO
    END DO

  END SUBROUTINE lri_int_ab

! *****************************************************************************

  SUBROUTINE lri_int_oo(sab,rab,fba,fbb,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: fba, fbb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_oo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iset, jset, lds, m1, m2, &
                                                maxco, maxl, ncoa, ncob, &
                                                nseta, nsetb, sgfa, sgfb, &
                                                stat, va, vb, wa, wb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: swork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

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

    CALL get_gto_basis_set(fba,maxco=va,maxl=wa)
    CALL get_gto_basis_set(fbb,maxco=vb,maxl=wb)
    maxco = MAX(va,vb)
    maxl = MAX(wa,wb)
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco),swork(lds,lds,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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

          ! calculate integrals
          CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                     lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                     rab,dab,sint,0,.FALSE.,swork,lds)

          CALL ocontract(sab(sgfa:m1,sgfb:m2),sint,sphi_a(:,sgfa:),sphi_b(:,sgfb:),&
               ncoa,ncob,nsgfa(iset),nsgfb(jset),error)
       END DO
    END DO
    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE lri_int_oo

! *****************************************************************************

  SUBROUTINE ricontract(abcint,sabc,ca,cb,cc,npa,npb,npc,nca,ncb,ncc,error)

    REAL(KIND=dp), DIMENSION(:, :, :)        :: abcint, sabc
    REAL(KIND=dp), DIMENSION(:, :)           :: ca, cb, cc
    INTEGER, INTENT(IN)                      :: npa, npb, npc, nca, ncb, ncc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ricontract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, m1, m2, m3, mca, mcb, &
                                                mcc, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cpc, cpp

    mca = SIZE(ca,1)
    mcb = SIZE(cb,1)
    mcc = SIZE(cc,1)

    m1 = SIZE(sabc,1)
    m2 = SIZE(sabc,2)
    m3 = SIZE(sabc,3)

    ALLOCATE(cpp(nca,m2,m3),cpc(nca,m2,ncc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL dgemm("T","N",nca,m2*m3,npa,1._dp,ca,mca,sabc,m1,0.0_dp,cpp,nca)
    CALL dgemm("N","N",nca*m2,ncc,npc,1._dp,cpp,nca*m2,cc,mcc,0.0_dp,cpc,nca*m2)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ic=1,ncc
       CALL dgemm("N","N",nca,ncb,npb,1._dp,cpc(1,1,ic),nca,cb,mcb,&
                  0.0_dp,abcint(1,1,ic),nca)
    END DO

    DEALLOCATE(cpc,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE ricontract

! *****************************************************************************

  SUBROUTINE ocontract(abint,sab,ca,cb,npa,npb,nca,ncb,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: abint
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab, ca, cb
    INTEGER, INTENT(IN)                      :: npa, npb, nca, ncb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ocontract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m1, m2, mca, mcb, nn, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cpp

    mca = SIZE(ca,1)
    mcb = SIZE(cb,1)

    m1 = SIZE(sab,1)
    m2 = SIZE(sab,2)

    nn = SIZE(abint,1)

    ALLOCATE(cpp(nca,m2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL dgemm("T","N",nca,m2,npa,1._dp,ca,mca,sab,m1,0.0_dp,cpp,nca)
    CALL dgemm("N","N",nca,ncb,npb,1._dp,cpp,nca,cb,mcb,0.0_dp,abint,nn)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE ocontract

! *****************************************************************************

END MODULE lri_environment_methods

