!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates forces for LRIGPW method
!>        lri : local resolution of the identity
!> \par History
!>      created Dorothea Golze [03.2014]
!> \authors Dorothea Golze  
! *****************************************************************************
MODULE lri_forces

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: &
       allocate_lri_force_components, deallocate_lri_force_components, &
       lri_density_type, lri_environment_type, lri_force_type, lri_int_type, &
       lri_kind_type, lri_list_type, lri_rhoab_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_forces'

  PUBLIC :: calculate_lri_forces

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief calculates the lri forces
!> \param lri_env ...
!> \param lri_density ...
!> \param qs_env ...
!> \param pmatrix density matrix
!> \param atomic_kind_set ...
!> \param use_virial ...
! *****************************************************************************
  SUBROUTINE calculate_lri_forces(lri_env,lri_density,qs_env,pmatrix,&
             atomic_kind_set,use_virial)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    LOGICAL                                  :: use_virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_lri_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ikind, ispin, &
                                                natom, nkind, nspin
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: v_dadr, v_dfdr
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_coef
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(atomic_kind, force, lri_coef, v_dadr, v_dfdr, virial)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN

      nkind = lri_env%lri_ints%nkind
      nspin = SIZE(pmatrix)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=natom)
      CALL get_qs_env(qs_env=qs_env,force=force,virial=virial)
 
      !calculate SUM_i integral(V*fbas_i)*davec/dR 
      CALL calculate_v_dadr(lri_env,lri_density,pmatrix,atomic_kind_set,&
                            use_virial,virial)

      DO ispin = 1, nspin

       lri_coef => lri_density%lri_coefs(ispin)%lri_kinds  
       
       DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom)
         DO iatom=1,natom
            v_dadr => lri_coef(ikind)%v_dadr(iatom,:)
            v_dfdr => lri_coef(ikind)%v_dfdr(iatom,:)

            force(ikind)%rho_lri_elec(:,iatom) = force(ikind)%rho_lri_elec(:,iatom)&
                                                 + v_dfdr(:) + v_dadr(:)

         END DO
       END DO
      ENDDO

    ENDIF
 
    CALL timestop(handle)

  END SUBROUTINE calculate_lri_forces

! *****************************************************************************
!> \brief calculates second term of derivative with respect to R, i.e.
!>        SUM_i integral(V * fbas_i)*davec/dR
!> \param lri_env ...
!> \param lri_density ...
!> \param pmatrix ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param virial ...
! *****************************************************************************
  SUBROUTINE calculate_v_dadr(lri_env,lri_density,pmatrix,atomic_kind_set,&
                              use_virial,virial)

    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: pmatrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(virial_type), POINTER               :: virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_v_dadr', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iac, iatom, ikind, ilist, ispin, &
      jatom, jkind, jneighbor, k, natom, nba, nbb, nfa, nfb, nkind, nlist, &
      nn, nneighbor, nspin, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: failure, found, trans
    REAL(KIND=dp), DIMENSION(3)              :: dcharge, dlambda, force_a, &
                                                force_b, nsdssn, nsdsst, &
                                                nsdt, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: st, v_dadra, v_dadrb, v_inta, &
                                                v_intb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: davec, dssn, dsst, dtvec, &
                                                pbij, sdssn, sdsst, sdt
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: ds
    TYPE(cp_dbcsr_type), POINTER             :: pmat
    TYPE(lri_force_type), POINTER            :: lri_force
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
    NULLIFY(atom_of_kind, lri_coef, lri_force, lrii, lri_rho, lrho, &
         nl_iterator, pbij, pmat, soo_list, v_dadra, v_dadrb, v_inta, &
         v_intb)
    NULLIFY(davec, ds, dssn, dsst, dtvec, sdssn, sdsst, sdt, st)

    IF ( ASSOCIATED(lri_env%soo_list) ) THEN
      soo_list => lri_env%soo_list

      nkind = lri_env%lri_ints%nkind
      nspin = SIZE(pmatrix)
 
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,natom=natom)
      ALLOCATE(atom_of_kind(natom),STAT=stat)
      CPASSERT(stat==0)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      DO ispin = 1, nspin

       pmat => pmatrix(ispin)%matrix
       lri_coef => lri_density%lri_coefs(ispin)%lri_kinds  
       lri_rho => lri_density%lri_rhos(ispin)%lri_list

       CALL neighbor_list_iterator_create(nl_iterator,soo_list)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,&
               jatom=jatom,nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,&
               r=rab)

          iac = ikind + nkind*(jkind - 1)

          IF(.NOT.ASSOCIATED(lri_env%lri_ints%lri_atom(iac)%lri_node)) CYCLE

          lrii => lri_env%lri_ints%lri_atom(iac)%lri_node(ilist)%lri_int(jneighbor)
          lrho => lri_rho%lri_atom(iac)%lri_node(ilist)%lri_rhoab(jneighbor)

          ! zero contribution for pairs aa or bb and periodc pairs aa' or bb'
          ! calculate forces for periodic pairs aa' and bb' only for virial
          IF(.NOT.lrii%calc_force_pair) CYCLE

          nfa = lrii%nfa
          nfb = lrii%nfb
          nba = lrii%nba
          nbb = lrii%nbb
          nn = nfa +  nfb

          CALL allocate_lri_force_components(lri_density%lri_force,nfa,nfb)
          lri_force => lri_density%lri_force
          ds        => lri_force%ds
          st        => lri_force%st
          dsst      => lri_force%dsst
          sdsst     => lri_force%sdsst
          dssn      => lri_force%dssn
          sdssn     => lri_force%sdssn
          sdt       => lri_force%sdt
          dtvec     => lri_force%dtvec
          davec     => lri_force%davec

          ! get the density matrix Pab
          NULLIFY (pbij)
          IF (iatom <= jatom) THEN
            CALL cp_dbcsr_get_block_p(matrix=pmat,row=iatom,col=jatom,block=pbij,found=found)
            trans = .FALSE.
          ELSE
            CALL cp_dbcsr_get_block_p(matrix=pmat,row=jatom,col=iatom,block=pbij,found=found)
            trans = .TRUE.
          END IF

          ! compute dtvec/dRa = SUM_ab Pab *d(a,b,x)/dRa 
          DO k = 1,3 
           IF (trans) THEN
              dcharge(k) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%dsoo(1:nba,1:nbb,k)) 
              DO i=1,nfa
                 dtvec(i,k) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%dabdaint(1:nba,1:nbb,i,k))
              END DO
              DO i=1,nfb
                 dtvec(nfa+i,k) = SUM(TRANSPOSE(pbij(1:nbb,1:nba))*lrii%dabbint(1:nba,1:nbb,i,k))
              END DO
           ELSE
              dcharge(k) = SUM(pbij(1:nba,1:nbb)*lrii%dsoo(1:nba,1:nbb,k))
              DO i=1,nfa
                 dtvec(i,k) = SUM(pbij(1:nba,1:nbb)*lrii%dabdaint(1:nba,1:nbb,i,k))
              END DO
              DO i=1,nfb
                 dtvec(nfa+i,k) = SUM(pbij(1:nba,1:nbb)*lrii%dabbint(1:nba,1:nbb,i,k))
              END DO
           END IF
          ENDDO

          ! construct derivative of S (overlap) matrix
          ds(1:nfa,1:nfa,1:3) = 0._dp 
          ds(1:nfa,nfa+1:nn,1:3) = lrii%dsab(1:nfa,1:nfb,1:3)
          DO k =1,3
           ds(nfa+1:nn,1:nfa,k) = TRANSPOSE(lrii%dsab(1:nfa,1:nfb,k))
          ENDDO
          ds(nfa+1:nn,nfa+1:nn,1:3) = 0._dp

          ! n*Sinv*dS*Sinv*t
          DO k =1,3
           st(1:nn)     = MATMUL(lrii%sinv(1:nn,1:nn),lrho%tvec(1:nn))
           dsst(1:nn,k) = MATMUL(ds(1:nn,1:nn,k),st(1:nn))
           sdsst(1:nn,k)= MATMUL(lrii%sinv(1:nn,1:nn),dsst(1:nn,k))
           nsdsst(k)    = SUM(sdsst(1:nn,k)*lrii%n(1:nn))
          ENDDO

          ! n*Sinv*dS*Sinv*n
          DO k =1,3
           dssn(1:nn,k) = MATMUL(ds(1:nn,1:nn,k),lrii%sn(1:nn))
           sdssn(1:nn,k)= MATMUL(lrii%sinv(1:nn,1:nn),dssn(1:nn,k))
           nsdssn(k)    = SUM(sdssn(1:nn,k)*lrii%n(1:nn))
          ENDDO

          ! Sinv*dtvec/dR_a 
          DO k =1,3
           sdt(1:nn,k) =  MATMUL(lrii%sinv(1:nn,1:nn),dtvec(1:nn,k))
           nsdt(k) = SUM(sdt(1:nn,k)*lrii%n(1:nn))
          ENDDO

          ! dlambda/dRa
          DO k =1,3
           dlambda(k) = (dcharge(k) + nsdsst(k) - nsdt(k))/lrii%nsn & 
                        + (lrho%charge - lrho%nst) * nsdssn(k)/(lrii%nsn*lrii%nsn)
          ENDDO
          
          ! davec/dRa, contains: davec_a/dRa (1..nfa) and davec_b/dRa (nfa+1..nn)
          DO k = 1,3 
           davec(1:nn,k) = -sdsst(1:nn,k) + sdt(1:nn,k)&
                           + dlambda(k)*lrii%sn(1:nn) - lrho%lambda * sdssn(1:nn,k)
          ENDDO

          atom_a = atom_of_kind(iatom)
          atom_b = atom_of_kind(jatom)

          v_dadra => lri_coef(ikind)%v_dadr(atom_a,:)
          v_dadrb => lri_coef(jkind)%v_dadr(atom_b,:)
          v_inta  => lri_coef(ikind)%v_int(atom_a,:) 
          v_intb  => lri_coef(jkind)%v_int(atom_b,:)
 
          ! sum over atom pairs
          DO k=1,3
           force_a(k) =   2.0_dp*SUM(davec(1:nfa,k)*v_inta(1:nfa))&
                        + 2.0_dp*SUM(davec(nfa+1:nn,k)*v_intb(1:nfb))
           v_dadra(k) = v_dadra(k) + force_a(k) 
           !minus: using that davec/dRa = - davec/dRb
           force_b(k) = - 2.0_dp*SUM(davec(nfa+1:nn,k)*v_intb(1:nfb))&
                        - 2.0_dp*SUM(davec(1:nfa,k)*v_inta(1:nfa))
           v_dadrb(k) = v_dadrb(k) + force_b(k)
          ENDDO

          ! contribution to virial
          IF(use_virial) THEN
            !periodic self-pairs aa' contribute only with factor 0.5
            IF(iatom == jatom) THEN
              CALL virial_pair_force (virial%pv_virial,0.5_dp,force_a,rab)
            ELSE 
              CALL virial_pair_force (virial%pv_virial,1.0_dp,force_a,rab)
            ENDIF
          ENDIF

          CALL deallocate_lri_force_components(lri_density%lri_force)

       END DO

       CALL neighbor_list_iterator_release(nl_iterator)

      ENDDO

      DEALLOCATE(atom_of_kind,STAT=stat)
      CPASSERT(stat==0)
 
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_v_dadr

END MODULE lri_forces
