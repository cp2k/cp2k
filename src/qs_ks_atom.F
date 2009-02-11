!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the Kohn-Sham matrix  contributions
!>      coming from local atomic densities
! *****************************************************************************
MODULE qs_ks_atom
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE sap_kind_types,                  ONLY: alist_type,&
                                             get_alist
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_atom'

  PUBLIC :: update_ks_atom

CONTAINS

! *****************************************************************************
!> \brief The correction to the KS matrix due to the GAPW local terms to the hartree and
!>      XC contributions is here added. The correspondig forces contribution are also calculated
!>      if required. Each sparse-matrix block A-B is corrected by all the atomic contributions
!>      centered on atoms C for which the triplet A-C-B exists (they are close enough)
!>      To this end special lists are used
!> \param qs_env qs enviroment, for the lists, the contraction coefficients and the
!>               pre calculated integrals of the potential with the atomic orbitals
!> \param ks KS matrix, sparse matrix
!> \param p density matrix, sparse matrix, needed only for the forces
!> \param forces switch for the calculation of the forces on atoms
!> \param tddft switch for TDDFT linear response
!> \param p_env perturbation theory environment
!> \param error cp_error_type
!> \par History
!>      created [MI]
!>      the loop over the spins is done internally [03-05,MI]
!>      Rewrite using new OCE matrices [08.02.09, jhu]        
! *****************************************************************************
  SUBROUTINE update_ks_atom(qs_env,ks,p,forces,tddft,p_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      INTENT(INOUT)                          :: ks
    TYPE(real_matrix_p_type), DIMENSION(:), &
      INTENT(IN)                             :: p
    LOGICAL, INTENT(IN)                      :: forces
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: p_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_ks_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia_kind, iab, iac, iat, iatom, ibc, ikind, ilist, &
      inode, ispin, istat, ja_kind, jatom, jkind, ka_kind, kac, katom, kbc, &
      kkind, lC_int, ldCPC, lPC, max_gau, max_nsgf, maxsoc, nat, natom, &
      nkind, nlist, nnode, nsetc, nsgfa, nsgfb, nsoctot, nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cell_b
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: my_tddft, paw_atom, use_virial
    REAL(dp), DIMENSION(3)                   :: rac, rbc
    REAL(dp), DIMENSION(3, 3)                :: force_tmp
    REAL(kind=dp)                            :: eps_cpc, factor, pmax
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dCPC_h, dCPC_s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: C_int_h, C_int_s, PC_h, PC_s
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b
    TYPE(alist_type), POINTER                :: alist, alist_ac, alist_bc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_h, mat_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom
    TYPE(rho_atom_type), POINTER             :: rho_at
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind_set,force,oce,para_env,rho_atom,sab_orb)
    NULLIFY(mat_h, mat_p ,orb_basis)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft
    factor = 1.0_dp
    IF (my_tddft) THEN
       IF (qs_env%dft_control%nspins == 1) THEN
          factor = 2.0_dp
       END IF
    END IF
    nspins = qs_env%dft_control%nspins

    CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           force=force,&
           oce=oce, &
           para_env=para_env,&
           rho_atom_set=rho_atom,&
           virial=virial,&
           sab_orb=sab_orb,error=error)

    IF (my_tddft) THEN
       rho_atom   => p_env%local_rho_set%rho_atom_set
    END IF

    eps_cpc = qs_env%dft_control%qs_control%gapw_control%eps_cpc

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)

    ALLOCATE (C_int_h(max_gau,max_nsgf,nspins),STAT=istat)
    ALLOCATE (C_int_s(max_gau,max_nsgf,nspins),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
                        "C_int_h,C_int_s",2*max_gau*max_nsgf*nspins*dp_size)
    lC_int = max_gau
    ALLOCATE(mat_h(nspins),mat_p(nspins), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
                        "mat_h,mat_p")
    DO ispin = 1,nspins
      NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
    END DO

    IF(forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
           "atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
      ALLOCATE(dCPC_h(max_gau,max_gau),STAT=istat)
      ALLOCATE(dCPC_s(max_gau,max_gau),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
            "dCPC_h,dCPC_s",2*(max_gau*max_gau)*dp_size)
      ldCPC = max_gau
      ALLOCATE(PC_h(max_nsgf,max_gau,nspins),STAT=istat)
      ALLOCATE(PC_s(max_nsgf,max_gau,nspins),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
            "PC_h,PC_s",2*(max_nsgf*max_gau*nspins)*dp_size)
      lPC = max_nsgf
      use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    END IF

    nkind = SIZE(atomic_kind_set,1)
    ! Collect the local potential contributions from all the processors
    DO ikind = 1,nkind
       NULLIFY(atomic_kind,atom_list)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,atom_list=atom_list,&
            paw_atom=paw_atom, natom=nat)
       IF(paw_atom) THEN
         DO iat = 1,nat
         iatom =  atom_list(iat)
            DO ispin = 1,nspins
              CALL  mp_sum(rho_atom(iatom)%ga_Vlocal_gb_h(ispin)%r_coef,&
                           para_env%group)
              CALL  mp_sum(rho_atom(iatom)%ga_Vlocal_gb_s(ispin)%r_coef,&
                           para_env%group)
            END DO
         END DO
       END IF
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis,&
                           paw_atom=paw_atom)
      CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                nsgf=nsgfa)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind, &
                             orb_basis_set=orb_basis)
        CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                nsgf=nsgfb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

        DO ilist=1,nlist

          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,neighbor=jatom,cell=cell_b)

            pmax = 0._dp
            DO ispin = 1,nspins
              NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
              IF (iatom <= jatom) THEN
                CALL get_block_node(matrix=ks(ispin)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    BLOCK=mat_h(ispin)%array)
              ELSE
                CALL get_block_node(matrix=ks(ispin)%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    BLOCK=mat_h(ispin)%array)
              END IF
              pmax = pmax+MAXVAL(ABS(mat_h(ispin)%array))

              IF(forces) THEN
                IF (iatom <= jatom) THEN
                  CALL get_block_node(matrix=p(ispin)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    BLOCK=mat_p(ispin)%array)
                ELSE
                  CALL get_block_node(matrix=p(ispin)%matrix,&
                                      block_row=jatom,&
                                      block_col=iatom,&
                                      BLOCK=mat_p(ispin)%array)
                END IF
              END IF
            END DO

            DO kkind=1,nkind

              atomic_kind => atomic_kind_set(kkind)

              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   orb_basis_set=orb_basis,&
                                   paw_atom=paw_atom)

              IF(.NOT. paw_atom)  CYCLE

              CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                      nset =nsetc,maxso=maxsoc)
              nsoctot = maxsoc*nsetc

              iac = ikind + nkind*(kkind - 1)
              ibc = jkind + nkind*(kkind - 1)

              IF (.NOT.ASSOCIATED(oce%intac(iac)%alist)) CYCLE
              IF (.NOT.ASSOCIATED(oce%intac(ibc)%alist)) CYCLE
              CALL get_alist(oce%intac(iac)%alist, alist_ac, iatom, error)
              CALL get_alist(oce%intac(ibc)%alist, alist_bc, jatom, error)
              DO kac=1,alist_ac%nclist
                DO kbc=1,alist_bc%nclist
                  IF ( alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom ) CYCLE
                  IF ( ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0) ) THEN
                    IF(pmax*alist_bc%clist(kbc)%maxac*alist_ac%clist(kac)%maxac < eps_cpc) CYCLE
                    katom = alist_ac%clist(kac)%catom
                    IF(iatom==katom .AND. ALL(alist_ac%clist(kac)%cell == 0)) THEN
                      C_coeff_hh_a => alist_ac%clist(kac)%achint
                      C_coeff_ss_a => alist_ac%clist(kac)%acint
                    ELSE
                      C_coeff_hh_a => alist_ac%clist(kac)%acint
                      C_coeff_ss_a => alist_ac%clist(kac)%acint
                    END IF
                    IF(jatom==katom .AND. ALL(alist_bc%clist(kbc)%cell == 0)) THEN
                      C_coeff_hh_b => alist_bc%clist(kbc)%achint
                      C_coeff_ss_b => alist_bc%clist(kbc)%acint
                    ELSE
                      C_coeff_hh_b => alist_bc%clist(kbc)%acint
                      C_coeff_ss_b => alist_bc%clist(kbc)%acint
                    END IF

                    rho_at => rho_atom(katom)

                    CALL add_vhxca_to_ks2(mat_h,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,C_coeff_ss_b,&
                         rho_at,nspins,iatom,nsgfa,jatom,nsgfb,katom,nsoctot,factor,C_int_h,C_int_s)

                    IF(forces)THEN
                      ia_kind = atom_of_kind(iatom)
                      ja_kind = atom_of_kind(jatom)
                      ka_kind = atom_of_kind(katom)
                      rho_at => rho_atom(katom)
                      force_tmp(1:3,1:3) = 0.0_dp
                      IF(iatom <= jatom) THEN
                        CALL add_vhxca_forces2(mat_p,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,C_coeff_ss_b,&
                             rho_at,force_tmp,nspins,iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                             dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC)
                        force(kkind)%vhxc_atom(1:3,ka_kind) = &
                                          force(kkind)%vhxc_atom(1:3,ka_kind) + &
                                          force_tmp(1:3,3)
                        force(ikind)%vhxc_atom(1:3,ia_kind) = &
                                         force(ikind)%vhxc_atom(1:3,ia_kind) + &
                                         force_tmp(1:3,1)
                        force(jkind)%vhxc_atom(1:3,ja_kind) = &
                                         force(jkind)%vhxc_atom(1:3,ja_kind) + &
                                         force_tmp(1:3,2)
                        IF (use_virial) THEN
                          rac = alist_ac%clist(kac)%rac
                          rbc = alist_bc%clist(kbc)%rac
                          CALL virial_pair_force(virial%pv_virial,1._dp,force_tmp(1:3,1),rac,error)
                          CALL virial_pair_force(virial%pv_virial,1._dp,force_tmp(1:3,2),rbc,error)
                        END IF
                      ELSE
                        CALL add_vhxca_forces2(mat_p,C_coeff_hh_b,C_coeff_hh_a,C_coeff_ss_b,C_coeff_ss_a,&
                             rho_at,force_tmp,nspins,jatom,nsgfb,iatom,nsgfa,katom,nsoctot,&
                             dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC)
                        force(kkind)%vhxc_atom(1:3,ka_kind) = &
                                          force(kkind)%vhxc_atom(1:3,ka_kind) + &
                                          force_tmp(1:3,3)
                        force(ikind)%vhxc_atom(1:3,ia_kind) = &
                                          force(ikind)%vhxc_atom(1:3,ia_kind) + &
                                          force_tmp(1:3,2)
                        force(jkind)%vhxc_atom(1:3,ja_kind) = &
                                          force(jkind)%vhxc_atom(1:3,ja_kind) + &
                                          force_tmp(1:3,1)
                        IF (use_virial) THEN
                          rac = alist_ac%clist(kac)%rac
                          rbc = alist_bc%clist(kbc)%rac
                          CALL virial_pair_force(virial%pv_virial,1._dp,force_tmp(1:3,2),rac,error)
                          CALL virial_pair_force(virial%pv_virial,1._dp,force_tmp(1:3,1),rbc,error)
                        END IF
                      END IF
                    END IF
                    EXIT ! search loop over jatom-katom list
                  END IF
                END DO
              END DO

            ENDDO    ! kkind

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO    ! inode
        ENDDO    ! ilist
      ENDDO    ! jkind
    ENDDO    ! ikind

    DEALLOCATE (C_int_h,C_int_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
                              "C_int_h,C_int_s")

    DO ispin = 1,nspins
      NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
    END DO
    DEALLOCATE(mat_h,mat_p, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
                              "mat_h,mat_p")

    IF(forces) THEN
       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
         "atom_of_kind")
       DEALLOCATE(DCPC_h,dCPC_s,PC_h,PC_s,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
         "DCPC_h,dCPC_s,PC_h,PC_s")
    END IF

    CALL timestop(handle)

  END SUBROUTINE update_ks_atom

! *****************************************************************************
  SUBROUTINE add_vhxca_to_ks2(mat_h,C_hh_a,C_hh_b,C_ss_a,C_ss_b,&
             rho_atom,nspins,ia,na,ja,nb,ka,nsp,factor,C_int_h,C_int_s)

    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_h
    REAL(KIND=dp), DIMENSION(:, :, :)        :: C_hh_a, C_hh_b, C_ss_a, C_ss_b
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      :: nspins, ia, na, ja, nb, ka, &
                                                nsp
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(:, :, :)        :: C_int_h, C_int_s

    CHARACTER(len=*), PARAMETER :: routineN = 'add_vhxca_to_ks2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, int_hard, int_soft
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s

    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,&
                      ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)

    DO ispin = 1,nspins
      int_hard => int_local_h(ispin)%r_coef
      int_soft => int_local_s(ispin)%r_coef
      h_block => mat_h(ispin)%array
      IF (ia <= ja) THEN
        CALL DGEMM('N','T',nsp,nb,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                   0.0_dp,C_int_h(1,1,ispin),SIZE(C_int_h,1))
        ! C_int_h(1:nsp,1:nb,ispin) = MATMUL(int_hard(1:nsp,1:nsp),TRANSPOSE(C_hh_b(1:nb,1:nsp,1)))
        CALL DGEMM('N','T',nsp,nb,nsp,1.0_dp,int_soft(1,1),SIZE(int_soft,1),C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                   0.0_dp,C_int_s(1,1,ispin),SIZE(C_int_s,1))
        ! C_int_s(1:nsp,1:nb,ispin) = MATMUL(int_soft(1:nsp,1:nsp),TRANSPOSE(C_ss_b(1:nb,1:nsp,1)))

        CALL DGEMM('N','N',na,nb,nsp, factor,C_hh_a(1,1,1),SIZE(C_hh_a,1),C_int_h(1,1,ispin),SIZE(C_int_h,1),&
                   1.0_dp,h_block(1,1),SIZE(h_block,1))
        ! h_block(1:na,1:nb) = h_block(1:na,1:nb) + factor*MATMUL(C_hh_a(1:na,1:nsp,1),C_int_h(1:nsp,1:nb,ispin))
        CALL DGEMM('N','N',na,nb,nsp,-factor,C_ss_a(1,1,1),SIZE(C_ss_a,1),C_int_s(1,1,ispin),SIZE(C_int_s,1),&
                   1.0_dp,h_block(1,1),SIZE(h_block,1))
        ! h_block(1:na,1:nb) = h_block(1:na,1:nb) - factor*MATMUL(C_ss_a(1:na,1:nsp,1),C_int_s(1:nsp,1:nb,ispin))
      ELSE
        int_hard => int_local_h(ispin)%r_coef
        int_soft => int_local_s(ispin)%r_coef

        CALL DGEMM('N','T',nsp,na,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                   0.0_dp,C_int_h(1,1,ispin),SIZE(C_int_h,1))
        ! C_int_h(1:nsp,1:na,ispin) = MATMUL(int_hard(1:nsp,1:nsp),TRANSPOSE(C_hh_a(1:na,1:nsp,1)))
        CALL DGEMM('N','T',nsp,na,nsp,1.0_dp,int_soft(1,1),SIZE(int_soft,1),C_ss_a(1,1,1),SIZE(C_ss_a,1),&
                   0.0_dp,C_int_s(1,1,ispin),SIZE(C_int_s,1))
        ! C_int_s(1:nsp,1:na,ispin) = MATMUL(int_soft(1:nsp,1:nsp),TRANSPOSE(C_ss_a(1:na,1:nsp,1)))

        CALL DGEMM('N','N',nb,na,nsp, factor,C_hh_b(1,1,1),SIZE(C_hh_b,1),C_int_h(1,1,ispin),SIZE(C_int_h,1),&
                   1.0_dp,h_block(1,1),SIZE(h_block,1))
        ! h_block(1:nb,1:na) = h_block(1:nb,1:na) + factor*MATMUL(C_hh_b(1:nb,1:nsp,1),C_int_h(1:nsp,1:na,ispin))
        CALL DGEMM('N','N',nb,na,nsp,-factor,C_ss_b(1,1,1),SIZE(C_ss_b,1),C_int_s(1,1,ispin),SIZE(C_int_s,1),&
                   1.0_dp,h_block(1,1),SIZE(h_block,1))
        ! h_block(1:nb,1:na) = h_block(1:nb,1:na) - factor*MATMUL(C_ss_b(1:nb,1:nsp,1),C_int_s(1:nsp,1:na,ispin))
      END IF
    END DO

  END SUBROUTINE add_vhxca_to_ks2

! *****************************************************************************
  SUBROUTINE add_vhxca_forces2(mat_p,C_hh_a,C_hh_b,C_ss_a,C_ss_b,&
             rho_atom,force,nspins,ia,na,ja,nb,ka,nsp,dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC)

    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_p
    REAL(KIND=dp), DIMENSION(:, :, :)        :: C_hh_a, C_hh_b, C_ss_a, C_ss_b
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(3, 3)                :: force
    INTEGER, INTENT(IN)                      :: nspins, ia, na, ja, nb, ka, &
                                                nsp
    REAL(KIND=dp), DIMENSION(:, :)           :: dCPC_h, dCPC_s
    INTEGER, INTENT(IN)                      :: ldCPC
    REAL(KIND=dp), DIMENSION(:, :, :)        :: PC_h, PC_s
    INTEGER, INTENT(IN)                      :: lPC

    CHARACTER(len=*), PARAMETER :: routineN = 'add_vhxca_forces2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, ispin
    REAL(dp), DIMENSION(:, :), POINTER       :: int_hard, int_soft
    REAL(KIND=dp)                            :: ieqj, trace
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s

!   *** factor 2 because only half of the pairs with ia =/ ja are  considered

    ieqj = 2.0_dp
    IF(ia == ja) ieqj = 1.0_dp

    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)

    DO ispin = 1,nspins
      p_block => mat_p(ispin)%array
      int_hard => int_local_h(ispin)%r_coef
      int_soft => int_local_s(ispin)%r_coef

      CALL DGEMM('N','N',na,nsp,nb,1.0_dp,p_block(1,1),SIZE(p_block,1),C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 0.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1))
      ! PC_h(1:na,1:nsp,ispin) = MATMUL(p_block(1:na,1:nb),C_hh_b(1:nb,1:nsp,1))
      CALL DGEMM('N','N',na,nsp,nb,1.0_dp,p_block(1,1),SIZE(p_block,1),C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                 0.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1))
      ! PC_s(1:na,1:nsp,ispin) = MATMUL(p_block(1:na,1:nb),C_ss_b(1:nb,1:nsp,1))

      DO dir = 2,4
        CALL DGEMM('T','N',nsp,nsp,na,1.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1),C_hh_a(1,1,dir),SIZE(C_hh_a,1),&
                   0.0_dp,dCPC_h(1,1),SIZE(dCPC_h,1))
        ! dCPC_h(1:nsp,1:nsp) = MATMUL(TRANSPOSE(PC_h(1:na,1:nsp,ispin)),C_hh_a(1:na,1:nsp,dir))
        trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)
        force(dir-1,3) = force(dir-1,3) + ieqj * trace
        force(dir-1,1) = force(dir-1,1) - ieqj * trace

        CALL DGEMM('T','N',nsp,nsp,na,1.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1),C_ss_a(1,1,dir),SIZE(C_ss_a,1),&
                   0.0_dp,dCPC_s(1,1),SIZE(dCPC_s,1))
        ! dCPC_s(1:nsp,1:nsp) = MATMUL(TRANSPOSE(PC_s(1:na,1:nsp,ispin)),C_ss_a(1:na,1:nsp,dir))
        trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)
        force(dir-1,3) = force(dir-1,3) -  ieqj * trace
        force(dir-1,1) = force(dir-1,1) +  ieqj * trace
      END DO

      ! j-k contributions

      CALL DGEMM('T','N',nb,nsp,na,1.0_dp,p_block(1,1),SIZE(p_block,1),C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                   0.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1))
      ! PC_h(1:nb,1:nsp,ispin) = MATMUL(TRANSPOSE(p_block(1:na,1:nb)),C_hh_a(1:na,1:nsp,1))
      CALL DGEMM('T','N',nb,nsp,na,1.0_dp,p_block(1,1),SIZE(p_block,1),C_ss_a(1,1,1),SIZE(C_ss_a,1),&
                   0.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1))
      ! PC_s(1:nb,1:nsp,ispin) = MATMUL(TRANSPOSE(p_block(1:na,1:nb)),C_ss_a(1:na,1:nsp,1))

      DO dir = 2,4
        CALL DGEMM('T','N',nsp,nsp,nb,1.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1),C_hh_b(1,1,dir),SIZE(C_hh_b,1),&
                   0.0_dp,dCPC_h(1,1),SIZE(dCPC_h,1))
        ! dCPC_h(1:nsp,1:nsp) = MATMUL(TRANSPOSE(PC_h(1:nb,1:nsp,ispin)),C_hh_b(1:nb,1:nsp,dir))
        trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)
        force(dir-1,3) = force(dir-1,3) +  ieqj * trace
        force(dir-1,2) = force(dir-1,2) -  ieqj * trace

        CALL DGEMM('T','N',nsp,nsp,nb,1.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1),C_ss_b(1,1,dir),SIZE(C_ss_b,1),&
                   0.0_dp,dCPC_s(1,1),SIZE(dCPC_s,1))
        ! dCPC_s(1:nsp,1:nsp) = MATMUL(TRANSPOSE(PC_s(1:nb,1:nsp,ispin)),C_ss_b(1:nb,1:nsp,dir))
        trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)
        force(dir-1,3) = force(dir-1,3) -  ieqj * trace
        force(dir-1,2) = force(dir-1,2) +  ieqj * trace
      END DO
    END DO

  END SUBROUTINE add_vhxca_forces2

END MODULE qs_ks_atom
