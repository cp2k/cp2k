!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the Kohn-Sham matrix  contributions
!>      coming from local atomic densities
! *****************************************************************************
MODULE qs_ks_atom
  USE ao_util,                         ONLY: trace_r_AxB
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_oce_methods,                  ONLY: prj_gather
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE sap_kind_types,                  ONLY: alist_post_align_blk,&
                                             alist_pre_align_blk,&
                                             alist_type,&
                                             get_alist
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      INTENT(INOUT)                          :: ks, p
    LOGICAL, INTENT(IN)                      :: forces
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: p_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_ks_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, ia_kind, iac, iat, iatom, ibc, ikind, ip, &
      ispin, istat, ja_kind, jatom, jkind, ka_kind, kac, katom, kbc, kkind, &
      lC_int, ldCPC, lPC, max_gau, max_nsgf, maxsoc, mepos, n_cont_a, &
      n_cont_b, nat, natom, nkind, nsetc, nsgfa, nsgfb, nsoctot, nspins, &
      num_pe
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cell_b
    INTEGER, DIMENSION(:), POINTER           :: atom_list, list_a, list_b
    LOGICAL                                  :: dista, distb, failure, found, &
                                                my_tddft, paw_atom, use_virial
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a_matrix, p_matrix
    REAL(dp), DIMENSION(3)                   :: rac, rbc
    REAL(dp), DIMENSION(3, 3)                :: force_tmp
    REAL(kind=dp)                            :: eps_cpc, factor
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: C_int_h, C_int_s, coc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dCPC_h, dCPC_s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: PC_h, PC_s
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b
    TYPE(alist_type), POINTER                :: alist_ac, alist_bc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_h, mat_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom
    TYPE(rho_atom_type), POINTER             :: rho_at
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(qs_kind_set,atomic_kind_set,force,oce,para_env,rho_atom,sab_orb)
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
           qs_kind_set=qs_kind_set,&
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

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set,maxsgf=max_nsgf,maxgtops=max_gau)

    ALLOCATE (C_int_h(max_gau*max_nsgf),C_int_s(max_gau*max_nsgf),coc(max_gau*max_gau),&
              a_matrix(max_gau,max_gau),p_matrix(max_nsgf,max_nsgf),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

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
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
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

    nkind = SIZE(qs_kind_set,1)
    ! Collect the local potential contributions from all the processors
    mepos = para_env%mepos
    num_pe = para_env%num_pe
    DO ikind = 1,nkind
       NULLIFY(atom_list)
       CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=atom_list, natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
       IF(paw_atom) THEN
         ! gather the atomic block integrals in a more compressed format
         bo = get_limit( nat, num_pe, mepos )
         DO iat = bo(1), bo(2)
            iatom = atom_list(iat)
            DO ispin = 1,nspins
              CALL prj_gather(rho_atom(iatom)%ga_Vlocal_gb_h(ispin)%r_coef,&
                              rho_atom(iatom)%cpc_h(ispin)%r_coef,qs_kind_set(ikind),error)
              CALL prj_gather(rho_atom(iatom)%ga_Vlocal_gb_s(ispin)%r_coef,&
                              rho_atom(iatom)%cpc_s(ispin)%r_coef,qs_kind_set(ikind),error)
            END DO
         END DO
         ! broadcast the CPC arrays to all processors (replicated data)
         DO ip=0,num_pe-1
           bo = get_limit( nat, num_pe, ip )
           DO iat = bo(1), bo(2)
             iatom = atom_list(iat)
             DO ispin = 1,nspins
               CALL mp_bcast(rho_atom(iatom)%cpc_h(ispin)%r_coef,ip,para_env%group)
               CALL mp_bcast(rho_atom(iatom)%cpc_s(ispin)%r_coef,ip,para_env%group)
             END DO
           END DO
         END DO
       END IF
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,moduleN,__LINE__,&
           "basis_set_list",nkind)
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
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,cell=cell_b)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       nsgfa = basis_set_a%nsgf
       nsgfb = basis_set_b%nsgf

       DO ispin = 1,nspins
         NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
         IF (iatom <= jatom) THEN
           CALL cp_dbcsr_get_block_p(matrix=ks(ispin)%matrix,&
                row=iatom,col=jatom,BLOCK=mat_h(ispin)%array,found=found)
         ELSE
           CALL cp_dbcsr_get_block_p(matrix=ks(ispin)%matrix,&
                row=jatom,col=iatom,BLOCK=mat_h(ispin)%array,found=found)
         END IF

         IF(forces) THEN
           IF (iatom <= jatom) THEN
             CALL cp_dbcsr_get_block_p(matrix=p(ispin)%matrix,&
                  row=iatom,col=jatom,BLOCK=mat_p(ispin)%array,found=found)
           ELSE
             CALL cp_dbcsr_get_block_p(matrix=p(ispin)%matrix,&
                  row=jatom,col=iatom,BLOCK=mat_p(ispin)%array,found=found)
           END IF
         END IF
       END DO

       DO kkind=1,nkind
         CALL get_qs_kind(qs_kind_set(kkind), orb_basis_set=orb_basis, paw_atom=paw_atom)

         IF(.NOT. paw_atom)  CYCLE

         CALL  get_gto_basis_set(gto_basis_set=orb_basis,nset=nsetc,maxso=maxsoc)

         iac = ikind + nkind*(kkind - 1)
         ibc = jkind + nkind*(kkind - 1)

         IF (.NOT.ASSOCIATED(oce%intac(iac)%alist)) CYCLE
         IF (.NOT.ASSOCIATED(oce%intac(ibc)%alist)) CYCLE
         CALL get_alist(oce%intac(iac), alist_ac, iatom, error)
         CALL get_alist(oce%intac(ibc), alist_bc, jatom, error)
         IF (.NOT.ASSOCIATED(alist_ac)) CYCLE
         IF (.NOT.ASSOCIATED(alist_bc)) CYCLE

         DO kac=1,alist_ac%nclist
           DO kbc=1,alist_bc%nclist
             IF ( alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom ) CYCLE
             IF ( ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0) ) THEN

               n_cont_a = alist_ac%clist(kac)%nsgf_cnt
               n_cont_b = alist_bc%clist(kbc)%nsgf_cnt
               IF(n_cont_a.EQ.0.OR.n_cont_b.EQ.0) CYCLE
               list_a => alist_ac%clist(kac)%sgf_list
               list_b => alist_bc%clist(kbc)%sgf_list

               katom = alist_ac%clist(kac)%catom
               IF(iatom==katom .AND. ALL(alist_ac%clist(kac)%cell == 0)) THEN
                 C_coeff_hh_a => alist_ac%clist(kac)%achint
                 C_coeff_ss_a => alist_ac%clist(kac)%acint
                 dista = .FALSE.
               ELSE
                 C_coeff_hh_a => alist_ac%clist(kac)%acint
                 C_coeff_ss_a => alist_ac%clist(kac)%acint
                 dista = .TRUE.
               END IF
               IF(jatom==katom .AND. ALL(alist_bc%clist(kbc)%cell == 0)) THEN
                 C_coeff_hh_b => alist_bc%clist(kbc)%achint
                 C_coeff_ss_b => alist_bc%clist(kbc)%acint
                 distb = .FALSE.
               ELSE
                 C_coeff_hh_b => alist_bc%clist(kbc)%acint
                 C_coeff_ss_b => alist_bc%clist(kbc)%acint
                 distb = .TRUE.
               END IF

               rho_at => rho_atom(katom)

               nsoctot = SIZE(C_coeff_ss_a,2)

               CALL add_vhxca_to_ks(mat_h,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,C_coeff_ss_b,&
                    rho_at,nspins,iatom,jatom,nsoctot,factor,&
                    list_a,n_cont_a,list_b,n_cont_b,C_int_h,C_int_s,a_matrix,dista,distb,coc,error)

               IF(forces)THEN
                 ia_kind = atom_of_kind(iatom)
                 ja_kind = atom_of_kind(jatom)
                 ka_kind = atom_of_kind(katom)
                 rho_at => rho_atom(katom)
                 force_tmp(1:3,1:3) = 0.0_dp
                 IF(iatom <= jatom) THEN
                    CALL add_vhxca_forces(mat_p,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,C_coeff_ss_b,&
                         rho_at,force_tmp,nspins,iatom,jatom,nsoctot,&
                         list_a,n_cont_a,list_b,n_cont_b,dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC,p_matrix,&
                         dista,distb,coc,error)

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
                    CALL add_vhxca_forces(mat_p,C_coeff_hh_b,C_coeff_hh_a,C_coeff_ss_b,C_coeff_ss_a,&
                         rho_at,force_tmp,nspins,jatom,iatom,nsoctot,&
                         list_b,n_cont_b,list_a,n_cont_a,dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC,p_matrix,&
                         distb,dista,coc,error)

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
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DO ispin = 1,nspins
      NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
    END DO

    DEALLOCATE(basis_set_list,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mat_h,mat_p,C_int_h,C_int_s,a_matrix,p_matrix,coc,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

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
!> \brief ...
!> \param mat_h ...
!> \param C_hh_a ...
!> \param C_hh_b ...
!> \param C_ss_a ...
!> \param C_ss_b ...
!> \param rho_atom ...
!> \param nspins ...
!> \param ia ...
!> \param ja ...
!> \param nsp ...
!> \param factor ...
!> \param lista ...
!> \param nconta ...
!> \param listb ...
!> \param ncontb ...
!> \param C_int_h ...
!> \param C_int_s ...
!> \param a_matrix ...
!> \param dista ...
!> \param distb ...
!> \param coc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE add_vhxca_to_ks(mat_h,C_hh_a,C_hh_b,C_ss_a,C_ss_b,&
             rho_atom,nspins,ia,ja,nsp,factor,lista,nconta,listb,ncontb,&
             C_int_h,C_int_s,a_matrix,dista,distb,coc,error)

    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_h
    REAL(KIND=dp), DIMENSION(:, :, :)        :: C_hh_a, C_hh_b, C_ss_a, C_ss_b
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      :: nspins, ia, ja, nsp
    REAL(KIND=dp), INTENT(IN)                :: factor
    INTEGER, DIMENSION(:), INTENT(IN)        :: lista
    INTEGER, INTENT(IN)                      :: nconta
    INTEGER, DIMENSION(:), INTENT(IN)        :: listb
    INTEGER, INTENT(IN)                      :: ncontb
    REAL(KIND=dp), DIMENSION(:)              :: C_int_h, C_int_s
    REAL(dp), DIMENSION(:, :)                :: a_matrix
    LOGICAL, INTENT(IN)                      :: dista, distb
    REAL(dp), DIMENSION(:)                   :: coc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_vhxca_to_ks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, j, k
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, int_hard, int_soft
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s

    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,cpc_h=int_local_h,cpc_s=int_local_s)

    DO ispin = 1,nspins
       h_block => mat_h(ispin)%array
       !
       int_hard => int_local_h(ispin)%r_coef
       int_soft => int_local_s(ispin)%r_coef
       !
       IF(ia <= ja) THEN
          IF (dista .AND. distb) THEN
            k=0
            DO i=1,nsp
              DO j=1,nsp
                k = k+1
                coc(k) = int_hard(j,i)-int_soft(j,i)
              END DO
            END DO
            CALL DGEMM('N','T',nsp,ncontb,nsp,1.0_dp,coc,nsp,C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 0.0_dp,C_int_h(1),nsp)
            CALL DGEMM('N','N',nconta,ncontb,nsp,factor,C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 C_int_h(1),nsp,0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSEIF (dista) THEN
            CALL DGEMM('N','T',nsp,ncontb,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),&
                 C_hh_b(1,1,1),SIZE(C_hh_b,1),0.0_dp,coc(1),nsp)
            CALL DGEMM('N','T',nsp,ncontb,nsp,-1.0_dp,int_soft(1,1),SIZE(int_soft,1),&
                 C_ss_b(1,1,1),SIZE(C_ss_b,1),1.0_dp,coc(1),nsp)
            CALL DGEMM('N','N',nconta,ncontb,nsp,factor,C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 coc(1),nsp,0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSEIF (distb) THEN
            CALL DGEMM('N','N',nconta,nsp,nsp,factor,C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 int_hard(1,1),SIZE(int_hard,1),0.0_dp,coc(1),nconta)
            CALL DGEMM('N','N',nconta,nsp,nsp,-factor,C_ss_a(1,1,1),SIZE(C_hh_a,1),&
                 int_soft(1,1),SIZE(int_hard,1),1.0_dp,coc(1),nconta)
            CALL DGEMM('N','T',nconta,ncontb,nsp,1.0_dp,coc(1),nconta,&
                 C_hh_b(1,1,1),SIZE(C_hh_b,1),0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSE
            CALL DGEMM('N','T',nsp,ncontb,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),&
                 C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 0.0_dp,C_int_h(1),nsp)
            CALL DGEMM('N','T',nsp,ncontb,nsp,1.0_dp,int_soft(1,1),SIZE(int_soft,1),&
                 C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                 0.0_dp,C_int_s(1),nsp)
            CALL DGEMM('N','N',nconta,ncontb,nsp,factor,C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 C_int_h(1),nsp,&
                 0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
            CALL DGEMM('N','N',nconta,ncontb,nsp,-factor,C_ss_a(1,1,1),SIZE(C_ss_a,1),&
                 C_int_s(1),nsp,&
                 1.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          END IF
          !
          CALL alist_post_align_blk(a_matrix,SIZE(a_matrix,1),h_block,SIZE(h_block,1),&
               lista,nconta,listb,ncontb)
       ELSE
          IF (dista .AND. distb) THEN
            k=0
            DO i=1,nsp
              DO j=1,nsp
                k = k+1
                coc(k) = int_hard(j,i)-int_soft(j,i)
              END DO
            END DO
            CALL DGEMM('N','T',nsp,nconta,nsp,1.0_dp,coc,nsp,C_hh_a(1,1,1),SIZE(C_hh_a,1),0.0_dp,C_int_h(1),nsp)
            CALL DGEMM('N','N',ncontb,nconta,nsp, factor,C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 C_int_h(1),nsp,0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSEIF (dista) THEN
            CALL DGEMM('N','N',ncontb,nsp,nsp,factor,C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 int_hard(1,1),SIZE(int_hard,1),0.0_dp,coc,ncontb)
            CALL DGEMM('N','N',ncontb,nsp,nsp,-factor,C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                 int_soft(1,1),SIZE(int_soft,1),1.0_dp,coc,ncontb)
            CALL DGEMM('N','T',ncontb,nconta,nsp,1.0_dp,coc,ncontb,&
                 C_hh_a(1,1,1),SIZE(C_hh_a,1),0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSEIF (distb) THEN
            CALL DGEMM('N','T',nsp,nconta,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),&
                 C_hh_a(1,1,1),SIZE(C_hh_a,1),0.0_dp,coc(1),nsp)
            CALL DGEMM('N','T',nsp,nconta,nsp,-1.0_dp,int_soft(1,1),SIZE(int_soft,1),&
                 C_ss_a(1,1,1),SIZE(C_ss_a,1),1.0_dp,coc(1),nsp)
            CALL DGEMM('N','N',ncontb,nconta,nsp, factor,C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 coc(1),nsp,0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          ELSE
            CALL DGEMM('N','T',nsp,nconta,nsp,1.0_dp,int_hard(1,1),SIZE(int_hard,1),&
                 C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 0.0_dp,C_int_h(1),nsp)
            CALL DGEMM('N','T',nsp,nconta,nsp,1.0_dp,int_soft(1,1),SIZE(int_soft,1),&
                 C_ss_a(1,1,1),SIZE(C_ss_a,1),&
                 0.0_dp,C_int_s(1),nsp)
            CALL DGEMM('N','N',ncontb,nconta,nsp, factor,C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 C_int_h(1),nsp,&
                 0.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
            CALL DGEMM('N','N',ncontb,nconta,nsp,-factor,C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                 C_int_s(1),nsp,&
                 1.0_dp,a_matrix(1,1),SIZE(a_matrix,1))
          END IF
          !
          CALL alist_post_align_blk(a_matrix,SIZE(a_matrix,1),h_block,SIZE(h_block,1),&
               listb,ncontb,lista,nconta)
       END IF
    END DO

  END SUBROUTINE add_vhxca_to_ks

! *****************************************************************************
!> \brief ...
!> \param mat_p ...
!> \param C_hh_a ...
!> \param C_hh_b ...
!> \param C_ss_a ...
!> \param C_ss_b ...
!> \param rho_atom ...
!> \param force ...
!> \param nspins ...
!> \param ia ...
!> \param ja ...
!> \param nsp ...
!> \param lista ...
!> \param nconta ...
!> \param listb ...
!> \param ncontb ...
!> \param dCPC_h ...
!> \param dCPC_s ...
!> \param ldCPC ...
!> \param PC_h ...
!> \param PC_s ...
!> \param lPC ...
!> \param p_matrix ...
!> \param dista ...
!> \param distb ...
!> \param coc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE add_vhxca_forces(mat_p,C_hh_a,C_hh_b,C_ss_a,C_ss_b,&
             rho_atom,force,nspins,ia,ja,nsp,lista,nconta,listb,ncontb,&
             dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC,p_matrix,dista,distb,coc,error)

    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: mat_p
    REAL(KIND=dp), DIMENSION(:, :, :)        :: C_hh_a, C_hh_b, C_ss_a, C_ss_b
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(3, 3)                :: force
    INTEGER, INTENT(IN)                      :: nspins, ia, ja, nsp
    INTEGER, DIMENSION(:), INTENT(IN)        :: lista
    INTEGER, INTENT(IN)                      :: nconta
    INTEGER, DIMENSION(:), INTENT(IN)        :: listb
    INTEGER, INTENT(IN)                      :: ncontb
    REAL(KIND=dp), DIMENSION(:, :)           :: dCPC_h, dCPC_s
    INTEGER, INTENT(IN)                      :: ldCPC
    REAL(KIND=dp), DIMENSION(:, :, :)        :: PC_h, PC_s
    INTEGER, INTENT(IN)                      :: lPC
    REAL(KIND=dp), DIMENSION(:, :)           :: p_matrix
    LOGICAL, INTENT(IN)                      :: dista, distb
    REAL(dp), DIMENSION(:)                   :: coc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_vhxca_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, ispin
    REAL(dp), DIMENSION(:, :), POINTER       :: int_hard, int_soft
    REAL(KIND=dp)                            :: ieqj, trace
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s

!   if(dista.and.distb) we could also here use ChPCh = CsPCs
!   *** factor 2 because only half of the pairs with ia =/ ja are  considered

    ieqj = 2.0_dp
    IF(ia == ja) ieqj = 1.0_dp

    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,cpc_h=int_local_h,cpc_s=int_local_s)

    DO ispin = 1,nspins
      p_block => mat_p(ispin)%array

      CALL alist_pre_align_blk(p_block,SIZE(p_block,1),p_matrix,SIZE(p_matrix,1),&
           lista,nconta,listb,ncontb)

      int_hard => int_local_h(ispin)%r_coef
      int_soft => int_local_s(ispin)%r_coef

      CALL DGEMM('N','N',nconta,nsp,ncontb,1.0_dp,p_matrix(1,1),SIZE(p_matrix,1),&
                 C_hh_b(1,1,1),SIZE(C_hh_b,1),&
                 0.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1))
      CALL DGEMM('N','N',nconta,nsp,ncontb,1.0_dp,p_matrix(1,1),SIZE(p_matrix,1),&
                 C_ss_b(1,1,1),SIZE(C_ss_b,1),&
                 0.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1))

      DO dir = 2,4
         CALL DGEMM('T','N',nsp,nsp,nconta,1.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1),&
                    C_hh_a(1,1,dir),SIZE(C_hh_a,1),&
                    0.0_dp,dCPC_h(1,1),SIZE(dCPC_h,1))
         trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)
         force(dir-1,3) = force(dir-1,3) + ieqj * trace
         force(dir-1,1) = force(dir-1,1) - ieqj * trace

         CALL DGEMM('T','N',nsp,nsp,nconta,1.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1),&
                    C_ss_a(1,1,dir),SIZE(C_ss_a,1),&
                    0.0_dp,dCPC_s(1,1),SIZE(dCPC_s,1))
         trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)
         force(dir-1,3) = force(dir-1,3) -  ieqj * trace
         force(dir-1,1) = force(dir-1,1) +  ieqj * trace
      END DO

      ! j-k contributions

      CALL DGEMM('T','N',ncontb,nsp,nconta,1.0_dp,p_matrix(1,1),SIZE(p_matrix,1),&
                 C_hh_a(1,1,1),SIZE(C_hh_a,1),&
                 0.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1))
      CALL DGEMM('T','N',ncontb,nsp,nconta,1.0_dp,p_matrix(1,1),SIZE(p_matrix,1),&
                 C_ss_a(1,1,1),SIZE(C_ss_a,1),&
                 0.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1))

      DO dir = 2,4
         CALL DGEMM('T','N',nsp,nsp,ncontb,1.0_dp,PC_h(1,1,ispin),SIZE(PC_h,1),&
                    C_hh_b(1,1,dir),SIZE(C_hh_b,1),&
                    0.0_dp,dCPC_h(1,1),SIZE(dCPC_h,1))
         trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)
         force(dir-1,3) = force(dir-1,3) +  ieqj * trace
         force(dir-1,2) = force(dir-1,2) -  ieqj * trace

         CALL DGEMM('T','N',nsp,nsp,ncontb,1.0_dp,PC_s(1,1,ispin),SIZE(PC_s,1),&
                    C_ss_b(1,1,dir),SIZE(C_ss_b,1),&
                    0.0_dp,dCPC_s(1,1),SIZE(dCPC_s,1))
         trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)
         force(dir-1,3) = force(dir-1,3) -  ieqj * trace
         force(dir-1,2) = force(dir-1,2) +  ieqj * trace
      END DO
    END DO

  END SUBROUTINE add_vhxca_forces

END MODULE qs_ks_atom
