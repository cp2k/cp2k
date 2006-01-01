!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-1004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_ks_atom [1.0] *
!!
!!   NAME
!!     qs_ks_atom
!!
!!   FUNCTION
!!     routines that build the Kohn-Sham matrix  contributions
!!     coming from local atomic densities
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!
!!
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_atom
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: tddfpt_singlet
  USE kinds,                           ONLY: dp, dp_size, int_size
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next,&
                                             reduced_3c_list_type
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_ks_atom'

  PUBLIC :: update_ks_atom

!******************************************************************************

CONTAINS

!******************************************************************************
!!
!!
!!
!!   FUNCTION
!!     The correction to the KS matrix due to the GAPW local terms to the hartree and
!!     XC contributions is here added. The correspondig forces contribution are also calculated
!!     if required. Each sparse-matrix block A-B is corrected by all the atomic contributions
!!     centered on atoms C for which the triplet A-C-B exists (they are close enough)
!!     To this end special lists are used
!!
!!   INPUTS
!!     qs_env : qs enviroment, for the lists, the contraction coefficients and the 
!!              pre calculated integrals of the potential with the atomic orbitals
!!      ks : KS matrix, sparse matrix
!!      p : density matrix, sparse matrix, needed only for the forces
!!      forces : switch for the calculation of the forces on atoms
!!      tddft : switch for TDDFT linear response
!!      p_env : perturbation theory environment
!!      error : cp_error_type
!!
!!   HISTOTY
!!     created [MI]
!!     the loop over the soins is done internally [03-05,MI]
!!
!!
!!
!******************************************************************************

  SUBROUTINE update_ks_atom(qs_env,ks,p,forces,tddft,p_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:),&
      INTENT(INOUT)                          :: ks
    TYPE(real_matrix_p_type), DIMENSION(:),&
      INTENT(IN)                             :: p
    LOGICAL, INTENT(IN)                      :: forces
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL   :: p_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine = "update_ks_atom"

    INTEGER :: handle, i, ia_kind, iab, iat, iatom, ikind, ilist, inode, &
               ispin, istat, ja_kind, jatom, jkind, ka_kind, katom, &
               kkind, knode, last_jatom, last_katom, lC_int, ldCPC, lPC,&
               max_gau, max_nsgf, &
               maxsoc, nat, natom, nkind, nlist, nnode, nnode_c, &
               nsetc, nsgfa, nsgfb, nsoctot, nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  ::  my_tddft, &
                                                 new_atom_b, new_atom_c,&
                                                 paw_atom
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:,:,:)                       :: C_int_h, C_int_s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dCPC_h, dCPC_s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:,:,:)                       :: PC_h, PC_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block
    REAL(dp), DIMENSION(3,3)                 :: force_tmp
    REAL(kind=dp)                            :: factor
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_2d_r_p_type), DIMENSION(:),&
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
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom
    TYPE(rho_atom_type), POINTER             :: rho_at

!******************************************************************************

    CALL timeset("update_ks_atom","I"," ",handle)

    NULLIFY(atomic_kind_set,force,oce,para_env,rho_atom,sab_orb)
    NULLIFY(mat_h, mat_p ,orb_basis, reduced_3c_oce)

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
           sab_orb=sab_orb)

    IF (my_tddft) THEN
       rho_atom   => p_env%local_rho_set%rho_atom_set
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)

    ALLOCATE (C_int_h(max_gau,max_nsgf,nspins),STAT=istat)
    ALLOCATE (C_int_s(max_gau,max_nsgf,nspins),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                        "C_int_h,C_int_s",2*max_gau*max_nsgf*nspins*dp_size)
    lC_int = max_gau
    ALLOCATE(mat_h(nspins),mat_p(nspins), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                        "mat_h,mat_p")
    DO ispin = 1,nspins
      NULLIFY(mat_h(ispin)%array,mat_p(ispin)%array)
    END DO


    IF(forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
           "atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
      ALLOCATE(dCPC_h(max_gau,max_gau),STAT=istat)
      ALLOCATE(dCPC_s(max_gau,max_gau),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
            "dCPC_h,dCPC_s",2*(max_gau*max_gau)*dp_size)
      ldCPC = max_gau
      ALLOCATE(PC_h(max_nsgf,max_gau,nspins),STAT=istat)
      ALLOCATE(PC_s(max_nsgf,max_gau,nspins),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
            "PC_h,PC_s",2*(max_nsgf*max_gau*nspins)*dp_size)
      lPC = max_nsgf
    END IF

    nkind = SIZE(atomic_kind_set,1)
    ! Collect the local potentialcontributions from all the processors
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

          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   reduced_3c_oce=reduced_3c_oce)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

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

              last_katom = 0

              nnode_c = reduced_3c_oce(kkind)%nnode
              DO knode=1,nnode_c

                katom = reduced_3c_oce(kkind)%index_atom(knode)

                IF (katom /= last_katom) THEN
                   new_atom_c = .TRUE.
                   last_katom = katom
                ELSE
                   new_atom_c = .FALSE.
                   CYCLE
                END IF

                rho_at => rho_atom(katom)
                CALL add_vhxca_to_ks(mat_h,oce,rho_at,nspins, &
                                    iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                                    factor,C_int_h,C_int_s)

                IF(forces)THEN
                  ia_kind = atom_of_kind(iatom)
                  ja_kind = atom_of_kind(jatom)
                  ka_kind = atom_of_kind(katom)
                  IF(iatom <= jatom) THEN
                    rho_at => rho_atom(katom)
                    force_tmp(1:3,1:3) = 0.0_dp 
                    CALL add_vhxca_forces(mat_p,oce,rho_at,force_tmp,&
                                    nspins,iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
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
                  ELSE
                    rho_at => rho_atom(katom)
                    force_tmp(1:3,1:3) = 0.0_dp 
                    CALL add_vhxca_forces(mat_p,oce,rho_at,force_tmp,&
                                    nspins,jatom,nsgfb,iatom,nsgfa,katom,nsoctot,&
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
                  END IF
                END IF

              END DO

            ENDDO    ! kkind

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO    ! inode
        ENDDO    ! ilist 
      ENDDO    ! jkind
    ENDDO    ! ikind

    DEALLOCATE (C_int_h,C_int_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                              "C_int_h,C_int_s")

    DEALLOCATE(mat_h,mat_p, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                              "mat_h,mat_p")

    IF(forces) THEN
       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
         "atom_of_kind")
       DEALLOCATE(DCPC_h,dCPC_s,PC_h,PC_s,STAT=istat)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE update_ks_atom
!******************************************************************************

  SUBROUTINE add_vhxca_to_ks(mat_h,oce,rho_atom,nspins,ia,na,ja,nb,&
                             ka,nsp,factor,C_int_h,C_int_s)

!   Arguments
    TYPE(cp_2d_r_p_type), DIMENSION(:),&
      POINTER                                :: mat_h
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      :: nspins, ia, na, ja, nb, ka, &
                                                nsp
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(:,:,:)          :: C_int_h,C_int_s

!   Parameters
    CHARACTER(LEN=*), PARAMETER              :: routine = "add_vhxca_to_ks"

!   Local Variables
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_ss , C_hh , h_block
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: int_hard, int_soft
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h,int_local_s
    INTEGER                                  :: is1, is2, ispin, istat

!******************************************************************************


    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,&
                      ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)


    IF (ia <= ja) THEN
      NULLIFY(C_hh)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_hh)

      NULLIFY(C_ss)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_ss)

      DO ispin = 1,nspins
        int_hard => int_local_h(ispin)%r_coef
        int_soft => int_local_s(ispin)%r_coef
        CALL dgemm("N","N",nsp,nb,nsp,&
                 1.0_dp,int_hard,nsp,&
                 C_hh,SIZE(C_hh,1),0.0_dp,&
                 C_int_h(1,1,ispin),SIZE(C_int_h,1))

        CALL dgemm("N","N",nsp,nb,nsp,&
                 1.0_dp,int_soft,nsp,&
                 C_ss,SIZE(C_ss,1),0.0_dp,&
                 C_int_s(1,1,ispin),SIZE(C_int_s,1))
      END DO

      NULLIFY(C_hh)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_hh)

      NULLIFY(C_ss)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_ss)

      DO ispin = 1,nspins
        h_block => mat_h(ispin)%array
        CALL dgemm("T","N",na,nb,nsp,&
                 factor,C_hh,SIZE(C_hh,1),&
                 C_int_h(1,1,ispin),SIZE(C_int_h,1),1.0_dp,&
                 h_block,SIZE(h_block,1))
       
        CALL dgemm("T","N",na,nb,nsp,&
                 -factor,C_ss,SIZE(C_ss,1),&
                 C_int_s(1,1,ispin),SIZE(C_int_s,1),1.0_dp,&
                 h_block,SIZE(h_block,1))
      END DO

    ELSE

      NULLIFY(C_hh)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_hh)

      NULLIFY(C_ss)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_ss)

      DO ispin = 1,nspins
        int_hard => int_local_h(ispin)%r_coef
        int_soft => int_local_s(ispin)%r_coef
        CALL dgemm("N","N",nsp,na,nsp,&
                 1.0_dp,int_hard,nsp,&
                 C_hh,SIZE(C_hh,1),0.0_dp,&
                 C_int_h(1,1,ispin),SIZE(C_int_h,1))

        CALL dgemm("N","N",nsp,na,nsp,&
                 1.0_dp,int_soft,nsp,&
                 C_ss,SIZE(C_ss,1),0.0_dp,&
                 C_int_s(1,1,ispin),SIZE(C_int_s,1))
      END DO

      NULLIFY(C_hh)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_hh)

      NULLIFY(C_ss)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_ss)

      DO ispin = 1,nspins
        h_block => mat_h(ispin)%array
 
        CALL dgemm("T","N",nb,na,nsp,&
                 factor,C_hh,SIZE(C_hh,1),&
                 C_int_h(1,1,ispin),SIZE(C_int_h,1),1.0_dp,&
                 h_block,SIZE(h_block,1))

        CALL dgemm("T","N",nb,na,nsp,&
                 -factor,C_ss,SIZE(C_ss,1),&
                 C_int_s(1,1,ispin),SIZE(C_int_s,1),1.0_dp,&
                 h_block,SIZE(h_block,1))
      END DO

    END IF

  END SUBROUTINE add_vhxca_to_ks

! *****************************************************************************

  SUBROUTINE add_vhxca_forces(mat_p,oce,rho_atom,force,nspins,&
                            ia,na,ja,nb,ka,nsp,dCPC_h,dCPC_s,ldCPC,PC_h,PC_s,lPC)

! *****************************************************************************

    TYPE(cp_2d_r_p_type), DIMENSION(:),&
      POINTER                                :: mat_p
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(3,3)                 :: force
    INTEGER, INTENT(IN)                      :: nspins, ia, na, ja, nb, ka, &
                                                nsp
    INTEGER, INTENT(IN)                      :: ldCPC, lPC
    REAL(KIND=dp), DIMENSION(:, :)           :: dCPC_h, dCPC_s
    REAL(KIND=dp), DIMENSION(:,:, :)         :: PC_h, PC_s

    CHARACTER(LEN=*), PARAMETER              :: routine = "add_vhxca_forces"

    INTEGER                                  :: dir, ispin
    REAL(KIND=dp)                            :: ieqj, trace
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_hh, C_ss
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h,int_local_s
    REAL(dp), DIMENSION(:,:), POINTER        :: int_hard, int_soft
    INTEGER                                  :: ig, jg, is1, is2, istat

! *****************************************************************************

!   *** factor 2 because only half of the pairs with ia =/ ja are  considered
    ieqj = 2.0_dp
    IF(ia == ja) ieqj = 1.0_dp

    NULLIFY(int_local_h,int_local_s,int_hard,int_soft)
    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)



    NULLIFY(C_hh)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_hh)

    NULLIFY(C_ss)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_ss)

    DO ispin = 1,nspins
      p_block => mat_p(ispin)%array
      CALL dgemm("N","T",na,nsp,nb,&
                 1.0_dp,p_block,SIZE(p_block,1),&
                 C_hh,SIZE(C_hh,1),0.0_dp,&
                 PC_h(1,1,ispin),na)

      CALL dgemm("N","T",na,nsp,nb,&
                 1.0_dp,p_block,SIZE(p_block,1),&
                 C_ss,SIZE(C_ss,1),0.0_dp,&
                 PC_s(1,1,ispin),lPC)
    END DO

    DO dir = 2,4
       NULLIFY(C_hh)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_hh)

       NULLIFY(C_ss)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_ss)


       DO ispin = 1,nspins
         int_hard => int_local_h(ispin)%r_coef
         int_soft => int_local_s(ispin)%r_coef

         CALL dgemm("T","T",nsp,nsp,na,&
                    1.0_dp,&
                    PC_h(1,1,ispin),na,C_hh,SIZE(C_hh,1),0.0_dp,&
                    dCPC_h,ldCPC)


         trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)

         force(dir-1,3) = force(dir-1,3) + ieqj * trace
         force(dir-1,1) = force(dir-1,1) - ieqj * trace


         CALL dgemm("T","T",nsp,nsp,na,&
                    1.0_dp,PC_s(1,1,ispin),lPC,&
                    C_ss,SIZE(C_ss,1),0.0_dp,&
                    dCPC_s,ldCPC)

         trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)

         force(dir-1,3) = force(dir-1,3) -  ieqj * trace
         force(dir-1,1) = force(dir-1,1) +  ieqj * trace
      END DO

    END DO

! j-k contributions


    NULLIFY(C_hh)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_hh)

    NULLIFY(C_ss)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=C_ss)

    DO ispin = 1, nspins
      p_block => mat_p(ispin)%array
      CALL dgemm("T","T",nb,nsp,na,&
                 1.0_dp,p_block,SIZE(p_block,1),&
                 C_hh,SIZE(C_hh,1),0.0_dp,&
                 PC_h(1,1,ispin),nb)

      CALL dgemm("T","T",nb,nsp,na,&
                 1.0_dp,p_block,SIZE(p_block,1),&
                 C_ss,SIZE(C_ss,1),0.0_dp,&
                 PC_s(1,1,ispin),lPC)
    END DO

    DO dir = 2,4
       NULLIFY(C_hh)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_hh)

       NULLIFY(C_ss)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=C_ss)

      DO ispin = 1,nspins
         int_hard => int_local_h(ispin)%r_coef
         int_soft => int_local_s(ispin)%r_coef

         CALL dgemm("T","T",nsp,nsp,nb,&
                    1.0_dp,PC_h(1,1,ispin),nb,&
                    C_hh,SIZE(C_hh,1),0.0_dp,&
                    dCPC_h,ldCPC)


         trace = trace_r_AxB(dCPC_h,ldCPC,int_hard,nsp,nsp,nsp)

         force(dir-1,3) = force(dir-1,3) +  ieqj * trace
         force(dir-1,2) = force(dir-1,2) -  ieqj * trace 


         CALL dgemm("T","T",nsp,nsp,nb,&
                    1.0_dp,PC_s(1,1,ispin),lPC,&
                    C_ss,SIZE(C_ss,1),0.0_dp,&
                    dCPC_s,ldCPC)

         trace = trace_r_AxB(dCPC_s,ldCPC,int_soft,nsp,nsp,nsp)

         force(dir-1,3) = force(dir-1,3) -  ieqj * trace
         force(dir-1,2) = force(dir-1,2) +  ieqj * trace
       END DO

    END DO

  END SUBROUTINE add_vhxca_forces

! *****************************************************************************

END MODULE qs_ks_atom
