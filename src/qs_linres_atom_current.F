!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief given the response wavefunctions obtained by the application
!>      of the (rxp), p, and ((dk-dl)xp) operators,
!>      here the current density vector (jx, jy, jz)
!>      is computed for the 3 directions of the magnetic field (Bx, By, Bz)
!> \par History
!>      created 02-2006 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_atom_current

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: current_gauge_atom,&
                                             current_gauge_r,&
                                             current_gauge_r_and_step_func
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_linres_op,                    ONLY: fac_vecp,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: allocate_jrho_atom_rad,&
                                             allocate_jrho_coeff,&
                                             current_env_type,&
                                             get_current_env,&
                                             jrho_atom_type,&
                                             set2zero_jrho_atom_rad
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next, reduced_3c_list_type
  USE qs_oce_methods,                  ONLY: set_up_oce
  USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                             allocate_vtriple,&
                                             build_reduced_3c_lists,&
                                             create_oce_set,&
                                             dist_list_4oce_type,&
                                             oce_matrix_type,&
                                             prepare_oce_coeff,&
                                             retrieve_sac_list
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: list_3c_jrho_atom, calculate_jrho_atom_coeff,calculate_jrho_atom_rad,&
       &    calculate_jrho_atom

  ! *** Global parameters (only in this module)
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_atom_current'

CONTAINS

! *****************************************************************************
!> \brief Calculate the expansion coefficients for the atomic terms
!>      of the current densitiy in GAPW
!> \par History
!>      07.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE calculate_jrho_atom_coeff(qs_env,current_env,mat_d0,mat_jp,mat_d1,mat_jp_rii,&
       &                               mat_jp_riii,iB,idir,error)
    !
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(current_env_type)                   :: current_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: mat_d0, mat_jp
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: mat_d1
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: mat_jp_rii, mat_jp_riii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_atom_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, iab, iat, iatom, idir2, ii, iii, ikind, ilist, &
      inode, ispin, istat, jatom, jkind, katom, kkind, kneighbor, last_jatom, &
      last_katom, max_gau, max_nsgf, maxso, maxsoc, mepos, nat, natom, nkind, &
      nlist, nnode, nnode_c, nset, nsetc, nsgfa, nsgfb, nsoctot, nspins, &
      num_pe
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, is_not_associated, &
                                                new_atom_b, new_atom_c, &
                                                paw_atom
    REAL(dp)                                 :: hard_radius_a, hard_radius_b, &
                                                hard_radius_c, rab(3), rbc(3)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: PC
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b, &
                                                tmp_coeff, zero_coeff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: mat_a, mat_b, mat_c, mat_d
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: a_block, b_block, c_block, &
                                                d_block, jp2_RARnu, jp_RARnu

    CALL timeset(routineN,handle)
    !
    failure =.FALSE.
    NULLIFY(atom_list,atomic_kind_set,dft_control,orb_basis_set,sab_orb,&
         &  jrho1_atom_set,oce,reduced_3c_oce,a_block,b_block,c_block,&
         &  d_block,jp_RARnu,jp2_RARnu,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,&
         &  C_coeff_ss_b,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
         &  sab_orb_neighbor_node,para_env,zero_coeff,tmp_coeff)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          dft_control=dft_control,&
         &          oce=oce,&
         &          sab_all=sab_orb,&
         &          para_env=para_env,&
         &          error=error)
    !
    CPPrecondition(ASSOCIATED(oce),cp_failure_level,routineP,error,failure)
    !
    CALL get_current_env(current_env=current_env,&
         &           jrho1_atom_set=jrho1_atom_set,&
         &           error=error)
    !
    CPPrecondition(ASSOCIATED(jrho1_atom_set),cp_failure_level,routineP,error,failure)
    !
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         &                   maxsgf=max_nsgf,&
         &                   maxgtops=max_gau)

    idir2 = 1
    IF(idir.NE.iB) THEN
       CALL set_vecp_rev(idir,iB,idir2)
    ENDIF
    CALL set_vecp(iB,ii,iii)
    !
    !
    ! Set pointers for the different gauge
    mat_a => mat_d0
    mat_b => mat_jp
    mat_c => mat_jp_rii
    mat_d => mat_jp_riii
    !
    ! Density-like matrices
    nkind  = SIZE(atomic_kind_set)
    natom  = SIZE(jrho1_atom_set)
    nspins = dft_control%nspins
    !
    ! Allocate some arrays
    ALLOCATE(a_block(nspins),b_block(nspins),c_block(nspins),d_block(nspins),&
         &   jp_RARnu(nspins),jp2_RARnu(nspins),PC(max_nsgf,max_gau),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    ! Reset CJC coefficients and local density arrays
    DO ikind = 1 ,nkind
       NULLIFY(atomic_kind, atom_list)
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            &               atom_list=atom_list,&
            &               natom=nat,&
            &               paw_atom=paw_atom)
       !
       ! Quick cycle if needed.
       IF(.NOT. paw_atom) CYCLE
       !
       ! Initialize the density matrix-like arrays.
       DO iat = 1,nat
          iatom = atom_list(iat)
          DO ispin = 1, nspins
             IF(ASSOCIATED(jrho1_atom_set(iatom)%cjc0_h(1)%r_coef)) THEN
                jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef = 0.0_dp
             ENDIF
          ENDDO ! ispin
       ENDDO ! iat
    ENDDO ! ikind
    !
    ! Three centers
    DO ikind = 1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            &               hard_radius=hard_radius_a,&
            &               orb_basis_set=orb_basis_set)
       IF(.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            &                 nsgf=nsgfa)

       DO jkind = 1,nkind
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               &               hard_radius=hard_radius_b,&
               &               orb_basis_set=orb_basis_set)

          IF(.NOT.ASSOCIATED(orb_basis_set)) CYCLE

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
               &                 nsgf=nsgfb)

          iab = ikind + nkind*(jkind - 1)

          IF(.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

          CALL get_neighbor_list_set(neighbor_list_set=&
               &                    sab_orb(iab)%neighbor_list_set,&
               &                    nlist=nlist)

          NULLIFY(sab_orb_neighbor_list)
          !
          ! Allocate temporary block to store the sum of the 3 blocks before contraction
          DO ispin = 1,nspins
             NULLIFY(jp_RARnu(ispin)%r_coef,jp2_RARnu(ispin)%r_coef)
             ALLOCATE(jp_RARnu(ispin)%r_coef(nsgfa,nsgfb),&
                  &   jp2_RARnu(ispin)%r_coef(nsgfa,nsgfb),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDDO

          DO ilist = 1,nlist
             IF(.NOT.ASSOCIATED(sab_orb_neighbor_list)) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             ENDIF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list

             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                  &                 atom=iatom,nnode=nnode)
             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

             DO inode = 1,nnode
                CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                     &                 neighbor=jatom,r=rab,&
                     &                 reduced_3c_oce=reduced_3c_oce)

                IF(jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                   CYCLE
                ENDIF
                !
                ! Take the block \mu\nu of jpab, jpab_ii and jpab_iii
                DO ispin = 1,nspins
                   NULLIFY(a_block(ispin)%r_coef)
                   NULLIFY(b_block(ispin)%r_coef)
                   NULLIFY(c_block(ispin)%r_coef)
                   NULLIFY(d_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_a(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=a_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_b(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=b_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_c(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=c_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_d(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=d_block(ispin)%r_coef)
                ENDDO
                !
                ! Loop over atoms
                DO kkind = 1,nkind

                   atomic_kind => atomic_kind_set(kkind)

                   CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        &               orb_basis_set=orb_basis_set,&
                        &               hard_radius=hard_radius_c,&
                        &               paw_atom=paw_atom)

                   CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                        &                 nset=nsetc,maxso=maxsoc)

                   nsoctot = maxsoc*nsetc
                   !
                   ! Quick cycle if needed.
                   IF(.NOT.paw_atom) CYCLE

                   last_katom = 0

                   nnode_c = reduced_3c_oce(kkind)%nnode
                   DO kneighbor = 1,nnode_c

                      katom = reduced_3c_oce(kkind)%index_atom(kneighbor)

                      rbc(1:3) = reduced_3c_oce(kkind)%rbc(1:3,kneighbor)

                      IF (katom /= last_katom) THEN
                         new_atom_c = .TRUE.
                         last_katom = katom
                      ELSE
                         new_atom_c = .FALSE.
                         CYCLE
                      END IF
                      !
                      IF(.NOT.ASSOCIATED(jrho1_atom_set(katom)%cjc0_h(1)%r_coef)) THEN
                         CALL allocate_jrho_coeff(jrho1_atom_set,katom,nsoctot,error)
                      ENDIF
                      !
                      ! Compute the modified Qai matrix as
                      ! mQai_\mu\nu = Qai_\mu\nu - Qbi_\mu\nu * (R_A-R_\nu)_ii
                      !             + Qci_\mu\nu * (R_A-R_\nu)_iii
                      DO ispin = 1,nspins
                         CALL DCOPY(nsgfa*nsgfb,b_block(ispin)%r_coef(1,1),1,&
                              &     jp_RARnu(ispin)%r_coef(1,1),1)
                         CALL DAXPY(nsgfa*nsgfb,-rbc(ii),d_block(ispin)%r_coef(1,1),1,&
                              &     jp_RARnu(ispin)%r_coef(1,1),1)
                         CALL DAXPY(nsgfa*nsgfb,rbc(iii),c_block(ispin)%r_coef(1,1),1,&
                              &     jp_RARnu(ispin)%r_coef(1,1),1)
                      ENDDO
                      !
                      ! Get the d_A's for the hard and soft densities.
                      NULLIFY(C_coeff_hh_a)
                      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                           &              block_row=katom,&
                           &              block_col=iatom,&
                           &              BLOCK=C_coeff_hh_a)
                      CPPrecondition(ASSOCIATED(C_coeff_hh_a),cp_failure_level,routineP,error,failure)
                      !
                      NULLIFY(C_coeff_ss_a)
                      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                           &              block_row=katom,&
                           &              block_col=iatom,&
                           &              BLOCK=C_coeff_ss_a)
                      CPPrecondition(ASSOCIATED(C_coeff_ss_a),cp_failure_level,routineP,error,failure)
                      !
                      ! Get the d_B's for the hard and soft densities.
                      NULLIFY(C_coeff_hh_b)
                      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                           &              block_row=katom,&
                           &              block_col=jatom,&
                           &              BLOCK=C_coeff_hh_b)
                      CPPrecondition(ASSOCIATED(C_coeff_hh_b),cp_failure_level,routineP,error,failure)
                      !
                      NULLIFY(C_coeff_ss_b)
                      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                           &              block_row=katom,&
                           &              block_col=jatom,&
                           &              BLOCK=C_coeff_ss_b)
                      CPPrecondition(ASSOCIATED(C_coeff_ss_b),cp_failure_level,routineP,error,failure)
                      !
                      DO ispin = 1,nspins
                         !------------------------------------------------------------------
                         ! P_\alpha\alpha'
                         ! hard
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,a_block(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc0_h(ispin)%r_coef(1,1),nsoctot)
                         !
                         ! soft
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &  1.0_dp,a_block(ispin)%r_coef(1,1),nsgfa,&
                                 &  C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                                 &  0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc0_s(ispin)%r_coef(1,1),nsoctot)
                         !------------------------------------------------------------------
                         ! mQai_\alpha\alpha'
                         ! hard
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,jp_RARnu(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_h(ispin)%r_coef(1,1),nsoctot)
                         !
                         ! soft
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,jp_RARnu(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_s(ispin)%r_coef(1,1),nsoctot)
                         !------------------------------------------------------------------
                         ! Qci_\alpha\alpha'
                         ! hard
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,c_block(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_ii_h(ispin)%r_coef(1,1),nsoctot)
                         !
                         ! soft
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,c_block(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_ii_s(ispin)%r_coef(1,1),nsoctot)
                         !------------------------------------------------------------------
                         ! Qbi_\alpha\alpha'
                         ! hard
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,d_block(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_h(ispin)%r_coef(1,1),nsoctot)
                         !
                         ! soft
                         CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                              &     1.0_dp,d_block(ispin)%r_coef(1,1),nsgfa,&
                              &     C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                              &     0.0_dp,PC(1,1),max_nsgf)
                         CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                              &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                              &     PC(1,1),max_nsgf,&
                              &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_s(ispin)%r_coef(1,1),nsoctot)
                         !------------------------------------------------------------------
                      ENDDO ! ispin
                      !
                   ENDDO ! kneighbor
                ENDDO ! kkind
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             ENDDO ! inode
             !
          ENDDO ! ilist
          DO ispin = 1,nspins
             DEALLOCATE(jp_RARnu(ispin)%r_coef,jp2_RARnu(ispin)%r_coef,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDDO
       ENDDO ! jkind
    ENDDO ! ikind
    !
    ! parallel sum up 
    DO ikind = 1 ,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            &               orb_basis_set=orb_basis_set,&
            &               atom_list=atom_list,&
            &               natom=nat,&
            &               paw_atom=paw_atom)
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            &                 nset=nsetc,maxso=maxsoc)
       nsoctot = maxsoc*nsetc
       !
       IF(.NOT. paw_atom) CYCLE
       !
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit(nat,num_pe,mepos)
       !
       ALLOCATE(zero_coeff(nsoctot,nsoctot),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       DO iat = 1,nat
          iatom = atom_list(iat)
          is_not_associated = .NOT.ASSOCIATED(jrho1_atom_set(iatom)%cjc0_h(1)%r_coef)
          !
          IF(iat.GE.bo(1).AND.iat.LE.bo(2).AND.is_not_associated) THEN
             CALL allocate_jrho_coeff(jrho1_atom_set,iatom,nsoctot,error)
          ENDIF
          !
          DO ispin = 1, nspins
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             !
             CALL mp_sum(tmp_coeff,para_env%group)
             tmp_coeff => jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
             !
             tmp_coeff => jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef
             IF(is_not_associated) THEN
                zero_coeff = 0.0_dp; tmp_coeff => zero_coeff
             ENDIF
             CALL mp_sum(tmp_coeff,para_env%group)
          ENDDO ! ispin
       ENDDO ! iat
       !
       DEALLOCATE(zero_coeff,STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       !
    ENDDO ! ikind
    !
    ! clean up
    DEALLOCATE(PC,a_block,b_block,c_block,d_block,jp_RARnu,jp2_RARnu,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE calculate_jrho_atom_coeff


  SUBROUTINE calculate_jrho_atom_rad(qs_env,current_env,IB,idir,error)
    !
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(current_env_type)                   :: current_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i1, i2, iat, iatom, icg, ikind, ipgf1, ipgf2, ir, &
      iset1, iset2, iso, iso1, iso1_first, iso1_last, iso2, iso2_first, &
      iso2_last, ispin, istat, l, l_iso, lmax12, lmax_expansion, lmin12, m1s, &
      m2s, max_iso_not0, max_nso, maxl, maxso, mepos, n1s, n2s, na, natom, &
      natom_tot, nkind, nr, nset, nspins, num_pe, size1, size2
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list, lmax, lmin, npgf
    LOGICAL                                  :: failure, paw_atom
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: is_set_to_0
    REAL(dp)                                 :: hard_radius
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2, gauge_h, gauge_s
    REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: cjc0_h_block, cjc0_s_block, &
      cjc_h_block, cjc_ii_h_block, cjc_ii_s_block, cjc_iii_h_block, &
      cjc_iii_s_block, cjc_s_block, dgg_1, gg, gg_lm1
    REAL(dp), DIMENSION(:, :), POINTER :: coeff, Fr_a_h, Fr_a_h_ii, &
      Fr_a_h_iii, Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, &
      Fr_b_h_iii, Fr_b_s, Fr_b_s_ii, Fr_b_s_iii, Fr_h, Fr_s, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz_asym
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom

!
!

    CALL timeset(routineN,handle)
    !
    failure =.FALSE.
    NULLIFY(atomic_kind_set,dft_control,para_env,&
         &  coeff,Fr_h,Fr_s,Fr_a_h,Fr_a_s,Fr_a_h_ii,Fr_a_s_ii,&
         &  Fr_a_h_iii,Fr_a_s_iii,Fr_b_h,Fr_b_s,Fr_b_h_ii,&
         &  Fr_b_s_ii,Fr_b_h_iii,Fr_b_s_iii,jrho1_atom_set,&
         &  jrho1_atom)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          dft_control=dft_control,&
         &          para_env=para_env,&
         &          error=error)
    !
    CALL get_current_env(current_env=current_env,&
         &           jrho1_atom_set=jrho1_atom_set,&
         &           error=error)
    !

    nkind = SIZE(atomic_kind_set)
    nspins = dft_control%nspins
    !
    natom_tot = SIZE(jrho1_atom_set,1)
    ALLOCATE(is_set_to_0(natom_tot,nspins),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    is_set_to_0(:,:) = .FALSE.

    !
    DO ikind = 1,nkind
       NULLIFY(atom_kind,atom_list,grid_atom,harmonics,orb_basis_set,&
            &  lmax,lmin,npgf,zet,grid_atom,harmonics,my_CG,my_CG_dxyz_asym)
       !
       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
            &               atom_list=atom_list,&
            &               grid_atom=grid_atom,&
            &               natom=natom,&
            &               paw_atom=paw_atom,&
            &               harmonics=harmonics,&
            &               hard_radius=hard_radius,&
            &               orb_basis_set=orb_basis_set)
       !
       ! Quick cycle if needed.
       IF(.NOT.paw_atom) CYCLE
       !
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            &                 lmax=lmax,lmin=lmin,&
            &                 maxl=maxl,npgf=npgf,&
            &                 nset=nset,zet=zet,&
            &                 maxso=maxso)
       !
       nr = grid_atom%nr
       na = grid_atom%ng_sphere
       max_iso_not0 = harmonics%max_iso_not0
       lmax_expansion = indso(1,max_iso_not0)
       !
       ! Distribute the atoms of this kind
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit(natom,num_pe,mepos)
       !
       my_CG => harmonics%my_CG
       my_CG_dxyz_asym => harmonics%my_CG_dxyz_asym
       !
       ! Allocate some arrays.
       max_nso = nsoset(maxl)
       ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),gg_lm1(nr,0:2*maxl),dgg_1(nr,0:2*maxl),&
            &   cjc0_h_block(max_nso,max_nso),cjc0_s_block(max_nso,max_nso),&
            &   cjc_h_block(max_nso,max_nso),cjc_s_block(max_nso,max_nso),&
            &   cjc_ii_h_block(max_nso,max_nso),cjc_ii_s_block(max_nso,max_nso),&
            &   cjc_iii_h_block(max_nso,max_nso),cjc_iii_s_block(max_nso,max_nso),&
            &   gauge_h(nr),gauge_s(nr),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       ! Compute the gauge
       SELECT CASE(current_env%gauge)
       CASE(current_gauge_r)
          ! d(r)=r
          gauge_h(1:nr) = grid_atom%rad(1:nr)
          gauge_s(1:nr) = grid_atom%rad(1:nr)
       CASE(current_gauge_r_and_step_func)
          ! step function
          gauge_h(1:nr) = 0e0_dp
          DO ir=1,nr
             IF(grid_atom%rad(ir).LE.hard_radius) THEN
                gauge_s(ir) = grid_atom%rad(ir)
             ELSE
                gauge_s(ir) = gauge_h(ir)
             ENDIF
          ENDDO
       CASE(current_gauge_atom)
          ! d(r)=A
          gauge_h(1:nr) = HUGE(0e0_dp)!0e0_dp
          gauge_s(1:nr) = HUGE(0e0_dp)!0e0_dp
       CASE DEFAULT
          CALL stop_program(routineP,"Unknown gauge, try again...")
       END SELECT
       !
       !
       m1s = 0
       DO iset1 = 1,nset
          n1s = nsoset(lmax(iset1))
          DO ipgf1  = 1,npgf(iset1)
             iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
             iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
             size1 = iso1_last - iso1_first + 1
             i1 = nsoset(lmin(iset1)-1)+1
             !
             g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))
             !
             m2s = 0
             DO iset2 = 1,nset
                n2s=nsoset(lmax(iset2))
                DO ipgf2 = 1,npgf(iset2)
                   iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
                   iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
                   size2 = iso2_last - iso2_first + 1
                   i2 = nsoset(lmin(iset2)-1)+1
                   !
                   g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))
                   !
                   lmin12 = lmin(iset1)+lmin(iset2)
                   lmax12 = lmax(iset1)+lmax(iset2)
                   !
                   gg = 0.0_dp
                   gg_lm1 = 0.0_dp
                   dgg_1 = 0.0_dp
                   !
                   ! Take only the terms of expansion for L < lmax_expansion
                   IF(lmin12.LE.lmax_expansion) THEN
                      !
                      IF(lmax12.GT.lmax_expansion) lmax12 = lmax_expansion
                      !
                      IF(lmin12==0) THEN
                         gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
                         gg_lm1(1:nr,lmin12) = 0.0_dp
                      ELSE
                         gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
                         gg_lm1(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12-1)*g1(1:nr)*g2(1:nr)
                      ENDIF
                      !
                      DO l=lmin12+1,lmax12
                         gg(1:nr,l) = grid_atom%rad(1:nr)*gg(1:nr,l-1)
                         gg_lm1(1:nr,l) = gg(1:nr,l-1)
                      ENDDO
                      !
                      DO l=lmin12,lmax12
                         dgg_1(1:nr,l)=2.0_dp*(zet(ipgf1,iset1)-zet(ipgf2,iset2))&
                              &              *gg(1:nr,l)*grid_atom%rad(1:nr)
                      ENDDO
                   ELSE
                      CYCLE
                   ENDIF ! lmin12
                   !
                   DO iat = bo(1),bo(2)
                      iatom = atom_list(iat)
                      !
                      DO ispin=1,nspins
                         !------------------------------------------------------------------
                         ! P_\alpha\alpha'
                         cjc0_h_block = HUGE(1.0_dp)
                         cjc0_s_block = HUGE(1.0_dp)
                         !
                         ! Hard term
                         coeff => jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef
                         cjc0_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !
                         ! Soft term
                         coeff => jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef
                         cjc0_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !------------------------------------------------------------------
                         ! mQai_\alpha\alpha'
                         cjc_h_block = HUGE(1.0_dp)
                         cjc_s_block = HUGE(1.0_dp)
                         !
                         ! Hard term
                         coeff => jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef
                         cjc_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !
                         ! Soft term
                         coeff => jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef
                         cjc_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !------------------------------------------------------------------
                         ! Qci_\alpha\alpha'
                         cjc_ii_h_block = HUGE(1.0_dp)
                         cjc_ii_s_block = HUGE(1.0_dp)
                         !
                         ! Hard term
                         coeff => jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef
                         cjc_ii_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &   coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !
                         ! Soft term
                         coeff => jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef
                         cjc_ii_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &   coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !------------------------------------------------------------------
                         ! Qbi_\alpha\alpha'
                         cjc_iii_h_block = HUGE(1.0_dp)
                         cjc_iii_s_block = HUGE(1.0_dp)
                         !
                         !
                         ! Hard term
                         coeff => jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef
                         cjc_iii_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &    coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !
                         ! Soft term
                         coeff => jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef
                         cjc_iii_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                              &    coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                         !------------------------------------------------------------------
                         !
                         ! Allocation radial functions
                         jrho1_atom => jrho1_atom_set(iatom)
                         IF(.NOT.ASSOCIATED(jrho1_atom%jrho_a_h(ispin)%r_coef)) THEN
                            CALL allocate_jrho_atom_rad(jrho1_atom,ispin,nr,na,max_iso_not0,error=error)
                            is_set_to_0(iatom,ispin) = .TRUE.
                         ELSE
                            IF(.NOT.is_set_to_0(iatom,ispin)) THEN
                               CALL set2zero_jrho_atom_rad(jrho1_atom,ispin,error=error)
                               is_set_to_0(iatom,ispin) = .TRUE.
                            ENDIF
                         ENDIF
                         !------------------------------------------------------------------
                         ! 
                         Fr_h => jrho1_atom%jrho_h(ispin)%r_coef
                         Fr_s => jrho1_atom%jrho_s(ispin)%r_coef
                         !------------------------------------------------------------------
                         ! 
                         Fr_a_h => jrho1_atom%jrho_a_h(ispin)%r_coef
                         Fr_a_s => jrho1_atom%jrho_a_s(ispin)%r_coef
                         Fr_b_h => jrho1_atom%jrho_b_h(ispin)%r_coef
                         Fr_b_s => jrho1_atom%jrho_b_s(ispin)%r_coef
                         !------------------------------------------------------------------
                         !
                         Fr_a_h_ii => jrho1_atom%jrho_a_h_ii(ispin)%r_coef
                         Fr_a_s_ii => jrho1_atom%jrho_a_s_ii(ispin)%r_coef
                         Fr_b_h_ii => jrho1_atom%jrho_b_h_ii(ispin)%r_coef
                         Fr_b_s_ii => jrho1_atom%jrho_b_s_ii(ispin)%r_coef
                         !------------------------------------------------------------------
                         !
                         Fr_a_h_iii => jrho1_atom%jrho_a_h_iii(ispin)%r_coef
                         Fr_a_s_iii => jrho1_atom%jrho_a_s_iii(ispin)%r_coef
                         Fr_b_h_iii => jrho1_atom%jrho_b_h_iii(ispin)%r_coef
                         Fr_b_s_iii => jrho1_atom%jrho_b_s_iii(ispin)%r_coef
                         !------------------------------------------------------------------
                         !
                         DO iso = 1,max_iso_not0
                            l_iso = indso(1,iso) ! not needed
                            !
                            DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                               !
                               iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                               iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                               IF(.NOT.(iso2>0.AND.iso1>0))THEN
                                  WRITE(*,*) 'iso1=',iso1
                                  WRITE(*,*) 'iso2=',iso2
                                  WRITE(*,*) '.... will stop!'
                               ENDIF
                               CPPrecondition(iso2>0.AND.iso1>0,cp_failure_level,routineP,error,failure)
                               !
                               l = indso(1,iso1) + indso(1,iso2)
                               IF(l.GT.lmax_expansion.OR.l.LT..0) THEN
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 1 l',l
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 1 lmax_expansion',lmax_expansion
                               ENDIF
                               CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
                               !------------------------------------------------------------------
                               ! P0
                               !
                               IF(current_env%gauge.EQ.current_gauge_atom) THEN
                                  ! Hard term
                                  Fr_h(1:nr,iso) = Fr_h(1:nr,iso) + &
                                       &           gg(1:nr,l)*cjc0_h_block(iso1,iso2)*&
                                       &           my_CG(iso1,iso2,iso)
                                  ! Soft term
                                  Fr_s(1:nr,iso) = Fr_s(1:nr,iso) + &
                                       &           gg(1:nr,l)*cjc0_s_block(iso1,iso2)*&
                                       &           my_CG(iso1,iso2,iso)
                               ELSE
                                  ! Hard term
                                  Fr_h(1:nr,iso) = Fr_h(1:nr,iso) + &
                                       &           gg(1:nr,l)*cjc0_h_block(iso1,iso2)*&
                                       &           my_CG(iso1,iso2,iso)*(grid_atom%rad(1:nr)-gauge_h(1:nr))
                                  ! Soft term
                                  Fr_s(1:nr,iso) = Fr_s(1:nr,iso) + &
                                       &           gg(1:nr,l)*cjc0_s_block(iso1,iso2)*&
                                       &           my_CG(iso1,iso2,iso)*(grid_atom%rad(1:nr)-gauge_s(1:nr))
                               ENDIF
                               !------------------------------------------------------------------
                               ! Rai
                               !
                               ! Hard term
                               Fr_a_h(1:nr,iso) = Fr_a_h(1:nr,iso) + &
                                    &             dgg_1(1:nr,l)*cjc_h_block(iso1,iso2)*&
                                    &             my_CG(iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_a_s(1:nr,iso) = Fr_a_s(1:nr,iso) + &
                                    &             dgg_1(1:nr,l)*cjc_s_block(iso1,iso2)*&
                                    &             my_CG(iso1,iso2,iso)
                               !------------------------------------------------------------------
                               ! Rci
                               !
                               IF(current_env%gauge.EQ.current_gauge_atom) THEN
                                  ! Hard term
                                  Fr_a_h_ii(1:nr,iso) = Fr_a_h_ii(1:nr,iso) + &
                                       &                dgg_1(1:nr,l)* &
                                       &                cjc_ii_h_block(iso1,iso2)* &
                                       &                my_CG(iso1,iso2,iso)
                                  ! Soft term
                                  Fr_a_s_ii(1:nr,iso) = Fr_a_s_ii(1:nr,iso) + &
                                       &                dgg_1(1:nr,l)* &
                                       &                cjc_ii_s_block(iso1,iso2)* &
                                       &                my_CG(iso1,iso2,iso)
                               ELSE
                                  ! Hard term
                                  Fr_a_h_ii(1:nr,iso) = Fr_a_h_ii(1:nr,iso) + &
                                       &                dgg_1(1:nr,l)*gauge_h(1:nr)* &
                                       &                cjc_ii_h_block(iso1,iso2)* &
                                       &                my_CG(iso1,iso2,iso)
                                  ! Soft term
                                  Fr_a_s_ii(1:nr,iso) = Fr_a_s_ii(1:nr,iso) + &
                                       &                dgg_1(1:nr,l)*gauge_s(1:nr)* &
                                       &                cjc_ii_s_block(iso1,iso2)* &
                                       &                my_CG(iso1,iso2,iso)
                               ENDIF
                               !------------------------------------------------------------------
                               ! Rbi
                               !
                               IF(current_env%gauge.EQ.current_gauge_atom) THEN
                                  ! Hard term
                                  Fr_a_h_iii(1:nr,iso) = Fr_a_h_iii(1:nr,iso) + &
                                       &                 dgg_1(1:nr,l)* &
                                       &                 cjc_iii_h_block(iso1,iso2)* &
                                       &                 my_CG(iso1,iso2,iso)
                                  ! Soft term
                                  Fr_a_s_iii(1:nr,iso) = Fr_a_s_iii(1:nr,iso) + &
                                       &                 dgg_1(1:nr,l)* &
                                       &                 cjc_iii_s_block(iso1,iso2)* &
                                       &                 my_CG(iso1,iso2,iso)
                               ELSE
                                  ! Hard term
                                  Fr_a_h_iii(1:nr,iso) = Fr_a_h_iii(1:nr,iso) + &
                                       &                 dgg_1(1:nr,l)*gauge_h(1:nr)* &
                                       &                 cjc_iii_h_block(iso1,iso2)* &
                                       &                 my_CG(iso1,iso2,iso)
                                  ! Soft term
                                  Fr_a_s_iii(1:nr,iso) = Fr_a_s_iii(1:nr,iso) + &
                                       &                 dgg_1(1:nr,l)*gauge_s(1:nr)* &
                                       &                 cjc_iii_s_block(iso1,iso2)* &
                                       &                 my_CG(iso1,iso2,iso)
                               ENDIF
                               !------------------------------------------------------------------
                            ENDDO !icg
                            !
                            !
                            DO icg = 1,harmonics%ncg_dxyz_asym(iso,iset1,iset2)
                               !
                               iso1 = harmonics%ind_cg_a_dxyz_asym(icg,iso,iset1,iset2)
                               iso2 = harmonics%ind_cg_b_dxyz_asym(icg,iso,iset1,iset2)
                               !
                               IF(.NOT.(iso2>0.AND.iso1>0))THEN
                                  WRITE(*,*) 'harmonics%ncg_dxyz_asym'
                                  WRITE(*,*) 'harmonics%ncg_dxyz_asym',harmonics%ncg_dxyz_asym(iso,iset1,iset2)
                                  WRITE(*,*) 'harmonics%ncg_dxyz',harmonics%ncg_dxyz(iso,iset1,iset2)
                                  WRITE(*,*) 'size(indso,1)=',SIZE(indso,1),' size(indso,2)=',SIZE(indso,2)
                                  WRITE(*,*) 'ubound(indso,2)=',UBOUND(indso,2)
                                  WRITE(*,*) 'iso1=',iso1,' iso2=',iso2,' iso=',iso,' icg=',icg
                                  WRITE(*,*) '.... will stop!'
                               ENDIF
                               CPPrecondition(iso2>0.AND.iso1>0,cp_failure_level,routineP,error,failure)
                               !
                               l = indso(1,iso1) + indso(1,iso2)
                               IF(l.GT.lmax_expansion) THEN
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 indso',LBOUND(indso,2),UBOUND(indso,2)
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 iso2',iso2
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 iso1',iso1
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 indso(1,iso1)',indso(1,iso1)
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 indso(1,iso2)',indso(1,iso2)
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 l',l
                                  WRITE(*,*) 'calculate_jrho_atom_rad: 2 lmax_expansion',lmax_expansion
                               ENDIF
                               CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
                               !------------------------------------------------------------------
                               ! Daij
                               !
                               ! Hard term
                               Fr_b_h(1:nr,iso) = Fr_b_h(1:nr,iso) + &
                                    &             gg_lm1(1:nr,l)*cjc_h_block(iso1,iso2)*&
                                    &             my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_b_s(1:nr,iso) = Fr_b_s(1:nr,iso) + &
                                    &             gg_lm1(1:nr,l)*cjc_s_block(iso1,iso2)*&
                                    &             my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               !
                               !------------------------------------------------------------------
                               ! Dcij
                               !
                               IF(current_env%gauge.EQ.current_gauge_atom) THEN
                                  ! Hard term
                                  Fr_b_h_ii(1:nr,iso) = Fr_b_h_ii(1:nr,iso) + &
                                       &                gg_lm1(1:nr,l)* &
                                       &                cjc_ii_h_block(iso1,iso2)* &
                                       &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                                  ! Soft term
                                  Fr_b_s_ii(1:nr,iso) = Fr_b_s_ii(1:nr,iso) + &
                                       &                gg_lm1(1:nr,l)* &
                                       &                cjc_ii_s_block(iso1,iso2)* &
                                       &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               ELSE
                                  ! Hard term
                                  Fr_b_h_ii(1:nr,iso) = Fr_b_h_ii(1:nr,iso) + &
                                       &                gg_lm1(1:nr,l)*gauge_h(1:nr)* &
                                       &                cjc_ii_h_block(iso1,iso2)* &
                                       &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                                  ! Soft term
                                  Fr_b_s_ii(1:nr,iso) = Fr_b_s_ii(1:nr,iso) + &
                                       &                gg_lm1(1:nr,l)*gauge_s(1:nr)* &
                                       &                cjc_ii_s_block(iso1,iso2)* &
                                       &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               ENDIF
                               !------------------------------------------------------------------
                               ! Dbij
                               !
                               IF(current_env%gauge.EQ.current_gauge_atom) THEN
                                  ! Hard term
                                  Fr_b_h_iii(1:nr,iso) = Fr_b_h_iii(1:nr,iso) + &
                                       &                 gg_lm1(1:nr,l)* &
                                       &                 cjc_iii_h_block(iso1,iso2)* &
                                       &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
                                  ! Soft term
                                  Fr_b_s_iii(1:nr,iso) = Fr_b_s_iii(1:nr,iso) + &
                                       &                 gg_lm1(1:nr,l)* &
                                       &                 cjc_iii_s_block(iso1,iso2)* &
                                       &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               ELSE
                                  ! Hard term
                                  Fr_b_h_iii(1:nr,iso) = Fr_b_h_iii(1:nr,iso) + &
                                       &                 gg_lm1(1:nr,l)*gauge_h(1:nr)* &
                                       &                 cjc_iii_h_block(iso1,iso2)* &
                                       &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
                                  ! Soft term
                                  Fr_b_s_iii(1:nr,iso) = Fr_b_s_iii(1:nr,iso) + &
                                       &                 gg_lm1(1:nr,l)*gauge_s(1:nr)* &
                                       &                 cjc_iii_s_block(iso1,iso2)* &
                                       &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               ENDIF
                               !------------------------------------------------------------------
                            ENDDO ! icg
                         ENDDO  ! iso
                         !
                      ENDDO ! ispin
                   ENDDO ! iat
                   !
                   !------------------------------------------------------------------
                   ! 
                ENDDO !ipgf2
                m2s = m2s+maxso
             ENDDO ! iset2
          ENDDO ! ipgf1
          m1s = m1s+maxso
       ENDDO ! iset1
       !
       DEALLOCATE(cjc0_h_block,cjc0_s_block,cjc_h_block,cjc_s_block,cjc_ii_h_block,cjc_ii_s_block,&
            &     cjc_iii_h_block,cjc_iii_s_block,g1,g2,gg,gg_lm1,dgg_1,gauge_h,gauge_s,&
            &     STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO ! ikind
    !
    !
    DEALLOCATE(is_set_to_0,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(handle)
    !
  END SUBROUTINE calculate_jrho_atom_rad


  SUBROUTINE calculate_jrho_atom_ang(jrho1_atom,jrho_h,jrho_s,grid_atom,&
       &                             harmonics,do_igaim,ratom,natm_gauge,&
       &                             hard_radius,iB,idir,ispin,error)
    !
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h, jrho_s
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    LOGICAL, INTENT(IN)                      :: do_igaim
    INTEGER, INTENT(IN)                      :: iB, idir, ispin, natm_gauge
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(dp), INTENT(IN) :: hard_radius, ratom(:,:)

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_ang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, idir2, iiB, iiiB, ir, &
                                                iso, max_iso_not0, na, nr, istat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rad_part, scale
    REAL(dp), DIMENSION(:, :), POINTER :: a, Fr_a_h, Fr_a_h_ii, Fr_a_h_iii, &
      Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, Fr_b_h_iii, Fr_b_s, &
      Fr_b_s_ii, Fr_b_s_iii, Fr_h, Fr_s, slm
    REAL(dp), DIMENSION(:), POINTER :: r
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: g
!
!
    failure = .FALSE.
    NULLIFY(Fr_h,Fr_s,Fr_a_h,Fr_a_s,Fr_a_h_ii,Fr_a_s_ii,Fr_a_h_iii,Fr_a_s_iii,&
         &  Fr_b_h,Fr_b_s,Fr_b_h_ii,Fr_b_s_ii,Fr_b_h_iii,Fr_b_s_iii,&
         &  a,slm)
    !
    CPPrecondition(ASSOCIATED(jrho_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)
    ! just to be sure...
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h_ii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s_ii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h_ii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s_ii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h_iii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s_iii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h_iii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s_iii),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_h_iii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_a_s_iii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_h_iii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_b_s_iii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    !
    !
    nr = grid_atom%nr
    na = grid_atom%ng_sphere
    max_iso_not0 = harmonics%max_iso_not0
    ALLOCATE(g(3,nr,na),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !------------------------------------------------------------------
    !
    Fr_h => jrho1_atom%jrho_h(ispin)%r_coef
    Fr_s => jrho1_atom%jrho_s(ispin)%r_coef
    !------------------------------------------------------------------
    !
    Fr_a_h => jrho1_atom%jrho_a_h(ispin)%r_coef !Rai
    Fr_a_s => jrho1_atom%jrho_a_s(ispin)%r_coef
    Fr_b_h => jrho1_atom%jrho_b_h(ispin)%r_coef !Daij
    Fr_b_s => jrho1_atom%jrho_b_s(ispin)%r_coef
    !------------------------------------------------------------------
    !
    Fr_a_h_ii => jrho1_atom%jrho_a_h_ii(ispin)%r_coef !Rci
    Fr_a_s_ii => jrho1_atom%jrho_a_s_ii(ispin)%r_coef
    Fr_b_h_ii => jrho1_atom%jrho_b_h_ii(ispin)%r_coef !Dcij
    Fr_b_s_ii => jrho1_atom%jrho_b_s_ii(ispin)%r_coef
    !------------------------------------------------------------------
    !
    Fr_a_h_iii => jrho1_atom%jrho_a_h_iii(ispin)%r_coef !Rbi
    Fr_a_s_iii => jrho1_atom%jrho_a_s_iii(ispin)%r_coef
    Fr_b_h_iii => jrho1_atom%jrho_b_h_iii(ispin)%r_coef !Dbij
    Fr_b_s_iii => jrho1_atom%jrho_b_s_iii(ispin)%r_coef
    !------------------------------------------------------------------
    !
    a   => harmonics%a
    slm => harmonics%slm
    r   => grid_atom%rad
    !
    CALL set_vecp(iB,iiB,iiiB)
    !
    scale = 0.0_dp
    idir2 = 1
    IF(idir.NE.iB) THEN
       CALL set_vecp_rev(idir,iB,idir2)
       scale = fac_vecp(idir,iB,idir2)
    ENDIF
    !
    ! Set the gauge
    CALL get_gauge()
    !
    DO ir = 1,nr
       DO iso = 1,max_iso_not0
          DO ia = 1,na
             IF(do_igaim) THEN
                !------------------------------------------------------------------
                ! Hard current density response
                ! radial(ia,ir) = (               aj(ia) * Rai(ir,iso) + Daij 
                !                  -  aii(ia) * ( aj(ia) * Rbi(ir,iso) + Dbij )
                !                  + aiii(ia) * ( aj(ia) * Rci(ir,iso) + Dcij ) 
                !                 ) * Ylm(ia)
                rad_part =                a(idir,ia)*Fr_a_h    (ir,iso)+Fr_b_h    (ir,iso) &
                     &   - g( iiB,ir,ia)*(a(idir,ia)*Fr_a_h_iii(ir,iso)+Fr_b_h_iii(ir,iso))&
                     &   + g(iiiB,ir,ia)*(a(idir,ia)*Fr_a_h_ii (ir,iso)+Fr_b_h_ii (ir,iso))&
                     &   + scale*(a(idir2,ia)*r(ir)-g(idir2,ir,ia))*Fr_h(ir,iso)
                !
                jrho_h(ir,ia) = jrho_h(ir,ia) + rad_part * slm(ia,iso)
                !------------------------------------------------------------------
                ! Soft current density response
                rad_part =                a(idir,ia)*Fr_a_s    (ir,iso)+Fr_b_s    (ir,iso) &
                     &   - g( iiB,ir,ia)*(a(idir,ia)*Fr_a_s_iii(ir,iso)+Fr_b_s_iii(ir,iso))&
                     &   + g(iiiB,ir,ia)*(a(idir,ia)*Fr_a_s_ii (ir,iso)+Fr_b_s_ii (ir,iso))&
                     &   + scale*(a(idir2,ia)*r(ir)-g(idir2,ir,ia))*Fr_s(ir,iso)
                !
                jrho_s(ir,ia) = jrho_s(ir,ia) + rad_part * slm(ia,iso)
                !------------------------------------------------------------------
             ELSE
                !------------------------------------------------------------------
                ! Hard current density response
                ! radial(ia,ir) = (               aj(ia) * Rai(ir,iso) + Daij 
                !                  -  aii(ia) * ( aj(ia) * Rbi(ir,iso) + Dbij )
                !                  + aiii(ia) * ( aj(ia) * Rci(ir,iso) + Dcij ) 
                !                 ) * Ylm(ia)
                rad_part =             a(idir,ia)*Fr_a_h    (ir,iso)+Fr_b_h    (ir,iso) &
                     &   - a( iiB,ia)*(a(idir,ia)*Fr_a_h_iii(ir,iso)+Fr_b_h_iii(ir,iso))&
                     &   + a(iiiB,ia)*(a(idir,ia)*Fr_a_h_ii (ir,iso)+Fr_b_h_ii (ir,iso))&
                     &   + scale*a(idir2,ia)*Fr_h(ir,iso)
                !
                jrho_h(ir,ia) = jrho_h(ir,ia) + rad_part * slm(ia,iso)
                !------------------------------------------------------------------
                ! Soft current density response
                rad_part =             a(idir,ia)*Fr_a_s    (ir,iso)+Fr_b_s    (ir,iso) &
                     &   - a( iiB,ia)*(a(idir,ia)*Fr_a_s_iii(ir,iso)+Fr_b_s_iii(ir,iso))&
                     &   + a(iiiB,ia)*(a(idir,ia)*Fr_a_s_ii (ir,iso)+Fr_b_s_ii (ir,iso))&
                     &   + scale*a(idir2,ia)*Fr_s(ir,iso)
                !
                jrho_s(ir,ia) = jrho_s(ir,ia) + rad_part * slm(ia,iso)
                !------------------------------------------------------------------
             ENDIF
          ENDDO ! ia
       ENDDO ! iso
    ENDDO ! ir
    !
    DEALLOCATE(g,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
  CONTAINS
    !
    SUBROUTINE get_gauge()
    INTEGER                                  :: iatom, ixyz, jatom
    REAL(dp)                                 :: ab, pa, pb, point(3), pra(3), &
                                                prb(3), res, tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buf

      ALLOCATE(buf(natm_gauge))
      DO ir = 1,nr
      DO ia = 1,na
         DO ixyz = 1,3
            g  (ixyz,ir,ia) = 0.0_dp
         ENDDO
         point(1) = r(ir) * a(1,ia)
         point(2) = r(ir) * a(2,ia)
         point(3) = r(ir) * a(3,ia)
         DO iatom = 1,natm_gauge
            buf(iatom) = 1.0_dp
            pra = point - ratom(:,iatom)
            pa  = SQRT(pra(1)**2 + pra(2)**2 + pra(3)**2 )
            DO jatom = 1,natm_gauge
               IF(iatom.EQ.jatom) CYCLE
               prb = point - ratom(:,jatom)
               pb  = SQRT( prb(1)**2 + prb(2)**2 + prb(3)**2 )
               ab  = SQRT( (pra(1)-prb(1))**2 + (pra(2)-prb(2))**2 + (pra(3)-prb(3))**2 )
               !
               tmp = ( pa - pb ) / ab
               tmp = 0.5_dp * ( 3.0_dp - tmp**2 ) * tmp
               tmp = 0.5_dp * ( 3.0_dp - tmp**2 ) * tmp
               tmp = 0.5_dp * ( 3.0_dp - tmp**2 ) * tmp
               buf(iatom) = buf(iatom) * 0.5_dp * ( 1.0_dp - tmp )
            ENDDO
         ENDDO
         DO ixyz = 1,3
            res = 0.0_dp
            DO iatom = 1,natm_gauge
               res = res + ratom(ixyz,iatom) * buf(iatom)
            ENDDO
            res = res / SUM(buf(1:natm_gauge))
            !
            g(ixyz,ir,ia) = res
         ENDDO
      ENDDO
      ENDDO
      DEALLOCATE(buf)
    END SUBROUTINE get_gauge
  END SUBROUTINE calculate_jrho_atom_ang

  SUBROUTINE calculate_jrho_atom(current_env,qs_env,iB,idir,error)
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, iatom, ikind, ispin, &
                                                istat, jatom, natm_gauge, &
                                                natm_tot, natom, nkind, nspins
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: do_igaim, failure, gapw, &
                                                paw_atom
    REAL(dp)                                 :: hard_radius, r(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ratom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.

    NULLIFY(para_env, dft_control)
    NULLIFY(jrho1_atom_set, grid_atom, harmonics)
    NULLIFY(atomic_kind_set, atom_kind, atom_list)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    cell=cell,&
                    para_env=para_env,error=error) 

    CALL get_current_env(current_env=current_env,&
                         jrho1_atom_set=jrho1_atom_set,&
                         error=error)

    do_igaim = .FALSE.
    IF(current_env%gauge.EQ.current_gauge_atom) do_igaim = .TRUE.

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    natm_tot = SIZE(particle_set)
    ALLOCATE(ratom(3,natm_tot),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (gapw) THEN
       DO ikind = 1,nkind
          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               harmonics=harmonics,grid_atom=grid_atom,&
                               natom=natom, hard_radius=hard_radius,&
                               paw_atom=paw_atom)
          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind
    
          bo = get_limit( natom, para_env%num_pe, para_env%mepos )

          DO iat = bo(1),bo(2)
             iatom = atom_list(iat)
             NULLIFY (jrho1_atom)
             jrho1_atom => jrho1_atom_set(iatom)

             natm_gauge = 0
             DO jatom = 1,natm_tot
                r(:) = pbc(particle_set(jatom)%r(:) - particle_set(iatom)%r(:),cell)
                IF(SQRT(SUM(r(:)**2)).LE.2.0_dp*hard_radius) THEN
                   natm_gauge = natm_gauge + 1
                   ratom(:,natm_gauge) = r(:)
                ENDIF
             ENDDO

             DO ispin = 1,nspins
                jrho1_atom%jrho_vec_rad_h(idir,ispin)%r_coef = 0.0_dp
                jrho1_atom%jrho_vec_rad_s(idir,ispin)%r_coef = 0.0_dp
                CALL calculate_jrho_atom_ang(jrho1_atom,&
                     jrho1_atom%jrho_vec_rad_h(idir,ispin)%r_coef,&
                     jrho1_atom%jrho_vec_rad_s(idir,ispin)%r_coef,&
                     grid_atom,harmonics,&
                     do_igaim,&
                     ratom,natm_gauge,hard_radius,iB,idir,ispin,error=error)
             END DO !ispin
          END DO !iat
       END DO !ikind
    END IF

    DEALLOCATE(ratom,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE calculate_jrho_atom 

! *****************************************************************************
!> \brief Prepare the 3 center lists used to cpmpute the atomic terms of the current density
!>      New lists need to be constructed because the density matrices
!>      are not symmetric, hence, the ij and the ji pairs  are not equivalent and both
!>      need to be considered
!> \par History
!>      07.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE list_3c_jrho_atom(qs_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'list_3c_jrho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, iab, iat, iatom, ikind, ilist, inode, istat, jatom, jkind, &
      last_jatom, nat, natom, nkind, nlist, nnode, nsgf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, new_atom_b, &
                                                oce_present
    REAL(dp)                                 :: rab(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    TYPE(dist_list_4oce_type), &
      DIMENSION(:), POINTER                  :: dist_list
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_oce, sbc_oce
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce

!
!

    failure = .FALSE.

    NULLIFY(sac_oce, sbc_oce, oce, &
         &  atomic_kind,atomic_kind_set,&
         &  sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
         &  sab_orb_neighbor_node,reduced_3c_oce,particle_set,&
         &  dist_list,distribution_2d)
    !
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    sab_all=sab_orb,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce,&
                    particle_set=particle_set,&
                    oce=oce,&
                    &error=error)

    oce_present = ASSOCIATED(sac_oce)

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(oce),cp_failure_level,routineP,error,failure)

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         &                   nsgf=nsgf)

    ALLOCATE(first_sgf(natom),last_sgf(natom),dist_list(natom),oceh(1),oces(1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1,natom
       dist_list(i)%ndist = 0
       NULLIFY(dist_list(i)%rdist)
    ENDDO
    NULLIFY(oceh(1)%block,oces(1)%block)

    CALL get_particle_set(particle_set=particle_set,&
         &                first_sgf=first_sgf,&
         &                last_sgf=last_sgf, error=error)
    !
    ! Create a new list
    CALL create_oce_set(oce,error=error)
    CALL allocate_oce_set(oce,natom,nkind,1,error=error)
    CALL prepare_oce_coeff(particle_set,oce,natom,nsgf,&
         &                 first_sgf,last_sgf,distribution_2d,&
         &                 error=error)
    CALL allocate_vtriple(oce%voce,nkind,natom,error)
    CPPrecondition(ASSOCIATED(oce%voce),cp_failure_level,routineP,error,failure)
    !
    !   *** Retrieve data of SAC_OCE and of SAC_3C neighbors of atom "iatom" ***
    DO ikind = 1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            &               natom=nat,&
            &               atom_list=atom_list)
       IF(oce_present) THEN
          DO iat = 1,nat
             iatom = atom_list(iat)
             CALL retrieve_sac_list(oce%voce,sac_oce,iatom,ikind,nkind,error=error)
          ENDDO
       ENDIF
    ENDDO

    DO ikind=1,nkind

       DO jkind=1,nkind

          iab = ikind + nkind*(jkind - 1)

          IF(.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

          CALL get_neighbor_list_set(neighbor_list_set=&
                                     sab_orb(iab)%neighbor_list_set,&
                                     nlist=nlist)
          NULLIFY( sab_orb_neighbor_list )

          DO ilist=1,nlist
             IF( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
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
                                       r=rab)

                IF(jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                END IF

                CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                       reduced_3c_oce=reduced_3c_oce)

!           *** Build reduced 3-center lists
                CALL build_reduced_3c_lists(atomic_kind_set, iatom, jatom,&
                     ikind, jkind, sbc_oce, oce%voce, reduced_3c_oce,&
                     .TRUE.,qs_env%dft_control%qs_control%eps_rho_rspace,error)
!
!            *** Calculate the coefficients of the one center expansions ***
                CALL set_up_oce(atomic_kind_set,reduced_3c_oce,oce,&
                              ikind,iatom,jkind,jatom,rab,&
                              oceh,oces,0,dist_list,error=error)

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             ENDDO  ! inode

          ENDDO  ! ilist

       ENDDO  ! jkind

    ENDDO  ! ikind
    !
    ! put the oce back to the env
    CALL set_qs_env(qs_env=qs_env,oce=oce,error=error)
    !
    DO i=1,natom
       IF(ASSOCIATED(dist_list(i)%rdist)) THEN
          DEALLOCATE(dist_list(i)%rdist,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDDO
    DEALLOCATE(first_sgf,last_sgf,dist_list,oceh,oces,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
  END SUBROUTINE list_3c_jrho_atom

END MODULE qs_linres_atom_current
