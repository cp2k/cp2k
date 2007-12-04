!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_nmr_current *
!!
!!   NAME
!!     qs_linres_nmr_current
!!
!!   FUNCTION
!!     given the response wavefunctions obtained by the application
!!     of the (rxp), p, and ((dk-dl)xp) operators,
!!     here the current density vector (jx, jy, jz)
!!     is computed for the 3 directions of the magnetic field (Bx, By, Bz)
!!
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 02-2006 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_atom_current

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
  USE qs_linres_nmr_op,                ONLY: fac_vecp,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: allocate_jrho_atom_rad,&
                                             get_nmr_env,&
                                             jrho_atom_type,&
                                             nmr_env_type,&
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
  USE input_constants,                 ONLY: nmr_gauge_r,&
                                             nmr_gauge_r_and_Bader,&
                                             nmr_gauge_r_and_step_function
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: calculate_jrho_atom_coeff,calculate_jrho_atom_rad,list_3c_jrho_atom,&
       &    shift_atom

  ! *** Global parameters (only in this module)
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_atom_current'

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_linres_nmr_atom_current/calculate_jrho_atom_coeff!!
!!   NAME
!!     calculate_jrho_atom_coeff
!!
!!   FUNCTION
!!     Calculate the expansion coefficients for the atomic terms
!!     of the current densitiy in GAPW
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [MI]
!!
!****************************************************************************

  SUBROUTINE calculate_jrho_atom_coeff(qs_env,nmr_env,mat_jp0,mat_jp,mat_jp_pii,&
       &                               mat_jp_piii,iB,idir,error)
    !
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: mat_jp0, mat_jp, mat_jp_pii, &
                                                mat_jp_piii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_atom_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iab, iat, iatom, ii, iii, ikind, ilist, inode, isgfa, &
      isgfb, ispin, istat, jatom, jkind, katom, kkind, kneighbor, last_jatom, &
      last_katom, max_gau, max_nsgf, maxsoc, nat, natom, nkind, nlist, nnode, &
      nnode_c, nsetc, nsgfa, nsgfb, nsoctot, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, new_atom_b, &
                                                new_atom_c, paw_atom
    REAL(dp)                                 :: rab(3), rbc(3)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: PC
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
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
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: jp0_block, jp_block, jp_ii_block, &
                                                jp_iii_block, jp_RARnu
    REAL(dp) :: hard_radius_a,hard_radius_b,hard_radius_c

    CALL timeset("calculate_jrho_atom_coeff","I"," ",handle)
    !
    failure =.FALSE.
    NULLIFY(atom_list,atomic_kind_set,dft_control,orb_basis_set,sab_orb,&
         &  jrho1_atom_set,oce,reduced_3c_oce,jp0_block,jp_block,jp_ii_block,&
         &  jp_ii_block,jp_RARnu,C_coeff_hh_a,C_coeff_hh_b,C_coeff_ss_a,&
         &  C_coeff_ss_b,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
         &  sab_orb_neighbor_node)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          dft_control=dft_control,&
         &          oce=oce,&
         &          sab_all=sab_orb,&
         &          error=error)
    !
    CPPrecondition(ASSOCIATED(oce),cp_failure_level,routineP,error,failure)
    !
    CALL get_nmr_env(nmr_env=nmr_env,&
         &           jrho1_atom_set=jrho1_atom_set,&
         &           error=error)
    !
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         &                   maxsgf=max_nsgf,&
         &                   maxgtops=max_gau)
    !
    ! Density-like matrices
    nkind  = SIZE(atomic_kind_set)
    natom  = SIZE(jrho1_atom_set)
    nspins = dft_control%nspins
    CALL set_vecp(iB,ii,iii)
    !
    ! Allocate some arrays
    ALLOCATE(jp0_block(nspins),jp_block(nspins),jp_ii_block(nspins),jp_iii_block(nspins),&
         &   jp_RARnu(nspins),PC(max_nsgf,max_gau),STAT=istat)
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
             jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef = 0.0_dp

             jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef = 0.0_dp
             jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef = 0.0_dp
             IF(nmr_env%store_current) THEN
                jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,iB,ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,iB,ispin)%r_coef = 0.0_dp
             ELSE
                jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef = 0.0_dp
                jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef = 0.0_dp
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
             NULLIFY(jp_RARnu(ispin)%r_coef)
             ALLOCATE(jp_RARnu(ispin)%r_coef(nsgfa,nsgfb),STAT=istat)
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
                   NULLIFY(jp0_block(ispin)%r_coef)
                   NULLIFY(jp_block(ispin)%r_coef)
                   NULLIFY(jp_ii_block(ispin)%r_coef)
                   NULLIFY(jp_iii_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_jp0(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=jp0_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_jp(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=jp_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_jp_pii(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=jp_ii_block(ispin)%r_coef)
                   CALL get_block_node(matrix=mat_jp_piii(ispin)%matrix,&
                        &              block_row=iatom,&
                        &              block_col=jatom,&
                        &              BLOCK=jp_iii_block(ispin)%r_coef)
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
                      ! check if the hard_radius overlap
                      !IF(hard_radius_b+hard_radius_c.GT.SQRT(rbc(1)**2+rbc(2)**2+rbc(3)**2)&
                      !     & .AND.SQRT(rbc(1)**2+rbc(2)**2+rbc(3)**2).GT.1.0E-10_dp) THEN
                      !   write(*,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
                      !   write(*,*) 'hard_radius_b',hard_radius_b,' hard_radius_c',hard_radius_c
                      !   write(*,*) '|rbc|',SQRT(rbc(1)**2+rbc(2)**2+rbc(3)**2)
                      !   write(*,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
                      !ENDIF
                      !
                      ! Compute the modified Qai matrix as
                      ! mQai_\mu\nu = Qai_\mu\nu - Qbi_\mu\nu * (R_A-R_\nu)_ii 
                      !             + Qci_\mu\nu * (R_A-R_\nu)_iii
                      DO ispin = 1,nspins
                         jp_RARnu(ispin)%r_coef=0.0_dp ! not needed
                      DO isgfb = 1,nsgfb
                      DO isgfa = 1,nsgfa
                         jp_RARnu(ispin)%r_coef(isgfa,isgfb) =  &
                              &  jp_block(ispin)%r_coef(isgfa,isgfb) &
                              & -jp_iii_block(ispin)%r_coef(isgfa,isgfb)*rbc(ii)&
                              & +jp_ii_block(ispin)%r_coef(isgfa,isgfb)*rbc(iii)
                      ENDDO
                      ENDDO
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
                      ! De-contract and copy into the final matrix.
                      ! The hard and soft density matrix-like terms are computed here.
                      IF(iatom==jatom) THEN
                         !
                         DO ispin = 1,nspins
                            !------------------------------------------------------------------
                            ! P_\alpha\alpha'
                            ! hard
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp0_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc0_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp0_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc0_s(ispin)%r_coef(1,1),nsoctot)
                            !------------------------------------------------------------------
                            ! mQai_\alpha\alpha'
                            ! hard
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_RARnu(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_RARnu(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_s(ispin)%r_coef(1,1),nsoctot)
                            !------------------------------------------------------------------
                            ! Qci_\alpha\alpha'
                            ! hard
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_ii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_ii_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_ii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_ii_s(ispin)%r_coef(1,1),nsoctot)
                            !------------------------------------------------------------------
                            ! Qbi_\alpha\alpha'
                            ! hard
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_iii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_iii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_s(ispin)%r_coef(1,1),nsoctot)
                            !------------------------------------------------------------------
                         ENDDO ! ispin
                         !
                      ELSE
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
                                 &     1.0_dp,jp0_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc0_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp0_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
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
                                 &     1.0_dp,jp_ii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_ii_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_ii_block(ispin)%r_coef(1,1),nsgfa,&
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
                                 &     1.0_dp,jp_iii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_hh_b(1,1),SIZE(C_coeff_hh_b,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_hh_a(1,1),SIZE(C_coeff_hh_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_h(ispin)%r_coef(1,1),nsoctot)
                            !
                            ! soft
                            CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 &     1.0_dp,jp_iii_block(ispin)%r_coef(1,1),nsgfa,&
                                 &     C_coeff_ss_b(1,1),SIZE(C_coeff_ss_b,1),&
                                 &     0.0_dp,PC(1,1),max_nsgf)
                            CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 &     1.0_dp,C_coeff_ss_a(1,1),SIZE(C_coeff_ss_a,1),&
                                 &     PC(1,1),max_nsgf,&
                                 &     1.0_dp,jrho1_atom_set(katom)%cjc_iii_s(ispin)%r_coef(1,1),nsoctot)
                            !------------------------------------------------------------------
                         ENDDO ! ispin
                      ENDIF ! iatom==jatom
                      !
                   ENDDO ! kneighbor
                ENDDO ! kkind
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             ENDDO ! inode
             !
          ENDDO ! ilist
          DO ispin = 1,nspins
             DEALLOCATE(jp_RARnu(ispin)%r_coef,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ENDDO
       ENDDO ! jkind
    ENDDO ! ikind
    !
    ! parallel sum up 
    DO iatom = 1,natom
       DO ispin = 1,nspins
          IF(ASSOCIATED(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)) THEN
             CALL mp_sum(jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef,qs_env%para_env%group)
             CALL mp_sum(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef,qs_env%para_env%group)
          ENDIF
       ENDDO
    ENDDO
    !
    ! clean up
    DEALLOCATE(PC,jp0_block,jp_block,jp_ii_block,jp_iii_block,jp_RARnu,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE calculate_jrho_atom_coeff

! *****************************************************************************

  SUBROUTINE calculate_jrho_atom_rad(qs_env,nmr_env,IB,idir,error)
    !
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i1, i2, iat, iatom, icg, ikind, ipgf1, ipgf2, iset1, &
      iset2, iso, iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, &
      ispin, istat, l, l_iso, lmax12, lmax_expansion, lmin12, m1s, m2s, &
      max_iso_not0, max_nso, maxl, maxso, mepos, n1s, n2s, natom, natom_tot, &
      nkind, nr, nset, nspins, num_pe, size1, size2
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list, lmax, lmin, npgf
    LOGICAL                                  :: failure, paw_atom
    LOGICAL, ALLOCATABLE, DIMENSION(:,:)     :: is_set_to_0
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2, gauge_h, gauge_s
    REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: cjc0_h_block,cjc0_s_block,cjc_h_block, cjc_ii_h_block, &
      cjc_ii_s_block, cjc_iii_h_block, cjc_iii_s_block, cjc_s_block, dgg_1, &
      dgg_2, gg, gg_lm1
    REAL(dp), DIMENSION(:, :), POINTER :: coeff, Fr_h, Fr_s, Fr_a_h, Fr_a_h_ii, &
      Fr_a_h_iii, Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, &
      Fr_b_h_iii, Fr_b_s, Fr_b_s_ii, Fr_b_s_iii, zet
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
    REAL(dp)                                 :: hard_radius,alpha
    INTEGER::ir
    !
    !
    CALL timeset("calculate_jrho_atom_rad","I"," ",handle)
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
    CALL get_nmr_env(nmr_env=nmr_env,&
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
       ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),gg_lm1(nr,0:2*maxl),&
            &   dgg_1(nr,0:2*maxl),dgg_2(nr,0:2*maxl),&
            &   cjc0_h_block(max_nso,max_nso),cjc0_s_block(max_nso,max_nso),&
            &   cjc_h_block(max_nso,max_nso),cjc_s_block(max_nso,max_nso),&
            &   cjc_ii_h_block(max_nso,max_nso),cjc_ii_s_block(max_nso,max_nso),&
            &   cjc_iii_h_block(max_nso,max_nso),cjc_iii_s_block(max_nso,max_nso),&
            &   gauge_h(nr),gauge_s(nr),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       ! Compute the gauge
       SELECT CASE(nmr_env%gauge)
       CASE(nmr_gauge_r)
          ! d(r)=r
          gauge_h(1:nr) = grid_atom%rad(1:nr)
          gauge_s(1:nr) = grid_atom%rad(1:nr)
       CASE(nmr_gauge_r_and_bader)
          ! d_h(r)=r*(1-exp(-alpha*r**4)))
          ! d_s(r)=r*   exp(-alpha*r**4)
          alpha = LOG(2.0_dp*hard_radius)/hard_radius**4
          gauge_h(1:nr) = grid_atom%rad(1:nr)*(1.0_dp-EXP(-alpha*grid_atom%rad2(1:nr)*grid_atom%rad2(1:nr)))
          gauge_s(1:nr) = grid_atom%rad(1:nr)*        EXP(-alpha*grid_atom%rad2(1:nr)*grid_atom%rad2(1:nr))
       CASE(nmr_gauge_r_and_step_function)
          ! step function
          gauge_h(1:nr) = 0d0
          do ir=1,nr
             if(grid_atom%rad(ir).LE.hard_radius) then
                gauge_s(ir) = grid_atom%rad(ir)
             else
                gauge_s(ir) = gauge_h(ir)
             endif
          enddo
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
                   dgg_2 = 0.0_dp
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
                         dgg_2(1:nr,l)=dgg_1(1:nr,l)*grid_atom%rad(1:nr)
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
                         cjc0_h_block = HUGE(1.0_dp)!0.0_dp
                         cjc0_s_block = HUGE(1.0_dp)!0.0_dp
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
                         cjc_h_block = HUGE(1.0_dp)!0.0_dp
                         cjc_s_block = HUGE(1.0_dp)!0.0_dp
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
                         cjc_ii_h_block = HUGE(1.0_dp)!0.0_dp
                         cjc_ii_s_block = HUGE(1.0_dp)!0.0_dp
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
                         cjc_iii_h_block = HUGE(1.0_dp)!0.0_dp
                         cjc_iii_s_block = HUGE(1.0_dp)!0.0_dp
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
                            CALL allocate_jrho_atom_rad(jrho1_atom,ispin,nr,max_iso_not0,error=error)
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
                               ! P
                               !
                               ! Hard term
                               Fr_h(1:nr,iso) = Fr_h(1:nr,iso) + &
                                    &           gg(1:nr,l)*cjc0_h_block(iso1,iso2)*&
                                    &           my_CG(iso1,iso2,iso)*(grid_atom%rad(1:nr)-gauge_h(1:nr))
                               !
                               ! Soft term
                               Fr_s(1:nr,iso) = Fr_s(1:nr,iso) + &
                                    &           gg(1:nr,l)*cjc0_s_block(iso1,iso2)*&
                                    &           my_CG(iso1,iso2,iso)*(grid_atom%rad(1:nr)-gauge_s(1:nr))
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
                               ! Hard term
                               Fr_a_h_ii(1:nr,iso) = Fr_a_h_ii(1:nr,iso) + &
                                    &                dgg_1(1:nr,l)*gauge_h(1:nr)*cjc_ii_h_block(iso1,iso2)*&
!                                    &                 dgg_2(1:nr,l)*cjc_ii_h_block(iso1,iso2)*&
                                    &                my_CG(iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_a_s_ii(1:nr,iso) = Fr_a_s_ii(1:nr,iso) + &
                                    &                dgg_1(1:nr,l)*gauge_s(1:nr)*cjc_ii_s_block(iso1,iso2)*&
!                                    &                 dgg_2(1:nr,l)*cjc_ii_s_block(iso1,iso2)*&
                                    &                my_CG(iso1,iso2,iso)
                               !------------------------------------------------------------------
                               ! Rbi
                               !
                               ! Hard term
                               Fr_a_h_iii(1:nr,iso) = Fr_a_h_iii(1:nr,iso) + &
                                    &                 dgg_1(1:nr,l)*gauge_h(1:nr)*cjc_iii_h_block(iso1,iso2)*&
!                                    &                 dgg_2(1:nr,l)*cjc_iii_h_block(iso1,iso2)*&
                                    &                 my_CG(iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_a_s_iii(1:nr,iso) = Fr_a_s_iii(1:nr,iso) + &
                                    &                 dgg_1(1:nr,l)*gauge_s(1:nr)*cjc_iii_s_block(iso1,iso2)*&
!                                    &                 dgg_2(1:nr,l)*cjc_iii_s_block(iso1,iso2)*&
                                    &                 my_CG(iso1,iso2,iso)
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
                               !------------------------------------------------------------------
                               ! Dcij
                               !
                               ! Hard term
                               Fr_b_h_ii(1:nr,iso) = Fr_b_h_ii(1:nr,iso) + &
                                    &                gg_lm1(1:nr,l)*gauge_h(1:nr)*cjc_ii_h_block(iso1,iso2)*&
!                                    &                gg(1:nr,l)*cjc_ii_h_block(iso1,iso2)*&
                                    &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_b_s_ii(1:nr,iso) = Fr_b_s_ii(1:nr,iso) + &
                                    &                gg_lm1(1:nr,l)*gauge_s(1:nr)*cjc_ii_s_block(iso1,iso2)*&
!                                       &                gg(1:nr,l)*cjc_ii_s_block(iso1,iso2)*&
                                    &                my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               !------------------------------------------------------------------
                               ! Dbij
                               !
                               ! Hard term
                               Fr_b_h_iii(1:nr,iso) = Fr_b_h_iii(1:nr,iso) + &
                                    &                 gg_lm1(1:nr,l)*gauge_h(1:nr)*cjc_iii_h_block(iso1,iso2)*&
!                                       &                 gg(1:nr,l)*cjc_iii_h_block(iso1,iso2)*&
                                    &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
                               !
                               ! Soft term
                               Fr_b_s_iii(1:nr,iso) = Fr_b_s_iii(1:nr,iso) + &
                                    &                 gg_lm1(1:nr,l)*gauge_s(1:nr)*cjc_iii_s_block(iso1,iso2)*&
!                                    &                 gg(1:nr,l)*cjc_iii_s_block(iso1,iso2)*&
                                    &                 my_CG_dxyz_asym(idir,iso1,iso2,iso)
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
            &     cjc_iii_h_block,cjc_iii_s_block,g1,g2,gg,gg_lm1,dgg_1,dgg_2,gauge_h,gauge_s,&
            &     STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO ! ikind
    !
    !
    DEALLOCATE(is_set_to_0,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE calculate_jrho_atom_rad

! *****************************************************************************

  SUBROUTINE shift_atom(qs_env,nmr_env,iB,idir,error)
    !
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shift_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, iat, iatom, idir2, idir3, ikind, ir, ispin, istat, &
      jatom, mepos, n_nics, na, natom, natom_tot, nkind, nr, nspins, num_pe, &
      output_unit
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, jatomisneighbor, &
                                                paw_atom, store_current
    REAL(dp) :: ddiff, distance, itegrated_jrho, jrho1_tot, r_iatom(3), &
      r_jatom(3), rdiff(3), rij(3), rvec(3), scale_fac, threshold
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: chemical_shift_loc_tmp
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h_grid, jrho_s_grid, &
                                                r_nics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift_loc, &
                                                chemical_shift_loc_nics
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!
!

    CALL timeset("shift_atom","I"," ",handle)
    !
    failure =.FALSE.
    NULLIFY(atomic_kind_set,cell,dft_control,para_env,particle_set,&
         &  chemical_shift_loc,chemical_shift_loc_nics,jrho1_atom_set,&
         &  jrho1_atom,r_nics,jrho_h_grid,jrho_s_grid,atom_kind,&
         &  atom_list,grid_atom,harmonics,logger)
    !
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_unit_nr(logger)
    !
    CALL get_qs_env(qs_env=qs_env,&
         &          atomic_kind_set=atomic_kind_set,&
         &          cell=cell,&
         &          dft_control=dft_control,&
         &          para_env=para_env,&
         &          particle_set=particle_set,&
         &          error=error)

    CALL get_nmr_env(nmr_env=nmr_env,&
         &           chemical_shift_loc=chemical_shift_loc, &
         &           chemical_shift_loc_nics=chemical_shift_loc_nics, &
         &           jrho1_atom_set=jrho1_atom_set,&
         &           n_nics=n_nics,&
         &           r_nics=r_nics,&
         &           error=error)
    !
    nkind = SIZE(atomic_kind_set,1)
    natom_tot = SIZE(particle_set,1)
    nspins = dft_control%nspins
    store_current = nmr_env%store_current
    itegrated_jrho = 0.0_dp
    !
    ALLOCATE(chemical_shift_loc_tmp(3,natom_tot),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    chemical_shift_loc_tmp = 0.0_dp
    !
    ! Loop over atoms to collocate the current on each atomic grid, JA
    ! Per each JA, loop over the points where the shift needs to be computed
    DO ikind = 1,nkind

       NULLIFY(atom_kind,atom_list,grid_atom,harmonics)
       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
            &               atom_list=atom_list,&
            &               natom=natom,&
            &               paw_atom=paw_atom,&
            &               harmonics=harmonics,&
            &               grid_atom=grid_atom)
       !
       na = grid_atom%ng_sphere
       nr = grid_atom%nr
       !
       ! Quick cycle if needed
       IF(.NOT.paw_atom) CYCLE
       !
       ! Distribute the atoms of this kind
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit( natom, num_pe, mepos )
       !
       DO iat = bo(1),bo(2)
          iatom = atom_list(iat)
          r_iatom(1:3) = pbc(particle_set(iatom)%r,cell)
          !
          NULLIFY(jrho1_atom,jrho_h_grid,jrho_s_grid)
          jrho1_atom => jrho1_atom_set(iatom)
          !
          DO ispin = 1,nspins
             
             IF(store_current) THEN
                jrho_h_grid => jrho1_atom%jrho_vec_rad_h(idir,iB,ispin)%r_coef
                jrho_s_grid => jrho1_atom%jrho_vec_rad_s(idir,iB,ispin)%r_coef
             ELSE
                jrho_h_grid => jrho1_atom%jrho_rad_h(ispin)%r_coef
                jrho_s_grid => jrho1_atom%jrho_rad_s(ispin)%r_coef
             ENDIF

             jrho_h_grid = 0.0_dp
             jrho_s_grid = 0.0_dp
             !
             ! Compute the local (hard/soft) current density response
             CALL calculate_jrho_atom_ang(jrho1_atom,jrho_h_grid,jrho_s_grid,&
                  &                       grid_atom,harmonics,iB,idir,ispin,&
                  &                       error=error)
             !
             ! loop over the atoms neighbors of iatom in terms of the current density
             ! for each compute the contribution to the shift coming from the
             ! local current density at iatom
             DO ia = 1,na
                DO ir = 1,nr
                   !
                   rvec(1) = grid_atom%rad(ir)*harmonics%a(1,ia)
                   rvec(2) = grid_atom%rad(ir)*harmonics%a(2,ia)
                   rvec(3) = grid_atom%rad(ir)*harmonics%a(3,ia)
                   !
                   !vw Can we compute that in calculate_jrho_atom_ang? wv!
                   jrho1_tot = (jrho_h_grid(ir,ia)-jrho_s_grid(ir,ia)) * grid_atom%weight(ia,ir)
                   itegrated_jrho = itegrated_jrho+jrho1_tot
                   !
                   ! local contribution to chi (NEEDS SOME MODIFICATIONS FOR PBC)
                   DO idir2 = 1,3
                      IF(idir2/=idir) THEN
                         CALL set_vecp_rev(idir,idir2,idir3)
                         scale_fac=fac_vecp(idir3,idir2,idir)
                         nmr_env%chi_tensor_loc(idir3,iB)=nmr_env%chi_tensor_loc(idir3,iB)+&
                              & scale_fac*(jrho1_tot)*nmr_env%chi_factor*&
                              & (rvec(idir2)+particle_set(iatom)%r(idir2))
                      ENDIF
                   ENDDO
                   !
                   !
                   DO jatom = 1,natom_tot
                      ! be sure you take the minimum image
                      rij(1:3) = pbc(particle_set(iatom)%r,particle_set(jatom)%r,cell)
                      distance = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                      ! the right thing to do is to build a list of neighbors, but for the time being
                      ! let's simply use a cutoff here

! my_minimum_image = .FALSE.
! IF(PRESENT(minimum_image)) THEN
!    my_minimum_image=minimum_image
!    Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
!    Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
!    Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
! END IF
!         IF(my_minimum_image) THEN
!           IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
!             sab_orb_neighbor_node => next(sab_orb_neighbor_node)
!             CYCLE
!           END IF
!         END IF

                      threshold = 60.0_dp
                      jatomisneighbor = (distance .LE. threshold)
                      IF(jatomisneighbor) THEN

                         rdiff(1) = rij(1) - rvec(1)
                         rdiff(2) = rij(2) - rvec(2)
                         rdiff(3) = rij(3) - rvec(3)
                         ddiff = SQRT( rdiff(1)*rdiff(1)+ rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3))
                         ! the component idir of J contributes to the components idir2/=idir of the shift
                         DO idir2 = 1,3

                            IF(idir2/=idir) THEN
                               CALL set_vecp_rev(idir,idir2,idir3)
                               scale_fac=fac_vecp(idir3,idir2,idir)

                               !IF(ddiff .GT. EPSILON(0.0_dp)) THEN
                               IF(ddiff.GT.1.0E-12_dp) THEN
                                  chemical_shift_loc_tmp(idir3,jatom) = &
                                       chemical_shift_loc_tmp(idir3,jatom) + &
                                       nmr_env%shift_factor_gapw*scale_fac*rdiff(idir2)&
                                       & * jrho1_tot/(ddiff*ddiff*ddiff)
                               ENDIF

                            ENDIF
                         ENDDO ! idir2
                         
                      ENDIF

                   ENDDO ! jatom

                   IF(nmr_env%do_nics) THEN

                      DO jatom = 1,n_nics
                         r_jatom(1:3) = r_nics(1:3,jatom)
                         rij(1:3) = pbc(particle_set(iatom)%r,r_jatom,cell)

                         distance = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
                         jatomisneighbor = (distance .LE. threshold)

                         IF(jatomisneighbor) THEN

                            rdiff(1) = rij(1) - rvec(1)
                            rdiff(2) = rij(2) - rvec(2)
                            rdiff(3) = rij(3) - rvec(3)
                            ddiff = SQRT( rdiff(1)*rdiff(1)+ rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3))
                            ! the component idir of J contributes to the components idir2/=idir of the shift
                            DO idir2 = 1,3

                               IF(idir2/=idir) THEN
                                  CALL set_vecp_rev(idir,idir2,idir3)
                                  scale_fac=fac_vecp(idir3,idir2,idir)

                                  IF(ddiff .GT. EPSILON(0.0_dp)) THEN
                                     chemical_shift_loc_nics(idir3,iB,jatom) = &
                                          chemical_shift_loc_nics(idir3,iB,jatom) + &
                                          scale_fac*rdiff(idir2)*jrho1_tot/(ddiff*ddiff*ddiff)
                                  ENDIF

                               ENDIF
                            ENDDO ! idir2
                            
                         ENDIF

                      ENDDO ! jatom

                   ENDIF ! do_nics
                ENDDO ! ir
             ENDDO ! ia
          ENDDO ! ispin
       ENDDO ! iat
    ENDDO ! ikind
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E24.16)') 'Integrated local j_'&
            &//ACHAR(idir+119)//ACHAR(iB+119)//'(r)=',itegrated_jrho
    ENDIF
    !
    CALL mp_sum(chemical_shift_loc_tmp,para_env%group)
    chemical_shift_loc(:,iB,:) = chemical_shift_loc(:,iB,:) + chemical_shift_loc_tmp(:,:)
    IF(nmr_env%do_nics) THEN
       !here need a chemical_shift_loc_nics_tmp... also for chi_loc
       CALL stop_program(routineP,"local nics needs some more work")
       CALL mp_sum(chemical_shift_loc_nics,para_env%group)
    ENDIF
    !
    DEALLOCATE(chemical_shift_loc_tmp,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(0.0_dp,handle)
    !
  END SUBROUTINE shift_atom

! *****************************************************************************

  SUBROUTINE calculate_jrho_atom_ang(jrho1_atom,jrho_h,jrho_s,grid_atom,&
       &                             harmonics,iB,idir,ispin,error)
    !
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h, jrho_s
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: iB, idir, ispin
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_ang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iiB, iiiB, idir2, ir, iso, &
                                                max_iso_not0, na, nr
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rad_part,scale
    REAL(dp), DIMENSION(:,:), POINTER :: a, Fr_h, Fr_s, Fr_a_h, Fr_a_h_ii, Fr_a_h_iii, &
         Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, Fr_b_h_iii, Fr_b_s, &
         Fr_b_s_ii, Fr_b_s_iii, slm
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
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_c_h_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho1_atom%jrho_c_s_ii(ispin)%r_coef),cp_failure_level,routineP,error,failure)
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
    DO ir = 1,nr
       DO iso = 1,max_iso_not0
          DO ia = 1,na
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
          ENDDO ! ia
       ENDDO ! iso
    ENDDO ! ir
    !
  END SUBROUTINE calculate_jrho_atom_ang

! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_atom_current/list_3c_jrho_atom
!!   NAME
!!    list_3c_jrho_atom
!!
!!   FUNCTION
!!     Prepare the 3 center lists used to cpmpute the atomic terms of the current density
!!     New lists need to be constructed because the density matrices
!!     are not symmetric, hence, the ij and the ji pairs  are not equivalent and both
!!     need to be considered
!!
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2006 created [MI]
!!
!****************************************************************************

  SUBROUTINE list_3c_jrho_atom(qs_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'list_3c_jrho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, iab, iat, iatom, ikind, ilist, inode, istat, jatom, jkind, &
      last_jatom, nat, natom, nkind, nlist, nnode, nsgf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, new_atom_b
    REAL(dp)                                 :: rab(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: C_coeff_hh, C_coeff_ss
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

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(sac_oce),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(sbc_oce),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(oce),cp_failure_level,routineP,error,failure)

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         &                   nsgf=nsgf)

    ALLOCATE(first_sgf(natom),last_sgf(natom),dist_list(natom),oceh(1),oces(1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineN,error,failure)

    DO i = 1,natom
       dist_list(i)%ndist = 0
       NULLIFY(dist_list(i)%rdist)
    ENDDO
    NULLIFY(oceh(1)%block,oces(1)%block)

    CALL get_particle_set(particle_set=particle_set,&
         &                first_sgf=first_sgf,&
         &                last_sgf=last_sgf)
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
       DO iat = 1,nat
          iatom = atom_list(iat)
          CALL retrieve_sac_list(oce%voce,sac_oce,iatom,ikind,nkind,error=error)
       ENDDO

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
    ! Scale the off diagonal projector by two. 
    ! This leads to the factor half in nmr_print_shift... NEED TO FIX THAT
    DO iatom=1,natom
       DO jatom=1,natom
          NULLIFY(C_coeff_hh,C_coeff_ss)
          CALL get_block_node(matrix=oce%hardo(1)%matrix,&
               &              block_row=iatom,&
               &              block_col=jatom,&
               &              BLOCK=C_coeff_hh)
          CALL get_block_node(matrix=oce%softo(1)%matrix,&
               &              block_row=iatom,&
               &              block_col=jatom,&
               &              BLOCK=C_coeff_ss)
          IF(iatom.NE.jatom) THEN
             IF(ASSOCIATED(C_coeff_hh))C_coeff_hh=2.0_dp*C_coeff_hh
             IF(ASSOCIATED(C_coeff_ss))C_coeff_ss=2.0_dp*C_coeff_ss
          ENDIF
       ENDDO
    ENDDO
    !
    ! put the oce back to the env
    CALL set_qs_env(qs_env=qs_env,oce=oce,error=error)
    !
    DO i=1,natom
       IF(ASSOCIATED(dist_list(i)%rdist)) THEN
          DEALLOCATE(dist_list(i)%rdist,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineN,error,failure)
       ENDIF
    ENDDO
    DEALLOCATE(first_sgf,last_sgf,dist_list,oceh,oces,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineN,error,failure)
    !
  END SUBROUTINE list_3c_jrho_atom

! *****************************************************************************
END MODULE qs_linres_nmr_atom_current
