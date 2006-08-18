!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005,2006  CP2K developers group                            !
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
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_linres_nmr_op,                ONLY: set_vecp
  USE qs_linres_types,                 ONLY: allocate_jrho_atom_rad,&
                                             get_nmr_env,&
                                             jrho_atom_type,&
                                             nmr_env_type,&
                                             set2zero_jrho_atom_rad
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next, reduced_3c_list_type
  USE qs_oce_types,                    ONLY: build_reduced_3c_lists,&
                                             deallocate_vtriple,&
                                             oce_matrix_type,&
                                             retrieve_sac_list
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE


  PRIVATE

 ! *** Public subroutines ***
  PUBLIC :: calculate_jrho_atom_coeff, calculate_jrho_atom_rad, list_3c_jrho_atom

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
!!     Calculate the expansio coefficients for the atomic terms 
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

  SUBROUTINE calculate_jrho_atom_coeff(qs_env,nmr_env,mat_jp,mat_jp_pii,&
             mat_jp_piii,iB,idir,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: mat_jp, mat_jp_pii, &
                                                mat_jp_piii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_atom_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iab, iat, iatom, ii, iii, ikind, ilist, inode, isgfa, &
      isgfb, ispin, istat, jatom, jkind, katom, kkind, kneighbor, last_jatom, &
      last_katom, max_gau, max_nsgf, maxsoc, nat, natom, nkind, nlist, nnode, &
      nnode_c, nsetc, nsgfa, nsgfb, nsoctot, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, new_atom_b, paw_atom
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
      POINTER                                :: jp_block, jp_ii_block, &
                                                jp_iii_block, jp_RARnu

! *****************************************************************************

    CALL timeset("calculate_jrho_atom_coeff","I"," ",handle)

    failure =.FALSE.
    NULLIFY (atom_list,atomic_kind_set,dft_control,orb_basis_set,sab_orb)
    NULLIFY (jrho1_atom_set,oce,reduced_3c_oce)
    NULLIFY (jp_block,jp_ii_block,jp_ii_block,jp_RARnu)
    NULLIFY (C_coeff_hh_a, C_coeff_hh_b,C_coeff_ss_a, C_coeff_ss_b)
    NULLIFY (sab_orb_neighbor_list,sab_orb_neighbor_list_local,sab_orb_neighbor_node)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    oce=oce,&
                    sab_all=sab_orb)

    CALL get_nmr_env(nmr_env=nmr_env,&
                     jrho1_atom_set=jrho1_atom_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=max_nsgf,&
                             maxgtops=max_gau)

   ! Density matrices

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(jrho1_atom_set)
    nspins = dft_control%nspins
    CALL set_vecp(iB,ii,iii)

    ALLOCATE(jp_block(nspins),jp_ii_block(nspins),jp_iii_block(nspins),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(jp_RARnu(nspins),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (PC(max_nsgf,max_gau),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

   ! Reset CJC coefficients and local density arrays
    DO ikind = 1 ,nkind
      NULLIFY(atomic_kind, atom_list)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,natom=nat,&
                           paw_atom=paw_atom)

      IF(.NOT. paw_atom) CYCLE
      DO iat = 1,nat
        iatom = atom_list(iat)
        DO ispin = 1, nspins

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
          END IF

        ENDDO  ! ispin
      ENDDO  ! iat
    END DO ! ikind


   ! Three centers 
    DO ikind = 1 ,nkind
      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)
      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             nsgf=nsgfa)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nsgf=nsgfb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

        ! Allocate temporary block to store the sum of the 3 blocks before contraction
        DO ispin = 1,nspins
          NULLIFY(jp_RARnu(ispin)%r_coef)
          ALLOCATE(jp_RARnu(ispin)%r_coef(nsgfa,nsgfb),STAT=istat) 
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO

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
                                   r=rab,&
                                   reduced_3c_oce=reduced_3c_oce)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

   ! Take the block \mu \nu of jpab, jpab_ii and jpab_iii, ! 
            DO ispin = 1,nspins
               NULLIFY(jp_block(ispin)%r_coef)
               NULLIFY(jp_ii_block(ispin)%r_coef)
               NULLIFY(jp_iii_block(ispin)%r_coef)
               CALL get_block_node(matrix=mat_jp(ispin)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=jp_block(ispin)%r_coef)
               CALL get_block_node(matrix=mat_jp_pii(ispin)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=jp_ii_block(ispin)%r_coef)
               CALL get_block_node(matrix=mat_jp_piii(ispin)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=jp_iii_block(ispin)%r_coef)
            ENDDO

   ! Loop over atoms
            DO kkind=1,nkind

              atomic_kind => atomic_kind_set(kkind)

              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   orb_basis_set=orb_basis_set,&
                                   paw_atom=paw_atom)

              CALL  get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                      nset =nsetc, maxso=maxsoc)

              nsoctot = maxsoc*nsetc

              IF(.NOT. paw_atom)  CYCLE

              last_katom = 0

              nnode_c = reduced_3c_oce(kkind)%nnode
              DO kneighbor = 1,nnode_c

                katom = reduced_3c_oce(kkind)%index_atom(kneighbor)
                rbc(1:3) = reduced_3c_oce(kkind)%rbc(1:3,kneighbor)

   ! multiply by (RA-R\nu)_iii and (RA-R\nu)_ii and sum
                DO ispin = 1,nspins
                DO isgfb = 1, nsgfb
                DO isgfa = 1, nsgfa
                  jp_RARnu(ispin)%r_coef(isgfa,isgfb) =  &
                    jp_block(ispin)%r_coef(isgfa,isgfb)  & 
                   -jp_iii_block(ispin)%r_coef(isgfa,isgfb)*rbc(ii)&
                   + jp_ii_block(ispin)%r_coef(isgfa,isgfb)*rbc(iii)
                END DO 
                END DO 
                END DO

                NULLIFY(C_coeff_hh_a)
                CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=C_coeff_hh_a)

                NULLIFY(C_coeff_ss_a)
                CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=C_coeff_ss_a)


   ! De-contract and copy into the final matrix
               IF(iatom==jatom) THEN

                 DO ispin = 1,nspins

                    !First Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_RARnu(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_h(ispin)%r_coef,nsoctot)
                    !First Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_RARnu(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_ss_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_s(ispin)%r_coef,nsoctot)

                    !Second Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_ii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_ii_h(ispin)%r_coef,nsoctot)
                    !Second Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_ii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_ss_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_ii_s(ispin)%r_coef,nsoctot)


                    !Third Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_iii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_iii_h(ispin)%r_coef,nsoctot)
                    !Third Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_iii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_a,SIZE(C_coeff_ss_a,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_iii_s(ispin)%r_coef,nsoctot)
     
                 END DO  ! ispin

               ELSE

                 NULLIFY(C_coeff_hh_b)
                 CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_hh_b)

                 NULLIFY(C_coeff_ss_b)
                 CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_ss_b)


                 DO ispin = 1,nspins

                    !First Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_RARnu(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_hh_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_h(ispin)%r_coef,nsoctot)
                    !First Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_RARnu(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_ss_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_s(ispin)%r_coef,nsoctot)

                    !Second Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_ii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_hh_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_ii_h(ispin)%r_coef,nsoctot)
                    !Second Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_ii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_ss_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_ii_s(ispin)%r_coef,nsoctot)


                    !Third Matrix hard
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_iii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_hh_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_iii_h(ispin)%r_coef,nsoctot)
                    !Third Matrix soft
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                               jp_iii_block(ispin)%r_coef,nsgfa,&
                               C_coeff_hh_b,SIZE(C_coeff_ss_b,1),&
                               0.0_dp,PC,max_nsgf)
                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                               SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                               jrho1_atom_set(katom)%cjc_iii_s(ispin)%r_coef,nsoctot)

                 END DO  ! ispin

               END IF


              END DO  ! kneighbor
            END DO  ! kkind
          sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO  ! inode

        END DO  ! ilist
        DO ispin = 1,nspins
          DEALLOCATE(jp_RARnu(ispin)%r_coef,STAT=istat) 
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
        END DO

      END DO ! jkind 
    END DO ! ikind 

    DO iatom = 1,natom

       DO ispin = 1,nspins
         IF(ASSOCIATED(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)) THEN
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef,qs_env%para_env%group)
         END IF
       ENDDO

    END DO

    DEALLOCATE (PC,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (jp_block,jp_ii_block,jp_iii_block,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (jp_RARnu,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_jrho_atom_coeff

! *****************************************************************************

  SUBROUTINE calculate_jrho_atom_rad(qs_env,nmr_env,IB,idir,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i1, i2, iat, iatom, icg, ikind, ipgf1, ipgf2, iset1, &
      iset2, iso, iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, &
      ispin, istat, l, l_iso, lmax12, lmax_expansion, lmin12, m1s, m2s, &
      max_iso_not0, max_nso, maxl, maxso, mepos, n1s, n2s, natom, nkind, nr, &
      nset, nspins, num_pe, size1, size2
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list, lmax, lmin, npgf
    LOGICAL                                  :: failure, paw_atom
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: cjc_h_block, cjc_ii_h_block, &
      cjc_ii_s_block, cjc_iii_h_block, cjc_iii_s_block, cjc_s_block, dgg_1, &
      dgg_2, gg, gg_lm1
    REAL(dp), DIMENSION(:, :), POINTER :: coeff, Fr_a_h, Fr_a_h_ii, &
      Fr_a_h_iii, Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, &
      Fr_b_h_iii, Fr_b_s, Fr_b_s_ii, Fr_b_s_iii, Fr_c_h_ii, Fr_c_h_iii, &
      Fr_c_s_ii, Fr_c_s_iii, zet
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

    CALL timeset("calculate_jrho_atom_rad","I"," ",handle)

    failure =.FALSE.
    NULLIFY (atomic_kind_set,dft_control,para_env)
    NULLIFY (coeff,Fr_a_h, Fr_a_s,Fr_a_h_ii, Fr_a_s_ii,Fr_a_h_iii, Fr_a_s_iii)
    NULLIFY (Fr_b_h, Fr_b_s,Fr_b_h_ii, Fr_b_s_ii,Fr_b_h_iii, Fr_b_s_iii)
    NULLIFY (Fr_c_h_ii, Fr_c_s_ii,Fr_c_h_iii, Fr_c_s_iii)
    NULLIFY (jrho1_atom_set,jrho1_atom)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    para_env=para_env)

    CALL get_nmr_env(nmr_env=nmr_env,&
                     jrho1_atom_set=jrho1_atom_set)


    nkind = SIZE(atomic_kind_set)
    nspins = dft_control%nspins

    DO ikind = 1,nkind
      NULLIFY (atom_kind,atom_list,grid_atom,harmonics,orb_basis_set)
      NULLIFY (lmax, lmin, npgf, zet) 
      NULLIFY (grid_atom,harmonics)
      NULLIFY (my_CG,my_CG_dxyz_asym)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                           grid_atom=grid_atom,natom=natom,&
                           paw_atom=paw_atom, harmonics=harmonics,&
                           orb_basis_set=orb_basis_set)

      IF(.NOT. paw_atom) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             lmax=lmax,lmin=lmin,&
                             maxl=maxl, npgf=npgf, nset=nset, zet=zet,&
                             maxso=maxso)

      max_iso_not0 = harmonics%max_iso_not0
      nr           = grid_atom%nr
      lmax_expansion = indso(1,max_iso_not0)
!     Distribute the atoms of this kind
      num_pe = para_env%num_pe
      mepos  = para_env%mepos
      bo = get_limit( natom, num_pe, mepos )

      my_CG           => harmonics%my_CG
      my_CG_dxyz_asym => harmonics%my_CG_dxyz_asym


      ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),gg_lm1(nr,0:2*maxl),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dgg_1(nr,0:2*maxl),dgg_2(nr,0:2*maxl),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

      max_nso = nsoset(maxl)
      ALLOCATE (cjc_h_block(max_nso,max_nso),cjc_s_block(max_nso,max_nso),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (cjc_ii_h_block(max_nso,max_nso),cjc_ii_s_block(max_nso,max_nso),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE (cjc_iii_h_block(max_nso,max_nso),cjc_iii_s_block(max_nso,max_nso),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

      m1s = 0
      DO iset1 = 1,nset
        n1s = nsoset(lmax(iset1))
        DO ipgf1  = 1,npgf(iset1)

          iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
          iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
          size1 = iso1_last - iso1_first + 1
          i1 = nsoset(lmin(iset1)-1)+1

          g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))

          m2s = 0
          DO iset2 = 1,nset
            n2s=nsoset(lmax(iset2))
            DO ipgf2 = 1,npgf(iset2)
              iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
              iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
              size2 = iso2_last - iso2_first + 1
              i2 = nsoset(lmin(iset2)-1)+1

              g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))

              lmin12 = lmin(iset1)+lmin(iset2)
              lmax12 = lmax(iset1)+lmax(iset2)

              gg = 0.0_dp
              gg_lm1 = 0.0_dp
              dgg_1 = 0.0_dp
              dgg_2 = 0.0_dp

              ! Take only the terms of expansion for L < lmax_expansion
              IF(lmin12 .LE. lmax_expansion) THEN

                IF(lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

                IF (lmin12 == 0) THEN
                  gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
                  gg_lm1(1:nr,lmin12) = 0.0_dp
                ELSE
                  gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
                  gg_lm1(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12-1)*g1(1:nr)*g2(1:nr)
                END IF

                DO l=lmin12+1,lmax12
                  gg(1:nr,l) = grid_atom%rad(1:nr)*gg(1:nr,l-1)
                  gg_lm1(1:nr,l) = gg(1:nr,l-1)
                  dgg_1(1:nr,l-1) = 2.0_dp*(zet(ipgf1,iset1) - zet(ipgf2,iset2))*&
                                    gg(1:nr,l)
                  IF(l-2>=lmin12) THEN
                    dgg_2(1:nr,l-2) = 2.0_dp*(zet(ipgf1,iset1)- zet(ipgf2,iset2))* & 
                                      gg(1:nr,l)
                  END IF
                END DO
                dgg_1(1:nr,lmax12)   = 2.0_dp*(zet(ipgf1,iset1)- zet(ipgf2,iset2))*&
                                       grid_atom%rad(1:nr) * gg(1:nr,lmax12)
                dgg_2(1:nr,lmax12-1) = 2.0_dp*(zet(ipgf1,iset1)- zet(ipgf2,iset2))*&
                                       grid_atom%rad(1:nr) * gg(1:nr,lmax12)
                dgg_2(1:nr,lmax12)   = dgg_2(1:nr,lmax12-1) * grid_atom%rad(1:nr)

              ELSE
                CYCLE
              END IF  ! lmin12

              DO iat = bo(1),bo(2)
                iatom = atom_list(iat)

                DO ispin=1,nspins

                  cjc_h_block  = 0.0_dp
                  cjc_s_block  = 0.0_dp
                  coeff => jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef
                  cjc_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                  coeff => jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef
                  cjc_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)

                  cjc_ii_h_block = 0.0_dp
                  cjc_ii_s_block = 0.0_dp
                  coeff => jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef
                  cjc_ii_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                  coeff => jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef
                  cjc_ii_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)

                  cjc_iii_h_block = 0.0_dp
                  cjc_iii_s_block = 0.0_dp
                  coeff => jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef
                  cjc_iii_h_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                  coeff => jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef
                  cjc_iii_s_block(i1:i1+size1-1,i2:i2+size2-1) = &
                             coeff(iso1_first:iso1_last,iso2_first:iso2_last)

                  jrho1_atom => jrho1_atom_set(iatom)
                  ! Allocation radial functions
                  IF(.NOT. ASSOCIATED(jrho1_atom_set(iatom)%jrho_a_h(ispin)%r_coef)) THEN
                    CALL allocate_jrho_atom_rad(jrho1_atom,ispin,nr,max_iso_not0, error=error)
                  ELSE
                    CALL set2zero_jrho_atom_rad(jrho1_atom,ispin,error=error)
                  END IF

                  Fr_a_h => jrho1_atom%jrho_a_h(ispin)%r_coef
                  Fr_a_s => jrho1_atom%jrho_a_s(ispin)%r_coef
                  Fr_b_h => jrho1_atom%jrho_b_h(ispin)%r_coef
                  Fr_b_s => jrho1_atom%jrho_b_s(ispin)%r_coef

                  Fr_a_h_ii => jrho1_atom%jrho_a_h_ii(ispin)%r_coef
                  Fr_a_s_ii => jrho1_atom%jrho_a_s_ii(ispin)%r_coef
                  Fr_b_h_ii => jrho1_atom%jrho_b_h_ii(ispin)%r_coef
                  Fr_b_s_ii => jrho1_atom%jrho_b_s_ii(ispin)%r_coef
                  Fr_c_h_ii => jrho1_atom%jrho_c_h_ii(ispin)%r_coef
                  Fr_c_s_ii => jrho1_atom%jrho_c_s_ii(ispin)%r_coef

                  Fr_a_h_iii => jrho1_atom%jrho_a_h_iii(ispin)%r_coef
                  Fr_a_s_iii => jrho1_atom%jrho_a_s_iii(ispin)%r_coef
                  Fr_b_h_iii => jrho1_atom%jrho_b_h_iii(ispin)%r_coef
                  Fr_b_s_iii => jrho1_atom%jrho_b_s_iii(ispin)%r_coef
                  Fr_c_h_iii => jrho1_atom%jrho_c_h_iii(ispin)%r_coef
                  Fr_c_s_iii => jrho1_atom%jrho_c_s_iii(ispin)%r_coef

                  DO iso = 1,max_iso_not0
                    l_iso = indso(1,iso)
                    DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                      iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                      iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)

                      l = indso(1,iso1) + indso(1,iso2)
                      CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)

                      Fr_a_h(1:nr,iso) = Fr_a_h(1:nr,iso) + &
                           dgg_1(1:nr,l)*cjc_h_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                      Fr_a_s(1:nr,iso) = Fr_a_s(1:nr,iso) + &
                           dgg_1(1:nr,l)*cjc_s_block(iso1,iso2)*my_CG(iso1,iso2,iso)

                      Fr_a_h_ii(1:nr,iso) = Fr_a_h_ii(1:nr,iso) + &
                           dgg_2(1:nr,l)*cjc_ii_h_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                      Fr_a_s_ii(1:nr,iso) = Fr_a_s_ii(1:nr,iso) + &
                           dgg_2(1:nr,l)*cjc_ii_s_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                     
                      Fr_c_h_ii(1:nr,iso) = Fr_c_h_ii(1:nr,iso) + &
                           gg(1:nr,l)*cjc_ii_h_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                      Fr_c_s_ii(1:nr,iso) = Fr_c_s_ii(1:nr,iso) + &
                           gg(1:nr,l)*cjc_ii_s_block(iso1,iso2)*my_CG(iso1,iso2,iso)

                      Fr_a_h_iii(1:nr,iso) = Fr_a_h_iii(1:nr,iso) + &
                           dgg_2(1:nr,l)*cjc_iii_h_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                      Fr_a_s_iii(1:nr,iso) = Fr_a_s_iii(1:nr,iso) + &
                           dgg_2(1:nr,l)*cjc_iii_s_block(iso1,iso2)*my_CG(iso1,iso2,iso)

                      Fr_c_h_iii(1:nr,iso) = Fr_c_h_iii(1:nr,iso) + &
                           gg(1:nr,l)*cjc_iii_h_block(iso1,iso2)*my_CG(iso1,iso2,iso)
                      Fr_c_s_iii(1:nr,iso) = Fr_c_s_iii(1:nr,iso) + &
                           gg(1:nr,l)*cjc_iii_s_block(iso1,iso2)*my_CG(iso1,iso2,iso)

                    END DO  !icg

                    DO icg = 1,harmonics%ncg_dxyz_asym(iso,iset1,iset2)
                      iso1 = harmonics%ind_cg_a_dxyz(icg,iso,iset1,iset2)
                      iso2 = harmonics%ind_cg_b_dxyz(icg,iso,iset1,iset2)

                      l = indso(1,iso1) + indso(1,iso2)


                      Fr_b_h(1:nr,iso) = Fr_b_h(1:nr,iso) + &
                                       gg_lm1(1:nr,l)*cjc_h_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)
                      Fr_b_s(1:nr,iso) = Fr_b_s(1:nr,iso) + &
                                       gg_lm1(1:nr,l)*cjc_s_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)

                      Fr_b_h_ii(1:nr,iso) = Fr_b_h_ii(1:nr,iso) + &
                                       gg(1:nr,l)*cjc_ii_h_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)
                      Fr_b_s_ii(1:nr,iso) = Fr_b_s_ii(1:nr,iso) + &
                                       gg(1:nr,l)*cjc_ii_s_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)

                      Fr_b_h_iii(1:nr,iso) = Fr_b_h_iii(1:nr,iso) + &
                                       gg(1:nr,l)*cjc_iii_h_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)
                      Fr_b_s_iii(1:nr,iso) = Fr_b_s_iii(1:nr,iso) + &
                                       gg(1:nr,l)*cjc_iii_s_block(iso1,iso2)*&
                                       my_CG_dxyz_asym(idir,iso1,iso2,iso)

                    END DO  ! icg

                  END DO  ! iso

                  NULLIFY(Fr_a_h,Fr_a_s,Fr_a_h_ii,Fr_a_s_ii,Fr_a_h_iii,Fr_a_s_iii)
                  NULLIFY(Fr_b_h,Fr_b_s,Fr_b_h_ii,Fr_b_s_ii,Fr_b_h_iii,Fr_b_s_iii)
                  NULLIFY(Fr_c_h_ii,Fr_c_s_ii,Fr_c_h_iii,Fr_c_s_iii)

                END DO  ! ispin

              END DO  ! iat

            END DO  !ipgf2
            m2s = m2s+maxso
          END DO  ! iset2
        END DO  ! ipgf1
        m1s = m1s+maxso
      END DO  ! iset1

      DEALLOCATE (cjc_h_block,cjc_s_block,STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE (cjc_ii_h_block,cjc_ii_s_block,STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE (cjc_iii_h_block,cjc_iii_s_block,STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(g1,g2,gg,gg_lm1,dgg_1,dgg_2,STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO  ! ikind

    CALL timestop(0.0_dp,handle)


  END SUBROUTINE calculate_jrho_atom_rad

! *****************************************************************************

  SUBROUTINE shift_atom(qs_env,nmr_env,iB,idir,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(IN)                      :: IB, idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'shift_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iat, iatom, ikind, &
                                                ispin, mepos, natom, nkind, &
                                                nspins, num_pe
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom, &
                                                store_current
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h_grid, jrho_s_grid
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom

    CALL timeset("shift_atom","I"," ",handle)

    failure =.FALSE.
    NULLIFY (atomic_kind_set,dft_control,para_env)
    NULLIFY (jrho1_atom_set,jrho1_atom)
    NULLIFY (jrho_h_grid, jrho_s_grid)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    para_env=para_env)

    CALL get_nmr_env(nmr_env=nmr_env,&
                     jrho1_atom_set=jrho1_atom_set)


    nkind = SIZE(atomic_kind_set)
    nspins = dft_control%nspins
    store_current = nmr_env%store_current

    ! Loop over atoms to collocate the current on each atomic grid, JA
    ! Per each JA, loop over the points where the shift needs to be computed

    DO ikind = 1,nkind

      NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                           natom=natom, paw_atom=paw_atom,&
                           harmonics=harmonics,grid_atom=grid_atom)

      IF (.NOT.paw_atom) CYCLE

      ! Distribute the atoms of this kind

      num_pe = para_env%num_pe
      mepos  = para_env%mepos
      bo = get_limit( natom, num_pe, mepos )

      DO iat = bo(1),bo(2) 
        iatom = atom_list(iat)

        NULLIFY (jrho1_atom,jrho_h_grid,jrho_s_grid)
        jrho1_atom => jrho1_atom_set(iatom)

        DO ispin = 1,nspins

          IF(store_current) THEN
            jrho_h_grid => jrho1_atom%jrho_vec_rad_h(idir,iB,ispin)%r_coef
            jrho_s_grid => jrho1_atom%jrho_vec_rad_s(idir,iB,ispin)%r_coef
          ELSE
            jrho_h_grid => jrho1_atom%jrho_rad_h(ispin)%r_coef
            jrho_s_grid => jrho1_atom%jrho_rad_s(ispin)%r_coef
          END IF

          jrho_h_grid = 0.0_dp
          jrho_s_grid = 0.0_dp
          CALL calculate_jrho_atom_ang(jrho1_atom,jrho_h_grid,jrho_s_grid,&
               grid_atom,harmonics,iB,idir,ispin,error=error)

       ! And now compute the shift somehow ......... if you can


        END DO 

      END DO  ! iat

    END DO  ! ikind

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE shift_atom

  SUBROUTINE calculate_jrho_atom_ang(jrho1_atom,jrho_h,jrho_s,grid_atom,&
             harmonics,iB,idir,ispin,error)

    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    REAL(dp), DIMENSION(:, :), POINTER       :: jrho_h, jrho_s
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: iB, idir, ispin
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_jrho_atom_ang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iiB, iiiB, ir, iso, &
                                                max_iso, na, nr
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rad_part
    REAL(dp), DIMENSION(:, :), POINTER :: Fr_a_h, Fr_a_h_ii, Fr_a_h_iii, &
      Fr_a_s, Fr_a_s_ii, Fr_a_s_iii, Fr_b_h, Fr_b_h_ii, Fr_b_h_iii, Fr_b_s, &
      Fr_b_s_ii, Fr_b_s_iii, Fr_c_h_ii, Fr_c_h_iii, Fr_c_s_ii, Fr_c_s_iii

    failure = .FALSE.
    NULLIFY (Fr_a_h, Fr_a_s,Fr_a_h_ii, Fr_a_s_ii,Fr_a_h_iii, Fr_a_s_iii)
    NULLIFY (Fr_b_h, Fr_b_s,Fr_b_h_ii, Fr_b_s_ii,Fr_b_h_iii, Fr_b_s_iii)
    NULLIFY (Fr_c_h_ii, Fr_c_s_ii,Fr_c_h_iii, Fr_c_s_iii)

    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(jrho_s),cp_failure_level,routineP,error,failure)

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      max_iso = harmonics%max_iso_not0

      Fr_a_h => jrho1_atom%jrho_a_h(ispin)%r_coef
      Fr_a_s => jrho1_atom%jrho_a_s(ispin)%r_coef
      Fr_b_h => jrho1_atom%jrho_b_h(ispin)%r_coef
      Fr_b_s => jrho1_atom%jrho_b_s(ispin)%r_coef
  
      Fr_a_h_ii => jrho1_atom%jrho_a_h_ii(ispin)%r_coef
      Fr_a_s_ii => jrho1_atom%jrho_a_s_ii(ispin)%r_coef
      Fr_b_h_ii => jrho1_atom%jrho_b_h_ii(ispin)%r_coef
      Fr_b_s_ii => jrho1_atom%jrho_b_s_ii(ispin)%r_coef
      Fr_c_h_ii => jrho1_atom%jrho_c_h_ii(ispin)%r_coef
      Fr_c_s_ii => jrho1_atom%jrho_c_s_ii(ispin)%r_coef
  
      Fr_a_h_iii => jrho1_atom%jrho_a_h_iii(ispin)%r_coef
      Fr_a_s_iii => jrho1_atom%jrho_a_s_iii(ispin)%r_coef
      Fr_b_h_iii => jrho1_atom%jrho_b_h_iii(ispin)%r_coef
      Fr_b_s_iii => jrho1_atom%jrho_b_s_iii(ispin)%r_coef
      Fr_c_h_iii => jrho1_atom%jrho_c_h_iii(ispin)%r_coef
      Fr_c_s_iii => jrho1_atom%jrho_c_s_iii(ispin)%r_coef

      CALL set_vecp(iB,iiB,iiiB)

      DO ir = 1,nr
        DO iso = 1,max_iso
          DO ia = 1,na

            rad_part = harmonics%a(idir,ia)*Fr_a_h(ir,iso)+Fr_b_h(ir,iso)-&
                       harmonics%a(iiB,ia)*&
                       (harmonics%a(idir,ia)*Fr_a_h_iii(ir,iso)+Fr_b_h_iii(ir,iso))+&
                       harmonics%a(iiiB,ia)*&
                       (harmonics%a(idir,ia)*Fr_a_h_ii(ir,iso)+Fr_b_h_ii(ir,iso))
            IF(idir==iiB) THEN
              rad_part = rad_part - Fr_c_h_iii(ir,iso)
            ELSEIF(idir==iiiB) THEN
              rad_part = rad_part + Fr_c_h_ii(ir,iso)
            END IF
                       
            jrho_h(ir,ia) = jrho_h(ir,ia) + rad_part * harmonics%slm(ia,iso)

            rad_part = harmonics%a(idir,ia)*Fr_a_s(ir,iso)+Fr_b_s(ir,iso)-&
                       harmonics%a(iiB,ia)*&
                       (harmonics%a(idir,ia)*Fr_a_s_iii(ir,iso)+Fr_b_s_iii(ir,iso))+&
                       harmonics%a(iiiB,ia)*&
                       (harmonics%a(idir,ia)*Fr_a_s_ii(ir,iso)+Fr_b_s_ii(ir,iso))
            IF(idir==iiB) THEN
              rad_part = rad_part - Fr_c_s_iii(ir,iso)
            ELSEIF(idir==iiiB) THEN
              rad_part = rad_part + Fr_c_s_ii(ir,iso)
            END IF
                       
            jrho_s(ir,ia) = jrho_s(ir,ia) + rad_part * harmonics%slm(ia,iso)

          END DO  ! ia
        END DO  ! iso
      
      END DO  ! ir 

  END SUBROUTINE calculate_jrho_atom_ang

! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_atom_current/list_3c_jrho_atom
!!   NAME
!!    list_3c_jrho_atom 
!!
!!   FUNCTION
!!     Prepare the 3 center lists used tocpmpute the atomic terms of the current density
!!     New lists need to be constructed because the density matrices
!!     are not symmetric, hence, the ij and the ji pairs  are not equivalens and bpth
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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'list_3c_jrho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iab, iat, iatom, ikind, &
                                                ilist, inode, jatom, jkind, &
                                                last_jatom, nat, natom, &
                                                nkind, nlist, nnode
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, new_atom_b
    REAL(dp)                                 :: rab(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_oce, sbc_oce
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce

! *****************************************************************************

    NULLIFY (sac_oce, sbc_oce, oce)
    NULLIFY (atomic_kind,atomic_kind_set)
    NULLIFY (sab_orb_neighbor_list,sab_orb_neighbor_list_local)
    NULLIFY (sab_orb_neighbor_node,reduced_3c_oce) 

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    sab_all=sab_orb,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce,&
                    oce=oce)

    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(sac_oce),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(sbc_oce),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(oce),cp_failure_level,routineP,error,failure)

    nkind = SIZE(atomic_kind_set)

!   *** Retrieve data of SAC_OCE and of SAC_3C neighbors of atom "iatom" ***

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=nat,atom_list=atom_list)
      DO iat = 1,nat
        iatom = atom_list(iat)
        CALL retrieve_sac_list(oce%voce,sac_oce,iatom,ikind,nkind,error=error)
      END DO

    END DO

    DO ikind=1,nkind

      DO jkind=1,nkind

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
                                   r=rab)

            IF (jatom /= last_jatom) THEN
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

! This should not be necessary
!            *** Calculate the coefficients of the one center expansions ***
!                CALL set_up_oce(atomic_kind_set,reduced_3c_oce,oce,&
!                              ikind,iatom,jkind,jatom,rab,&
!                              oceh,oces,nder,dist_list,error=error)
!

          END DO  ! inode

        END DO  ! ilist

      END DO  ! jkind

    END DO  ! ikind

!   the vlist can be deallocated already here
    CALL deallocate_vtriple(oce%voce)

  END SUBROUTINE list_3c_jrho_atom

! *****************************************************************************
END MODULE qs_linres_nmr_atom_current
