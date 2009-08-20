!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE dkh_main 
  USE ai_kinetic,                      ONLY: kinetic
  USE ai_overlap_new,                  ONLY: overlap
  USE ai_overlap_ppl,                  ONLY: overlap_ppl
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_dbcsr_operations,             ONLY: &
       dbcsr_add_block_node, dbcsr_allocate_matrix, &
       dbcsr_allocate_matrix_set, dbcsr_deallocate_matrix_set, dbcsr_from_sm, &
       dbcsr_get_block_node, dbcsr_obj_p_type, &
       dbcsr_replicate_matrix_structure, sm_from_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_syrk,&
                                             cp_fm_transpose,&
                                             cp_fm_triangular_multiply,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm
  USE dbcsr_operations,                ONLY: dbcsr_add
  USE dbcsr_util,                      ONLY: dbcsr_finalize
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: rel_pot_erfc,&
                                             rel_pot_full,&
                                             rel_trans_atom,&
                                             rel_trans_full
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE physcon,                         ONLY: a_fine
  USE qs_all_potential,                ONLY: all_integrals,&
                                             vall_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_neighbor_list_types,          ONLY: &
       find_neighbor_list, first_list, first_node, get_neighbor_list, &
       get_neighbor_list_set, get_neighbor_node, neighbor_list_set_p_type, &
       neighbor_list_type, neighbor_node_type, next, reduced_3c_list_type
  USE qs_oce_types,                    ONLY: allocate_vtriple,&
                                             deallocate_vtriple,&
                                             oce_matrix_type,&
                                             retrieve_sac_list,&
                                             vtriple_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  TYPE vppl_type
     REAL(KIND = dp)                          :: alpha_ppl,ppl_radius
     REAL(KIND = dp), DIMENSION(:), POINTER   :: cexp_ppl,r2
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: r
     INTEGER, DIMENSION(:), POINTER           :: neighbor
  END TYPE vppl_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dkh_main'


  PUBLIC                               :: dkh_mol_integrals,&
                                          dkh_atom_transformation

CONTAINS

! *****************************************************************************
!> \brief Calculation of all integrals needed for DKH transformation 
!> \author Creation (11.03.2002,MK)
!>      Non-redundant calculation of the non-local part of the GTH PP (22.05.2003,MK)
!>      New parallelization scheme (27.06.2003,MK) 
!>      contracted to compute only uncontracted integrals for DKH trafo (17.10.2008,JT)
! *****************************************************************************
  SUBROUTINE dkh_mol_integrals(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dkh_mol_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iac, iat, iatom, iblock_col, iblock_row, ico, &
      icol, ikind, ilist, inode, irow, iset, ithread, j, jatom, jco, jkind, &
      jset, k, kkind, kneighbor, last_jatom, ldai, ldsab, maxblock, maxco, &
      maxdco, maxder, maxl, maxlgto, maxlppl, maxlppnl, maxppnl, maxsgf, n, &
      nat, natom, ncoa, ncoa_sum, ncob, nder, neighbor_list_id, nkind, nlist, &
      nneighbor, nnode, nrow, nseta, nsetb, nsgf, nthread, pgfa, pgfb, sgfa, &
      sgfb, stat, tot_npgfa, tot_npgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_pgf, first_sgf, &
                                                last_pgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                lb_max, lb_min, npgfa, npgfb, &
                                                nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: all_potential_present, dkh_erfc, dkh_full, dkh_loc_atom, &
      failure, gth_potential_present, kinetic_m, new_atom_b, oce_present, &
      ppl_present, return_s_derivatives, return_t_derivatives, rho0_present
    REAL(KIND=dp)                            :: dab, dac, dbc, f0, rab2, &
                                                rac2, rbc2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: rab, rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: hab, hp_block, hpgf_block, &
      loc_hp_block, loc_n_block, loc_pVp_block, loc_s_block, loc_t_block, &
      n_block, nab, pVp_block, pVpab, rpgfa, rpgfb, s_block, sab, spgf_block, &
      sphi_a, sphi_b, t_block, tab, tpgf_block, work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: habt, nabt, pabt, pVpabt, &
                                                sabt, sdab, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: sdabt
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER                :: matrix_pVp_full, &
                                                matrix_s_full, matrix_t_full, &
                                                matrix_v_full
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), POINTER :: matrix_h, matrix_hpgf, &
      matrix_pVppgf, matrix_s, matrix_spgf, matrix_t, matrix_tpgf, matrix_vnuc
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_3c, sac_ppl, &
                                                sbc_3c
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local, &
                                                sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sac_ppl_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_h_sm, &
      matrix_hpgf_sm, matrix_pVppgf_sm, matrix_s_sm, matrix_spgf_sm, &
      matrix_t_sm, matrix_tpgf_sm
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce, &
                                                reduced_3c_rho0
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: input
    TYPE(vall_type), DIMENSION(:), POINTER   :: vall
    TYPE(vall_type), DIMENSION(:, :), &
      POINTER                                :: vallt
    TYPE(virial_type), POINTER               :: virial
    TYPE(vppl_type), DIMENSION(:), POINTER   :: vppl
    TYPE(vppl_type), DIMENSION(:, :), &
      POINTER                                :: vpplt
    TYPE(vtriple_type), DIMENSION(:, :), &
      POINTER                                :: VH_3c_list

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY (atomic_kind_set)
    NULLIFY (distribution_2d)
    NULLIFY (matrix_h_sm)
    NULLIFY (habt)
    NULLIFY (matrix_hpgf_sm)               
    NULLIFY (nabt)
    NULLIFY (pVpabt)
    NULLIFY (VH_3c_list)
    NULLIFY (pabt)
    NULLIFY (particle_set)
    NULLIFY (matrix_pVppgf_sm)              
    NULLIFY (matrix_s_sm)
    NULLIFY (matrix_s)
    NULLIFY (matrix_h)
    NULLIFY (matrix_spgf)
    NULLIFY (matrix_tpgf)
    NULLIFY (matrix_hpgf)
    NULLIFY (matrix_vnuc)
    NULLIFY (matrix_pVppgf)
    NULLIFY (matrix_spgf_sm)                 
    NULLIFY (sab_orb)
    NULLIFY (sabt)
    NULLIFY (sac_ppl)
    NULLIFY (sac_3c)
    NULLIFY (sbc_3c)
    NULLIFY (scf_control)
    NULLIFY (sdabt)
    NULLIFY (matrix_t_sm)
    NULLIFY (matrix_t)
    NULLIFY (matrix_tpgf_sm)                  
    NULLIFY (matrix_vnuc)                    
    NULLIFY (workt)
    NULLIFY (para_env)
    NULLIFY (atom_list)

    ! Default
    oce_present = .FALSE.
    rho0_present = .FALSE.
    kinetic_m = .FALSE.
    return_s_derivatives = .FALSE.

    IF(qs_env%rel_control%rel_potential==rel_pot_erfc)THEN
       dkh_erfc=.TRUE.
    ELSE
       dkh_erfc=.FALSE.
    END IF
    IF(qs_env%rel_control%rel_DKH_transformation==rel_trans_atom)THEN
       dkh_loc_atom = .TRUE.
    ELSE
       dkh_loc_atom = .FALSE.
    END IF
    IF(qs_env%rel_control%rel_DKH_transformation==rel_trans_full)THEN
       dkh_full = .TRUE.
    ELSE
       dkh_full = .FALSE.
    END IF

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         matrix_h=matrix_h_sm,&
         kinetic=matrix_t_sm,&
         oce=oce, &
         particle_set=particle_set,&
         matrix_s=matrix_s_sm,&
         sab_orb=sab_orb,&
         sac_ppl=sac_ppl,&
         VH_3c_list=VH_3c_list,&
         sac_3c=sac_3c,&
         sbc_3c=sbc_3c,&
         neighbor_list_id=neighbor_list_id,&
         scf_control=scf_control,&
         distribution_2d=distribution_2d,&
         virial=virial,&
         input=input,para_env=para_env,error=error)

       CALL dbcsr_allocate_matrix_set( matrix_s, SIZE(matrix_s_sm), error )
       CALL dbcsr_allocate_matrix_set( matrix_h, SIZE(matrix_h_sm), error )
       IF(kinetic_m) CALL dbcsr_allocate_matrix_set( matrix_t, SIZE(matrix_t_sm), error )

       DO i = 1,SIZE(matrix_s)
          ALLOCATE(matrix_s(i)%matrix)!super dirty
          ALLOCATE(matrix_h(i)%matrix)!super dirty
          CALL dbcsr_from_sm(matrix_s(i)%matrix, matrix_s_sm(i)%matrix, error)
          CALL dbcsr_from_sm(matrix_h(i)%matrix, matrix_h_sm(i)%matrix, error)
          IF(kinetic_m) THEN
             ALLOCATE(matrix_t(i)%matrix)!super dirty
             CALL dbcsr_from_sm(matrix_t(i)%matrix, matrix_t_sm(i)%matrix, error)
          ENDIF
       ENDDO


    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    nder = 0 
    return_t_derivatives = .FALSE.
    maxder = ncoset(nder)

    !!  ********************************************************************
    !!  * Allocate the basic matrices matrix_spgf, matrix_s and construct
    !!  * the full_matrix type needed for the DKH transformation 
    !!  ******************************************************************** 

    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (first_pgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_pgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
         first_pgf=first_pgf,&          
         last_pgf=last_pgf,&           
         first_sgf=first_sgf,&
         last_sgf=last_sgf,&
         error = error )
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = qs_env%blacs_env,&
         nrow_global = last_pgf(natom),&
         ncol_global = last_pgf(natom),&
         error = error)

    !   assert here that the S-matrix is of the type distribution_2d if this one is non-null
    !   we assert that the sparsity of s is determined by the neighbor list that is used to build it
    !   is that right ? i.e. this implies that even if forces are computed the sparsity of s remains the same.
    !   this seems necessary anyway.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    CALL dbcsr_allocate_matrix_set(matrix_spgf,1,error=error)
    CALL dbcsr_allocate_matrix(matrix=matrix_spgf(1)%matrix,&
         nblock_row=natom,&
         nblock_col=natom,&
         nrow=nsgf,& 
         ncol=nsgf,&
         first_row=first_pgf,&                       
         last_row=last_pgf,&
         first_col=first_pgf,&
         last_col=last_pgf,&
         matrix_name="OVERLAP MATRIX PGF",&
         matrix_symmetry="symmetric",&
         sparsity_id=neighbor_list_id, &
         distribution_2d=distribution_2d,&
         error=error)

    DEALLOCATE (first_pgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_pgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    !   *********************************************************************
    !   *** Construct all other matrices ( uncontracted GTO ) ***************
    !   *********************************************************************

    CALL dbcsr_allocate_matrix_set(matrix_tpgf,1,error=error)             
    CALL dbcsr_replicate_matrix_structure(source=matrix_spgf(1)%matrix,&
         TARGET=matrix_tpgf(1)%matrix,&
         target_name="KINETIC ENERGY MATRIX PGF",&
         error=error)


    CALL dbcsr_allocate_matrix_set(matrix_hpgf,1,error=error)
    CALL dbcsr_replicate_matrix_structure(source=matrix_spgf(1)%matrix,&
         TARGET=matrix_hpgf(1)%matrix,&
         target_name="CORE HAMILTONIAN MATRIX PGF",&
         error=error)


    CALL dbcsr_allocate_matrix_set(matrix_vnuc,1,error=error)                                   
    CALL dbcsr_replicate_matrix_structure(source=matrix_spgf(1)%matrix,&            
         TARGET=matrix_vnuc(1)%matrix,&            
         target_name="NUCLEAR POTENTIAL ENERGY MATRIX",&
         error=error)


    CALL dbcsr_allocate_matrix_set(matrix_pVppgf,1,error=error)
    CALL dbcsr_replicate_matrix_structure(source=matrix_spgf(1)%matrix,&
         TARGET=matrix_pVppgf(1)%matrix,&
         target_name="pVp MATRIX PGF",&
         error=error)

    !   *********************************************************************
    !   ***  Allocate work storage                            ***************
    !   *********************************************************************

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present,&
         all_potential_present=all_potential_present, &
         maxco=maxco,&
         maxlgto=maxlgto,&
         maxlppl=maxlppl,&
         maxlppnl=maxlppnl,&
         maxppnl=maxppnl,&
         maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl,maxlppnl)
    CALL init_orbital_pointers(maxl+nder+1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxdco,&
         maxder=MAX(1,nder))

    ldsab = MAX(maxco,ncoset(maxlppnl),maxsgf,maxppnl)
    maxblock = MAX(SIZE(matrix_spgf),maxder)
    nthread = 1
    !   $  nthread = omp_get_max_threads()
    ldai = ncoset(maxl+nder+1)

    CALL reallocate(habt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    CALL reallocate(sabt,1,ldsab,1,ldsab*maxblock,0,nthread-1)
    CALL reallocate(sdabt,1,maxdco,1,maxco,1,4,0,nthread-1)
    CALL reallocate(workt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(pVpabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(nabt,1,maxco,1,maxco,0,nthread-1)

    IF(qs_env%dft_control%qs_control%gapw .OR. qs_env%dft_control%qs_control%gapw_xc) THEN
       rho0_present = ASSOCIATED(sac_3c)
       IF(rho0_present) THEN
          CALL allocate_vtriple(VH_3c_list,nkind,natom,error=error)
       END IF
    ENDIF

    !   *********************************************************************
    !   ***  Load GTH pseudo potential data (local part -> PPL)  ************
    !   *********************************************************************

    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN
       ALLOCATE (vpplt(nkind,0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,gth_potential=gth_potential)

          DO j=0,nthread-1
             NULLIFY (vpplt(kkind,j)%cexp_ppl)
             IF (ASSOCIATED(gth_potential)) THEN
                CALL get_potential(potential=gth_potential,&
                     alpha_ppl=vpplt(kkind,j)%alpha_ppl,&
                     cexp_ppl=vpplt(kkind,j)%cexp_ppl,&
                     ppl_radius=vpplt(kkind,j)%ppl_radius)
             END IF
             NULLIFY (vpplt(kkind,j)%r2)
             NULLIFY (vpplt(kkind,j)%neighbor)
             NULLIFY (vpplt(kkind,j)%r)
          END DO
       END DO
    END IF

    !   *********************************************************************
    !   ***  Load all potential data                             ************
    !   *********************************************************************

    IF (all_potential_present) THEN
       ALLOCATE(vallt(nkind,0:nthread-1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO kkind=1,nkind
          atomic_kind => atomic_kind_set(kkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               all_potential=all_potential)
          DO j=0,nthread-1
             NULLIFY(vallt(kkind,j)%neighbor)
             NULLIFY(vallt(kkind,j)%force_c)
             IF (ASSOCIATED(all_potential)) THEN
                CALL get_potential(potential=all_potential,&
                     alpha_core_charge=vallt(kkind,j)%alpha_c,&
                     zeff=vallt(kkind,j)%zeta_c,&
                     ccore_charge=vallt(kkind,j)%core_charge,&
                     core_charge_radius=vallt(kkind,j)%core_radius)
             ELSE
                vallt(kkind,j)%alpha_c = 0.0_dp
                vallt(kkind,j)%zeta_c = 0.0_dp
                vallt(kkind,j)%core_charge = 0.0_dp
                vallt(kkind,j)%core_radius = 0.0_dp
             END IF
          END DO
       END DO
    END IF

    DO ikind = 1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,atom_list=atom_list)
       DO iat = 1,nat
          iatom = atom_list(iat)
          IF(rho0_present) THEN
             CALL retrieve_sac_list(VH_3c_list,sac_3c,iatom,ikind,nkind,error=error)
          ENDIF
       END DO
    END DO

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            first_sgf=first_sgfa,&
            lmax=la_max,&
            lmin=la_min,&
            nco_sum=ncoa_sum,&
            npgf=npgfa,&
            nset=nseta,&
            nsgf_set=nsgfa,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi=sphi_a,&
            zet=zeta)

       DO jkind=1,nkind
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)

          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
          IF(dkh_loc_atom)THEN
             IF(ikind/=jkind) CYCLE 
          END IF
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
               first_sgf=first_sgfb,&
               lmax=lb_max,&
               lmin=lb_min,&
               npgf=npgfb,&
               nset=nsetb,&
               nsgf_set=nsgfb,&
               pgf_radius=rpgfb,&
               set_radius=set_radius_b,&
               sphi=sphi_b,&
               zet=zetb)

          iab = ikind + nkind*(jkind - 1)

          IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
          CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)

          NULLIFY ( sab_orb_neighbor_list )
          ithread = 0
          ! !$      ithread = omp_get_thread_num()
          hab => habt(:,:,ithread)
          sab => sabt(:,:,ithread)
          sdab => sdabt(:,:,:,ithread)
          work => workt(:,:,ithread)
          pVpab => pVpabt(:,:,ithread)
          nab => nabt(:,:,ithread)

          IF (ppl_present) THEN
             vppl => vpplt(:,ithread)
          END IF
          IF (all_potential_present) THEN
             vall => vallt(:,ithread)
          END IF
          ALLOCATE (ai_work(ldai,ldai,MAX(1,ncoset(maxlppl)),ncoset(nder+1)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,atom=iatom,nnode=nnode)

             IF (ppl_present) THEN
                DO kkind=1,nkind
                   iac = ikind + nkind*(kkind - 1)
                   IF (.NOT.ASSOCIATED(sac_ppl(iac)%neighbor_list_set)) CYCLE
                   sac_ppl_neighbor_list => find_neighbor_list(sac_ppl(iac)%neighbor_list_set,atom=iatom)

                   CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,nnode=nneighbor)

                   ALLOCATE (vppl(kkind)%r2(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%neighbor(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%r(3,nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                   sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)

                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                           neighbor=vppl(kkind)%neighbor(kneighbor),&
                           r=vppl(kkind)%r(:,kneighbor))
                      vppl(kkind)%r2(kneighbor) =&
                           vppl(kkind)%r(1,kneighbor)*vppl(kkind)%r(1,kneighbor) +&
                           vppl(kkind)%r(2,kneighbor)*vppl(kkind)%r(2,kneighbor) +&
                           vppl(kkind)%r(3,kneighbor)*vppl(kkind)%r(3,kneighbor)
                      sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
                   END DO
                END DO
             END IF

             IF(all_potential_present) THEN
                DO kkind = 1,nkind
                   IF(ASSOCIATED(VH_3c_list(kkind,iatom)%neighbor)) &
                        vall(kkind)%neighbor =>  VH_3c_list(kkind,iatom)%neighbor
                END DO
             END IF

             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

             DO inode=1,nnode
                CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)
                IF(rho0_present) THEN
                   CALL get_neighbor_node(sab_orb_neighbor_node,reduced_3c_rho0=reduced_3c_rho0)
                END IF

                IF(dkh_loc_atom)THEN
                   IF(iatom /= jatom) THEN
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      IF (jatom /= last_jatom) THEN
                         new_atom_b = .TRUE.
                         last_jatom = jatom
                      ELSE
                         new_atom_b = .FALSE.
                      END IF
                      CYCLE
                   END IF
                   rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                END IF

                IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                END IF

                ! *** Use the symmetry of the first derivatives ***
                IF (iatom == jatom) THEN
                   f0 = 1.0_dp
                ELSE
                   f0 = 2.0_dp
                END IF

                IF (new_atom_b) THEN
                   IF (iatom <= jatom) THEN
                      irow = iatom
                      icol = jatom
                   ELSE
                      irow = jatom
                      icol = iatom
                   END IF

                   NULLIFY (s_block)
                   CALL dbcsr_add_block_node(matrix=matrix_spgf(1)%matrix,&
                        block_row=irow,&                 
                        block_col=icol,&                 
                        block=s_block,&                  
                        error=error)                     

                   NULLIFY(t_block)
                   CALL dbcsr_add_block_node(matrix=matrix_tpgf(1)%matrix,&
                        block_row=irow,&
                        block_col=icol,&
                        block=t_block,&
                        error=error)

                   NULLIFY (n_block)
                   CALL dbcsr_add_block_node(matrix=matrix_vnuc(1)%matrix,&
                        block_row=irow,&
                        block_col=icol,&
                        block=n_block,&
                        error=error)

                   NULLIFY (hp_block)
                   CALL dbcsr_add_block_node(matrix=matrix_hpgf(1)%matrix,&
                        block_row=irow,&
                        block_col=icol,&
                        block=hp_block,&
                        error=error)

                   NULLIFY (pVp_block)
                   CALL dbcsr_add_block_node(matrix=matrix_pVppgf(1)%matrix,&
                        block_row=irow,&
                        block_col=icol,&
                        block=pVp_block,&
                        error=error)
                END IF

                IF(dkh_loc_atom)THEN
                   n = SIZE(s_block,1)
                   ALLOCATE(loc_s_block(n,n))
                   ALLOCATE(loc_t_block(n,n))
                   ALLOCATE(loc_n_block(n,n))
                   ALLOCATE(loc_hp_block(n,n))
                   ALLOCATE(loc_pVp_block(n,n))
                   loc_s_block(:,:)=0.0_dp
                   loc_t_block(:,:)=0.0_dp
                   loc_n_block(:,:)=0.0_dp
                   loc_hp_block(:,:)=0.0_dp
                   loc_pVp_block(:,:)=0.0_dp
                END IF

                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)
                nrow = 1

                DO iset=1,nseta
                   pgfa=1
                   DO i=2,iset
                      DO k = la_min(i-1),la_max(i-1)
                         pgfa=pgfa+npgfa(i-1)*(k+1)*(k+2)/2
                      END DO
                   END DO

                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   tot_npgfa = 0
                   DO k = la_min(iset),la_max(iset)
                      tot_npgfa = (k+1)*(k+2)/2 + tot_npgfa
                   END DO
                   DO jset=1,nsetb
                      pgfb=1
                      DO i=2,jset
                         DO k = lb_min(i-1),lb_max(i-1)
                            pgfb=pgfb+npgfb(i-1)*(k+1)*(k+2)/2
                         END DO
                      END DO

                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      tot_npgfb = 0
                      DO k = lb_min(jset),lb_max(jset)
                         tot_npgfb = (k+1)*(k+2)/2 + tot_npgfb 
                      END DO

                      IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN
                         ! *** Calculate the primitive overlap integrals ***
                         CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                              rpgfa(:,iset),zeta(:,iset),&
                              lb_max(jset),lb_min(jset),npgfb(jset),&
                              rpgfb(:,jset),zetb(:,jset),&
                              rab,dab,sab,nder,return_s_derivatives,&
                              ai_work,ldai,sdab)
                         IF(dkh_loc_atom.and.rab2 == 0.0_dp) THEN
                            CALL store_integral(iatom,jatom,sab,loc_s_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         ELSE
                            CALL store_integral(iatom,jatom,sab,s_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         END IF
                         ! *** Calculate the primitive kinetic energy integrals ***
                         CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),&
                              rpgfa(:,iset),zeta(:,iset),&
                              lb_max(jset),lb_min(jset),npgfb(jset),&
                              rpgfb(:,jset),zetb(:,jset),&
                              rab,dab,sdab,hab,nder,return_t_derivatives,&
                              ai_work,ldai)
                         IF(dkh_loc_atom.and.rab2 == 0.0_dp) THEN
                            CALL store_integral(iatom,jatom,hab,loc_t_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         ELSE
                            CALL store_integral(iatom,jatom,hab,t_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         END IF
                         ! *** Loop over the GTH pseudo potential atoms (local part) ***
                         IF (ppl_present) THEN
                            hab=0.0_dp
                            pVpab=0.0_dp
                            DO kkind=1,nkind
                               IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE
                               DO kneighbor=1,SIZE(vppl(kkind)%neighbor)
                                  rac(:) = vppl(kkind)%r(:,kneighbor)
                                  rac2 = vppl(kkind)%r2(kneighbor)
                                  dac = SQRT(rac2)

                                  IF (set_radius_a(iset) + vppl(kkind)%ppl_radius < dac) CYCLE
                                  rbc(:) = rac(:) - rab(:)
                                  rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                                  dbc = SQRT(rbc2)
                                  IF (set_radius_b(jset) + vppl(kkind)%ppl_radius < dbc) CYCLE

                                  ! *** Calculate the GTH pseudo potential forces ***
                                  CALL overlap_ppl(&
                                       la_max(iset),la_min(iset),npgfa(iset),&
                                       rpgfa(:,iset),zeta(:,iset),&
                                       lb_max(jset),lb_min(jset),npgfb(jset),&
                                       rpgfb(:,jset),zetb(:,jset),&
                                       vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                                       vppl(kkind)%ppl_radius,&
                                       rab,dab,rac,dac,rbc,dbc,&
                                       hab,0,0,.FALSE.,ai_work,pVpab=pVpab,dkh_erfc=dkh_erfc)
                               END DO
                            END DO
                            IF(dkh_loc_atom.and.rab2 == 0.0_dp) THEN
                               CALL store_integral(iatom,jatom,hab,hp_block,pgfa,pgfb,npgfa(iset),&
                                    tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                    ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                    ncoset(lb_min(jset)-1)) 
                               CALL store_integral(iatom,jatom,pVpab,pVp_block,pgfa,pgfb,npgfa(iset),&
                                    tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                    ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                    ncoset(lb_min(jset)-1))
                            ELSE
                               CALL store_integral(iatom,jatom,hab,hp_block,pgfa,pgfb,npgfa(iset),&
                                    tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                    ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                    ncoset(lb_min(jset)-1)) 
                               CALL store_integral(iatom,jatom,pVpab,pVp_block,pgfa,pgfb,npgfa(iset),&
                                    tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                    ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                    ncoset(lb_min(jset)-1))
                            END IF
                         END IF
                      ELSE
                         DO jco=1,ncob
                            DO ico=1,ncoa
                               hab(ico,jco) = 0.0_dp
                            END DO
                         END DO
                      END IF

                      IF (all_potential_present ) THEN
                         CALL all_integrals(la_max(iset),la_min(iset),npgfa(iset),&
                              rpgfa(:,iset),zeta(:,iset),&
                              lb_max(jset),lb_min(jset),npgfb(jset),&
                              rpgfb(:,jset),zetb(:,jset),&
                              rab,dab,jkind,nkind,sbc_3c,&
                              reduced_3c_rho0,vall,&
                              hab,nder,nab=nab,pVpab=pVpab,dkh_erfc=dkh_erfc,error=error)
                         IF(dkh_loc_atom.and.rab2 == 0.0_dp) THEN
                            CALL store_integral(iatom,jatom,nab,loc_n_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                            CALL store_integral(iatom,jatom,hab,loc_hp_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1)) 
                            CALL store_integral(iatom,jatom,pVpab,loc_pVp_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         ELSE
                            CALL store_integral(iatom,jatom,nab,n_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                            CALL store_integral(iatom,jatom,hab,hp_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1)) 
                            CALL store_integral(iatom,jatom,pVpab,pVp_block,pgfa,pgfb,npgfa(iset),&
                                 tot_npgfa,npgfb(jset),tot_npgfb, ncoset(la_max(iset)),&
                                 ncoset(la_min(iset)-1),ncoset(lb_max(jset)),&
                                 ncoset(lb_min(jset)-1))
                         END IF
                         DO kkind = 1,nkind
                            IF (.NOT. ASSOCIATED(vall(kkind)%force_c)) CYCLE
                            DEALLOCATE(vall(kkind)%force_c,STAT=stat)
                         END DO
                      END IF
                   END DO
                   nrow = nrow + ncoa
                END DO
                IF(dkh_loc_atom.and.rab2 == 0.0_dp) THEN
                   IF(dkh_erfc)THEN
                      loc_hp_block(:,:) = loc_hp_block(:,:) + loc_n_block (:,:) 
                   END IF

                   CALL dkh_atom_transformation(loc_s_block,loc_hp_block,loc_t_block,loc_pVp_block,n,&
                        qs_env%rel_control%rel_dkh_order) 
                   loc_hp_block(:,:) = loc_hp_block(:,:) + loc_t_block (:,:)
                   IF(.not.dkh_erfc)THEN
                      loc_hp_block(:,:) = loc_hp_block(:,:) + loc_n_block (:,:) 
                   END IF
                   hp_block(:,:) = hp_block(:,:) + loc_hp_block(:,:)
                   t_block(:,:) = t_block(:,:) + loc_t_block(:,:)
                   s_block(:,:) = s_block(:,:) + loc_s_block(:,:)
                ENDIF

                IF(dkh_loc_atom) THEN
                   DEALLOCATE(loc_s_block)
                   DEALLOCATE(loc_t_block)
                   DEALLOCATE(loc_n_block)
                   DEALLOCATE(loc_pVp_block)
                   DEALLOCATE(loc_hp_block)
                END IF

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             END DO
             IF (ppl_present) THEN
                DO kkind=1,nkind
                   IF (ASSOCIATED(vppl(kkind)%r2)) THEN
                      DEALLOCATE (vppl(kkind)%r2,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      DEALLOCATE (vppl(kkind)%neighbor,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      DEALLOCATE (vppl(kkind)%r,STAT=stat)
                      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   END IF
                END DO
             END IF

             IF (all_potential_present) THEN
                DO kkind=1,nkind
                   IF (ASSOCIATED(vall(kkind)%neighbor)) THEN
                      NULLIFY (vall(kkind)%neighbor)
                      NULLIFY (vall(kkind)%force_c)
                   END IF
                END DO
             END IF



          END DO


          IF (ALLOCATED(ai_work)) THEN
             DEALLOCATE (ai_work,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END DO

    END DO


    CALL dbcsr_finalize(matrix_spgf(1)%matrix, error=error)
    CALL dbcsr_finalize(matrix_tpgf(1)%matrix, error=error)
    CALL dbcsr_finalize(matrix_vnuc(1)%matrix, error=error)
    CALL dbcsr_finalize(matrix_hpgf(1)%matrix, error=error)
    CALL dbcsr_finalize(matrix_pVppgf(1)%matrix, error=error)


    IF(dkh_full) THEN
       IF(qs_env%rel_control%rel_potential.eq.rel_pot_erfc)THEN
          CALL dbcsr_add(matrix_hpgf(1)%matrix,matrix_vnuc(1)%matrix,&
                         alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
       END IF

       CALL cp_fm_create(matrix_t_full, matrix_full, error=error)
       CALL cp_fm_create(matrix_s_full, matrix_full, error=error)
       CALL cp_fm_create(matrix_v_full, matrix_full, error=error)
       CALL cp_fm_create(matrix_pVp_full, matrix_full, error=error)

    !super dirty>
    CALL allocate_matrix_set(matrix_spgf_sm,1,error=error)
    CALL allocate_matrix_set(matrix_hpgf_sm,1,error=error)
    CALL allocate_matrix_set(matrix_tpgf_sm,1,error=error)
    CALL allocate_matrix_set(matrix_pVppgf_sm,1,error=error)
    CALL sm_from_dbcsr ( matrix_spgf_sm(1)%matrix, matrix_spgf(1)%matrix, distribution_2d, error=error )
    CALL sm_from_dbcsr ( matrix_hpgf_sm(1)%matrix, matrix_hpgf(1)%matrix, distribution_2d, error=error )
    CALL sm_from_dbcsr ( matrix_tpgf_sm(1)%matrix, matrix_tpgf(1)%matrix, distribution_2d, error=error )
    CALL sm_from_dbcsr ( matrix_pVppgf_sm(1)%matrix, matrix_pVppgf(1)%matrix, distribution_2d, error=error )
    !super dirty<

       CALL copy_sm_to_fm(matrix_spgf_sm(1)%matrix,matrix_s_full,error=error)
       CALL copy_sm_to_fm(matrix_hpgf_sm(1)%matrix,matrix_v_full,error=error)
       CALL copy_sm_to_fm(matrix_tpgf_sm(1)%matrix,matrix_t_full,error=error)
       CALL copy_sm_to_fm(matrix_pVppgf_sm(1)%matrix,matrix_pVp_full,error=error)

       CALL dkh_full_transformation(qs_env,matrix_s_full,matrix_v_full,matrix_t_full,&
            matrix_pVp_full, matrix_s_full%matrix_struct%nrow_global,&
            qs_env%rel_control%rel_dkh_order,error)


       CALL copy_fm_to_sm(matrix_s_full,matrix_spgf_sm(1)%matrix)
       CALL copy_fm_to_sm(matrix_v_full,matrix_hpgf_sm(1)%matrix)
       CALL copy_fm_to_sm(matrix_t_full,matrix_tpgf_sm(1)%matrix)

    !super dirty>
    CALL dbcsr_deallocate_matrix_set(matrix_tpgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_hpgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_spgf,error)
    CALL dbcsr_allocate_matrix_set(matrix_tpgf,1,error=error)
    CALL dbcsr_allocate_matrix_set(matrix_hpgf,1,error=error)
    CALL dbcsr_allocate_matrix_set(matrix_spgf,1,error=error)
    ALLOCATE(matrix_spgf(1)%matrix,matrix_hpgf(1)%matrix,matrix_tpgf(1)%matrix)
     CALL dbcsr_from_sm(matrix_spgf(1)%matrix,matrix_spgf_sm(1)%matrix,error)
     CALL dbcsr_from_sm(matrix_hpgf(1)%matrix,matrix_hpgf_sm(1)%matrix,error)
     CALL dbcsr_from_sm(matrix_tpgf(1)%matrix,matrix_tpgf_sm(1)%matrix,error)
     CALL deallocate_matrix_set(matrix_spgf_sm  ,error)
     CALL deallocate_matrix_set(matrix_hpgf_sm  ,error)
     CALL deallocate_matrix_set(matrix_tpgf_sm  ,error)
     CALL deallocate_matrix_set(matrix_pVppgf_sm,error)
    !super dirty<

       CALL dbcsr_add(matrix_hpgf(1)%matrix,matrix_tpgf(1)%matrix,&
                      alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)

       IF(qs_env%rel_control%rel_potential.eq.rel_pot_full)THEN
          CALL dbcsr_add(matrix_hpgf(1)%matrix,matrix_vnuc(1)%matrix,&
                         alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
       END IF

       CALL cp_fm_release(matrix_t_full, error=error)
       CALL cp_fm_release(matrix_s_full, error=error)
       CALL cp_fm_release(matrix_v_full, error=error)
       CALL cp_fm_release(matrix_pVp_full, error=error)
    END IF

    CALL cp_fm_struct_release(matrix_full, error=error)

    ! Transfer the sparse matrices with primitive gaussian to the sparse matrices
    ! with contracted gaussians   

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            first_sgf=first_sgfa,&
            lmax=la_max,&
            lmin=la_min,&
            nco_sum=ncoa_sum,&
            npgf=npgfa,&
            nset=nseta,&
            nsgf_set=nsgfa,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi = sphi_a,&
            zet=zeta)

       DO jkind=1,nkind
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
          IF(dkh_loc_atom)THEN
             IF(ikind/=jkind) CYCLE 
          END IF

          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
               first_sgf=first_sgfb,&
               lmax=lb_max,&
               lmin=lb_min,&
               npgf=npgfb,&
               nset=nsetb,&
               nsgf_set=nsgfb,&
               pgf_radius=rpgfb,&
               set_radius=set_radius_b,&
               sphi = sphi_b,&
               zet=zetb)

          iab = ikind + nkind*(jkind - 1)

          IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
          CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)

          NULLIFY ( sab_orb_neighbor_list )

          DO ilist=1,nlist  
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
                sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,atom=iatom,nnode=nnode)

             ! *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***

             IF (ppl_present) THEN
                DO kkind=1,nkind
                   iac = ikind + nkind*(kkind - 1)
                   IF (.NOT.ASSOCIATED(sac_ppl(iac)%neighbor_list_set)) CYCLE
                   sac_ppl_neighbor_list => find_neighbor_list(sac_ppl(iac)%neighbor_list_set,atom=iatom)

                   CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,nnode=nneighbor)

                   ALLOCATE (vppl(kkind)%r2(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%neighbor(nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE (vppl(kkind)%r(3,nneighbor),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                   sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)

                   DO kneighbor=1,nneighbor
                      CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                           neighbor=vppl(kkind)%neighbor(kneighbor),&
                           r=vppl(kkind)%r(:,kneighbor))
                      vppl(kkind)%r2(kneighbor) =&
                           vppl(kkind)%r(1,kneighbor)*vppl(kkind)%r(1,kneighbor) +&
                           vppl(kkind)%r(2,kneighbor)*vppl(kkind)%r(2,kneighbor) +&
                           vppl(kkind)%r(3,kneighbor)*vppl(kkind)%r(3,kneighbor)
                      sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
                   END DO
                END DO
             END IF
             IF(all_potential_present) THEN
                DO kkind = 1,nkind
                   IF(ASSOCIATED(VH_3c_list(kkind,iatom)%neighbor)) &
                        vall(kkind)%neighbor =>  VH_3c_list(kkind,iatom)%neighbor
                END DO
             END IF

             last_jatom = 0
             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

             DO inode=1,nnode
                CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)
                atomic_kind => atomic_kind_set(ikind) !localization?

                IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                ELSE
                   new_atom_b = .FALSE.
                END IF

                IF(rho0_present) THEN
                   CALL get_neighbor_node(sab_orb_neighbor_node,reduced_3c_rho0=reduced_3c_rho0)
                END IF
                IF(oce_present) THEN
                   CALL get_neighbor_node(sab_orb_neighbor_node,reduced_3c_oce=reduced_3c_oce)
                END IF

                ! *** Use the symmetry of the first derivatives ***
                IF (iatom == jatom) THEN
                   f0 = 1.0_dp
                ELSE
                   f0 = 2.0_dp
                END IF

                ! *** Create matrix blocks for a new matrix block column ***
                IF (new_atom_b) THEN
                   IF (iatom <= jatom) THEN
                      iblock_row = iatom
                      iblock_col = jatom
                   ELSE
                      iblock_row = jatom
                      iblock_col = iatom
                   END IF

                   NULLIFY(spgf_block)
                   CALL dbcsr_get_block_node(matrix=matrix_spgf(1)%matrix,&
                        block_row=iblock_row,&
                        block_col=iblock_col,&
                        BLOCK=spgf_block)

                   NULLIFY(hpgf_block)
                   CALL dbcsr_get_block_node(matrix=matrix_hpgf(1)%matrix,&
                        block_row=iblock_row,&
                        block_col=iblock_col,&
                        BLOCK=hpgf_block)

                   NULLIFY(tpgf_block)
                   CALL dbcsr_get_block_node(matrix=matrix_tpgf(1)%matrix,&
                        block_row=iblock_row,&
                        block_col=iblock_col,&
                        BLOCK=tpgf_block)

                   NULLIFY(s_block)
                   CALL dbcsr_get_block_node(matrix=matrix_s(1)%matrix,&
                        block_row=iblock_row,&
                        block_col=iblock_col,&
                        BLOCK=s_block)

                   NULLIFY(hp_block)
                   CALL dbcsr_get_block_node(matrix=matrix_h(1)%matrix,&
                        block_row=iblock_row,&
                        block_col=iblock_col,&
                        BLOCK=hp_block)

                   IF(kinetic_m)THEN
                      NULLIFY(t_block)
                      CALL dbcsr_get_block_node(matrix=matrix_t(1)%matrix,&
                           block_row=iblock_row,&
                           block_col=iblock_col,&
                           BLOCK=t_block)
                   END IF
                END IF

                IF(dkh_loc_atom)THEN
                   IF(iatom/=jatom)THEN
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      CYCLE
                   END IF
                END IF

                DO iset=1,nseta
                   pgfa = 1
                   DO i=2,iset
                      DO k = la_min(i-1),la_max(i-1)
                         pgfa=pgfa+npgfa(i-1)*(k+1)*(k+2)/2
                      END DO
                   END DO
                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   tot_npgfa = 0
                   DO k = la_min(iset),la_max(iset)
                      tot_npgfa = (k+1)*(k+2)/2 + tot_npgfa
                   END DO

                   DO jset=1,nsetb
                      pgfb = 1
                      DO i=2,jset
                         DO k = lb_min(i-1),lb_max(i-1)
                            pgfb=pgfb+npgfb(i-1)*(k+1)*(k+2)/2
                         END DO
                      END DO

                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      tot_npgfb = 0
                      DO k = lb_min(jset),lb_max(jset)
                         tot_npgfb = (k+1)*(k+2)/2 + tot_npgfb 
                      END DO

                      sab = 0.0_dp
                      CALL present_integral(iatom,jatom,sab,spgf_block,pgfa,pgfb,&
                           npgfa(iset),tot_npgfa,npgfb(jset),tot_npgfb,&
                           ncoset(la_max(iset)),ncoset(la_min(iset)-1),&
                           ncoset(lb_max(jset)),ncoset(lb_min(jset)-1)) 

                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,sab(1,1),SIZE(sab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                         CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              work(1,1),SIZE(work,1),&
                              0.0_dp,s_block(sgfa,sgfb),&
                              SIZE(s_block,1))
                      ELSE
                         CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              0.0_dp,s_block(sgfb,sgfa),&
                              SIZE(s_block,1))
                      END IF

                      IF(kinetic_m)THEN   
                         tab = 0.0_dp
                         CALL present_integral(iatom,jatom,tab,tpgf_block,pgfa,pgfb,&
                              npgfa(iset),tot_npgfa,npgfb(jset),tot_npgfb,&
                              ncoset(la_max(iset)),ncoset(la_min(iset)-1),&
                              ncoset(lb_max(jset)),ncoset(lb_min(jset)-1)) 

                         CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                              1.0_dp,hab(1,1),SIZE(hab,1),&
                              sphi_b(1,sgfb),SIZE(sphi_b,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                         IF (iatom <= jatom) THEN
                            CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 0.0_dp,t_block(sgfa,sgfb),&
                                 SIZE(t_block,1))
                         ELSE
                            CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 1.0_dp,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 0.0_dp,t_block(sgfb,sgfa),&
                                 SIZE(t_block,1))
                         END IF
                      END IF

                      hab = 0.0_dp
                      CALL present_integral(iatom,jatom,hab,hpgf_block,pgfa,pgfb,&
                           npgfa(iset),tot_npgfa,npgfb(jset),tot_npgfb,&
                           ncoset(la_max(iset)),ncoset(la_min(iset)-1),&
                           ncoset(lb_max(jset)),ncoset(lb_min(jset)-1)) 

                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                         CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              work(1,1),SIZE(work,1),&
                              0.0_dp,hp_block(sgfa,sgfb),SIZE(hp_block,1))
                      ELSE
                         CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              0.0_dp,hp_block(sgfb,sgfa),SIZE(hp_block,1))
                      END IF
                   END DO
                END DO
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
             END DO
          END DO
          IF (ppl_present) THEN
             DO kkind=1,nkind
                IF (ASSOCIATED(vppl(kkind)%r2)) THEN
                   DEALLOCATE (vppl(kkind)%r2,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   DEALLOCATE (vppl(kkind)%neighbor,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   DEALLOCATE (vppl(kkind)%r,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
             END DO
          END IF

          IF (all_potential_present) THEN
             DO kkind=1,nkind
                IF (ASSOCIATED(vall(kkind)%neighbor)) THEN
                   NULLIFY (vall(kkind)%neighbor)
                   NULLIFY (vall(kkind)%force_c)
                END IF
             END DO
          END IF
       END DO
    END DO
    DO kkind=1,nkind
       IF (ALLOCATED(ai_work)) THEN
          DEALLOCATE (ai_work,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    CALL sm_from_dbcsr ( matrix_s_sm(1)%matrix, matrix_s(1)%matrix, distribution_2d, error=error )
    CALL sm_from_dbcsr ( matrix_h_sm(1)%matrix, matrix_h(1)%matrix, distribution_2d, error=error )
    IF(kinetic_m) CALL sm_from_dbcsr ( matrix_t_sm(1)%matrix, matrix_t(1)%matrix, distribution_2d, error=error )

    CALL set_qs_env(qs_env=qs_env,&
         matrix_s=matrix_s_sm,&
         matrix_h=matrix_h_sm,&
         error=error)

    IF(kinetic_m)THEN
       CALL set_qs_env(qs_env=qs_env,&
            kinetic=matrix_t_sm,&
            error=error)
    END IF

    ! *** Release work storage ***
    IF (ppl_present) THEN
       DO kkind=1,nkind
          DO j=0,nthread-1
             NULLIFY (vpplt(kkind,j)%cexp_ppl)
          END DO
       END DO
       DEALLOCATE (vpplt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (all_potential_present) THEN
       DO kkind=1,nkind
          DO j =0,nthread-1
             NULLIFY (vallt(kkind,j)%neighbor)
             NULLIFY (vallt(kkind,j)%force_c)
          END DO
       END DO
       DEALLOCATE (vallt,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF(ASSOCIATED(VH_3c_list))THEN
       CALL deallocate_vtriple(VH_3c_list,error=error)
    END IF


    DEALLOCATE (habt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sdabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (workt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (pabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (pVpabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (nabt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! deallocate temporary DKH matrices
    CALL dbcsr_deallocate_matrix_set(matrix_tpgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_hpgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_spgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_vnuc,error)
    CALL dbcsr_deallocate_matrix_set(matrix_pVppgf,error)
    CALL dbcsr_deallocate_matrix_set(matrix_s,error)
    CALL dbcsr_deallocate_matrix_set(matrix_h,error)
    CALL dbcsr_deallocate_matrix_set(matrix_t,error)

    CALL timestop(handle)

  END SUBROUTINE dkh_mol_integrals

! *****************************************************************************
  SUBROUTINE store_integral2(iatom,jatom,small_block,block,pgfa,pgfb,npgfa,tot_npgfa,npgfb,tot_npgfb,&
       ncoseta,ncosetam,ncosetb,ncosetbm)

    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: small_block
    REAL(KIND=dp), DIMENSION(:, :)           :: block
    INTEGER, INTENT(IN)                      :: pgfa, pgfb, npgfa, tot_npgfa, &
                                                npgfb, tot_npgfb, ncoseta, &
                                                ncosetam, ncosetb, ncosetbm

    INTEGER                                  :: i, i2, j, j2

    IF (iatom <= jatom) THEN
       DO i   = 1,npgfa
          DO i2   = 1,tot_npgfa
             DO j   = 1,npgfb
                DO j2 = 1,tot_npgfb
                   block(pgfa-1+(i-1)*tot_npgfa+i2,pgfb-1+(j-1)*tot_npgfb+j2)=&
                        small_block((i-1)*ncoseta+ncosetam+i2,&
                        (j-1)*ncosetb+ncosetbm+j2)+block(pgfa-1+(i-1)*tot_npgfa+i2,pgfb-1+(j-1)*tot_npgfb+j2)
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO i = 1,npgfa
          DO i2 = 1,tot_npgfa
             DO j = 1,npgfb
                DO j2 = 1,tot_npgfb
                   block(pgfb-1+(j-1)*tot_npgfb+j2,pgfa-1+(i-1)*tot_npgfa+i2)=&
                        small_block((i-1)*ncoseta+ncosetam+i2,&
                        (j-1)*ncosetb+ncosetbm+j2)+block(pgfb-1+(j-1)*tot_npgfb+j2,pgfa-1+(i-1)*tot_npgfa+i2)
                END DO
             END DO
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE store_integral2

! *****************************************************************************
  SUBROUTINE store_integral(iatom,jatom,small_block,block,pgfa,pgfb,npgfa,tot_npgfa,npgfb,tot_npgfb,&
       ncoseta,ncosetam,ncosetb,ncosetbm)

    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: small_block, block
    INTEGER, INTENT(IN)                      :: pgfa, pgfb, npgfa, tot_npgfa, &
                                                npgfb, tot_npgfb, ncoseta, &
                                                ncosetam, ncosetb, ncosetbm

    INTEGER                                  :: i, i2, j, j2

    IF (iatom <= jatom) THEN
       DO i   = 1,npgfa
          DO i2   = 1,tot_npgfa
             DO j   = 1,npgfb
                DO j2 = 1,tot_npgfb
                   block(pgfa-1+(i-1)*tot_npgfa+i2,pgfb-1+(j-1)*tot_npgfb+j2)=&
                        small_block((i-1)*ncoseta+ncosetam+i2,&
                        (j-1)*ncosetb+ncosetbm+j2)+block(pgfa-1+(i-1)*tot_npgfa+i2,pgfb-1+(j-1)*tot_npgfb+j2)
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO i = 1,npgfa
          DO i2 = 1,tot_npgfa
             DO j = 1,npgfb
                DO j2 = 1,tot_npgfb
                   block(pgfb-1+(j-1)*tot_npgfb+j2,pgfa-1+(i-1)*tot_npgfa+i2)=&
                        small_block((i-1)*ncoseta+ncosetam+i2,&
                        (j-1)*ncosetb+ncosetbm+j2)+block(pgfb-1+(j-1)*tot_npgfb+j2,pgfa-1+(i-1)*tot_npgfa+i2)
                END DO
             END DO
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE store_integral

! *****************************************************************************
  SUBROUTINE present_integral(iatom,jatom,small_block,block,pgfa,pgfb,npgfa,tot_npgfa,npgfb,tot_npgfb,&
       ncoseta,ncosetam,ncosetb,ncosetbm)

    INTEGER, INTENT(IN)                      :: iatom, jatom
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: small_block, block
    INTEGER, INTENT(IN)                      :: pgfa, pgfb, npgfa, tot_npgfa, &
                                                npgfb, tot_npgfb, ncoseta, &
                                                ncosetam, ncosetb, ncosetbm

    INTEGER                                  :: i, i2, j, j2

    IF (iatom <= jatom) THEN
       DO i = 1,npgfa
          DO i2 = 1,tot_npgfa
             DO j = 1,npgfb
                DO j2 = 1,tot_npgfb
                   small_block((i-1)*ncoseta+ncosetam+i2, (j-1)*ncosetb+ncosetbm+j2)=&
                        block(pgfa-1+(i-1)*tot_npgfa+i2,pgfb-1+(j-1)*tot_npgfb+j2)
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO i = 1,npgfa
          DO i2 = 1,tot_npgfa
             DO j = 1,npgfb
                DO j2 = 1,tot_npgfb
                   small_block((i-1)*ncoseta+ncosetam+i2, (j-1)*ncosetb+ncosetbm+j2)=&
                        block(pgfb-1+(j-1)*tot_npgfb+j2,pgfa-1+(i-1)*tot_npgfa+i2)
                END DO
             END DO
          END DO
       END DO
    END IF

    RETURN

  END SUBROUTINE present_integral

! *****************************************************************************
!> \brief 2th order DKH calculations                                 
!>                                                                     
!> \author
!>     Jens Thar, Barbara Kirchner: Uni Bonn (09/2006)                 
!>     Markus Reiher: ETH Zurich (09/2006)                             
!>                                                                     
!> \par Literature                                              
!>  M. Reiher, A. Wolf, J. Chem. Phys. 121 (2004) 10944-10956          
!>  A. Wolf, M. Reiher, B. A. Hess, J. Chem. Phys. 117 (2002) 9215-9226
!>                                                                     
!>\par Note
!>      based on subroutines for DKH1 to DKH5 by                       
!>       A. Wolf, M. Reiher, B. A. Hess                                
!>                                                                      
!>  INPUT:                                                              
!>    qs_env (:)        The quickstep environment                       
!>    n                 Number of primitive gaussians                   
!>    matrix_s    (:,:) overlap matrix                                  
!>    matrix_pVp  (:,:) pVp matrix                                      
!>                                                                      
!>  IN_OUT:                                                             
!>    matrix_v    (:,:) input: nonrelativistic potential energy matrix  
!>                      output: (ev1+ev2)                               
!>    matrix_t    (:,:) input: kinetic energy matrix                    
!>                      output: kinetic part of hamiltonian              
!>                      in position space                               
!>                                                                      
!>  INTERNAL                                                            
!>    sinv (:,:) inverted, orthogonalized overlap matrix                
!>    ev0t (:)   DKH-even0 matrix in T-basis                            
!>    e    (:)   e=SQRT(p^2c^2+c^4)                                     
!>    eig  (:,:) eigenvectors of sinv' h sinv                           
!>    tt   (:)   eigenvalues of sinv' h sinv                            
!>    revt (:,:) reverse transformation matrix T-basis -> position space
!>    aa   (:)   kinematical factors f. DKH SQRT((c^2+e(i))/(2.0*e(i))) 
!>    rr   (:)   kinematical factors f. DKH c/(c^2+e(i))                
!>    vt   (:,:) non relativistic potential matrix in T-basis           
!>    pvpt (:,:) pvp integral matrix in T-basis                          
!>    ev1t (:,:) DKH-even1 matrix in T-basis                            
!>    evt2 (:,:) DKH-even2 matrix in T-basis                            
!>    ev1  (:,:) DKH-even1 matrix in position space                     
!>    ev2  (:,:) DKH-even2 matrix in position space                     
!>    ove (:,:) scratch                                                 
!>    aux (:,:) scratch                                                 
!>    velit  velocity of light 137 a.u.                                 
!>    prea   prefactor, 1/137^2                                         
!>    con2   prefactor, 2/137^2                                         
!>    con    prefactor, 137^2                                           
! *****************************************************************************
  SUBROUTINE DKH_full_transformation (qs_env,matrix_s,matrix_v,matrix_t,matrix_pVp,n,dkh_order,error) 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: matrix_s, matrix_v, matrix_t, &
                                                matrix_pVp
    INTEGER, INTENT(IN)                      :: n, dkh_order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DKH_full_transformation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    REAL(KIND=dp)                            :: velit
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aa, e, ev0t, rr, tt
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER :: matrix_aux, matrix_aux2, matrix_eig, &
      matrix_ev1, matrix_ev2, matrix_ev3, matrix_ev4, matrix_pe1p, &
      matrix_rev, matrix_se, matrix_sinv

    CALL timeset(routineN,handle)
    !-----------------------------------------------------------------------
    !     Define velocity of light
    !-----------------------------------------------------------------------

    !     velit = 137.0359895_dp   ! 1/a_fine  physcon
    velit = 1._dp/a_fine 

    !-----------------------------------------------------------------------
    !     Construct the matrix structure
    !-----------------------------------------------------------------------

    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = qs_env%blacs_env,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    !-----------------------------------------------------------------------
    !     Allocate some matrices
    !-----------------------------------------------------------------------

    ALLOCATE(e(n))
    ALLOCATE(aa(n))
    ALLOCATE(rr(n))
    ALLOCATE(tt(n))
    ALLOCATE(ev0t(n))

    NULLIFY(matrix_eig)
    NULLIFY(matrix_aux)
    NULLIFY(matrix_aux2)
    NULLIFY(matrix_rev)
    NULLIFY(matrix_se)
    NULLIFY(matrix_ev1)
    NULLIFY(matrix_ev2)
    NULLIFY(matrix_ev3)
    NULLIFY(matrix_ev4)
    NULLIFY(matrix_sinv)
    NULLIFY(matrix_pe1p)

    CALL cp_fm_create(matrix_eig,  matrix_full, error=error)
    CALL cp_fm_create(matrix_aux,  matrix_full, error=error)
    CALL cp_fm_create(matrix_aux2,  matrix_full, error=error)
    CALL cp_fm_create(matrix_rev,  matrix_full, error=error)
    CALL cp_fm_create(matrix_se,   matrix_full, error=error)
    CALL cp_fm_create(matrix_ev1,  matrix_full, error=error)
    CALL cp_fm_create(matrix_ev2,  matrix_full, error=error)
    CALL cp_fm_create(matrix_sinv, matrix_full, error=error)
    CALL cp_fm_create(matrix_ev3,  matrix_full, error=error)
    CALL cp_fm_create(matrix_ev4,  matrix_full, error=error)
    CALL cp_fm_create(matrix_pe1p,  matrix_full, error=error)


    !-----------------------------------------------------------------------
    !     Now with Cholesky decomposition
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_s,matrix_sinv,error)
    CALL cp_fm_cholesky_decompose(matrix_sinv,n,error=error)

    !-----------------------------------------------------------------------
    !     Calculate matrix representation from nonrelativistic T matrix
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_t,matrix_sinv, error)  
    CALL cp_fm_syevd(matrix_t,matrix_eig,tt,error)

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of Hamiltonian in T-basis 
    !-----------------------------------------------------------------------

    CALL kintegral (n,ev0t,tt,e,velit)       

    !-----------------------------------------------------------------------
    !     Calculate reverse transformation matrix revt
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_eig,matrix_rev,error)
    CALL cp_fm_triangular_multiply(matrix_sinv,matrix_rev,transpose_tr=.TRUE.,error=error)

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of the Hamiltonian  
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_rev,matrix_aux,error)  
    CALL cp_fm_column_scale(matrix_aux,ev0t) 
    CALL cp_fm_gemm("N","T",n,n,n,1.0_dp,matrix_rev,matrix_aux,0.0_dp,matrix_t,error)

    !-----------------------------------------------------------------------
    !     Calculate kinematical factors for DKH
    !     only vectors present - will be done on every CPU
    !-----------------------------------------------------------------------

    DO i=1,n
       aa(i)=SQRT((velit*velit+e(i)) / (2.0_dp*e(i)))            
       rr(i)=SQRT(velit*velit)/(velit*velit+e(i))                        
    END DO

    !----------------------------------------------------------------------- 
    !     Transform v integrals to T-basis (v -> v(t))
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_v,matrix_sinv, error)  
    CALL cp_fm_upper_to_full(matrix_v,matrix_aux,error)
    CALL cp_fm_gemm("T","N",n,n,n,1.0_dp,matrix_eig,matrix_v,0.0_dp,matrix_aux,error) 
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_eig,0.0_dp,matrix_v,error)

    !-----------------------------------------------------------------------
    !     Transform pVp integrals to T-basis (pVp -> pVp(t))
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_pVp,matrix_sinv, error)  
    CALL cp_fm_upper_to_full(matrix_pVp,matrix_aux,error)
    CALL cp_fm_gemm("T","N",n,n,n,1.0_dp,matrix_eig,matrix_pVp,0.0_dp,matrix_aux,error)
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_eig,0.0_dp,matrix_pVp,error)

    !-----------------------------------------------------------------------
    !     Calculate even1 in T-basis
    !-----------------------------------------------------------------------

    CALL even1(matrix_ev1,matrix_v,matrix_pvp,aa,rr,matrix_aux,matrix_aux2,error)

    !-----------------------------------------------------------------------
    !     Calculate even2 in T-basis
    !-----------------------------------------------------------------------

    CALL even2c (n,matrix_ev2,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !-----------------------------------------------------------------------
    !     Calculate even3 in T-basis, only if requested
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.3) THEN       
       CALL peven1p(n,matrix_pe1p,matrix_v,matrix_pvp,matrix_aux,matrix_aux2,aa,rr,tt,error)
       CALL even3b(n,matrix_ev3,matrix_ev1,matrix_pe1p,matrix_v,matrix_pvp,aa,rr,tt,e,matrix_aux,error)

       !-----------------------------------------------------------------------
       !     Transform even3 back to position space
       !-----------------------------------------------------------------------

       CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev3,0.0_dp,matrix_aux,error)
       CALL cp_fm_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev3,error)

       !-----------------------------------------------------------------------
       !     Calculate even4 in T-basis, only if requested
       !-----------------------------------------------------------------------

       IF (dkh_order.ge.4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="DKH order greater than 3 not yet available", &
               error=error, error_level=cp_failure_level)
          !          CALL even4a(n,matrix_ev4%local_data,matrix_ev2%local_data,matrix_pe1p%local_data,matrix_v%local_data,&
          !                      matrix_pvp%local_data,aa,rr,tt,e)

          !-----------------------------------------------------------------------
          !     Transform even4 back to position space
          !-----------------------------------------------------------------------

          !        CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev4,0.0_dp,matrix_aux,error)
          !        CALL cp_fm_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev4,error)

       END IF
    END IF

    !----------------------------------------------------------------------
    !     Transform even1 back to position space
    !----------------------------------------------------------------------

    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev1,0.0_dp,matrix_aux,error)
    CALL cp_fm_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev1,error)

    !-----------------------------------------------------------------------
    !     Transform even2 back to position space
    !-----------------------------------------------------------------------

    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev2,0.0_dp,matrix_aux,error)
    CALL cp_fm_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev2,error)


    !-----------------------------------------------------------------------
    !     Calculate v in position space
    !-----------------------------------------------------------------------
    !
    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev1,1.0_dp,matrix_ev2,error)
    CALL cp_fm_upper_to_full(matrix_ev1,matrix_aux,error)    
    CALL cp_fm_to_fm(matrix_ev1,matrix_v,error)
    IF(dkh_order.ge.3) THEN
       CALL cp_fm_scale_and_add(1.0_dp,matrix_v,1.0_dp,matrix_ev3,error)
       IF(dkh_order.ge.4) THEN
          CALL cp_fm_scale_and_add(1.0_dp,matrix_v,1.0_dp,matrix_ev4,error)
       END IF
    END IF

    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_eig,  error=error)
    CALL cp_fm_release(matrix_aux,  error=error)
    CALL cp_fm_release(matrix_aux2, error=error)
    CALL cp_fm_release(matrix_rev,  error=error)
    CALL cp_fm_release(matrix_se,   error=error)
    CALL cp_fm_release(matrix_ev1,  error=error)
    CALL cp_fm_release(matrix_ev2,  error=error)
    CALL cp_fm_release(matrix_sinv, error=error)
    CALL cp_fm_release(matrix_ev3,  error=error)
    CALL cp_fm_release(matrix_ev4,  error=error)
    CALL cp_fm_release(matrix_pe1p, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)

    DEALLOCATE(ev0t,e,aa,rr,tt)

    CALL timestop(handle)

  END SUBROUTINE DKH_full_transformation

! *****************************************************************************
  SUBROUTINE kintegral (n,ev0t,tt,e,velit)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ev0t
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tt
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: e
    REAL(KIND=dp), INTENT(IN)                :: velit

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: con, con2, prea, ratio, tv1, &
                                                tv2, tv3, tv4

    prea=1/(velit*velit)
    con2=prea+prea
    con=1.0_dp/prea

    DO i=1,n
       IF (tt(i).LT.0.0_dp) THEN
          WRITE(*,*) ' dkh_main.F | tt(',i,') = ',tt(i)
       END IF

       !       If T is sufficiently small, use series expansion to avoid
       !       cancellation, otherwise calculate SQRT directly

       ev0t(i)=tt(i)
       ratio=tt(i)/velit
       IF (ratio.LE.0.02_dp) THEN
          tv1=tt(i)
          tv2=-tv1*tt(i)*prea*0.5_dp
          tv3=-tv2*tt(i)*prea
          tv4=-tv3*tt(i)*prea*1.25_dp
          ev0t(i)=tv1+tv2+tv3+tv4
       ELSE
          ev0t(i)=con*(SQRT(1.0_dp+con2*tt(i))-1.0_dp)   
       END IF
       e(i)=ev0t(i)+con           
    END DO

    RETURN
  END SUBROUTINE kintegral

! *****************************************************************************
  SUBROUTINE even1(matrix_ev1,matrix_v,matrix_pvp,aa,rr,matrix_aux,matrix_aux2,error)
    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     ev1  out  even1 output matrix                                    -
    !     v    in   potential matrix v in T-space                          -
    !     pvp  in   pvp matrix in T-space                                  -
    !     aa   in   A-factors (diagonal)                                   -
    !     rr   in   R-factors (diagonal)                                   -
    !                                                                      -
    !-----------------------------------------------------------------------
    TYPE(cp_fm_type), POINTER                :: matrix_ev1, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux, matrix_aux2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL cp_fm_to_fm(matrix_v,matrix_aux,error)
    CALL cp_fm_column_scale(matrix_aux,aa)
    CALL cp_fm_transpose(matrix_aux,matrix_ev1,error)
    CALL cp_fm_column_scale(matrix_ev1,aa)

    CALL cp_fm_to_fm(matrix_pVp,matrix_aux,error) 
    CALL cp_fm_column_scale(matrix_aux,aa)
    CALL cp_fm_column_scale(matrix_aux,rr)
    CALL cp_fm_transpose(matrix_aux,matrix_aux2,error)
    CALL cp_fm_column_scale(matrix_aux2,aa)
    CALL cp_fm_column_scale(matrix_aux2,rr)

    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev1,1.0_dp,matrix_aux2,error)

    RETURN
  END SUBROUTINE even1

! *****************************************************************************
  SUBROUTINE peven1p(n,matrix_pe1p,matrix_v,matrix_pvp,matrix_aux,matrix_aux2,aa,rr,tt,error)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n      in   dimension of matrices                                -
    !     pev1tp out  peven1p output matrix                                -
    !     vt     in   potential matrix v in T-space                        -
    !     pvpt   in   pvp matrix in T-space                                -
    !     aa     in   A-factors (diagonal)                                 -
    !     rr     in   R-factors (diagonal)                                 -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_pe1p, matrix_v, &
                                                matrix_pvp, matrix_aux, &
                                                matrix_aux2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    REAL(KIND=dp), DIMENSION(n)              :: vec_ar, vec_arrt
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    DO i=1,n
       vec_ar(i)=aa(i)*rr(i)
       vec_arrt(i)=vec_ar(i)*rr(i)*tt(i)
    END DO

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    CALL cp_fm_get_info(matrix_v, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = vec_arrt(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_aux2,error)
    CALL cp_fm_schur_product(matrix_v,matrix_aux,matrix_pe1p,error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = vec_ar(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_aux2,error)
    CALL cp_fm_schur_product(matrix_pvp,matrix_aux,matrix_aux2,error)

    CALL cp_fm_scale_and_add(4.0_dp,matrix_pe1p,1.0_dp,matrix_aux2,error)  

    CALL cp_fm_release(vec_a,error=error)
    CALL cp_fm_struct_release(vec_full,error=error)

    RETURN
  END SUBROUTINE peven1p

  !-----------------------------------------------------------------------

  SUBROUTINE even2c (n,matrix_ev2,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified: 20.02.2002 - DKH2                 *
    !                                                                      *
    !     2nd order DK-approximation ( original DK-transformation with     *
    !                                       U = SQRT(1+W^2) + W        )   *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (6.2.2002)                                         *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !     2008       Jens Thar: transfer to CP2K                           *
    !                                                                      *
    !     ev2 = 1/2 [W1,O1]                                                *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n           in   Dimension of matrices                           *
    !     ev2         out  even2 output matrix = final result              *
    !     vv          in   potential v                                     *
    !     gg          in   pvp                                             *
    !     aa          in   A-Factors (DIAGONAL)                            *
    !     rr          in   R-Factors (DIAGONAL)                            *
    !     tt          in   Nonrel. kinetic Energy (DIAGONAL)               *
    !     e           in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)   *
    !     matrix_v    symmetric (n x n)-matrix containing (A V A)          *
    !     matrix_pvp  symmetric (n x n)-matrix containing (A P V P A)      *
    !     vh          symmetric (n x n)-matrix containing (A V~ A)         *
    !     pvph        symmetric (n x n)-matrix containing (A P V~ P A)     *
    !     w1o1        W1*O1 (2-component form)                             *
    !     o1w1        O1*W1 (2-component form)                             *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_ev2, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER                :: matrix_apVpa, matrix_apVVpa, &
                                                matrix_aux2, matrix_ava, &
                                                matrix_avva

!     result  intermediate result of even2-calculation
!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH2                      
!-----------------------------------------------------------------------

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    NULLIFY(matrix_aux2)
    NULLIFY(matrix_ava)
    NULLIFY(matrix_avva)
    NULLIFY(matrix_apVpa)
    NULLIFY(matrix_apVVpa)

    CALL cp_fm_create(matrix_aux2,    matrix_full, error=error)
    CALL cp_fm_create(matrix_ava,     matrix_full, error=error)
    CALL cp_fm_create(matrix_avva,    matrix_full, error=error)
    CALL cp_fm_create(matrix_apVpa,   matrix_full, error=error)
    CALL cp_fm_create(matrix_apVVpa,  matrix_full, error=error)

    CALL cp_fm_to_fm(matrix_v,matrix_ava,error)
    CALL cp_fm_to_fm(matrix_v,matrix_avva,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVpa,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVVpa,error)

    !  Calculate  v = A V A:

    CALL mat_axa(matrix_v,matrix_ava,n,aa,matrix_aux,error)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra(matrix_pVp,matrix_apVpa,n,aa,rr,matrix_aux,error)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h(matrix_v,matrix_avva,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_avva,matrix_aux2,error)
    CALL mat_axa(matrix_aux2,matrix_avva,n,aa,matrix_aux,error)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h(matrix_pVp,matrix_apVVpa,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_apVVpa,matrix_aux2,error)
    CALL mat_arxra(matrix_aux2,matrix_apVVpa,n,aa,rr,matrix_aux,error)

    !  Calculate w1o1:

    CALL cp_fm_gemm("N","N",n,n,n,-1.0_dp,matrix_apVVpa,matrix_ava,0.0_dp,matrix_aux2,error)
    CALL mat_muld(matrix_aux2,matrix_apVVpa,matrix_apVpa,n,  1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL mat_mulm(matrix_aux2,matrix_avva,  matrix_ava,n,    1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL cp_fm_gemm("N","N",n,n,n,-1.0_dp,matrix_avva,matrix_apVpa,1.0_dp,matrix_aux2,error)

    !  Calculate o1w1 (already stored in ev2):

    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_apVpa,matrix_avva,0.0_dp,matrix_ev2,error)
    CALL mat_muld(matrix_ev2,matrix_apVpa,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL mat_mulm(matrix_ev2,matrix_ava,  matrix_avva,n,    -1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_ava,matrix_apVVpa,1.0_dp,matrix_ev2,error)

    !-----------------------------------------------------------------------
    !     2.   1/2 [W1,O1] = 1/2 W1O1 -  1/2 O1W1  
    !-----------------------------------------------------------------------

    CALL cp_fm_scale_and_add(-0.5_dp,matrix_ev2,0.5_dp,matrix_aux2,error)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!  
    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_aux2,   error=error)
    CALL cp_fm_release(matrix_ava,    error=error)
    CALL cp_fm_release(matrix_avva,   error=error)
    CALL cp_fm_release(matrix_apVpa,  error=error)
    CALL cp_fm_release(matrix_apVVpa, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)     

    WRITE (*,*) "CAW:  DKH2 with even2c (Alex)"
    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even2c

  !-----------------------------------------------------------------------

  SUBROUTINE even3b (n,matrix_ev3,matrix_ev1,matrix_pe1p,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified:  20.2.2002 - DKH3                 *
    !                                                                      *
    !     3rd order DK-approximation (generalised DK-transformation)       *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (7.2.2002)                                         *
    !                                                                      *
    !     ev3 = 1/2 [W1,[W1,E1]]                                           *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev3     out  even3 output matrix = final result                  *
    !     e1      in   E1 = even1-operator                                 *
    !     pe1p    in   pE1p                                                *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation
    !     vh      symmetric (n x n)-matrix containing (A V~ A)
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)
    !     e1      E1
    !     pe1p    pE1p
    !     w1w1    (W1)^2
    !     w1e1w1  W1*E1*W1
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,2)
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_ev3, matrix_ev1, &
                                                matrix_pe1p, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER                :: matrix_apVVpa, matrix_aux2, &
                                                matrix_avva, matrix_w1e1w1, &
                                                matrix_w1w1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH3
!-----------------------------------------------------------------------

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    NULLIFY(matrix_aux2)
    NULLIFY(matrix_w1w1)
    NULLIFY(matrix_w1e1w1)
    NULLIFY(matrix_avva)
    NULLIFY(matrix_apVVpa)

    CALL cp_fm_create(matrix_aux2,    matrix_full, error=error)
    CALL cp_fm_create(matrix_w1w1,    matrix_full, error=error)
    CALL cp_fm_create(matrix_w1e1w1,  matrix_full, error=error)
    CALL cp_fm_create(matrix_avva,    matrix_full, error=error)
    CALL cp_fm_create(matrix_apVVpa,  matrix_full, error=error)

    CALL cp_fm_to_fm(matrix_v,matrix_avva,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVVpa,error)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h(matrix_v,matrix_avva,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_avva,matrix_aux2,error)
    CALL mat_axa(matrix_aux2,matrix_avva,n,aa,matrix_aux,error)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h(matrix_pVp,matrix_apVVpa,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_apVVpa,matrix_aux2,error)
    CALL mat_arxra(matrix_aux2,matrix_apVVpa,n,aa,rr,matrix_aux,error)

    !  Calculate w1w1:

    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_apVVpa,matrix_avva,0.0_dp,matrix_w1w1,error)
    CALL mat_muld(matrix_w1w1,matrix_apVVpa,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_aux2,error)
    CALL mat_mulm(matrix_w1w1,matrix_avva,  matrix_avva,n,    -1.0_dp,1.0_dp,tt,rr,matrix_aux2,error)
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_avva,matrix_apVVpa,1.0_dp,matrix_w1w1,error)

    !  Calculate w1e1w1: (warning: ev3 is scratch array)

    CALL mat_muld(matrix_aux,matrix_apVVpa,matrix_pe1p,n,  1.0_dp,0.0_dp,tt,rr,matrix_aux2,error)
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_avva,matrix_pe1p,0.0_dp,matrix_aux2,error)
    CALL cp_fm_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_avva,0.0_dp,matrix_w1e1w1,error)
    CALL mat_muld(matrix_w1e1w1,matrix_aux,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_ev3,error)
    CALL cp_fm_gemm("N","N",n,n,n,-1.0_dp,matrix_aux2,matrix_avva,1.0_dp,matrix_w1e1w1,error)
    CALL mat_muld(matrix_w1e1w1,matrix_aux2,matrix_apVVpa,n,  1.0_dp,1.0_dp,tt,rr,matrix_ev3,error)

    !-----------------------------------------------------------------------
    !     2.   ev3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
    !-----------------------------------------------------------------------

    CALL cp_fm_gemm("N","N",n,n,n,0.5_dp,matrix_w1w1,matrix_ev1,0.0_dp,matrix_ev3,error)
    CALL cp_fm_gemm("N","N",n,n,n,0.5_dp,matrix_ev1,matrix_w1w1,1.0_dp,matrix_ev3,error)
    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev3,-1.0_dp,matrix_w1e1w1,error)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_aux2,   error=error)
    CALL cp_fm_release(matrix_avva,   error=error)
    CALL cp_fm_release(matrix_apVVpa, error=error)
    CALL cp_fm_release(matrix_w1w1,   error=error)
    CALL cp_fm_release(matrix_w1e1w1, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)

    WRITE (*,*) "CAW:  DKH3 with even3b (Alex)"
    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even3b

  !-----------------------------------------------------------------------

  SUBROUTINE even4a (n,ev4,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf,   last modified: 25.02.2002   --   DKH4          *
    !                                                                      *
    !     4th order DK-approximation (scalar = spin-free)                  *
    !                                                                      *
    !     Version: 1.2  (25.2.2002) :  Elegant (short) way of calculation  *
    !                                  included                            *
    !              1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (8.2.2002)                                         *
    !                                                                      *
    !     ev4  =  1/2 [W2,[W1,E1]] + 1/8 [W1,[W1,[W1,O1]]]  =              *
    !                                                                      *
    !          =      sum_1        +         sum_2                         *
    !                                                                      *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                (not working yet)                                     *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev4     out  even4 output matrix = final result                  *
    !     e1     in   E1                                                   *
    !     pe1p   in   p(E1)p                                               *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     v       symmetric (n x n)-matrix containing (A V A)              *
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     w1w1    (W1)^2                                                   *
    !     w1o1    W1*O1      (2-component formulation)                     *
    !     o1w1    O1*W1      (2-component formulation)                     *
    !     e1      symmetric (n x n)-matrix containing E1                   *
    !     pe1p    symmetric (n x n)-matrix containing p(E1)p               *
    !     sum_i   2 addends defined above  (i=1,2)                         *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,..,4)            *
    !     scrh_i  temp. (n x n)-scr.-mat. with energy-denom. (i=1,..,4)    *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev4
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: o1w1, pvp, pvph, scr_1, &
      scr_2, scr_3, scr_4, scrh_1, scrh_2, scrh_3, scrh_4, sum_1, sum_2, v, &
      vh, w1o1, w1w1

!C-----------------------------------------------------------------------
!C     1.   General Structures and Patterns for DKH4
!C-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev4=0.0_dp
    !  Calculate  v = A V A:

    !     CALL mat_axa(v,n,aa)

    !  Calculate  pvp = A P V P A:

    !     CALL mat_arxra(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    !     CALL mat_1_over_h(vh,n,e)
    !     CALL mat_axa(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    !     CALL mat_1_over_h(pvph,n,e)
    !     CALL mat_arxra(pvph,n,aa,rr)


    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    w1w1 = 0.0_dp
    ALLOCATE(w1o1(n,n))
    w1o1 = 0.0_dp
    ALLOCATE(o1w1(n,n))
    o1w1 = 0.0_dp
    ALLOCATE(sum_1(n,n))
    sum_1 = 0.0_dp
    ALLOCATE(sum_2(n,n))
    sum_2 = 0.0_dp
    ALLOCATE(scr_1(n,n))
    scr_1 = 0.0_dp
    ALLOCATE(scr_2(n,n))
    scr_2 = 0.0_dp
    ALLOCATE(scr_3(n,n))
    scr_3 = 0.0_dp
    ALLOCATE(scr_4(n,n))
    scr_4 = 0.0_dp
    ALLOCATE(scrh_1(n,n))
    scrh_1 = 0.0_dp
    ALLOCATE(scrh_2(n,n))
    scrh_2 = 0.0_dp
    ALLOCATE(scrh_3(n,n))
    scrh_3 = 0.0_dp
    ALLOCATE(scrh_4(n,n))
    scrh_4 = 0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    !      CALL mat_muld(w1w1,pvph,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(w1w1,vh,  vh,n,   -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    !      CALL mat_muld(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    !      CALL mat_muld(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !   2. sum_1 = 1/2 [W2,[W1,E1]] = 1/2 (W2W1E1 - W2E1W1 - W1E1W2 + E1W1W2)
    !-----------------------------------------------------------------------

    !  scr_i & scrh_i  for steps 2a (W2W1E1)  and 2b (W2E1W1):

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scr_1,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    !      CALL mat_muld(scr_4, pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)

    !      CALL mat_muld(scrh_1,pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    !      CALL mat_1_over_h(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scrh_2,n)
    !      CALL mat_1_over_h(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scrh_3,n)
    !      CALL mat_1_over_h(scrh_3,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scrh_4,n)
    !      CALL mat_1_over_h(scrh_4,n,e)

    !  2a)  sum_1 = 1/2 W2W1E1               ( [1]-[8] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_1,n,scr_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_1,scr_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_2,n,scr_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_2,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_3,n,scr_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_3,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(sum_1,scrh_4,scr_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_2,n,1.0_dp,sum_1,n)


    !  2b)  sum_1 = - 1/2 W2E1W1 (+ sum_1)   ( [9]-[16] )

    !      CALL mat_muld(sum_1,scrh_1,scr_3,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_1,scr_4,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_2,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_2,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_3,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_3,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_4,n,scr_4,n,1.0_dp,sum_1,n)


    !  scr_i & scrh_i  for steps 2c (W1E1W2)  and 2d (E1W1W2):

    !      CALL mat_muld(scr_1, pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scr_3,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scr_4,n)

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scrh_1,n)
    !      CALL mat_1_over_h(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scrh_2,n)
    !      CALL mat_1_over_h(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    !      CALL mat_1_over_h(scrh_3,n,e)
    !      CALL mat_muld(scrh_4,pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)
    !      CALL mat_1_over_h(scrh_4,n,e)

    !  2c)  sum_1 = - 1/2 W1E1W2 (+ sum_1)   ( [17]-[24] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scr_1,n,scrh_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_1,scrh_2,n,-0.5_dp,1.0_dp,tt,rr) 
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_2,n,scrh_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_2,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_1,scrh_3,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_1,scrh_4,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_2,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_2,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)

    !  2d)  sum_1 = 1/2 E1W1W2 (+ sum_1)     ( [25]-[32] )

    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_3,n,scrh_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_3,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(sum_1,scr_4,scrh_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_2,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_3,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_3,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scr_4,n,scrh_4,n,1.0_dp,sum_1,n)


    !-----------------------------------------------------------------------
    !   3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =
    !
    !             = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.125_dp,w1w1,n,w1o1,n,0.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.375_dp,w1w1,n,o1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,0.375_dp,w1o1,n,w1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.125_dp,o1w1,n,w1w1,n,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   4.  result = sum_1 + sum_2
    !-----------------------------------------------------------------------

    CALL mat_add(ev4,1.0_dp,sum_1,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   5. Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(v,pvp,vh,pvph,w1w1,w1o1,o1w1,sum_1,sum_2)
    DEALLOCATE(scr_1,scr_2,scr_3,scr_4,scrh_1,scrh_2,scrh_3,scrh_4)

    WRITE (*,*) "CAW:  DKH4 with even4a (Alex)"
    WRITE (*,*) "JT:   Now available in cp2k" 

    RETURN
  END SUBROUTINE even4a

  !-----------------------------------------------------------------------
  !                                                                      -
  !     Matrix routines for DKH-procedure                                -
  !     Alexander Wolf                                                   -
  !     modifed: Jens Thar: Mem manager deleted                          -
  !     This file contains the                                           -
  !      following subroutines:                                          -     
  !                                 1. mat_1_over_h                      -
  !                                 2. mat_axa                           -
  !                                 3. mat_arxra                         -
  !                                 4. mat_mulm                          -
  !                                 5. mat_muld                          -
  !                                 6. mat_add                           -
  !                                                                      -
  !-----------------------------------------------------------------------

  SUBROUTINE mat_1_over_h (matrix_p,matrix_pp,n,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !   2. SR mat_1_over_h: Transform matrix p into matrix p/(e(i)+e(j))   *
    !                                                                      *
    !   p    in  REAL(:,:) :   matrix p                                    *
    !   e    in  REAL(:)   :   rel. energy (diagonal)                      *
    !   n    in  INTEGER                                                   *
    !                                                                      *
    !***********************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_pp
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, j, ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices

    CALL cp_fm_get_info(matrix_aux, nrow_local=nrow_local, ncol_local=ncol_local,&
         row_indices=row_indices,col_indices=col_indices,error=error)

    DO j=1,ncol_local
       DO i=1,nrow_local
          matrix_aux%local_data(i,j)=1/(e(row_indices(i))+e(col_indices(j)))
       ENDDO
    ENDDO


    CALL cp_fm_schur_product(matrix_p,matrix_aux,matrix_pp,error)

    RETURN

  END SUBROUTINE mat_1_over_h
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_axa (matrix_x,matrix_axa,n,a,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   3. SR mat_axa: Transform matrix p into matrix  a*p*a               *
    !C                                                                      *
    !C   p    in  REAL(:,:):   matrix p                                     *
    !C   a    in  REAL(:)  :   A-factors (diagonal)                         *
    !CJT n    in  INTEGER  :   dimension of matrix p                        *  
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_axa
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    CALL cp_fm_struct_get(matrix_x%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    CALL cp_fm_get_info(matrix_x, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = a(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_axa,error)
    CALL cp_fm_schur_product(matrix_x,matrix_aux,matrix_axa,error)

    !     DO i=1,n
    !       DO j=1,n
    !          p(i,j)=p(i,j)*a(i)*a(j)
    !       ENDDO
    !     ENDDO


    CALL cp_fm_release(vec_a,error)
    CALL cp_fm_struct_release(vec_full,error)

    RETURN
  END SUBROUTINE mat_axa
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_arxra (matrix_x,matrix_axa,n,a,r,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   4. SR mat_arxra: Transform matrix p into matrix  a*r*p*r*a         *
    !C                                                                      *
    !C   p    in  REAL(:,:) :   matrix p                                    *
    !C   a    in  REAL(:)   :   A-factors (diagonal)                        *
    !C   r    in  REAL(:)   :   R-factors (diagonal)                        *
    !C   n    in  INTEGER   :   dimension of matrix p                       *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_axa
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, r
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    CALL cp_fm_struct_get(matrix_x%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    CALL cp_fm_get_info(matrix_aux, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)


    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = a(row_indices(i))*r(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_axa,error)
    CALL cp_fm_schur_product(matrix_x,matrix_aux,matrix_axa,error)


    CALL cp_fm_release(vec_a,error)
    CALL cp_fm_struct_release(vec_full,error)

    RETURN
  END SUBROUTINE mat_arxra
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_mulm (matrix_p,matrix_q,matrix_r,n,alpha,beta,t,rr,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   5. SR mat_mulm:  Multiply matrices according to:                   *
    !C                                                                      *
    !C                      p = alpha*q*(..P^2..)*r + beta*p                *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,.):   matrix r                                  *
    !C   n      in   INTEGER  :   dimension n of matrices                   *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_q, matrix_r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(n)              :: vec

    CALL cp_fm_to_fm(matrix_q,matrix_aux,error)

    DO i=1,n
       vec(i)=2.0_dp*t(i)*rr(i)*rr(i)
    END DO
    CALL cp_fm_column_scale(matrix_aux,vec)

    CALL cp_fm_gemm("N","N",n,n,n,alpha,matrix_aux,matrix_r,beta,matrix_p,error)
    RETURN

  END SUBROUTINE mat_mulm
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_muld (matrix_p,matrix_q,matrix_r,n,alpha,beta,t,rr,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   16. SR mat_muld:  Multiply matrices according to:                  *
    !C                                                                      *
    !C                      p = alpha*q*(..1/P^2..)*r + beta*p              *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,:):   matrix r                                  *
    !C   n      in   INTEGER  :   Dimension of all matrices                 *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_q, matrix_r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(n)              :: vec

    CALL cp_fm_to_fm(matrix_q,matrix_aux,error)

    DO i=1,n
       vec(i)=0.5_dp/(t(i)*rr(i)*rr(i))
    END DO

    CALL cp_fm_column_scale(matrix_aux,vec)

    CALL cp_fm_gemm("N","N",n,n,n,alpha,matrix_aux,matrix_r,beta,matrix_p,error)

    RETURN

  END SUBROUTINE mat_muld
  !C-----------------------------------------------------------------------


  SUBROUTINE DKH_atom_transformation (s,v,h,pVp,n,dkh_order) 

    !-----------------------------------------------------------------------
    !                                                                      *
    !  INPUT:                                                              *
    !    n          Number of primitive gaussians                          *
    !    s    (:,:) overlap matrix                                         *
    !    pVp  (:,:) pVp matrix                                             *
    !                                                                      *
    !  IN_OUT:                                                             *
    !    v    (:,:) input: nonrelativistic potential energy matrix         *
    !               output: (ev1+ev2)                                      *
    !    h    (:,:) input: kinetic energy matrix                           *
    !               output: kinetic part of hamiltonian in position space  *
    !                                                                      *
    !  INTERNAL                                                            *
    !    sinv (:,:) inverted, orthogonalized overlap matrix                *
    !    ev0t (:)   DKH-even0 matrix in T-basis                            *
    !    e    (:)   e=SQRT(p^2c^2+c^4)                                     *
    !    eig  (:,:) eigenvectors of sinv' h sinv                           *
    !    tt   (:)   eigenvalues of sinv' h sinv                            *
    !    revt (:,:) reverse transformation matrix T-basis -> position space*
    !    aa   (:)   kinematical factors f. DKH SQRT((c^2+e(i))/(2.0*e(i))) *
    !    rr   (:)   kinematical factors f. DKH c/(c^2+e(i))                *
    !    vt   (:,:) non relativistic potential matrix in T-basis           *
    !    pvpt (:,:) pvp integral matrix in T-basis                         * 
    !    ev1t (:,:) DKH-even1 matrix in T-basis                            *
    !    evt2 (:,:) DKH-even2 matrix in T-basis                            *
    !    ev1  (:,:) DKH-even1 matrix in position space                     *
    !    ev2  (:,:) DKH-even2 matrix in position space                     *
    !    ove (:,:) scratch                                                 *
    !    aux (:,:) scratch                                                 *
    !    velit  velocity of light 137 a.u.                                 *
    !    prea   prefactor, 1/137^2                                         *
    !    con2   prefactor, 2/137^2                                         *
    !    con    prefactor, 137^2                                           *
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: s, v, h, pVp
    INTEGER, INTENT(IN)                      :: n, dkh_order

    INTEGER                                  :: i, j, k
    REAL(KIND=dp)                            :: velit
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aa, e, ev0t, rr, tt
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: aux, eig, ev1, ev1t, ev2, &
                                                ev2t, ev3, ev3t, ev4, ev4t, &
                                                ove, pev1tp, pVpt, revt, &
                                                sinv, vt

    IF (dkh_order.lt.0) RETURN     

    !CAW  pp: p^2-values (in momentum-space), stored as matrix!!
    !-----------------------------------------------------------------------
    !     Define velocity of light
    !-----------------------------------------------------------------------
    !      velit = 137.0359895_dp

    velit = 1._dp/a_fine

    !-----------------------------------------------------------------------
    !     Allocate some matrices
    !-----------------------------------------------------------------------

    ALLOCATE(eig(n,n))
    ALLOCATE(sinv(n,n))
    ALLOCATE(revt(n,n))
    ALLOCATE(aux(n,n))
    ALLOCATE(ove(n,n)) 
    ALLOCATE(ev0t(n))
    ALLOCATE(e(n))
    ALLOCATE(aa(n))
    ALLOCATE(rr(n))
    ALLOCATE(tt(n))
    ALLOCATE(ev1t(n,n))
    ALLOCATE(ev2t(n,n))
    ALLOCATE(ev3t(n,n))
    ALLOCATE(ev4t(n,n))
    ALLOCATE(vt(n,n))
    ALLOCATE(pVpt(n,n))
    ALLOCATE(pev1tp(n,n)) 
    ALLOCATE(ev1(n,n))
    ALLOCATE(ev2(n,n))
    ALLOCATE(ev3(n,n))
    ALLOCATE(ev4(n,n))

    !-----------------------------------------------------------------------
    !     Schmidt-orthogonalize overlap matrix
    !-----------------------------------------------------------------------

    CALL sog (n,s,sinv)      

    !-----------------------------------------------------------------------
    !     Calculate matrix representation from nonrelativistic T matrix
    !-----------------------------------------------------------------------

    CALL dkh_diag ( h,n,eig,tt,sinv,aux,0 )

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of Hamiltonian in T-basis 
    !-----------------------------------------------------------------------

    CALL kintegral_a (n,ev0t,tt,e,velit)

    !-----------------------------------------------------------------------
    !     Calculate reverse transformation matrix revt
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,1.0_dp,sinv,n,eig,n,0.0_dp,aux,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,s,n,aux,n,0.0_dp,revt,n)      

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of the Hamiltonian  
    !-----------------------------------------------------------------------

    h = 0.0_dp
    DO i=1,n
       DO j=1,i
          DO k=1,n
             h(i,j)=h(i,j)+revt(i,k)*revt(j,k)*ev0t(k)
             h(j,i)=h(i,j)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------------------
    !     Calculate kinematical factors for DKH
    !-----------------------------------------------------------------------

    DO i=1,n
       aa(i)=SQRT((velit*velit+e(i)) / (2.0_dp*e(i)))            
       rr(i)=SQRT(velit*velit)/(velit*velit+e(i))                        
    END DO

    !----------------------------------------------------------------------- 
    !     Transform v integrals to T-basis (v -> vt)
    !-----------------------------------------------------------------------

    CALL trsm(v,sinv,ove,n,aux)
    CALL trsm(ove,eig,vt,n,aux)

    !-----------------------------------------------------------------------
    !     Transform pVp integrals to T-basis (pVp -> pVpt)
    !-----------------------------------------------------------------------

    CALL trsm(pVp,sinv,ove,n,aux)
    CALL trsm(ove,eig,pVpt,n,aux)

    !-----------------------------------------------------------------------
    !     Calculate even1 in T-basis
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.1) THEN       
       CALL even1_a(n,ev1t,vt,pvpt,aa,rr)

       !----------------------------------------------------------------------
       !     Transform even1 back to position space
       !----------------------------------------------------------------------

       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev1t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev1,n)  
    END IF

    !-----------------------------------------------------------------------
    !     Calculate even2 in T-basis
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.2) THEN       
       CALL even2c_a (n,ev2t,vt,pvpt,aa,rr,tt,e)

       !-----------------------------------------------------------------------
       !     Transform even2 back to position space
       !-----------------------------------------------------------------------

       aux=0.0_dp
       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev2t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev2,n)
    END IF

    !-----------------------------------------------------------------------
    !     Calculate even3 in T-basis, only if requested
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.3) THEN       
       CALL peven1p_a(n,pev1tp,vt,pvpt,aa,rr,tt)
       CALL even3b_a(n,ev3t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)

       !-----------------------------------------------------------------------
       !     Transform even3 back to position space
       !-----------------------------------------------------------------------
       aux=0.0_dp
       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev3t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev3,n)

       !-----------------------------------------------------------------------
       !     Calculate even4 in T-basis, only if requested
       !-----------------------------------------------------------------------

       IF (dkh_order.ge.4) THEN
          CALL even4a_a(n,ev4t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)

          !-----------------------------------------------------------------------
          !     Transform even4 back to position space
          !-----------------------------------------------------------------------
          aux=0.0_dp
          CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev4t,n,0.0_dp,aux,n)
          CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev4,n)
       END IF
    END IF

    IF (dkh_order.ge.4) THEN       
       STOP "DKH 4"
    END IF
    !-----------------------------------------------------------------------
    !     Calculate v in position space
    !-----------------------------------------------------------------------

    IF(dkh_order.ge.1) THEN
       CALL mat_add2(v,0.0_dp,1.0_dp,ev1,n)
    END IF
    IF(dkh_order.ge.2) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev2,n)
    END IF
    IF(dkh_order.ge.3) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev3,n)
    END IF
    IF(dkh_order.ge.4) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev4,n)
    END IF

    !-----------------------------------------------------------------------

    DEALLOCATE(eig,sinv,revt,ove,aux,vt,pVpt,ev1,ev2,ev3,ev4,ev1t,ev2t,ev3t,ev4t,pev1tp) 
    DEALLOCATE(ev0t,e,aa,rr,tt)

    RETURN
  END SUBROUTINE dkh_atom_transformation

  !-----------------------------------------------------------------------


  SUBROUTINE kintegral_a (n,ev0t,tt,e,velit)

    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ev0t
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tt
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: e
    REAL(KIND=dp), INTENT(IN)                :: velit

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: con, con2, prea, ratio, tv1, &
                                                tv2, tv3, tv4

!-----------------------------------------------------------------------

    DO i=1,n
       IF (tt(i).LT.0.0_dp) THEN
          WRITE(*,*) ' dkh_main.F | tt(',i,') = ',tt(i)
       END IF

       !       Calculate some constants

       prea=1/(velit*velit)
       con2=prea+prea
       con=1.0_dp/prea

       !       If T is sufficiently small, use series expansion to avoid
       !       cancellation, otherwise calculate SQRT directly

       ev0t(i)=tt(i)
       ratio=tt(i)/velit
       IF (ratio.LE.0.02_dp) THEN
          tv1=tt(i)
          tv2=-tv1*tt(i)*prea/2.0_dp
          tv3=-tv2*tt(i)*prea
          tv4=-tv3*tt(i)*prea*1.25_dp
          ev0t(i)=tv1+tv2+tv3+tv4
       ELSE
          ev0t(i)=con*(SQRT(1.0_dp+con2*tt(i))-1.0_dp)   
       END IF
       e(i)=ev0t(i)+con           
    END DO

    RETURN
  END SUBROUTINE kintegral_a


  !-----------------------------------------------------------------------

  SUBROUTINE even1_a(n,ev1t,vt,pvpt,aa,rr)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n    in   dimension of matrices                                  -
    !     ev1t out  even1 output matrix                                    -
    !     vt   in   potential matrix v in T-space                          -
    !     pvpt in   pvp matrix in T-space                                  -
    !     aa   in   A-factors (diagonal)                                   -
    !     rr   in   R-factors (diagonal)                                   -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev1t
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vt, pvpt
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr

    INTEGER                                  :: i, j

!-----------------------------------------------------------------------

    DO i=1,n
       DO j=1,i
          ev1t(i,j)=vt(i,j)*aa(i)*aa(j)+pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          ev1t(j,i)=ev1t(i,j)
       END DO
    END DO

    RETURN
  END SUBROUTINE even1_a

  !-----------------------------------------------------------------------

  SUBROUTINE peven1p_a(n,pev1tp,vt,pvpt,aa,rr,tt)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n      in   dimension of matrices                                -
    !     pev1tp out  peven1p output matrix                                -
    !     vt     in   potential matrix v in T-space                        -
    !     pvpt   in   pvp matrix in T-space                                -
    !     aa     in   A-factors (diagonal)                                 -
    !     rr     in   R-factors (diagonal)                                 -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pev1tp
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vt, pvpt
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt

    INTEGER                                  :: i, j

!-----------------------------------------------------------------------

    DO i=1,n
       DO j=1,i
          pev1tp(i,j)=4.0_dp*vt(i,j)*aa(i)*aa(j)*rr(i)*rr(i)*rr(j)*rr(j)*tt(i)*tt(j)+&
               pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          pev1tp(j,i)=pev1tp(i,j)
       END DO
    END DO

    RETURN
  END SUBROUTINE peven1p_a

  !-----------------------------------------------------------------------

  SUBROUTINE even2c_a (n,ev2,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified: 20.02.2002 - DKH2                 *
    !                                                                      *
    !     2nd order DK-approximation ( original DK-transformation with     *
    !                                       U = SQRT(1+W^2) + W        )   *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (6.2.2002)                                         *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !     ev2 = 1/2 [W1,O1]                                                *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev2     out  even2 output matrix = final result                  *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation                 *
    !     v       symmetric (n x n)-matrix containing (A V A)              * 
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         * 
    !     w1o1    W1*O1 (2-component form)                                 * 
    !     o1w1    O1*W1 (2-component form)                                 *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: o1w1, pvp, pvph, v, vh, w1o1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH2                      
!-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev2=0.0_dp
    !  Calculate  v = A V A:

    CALL mat_axa_a(v,n,aa)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra_a(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)

    !  Create/Initialize necessary matrices:
    ALLOCATE(w1o1(n,n))
    ALLOCATE(o1w1(n,n))
    w1o1=0.0_dp
    o1w1=0.0_dp

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    CALL mat_muld_a(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    CALL mat_muld_a(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)
    !  Calculate in symmetric pakets

    !-----------------------------------------------------------------------
    !     2.   1/2 [W1,O1] = 1/2 W1O1 -  1/2 O1W1  
    !-----------------------------------------------------------------------

    CALL mat_add (ev2,0.5_dp,w1o1,-0.5_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!  
    !-----------------------------------------------------------------------

    DEALLOCATE(v,vh,pvp,pvph,w1o1,o1w1)

    WRITE (*,*) "CAW:  DKH2 with even2c (Alex)"
    WRITE (*,*) "!JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even2c_a

  !-----------------------------------------------------------------------

  SUBROUTINE even3b_a (n,ev3,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified:  20.2.2002 - DKH3                 *
    !                                                                      *
    !     3rd order DK-approximation (generalised DK-transformation)       *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (7.2.2002)                                         *
    !                                                                      *
    !     ev3 = 1/2 [W1,[W1,E1]]                                           *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev3     out  even3 output matrix = final result                  *
    !     e1      in   E1 = even1-operator                                 *
    !     pe1p    in   pE1p                                                *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation                 * 
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     e1      E1                                                       * 
    !     pe1p    pE1p                                                     *
    !     w1w1    (W1)^2                                                   *
    !     w1e1w1  W1*E1*W1                                                 *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,2)               * 
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev3
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pvph, scr_1, scr_2, vh, &
                                                w1e1w1, w1w1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH3
!-----------------------------------------------------------------------

    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    vh=0.0_dp
    pVph=0.0_dp
    vh(1:n,1:n)=vv(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev3=0.0_dp

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)

    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    ALLOCATE(w1e1w1(n,n))
    ALLOCATE(scr_1(n,n))
    ALLOCATE(scr_2(n,n))
    w1w1=0.0_dp
    w1e1w1=0.0_dp
    scr_1=0.0_dp
    scr_2=0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    CALL mat_muld_a(w1w1,pvph,pvph,n,-1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1w1,vh,  vh,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1e1w1:
    CALL mat_muld_a(scr_1 ,pvph ,pe1p,n, 1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,scr_1,n,vh,n,0.0_dp,w1e1w1,n)
    CALL mat_muld_a(w1e1w1,scr_1,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,scr_2,n,vh,n,1.0_dp,w1e1w1,n)
    CALL mat_muld_a(w1e1w1,scr_2,pvph,n, 1.0_dp,1.0_dp,tt,rr)


    !-----------------------------------------------------------------------
    !     2.   ev3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.5_dp,w1w1,n,e1,n,0.0_dp,ev3,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,e1,n,w1w1,n,1.0_dp,ev3,n)
    CALL mat_add2 (ev3,1.0_dp,-1.0_dp,w1e1w1,n)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(vh,pvph,w1w1,w1e1w1,scr_1,scr_2)

    WRITE (*,*) "CAW:  DKH3 with even3b (Alex)"
    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even3b_a

  !-----------------------------------------------------------------------

  SUBROUTINE even4a_a (n,ev4,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf,   last modified: 25.02.2002   --   DKH4          *
    !                                                                      *
    !     4th order DK-approximation (scalar = spin-free)                  *
    !                                                                      *
    !     Version: 1.2  (25.2.2002) :  Elegant (short) way of calculation  *
    !                                  included                            *
    !              1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (8.2.2002)                                         *
    !                                                                      *
    !     ev4  =  1/2 [W2,[W1,E1]] + 1/8 [W1,[W1,[W1,O1]]]  =              *
    !                                                                      *
    !          =      sum_1        +         sum_2                         *
    !                                                                      *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev4     out  even4 output matrix = final result                  *
    !     e1     in   E1                                                   *
    !     pe1p   in   p(E1)p                                               *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     v       symmetric (n x n)-matrix containing (A V A)              *
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     w1w1    (W1)^2                                                   *
    !     w1o1    W1*O1      (2-component formulation)                     *
    !     o1w1    O1*W1      (2-component formulation)                     *
    !     e1      symmetric (n x n)-matrix containing E1                   *
    !     pe1p    symmetric (n x n)-matrix containing p(E1)p               *
    !     sum_i   2 addends defined above  (i=1,2)                         *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,..,4)            *
    !     scrh_i  temp. (n x n)-scr.-mat. with energy-denom. (i=1,..,4)    *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev4
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: o1w1, pvp, pvph, scr_1, &
      scr_2, scr_3, scr_4, scrh_1, scrh_2, scrh_3, scrh_4, sum_1, sum_2, v, &
      vh, w1o1, w1w1

!C-----------------------------------------------------------------------
!C     1.   General Structures and Patterns for DKH4
!C-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev4=0.0_dp
    !  Calculate  v = A V A:

    CALL mat_axa_a(v,n,aa)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra_a(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)


    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    w1w1 = 0.0_dp
    ALLOCATE(w1o1(n,n))
    w1o1 = 0.0_dp
    ALLOCATE(o1w1(n,n))
    o1w1 = 0.0_dp
    ALLOCATE(sum_1(n,n))
    sum_1 = 0.0_dp
    ALLOCATE(sum_2(n,n))
    sum_2 = 0.0_dp
    ALLOCATE(scr_1(n,n))
    scr_1 = 0.0_dp
    ALLOCATE(scr_2(n,n))
    scr_2 = 0.0_dp
    ALLOCATE(scr_3(n,n))
    scr_3 = 0.0_dp
    ALLOCATE(scr_4(n,n))
    scr_4 = 0.0_dp
    ALLOCATE(scrh_1(n,n))
    scrh_1 = 0.0_dp
    ALLOCATE(scrh_2(n,n))
    scrh_2 = 0.0_dp
    ALLOCATE(scrh_3(n,n))
    scrh_3 = 0.0_dp
    ALLOCATE(scrh_4(n,n))
    scrh_4 = 0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    CALL mat_muld_a(w1w1,pvph,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1w1,vh,  vh,n,   -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    CALL mat_muld_a(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    CALL mat_muld_a(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !   2. sum_1 = 1/2 [W2,[W1,E1]] = 1/2 (W2W1E1 - W2E1W1 - W1E1W2 + E1W1W2)
    !-----------------------------------------------------------------------

    !  scr_i & scrh_i  for steps 2a (W2W1E1)  and 2b (W2E1W1):

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scr_1,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    CALL mat_muld_a(scr_4, pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)

    CALL mat_muld_a(scrh_1,pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL mat_1_over_h_a(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scrh_2,n)
    CALL mat_1_over_h_a(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scrh_3,n)
    CALL mat_1_over_h_a(scrh_3,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scrh_4,n)
    CALL mat_1_over_h_a(scrh_4,n,e)

    !  2a)  sum_1 = 1/2 W2W1E1               ( [1]-[8] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_1,n,scr_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_1,scr_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_2,n,scr_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_2,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_3,n,scr_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_3,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(sum_1,scrh_4,scr_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_2,n,1.0_dp,sum_1,n)


    !  2b)  sum_1 = - 1/2 W2E1W1 (+ sum_1)   ( [9]-[16] )

    CALL mat_muld_a(sum_1,scrh_1,scr_3,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_1,scr_4,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_2,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_2,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_3,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_3,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_4,n,scr_4,n,1.0_dp,sum_1,n)


    !  scr_i & scrh_i  for steps 2c (W1E1W2)  and 2d (E1W1W2):

    CALL mat_muld_a(scr_1, pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scr_3,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scr_4,n)

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scrh_1,n)
    CALL mat_1_over_h_a(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scrh_2,n)
    CALL mat_1_over_h_a(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    CALL mat_1_over_h_a(scrh_3,n,e)
    CALL mat_muld_a(scrh_4,pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)
    CALL mat_1_over_h_a(scrh_4,n,e)

    !  2c)  sum_1 = - 1/2 W1E1W2 (+ sum_1)   ( [17]-[24] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scr_1,n,scrh_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_1,scrh_2,n,-0.5_dp,1.0_dp,tt,rr) 
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_2,n,scrh_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_2,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_1,scrh_3,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_1,scrh_4,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_2,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_2,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)

    !  2d)  sum_1 = 1/2 E1W1W2 (+ sum_1)     ( [25]-[32] )

    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_3,n,scrh_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_3,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(sum_1,scr_4,scrh_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_2,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_3,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_3,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scr_4,n,scrh_4,n,1.0_dp,sum_1,n)


    !-----------------------------------------------------------------------
    !   3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =
    !
    !             = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.125_dp,w1w1,n,w1o1,n,0.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.375_dp,w1w1,n,o1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,0.375_dp,w1o1,n,w1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.125_dp,o1w1,n,w1w1,n,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   4.  result = sum_1 + sum_2
    !-----------------------------------------------------------------------

    CALL mat_add(ev4,1.0_dp,sum_1,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   5. Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(v,pvp,vh,pvph,w1w1,w1o1,o1w1,sum_1,sum_2)
    DEALLOCATE(scr_1,scr_2,scr_3,scr_4,scrh_1,scrh_2,scrh_3,scrh_4)

    WRITE (*,*) "CAW:  DKH4 with even4a (Alex)"
    WRITE (*,*) "JT:   Now available in cp2k" 

    RETURN
  END SUBROUTINE even4a_a

  !-----------------------------------------------------------------------
  !                                                                      -
  !     Matrix routines for DKH-procedure                                -
  !     Alexander Wolf                                                   -
  !     modifed: Jens Thar: Mem manager deleted                          -
  !     This file contains the                                           -
  !      following subroutines:                                          -     
  !                                 1. mat_1_over_h                      -
  !                                 2. mat_axa                           -
  !                                 3. mat_arxra                         -
  !                                 4. mat_mulm                          -
  !                                 5. mat_muld                          -
  !                                 6. mat_add                           -
  !                                                                      -
  !-----------------------------------------------------------------------

  SUBROUTINE mat_1_over_h_a (p,n,e)

    !***********************************************************************
    !                                                                      *
    !   2. SR mat_1_over_h: Transform matrix p into matrix p/(e(i)+e(j))   *
    !                                                                      *
    !   p    in  REAL(:,:) :   matrix p                                    *
    !   e    in  REAL(:)   :   rel. energy (diagonal)                      *
    !   n    in  INTEGER                                                   *
    !                                                                      *
    !***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)/(e(i)+e(j))
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_1_over_h_a
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_axa_a (p,n,a)

    !C***********************************************************************
    !C                                                                      *
    !C   3. SR mat_axa: Transform matrix p into matrix  a*p*a               *
    !C                                                                      *
    !C   p    in  REAL(:,:):   matrix p                                     *
    !C   a    in  REAL(:)  :   A-factors (diagonal)                         *
    !CJT n    in  INTEGER  :   dimension of matrix p                        *  
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)*a(i)*a(j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_axa_a
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_arxra_a (p,n,a,r)

    !C***********************************************************************
    !C                                                                      *
    !C   4. SR mat_arxra: Transform matrix p into matrix  a*r*p*r*a         *
    !C                                                                      *
    !C   p    in  REAL(:,:) :   matrix p                                    *
    !C   a    in  REAL(:)   :   A-factors (diagonal)                        *
    !C   r    in  REAL(:)   :   R-factors (diagonal)                        *
    !C   n    in  INTEGER   :   dimension of matrix p                       *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, r

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)*a(i)*a(j)*r(i)*r(j)   
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_arxra_a
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_mulm_a (p,q,r,n,alpha,beta,t,rr)

    !C***********************************************************************
    !C                                                                      *
    !C   5. SR mat_mulm:  Multiply matrices according to:                   *
    !C                                                                      *
    !C                      p = alpha*q*(..P^2..)*r + beta*p                *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,.):   matrix r                                  *
    !C   n      in   INTEGER  :   dimension n of matrices                   *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p, q, r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr

    INTEGER                                  :: i, j
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: qtemp

    ALLOCATE(qtemp(n,n))

    DO i=1,n
       DO j=1,n
          qtemp(i,j)=q(i,j)*2.0_dp*t(j)*rr(j)*rr(j)
       ENDDO
    ENDDO

    CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)
    RETURN
    DEALLOCATE(qtemp)
  END SUBROUTINE mat_mulm_a
  !C-----------------------------------------------------------------------



  SUBROUTINE mat_muld_a (p,q,r,n,alpha,beta,t,rr)

    !C***********************************************************************
    !C                                                                      *
    !C   16. SR mat_muld:  Multiply matrices according to:                  *
    !C                                                                      *
    !C                      p = alpha*q*(..1/P^2..)*r + beta*p              *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,:):   matrix r                                  *
    !C   n      in   INTEGER  :   Dimension of all matrices                 *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p, q, r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr

    INTEGER                                  :: i, j
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: qtemp

    ALLOCATE(qtemp(n,n)) 

    DO i=1,n
       DO j=1,n
          qtemp(i,j)=q(i,j)*0.5_dp/(t(j)*rr(j)*rr(j))
       ENDDO
    ENDDO

    CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)

    RETURN
    DEALLOCATE(qtemp) 
  END SUBROUTINE mat_muld_a
  !C-----------------------------------------------------------------------

  SUBROUTINE mat_add2 (p,alpha,beta,r,n)

    !C***********************************************************************
    !C                                                                      *
    !C   19. SR mat_add:  Add two matrices of the same size according to:   *
    !C                                                                      *
    !C                            p = alpha*p + beta*r                      *
    !C                                                                      *
    !C                    and store them in the first                       *
    !C   p      out  REAL(:,:)  :   matrix p                                *
    !C   r      in   REAL(:,:)  :   matrix r                                *
    !C   alpha  in   REAL(dp)                                               *
    !C   beta   in   REAL(dp)                                               *
    !C                                                                      *
    !C   Matrix p must already exist before calling this SR!!               *
    !C                                                                      *
    !C  [written by: Alexander Wolf,  20.2.2002,  v1.0]                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j

!C  Add matrices:

    DO i=1,n
       DO j=1,n
          p(i,j) = alpha*p(i,j) + beta*r(i,j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_add2

  !CAW---------------------------------------------------------------------


  SUBROUTINE mat_add (p,alpha,q,beta,r,n)

    !C***********************************************************************
    !C                                                                      *
    !C   19. SR mat_add:  Add two matrices of the same size according to:   *
    !C                                                                      *
    !C                            p = alpha*q + beta*r                      *
    !C                                                                      *
    !C   p      out  REAL(:,:)  :   matrix p                                *
    !C   q      in   REAL(:,:)  :   matrix q                                *
    !C   r      in   REAL(:,:)  :   matrix r                                *
    !C   alpha  in   REAL(dp)                                               *
    !C   beta   in   REAL(dp)                                               *
    !C                                                                      *
    !C   Matrix p must already exist before calling this SR!!               *
    !C                                                                      *
    !C  [written by: Alexander Wolf,  20.2.2002,  v1.0]                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: p
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: q
    REAL(KIND=dp), INTENT(IN)                :: beta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j

!C  Add matrices:

    DO i=1,n
       DO j=1,n
          p(i,j) = alpha*q(i,j) + beta*r(i,j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_add

  !CAW---------------------------------------------------------------------

  SUBROUTINE TRSM ( W,B,C,N,H)


    REAL(KIND=dp), DIMENSION(:, :)           :: W, B, C
    INTEGER                                  :: N
    REAL(KIND=dp), DIMENSION(:, :)           :: H

    INTEGER                                  :: I, IJ, J, K, L

!C
!C     TRANSFORM SYMMETRIC matrix A by UNITARY TRANSFORMATION
!C     IN B. RESULT IS IN C
!C
!CAW      C = B^{dagger} * A * B

    IJ=0
    DO I=1,N
       DO J=1,I
          IJ=IJ+1
          C(I,J)=0.0_dp
          C(J,I)=0.0_dp
          H(I,J)=0.0_dp
          H(J,I)=0.0_dp
       END DO
    END DO
    DO I=1,N
       DO L=1,N
          DO K=1,N
             H(I,L)=B(K,I)*W(K,L)+H(I,L)
          END DO
       END DO
    END DO

    IJ=0
    DO I=1,N
       DO J=1,I
          IJ=IJ+1
          DO L=1,N
             C(I,J)=H(I,L)*B(L,J)+C(I,J)
             C(J,I)=C(I,J)
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE TRSM

  !JT----------------------------------------------------------------------

  SUBROUTINE dkh_diag (matrix_t_pgf,n,eig,ew,matrix_sinv_pgf,aux,ic)


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: matrix_t_pgf
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: eig
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: ew
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix_sinv_pgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: aux
    INTEGER                                  :: ic

    INTEGER                                  :: n2

    eig = 0.0_dp
    aux = 0.0_dp

    CALL dgemm("N","N",n,n,n,1.0_dp,matrix_t_pgf,n,matrix_sinv_pgf,n,0.0_dp,eig,n)

    aux = 0.0_dp

    CALL dgemm("T","N",n,n,n,1.0_dp,matrix_sinv_pgf,n,eig,n,0.0_dp,aux,n)

    n2=3*n-1

    CALL JACOB2 ( AUX,EIG,EW,N,IC )
    RETURN
  END SUBROUTINE dkh_diag

  !CAW---------------------------------------------------------------------

  SUBROUTINE JACOB2 ( sogt,eigv,eigw,n,ic )


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: eigw
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(OUT)                            :: eigv
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(INOUT)                          :: sogt
    INTEGER, INTENT(IN)                      :: ic

    INTEGER                                  :: i, il, im, ind, j, k, l, ll, &
                                                m, mm
    REAL(KIND=dp)                            :: cost, cost2, ext_norm, sincs, &
                                                sint, sint2, thr, thr_min, &
                                                tol, u1, x, xy, y

    tol=1.0E-15
    ext_norm=0.0_dp
    u1=REAL(n)
    DO i=1,n
       eigv(i,i)=1.0_dp
       eigw(i)=sogt(i,i)
       DO j=1,i
          IF(i.ne.j) THEN
             eigv(i,j)=0.0_dp
             eigv(j,i)=0.0_dp
             ext_norm=ext_norm+sogt(i,j)*sogt(i,j)
          END IF
       END DO
    END DO

    IF (ext_norm.gt.0.0_dp) THEN 
       ext_norm=SQRT(2.0_dp*ext_norm)
       thr_min=ext_norm*tol/u1    
       ind=0
       thr=ext_norm

       DO
          thr=thr/u1              
          DO
             l=1                       
             DO                   
                m=l+1
                DO                
                   IF ((ABS(sogt(m,l))-thr).ge.0.0_dp) THEN     
                      ind=1
                      x=0.5_dp*(eigw(l)-eigw(m))                       
                      y=-sogt(m,l)/SQRT(sogt(m,l)*sogt(m,l)+x*x)   
                      IF (x.lt.0.0_dp) y=-y                      

                      IF (y.gt.1.0_dp) y=1.0_dp                
                      IF (y.lt.-1.0_dp) y=-1.0_dp                      
                      xy=1.0_dp-y*y                                    
                      sint=y/SQRT(2.0_dp*(1.0_dp+SQRT(xy))) 
                      sint2=sint*sint
                      cost2=1.0_dp-sint2
                      cost=SQRT(cost2)
                      sincs=sint*cost

                      DO i=1,n                               
                         IF((i-m).ne.0) THEN
                            IF ((i-m).lt.0) THEN
                               im=m
                               mm=i
                            ELSE
                               im=i
                               mm=m
                            END IF
                            IF ((i-l).ne.0) THEN
                               IF ((i-l).lt.0) THEN
                                  il=l
                                  ll=i
                               ELSE
                                  il=i
                                  ll=l
                               END IF
                               x=sogt(il,ll)*cost-sogt(im,mm)*sint         
                               sogt(im,mm)=sogt(il,ll)*sint+sogt(im,mm)*cost 
                               sogt(il,ll)=x  
                            END IF
                         END IF

                         x=eigv(i,l)*cost-eigv(i,m)*sint                  
                         eigv(i,m)=eigv(i,l)*sint+eigv(i,m)*cost
                         eigv(i,l)=x
                      END DO

                      x=2.0_dp*sogt(m,l)*sincs                          
                      y=eigw(l)*cost2+eigw(m)*sint2-x
                      x=eigw(l)*sint2+eigw(m)*cost2+x
                      sogt(m,l)=(eigw(l)-eigw(m))*sincs+sogt(m,l)*(cost2-sint2)
                      eigw(l)=y
                      eigw(m)=x
                   END IF
                   IF ((m-n).eq.0) EXIT
                   m=m+1
                END DO
                IF ((l-m+1).eq.0) EXIT
                l=l+1
             END DO
             IF((ind-1).ne.0.0_dp) EXIT
             ind=0
          END DO
          IF ((thr-thr_min).le.0.0_dp) EXIT
       END DO
    END IF

    IF (ic.ne.0) THEN
       DO i=1,n
          DO j=1,n
             IF ((eigw(i)-eigw(j)).gt.0.0_dp) THEN
                x=eigw(i)
                eigw(i)=eigw(j)
                eigw(j)=x
                DO k=1,n
                   y=eigv(k,i)
                   eigv(k,i)=eigv(k,j)
                   eigv(k,j)=y
                END DO
             END IF
          END DO
       END DO

    END IF
    RETURN
  END SUBROUTINE JACOB2

  !JT---------------------------------------------------------------------

  SUBROUTINE SOG (n,matrix_s_pgf,matrix_sinv_pgf)


    INTEGER                                  :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix_s_pgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: matrix_sinv_pgf

    INTEGER                                  :: i, j, jn, k
    REAL(KIND=dp)                            :: diag_s, row_sum, scalar
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: a
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: g

!
!     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
!     ORTHOGONALIZED BASIS.
!     sinv-1*matrix_s_pgf*sinv = "orthogonal matrix"
!     n              dimension of matrices
!     matrix_s_pgf   original overlap matrix
!     matrix_sinv_pgf new overlap matrix 
!     g              scratch
!     a              scratch
!

    ALLOCATE(a(n))
    ALLOCATE(g(n,n))

    DO jn=1,n                                                         
       diag_s = matrix_s_pgf(jn,jn)
       g(jn,jn)=1.0_dp

       IF(jn.ne.1) THEN                                               
          DO j=1,jn-1                                                 
             scalar=0.0_dp
             DO i=1,j                                                 
                scalar=scalar+matrix_s_pgf(i,jn)*g(i,j)               
             END DO
             diag_s=diag_s-scalar*scalar                            
             a(j) = scalar
          END DO

          DO j=1,jn-1                                             
             row_sum=0.0_dp
             DO k=j,jn-1                                          
                row_sum=row_sum+a(k)*g(j,k)
             END DO
             g(j,jn)=-row_sum                                    
          END DO
       END IF

       diag_s=1.0_dp/SQRT(diag_s)                          
       DO i=1,jn
          g(i,jn)=g(i,jn)*diag_s                            
       END DO
    END DO

    DO j=1,n
       DO i=1,j
          matrix_sinv_pgf(j,i)=0.0_dp 
          matrix_sinv_pgf(i,j)=g(i,j)
       END DO
    END DO
    DEALLOCATE(a,g)

    RETURN

  END SUBROUTINE SOG

END MODULE dkh_main
