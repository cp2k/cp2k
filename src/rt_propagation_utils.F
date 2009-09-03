!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_utils


  
  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_add_block_node,&
                                             dbcsr_allocate_matrix,&
                                             dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_obj_p_type,&
                                             dbcsr_replicate_matrix_struct,&
                                             sm_from_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose,&
                                             cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE dbcsr_util,                      ONLY: dbcsr_finalize
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: use_restart_wfn,&
                                             use_rt_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             read_rt_mos_from_restart
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE scp_environment_types,           ONLY: scp_environment_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: calc_S_derivs,&
          get_restart_wfn,&
          calculate_dp_imag,&
          calc_c_mat_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'

!***
  
  CONTAINS



! *****************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD. 
!>      
!> \param qs_env the qs environment 
!> \param error 
!>
!> \par History
!>      02.2009 created [Manuel Guidon]
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE calc_S_derivs(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_S_derivs', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=80)                        :: name
    INTEGER :: atom_a, atom_b, first_col, handle, i, iab, iatom, icol, &
      icol_atom, icol_global, icol_local, ikind, ilist, inode, irow, &
      irow_atom, irow_global, irow_local, iset, j, jatom, jkind, jset, &
      last_jatom, ldai, ldsab, m, maxblock, maxco, maxcoa, maxdco, maxdcoa, &
      maxder, maxl, maxlgto, maxlgtoa, maxlppl, maxlppla, maxsgf, maxsgfa, n, &
      natom, ncoa, ncob, ncol_local, nder, neighbor_list_id, nkind, nlist, &
      nnode, nrow, nrow_local, nseta, nsetb, nsgf, nthread, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, first_sgf, &
                                                last_sgf
    INTEGER, DIMENSION(6, 2)                 :: c_map_mat
    INTEGER, DIMENSION(:), POINTER :: col_atom_index, col_indices, la_max, &
      la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, row_atom_index, &
      row_indices
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, new_atom_b, &
                                                return_s_derivatives
    REAL(KIND=dp)                            :: dab, f, rab2, scale
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data, rpgfa, rpgfb, &
                                                sab, sphi_a, sphi_b, work, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: sint
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: sintt
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: C_mat
    TYPE(cp_fm_type), POINTER                :: B_mat, tmp_mat, tmp_mat2
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: s_derivs
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_derivs_sm
    TYPE(rt_prop_type), POINTER              :: rtp

    failure = .FALSE.
    CALL timeset(routineN,handle)

    return_s_derivatives=.TRUE.

    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY(s_derivs_sm)
    NULLIFY (particle_set)
    NULLIFY (sab)
    NULLIFY(rtp)
    NULLIFY (work)
    NULLIFY(distribution_2d)
    NULLIFY(s_derivs)

    CALL get_qs_env(qs_env=qs_env,&
                    rtp=rtp,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    distribution_2d=distribution_2d,&
                    sab_orb=sab_orb,&
                    force=force,&
                    error=error)

    CALL get_rtp(rtp=rtp,B_mat=B_mat,C_mat=C_mat,error=error)
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    nder = 2
    maxder = ncoset(nder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)

    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
         error=error)

    ! assert here that the S-matrix is of the type distribution_2d if this one is non-null
    ! we assert that the sparsity of s is determined by the neighbor list that is used to build it
    ! is that right ? i.e. this implies that even if forces are computed the sparsity of s remains the same.
    ! this seems necessary anyway.
    CALL dbcsr_allocate_matrix_set(s_derivs,maxder-1,error=error)
    CALL dbcsr_allocate_matrix(matrix=s_derivs(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=first_sgf,&
                         last_row=last_sgf,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name="X DERIVATIVE",&
                         matrix_symmetry="antisymmetric",&
                         sparsity_id=neighbor_list_id, &
                         distribution_2d=distribution_2d,error=error)

    

    DO m=2,3
       cgfsym = cgf_symbol(1,indco(1:3,m+1))
       name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE OVERLAP MATRIX "//&
            "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL dbcsr_replicate_matrix_struct(source=s_derivs(1)%matrix,&
                                      TARGET=s_derivs(m)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="antisymmetric",error=error)
   END DO


    DO m=4,maxder-1
       cgfsym = cgf_symbol(1,indco(1:3,m+1))
       name = TRIM(cgfsym(4:))//" SECOND DERIVATIVE OF THE OVERLAP MATRIX "//&
            "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL dbcsr_replicate_matrix_struct(source=s_derivs(1)%matrix,&
                                      TARGET=s_derivs(m)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="symmetric",error=error)
   END DO
    
    ! *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxcoa,&
                             maxlgto=maxlgtoa,&
                             maxlppl=maxlppla,&
                             maxsgf=maxsgfa)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdcoa,&
                             maxder=MAX(1,nder))

    nder=2
    maxdco = maxdcoa
    maxco = maxcoa
    maxlgto = maxlgtoa
    maxlppl = maxlppla
    maxsgf = maxsgfa
    maxl = maxlgto
    maxder =  ncoset(nder)
    nthread = 1

    CALL init_orbital_pointers(maxl+nder+1)


    ALLOCATE (sintt(SIZE(s_derivs),0:nthread-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO m=1,SIZE(sintt,1)
       DO j=0,nthread-1
          NULLIFY (sintt(m,j)%block)
       END DO
    END DO




    ldsab = MAX(maxco,maxsgf)
    maxblock = MAX(4,maxder) 
    ldai = ncoset(maxl+nder+1)

    ALLOCATE(sab(ldsab,ldsab*maxblock))
    ALLOCATE(work(ldsab,ldsab*maxder))


    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
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

        sint => sintt(:,0)

        NULLIFY ( sab_orb_neighbor_list )
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

          atom_a = atom_of_kind(iatom)

          last_jatom = 0
          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF



            atom_b = atom_of_kind(jatom)

            IF (new_atom_b) THEN
               IF (iatom <= jatom) THEN
                  irow = iatom
                  icol = jatom
               ELSE
                  irow = jatom
                  icol = iatom
               END IF
              DO m=1,SIZE(s_derivs)
                 NULLIFY (sint(m)%block)
                 CALL dbcsr_add_block_node(matrix=s_derivs(m)%matrix,&
                      block_row=irow,&
                      block_col=icol,&
                      BLOCK=sint(m)%block,error=error)
              END DO
            END IF

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)
            nrow = 1
            DO iset=1,nseta
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              DO jset=1,nsetb
                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)
                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN
!                  IF ( iatom /= jatom ) THEN

                    ! *** Calculate the primitive overlap integrals ***
                    ! *** and the corresponding force contribution  ***
                    CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai)         

                    ! *** Contraction step (overlap matrix and its derivatives) ***
                         DO m=1,SIZE(s_derivs)
                            IF ((m > 1).AND.(.NOT.return_s_derivatives)) CYCLE
                            first_col = (m)*SIZE(sab,1) + 1
                            CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_dp,sab(1,first_col),SIZE(sab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                            IF (iatom <= jatom) THEN
                               CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    work(1,1),SIZE(work,1),&
                                    1.0_dp,sint(m)%block(sgfa,sgfb),&
                                    SIZE(sint(m)%block,1))
                            ELSE
                               IF (m < 4) THEN
                                  f = -1.0_dp
                               ELSE
                                  f = 1.0_dp
                               END IF
                               CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                    f,work(1,1),SIZE(work,1),&
                                    sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    1.0_dp,sint(m)%block(sgfb,sgfa),&
                                    SIZE(sint(m)%block,1))
                            END IF
                         END DO
!                      END IF
                END IF
              END DO
              nrow = nrow + ncoa
            END DO
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO
        END DO
        DEALLOCATE (ai_work,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO
    END DO

    CALL allocate_matrix_set(s_derivs_sm,SIZE(s_derivs),error=error)
    DO i = 1,SIZE(s_derivs)
       CALL dbcsr_finalize(s_derivs(i)%matrix,error=error)
       CALL sm_from_dbcsr ( s_derivs_sm(i)%matrix, s_derivs(i)%matrix, distribution_2d, error=error )
    ENDDO
    CALL dbcsr_deallocate_matrix_set ( s_derivs, error=error )
 
    DEALLOCATE (sab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sintt,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    ! here we have the S derivativs needed, lets calculate the B and C matrix

    CALL cp_fm_create(tmp_mat,&
                 matrix_struct=rtp%exp_H_new(1)%matrix%matrix_struct,&
                 name="TMP",&
                 error=error)
    CALL cp_fm_create(tmp_mat2,&
                 matrix_struct=rtp%exp_H_new(1)%matrix%matrix_struct,&
                 name="TMP2",&
                 error=error)

    CALL cp_fm_get_info(B_mat, nrow_local=nrow_local, ncol_local=ncol_local, &
          row_indices=row_indices, col_indices=col_indices,&
          nrow_global=irow_global,ncol_global=icol_global,&
          local_data=local_data,error=error)
    local_data=0._dp
    ALLOCATE(row_atom_index(SIZE(row_indices)))
    DO irow_local=1,nrow_local
       irow_global=row_indices(irow_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=irow_global .AND. irow_global <= last_sgf(iatom)) EXIT
       ENDDO
       row_atom_index(irow_local)=iatom
    ENDDO 

    ALLOCATE(col_atom_index(SIZE(col_indices)))
    DO icol_local=1,ncol_local
       icol_global=col_indices(icol_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=icol_global .AND. icol_global <= last_sgf(iatom)) EXIT
       ENDDO
       col_atom_index(icol_local)=iatom
    ENDDO 

    !calculate scalar product v(Rb)*<alpha|d/dRb beta> (B_mat), and store the first derivatives

    CALL cp_fm_set_all(B_mat,0.0_dp,0.0_dp,error)

    DO m=1,3

       CALL copy_sm_to_fm(s_derivs_sm(m)%matrix,tmp_mat2,error=error)
       CALL cp_fm_upper_to_full(tmp_mat2,tmp_mat,error)
       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             irow_global=row_indices(irow_local)
             irow_atom=row_atom_index(irow_local)
             icol_global=col_indices(icol_local)
             icol_atom=col_atom_index(icol_local)
             !signs are right, since s_derivs is <d/dR alpha|beta>
             IF (irow_global>icol_global) THEN
                scale=1.0_dp 
             ELSE
                scale=-1.0_dp 
             ENDIF
             local_data(irow_local,icol_local)=local_data(irow_local,icol_local)+&
                  scale*(particle_set(icol_atom)%v(m))*tmp_mat2%local_data(irow_local,icol_local)
             rtp%S_der(m)%matrix%local_data(irow_local,icol_local)=scale*tmp_mat2%local_data(irow_local,icol_local)
          ENDDO
       ENDDO
    END DO
    

    !calculate C matrix: v(Rb)*<d/dRa alpha| d/dRb beta>
    c_map_mat=0
    n=0
    DO j=1,3
       DO m=j,3
          n=n+1 
          c_map_mat(n,1)=j
          IF(m==j)CYCLE
          c_map_mat(n,2)=m
       END DO
    END DO

    DO m=1,3
       CALL cp_fm_set_all(C_mat(m)%matrix,0.0_dp,0.0_dp,error)
    END DO
    
    DO m=1,6
       CALL copy_sm_to_fm(s_derivs_sm(m+3)%matrix,tmp_mat2,error=error)
       CALL cp_fm_upper_to_full(tmp_mat2,tmp_mat,error)
       DO j=1,2
          IF(c_map_mat(m,j)==0)CYCLE
          CALL cp_fm_scale_and_add(1.0_dp,C_mat(c_map_mat(m,j))%matrix,-1.0_dp,tmp_mat2,error)
       END DO
    END DO
    DO m=1,3
       NULLIFY(local_data)
       CALL cp_fm_get_info(C_mat(m)%matrix,local_data=local_data,error=error)
       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             irow_global=row_indices(irow_local)
             irow_atom=row_atom_index(irow_local)
             icol_global=col_indices(icol_local)
             icol_atom=col_atom_index(icol_local)

             local_data(irow_local,icol_local)=local_data(irow_local,icol_local)*&
                  particle_set(irow_atom)%v(m)
          ENDDO
       ENDDO
    END DO

    CALL cp_fm_release(tmp_mat,error)
    CALL cp_fm_release(tmp_mat2,error)
    CALL deallocate_matrix_set(s_derivs_sm,error)
    DEALLOCATE(row_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(col_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
    
  END SUBROUTINE calc_S_derivs


! *****************************************************************************
!> \brief reads the restart file. At the moment only SCF (means only real)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************


  SUBROUTINE get_restart_wfn(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_restart_wfn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id_nr, ispin, nspin, &
                                                qs_env_id
    LOGICAL                                  :: scp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv, s_sparse
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(atomic_kind_set,mo_array,particle_set)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, matrix_s=s_sparse,&
         input=input,scp_env=scp_env,&
         id_nr=qs_env_id, dft_control=dft_control,&
         error=error)

    id_nr=0
    nspin=SIZE(mo_array)
    p_rmpv => qs_env%rho%rho_ao
    scp = dft_control%scp
    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)

    SELECT CASE(dft_control%rtp_control%initial_wfn)
    CASE(use_restart_wfn)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            scp=scp, scp_env=scp_env, error=error)
       DO i=1,SIZE(qs_env%mos)
          CALL cp_fm_to_fm(qs_env%mos(i)%mo_set%mo_coeff,qs_env%rtp%mos_old(2*i-1)%matrix,error)
          CALL cp_fm_set_all(qs_env%rtp%mos_new(2*i)%matrix,zero,zero,error)
       END DO
    CASE(use_rt_restart)
       CALL read_rt_mos_from_restart(mo_array,qs_env%rtp%mos_old,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr,dft_control%multiplicity,dft_section,&
            scp, scp_env, error)
    END SELECT
    

    DO ispin=1,nspin
       CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
            p_rmpv(ispin)%matrix,error=error)
    ENDDO
    

  END SUBROUTINE get_restart_wfn

! *****************************************************************************
!> \brief calculats the imaginary part of the density matrix
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calculate_dp_imag(dp_imag,mos_new,mo_set,error)
    TYPE(cp_fm_type), POINTER                :: dp_imag
    TYPE(cp_fm_p_type), DIMENSION(2)         :: mos_new
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dp_imag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp

    CALL timeset(routineN,handle)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp

    ELSE
       alpha=mo_set%maxocc
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
    ENDIF
    CALL get_mo_set(mo_set, nao=nao,nmo=nmo)

    CALL cp_fm_gemm("N","T",nao,nao,nmo,alpha,mos_new(1)%matrix,fm_tmp,0.0_dp,dp_imag,error)

    CALL cp_fm_gemm("N","T",nao,nao,nmo,-alpha,fm_tmp,mos_new(1)%matrix,1.0_dp,dp_imag,error)

    CALL cp_fm_release(fm_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_dp_imag

! *****************************************************************************
!> \brief calculats the three additional force contribution needed in EMD
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der 
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_c_mat_force(qs_env,dp_imag,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: dp_imag
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol_atom, icol_global, icol_local, ikind, &
      irow_atom, irow_global, irow_local, ispin, kind_atom, nao, natom, &
      ncol_local, nmo, nrow_local, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, col_atom_index, &
                                                first_sgf, kind_of, last_sgf, &
                                                row_atom_index
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, mos_new, S_der, SinvH
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: B_mat, p, S_inv, SinvB, &
                                                tmp_mat, tmp_mat1, tmp_mat2, &
                                                tmp_mat3, tmp_mo1, tmp_mo2
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(rtp,particle_set,atomic_kind_set,mos)
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,particle_set=particle_set,&
                 atomic_kind_set=atomic_kind_set,mos=mos,force=force,error=error)

    CALL get_rtp(rtp=rtp,C_mat=C_mat,error=error)

    natom = SIZE(particle_set)
    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
         error=error)

    CALL cp_fm_get_info(C_mat(1)%matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
          row_indices=row_indices, col_indices=col_indices, &
          error=error)

    CALL cp_fm_create(tmp_mat,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat1,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat2,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat3,C_mat(1)%matrix%matrix_struct,error=error)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(row_atom_index(SIZE(row_indices)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO irow_local=1,nrow_local
       irow_global=row_indices(irow_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=irow_global .AND. irow_global <= last_sgf(iatom)) EXIT
       ENDDO
       row_atom_index(irow_local)=iatom
    ENDDO 

    ALLOCATE(col_atom_index(SIZE(col_indices)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO icol_local=1,ncol_local
       icol_global=col_indices(icol_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=icol_global .AND. icol_global <= last_sgf(iatom)) EXIT
       ENDDO
       col_atom_index(icol_local)=iatom
    ENDDO 

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind,kind_of=kind_of)


    CALL cp_fm_create(dp_imag,&
         matrix_struct=qs_env%rtp%B_mat%matrix_struct,&
         name="dp_imag",&
         error=error)


    DO ispin=1,SIZE(mos)

       CALL calculate_dp_imag(dp_imag,qs_env%rtp%mos_new(2*ispin-1:2*ispin),mos(ispin)%mo_set,error)

       alpha=mos(ispin)%mo_set%maxocc

       ! calculate B_mat force 

       CALL get_rtp(rtp=rtp,B_mat=B_mat,S_der=S_der,S_inv=S_inv,SinvH=SinvH,SinvB=SinvB,mos_new=mos_new,error=error)
       CALL cp_fm_get_info(mos_new(2*ispin-1)%matrix,nrow_global=nao,ncol_global=nmo,error=error)

       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,ncol_global=nao,&
            para_env=mos_new(2*ispin-1)%matrix%matrix_struct%para_env, &
            context=mos_new(2*ispin-1)%matrix%matrix_struct%context,error=error)

       CALL cp_fm_create(tmp_mo1,tmp_fm_struct,error=error)
       CALL cp_fm_create(tmp_mo2,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(2*ispin)%matrix,SinvB,0.0_dp,tmp_mo1,error=error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(2*ispin-1)%matrix,SinvB,0.0_dp,tmp_mo2,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,1.0_dp,mos_new(2*ispin-1)%matrix,tmp_mo1,0.0_dp,tmp_mat1,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,-1.0_dp,mos_new(2*ispin)%matrix,tmp_mo2,1.0_dp,tmp_mat1,error=error)
       
       CALL cp_fm_transpose(tmp_mat1,tmp_mat,error)
       CALL cp_fm_scale_and_add(0.0_dp,tmp_mat2,1.0_dp,tmp_mat,error)
       CALL cp_fm_scale_and_add(0.0_dp,tmp_mat3,1.0_dp,tmp_mat1,error)

       !calculate the W matrix analogue for EMD

       CALL cp_fm_create(p,C_mat(1)%matrix%matrix_struct,error=error)
       CALL copy_sm_to_fm(qs_env%rho%rho_ao(ispin)%matrix,p,error=error)
       CALL cp_fm_upper_to_full(p,tmp_mat,error)

       
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(2*ispin-1)%matrix,SinvH(ispin)%matrix,0.0_dp,tmp_mo1,error=error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(2*ispin)%matrix,SinvH(ispin)%matrix,0.0_dp,tmp_mo2,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,1.0_dp,mos_new(2*ispin-1)%matrix,tmp_mo1,0.0_dp,tmp_mat1,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,1.0_dp,mos_new(2*ispin)%matrix,tmp_mo2,1.0_dp,tmp_mat1,error=error)
        
!      CALL cp_fm_gemm("T","T",nao,nao,nao,1.0_dp,p,SinvH(ispin)%matrix,0.0_dp,tmp_mat1,error=error)
       CALL cp_fm_transpose(tmp_mat1,tmp_mat,error)

       CALL cp_fm_scale_and_add(1.0_dp,tmp_mat2,1.0_dp,tmp_mat,error)
       CALL cp_fm_scale_and_add(1.0_dp,tmp_mat3,1.0_dp,tmp_mat1,error)

       !C_mat forces

       CALL cp_fm_transpose(dp_imag,tmp_mat,error)
       CALL cp_fm_get_info(tmp_mat,local_data=local_data,error=error)

       !This loop evaluates all forces, and puts them in the right place
       DO i=1,3
          CALL cp_fm_transpose(C_mat(i)%matrix,tmp_mat1,error)

          DO icol_local=1,ncol_local
             DO irow_local=1,nrow_local

                icol_atom=col_atom_index(icol_local)
                irow_atom=row_atom_index(irow_local)
                ikind=kind_of(icol_atom)
                kind_atom=atom_of_kind(icol_atom)

                force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
                     tmp_mat2%local_data(irow_local,icol_local)* S_der(i)%matrix%local_data(irow_local,icol_local)-&
                     local_data(irow_local,icol_local)* C_mat(i)%matrix%local_data(irow_local,icol_local)

                ikind=kind_of(irow_atom)
                kind_atom=atom_of_kind(irow_atom)

                force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)-&
                     tmp_mat3%local_data(irow_local,icol_local)* S_der(i)%matrix%local_data(irow_local,icol_local)+ &
                     local_data(irow_local,icol_local)* tmp_mat1%local_data(irow_local,icol_local)             
             END DO
          END DO
       END DO

       ! recall QS forces, at this point have the other sign.
       CALL cp_fm_release(p,error)
       CALL cp_fm_release(tmp_mo1,error)
       CALL cp_fm_release(tmp_mo2,error)
       
    END DO
    ! recall QS forces, at this point have the other sign.
    DO i=1,SIZE(force)
       force(i)%ehrenfest(:,:)=- force(i)%ehrenfest(:,:)    
    END DO

    CALL cp_fm_release(tmp_mat,error)
    CALL cp_fm_release(tmp_mat1,error)
    CALL cp_fm_release(tmp_mat2,error)
    CALL cp_fm_release(tmp_mat3,error)
    CALL cp_fm_release(dp_imag,error)

    DEALLOCATE(row_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(col_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE calc_c_mat_force

END MODULE rt_propagation_utils
