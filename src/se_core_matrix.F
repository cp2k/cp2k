!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Hamiltonian integral matrix <a|H|b> for
!>      semi-empirical methods
!> \author JGH
! *****************************************************************************
MODULE se_core_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_from_sm,&
                                             dbcsr_obj_p_type
  USE cp_dbcsr_output,                 ONLY: dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_init
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_distribute,&
                                             dbcsr_get_block_diag,&
                                             dbcsr_get_block_p,&
                                             dbcsr_replicate_all,&
                                             dbcsr_set,&
                                             dbcsr_sum_replicated
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pdg, &
       do_method_pm3, do_method_pm6, do_method_pnnl, do_method_rm1, &
       use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_core_matrix'

  PUBLIC :: build_se_core_matrix

CONTAINS

! *****************************************************************************
  SUBROUTINE build_se_core_matrix(qs_env,para_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_se_core_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, i, iab, iatom, icol, icor, ikind, ilist, &
      inode, irow, istat, itype, iw, j, jatom, jkind, natom, natorb_a, &
      natorb_b, nkind, nlist, nnode
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, found, tr, &
                                                use_virial
    REAL(KIND=dp)                            :: delta, dr, econst, eheat, &
                                                eisol, enuclear, kh, udd, &
                                                uff, upp, uss
    REAL(KIND=dp), DIMENSION(16)             :: ha, hb, ua
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta_a, beta_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsmat, h_block, h_blocka, &
                                                pabmat, pamat, s_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_obj), POINTER                 :: diagmat_h, diagmat_p
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p_sm
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

     NULLIFY ( matrix_p_sm,rho,force,atomic_kind_set,sab_orb,&
               diagmat_h,diagmat_p,particle_set, matrix_p )

     ! calculate overlap matrix
     IF(calculate_forces) THEN
        CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, error=error)
        CALL build_overlap_matrix(qs_env,para_env,nderivative=1, matrix_s=matrix_s,&
                              matrix_name="OVERLAP",&
                              basis_set_id_a=use_orb_basis_set,&
                              basis_set_id_b=use_orb_basis_set, &
                              neighbor_list_sab=qs_env%sab_orb,&
                              error=error)
        CALL set_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
        CALL get_qs_env(qs_env=qs_env, virial=virial, error=error)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
     ELSE
        CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, error=error)
        CALL build_overlap_matrix(qs_env,para_env,matrix_s=matrix_s,&
                              matrix_name="OVERLAP",&
                              basis_set_id_a=use_orb_basis_set,&
                              basis_set_id_b=use_orb_basis_set, &
                              neighbor_list_sab=qs_env%sab_orb,&
                              error=error)
        CALL set_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
        use_virial = .FALSE.
     END IF

     CALL get_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h,&
                     matrix_s=matrix_s,&
                     atomic_kind_set=atomic_kind_set,&
                     dft_control=dft_control,error=error)

     IF(calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       particle_set=particle_set,&
                       rho=rho,&
                       force=force,error=error)
       matrix_p_sm => rho%rho_ao

       NULLIFY(matrix_p)!sm->dbcsr
       CALL dbcsr_allocate_matrix_set(matrix_p,SIZE(matrix_p_sm),error)!sm->dbcsr
       DO i=1,SIZE(matrix_p)!sm->dbcsr
          ALLOCATE(matrix_p(i)%matrix)
          CALL dbcsr_from_sm( matrix_p(i)%matrix, matrix_p_sm(i)%matrix, error=error )!sm->dbcsr
       ENDDO!sm->dbcsr

       IF (SIZE(matrix_p) == 2) THEN
          CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, alpha_scalar=1.0_dp, beta_scalar=1.0_dp, error=error)
       END IF
       natom = SIZE (particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       delta = dft_control%qs_control%se_control%delta
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)
       ALLOCATE(diagmat_p)!sm->dbcsr
       CALL dbcsr_init (diagmat_p)
       CALL dbcsr_get_block_diag(matrix_p(1)%matrix, diagmat_p, error=error)
       CALL dbcsr_replicate_all(diagmat_p)
     END IF

     ! Allocate the core Hamiltonian matrix
     CALL dbcsr_allocate_matrix_set(matrix_h,1,error=error)
     ALLOCATE(matrix_h(1)%matrix)
     CALL dbcsr_init(matrix_h(1)%matrix)
     CALL dbcsr_copy(matrix_h(1)%matrix,matrix_s(1)%matrix,"CORE HAMILTONIAN MATRIX",error=error)
     CALL dbcsr_set(matrix_h(1)%matrix,0.0_dp,error=error)
     !CALL dbcsr_replicate_matrix_struct(source=matrix_s(1)%matrix,&
     !               TARGET=matrix_h(1)%matrix,&
     !               target_name="CORE HAMILTONIAN MATRIX",error=error)

     ! Allocate a diagonal block matrix
     ALLOCATE(diagmat_h)!sm->dbcsr
     CALL dbcsr_init(diagmat_h)
     CALL dbcsr_get_block_diag(matrix_s(1)%matrix, diagmat_h, error=error)
     CALL dbcsr_set(diagmat_h, 0.0_dp, error=error)
     CALL dbcsr_replicate_all(diagmat_h)

     CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,error=error)

     ! kh might be set in qs_control
     itype = get_se_type(dft_control%qs_control%method_id)
     kh = 0.5_dp

     nkind = SIZE(atomic_kind_set)

     enuclear = 0.0_dp
     econst   = 0.0_dp
     DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,&
                             se_parameter=se_kind_a)
        CALL get_se_param(se_kind_a,&
                          defined=defined,&
                          natorb=natorb_a,&
                          beta=beta_a,&
                          uss=uss,upp=upp,udd=udd,uff=uff,&
                          eisol=eisol,eheat=eheat)

        IF (.NOT.defined .OR. natorb_a < 1) CYCLE

        econst=econst-(eisol-eheat)*REAL(natom,dp)

        ha( 1   ) = beta_a(0)
        ha( 2: 4) = beta_a(1)
        ha( 5: 9) = beta_a(2)
        ha(10:16) = beta_a(3)

        ua( 1   ) = uss
        ua( 2: 4) = upp
        ua( 5: 9) = udd
        ua(10:16) = uff

        DO jkind=1,nkind
           atomic_kind => atomic_kind_set(jkind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                se_parameter=se_kind_b)
           CALL get_se_param(se_kind_b,&
                             defined=defined,&
                             natorb=natorb_b,&
                             beta=beta_b)
           IF (.NOT.defined .OR. natorb_b < 1) CYCLE

           hb( 1   ) = beta_b(0)
           hb( 2: 4) = beta_b(1)
           hb( 5: 9) = beta_b(2)
           hb(10:16) = beta_b(3)

           iab = ikind + nkind*(jkind - 1)
           IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(iab)%neighbor_list_set,&
                                      nlist=nlist)
           NULLIFY(sab_orb_neighbor_list)
           DO ilist=1,nlist
             IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
               sab_orb_neighbor_list => &
                          first_list(sab_orb(iab)%neighbor_list_set)
             ELSE
               sab_orb_neighbor_list => next(sab_orb_neighbor_list)
             END IF
             sab_orb_neighbor_list_local => sab_orb_neighbor_list
             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                    atom=iatom,&
                                    nnode=nnode)

             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

             SELECT CASE (dft_control%qs_control%method_id)
               CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,&
                     do_method_pm3,do_method_pm6,do_method_mndod, do_method_pnnl)
                 NULLIFY(h_blocka)
                 CALL dbcsr_get_block_p(diagmat_h,&
                      iatom,iatom,h_blocka,tr,found)
                 CPPostcondition(ASSOCIATED(h_blocka),cp_failure_level,routineP,error,failure)
                 IF(calculate_forces) THEN
                   CALL dbcsr_get_block_p(diagmat_p,&
                        iatom,iatom,pamat,tr,found)
                   CPPostcondition(ASSOCIATED(pamat),cp_failure_level,routineP,error,failure)
                 END IF
             END SELECT
             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                      neighbor=jatom,r=rij)
               dr = SUM(rij(:)**2)
               IF (iatom == jatom .AND. dr < rij_threshold) THEN

                 SELECT CASE (dft_control%qs_control%method_id)
                    CASE DEFAULT
                      CALL stop_program(routineN,moduleN,__LINE__,&
                                      "Method not available",para_env)
                    CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,&
                          do_method_pm3,do_method_pm6,do_method_mndod, do_method_pnnl)
                      DO i=1,SIZE(h_blocka,1)
                        h_blocka(i,i)=h_blocka(i,i)+ua(i)
                      END DO
                 END SELECT

               ELSE
                 IF (iatom <= jatom) THEN
                   irow = iatom
                   icol = jatom
                 ELSE
                   irow = jatom
                   icol = iatom
                 END IF
                 NULLIFY (h_block)
                 CALL dbcsr_get_block_p(matrix_h(1)%matrix,&
                      irow,icol,h_block,tr,found)
                 CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)
                 ! two-centre one-electron term
                 NULLIFY(s_block)
                 CALL dbcsr_get_block_p(matrix_s(1)%matrix,&
                      irow,icol,s_block,tr,found)
                 CPPostcondition(ASSOCIATED(s_block),cp_failure_level,routineP,error,failure)
                 IF ( irow == iatom ) THEN
                   DO i=1,SIZE(h_block,1)
                     DO j=1,SIZE(h_block,2)
                       h_block(i,j)=h_block(i,j)+kh*(ha(i)+hb(j))*s_block(i,j)
                     END DO
                   END DO
                 ELSE
                   DO i=1,SIZE(h_block,1)
                     DO j=1,SIZE(h_block,2)
                       h_block(i,j)=h_block(i,j)+kh*(ha(j)+hb(i))*s_block(i,j)
                     END DO
                   END DO
                 END IF
                 IF(calculate_forces) THEN
                   atom_a = atom_of_kind(iatom)
                   atom_b = atom_of_kind(jatom)
                   CALL dbcsr_get_block_p(matrix_p(1)%matrix,&
                        irow,icol,pabmat,tr,found)
                   CPPostcondition(ASSOCIATED(pabmat),cp_failure_level,routineP,error,failure)
                   DO icor=1,3
                     force_ab(icor) = 0._dp
                     CALL dbcsr_get_block_p(matrix_s(icor+1)%matrix,&
                          irow,icol,dsmat,tr,found)
                     CPPostcondition(ASSOCIATED(dsmat),cp_failure_level,routineP,error,failure)
                     dsmat=2._dp*kh*dsmat*pabmat
                     IF ( irow == iatom ) THEN
                        DO i=1,SIZE(h_block,1)
                          DO j=1,SIZE(h_block,2)
                            force_ab(icor)=force_ab(icor)+(ha(i)+hb(j))*dsmat(i,j)
                          END DO
                        END DO
                      ELSE
                        DO i=1,SIZE(h_block,1)
                          DO j=1,SIZE(h_block,2)
                            force_ab(icor)=force_ab(icor)+(ha(j)+hb(i))*dsmat(i,j)
                          END DO
                        END DO
                      END IF
                   END DO
                 END IF

               END IF

               IF(calculate_forces .AND. (iatom/=jatom .OR. dr > rij_threshold)) THEN
                 IF ( irow == iatom ) force_ab = -force_ab
                 force(ikind)%all_potential(:,atom_a) =&
                     force(ikind)%all_potential(:,atom_a) - force_ab(:)
                 force(jkind)%all_potential(:,atom_b) =&
                     force(jkind)%all_potential(:,atom_b) + force_ab(:)
                 IF (use_virial) THEN
                   CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
                 END IF
               END IF

               sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     CALL dbcsr_sum_replicated(diagmat_h)
     CALL dbcsr_distribute(diagmat_h)
     CALL dbcsr_add(matrix_h(1)%matrix, diagmat_h,error=error)
     CALL set_qs_env(qs_env=qs_env,matrix_h=matrix_h,error=error)

     qs_env%energy%core_self    = econst

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL dbcsr_write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,para_env,&
             scale=evolt,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
     END IF

     IF(calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
       END IF
       DEALLOCATE(atom_of_kind)
       CALL dbcsr_deallocate_matrix(diagmat_p,error=error)
     END IF

     CALL dbcsr_deallocate_matrix(diagmat_h,error=error)
     CALL dbcsr_deallocate_matrix_set(matrix_p,error=error)

  END SUBROUTINE build_se_core_matrix

END MODULE se_core_matrix

