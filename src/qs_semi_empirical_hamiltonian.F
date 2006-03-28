!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004-2006  CP2K developers group
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_semi_empirical_hamiltonian [1.1] *
!!
!!   FUNCTION
!!     Calculation of the Hamiltonian integral matrix <a|H|b> for 
!!     semi-empirical methods
!!
!!   AUTHOR
!!     JGH
!!
!!***
!******************************************************************************

MODULE qs_semi_empirical_hamiltonian

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE nddo_methods,                    ONLY: build_nddo_fock_matrix
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_overlap,                      ONLY: build_overlap_matrix,&
                                             write_sparse_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_param,&
                                             get_se_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             cp_sm_scale_and_add,&
                                             create_dbm_from_sm,&
                                             deallocate_dbm_from_sm,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_semi_empirical_hamiltonian'

  PUBLIC :: build_se_core_matrix, build_se_fock_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_se_core_matrix(qs_env,globenv,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_se_core_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, i, iab, iatom, icol, icor, ikind, ilist, &
      inode, irow, istat, itype, iw, j, jatom, jkind, natom, natorb_a, &
      natorb_b, nkind, nlist, nnode
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure
    REAL(KIND=dp)                            :: delta, econst, eheat, eisol, &
                                                enuclear, kh, range, rc, udd, &
                                                uff, upp, uss, dr
    REAL(KIND=dp), DIMENSION(16)             :: ha, hb, ua
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta_a, beta_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsmat, h_block, h_blocka, &
                                                pabmat, pamat, s_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
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
      POINTER                                :: matrix_h, matrix_p, matrix_s
    TYPE(real_matrix_type), POINTER          :: diagmat_h, diagmat_p
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

     NULLIFY ( matrix_h,matrix_s,matrix_p,rho,force,atomic_kind_set,sab_orb,&
               diagmat_h,diagmat_p,particle_set )

     ! calculate overlap matrix
     IF(calculate_forces) THEN
       CALL build_overlap_matrix(qs_env=qs_env,globenv=globenv,nderivative=1,error=error)
     ELSE
       CALL build_overlap_matrix(qs_env=qs_env,globenv=globenv,error=error)
     END IF

     CALL get_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h,&
                     matrix_s=matrix_s,&
                     atomic_kind_set=atomic_kind_set,&
                     dft_control=dft_control,&
                     para_env=para_env)

     IF(calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       particle_set=particle_set,&
                       rho=rho,&
                       force=force)
       matrix_p => rho%rho_ao
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,&
               matrix_p(2)%matrix,1.0_dp)
       END IF
       natom = SIZE (particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       delta = dft_control%qs_control%se_control%delta
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)
       CALL create_dbm_from_sm(diagmat_p,matrix_p(1)%matrix,&
                               "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env)
     END IF

     ! Allocate the core Hamiltonian matrix
     CALL allocate_matrix_set(matrix_h,1)
     CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                    TARGET=matrix_h(1)%matrix,&
                    target_name="CORE HAMILTONIAN MATRIX")

     ! Allocate a diagonal block matrix
     CALL create_dbm_from_sm(diagmat_h,matrix_s(1)%matrix,&
                             "DIAGONAL BLOCK H MATRIX",.FALSE.)

     CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb)

     ! kh might be set in qs_control
     kh=1._dp
     IF(dft_control%qs_control%method == "EHT") kh=1.75_dp
     itype = get_se_type(dft_control%qs_control%method)
     kh = 0.5_dp * kh

     nkind = SIZE(atomic_kind_set)

     enuclear = 0._dp
     econst = 0._dp

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

        ha(1) = beta_a(0)
        ha(2:4) = beta_a(1)
        ha(5:9) = beta_a(2)
        ha(10:16) = beta_a(3)
        ua(1) = uss
        ua(2:4) = upp
        ua(5:9) = udd
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

           hb(1) = beta_b(0)
           hb(2:4) = beta_b(1)
           hb(5:9) = beta_b(2)
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

             SELECT CASE (dft_control%qs_control%method)
               CASE ("EHT")
               CASE ("MNDO","AM1","PM3","PDG")
                 NULLIFY(h_blocka)
                 CALL get_block_node(matrix=diagmat_h,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     BLOCK=h_blocka)
                 CPPostcondition(ASSOCIATED(h_blocka),cp_failure_level,routineP,error,failure)
                 IF(calculate_forces) THEN
                   CALL get_block_node(matrix=diagmat_p,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pamat)
                   CPPostcondition(ASSOCIATED(pamat),cp_failure_level,routineP,error,failure)
                 END IF
             END SELECT

             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                      neighbor=jatom,r=rij)

               dr = SUM(rij(:)**2)
               IF (iatom == jatom .AND. dr < 0.001_dp) THEN

                 SELECT CASE (dft_control%qs_control%method)
                    CASE DEFAULT
                      CALL stop_program(routineN,moduleN,__LINE__,&
                                      "Method not available",globenv)
                    CASE ("EHT")
                      DO i=1,SIZE(h_blocka,1)
                        h_blocka(i,i)=h_blocka(i,i)+ha(i)
                      END DO
                    CASE ("AM1","PM3","MNDO","PDG")
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
                 CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     BLOCK=h_block)
                 CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

                 ! two-centre one-electron term
                 NULLIFY(s_block)
                 CALL get_block_node(matrix=matrix_s(1)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     BLOCK=s_block)
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
                   CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=pabmat)
                   CPPostcondition(ASSOCIATED(pabmat),cp_failure_level,routineP,error,failure)
                   DO icor=1,3
                     force_ab(icor) = 0._dp
                     CALL get_block_node(matrix=matrix_s(icor+1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         BLOCK=dsmat)
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

               IF(calculate_forces .AND. (iatom/=jatom .OR. dr > 0.001_dp)) THEN
                 IF ( irow == iatom ) force_ab = -force_ab
                 force(ikind)%all_potential(:,atom_a) =&
                     force(ikind)%all_potential(:,atom_a) - force_ab(:)
                 force(jkind)%all_potential(:,atom_b) =&
                     force(jkind)%all_potential(:,atom_b) + force_ab(:)
               END IF

               sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     CALL deallocate_dbm_from_sm(diagmat_h,matrix_h(1)%matrix,.TRUE.,para_env)

     CALL set_qs_env(qs_env=qs_env,matrix_h=matrix_h)

     qs_env%energy%core_self    = econst

     IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,globenv,&
             scale=evolt,output_unit=iw)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", error=error)
     END IF

     IF(calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix,1.0_dp,&
               matrix_p(2)%matrix,-1.0_dp)
       END IF
       DEALLOCATE(atom_of_kind)
       CALL deallocate_dbm_from_sm(diagmat_p,matrix_p(1)%matrix,.FALSE.)
     END IF

  END SUBROUTINE build_se_core_matrix

! *****************************************************************************

  SUBROUTINE build_se_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_fock_matrix', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control

!   ----------------------------------------------------------------------------

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control)

    SELECT CASE (dft_control%qs_control%method)
      CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE ("AM1","PM3","MNDO","PDG")
        CALL build_nddo_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)
    END SELECT

  END SUBROUTINE build_se_fock_matrix

! *****************************************************************************

END MODULE qs_semi_empirical_hamiltonian

! *****************************************************************************
