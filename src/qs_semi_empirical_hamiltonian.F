!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE simulation_cell,                 ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_overlap,                      ONLY: build_overlap_matrix,&
                                             write_sparse_matrix
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             get_qs_env,&
                                             set_qs_env
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_integrals,        ONLY: rotint,rotnuc,drotint,drotnuc
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_param
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             cp_sm_scale_and_add,&
                                             transfer_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "qs_semi_empirical_hamiltonian"

  PUBLIC :: build_se_core_matrix, build_se_fock_matrix

  INTERFACE fock1
    MODULE PROCEDURE rfock1, ufock1
  END INTERFACE

  INTERFACE fock2
    MODULE PROCEDURE rfock2, ufock2
  END INTERFACE

  INTERFACE dfock2
    MODULE PROCEDURE rdfock2, udfock2
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_se_core_matrix(qs_env,globenv,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "build_se_core_matrix",&
      routineP = moduleN//':'//routineN

    TYPE(real_matrix_p_type), DIMENSION(:), &
                               POINTER       :: matrix_h,matrix_s,matrix_p
    TYPE(atomic_kind_type), DIMENSION(:), &
                               POINTER       :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(neighbor_list_set_p_type), &
                       DIMENSION(:), POINTER :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list,&
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_force_type), DIMENSION(:), &
                          POINTER            :: force
    TYPE(particle_type), DIMENSION(:), &
                          POINTER            :: particle_set

    INTEGER                                  :: i,j,k,ikind,jkind,nkind
    INTEGER                                  :: iab,ilist,nlist,iatom,jatom,&
                                                inode,nnode,istat,&
                                                irow,icol,natorb_a,natorb_b
    LOGICAL                                  :: defined,failure
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_block,s_block
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_blocka,h_blockb
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: pamat,pbmat,pabmat,dsmat
    REAL(KIND = dp), DIMENSION(:), POINTER   :: beta_a,beta_b
    REAL(KIND = dp), DIMENSION(16)           :: ha,hb,ua
    INTEGER                                  :: za,zb,itype,natom,icor,&
                                                atom_a,atom_b
    REAL(KIND = dp)                          :: kh,enuc,uss,upp,udd,uff
    REAL(KIND = dp), DIMENSION(3)            :: rij,denuc,force_ab
    REAL(KIND = dp), DIMENSION(10)           :: e1b,e2a,pvec
    REAL(KIND = dp), DIMENSION(10,3)         :: de1b,de2a
    REAL(KIND = dp)                          :: enuclear,eisol,eheat,econst,&
                                                delta
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind

!   ---------------------------------------------------------------------------

     NULLIFY ( matrix_h,matrix_s,matrix_p,rho,force,atomic_kind_set,sab_orb,&
               particle_set )

     ! calculate overlap matrix
     IF(calculate_forces) THEN
       CALL build_overlap_matrix(qs_env=qs_env,globenv=globenv,nderivative=1)
     ELSE
       CALL build_overlap_matrix(qs_env=qs_env,globenv=globenv)
     END IF

     CALL get_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h,&
                     matrix_s=matrix_s,&
                     atomic_kind_set=atomic_kind_set,&
                     dft_control=dft_control)

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
     END IF

     ! Allocate the core Hamiltonian matrix
     CALL allocate_matrix_set(matrix_h,1)
     CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                    target=matrix_h(1)%matrix,&
                    target_name="CORE HAMILTONIAN MATRIX")

     CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb)

     ! kh might be set in qs_control
     SELECT CASE (dft_control%qs_control%method)
       CASE DEFAULT
         kh=1._dp
         itype=0
       CASE ("EHT")
         kh=1.75_dp
         itype=0
       CASE ("AM1")
         kh=1._dp
         itype=2
       CASE ("PM3")
         kh=1._dp
         itype=3
       CASE ("MNDO")
         kh=1._dp
         itype=1
     END SELECT
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
                          zeff=za,&
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
                             zeff=zb,&
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
               CASE ("MNDO","AM1","PM3")
                 NULLIFY(h_blocka)
                 CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=h_blocka)
                 CPPostcondition(ASSOCIATED(h_blocka),cp_failure_level,routineP,error,failure)
                 IF(calculate_forces) THEN
                   CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       block=pamat)
                   CPPostcondition(ASSOCIATED(pamat),cp_failure_level,routineP,error,failure)
                 END IF
             END SELECT

             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                      neighbor=jatom,r=rij)

               IF (iatom == jatom) THEN

                 SELECT CASE (dft_control%qs_control%method)
                    CASE DEFAULT
                      CALL stop_program(routineN,moduleN,__LINE__,&
                                      "Method not available",globenv)
                    CASE ("EHT")
                      DO i=1,SIZE(h_blocka,1)
                        h_blocka(i,i)=h_blocka(i,i)+ha(i)
                      END DO
                    CASE ("AM1","PM3","MNDO")
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
                                     block=h_block)
                 CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

                 ! two-centre one-electron term
                 NULLIFY(s_block)
                 CALL get_block_node(matrix=matrix_s(1)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=s_block)
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
                                       block=pabmat)
                   CPPostcondition(ASSOCIATED(pabmat),cp_failure_level,routineP,error,failure)
                   DO icor=1,3
                     force_ab(icor) = 0._dp
                     CALL get_block_node(matrix=matrix_s(icor+1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=dsmat)
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

                 SELECT CASE (dft_control%qs_control%method)
                   CASE ("EHT")
                   CASE ("MNDO","AM1","PM3")
                     IF ( irow == iatom ) THEN
                       CALL rotnuc (se_kind_a,se_kind_b,rij,e1b,e2a,enuc,itype)
                     ELSE
                       CALL rotnuc (se_kind_b,se_kind_a,-rij,e2a,e1b,enuc,itype)
                     END IF
                     enuclear = enuclear + enuc
                     ! one-centre one-electron terms
                     NULLIFY(h_blockb)
                     CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=h_blockb)
                     CPPostcondition(ASSOCIATED(h_blockb),cp_failure_level,routineP,error,failure)
                     ! contribution to iatom
                     IF ( SIZE(h_blocka,1) == 1 ) THEN
                       h_blocka(1,1)=h_blocka(1,1)+e1b(1)
                     ELSE
                       CALL expand ( h_blocka, e1b )
                     END IF
                     ! contribution to jatom
                     IF ( SIZE(h_blockb,1) == 1 ) THEN
                       h_blockb(1,1)=h_blockb(1,1)+e2a(1)
                     ELSE
                       CALL expand ( h_blockb, e2a )
                     END IF
                     IF(calculate_forces) THEN
                       CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                           block_row=jatom,&
                                           block_col=jatom,&
                                           block=pbmat)
                       CPPostcondition(ASSOCIATED(pbmat),cp_failure_level,routineP,error,failure)
                       IF ( irow == iatom ) THEN
                         CALL drotnuc (se_kind_a,se_kind_b,rij,de1b,de2a,&
                                       denuc,itype,delta)
                       ELSE
                         CALL drotnuc (se_kind_b,se_kind_a,-rij,de2a,de1b,&
                                       denuc,itype,delta)
                       END IF
                       force_ab(1:3)=force_ab(1:3)-denuc(1:3)
                       IF ( SIZE(h_blocka,1) == 1 ) THEN
                         force_ab(1:3)=force_ab(1:3)-de1b(1,1:3)*pamat(1,1)
                       ELSE
                         CALL contract ( pamat, pvec )
                         DO icor=1,3
                           force_ab(icor)=force_ab(icor)-&
                                          SUM(de1b(1:10,icor)*pvec(1:10))
                         END DO
                       END IF
                       IF ( SIZE(h_blockb,1) == 1 ) THEN
                         force_ab(1:3)=force_ab(1:3)-de2a(1,1:3)*pbmat(1,1)
                       ELSE
                         CALL contract ( pbmat, pvec )
                         DO icor=1,3
                           force_ab(icor)=force_ab(icor)-&
                                          SUM(de2a(1:10,icor)*pvec(1:10))
                         END DO
                       END IF
                     END IF

                 END SELECT

               END IF

               IF(calculate_forces .AND. iatom/=jatom) THEN
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

     CALL set_qs_env(qs_env=qs_env,matrix_h=matrix_h)

     qs_env%energy%core_overlap = enuclear
     qs_env%energy%core_self    = econst

     IF (globenv%print%core_hamiltonian_matrix) THEN
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,globenv,&
             scale=evolt)
     END IF

     !deallocate overlap matrix for orthogonal methods
     IF ( dft_control%qs_control%se_control%orthogonal_basis ) THEN
       !
     END IF

     IF(calculate_forces) THEN
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_sm_scale_and_add(matrix_p(1)%matrix,1.0_dp,&
               matrix_p(2)%matrix,-1.0_dp)
       END IF
     END IF

  END SUBROUTINE build_se_core_matrix

! *****************************************************************************

  SUBROUTINE build_se_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
        POINTER                              :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), INTENT(inout), &
      TARGET                                 :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_fock', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure,ionode,defined
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: handle,nspins,output_unit,&
                                                ispin,ikind,nkind,jkind
    INTEGER                                  :: za,natorb_a,zb,natorb_b,&
                                                ilist,nlist,iatom,jatom,&
                                                inode,nnode,iab,irow,icol
    REAL ( dp ), DIMENSION (3)               :: rij
    REAL ( dp )                              :: ecoul
    REAL ( dp )                              :: gss,gsp,gpp,gp2,hsp
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h,matrix_p
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(atomic_kind_type), DIMENSION(:), &
                               POINTER       :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(neighbor_list_set_p_type), &
                       DIMENSION(:), POINTER :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list,&
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(qs_force_type), DIMENSION(:), &
                          POINTER            :: force
    TYPE(particle_type), DIMENSION(:), &
                          POINTER            :: particle_set

    REAL(KIND = dp), DIMENSION(:,:), POINTER :: p_block_a,p_block_b,&
                                                pa_block_a,pa_block_b,&
                                                pb_block_a,pb_block_b,&
                                                ks_block_a,ks_block_b,&
                                                ksa_block_a,ksa_block_b,&
                                                ksb_block_a,ksb_block_b
    INTEGER                                  :: atom_a,atom_b,natom,istat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    REAL(KIND = dp)                          :: delta
    REAL(KIND = dp), DIMENSION(3)            :: force_ab

!   ----------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,force,particle_set,logger)

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         cell=cell,&
         para_env=para_env)

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       DO ispin=1,nspins
          ! copy the core matrix into the fock matrix
          CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix)
       END DO

       matrix_p => rho%rho_ao

       CALL get_qs_env(qs_env=qs_env,&
                       sab_orb=sab_orb,&
                       atomic_kind_set=atomic_kind_set)

       nkind = SIZE(atomic_kind_set)

       IF(calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env,&
                         particle_set=particle_set,&
                         force=force)
         natom = SIZE (particle_set)
         ALLOCATE (atom_of_kind(natom),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         delta = dft_control%qs_control%se_control%delta
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                  atom_of_kind=atom_of_kind)
       END IF

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,&
                            zeff=za,&
                            defined=defined,&
                            natorb=natorb_a,&
                            gss=gss,gsp=gsp,gpp=gpp,gp2=gp2,hsp=hsp)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,&
                               zeff=zb,&
                               defined=defined,&
                               natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                                        sab_orb(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_orb_neighbor_list)

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

               CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
               IF ( nspins >= 1 ) THEN
                 CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=pa_block_a)
                 CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=ksa_block_a)
                 CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
               END IF
               IF ( nspins >= 2 ) THEN
                 CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=pa_block_b)
                 CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=ksa_block_b)
                 CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
               END IF

               DO inode=1,nnode
                 CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                        neighbor=jatom,r=rij)

                 IF ( iatom == jatom ) THEN
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins == 1 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,pa_block_a,ksa_block_a,error)
                   ELSE IF ( nspins == 2 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,pa_block_a,ksa_block_a,&
                                pa_block_b,ksa_block_b,error)
                   END IF

                 ELSE

                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins >= 1 ) THEN
                     CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=pb_block_a)
                     CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=ksb_block_a)
                     CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                   END IF
                   IF ( nspins >= 2 ) THEN
                     CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=pb_block_b)
                     CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=ksb_block_b)
                     CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                   END IF

                   IF (iatom <= jatom) THEN
                     irow = iatom
                     icol = jatom
                   ELSE
                     irow = jatom
                     icol = iatom
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins >= 1 ) THEN
                     CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=ks_block_a)
                     CPPostcondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=p_block_a)
                     CPPostcondition(ASSOCIATED(p_block_a),cp_failure_level,routineP,error,failure)
                   END IF
                   IF ( nspins >= 2 ) THEN
                     CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=ks_block_b)
                     CPPostcondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=p_block_b)
                     CPPostcondition(ASSOCIATED(p_block_b),cp_failure_level,routineP,error,failure)
                   END IF

                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins == 1 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,p_block_a,&
                                  ksa_block_a,ksb_block_a,ks_block_a,error)
                     ELSE
                       CALL fock2(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,p_block_a,&
                                  ksb_block_a,ksa_block_a,ks_block_a,error)
                     ENDIF
                   ELSE IF ( nspins == 2 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,p_block_a,&
                                  ksa_block_a,ksb_block_a,ks_block_a,&
                                  pa_block_b,pb_block_b,p_block_b,&
                                  ksa_block_b,ksb_block_b,ks_block_b,error)
                     ELSE
                       CALL fock2(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,p_block_a,&
                                  ksb_block_a,ksa_block_a,ks_block_a,&
                                  pb_block_b,pa_block_b,p_block_b,&
                                  ksb_block_b,ksa_block_b,ks_block_b,error)
                     ENDIF
                   END IF

                   IF(calculate_forces) THEN
                     force_ab(1:3) = 0._dp
                     CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                     IF ( nspins == 1 ) THEN
                       IF ( irow == iatom ) THEN
                         CALL dfock2(se_kind_a,se_kind_b,rij,pa_block_a,&
                              pb_block_a,p_block_a,force_ab,delta,error)
                       ELSE
                         CALL dfock2(se_kind_b,se_kind_a,-rij,pb_block_a,&
                              pa_block_a,p_block_a,force_ab,delta,error)
                       ENDIF
                     ELSE IF ( nspins == 2 ) THEN
                       IF ( irow == iatom ) THEN
                         CALL dfock2(se_kind_a,se_kind_b,rij,pa_block_a,&
                              pb_block_a,p_block_a,pa_block_b,pb_block_b,&
                              p_block_b,force_ab,delta,error)
                       ELSE
                         CALL dfock2(se_kind_b,se_kind_a,-rij,pb_block_a,&
                              pa_block_a,p_block_a,pb_block_b,pa_block_b,&
                              p_block_b,force_ab,delta,error)
                       ENDIF
                     END IF
                     atom_a = atom_of_kind(iatom)
                     atom_b = atom_of_kind(jatom)
                     IF ( irow /= iatom ) force_ab = -force_ab
                     force(ikind)%rho_elec(:,atom_a) =&
                       force(ikind)%rho_elec(:,atom_a) - force_ab(:)
                     force(jkind)%rho_elec(:,atom_b) =&
                       force(jkind)%rho_elec(:,atom_b) + force_ab(:)
                   END IF

                 END IF

                 sab_orb_neighbor_node => next(sab_orb_neighbor_node)

               END DO ! inode => jatom(atom B)

             END DO ! ilist => iatom(atom A)

          END DO ! jkind

       END DO ! ikind

       IF (qs_env%qmmm) THEN
       END IF

       energy%hartree=-energy%core
       DO ispin=1,nspins
          CALL calculate_ecore(ks_matrix(ispin)%matrix,&
                               matrix_p(ispin)%matrix,ecoul,para_env)
          energy%hartree=energy%hartree+ecoul
       END DO

       energy%mulliken = 0._dp
       energy%exc = 0._dp

       energy%total = energy%core + energy%core_overlap + 0.5_dp*energy%hartree + &
            energy%qmmm_el + energy%mulliken

       IF (ionode.AND.logger%print_keys%scf_energies) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Two-electron integral energy:                  ",energy%hartree
       END IF
       IF (qs_env%qmmm) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
       END IF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_se_fock_matrix

! *****************************************************************************

  SUBROUTINE rfock1 (gss,gsp,gpp,gp2,hsp,pmat,fmat,error)
    REAL(dp), INTENT (IN)                    :: gss,gsp,gpp,gp2,hsp
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pmat
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fmat
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'rfock1', &
      routineP = moduleN//':'//routineN
    INTEGER     :: n

    n = SIZE ( pmat,1 )

    CPPrecondition(n<=4,cp_failure_level,routineP,error,failure)
    IF ( n == 1 ) THEN
      fmat(1,1) = fmat(1,1) + 0.5_dp * pmat(1,1)*gss
    END IF
    IF ( n == 4 ) THEN
      fmat(1,1) = fmat(1,1) + 0.5_dp * pmat(1,1)*gss + &
         (pmat(2,2)+pmat(3,3)+pmat(4,4))*(gsp-0.5_dp*hsp)
      fmat(2,2) = fmat(2,2) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(2,2)*gpp + &
                              (pmat(3,3)+pmat(4,4))*gp2 - &
                              0.25_dp * (pmat(3,3)+pmat(4,4))*(gpp-gp2)
      fmat(3,3) = fmat(3,3) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(3,3)*gpp + &
                              (pmat(2,2)+pmat(4,4))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(4,4))*(gpp-gp2)
      fmat(4,4) = fmat(4,4) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(4,4)*gpp + &
                              (pmat(2,2)+pmat(3,3))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(3,3))*(gpp-gp2)

      fmat(1,2) = fmat(1,2) + 0.5_dp * pmat(1,2)*(3._dp*hsp-gsp)
      fmat(2,1) = fmat(1,2)
      fmat(1,3) = fmat(1,3) + 0.5_dp * pmat(1,3)*(3._dp*hsp-gsp)
      fmat(3,1) = fmat(1,3)
      fmat(1,4) = fmat(1,4) + 0.5_dp * pmat(1,4)*(3._dp*hsp-gsp)
      fmat(4,1) = fmat(1,4)

      fmat(2,3) = fmat(2,3) + pmat(2,3)*(gpp-gp2) - 0.25_dp * pmat(2,3)*(gpp+gp2)
      fmat(3,2) = fmat(2,3)
      fmat(2,4) = fmat(2,4) + pmat(2,4)*(gpp-gp2) - 0.25_dp * pmat(2,4)*(gpp+gp2)
      fmat(4,2) = fmat(2,4)
      fmat(3,4) = fmat(3,4) + pmat(3,4)*(gpp-gp2) - 0.25_dp * pmat(3,4)*(gpp+gp2)
      fmat(4,3) = fmat(3,4)
    END IF

  END SUBROUTINE rfock1

  SUBROUTINE ufock1 (gss,gsp,gpp,gp2,hsp,pmata,fmata,pmatb,fmatb,error)
    REAL(dp), INTENT (IN)                    :: gss,gsp,gpp,gp2,hsp
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pmata,pmatb
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fmata,fmatb
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'ufock1', &
      routineP = moduleN//':'//routineN
    INTEGER     :: n
    REAL(dp)    :: pss,pxx,pyy,pzz

    n = SIZE ( pmata,1 )

    CPPrecondition(n<=4,cp_failure_level,routineP,error,failure)
    IF ( n == 1 ) THEN
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss
    END IF
    IF ( n == 4 ) THEN
      pss = pmata(1,1)+pmatb(1,1)
      pxx = pmata(2,2)+pmatb(2,2)
      pyy = pmata(3,3)+pmatb(3,3)
      pzz = pmata(4,4)+pmatb(4,4)
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss + &
         (pxx+pyy+pzz)*gsp - (pmata(2,2)+pmata(3,3)+pmata(4,4))*hsp
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss + &
         (pxx+pyy+pzz)*gsp - (pmatb(2,2)+pmatb(3,3)+pmatb(4,4))*hsp
      fmata(2,2) = fmata(2,2) + pmatb(2,2)*gpp - pmata(1,1)*hsp + pss*gsp + (pyy+pzz)*gp2 -&
                    0.5_dp*(pmata(3,3)+pmata(4,4))*(gpp-gp2)
      fmata(3,3) = fmata(3,3) + pmatb(3,3)*gpp - pmata(1,1)*hsp + pss*gsp + (pxx+pzz)*gp2 -&
                    0.5_dp*(pmata(2,2)+pmata(4,4))*(gpp-gp2)
      fmata(4,4) = fmata(4,4) + pmatb(4,4)*gpp - pmata(1,1)*hsp + pss*gsp + (pxx+pyy)*gp2 -&
                    0.5_dp*(pmata(2,2)+pmata(3,3))*(gpp-gp2)
      fmatb(2,2) = fmatb(2,2) + pmata(2,2)*gpp - pmatb(1,1)*hsp + pss*gsp + (pyy+pzz)*gp2 -&
                    0.5_dp*(pmatb(3,3)+pmatb(4,4))*(gpp-gp2)
      fmatb(3,3) = fmatb(3,3) + pmata(3,3)*gpp - pmatb(1,1)*hsp + pss*gsp + (pxx+pzz)*gp2 -&
                    0.5_dp*(pmatb(2,2)+pmatb(4,4))*(gpp-gp2)
      fmatb(4,4) = fmatb(4,4) + pmata(4,4)*gpp - pmatb(1,1)*hsp + pss*gsp + (pxx+pyy)*gp2 -&
                    0.5_dp*(pmatb(2,2)+pmatb(3,3))*(gpp-gp2)

      fmata(1,2) = fmata(1,2) + pmata(1,2)*(hsp-gsp) + 2._dp*pmatb(1,2)*hsp
      fmata(2,1) = fmata(1,2)
      fmata(1,3) = fmata(1,3) + pmata(1,3)*(hsp-gsp) + 2._dp*pmatb(1,3)*hsp
      fmata(3,1) = fmata(1,3)
      fmata(1,4) = fmata(1,4) + pmata(1,4)*(hsp-gsp) + 2._dp*pmatb(1,4)*hsp
      fmata(4,1) = fmata(1,4)
      fmatb(1,2) = fmatb(1,2) + pmatb(1,2)*(hsp-gsp) + 2._dp*pmata(1,2)*hsp
      fmatb(2,1) = fmatb(1,2)
      fmatb(1,3) = fmatb(1,3) + pmatb(1,3)*(hsp-gsp) + 2._dp*pmata(1,3)*hsp
      fmatb(3,1) = fmatb(1,3)
      fmatb(1,4) = fmatb(1,4) + pmatb(1,4)*(hsp-gsp) + 2._dp*pmata(1,4)*hsp
      fmatb(4,1) = fmatb(1,4)

      fmata(2,3) = fmata(2,3) + (pmata(2,3)+pmatb(2,3))*(gpp-gp2) &
                              - 0.5_dp * pmata(2,3)*(gpp+gp2)
      fmata(3,2) = fmata(2,3)
      fmata(2,4) = fmata(2,4) + (pmata(2,4)+pmatb(2,4))*(gpp-gp2) &
                              - 0.5_dp * pmata(2,4)*(gpp+gp2)
      fmata(4,2) = fmata(2,4)
      fmata(3,4) = fmata(3,4) + (pmata(3,4)+pmatb(3,4))*(gpp-gp2) &
                              - 0.5_dp * pmata(3,4)*(gpp+gp2)
      fmata(4,3) = fmata(3,4)
      fmatb(2,3) = fmatb(2,3) + (pmata(2,3)+pmatb(2,3))*(gpp-gp2) &
                              - 0.5_dp * pmatb(2,3)*(gpp+gp2)
      fmatb(3,2) = fmatb(2,3)
      fmatb(2,4) = fmatb(2,4) + (pmata(2,4)+pmatb(2,4))*(gpp-gp2) &
                              - 0.5_dp * pmatb(2,4)*(gpp+gp2)
      fmatb(4,2) = fmatb(2,4)
      fmatb(3,4) = fmatb(3,4) + (pmata(3,4)+pmatb(3,4))*(gpp-gp2) &
                              - 0.5_dp * pmatb(3,4)*(gpp+gp2)
      fmatb(4,3) = fmatb(3,4)

    END IF
   
  END SUBROUTINE ufock1

! *****************************************************************************

  SUBROUTINE rfock2(sepa,sepb,rij,pa,pb,pab,fa,fb,fab,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa,pb,pab
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa,fb,fab
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2', &
      routineP = moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(100)          :: wint
    REAL(KIND = dp), DIMENSION(10)           :: ppa,ppb
    REAL(KIND = dp), DIMENSION(10)           :: ca,cb
    REAL(KIND = dp), DIMENSION(4,4)          :: kab
    INTEGER                                  :: na,nb,i,ii

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fa(1,1) = fa(1,1) + pb(1,1)*wint(1)
      fb(1,1) = fb(1,1) + pa(1,1)*wint(1)
      fab(1,1) = fab(1,1) - 0.5_dp * pab(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      fa(1,1) = fa(1,1) + SUM ( ppb(1:10)*wint(1:10) )
      ppa(1:10) = pa(1,1)*wint(1:10)
      CALL expand( fb, ppa )
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(1,4) + wint(4)*pab(1,2) + wint(7)*pab(1,3) )
      fab(1,2) = fab(1,2) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(1,4) + wint(6)*pab(1,2) + wint(9)*pab(1,3) )
      fab(1,3) = fab(1,3) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(1,4) + wint(9)*pab(1,2) + wint(10)*pab(1,3) )
      fab(1,4) = fab(1,4) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(1,4) + wint(5)*pab(1,2) + wint(8)*pab(1,3) )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      fb(1,1) = fb(1,1) + SUM ( ppa(1:10)*wint(1:10) )
      ppb(1:10) = pb(1,1)*wint(1:10)
      CALL expand( fa, ppb )
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(4,1) + wint(4)*pab(2,1) + wint(7)*pab(3,1) )
      fab(2,1) = fab(2,1) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(4,1) + wint(6)*pab(2,1) + wint(9)*pab(3,1) )
      fab(3,1) = fab(3,1) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(4,1) + wint(9)*pab(2,1) + wint(10)*pab(3,1) )
      fab(4,1) = fab(4,1) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(4,1) + wint(5)*pab(2,1) + wint(8)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa, ppa )
      CALL contract( pb, ppb )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*ppb(1:10) )
        cb(i) = SUM ( wint(i:90+i:10)*ppa(1:10) )
      END DO
      CALL expand( fa, ca )
      CALL expand( fb, cb )
      CALL exchange_sp ( kab, pab, wint )
      fab(1:4,1:4) = fab(1:4,1:4) - 0.5_dp*kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rfock2

  SUBROUTINE ufock2(sepa,sepb,rij,pa_a,pb_a,pab_a,fa_a,fb_a,fab_a,&
                    pa_b,pb_b,pab_b,fa_b,fb_b,fab_b,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_a,pb_a,pab_a
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_b,pb_b,pab_b
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa_a,fb_a,fab_a
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa_b,fb_b,fab_b
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND = dp)                          :: ess
    REAL(KIND = dp), DIMENSION(100)          :: wint

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2', &
      routineP = moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(10)           :: ppa_a,ppa_b,ppb_a,ppb_b
    REAL(KIND = dp), DIMENSION(10)           :: ca,cb
    REAL(KIND = dp), DIMENSION(4,4)          :: kab
    INTEGER                                  :: na,nb,i,ii

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fa_a(1,1) = fa_a(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_a(1,1) = fb_a(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
      fab_a(1,1) = fab_a(1,1) - pab_a(1,1)*wint(1)
      fa_b(1,1) = fa_b(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_b(1,1) = fb_b(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
      fab_b(1,1) = fab_b(1,1) - pab_b(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * wint(1:10) )
      fa_a(1,1) = fa_a(1,1) + ess
      fa_b(1,1) = fa_b(1,1) + ess
      ppa_a(1:10) = (pa_a(1,1)+pa_b(1,1))*wint(1:10)
      CALL expand( fb_a, ppa_a )
      CALL expand( fb_b, ppa_a )
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(1,4) - wint(4)*pab_a(1,2) - wint(7)*pab_a(1,3)
      fab_a(1,2) = fab_a(1,2) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(1,4) - wint(6)*pab_a(1,2) - wint(9)*pab_a(1,3)
      fab_a(1,3) = fab_a(1,3) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(1,4) - wint(9)*pab_a(1,2) - wint(10)*pab_a(1,3)
      fab_a(1,4) = fab_a(1,4) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(1,4) - wint(5)*pab_a(1,2) - wint(8)*pab_a(1,3)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(1,4) - wint(4)*pab_b(1,2) - wint(7)*pab_b(1,3)
      fab_b(1,2) = fab_b(1,2) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(1,4) - wint(6)*pab_b(1,2) - wint(9)*pab_b(1,3)
      fab_b(1,3) = fab_b(1,3) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(1,4) - wint(9)*pab_b(1,2) - wint(10)*pab_b(1,3)
      fab_b(1,4) = fab_b(1,4) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(1,4) - wint(5)*pab_b(1,2) - wint(8)*pab_b(1,3)
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * wint(1:10) )
      fb_a(1,1) = fb_a(1,1) + ess
      fb_b(1,1) = fb_b(1,1) + ess
      ppb_a(1:10) = (pb_a(1,1)+pb_b(1,1))*wint(1:10)
      CALL expand( fa_a, ppb_a )
      CALL expand( fa_b, ppb_a )
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(4,1) - wint(4)*pab_a(2,1) - wint(7)*pab_a(3,1)
      fab_a(2,1) = fab_a(2,1) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(4,1) - wint(6)*pab_a(2,1) - wint(9)*pab_a(3,1)
      fab_a(3,1) = fab_a(3,1) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(4,1) - wint(9)*pab_a(2,1) - wint(10)*pab_a(3,1)
      fab_a(4,1) = fab_a(4,1) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(4,1) - wint(5)*pab_a(2,1) - wint(8)*pab_a(3,1)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(4,1) - wint(4)*pab_b(2,1) - wint(7)*pab_b(3,1)
      fab_b(2,1) = fab_b(2,1) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(4,1) - wint(6)*pab_b(2,1) - wint(9)*pab_b(3,1)
      fab_b(3,1) = fab_b(3,1) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(4,1) - wint(9)*pab_b(2,1) - wint(10)*pab_b(3,1)
      fab_b(4,1) = fab_b(4,1) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(4,1) - wint(5)*pab_b(2,1) - wint(8)*pab_b(3,1)
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*(ppb_a(1:10)+ppb_b(1:10)) )
        cb(i) = SUM ( wint(i:90+i:10)*(ppa_a(1:10)+ppa_b(1:10)) )
      END DO
      CALL expand( fa_a, ca )
      CALL expand( fa_b, ca )
      CALL expand( fb_a, cb )
      CALL expand( fb_b, cb )
      CALL exchange_sp ( kab, pab_a, wint )
      fab_a(1:4,1:4) = fab_a(1:4,1:4) - kab(1:4,1:4)
      CALL exchange_sp ( kab, pab_b, wint )
      fab_b(1:4,1:4) = fab_b(1:4,1:4) - kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE ufock2

! *****************************************************************************

  SUBROUTINE rdfock2(sepa,sepb,rij,pa,pb,pab,force,delta,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa,pb,pab
    REAL(dp), DIMENSION(:), INTENT (INOUT)   :: force
    REAL(dp), INTENT (IN)                    :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'rdfock2', &
      routineP = moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(100,3)        :: dwint
    REAL(KIND = dp), DIMENSION(10)           :: ppa,ppb
    REAL(KIND = dp), DIMENSION(10)           :: ca
    REAL(KIND = dp), DIMENSION(4,4)          :: kab
    INTEGER                                  :: na,nb,i,ii,j

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) + pa(1,1)*pb(1,1)*dwint(1,:)
      force(:) = force(:) - 0.5_dp*pab(1,1)*pab(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      DO i=1,3
        force(i) = force(i) + pa(1,1) * SUM ( ppb(1:10)*dwint(1:10,i) )
      END DO
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(1,4) + dwint(4,:)*pab(1,2) + dwint(7,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,2)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(1,4) + dwint(6,:)*pab(1,2) + dwint(9,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,3)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(1,4) + dwint(9,:)*pab(1,2) + dwint(10,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,4)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(1,4) + dwint(5,:)*pab(1,2) + dwint(8,:)*pab(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      DO i=1,3
        force(i) = force(i) + pb(1,1) * SUM ( ppa(1:10)*dwint(1:10,i) )
      END DO
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(4,1) + dwint(4,:)*pab(2,1) + dwint(7,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(2,1)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(4,1) + dwint(6,:)*pab(2,1) + dwint(9,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(3,1)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(4,1) + dwint(9,:)*pab(2,1) + dwint(10,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(4,1)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(4,1) + dwint(5,:)*pab(2,1) + dwint(8,:)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      CALL contract( pa, ppa )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*ppb(1:10) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*ppa(1:10) )
        CALL exchange_sp ( kab, pab, dwint(:,j) )
        force(j) = force(j) - 0.5_dp * SUM(pab*kab)
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rdfock2

  SUBROUTINE udfock2(sepa,sepb,rij,pa_a,pb_a,pab_a,pa_b,pb_b,pab_b,&
             force,delta,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_a,pb_a,pab_a
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_b,pb_b,pab_b
    REAL(dp), DIMENSION(:), INTENT (INOUT)   :: force
    REAL(dp), INTENT (IN)                    :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'udfock2', &
      routineP = moduleN//':'//routineN
    REAL(KIND = dp), DIMENSION(100,3)        :: dwint
    REAL(KIND = dp), DIMENSION(10)           :: ppa_a,ppa_b,ppb_a,ppb_b
    REAL(KIND = dp), DIMENSION(10)           :: ca
    REAL(KIND = dp), DIMENSION(4,4)          :: kab
    INTEGER                                  :: na,nb,i,ii,j
    REAL(dp)                                 :: pta,ptb,ess

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      pta = pa_a(1,1)+pa_b(1,1)
      ptb = pb_a(1,1)+pb_b(1,1)
      force(:) = force(:) + pta*ptb*dwint(1,:)
      force(:) = force(:) - pab_a(1,1)*pab_a(1,1)*dwint(1,:)
      force(:) = force(:) - pab_b(1,1)*pab_b(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,3
        ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pa_a(1,1)+pa_b(1,1)) * ess
      END DO
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(1,4) + dwint(4,:)*pab_a(1,2) + dwint(7,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,2)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(1,4) + dwint(6,:)*pab_a(1,2) + dwint(9,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,3)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(1,4) + dwint(9,:)*pab_a(1,2) + dwint(10,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,4)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(1,4) + dwint(5,:)*pab_a(1,2) + dwint(8,:)*pab_a(1,3) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(1,4) + dwint(4,:)*pab_b(1,2) + dwint(7,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,2)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(1,4) + dwint(6,:)*pab_b(1,2) + dwint(9,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,3)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(1,4) + dwint(9,:)*pab_b(1,2) + dwint(10,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,4)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(1,4) + dwint(5,:)*pab_b(1,2) + dwint(8,:)*pab_b(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      DO i=1,3
        ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pb_a(1,1)+pb_b(1,1)) * ess
      END DO
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(4,1) + dwint(4,:)*pab_a(2,1) + dwint(7,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(2,1)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(4,1) + dwint(6,:)*pab_a(2,1) + dwint(9,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(3,1)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(4,1) + dwint(9,:)*pab_a(2,1) + dwint(10,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(4,1)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(4,1) + dwint(5,:)*pab_a(2,1) + dwint(8,:)*pab_a(3,1) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(4,1) + dwint(4,:)*pab_b(2,1) + dwint(7,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(2,1)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(4,1) + dwint(6,:)*pab_b(2,1) + dwint(9,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(3,1)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(4,1) + dwint(9,:)*pab_b(2,1) + dwint(10,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(4,1)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(4,1) + dwint(5,:)*pab_b(2,1) + dwint(8,:)*pab_b(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*(ppb_a(1:10)+ppb_b(1:10)) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*(ppa_a(1:10)+ppa_b(1:10)) )
        CALL exchange_sp ( kab, pab_a, dwint(:,j) )
        force(j) = force(j) - SUM(pab_a*kab)
        CALL exchange_sp ( kab, pab_b, dwint(:,j) )
        force(j) = force(j) - SUM(pab_b*kab)
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE udfock2

! *****************************************************************************

  SUBROUTINE contract ( amat, avec )
    REAL ( dp ), INTENT ( IN ), DIMENSION (:,:)    :: amat
    REAL ( dp ), INTENT ( OUT ), DIMENSION (:)     :: avec

    avec( 1) = amat(1,1) 
    avec( 2) = amat(4,1) + amat(1,4)
    avec( 3) = amat(4,4) 
    avec( 4) = amat(2,1) + amat(1,2)
    avec( 5) = amat(2,4) + amat(4,2)
    avec( 6) = amat(2,2) 
    avec( 7) = amat(3,1) + amat(1,3)
    avec( 8) = amat(3,4) + amat(4,3)
    avec( 9) = amat(3,2) + amat(2,3)
    avec(10) = amat(3,3) 

  END SUBROUTINE contract 

  SUBROUTINE expand ( amat, avec )
    REAL ( dp ), INTENT ( OUT ), DIMENSION (:,:)  :: amat
    REAL ( dp ), INTENT ( IN ), DIMENSION (:)     :: avec

    amat(1,1) = amat(1,1) + avec( 1)
    amat(4,1) = amat(4,1) + avec( 2)
    amat(1,4) = amat(1,4) + avec( 2)
    amat(4,4) = amat(4,4) + avec( 3)
    amat(2,1) = amat(2,1) + avec( 4)
    amat(1,2) = amat(1,2) + avec( 4)
    amat(2,4) = amat(2,4) + avec( 5)
    amat(4,2) = amat(4,2) + avec( 5)
    amat(2,2) = amat(2,2) + avec( 6)
    amat(3,1) = amat(3,1) + avec( 7)
    amat(1,3) = amat(1,3) + avec( 7)
    amat(3,4) = amat(3,4) + avec( 8)
    amat(4,3) = amat(4,3) + avec( 8)
    amat(3,2) = amat(3,2) + avec( 9)
    amat(2,3) = amat(2,3) + avec( 9)
    amat(3,3) = amat(3,3) + avec(10)

  END SUBROUTINE expand

  SUBROUTINE exchange_sp ( kab, pab, w )
    REAL ( dp ), INTENT ( OUT ), DIMENSION (:,:)  :: kab
    REAL ( dp ), INTENT ( IN ), DIMENSION (:,:)   :: pab
    REAL ( dp ), INTENT ( IN ), DIMENSION (:)     :: w
    REAL ( dp ), DIMENSION (16)                   :: pk

    pk( 1) = pab(1,1)
    pk( 2) = pab(1,4)
    pk( 3) = pab(1,2)
    pk( 4) = pab(1,3)
    pk( 5) = pab(4,1)
    pk( 6) = pab(4,4)
    pk( 7) = pab(4,2)
    pk( 8) = pab(4,3)
    pk( 9) = pab(2,1)
    pk(10) = pab(2,4)
    pk(11) = pab(2,2)
    pk(12) = pab(2,3)
    pk(13) = pab(3,1)
    pk(14) = pab(3,4)
    pk(15) = pab(3,2)
    pk(16) = pab(3,3)

    kab(1,1)= &
      +PK( 1)*W(  1)+PK( 2)*W(  2)+PK( 3)*W(  4)+PK( 4)*W(  7) &
      +PK( 5)*W( 11)+PK( 6)*W( 12)+PK( 7)*W( 14)+PK( 8)*W( 17) &
      +PK( 9)*W( 31)+PK(10)*W( 32)+PK(11)*W( 34)+PK(12)*W( 37) &
      +PK(13)*W( 61)+PK(14)*W( 62)+PK(15)*W( 64)+PK(16)*W( 67)
    kab(1,4)= &
      +PK( 1)*W(  2)+PK( 2)*W(  3)+PK( 3)*W(  5)+PK( 4)*W(  8) &
      +PK( 5)*W( 12)+PK( 6)*W( 13)+PK( 7)*W( 15)+PK( 8)*W( 18) &
      +PK( 9)*W( 32)+PK(10)*W( 33)+PK(11)*W( 35)+PK(12)*W( 38) &
      +PK(13)*W( 62)+PK(14)*W( 63)+PK(15)*W( 65)+PK(16)*W( 68)
    kab(1,2)= &
      +PK( 1)*W(  4)+PK( 2)*W(  5)+PK( 3)*W(  6)+PK( 4)*W(  9) &
      +PK( 5)*W( 14)+PK( 6)*W( 15)+PK( 7)*W( 16)+PK( 8)*W( 19) &
      +PK( 9)*W( 34)+PK(10)*W( 35)+PK(11)*W( 36)+PK(12)*W( 39) &
      +PK(13)*W( 64)+PK(14)*W( 65)+PK(15)*W( 66)+PK(16)*W( 69)
    kab(1,3)= &
      +PK( 1)*W(  7)+PK( 2)*W(  8)+PK( 3)*W(  9)+PK( 4)*W( 10) &
      +PK( 5)*W( 17)+PK( 6)*W( 18)+PK( 7)*W( 19)+PK( 8)*W( 20) &
      +PK( 9)*W( 37)+PK(10)*W( 38)+PK(11)*W( 39)+PK(12)*W( 40) &
      +PK(13)*W( 67)+PK(14)*W( 68)+PK(15)*W( 69)+PK(16)*W( 70)
    kab(4,1)= &
      +PK( 1)*W( 11)+PK( 2)*W( 12)+PK( 3)*W( 14)+PK( 4)*W( 17) &
      +PK( 5)*W( 21)+PK( 6)*W( 22)+PK( 7)*W( 24)+PK( 8)*W( 27) &
      +PK( 9)*W( 41)+PK(10)*W( 42)+PK(11)*W( 44)+PK(12)*W( 47) &
      +PK(13)*W( 71)+PK(14)*W( 72)+PK(15)*W( 74)+PK(16)*W( 77)
    kab(4,4)= &
      +PK( 1)*W( 12)+PK( 2)*W( 13)+PK( 3)*W( 15)+PK( 4)*W( 18) &
      +PK( 5)*W( 22)+PK( 6)*W( 23)+PK( 7)*W( 25)+PK( 8)*W( 28) &
      +PK( 9)*W( 42)+PK(10)*W( 43)+PK(11)*W( 45)+PK(12)*W( 48) &
      +PK(13)*W( 72)+PK(14)*W( 73)+PK(15)*W( 75)+PK(16)*W( 78)
    kab(4,2)= &
      +PK( 1)*W( 14)+PK( 2)*W( 15)+PK( 3)*W( 16)+PK( 4)*W( 19) &
      +PK( 5)*W( 24)+PK( 6)*W( 25)+PK( 7)*W( 26)+PK( 8)*W( 29) &
      +PK( 9)*W( 44)+PK(10)*W( 45)+PK(11)*W( 46)+PK(12)*W( 49) &
      +PK(13)*W( 74)+PK(14)*W( 75)+PK(15)*W( 76)+PK(16)*W( 79)
    kab(4,3)= &
      +PK( 1)*W( 17)+PK( 2)*W( 18)+PK( 3)*W( 19)+PK( 4)*W( 20) &
      +PK( 5)*W( 27)+PK( 6)*W( 28)+PK( 7)*W( 29)+PK( 8)*W( 30) &
      +PK( 9)*W( 47)+PK(10)*W( 48)+PK(11)*W( 49)+PK(12)*W( 50) &
      +PK(13)*W( 77)+PK(14)*W( 78)+PK(15)*W( 79)+PK(16)*W( 80)
    kab(2,1)= &
      +PK( 1)*W( 31)+PK( 2)*W( 32)+PK( 3)*W( 34)+PK( 4)*W( 37) &
      +PK( 5)*W( 41)+PK( 6)*W( 42)+PK( 7)*W( 44)+PK( 8)*W( 47) &
      +PK( 9)*W( 51)+PK(10)*W( 52)+PK(11)*W( 54)+PK(12)*W( 57) &
      +PK(13)*W( 81)+PK(14)*W( 82)+PK(15)*W( 84)+PK(16)*W( 87)
    kab(2,4)= &
      +PK( 1)*W( 32)+PK( 2)*W( 33)+PK( 3)*W( 35)+PK( 4)*W( 38) &
      +PK( 5)*W( 42)+PK( 6)*W( 43)+PK( 7)*W( 45)+PK( 8)*W( 48) &
      +PK( 9)*W( 52)+PK(10)*W( 53)+PK(11)*W( 55)+PK(12)*W( 58) &
      +PK(13)*W( 82)+PK(14)*W( 83)+PK(15)*W( 85)+PK(16)*W( 88)
    kab(2,2)= &
      +PK( 1)*W( 34)+PK( 2)*W( 35)+PK( 3)*W( 36)+PK( 4)*W( 39) &
      +PK( 5)*W( 44)+PK( 6)*W( 45)+PK( 7)*W( 46)+PK( 8)*W( 49) &
      +PK( 9)*W( 54)+PK(10)*W( 55)+PK(11)*W( 56)+PK(12)*W( 59) &
      +PK(13)*W( 84)+PK(14)*W( 85)+PK(15)*W( 86)+PK(16)*W( 89)
    kab(2,3)= &
      +PK( 1)*W( 37)+PK( 2)*W( 38)+PK( 3)*W( 39)+PK( 4)*W( 40) &
      +PK( 5)*W( 47)+PK( 6)*W( 48)+PK( 7)*W( 49)+PK( 8)*W( 50) &
      +PK( 9)*W( 57)+PK(10)*W( 58)+PK(11)*W( 59)+PK(12)*W( 60) &
      +PK(13)*W( 87)+PK(14)*W( 88)+PK(15)*W( 89)+PK(16)*W( 90)
    kab(3,1)= &
      +PK( 1)*W( 61)+PK( 2)*W( 62)+PK( 3)*W( 64)+PK( 4)*W( 67) &
      +PK( 5)*W( 71)+PK( 6)*W( 72)+PK( 7)*W( 74)+PK( 8)*W( 77) &
      +PK( 9)*W( 81)+PK(10)*W( 82)+PK(11)*W( 84)+PK(12)*W( 87) &
      +PK(13)*W( 91)+PK(14)*W( 92)+PK(15)*W( 94)+PK(16)*W( 97)
    kab(3,4)= &
      +PK( 1)*W( 62)+PK( 2)*W( 63)+PK( 3)*W( 65)+PK( 4)*W( 68) &
      +PK( 5)*W( 72)+PK( 6)*W( 73)+PK( 7)*W( 75)+PK( 8)*W( 78) &
      +PK( 9)*W( 82)+PK(10)*W( 83)+PK(11)*W( 85)+PK(12)*W( 88) &
      +PK(13)*W( 92)+PK(14)*W( 93)+PK(15)*W( 95)+PK(16)*W( 98)
    kab(3,2)= &
      +PK( 1)*W( 64)+PK( 2)*W( 65)+PK( 3)*W( 66)+PK( 4)*W( 69) &
      +PK( 5)*W( 74)+PK( 6)*W( 75)+PK( 7)*W( 76)+PK( 8)*W( 79) &
      +PK( 9)*W( 84)+PK(10)*W( 85)+PK(11)*W( 86)+PK(12)*W( 89) &
      +PK(13)*W( 94)+PK(14)*W( 95)+PK(15)*W( 96)+PK(16)*W( 99)
    kab(3,3)= &
      +PK( 1)*W( 67)+PK( 2)*W( 68)+PK( 3)*W( 69)+PK( 4)*W( 70) &
      +PK( 5)*W( 77)+PK( 6)*W( 78)+PK( 7)*W( 79)+PK( 8)*W( 80) &
      +PK( 9)*W( 87)+PK(10)*W( 88)+PK(11)*W( 89)+PK(12)*W( 90) &
      +PK(13)*W( 97)+PK(14)*W( 98)+PK(15)*W( 99)+PK(16)*W(100)

  END SUBROUTINE exchange_sp

! *****************************************************************************

END MODULE qs_semi_empirical_hamiltonian

! *****************************************************************************
