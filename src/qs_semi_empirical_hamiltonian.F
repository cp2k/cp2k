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
  USE physcon,                         ONLY: evolt
  USE qs_overlap,                      ONLY: build_overlap_matrix,&
                                             write_sparse_matrix
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             get_qs_env,&
                                             set_qs_env
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_energy_types,                 ONLY: qs_energy_type
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
  USE semi_empirical_integrals,        ONLY: rotint
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_param
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
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
                               POINTER       :: matrix_h,matrix_s
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

    INTEGER                                  :: i,j,k,ikind,jkind,nkind
    INTEGER                                  :: iab,ilist,nlist,iatom,jatom,&
                                                inode,nnode,&
                                                irow,icol,natorb_a,natorb_b
    LOGICAL                                  :: defined,failure
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_block,s_block
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_blocka,h_blockb
    REAL(KIND = dp), DIMENSION(:), POINTER   :: beta_a,beta_b
    REAL(KIND = dp), DIMENSION(16)           :: ha,hb,ua
    INTEGER                                  :: za,zb,itype,natom
    REAL(KIND = dp)                          :: kh,enuc,uss,upp,udd,uff
    REAL(KIND = dp), DIMENSION(3)            :: rij
    REAL(KIND = dp), DIMENSION(100)          :: wint
    REAL(KIND = dp), DIMENSION(10)           :: e1b,e2a
    REAL(KIND = dp)                          :: enuclear,eisol,eheat,econst

!   ---------------------------------------------------------------------------

     NULLIFY ( matrix_h,matrix_s,atomic_kind_set,sab_orb )

     ! calculate overlap matrix
     CALL build_overlap_matrix(qs_env=qs_env,globenv=globenv)

     CALL get_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h,&
                     matrix_s=matrix_s)

     ! Allocate the core Hamiltonian matrix
     CALL allocate_matrix_set(matrix_h,1)
     CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                     target=matrix_h(1)%matrix,&
                                     target_name="CORE HAMILTONIAN MATRIX")

     CALL get_qs_env(qs_env=qs_env,&
                     sab_orb=sab_orb,&
                     atomic_kind_set=atomic_kind_set,&
                     dft_control=dft_control)

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
               sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
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
                 IF (.NOT.ASSOCIATED(h_blocka)) THEN
                   CALL stop_program(routineN,moduleN,__LINE__,&
                                    "A hmatrix block is missing",globenv)
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
                 CPPrecondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

                 ! two-centre one-electron term
                 NULLIFY(s_block)
                 CALL get_block_node(matrix=matrix_s(1)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=s_block)
                 CPPrecondition(ASSOCIATED(s_block),cp_failure_level,routineP,error,failure)
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

                 SELECT CASE (dft_control%qs_control%method)
                   CASE ("EHT")
                   CASE ("MNDO","AM1","PM3")
                     IF ( irow == iatom ) THEN
                       CALL rotint (se_kind_a,se_kind_b,rij,wint,e1b,e2a,&
                                    enuc,itype)
                     ELSE
                       CALL rotint (se_kind_b,se_kind_a,-rij,wint,e2a,e1b,&
                                    enuc,itype)
                     END IF
                     enuclear = enuclear + enuc
                     ! one-centre one-electron terms
                     NULLIFY(h_blockb)
                     CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=h_blockb)
                     CPPrecondition(ASSOCIATED(h_blockb),cp_failure_level,routineP,error,failure)
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

                 END SELECT

               END IF

               sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     CALL set_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h)

     qs_env%energy%core_overlap = enuclear
     qs_env%energy%core_self    = econst

     IF (globenv%print%core_hamiltonian_matrix) THEN
       CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,globenv,&
             scale=evolt)
     END IF

     !deallocate overlap matrix and initialize with unit matrix 
     !for orthogonal methods
     IF ( dft_control%qs_control%se_control%orthogonal_basis ) THEN
       !
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
                                                ispin,itype,ikind,nkind,jkind
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
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: p_block_a,p_block_b,&
                                                pa_block_a,pa_block_b,&
                                                pb_block_a,pb_block_b,&
                                                ks_block_a,ks_block_b,&
                                                ksa_block_a,ksa_block_b,&
                                                ksb_block_a,ksb_block_b

!   ----------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,logger)

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

       SELECT CASE (dft_control%qs_control%method)
         CASE DEFAULT
           itype=0
         CASE ("EHT")
           itype=0
         CASE ("AM1")
           itype=2
         CASE ("PM3")
           itype=3
         CASE ("MNDO")
           itype=1
       END SELECT

       nkind = SIZE(atomic_kind_set)

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

               IF ( nspins >= 1 ) THEN
                 CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=pa_block_a)
                 CPPrecondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=ksa_block_a)
                 CPPrecondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
               END IF
               IF ( nspins >= 2 ) THEN
                 CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=pa_block_b)
                 CPPrecondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=ksa_block_b)
                 CPPrecondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
               END IF
               CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)

               DO inode=1,nnode
                 CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                        neighbor=jatom,r=rij)

                 IF ( iatom == jatom ) THEN
                   IF ( nspins == 1 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,pa_block_a,ksa_block_a,error)
                   ELSE IF ( nspins == 2 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,pa_block_a,ksa_block_a,&
                                pa_block_b,ksa_block_b,error)
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)

                 ELSE

                   IF ( nspins >= 1 ) THEN
                     CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=pb_block_a)
                     CPPrecondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=ksb_block_a)
                     CPPrecondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                   END IF
                   IF ( nspins >= 2 ) THEN
                     CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=pb_block_b)
                     CPPrecondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         block=ksb_block_b)
                     CPPrecondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)

                   IF (iatom <= jatom) THEN
                     irow = iatom
                     icol = jatom
                   ELSE
                     irow = jatom
                     icol = iatom
                   END IF
                   IF ( nspins >= 1 ) THEN
                     CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=ks_block_a)
                     CPPrecondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=p_block_a)
                     CPPrecondition(ASSOCIATED(p_block_a),cp_failure_level,routineP,error,failure)
                   END IF
                   IF ( nspins >= 2 ) THEN
                     CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=ks_block_b)
                     CPPrecondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                         block_row=irow,&
                                         block_col=icol,&
                                         block=p_block_b)
                     CPPrecondition(ASSOCIATED(p_block_b),cp_failure_level,routineP,error,failure)
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)

                   IF ( nspins == 1 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,p_block_a,&
                                  ksa_block_a,ksb_block_a,ks_block_a,itype,error)
                     ELSE
                       CALL fock2(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,p_block_a,&
                                  ksb_block_a,ksa_block_a,ks_block_a,itype,error)
                     ENDIF
                   ELSE IF ( nspins == 2 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,p_block_a,&
                                  ksa_block_a,ksb_block_a,ks_block_a,&
                                  pa_block_b,pb_block_b,p_block_b,&
                                  ksa_block_b,ksb_block_b,ks_block_b,itype,error)
                     ELSE
                       CALL fock2(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,p_block_a,&
                                  ksb_block_a,ksa_block_a,ks_block_a,&
                                  pb_block_a,pa_block_a,p_block_a,&
                                  ksb_block_a,ksa_block_a,ks_block_a,itype,error)
                     ENDIF
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
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
                              0.25_dp * (pmat(3,3)+pmat(4,4))*(GPP-GP2)
      fmat(3,3) = fmat(3,3) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(3,3)*gpp + &
                              (pmat(2,2)+pmat(4,4))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(4,4))*(GPP-GP2)
      fmat(4,4) = fmat(4,4) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(4,4)*gpp + &
                              (pmat(2,2)+pmat(3,3))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(3,3))*(GPP-GP2)

      fmat(1,2) = fmat(1,2) + 0.5_dp * pmat(1,2)*(3._dp*hsp-gsp)
      fmat(2,1) = fmat(1,2)
      fmat(1,3) = fmat(1,3) + 0.5_dp * pmat(1,3)*(3._dp*hsp-gsp)
      fmat(3,1) = fmat(1,3)
      fmat(1,4) = fmat(1,4) + 0.5_dp * pmat(1,4)*(3._dp*hsp-gsp)
      fmat(4,1) = fmat(1,4)

      fmat(2,3) = fmat(2,3) - 0.5_dp * pmat(2,3)*gp2
      fmat(3,2) = fmat(2,3)
      fmat(2,4) = fmat(2,4) - 0.5_dp * pmat(2,4)*gp2
      fmat(4,2) = fmat(2,4)
      fmat(3,4) = fmat(3,4) - 0.5_dp * pmat(3,4)*gp2
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

    n = SIZE ( pmata,1 )

    CPPrecondition(n<=4,cp_failure_level,routineP,error,failure)
    IF ( n == 1 ) THEN
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss
    END IF
    IF ( n == 4 ) THEN
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss + &
         (pmata(2,2)+pmata(3,3)+pmata(4,4))*(hsp-gsp) + &
         (pmatb(2,2)+pmatb(3,3)+pmatb(4,4))*hsp
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss + &
         (pmatb(2,2)+pmatb(3,3)+pmatb(4,4))*(hsp-gsp) + &
         (pmata(2,2)+pmata(3,3)+pmata(4,4))*hsp
      fmata(2,2) = fmata(2,2) + (pmata(1,1)+pmatb(1,1))*hsp - pmata(1,1)*gsp +&
         pmatb(2,2)*gpp - (pmata(3,3)+pmata(4,4))*gp2                          
      fmata(3,3) = fmata(3,3) + (pmata(1,1)+pmatb(1,1))*hsp - pmata(1,1)*gsp +&
         pmatb(3,3)*gpp - (pmata(2,2)+pmata(4,4))*gp2                          
      fmata(4,4) = fmata(4,4) + (pmata(1,1)+pmatb(1,1))*hsp - pmata(1,1)*gsp +&
         pmatb(4,4)*gpp - (pmata(2,2)+pmata(3,3))*gp2                          
      fmatb(2,2) = fmatb(2,2) + (pmata(1,1)+pmatb(1,1))*hsp - pmatb(1,1)*gsp +&
         pmata(2,2)*gpp - (pmatb(3,3)+pmatb(4,4))*gp2                          
      fmatb(3,3) = fmatb(3,3) + (pmata(1,1)+pmatb(1,1))*hsp - pmatb(1,1)*gsp +&
         pmata(3,3)*gpp - (pmatb(2,2)+pmatb(4,4))*gp2                          
      fmatb(4,4) = fmatb(4,4) + (pmata(1,1)+pmatb(1,1))*hsp - pmatb(1,1)*gsp +&
         pmata(4,4)*gpp - (pmatb(2,2)+pmatb(3,3))*gp2                          
      fmata(1,2) = fmata(1,2) + (pmata(1,2)+2._dp*pmatb(1,2))*gsp - pmata(1,2)*hsp
      fmata(2,1) = fmata(1,2)
      fmata(1,3) = fmata(1,3) + (pmata(1,3)+2._dp*pmatb(1,3))*gsp - pmata(1,3)*hsp
      fmata(3,1) = fmata(1,3)
      fmata(1,4) = fmata(1,4) + (pmata(1,4)+2._dp*pmatb(1,4))*gsp - pmata(1,4)*hsp
      fmata(4,1) = fmata(1,4)
      fmatb(1,2) = fmatb(1,2) + (pmatb(1,2)+2._dp*pmata(1,2))*gsp - pmatb(1,2)*hsp
      fmatb(2,1) = fmatb(1,2)
      fmatb(1,3) = fmatb(1,3) + (pmatb(1,3)+2._dp*pmata(1,3))*gsp - pmatb(1,3)*hsp
      fmatb(3,1) = fmatb(1,3)
      fmatb(1,4) = fmatb(1,4) + (pmatb(1,4)+2._dp*pmata(1,4))*gsp - pmatb(1,4)*hsp
      fmatb(4,1) = fmatb(1,4)
      fmata(2,3) = fmata(2,3) + (pmata(2,3)+2._dp*pmatb(2,3))*gp2
      fmata(3,2) = fmata(2,3)
      fmata(2,4) = fmata(2,4) + (pmata(2,4)+2._dp*pmatb(2,4))*gp2
      fmata(4,2) = fmata(2,4)
      fmata(3,4) = fmata(3,4) + (pmata(3,4)+2._dp*pmatb(3,4))*gp2
      fmata(4,3) = fmata(3,4)
      fmatb(2,3) = fmatb(2,3) + (pmatb(2,3)+2._dp*pmata(2,3))*gp2
      fmatb(3,2) = fmatb(2,3)
      fmatb(2,4) = fmatb(2,4) + (pmatb(2,4)+2._dp*pmata(2,4))*gp2
      fmatb(4,2) = fmatb(2,4)
      fmatb(3,4) = fmatb(3,4) + (pmatb(3,4)+2._dp*pmata(3,4))*gp2
      fmatb(4,3) = fmatb(3,4)
    END IF
   
  END SUBROUTINE ufock1

! *****************************************************************************
  SUBROUTINE rfock2(sepa,sepb,rij,pa,pb,pab,fa,fb,fab,itype,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa,pb,pab
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa,fb,fab
    INTEGER, INTENT (IN)                     :: itype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2', &
      routineP = moduleN//':'//routineN
    REAL(KIND = dp)                          :: enuc
    REAL(KIND = dp), DIMENSION(100)          :: wint
    REAL(KIND = dp), DIMENSION(10)           :: e1b,e2a
    REAL(KIND = dp), DIMENSION(10)           :: ppa,ppb
    REAL(KIND = dp), DIMENSION(10)           :: ca,cb
    REAL(KIND = dp), DIMENSION(4,4)          :: kab
    INTEGER                                  :: na,nb,i,ii

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL rotint (sepa,sepb,rij,wint,e1b,e2a,enuc,itype)

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
                    pa_b,pb_b,pab_b,fa_b,fb_b,fab_b,itype,error)
    TYPE(semi_empirical_type), INTENT (IN)   :: sepa,sepb
    REAL(dp), DIMENSION(:), INTENT (IN)      :: rij
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_a,pb_a,pab_a
    REAL(dp), DIMENSION(:,:), INTENT (IN)    :: pa_b,pb_b,pab_b
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa_a,fb_a,fab_a
    REAL(dp), DIMENSION(:,:), INTENT (INOUT) :: fa_b,fb_b,fab_b
    INTEGER, INTENT (IN)                     :: itype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND = dp)                          :: enuc
    REAL(KIND = dp), DIMENSION(100)          :: wint
    REAL(KIND = dp), DIMENSION(10)           :: e1b,e2a

    LOGICAL                                  :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: na,nb

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL rotint (sepa,sepb,rij,wint,e1b,e2a,enuc,itype)

  END SUBROUTINE ufock2

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
