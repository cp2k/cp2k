!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
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
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: evolt
  USE qs_overlap,                      ONLY: build_overlap_matrix,&
                                             write_sparse_matrix
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             get_qs_env,&
                                             set_qs_env
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE semi_empirical_integrals,        ONLY: rotint
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_param
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_semi_empirical_hamiltonian"

  PUBLIC :: build_se_core_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_se_core_matrix(qs_env,globenv,calculate_forces)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_se_core_matrix"

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
                                                inode,nnode,last_jatom,&
                                                irow,icol,natorb_a,natorb_b
    LOGICAL                                  :: new_atom_b,defined
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_block,s_block
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_blocka,h_blockb
    REAL(KIND = dp), DIMENSION(:), POINTER   :: beta_a,beta_b
    REAL(KIND = dp), DIMENSION(16)           :: ha,hb,ua
    INTEGER                                  :: za,zb,kr,itype,natom
    REAL(KIND = dp)                          :: kh,enuc,uss,upp,udd,uff
    REAL(KIND = dp), DIMENSION(3)            :: rij
    REAL(KIND = dp), DIMENSION(72)           :: wint
    REAL(KIND = dp), DIMENSION(10)           :: e1b,e2a

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

     DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             se_parameter=se_kind_a)
        CALL get_se_param(se_kind_a,&
                          zeff=za,&
                          defined=defined,&
                          natorb=natorb_a,&
                          beta=beta_a,&
                          uss=uss,upp=upp,udd=udd,uff=uff)

        IF (.NOT.defined .OR. natorb_a < 1) CYCLE

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
             last_jatom = 0

             SELECT CASE (dft_control%qs_control%method)
               CASE ("EHT")
               CASE ("MNDO","AM1","PM3")
                 NULLIFY(h_blocka)
                 CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     block=h_blocka)
                 IF (.NOT.ASSOCIATED(h_blocka)) THEN
                   CALL stop_program(routine_name,module_name,__LINE__,&
                                    "A hmatrix block is missing",globenv)
                 END IF
             END SELECT

             DO inode=1,nnode
               CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                      neighbor=jatom,r=rij)

               IF (iatom == jatom) THEN

                 SELECT CASE (dft_control%qs_control%method)
                    CASE DEFAULT
                      CALL stop_program(routine_name,module_name,__LINE__,&
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

                 IF (jatom /= last_jatom) THEN
                   new_atom_b = .TRUE.
                   last_jatom = jatom
                 ELSE
                   new_atom_b = .FALSE.
                 END IF
                 IF (new_atom_b) THEN
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
                   IF (.NOT.ASSOCIATED(h_blocka)) THEN
                     CALL stop_program(routine_name,module_name,__LINE__,&
                                      "A hmatrix block is missing",globenv)
                   END IF
                 END IF

                 ! two-centre one-electron term
                 NULLIFY(s_block)
                 CALL get_block_node(matrix=matrix_s(1)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=s_block)
                 IF (.NOT.ASSOCIATED(s_block)) THEN
                   CALL stop_program(routine_name,module_name,__LINE__,&
                                     "A overlap matrix block is missing",&
                                     globenv)
                 END IF
                 DO i=1,SIZE(h_block,1)
                   DO j=1,SIZE(h_block,2)
                     h_block(i,j)=h_block(i,j)+kh*(ha(i)+hb(j))*s_block(i,j)
                   END DO
                 END DO
                 SELECT CASE (dft_control%qs_control%method)
                   CASE ("EHT")
                   CASE ("MNDO","AM1","PM3")
                   CALL rotint (se_kind_a,se_kind_b,rij,wint,kr,e1b,e2a,&
                                enuc,itype)
                   ! one-centre one-electron terms
                   NULLIFY(h_blockb)
                   CALL get_block_node(matrix=matrix_h(1)%matrix,&
                                       block_row=jatom,&
                                       block_col=jatom,&
                                       block=h_blockb)
                   IF (.NOT.ASSOCIATED(h_blockb)) THEN
                     CALL stop_program(routine_name,module_name,__LINE__,&
                                      "A hmatrix block is missing",globenv)
                   END IF
                   ! contribution to iatom
                   k=0
                   DO i=1,SIZE(h_blocka,1)
                     DO j=1,i-1
                       k=k+1
                       h_blocka(j,i)=h_blocka(j,i)+e1b(k)
                       h_blocka(i,j)=h_blocka(i,j)+e1b(k)
                     END DO
                     k=k+1
                     h_blocka(i,i)=h_blocka(i,i)+e1b(k)
                   END DO
                   ! contribution to jatom
                   k=0
                   DO i=1,SIZE(h_blockb,1)
                     DO j=1,i-1
                       k=k+1
                       h_blockb(j,i)=h_blockb(j,i)+e2a(k)
                       h_blockb(i,j)=h_blockb(i,j)+e2a(k)
                     END DO
                     k=k+1
                     h_blockb(i,i)=h_blockb(i,i)+e2a(k)
                   END DO

                 END SELECT

               END IF

               sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO ! inode => jatom(atom B)

           END DO ! ilist => iatom(atom A)

        END DO ! jkind

     END DO ! ikind

     CALL set_qs_env(qs_env=qs_env,&
                     matrix_h=matrix_h)

!deb
     CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,globenv,scale=evolt)
!deb

!deb
     WRITE(6,*) "FINISH HCORE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
!deb

  END SUBROUTINE build_se_core_matrix

! *****************************************************************************

END MODULE qs_semi_empirical_hamiltonian

! *****************************************************************************
