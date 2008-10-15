!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      created 07.2005
!> \author MI (07.2005)
! *****************************************************************************
MODULE qs_operators_ao

  USE ai_angmom,                       ONLY: angmom
  USE ai_moments,                      ONLY: diff_momop,&
                                             diffop,&
                                             moment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_operators_ao'

! *** Public subroutines ***

  PUBLIC :: p_xyz_ao, rmc_x_p_xyz_ao, rRc_xyz_ao, rRc_xyz_der_ao, set_up_op_sm

CONTAINS

! *****************************************************************************
!> \brief Calculation of the components of the dipole operator in the velocity form
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \param op_sm matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param minimum_image take into account only the first neighbors in the lists
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE p_xyz_ao(op_sm,qs_env,minimum_image,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, inode, irow, iset, &
      istat, ithread, j, jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, &
      ldwork, maxl, ncoa, ncob, nkind, nlist, nnode, nseta, nsetb, nthread, &
      sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, my_minimum_image, &
                                                new_atom_b
    REAL(KIND=dp)                            :: alpha, dab, Lxo2, Lyo2, Lzo2, &
                                                rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: difab, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: difabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip, op_dip2
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_dipt, op_dipt2
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (difab,difabt,op_dip,op_dipt,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell=cell,particle_set=particle_set,&
         sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF

    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()

    ALLOCATE(difabt(ldab,ldab,3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    difabt(1:ldab,1:ldab,1:3,0:nthread-1) = 0.0_dp
    ALLOCATE(workt(ldwork,ldwork,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    workt(1:ldwork,1:ldwork,0:nthread-1) = 0.0_dp
    ALLOCATE (op_dipt(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_dipt(i,j)%block)
      END DO
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

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
      ldsa = SIZE(sphi_a,1)
      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

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
        ldsb = SIZE(sphi_b,1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,difab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_dip,op_dip2,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local,alpha)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,difabt,qs_env)&
!$OMP shared(op_dipt,op_dipt2)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2,Lyo2,Lzo2)
        ithread = 0
!$      ithread = omp_get_thread_num()

        difab => difabt(:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_dip => op_dipt(:,ithread)

!$OMP do
        DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0
          alpha = 1.0_dp

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

            IF(my_minimum_image) THEN
              IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            END IF

            rb = rab + ra

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
                alpha = 1.0_dp
              ELSE
                irow = jatom
                icol = iatom
                IF(op_sm(1)%matrix%symmetry=="antisymmetric") THEN
                  alpha = -1.0_dp
                END IF
              END IF

              DO i = 1,3
                 NULLIFY(op_dip(i)%block)
                 CALL get_block_node(matrix=op_sm(i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip(i)%block)
              END DO
            END IF  ! new_atom_b
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!            *** Calculate the primitive overlap integrals ***
                   CALL diffop(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                        zetb(:,jset),rpgfb(:,jset),lb_min(jset),rab,difab)

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        alpha,difab(1,1,1),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(1)%block(sgfa,sgfb),&
                          SIZE(op_dip(1)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(1)%block(sgfb,sgfa),&
                          SIZE(op_dip(1)%block,1))
                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        alpha,difab(1,1,2),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(2)%block(sgfa,sgfb),&
                          SIZE(op_dip(2)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(2)%block(sgfb,sgfa),&
                          SIZE(op_dip(2)%block,1))
                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        alpha,difab(1,1,3),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(3)%block(sgfa,sgfb),&
                          SIZE(op_dip(3)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(3)%block(sgfb,sgfa),&
                         SIZE(op_dip(3)%block,1))
                  END IF
                END IF  !  >= dab

              END DO ! jset

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind

    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_dipt(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_dipt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(difabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE p_xyz_ao

! *****************************************************************************
!> \brief Calculation of the components of the (r-c)xp operator
!>      where p is the gradient and c is the center of the primitive
!>      function multiplied from the .
!>      The resulting sparse matrix is not symmetric, therefore all_pair
!>      lists are used, which contain the either ab and ba
!>      The center c will be the position of
!>      iatom and the position of jatom respectively.
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \param op_sm matrix representation of the (r-c)xp operator
!>               in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param minimum_image take into account only the first neighbors in the lists
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE rmc_x_p_xyz_ao(op_sm,qs_env,minimum_image,wancen,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: wancen
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rmc_x_p_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, inode, irow, iset, &
      istat, ithread, j, jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, &
      ldwork, maxl, ncoa, ncob, nkind, nlist, nnode, nrow, nseta, nsetb, &
      nthread, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, my_minimum_image, &
                                                new_atom_b
    REAL(KIND=dp)                            :: dab, Lxo2, Lyo2, Lzo2, rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc, rc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: angab, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: angabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_angij
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_angijt
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (angab,angabt,op_angij,op_angijt,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell=cell,particle_set=particle_set,&
         sab_all=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF

    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()

    ALLOCATE (angabt(ldab,ldab,3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    angabt(1:ldab,1:ldab,1:3,0:nthread-1)=0.0_dp
    ALLOCATE (workt(ldwork,ldwork,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    workt(1:ldwork,1:ldwork,0:nthread-1)=0.0_dp
    ALLOCATE (op_angijt(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_angijt(i,j)%block)
      END DO
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

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
      ldsa = SIZE(sphi_a,1)
      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

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
        ldsb = SIZE(sphi_b,1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,angab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_angij,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP private(rc,rac,rbc)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,angabt,qs_env)&
!$OMP shared(op_angijt)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2,Lyo2,Lzo2)&
!$OMP shared(wancen)
        ithread = 0
!$      ithread = omp_get_thread_num()

        angab => angabt(:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_angij => op_angijt(:,ithread)

!$OMP do
       DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

            IF(my_minimum_image) THEN
              IF(ABS(rab(1)) > ABS(Lxo2) .OR. ABS(rab(2)) > ABS(Lyo2) .OR. ABS(rab(3)) > ABS(Lzo2)) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            END IF

            rb = rab + ra

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

!dbg
!   write(*,*) 'atoms ', iatom, jatom, rab(1:3)
!dbg
            IF (new_atom_b) THEN
            ! iatom , jatom
              DO i = 1,3
                 NULLIFY(op_angij(i)%block)
                 irow = iatom
                 icol = jatom
                 CALL get_block_node(matrix=op_sm(i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_angij(i)%block)
              END DO
            END IF  !new_atom_b
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!            *** Calculate the primitive overlap integrals (iatom jatom)***
                   IF(PRESENT(wancen)) THEN
                      rc = wancen
                      rac = pbc(rc,ra,cell)
                      rbc = rac + rab
                   ELSE
                      rc(1:3) = rb(1:3)
                      rac(1:3) = -rab(1:3)
                      rbc(1:3) = 0.0_dp
                   END IF

!                  rc(1:3) = 0.5_dp
!                  rac = pbc(rc,ra,cell)
!!                  rac(1:3) = ra(1:3)-0.5_dp
!                  rbc = pbc(rc,rb,cell)
!!                  rbc(1:3) = rb(1:3)-0.5_dp

                  CALL angmom(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                        zetb(:,jset),rpgfb(:,jset),lb_min(jset),rac,rbc,angab)
!            *** Contraction ***
                  DO i = 1,3
                     CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                          1.0_dp,angab(1,1,i),ldab,sphi_b(1,sgfb),ldsb,&
                          0.0_dp,work(1,1),ldwork)
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_angij(i)%block(sgfa,sgfb),&
                          SIZE(op_angij(i)%block,1))
                  END DO

                END IF  !  >= dab

              END DO ! jset

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A
!$OMP end parallel

      END DO ! jkind

    END DO ! ikind

    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_angijt(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_angijt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(angabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE rmc_x_p_xyz_ao

! *****************************************************************************
!> \brief Calculation of the components of the dipole operator in the length form
!>      by taking the relative position operator r-Rc, with respect a reference point Rc
!>      Probably it does not work for PBC, or maybe yes if the wfn are
!>      sufficiently localized
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \param op_sm matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param rc reference vector position
!> \param order maximum order of the momentum, for the doipole order = 1
!> \param minimum_image take into account only the first neighbors in the lists
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE rRc_xyz_ao(op_sm,qs_env,rc,order,minimum_image,soft,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: Rc(3)
    INTEGER, INTENT(IN)                      :: order
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image, soft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rRc_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, imom, inode, irow, &
      iset, istat, ithread, j, jatom, jkind, jset, last_jatom, ldab, ldsa, &
      ldsb, ldwork, M_dim, maxl, ncoa, ncob, nkind, nlist, nnode, nrow, &
      nseta, nsetb, nthread, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, my_minimum_image, &
                                                my_soft, new_atom_b
    REAL(KIND=dp)                            :: dab, Lxo2, Lyo2, Lzo2, rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mab, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: mabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_dipt
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (mab,mabt,op_dip,op_dipt,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell=cell,particle_set=particle_set,sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF
    my_soft = .FALSE.
    IF(PRESENT(soft)) THEN
      my_soft = soft
    END IF

    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()

    M_dim = ncoset(order)-1
    CPPostcondition(M_dim<=SIZE(op_sm,1),cp_failure_level,routineP,error,failure)

    ALLOCATE(mabt(ldab,ldab,M_dim,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mabt(1:ldab,1:ldab,1:M_dim,0:nthread-1) = 0.0_dp
    ALLOCATE(workt(ldwork,ldwork,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    workt(1:ldwork,1:ldwork,0:nthread-1) = 0.0_dp
    ALLOCATE (op_dipt(M_dim,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imom = 1, M_dim
      DO j=0,nthread-1
        NULLIFY (op_dipt(imom,j)%block)
      END DO
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb=my_soft,&
                           orb_basis_set=orb_basis_set)

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
      ldsa = SIZE(sphi_a,1)
      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             softb = my_soft,&
                             orb_basis_set=orb_basis_set)

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
        ldsb = SIZE(sphi_b,1)

        iab = ikind + nkind*(jkind - 1)
        ldsb = SIZE(sphi_b,1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,mab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb,rac,rbc)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_dip,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local,imom)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,mabt,qs_env)&
!$OMP shared(op_dipt,order,M_dim,rc)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2,Lyo2,Lzo2)
        ithread = 0
!$      ithread = omp_get_thread_num()

        mab => mabt(:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_dip => op_dipt(:,ithread)

!$OMP do
        DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

            IF(my_minimum_image) THEN
              IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            END IF

            rb = rab + ra

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

              DO imom = 1, M_dim
                 NULLIFY(op_dip(imom)%block)
                 CALL get_block_node(matrix=op_sm(imom)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip(imom)%block)
              END DO   ! imom
            END IF  ! new_atom_b

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            nrow = 0
            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

                  rac = pbc(rc,ra,cell)
                  rbc = pbc(rc,rb,cell)

!            *** Calculate the primitive overlap integrals ***
                  CALL moment(la_max(iset),npgfa(iset),zeta(:,iset),&
                    rpgfa(:,iset),la_min(iset),&
                    lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),&
                    lb_min(jset), order,rac,rbc,mab)

                  DO imom = 1,M_dim
!                 *** Contraction ***
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,mab(1,1,imom),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                    IF(iatom<=jatom) THEN
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(imom)%block(sgfa,sgfb),&
                          SIZE(op_dip(imom)%block,1))
                    ELSE
                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(imom)%block(sgfb,sgfa),&
                          SIZE(op_dip(imom)%block,1))
                    END IF

                  END DO  ! imom
                END IF  !  >= dab

              END DO ! jset

              nrow = nrow + ncoa

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind
    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_dipt(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_dipt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(mabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE   rRc_xyz_ao

! *****************************************************************************
!> \brief Calculation of the  multipole operators integrals
!>      and of its derivatives of the type
!>      [\mu | op | d(\nu)/dR(\nu)]-[d(\mu)/dR(\mu)| op | \nu]
!>      by taking the relative position operator r-Rc, with respect a reference point Rc
!>      The derivative are with respect to the primitive position,
!>      The multipole operator is symmetric and if it does not depend on R(\mu) or R(\nu)
!>      therefore  [\mu | op | d(\nu)/dR(\nu)] = -[d(\mu)/dR(\mu)| op | \nu]
!>        [\mu|op|d(\nu)/dR]-[d(\mu)/dR|op|\nu]=2[\mu|op|d(\nu)/dR]
!>      When it is not the case a correction term is needed
!> 
!>     The momentum operator [\mu|M|\nu] is symmetric, the number of components is
!>     detrermined by the order: 3 for order 1 (x,y,x), 9 for order 2(xx,xy,xz,yy,yz,zz)
!>     The derivative of the type [\mu | op | d(\nu)/dR_i(\nu)], where
!>     i indicates the cartesian direction, is antisymmetric only when
!>     the no component M =(r_i) or (r_i r_j) is in the same cartesian
!>     direction of the derivative,  indeed
!>   d([\mu|M|\nu])/dr_i = [d(\mu)/dr_i|M|\nu] + [\mu|M|d(\nu)/dr_i] + [\mu |d(M)/dr_i|\nu]
!>   d([\mu|M|\nu])/dr_i = -[d(\mu)/dR_i(\mu)|M|\nu] -[\mu|M|d(\nu)/dR_i(\nu)] + [\mu |d(M)/dr_i|\nu]
!>     Therefore we cannot use an antisymmetric matrix
!> 
!>     The same holds for the derivative with respect to the electronic position r
!>     taking into account that [\mu|op|d(\nu)/dR] = -[\mu|op|d(\nu)/dr]
!> \param op_sm matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param rc reference vector position
!> \param order maximum order of the momentum, for the dipole order = 1
!> \param minimum_image take into account only the first neighbors in the lists
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      Probably it does not work for PBC, or maybe yes if the wfn are
!>      sufficiently localized
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \par History
!>      03.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE rRc_xyz_der_ao(op_sm,op_sm_der,qs_env,rc,order,minimum_image,soft,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_sm_der
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: Rc(3)
    INTEGER, INTENT(IN)                      :: order
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image, soft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rRc_xyz_der_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, idir, ikind, ilist, imom, inode, &
      ipgf, irow, iset, istat, ithread, j, jatom, jkind, jpgf, jset, &
      last_jatom, lda_min, ldab, ldb_min, ldsa, ldsb, ldwork, M_dim, maxl, &
      na, nb, ncoa, ncob, nda, ndb, nkind, nlist, nnode, nseta, nsetb, &
      nthread, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, my_minimum_image, &
                                                my_soft, new_atom_b
    REAL(KIND=dp)                            :: alpha, alpha_der, dab, Lxo2, &
                                                Lyo2, Lzo2, rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mab, mab_tmp, workt
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: difmab, mabt
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: difmabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_dip_der, op_dipt
    TYPE(block_p_type), DIMENSION(:, :, :), &
      POINTER                                :: op_dip_dert
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!dbg

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(op_sm),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(op_sm_der),cp_failure_level,routineP,error,failure)
    IF(.NOT.op_sm_der(1,1)%matrix%symmetry=="none") THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (difmab,difmabt,mab,mab_tmp,mabt)
    NULLIFY (op_dip,op_dip_der,op_dipt,op_dip_dert,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell=cell,particle_set=particle_set,&
         sab_all=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF
    my_soft = .FALSE.
    IF(PRESENT(soft)) THEN
      my_soft = soft
    END IF

    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()

    M_dim = ncoset(order)-1
    CPPostcondition(M_dim<=SIZE(op_sm,1),cp_failure_level,routineP,error,failure)

    ALLOCATE(mabt(ldab,ldab,M_dim,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mabt(1:ldab,1:ldab,1:M_dim,0:nthread-1) = 0.0_dp
    ALLOCATE(difmabt(ldab,ldab,M_dim,3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    difmabt(1:ldab,1:ldab,1:M_dim,1:3,0:nthread-1) = 0.0_dp

    ALLOCATE(workt(ldwork,ldwork,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    workt(1:ldwork,1:ldwork,0:nthread-1) = 0.0_dp
    ALLOCATE (op_dipt(M_dim,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_dip_dert(M_dim,3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imom = 1, M_dim
      DO j=0,nthread-1
        NULLIFY (op_dipt(imom,j)%block)
        DO i =1,3
          NULLIFY (op_dip_dert(imom,i,j)%block)
        END DO
      END DO
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb=my_soft,&
                           orb_basis_set=orb_basis_set)

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
      ldsa = SIZE(sphi_a,1)
      !
      ! protect if no soft basis set
      IF(ldsa.EQ.0) CYCLE
      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             softb = my_soft,&
                             orb_basis_set=orb_basis_set)

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
        ldsb = SIZE(sphi_b,1)
        !
        ! protect if no soft basis set
        IF(ldsb.EQ.0) CYCLE
        
        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,mab,difmab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb,rac,rbc)&
!$OMP private(iset,jset,ncoa,ncob,nnode,inode,irow,icol,imom,idir)&
!$OMP private(lda_min,mab_tmp,na,nda,nb,ndb,ipgf,jpgf,i,j)&
!$OMP private(op_dip,op_dip_der,rab,rab2,dab)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local,alpha,failure,ldb_min)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,mabt,difmabt,qs_env)&
!$OMP shared(op_dipt,op_dip_dert,order,M_dim,rc)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm,op_sm_der)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2,Lyo2,Lzo2,alpha_der,error)
        ithread = 0
!$      ithread = omp_get_thread_num()

        mab => mabt(:,:,:,ithread)
        difmab => difmabt(:,:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_dip => op_dipt(:,ithread)
        op_dip_der => op_dip_dert(:,:,ithread)

!$OMP do
        DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

            IF(my_minimum_image) THEN
              IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            END IF

            rb = rab + ra

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            IF (new_atom_b) THEN
                irow = iatom
                icol = jatom
                alpha_der = 2.0_dp

              DO imom = 1, M_dim
                 NULLIFY(op_dip(imom)%block)
                 CALL get_block_node(matrix=op_sm(imom)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip(imom)%block)
                 DO idir = 1,3
                   NULLIFY(op_dip_der(imom,idir)%block)
                   CALL get_block_node(matrix=op_sm_der(imom,idir)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip_der(imom,idir)%block)
                 END DO  ! idir
              END DO   ! imom
            END IF  ! new_atom_b

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

                   rac = pbc(rc,ra,cell)
                   rbc = rac + rab
!                  rac = pbc(rc,ra,cell)
!                  rbc = pbc(rc,rb,cell)

                  ALLOCATE(mab_tmp(npgfa(iset)*ncoset(la_max(iset)+1),&
                           npgfb(jset)*ncoset(lb_max(jset)+1),ncoset(order)-1),STAT=istat)
                  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

                  lda_min = MAX ( 0, la_min(iset)-1)
                  ldb_min = MAX ( 0, lb_min(jset)-1)
!            *** Calculate the primitive overlap integrals ***
                  CALL moment(la_max(iset)+1,npgfa(iset),zeta(:,iset),&
                       rpgfa(:,iset),lda_min,&
                       lb_max(jset)+1,npgfb(jset),zetb(:,jset),rpgfb(:,jset),&
                       ldb_min, order,rac,rbc,mab_tmp)

!            *** Calculate the derivatives
                  CALL diff_momop(la_max(iset),npgfa(iset),zeta(:,iset),&
                       rpgfa(:,iset),la_min(iset),lb_max(jset),npgfb(jset),&
                       zetb(:,jset),rpgfb(:,jset),lb_min(jset),order,rac,rbc,&
                       difmab,mab_ext=mab_tmp)

! Contract and copy in the sparse matrix
                  mab =0.0_dp
                  DO imom = 1,M_dim
                    na = 0
                    nda = 0
                    DO ipgf=1,npgfa(iset)
                      nb =0
                      ndb = 0
                      DO jpgf = 1,npgfb(jset)
                        DO j = 1,ncoset(lb_max(jset))
                          DO i = 1,ncoset(la_max(iset))
                            mab(i+na,j+nb,imom) = mab_tmp(i+nda,j+ndb,imom)
                          END DO  ! i
                        END DO  ! j
                        nb = nb + ncoset(lb_max(jset))
                        ndb = ndb + ncoset(lb_max(jset)+1)
                      END DO  ! jpgf
                      na = na + ncoset(la_max(iset))
                      nda = nda + ncoset(la_max(iset)+1)
                    END DO  ! ipgf

!                 *** Contraction ***
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,mab(1,1,imom),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(imom)%block(sgfa,sgfb),&
                          SIZE(op_dip(imom)%block,1))
!dbg
!   num1 =500+10*imom
!   write(num1,*) ' imom ', imom,' iset ' , iset, nsgfa(iset), sgfa, ' jset ', jset,  nsgfb(jset), sgfb
!   do i = sgfa,sgfa+nsgfa(iset)-1
!      write(num1,'(I4,50f12.6)') i, op_dip(imom)%block(i,sgfb:sgfb+nsgfb(jset)-1)
!   end do
!dbg

                    alpha = -1.0_dp!-alpha_der
                    DO idir = 1,3
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                          alpha,difmab(1,1,imom,idir),ldab,sphi_b(1,sgfb),ldsb,&
                          0.0_dp,work(1,1),ldwork)
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                            1.0_dp,sphi_a(1,sgfa),ldsa,&
                            work(1,1),ldwork,&
                            1.0_dp,op_dip_der(imom,idir)%block(sgfa,sgfb),&
                            SIZE(op_dip_der(imom,idir)%block,1))
!dbg
!   num1 =1000+100*imom+10*idir+2
!   write(num1,*) ' imom ', imom, ' der. ', idir, ' iset '  , iset, nsgfa(iset), &
!        sgfa, ' jset ', jset,  nsgfb(jset), sgfb
!   do i = sgfa,sgfa+nsgfa(iset)-1
!      write(num1,'(I4,50f12.6)') i, op_dip_der(imom,idir)%block(i,sgfb:sgfb+nsgfb(jset)-1)
!   end do
!dbg

                    END DO  ! idir

                  END DO  ! imom

                  DEALLOCATE(mab_tmp,STAT=istat)
                  CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
                END IF  !  >= dab

              END DO ! jset

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind
    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_dipt(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_dipt,op_dip_dert, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(mabt,difmabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE   rRc_xyz_der_ao

! *****************************************************************************
!> \brief allocate a not symetric sparse matrix
!> \note
!>      I am not sure this can really work
!>      How is treated a non symmetric sparse matrix? How is it distributed
!>      This routine can work only if the element ij and ji are on the same processor
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE set_up_op_sm(matrix,qs_env,symmetry,name,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(LEN=*)                         :: symmetry, name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_up_op_sm', &
      routineP = moduleN//':'//routineN

    INTEGER :: iab, iatom, icol, ikind, ilist, inode, irow, istat, jatom, &
      jkind, last_jatom, natom, neighbor_list_id, nkind, nlist, nnode, nsgf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    LOGICAL                                  :: failure, new_atom_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.

    NULLIFY(particle_set,distribution_2d,atomic_kind_set,sab_orb)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,&
         distribution_2d=distribution_2d,&
         neighbor_list_id=neighbor_list_id,error=error)

    IF(symmetry=="none") THEN
      CALL get_qs_env(qs_env=qs_env,&
           sab_all=sab_orb,error=error)
    ELSE
      CALL get_qs_env(qs_env=qs_env,&
           sab_orb=sab_orb,error=error)
    END IF

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    CALL allocate_matrix(matrix=matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=first_sgf,&
                         last_row=last_sgf,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name=name,&
                         matrix_symmetry=symmetry,&
                         sparsity_id=neighbor_list_id, &
                         distribution_2d=distribution_2d,error=error)

    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

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
                                   neighbor=jatom)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            NULLIFY(block)
            IF(symmetry=="none") THEN
              CALL add_block_node(matrix=matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=block,error=error)
            ELSE
              IF(iatom<=jatom) THEN
                irow = iatom
                icol = jatom
              ELSE
                irow = jatom
                icol = iatom
              END IF
              CALL add_block_node(matrix=matrix,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  block=block,error=error)
            END IF
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO  ! inode

        END DO  ! ilist
      END DO  ! jkind
    END DO  ! ikind

  END SUBROUTINE set_up_op_sm

END MODULE qs_operators_ao
