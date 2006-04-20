!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/qs_operators_ao   [1.0] *
!!
!!   NAME
!!      qs_operators_ao
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MI (07.2005)
!!
!!   MODIFICATION HISTORY
!!     created 07.2005
!!
!!   SOURCE
!******************************************************************************

MODULE qs_operators_ao

  USE ai_angmom,                       ONLY: angmom
  USE ai_moments,                      ONLY: diffop,&
                                             moment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
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
       block_p_type, first_list, first_node, get_neighbor_list, &
       get_neighbor_list_set, get_neighbor_node, neighbor_list_set_p_type, &
       neighbor_list_type, neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             add_1d_block_node,&
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

  PUBLIC :: p_xyz_ao, rmc_x_p_xyz_ao, rRc_xyz_ao, set_up_op_sm

!****************************************************************************

CONTAINS

!***************************************************************************
!!****f* qs_operators_ao/p_xyz_ao [1.0] *
!!
!!   NAME
!!     p_xyz_op 
!!
!!   FUNCTION
!!     Calculation of the components of the dipole operator in the velocity form
!!     The elements of the  sparse matrices are the integrals in the
!!     basis functions
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - op_sm: matrix representation of the p operator
!!              calculated in terms of the contracted basis functions
!!     - qs_env: enviroment for the lists and the basis sets
!!     - minimum_image : take into account only the first neighbors in the lists
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE p_xyz_ao(op_sm,qs_env,minimum_image,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_xyz_ao', &
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

    CALL timeset(routineN,"I"," ",handle)

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
         sab_all=sab_orb)

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
    ALLOCATE (op_dipt2(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_dipt(i,j)%block)
        NULLIFY (op_dipt2(i,j)%block)
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
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_dip,op_dip2,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,difabt,qs_env)&
!$OMP shared(op_dipt,op_dipt2)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2)
        ithread = 0
!$      ithread = omp_get_thread_num()

        difab => difabt(:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_dip => op_dipt(:,ithread)
        op_dip2 => op_dipt2(:,ithread)

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
!              IF (iatom <= jatom) THEN
!                irow = iatom
!                icol = jatom
!              ELSE
!                irow = jatom
!                icol = iatom
!              END IF
!dbg
!      write(*,*) 'atoms ', iatom, jatom
!dbg
              DO i = 1,3
                 NULLIFY(op_dip(i)%block)
                 irow = iatom
                 icol = jatom
                 CALL get_block_node(matrix=op_sm(i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip(i)%block)
!                 IF(iatom/=jatom) THEN
!                   NULLIFY(op_dip2(i)%block)
!                   irow = jatom
!                   icol = iatom
!                   CALL get_block_node(matrix=op_sm(i)%matrix,&
!                                     block_row=irow,&
!                                     block_col=icol,&
!                                     block=op_dip2(i)%block)
!                 END IF
              END DO
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

!dbg
!      write(*,*) 'sets ', iset, jset, la_max(iset), lb_max(jset)
!dbg
!            *** Calculate the primitive overlap integrals ***
                   CALL diffop(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                        zetb(:,jset),rpgfb(:,jset),lb_min(jset),rab,difab)
!dbg
!   write(*,*) 'primitives straight', iatom, jatom, ncoa,ncob
!   do i = 1,ncoa
!    write(*,'(i4,23f10.5)') i,( difab(i,j,1), j=1,ncob)
!   end do
!dbg
!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,difab(1,1,1),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
!                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(1)%block(sgfa,sgfb),&
                          SIZE(op_dip(1)%block,1))
!                   ELSE
!                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
!                          1.0_dp,work(1,1),ldwork,&
!                          sphi_a(1,sgfa),ldsa,&
!                          1.0_dp,op_dip(1)%block(sgfb,sgfa),&
!                          SIZE(op_dip(1)%block,1))
!                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,difab(1,1,2),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
!                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(2)%block(sgfa,sgfb),&
                          SIZE(op_dip(2)%block,1))
!                   ELSE
!                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
!                          1.0_dp,work(1,1),ldwork,&
!                          sphi_a(1,sgfa),ldsa,&
!                          1.0_dp,op_dip(2)%block(sgfb,sgfa),&
!                          SIZE(op_dip(2)%block,1))
!                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                        1.0_dp,difab(1,1,3),ldab,sphi_b(1,sgfb),ldsb,&
                        0.0_dp,work(1,1),ldwork)
!                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(3)%block(sgfa,sgfb),&
                          SIZE(op_dip(3)%block,1))
!                   ELSE
!                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
!                          1.0_dp,work(1,1),ldwork,&
!                          sphi_a(1,sgfa),ldsa,&
!                          1.0_dp,op_dip(3)%block(sgfb,sgfa),&
!                          SIZE(op_dip(3)%block,1))
!                   END IF
!dbg
!                IF(iatom == jatom) GOTO 100
!            *** Calculate the primitive overlap integrals ***
!                   CALL diffop(lb_max(jset),npgfb(jset),zetb(:,jset),&
!                        rpgfb(:,jset),lb_min(jset), la_max(iset),npgfa(iset),&
!                        zeta(:,iset),rpgfa(:,iset),la_min(iset),-rab,difab)
!            *** Contraction ***
!                   CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
!                        1.0_dp,difab(1,1,1),ldab,sphi_a(1,sgfa),ldsa,&
!                        0.0_dp,work(1,1),ldwork)
!                   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
!                        1.0_dp,sphi_b(1,sgfb),ldsb,&
!                        work(1,1),ldwork,&
!                        1.0_dp,op_dip2(1)%block(sgfb,sgfa),&
!                        SIZE(op_dip2(1)%block,1))
!            *** Contraction ***
!                   CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
!                        1.0_dp,difab(1,1,2),ldab,sphi_a(1,sgfa),ldsa,&
!                        0.0_dp,work(1,1),ldwork)
!                   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
!                        1.0_dp,sphi_b(1,sgfb),ldsb,&
!                        work(1,1),ldwork,&
!                        1.0_dp,op_dip2(2)%block(sgfb,sgfa),&
!                        SIZE(op_dip2(2)%block,1))
!!            *** Contraction ***
!                   CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
!                        1.0_dp,difab(1,1,3),ldab,sphi_a(1,sgfa),ldsa,&
!                        0.0_dp,work(1,1),ldwork)
!                   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
!                        1.0_dp,sphi_b(1,sgfb),ldsb,&
!                        work(1,1),ldwork,&
!                        1.0_dp,op_dip2(3)%block(sgfb,sgfa),&
!                        SIZE(op_dip2(3)%block,1))

!dbg
!   write(*,*) 'primitives reverse', jatom, iatom, ncob,ncoa
!   do i = 1,ncob
!    write(*,'(i4,23f10.5)') i, (difab(i,j,1), j=1,ncoa)
!   end do
!dbg
! 100 CONTINUE
!dbg
                END IF  !  >= dab

              END DO ! jset

              nrow = nrow + ncoa

            END DO ! iset

!dbg
!   write(*,*) 'contracted straight', iatom, jatom, ldsa, ldsb
!   do i = 1,size(op_dip(1)%block,1)
!    write(*,'(i4,23f10.5)') i, (op_dip(1)%block(i,j), j=1,size(op_dip(1)%block,2))
!   end do
!dbg
!dbg
! IF(jatom/=iatom) THEN
!   write(*,*) 'contracted reverse', jatom, iatom, ldsb,ldsa
!   do i = 1,size(op_dip2(1)%block,1)
!    write(*,'(i4,23f10.5)') i, (op_dip2(1)%block(i,j), j=1,size(op_dip2(1)%block,2))
!   end do
! END IF
!dbg
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind
    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_dipt(i,j)%block)
           NULLIFY(op_dipt2(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_dipt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(op_dipt2, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(difabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
 
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE p_xyz_ao


!***************************************************************************
!!****f* qs_operators_ao/rmc_x_p_xyz_ao [1.0] *
!!
!!   NAME
!!     rmc_x_p_xyz_op 
!!
!!   FUNCTION
!!     Calculation of the components of the (r-c)xp operator 
!!     where p is in the velocity form and c is the center of the primitive
!!     function multiplied from the left.
!!     The resulting sparse matrix is not symmetric, therefore all_couples
!!     lists are used, which cntain the couples ab and ba
!!     The center c will be the position of
!!     iatom and the position of jatom respectively.
!!     The elements of the  sparse matrices are the integrals in the
!!     basis functions
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - op_sm: matrix representation of the (r-c)xp operator 
!!              in terms of the contracted basis functions
!!     - qs_env: enviroment for the lists and the basis sets
!!     - minimum_image : take into account only the first neighbors in the lists
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE rmc_x_p_xyz_ao(op_sm,qs_env,minimum_image,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

    CALL timeset(routineN,"I"," ",handle)

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
         sab_all=sab_orb)

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
!$OMP shared(my_minimum_image,Lxo2)
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

            IF (new_atom_b) THEN
            ! iatom , jatom and  jatom, iatom
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


            nrow = 0
            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!            *** Calculate the primitive overlap integrals (iatom jatom)***
                  rc(1:3) = ra(1:3)
                  rac(1:3) = 0.0_dp
                  rbc(1:3) = rab(1:3)
!                  rac(1:3) = ra(1:3)
!                  rbc(1:3) = rb(1:3)
                  CALL angmom(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                        zetb(:,jset),rpgfb(:,jset),lb_min(jset),rac,rbc,angab)
!dbg
!   write(*,*) 'primitives straight', iatom, jatom, ncoa,ncob
!    write(*,*) 'set ', iset, jset, la_max(iset), lb_max(jset)
!   do i = 1,ncoa
!    write(*,'(i4,23f10.5)') i,(angab(i,j,1), j=1,ncob)
!   end do
!dbg
!   
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

              nrow = nrow + ncoa

            END DO ! iset
!dbg
!   write(*,*) 'contracted straight', iatom, jatom, ldsa, ldsb
!   do i = 1,size(op_angij(1)%block,1)
!    write(*,'(i4,23f10.5)') i, (op_angij(1)%block(i,j), j=1,size(op_angij(1)%block,2))
!   end do
!dbg


            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A
!$OMP end parallel

      END DO ! jkind

    END DO ! ikind
 !dbg
! stop 'integrals'
!dbg
    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_angijt(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_angijt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(angabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE rmc_x_p_xyz_ao

!***************************************************************************
!!****f* qs_operators_ao/rRc_xyz_ao [1.0] *
!!
!!   NAME
!!     rRc_xyz_op 
!!
!!   FUNCTION
!!     Calculation of the components of the dipole operator in the length form
!!     by taking the relative position operator r-Rc, with respect a reference point Rc
!!     Probably it does not work for PBC, or maybe yes if the wfn are 
!!     sufficiently localized
!!     The elements of the  sparse matrices are the integrals in the
!!     basis functions
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - op_sm: matrix representation of the p operator
!!              calculated in terms of the contracted basis functions
!!     - qs_env: enviroment for the lists and the basis sets
!!     - rc: reference vector position
!!     - order: maximum order of the momentum, for the doipole order = 1 
!!     - minimum_image : take into account only the first neighbors in the lists
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE rRc_xyz_ao(op_sm,qs_env,rc,order,minimum_image,soft,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: Rc(3)
    INTEGER, INTENT(IN)                      :: order
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image, soft
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rRc_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, imom, inode, irow, iset, &
      istat, ithread, j, jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, &
      ldwork, M_dim, maxl, ncoa, ncob, nkind, nlist, nnode, nrow, nseta, nsetb, &
      nthread, sgfa, sgfb
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

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (mab,mabt,op_dip,op_dipt,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell=cell,particle_set=particle_set,sab_orb=sab_orb)

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
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_dip,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,mabt,qs_env)&
!$OMP shared(op_dipt)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm)&
!$OMP shared(atomic_kind_set,particle_set,cell)&
!$OMP shared(my_minimum_image,Lxo2)
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
 
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE   rRc_xyz_ao

!******************************************************************************
!!****f* qs_operators_ao/set_up_op_sm
!!
!!   NAME
!!     set_up_op_sm
!!
!!   FUNCTION
!!     allocate a not symetric sparse matrix
!!
!!   NOTE
!!     I am not sure this can really work
!!     How is treated a non symmetric sparse matrix? How is it distributed
!!     This routine can work only if the element ij and ji are on the same processor
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE set_up_op_sm(matrix,qs_env,symmetry,name,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(LEN=*)                         :: symmetry, name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_up_op_sm', &
      routineP = moduleN//':'//routineN

    INTEGER :: iab, iatom, icol, ikind, ilist, inode, irow, istat, jatom, jkind, &
      last_jatom, natom, neighbor_list_id, nkind, nlist, nnode, nsgf
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

!******************************************************************************

    failure = .FALSE.

    NULLIFY(particle_set,distribution_2d,atomic_kind_set,sab_orb)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,&
         distribution_2d=distribution_2d,&
         neighbor_list_id=neighbor_list_id)

    IF(symmetry=="none") THEN
      CALL get_qs_env(qs_env=qs_env,&
           sab_all=sab_orb)
    ELSE
      CALL get_qs_env(qs_env=qs_env,&
           sab_all=sab_orb)
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
                         distribution_2d=distribution_2d)

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
                                  block=block)
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
                                  block=block)
            END IF
!            CALL add_1d_block_node(matrix=matrix,&
!                                block_row=iatom,&
!                                block_col=jatom,&
!                                block=block)
!
!            IF(iatom/=jatom) THEN
!              NULLIFY(block)
!              CALL add_1d_block_node(matrix=matrix,&
!                                  block_row=jatom,&
!                                  block_col=iatom,&
!                                  block=block)
!            END IF
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO  ! inode

        END DO  ! ilist
      END DO  ! jkind
    END DO  ! ikind

  END SUBROUTINE set_up_op_sm

END MODULE qs_operators_ao
