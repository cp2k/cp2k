!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_moments [1.0] *
!!
!!   NAME
!!     qs_moments
!!
!!   FUNCTION
!!     Calculates the moment integrals <a|r^m|b>
!!
!!   AUTHOR
!!     JGH (20.07.2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_moments

  USE ai_moments,                      ONLY: moment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_moments'

! *** Public subroutines ***

  PUBLIC :: build_moment_matrix

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_moment_matrix(qs_env,moments,nmoments,ref_point,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: moments
    INTEGER, INTENT(IN)                      :: nmoments
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: ref_point
    TYPE(cp_error_type), INTENT(INOUT)  &
                                             :: error

    TYPE block_p_type
      REAL(KIND = dp), DIMENSION(:,:), POINTER :: block
    END TYPE block_p_type

    CHARACTER(LEN=*), PARAMETER :: routineN = "build_moment_matrix", &
      routineP = moduleN//':'//routineN

    REAL(KIND = dp), DIMENSION(3)                         :: rab, rac, rbc, ra, rb, rc
    REAL(KIND = dp), DIMENSION(:), POINTER                :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER                        :: la_max,la_min,lb_max,&
                                                             lb_min,npgfa,npgfb,&
                                                             nsgfa,nsgfb
    REAL(KIND = dp), DIMENSION(:,:), POINTER              :: rpgfa,rpgfb,sphi_a,&
                                                             sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER                      :: first_sgfa,first_sgfb

    TYPE(block_p_type), DIMENSION(:), ALLOCATABLE         :: mint
    REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE          :: work
    REAL(KIND = dp), DIMENSION(:,:,:), ALLOCATABLE        :: mab

    INTEGER                                               :: nm, i, nkind, &
       natom, ikind, jkind,ilist, nlist, inode, nnode, iatom, jatom, iab, &
       last_jatom, nseta, nsetb, iset, jset, ncoa, ncob, sgfa, sgfb, istat, &
       maxsgf, maxco, icol, irow
    TYPE(cell_type), POINTER                              :: cell
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER                       :: atomic_kind
    TYPE(neighbor_list_set_type), POINTER                 :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER                     :: neighbor_list
    TYPE(neighbor_node_type), POINTER                     :: neighbor_node
    TYPE(gto_basis_set_type), POINTER                     :: orb_basis_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: matrix_s
    INTEGER                                               :: handle
    REAL(KIND = dp)                                       :: rab2, dab
    LOGICAL                                               :: failure, new_atom_b

! ------------------------------------------------------------------------------

    IF ( nmoments < 1 ) RETURN

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    NULLIFY (matrix_s)
    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s,error=error)

    nm = (6 + 11*nmoments + 6*nmoments**2 + nmoments**3)/6 - 1
    CALL allocate_matrix_set ( moments, nm ,error=error)

    DO i=1,nm
      CALL replicate_matrix(source=matrix_s(1)%matrix,&
                            TARGET=moments(i)%matrix,&
                            target_name="Moments",error=error)
    END DO

    NULLIFY (atomic_kind_set,particle_set,sab_orb,cell)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell,&
                    sab_orb=sab_orb,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf)

    ALLOCATE (mab(maxco,maxco,nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mab(:,:,:) = 0.0_dp

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(:,:) = 0.0_dp

    ALLOCATE (mint(nm),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
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

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        neighbor_list_set => sab_orb(iab)%neighbor_list_set

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          last_jatom = 0

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=jatom)

            ! this guarentees minimum image convention
            ! anything else would not make sense
            IF (jatom == last_jatom) CYCLE

            last_jatom = jatom

            IF (iatom <= jatom) THEN
              irow = iatom
              icol = jatom
            ELSE
              irow = jatom
              icol = iatom
            END IF

            DO i=1,nm
              NULLIFY (mint(i)%block)
              CALL get_block_node(matrix=moments(i)%matrix,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  BLOCK=mint(i)%block)
              mint(i)%block = 0._dp
            END DO

            ra(:) = pbc(particle_set(iatom)%r(:),cell)
            rb(:) = pbc(particle_set(jatom)%r(:),cell)
            rc(:) = 0.5_dp * (ref_point(:,iatom) + ref_point(:,jatom))
            rab(:) = ra(:) - rb(:)
            rab(:) = pbc(rab,cell)
            rac(:) = ra(:) - rc(:)
            rac(:) = pbc(rac,cell)
            rbc(:) = rb(:) - rc(:)
            rbc(:) = pbc(rbc,cell)
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

!               *** Calculate the primitive integrals ***

                CALL moment(la_max(iset),npgfa(iset),zeta(:,iset),&
                            rpgfa(:,iset),la_min(iset),&
                            lb_max(jset),npgfb(jset),zetb(:,jset),&
                            rpgfb(:,jset),lb_min(jset),&
                            nmoments,rac,rbc,mab)

!               *** Contraction step ***

                DO i=1,nm

                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                             1.0_dp,mab(1,1,i),SIZE(mab,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_dp,work(1,1),SIZE(work,1))

                  IF (iatom <= jatom) THEN

                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               work(1,1),SIZE(work,1),&
                               1.0_dp,mint(i)%block(sgfa,sgfb),&
                               SIZE(mint(i)%block,1))

                  ELSE

                    CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                               1.0_dp,work(1,1),SIZE(work,1),&
                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               1.0_dp,mint(i)%block(sgfb,sgfa),&
                               SIZE(mint(i)%block,1))

                  END IF

                END DO


              END DO
            END DO

            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (mab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i=1,nm
      NULLIFY (mint(i)%block)
    END DO
    DEALLOCATE (mint,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_moment_matrix

!***************************************************************************

END MODULE qs_moments

! *****************************************************************************
