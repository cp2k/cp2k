!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Distribution of the spin orbit integral matrix. 
!> \par History
!> \author VW (27.02.2009)
! *****************************************************************************
MODULE qs_spin_orbit
  USE ai_spin_orbit,                   ONLY: pso
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_dbcsr_operations,             ONLY: dbcsr_obj_p_type
  USE cp_dbcsr_output,                 ONLY: dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_spin_orbit'

! *** Public subroutines ***

  PUBLIC :: build_pso_matrix

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the paramagnetic spin orbit matrix over
!>          Cartesian Gaussian functions.
!> \author  VW
!> \date    27.02.2009
!> \version 1.0
! *****************************************************************************

  SUBROUTINE build_pso_matrix(qs_env,matrix_so,rc,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_so
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_pso_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, inode, irow, iset, &
      istat, iw, jatom, jkind, jset, last_jatom, ldai, maxco, maxlgto, &
      maxsgf, natom, ncoa, ncob, neighbor_list_id, nkind, nlist, nnode, &
      nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: found, new_atom_b, tr
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rr_work, soab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: soint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)

    NULLIFY(cell,sab_orb,atomic_kind_set,particle_set,para_env)
    NULLIFY(logger)

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    para_env=para_env,&
                    sab_orb=sab_orb,&
                    cell=cell,&
                    distribution_2d=distribution_2d,&
                    error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)


!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxsgf=maxsgf)

    ldai = ncoset(maxlgto+1)
    CALL init_orbital_pointers(ldai)

    ALLOCATE(rr_work(0:2*maxlgto+2,ldai,ldai),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ff", (2*maxlgto+2)*dp_size)

    ALLOCATE (soab(maxco,maxco,3),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "soab",maxco*maxco*3*dp_size)

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",maxco*maxsgf*dp_size)

    ALLOCATE (soint(3),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "soint",3*int_size)

    rr_work(:,:,:) = 0.0_dp
    soab(:,:,:) = 0.0_dp
    work(:,:) = 0.0_dp

    DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)

       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=basis_set_a,&
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

          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_b)

          IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE

          CALL get_gto_basis_set(gto_basis_set=basis_set_b,&
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

          CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,nlist=nlist)

          neighbor_list => first_list(neighbor_list_set)

          DO ilist=1,nlist

             CALL get_neighbor_list(neighbor_list=neighbor_list,atom=iatom,nnode=nnode)

             ra = pbc(particle_set(iatom)%r,cell)

             last_jatom = 0

             neighbor_node => first_node(neighbor_list)

             DO inode=1,nnode

                CALL get_neighbor_node(neighbor_node=neighbor_node,neighbor=jatom,r=rab(:))

                rb = rab + ra
                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)
                rac = pbc(ra,rc,cell)
                rbc = rac - rab

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

                   DO i=1,3
                      NULLIFY(soint(i)%block)
                      CALL dbcsr_get_block_p(matrix=matrix_so(i)%matrix,&
                           row=irow,col=icol,BLOCK=soint(i)%block,tr=tr,found=found)
                   ENDDO
                END IF

                DO iset=1,nseta

                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)

                   DO jset=1,nsetb

                      IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)

                      ! *** Calculate the primitive fermi contact integrals ***

                      CALL pso(la_max(iset),la_min(iset),npgfa(iset),&
                               rpgfa(:,iset),zeta(:,iset),&
                               lb_max(jset),lb_min(jset),npgfb(jset),&
                               rpgfb(:,jset),zetb(:,jset),&
                               rac,rbc,rab,soab,SIZE(rr_work,1),SIZE(rr_work,2),rr_work,error)

                      ! *** Contraction step ***

                      DO i=1,3

                         CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                    1.0_dp,soab(1,1,i),SIZE(soab,1),&
                                    sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))

                         IF (iatom <= jatom) THEN

                            CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                       1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                       work(1,1),SIZE(work,1),&
                                       1.0_dp,soint(i)%block(sgfa,sgfb),&
                                       SIZE(soint(i)%block,1))

                         ELSE

                            CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                      -1.0_dp,work(1,1),SIZE(work,1),&
                                       sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                       1.0_dp,soint(i)%block(sgfb,sgfa),&
                                       SIZE(soint(i)%block,1))
                         ENDIF

                      ENDDO

                   ENDDO

                ENDDO

                neighbor_node => next(neighbor_node)

             ENDDO

             neighbor_list => next(neighbor_list)

          ENDDO

       ENDDO

    ENDDO

    ! *** Release work storage ***

    DEALLOCATE (soab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "soab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

    DEALLOCATE (soint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "soint")

!   *** Print the spin orbit matrix, if requested ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/PSO",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/PSO",&
            extension=".Log",error=error)
       CALL dbcsr_write_sparse_matrix(matrix_so(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL dbcsr_write_sparse_matrix(matrix_so(2)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL dbcsr_write_sparse_matrix(matrix_so(3)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/PSO", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_pso_matrix

! *****************************************************************************

END MODULE qs_spin_orbit

