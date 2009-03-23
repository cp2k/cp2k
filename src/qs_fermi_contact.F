!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Distribution of the Fermi contact integral matrix. 
!> \par History
!> \author VW (27.02.2009)
! *****************************************************************************
MODULE qs_fermi_contact

  USE ai_fermi_contact,                ONLY: fermi_contact
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
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
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fermi_contact'

! *** Public subroutines ***

  PUBLIC :: build_fermi_contact_matrix

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the Fermi contact matrix over Cartesian 
!>          Gaussian functions.
!> \author  VW
!> \date    27.02.2009
!> \version 1.0
! *****************************************************************************

  SUBROUTINE build_fermi_contact_matrix(qs_env,matrix_fc,rc,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_fc
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fermi_contact_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iab, iatom, icol, ikind, ilist, inode, irow, iset, &
      istat, iw, jatom, jkind, jset, last_jatom, ldai, ldfc, maxco, maxlgto, &
      maxsgf, natom, ncoa, ncob, nkind, nlist, nnode, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: new_atom_b
    REAL(KIND=dp)                            :: dab, dac, dbc, max_radius_a, &
                                                max_radius_b, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fcab, work
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: fcint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
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
                    para_env=para_env,&
                    sab_orb=sab_orb,&
                    cell=cell,&
                    error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ! *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxsgf=maxsgf)

    ldai = ncoset(maxlgto)
    CALL init_orbital_pointers(ldai)

    ldfc = maxco
    ALLOCATE (fcab(ldfc,ldfc),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "fcab",ldfc*ldfc*dp_size)
    fcab(:,:) = 0.0_dp

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",maxco*maxsgf*dp_size)
    work(:,:) = 0.0_dp

    ALLOCATE (fcint(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "fcint",1*int_size)
    NULLIFY(fcint(1)%block)

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

       max_radius_a = MAXVAL(set_radius_a(1:nseta))

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

          max_radius_b = MAXVAL(set_radius_b(1:nsetb))

          iab = ikind + nkind*(jkind - 1)

          IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

          neighbor_list_set => sab_orb(iab)%neighbor_list_set

          CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,nlist=nlist)

          neighbor_list => first_list(neighbor_list_set)

          DO ilist=1,nlist

             CALL get_neighbor_list(neighbor_list=neighbor_list,atom=iatom,nnode=nnode)

             ra = pbc(particle_set(iatom)%r,cell)

             rac = rc - ra
             dac = SQRT(rac(1)**2 + rac(2)**2 + rac(3)**2)

             last_jatom = 0

             neighbor_node => first_node(neighbor_list)

             DO inode=1,nnode

                CALL get_neighbor_node(neighbor_node=neighbor_node,neighbor=jatom,r=rab(:))

                rb = rab + ra
                rbc = rc - rb
                dbc = SQRT(rbc(1)**2 + rbc(2)**2 + rbc(3)**2)

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

                   NULLIFY (fcint(1)%block)
                   CALL get_block_node(matrix=matrix_fc(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=fcint(1)%block)
                ENDIF

                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)

                DO iset=1,nseta

                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)

                   DO jset=1,nsetb

                      IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)

                      ! *** Calculate the primitive fermi contact integrals ***

                      CALL fermi_contact(la_max(iset),la_min(iset),npgfa(iset),&
                                         rpgfa(:,iset),zeta(:,iset),&
                                         lb_max(jset),lb_min(jset),npgfb(jset),&
                                         rpgfb(:,jset),zetb(:,jset),&
                                         rac,rbc,dab,fcab,SIZE(fcab,1))

                      ! *** Contraction step ***

                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_dp,fcab(1,1),SIZE(fcab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))

                      IF (iatom <= jatom) THEN

                         CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    work(1,1),SIZE(work,1),&
                                    1.0_dp,fcint(1)%block(sgfa,sgfb),&
                                    SIZE(fcint(1)%block,1))

                      ELSE

                         CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                    1.0_dp,work(1,1),SIZE(work,1),&
                                    sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    1.0_dp,fcint(1)%block(sgfb,sgfa),&
                                    SIZE(fcint(1)%block,1))

                      ENDIF

                   ENDDO

                ENDDO

                neighbor_node => next(neighbor_node)

             ENDDO

             neighbor_list => next(neighbor_list)

          END DO

       END DO

    END DO

    ! *** Release work storage ***

    DEALLOCATE (fcab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "fcab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

    NULLIFY (fcint(1)%block)
    DEALLOCATE (fcint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "fcint")

!   *** Print the Fermi contact matrix, if requested ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/FERMI_CONTACT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/FERMI_CONTACT",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_fc(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/FERMI_CONTACT", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_fermi_contact_matrix

! *****************************************************************************

END MODULE qs_fermi_contact

