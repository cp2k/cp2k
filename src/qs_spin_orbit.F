!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
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
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_so
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_pso_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, iw, &
      jatom, jkind, jset, last_jatom, ldai, maxco, maxlgto, maxsgf, natom, &
      ncoa, ncob, neighbor_list_id, nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: found, new_atom_b
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
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
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

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "basis_set_list",nkind)
    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       IF(inode==1) last_jatom = 0

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
             CALL cp_dbcsr_get_block_p(matrix=matrix_so(i)%matrix,&
                  row=irow,col=icol,BLOCK=soint(i)%block,found=found)
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

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    ! *** Release work storage ***

    DEALLOCATE (basis_set_list,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "basis_set_list")

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
       CALL cp_dbcsr_write_sparse_matrix(matrix_so(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_so(2)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_so(3)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/PSO", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_pso_matrix

! *****************************************************************************

END MODULE qs_spin_orbit

