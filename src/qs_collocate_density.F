!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2002 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_collocate_density [1.0] *
!!
!!   NAME
!!     qs_collocate_density
!!
!!   FUNCTION
!!     Calculate the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!     1) Joost VandeVondele (01.2002)
!!
!!   MODIFICATION HISTORY
!!     - rewrote collocate for increased accuracy and speed
!!     - introduced the PGI hack for increased speed with that compiler
!!       (22.02.02)
!!     - Added Multiple Grid feature
!!     - new way to get over the grid (01.03.02)
!!     - removed timing calls since they were getting expensive
!!     - Updated with the new QS data structures (09.04.02,MK)
!!     - introduction of the real space grid type ( prelim. version JVdV 05.02)
!!     - parallel FFT (JGH 22.05.02)
!!     - multigrid arrays independent from density (JGH 30.08.02)
!!     - old density stored in g space (JGH 30.08.02)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_collocate_density

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE atomic_kind_types,      ONLY: atomic_kind_type,&
                                    get_atomic_kind,&
                                    get_atomic_kind_set
  USE basis_set_types,        ONLY: get_gto_basis_set,&
                                    gto_basis_set_type
  USE coefficient_types,      ONLY: coeff_allocate, coeff_deallocate, &
                                    coeff_add,coeff_sumup,&
                                    coeff_copy,&
                                    coeff_transform_space,&
                                    coeff_type,&
                                    coeff_zero
  USE qs_blacs,               ONLY: blacs_get_element,cp_full_matrix_type, &
                                    cp_fm_get_info
  USE realspace_grid_types,   ONLY: realspace_grid_type, rs_grid_allocate,  &
                                    rs_grid_deallocate, &
                                    rs_grid_setup, rs_pw_transfer, rs_grid_zero
  USE cube_utils,             ONLY: cube_info_type,&
                                    return_cube
  USE dft_types,              ONLY: dft_control_type
  USE gaussian_gridlevels,    ONLY: gaussian_gridlevel,&
                                    gridlevel_info_type
  USE global_types,           ONLY: global_environment_type
  USE l_utils,                ONLY: l_info_type,&
                                    return_l_info
  USE mathconstants,          ONLY: pi,twopi
  USE message_passing,        ONLY: mp_sum
  USE orbital_pointers,       ONLY: coset,indco,ncoset,nso
  USE particle_types,         ONLY: particle_type
  USE pw_types,               ONLY: COMPLEXDATA1D,COMPLEXDATA3D,&
                                    REALDATA3D,REALSPACE,&
                                    RECIPROCALSPACE
  USE qs_environment_types,   ONLY: get_qs_env,&
                                    qs_environment_type
  USE qs_interactions,        ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types, ONLY: first_list,&
                                    first_node,&
                                    get_neighbor_list,&
                                    get_neighbor_node,&
                                    neighbor_list_set_type,&
                                    neighbor_list_set_p_type,&
                                    neighbor_list_type,&
                                    neighbor_node_type,&
                                    next
  USE simulation_cell,        ONLY: cell_type,&
                                    pbc
  USE sparse_matrix_types,    ONLY: add_block_node,&
                                    allocate_matrix,&
                                    copy_sparse_to_full_matrix,&
                                    deallocate_matrix,&
                                    get_block_node,&
                                    get_matrix_info,&
                                    real_matrix_p_type,&
                                    real_matrix_type
  USE termination,            ONLY: stop_memory,&
                                    stop_program
  USE timings,                ONLY: timeset,&
                                    timestop

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: calculate_rho_core,&
            calculate_rho_elec,&
            calculate_wavefunction

! *** Public functions ***

  PUBLIC :: calculate_total_rho

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_rho_core(rho_core,cube_info,l_info,total_rho,qs_env,globenv)

    TYPE(coeff_type), INTENT(OUT)             :: rho_core
    TYPE(cube_info_type), INTENT(IN)          :: cube_info
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(l_info_type), INTENT(IN)             :: l_info
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    REAL(wp), INTENT(OUT)                     :: total_rho

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    TYPE(dft_control_type), POINTER :: dft_control

    REAL(wp) :: alpha,eps_rho_rspace
    INTEGER  :: atom_a,handle,iatom,ikind,natom

    REAL(wp), DIMENSION(3)   :: ra
    REAL(wp), DIMENSION(1,1) :: pab

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    TYPE(realspace_grid_type)                     :: rs_rho
    INTEGER, DIMENSION(3)                         :: nsmax
    TYPE(coeff_type)                              :: rhoc_r

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho_core","I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set)

    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=100
    call rs_grid_setup(rs_rho,rho_core%pw%pw_grid,nsmax)
    call rs_grid_allocate(rs_rho)
    call rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      DO iatom=1,natom

        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)

!MK May be for GAPW the analytic FT is useful to avoid a to high
!MK cutoff due to a big exponent of the core charge distribution
!MK     CALL collocate_pgf_product_gspace(0,alpha,0,0,0.0_wp,0,ra,&
!MK                                       (/0.0_wp,0.0_wp,0.0_wp/),&
!MK                                       0.0_wp,-1.0_wp,pab,rho%pw)

        IF(rs_rho%group_head) &
        CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_wp,0,ra,&
                                          (/0.0_wp,0.0_wp,0.0_wp/),&
                                          0.0_wp,-1.0_wp,pab,rs_rho,&
                                          cube_info,l_info,eps_rho_rspace)

      END DO

    END DO

    CALL coeff_allocate(rhoc_r,rho_core%pw%pw_grid,REALDATA3D)
    rhoc_r%pw%in_space = REALSPACE

    call rs_pw_transfer(rs_rho,rhoc_r%pw,"FORWARD")
    call rs_grid_deallocate(rs_rho)

    total_rho = calculate_total_rho(rhoc_r)

    CALL coeff_transform_space(rhoc_r,rho_core)

    CALL coeff_deallocate(rhoc_r)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_rho_core

! *****************************************************************************

  SUBROUTINE calculate_rho_elec(p,p_old,rebuild,rho,rho_gspace,rho_gspace_old,&
                                mgrid_rspace,mgrid_gspace,&
                                gridlevel_info,cube_info,l_info,total_rho,&
                                qs_env,globenv)

    TYPE(coeff_type), INTENT(INOUT)             :: rho,rho_gspace
    TYPE(coeff_type), INTENT(INOUT)             :: rho_gspace_old
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(qs_environment_type), INTENT(IN)       :: qs_env
    TYPE(real_matrix_type), POINTER             :: p,p_old
    REAL(wp), INTENT(OUT)                       :: total_rho
    LOGICAL, INTENT(IN)                         :: rebuild

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_rho_elec (MODULE qs_collocate_density)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    REAL(wp) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ineighbor,ipgf,iset,ishell,istat,j,jatom,jco,jkind,jpgf,jset,&
                jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,ncob,&
                nkind,nneighbor,nseta,nsetb,sgfa,sgfb
    LOGICAL  :: skip_block

    REAL(wp), DIMENSION(3) :: ra,rab,rb

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pab,pab_old,work

    REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(wp), DIMENSION(:,:), POINTER :: p_block,p_old_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    INTEGER, DIMENSION(3)             :: nsmax
    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER  :: rs_rho

!   ---------------------------------------------------------------------------



    allocate(rs_rho(gridlevel_info%ngrid_levels),STAT=istat)
    if (istat.ne.0) then
       CALL stop_program("calculate_rho_elec","rs_rho")
    endif
    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=100

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO
    CALL timeset("calculate_rho_elec","I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,maxco),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*maxco*wp_size)

    ALLOCATE (work(maxco,maxsgf_set),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*maxsgf_set*wp_size)

    IF (.NOT.rebuild) THEN
      ALLOCATE (pab_old(maxco,maxco),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pab_old",maxco*maxco*wp_size)
    END IF


    nkind = SIZE(atomic_kind_set)

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

      DO jkind=ikind,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               kind_radius=kind_radius_b,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        ab = ikind + jkind*(jkind - 1)/2

        sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)

        DO WHILE (ASSOCIATED(sab_orb_neighbor_list))

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO WHILE (ASSOCIATED(sab_orb_neighbor_node))

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=p,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=p_block)
            ELSE
              CALL get_block_node(matrix=p,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=p_block)
            END IF

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

            IF (.NOT.ASSOCIATED(p_block)) THEN
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            skip_block = .TRUE.

            IF (rebuild) THEN
              loop1: DO j=1,SIZE(p_block,2)
                DO i=1,SIZE(p_block,1)
                  IF (ABS(p_block(i,j)) > 1.0E-14_wp) THEN
                    skip_block = .FALSE.
                    EXIT loop1
                  END IF
                END DO
              END DO loop1
              IF (skip_block) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            ELSE
              IF (iatom <= jatom) THEN
                CALL get_block_node(matrix=p_old,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=p_old_block)
              ELSE
                CALL get_block_node(matrix=p_old,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=p_old_block)
              END IF
              loop2: DO j=1,SIZE(p_block,2)
                DO i=1,SIZE(p_block,1)
                  IF (ABS(p_block(i,j) - p_old_block(i,j)) > 1.0E-14_wp) THEN
                    skip_block = .FALSE.
                    EXIT loop2
                  END IF
                END DO
              END DO loop2
              IF (skip_block) THEN
                sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                CYCLE
              END IF
            ENDIF

            rb(:) = ra(:) + rab(:)
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              IF (iatom == jatom) THEN
                first_setb = iset
              ELSE
                first_setb = 1
              END IF

              DO jset=first_setb,nsetb

                IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

!               *** Decontract the density matrix set block ***

                IF (iatom <= jatom) THEN
                  CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                             1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             p_block(sgfa,sgfb),SIZE(p_block,1),&
                             0.0_wp,work(1,1),SIZE(work,1))
                  CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_wp,pab(1,1),SIZE(pab,1))
                  IF (.NOT.rebuild) THEN
                    CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                               1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               p_old_block(sgfa,sgfb),SIZE(p_old_block,1),&
                               0.0_wp,work(1,1),SIZE(work,1))
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab_old(1,1),SIZE(pab_old,1))
                    DO jco=1,ncob
                      DO ico=1,ncoa
                        pab(ico,jco) = pab(ico,jco) - pab_old(ico,jco)
                      END DO
                    END DO
                  END IF
                ELSE
                  CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                             1.0_wp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             p_block(sgfb,sgfa),SIZE(p_block,1),&
                             0.0_wp,work(1,1),SIZE(work,1))
                  CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             0.0_wp,pab(1,1),SIZE(pab,1))
                  IF (.NOT.rebuild) THEN
                    CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                               1.0_wp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               p_old_block(sgfb,sgfa),SIZE(p_old_block,1),&
                               0.0_wp,work(1,1),SIZE(work,1))
                    CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               0.0_wp,pab_old(1,1),SIZE(pab_old,1))
                    DO jco=1,ncoa
                      DO ico=1,ncob
                        pab(ico,jco) = pab(ico,jco) - pab_old(ico,jco)
                      END DO
                    END DO
                  END IF
                END IF

                DO ipgf=1,npgfa(iset)

                  IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                  na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                  na2 = ipgf*ncoset(la_max(iset))

                  IF ((iatom == jatom).AND.(iset == jset)) THEN
                    first_pgfb = ipgf
                  ELSE
                    first_pgfb = 1
                  END IF

                  DO jpgf=first_pgfb,npgfb(jset)

                    IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                    nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                    nb2 = jpgf*ncoset(lb_max(jset))

                    IF ((iatom == jatom).AND.&
                        (iset == jset).AND.&
                        (ipgf == jpgf)) THEN
                      scale = 1.0_wp
                    ELSE
                      scale = 2.0_wp
                    END IF

                    zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                    igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

                    IF (iatom <= jatom) THEN
                      CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                        ra,rab,rab2,scale,pab(na1:na2,nb1:nb2),&
                        !    rho(igrid_level)%pw,cube_info(igrid_level),&
                        rs_rho(igrid_level),cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
                    ELSE
                      CALL collocate_pgf_product_rspace(&
                        lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        rb,-rab,rab2,scale,pab(nb1:nb2,na1:na2),&
                        !   rho(igrid_level)%pw,cube_info(igrid_level),&
                        rs_rho(igrid_level),cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
                    END IF

                  END DO

                END DO

              END DO

            END DO

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_pw_transfer(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw,"FORWARD")
       call rs_grid_deallocate(rs_rho(igrid_level))
    ENDDO
    deallocate(rs_rho)

!   *** Get the density on the auxillary basis cutoff grid ***

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),mgrid_gspace(igrid_level))
      CALL coeff_add(rho_gspace,mgrid_gspace(igrid_level),rho_gspace)
    END DO

!   *** Save real-space density ***
    IF (.NOT. rebuild) THEN
      CALL coeff_add(rho_gspace,rho_gspace_old,rho_gspace)
    END IF
    CALL coeff_copy(rho_gspace,rho_gspace_old)

!   *** The latter is needed for a consistent density in real and g space ***

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    IF (.NOT.rebuild) THEN
      DEALLOCATE (pab_old,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pab_old")
    END IF

    total_rho = calculate_total_rho(rho)
    CALL timestop(0.0_wp,handle)


  END SUBROUTINE calculate_rho_elec

! *****************************************************************************
! modified calculate_rho_elec, should write the wavefunction represented by 
! the vector eigenvector on the grid. Is not efficient, in the sence that 
! it's presumably dominated by the FFT and the rs->pw and back routines
! especially in parrallel things are not efficient. 
! Currently it's only meant to provide a way to plot an MO
! *****************************************************************************
  SUBROUTINE calculate_wavefunction(mo_vectors,ivector,rho,rho_gspace,&
                                    mgrid_rspace,mgrid_gspace,&
                                    gridlevel_info,cube_info,l_info,&
                                    qs_env,globenv)

    TYPE(cp_full_matrix_type), POINTER :: mo_vectors
    INTEGER :: ivector

    TYPE(coeff_type), INTENT(INOUT)             :: rho,rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(global_environment_type), INTENT(IN)   :: globenv
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(qs_environment_type), INTENT(IN)       :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_wavefunction (MODULE qs_collocate_density)"

!   *** Local variables ***
    REAL(wp), DIMENSION(:), POINTER   :: eigenvector   

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    REAL(wp) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ineighbor,ipgf,iset,ishell,istat,j,jatom,jco,jkind,jpgf,jset,&
                jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,ncob,&
                nkind,nneighbor,nseta,nsetb,sgfa,sgfb
    LOGICAL  :: skip_block

    REAL(wp), DIMENSION(3) :: ra,rab,rb

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pab,work

    REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(wp), DIMENSION(:,:), POINTER :: p_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    INTEGER, DIMENSION(3)             :: nsmax
    TYPE(REALSPACE_GRID_TYPE), DIMENSION(:), POINTER  :: rs_rho
    INTEGER :: offset,natom,nao,ncol_global

!   ---------------------------------------------------------------------------
    CALL timeset("calculate_wavefunction","I","",handle)

    CALL cp_fm_get_info(matrix=mo_vectors,nrow_global=nao,ncol_global=ncol_global) 

    allocate(eigenvector(nao))
    ! only some processors have non-zero elements. The sum of eigenvector on 
    ! all cpus is the real eigenvector. combined with the rs->pw this will give the
    ! real density
    eigenvector(:)=0.0_wp
    DO i=1,nao 
       CALL blacs_get_element(mo_vectors,i,ivector,eigenvector(i))
    ENDDO

    allocate(rs_rho(gridlevel_info%ngrid_levels),STAT=istat)
    if (istat.ne.0) then
       CALL stop_program("calculate_rho_elec","rs_rho")
    endif
    ! hard coded for the time being
    ! gives the minimum thickness of a slice (or related quantity)
    ! is not yet used in rs_grid_setup
    nsmax(:)=1000

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_grid_setup(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw%pw_grid, & 
                          nsmax)
       call rs_grid_allocate(rs_rho(igrid_level))
       call rs_grid_zero(rs_rho(igrid_level))
    ENDDO

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*wp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*wp_size)

    offset=0

    DO iatom=1,natom


      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               orb_basis_set=orb_basis_set)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      ra(:) = pbc(particle_set(iatom)%r,cell)
      rb(:) = 0.0_wp
      rab2  = 0.0_wp
      dab   = 0.0_wp


      DO iset=1,nseta


         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         do i=1,nsgfa(iset)
            work(i,1)=eigenvector(offset+i)
         enddo

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_wp,pab(1,1),SIZE(pab,1))

         DO ipgf=1,npgfa(iset)


            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))


            scale = 1.0_wp
            zetp = zeta(ipgf,iset)
            igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

            CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_wp,0,&
                        ra,rab,rab2,scale,pab(na1:na2,1:1),&
                        rs_rho(igrid_level),cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       call rs_pw_transfer(rs_rho(igrid_level),mgrid_rspace(igrid_level)%pw,"FORWARD")
       call rs_grid_deallocate(rs_rho(igrid_level))
    ENDDO

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),mgrid_gspace(igrid_level))
      CALL coeff_add(rho_gspace,mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    deallocate(rs_rho)
    deallocate(eigenvector)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_wavefunction
! *****************************************************************************

  FUNCTION calculate_total_rho(rho) RESULT(total_rho)

    TYPE(coeff_type), INTENT(IN), TARGET :: rho

    REAL(wp) :: total_rho

!   *** Local variables ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION calculate_total_rho (MODULE qs_collocate_density)"

!   *** Local variables ***

    INTEGER, DIMENSION(:), POINTER :: lb_grid

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho%pw%cc3d)) THEN
      lb_grid => rho%pw%pw_grid%bounds(1,:)
      total_rho = -rho%pw%pw_grid%vol*rho%pw%cc3d(lb_grid(1),&
                                                  lb_grid(2),&
                                                  lb_grid(3))
    ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
      total_rho = -rho%pw%pw_grid%dvol*SUM(rho%pw%cr3d)
    ELSE IF (ASSOCIATED(rho%pw%cc)) THEN
      IF ( rho%pw%pw_grid%have_g0 ) THEN
        total_rho = -rho%pw%pw_grid%vol*rho%pw%cc(1)
      ELSE
        total_rho = 0._wp
      END IF
    ELSE
      CALL stop_program(routine,"No density coefficients available")
    END IF
    IF (rho%pw%pw_grid%para%mode.ne.0) THEN
       CALL mp_sum(total_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_rho

! *****************************************************************************
!MK
!MK  SUBROUTINE collocate_pgf_product_gspace(la_max,zeta,la_min,&
!MK                                          lb_max,zetb,lb_min,&
!MK                                          ra,rab,rab2,scale,pab,pw)
!MK
!MK    TYPE(pw_type), TARGET, INTENT(INOUT) :: pw
!MK    REAL(wp), INTENT(IN)                 :: rab2,scale,zeta,zetb
!MK    INTEGER, INTENT(IN)                  :: la_max,la_min,lb_max,lb_min
!MK    REAL(wp), DIMENSION(3), INTENT(IN)   :: ra,rab
!MK
!MK    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab
!MK
!MK!   *** Local parameters ***
!MK
!MK    CHARACTER(LEN=*), PARAMETER :: routine =&
!MK      "SUBROUTINE collocate_pgf_product_gspace (MODULE qs_collocate_density)"
!MK
!MK!   *** Local variables ***
!MK
!MK    COMPLEX(wp) :: yz,z
!MK    REAL(wp)    :: f,fa,fb,pij,prefactor,rzetp,twozetp,zetp
!MK    INTEGER     :: ax,ay,az,bx,by,bz,handle,ico,i,ig,ig2,istat,jco,jg,kg,la,&
!MK                   lb,lb_cube_min,ub_cube_max
!MK
!MK    COMPLEX(wp), DIMENSION(3) :: phasefactor
!MK    REAL(wp), DIMENSION(3)    :: dg,expfactor,fap,fbp,rap,rbp,rp
!MK    INTEGER, DIMENSION(3)     :: lb_cube,ub_cube
!MK
!MK    COMPLEX(wp), DIMENSION(:), ALLOCATABLE :: rag,rbg
!MK    REAL(wp), DIMENSION(:), ALLOCATABLE    :: g
!MK
!MK    COMPLEX(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: cubeaxis
!MK
!MK!   ---------------------------------------------------------------------------
!MK
!MK    dg(:) = twopi/(pw%pw_grid%npts(:)*pw%pw_grid%dr(:))
!MK
!MK    zetp = zeta + zetb
!MK    rzetp = 1.0_wp/zetp
!MK    f = zetb*rzetp
!MK    rap(:) = f*rab(:)
!MK    rbp(:) = rap(:) - rab(:)
!MK    rp(:) = ra(:) + rap(:)
!MK    twozetp = 2.0_wp*zetp
!MK    fap(:) = twozetp*rap(:)
!MK    fbp(:) = twozetp*rbp(:)
!MK
!MK    prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
!MK    phasefactor(:) = EXP(CMPLX(0.0_wp,-rp(:)*dg(:),wp))
!MK    expfactor(:) = EXP(-0.25*rzetp*dg(:)*dg(:))
!MK
!MK    lb_cube(:) = pw%pw_grid%bounds(1,:)
!MK    ub_cube(:) = pw%pw_grid%bounds(2,:)
!MK
!MK    lb_cube_min = MINVAL(lb_cube(:))
!MK    ub_cube_max = MAXVAL(ub_cube(:))
!MK
!MK    ALLOCATE (cubeaxis(lb_cube_min:Ub_cube_max,3,0:la_max,0:lb_max),STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"cubeaxis",0)
!MK
!MK    ALLOCATE (g(lb_cube_min:ub_cube_max),STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"g",0)
!MK
!MK    ALLOCATE (rag(lb_cube_min:ub_cube_max),STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"rag",0)
!MK
!MK    ALLOCATE (rbg(lb_cube_min:ub_cube_max),STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"rbg",0)
!MK
!MK    DO i=1,3
!MK
!MK      DO ig=lb_cube(i),ub_cube(i)
!MK        ig2 = ig*ig
!MK        cubeaxis(ig,i,0,0) = expfactor(i)**ig2*phasefactor(i)**ig
!MK      END DO
!MK
!MK      IF (la_max > 0) THEN
!MK        DO ig=lb_cube(i),ub_cube(i)
!MK          g(ig) = REAL(ig,wp)*dg(i)
!MK          rag(ig) = CMPLX(fap(i),-g(ig),wp)
!MK          cubeaxis(ig,i,1,0) = rag(ig)*cubeaxis(ig,i,0,0)
!MK        END DO
!MK        DO la=2,la_max
!MK          fa = REAL(la-1,wp)*twozetp
!MK          DO ig=lb_cube(i),ub_cube(i)
!MK            cubeaxis(ig,i,la,0) = rag(ig)*cubeaxis(ig,i,la-1,0) +&
!MK                                  fa*cubeaxis(ig,i,la-2,0)
!MK          END DO
!MK        END DO
!MK        IF (lb_max > 0) THEN
!MK          fa = twozetp
!MK          DO ig=lb_cube(i),ub_cube(i)
!MK            rbg(ig) = CMPLX(fbp(i),-g(ig),wp)
!MK            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
!MK            cubeaxis(ig,i,1,1) = rbg(ig)*cubeaxis(ig,i,1,0) +&
!MK                                 fa*cubeaxis(ig,i,0,0)
!MK          END DO
!MK          DO lb=2,lb_max
!MK            fb = REAL(lb-1,wp)*twozetp
!MK            DO ig=lb_cube(i),ub_cube(i)
!MK              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
!MK                                    fb*cubeaxis(ig,i,0,lb-2)
!MK              cubeaxis(ig,i,1,lb) = rbg(ig)*cubeaxis(ig,i,1,lb-1) +&
!MK                                    fb*cubeaxis(ig,i,1,lb-2) +&
!MK                                    fa*cubeaxis(ig,i,0,lb-1)
!MK            END DO
!MK          END DO
!MK          DO la=2,la_max
!MK            fa = REAL(la,wp)*twozetp
!MK            DO lb=1,lb_max
!MK              fb = REAL(lb-1,wp)*twozetp
!MK              DO ig=lb_cube(i),ub_cube(i)
!MK                cubeaxis(ig,i,la,lb) = rbg(ig)*cubeaxis(ig,i,la,lb-1) +&
!MK                                       fb*cubeaxis(ig,i,la,lb-2) +&
!MK                                       fa*cubeaxis(ig,i,la-1,lb-1)
!MK              END DO
!MK            END DO
!MK          END DO
!MK        END IF
!MK      ELSE
!MK        IF (lb_max > 0) THEN
!MK          DO ig=lb_cube(i),ub_cube(i)
!MK            g(ig) = REAL(ig,wp)*dg(i)
!MK            rbg(ig) = CMPLX(fbp(i),-g(ig),wp)
!MK            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
!MK          END DO
!MK          DO lb=2,lb_max
!MK            fb = REAL(lb-1,wp)*twozetp
!MK            DO ig=lb_cube(i),ub_cube(i)
!MK              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
!MK                                    fb*cubeaxis(ig,i,0,lb-2)
!MK            END DO
!MK          END DO
!MK        END IF
!MK      END IF
!MK
!MK    END DO
!MK
!MK    DO la=0,la_max
!MK      DO lb=0,lb_max
!MK        IF (la + lb == 0) CYCLE
!MK        fa = (1.0_wp/twozetp)**(la + lb)
!MK        DO i=1,3
!MK          DO ig=lb_cube(i),ub_cube(i)
!MK            cubeaxis(ig,i,la,lb) = fa*cubeaxis(ig,i,la,lb)
!MK          END DO
!MK        END DO
!MK      END DO
!MK    END DO
!MK
!MK!   *** Add the current primitive Gaussian function product to grid ***
!MK
!MK    DO ico=ncoset(la_min-1)+1,ncoset(la_max)
!MK
!MK      ax = indco(1,ico)
!MK      ay = indco(2,ico)
!MK      az = indco(3,ico)
!MK
!MK      DO jco=ncoset(lb_min-1)+1,ncoset(lb_max)
!MK
!MK        pij = prefactor*pab(ico,jco)
!MK
!MK        IF (ABS(pij) < eps_rho_gspace) CYCLE
!MK
!MK        bx = indco(1,jco)
!MK        by = indco(2,jco)
!MK        bz = indco(3,jco)
!MK
!MK        DO kg=lb_cube(3),ub_cube(3)
!MK          z = pij*cubeaxis(kg,3,az,bz)
!MK          IF (ABS(z) < eps_rho_gspace) CYCLE
!MK          DO jg=lb_cube(2),ub_cube(2)
!MK            yz = cubeaxis(jg,2,ay,by)*z
!MK            IF (ABS(yz) < eps_rho_gspace) CYCLE
!MK            DO ig=lb_cube(1),ub_cube(1)
!MK              pw%cc3d(ig,jg,kg) = pw%cc3d(ig,jg,kg) + cubeaxis(ig,1,ax,bx)*yz
!MK            END DO
!MK          END DO
!MK        END DO
!MK
!MK      END DO
!MK
!MK    END DO
!MK
!MK    DEALLOCATE (cubeaxis,STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"cubeaxis")
!MK
!MK    DEALLOCATE (g,STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"g")
!MK
!MK    DEALLOCATE (rag,STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"rag")
!MK
!MK    DEALLOCATE (rbg,STAT=istat)
!MK    IF (istat /= 0) CALL stop_memory(routine,"rbg")
!MK
!MK  END SUBROUTINE collocate_pgf_product_gspace
!MK
! *****************************************************************************

  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,&
                                          rsgrid,cube_info,l_info,eps_rho_rspace)

    ! TYPE(pw_type), INTENT(INOUT), TARGET :: pw
    TYPE(realspace_grid_type), INTENT(INOUT), TARGET :: rsgrid
    TYPE(cube_info_type), INTENT(IN)     :: cube_info
    TYPE(l_info_type), INTENT(IN)        :: l_info
    REAL(wp), INTENT(IN)                 :: eps_rho_rspace,rab2,scale,zeta,zetb
    INTEGER, INTENT(IN)                  :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN)   :: ra,rab

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab

!   *** Local variables ***

    REAL(wp) :: f,prefactor,radius,zetp
    INTEGER  :: i,ico,ig,j,jco,jg,k,kg

    REAL(wp), DIMENSION(3) :: dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(wp), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER :: lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER :: coef_max,lxa,lxyz,l,gridbounds(2,3)
    INTEGER, pointer,dimension(:) :: ly_max,lz_max
    REAL(wp) :: a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,&
                rpg,s,pg

    REAL(wp), pointer :: dpx(:,:)
    REAL(wp), pointer :: pzyx(:), alpha(:,:)
    integer,  pointer :: ipzyx(:,:,:,:,:,:)
    REAL(wp), pointer :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    INTEGER, pointer, dimension(:)  :: sphere_bounds
    REAL(wp)  :: xa,xb,xap,xbp,polypart,cutoff
    INTEGER  :: start,length,offset,cmax
    integer  handle

!   ---------------------------------------------------------------------------


    coef_max=la_max+lb_max+1

    zetp      = zeta + zetb
    f         = zetb/zetp
    prefactor = scale*EXP(-zeta*f*rab2)
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)  
    rb(:)     = ra(:)+rab(:)
    cutoff    = 0.0_wp

    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp,&
                                    zetp,eps_rho_rspace,prefactor,cutoff)

    IF (radius .eq. 0.0_wp ) THEN
      RETURN
    END IF

!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    !dr(:) = pw%pw_grid%dr(:)
    !ng(:) = pw%pw_grid%npts(:)
    !grid => pw%cr3d(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min,la_max,lb_min,lb_max,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** position of the gaussian product

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),wp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    DO i=1,3
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         if (start+length.ge.ub_cube(i)) EXIT
         start=start+length+1
        END DO
    ENDDO

!   *** initialise the p terms and loop logic
    lxyz=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             pzyx(lxyz)=prefactor*pab(ico,jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

 
!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,wp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max
       zbp=1.0_wp
       DO lzb=0,lb_max
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO 
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,wp)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max
       ybp=1.0_wp
       DO lyb=0,lb_max
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

 
!   *** make the alpha matrix ***
    alpha(:,:)=0.0_wp
    lx=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       lx=lx+1
       binomial_k_lxa=1.0_wp
       a=1.0_wp
       DO k=0,lxa
        binomial_l_lxb=1.0_wp
        b=1.0_wp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*real(lxb-l,wp)/real(l+1,wp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*real(lxa-k,wp)/real(k+1,wp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO

    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,wp)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

!   *** do the loop over the grid
    CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax), &
                        polz(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
                        lxyz_max,coef_max, &
                        cmax,ly_max(1),lz_max(1),gridbounds(1,1),&
                        map(-cmax,1),sphere_bounds(1))

  END SUBROUTINE collocate_pgf_product_rspace

! *****************************************************************************

END MODULE qs_collocate_density
