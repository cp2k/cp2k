!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003 CP2K developers group                           !
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
!!     - distributed real space code (JGH 17.07.03)
!!     - refactoring and new loop ordering (JGH 23.11.03)
!!     - OpenMP parallelization (JGH 03.12.03)
!!     - Modified to compute tau (Joost 12.03)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_collocate_density
! *****************************************************************************
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE cp_control_types,                       ONLY: dft_control_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dp,&
                                             dp_size
  USE l_utils,                         ONLY: l_info_type,&
                                             return_l_info
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_sum,&
                                             mp_sync
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset,&
                                             nso
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pool_init_coeff,&
                                             pw_pool_give_back_coeff
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_p_type,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             real_matrix_type,&
                                             allocate_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='qs_collocate_density'
! *** Public subroutines ***

  PUBLIC :: calculate_rho_core,&
            calculate_rho_elec,&
            calculate_wavefunction,&
            collocate_pgf_product_rspace

! *** Public functions ***

  PUBLIC :: calculate_total_rho,calculate_total_abs_rho

! *** Public type ***

  TYPE lgrid_type
     INTEGER :: ldim
     REAL(dp), DIMENSION(:), POINTER :: r
  END TYPE lgrid_type

  PUBLIC :: lgrid_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_rho_core(rho_core,cube_info,l_info,total_rho,qs_env,error)

    TYPE(coeff_type), INTENT(INOUT)                :: rho_core
    TYPE(cube_info_type), INTENT(IN)             :: cube_info
    TYPE(l_info_type), INTENT(IN)                :: l_info
    TYPE(qs_environment_type), INTENT(IN)        :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    REAL(KIND = dp), INTENT(OUT)                        :: total_rho

!   *** Local variables ***

    logical :: failure
    CHARACTER(len=*), PARAMETER :: routineN='calculate_rho_core',&
         routineP=moduleN//':'//routineN
    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: cell
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(cp_para_env_type), POINTER :: para_env

    REAL(KIND = dp) :: alpha,eps_rho_rspace
    INTEGER  :: atom_a,handle,iatom,ikind,natom,handled
    INTEGER  :: bo(2), npme, dir, j, ierr, ncurr

    REAL(KIND = dp), DIMENSION(3)   :: ra
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: pab
    INTEGER, DIMENSION ( :, : ), POINTER :: tasks
    INTEGER  :: ithread, omp_get_thread_num

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    TYPE(realspace_grid_type), POINTER            :: rs_rho
    TYPE(coeff_type)                              :: rhoc_r
    TYPE(pw_env_type), POINTER                    :: pw_env
    TYPE(cp_rs_pool_type), POINTER                :: auxbas_rs_pool
    TYPE(pw_pool_type), POINTER                   :: auxbas_pw_pool

!   ---------------------------------------------------------------------------
    CALL timeset("calculate_rho_core","I","",handle)

    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_pool,auxbas_pw_pool,tasks)
    ALLOCATE(pab(1,1),STAT=ierr)
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab", 1 )

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env)
    CALL pw_env_get(pw_env,auxbas_rs_pool=auxbas_rs_pool,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)

    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_rho)
    CALL rs_grid_zero(rs_rho)

    ncurr = -1

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      bo = get_limit ( natom, para_env%num_pe, para_env%mepos )

      npme = bo(2) - bo(1) + 1
      IF ( ncurr < npme ) THEN
        CALL reallocate ( tasks, 1, 2, 1, npme )
        ncurr = npme
      ELSE
        tasks = 0
      END IF

      dir = rs_rho%direction

!$OMP parallel do &
!$OMP default(none) &
!$OMP private(iatom,j,atom_a,ra) &
!$OMP shared(bo,tasks,dir,rs_rho,atom_list,particle_set,cell)
      DO iatom=bo(1),bo(2)

        j = iatom - bo(1) + 1
        tasks ( 1, j ) = iatom
        IF ( dir /= 0) THEN
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          tasks ( 2, j ) = FLOOR(ra(dir)/rs_rho%dr(dir))
          tasks ( 2, j ) = MODULO ( tasks ( 2, j ), rs_rho%npts(dir) )
          tasks ( 2, j ) = tasks ( 2, j ) + rs_rho%lb(dir)
        END IF

      END DO

      CALL rs_get_my_tasks ( rs_rho, tasks, npme )

      ithread=0
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,iatom,atom_a,ra,ithread) &
!$OMP shared(npme,tasks,atom_list,particle_set,cell,pab,rs_rho) &
!$OMP shared(cube_info,l_info,eps_rho_rspace,alpha)
      DO j=1,npme

        iatom = tasks(1,j)
        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)

!$      ithread=omp_get_thread_num()
        CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
             (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,-1.0_dp,pab,0,0,rs_rho,&
             cube_info,l_info,eps_rho_rspace,ithread=ithread)

      END DO
    
    END DO

    DEALLOCATE ( tasks, STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "tasks" )
    DEALLOCATE ( pab, STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab" )

    CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rhoc_r, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)

    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,"FORWARD")
    CALL rs_pool_give_back_rs(auxbas_rs_pool,rs_rho,error=error)

    total_rho = calculate_total_rho(rhoc_r)

    CALL coeff_transform_space(rhoc_r,rho_core)

    CALL pw_pool_give_back_coeff(auxbas_pw_pool, rhoc_r, error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_core

! *****************************************************************************
! both rho,rho_gspace are the new rho, also rho_gspace_old is overwritten with
! the new rho
! *****************************************************************************

  SUBROUTINE calculate_rho_elec(p,p_old,rebuild,rho,rho_gspace,rho_gspace_old,&
                                mgrid_rspace,mgrid_gspace,&
                                gridlevel_info,cube_info,l_info,total_rho,&
                                qs_env, compute_tau, error)

    TYPE(coeff_type), INTENT(INOUT)             :: rho,rho_gspace
    TYPE(coeff_type), INTENT(INOUT)             :: rho_gspace_old
    TYPE(coeff_type), DIMENSION(:), POINTER     :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(gridlevel_info_type), INTENT(IN)       :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)               :: l_info
    TYPE(qs_environment_type), INTENT(IN)       :: qs_env
    TYPE(real_matrix_type), POINTER             :: p,p_old
    TYPE(cp_error_type), INTENT(inout),OPTIONAL :: error
    REAL(KIND = dp), INTENT(OUT)                :: total_rho
    LOGICAL, INTENT(IN)                         :: rebuild
    LOGICAL, INTENT(IN), OPTIONAL               :: compute_tau

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
    TYPE(real_matrix_type), POINTER   :: deltap

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(realspace_grid_p_type), DIMENSION(:), POINTER    :: rs_rho
    TYPE(pw_env_type), POINTER                            :: pw_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), POINTER        :: rs_pools
    TYPE(cp_para_env_type), POINTER                       :: para_env

    REAL(KIND = dp) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ilist,ineighbor,inode,ipgf,iset,ishell,istat,j,jatom,jco,jkind,&
                jpgf,jset,jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,&
                ncob,nkind,nlist,nneighbor,nnode,nseta,nsetb,sgfa,sgfb,&
                ijatoms,ijsets,k
    INTEGER  :: handled,brow,bcol,pdim1,pdim2,n,dir,tp,itask,npme,nmax,maxsgf
    LOGICAL  :: reread

    REAL(KIND = dp), DIMENSION(3) :: ra,rab,rb,rp

    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: pab,work,pblock
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: dist_ab,pabt,workt
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: dab_local

    REAL(KIND = dp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: p_block,p_old_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb

    INTEGER, DIMENSION(:,:,:), POINTER  :: tasks
    INTEGER, DIMENSION(:,:), POINTER    :: tasks_local, ival, latom
    INTEGER, DIMENSION(:), POINTER      :: ntasks
    INTEGER, PARAMETER                  :: max_tasks = 2000, add_tasks = 1000
    INTEGER                             :: curr_tasks, natom_pairs
    INTEGER, DIMENSION(:,:), POINTER    :: asets, atasks
    INTEGER                             :: ithread, nthread, nalimit(2)
    INTEGER                             :: omp_get_thread_num, omp_get_max_threads
    INTEGER, DIMENSION(:), POINTER      :: ub, lb
    TYPE(lgrid_type)                    :: lgrid
    LOGICAL                             :: my_compute_tau
!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,deltap,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,pblock,dist_ab,dab_local,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,p_block,p_old_block,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         tasks,tasks_local,ival,latom,ntasks,asets,atasks,pabt,workt)

    ! by default, do not compute the kinetic energy density (tau)
    ! if compute_tau, all grids referening to rho are actually tau
    IF (PRESENT(compute_tau)) THEN 
       my_compute_tau = compute_tau
    ELSE
       my_compute_tau = .FALSE.
    ENDIF

    IF (my_compute_tau) THEN
       call timeset("calculate_rho_tau","I","",handle)
    ELSE
       call timeset("calculate_rho_elec","I","",handle)
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    para_env=para_env,&
                    pw_env=pw_env)
    CALL pw_env_get(pw_env, rs_pools=rs_pools,error=error)

    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf,&
                             maxsgf_set=maxsgf_set)

    IF ( nthread > 1 ) THEN
      n=0
      DO igrid_level = 1,gridlevel_info%ngrid_levels
        n = MAX(n,rs_rho(igrid_level)%rs_grid%ngpts_local)
      END DO
      n = n*nthread
      CALL reallocate(lgrid%r,1,n)
    END IF

    nkind = SIZE(atomic_kind_set)

    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(dist_ab,1,3,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(latom,1,2,1,max_tasks)
    CALL reallocate(dab_local,1,3,1,max_tasks)
    CALL reallocate(atasks,1,2,1,max_tasks)
    CALL reallocate(asets,1,2,1,max_tasks)
    curr_tasks = max_tasks

!   *** Initialize working density matrix ***

    CALL allocate_matrix(matrix=deltap,&
                         nrow=p%nrow,&
                         ncol=p%ncol,&
                         nblock_row=p%nblock_row,&
                         nblock_col=p%nblock_col,&
                         first_row=p%first_row(:),&
                         last_row=p%last_row(:),&
                         first_col=p%first_col(:),&
                         last_col=p%last_col(:),&
                         matrix_name="DeltaP",&
                         matrix_symmetry=p%symmetry)

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

        ab = ikind + nkind*(jkind - 1)

        IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN

           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(ab)%neighbor_list_set,&
                                      nlist=nlist)
           sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
        ELSE
           nlist=0
        END IF

        ntasks = 0
        tasks = 0

        DO ilist = 1, nlist

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom,nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO inode = 1, nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

             IF (iatom <= jatom) THEN
                brow = iatom
                bcol = jatom
             ELSE
                brow = jatom
                bcol = iatom
             END IF

             CALL get_block_node(matrix=p,&
                                 block_row=brow,&
                                 block_col=bcol,&
                                 BLOCK=p_block)

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

             IF (.NOT.ASSOCIATED(p_block)) THEN
               sab_orb_neighbor_node => next(sab_orb_neighbor_node)
               CYCLE
             END IF

             pdim1 = SIZE(p_block,1)
             pdim2 = SIZE(p_block,2)

             NULLIFY ( pblock )
             CALL add_block_node ( deltap, brow, bcol, pblock )

             IF (rebuild) THEN
               pblock = p_block
             ELSE
               CALL get_block_node(matrix=p_old,&
                                   block_row=brow,&
                                   block_col=bcol,&
                                   BLOCK=p_old_block)
               IF (.NOT.ASSOCIATED(p_old_block)) &
                  CALL stop_program(routine,"p_old_block not associated")
               pblock = p_block - p_old_block
             END IF
   
             IF ( ALL ( 100.0_dp*ABS(pblock) < eps_rho_rspace) ) THEN
               sab_orb_neighbor_node => next(sab_orb_neighbor_node)
               CYCLE
             END IF

             rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
             dab = SQRT(rab2)

             DO iset=1,nseta

               IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE
   
               IF (iatom == jatom) THEN
                 first_setb = iset
               ELSE
                 first_setb = 1
               END IF
 
               DO jset=first_setb,nsetb
 
                 IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
 
                 DO ipgf=1,npgfa(iset)

                   IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                   IF ((iatom == jatom).AND.(iset == jset)) THEN
                     first_pgfb = ipgf
                   ELSE
                     first_pgfb = 1
                   END IF

                   DO jpgf=first_pgfb,npgfb(jset)

                     IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                     zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                     igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

                     ntasks(igrid_level) = ntasks(igrid_level) + 1
                     n = ntasks(igrid_level)
                     IF ( n > curr_tasks ) THEN
                       curr_tasks = curr_tasks + add_tasks
                       CALL reallocate(tasks,1,8,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                       CALL reallocate(dist_ab,1,3,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                     END IF
                     dir = rs_rho(igrid_level)%rs_grid%direction
                     tasks (1,n,igrid_level) = n
                     IF ( dir /= 0) THEN
                       rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                       rp(:) = pbc(rp,cell)
                       tp = FLOOR(rp(dir)/rs_rho(igrid_level)%rs_grid%dr(dir))
                       tp = MODULO ( tp, rs_rho(igrid_level)%rs_grid%npts(dir) )
                       tasks (2,n,igrid_level) = tp + rs_rho(igrid_level)%rs_grid%lb(dir)
                     END IF
                     tasks (3,n,igrid_level) = iatom
                     tasks (4,n,igrid_level) = jatom
                     tasks (5,n,igrid_level) = iset
                     tasks (6,n,igrid_level) = jset
                     tasks (7,n,igrid_level) = ipgf
                     tasks (8,n,igrid_level) = jpgf
                     dist_ab (1,n,igrid_level) = rab(1)
                     dist_ab (2,n,igrid_level) = rab(2)
                     dist_ab (3,n,igrid_level) = rab(3)

                   END DO

                 END DO

               END DO

             END DO

             sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

        DO igrid_level = 1, gridlevel_info%ngrid_levels
          n = ntasks ( igrid_level )
          IF ( n > SIZE ( tasks_local, 2 ) ) &
            CALL reallocate(tasks_local,1,2,1,n)
          IF ( n > SIZE ( ival, 2 ) ) &
            CALL reallocate(ival,1,6,1,n)
          IF ( n > SIZE ( latom, 2 ) ) &
            CALL reallocate(latom,1,2,1,n)
          IF ( n > SIZE ( dab_local, 2 ) ) &
            CALL reallocate(dab_local,1,3,1,n)

!$OMP parallel do private(i)
          DO i=1,n
            tasks_local(1,i) = tasks(1,i,igrid_level)
            tasks_local(2,i) = tasks(2,i,igrid_level)
            latom(1,i) = tasks(3,i,igrid_level)
            latom(2,i) = tasks(4,i,igrid_level)
            ival(1,i) = tasks(3,i,igrid_level)
            ival(2,i) = tasks(4,i,igrid_level)
            ival(3,i) = tasks(5,i,igrid_level)
            ival(4,i) = tasks(6,i,igrid_level)
            ival(5,i) = tasks(7,i,igrid_level)
            ival(6,i) = tasks(8,i,igrid_level)
            dab_local(1,i) = dist_ab(1,i,igrid_level)
            dab_local(2,i) = dist_ab(2,i,igrid_level)
            dab_local(3,i) = dist_ab(3,i,igrid_level)
          END DO
!$OMP parallel do private(i)
          DO i=n+1,SIZE(tasks_local,2)
            tasks_local(1,i) = 0
            tasks_local(2,i) = 0
          END DO
          npme = 0
          CALL rs_get_my_tasks ( rs_rho(igrid_level)%rs_grid, tasks_local, &
               npme, ival=ival, rval=dab_local, pmat=deltap, pcor=latom )
          CALL rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )
 
          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            lb => rs_rho(igrid_level)%rs_grid%lb_local
            ub => rs_rho(igrid_level)%rs_grid%ub_local
            lgrid%ldim = rs_rho(igrid_level)%rs_grid%ngpts_local
!$OMP parallel private(ithread,n)
!$          ithread = omp_get_thread_num()
            n = ithread*lgrid%ldim + 1
            CALL dcopy(lgrid%ldim,0._dp,0,lgrid%r(n),1)
!$OMP end parallel
          END IF
!$OMP parallel &
!$OMP default(none) &
!$OMP private(ijatoms,ithread,itask,iatom,jatom,ra,brow,bcol,p_block) &
!$OMP private(ijsets,iset,jset,ncoa,ncob,sgfa,sgfb) &
!$OMP private(work,pab,istat,nalimit) &
!$OMP private(rb,rab,rab2,ipgf,jpgf,na1,na2,nb1,nb2,scale) &
!$OMP shared(maxco,maxsgf_set,natom_pairs,atasks,asets,ival,particle_set,cell) &
!$OMP shared(deltap,npgfa,npgfb,ncoset,la_max,lb_max,first_sgfa,first_sgfb) &
!$OMP shared(nsgfa,nsgfb,sphi_a,sphi_b,dab_local,la_min,lb_min,zeta,zetb) &
!$OMP shared(rs_rho,igrid_level,cube_info,l_info,eps_rho_rspace,lgrid,nthread) &
!$OMP shared(workt,pabt,my_compute_tau)
          ithread = 0
!$        ithread = omp_get_thread_num()
          nalimit = get_limit(natom_pairs,nthread,ithread)
          pab => pabt(:,:,ithread)
          work => workt(:,:,ithread)
!         DO ijatoms = nalimit(1),nalimit(2)
!$OMP do
          DO ijatoms = 1,natom_pairs
            itask = atasks(1,asets(1,ijatoms))
            iatom  = ival (1,itask)
            jatom  = ival (2,itask)
            ra(:) = pbc(particle_set(iatom)%r,cell)
            IF (iatom <= jatom) THEN
              brow = iatom
              bcol = jatom
            ELSE
              brow = jatom
              bcol = iatom
            END IF
            CALL get_block_node(matrix=deltap,&
                                block_row=brow,&
                                block_col=bcol,&
                                BLOCK=p_block)
            IF (.NOT.ASSOCIATED(p_block)) &
               CALL stop_program(routine,"p_block not associated in deltap")
            DO ijsets = asets(1,ijatoms), asets(2,ijatoms)
              itask = atasks(1,ijsets)
              iset   = ival (3,itask)
              jset   = ival (4,itask)
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              ncob = npgfb(jset)*ncoset(lb_max(jset))
              sgfb = first_sgfb(1,jset)
              IF (iatom <= jatom) THEN
                CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                           1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           p_block(sgfa,sgfb),SIZE(p_block,1),&
                           0.0_dp,work(1,1),maxco)
                CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                           1.0_dp,work(1,1),maxco,&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,pab(1,1),maxco)
              ELSE
                CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                           1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           p_block(sgfb,sgfa),SIZE(p_block,1),&
                           0.0_dp,work(1,1),maxco)
                CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                           1.0_dp,work(1,1),maxco,&
                           sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           0.0_dp,pab(1,1),maxco)
              END IF
              DO itask = atasks(1,ijsets),atasks(2,ijsets)
                rab(:) = dab_local (:,itask)
                rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                rb(:) = ra(:) + rab(:)
                ipgf   = ival (5,itask)
                jpgf   = ival (6,itask)
                na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                na2 = ipgf*ncoset(la_max(iset))
                nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                nb2 = jpgf*ncoset(lb_max(jset))
 
                IF ((iatom == jatom).AND.&
                    (iset == jset).AND.&
                    (ipgf == jpgf)) THEN
                  scale = 1.0_dp
                ELSE
                  scale = 2.0_dp
                END IF
 
                IF ( nthread > 1 ) THEN
                  IF (iatom <= jatom) THEN
                    CALL collocate_pgf_product_rspace(&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread,compute_tau=my_compute_tau)
                  ELSE
                    CALL collocate_pgf_product_rspace(&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread,compute_tau=my_compute_tau)
                  END IF
                ELSE
                  IF (iatom <= jatom) THEN
                    CALL collocate_pgf_product_rspace(&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,compute_tau=my_compute_tau)
                  ELSE
                    CALL collocate_pgf_product_rspace(&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,compute_tau=my_compute_tau)
                  END IF
                END IF
 
              END DO
 
            END DO
 
          END DO
!$OMP end parallel
          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            n = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)
            DO i=1,nthread
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,k) &
!$OMP shared(i,lb,ub,lgrid,rs_rho,n,igrid_level)
              DO j=lb(3),ub(3)
                k = lgrid%ldim*(i-1) + n*(j-lb(3)) + 1
                CALL daxpy (n,1._dp,lgrid%r(k),1,&
                  rs_rho(igrid_level)%rs_grid%r(lb(1),lb(2),j),1)
              END DO
            END DO
          END IF

        END DO

      END DO

    END DO

!   *** Release work storage ***

    CALL deallocate_matrix ( deltap )

    IF ( nthread > 1 ) THEN
      DEALLOCATE (lgrid%r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"lgrid%r")
    END IF

    DEALLOCATE (pabt,workt,ntasks,tasks,tasks_local,ival,latom,&
        dist_ab,dab_local,asets,atasks,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pabt,workt,ntasks,"//&
        "tasks,tasks_local,ival,latom,dist_ab,dab_local,asets,atasks")

    IF (gridlevel_info%ngrid_levels==1 .AND. .NOT.rho%pw%pw_grid%spherical&
         .AND. rebuild) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho%pw,"FORWARD")
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)
!FM copy also in mgrid_rspace?
       CALL coeff_transform_space(rho,rho_gspace)
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,"FORWARD")
       ENDDO

       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

       !   *** Get the density on the auxillary basis cutoff grid ***

       CALL coeff_zero(rho_gspace)
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL coeff_transform_space(mgrid_rspace(igrid_level),&
               mgrid_gspace(igrid_level))
          CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
       END DO

       !   *** Save real-space density ***
       IF (.NOT. rebuild) THEN
          CALL coeff_sumup(rho_gspace_old,rho_gspace)
       END IF

       CALL coeff_transform_space(rho_gspace,rho)

    END IF
    
    CALL coeff_copy(rho_gspace,rho_gspace_old)
    total_rho = calculate_total_rho(rho)

    CALL timestop(0.0_dp,handle)

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
                                    qs_env,error)

    TYPE(cp_fm_type), POINTER :: mo_vectors
    INTEGER :: ivector

    TYPE(coeff_type), INTENT(INOUT)              :: rho,rho_gspace
    TYPE(coeff_type), DIMENSION(:), POINTER      :: mgrid_rspace,mgrid_gspace
    TYPE(cube_info_type), DIMENSION(:), POINTER  :: cube_info
    TYPE(gridlevel_info_type), INTENT(IN)        :: gridlevel_info
    TYPE(l_info_type), INTENT(IN)                :: l_info
    TYPE(qs_environment_type), INTENT(IN)        :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_wavefunction (MODULE qs_collocate_density)"

!   *** Local variables ***
    REAL(KIND = dp), DIMENSION(:), POINTER   :: eigenvector   

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(dft_control_type), POINTER   :: dft_control
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    REAL(KIND = dp) :: dab,eps_rho_rspace,kind_radius_b,rab2,scale,zetp
    INTEGER  :: ab,first_pgfb,first_setb,handle,i,iatom,ico,igrid_level,ikind,&
                ineighbor,ipgf,iset,ishell,istat,j,jatom,jco,jkind,jpgf,jset,&
                jshell,katom,maxco,maxsgf_set,na1,na2,nb1,nb2,ncoa,ncob,&
                nkind,nneighbor,nseta,nsetb,sgfa,sgfb

    REAL(KIND = dp), DIMENSION(3) :: ra,rab,rb

    REAL(KIND = dp), DIMENSION(:,:), POINTER :: pab,work

    REAL(KIND = dp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: p_block,rpgfa,rpgfb,&
                                         sphi_a,sphi_b,zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    TYPE(REALSPACE_GRID_P_TYPE), DIMENSION(:), POINTER  :: rs_rho
    INTEGER :: offset,natom,nao,ncol_global
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), POINTER :: rs_pools
    LOGICAL :: local

!   ---------------------------------------------------------------------------
    CALL timeset("calculate_wavefunction","I","",handle)

    CALL cp_fm_get_info(matrix=mo_vectors,nrow_global=nao,ncol_global=ncol_global) 

    NULLIFY(eigenvector, atomic_kind,cell,dft_control, orb_basis_set,&
         sab_orb_neighbor_list,&
         sab_orb_neighbor_node, atomic_kind_set, sab_orb, particle_set,&
         pab,work,set_radius_a, set_radius_b, la_max, la_min,&
         lb_max, npgfa, npgfb, nsgfa, nsgfb, p_block,&
         rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb, first_sgfa, first_sgfb,&
         rs_rho,pw_env,rs_pools)

    ALLOCATE(eigenvector(nao), stat=istat)
    IF (istat.NE.0) THEN
       CALL stop_program("calculate_wavefunction","eigenvector")
    ENDIF
    ! only some processors have non-zero elements. The sum of eigenvectors on 
    ! all cpus is the real eigenvector. combined with the rs->pw this will give the
    ! real density
    DO i=1,nao 
       CALL cp_fm_get_element(mo_vectors,i,ivector,eigenvector(i),local)
       IF (.NOT. local) eigenvector(i)=0.0_dp ! some more intelligent scheme is needed
    ENDDO

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    pw_env=pw_env)
    CALL pw_env_get(pw_env, rs_pools=rs_pools,error=error)

    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*dp_size)

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
      rab(:) = 0.0_dp
      rab2  = 0.0_dp
      dab   = 0.0_dp

      DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO i=1,nsgfa(iset)
            work(i,1)=eigenvector(offset+i)
         ENDDO

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dp,pab(1,1),SIZE(pab,1))

         DO ipgf=1,npgfa(iset)

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))

            scale = 1.0_dp
            zetp = zeta(ipgf,iset)
            igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

            CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_dp,0,&
                        ra,rab,rab2,scale,pab,na1-1,0,&
                        rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
            mgrid_rspace(igrid_level)%pw,"FORWARD")
    ENDDO
    
    CALL rs_pools_give_back_rs_vect(rs_pools,rs_rho,error=error)

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),&
           mgrid_gspace(igrid_level))
      CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    DEALLOCATE(eigenvector)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_wavefunction
! *****************************************************************************

  FUNCTION calculate_total_rho(rho) RESULT(total_rho)

    TYPE(coeff_type), INTENT(IN), TARGET :: rho

    REAL(KIND = dp) :: total_rho

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
        total_rho = 0._dp
      END IF
    ELSE
      CALL stop_program(routine,"No density coefficients available")
    END IF
    IF (rho%pw%pw_grid%para%mode.NE.0) THEN
       CALL mp_sum(total_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_rho

! *****************************************************************************

  FUNCTION calculate_total_abs_rho(rho) RESULT(total_abs_rho)

    TYPE(coeff_type), INTENT(IN), TARGET :: rho

    REAL(KIND = dp) :: total_abs_rho

!   *** Local variables ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "FUNCTION calculate_total_rho (MODULE qs_collocate_density)"

!   *** Local variables ***

    INTEGER, DIMENSION(:), POINTER :: lb_grid

!   ---------------------------------------------------------------------------
    total_abs_rho = 0.0_dp
    IF (ASSOCIATED(rho%pw%cr3d)) THEN
      total_abs_rho = rho%pw%pw_grid%dvol*SUM(ABS(rho%pw%cr3d))
    ELSE
      CALL stop_program(routine,"Need density coefficients in real space !")
    END IF
    IF (rho%pw%pw_grid%para%mode.NE.0) THEN
       CALL mp_sum(total_abs_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_abs_rho

! *****************************************************************************
  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,o1,o2,&
                                          rsgrid,cube_info,l_info,&
                                          eps_rho_rspace,lgrid,ithread,compute_tau)

    TYPE(realspace_grid_type), INTENT(INOUT), TARGET :: rsgrid
    TYPE(cube_info_type), INTENT(IN)  :: cube_info
    TYPE(l_info_type), INTENT(IN)     :: l_info
    REAL(KIND = dp), INTENT(IN)       :: eps_rho_rspace,rab2,scale,zeta,zetb
    INTEGER, INTENT(IN)               :: la_max,la_min,lb_max,lb_min
    REAL(KIND = dp), DIMENSION(3), INTENT(IN)   :: ra,rab
    TYPE(lgrid_type),OPTIONAL         :: lgrid
    INTEGER, INTENT(IN), OPTIONAL     :: ithread
    LOGICAL, INTENT(IN), OPTIONAL     :: compute_tau

    REAL(KIND = dp), DIMENSION(:,:), POINTER    :: pab
    INTEGER                              :: o1,o2

!   *** Local variables ***

    REAL(KIND = dp) :: f,prefactor,radius,zetp
    INTEGER  :: i,ico,ig,j,jco,jg,k,kg,ico_l,jco_l

    REAL(KIND = dp), DIMENSION(3) :: dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER :: lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER :: coef_max,lxa,lxyz,l,gridbounds(2,3)
    INTEGER, POINTER,DIMENSION(:) :: ly_max,lz_max
    REAL(KIND = dp) :: a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,&
                rpg,s,pg

    REAL(KIND = dp), POINTER :: dpx(:,:)
    REAL(KIND = dp), POINTER :: pzyx(:), alpha(:,:)
    INTEGER,  POINTER :: ipzyx(:,:,:,:,:,:)
    REAL(KIND = dp), POINTER :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    INTEGER, POINTER, DIMENSION(:)  :: sphere_bounds
    REAL(KIND = dp)  :: xa,xb,xap,xbp,polypart,cutoff
    INTEGER  :: start,length,offset,cmax
    INTEGER  :: handle
    INTEGER  :: ithread_l
    LOGICAL  :: my_compute_tau
    INTEGER  :: la_min_local,la_max_local, lb_min_local, lb_max_local,o1_local,o2_local
    REAL(KIND = dp), DIMENSION(:,:), POINTER    :: pab_local

!   ---------------------------------------------------------------------------

    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF

    IF (PRESENT(compute_tau)) THEN
       my_compute_tau=compute_tau
    ELSE
       my_compute_tau=.FALSE.
    ENDIF

    zetp      = zeta + zetb
    f         = zetb/zetp
    prefactor = scale*EXP(-zeta*f*rab2)
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)  
    rb(:)     = ra(:)+rab(:)
    cutoff    = 0.0_dp
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                    zetp,eps_rho_rspace,prefactor,cutoff)

    IF (radius .EQ. 0.0_dp ) THEN
      RETURN
    END IF

    ! it's a choice to compute lX_min/max, pab here, 
    ! this way we get the same radius as we use for the corresponding density
    IF (my_compute_tau) THEN
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
        ! (ddx pgf_a ) (ddx pgf_b) = (a pgf_{a-1} - 2*zeta*pgf_{a+1}*(b pgf_{x-1} - 2 * zetb * pgf_{b+1})
        ! cleaner would possibly be to touch pzyx directly (avoiding the following allocate)
        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                     ! this element of pab results in 12 elements of pab_local
                     ico=coset(lxa,lya,lza)
                     jco=coset(lxb,lyb,lzb)
                     ! x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)

                     ico_l=coset(MAX(lxa-1,0),lya,lza)
                     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxa* lxb*pab(o1+ico,o2+jco)
                     ico_l=coset(MAX(lxa-1,0),lya,lza)
                     jco_l=coset(   (lxb+1  ),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lxa*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(   (lxa+1  ),lya,lza)
                     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lxb*pab(o1+ico,o2+jco)
                     ico_l=coset(   (lxa+1  ),lya,lza)
                     jco_l=coset(   (lxb+1  ),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

                     ! y 

                     ico_l=coset(lxa,MAX(lya-1,0),lza)
                     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lya* lyb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,MAX(lya-1,0),lza)
                     jco_l=coset(lxb,   (lyb+1  ),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lya*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,   (lya+1  ),lza)
                     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lyb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,   (lya+1  ),lza)
                     jco_l=coset(lxb,   (lyb+1  ),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

                     ! z 

                     ico_l=coset(lxa,lya,MAX(lza-1,0))
                     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lza* lzb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,MAX(lza-1,0))
                     jco_l=coset(lxb,lyb,   (lzb+1  ))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lza*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,   (lza+1  ))
                     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lzb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,   (lza+1  ))
                     jco_l=coset(lxb,lyb,   (lzb+1  ))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

              ENDDO
              ENDDO
           ENDDO
           ENDDO
        ENDDO
        ENDDO
        o1_local=0
        o2_local=0
        pab_local=pab_local * 0.5_dp
    ELSE
        la_max_local=la_max
        la_min_local=la_min
        lb_max_local=lb_max
        lb_min_local=lb_min
        pab_local => pab
        o1_local=o1
        o2_local=o2
    ENDIF

    coef_max=la_max_local+lb_max_local+1
!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min_local,la_max_local,lb_min_local,lb_max_local,ithread_l,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is 1,1,1 on the grid
!

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
    DO i=1,3
      IF ( rsgrid % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i),ng(i))+rsgrid%lb(i)
        offset=offset-rsgrid%lb_local(i)+1
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

!   *** initialise the p terms and loop logic
    lxyz=0
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
          DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
          DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             pzyx(lxyz)=prefactor*pab_local(o1_local+ico,o2_local+jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
 
!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,dp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max_local
       zbp=1.0_dp
       DO lzb=0,lb_max_local
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max_local
      DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
          DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
          DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
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
      rpg = REAL(ig,dp)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max_local
       ybp=1.0_dp
       DO lyb=0,lb_max_local
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max_local
      DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO
 
!   *** make the alpha matrix ***
    alpha(:,:)=0.0_dp
    lx=0
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       lx=lx+1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO

    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,dp)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

    IF ( PRESENT ( lgrid ) ) THEN
      ig = lgrid%ldim * ithread_l + 1
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),lgrid%r(ig),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    ELSE
!   *** do the loop over the grid
!   notice this is not the same as critical or so, since we may have 
!   several different rsgrids that all use the same function
!   I guess we need a flush (of the grid) to guarantee that we are working 
!   on an uptodate copy of the grid
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    END IF

    IF (my_compute_tau) THEN
       DEALLOCATE(pab_local)
    ENDIF

  END SUBROUTINE collocate_pgf_product_rspace

! *****************************************************************************

END MODULE qs_collocate_density
