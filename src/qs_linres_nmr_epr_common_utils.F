!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief given the response wavefunctions obtained by the application
!>      of the (rxp), p, and ((dk-dl)xp) operators,
!>      here the current density vector (jx, jy, jz)
!>      is computed for the 3 directions of the magnetic field (Bx, By, Bz)
!> \par History
!>      created 02-2006 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_nmr_epr_common_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pools_create_rs_vect
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: gaussi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace,&
                                             density_rs2pw,&
                                             lgrid_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_epr_op,                ONLY: set_vecp
  USE qs_modify_pab_block,             ONLY: FUNC_ADBmDAB,&
                                             FUNC_ARDBmDARB
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             rs_find_node,&
                                             rs_grid_zero
  USE realspace_task_selection,        ONLY: distribute_matrix,&
                                             int2pair,&
                                             pair2int,&
                                             rs_get_my_tasks
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: calculate_jrho_resp,mult_G_ov_G2_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_epr_common_utils'

CONTAINS

! *****************************************************************************
!> \brief Calculation of the idir component of the response current density
!>       in the presence of a constant magnetic field in direction iB
!>       the current density is collocated on the pw grid in real space
!> \note
!>       The collocate is done in three parts, one for each density matrix
!>       In all cases the density matrices and therefore the collocation
!>       are not symmetric, that means that all the pairs (ab and ba) have
!>       to be considered separately
!> 
!>       mat_jp_{\mu\nu} is multiplied by
!>           f_{\mu\nu} = \phi_{\mu} (d\phi_{\nu}/dr)_{idir} -
!>                        (d\phi_{\mu}/dr)_{idir} \phi_{\nu}
!> 
!>       mat_jp_rii_{\mu\nu} is multiplied by
!>           f_{\mu\nu} = \phi_{\mu} (r - R_{\nu})_{iiiB} (d\phi_{\nu}/dr)_{idir} -
!>                        (d\phi_{\mu}/dr)_{idir} (r - R_{\nu})_{iiiB} \phi_{\nu} +
!>                         \phi_{\mu} \phi_{\nu}  (last term only if iiiB=idir)
!> 
!>       mat_jp_riii_{\mu\nu} is multiplied by
!>                             (be careful: change in sign with respect to previous)
!>           f_{\mu\nu} = -\phi_{\mu} (r - R_{\nu})_{iiB} (d\phi_{\nu}/dr)_{idir} +
!>                        (d\phi_{\mu}/dr)_{idir} (r - R_{\nu})_{iiB} \phi_{\nu} -
!>                         \phi_{\mu} \phi_{\nu}  (last term only if iiB=idir)
!> 
!>       All the terms sum up to the same grid
! *****************************************************************************
  SUBROUTINE calculate_jrho_resp(mat_jp,mat_jp_rii,mat_jp_riii,iB,idir,&
       rho_rs, rho_gs, qs_env, soft_valid, error)

    TYPE(real_matrix_type), POINTER          :: mat_jp, mat_jp_rii, &
                                                mat_jp_riii
    INTEGER, INTENT(IN)                      :: iB, idir
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_rs, rho_gs
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_resp', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, bcol, brow, curr_tasks, dest, dir, handle, i, iatom, &
      iatom_old, igrid_level, iiB, iiiB, ikind, ikind_old, ilist, inode, &
      ipgf, iset, iset_old, istat, itask, ithread, jatom, jatom_old, jkind, &
      jkind_old, jpgf, jset, jset_old, maxco, maxpgf, maxset, maxsgf, &
      maxsgf_set, na1, na2, natom, nb1, nb2, ncoa, ncob, nkind, nlist, nnode, &
      nseta, nsetb, ntasks, nthread, sgfa, sgfb, stat, tp(3)
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, distributed_rs_grids, failure, &
      map_consistent, my_minimum_image, my_soft
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
                                                kind_radius_a, kind_radius_b, &
                                                Lxo2, Lyo2, Lzo2, rab2, &
                                                scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: dist_ab, jp_block, &
      jp_block_rii, jp_block_riii, jpab, jpab_ii, jpab_iii, jpblock, &
      jpblock_rii, jpblock_riii, rpgfa, rpgfb, sphi_a, sphi_b, work, zeta, &
      zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: jpabt, jpabt_ii, jpabt_iii, &
                                                workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(real_matrix_type), POINTER          :: deltajp, deltajp_rii, &
                                                deltajp_riii
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(section_vals_type), POINTER         :: input, interp_section

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,dist_ab,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         dist_ab,tasks,workt)
    NULLIFY(deltajp,deltajp_rii,deltajp_riii)
    NULLIFY(jp_block,jp_block_rii,jp_block_riii)
    NULLIFY(jpblock,jpblock_rii,jpblock_riii)
    NULLIFY(jpabt,jpabt_ii,jpabt_iii)
    NULLIFY(lgrid%r)

    !    debug_count=debug_count+1

    CALL timeset(routineN,"I"," ",handle)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         sab_all=sab_orb,&
         para_env=para_env,&
         input=input,&
         pw_env=pw_env,error=error)

    ! Component of appearing in the vector product rxp, iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    !    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)

    !   Check that the neighbor list with all the pairs is associated
    CPPrecondition(ASSOCIATED(sab_orb),cp_failure_level,routineP,error,failure)
    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)

    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( .NOT. ALL (rs_rho(igrid_level)%rs_grid%perd == 1 ) ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1

    !   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf=maxsgf,&
         maxsgf_set=maxsgf_set)
    my_minimum_image = .TRUE.
    !    IF(PRESENT(minimum_image)) THEN
    !       my_minimum_image=minimum_image
    Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
    Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
    Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    !    END IF

    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    ntasks = 0
    IF ( nthread > 1 ) THEN
       DO igrid_level = 1,gridlevel_info%ngrid_levels
          ntasks = MAX(ntasks,rs_rho(igrid_level)%rs_grid%ngpts_local)
       END DO
       ntasks = ntasks*nthread
       CALL reallocate(lgrid%r,1,ntasks)
    END IF

    nkind = SIZE(atomic_kind_set)

    CALL reallocate(jpabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_ii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(jpabt_iii,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(tasks,1,5,1,max_tasks)
    CALL reallocate(dist_ab,1,3,1,max_tasks)

    tasks = 0
    ntasks = 0
    curr_tasks = SIZE(tasks,2)

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    !   *** Initialize working density matrix ***

    ! distributed rs grids require a matrix that will be changed (rs_get_my_tasks)
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       CALL allocate_matrix(matrix=deltajp,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
       CALL allocate_matrix(matrix=deltajp_rii,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
       CALL allocate_matrix(matrix=deltajp_riii,&
            nrow=mat_jp%nrow,&
            ncol=mat_jp%ncol,&
            nblock_row=mat_jp%nblock_row,&
            nblock_col=mat_jp%nblock_col,&
            first_row=mat_jp%first_row(:),&
            last_row=mat_jp%last_row(:),&
            first_col=mat_jp%first_col(:),&
            last_col=mat_jp%last_col(:),&
            matrix_name="DeltaP",&
            sparsity_id=-1, &   ! basically unknown sparsity in parallel
            matrix_symmetry=mat_jp%symmetry,error=error)
    ELSE
       deltajp=>mat_jp
       deltajp_rii=>mat_jp_rii
       deltajp_riii=>mat_jp_riii
    ENDIF

    loop_ikind: DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_soft, &
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            first_sgf=first_sgfa,&
            kind_radius=kind_radius_a,&
            lmax=la_max,&
            lmin=la_min,&
            npgf=npgfa,&
            nset=nseta,&
            nsgf_set=nsgfa,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi=sphi_a,&
            zet=zeta)

       loop_jkind: DO jkind=1,nkind

          atomic_kind => atomic_kind_set(jkind)

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               softb = my_soft, &
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

          loop_ilist:DO ilist = 1, nlist

             CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                  atom=iatom,nnode=nnode)

             ra(:) = pbc(particle_set(iatom)%r,cell)

             sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

             loop_inode: DO inode = 1, nnode

                CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                     neighbor=jatom,&
                     r=rab(:))

                IF(my_minimum_image) THEN
                   IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) THEN
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      CYCLE
                   END IF
                END IF

                brow = iatom
                bcol = jatom

                CALL get_block_node(matrix=mat_jp,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block)
                CALL get_block_node(matrix=mat_jp_rii,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block_rii)
                CALL get_block_node(matrix=mat_jp_riii,&
                     block_row=brow,&
                     block_col=bcol,&
                     BLOCK=jp_block_riii)

                IF (.NOT.ASSOCIATED(jp_block)) THEN
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                   CYCLE
                END IF

                IF (distributed_rs_grids) THEN
                   NULLIFY (jpblock,jpblock_rii,jp_block_riii )
                   CALL add_block_node ( deltajp, brow, bcol, jpblock ,error=error)
                   jpblock = jp_block
                   CALL add_block_node ( deltajp_rii, brow, bcol, jpblock_rii ,error=error)
                   jpblock_rii = jp_block_rii
                   CALL add_block_node ( deltajp_riii, brow, bcol, jpblock_riii ,error=error)
                   jpblock_riii = jp_block_riii
                ELSE
                   jpblock => jp_block
                   jpblock_rii => jp_block_rii
                   jpblock_riii => jp_block_riii
                ENDIF

                IF (.NOT. map_consistent) THEN
                   IF ( ALL ( 100.0_dp*ABS(jpblock)      < eps_rho_rspace) .AND. & 
                        ALL ( 100.0_dp*ABS(jpblock_rii)  < eps_rho_rspace) .AND. & 
                        ALL ( 100.0_dp*ABS(jpblock_riii) < eps_rho_rspace) ) THEN 
                      sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                      CYCLE
                   END IF
                END IF

                rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                dab = SQRT(rab2)

                loop_iset : DO iset=1,nseta
                   IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE

                   loop_jset: DO jset = 1,nsetb
                      IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                      loop_ipgf: DO ipgf=1,npgfa(iset)
                         IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                         loop_jpgf: DO jpgf=1,npgfb(jset)
                            IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                            zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                            IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                               igrid_level = 1
                            ELSE
                               igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                            ENDIF

                            ntasks = ntasks + 1
                            IF ( ntasks > curr_tasks ) THEN
                               curr_tasks = curr_tasks*mult_tasks
                               CALL reallocate(tasks,1,5,1,curr_tasks)
                               CALL reallocate(dist_ab,1,3,1,curr_tasks)
                            END IF

                            rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                            rp(:) = pbc(rp,cell)
                            DO dir = 1,3
                               tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),rp)*rs_rho(igrid_level)%rs_grid%npts(dir))
                               tp(dir) = MODULO ( tp(dir), rs_rho(igrid_level)%rs_grid%npts(dir) )
                               tp(dir) = tp(dir) + rs_rho(igrid_level)%rs_grid%lb(dir)
                            END DO

                            IF(rs_rho(igrid_level)%rs_grid%distributed) THEN
                               CALL rs_find_node(rs_rho(igrid_level)%rs_grid,tp,dest)
                               tasks(4,ntasks) = 1
                            ELSE
                               dest = rs_rho(igrid_level)%rs_grid % my_pos
                               tasks(4,ntasks) = 0
                            ENDIF

                            tasks (1,ntasks) = dest
                            tasks (2,ntasks) = rs_rho(igrid_level)%rs_grid % my_pos
                            CALL pair2int(tasks(3,ntasks),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

                            tasks(5,ntasks) = 1 

                            dist_ab (:,ntasks) = rab(:)

                         END DO loop_jpgf

                      END DO loop_ipgf

                   END DO loop_jset

                END DO loop_iset

                sab_orb_neighbor_node => next(sab_orb_neighbor_node)

             END DO loop_inode

             sab_orb_neighbor_list => next(sab_orb_neighbor_list)

          END DO loop_ilist

       END DO loop_jkind

    END DO loop_ikind

    ! sorts / redistributes the task list
    CALL rs_get_my_tasks ( rs_rho(1)%rs_grid, distributed_rs_grids, ntasks, natom,&
            maxset, maxpgf, tasks, dist_ab, atom_pair_send, atom_pair_recv,&
            .FALSE. , error)

    IF (distributed_rs_grids) THEN
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp_rii, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
        CALL distribute_matrix (rs_rho(1)%rs_grid, deltajp_riii, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
    ENDIF

    ithread = 0
    jpab => jpabt(:,:,ithread)
    jpab_ii => jpabt_ii(:,:,ithread)
    jpab_iii => jpabt_iii(:,:,ithread)
    work => workt(:,:,ithread)

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1 
    ikind_old = -1 ; jkind_old = -1

    loop_tasks: DO itask = 1,ntasks

       CALL int2pair(tasks(3,itask),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN 

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          IF(iatom .NE. iatom_old ) ra(:) = pbc(particle_set(iatom)%r,cell)

          brow = iatom
          bcol = jatom

          IF (ikind .NE. ikind_old ) THEN
             CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

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
          ENDIF

          IF (jkind .NE. jkind_old ) THEN

             CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

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

          ENDIF

          CALL get_block_node(matrix=deltajp,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block)
          CALL get_block_node(matrix=deltajp_rii,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block_rii)
          CALL get_block_node(matrix=deltajp_riii,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=jp_block_riii)
          IF (.NOT.ASSOCIATED(jp_block)) &
               CALL stop_program(routineP,"p_block not associated in deltap")

          iatom_old = iatom
          jatom_old = jatom
          ikind_old = ikind
          jkind_old = jkind

          atom_pair_changed = .TRUE.

       ELSE

          atom_pair_changed = .FALSE.

       ENDIF

       IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          ! Decontraction step for the selected blocks of the 3 density matrices
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block(sgfa,sgfb),SIZE(jp_block,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab(1,1),maxco)
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_rii(sgfa,sgfb),SIZE(jp_block_rii,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_ii(1,1),maxco)
          CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               jp_block_riii(sgfa,sgfb),SIZE(jp_block_riii,1),&
               0.0_dp,work(1,1),maxco)
          CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
               1.0_dp,work(1,1),maxco,&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,jpab_iii(1,1),maxco)

          iset_old = iset
          jset_old = jset

       ENDIF

       rab(:) = dist_ab (:,itask)
       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       rb(:) = ra(:) + rab(:)
       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       scale = 1.0_dp

       ! Four calls to the general collocate density, to multply the correct function
       ! to each density matrix

       ! here the decontracted mat_jp_{ab} is multiplied by
       !     f_{ab} = g_{a} (dg_{b}/dr)_{idir} - (dg_{a}/dr)_{idir} g_{b}
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ADBmDAB,&
            idir=idir,&
            map_consistent=map_consistent,error=error)
       ! here the decontracted mat_jp_rii{ab} is multiplied by
       !     f_{ab} = g_{a} (r - R_{b})_{iiB} (dg_{b}/dr)_{idir} -
       !             (dg_{a}/dr)_{idir} (r - R_{b})_{iiB} g_{b}
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab_ii,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ARDBmDARB,&
            idir=idir,ir=iiiB,&
            map_consistent=map_consistent,error=error)
       ! here the decontracted mat_jp_riii{ab} is multiplied by
       !     f_{ab} = -g_{a} (r - R_{b})_{iiB} (dg_{b}/dr)_{idir} +
       !             (dg_{a}/dr)_{idir} (r - R_{b})_{iiB} g_{b}
       scale = -1.0_dp
       CALL collocate_pgf_product_rspace(la_max(iset),zeta(ipgf,iset),&
            la_min(iset),lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
            ra,rab,rab2,scale,jpab_iii,na1-1,nb1-1,&
            rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
            eps_rho_rspace,&
            ga_gb_function=FUNC_ARDBmDARB,&
            idir=idir,ir=iiB,&
            map_consistent=map_consistent,error=error)

    END DO loop_tasks

    !   *** Release work storage ***

    IF (distributed_rs_grids) THEN
       CALL deallocate_matrix ( deltajp ,error=error)
       CALL deallocate_matrix ( deltajp_rii ,error=error)
       CALL deallocate_matrix ( deltajp_riii ,error=error)
    END IF

    IF ( nthread > 1 ) THEN
       DEALLOCATE (lgrid%r,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"lgrid%r")
    END IF

    DEALLOCATE (jpabt,jpabt_ii,jpabt_iii,workt,tasks,dist_ab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"jpabt,workt,tasks,dist_ab")

    IF (distributed_rs_grids) THEN
       DEALLOCATE(atom_pair_send,STAT=stat)
       IF (stat/=0) CALL stop_memory("collocate_density","atom_pair_send")        
       DEALLOCATE(atom_pair_recv,STAT=stat)
       IF (stat/=0) CALL stop_memory("collocate_density","atom_pair_recv")
    ENDIF

    CALL density_rs2pw(pw_env,rs_rho,rho_rs,rho_gs,&
         interp_section=interp_section,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_jrho_resp

! *****************************************************************************
!> \brief Given the current density on the PW grid in reciprcal space
!>       (obtained by FFT), calculate the integral
!>         \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!>       which in reciprcal space reads  (for G/=0)
!>          i G/|G|^2 x J(G)
!> \note
!>      The G=0 component is not comnputed here, but can be evaluated
!>      through the susceptibility and added to the shift in a second time
!> 
!>      This method would not work for a non periodic system
!>      It should be generalized like the calculation of Hartree
!> \author MI
! *****************************************************************************
  SUBROUTINE  mult_G_ov_G2_grid(cell,pw_pool,rho_gspace,funcG_times_rho,idir,my_chi,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_p_type), POINTER                 :: rho_gspace
    TYPE(pw_p_type)                          :: funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
    REAL(dp), INTENT(IN)                     :: my_chi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ig, ng
    LOGICAL                                  :: failure
    REAL(dp)                                 :: g2
    TYPE(pw_grid_type), POINTER              :: grid
    CHARACTER(len=*), PARAMETER :: routineN = 'mult_G_ov_G2_grid', &
      routineP = moduleN//':'//routineN

    TYPE(pw_type), POINTER                   :: frho, influence_fn

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)

    CALL pw_pool_create_pw ( pw_pool, influence_fn,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE ,error=error)

    grid => influence_fn % pw_grid
    DO ig = grid % first_gne0, grid % ngpts_cut_local
       g2 = grid % gsq ( ig )
       influence_fn%cc(ig) = gaussi * grid % g(idir,ig)/g2
    END DO  ! ig
    IF ( grid % have_g0 )  influence_fn%cc ( 1 ) = 0.0_dp

    frho => funcG_times_rho%pw
    CALL pw_transfer (rho_gspace%pw,frho)

    ng = SIZE(grid % gsq)
    frho%cc(1:ng) = frho%cc(1:ng)*influence_fn % cc  ( 1 : ng )
    IF ( grid % have_g0 ) frho%cc(1) = my_chi

    CALL pw_pool_give_back_pw(pw_pool,influence_fn,&
         accept_non_compatible=.TRUE.,error=error)

  END SUBROUTINE  mult_G_ov_G2_grid

END MODULE qs_linres_nmr_epr_common_utils

