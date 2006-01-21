!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003, 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!


!!****** cp2k/kg_gpw_collocate_den [1.0] *
!!
!!   NAME
!!     kg_gpw_collocate_den
!!
!!   FUNCTION
!!     Calculate the plane wave density for each molecule independently
!!     by collocating the primitive Gaussian functions on the grids 
!!     constructed in the molecolar box. The molecule is centered in the
!!     small box e no PBC are used.
!!     Each molecule is handled by one processor entirely, therefore
!!     no further distribution is allowed, unless OMP is used (if it works)
!!
!!   AUTHOR
!!     MI (05.01.2005)
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_collocate_den

! *****************************************************************************

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE l_utils,                         ONLY: l_info_type,&
                                             return_l_info
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_restrict_l,&
                                             pw_prolongate_l,&
                                             pw_sumup
  USE qs_collocate_density,            ONLY: calculate_total_rho,&
                                             collocate_pgf_product_rspace,&
                                             lgrid_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_grid_zero,&
                                             rs_pw_transfer,&
                                             rs2pw,&
                                             pw2rs
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: module_name='kg_gpw_collocate_den'

  INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                              max_tasks = 2000
  REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

! *** Public functions ***

  PUBLIC :: calculate_rho_mol, integrate_mol_potential

! *****************************************************************************

CONTAINS

!****f* kg_gpw_collocate_den/calculate_rho_mol
!!
!!   NAME
!!     calculate_rho_mol
!!
!!   FUNCTION
!!     This function is called within a loop over all the molecules
!!     to calculate the molecular density on the molecular grid only. 
!!     The molecular densities are used by KG_GPW to calculate the Kinetic energy
!!     correction by subtracting the sum of the kinetic energies calculated
!!     for each molecule independently
!!
!!   NOTES
!!     In KG_GPW a molecular pw environment is created for each molecule kind
!!     which is based on a molecular cell box  that does not use PBC
!!     The molecule is located in the center of the bos and the coordinates passed
!!     to this routine are internal coordinates.
!!     At this level, it should not be necessary to have the atomic positions in the global box 
!!
!!   INPUTS
!!     - qs_env: the qs environment
!!     - matrix_p: global density matrix
!!     - rho_r: molecular density on the molecular real space mesh
!!     - rho_g: molecular density on the molecular reciprocal space mesh
!!     - total_rho: integral of the molecular density (output)
!!     - atom: global index of the atoms in the molecule (input)
!!     - kind: index of the kind of the atoms in the molecule (input)
!!     - ratom: internal coordinates of the atoms in the box of the molecule (input)
!!     - compute_tau: logical to control if tau is required
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!    
!!
!!*** **************************************************************************

  SUBROUTINE calculate_rho_mol(qs_env,matrix_p,rho_r,rho_g,total_rho,pw_env,&
                               atom,kind,ratom,compute_tau,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_type), POINTER          :: matrix_p
    TYPE(coeff_type), INTENT(INOUT)          :: rho_r, rho_g
    REAL(dp), INTENT(OUT)                    :: total_rho
    TYPE(pw_env_type), POINTER               :: pw_env
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom,kind
    REAL(dp), DIMENSION(:,:), INTENT(IN)     :: ratom
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_rho_mol',&
      routineP = module_name//':'//routineN

    INTEGER :: brow, bcol, curr_tasks, first_pgfb, first_setb, handle, i, iat, &
               iatom, igrid_level, ikind, ikind_old, ipgf, itask, ithread, iset, istat,&
               j, jat, jatom, jkind, jkind_old, jpgf, jset, k, maxco, maxsgf, &
               maxsgf_set, n, na1, na2, nat_mol, nb1, nb2, ncoa, ncob, nthread,&
               nseta, nsetb, num_tasks, sgfa, sgfb, igridlevel
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, ntasks, ub
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, &
                                                first_sgfb, ival, latom, &
                                                tasks_local
    INTEGER, DIMENSION(:, :, :), POINTER     :: tasks
    LOGICAL                                  :: failure, map_consistent, &
                                                my_compute_tau, distributed_rs_grids
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
                                                kind_radius_b, rab2, scale, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, pab, &
                                                rpgfa, rpgfb, sphi_a, &
                                                sphi_b, work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pabt, workt

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace, &
                                                mgrid_temp_rspace
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b
    TYPE(l_info_type)                        :: l_info
    TYPE(lgrid_type)                         :: lgrid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho

    INTEGER                                  :: omp_get_max_threads,&
                                                omp_get_thread_num
!   ---------------------------------------------------------------------------

    failure=.FALSE.

    ! by default, do not compute the kinetic energy density (tau)
    ! if compute_tau, all grids referening to rho are actually tau
    IF (PRESENT(compute_tau)) THEN 
       my_compute_tau = compute_tau
    ELSE
       my_compute_tau = .FALSE.
    ENDIF

    IF (my_compute_tau) THEN
       CALL timeset("calculate_rho_tau","I","",handle)
    ELSE
       CALL timeset("calculate_rho_elec","I","",handle)
    ENDIF

    NULLIFY(atomic_kind_set,dft_control)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control)

    ! *** assign from pw_env
    NULLIFY(gridlevel_info,cube_info)
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info
    l_info=pw_env%l_info

    ! *** set up the pw multi-grids 
    NULLIFY(rs_pools,pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    NULLIFY(mgrid_rspace, mgrid_temp_rspace)
    ALLOCATE(mgrid_rspace(SIZE(pw_pools)) ,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)

    IF (dft_control % qs_control % realspace_mgrids) THEN
        ALLOCATE(mgrid_temp_rspace(SIZE(pw_pools)) ,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pools_init_coeffs(pw_pools,mgrid_temp_rspace,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)
    ELSE
        ALLOCATE(mgrid_gspace(SIZE(pw_pools)) ,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE, error=error)
    ENDIF

    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    NULLIFY(rs_rho)
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, force_env_section=qs_env%input, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf,&
                             maxsgf_set=maxsgf_set)

    IF ( nthread > 1 ) THEN
      NULLIFY(lgrid%r)
      n=0
      DO igrid_level = 1,gridlevel_info%ngrid_levels
        n = MAX(n,rs_rho(igrid_level)%rs_grid%ngpts_local)
      END DO
      n = n*nthread
      CALL reallocate(lgrid%r,1,n)
    END IF

    NULLIFY(pabt,workt,ntasks,tasks,tasks_local,ival,latom)
    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(latom,1,2,1,max_tasks)
    curr_tasks = max_tasks

    nat_mol = SIZE(atom,1)

    ikind_old = 0
    jkind_old = 0
    DO iat = 1,nat_mol
       ikind = kind(iat)
       iatom = atom(iat)
       ra(1:3) = ratom(1:3,iat)

       IF(ikind /= ikind_old) THEN 
         NULLIFY(atomic_kind,orb_basis_a,la_max,la_min,npgfa,nsgfa,&
                 rpgfa,set_radius_a,sphi_a,zeta)
         atomic_kind => atomic_kind_set(ikind)

         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              orb_basis_set=orb_basis_a)

         IF (.NOT.ASSOCIATED(orb_basis_a)) CYCLE
         CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
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
         ikind_old = ikind 
       END IF
       IF (.NOT.ASSOCIATED(orb_basis_a)) CYCLE

       DO jat = iat, nat_mol
          jkind = kind(jat)
          jatom = atom(jat)
          rb(1:3) = ratom(1:3,jat)

          IF(jkind /= jkind_old) THEN
            NULLIFY(atomic_kind,orb_basis_b,lb_max,lb_min,npgfb,nsgfb,&
                 rpgfb,set_radius_b,sphi_b,zetb)

            atomic_kind => atomic_kind_set(jkind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              orb_basis_set=orb_basis_b)

            IF (.NOT.ASSOCIATED(orb_basis_b)) CYCLE

            CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
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
            jkind_old = jkind
          END IF
          IF (.NOT.ASSOCIATED(orb_basis_b)) CYCLE

          ntasks = 0
          tasks = 0

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          ! bad, should do better loop ordering XXXXXXXXXX
          NULLIFY(p_block)
          CALL get_block_node(matrix=matrix_p,&
                              block_row=brow,&
                              block_col=bcol,&
                              BLOCK=p_block)
          IF (.NOT.ASSOCIATED(p_block)) CYCLE

          IF (.NOT. map_consistent) THEN
             IF ( ALL ( 100.0_dp*ABS(p_block) < eps_rho_rspace) ) CYCLE
          END IF

          rab2 = 0.0_dp
          DO i = 1,3
            rab(i) = rb(i)-ra(i)
            rab2 = rab2 + rab(i)*rab(i)
          END DO
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

                  IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                     igrid_level = 1
                  ELSE
                     igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                  ENDIF

                  ntasks(igrid_level) = ntasks(igrid_level) + 1
                  n = ntasks(igrid_level)
                  IF ( n > curr_tasks ) THEN
                    curr_tasks = curr_tasks*mult_tasks
                    CALL reallocate(tasks,1,8,1,curr_tasks,&
                                    1,gridlevel_info%ngrid_levels)
                  END IF

                  tasks (1,n,igrid_level) = n
                  tasks (3,n,igrid_level) = iatom
                  tasks (4,n,igrid_level) = jatom
                  tasks (5,n,igrid_level) = iset
                  tasks (6,n,igrid_level) = jset
                  tasks (7,n,igrid_level) = ipgf
                  tasks (8,n,igrid_level) = jpgf

                END DO  ! jpgf
              END DO  ! ipgf
            END DO  ! jset
          END DO  ! iset

          DO igrid_level = 1, gridlevel_info%ngrid_levels
            num_tasks = ntasks ( igrid_level )
            IF ( num_tasks  > SIZE ( tasks_local, 2 ) ) &
              CALL reallocate(tasks_local,1,2,1,num_tasks )
            IF ( num_tasks  > SIZE ( ival, 2 ) ) &
              CALL reallocate(ival,1,6,1,num_tasks )
            IF ( num_tasks  > SIZE ( latom, 2 ) ) &
              CALL reallocate(latom,1,2,1,num_tasks )

!$OMP parallel do private(i)
            DO i=1,num_tasks
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
            END DO

!$OMP parallel do private(i)
            DO i=num_tasks+1,SIZE(tasks_local,2)
              tasks_local(1,i) = 0
              tasks_local(2,i) = 0
            END DO

            IF ( nthread > 1 ) THEN
              lb => rs_rho(igrid_level)%rs_grid%lb_local
              ub => rs_rho(igrid_level)%rs_grid%ub_local
              lgrid%ldim = rs_rho(igrid_level)%rs_grid%ngpts_local
!$OMP parallel private(ithread,n)
!$            ithread = omp_get_thread_num()
              n = ithread*lgrid%ldim + 1
              CALL dcopy(lgrid%ldim,0._dp,0,lgrid%r(n),1)
!$OMP end parallel
            END IF

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,itask,iset,jset,ncoa,ncob,sgfa,sgfb) &
!$OMP private(work,pab,istat,ipgf,jpgf,na1,na2,nb1,nb2,scale) &
!$OMP shared(iatom,jatom,ra,rb,rab,rab2,brow,bcol,p_block) &
!$OMP shared(maxco,maxsgf_set,ival,num_tasks) &
!$OMP shared(npgfa,npgfb,ncoset,la_max,lb_max,first_sgfa,first_sgfb) &
!$OMP shared(nsgfa,nsgfb,sphi_a,sphi_b,la_min,lb_min,zeta,zetb) &
!$OMP shared(rs_rho,igrid_level,cube_info,l_info,eps_rho_rspace,lgrid,nthread) &
!$OMP shared(workt,pabt,my_compute_tau,map_consistent)
            ithread = 0
!$          ithread = omp_get_thread_num()
            pab => pabt(:,:,ithread)
            work => workt(:,:,ithread)
!$OMP do 
            DO itask = 1,num_tasks 
 
              IF (.NOT.ASSOCIATED(p_block)) &
                  CALL stop_program(routineP,"p_block not associated in matrixp") 
              iset = ival (3,itask)
              jset = ival (4,itask)
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
                       l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread, &
                       compute_tau=my_compute_tau,map_consistent=map_consistent)
                ELSE
                  CALL collocate_pgf_product_rspace(&
                       lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                       la_max(iset),zeta(ipgf,iset),la_min(iset),&
                       rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                       rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                       l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread, &
                       compute_tau=my_compute_tau,map_consistent=map_consistent)
                END IF
              ELSE
                IF (iatom <= jatom) THEN
                  CALL collocate_pgf_product_rspace(&
                       la_max(iset),zeta(ipgf,iset),la_min(iset),&
                       lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                       ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                       rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                       l_info,eps_rho_rspace,compute_tau=my_compute_tau, &
                       map_consistent=map_consistent)
                ELSE
                  CALL collocate_pgf_product_rspace(&
                       lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                       la_max(iset),zeta(ipgf,iset),la_min(iset),&
                       rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                       rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                       l_info,eps_rho_rspace,compute_tau=my_compute_tau, &
                       map_consistent=map_consistent)
                END IF
              END IF
            END DO  ! itask

!$OMP end parallel
            IF ( nthread > 1 ) THEN
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
          END DO  ! igrid_level

       END DO  ! jat
    END DO  ! iat

    IF ( nthread > 1 ) THEN
      DEALLOCATE (lgrid%r,STAT=istat)
      CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error) 
    END IF

    DEALLOCATE (pabt,workt,ntasks,tasks,tasks_local,ival,latom,STAT=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho_r%pw,rs2pw)
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)
       CALL coeff_transform_space(rho_r,rho_g)
       IF (rho_r%pw%pw_grid%spherical) THEN ! rho_g = rho_r
          CALL coeff_transform_space(rho_g,rho_r)
       ENDIF
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw)
       ENDDO
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

       ! we want both rho_r and rho_g 
       IF (dft_control % qs_control % realspace_mgrids) THEN 
           CALL pw_copy(mgrid_rspace(gridlevel_info%ngrid_levels)%pw,mgrid_temp_rspace(gridlevel_info%ngrid_levels)%pw)
           DO igridlevel=gridlevel_info%ngrid_levels,2,-1
              ! prologate to the next grid and addup
              CALL pw_prolongate_l(mgrid_temp_rspace(igridlevel)%pw,mgrid_temp_rspace(igridlevel-1)%pw)
              ! add the next grid to the prolongated grid
              CALL pw_sumup(mgrid_rspace(igridlevel-1)%pw,mgrid_temp_rspace(igridlevel-1)%pw)
           ENDDO 
           CALL pw_copy(mgrid_temp_rspace(1)%pw,rho_r%pw)
           CALL coeff_transform_space(rho_r,rho_g)
        ELSE
           CALL coeff_zero(rho_g)
           DO igrid_level=1,gridlevel_info%ngrid_levels
              CALL coeff_transform_space(mgrid_rspace(igrid_level),&
                   mgrid_gspace(igrid_level))
              CALL coeff_sumup(mgrid_gspace(igrid_level),rho_g)
           END DO
           CALL coeff_transform_space(rho_g,rho_r)
       ENDIF
    END IF
    
    total_rho = calculate_total_rho(rho_r)

    ! *** give back the pw multi-grids
    IF (dft_control % qs_control % realspace_mgrids) THEN
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_temp_rspace,&
                                      error=error)
       DEALLOCATE(mgrid_temp_rspace,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
                                      error=error)
       DEALLOCATE(mgrid_gspace,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
         error=error)
    DEALLOCATE(mgrid_rspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)


  END SUBROUTINE calculate_rho_mol

!****f* kg_gpw_collocate_den/integrate_mol_potential
!!
!!   NAME
!!     integrate_mol_potential
!!
!!   FUNCTION
!!     This function is called within a loop over all the molecules
!!     to calculate  integrate a given potential defined on the molecular grid only. 
!!     This is used by KG_GPW to calculate the Kinetic energy
!!     correction as sum over the molecular contributions to the kinetic energy calculated
!!     through the kinetic energy functional of the single molecular density,
!!     for each molecule independently.
!!
!!   NOTES
!!     In KG_GPW a molecular pw environment is created for each molecule kind
!!     The corresponding molecular cell box  does not use PBC and is smaller than the global box.
!!     The molecule is located in the center of the box and the coordinates passed
!!     to this routine are internal coordinates.
!!     At this level, it should not be necessary to have the atomic positions in the global box 
!!
!!   INPUTS
!!     - qs_env: the qs environment
!!     - vxc_mol: potential on the radial grid as calculated from the molecular density
!!     - matrix_p: global density matrix
!!     - matrix_h: global ks matrix
!!     - pw_env: pw environment for this molecule kind
!!     - atom: global index of the atoms in the molecule (input)
!!     - kind: index of the kind of the atoms in the molecule (input)
!!     - ratom: internal coordinates of the atoms in the box of the molecule (input)
!!     - forces_mol: forces on the atoms of the molecule, coming from this additional potential
!!     - compute_tau: logical to control if tau is required
!!     - error: variable to control error logging, stopping, see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!    
!!
!!*** **************************************************************************

  SUBROUTINE integrate_mol_potential(qs_env,vxc_mol,matrix_p,&
                                    matrix_h,pw_env,atom, kind, ratom,&
                                    forces_mol,compute_tau,error)
 
    TYPE(qs_environment_type), POINTER       :: qs_env
     TYPE(coeff_type), INTENT(IN)            :: vxc_mol
    TYPE(real_matrix_p_type), INTENT(IN)     :: matrix_p
    TYPE(real_matrix_p_type), INTENT(INOUT)  :: matrix_h
    TYPE(pw_env_type), POINTER               :: pw_env
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom,kind
    REAL(dp), DIMENSION(:,:), INTENT(IN)     :: ratom
    REAL(dp), DIMENSION(:,:), INTENT(INOUT),&
      OPTIONAL                               :: forces_mol
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_mol_potential',&
      routineP = module_name//':'//routineN


    INTEGER :: auxbas_grid, bcol, brow, curr_tasks,  &
               handle, i, iat, iatom, igrid_level,&
               ijsets, ikind, ikind_old, ipgf, iset, istat, itask, ithread, jat, jatom, &
               jkind, jkind_old, jpgf, jset, maxco, maxsgf, maxsgf_set,&
               n, na1, na2, nat_mol, nb1, nb2, ncoa, ncob, &
               npme, nset_pairs, nseta, nsetb, nthread,&
               num_tasks, omp_get_max_threads, omp_get_thread_num, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, ntasks
    INTEGER, DIMENSION(:, :), POINTER        :: atasks, first_sgfa, &
                                                first_sgfb, ival, &
                                                tasks_local
    INTEGER, DIMENSION(:, :, :), POINTER     :: tasks
    LOGICAL                                  :: calculate_forces, failure, &
                                                map_consistent, my_compute_tau
    REAL(dp)                                 :: dab, eps_gvg_rspace,&
                                                kind_radius_b, rab2, scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, hab , p_block, pab, &
                                                rpgfa, rpgfb, sphi_a, &
                                                sphi_b, work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: habt, pabt, workt

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace, &
                                                mgrid_temp_rspace
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b
    TYPE(l_info_type)                        :: l_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v

!   ---------------------------------------------------------------------------

    failure=.FALSE.

    ! by default, do not compute the kinetic energy density (tau)
    ! if compute_tau, all grids referening to rho are actually tau
    IF (PRESENT(compute_tau)) THEN 
       my_compute_tau = compute_tau
    ELSE
       my_compute_tau = .FALSE.
    ENDIF

    IF (my_compute_tau) THEN
       CALL timeset("calculate_rho_tau","I","",handle)
    ELSE
       CALL timeset("calculate_rho_elec","I","",handle)
    ENDIF

    calculate_forces = .FALSE.
    IF(PRESENT(forces_mol)) calculate_forces = .TRUE.

    NULLIFY(atomic_kind_set,dft_control)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control)

    ! *** assign from pw_env for this molecule kind
    NULLIFY(gridlevel_info,cube_info)
    auxbas_grid=pw_env%auxbas_grid
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info
    l_info=pw_env%l_info

    ! *** set up the pw multi-grids
    NULLIFY(rs_pools,pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    NULLIFY(mgrid_rspace, mgrid_temp_rspace)
    ALLOCATE(mgrid_rspace(SIZE(pw_pools)) ,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)

    ! use either realspace or fft techniques to get the potential on the rs multigrids
    IF ( dft_control % qs_control % realspace_mgrids ) THEN
       ALLOCATE(mgrid_temp_rspace(SIZE(pw_pools)),stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_temp_rspace,&
                                 use_data = REALDATA3D,&
                                 in_space = REALSPACE, error=error)
       CALL coeff_copy(vxc_mol,mgrid_temp_rspace(1))
       CALL coeff_copy(vxc_mol,mgrid_rspace(1))
       DO igrid_level=2,gridlevel_info%ngrid_levels
          CALL pw_restrict_l(mgrid_temp_rspace(igrid_level-1)%pw,mgrid_rspace(igrid_level)%pw)
          CALL coeff_copy(mgrid_rspace(igrid_level),mgrid_temp_rspace(igrid_level))
       ENDDO
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_temp_rspace,&
                                      error=error)
       DEALLOCATE(mgrid_temp_rspace,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ELSE
       ALLOCATE(mgrid_gspace(SIZE(pw_pools)),stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
                                 use_data = COMPLEXDATA1D,&
                                 in_space = RECIPROCALSPACE, error=error)
       CALL coeff_transform_space(vxc_mol,mgrid_gspace(auxbas_grid))
       DO igrid_level=1,gridlevel_info%ngrid_levels
         IF ( igrid_level /= auxbas_grid ) THEN
              CALL coeff_copy(mgrid_gspace(auxbas_grid),mgrid_gspace(igrid_level))
              CALL coeff_transform_space(mgrid_gspace(igrid_level),&
                                             mgrid_rspace(igrid_level))
         ELSE
              IF (mgrid_gspace(auxbas_grid)%pw%pw_grid%spherical) THEN
                  CALL coeff_transform_space(mgrid_gspace(auxbas_grid),&
                                             mgrid_rspace(auxbas_grid))
              ELSE ! fft forward + backward should be identical
                  CALL coeff_copy(vxc_mol,mgrid_rspace(auxbas_grid))
              ENDIF
         ENDIF
         ! *** Multiply by the grid volume element ratio ***
         IF ( igrid_level /= auxbas_grid ) THEN
            scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                    mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
            mgrid_rspace(igrid_level)%pw%cr3d = &
                                      scale*mgrid_rspace(igrid_level)%pw%cr3d
         END IF
       END DO
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
                                      error=error)
       DEALLOCATE(mgrid_gspace,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL rs_pools_create_rs_vect(rs_pools, rs_v, force_env_section=qs_env%input, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_v(igrid_level)%rs_grid,&
                           mgrid_rspace(igrid_level)%pw,pw2rs)
    ENDDO

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
         error=error)
    DEALLOCATE(mgrid_rspace,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    map_consistent=dft_control%qs_control%map_consistent
    IF (map_consistent) THEN
      eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace ! needs to be consistent with rho_rspace
    ELSE
      eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
    ENDIF

    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf,&
                             maxsgf_set=maxsgf_set)

    NULLIFY ( pabt, habt, workt )
    CALL reallocate(habt,1,maxco,1,maxco,0,nthread)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread)
    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread)

    NULLIFY(atasks,ntasks,tasks,tasks_local,ival)

    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(atasks,1,2,1,max_tasks)
    curr_tasks = max_tasks

    nat_mol = SIZE(atom,1)

    ikind_old = 0
    jkind_old = 0
    DO iat = 1,nat_mol
       ikind = kind(iat)
       iatom = atom(iat)
       ra(1:3) = ratom(1:3,iat)

       IF(ikind /= ikind_old) THEN
         NULLIFY(atomic_kind,orb_basis_a,la_max,la_min,npgfa,nsgfa,&
                 rpgfa,set_radius_a,sphi_a,zeta)
         atomic_kind => atomic_kind_set(ikind)

         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              orb_basis_set=orb_basis_a)

         IF (.NOT.ASSOCIATED(orb_basis_a)) CYCLE
         CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
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
         ikind_old = ikind
       END IF
       IF (.NOT.ASSOCIATED(orb_basis_a)) CYCLE

       DO jat = iat, nat_mol
          jkind = kind(jat)
          jatom = atom(jat)
          rb(1:3) = ratom(1:3,jat)

          IF(jkind /= jkind_old) THEN
            NULLIFY(atomic_kind,orb_basis_b,lb_max,lb_min,npgfb,nsgfb,&
                 rpgfb,set_radius_b,sphi_b,zetb)

            atomic_kind => atomic_kind_set(jkind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                              orb_basis_set=orb_basis_b)

            IF (.NOT.ASSOCIATED(orb_basis_b)) CYCLE

            CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
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
            jkind_old = jkind
          END IF
          IF (.NOT.ASSOCIATED(orb_basis_b)) CYCLE

          ntasks = 0
          tasks = 0

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          NULLIFY(h_block)
          CALL get_block_node(matrix=matrix_h%matrix,&
                              block_row=brow,&
                              block_col=bcol,&
                              BLOCK=h_block)
          CPPrecondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

          NULLIFY(p_block)
          CALL get_block_node(matrix=matrix_p%matrix,&
                              block_row=brow,&
                              block_col=bcol,&
                              BLOCK=p_block)
          CPPrecondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)

          rab2 = 0.0_dp
          DO i = 1,3
            rab(i) = rb(i)-ra(i)
            rab2 = rab2 + rab(i)*rab(i)
          END DO
          dab = SQRT(rab2)

          DO iset=1,nseta

            IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE


            DO jset=1,nsetb

              IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

              DO ipgf=1,npgfa(iset)

                IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE


                DO jpgf=1,npgfb(jset)

                  IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                  zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                  IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                     igrid_level = 1
                  ELSE
                     igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                  ENDIF

                  ntasks(igrid_level) = ntasks(igrid_level) + 1
                  n = ntasks(igrid_level)
                  IF ( n > curr_tasks ) THEN
                    curr_tasks = curr_tasks*mult_tasks
                    CALL reallocate(tasks,1,8,1,curr_tasks,&
                                    1,gridlevel_info%ngrid_levels)
                  END IF

                  tasks (1,n,igrid_level) = n
                  tasks (3,n,igrid_level) = iatom
                  tasks (4,n,igrid_level) = jatom
                  tasks (5,n,igrid_level) = iset
                  tasks (6,n,igrid_level) = jset
                  tasks (7,n,igrid_level) = ipgf
                  tasks (8,n,igrid_level) = jpgf

                END DO  ! jpgf
              END DO  ! ipgf
            END DO  ! jset
          END DO  ! iset

          DO igrid_level = 1, gridlevel_info%ngrid_levels
            num_tasks = ntasks ( igrid_level )
            IF ( num_tasks  > SIZE ( tasks_local, 2 ) ) &
              CALL reallocate(tasks_local,1,2,1,num_tasks)
            IF ( num_tasks > SIZE ( ival, 2 ) ) &
              CALL reallocate(ival,1,6,1,num_tasks)

!$OMP parallel do private(i)
            DO i=1,num_tasks
              tasks_local(1,i) = tasks(1,i,igrid_level)
              tasks_local(2,i) = tasks(2,i,igrid_level)
              ival(1,i) = tasks(3,i,igrid_level)
              ival(2,i) = tasks(4,i,igrid_level)
              ival(3,i) = tasks(5,i,igrid_level)
              ival(4,i) = tasks(6,i,igrid_level)
              ival(5,i) = tasks(7,i,igrid_level)
              ival(6,i) = tasks(8,i,igrid_level)
            END DO

!$OMP parallel do private(i)
            DO i=num_tasks+1,SIZE(tasks_local,2)
              tasks_local(1,i) = 0
              tasks_local(2,i) = 0
            END DO

! fully replicated grids, each processor can process all its tasks
           ! get number of tasks available locally
            npme = SIZE ( tasks_local, 2 )
            DO i = 1, SIZE ( tasks_local, 2)
              IF ( tasks_local ( 1, i ) <= 0 ) THEN
                 npme = i - 1
                 EXIT
              END IF
            END DO

            CALL  pair_get_loop_vars( npme, ival, nset_pairs, atasks )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,itask,iset,jset,ncoa,ncob,sgfa,sgfb) &
!$OMP private(work,hab,pab,istat,ipgf,jpgf,na1,na2,nb1,nb2,scale) &
!$OMP private(force_a,force_b) &
!$OMP shared(iatom,jatom,ra,rb,rab,rab2,brow,bcol,p_block,h_block) &
!$OMP shared(maxco,maxsgf_set,ival,num_tasks,atasks,calculate_forces) &
!$OMP shared(npgfa,npgfb,ncoset,la_max,lb_max,first_sgfa,first_sgfb) &
!$OMP shared(nsgfa,nsgfb,sphi_a,sphi_b,la_min,lb_min,eps_gvg_rspace,zeta,zetb) &
!$OMP shared(rs_v,igrid_level,cube_info,l_info,nthread,nset_pairs) &
!$OMP shared(workt,habt,pabt,my_compute_tau,map_consistent,forces_mol,iat,jat)
            ithread = 0
!$          ithread = omp_get_thread_num()
            hab => habt(:,:,ithread)
            pab => pabt(:,:,ithread)
            work => workt(:,:,ithread)

!$OMP do
            DO ijsets = 1,nset_pairs
              IF (calculate_forces) THEN
                force_a(:) = 0.0_dp
                force_b(:) = 0.0_dp
              END IF

              itask = atasks(1,ijsets)
              iset = ival (3,itask)
              jset = ival (4,itask)
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              ncob = npgfb(jset)*ncoset(lb_max(jset))
              sgfb = first_sgfb(1,jset)

              IF(calculate_forces) THEN
                IF (iatom <= jatom) THEN
                   CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              p_block(sgfa,sgfb),SIZE(p_block,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                   CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_b(1,sgfb),SIZE(sphi_b,1),&
                              0.0_dp,pab(1,1),SIZE(pab,1))
                ELSE
                   CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                              1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                              p_block(sgfb,sgfa),SIZE(p_block,1),&
                              0.0_dp,work(1,1),SIZE(work,1))
                   CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                              1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              0.0_dp,pab(1,1),SIZE(pab,1))
               END IF
              ELSE
                pab = 0._dp
              END IF
              hab = 0._dp
              DO itask = atasks(1,ijsets),atasks(2,ijsets)
                ipgf   = ival (5,itask)
                jpgf   = ival (6,itask)
                na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                na2 = ipgf*ncoset(la_max(iset))
                nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                nb2 = jpgf*ncoset(lb_max(jset))

                IF(calculate_forces) THEN
                  IF (iatom <= jatom) THEN
                     CALL integrate_pgf_product_rspace(&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          ra,rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,pab=pab,o1=na1-1,o2=nb1-1, &
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=.TRUE.,&
                          force_a=force_a,force_b=force_b,ithread=ithread,&
                          compute_tau=my_compute_tau,map_consistent=map_consistent)
                   ELSE
                     CALL integrate_pgf_product_rspace(&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          rb,-rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,pab=pab,o1=nb1-1,o2=na1-1, &
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=.TRUE.,&
                          force_a=force_b,force_b=force_a,ithread=ithread,&
                          compute_tau=my_compute_tau,map_consistent=map_consistent)
                   END IF
                ELSE
                   IF (iatom <= jatom) THEN
                     CALL integrate_pgf_product_rspace(&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          ra,rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,o1=na1-1,o2=nb1-1,&
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=.FALSE.,&
                          ithread=ithread,&
                          compute_tau=my_compute_tau,map_consistent=map_consistent)
                   ELSE
                     CALL integrate_pgf_product_rspace(&
                          lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                          la_max(iset),zeta(ipgf,iset),la_min(iset),&
                          rb,-rab,rab2,rs_v(igrid_level)%rs_grid,&
                          cube_info(igrid_level),l_info,&
                          hab,o1=nb1-1,o2=na1-1,&
                          eps_gvg_rspace=eps_gvg_rspace,&
                          calculate_forces=.FALSE.,&
                          ithread=ithread, &
                          compute_tau=my_compute_tau,map_consistent=map_consistent)
                   END IF
                END IF
              END DO ! itask

              IF (iatom <= jatom) THEN
                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                           -1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
              ELSE
                CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
                           -1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
              END IF

              IF (calculate_forces) THEN
!$OMP critical (qs_integrate_force)
                 forces_mol(:,iat)= forces_mol(:,iat)+ 2.0_dp*force_a(:)
                 IF (iatom /= jatom) THEN
                    forces_mol(:,jat) = forces_mol(:,jat) + 2.0_dp*force_b(:)
                 END IF
!$OMP end critical (qs_integrate_force)
              END IF

            END DO  ! ijsets

!$OMP end parallel

          END DO  ! igrid_level
       END DO  ! jat
    END DO  ! iat

    DEALLOCATE (habt,workt,ntasks,tasks,tasks_local,ival,&
                pabt,atasks,STAT=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    CALL rs_pools_give_back_rs_vect(rs_pools, rs_v)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_mol_potential

!***************************************************************************

  SUBROUTINE pair_get_loop_vars ( npme, ival, nset_pairs, atasks )

    INTEGER, INTENT(IN)                      :: npme
    INTEGER, INTENT(OUT)                     :: nset_pairs
    INTEGER, DIMENSION(6, npme), INTENT(IN)  :: ival
    INTEGER, DIMENSION(:, :), POINTER        :: atasks

    INTEGER                                  :: iset, &
                                                iset_old, itask,  &
                                                jset, jset_old

     IF(SIZE(atasks,2) < npme) CALL reallocate(atasks,1,2,1,npme)

     nset_pairs = 0
     iset_old = 0
     jset_old = 0
     DO itask = 1,npme
       iset = ival(3,itask)
       jset = ival(4,itask)
       IF ( iset /= iset_old .OR. jset /= jset_old ) THEN
         IF(nset_pairs>0) atasks(2,nset_pairs) = itask - 1
         nset_pairs = nset_pairs + 1
         atasks(1,nset_pairs) = itask
         iset_old = iset
         jset_old = jset
        END IF
     END DO  ! itask 
     IF(nset_pairs>0) atasks(2,nset_pairs) = npme

  END SUBROUTINE pair_get_loop_vars

!***************************************************************************
END MODULE kg_gpw_collocate_den
