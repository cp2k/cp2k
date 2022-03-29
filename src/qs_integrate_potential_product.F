!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!>      JGH [26.06.15] modification to allow for k-points
!> \author Matthias Krack (03.04.2001)
! **************************************************************************************************
MODULE qs_integrate_potential_product
   USE admm_types,                      ONLY: admm_type
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE grid_api,                        ONLY: grid_integrate_task_list,&
                                              integrate_pgf_product
   USE input_constants,                 ONLY: do_admm_exch_scaling_merlot
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              rs_grid_release,&
                                              rs_grid_retain
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_methods,               ONLY: rs_copy_to_buffer,&
                                              rs_copy_to_matrices,&
                                              rs_distribute_matrix,&
                                              rs_gather_matrices,&
                                              rs_scatter_matrices
   USE task_list_types,                 ONLY: atom_pair_type,&
                                              task_list_type,&
                                              task_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_product'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_v_rspace
   PUBLIC :: integrate_v_dbasis

CONTAINS

! **************************************************************************************************
!> \brief Integrate a potential v_rspace over the derivatives of the basis functions
!>         < da/dR | V | b > + < a | V | db/dR >
!>        Adapted from the old version of integrate_v_rspace (ED)
!> \param v_rspace ...
!> \param matrix_vhxc_dbasis ...
!> \param matrix_p ...
!> \param qs_env ...
!> \param lambda The atom index.
! **************************************************************************************************
   SUBROUTINE integrate_v_dbasis(v_rspace, matrix_vhxc_dbasis, matrix_p, qs_env, lambda)
      TYPE(pw_p_type)                                    :: v_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vhxc_dbasis
      TYPE(dbcsr_type), POINTER                          :: matrix_p
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: lambda

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_dbasis'

      INTEGER :: bcol, brow, group, handle, i, iatom, igrid_level, ikind, ikind_old, ilevel, img, &
         ipair, ipgf, ipgf_new, iset, iset_new, iset_old, itask, ithread, jatom, jkind, jkind_old, &
         jpgf, jpgf_new, jset, jset_new, jset_old, maxco, maxpgf, maxset, maxsgf_set, na1, na2, &
         natom, nb1, nb2, ncoa, ncob, nimages, nkind, nseta, nsetb, nthread, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: block_touched
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL :: atom_pair_changed, atom_pair_done, dh_duplicated, distributed_grids, found, &
         my_compute_tau, my_gapw, new_set_pair_coming, pab_required, scatter, use_subpatch
      REAL(KIND=dp)                                      :: eps_rho_rspace, f, prefactor, radius, &
                                                            scalef, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra, rab, rab_inv, rb, &
                                                            rp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, hab, p_block, pab, rpgfa, &
                                                            rpgfb, sphi_a, sphi_b, work, zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: habt, hadb, hdab, pabt, workt
      REAL(kind=dp), DIMENSION(:, :, :, :), POINTER      :: hadbt, hdabt
      TYPE(atom_pair_type), DIMENSION(:), POINTER        :: atom_pair_recv, atom_pair_send
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: vhxc_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(task_type), DIMENSION(:), POINTER             :: tasks

      CALL timeset(routineN, handle)
      NULLIFY (pw_env)

      ! get the task lists
      my_gapw = .FALSE.
      CALL get_qs_env(qs_env=qs_env, &
                      task_list=task_list, &
                      task_list_soft=task_list_soft)
      CPASSERT(ASSOCIATED(task_list))

      ! the information on the grids is provided through pw_env
      ! pw_env has to be the parent env for the potential grid (input)
      ! there is an option to provide an external grid
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)

      ! *** assign from pw_env
      gridlevel_info => pw_env%gridlevel_info

      ! get all the general information on the system we are working on
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      natom=natom, &
                      dft_control=dft_control, &
                      particle_set=particle_set)

      ! *** set up the rs multi-grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho)
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         distributed_grids = rs_rho(igrid_level)%rs_grid%desc%distributed
      END DO
      ! get mpi group from rs_rho
      group = rs_rho(1)%rs_grid%desc%group

      ! transform the potential on the rs_multigrids
      CALL potential_pw2rs(rs_rho, v_rspace, pw_env)

      nkind = SIZE(qs_kind_set)

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type="ORB")

      ! short cuts to task list variables
      tasks => task_list%tasks
      atom_pair_send => task_list%atom_pair_send
      atom_pair_recv => task_list%atom_pair_recv

      ! needs to be consistent with rho_rspace
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      !   *** Initialize working density matrix ***
      ! distributed rs grids require a matrix that will be changed
      ! whereas this is not the case for replicated grids
      ALLOCATE (deltap(dft_control%nimages))
      IF (distributed_grids) THEN
         ! this matrix has no strict sparsity pattern in parallel
         ! deltap%sparsity_id=-1
         CALL dbcsr_copy(deltap(1)%matrix, matrix_p, name="DeltaP")
      ELSE
         deltap(1)%matrix => matrix_p
      END IF
      nthread = 1
!$    nthread = omp_get_max_threads()

      !   *** Allocate work storage ***
      NULLIFY (pabt, habt, workt)
      ALLOCATE (habt(maxco, maxco, 0:nthread))
      ALLOCATE (workt(maxco, maxsgf_set, 0:nthread))
      ALLOCATE (hdabt(3, maxco, maxco, 0:nthread))
      ALLOCATE (hadbt(3, maxco, maxco, 0:nthread))
      ALLOCATE (pabt(maxco, maxco, 0:nthread))

      IF (distributed_grids) THEN
         CALL rs_distribute_matrix(rs_descs, deltap, atom_pair_send, atom_pair_recv, &
                                   nimages, scatter=.TRUE.)
      END IF

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(workt,habt,hdabt,hadbt,pabt,tasks,particle_set,natom,maxset), &
!$OMP SHARED(maxpgf,my_gapw,matrix_vhxc_dbasis,deltap), &
!$OMP SHARED(pab_required,ncoset,rs_rho,my_compute_tau), &
!$OMP SHARED(eps_rho_rspace,force,cell), &
!$OMP SHARED(gridlevel_info,task_list,block_touched,nthread,qs_kind_set), &
!$OMP SHARED(nimages,lambda, dh_duplicated), &
!$OMP PRIVATE(ithread,work,hab,hdab,hadb,pab,iset_old,jset_old), &
!$OMP PRIVATE(ikind_old,jkind_old,iatom,jatom,iset,jset,ikind,jkind,ilevel,ipgf,jpgf), &
!$OMP PRIVATE(img,brow,bcol,orb_basis_set,first_sgfa,la_max,la_min,npgfa,nseta,nsgfa), &
!$OMP PRIVATE(rpgfa,set_radius_a,sphi_a,zeta,first_sgfb,lb_max,lb_min,npgfb), &
!$OMP PRIVATE(nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,found), &
!$OMP PRIVATE(force_a,force_b,my_virial_a,my_virial_b,atom_pair_changed,h_block, vhxc_block), &
!$OMP PRIVATE(p_block,ncoa,sgfa,ncob,sgfb,rab,ra,rb,rp,zetp,f,prefactor,radius,igrid_level), &
!$OMP PRIVATE(na1,na2,nb1,nb2,use_subpatch,rab_inv,new_set_pair_coming,atom_pair_done), &
!$OMP PRIVATE(iset_new,jset_new,ipgf_new,jpgf_new,scalef), &
!$OMP PRIVATE(itask)

      IF (.NOT. ALLOCATED(vhxc_block)) ALLOCATE (vhxc_block(3))

      ithread = 0
!$    ithread = omp_get_thread_num()
      work => workt(:, :, ithread)
      hab => habt(:, :, ithread)
      pab => pabt(:, :, ithread)
      hdab => hdabt(:, :, :, ithread)
      hadb => hadbt(:, :, :, ithread)

      iset_old = -1; jset_old = -1
      ikind_old = -1; jkind_old = -1

      ! Here we loop over gridlevels first, finalising the matrix after each grid level is
      ! completed.  On each grid level, we loop over atom pairs, which will only access
      ! a single block of each matrix, so with OpenMP, each matrix block is only touched
      ! by a single thread for each grid level
      loop_gridlevels: DO igrid_level = 1, gridlevel_info%ngrid_levels
!$OMP BARRIER
!$OMP DO schedule (dynamic, MAX(1,task_list%npairs(igrid_level)/(nthread*50)))
         loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
         loop_tasks: DO itask = task_list%taskstart(ipair, igrid_level), task_list%taskstop(ipair, igrid_level)
            ilevel = tasks(itask)%grid_level
            img = tasks(itask)%image
            iatom = tasks(itask)%iatom
            jatom = tasks(itask)%jatom
            iset = tasks(itask)%iset
            jset = tasks(itask)%jset
            ipgf = tasks(itask)%ipgf
            jpgf = tasks(itask)%jpgf

            ! At the start of a block of tasks, get atom data (and kind data, if needed)
            IF (itask .EQ. task_list%taskstart(ipair, igrid_level)) THEN

               ikind = particle_set(iatom)%atomic_kind%kind_number
               jkind = particle_set(jatom)%atomic_kind%kind_number

               ra(:) = pbc(particle_set(iatom)%r, cell)

               IF (iatom <= jatom) THEN
                  brow = iatom
                  bcol = jatom
               ELSE
                  brow = jatom
                  bcol = iatom
               END IF

               IF (ikind .NE. ikind_old) THEN
                  CALL get_qs_kind(qs_kind_set(ikind), &
                                   softb=my_gapw, &
                                   basis_set=orb_basis_set, basis_type="ORB")

                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfa, &
                                         lmax=la_max, &
                                         lmin=la_min, &
                                         npgf=npgfa, &
                                         nset=nseta, &
                                         nsgf_set=nsgfa, &
                                         pgf_radius=rpgfa, &
                                         set_radius=set_radius_a, &
                                         sphi=sphi_a, &
                                         zet=zeta)
               END IF

               IF (jkind .NE. jkind_old) THEN
                  CALL get_qs_kind(qs_kind_set(jkind), &
                                   softb=my_gapw, &
                                   basis_set=orb_basis_set, basis_type="ORB")
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         first_sgf=first_sgfb, &
                                         lmax=lb_max, &
                                         lmin=lb_min, &
                                         npgf=npgfb, &
                                         nset=nsetb, &
                                         nsgf_set=nsgfb, &
                                         pgf_radius=rpgfb, &
                                         set_radius=set_radius_b, &
                                         sphi=sphi_b, &
                                         zet=zetb)

               END IF

               DO i = 1, 3
                  NULLIFY (vhxc_block(i)%block)
                  CALL dbcsr_get_block_p(matrix_vhxc_dbasis(i)%matrix, brow, bcol, vhxc_block(i)%block, found)
                  CPASSERT(found)
               END DO

               CALL dbcsr_get_block_p(matrix=deltap(img)%matrix, &
                                      row=brow, col=bcol, BLOCK=p_block, found=found)
               CPASSERT(found)

               ikind_old = ikind
               jkind_old = jkind

               atom_pair_changed = .TRUE.

            ELSE

               atom_pair_changed = .FALSE.

            END IF

            IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               IF (iatom <= jatom) THEN
                  work(1:ncoa, 1:nsgfb(jset)) = MATMUL(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1), &
                                                       p_block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1))
                  pab(1:ncoa, 1:ncob) = MATMUL(work(1:ncoa, 1:nsgfb(jset)), TRANSPOSE(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1)))
               ELSE
                  work(1:ncob, 1:nsgfa(iset)) = MATMUL(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1), &
                                                       p_block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1))
                  pab(1:ncob, 1:ncoa) = MATMUL(work(1:ncob, 1:nsgfa(iset)), TRANSPOSE(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1)))
               END IF

               IF (iatom <= jatom) THEN
                  hab(1:ncoa, 1:ncob) = 0._dp
                  hdab(:, 1:ncoa, 1:ncob) = 0._dp
                  hadb(:, 1:ncoa, 1:ncob) = 0._dp
               ELSE
                  hab(1:ncob, 1:ncoa) = 0._dp
                  hdab(:, 1:ncob, 1:ncoa) = 0._dp
                  hadb(:, 1:ncob, 1:ncoa) = 0._dp
               END IF

               iset_old = iset
               jset_old = jset

            END IF

            rab = tasks(itask)%rab
            rb(:) = ra(:) + rab(:)
            zetp = zeta(ipgf, iset) + zetb(jpgf, jset)

            f = zetb(jpgf, jset)/zetp
            rp(:) = ra(:) + f*rab(:)
            prefactor = EXP(-zeta(ipgf, iset)*f*DOT_PRODUCT(rab, rab))
            radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                              lb_min=lb_min(jset), lb_max=lb_max(jset), &
                                              ra=ra, rb=rb, rp=rp, &
                                              zetp=zetp, eps=eps_rho_rspace, &
                                              prefactor=prefactor, cutoff=1.0_dp)

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))
            nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
            nb2 = jpgf*ncoset(lb_max(jset))

            ! check whether we need to use fawzi's generalised collocation scheme
            IF (rs_rho(igrid_level)%rs_grid%desc%distributed) THEN
               !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
               IF (tasks(itask)%dist_type .EQ. 2) THEN
                  use_subpatch = .TRUE.
               ELSE
                  use_subpatch = .FALSE.
               END IF
            ELSE
               use_subpatch = .FALSE.
            END IF

            IF (iatom <= jatom) THEN
               IF (iatom == lambda) &
                  CALL integrate_pgf_product( &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  ra, rab, rs_rho(igrid_level)%rs_grid, &
                  hab, o1=na1 - 1, o2=nb1 - 1, &
                  radius=radius, &
                  calculate_forces=.TRUE., &
                  compute_tau=.FALSE., &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern, &
                  hdab=hdab, pab=pab)
               IF (jatom == lambda) &
                  CALL integrate_pgf_product( &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  ra, rab, rs_rho(igrid_level)%rs_grid, &
                  hab, o1=na1 - 1, o2=nb1 - 1, &
                  radius=radius, &
                  calculate_forces=.TRUE., &
                  compute_tau=.FALSE., &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern, &
                  hadb=hadb, pab=pab)
            ELSE
               rab_inv = -rab
               IF (iatom == lambda) &
                  CALL integrate_pgf_product( &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  rb, rab_inv, rs_rho(igrid_level)%rs_grid, &
                  hab, o1=nb1 - 1, o2=na1 - 1, &
                  radius=radius, &
                  calculate_forces=.TRUE., &
                  force_a=force_b, force_b=force_a, &
                  compute_tau=.FALSE., &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern, &
                  hadb=hadb, pab=pab)
               IF (jatom == lambda) &
                  CALL integrate_pgf_product( &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  rb, rab_inv, rs_rho(igrid_level)%rs_grid, &
                  hab, o1=nb1 - 1, o2=na1 - 1, &
                  radius=radius, &
                  calculate_forces=.TRUE., &
                  force_a=force_b, force_b=force_a, &
                  compute_tau=.FALSE., &
                  use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern, &
                  hdab=hdab, pab=pab)
            END IF

            new_set_pair_coming = .FALSE.
            atom_pair_done = .FALSE.
            IF (itask < task_list%taskstop(ipair, igrid_level)) THEN
               ilevel = tasks(itask + 1)%grid_level
               img = tasks(itask + 1)%image
               iatom = tasks(itask + 1)%iatom
               jatom = tasks(itask + 1)%jatom
               iset_new = tasks(itask + 1)%iset
               jset_new = tasks(itask + 1)%jset
               ipgf_new = tasks(itask + 1)%ipgf
               jpgf_new = tasks(itask + 1)%jpgf
               IF (iset_new .NE. iset .OR. jset_new .NE. jset) THEN
                  new_set_pair_coming = .TRUE.
               END IF
            ELSE
               ! do not forget the last block
               new_set_pair_coming = .TRUE.
               atom_pair_done = .TRUE.
            END IF

            IF (new_set_pair_coming) THEN

               DO i = 1, 3
                  hdab(i, :, :) = hdab(i, :, :) + hadb(i, :, :)
                  IF (iatom <= jatom) THEN
                     work(1:ncoa, 1:nsgfb(jset)) = MATMUL(hdab(i, 1:ncoa, 1:ncob), sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1))
                     vhxc_block(i)%block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) = &
                        vhxc_block(i)%block(sgfa:sgfa + nsgfa(iset) - 1, sgfb:sgfb + nsgfb(jset) - 1) + &
                        MATMUL(TRANSPOSE(sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1)), work(1:ncoa, 1:nsgfb(jset)))
                  ELSE
                     work(1:ncob, 1:nsgfa(iset)) = MATMUL(hdab(i, 1:ncob, 1:ncoa), sphi_a(1:ncoa, sgfa:sgfa + nsgfa(iset) - 1))
                     vhxc_block(i)%block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) = &
                        vhxc_block(i)%block(sgfb:sgfb + nsgfb(jset) - 1, sgfa:sgfa + nsgfa(iset) - 1) + &
                        MATMUL(TRANSPOSE(sphi_b(1:ncob, sgfb:sgfb + nsgfb(jset) - 1)), work(1:ncob, 1:nsgfa(iset)))
                  END IF
               END DO
            END IF  ! new_set_pair_coming

         END DO loop_tasks
         END DO loop_pairs
!$OMP END DO

         DO i = 1, 3
            CALL dbcsr_finalize(matrix_vhxc_dbasis(i)%matrix)
         END DO

      END DO loop_gridlevels

!$OMP END PARALLEL

      IF (distributed_grids) THEN
         ! Reconstruct H matrix if using distributed RS grids
         ! note send and recv direction reversed WRT collocate
         scatter = .FALSE.
         CALL rs_distribute_matrix(rs_descs, matrix_vhxc_dbasis, atom_pair_recv, atom_pair_send, &
                                   dft_control%nimages, scatter=.FALSE.)
      END IF

      IF (distributed_grids) THEN
         CALL dbcsr_deallocate_matrix_set(deltap)
      ELSE
         DO img = 1, dft_control%nimages
            NULLIFY (deltap(img)%matrix)
         END DO
         DEALLOCATE (deltap)
      END IF

      DEALLOCATE (pabt, habt, workt, hdabt, hadbt)

      IF (ASSOCIATED(rs_rho)) THEN
         DO i = 1, SIZE(rs_rho)
            CALL rs_grid_release(rs_rho(i)%rs_grid)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE integrate_v_dbasis

! **************************************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \param v_rspace ...
!> \param hmat ...
!> \param hmat_kp ...
!> \param pmat ...
!> \param pmat_kp ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param force_adm whether force of in aux. dens. matrix is calculated
!> \param ispin ...
!> \param compute_tau ...
!> \param gapw ...
!> \param basis_type ...
!> \param pw_env_external ...
!> \param task_list_external ...
!> \par History
!>      IAB (29-Apr-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      Some refactoring, get priorities for options correct (JGH, 04.2014)
!>      Added options to allow for k-points
!>      For a smooth transition we allow for old and new (vector) matrices (hmat, pmat) (JGH, 06.2015)
!> \note
!>     integrates a given potential (or other object on a real
!>     space grid) = v_rspace using a multi grid technique (mgrid_*)
!>     over the basis set producing a number for every element of h
!>     (should have the same sparsity structure of S)
!>     additional screening is available using the magnitude of the
!>     elements in p (? I'm not sure this is a very good idea)
!>     this argument is optional
!>     derivatives of these matrix elements with respect to the ionic
!>     coordinates can be computed as well
! **************************************************************************************************
   SUBROUTINE integrate_v_rspace(v_rspace, hmat, hmat_kp, pmat, pmat_kp, &
                                 qs_env, calculate_forces, force_adm, ispin, &
                                 compute_tau, gapw, basis_type, pw_env_external, task_list_external)

      TYPE(pw_p_type)                                    :: v_rspace
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: hmat_kp
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: pmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pmat_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_adm
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, gapw
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace'

      CHARACTER(len=default_string_length)               :: my_basis_type
      INTEGER                                            :: atom_a, group, handle, i, iatom, &
                                                            igrid_level, ikind, img, maxco, &
                                                            maxsgf_set, natoms, nimages, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: calculate_virial, distributed_grids, &
                                                            do_kp, my_compute_tau, my_force_adm, &
                                                            my_gapw, pab_required
      REAL(KIND=dp)                                      :: admm_scal_fac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: forces_array
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_matrix
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap, dhmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! we test here if the provided operator matrices are consistent
      CPASSERT(PRESENT(hmat) .OR. PRESENT(hmat_kp))
      do_kp = .FALSE.
      IF (PRESENT(hmat_kp)) do_kp = .TRUE.
      IF (PRESENT(pmat)) THEN
         CPASSERT(PRESENT(hmat))
      ELSE IF (PRESENT(pmat_kp)) THEN
         CPASSERT(PRESENT(hmat_kp))
      END IF

      NULLIFY (pw_env, admm_env)

      ! this routine works in two modes:
      ! normal mode : <a| V | b>
      ! tau mode    : < nabla a| V | nabla b>
      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau

      my_force_adm = .FALSE.
      IF (PRESENT(force_adm)) my_force_adm = force_adm

      ! this sets the basis set to be used. GAPW(==soft basis) overwrites basis_type
      ! default is "ORB"
      my_gapw = .FALSE.
      IF (PRESENT(gapw)) my_gapw = gapw
      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      ! get the task lists
      ! task lists have to be in sync with basis sets
      ! there is an option to provide the task list from outside (not through qs_env)
      ! outside option has highest priority
      CALL get_qs_env(qs_env=qs_env, &
                      task_list=task_list, &
                      task_list_soft=task_list_soft)
      IF (.NOT. my_basis_type == "ORB") THEN
         CPASSERT(PRESENT(task_list_external))
      END IF
      IF (my_gapw) task_list => task_list_soft
      IF (PRESENT(task_list_external)) task_list => task_list_external
      CPASSERT(ASSOCIATED(task_list))

      ! the information on the grids is provided through pw_env
      ! pw_env has to be the parent env for the potential grid (input)
      ! there is an option to provide an external grid
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external

      ! get all the general information on the system we are working on
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      natom=natoms, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      force=force, &
                      virial=virial)

      admm_scal_fac = 1.0_dp
      IF (my_force_adm) THEN
         CALL get_qs_env(qs_env=qs_env, admm_env=admm_env)
         ! Calculate bare scaling of force according to Merlot, 1. IF: ADMMP, 2. IF: ADMMS,
         IF ((.NOT. admm_env%charge_constrain) .AND. &
             (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = admm_env%gsi(ispin)**2
         ELSE IF (admm_env%charge_constrain .AND. &
                  (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)
         END IF
      END IF

      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO

      ! get mpi group from rs_v
      group = rs_v(1)%rs_grid%desc%group

      ! assign from pw_env
      gridlevel_info => pw_env%gridlevel_info

      ! transform the potential on the rs_multigrids
      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         CPASSERT(do_kp)
      END IF
      nkind = SIZE(qs_kind_set)
      calculate_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer) .AND. calculate_forces
      pab_required = (PRESENT(pmat) .OR. PRESENT(pmat_kp)) .AND. calculate_forces

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)

      distributed_grids = .FALSE.
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF (rs_v(igrid_level)%rs_grid%desc%distributed) THEN
            distributed_grids = .TRUE.
         END IF
      END DO

      ALLOCATE (forces_array(3, natoms))

      IF (pab_required) THEN
         ! initialize the working pmat structures
         ALLOCATE (deltap(nimages))
         IF (do_kp) THEN
            DO img = 1, nimages
               deltap(img)%matrix => pmat_kp(img)%matrix
            END DO
         ELSE
            deltap(1)%matrix => pmat%matrix
         END IF

         ! Distribute matrix blocks.
         IF (distributed_grids) THEN
            CALL rs_scatter_matrices(deltap, task_list%pab_buffer, task_list, group)
         ELSE
            CALL rs_copy_to_buffer(deltap, task_list%pab_buffer, task_list)
         END IF
         DEALLOCATE (deltap)
      END IF

      ! Map all tasks from the grids
      CALL grid_integrate_task_list(task_list=task_list%grid_task_list, &
                                    compute_tau=my_compute_tau, &
                                    calculate_forces=calculate_forces, &
                                    calculate_virial=calculate_virial, &
                                    pab_blocks=task_list%pab_buffer, &
                                    rs_grids=rs_v, &
                                    hab_blocks=task_list%hab_buffer, &
                                    forces=forces_array, &
                                    virial=virial_matrix)

      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natoms), kind_of(natoms))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
!$OMP PARALLEL DO DEFAULT(NONE)  PRIVATE(atom_a, ikind) &
!$OMP             SHARED(natoms, force, forces_array, atom_of_kind, kind_of, admm_scal_fac)
         DO iatom = 1, natoms
            atom_a = atom_of_kind(iatom)
            ikind = kind_of(iatom)
            force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + admm_scal_fac*forces_array(:, iatom)
         END DO
!$OMP END PARALLEL DO
         DEALLOCATE (atom_of_kind, kind_of)
      END IF

      IF (calculate_virial) THEN
         virial%pv_virial = virial%pv_virial + admm_scal_fac*virial_matrix
      END IF

      ! Gather all matrix images into a single array.
      ALLOCATE (dhmat(nimages))
      IF (PRESENT(hmat_kp)) THEN
         CPASSERT(.NOT. PRESENT(hmat))
         DO img = 1, nimages
            dhmat(img)%matrix => hmat_kp(img)%matrix
         END DO
      ELSE
         CPASSERT(PRESENT(hmat) .AND. nimages == 1)
         dhmat(1)%matrix => hmat%matrix
      END IF

      ! Distribute matrix blocks.
      IF (distributed_grids) THEN
         CALL rs_gather_matrices(task_list%hab_buffer, dhmat, task_list, group)
      ELSE
         CALL rs_copy_to_matrices(task_list%hab_buffer, dhmat, task_list)
      END IF
      DEALLOCATE (dhmat)

      IF (ASSOCIATED(rs_v)) THEN
         DO i = 1, SIZE(rs_v)
            CALL rs_grid_release(rs_v(i)%rs_grid)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE integrate_v_rspace

END MODULE qs_integrate_potential_product
