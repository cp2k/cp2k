!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb contributions in DFTB
!> \author JGH
! *****************************************************************************
MODULE qs_dftb_coulomb

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_multiply_local
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dgs,                             ONLY: dg_sum_patch,&
                                             dg_sum_patch_force_1d,&
                                             dg_sum_patch_force_3d
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erfc
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE f77_blas
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pme_tools,                       ONLY: get_center,&
                                             set_list
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: get_pw_grid_info
  USE pw_methods,                      ONLY: pw_integral_a2b,&
                                             pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_rebuild,&
                                             pw_poisson_solve
  USE pw_poisson_types,                ONLY: greens_fn_type,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_p_type, &
       realspace_grid_type, rs2pw, rs_grid_create, rs_grid_release, &
       rs_grid_set_box, rs_grid_zero, rs_pw_transfer
  USE spme,                            ONLY: get_patch
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_coulomb'

  PUBLIC :: build_dftb_coulomb

CONTAINS

! *****************************************************************************
  SUBROUTINE build_dftb_coulomb(ks_env,qs_env,ks_matrix,rho,mcharge,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(dp), DIMENSION(:)                   :: mcharge
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_i, atom_j, blk, &
                                                ewald_type, handle, i, ia, &
                                                iatom, ikind, istat, jatom, &
                                                jkind, natom, nmat
    INTEGER, DIMENSION(3)                    :: periodic
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    LOGICAL                                  :: do_ewald, failure, found, tr, &
                                                use_virial
    REAL(KIND=dp)                            :: alpha, deth, dr, fi, gmij
    REAL(KIND=dp), DIMENSION(3)              :: fij, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsblock, gblock, gmcharge, &
                                                ksblock, ksblock_2, mcin, &
                                                mcout, pblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: gamma_matrix, matrix_p, &
                                                matrix_s
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: n_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(gamma_matrix, matrix_p, matrix_s)

    failure=.FALSE.
    use_virial=.FALSE.

    IF ( calculate_forces ) THEN
      nmat = 4
    ELSE
      nmat = 1
    END IF

    natom = SIZE(mcharge)
    ALLOCATE(gmcharge(natom,nmat),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    gmcharge = 0._dp

    CALL get_qs_env(qs_env=qs_env,gamma_matrix=gamma_matrix,&
                    particle_set=particle_set,cell=cell,&
                    virial=virial,error=error)

    IF ( calculate_forces ) THEN
      use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    END IF

    ALLOCATE(mcin(natom,1),mcout(natom,1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mcin(:,1) = mcharge(:)

    DO i=1,nmat
      mcout = 0._dp
      CALL cp_dbcsr_multiply_local(gamma_matrix(i)%matrix,mcin,mcout,ncol=1,error=error)
      gmcharge(:,i) = mcout(:,1)
    END DO

    DEALLOCATE(mcin,mcout,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF ( calculate_forces .AND. use_virial ) THEN
       CALL cp_dbcsr_iterator_start(iter, gamma_matrix(1)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          NULLIFY(gblock)
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, gblock, tr, blk)
          DO i=1,3
             NULLIFY(gblock)
             CALL cp_dbcsr_get_block_p(matrix=gamma_matrix(i+1)%matrix,&
                  row=iatom,col=jatom,block=gblock,tr=tr,found=found)
             fij(i) = gblock(1,1)*mcharge(iatom)*mcharge(jatom)
          END DO
          rij = particle_set(iatom)%r - particle_set(jatom)%r
          rij = pbc(rij,cell)
          CALL virial_pair_force ( virial%pv_virial, 1._dp, fij, rij, error)
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    END IF

    ! 1/R contribution
    do_ewald = qs_env%dft_control%qs_control%dftb_control%do_ewald
    IF (do_ewald) THEN
      ! Ewald sum
      NULLIFY(ewald_env,ewald_pw)
      CALL get_qs_env(qs_env=qs_env,&
                      ewald_env=ewald_env,ewald_pw=ewald_pw,&
                      error=error)
      CALL get_cell(cell=cell,periodic=periodic,deth=deth)
      CALL ewald_env_get(ewald_env,alpha=alpha,ewald_type=ewald_type,error=error)
      CALL get_qs_env(qs_env=qs_env,sac_3c=n_list,error=error)
      CALL dftb_ewald_overlap(gmcharge,mcharge,alpha,n_list,&
                              virial,use_virial,error=error)
      SELECT CASE(ewald_type)
      CASE DEFAULT
        CALL cp_unimplemented_error(routineP,"Invalid Ewald type",error)
      CASE(do_ewald_none)
        CPErrorMessage(cp_failure_level,routineP,"Not allowed with DFTB",error)
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE(do_ewald_ewald)
        CALL cp_unimplemented_error(routineP,&
             "Standard Ewald not implemented in DFTB",error)
      CASE(do_ewald_pme)
        CALL cp_unimplemented_error(routineP,&
             "PME not implemented in DFTB",error)
      CASE(do_ewald_spme)
        CALL dftb_spme_evaluate (ewald_env,ewald_pw,particle_set,cell,&
          gmcharge,mcharge,calculate_forces,virial,use_virial,error)
      END SELECT
    ELSE
      ! direct sum
      CALL get_qs_env(qs_env=qs_env,&
                      local_particles=local_particles,error=error)
      DO ikind=1,SIZE(local_particles%n_el)
        DO ia=1,local_particles%n_el(ikind)
          iatom=local_particles%list(ikind)%array(ia)
          DO jatom=1,iatom-1
            rij = particle_set(iatom)%r - particle_set(jatom)%r
            rij = pbc(rij,cell)
            dr = SQRT(SUM(rij(:)**2))
            gmcharge(iatom,1)=gmcharge(iatom,1)+mcharge(jatom)/dr
            gmcharge(jatom,1)=gmcharge(jatom,1)+mcharge(iatom)/dr
            DO i=2,nmat
              gmcharge(iatom,i)=gmcharge(iatom,i)+rij(i-1)*mcharge(jatom)/dr**3
              gmcharge(jatom,i)=gmcharge(jatom,i)-rij(i-1)*mcharge(iatom)/dr**3
            END DO
          END DO
        END DO
      END DO
      CPPostcondition(.NOT.use_virial,cp_failure_level,routineP,error,failure)
    END IF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    force=force,para_env=para_env,error=error)
    CALL mp_sum(gmcharge(:,1),para_env%group)

    IF (do_ewald) THEN
      ! add self charge interaction and background charge contribution
      gmcharge(:,1) = gmcharge(:,1) - 2._dp*alpha*oorootpi*mcharge(:)
      IF ( ANY (periodic(:)==1) ) THEN
        gmcharge(:,1) = gmcharge(:,1) - pi / alpha**2 / deth
      END IF
    END IF

    energy%hartree = energy%hartree + 0.5_dp*SUM(mcharge(:)*gmcharge(:,1))
    IF ( calculate_forces ) THEN
      ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of,&
                               atom_of_kind=atom_of_kind)

      gmcharge(:,2)=gmcharge(:,2)*mcharge(:)
      gmcharge(:,3)=gmcharge(:,3)*mcharge(:)
      gmcharge(:,4)=gmcharge(:,4)*mcharge(:)
      DO iatom=1,natom
        ikind  = kind_of(iatom)
        atom_i = atom_of_kind(iatom)
        force(ikind)%rho_elec(1,atom_i) =&
            force(ikind)%rho_elec(1,atom_i) - gmcharge(iatom,2)
        force(ikind)%rho_elec(2,atom_i) =&
            force(ikind)%rho_elec(2,atom_i) - gmcharge(iatom,3)
        force(ikind)%rho_elec(3,atom_i) =&
            force(ikind)%rho_elec(3,atom_i) - gmcharge(iatom,4)
      END DO
    END IF

    IF ( .NOT. just_energy ) THEN
      CALL get_qs_env(qs_env=qs_env,&
                      matrix_s=matrix_s,error=error)
      matrix_p => rho%rho_ao

      IF ( calculate_forces .AND. SIZE(matrix_p) == 2) THEN
        CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
             alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
      END IF

      CALL cp_dbcsr_iterator_start(iter, ks_matrix(1)%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, ksblock, tr, blk)
         NULLIFY(sblock,ksblock_2)
         IF (SIZE(ks_matrix,1)>1) CALL cp_dbcsr_get_block_p(matrix=ks_matrix(2)%matrix,&
              row=iatom,col=jatom,block=ksblock_2,tr=tr,found=found)
         CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix,&
              row=iatom,col=jatom,block=sblock,tr=tr,found=found)
         gmij = 0.5_dp*(gmcharge(iatom,1)+gmcharge(jatom,1))
         ksblock = ksblock - gmij*sblock
         IF (SIZE(ks_matrix,1)>1) ksblock_2 = ksblock_2 - gmij*sblock
         IF ( calculate_forces ) THEN
            ikind  = kind_of(iatom)
            atom_i = atom_of_kind(iatom)
            jkind  = kind_of(jatom)
            atom_j = atom_of_kind(jatom)
            NULLIFY(pblock)
            CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
                 row=iatom,col=jatom,block=pblock,tr=tr,found=found)
            DO i=1,3
               NULLIFY(dsblock)
               CALL cp_dbcsr_get_block_p(matrix=matrix_s(1+i)%matrix,&
                    row=iatom,col=jatom,block=dsblock,tr=tr,found=found)
               fi = -gmij*SUM(pblock*dsblock)
               force(ikind)%rho_elec(i,atom_i) =&
                    force(ikind)%rho_elec(i,atom_i) + fi
               force(jkind)%rho_elec(i,atom_j) =&
                    force(jkind)%rho_elec(i,atom_j) - fi
               fij(i) = fi
            END DO
            IF (use_virial) THEN
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               rij = pbc(rij,cell)
               CALL virial_pair_force ( virial%pv_virial, -1._dp, fij, rij, error)
            END IF
         END IF
      END DO
      CALL cp_dbcsr_iterator_stop(iter)
      IF ( calculate_forces .AND. SIZE(matrix_p) == 2) THEN
         CALL cp_dbcsr_add(matrix_p(1)%matrix,matrix_p(2)%matrix,&
                        alpha_scalar=1.0_dp,beta_scalar=-1.0_dp,error=error)
      END IF
    END IF

    IF ( calculate_forces ) THEN
      DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(gmcharge,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE build_dftb_coulomb

! *****************************************************************************
  SUBROUTINE dftb_spme_evaluate (ewald_env,ewald_pw,particle_set,box,&
       gmcharge,mcharge,calculate_forces,virial,use_virial,error)

    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: box
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: gmcharge
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(inout)                          :: mcharge
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(virial_type), POINTER               :: virial
    LOGICAL, INTENT(in)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_spme_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, ipart, j, &
                                                n, npart, o_spline, p1, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: center
    INTEGER, DIMENSION(3)                    :: npts
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, dvols, fat(3), ffa, &
                                                fint, vgc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rhos
    REAL(KIND=dp), DIMENSION(3, 3)           :: f_stress, h_stress
    TYPE(greens_fn_type), POINTER            :: green
    TYPE(pw_grid_type), POINTER              :: grid_spme
    TYPE(pw_p_type), DIMENSION(3)            :: dphi_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_type), POINTER                   :: phi_g, rhob_g, rhob_r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    TYPE(realspace_grid_type), POINTER       :: rden, rpot

!------------------------------------------------------------------------------
!-------------- INITIALISATION ---------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL ewald_env_get ( ewald_env, alpha=alpha, o_spline = o_spline, group = group, &
         error=error)
    NULLIFY(green, poisson_env, pw_pool)
    CALL ewald_pw_get ( ewald_pw, pw_big_pool=pw_pool, rs_desc=rs_desc, &
         poisson_env=poisson_env)
    CALL pw_poisson_rebuild(poisson_env,error=error)
    green => poisson_env%green_fft
    grid_spme => pw_pool % pw_grid

    CALL get_pw_grid_info(grid_spme,dvol=dvols,npts=npts,error=error)

    npart = SIZE ( particle_set )

    n = o_spline
    ALLOCATE ( rhos ( n, n, n ), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL rs_grid_create(rden, rs_desc, error=error)
    CALL rs_grid_set_box ( grid_spme, rs=rden, error=error )
    CALL rs_grid_zero ( rden )

    ALLOCATE ( center ( 3, npart ), STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_center ( particle_set, box, center, npts, n )

    !-------------- DENSITY CALCULATION ----------------
    ipart = 0
    DO
       CALL set_list ( particle_set, npart, center, p1, rden, ipart)
       IF ( p1 == 0 ) EXIT

       ! calculate function on small boxes
       CALL get_patch ( particle_set, box, green, npts, p1, rhos, is_core=.FALSE.,&
            is_shell=.FALSE., unit_charge=.TRUE.)
       rhos = rhos * mcharge(p1)

       ! add boxes to real space grid (big box)
       CALL dg_sum_patch ( rden, rhos, center(:,p1) )
    END DO

    CALL pw_pool_create_pw ( pw_pool, rhob_r, use_data = REALDATA3D,&
         in_space = REALSPACE, error=error)

    CALL rs_pw_transfer ( rden, rhob_r, rs2pw, error=error)

    ! transform density to G space and add charge function
    CALL pw_pool_create_pw ( pw_pool, rhob_g, use_data = COMPLEXDATA1D, &
         in_space = RECIPROCALSPACE ,error=error)
    CALL pw_transfer ( rhob_r, rhob_g, error=error)
    ! update charge function
    rhob_g % cc = rhob_g % cc * green % p3m_charge % cr

    !-------------- ELECTROSTATIC CALCULATION -----------

    ! allocate intermediate arrays
    DO i = 1, 3
       NULLIFY(dphi_g(i)%pw)
       CALL pw_pool_create_pw ( pw_pool, dphi_g ( i )%pw,  &
            use_data = COMPLEXDATA1D,in_space = RECIPROCALSPACE, error=error)
    END DO
    CALL pw_pool_create_pw ( pw_pool, phi_g,  &
         use_data = COMPLEXDATA1D, in_space = RECIPROCALSPACE, error=error)
    IF ( use_virial ) THEN
       CALL pw_poisson_solve ( poisson_env, rhob_g, vgc, phi_g, dphi_g, &
            h_stress=h_stress, error=error)
    ELSE
       CALL pw_poisson_solve ( poisson_env, rhob_g, vgc, phi_g, dphi_g, &
            error=error)
    END IF

    CALL pw_pool_give_back_pw ( pw_pool, rhob_g ,error=error)

    CALL rs_grid_create(rpot, rs_desc, error=error )
    CALL rs_grid_set_box ( grid_spme, rs=rpot, error=error )
    CALL rs_grid_zero ( rpot )
    phi_g%cc = phi_g%cc * green%p3m_charge%cr
    CALL pw_transfer ( phi_g, rhob_r, error=error)
    CALL pw_pool_give_back_pw ( pw_pool, phi_g ,error=error)
    CALL rs_pw_transfer ( rpot, rhob_r, pw2rs, error=error)

    !---------- END OF ELECTROSTATIC CALCULATION --------

    !------------- STRESS TENSOR CALCULATION ------------

    IF ( use_virial ) THEN
       DO i = 1, 3
          DO j = i, 3
             f_stress(i,j) = pw_integral_a2b(dphi_g(i)%pw,dphi_g(j)%pw,error=error)
             f_stress(j,i) = f_stress(i,j)
          END DO
       END DO
       ffa = ( 1.0_dp / fourpi )  * ( 0.5_dp / alpha ) ** 2
       virial%pv_virial = virial%pv_virial - ffa * f_stress + h_stress
    END IF

    !--------END OF STRESS TENSOR CALCULATION -----------

    IF ( calculate_forces ) THEN
       ! move derivative of potential to real space grid and
       ! multiply by charge function in g-space
       ALLOCATE ( drpot(1:3), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, 3
          CALL rs_grid_create(drpot(i)%rs_grid, rs_desc, error=error)
          CALL rs_grid_set_box ( grid_spme, rs=drpot(i)%rs_grid, error=error )
          dphi_g ( i ) % pw % cc = dphi_g ( i ) % pw % cc * green % p3m_charge % cr
          CALL pw_transfer ( dphi_g ( i )%pw, rhob_r, error=error)
          CALL pw_pool_give_back_pw ( pw_pool, dphi_g ( i )%pw ,error=error)
          CALL rs_pw_transfer ( drpot ( i ) % rs_grid, rhob_r, pw2rs, error=error)
       END DO
    ELSE
       DO i = 1, 3
          CALL pw_pool_give_back_pw ( pw_pool, dphi_g ( i )%pw ,error=error)
       END DO
    END IF
    CALL pw_pool_give_back_pw ( pw_pool, rhob_r ,error=error)

    !----------------- FORCE CALCULATION ----------------

    ipart = 0
    DO

       CALL set_list ( particle_set, npart, center, p1, rden, ipart )
       IF ( p1 == 0 ) EXIT

       ! calculate function on small boxes
       CALL get_patch ( particle_set, box, green, npts, p1, rhos, is_core=.FALSE.,&
            is_shell=.FALSE., unit_charge=.TRUE.)

       CALL dg_sum_patch_force_1d ( rpot, rhos, center(:,p1), fint )
       gmcharge(p1,1) = gmcharge(p1,1) + fint * dvols

       IF ( calculate_forces ) THEN
          CALL dg_sum_patch_force_3d ( drpot, rhos, center(:,p1), fat )
          gmcharge(p1,2) = gmcharge(p1,2) - fat(1) * dvols
          gmcharge(p1,3) = gmcharge(p1,3) - fat(2) * dvols
          gmcharge(p1,4) = gmcharge(p1,4) - fat(3) * dvols
       END IF

    END DO

    !--------------END OF FORCE CALCULATION -------------

    !------------------CLEANING UP ----------------------

    CALL rs_grid_release(rden, error=error)
    CALL rs_grid_release(rpot, error=error)
    IF ( calculate_forces ) THEN
       IF (ASSOCIATED(drpot)) THEN
         DO i = 1, 3
           CALL rs_grid_release(drpot(i)%rs_grid, error=error)
         END DO
         DEALLOCATE ( drpot, STAT = stat )
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    DEALLOCATE ( rhos, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( center, STAT = stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE dftb_spme_evaluate

! *****************************************************************************
  SUBROUTINE dftb_ewald_overlap(gmcharge,mcharge,alpha,n_list,&
                                virial,use_virial,error)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: gmcharge
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: mcharge
    REAL(KIND=dp), INTENT(in)                :: alpha
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: n_list
    TYPE(virial_type), POINTER               :: virial
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dftb_ewald_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iab, iatom, ilist, &
                                                inode, jatom, nlist, nmat, &
                                                nnode
    REAL(KIND=dp)                            :: dfr, dr, fr, pfr, rij(3)
    TYPE(neighbor_list_type), POINTER        :: n_list_ab, n_list_local
    TYPE(neighbor_node_type), POINTER        :: n_list_node

    CALL timeset(routineN,handle)

    nmat = SIZE(gmcharge,2)

    DO iab = 1, SIZE(n_list)
      IF (.NOT.ASSOCIATED(n_list(iab)%neighbor_list_set)) CYCLE
      CALL get_neighbor_list_set(neighbor_list_set=&
                                 n_list(iab)%neighbor_list_set,&
                                 nlist=nlist)
      NULLIFY(n_list_ab)

      DO ilist=1,nlist
        IF ( .NOT. ASSOCIATED(n_list_ab) ) THEN
          n_list_ab => first_list(n_list(iab)%neighbor_list_set)
        ELSE
          n_list_ab => next(n_list_ab)
        END IF
        n_list_local => n_list_ab
        CALL get_neighbor_list(neighbor_list=n_list_local,&
                               atom=iatom,nnode=nnode)

        n_list_node => first_node(n_list_local)

        DO inode=1,nnode
          CALL get_neighbor_node(neighbor_node=n_list_node,&
                                 neighbor=jatom,r=rij)

          dr = SQRT(SUM(rij(:)**2))
          IF (dr > 1.e-10) THEN
            fr = erfc(alpha*dr)/dr
            gmcharge(iatom,1)=gmcharge(iatom,1)+mcharge(jatom)*fr
            gmcharge(jatom,1)=gmcharge(jatom,1)+mcharge(iatom)*fr
            IF ( nmat > 1 ) THEN
              dfr = -2._dp*alpha*EXP(-alpha*alpha*dr*dr)*oorootpi/dr-fr/dr
              dfr = -dfr/dr
              DO i=2,nmat
                gmcharge(iatom,i)=gmcharge(iatom,i)-rij(i-1)*mcharge(jatom)*dfr
                gmcharge(jatom,i)=gmcharge(jatom,i)+rij(i-1)*mcharge(iatom)*dfr
              END DO
            END IF
            IF ( use_virial ) THEN
              pfr = -dfr*mcharge(iatom)*mcharge(jatom)
              CALL virial_pair_force ( virial%pv_virial, -pfr, rij, rij, error)
            END IF
          END IF

          n_list_node => next(n_list_node)

        END DO ! inode => jatom(atom B)

      END DO ! ilist => iatom(atom A)

    END DO

    CALL timestop(handle)

  END SUBROUTINE dftb_ewald_overlap

END MODULE qs_dftb_coulomb

